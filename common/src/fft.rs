use ark_ff::{Field, PrimeField, Zero};
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_poly::univariate::DensePolynomial;
use num_bigint::BigUint;
use num_integer::Integer;

struct CooleyTukeyDomain<F: PrimeField> {
    n: usize,
    n1: usize,
    n2: usize,
    g1: F,
    g2: F,
    twiddles: Vec<Vec<F>>,
}

impl<F: PrimeField> CooleyTukeyDomain<F> {
    fn new(g: F, n1: usize, n2: usize) -> Self {
        let n = n1 * n2;
        assert!(g.pow([n as u64]).is_one());
        let g2 = g.pow([n1 as u64]);
        let g1 = g.pow([n2 as u64]);
        assert!(g1.pow([n1 as u64]).is_one());
        assert!(g2.pow([n2 as u64]).is_one());

        let mut twiddles = vec![vec![F::zero(); n1]; n2];
        for i2 in 0..n2 {
            let mut inner = vec![F::zero(); n1];
            for k1 in 0..n1 {
                inner[k1] = g.pow([(k1 * i2) as u64]);
            }
            twiddles[i2] = inner;
        }

        Self {
            n,
            n1,
            n2,
            g1,
            g2,
            twiddles,
        }
    }

    fn fft(&self, coeffs: &[F]) -> Vec<F> {
        let n = self.n;
        let n1 = self.n1;
        let n2 = self.n2;
        let wn1 = self.g2;
        let wn2 = self.g1;
        assert_eq!(coeffs.len(), n);

        let mut inner_dfts = vec![vec![F::zero(); n2]; n1];
        for k1 in 0..n1 {
            let inner_coeffs: Vec<_> = coeffs.iter().cloned().skip(k1).step_by(n1).collect();
            assert_eq!(inner_coeffs.len(), n2);
            inner_dfts[k1] = dft(&wn1, &inner_coeffs);
        }

        let mut res = vec![F::zero(); n];
        for i2 in 0..n2 {
            let mut outer_coeffs = vec![F::zero(); n1];
            for k1 in 0..n1 {
                outer_coeffs[k1] = self.twiddles[i2][k1] * inner_dfts[k1][i2];
            }
            let outer_dft = dft(&wn2, &outer_coeffs);
            for i1 in 0..n1 {
                res[i1 * n2 + i2] = outer_dft[i1];
            }
        }
        res
    }
}

fn dft<F: Field>(w: &F, coeffs: &[F]) -> Vec<F> {
    let n = coeffs.len();
    let poly = DensePolynomial::from_coefficients_slice(coeffs);
    let mut res = vec![F::zero(); n];
    res[0] = poly.evaluate(&F::one());
    res[1] = poly.evaluate(w);
    let mut wi = *w;
    for i in 2..n {
        wi *= w;
        res[i] = poly.evaluate(&wi);
    }
    res
}

fn gen<F: PrimeField>(n: usize) -> F {
    let field_size: BigUint = F::MODULUS.into();
    let multiplicative_group_order = field_size - 1u8;
    let (cofactor, rem) = multiplicative_group_order.div_rem(&n.into());
    assert!(rem.is_zero());
    let cofactor: F::BigInt = cofactor.try_into().unwrap();
    let gen = F::GENERATOR.pow(cofactor);
    gen
}

#[cfg(test)]
mod tests {
    use ark_std::{test_rng, UniformRand};

    use super::*;

    #[test]
    fn test_cooley_tukey() {
        use ark_bls12_381::Fq;
        let rng = &mut test_rng();

        let n1 = 23;
        let n2 = 47;
        let n = n1 * n2;

        let g = gen::<Fq>(n);

        let domain = CooleyTukeyDomain::new(g, n1, n2);
        let x: Vec<_> = (0..n).map(|_| Fq::rand(rng)).collect();
        let fft = domain.fft(&x);
        let dft = dft(&g, &x);
        assert_eq!(fft, dft);
    }
}
