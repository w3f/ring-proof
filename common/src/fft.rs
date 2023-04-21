use ark_ff::{Field, PrimeField, Zero};
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_poly::univariate::DensePolynomial;
use num_bigint::BigUint;
use num_integer::Integer;

trait FftDomain<F: Field> {
    fn fft(&self, coeffs: &[F]) -> Vec<F>;
    fn n(&self) -> usize;
    fn w(&self) -> F;
}

struct NaiveDomain<F: Field> {
    n: usize,
    w: F,
}

impl<F: Field> NaiveDomain<F> {
    fn new(w: F, n: usize) -> Self {
        assert!(w.pow([n as u64]).is_one());
        Self {
            n,
            w,
        }
    }
}

impl<F: Field> FftDomain<F> for NaiveDomain<F> {
    fn fft(&self, coeffs: &[F]) -> Vec<F> {
        assert_eq!(coeffs.len(), self.n);
        let poly = DensePolynomial::from_coefficients_slice(coeffs);
        let mut res = vec![F::zero(); self.n];
        res[0] = poly.evaluate(&F::one());
        res[1] = poly.evaluate(&self.w);
        let mut wi = self.w;
        for i in 2..self.n {
            wi *= self.w;
            res[i] = poly.evaluate(&wi);
        }
        res
    }

    fn n(&self) -> usize {
        self.n
    }

    fn w(&self) -> F {
        self.w
    }
}

struct CooleyTukeyDomain<F: Field, D: FftDomain<F>> {
    n: usize,
    w: F,
    d1: D,
    d2: D,
    twiddles: Vec<Vec<F>>,
}

impl<F: Field, D: FftDomain<F>> CooleyTukeyDomain<F, D> {
    fn new(w: F, d1: D, d2: D) -> Self {
        let n = d1.n() * d2.n();
        assert!(w.pow([n as u64]).is_one());
        let mut twiddles = vec![vec![F::zero(); d1.n()]; d2.n()];
        for i2 in 0..d2.n() {
            let mut inner = vec![F::zero(); d1.n()];
            for k1 in 0..d1.n() {
                inner[k1] = w.pow([(k1 * i2) as u64]);
            }
            twiddles[i2] = inner;
        }

        Self {
            n,
            w,
            d1,
            d2,
            twiddles,
        }
    }
}

impl<F: Field, D: FftDomain<F>> FftDomain<F> for CooleyTukeyDomain<F, D> {
    fn fft(&self, coeffs: &[F]) -> Vec<F> {
        let n = self.n;
        let n1 = self.d1.n();
        let n2 = self.d2.n();
        assert_eq!(coeffs.len(), n);

        let mut inner_dfts = vec![vec![F::zero(); n2]; n1];
        for k1 in 0..n1 {
            let inner_coeffs: Vec<_> = coeffs.iter().cloned().skip(k1).step_by(n1).collect();
            assert_eq!(inner_coeffs.len(), n2);
            inner_dfts[k1] = self.d2.fft(&inner_coeffs);
        }

        let mut res = vec![F::zero(); n];
        for i2 in 0..n2 {
            let mut outer_coeffs = vec![F::zero(); n1];
            for k1 in 0..n1 {
                outer_coeffs[k1] = self.twiddles[i2][k1] * inner_dfts[k1][i2];
            }
            let outer_dft = self.d1.fft(&outer_coeffs);
            for i1 in 0..n1 {
                res[i1 * n2 + i2] = outer_dft[i1];
            }
        }
        res
    }

    fn n(&self) -> usize {
        self.n
    }

    fn w(&self) -> F {
        self.w
    }
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
    use ark_std::{end_timer, start_timer, test_rng, UniformRand};

    use super::*;

    #[test]
    fn test_cooley_tukey() {
        use ark_bls12_381::Fq;
        let rng = &mut test_rng();

        let n1 = 23;
        let n2 = 47;
        let n = n1 * n2;

        let w = gen::<Fq>(n);
        let w1 = w.pow([n2 as u64]);
        let w2 = w.pow([n1 as u64]);
        let d1 = NaiveDomain::new(w1, n1);
        let d2 = NaiveDomain::new(w2, n2);
        let ctd = CooleyTukeyDomain::new(w, d1, d2);
        let nd = NaiveDomain::new(w, n);

        let coeffs: Vec<_> = (0..n).map(|_| Fq::rand(rng)).collect();

        let t_fft = start_timer!(|| "fft");
        let fft = ctd.fft(&coeffs);
        end_timer!(t_fft);

        let t_dft = start_timer!(|| "dft");
        let dft = nd.fft(&coeffs);
        end_timer!(t_dft);
        assert_eq!(fft, dft);
    }
}
