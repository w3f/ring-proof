#![cfg_attr(not(feature = "std"), no_std)]

use ark_ff::{FftField, PrimeField};
use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain, Polynomial};
use ark_poly::univariate::DensePolynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{vec, vec::Vec};
use fflonk::pcs::{Commitment, PCS};

pub mod gadgets;
pub mod test_helpers;
pub mod piop;
pub mod prover;
pub mod verifier;
pub mod transcript;
pub mod domain;

pub trait Column<F: FftField> {
    fn domain(&self) -> GeneralEvaluationDomain<F>;
    fn domain_4x(&self) -> GeneralEvaluationDomain<F>;
    fn as_poly(&self) -> &DensePolynomial<F>;
    fn size(&self) -> usize {
        self.domain().size()
    }
    fn evaluate(&self, z: &F) -> F {
        self.as_poly().evaluate(z)
    }
}

#[derive(Clone)]
pub struct FieldColumn<F: FftField> {
    // actual (constrained) len of the input in evaluation form
    len: usize,
    poly: DensePolynomial<F>,
    evals: Evaluations<F>,
    evals_4x: Evaluations<F>,
}

impl<F: FftField> FieldColumn<F> {
    pub fn shifted_4x(&self) -> Evaluations<F> {
        let mut evals_4x = self.evals_4x.evals.clone();
        evals_4x.rotate_left(4);
        Evaluations::from_vec_and_domain(evals_4x, self.domain_4x())
    }

    pub fn vals(&self) -> &[F] {
        &self.evals.evals[..self.len]
    }
}

impl<F: FftField> Column<F> for FieldColumn<F> {
    fn domain(&self) -> GeneralEvaluationDomain<F> {
        self.evals.domain()
    }

    fn domain_4x(&self) -> GeneralEvaluationDomain<F> {
        self.evals_4x.domain()
    }

    fn as_poly(&self) -> &DensePolynomial<F> {
        &self.poly
    }
}

pub fn const_evals<F: FftField>(c: F, domain: GeneralEvaluationDomain<F>) -> Evaluations<F> {
    Evaluations::from_vec_and_domain(vec![c; domain.size()], domain)
}


pub trait ColumnsEvaluated<F: PrimeField>: CanonicalSerialize + CanonicalDeserialize {
    fn to_vec(self) -> Vec<F>;
}

pub trait ColumnsCommited<F: PrimeField, C: Commitment<F>>: CanonicalSerialize + CanonicalDeserialize {
    fn to_vec(self) -> Vec<C>;
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<F, CS, Commitments, Evaluations>
    where
        F: PrimeField,
        CS: PCS<F>,
        Commitments: ColumnsCommited<F, CS::C>,
        Evaluations: ColumnsEvaluated<F>, {
    pub column_commitments: Commitments,
    pub columns_at_zeta: Evaluations,
    pub quotient_commitment: CS::C,
    pub lin_at_zeta_omega: F,
    pub agg_at_zeta_proof: CS::Proof,
    pub lin_at_zeta_omega_proof: CS::Proof,
}


#[cfg(test)]
mod tests {
    use ark_ed_on_bls12_381_bandersnatch::Fq;
    use ark_ff::{FftField, Field, One, Zero};
    use ark_poly::{Evaluations, GeneralEvaluationDomain, EvaluationDomain, DenseUVPolynomial, Polynomial};
    use ark_poly::univariate::DensePolynomial;
    use ark_std::{test_rng, UniformRand};

    // Divides a polynomial presented in evaluation form over the coset of a domain by the vanishing polynomial of the 2 times smaller domain.
    fn divide_by_vanishing_poly_on_coset<F: FftField>(domain: GeneralEvaluationDomain<F>, evals: &[F]) -> Vec<F> {
        assert_eq!(2 * domain.size(), evals.len());
        // Evaluations of the vanishing polynomial of a domain at the points of the same domain are all the same = g^N - 1.
        let z_coset_eval = domain.evaluate_vanishing_polynomial(F::GENERATOR).inverse().unwrap();
        evals.iter().step_by(2).map(|e| z_coset_eval * e).collect()
    }

    #[test]
    fn divide_by_z_on_coset() {
        let rng = &mut test_rng();

        let log_n = 2;
        let n = 2usize.pow(log_n);
        let domain = GeneralEvaluationDomain::<Fq>::new(n).unwrap();

        let divisor = domain.vanishing_polynomial();
        let quotient = DensePolynomial::<Fq>::rand(n - 1, rng);
        let poly = &quotient * &divisor.into();
        let (q, r) = poly.divide_by_vanishing_poly(domain).unwrap();
        assert_eq!(q, quotient);
        assert!(r.is_zero());

        let domain_2x = GeneralEvaluationDomain::<Fq>::new(2 * n).unwrap();
        let evals_on_coset_2x = domain_2x.coset_fft(&poly.coeffs);

        // domain.divide_by_vanishing_poly_on_coset_in_place(&mut evals_on_coset_2x);
        // let mut quotient_on_coset_2x = domain_2x.coset_fft(&quotient.coeffs);
        // assert_eq!(quotient_on_coset_2x, evals_on_coset_2x);

        // The code above doesn't works when the vanishing polynomial of the smaller domain is evaluated at the points of a larger domain.

        let quotient_evals = divide_by_vanishing_poly_on_coset(domain, &evals_on_coset_2x);
        let quotient_coeffs = domain.coset_ifft(&quotient_evals);
        assert_eq!(DensePolynomial::from_coefficients_vec(quotient_coeffs), quotient);
    }

    #[test]
    fn secret_bits() {
        let rng = &mut test_rng();

        let log_n = 2;
        let n = 2usize.pow(log_n);
        let domain = GeneralEvaluationDomain::<Fq>::new(n).unwrap();
        let domain_x2 = GeneralEvaluationDomain::<Fq>::new(2 * domain.size()).unwrap();
        let blinding_size = 0;
        let input_size = domain.size() - blinding_size;
        let input: Vec<Fq> = (0..input_size).map(|_| if bool::rand(rng) { Fq::one() } else { Fq::zero() }).collect();
        let blinding = (0..blinding_size).map(|_| Fq::rand(rng)).collect();

        let blinded_input = [input, blinding].concat();

        let column = Evaluations::from_vec_and_domain(blinded_input, domain);
        let column_poly = column.interpolate();
        assert_eq!(column_poly.degree(), n - 1);
        let column_on_coset_x2 = domain_x2.coset_fft(&column_poly.coeffs);
        let column_on_coset_x2 = Evaluations::from_vec_and_domain(column_on_coset_x2, domain_x2);

        // constraint
        let ones = vec![Fq::one(); domain_x2.size()];
        let ones = Evaluations::from_vec_and_domain(ones, domain_x2);
        let constraint_evals_on_coset_x2 = &(&ones - &column_on_coset_x2) * &column_on_coset_x2;
        let constraint_coeffs = domain_x2.coset_ifft(&constraint_evals_on_coset_x2.evals);
        let constraint_poly = DensePolynomial::from_coefficients_vec(constraint_coeffs);

        let (expected_quotient, remainder) = constraint_poly.divide_by_vanishing_poly(domain).unwrap();
        assert!(remainder.is_zero());

        let quotient_evals = divide_by_vanishing_poly_on_coset(domain, &constraint_evals_on_coset_x2.evals);
        let quotient_coeffs = domain.coset_ifft(&quotient_evals);
        let quotient_poly = DensePolynomial::from_coefficients_vec(quotient_coeffs);

        assert_eq!(quotient_poly, expected_quotient);
    }
}