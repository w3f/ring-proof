use ark_ff::{FftField, Field, PrimeField};
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Evaluations, GeneralEvaluationDomain as Domain, Polynomial};
use ark_poly::univariate::DensePolynomial;
use fflonk::pcs::{Commitment, PCS};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_serialize::{Read, SerializationError, Write};


pub mod gadgets;
pub mod test_helpers;
pub mod piop;
pub mod prover;
pub mod verifier;
pub mod transcript;
pub mod setup;

pub trait Column<F: FftField> {
    fn domain(&self) -> Domain<F>;
    fn domain_4x(&self) -> Domain<F>;
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
    poly: DensePolynomial<F>,
    evals: Evaluations<F>,
    evals_4x: Evaluations<F>,
}

impl<F: FftField> FieldColumn<F> {
    pub fn init(vals: Vec<F>) -> Self {
        let n = vals.len();
        let domain = Domain::new(n)
            .expect("field is not smooth enough to construct domain");
        let evals = Evaluations::from_vec_and_domain(vals, domain);
        let poly = evals.interpolate_by_ref();
        let evals_4x = Self::amplify(&poly, n);
        Self { poly, evals, evals_4x }
    }

    pub fn from_poly(poly: DensePolynomial<F>, n: usize) -> Self {
        assert!(poly.degree() < n);
        let evals_4x = Self::amplify(&poly, n);
        let domain = Domain::new(n)
            .expect("field is not smooth enough to construct domain");
        let evals = evals_4x.evals.iter().step_by(4).cloned().collect();
        let evals = Evaluations::from_vec_and_domain(evals, domain);
        assert_eq!(poly, evals.interpolate_by_ref());
        Self { poly, evals, evals_4x }
    }

    // Amplifies the number of the evaluations of the polynomial so it can be multiplied in linear time.
    fn amplify(poly: &DensePolynomial<F>, n: usize) -> Evaluations<F, Domain<F>> {
        let domain_4x = Domain::new(4 * n)
            .expect("field is not smooth enough to construct domain");
        let evals_4x = poly.evaluate_over_domain_by_ref(domain_4x);
        evals_4x
    }

    pub fn shifted_4x(&self) -> Evaluations<F> {
        let mut evals_4x = self.evals_4x.evals.clone();
        evals_4x.rotate_left(4);
        Evaluations::from_vec_and_domain(evals_4x, self.domain_4x())
    }
}

impl<F: FftField> Column<F> for FieldColumn<F> {
    fn domain(&self) -> Domain<F> {
        self.evals.domain()
    }

    fn domain_4x(&self) -> Domain<F> {
        self.evals_4x.domain()
    }

    fn as_poly(&self) -> &DensePolynomial<F> {
        &self.poly
    }
}

pub fn const_evals<F: FftField>(c: F, domain: Domain<F>) -> Evaluations<F> {
    Evaluations::from_vec_and_domain(vec![c; domain.size()], domain)
}

fn l_i<F: FftField>(i: usize, n: usize) -> Vec<F> {
    let mut l_i = vec![F::zero(); n];
    l_i[i] = F::one();
    l_i
}

pub fn l_first<F: FftField>(n: usize) -> Vec<F> {
    l_i(0, n)
}

pub fn l_last<F: FftField>(n: usize) -> Vec<F> {
    l_i(n - 1, n)
}

pub fn not_last<F: FftField>(domain: Domain<F>) -> DensePolynomial<F> {
    let x = &DensePolynomial::from_coefficients_slice(&[F::zero(), F::one()]);
    let w_last = domain.group_gen().pow(&[domain.size() as u64 - 1]);
    let w_last = &DensePolynomial::from_coefficients_slice(&[w_last]);
    x - w_last
}

pub trait ColumnsEvaluated<F: PrimeField>: CanonicalSerialize + CanonicalDeserialize {
    fn to_vec(self) -> Vec<F>;
}

pub trait ColumnsCommited<F: PrimeField, C: Commitment<F>>: CanonicalSerialize + CanonicalDeserialize {
    fn to_vec(self) -> Vec<C>;
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
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
