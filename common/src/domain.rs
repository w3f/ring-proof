use ark_ff::{FftField, Zero};
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Evaluations, GeneralEvaluationDomain, Polynomial};
use ark_poly::univariate::DensePolynomial;
use crate::FieldColumn;

// Domains for performing calculations with constraint polynomials of degree up to 4.
struct Domains<F: FftField> {
    x1: GeneralEvaluationDomain<F>,
    x4: GeneralEvaluationDomain<F>,
}

impl<F: FftField> Domains<F> {
    fn new(n: usize) -> Self {
        let x1 = GeneralEvaluationDomain::<F>::new(n).unwrap_or_else(|| panic!("No domain of size {}", n));
        let x4 = GeneralEvaluationDomain::<F>::new(4 * n).unwrap_or_else(|| panic!("No domain of size {}", 4 * n));
        Self { x1, x4 }
    }

    fn column_from_evals(&self, evals: Vec<F>, len: usize) -> FieldColumn<F> {
        assert_eq!(evals.len(), self.x1.size());
        let evals = Evaluations::from_vec_and_domain(evals, self.x1);
        let poly = evals.interpolate_by_ref();
        let evals_4x = poly.evaluate_over_domain_by_ref(self.x4);
        FieldColumn { len, poly, evals, evals_4x }
    }

    fn column_from_poly(&self, poly: DensePolynomial<F>, len: usize) -> FieldColumn<F> {
        assert!(poly.degree() < self.x1.size());
        let evals_4x = self.amplify(&poly);
        let evals = evals_4x.evals.iter().step_by(4).cloned().collect();
        let evals = Evaluations::from_vec_and_domain(evals, self.x1);
        FieldColumn { len, poly, evals, evals_4x }
    }

    // Amplifies the number of the evaluations of the polynomial so it can be multiplied in linear time.
    fn amplify(&self, poly: &DensePolynomial<F>) -> Evaluations<F> {
        poly.evaluate_over_domain_by_ref(self.x4)
    }
}

pub struct Domain<F: FftField> {
    domains: Domains<F>,
    pub capacity: usize,
    pub not_last_row: FieldColumn<F>,
    pub l_first: FieldColumn<F>,
    pub l_last: FieldColumn<F>,
}

impl<F: FftField> Domain<F> {
    pub fn new(n: usize) -> Self {
        let domains = Domains::new(n);
        let capacity = domains.x1.size();
        let not_last_row = not_last(domains.x1);
        let not_last_row = domains.column_from_poly(not_last_row, capacity);
        let l_first = domains.column_from_evals(l_i(0, capacity), capacity);
        let l_last = domains.column_from_evals(l_i(capacity - 1, capacity), capacity);
        Self {
            domains,
            capacity,
            not_last_row,
            l_first,
            l_last,
        }
    }

    pub(crate) fn divide_by_vanishing_poly<>(
        &self,
        poly: &DensePolynomial<F>,
    ) -> DensePolynomial<F> {
        let (quotient, remainder) = poly.divide_by_vanishing_poly(self.domains.x1).unwrap(); //TODO error-handling
        assert!(remainder.is_zero()); //TODO error-handling
        quotient
    }

    pub fn column(&self, mut evals: Vec<F>) -> FieldColumn<F> {
        let len = evals.len();
        assert!(len <= self.capacity);
        evals.resize(self.capacity, F::zero());
        self.domains.column_from_evals(evals, len)
    }
}

fn l_i<F: FftField>(i: usize, n: usize) -> Vec<F> {
    let mut l_i = vec![F::zero(); n];
    l_i[i] = F::one();
    l_i
}

pub fn not_last<F: FftField>(domain: GeneralEvaluationDomain<F>) -> DensePolynomial<F> {
    let x = &DensePolynomial::from_coefficients_slice(&[F::zero(), F::one()]);
    let w_last = domain.group_gen().pow(&[domain.size() as u64 - 1]);
    let w_last = &DensePolynomial::from_coefficients_slice(&[w_last]);
    x - w_last
}