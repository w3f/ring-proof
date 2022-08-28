use ark_ff::{FftField, Zero};
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Evaluations, GeneralEvaluationDomain, Polynomial};
use ark_poly::univariate::DensePolynomial;
use ark_std::test_rng;
use crate::FieldColumn;

const ZK_ROWS: usize = 3;

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
    hiding: bool,
    pub capacity: usize,
    pub not_last_row: FieldColumn<F>,
    pub l_first: FieldColumn<F>,
    pub l_last: FieldColumn<F>,
    zk_rows_vanishing_poly: Option<DensePolynomial<F>>,
}

impl<F: FftField> Domain<F> {
    pub fn new(n: usize, hiding: bool) -> Self {
        let domains = Domains::new(n);
        let size = domains.x1.size();
        let capacity = if hiding { size - ZK_ROWS } else { size };
        let last_row_index = capacity - 1;

        let l_first = l_i(0, size);
        let l_first = domains.column_from_evals(l_first, capacity);
        let l_last = l_i(last_row_index, size);
        let l_last = domains.column_from_evals(l_last, capacity);
        let not_last_row = vanishes_on_row(last_row_index, domains.x1);
        let not_last_row = domains.column_from_poly(not_last_row, capacity);

        let zk_rows_vanishing_poly = hiding.then(|| vanishes_on_last_3_rows(domains.x1));

        Self {
            domains,
            hiding,
            capacity,
            not_last_row,
            l_first,
            l_last,
            zk_rows_vanishing_poly,
        }
    }

    pub(crate) fn divide_by_vanishing_poly<>(
        &self,
        poly: &DensePolynomial<F>,
    ) -> DensePolynomial<F> {
        let (quotient, remainder) = if self.hiding {
            let exclude_zk_rows = poly * self.zk_rows_vanishing_poly.as_ref().unwrap();
            exclude_zk_rows.divide_by_vanishing_poly(self.domains.x1).unwrap() //TODO error-handling
        } else {
            poly.divide_by_vanishing_poly(self.domains.x1).unwrap() //TODO error-handling
        };
        assert!(remainder.is_zero()); //TODO error-handling
        quotient
    }

    pub fn column(&self, mut evals: Vec<F>) -> FieldColumn<F> {
        let len = evals.len();
        assert!(len <= self.capacity);
        evals.resize(self.capacity, F::zero());
        if self.hiding {
            evals.resize_with(self.domains.x1.size(), || F::rand(&mut test_rng())); //TODO
        }
        self.domains.column_from_evals(evals, len)
    }

    // public column
    pub fn selector(&self, mut evals: Vec<F>) -> FieldColumn<F> {
        let len = evals.len();
        assert!(len <= self.capacity);
        evals.resize(self.domains.x1.size(), F::zero());
        self.domains.column_from_evals(evals, len)
    }

    pub fn omega(&self) -> F {
        self.domains.x1.group_gen()
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

// (x - w^i)
pub fn vanishes_on_row<F: FftField>(i: usize, domain: GeneralEvaluationDomain<F>) -> DensePolynomial<F> {
    assert!(i < domain.size());
    let w = domain.group_gen();
    let wi = w.pow(&[i as u64]);
    let wi = DensePolynomial::from_coefficients_slice(&[wi]);
    let x = DensePolynomial::from_coefficients_slice(&[F::zero(), F::one()]);
    &x - &wi
}

/// The polynomial
/// (x - w^{n - 3}) * (x - w^{n - 2}) * (x - w^{n - 1})
pub fn vanishes_on_last_3_rows<F: FftField>(domain: GeneralEvaluationDomain<F>) -> DensePolynomial<F> {
    let w = domain.group_gen();
    let n3 = (domain.size() - ZK_ROWS) as u64;
    let w3 = w.pow(&[n3]);
    let w2 = w3 * w;
    let w1 = w2 * w;
    assert_eq!(w1, domain.group_gen_inv());
    let x = DensePolynomial::from_coefficients_slice(&[F::zero(), F::one()]); // X
    let c = |a: F| DensePolynomial::from_coefficients_slice(&[a]);
    &(&(&x - &c(w3)) * &(&x - &c(w2))) * &(&x - &c(w1))
}
