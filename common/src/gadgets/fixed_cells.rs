use ark_ff::{FftField, Field, Zero};
use ark_poly::Evaluations;
use ark_poly::univariate::DensePolynomial;

use crate::{Column, const_evals, FieldColumn};
use crate::gadgets::VerifierGadget;

pub struct FixedCells<F: FftField> {
    col: FieldColumn<F>,
    l_first: FieldColumn<F>,
    l_last: FieldColumn<F>,
}


pub struct FixedCellsValues<F: Field> {
    pub col: F,
    pub col_first: F,
    pub col_last: F,
    pub l_first: F,
    pub l_last: F,
}


impl<F: FftField> FixedCells<F> {
    pub fn init(col: FieldColumn<F>, l_first: FieldColumn<F>, l_last: FieldColumn<F>) -> Self {
        let n = col.size();
        assert_eq!(l_first.size(), n);
        assert_eq!(l_last.size(), n);
        Self { col, l_first, l_last }
    }

    pub fn constraints(&self) -> Vec<Evaluations<F>> {
        let col = &self.col;
        let domain = col.domain_4x();
        let first = &const_evals(col.evals.evals[0], domain);
        let last = &const_evals(col.evals.evals[col.size() - 1], domain);
        let col = &self.col.evals_4x;
        let l_first = &self.l_first.evals_4x;
        let l_last = &self.l_last.evals_4x;
        let c = &(l_first * &(col - first)) + &(l_last * &(col - last));
        vec![c]
    }

    pub fn constraints_linearized(&self, _z: &F) -> Vec<DensePolynomial<F>> {
        vec![DensePolynomial::zero()]
    }
}


impl<F: Field> VerifierGadget<F> for FixedCellsValues<F> {
    fn evaluate_constraints_main(&self) -> Vec<F> {
        let c = (self.col - self.col_first) * self.l_first + (self.col - self.col_last) * self.l_last;
        vec![c]
    }
}