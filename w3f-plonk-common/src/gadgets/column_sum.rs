use ark_ff::{FftField, Field};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{Evaluations, GeneralEvaluationDomain};
use ark_std::rc::Rc;
use ark_std::{vec, vec::Vec};

use crate::domain::Domain;
use crate::gadgets::{ProverGadget, VerifierGadget};
use crate::{Column, FieldColumn};

pub struct ColumnSumPolys<F: FftField> {
    /// Input column.
    /// Should have length `n-1`, where `n` is the domain "capacity" (domain.size - ZK_ROWS):
    /// `col[0], ..., col[n-2]`
    pub col: Rc<FieldColumn<F>>,
    /// Partial sums of `col`: `acc[0] = 0, acc[i] = col[0] + ... + col[i-1], i = 1,...,n-1`
    pub acc: Rc<FieldColumn<F>>,
    pub not_last: FieldColumn<F>,
}

pub struct ColumnSumEvals<F: Field> {
    pub col: F,
    pub acc: F,
    pub not_last: F,
}

impl<F: FftField> ColumnSumPolys<F> {
    pub fn init(col: Rc<FieldColumn<F>>, domain: &Domain<F>) -> Self {
        assert_eq!(col.len, domain.capacity - 1); // last element is not constrained
        let partial_sums = Self::partial_sums(col.vals());
        let mut acc = vec![F::zero()];
        acc.extend(partial_sums);
        let acc = domain.private_column(acc);
        let acc = Rc::new(acc);
        Self {
            col,
            acc,
            not_last: domain.not_last_row.clone(),
        }
    }

    /// Returns `col[0], col[0] + col[1], ..., col[0] + col[1] + ... + col[n-1]`.
    fn partial_sums(col: &[F]) -> Vec<F> {
        col.iter()
            .scan(F::zero(), |state, &x| {
                *state += x;
                Some(*state)
            })
            .collect()
    }
}

impl<F: FftField> ProverGadget<F> for ColumnSumPolys<F> {
    fn witness_columns(&self) -> Vec<DensePolynomial<F>> {
        vec![self.acc.poly.clone()]
    }

    fn constraints(&self) -> Vec<Evaluations<F>> {
        let col = &self.col.evals_4x;
        let acc = &self.acc.evals_4x;
        let acc_shifted = &self.acc.shifted_4x();
        let not_last = &self.not_last.evals_4x;
        let c = &(&(acc_shifted - acc) - col) * not_last;
        vec![c]
    }

    fn constraints_linearized(&self, z: &F) -> Vec<DensePolynomial<F>> {
        let c = &self.acc.poly * self.not_last.evaluate(z);
        vec![c]
    }

    fn domain(&self) -> GeneralEvaluationDomain<F> {
        self.col.evals.domain()
    }
}

impl<F: Field> VerifierGadget<F> for ColumnSumEvals<F> {
    fn evaluate_constraints_main(&self) -> Vec<F> {
        let c = (-self.acc - self.col) * self.not_last;
        vec![c]
    }
}

#[cfg(test)]
mod tests {
    use ark_ed_on_bls12_381_bandersnatch::Fq;
    use ark_ff::Zero;
    use ark_poly::Polynomial;
    use ark_std::test_rng;

    use crate::domain::Domain;
    use crate::test_helpers::random_vec;

    use super::*;

    fn _test_column_sum_gadget(hiding: bool) {
        let rng = &mut test_rng();

        let log_n = 10;
        let n = 2usize.pow(log_n);
        let domain = Domain::new(n, hiding);

        let col = random_vec(domain.capacity - 1, rng);
        let sum = col.iter().sum();
        let col = Rc::new(domain.private_column(col));

        let gadget = ColumnSumPolys::<Fq>::init(col, &domain);

        let acc = &gadget.acc.evals.evals;
        assert!(acc[0].is_zero());
        assert_eq!(acc[domain.capacity - 1], sum);

        let constraint_poly = gadget.constraints()[0].interpolate_by_ref();

        assert_eq!(constraint_poly.degree(), n);
        domain.divide_by_vanishing_poly(&constraint_poly);
    }

    #[test]
    fn test_column_sum_gadget() {
        _test_column_sum_gadget(false);
        _test_column_sum_gadget(true);
    }
}
