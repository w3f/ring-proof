use ark_ff::{FftField, Field};
use ark_poly::{Evaluations, GeneralEvaluationDomain};
use ark_poly::univariate::DensePolynomial;

pub mod booleanity;
pub mod inner_prod_pub;
pub mod sw_cond_add;
pub mod fixed_cells;

pub trait ProverGadget<F: FftField> {
    // Columns populated by the gadget.
    fn witness_columns(&self) -> Vec<DensePolynomial<F>>;

    // Constraint polynomials.
    fn constraints(&self) -> Vec<Evaluations<F>>;

    // 'Linearized' parts of the constraint polynomials.
    fn constraints_linearized(&self, zeta: &F) -> Vec<DensePolynomial<F>>;

    // Subgroup over which the columns are defined.
    fn domain(&self) -> GeneralEvaluationDomain<F>;
}

pub trait VerifierGadget<F: Field> {
    fn evaluate_constraints_main(&self) -> Vec<F>;
}

#[cfg(test)]
mod tests {
    use ark_ff::Zero;
    use ark_poly::EvaluationDomain;
    use ark_std::test_rng;

    use super::*;

    pub fn test_gadget<F: FftField, G: ProverGadget<F>>(gadget: G) {
        let domain = gadget.domain();
        let cs_polys: Vec<DensePolynomial<F>> = gadget.constraints().into_iter()
            .map(|c| c.interpolate())
            .collect();

        let mut rs = cs_polys.iter()
            .map(|c| c.divide_by_vanishing_poly(domain).unwrap().1);
        assert!(rs.all(|r| r.is_zero()));

        let _zeta = F::rand(&mut test_rng());
        let _omega = domain.group_gen();

        // let evaluated = gadget.evaluate_assignment(&zeta);
        // let cs_evaluated_main = evaluated.evaluate_constraints_main();
        // assert_eq!(cs_evaluated_main.len(), cs_polys.len());
        // let constraints_linearized = gadget.constraints_linearized(&evaluated);
        // assert_eq!(constraints_linearized.len(), cs_polys.len());
        //
        // for ((c_full, c_ev_main), c_lin) in cs_polys.iter()
        //     .zip(cs_evaluated_main)
        //     .zip(constraints_linearized) {
        //     assert_eq!(c_full.evaluate(&zeta), c_ev_main + c_lin.evaluate(&(zeta * omega)));
        // }
    }
}