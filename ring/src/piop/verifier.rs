use ark_ec::short_weierstrass::SWCurveConfig;
use ark_ff::PrimeField;
use ark_poly::{Evaluations, Polynomial};
use fflonk::pcs::Commitment;
use common::domain::EvaluatedDomain;
use common::gadgets::booleanity::BooleanityValues;
use common::gadgets::fixed_cells::FixedCellsValues;
use common::gadgets::inner_prod::InnerProdValues;
use common::gadgets::sw_cond_add::CondAddValues;
use common::gadgets::VerifierGadget;
use common::piop::VerifierPiop;
use crate::piop::{RingCommitments, RingEvaluations};
use crate::piop::params::PiopParams;

pub struct PiopVerifier<F: PrimeField, C: Commitment<F>> {
    domain_evals: EvaluatedDomain<F>,
    points: [C; 2],
    columns: RingCommitments<F, C>,
    evals: RingEvaluations<F>,
    inner_prod: InnerProdValues<F>,
    cond_add: CondAddValues<F>,
    booleanity: BooleanityValues<F>,
    cond_add_acc_x: FixedCellsValues<F>,
    cond_add_acc_y: FixedCellsValues<F>,
    inner_prod_acc: FixedCellsValues<F>,
}

impl<F: PrimeField, C: Commitment<F>> PiopVerifier<F, C> {
    pub fn init<Curve: SWCurveConfig<BaseField=F>>(
        piop_params: &PiopParams<F, Curve>,
        domain_evals: EvaluatedDomain<F>,
        points: &[C; 2],
        columns: RingCommitments<F, C>,
        evals: RingEvaluations<F>,
        init: (F, F),
        result: (F, F),
        zeta: F,
    ) -> Self {
        let keyset_part_selector = piop_params.keyset_part_selector();
        let keyset_part_selector = Evaluations::from_vec_and_domain(keyset_part_selector, domain_evals.domain);
        let keyset_part_selector_at_zeta = keyset_part_selector.interpolate().evaluate(&zeta);

        let cond_add = CondAddValues {
            bitmask: evals.bits,
            points: (evals.points[0], evals.points[1]),
            not_last: domain_evals.not_last_row,
            acc: (evals.cond_add_acc[0], evals.cond_add_acc[1]),
        };

        let inner_prod = InnerProdValues {
            a: keyset_part_selector_at_zeta,
            b: evals.bits,
            not_last: domain_evals.not_last_row,
            acc: evals.inn_prod_acc,
        };

        let booleanity = BooleanityValues {
            bits: evals.bits,
        };

        let cond_add_acc_x = FixedCellsValues {
            col: evals.cond_add_acc[0],
            col_first: init.0,
            col_last: result.0,
            l_first: domain_evals.l_first,
            l_last: domain_evals.l_last,
        };

        let cond_add_acc_y = FixedCellsValues {
            col: evals.cond_add_acc[1],
            col_first: init.1,
            col_last: result.1,
            l_first: domain_evals.l_first,
            l_last: domain_evals.l_last,
        };

        let inner_prod_acc = FixedCellsValues {
            col: evals.inn_prod_acc,
            col_first: F::zero(),
            col_last: F::one(),
            l_first: domain_evals.l_first,
            l_last: domain_evals.l_last,
        };

        Self {
            domain_evals,
            points: points.clone(),
            columns,
            evals,
            inner_prod,
            cond_add,
            booleanity,
            cond_add_acc_x,
            cond_add_acc_y,
            inner_prod_acc,
        }
    }
}

impl<F: PrimeField, C: Commitment<F>> VerifierPiop<F, C> for PiopVerifier<F, C> {
    const N_CONSTRAINTS: usize = 7;
    const N_COLUMNS: usize = 6;

    fn precommitted_columns(&self) -> Vec<C> {
        self.points.to_vec()
    }

    fn evaluate_constraints_main(&self) -> Vec<F> {
        vec![
            self.inner_prod.evaluate_constraints_main(),
            self.cond_add.evaluate_constraints_main(),
            self.booleanity.evaluate_constraints_main(),
            self.cond_add_acc_x.evaluate_constraints_main(),
            self.cond_add_acc_y.evaluate_constraints_main(),
            self.inner_prod_acc.evaluate_constraints_main(),
        ].concat()
    }

    fn constraint_polynomials_linearized_commitments(&self) -> Vec<C> {
        let inner_prod_acc = self.columns.inn_prod_acc.mul(self.inner_prod.not_last);
        let acc_x = &self.columns.cond_add_acc[0];
        let acc_y = &self.columns.cond_add_acc[1];
        // let (c_acc_x_1, c_acc_y_1) = self.cond_add.acc_coeffs_1();
        // let (c_acc_x_2, c_acc_y_2) = self.cond_add.acc_coeffs_2();
        // vec![
        //     vec![(inner_prod_acc, F::one())],
        //     vec![(cond_add_acc_x.clone(), c_acc_x_1), (cond_add_acc_y.clone(), c_acc_y_1)],
        //     vec![(cond_add_acc_x, c_acc_x_2), (cond_add_acc_y, c_acc_y_2)],
        // ]

        let (c_acc_x, c_acc_y) = self.cond_add.acc_coeffs_1();
        let c1_lin = acc_x.mul(c_acc_x) + acc_y.mul(c_acc_y);

        let (c_acc_x, c_acc_y) = self.cond_add.acc_coeffs_2();
        let c2_lin = acc_x.mul(c_acc_x) + acc_y.mul(c_acc_y);

        vec![
            inner_prod_acc,
            c1_lin,
            c2_lin,
        ]
    }

    fn domain_evaluated(&self) -> &EvaluatedDomain<F> {
        &self.domain_evals
    }
}
