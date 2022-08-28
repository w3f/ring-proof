use ark_ec::AffineCurve;
use ark_ff::PrimeField;
use fflonk::pcs::Commitment;
use common::gadgets::booleanity::BooleanityValues;
use common::gadgets::fixed_cells::FixedCellsValues;
use common::gadgets::inner_prod_pub::InnerProdValues;
use common::gadgets::sw_cond_add::{CondAdd, CondAddValues};
use common::gadgets::VerifierGadget;
use common::piop::VerifierPiop;
use crate::piop::{RingCommitments, RingEvaluations, SelectorsValues};

pub struct PiopVerifier<F: PrimeField, C: Commitment<F>> {
    selectors: SelectorsValues<F>,
    points: [C; 2],
    columns: RingCommitments<F, C>,
    evals: RingEvaluations<F>,
    inner_prod: InnerProdValues<F>,
    cond_add: CondAddValues<F>,
    booleanity: BooleanityValues<F>,
    fixed_cells_acc_x: FixedCellsValues<F>,
    fixed_cells_acc_y: FixedCellsValues<F>,
}

impl<F: PrimeField, C: Commitment<F>> PiopVerifier<F, C> {
    pub fn init(points: &[C; 2],
                columns: RingCommitments<F, C>,
                evals: RingEvaluations<F>,
                selectors: SelectorsValues<F>,
                init: (F, F),
                result: (F, F),
    ) -> Self {

        let cond_add = CondAddValues {
            bitmask: evals.bits,
            points: (evals.points[0], evals.points[1]),
            not_last: selectors.not_last,
            acc: (evals.cond_add_acc[0], evals.cond_add_acc[1]),
        };

        let inner_prod = InnerProdValues {
            a: selectors.ring_selector,
            b: evals.bits,
            l_last: selectors.l_last, //TODO: can be done in O(1)
            acc: evals.inn_prod_acc,
            inner_prod: F::one(),
        };

        let booleanity = BooleanityValues {
            bits: evals.bits,
        };

        let fixed_cells_acc_x = FixedCellsValues {
            col: evals.cond_add_acc[0],
            col_first: init.0,
            col_last: result.0,
            l_first: selectors.l_first,
            l_last: selectors.l_last,
        };

        let fixed_cells_acc_y = FixedCellsValues {
            col: evals.cond_add_acc[1],
            col_first: init.1,
            col_last: result.1,
            l_first: selectors.l_first,
            l_last: selectors.l_last,
        };

        Self {
            selectors,
            points: points.clone(),
            columns,
            evals,
            inner_prod,
            cond_add,
            booleanity,
            fixed_cells_acc_x,
            fixed_cells_acc_y,
        }
    }
}

impl<F: PrimeField, C: Commitment<F>> VerifierPiop<F, C> for PiopVerifier<F, C> {
    const N_CONSTRAINTS: usize = 6;
    const N_COLUMNS: usize = 6;

    fn precommitted_columns(&self) -> Vec<C> {
        self.points.to_vec()
    }

    fn evaluate_constraints_main(&self) -> Vec<F> {
        vec![
            self.inner_prod.evaluate_constraints_main(),
            self.cond_add.evaluate_constraints_main(),
            self.booleanity.evaluate_constraints_main(),
            self.fixed_cells_acc_x.evaluate_constraints_main(),
            self.fixed_cells_acc_y.evaluate_constraints_main(),
        ].concat()
    }

    fn constraint_polynomials_linearized_commitments(&self) -> Vec<C> {
        let inner_prod_acc = self.columns.inn_prod_acc.clone();
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

    fn get_n(&self) -> (usize, F) {
        (self.selectors.n, self.selectors.omega)
    }
}
