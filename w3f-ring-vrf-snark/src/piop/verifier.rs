use ark_ec::twisted_edwards::{Affine, TECurveConfig};
use ark_ec::AffineRepr;
use ark_ff::PrimeField;
use ark_std::marker::PhantomData;
use ark_std::{vec, vec::Vec};
use w3f_pcs::pcs::Commitment;

use w3f_plonk_common::domain::EvaluatedDomain;
use w3f_plonk_common::gadgets::booleanity::BooleanityValues;
use w3f_plonk_common::gadgets::ec::CondAddValues;
use w3f_plonk_common::gadgets::fixed_cells::FixedCellsValues;
use w3f_plonk_common::gadgets::inner_prod::InnerProdValues;
use w3f_plonk_common::gadgets::VerifierGadget;
use w3f_plonk_common::piop::VerifierPiop;

use crate::piop::{FixedColumnsCommitted, RingCommitments};
use crate::RingEvaluations;

pub struct PiopVerifier<F: PrimeField, C: Commitment<F>, P: AffineRepr<BaseField = F>> {
    domain_evals: EvaluatedDomain<F>,
    fixed_columns_committed: FixedColumnsCommitted<F, C>,
    witness_columns_committed: RingCommitments<F, C>,
    // Gadget verifiers:
    booleanity_of_signer_index: BooleanityValues<F>,
    booleanity_of_secret_key_bits: BooleanityValues<F>,

    sole_signer_inner_prod: InnerProdValues<F>,
    sole_signer_inner_prod_acc: FixedCellsValues<F>,

    cond_add_pubkey: CondAddValues<F, P>,
    cond_add_pubkey_acc_x: FixedCellsValues<F>,
    cond_add_pubkey_acc_y: FixedCellsValues<F>,

    //cond_add_gen_multiples: CondAddValuesT,
    cond_add_gen_multiples_acc_x: FixedCellsValues<F>,
    cond_add_gen_multiples_acc_y: FixedCellsValues<F>,

    //cond_add_vrfout: CondAddValuesT,
    cond_add_vrfout_acc_x: FixedCellsValues<F>,
    cond_add_vrfout_acc_y: FixedCellsValues<F>,
}

impl<F: PrimeField, C: Commitment<F>, P: AffineRepr<BaseField = F>> PiopVerifier<F, C, P> {
    pub fn init(
        domain_evals: EvaluatedDomain<F>,
        fixed_columns_committed: FixedColumnsCommitted<F, C>,
        witness_columns_committed: RingCommitments<F, C>,
        all_columns_evaluated: RingEvaluations<F>,
        init: (F, F),
        result: (F, F),
    ) -> Self {
        let cond_add_pubkey = CondAddValues {
            bitmask: all_columns_evaluated.signer_index,
            points: (
                all_columns_evaluated.pubkey_points[0],
                all_columns_evaluated.pubkey_points[1],
            ),
            not_last: domain_evals.not_last_row,
            acc: (
                all_columns_evaluated.cond_add_pubkey_acc[0],
                all_columns_evaluated.cond_add_pubkey_acc[1],
            ),
            _phantom: PhantomData,
        };

        // TODO: These all need to be added
        // let cond_add_gen_multiples = CondAddValuesT::init(
        //     all_columns_evaluated.signer_secret_key_bits,
        //     (
        //         all_columns_evaluated.gen_multiples_points[0],
        //         all_columns_evaluated.gen_multiples_points[1],
        //     ),
        //     domain_evals.not_last_row,
        //     (
        //         all_columns_evaluated.cond_add_gen_multiples_acc[0],
        //         all_columns_evaluated.cond_add_gen_multiples_acc[1],
        //     ),
        // );

        // let cond_add_vrfout = CondAddValuesT::init(
        //     all_columns_evaluated.signer_secret_key_bits,
        //     (
        //         all_columns_evaluated.vrfout_points[0],
        //         all_columns_evaluated.vrfout_points[1],
        //     ),
        //     domain_evals.not_last_row,
        //     (
        //         all_columns_evaluated.cond_add_vrfout_acc[0],
        //         all_columns_evaluated.cond_add_vrfout_acc[1],
        //     ),
        // );

        let sole_signer_inner_prod = InnerProdValues {
            a: all_columns_evaluated.ring_selector,
            b: all_columns_evaluated.signer_index,
            not_last: domain_evals.not_last_row,
            acc: all_columns_evaluated.sole_signer_inn_prod_acc,
        };

        let booleanity_of_signer_index = BooleanityValues {
            bits: all_columns_evaluated.signer_index,
        };

        let booleanity_of_secret_key_bits = BooleanityValues {
            bits: all_columns_evaluated.signer_secret_key_bits,
        };

        let cond_add_pubkey_acc_x = FixedCellsValues {
            col: all_columns_evaluated.cond_add_pubkey_acc[0],
            col_first: init.0,
            col_last: result.0,
            l_first: domain_evals.l_first,
            l_last: domain_evals.l_last,
        };

        let cond_add_pubkey_acc_y = FixedCellsValues {
            col: all_columns_evaluated.cond_add_pubkey_acc[1],
            col_first: init.1,
            col_last: result.1,
            l_first: domain_evals.l_first,
            l_last: domain_evals.l_last,
        };

        let cond_add_gen_multiples_acc_x = FixedCellsValues {
            col: all_columns_evaluated.cond_add_gen_multiples_acc[0],
            col_first: init.0,
            col_last: result.0,
            l_first: domain_evals.l_first,
            l_last: domain_evals.l_last,
        };

        let cond_add_gen_multiples_acc_y = FixedCellsValues {
            col: all_columns_evaluated.cond_add_gen_multiples_acc[1],
            col_first: init.1,
            col_last: result.1,
            l_first: domain_evals.l_first,
            l_last: domain_evals.l_last,
        };

        let cond_add_vrfout_acc_x = FixedCellsValues {
            col: all_columns_evaluated.cond_add_vrfout_acc[0],
            col_first: init.0,
            col_last: result.0,
            l_first: domain_evals.l_first,
            l_last: domain_evals.l_last,
        };

        let cond_add_vrfout_acc_y = FixedCellsValues {
            col: all_columns_evaluated.cond_add_vrfout_acc[1],
            col_first: init.1,
            col_last: result.1,
            l_first: domain_evals.l_first,
            l_last: domain_evals.l_last,
        };

        let sole_signer_inner_prod_acc = FixedCellsValues {
            col: all_columns_evaluated.sole_signer_inn_prod_acc,
            col_first: F::zero(),
            col_last: F::one(),
            l_first: domain_evals.l_first,
            l_last: domain_evals.l_last,
        };

        Self {
            domain_evals,
            fixed_columns_committed,
            witness_columns_committed,
            cond_add_pubkey,

            sole_signer_inner_prod,

            booleanity_of_signer_index,
            booleanity_of_secret_key_bits,

            cond_add_pubkey_acc_x,
            cond_add_pubkey_acc_y,

            cond_add_gen_multiples_acc_x,
            cond_add_gen_multiples_acc_y,

            cond_add_vrfout_acc_x,
            cond_add_vrfout_acc_y,

            sole_signer_inner_prod_acc,
        }
    }
}

impl<F: PrimeField, C: Commitment<F>, Jubjub: TECurveConfig<BaseField = F>> VerifierPiop<F, C>
    for PiopVerifier<F, C, Affine<Jubjub>>
{
    const N_CONSTRAINTS: usize = 7;
    const N_COLUMNS: usize = 7;

    fn precommitted_columns(&self) -> Vec<C> {
        self.fixed_columns_committed.as_vec()
    }

    fn evaluate_constraints_main(&self) -> Vec<F> {
        [
            self.sole_signer_inner_prod.evaluate_constraints_main(),
            self.cond_add_pubkey.evaluate_constraints_main(),
            self.booleanity_of_signer_index.evaluate_constraints_main(),
            self.booleanity_of_secret_key_bits
                .evaluate_constraints_main(),
            self.cond_add_pubkey_acc_x.evaluate_constraints_main(),
            self.cond_add_pubkey_acc_y.evaluate_constraints_main(),
            self.cond_add_gen_multiples_acc_x
                .evaluate_constraints_main(),
            self.cond_add_gen_multiples_acc_y
                .evaluate_constraints_main(),
            self.cond_add_vrfout_acc_x.evaluate_constraints_main(),
            self.cond_add_vrfout_acc_y.evaluate_constraints_main(),
            self.sole_signer_inner_prod_acc.evaluate_constraints_main(),
        ]
        .concat()
    }

    fn constraint_polynomials_linearized_commitments(&self) -> Vec<C> {
        let inner_prod_acc = self
            .witness_columns_committed
            .sole_signer_inn_prod_acc
            .mul(self.sole_signer_inner_prod.not_last);
        let pubkey_acc_x = &self.witness_columns_committed.cond_add_pubkey_acc[0];
        let pubkey_acc_y = &self.witness_columns_committed.cond_add_pubkey_acc[1];

        let (c_acc_x, c_acc_y) = self.cond_add_pubkey.acc_coeffs_1();
        let c1_lin = pubkey_acc_x.mul(c_acc_x) + pubkey_acc_y.mul(c_acc_y);

        let (c_acc_x, c_acc_y) = self.cond_add_pubkey.acc_coeffs_2();
        let c2_lin = pubkey_acc_x.mul(c_acc_x) + pubkey_acc_y.mul(c_acc_y);

        vec![inner_prod_acc, c1_lin, c2_lin]
    }

    fn domain_evaluated(&self) -> &EvaluatedDomain<F> {
        &self.domain_evals
    }
}
