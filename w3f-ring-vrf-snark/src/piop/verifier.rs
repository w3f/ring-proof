use ark_ec::twisted_edwards::{Affine, TECurveConfig};
use ark_ec::AffineRepr;
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_std::marker::PhantomData;
use ark_std::{vec, vec::Vec};
use w3f_pcs::pcs::Commitment;

use w3f_plonk_common::domain::EvaluatedDomain;
use w3f_plonk_common::gadgets::booleanity::BooleanityValues;
use w3f_plonk_common::gadgets::ec::te_doubling::DoublingValues;
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
    // booleanity_of_signer_index: BooleanityValues<F>,
    booleanity_of_secret_key_bits: BooleanityValues<F>,
    pk_from_sk: CondAddValues<F, P>,
    pk_from_sk_x: FixedCellsValues<F>,
    pk_from_sk_y: FixedCellsValues<F>,
    doublings_of_vrf_in: DoublingValues<F, P>,
    // powers_of_in: PowersOfTwoMultipleValuesTE<F, P>, //TODO

    // sole_signer_inner_prod: InnerProdValues<F>,
    // sole_signer_inner_prod_acc: FixedCellsValues<F>,

    // cond_add_pubkey: CondAddValues<F, P>,
    // cond_add_pubkey_acc_x: FixedCellsValues<F>,
    // cond_add_pubkey_acc_y: FixedCellsValues<F>,

    // vrf_out: CondAddValues<F, P>,
    // vrf_out_x: FixedCellsValues<F>,
    // vrf_out_y: FixedCellsValues<F>,
}

impl<F: PrimeField, C: Commitment<F>, P: AffineRepr<BaseField = F>> PiopVerifier<F, C, P> {
    pub fn init(
        domain_evals: EvaluatedDomain<F>,
        fixed_columns_committed: FixedColumnsCommitted<F, C>,
        witness_columns_committed: RingCommitments<F, C>,
        all_columns_evaluated: RingEvaluations<F>,
        seed: (F, F),
        vrf_in: (F, F),
        vrf_out: (F, F),
        signer_pk: (F, F), //TODO REMOVE CRITICAL111
    ) -> Self {
        // C1: `signer_index` column is boolean
        // let booleanity_of_signer_index = BooleanityValues {
        //     bits: all_columns_evaluated.signer_index,
        // };
        // C2 + С3: <signer_index, ring_selector> == 1
        // let sole_signer_inner_prod = InnerProdValues {
        //     a: all_columns_evaluated.ring_selector,
        //     b: all_columns_evaluated.signer_index,
        //     not_last: domain_evals.not_last_row,
        //     acc: all_columns_evaluated.k_is_one_bit,
        // };
        // let inner_prod_acc = FixedCellsValues {
        //     col: all_columns_evaluated.k_is_one_bit,
        //     col_first: F::zero(),
        //     col_last: F::one(),
        //     l_first: domain_evals.l_first,
        //     l_last: domain_evals.l_last,
        // };

        // C4: PK_k := <k, R>
        // let pk_from_k_gadget = CondAddValues {
        //     bitmask: all_columns_evaluated.signer_index,
        //     points: (
        //         all_columns_evaluated.pks[0],
        //         all_columns_evaluated.pks[1],
        //     ),
        //     not_last: domain_evals.not_last_row,
        //     acc: (
        //         all_columns_evaluated.pk_from_k_acc[0],
        //         all_columns_evaluated.pk_from_k_acc[1],
        //     ),
        //     _phantom: PhantomData,
        // };

        // C5: `sk` is boolean
        let booleanity_of_secret_key_bits = BooleanityValues {
            bits: all_columns_evaluated.sk_bits,
        };

        // C6: `PK_sk := sk.G`
        let pk_from_sk = CondAddValues {
            bitmask: all_columns_evaluated.sk_bits,
            points: (
                all_columns_evaluated.doublings_of_g[0],
                all_columns_evaluated.doublings_of_g[1],
            ),
            not_last: domain_evals.not_last_row,
            acc: (
                all_columns_evaluated.pk_from_sk[0],
                all_columns_evaluated.pk_from_sk[1],
            ),
            _phantom: Default::default(),
        };

        // C7 + C8: `PK_sk == PK`
        let pk_from_sk_val_x = FixedCellsValues {
            col: all_columns_evaluated.pk_from_sk[0],
            col_first: seed.0,
            col_last: signer_pk.0,
            l_first: domain_evals.l_first,
            l_last: domain_evals.l_last,
        };
        let pk_from_sk_val_y = FixedCellsValues {
            col: all_columns_evaluated.pk_from_sk[1],
            col_first: seed.1,
            col_last: signer_pk.1,
            l_first: domain_evals.l_first,
            l_last: domain_evals.l_last,
        };

        let doublings_of_vrf_in = DoublingValues {
            doublings: (
                all_columns_evaluated.doublings_of_vrf_in[0],
                all_columns_evaluated.doublings_of_vrf_in[1],
            ),
            not_last: domain_evals.not_last_row,
            _phantom: Default::default(),
        };
        // // C9 + C10: `PK_k == PK`
        // let pk_from_k_val_x = FixedCellsValues {
        //     col: all_columns_evaluated.pk_from_k_acc[0],
        //     col_first: seed.0,
        //     col_last: signer_pk.0,
        //     l_first: domain_evals.l_first,
        //     l_last: domain_evals.l_last,
        // };
        // let pk_from_k_val_y = FixedCellsValues {
        //     col: all_columns_evaluated.pk_from_k_acc[1],
        //     col_first: seed.1,
        //     col_last: signer_pk.1,
        //     l_first: domain_evals.l_first,
        //     l_last: domain_evals.l_last,
        // };

        // // C11: 2-adic powers of the `vrf_input` //TODO
        // let power_of_in = all_columns_evaluated.powers_of_in;
        // // let power_of_in_gadget = PowersOfTwoMultipleValuesTE::init(vrf_in, domain_evals.not_last_row, (power_of_in[0], power_of_in[1]));
        //
        // // C12: `snark_out := <sk, powers_of_in>`
        // let snark_vrf_out = CondAddValues {
        //     bitmask: all_columns_evaluated.signer_sk, // sk
        //     points: (
        //         power_of_in[0],
        //         power_of_in[1],
        //     ),
        //     not_last: domain_evals.not_last_row,
        //     acc: (
        //         all_columns_evaluated.vrf_out_acc[0],
        //         all_columns_evaluated.vrf_out_acc[1],
        //     ),
        //     _phantom: PhantomData,
        // };
        // C13+14: `snark_out == vrf_out`
        // let snark_vrf_out_x = FixedCellsValues {
        //     col: all_columns_evaluated.vrf_out_acc[0],
        //     col_first: seed.0,
        //     col_last: vrf_out.0,
        //     l_first: domain_evals.l_first,
        //     l_last: domain_evals.l_last,
        // };
        // let snark_vrf_out_y = FixedCellsValues {
        //     col: all_columns_evaluated.vrf_out_acc[1],
        //     col_first: seed.1,
        //     col_last: vrf_out.1,
        //     l_first: domain_evals.l_first,
        //     l_last: domain_evals.l_last,
        // };

        Self {
            domain_evals,
            fixed_columns_committed,
            witness_columns_committed,
            // cond_add_pubkey: pk_from_k_gadget,

            // sole_signer_inner_prod,

            // booleanity_of_signer_index,
            booleanity_of_secret_key_bits,

            // cond_add_pubkey_acc_x: pk_from_k_val_x
            // cond_add_pubkey_acc_y: pk_from_k_val_y
            pk_from_sk: pk_from_sk,
            pk_from_sk_x: pk_from_sk_val_x,
            pk_from_sk_y: pk_from_sk_val_y,

            // vrf_out: snark_vrf_out,
            // vrf_out_x: snark_vrf_out_x,
            // vrf_out_y: snark_vrf_out_y,

            // sole_signer_inner_prod_acc: inner_prod_acc,
            // powers_of_in: power_of_in_gadget,
            doublings_of_vrf_in,
        }
    }
}

impl<F: PrimeField, C: Commitment<F>, Jubjub: TECurveConfig<BaseField = F>> VerifierPiop<F, C>
    for PiopVerifier<F, C, Affine<Jubjub>>
{
    const N_CONSTRAINTS: usize = 5;
    const N_COLUMNS: usize = 7;

    fn precommitted_columns(&self) -> Vec<C> {
        self.fixed_columns_committed.as_vec()
    }

    fn evaluate_constraints_main(&self) -> Vec<F> {
        [
            self.booleanity_of_secret_key_bits
                .evaluate_constraints_main(),
            self.pk_from_sk.evaluate_constraints_main(),
            self.doublings_of_vrf_in.evaluate_constraints_main(),
            // self.sole_signer_inner_prod.evaluate_constraints_main(),
            // self.cond_add_pubkey.evaluate_constraints_main(),
            // self.booleanity_of_signer_index.evaluate_constraints_main(),

            // self.cond_add_pubkey_acc_x.evaluate_constraints_main(),
            // self.cond_add_pubkey_acc_y.evaluate_constraints_main(),
            // self.pk_from_sk_x
            //     .evaluate_constraints_main(),
            // self.pk_from_sk_y
            //     .evaluate_constraints_main(),
            // self.vrf_out_x.evaluate_constraints_main(),
            // self.vrf_out_y.evaluate_constraints_main(),
            // self.sole_signer_inner_prod_acc.evaluate_constraints_main(),
        ]
        .concat()
    }

    fn constraint_polynomials_linearized_commitments(&self) -> Vec<C> {
        let acc_x = &self.witness_columns_committed.pk_from_sk[0];
        let acc_y = &self.witness_columns_committed.pk_from_sk[1];

        let (c_acc_x, c_acc_y) = self.pk_from_sk.acc_coeffs_1();
        let c1_lin = acc_x.mul(c_acc_x) + acc_y.mul(c_acc_y);

        let (c_acc_x, c_acc_y) = self.pk_from_sk.acc_coeffs_2();
        let c2_lin = acc_x.mul(c_acc_x) + acc_y.mul(c_acc_y);

        let doublings_of_in_x = self.witness_columns_committed.doublings_of_vrf_in[0].clone();
        let doublings_of_in_y = self.witness_columns_committed.doublings_of_vrf_in[1].clone();
        let c3c4 = self
            .doublings_of_vrf_in
            .zeta_omega_poly_commitment(doublings_of_in_x, doublings_of_in_y);

        vec![c1_lin, c2_lin, c3c4[0].clone(), c3c4[1].clone()]
    }

    fn domain_evaluated(&self) -> &EvaluatedDomain<F> {
        &self.domain_evals
    }
}
