use ark_ec::twisted_edwards::{Affine, TECurveConfig};
use ark_ec::AffineRepr;
use ark_ff::PrimeField;
use ark_std::{vec, vec::Vec};
use w3f_pcs::pcs::Commitment;

use w3f_plonk_common::domain::EvaluatedDomain;
use w3f_plonk_common::gadgets::booleanity::BooleanityValues;
use w3f_plonk_common::gadgets::ec::te_doubling::DoublingValues;
use w3f_plonk_common::gadgets::ec::CondAddValues;
use w3f_plonk_common::gadgets::fixed_cells::FixedCellsValues;
use w3f_plonk_common::gadgets::VerifierGadget;
use w3f_plonk_common::piop::VerifierPiop;

use crate::piop::{FixedColumnsCommitted, RingCommitments};
use crate::RingEvaluations;

pub struct PiopVerifier<F: PrimeField, C: Commitment<F>, P: AffineRepr<BaseField = F>> {
    domain_evals: EvaluatedDomain<F>,
    // columns
    fixed_columns_committed: FixedColumnsCommitted<F, C>,
    witness_columns_committed: RingCommitments<F, C>,
    // gadgets
    sk_bits_bool: BooleanityValues<F>,
    pk_from_sk: CondAddValues<F, P>,
    doublings_of_in_gadget: DoublingValues<F, P>,
    out_from_in: CondAddValues<F, P>,
    out_from_in_x: FixedCellsValues<F>,
    out_from_in_y: FixedCellsValues<F>,
}

impl<F: PrimeField, C: Commitment<F>, P: AffineRepr<BaseField = F>> PiopVerifier<F, C, P> {
    pub fn init(
        domain_evals: EvaluatedDomain<F>,
        fixed_columns_committed: FixedColumnsCommitted<F, C>,
        witness_columns_committed: RingCommitments<F, C>,
        all_columns_evaluated: RingEvaluations<F>,
        seed: (F, F),
        _vrf_in: (F, F), //TODO
        vrf_out: (F, F),
    ) -> Self {
        let sk_bits_bool = BooleanityValues {
            bits: all_columns_evaluated.sk_bits,
        };

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

        let doublings_of_in_gadget = DoublingValues {
            doublings: (
                all_columns_evaluated.doublings_of_in[0],
                all_columns_evaluated.doublings_of_in[1],
            ),
            not_last: domain_evals.not_last_row,
            _phantom: Default::default(),
        };
        let out_from_in = CondAddValues {
            bitmask: all_columns_evaluated.sk_bits,
            points: (
                all_columns_evaluated.doublings_of_in[0],
                all_columns_evaluated.doublings_of_in[1],
            ),
            not_last: domain_evals.not_last_row,
            acc: (
                all_columns_evaluated.out_from_in[0],
                all_columns_evaluated.out_from_in[1],
            ),
            _phantom: Default::default(),
        };
        let out_from_in_x = FixedCellsValues {
            col: all_columns_evaluated.out_from_in[0],
            col_first: seed.0,
            col_last: vrf_out.0,
            l_first: domain_evals.l_first,
            l_last: domain_evals.l_last,
        };
        let out_from_in_y = FixedCellsValues {
            col: all_columns_evaluated.out_from_in[1],
            col_first: seed.1,
            col_last: vrf_out.1,
            l_first: domain_evals.l_first,
            l_last: domain_evals.l_last,
        };

        Self {
            domain_evals,
            fixed_columns_committed,
            witness_columns_committed,
            sk_bits_bool,
            pk_from_sk,
            doublings_of_in_gadget,
            out_from_in,
            out_from_in_x,
            out_from_in_y,
        }
    }
}

impl<F: PrimeField, C: Commitment<F>, Jubjub: TECurveConfig<BaseField = F>> VerifierPiop<F, C>
    for PiopVerifier<F, C, Affine<Jubjub>>
{
    const N_CONSTRAINTS: usize = 9;
    const N_COLUMNS: usize = 9;

    fn precommitted_columns(&self) -> Vec<C> {
        self.fixed_columns_committed.as_vec()
    }

    fn evaluate_constraints_main(&self) -> Vec<F> {
        [
            self.sk_bits_bool.evaluate_constraints_main(),
            self.pk_from_sk.evaluate_constraints_main(),
            self.doublings_of_in_gadget.evaluate_constraints_main(),
            self.out_from_in.evaluate_constraints_main(),
            self.out_from_in_x.evaluate_constraints_main(),
            self.out_from_in_y.evaluate_constraints_main(),
        ]
        .concat()
    }

    fn constraint_polynomials_linearized_commitments(&self) -> Vec<C> {
        let pk_x = &self.witness_columns_committed.pk_from_sk[0];
        let pk_y = &self.witness_columns_committed.pk_from_sk[1];
        let (pk_x_coeff, pk_y_coeff) = self.pk_from_sk.acc_coeffs_1();
        let pk_from_sk_c1_lin = pk_x.mul(pk_x_coeff) + pk_y.mul(pk_y_coeff);
        let (pk_x_coeff, pk_y_coeff) = self.pk_from_sk.acc_coeffs_2();
        let pk_from_sk_c2_lin = pk_x.mul(pk_x_coeff) + pk_y.mul(pk_y_coeff);

        let doublings_of_in_x = self.witness_columns_committed.doublings_of_in[0].clone();
        let doublings_of_in_y = self.witness_columns_committed.doublings_of_in[1].clone();
        let doublings_of_in_lin = self
            .doublings_of_in_gadget
            .zeta_omega_poly_commitment(doublings_of_in_x, doublings_of_in_y);

        let out_x = &self.witness_columns_committed.out_from_in[0];
        let out_y = &self.witness_columns_committed.out_from_in[1];
        let (out_x_coeff, out_y_coeff) = self.out_from_in.acc_coeffs_1();
        let out_from_in_c1_lin = out_x.mul(out_x_coeff) + out_y.mul(out_y_coeff);
        let (out_x_coeff, out_y_coeff) = self.out_from_in.acc_coeffs_2();
        let out_from_in_c2_lin = out_x.mul(out_x_coeff) + out_y.mul(out_y_coeff);

        vec![
            pk_from_sk_c1_lin,
            pk_from_sk_c2_lin,
            doublings_of_in_lin[0].clone(),
            doublings_of_in_lin[1].clone(),
            out_from_in_c1_lin,
            out_from_in_c2_lin,
        ]
    }

    fn domain_evaluated(&self) -> &EvaluatedDomain<F> {
        &self.domain_evals
    }
}
