use ark_ff::PrimeField;
use ark_poly::GeneralEvaluationDomain;
use ark_serialize::CanonicalSerialize;
use ark_std::vec::Vec;
use w3f_pcs::pcs::{PcsParams, PCS};
use rand_core::RngCore;

use crate::{ColumnsCommited, ColumnsEvaluated};

pub trait PlonkTranscript<F: PrimeField, CS: PCS<F>>: Clone {
    fn add_protocol_params(
        &mut self,
        domain: &GeneralEvaluationDomain<F>,
        pcs_raw_vk: &<CS::Params as PcsParams>::RVK,
    ) {
        self._add_serializable(b"domain", domain);
        self._add_serializable(b"pcs_raw_vk", pcs_raw_vk);
    }

    fn add_precommitted_cols(&mut self, precommitted_cols: &[CS::C; 2]) {
        self._add_serializable(b"precommitted_cols", precommitted_cols);
    }

    fn add_instance(&mut self, instance: &impl CanonicalSerialize) {
        self._add_serializable(b"instance", instance);
    }

    fn add_committed_cols(&mut self, committed_cols: &impl ColumnsCommited<F, CS::C>) {
        self._add_serializable(b"committed_cols", committed_cols);
    }

    fn get_constraints_aggregation_coeffs(&mut self, n: usize) -> Vec<F> {
        self._128_bit_coeffs(b"constraints_aggregation", n)
    }

    fn add_quotient_commitment(&mut self, point: &CS::C) {
        self._add_serializable(b"quotient", point);
    }

    fn add_kzg_proofs(&mut self, in_zeta: &CS::Proof, in_zeta_omega: &CS::Proof) {
        self._add_serializable(b"kzg_proof_zeta", in_zeta);
        self._add_serializable(b"kzg_proof_zeta_omega", in_zeta_omega);
    }

    fn get_evaluation_point(&mut self) -> F {
        self._128_bit_point(b"evaluation_point")
    }

    fn add_evaluations(&mut self, evals: &impl ColumnsEvaluated<F>, r_at_zeta_omega: &F) {
        self._add_serializable(b"register_evaluations", evals);
        self._add_serializable(b"shifted_linearization_evaluation", r_at_zeta_omega);
    }

    fn get_kzg_aggregation_challenges(&mut self, n: usize) -> Vec<F> {
        self._128_bit_coeffs(b"kzg_aggregation", n)
    }

    fn _128_bit_point(&mut self, label: &'static [u8]) -> F;

    fn _128_bit_coeffs(&mut self, label: &'static [u8], n: usize) -> Vec<F> {
        (0..n).map(|_| self._128_bit_point(label)).collect()
    }

    fn _add_serializable(&mut self, label: &'static [u8], message: &impl CanonicalSerialize);

    fn to_rng(self) -> impl RngCore;
}
