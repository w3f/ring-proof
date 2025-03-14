use ark_ff::{Field, PrimeField};
use ark_serialize::CanonicalSerialize;
use ark_std::rand::Rng;
use ark_std::{vec, vec::Vec};
use rand_core::RngCore;
use w3f_pcs::pcs::{Commitment, PcsParams, PCS};

use crate::piop::VerifierPiop;
use crate::transcript::PlonkTranscript;
use crate::{ColumnsCommited, ColumnsEvaluated, Proof};

pub struct PlonkVerifier<F: PrimeField, CS: PCS<F>, T: PlonkTranscript<F, CS>> {
    // Polynomial commitment scheme verifier's key.
    pcs_vk: CS::VK,
    // Transcript,
    // initialized with the public parameters and the commitments to the precommitted columns.
    transcript_prelude: T,
}

impl<F: PrimeField, CS: PCS<F>, T: PlonkTranscript<F, CS>> PlonkVerifier<F, CS, T> {
    pub fn init(
        pcs_vk: <CS::Params as PcsParams>::VK,
        verifier_key: &impl CanonicalSerialize,
        empty_transcript: T,
    ) -> Self {
        let mut transcript_prelude = empty_transcript;
        transcript_prelude._add_serializable(b"vk", verifier_key);

        Self {
            pcs_vk,
            transcript_prelude,
        }
    }

    pub fn verify<Piop, Commitments, Evaluations, R: Rng>(
        &self,
        piop: Piop,
        proof: Proof<F, CS, Commitments, Evaluations>,
        challenges: Challenges<F>,
        rng: &mut R,
    ) -> bool
    where
        Piop: VerifierPiop<F, CS::C>,
        Commitments: ColumnsCommited<F, CS::C>,
        Evaluations: ColumnsEvaluated<F>,
    {
        let q_zeta = piop.evaluate_q_at_zeta(&challenges.alphas, proof.lin_at_zeta_omega);

        let mut columns = [
            piop.precommitted_columns(),
            proof.column_commitments.to_vec(),
        ]
        .concat();
        columns.push(proof.quotient_commitment.clone());

        let mut columns_at_zeta = proof.columns_at_zeta.to_vec();
        columns_at_zeta.push(q_zeta);

        let agg_comm = CS::C::combine(&challenges.nus, &columns);
        let agg_at_zeta = columns_at_zeta
            .into_iter()
            .zip(challenges.nus.iter())
            .map(|(y, r)| y * r)
            .sum();

        let lin_comm = piop.lin_poly_commitment(&challenges.alphas);

        let zeta = challenges.zeta;
        let zeta_omega = zeta * piop.domain_evaluated().omega();

        CS::batch_verify(
            &self.pcs_vk,
            vec![agg_comm, lin_comm],
            vec![zeta, zeta_omega],
            vec![agg_at_zeta, proof.lin_at_zeta_omega],
            vec![proof.agg_at_zeta_proof, proof.lin_at_zeta_omega_proof],
            rng,
        )
        .is_ok()
    }

    pub fn restore_challenges<Commitments, Evaluations>(
        &self,
        instance: &impl CanonicalSerialize,
        proof: &Proof<F, CS, Commitments, Evaluations>,
        n_polys: usize,
        n_constraints: usize,
    ) -> (Challenges<F>, impl RngCore)
    where
        Commitments: ColumnsCommited<F, CS::C>,
        Evaluations: ColumnsEvaluated<F>,
    {
        let mut transcript = self.transcript_prelude.clone();
        transcript.add_instance(instance);
        transcript.add_committed_cols(&proof.column_commitments);
        // let r = transcript.get_bitmask_aggregation_challenge();
        // transcript.append_2nd_round_register_commitments(&proof.additional_commitments);
        let alphas = transcript.get_constraints_aggregation_coeffs(n_constraints);
        transcript.add_quotient_commitment(&proof.quotient_commitment);
        let zeta = transcript.get_evaluation_point();
        transcript.add_evaluations(&proof.columns_at_zeta, &proof.lin_at_zeta_omega);
        let nus = transcript.get_kzg_aggregation_challenges(n_polys);
        transcript.add_kzg_proofs(&proof.agg_at_zeta_proof, &proof.lin_at_zeta_omega_proof);
        let challenges = Challenges { alphas, zeta, nus };
        (challenges, transcript.to_rng())
    }
}

pub struct Challenges<F: Field> {
    pub alphas: Vec<F>,
    pub zeta: F,
    pub nus: Vec<F>,
}
