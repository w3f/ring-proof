use ark_ff::{Field, PrimeField};
use ark_serialize::CanonicalSerialize;
use ark_std::{vec, vec::Vec};
use ark_std::rand::Rng;
use fflonk::pcs::{Commitment, PCS, PcsParams};
use rand_chacha::ChaCha20Rng;

use crate::{ColumnsCommited, ColumnsEvaluated, Proof};
use crate::piop::VerifierPiop;
use crate::transcript::Transcript;

pub struct PlonkVerifier<F: PrimeField, CS: PCS<F>, T: Transcript<F, CS>> {
    // Polynomial commitment scheme verifier's key.
    pcs_vk: CS::VK,
    // Transcript,
    // initialized with the public parameters and the commitments to the precommitted columns.
    transcript_prelude: T,
}

impl<F: PrimeField, CS: PCS<F>, T: Transcript<F, CS>> PlonkVerifier<F, CS, T> {
    pub fn init(pcs_vk: <CS::Params as PcsParams>::VK,
                verifier_key: &impl CanonicalSerialize,
                empty_transcript: T) -> Self {
        let mut transcript_prelude = empty_transcript;
        transcript_prelude._add_serializable(b"vk", verifier_key);

        Self {
            pcs_vk,
            transcript_prelude,
        }
    }

    //~ #### Plonk.Verify\
    //~ **Inputs**:\
    //~   - $Piop$: an object of Piop type\
    //~   - $Proof$: a proof tuple as defined in ???\
    //~   - $Challenges: ([\alpha_1,...,\alpha_n, \zeta, [\nu_1,..,nu_n])$ A Plonk Verifier challenge defined in ???\
    //~   - $H$: R Random oracle\
    //~ **Output**:\
    //~   - A boolean value indicating if the $Proof$ represents a correct proof\
    //~
    //~ *****  
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
        //~ $C \leftarrow EvaluateConsttrain(Poip)$  
        //~ $E \leftarrow \sum_i^n \alpha[i] * C[i]$  
        //~ $D \leftarrow DomainEvaluated(Piop)$  
        let eval: F = piop.evaluate_constraints_main().iter().zip(challenges.alphas.iter()).map(|(c, alpha)| *alpha * c).sum();
        let zeta = challenges.zeta;
        let domain_evaluated = piop.domain_evaluated();

        //~ $q_{zeta} \leftarrow \frac{D}{\Omega(C + Proof.x)}$  
        let q_zeta = domain_evaluated.divide_by_vanishing_poly_in_zeta(eval + proof.lin_at_zeta_omega);

        let mut columns = [
            piop.precommitted_columns(),
            proof.column_commitments.to_vec(),
        ].concat();
        columns.push(proof.quotient_commitment.clone());

        let mut columns_at_zeta = proof.columns_at_zeta.to_vec();
        columns_at_zeta.push(q_zeta);

        let cl = CS::C::combine(&challenges.nus, &columns);
        let agg_y = columns_at_zeta.into_iter().zip(challenges.nus.iter()).map(|(y, r)| y * r).sum();

        let lin_pices = piop.constraint_polynomials_linearized_commitments();
        let lin_comm = CS::C::combine(&challenges.alphas[..3], &lin_pices);

        let zeta_omega = zeta * domain_evaluated.omega();
        
        //~ **return** $BatchVerify()$  
        CS::batch_verify(&self.pcs_vk, vec![cl, lin_comm], vec![challenges.zeta, zeta_omega], vec![agg_y, proof.lin_at_zeta_omega], vec![proof.agg_at_zeta_proof, proof.lin_at_zeta_omega_proof], rng)
        //~      
        //~ *****    

    }

    pub fn restore_challenges<Commitments, Evaluations>(
        &self,
        instance: &impl CanonicalSerialize,
        proof: &Proof<F, CS, Commitments, Evaluations>,
        n_polys: usize,
        n_constraints: usize,
    ) -> (Challenges<F>, ChaCha20Rng)
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
        let challenges = Challenges {
            alphas,
            zeta,
            nus,
        };
        (challenges, transcript.to_rng())
    }
}


//~ ### Challenge  
//~
//~ **Definition**: *Plonk verifier challange* is defined as triple:
//~ $$([\alpha_1,...,\alpha_n, \zeta, [\nu_1,..,nu_n])$$
//~ where \alpha_i, zeta and nu_i are all elements of Plonk Scalar Field.
//~
pub struct Challenges<F: Field> {
    pub alphas: Vec<F>,
    pub zeta: F,
    pub nus: Vec<F>,
}

