use ark_ec::twisted_edwards::{Affine, TECurveConfig};
use ark_ec::CurveGroup;
use ark_ff::PrimeField;
use fflonk::pcs::{RawVerifierKey, PCS};

use common::domain::EvaluatedDomain;
use common::piop::VerifierPiop;
use common::transcript::PlonkTranscript;
use common::verifier::PlonkVerifier;

use crate::piop::params::PiopParams;
use crate::piop::{FixedColumnsCommitted, PiopVerifier, VerifierKey};
use crate::RingProof;

pub struct RingVerifier<F, CS, Jubjub, T>
where
    F: PrimeField,
    CS: PCS<F>,
    Jubjub: TECurveConfig<BaseField = F>,
    T: PlonkTranscript<F, CS>,
{
    piop_params: PiopParams<F, Jubjub>,
    fixed_columns_committed: FixedColumnsCommitted<F, CS::C>,
    plonk_verifier: PlonkVerifier<F, CS, T>,
}

impl<F, CS, Jubjub, T> RingVerifier<F, CS, Jubjub, T>
where
    F: PrimeField,
    CS: PCS<F>,
    Jubjub: TECurveConfig<BaseField = F>,
    T: PlonkTranscript<F, CS>,
{
    pub fn init(
        verifier_key: VerifierKey<F, CS>,
        piop_params: PiopParams<F, Jubjub>,
        empty_transcript: T,
    ) -> Self {
        let pcs_vk = verifier_key.pcs_raw_vk.prepare();
        let plonk_verifier = PlonkVerifier::init(pcs_vk, &verifier_key, empty_transcript);
        Self {
            piop_params,
            fixed_columns_committed: verifier_key.fixed_columns_committed,
            plonk_verifier,
        }
    }

    pub fn verify_ring_proof(&self, proof: RingProof<F, CS>, result: Affine<Jubjub>) -> bool {
        let (challenges, mut rng) = self.plonk_verifier.restore_challenges(
            &result,
            &proof,
            // '1' accounts for the quotient polynomial that is aggregated together with the columns
            PiopVerifier::<F, CS::C, Affine<Jubjub>>::N_COLUMNS + 1,
            PiopVerifier::<F, CS::C, Affine<Jubjub>>::N_CONSTRAINTS,
        );
        let seed = self.piop_params.seed;
        let seed_plus_result = (seed + result).into_affine();
        let domain_eval = EvaluatedDomain::new(
            self.piop_params.domain.domain(),
            challenges.zeta,
            self.piop_params.domain.hiding,
        );

        let piop = PiopVerifier::<_, _, Affine<Jubjub>>::init(
            domain_eval,
            self.fixed_columns_committed.clone(),
            proof.column_commitments.clone(),
            proof.columns_at_zeta.clone(),
            (seed.x, seed.y),
            (seed_plus_result.x, seed_plus_result.y),
        );

        self.plonk_verifier
            .verify(piop, proof, challenges, &mut rng)
    }

    pub fn piop_params(&self) -> &PiopParams<F, Jubjub> {
        &self.piop_params
    }
}
