use ark_ec::AffineRepr;
use ark_ec::CurveGroup;
use ark_ff::PrimeField;
use common::gadgets::cond_add::CondAddValuesFor;
use common::gadgets::cond_add::{AffineCondAdd, CondAdd};
use fflonk::pcs::{RawVerifierKey, PCS};

use common::domain::EvaluatedDomain;
use common::piop::VerifierPiop;
use common::transcript::PlonkTranscript;
use common::verifier::PlonkVerifier;

use crate::piop::{params::PiopParams, FixedColumnsCommitted, PiopVerifier, VerifierKey};
use crate::ArkTranscript;
use crate::RingProof;

pub struct RingVerifier<
    F: PrimeField,
    CS: PCS<F>,
    P: AffineRepr<BaseField = F>,
    T: PlonkTranscript<F, CS> = ArkTranscript,
> {
    piop_params: PiopParams<F, P>,
    fixed_columns_committed: FixedColumnsCommitted<F, CS::C>,
    plonk_verifier: PlonkVerifier<F, CS, T>,
}

impl<F: PrimeField, CS: PCS<F>, P: AffineRepr<BaseField = F>, T: PlonkTranscript<F, CS>>
    RingVerifier<F, CS, P, T>
where
    P: AffineCondAdd,
{
    pub fn init(
        verifier_key: VerifierKey<F, CS>,
        piop_params: PiopParams<F, P>,
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

    pub fn verify(&self, proof: RingProof<F, CS>, result: P) -> bool {
        let (challenges, mut rng) = self.plonk_verifier.restore_challenges(
            &result,
            &proof,
            // '1' accounts for the quotient polynomial that is aggregated together with the columns
            PiopVerifier::<F, CS::C, <P::CondAddT as CondAdd<F, P>>::Values>::N_COLUMNS + 1,
            PiopVerifier::<F, CS::C, <P::CondAddT as CondAdd<F, P>>::Values>::N_CONSTRAINTS,
        );
        let seed = self.piop_params.seed;
        let seed_plus_result = (seed + result).into_affine();
        let domain_eval = EvaluatedDomain::new(
            self.piop_params.domain.domain(),
            challenges.zeta,
            self.piop_params.domain.hiding,
        );

        let piop: PiopVerifier<F, <CS as PCS<F>>::C, CondAddValuesFor<P>> = PiopVerifier::init(
            domain_eval,
            self.fixed_columns_committed.clone(),
            proof.column_commitments.clone(),
            proof.columns_at_zeta.clone(),
            (*seed.x().unwrap(), *seed.y().unwrap()),
            (
                *seed_plus_result.x().unwrap(),
                *seed_plus_result.y().unwrap(),
            ),
        );

        self.plonk_verifier
            .verify(piop, proof, challenges, &mut rng)
    }

    pub fn piop_params(&self) -> &PiopParams<F, P> {
        &self.piop_params
    }
}
