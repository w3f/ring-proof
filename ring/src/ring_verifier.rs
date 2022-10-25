use std::marker::PhantomData;
use ark_ec::CurveGroup;
use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::PrimeField;
use fflonk::pcs::{PCS, PcsParams};
use common::domain::{Domain, EvaluatedDomain};

use common::gadgets::sw_cond_add::CondAdd;
use common::piop::VerifierPiop;
use common::verifier::PlonkVerifier;
use crate::piop::params::PiopParams;

use crate::piop::{FixedColumnsCommitted, PiopVerifier};
use crate::RingProof;

pub struct RingVerifier<F: PrimeField, CS: PCS<F>, Curve: SWCurveConfig<BaseField=F>> {
    domain: Domain<F>,
    fixed_cols: FixedColumnsCommitted<F, CS::C>,
    plonk_verifier: PlonkVerifier<F, CS, merlin::Transcript>,
    phantom: PhantomData<Curve>,

}

impl<F: PrimeField, CS: PCS<F>, Curve: SWCurveConfig<BaseField=F>> RingVerifier<F, CS, Curve> {
    pub fn init(raw_vk: &<CS::Params as PcsParams>::RVK,
                domain: Domain<F>,
                fixed_cols: FixedColumnsCommitted<F, CS::C>,
                empty_transcript: merlin::Transcript,
    ) -> Self {
        let plonk_verifier = PlonkVerifier::init(raw_vk, domain.domain(), fixed_cols.points_committed.clone(), empty_transcript);

        Self {
            domain,
            fixed_cols,
            plonk_verifier,
            phantom: PhantomData, //TODO

        }
    }

    pub fn verify_ring_proof(&self, proof: RingProof<F, CS>, result: Affine<Curve>) -> bool {
        let challenges = self.plonk_verifier.restore_challenges(
            &result,
            &proof,
            // '1' accounts for the quotient polynomial that is aggregated together with the columns
            PiopVerifier::<F, CS::C>::N_COLUMNS + 1,
            PiopVerifier::<F, CS::C>::N_CONSTRAINTS,
        );
        let init = CondAdd::<F, Affine<Curve>>::point_in_g1_complement();
        let init_plus_result = (init + result).into_affine();
        let domain_eval = EvaluatedDomain::new(self.domain.domain(), challenges.zeta, self.domain.hiding);

        let piop = PiopVerifier::init::<Curve>(
            domain_eval,
            &self.fixed_cols,
            proof.column_commitments.clone(),
            proof.columns_at_zeta.clone(),
            (init.x, init.y),
            (init_plus_result.x, init_plus_result.y),
        );
        self.plonk_verifier.verify(piop, proof, challenges)
    }
}

