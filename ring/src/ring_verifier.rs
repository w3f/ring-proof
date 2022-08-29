use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::PrimeField;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use fflonk::pcs::{PCS, PcsParams};
use common::domain::{Domain, EvaluatedDomain};

use common::gadgets::sw_cond_add::CondAdd;
use common::piop::VerifierPiop;
use common::verifier::PlonkVerifier;
use crate::piop::params::PiopParams;

use crate::piop::PiopVerifier;
use crate::piop::SelectorColumns;
use crate::RingProof;

pub struct RingVerifier<F: PrimeField, CS: PCS<F>, Curve: SWCurveConfig<BaseField=F>> {
    piop_params: PiopParams<F, Curve>,
    points_comm: [CS::C; 2],
    domain_size: usize,
    keyset_size: usize,

    plonk_verifier: PlonkVerifier<F, CS, merlin::Transcript>,
}

impl<F: PrimeField, CS: PCS<F>, Curve: SWCurveConfig<BaseField=F>> RingVerifier<F, CS, Curve> {
    pub fn init(raw_vk: &<CS::Params as PcsParams>::RVK,
                piop_params: PiopParams<F, Curve>,
                points_comm: [CS::C; 2],
                domain_size: usize,
                keyset_size: usize,
                empty_transcript: merlin::Transcript,
    ) -> Self {
        let domain = piop_params.domain.domain();
        let plonk_verifier = PlonkVerifier::init(raw_vk, domain, points_comm.clone(), empty_transcript);

        Self {
            piop_params,
            points_comm,
            domain_size,
            keyset_size,
            plonk_verifier,
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
        let init_plus_result = init + result;
        let domain_eval = EvaluatedDomain::new(self.piop_params.domain.domain(), challenges.zeta, self.piop_params.domain.hiding);

        let piop = PiopVerifier::init(
            &self.piop_params,
            domain_eval,
            &self.points_comm,
            proof.column_commitments.clone(),
            proof.columns_at_zeta.clone(),
            (init.x, init.y),
            (init_plus_result.x, init_plus_result.y),
            challenges.zeta,
        );
        self.plonk_verifier.verify(piop, proof, challenges)
    }
}

