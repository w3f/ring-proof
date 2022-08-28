use ark_ec::AffineCurve;
use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::PrimeField;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use fflonk::pcs::{PCS, PcsParams};
use common::domain::Domain;

use common::gadgets::sw_cond_add::CondAdd;
use common::piop::VerifierPiop;
use common::verifier::PlonkVerifier;

use crate::piop::PiopVerifier;
use crate::piop::SelectorColumns;
use crate::RingProof;

pub struct RingVerifier<F: PrimeField, CS: PCS<F>> {
    points_comm: [CS::C; 2],
    domain_size: usize,
    keyset_size: usize,

    plonk_verifier: PlonkVerifier<F, CS, merlin::Transcript>,
}

impl<F: PrimeField, CS: PCS<F>> RingVerifier<F, CS> {
    pub fn init(raw_vk: &<CS::Params as PcsParams>::RVK,
                points_comm: [CS::C; 2],
                domain_size: usize,
                keyset_size: usize,
                empty_transcript: merlin::Transcript,
    ) -> Self {
        let domain = GeneralEvaluationDomain::new(domain_size).unwrap();
        let plonk_verifier = PlonkVerifier::init(raw_vk, domain, points_comm.clone(), empty_transcript);

        Self {
            points_comm,
            domain_size,
            keyset_size,
            plonk_verifier,
        }
    }

    pub fn verify_ring_proof<Curve: SWCurveConfig<BaseField=F>>(&self, proof: RingProof<F, CS>, result: Affine<Curve>) -> bool {
        let challenges = self.plonk_verifier.restore_challenges(
            &result,
            &proof,
            // '1' accounts for the quotient polynomial that is aggregated together with the columns
            PiopVerifier::<F, CS::C>::N_COLUMNS + 1,
            PiopVerifier::<F, CS::C>::N_CONSTRAINTS,
        );
        let domain = Domain::new(self.domain_size, false);
        let selectors = SelectorColumns::init(&domain, self.keyset_size);
        let selectors_at_zeta = selectors.evaluate(&challenges.zeta);
        let init = CondAdd::<F, Affine<Curve>>::point_in_g1_complement();
        let init_plus_result = init + result;
        let piop = PiopVerifier::init(&self.points_comm,
                                      proof.column_commitments.clone(),
                                      proof.columns_at_zeta.clone(),
                                      selectors_at_zeta,
                                      (init.x, init.y),
                                      (init_plus_result.x, init_plus_result.y),
        );
        self.plonk_verifier.verify(piop, proof, challenges)
    }
}

