use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::PrimeField;
use ark_poly::EvaluationDomain;
use ark_std::test_rng;
use fflonk::pcs::PCS;
use common::domain::Domain;

use common::gadgets::sw_cond_add::AffineColumn;
use common::piop::ProverPiop;
use common::prover::PlonkProver;
use common::setup::Setup;

use crate::piop::params::PiopParams;
use crate::piop::PiopProver;
use crate::piop::SelectorColumns;
use crate::RingProof;

pub struct RingProver<F: PrimeField, CS: PCS<F>, Curve: SWCurveConfig<BaseField=F>> {
    piop_params: PiopParams<F, Curve>,
    selectors: SelectorColumns<F>,
    points: AffineColumn<F, Affine<Curve>>,
    k: usize,

    plonk_prover: PlonkProver<F, CS, merlin::Transcript>,
}


impl<F: PrimeField, CS: PCS<F>, Curve: SWCurveConfig<BaseField=F>> RingProver<F, CS, Curve> {
    pub fn init(setup: Setup<F, CS>,
                piop_params: PiopParams<F, Curve>,
                keys: Vec<Affine<Curve>>,
                k: usize,
                empty_transcript: merlin::Transcript,
    ) -> Self {
        let domain = Domain::new(piop_params.domain.size());
        let selectors = SelectorColumns::init(&domain, piop_params.keyset_part_size);
        let points = PiopProver::keyset_column(&domain, &piop_params, &keys);
        let points_comm = [setup.commit_to_column(&points.xs), setup.commit_to_column(&points.ys)];

        let plonk_prover = PlonkProver::init(setup, &points_comm, empty_transcript);

        Self {
            piop_params,
            selectors,
            points,
            k,
            plonk_prover,
        }
    }


    pub fn prove(&self, t: Curve::ScalarField) -> RingProof<F, CS> {
        let domain = Domain::new(self.piop_params.domain.size());
        let piop = PiopProver::init(&domain, &self.piop_params, self.selectors.clone(), self.points.clone(), self.k, t);
        self.plonk_prover.prove(piop)
    }
}

