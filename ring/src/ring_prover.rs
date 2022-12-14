use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::PrimeField;
use fflonk::pcs::PCS;

use common::gadgets::sw_cond_add::AffineColumn;
use common::prover::PlonkProver;
use common::setup::Setup;

use crate::piop::params::PiopParams;
use crate::piop::PiopProver;
use crate::RingProof;

pub struct RingProver<F: PrimeField, CS: PCS<F>, Curve: SWCurveConfig<BaseField=F>> {
    piop_params: PiopParams<F, Curve>,
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
        let points = PiopProver::keyset_column(&piop_params, &keys);
        let points_comm = [setup.commit_to_column(&points.xs), setup.commit_to_column(&points.ys)];

        let plonk_prover = PlonkProver::init(setup, &points_comm, empty_transcript);

        Self {
            piop_params,
            points,
            k,
            plonk_prover,
        }
    }


    pub fn prove(&self, t: Curve::ScalarField) -> RingProof<F, CS> {
        let piop = PiopProver::build(&self.piop_params, self.points.clone(), self.k, t);
        self.plonk_prover.prove(piop)
    }
}

