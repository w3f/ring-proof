use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::PrimeField;
use ark_std::Zero;
use fflonk::pcs::{PCS, PcsParams};

use common::gadgets::sw_cond_add::AffineColumn;
use common::piop::ProverPiop;
use common::prover::PlonkProver;
use common::setup::Setup;

use crate::piop::params::PiopParams;
use crate::piop::{FixedColumns, Indexer, PiopProver};
use crate::{RingCommitments, RingProof};


pub struct RingProver<F: PrimeField, CS: PCS<F>, Curve: SWCurveConfig<BaseField=F>> {
    piop_params: PiopParams<F, Curve>,
    fixed_columns: FixedColumns<F, Curve>,
    k: usize,
    plonk_prover: PlonkProver<F, CS, merlin::Transcript>,
    pub cols: RingCommitments<F, CS::C>,
}


impl<F: PrimeField, CS: PCS<F>, Curve: SWCurveConfig<BaseField=F>> RingProver<F, CS, Curve> {
    pub fn init(setup: Setup<F, CS>,
                piop_params: PiopParams<F, Curve>,
                index: Indexer<F, Curve, CS>,
                k: usize,
                empty_transcript: merlin::Transcript,
    ) -> Self {
        let fixed_columns  = index.fixed_columns;

        let piop = PiopProver::build(&piop_params, &fixed_columns, k, Curve::ScalarField::zero());
        let cols = piop.committed_columns(|p| CS::commit(&setup.pcs_params.ck(), p));
        let plonk_prover = PlonkProver::init(setup, &index.fixed_columns_committed.points_committed, empty_transcript);
        Self {
            piop_params,
            fixed_columns,
            k,
            plonk_prover,
            cols,
        }
    }


    pub fn prove(&self, t: Curve::ScalarField) -> RingProof<F, CS> {
        let piop = PiopProver::build(&self.piop_params, &self.fixed_columns, self.k, t);
        self.plonk_prover.prove(piop)
    }
}

