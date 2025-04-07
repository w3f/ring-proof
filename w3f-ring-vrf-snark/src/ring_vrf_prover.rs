use ark_ec::twisted_edwards::{Affine, TECurveConfig};
use ark_ff::PrimeField;
use w3f_pcs::pcs::PCS;
use w3f_plonk_common::piop::ProverPiop;
use w3f_plonk_common::prover::PlonkProver;
use w3f_plonk_common::transcript::PlonkTranscript;

use crate::piop::params::PiopParams;
use crate::piop::{FixedColumns, PiopProver, ProverKey};
use crate::{ArkTranscript, RingVrfProof};

pub struct RingVrfProver<F, CS, Curve, T = ArkTranscript>
where
    F: PrimeField,
    CS: PCS<F>,
    Curve: TECurveConfig<BaseField = F>,
    T: PlonkTranscript<F, CS>,
{
    piop_params: PiopParams<F, Curve>,
    fixed_columns: FixedColumns<F, Affine<Curve>>,
    pk_index: usize,
    sk: Curve::ScalarField,
    plonk_prover: PlonkProver<F, CS, T>,
}

impl<F, CS, Curve, T> RingVrfProver<F, CS, Curve, T>
where
    F: PrimeField,
    CS: PCS<F>,
    Curve: TECurveConfig<BaseField = F>,
    T: PlonkTranscript<F, CS>,
{
    pub fn init(
        prover_key: ProverKey<F, CS, Affine<Curve>>,
        piop_params: PiopParams<F, Curve>,
        pk_index: usize,
        sk: Curve::ScalarField,
        empty_transcript: T,
    ) -> Self {
        let ProverKey {
            pcs_ck,
            fixed_columns,
            verifier_key,
        } = prover_key;

        let plonk_prover = PlonkProver::init(pcs_ck, verifier_key, empty_transcript);

        Self {
            piop_params,
            fixed_columns,
            pk_index,
            sk,
            plonk_prover,
        }
    }

    pub fn prove(&self, vrf_in: Affine<Curve>) -> (Affine<Curve>, RingVrfProof<F, CS>) {
        let piop: PiopProver<F, Curve> = PiopProver::build(
            &self.piop_params,
            self.fixed_columns.clone(),
            self.pk_index,
            self.sk,
            vrf_in,
        );
        let vrf_out = <PiopProver<F, Curve> as ProverPiop<F, CS::C>>::result(&piop);
        let proof = self.plonk_prover.prove(piop);
        (vrf_out, proof)
    }

    pub fn piop_params(&self) -> &PiopParams<F, Curve> {
        &self.piop_params
    }
}
