use ark_ec::AffineRepr;
use ark_ff::PrimeField;
use fflonk::pcs::PCS;

use common::gadgets::cond_add::AffineCondAdd;
use common::prover::PlonkProver;
use common::transcript::PlonkTranscript;

use crate::piop::{params::PiopParams, FixedColumns, PiopProver, ProverKey};
use crate::{ArkTranscript, RingProof};

pub struct RingProver<
    F: PrimeField,
    CS: PCS<F>,
    P: AffineRepr<BaseField = F>,
    T: PlonkTranscript<F, CS> = ArkTranscript,
> {
    piop_params: PiopParams<F, P>,
    fixed_columns: FixedColumns<F, P>,
    k: usize,
    plonk_prover: PlonkProver<F, CS, T>,
}

impl<F: PrimeField, CS: PCS<F>, P, T: PlonkTranscript<F, CS>> RingProver<F, CS, P, T>
where
    P: AffineCondAdd<BaseField = F>,
{
    pub fn init(
        prover_key: ProverKey<F, CS, P>,
        piop_params: PiopParams<F, P>,
        k: usize,
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
            k,
            plonk_prover,
        }
    }

    pub fn prove(&self, secret_key: P::ScalarField, vrf_input: P) -> RingProof<F, CS> {
        let piop: PiopProver<F, P, P::CondAddT> =
            PiopProver::build(&self.piop_params, self.fixed_columns.clone(), self.k, secret_key, vrf_input);
        self.plonk_prover.prove(piop)
    }

    pub fn piop_params(&self) -> &PiopParams<F, P> {
        &self.piop_params
    }
}
