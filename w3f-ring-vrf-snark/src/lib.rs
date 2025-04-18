#![cfg_attr(not(feature = "std"), no_std)]

use crate::piop::{RingCommitments, RingEvaluations};
use ark_ff::PrimeField;
use ark_serialize::CanonicalSerialize;
use ark_std::rand::RngCore;
use w3f_pcs::pcs::PCS;
use w3f_plonk_common::transcript::PlonkTranscript;
use w3f_plonk_common::Proof;

pub use crate::piop::{params::PiopParams, FixedColumnsCommitted, ProverKey, VerifierKey};
pub use piop::index;
pub use w3f_pcs::pcs;
pub use w3f_plonk_common::domain::Domain;

mod piop;
pub mod ring_vrf_prover;
pub mod ring_vrf_verifier;

pub type RingVrfProof<F, CS> =
    Proof<F, CS, RingCommitments<F, <CS as PCS<F>>::C>, RingEvaluations<F>>;

#[derive(Clone)]
pub struct ArkTranscript(ark_transcript::Transcript);

impl<F: PrimeField, CS: PCS<F>> PlonkTranscript<F, CS> for ArkTranscript {
    fn _128_bit_point(&mut self, label: &'static [u8]) -> F {
        self.0.challenge(label).read_reduce()
    }

    fn _add_serializable(&mut self, label: &'static [u8], message: &impl CanonicalSerialize) {
        self.0.label(label);
        self.0.append(message);
    }

    fn to_rng(mut self) -> impl RngCore {
        self.0.challenge(b"transcript_rng")
    }
}

impl ArkTranscript {
    pub fn new(label: &'static [u8]) -> Self {
        Self(ark_transcript::Transcript::new_labeled(label))
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Bls12_381;
    use ark_ec::CurveGroup;
    use ark_ed_on_bls12_381_bandersnatch::{BandersnatchConfig, EdwardsAffine, Fq, Fr};
    use ark_std::ops::Mul;
    use ark_std::rand::Rng;
    use ark_std::{end_timer, start_timer, test_rng, UniformRand};
    use w3f_pcs::pcs::kzg::KZG;
    use w3f_plonk_common::test_helpers::random_vec;

    use crate::ring_vrf_prover::RingVrfProver;
    use crate::ring_vrf_verifier::RingVrfVerifier;

    use super::*;

    fn _test_ring_proof<CS: PCS<Fq>>(domain_size: usize, keyset_size: usize) {
        let rng = &mut test_rng();

        let (pcs_params, piop_params) = setup::<_, CS>(rng, domain_size);

        let mut pks = random_vec::<EdwardsAffine, _>(keyset_size, rng);
        let pk_index = rng.gen_range(0..keyset_size); // prover's secret index
        let sk = Fr::rand(rng); // prover's secret scalar
        let pk = piop_params.g.mul(sk).into_affine();
        pks[pk_index] = pk;

        let (prover_key, verifier_key) = index::<_, CS, _>(&pcs_params, &piop_params, &pks);

        let vrf_in = EdwardsAffine::rand(rng);
        let vrf_out = vrf_in.mul(sk).into_affine();

        let ring_prover = RingVrfProver::init(
            prover_key,
            piop_params.clone(),
            pk_index,
            sk,
            ArkTranscript::new(b"w3f-ring-vrf-snark-test"),
        );
        let t_prove = start_timer!(|| "Prove");
        let (res, proof) = ring_prover.prove(vrf_in);
        end_timer!(t_prove);
        assert_eq!(res, vrf_out);

        let ring_verifier = RingVrfVerifier::init(
            verifier_key,
            piop_params,
            ArkTranscript::new(b"w3f-ring-vrf-snark-test"),
        );
        let t_verify = start_timer!(|| "Verify");
        let res = ring_verifier.verify(vrf_in, vrf_out, proof);
        end_timer!(t_verify);
        assert!(res);
    }

    // #[test]
    // fn test_lagrangian_commitment() {
    //     let rng = &mut test_rng();
    //
    //     let domain_size = 2usize.pow(9);
    //
    //     let (pcs_params, piop_params) = setup::<_, KZG<Bls12_381>>(rng, domain_size);
    //     let ring_builder_key = RingBuilderKey::from_srs(&pcs_params, domain_size);
    //
    //     let max_keyset_size = piop_params.keyset_part_size;
    //     let keyset_size: usize = rng.gen_range(0..max_keyset_size);
    //     let pks = random_vec::<EdwardsAffine, _>(keyset_size, rng);
    //
    //     let (_, verifier_key) = index::<_, KZG<Bls12_381>, _>(&pcs_params, &piop_params, &pks);
    //
    //     let ring = Ring::<_, Bls12_381, _>::with_keys(&piop_params, &pks, &ring_builder_key);
    //
    //     let fixed_columns_committed = FixedColumnsCommitted::from_ring(&ring);
    //     assert_eq!(
    //         fixed_columns_committed,
    //         verifier_key.fixed_columns_committed
    //     );
    // }

    fn setup<R: Rng, CS: PCS<Fq>>(
        rng: &mut R,
        domain_size: usize,
    ) -> (CS::Params, PiopParams<Fq, BandersnatchConfig>) {
        let setup_degree = 3 * domain_size;
        let pcs_params = CS::setup(setup_degree, rng);
        let domain = Domain::new(domain_size, true);
        let seed = EdwardsAffine::rand(rng);
        let padding = EdwardsAffine::rand(rng);
        let piop_params = PiopParams::setup(domain, seed, padding);
        (pcs_params, piop_params)
    }

    #[test]
    fn test_ring_proof_kzg() {
        _test_ring_proof::<KZG<Bls12_381>>(2usize.pow(9), 500);
    }

    #[test]
    fn test_ring_proof_id() {
        _test_ring_proof::<pcs::IdentityCommitment>(2usize.pow(10), 2usize.pow(10) - 4); // no padding
        _test_ring_proof::<pcs::IdentityCommitment>(2usize.pow(9), 253);
        _test_ring_proof::<pcs::IdentityCommitment>(2usize.pow(9), 1);
    }
}
