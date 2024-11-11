#![cfg_attr(not(feature = "std"), no_std)]

use ark_ec::AffineRepr;
use ark_ff::PrimeField;
use ark_serialize::CanonicalSerialize;
use ark_std::rand::RngCore;
use fflonk::pcs::PCS;

pub use common::domain::Domain;
pub use common::gadgets::cond_add::AffineCondAdd;

use common::Proof;
pub use piop::index;

pub use crate::piop::{params::PiopParams, FixedColumnsCommitted, ProverKey, VerifierKey};
use crate::piop::{RingCommitments, RingEvaluations};

mod piop;
pub mod ring;
pub mod ring_prover;
pub mod ring_verifier;

pub type RingProof<F, CS> = Proof<F, CS, RingCommitments<F, <CS as PCS<F>>::C>, RingEvaluations<F>>;

/// Polynomial Commitment Schemes.
pub use fflonk::pcs;

#[derive(Clone)]
pub struct ArkTranscript(ark_transcript::Transcript);

impl<F: PrimeField, CS: PCS<F>> common::transcript::PlonkTranscript<F, CS> for ArkTranscript {
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

/// Simple try and increment hash to curve utility.
///
/// This can be used to generage a point within the prime order subgroup
/// given a seed message.
pub fn hash_to_curve<F: PrimeField, P: AffineRepr<BaseField = F>>(message: &[u8]) -> P {
    let mut t = ark_transcript::Transcript::new_blank();
    t.append_slice(message);
    let mut rng = t.challenge(b"rng");
    P::rand(&mut rng)
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Bls12_381;
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ed_on_bls12_381_bandersnatch::{EdwardsAffine, Fq, Fr, SWAffine};
    use ark_std::rand::Rng;
    use ark_std::{end_timer, start_timer, test_rng, UniformRand};
    use fflonk::pcs::kzg::KZG;

    use common::test_helpers::random_vec;

    use crate::piop::FixedColumnsCommitted;
    use crate::ring::{Ring, RingBuilderKey};
    use crate::ring_prover::RingProver;
    use crate::ring_verifier::RingVerifier;
    use common::gadgets::cond_add::AffineCondAdd;

    #[cfg(feature = "intensive-benchmarking")]
    use std::hint::black_box;

    use super::*;

    fn _test_ring_proof<CS, P>(domain_size: usize, _repeat: usize)
    where
        CS: PCS<Fq>,
        P: AffineCondAdd<BaseField = Fq, ScalarField = Fr>,
    {
        let rng = &mut test_rng();

        let (pcs_params, piop_params) = setup::<_, CS, P>(rng, domain_size);

        let max_keyset_size = piop_params.keyset_part_size;
        let keyset_size: usize = rng.gen_range(0..max_keyset_size);
        let pks = random_vec::<P, _>(keyset_size, rng);
        let k = rng.gen_range(0..keyset_size); // prover's secret index
        let pk = pks[k];

        let (prover_key, verifier_key) = index::<_, CS, _>(&pcs_params, &piop_params, &pks);

        // PROOF generation
        let secret = Fr::rand(rng); // prover's secret scalar
        let result = piop_params.h.mul(secret) + pk;
        let ring_prover = RingProver::init(
            prover_key,
            piop_params.clone(),
            k,
            ArkTranscript::new(b"ring-vrf-test"),
        );
        let t_prove = start_timer!(|| "Prove");

        #[cfg(feature = "intensive-benchmarking")]
        let mut proofs: Vec<RingProof<_, CS>> = vec![];

        #[cfg(feature = "intensive-benchmarking")]
        for _ in 0.._repeat - 1 {
            let proof = black_box(ring_prover.prove(secret));
            proofs.push(proof);
        }

        let proof = ring_prover.prove(secret);
        end_timer!(t_prove);

        let ring_verifier = RingVerifier::init(
            verifier_key,
            piop_params,
            ArkTranscript::new(b"ring-vrf-test"),
        );
        let t_verify = start_timer!(|| "Verify");

        #[cfg(feature = "intensive-benchmarking")]
        for _ in 0.._repeat - 1 {
            let proof = proofs.pop().unwrap();
            black_box(ring_verifier.verify(proof, result.into_affine()));
        }
        let res = ring_verifier.verify(proof, result.into_affine());
        end_timer!(t_verify);
        assert!(res);
    }

    fn _test_lagrangian_commitment<P: AffineRepr<BaseField = Fq>>() {
        let rng = &mut test_rng();

        let domain_size = 2usize.pow(9);

        let (pcs_params, piop_params) = setup::<_, KZG<Bls12_381>, P>(rng, domain_size);
        let ring_builder_key = RingBuilderKey::from_srs(&pcs_params, domain_size);

        let max_keyset_size = piop_params.keyset_part_size;
        let keyset_size: usize = rng.gen_range(0..max_keyset_size);
        let pks = random_vec::<P, _>(keyset_size, rng);

        let (_, verifier_key) = index::<_, KZG<Bls12_381>, _>(&pcs_params, &piop_params, &pks);

        let ring = Ring::<_, Bls12_381, _>::with_keys(&piop_params, &pks, &ring_builder_key);

        let fixed_columns_committed = FixedColumnsCommitted::from_ring(&ring);
        assert_eq!(
            fixed_columns_committed,
            verifier_key.fixed_columns_committed
        );
    }

    fn setup<R: Rng, CS: PCS<Fq>, P: AffineRepr<BaseField = Fq>>(
        rng: &mut R,
        domain_size: usize,
    ) -> (CS::Params, PiopParams<Fq, P>) {
        let setup_degree = 3 * domain_size;
        let pcs_params = CS::setup(setup_degree, rng);

        let domain = Domain::new(domain_size, true);
        let h = P::rand(rng);
        let seed = P::rand(rng);
        let pad = P::rand(rng);
        let piop_params = PiopParams::setup(domain, h, seed, pad);

        (pcs_params, piop_params)
    }

    #[test]
    fn test_lagrangian_commitment_sw() {
        _test_lagrangian_commitment::<SWAffine>();
    }

    #[test]
    fn test_lagrangian_commitment_te() {
        _test_lagrangian_commitment::<SWAffine>();
    }

    #[test]
    fn test_ring_proof_kzg_sw() {
        _test_ring_proof::<KZG<Bls12_381>, SWAffine>(2usize.pow(10), 1);
    }

    #[test]
    fn test_ring_proof_kzg_te() {
        _test_ring_proof::<KZG<Bls12_381>, EdwardsAffine>(2usize.pow(10), 1);
    }

    #[test]
    fn test_ring_proof_id_sw() {
        _test_ring_proof::<fflonk::pcs::IdentityCommitment, SWAffine>(2usize.pow(10), 1);
    }

    #[test]
    fn test_ring_proof_id_te() {
        _test_ring_proof::<fflonk::pcs::IdentityCommitment, EdwardsAffine>(2usize.pow(10), 1);
    }

    #[cfg(feature = "intensive-benchmarking")]
    #[test]
    fn test_16k_ring_10_proofs_kzg_sw() {
        _test_ring_proof::<KZG<Bls12_381>, SWAffine>(2usize.pow(14), 10);
    }

    #[cfg(feature = "intensive-benchmarking")]
    #[test]
    fn test_16k_ring_10_proofs_kzg_te() {
        _test_ring_proof::<KZG<Bls12_381>, EdwardsAffine>(2usize.pow(14), 10);
    }
}
