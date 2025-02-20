#![cfg_attr(not(feature = "std"), no_std)]

use ark_ff::PrimeField;
use ark_serialize::CanonicalSerialize;
use ark_std::rand::RngCore;
use w3f_pcs::pcs::PCS;

pub use piop::index;
pub use w3f_plonk_common::domain::Domain;
use w3f_plonk_common::Proof;

pub use crate::piop::{params::PiopParams, FixedColumnsCommitted, ProverKey, VerifierKey};
use crate::piop::{RingCommitments, RingEvaluations};

mod piop;
pub mod ring;
pub mod ring_prover;
pub mod ring_verifier;

pub type RingProof<F, CS> = Proof<F, CS, RingCommitments<F, <CS as PCS<F>>::C>, RingEvaluations<F>>;

/// Polynomial Commitment Schemes.
pub use w3f_pcs::pcs;

#[derive(Clone)]
pub struct ArkTranscript(ark_transcript::Transcript);

impl<F: PrimeField, CS: PCS<F>> w3f_plonk_common::transcript::PlonkTranscript<F, CS>
    for ArkTranscript
{
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
    use ark_ec::twisted_edwards::{Affine, TECurveConfig};
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ed_on_bls12_381_bandersnatch::{BandersnatchConfig, EdwardsAffine, Fq, Fr};
    use ark_std::ops::Mul;
    use ark_std::rand::Rng;
    use ark_std::{end_timer, start_timer, test_rng, UniformRand};
    use w3f_pcs::pcs::kzg::KZG;

    use w3f_plonk_common::test_helpers::random_vec;

    use crate::piop::FixedColumnsCommitted;
    use crate::ring::{Ring, RingBuilderKey};
    use crate::ring_prover::RingProver;
    use crate::ring_verifier::RingVerifier;

    use super::*;

    // Try and increment hash to curve.
    fn hash_to_curve<F: PrimeField, Curve: TECurveConfig<BaseField = F>>(
        message: &[u8],
    ) -> Affine<Curve> {
        use blake2::Digest;
        let mut seed = message.to_vec();
        let cnt_offset = seed.len();
        seed.push(0);
        loop {
            let hash: [u8; 64] = blake2::Blake2b::digest(&seed[..]).into();
            let x = F::from_le_bytes_mod_order(&hash);
            if let Some(point) = Affine::<Curve>::get_point_from_y_unchecked(x, false) {
                let point = point.clear_cofactor();
                assert!(point.is_in_correct_subgroup_assuming_on_curve());
                return point;
            }
            seed[cnt_offset] += 1;
        }
    }

    pub(crate) fn padding_point<F: PrimeField, Curve: TECurveConfig<BaseField = F>>(
    ) -> Affine<Curve> {
        hash_to_curve(b"/w3f/w3f-ring-proof/padding")
    }

    fn setup<R: Rng, CS: PCS<Fq>>(
        rng: &mut R,
        domain_size: usize,
    ) -> (CS::Params, PiopParams<Fq, BandersnatchConfig>) {
        let setup_degree = 3 * domain_size;
        let pcs_params = CS::setup(setup_degree, rng);

        let domain = Domain::new(domain_size, true);
        let h = EdwardsAffine::rand(rng);
        let seed = EdwardsAffine::rand(rng);
        let pad = padding_point();
        let piop_params = PiopParams::setup(domain, h, seed, pad);

        (pcs_params, piop_params)
    }

    fn _test_ring_proof<CS: PCS<Fq>>(domain_size: usize) {
        let rng = &mut test_rng();

        let (pcs_params, piop_params) = setup::<_, CS>(rng, domain_size);

        let max_keyset_size = piop_params.keyset_part_size;
        let keyset_size: usize = rng.gen_range(0..max_keyset_size);
        let pks = random_vec::<EdwardsAffine, _>(keyset_size, rng);
        let k = rng.gen_range(0..keyset_size); // prover's secret index
        let pk = pks[k].clone();

        let (prover_key, verifier_key) = index::<_, CS, _>(&pcs_params, &piop_params, &pks);

        // PROOF generation
        let secret = Fr::rand(rng); // prover's secret scalar
        let result = piop_params.h.mul(secret) + pk;
        let ring_prover = RingProver::init(
            prover_key,
            piop_params.clone(),
            k,
            ArkTranscript::new(b"w3f-ring-proof-test"),
        );
        let t_prove = start_timer!(|| "Prove");
        let proof = ring_prover.prove(secret);
        end_timer!(t_prove);

        let ring_verifier = RingVerifier::init(
            verifier_key,
            piop_params,
            ArkTranscript::new(b"w3f-ring-proof-test"),
        );
        let t_verify = start_timer!(|| "Verify");
        let res = ring_verifier.verify(proof, result.into_affine());
        end_timer!(t_verify);
        assert!(res);
    }

    #[test]
    fn test_lagrangian_commitment() {
        let rng = &mut test_rng();

        let domain_size = 2usize.pow(9);

        let (pcs_params, piop_params) = setup::<_, KZG<Bls12_381>>(rng, domain_size);
        let ring_builder_key = RingBuilderKey::from_srs(&pcs_params, domain_size);

        let max_keyset_size = piop_params.keyset_part_size;
        let keyset_size: usize = rng.gen_range(0..max_keyset_size);
        let pks = random_vec::<EdwardsAffine, _>(keyset_size, rng);

        let (_, verifier_key) = index::<_, KZG<Bls12_381>, _>(&pcs_params, &piop_params, &pks);

        let ring = Ring::<_, Bls12_381, _>::with_keys(&piop_params, &pks, &ring_builder_key);

        let fixed_columns_committed = FixedColumnsCommitted::from_ring(&ring);
        assert_eq!(
            fixed_columns_committed,
            verifier_key.fixed_columns_committed
        );
    }

    #[test]
    fn test_ring_proof_kzg() {
        _test_ring_proof::<KZG<Bls12_381>>(2usize.pow(10));
    }

    #[test]
    fn test_ring_proof_id() {
        _test_ring_proof::<w3f_pcs::pcs::IdentityCommitment>(2usize.pow(10));
    }
}
