#![cfg_attr(not(feature = "std"), no_std)]

use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ec::AffineRepr;
use ark_ff::{Field, One, Zero};
use ark_std::rand;
use fflonk::pcs::PCS;

pub use common::domain::Domain;
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

/// Transcript for `RingProver` and `RingVerifier` construction.
pub use merlin::Transcript;

// Calling the method for a prime-order curve results in an infinite loop.
pub fn find_complement_point<Curve: SWCurveConfig>() -> Affine<Curve> {
    let mut x = Curve::BaseField::zero();
    loop {
        let p = Affine::<Curve>::get_point_from_x_unchecked(x, false);
        if p.is_some() && !p.unwrap().is_in_correct_subgroup_assuming_on_curve() {
            return p.unwrap();
        }
        x = x + Curve::BaseField::one()
    }
}

pub fn find_random_point<F: Field, P: AffineRepr<BaseField = F>>() -> P {
    let mut x: u8 = 0;
    loop {
        let p = P::from_random_bytes(&[x]);
        if p.is_some() && !p.unwrap().is_zero() {
            // && !p.unwrap().is_in_correct_subgroup_assuming_on_curve() {
            return p.unwrap().clear_cofactor();
        }
        x = x + 1;
    }
}

// TODO: switch to better hash to curve when available
pub fn hash_to_curve<A: AffineRepr>(message: &[u8]) -> A {
    use ark_std::rand::SeedableRng;
    use blake2::Digest;

    let seed = blake2::Blake2s::digest(message);
    let rng = &mut rand::rngs::StdRng::from_seed(seed.into());
    A::rand(rng)
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Bls12_381;
    use ark_ec::CurveGroup;
    use ark_ed_on_bls12_381_bandersnatch::{BandersnatchConfig, EdwardsAffine, Fq, Fr, SWAffine};
    use ark_ff::MontFp;
    use ark_std::rand::Rng;
    use ark_std::{end_timer, start_timer, test_rng, UniformRand};
    use fflonk::pcs::kzg::KZG;
    use merlin::Transcript;

    use common::test_helpers::random_vec;

    use crate::piop::FixedColumnsCommitted;
    use crate::ring::{Ring, RingBuilderKey};
    use crate::ring_prover::RingProver;
    use crate::ring_verifier::RingVerifier;
    use common::gadgets::cond_add::CondAdd;
    use common::gadgets::sw_cond_add::SwCondAdd;
    use common::gadgets::te_cond_add::TeCondAdd;
    use common::gadgets::ProverGadget;
    use common::gadgets::VerifierGadget;

    use super::*;

    fn _test_ring_proof<
        CS: PCS<Fq>,
        P: AffineRepr<BaseField = Fq, ScalarField = Fr>,
        CondAddT: CondAdd<Fq, P> + ProverGadget<Fq>,
    >(
        domain_size: usize,
    ) where
        CondAddT::CondAddValT: VerifierGadget<Fq>,
    {
        let rng = &mut test_rng();

        let (pcs_params, piop_params) = setup::<_, CS, P>(rng, domain_size);

        let max_keyset_size = piop_params.keyset_part_size;
        let keyset_size: usize = rng.gen_range(0..max_keyset_size);
        let pks = random_vec::<P, _>(keyset_size, rng);
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
            Transcript::new(b"ring-vrf-test"),
        );
        let t_prove = start_timer!(|| "Prove");
        let proof = ring_prover.prove::<CondAddT>(secret);
        end_timer!(t_prove);

        let ring_verifier =
            RingVerifier::init(verifier_key, piop_params, Transcript::new(b"ring-vrf-test"));
        let t_verify = start_timer!(|| "Verify");
        let res =
            ring_verifier.verify_ring_proof::<CondAddT::CondAddValT>(proof, result.into_affine());
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
        let seed = find_random_point::<Fq, P>();
        let piop_params = PiopParams::setup(domain, h, seed);

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
    fn test_complement_point() {
        let p = find_complement_point::<BandersnatchConfig>();
        assert!(p.is_on_curve());
        assert!(!p.is_in_correct_subgroup_assuming_on_curve());
        assert_eq!(
            p,
            SWAffine::new_unchecked(
                MontFp!("0"),
                MontFp!(
                    "11982629110561008531870698410380659621661946968466267969586599013782997959645"
                )
            )
        )
    }

    #[test]
    fn test_ring_proof_kzg_sw() {
        _test_ring_proof::<KZG<Bls12_381>, SWAffine, SwCondAdd<Fq, SWAffine>>(2usize.pow(10));
    }

    #[test]
    fn test_ring_proof_kzg_te() {
        _test_ring_proof::<KZG<Bls12_381>, EdwardsAffine, TeCondAdd<Fq, EdwardsAffine>>(
            2usize.pow(10),
        );
    }

    #[test]
    fn test_ring_proof_id_sw() {
        _test_ring_proof::<fflonk::pcs::IdentityCommitment, SWAffine, SwCondAdd<Fq, SWAffine>>(
            2usize.pow(10),
        );
    }

    #[test]
    fn test_ring_proof_id_te() {
        _test_ring_proof::<
            fflonk::pcs::IdentityCommitment,
            EdwardsAffine,
            TeCondAdd<Fq, EdwardsAffine>,
        >(2usize.pow(10));
    }
}
