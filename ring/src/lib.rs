use ark_ec::pairing::Pairing;
use fflonk::pcs::PCS;

use common::Proof;

use crate::piop::{RingCommitments, RingEvaluations};

mod piop;
pub mod ring_prover;
pub mod ring_verifier;

type RingProof<F, CS> = Proof<F, CS, RingCommitments<F, <CS as PCS<F>>::C>, RingEvaluations<F>>;

#[cfg(test)]
mod tests {
    use std::ops::Mul;
    use ark_bls12_381::Bls12_381;

    use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
    use ark_ec::pairing::Pairing;
    use ark_ed_on_bls12_381_bandersnatch::{Fq, Fr, SWAffine};
    use ark_poly::EvaluationDomain;
    use ark_std::{end_timer, start_timer, test_rng, UniformRand};
    use ark_std::rand::Rng;
    use fflonk::pcs::kzg::KZG;
    use fflonk::pcs::kzg::params::KzgCommitterKey;
    use fflonk::pcs::PCS;
    use fflonk::pcs::PcsParams;
    use merlin::Transcript;

    use common::domain::Domain;
    use common::setup::Setup;
    use common::test_helpers::*;

    use crate::piop::params::PiopParams;
    use crate::piop::{Indexer, PiopProver};
    use crate::ring_prover::RingProver;
    use crate::ring_verifier::RingVerifier;

    pub fn lagrangian_basis<E: Pairing, D: EvaluationDomain<E::ScalarField>>(powers: &[E::G1Affine], domain: D) -> Vec<E::G1Affine> {
        assert_eq!(powers.len(), domain.size());
        let mut powers_in_g1_proj: Vec<_> = powers.iter()
            .map(|p| p.into_group())
            .collect();
        domain.ifft_in_place(&mut powers_in_g1_proj);
        E::G1::normalize_batch(&powers_in_g1_proj)
    }

    fn _test_ring_proof(domain_size: usize, hiding: bool) {
        let rng = &mut test_rng();

        // SETUP per curve and domain
        let domain = Domain::new(domain_size, hiding);
        let piop_params = PiopParams::setup(domain.clone(), &mut test_rng());
        let setup = Setup::<Fq, KZG<Bls12_381>>::generate(domain_size, rng);


        let max_keyset_size = piop_params.keyset_part_size;
        let keyset_size: usize = rng.gen_range(0..max_keyset_size);
        let pks = random_vec::<SWAffine, _>(keyset_size, rng);
        let k = rng.gen_range(0..keyset_size); // prover's secret index
        let pk = &pks[k];

        let index = Indexer::new(&setup, &piop_params, &pks);
        let comm = index.fixed_columns_committed.clone();
        let vk = &setup.pcs_params.raw_vk();

        let powers = &setup.pcs_params.powers_in_g1[..domain_size];
        let lag = lagrangian_basis::<Bls12_381, _>(powers, domain.domain());

        // PROOF generation
        let secret = Fr::rand(rng); // prover's secret scalar
        let result = piop_params.h.mul(secret) + pk;
        let ring_prover = RingProver::init(setup, piop_params, index, k, Transcript::new(b"ring-vrf-test"));

        let cols = ring_prover.cols;
        assert_eq!(cols.bits.0, lag[k]);
        // let t_prove = start_timer!(|| "Prove");
        // let proof = ring_prover.prove(secret);
        // end_timer!(t_prove);



        // let ring_verifier = RingVerifier::init(vk, domain.clone(), comm, Transcript::new(b"ring-vrf-test"));
        // let t_verify = start_timer!(|| "Verify");
        // let res = ring_verifier.verify_ring_proof(proof, result.into_affine());
        // end_timer!(t_verify);
        //
        // assert!(res);
    }

    #[test]
    fn test_ring_proof() {
        // _test_ring_proof::<fflonk::pcs::kzg::KZG<ark_bls12_381::Bls12_381>>(2usize.pow(12));
        // _test_ring_proof::<fflonk::pcs::IdentityCommitment>(2usize.pow(10), false);
        // _test_ring_proof::<fflonk::pcs::IdentityCommitment>(2usize.pow(10), true);
        _test_ring_proof(2usize.pow(10), false);
    }
}
