use fflonk::pcs::PCS;

use common::Proof;

use crate::piop::{RingCommitments, RingEvaluations};

mod piop;
pub mod ring_prover;
pub mod ring_verifier;
mod vrf;

type RingProof<F, CS> = Proof<F, CS, RingCommitments<F, <CS as PCS<F>>::C>, RingEvaluations<F>>;

#[cfg(test)]
mod tests {
    use std::ops::Mul;

    use ark_ec::CurveGroup;
    use ark_ed_on_bls12_381_bandersnatch::{Fq, Fr, SWAffine};
    use ark_std::{end_timer, start_timer, test_rng, UniformRand};
    use ark_std::rand::Rng;
    use fflonk::pcs::PCS;
    use fflonk::pcs::PcsParams;
    use merlin::Transcript;

    use common::domain::Domain;
    use common::setup::Setup;
    use common::test_helpers::*;

    use crate::piop::params::PiopParams;
    use crate::piop::PiopProver;
    use crate::ring_prover::RingProver;
    use crate::ring_verifier::RingVerifier;

    fn _test_ring_proof<CS: PCS<Fq>>(domain_size: usize, hiding: bool) {
        let rng = &mut test_rng();

        // SETUP per curve and domain
        let domain = Domain::new(domain_size, hiding);
        let piop_params = PiopParams::setup(domain.clone(), &mut test_rng());
        let piop_params2 = PiopParams::setup(domain, &mut test_rng());
        assert_eq!(piop_params.h, piop_params2.h);
        let setup = Setup::<Fq, CS>::generate(domain_size, rng);

        let max_keyset_size = piop_params.keyset_part_size;
        let keyset_size: usize = rng.gen_range(0..max_keyset_size);
        let pks = random_vec::<SWAffine, _>(keyset_size, rng);
        let k = rng.gen_range(0..keyset_size); // prover's secret index
        let pk = &pks[k];

        let points = PiopProver::keyset_column( &piop_params, &pks);
        let points_comm = [setup.commit_to_column(&points.xs), setup.commit_to_column(&points.ys)];
        let vk = &setup.pcs_params.raw_vk();

        // PROOF generation
        let secret = Fr::rand(rng); // prover's secret scalar
        let result = piop_params.h.mul(secret) + pk;
        let ring_prover = RingProver::init(setup, piop_params, pks, k, Transcript::new(b"ring-vrf-test"));

        let t_prove = start_timer!(|| "Prove");
        let proof = ring_prover.prove(secret);
        end_timer!(t_prove);


        let ring_verifier = RingVerifier::init(vk, piop_params2, points_comm, domain_size, max_keyset_size, Transcript::new(b"ring-vrf-test"));
        let t_verify = start_timer!(|| "Verify");
        let res = ring_verifier.verify_ring_proof(proof, result.into_affine());
        end_timer!(t_verify);

        assert!(res);
    }

    #[test]
    fn test_ring_proof() {
        // _test_ring_proof::<fflonk::pcs::kzg::KZG<ark_bls12_381::Bls12_381>>(2usize.pow(12));
        _test_ring_proof::<fflonk::pcs::IdentityCommitment>(2usize.pow(10), false);
        _test_ring_proof::<fflonk::pcs::IdentityCommitment>(2usize.pow(10), true);
    }
}
