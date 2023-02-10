use std::marker::PhantomData;

use ark_ec::AffineRepr;
use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use fflonk::pcs::{Commitment, PCS, PcsParams};

use common::{Column, FieldColumn, Proof};
use common::gadgets::sw_cond_add::AffineColumn;

use crate::piop::{RingCommitments, RingEvaluations};
use crate::piop::params::PiopParams;

mod piop;
pub mod ring_prover;
pub mod ring_verifier;

type RingProof<F, CS> = Proof<F, CS, RingCommitments<F, <CS as PCS<F>>::C>, RingEvaluations<F>>;

// Columns commitment to which the verifier knows (or trusts).
pub struct FixedColumns<F: PrimeField, G: AffineRepr<BaseField=F>> {
    // Public keys of the ring participants in order,
    // followed by the powers-of-2 multiples of the second Pedersen base.
    // pk_1, ..., pk_n, H, 2H, 4H, ..., 2^sH
    // 1          n                     n+s+1
    points: AffineColumn<F, G>,
    // Binary column that highlights which rows of the table correspond to the ring.
    // 1, 1, ..., 1, 0, 0, ..., 0
    // 1          n
    ring_selector: FieldColumn<F>,
}

// Commitments to the fixed columns (see above).
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct FixedColumnsCommitted<F: PrimeField, C: Commitment<F>> {
    points: [C; 2],
    ring_selector: C,
    phantom: PhantomData<F>,
}

impl<F: PrimeField, G: AffineRepr<BaseField=F>> FixedColumns<F, G> {
    fn commit<CS: PCS<F>>(&self, ck: &CS::CK) -> FixedColumnsCommitted<F, CS::C> {
        let points = [
            CS::commit(ck, self.points.xs.as_poly()),
            CS::commit(ck, self.points.ys.as_poly()),
        ];
        let ring_selector = CS::commit(ck, self.ring_selector.as_poly());
        FixedColumnsCommitted { points, ring_selector, phantom: Default::default() }
    }
}

pub struct ProverKey<F: PrimeField, CS: PCS<F>, G: AffineRepr<BaseField=F>> {
    pcs_ck: CS::CK,
    fixed_columns: FixedColumns<F, G>,
    verifier_key: VerifierKey<F, CS>, // used in the Fiat-Shamir transform
}


#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerifierKey<F: PrimeField, CS: PCS<F>> {
    pcs_raw_vk: <CS::Params as PcsParams>::RVK,
    fixed_columns_committed: FixedColumnsCommitted<F, CS::C>,
}

pub fn index<F: PrimeField, CS: PCS<F>, Curve: SWCurveConfig<BaseField=F>>(
    pcs_params: CS::Params,
    piop_params: &PiopParams<F, Curve>,
    keys: Vec<Affine<Curve>>,
) -> (ProverKey<F, CS, Affine<Curve>>, VerifierKey<F, CS>) {
    let pcs_ck = pcs_params.ck();
    let pcs_raw_vk = pcs_params.raw_vk();
    let fixed_columns = piop_params.fixed_columns(&keys);
    let fixed_columns_committed = fixed_columns.commit::<CS>(&pcs_ck);
    let verifier_key = VerifierKey {
        pcs_raw_vk: pcs_raw_vk.clone(),
        fixed_columns_committed: fixed_columns_committed.clone()
    };
    let prover_key = ProverKey { pcs_ck, fixed_columns, verifier_key };
    let verifier_key = VerifierKey { pcs_raw_vk, fixed_columns_committed };
    (prover_key, verifier_key)
}

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
