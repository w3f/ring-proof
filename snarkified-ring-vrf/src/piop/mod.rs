use ark_ec::pairing::Pairing;
use ark_ec::{AdditiveGroup, AffineRepr, CurveGroup};
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::marker::PhantomData;
use ark_std::{vec, vec::Vec};
use fflonk::pcs::kzg::commitment::KzgCommitment;
use fflonk::pcs::kzg::params::RawKzgVerifierKey;
use fflonk::pcs::kzg::KZG;
use fflonk::pcs::{Commitment, PcsParams, PCS};

use common::AffineColumn;
use common::{Column, ColumnsCommited, ColumnsEvaluated, FieldColumn};
pub(crate) use prover::PiopProver;
pub(crate) use verifier::PiopVerifier;

use crate::ring_vrf::Ring;
use crate::PiopParams;

pub mod params;
mod prover;
mod verifier;

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct RingCommitments<F: PrimeField, C: Commitment<F>> {
    pub(crate) signer_index: C,
    pub(crate) signer_secret_key_bits: C,
    pub(crate) sole_signer_inn_prod_acc: C,
    pub(crate) cond_add_pubkey_acc: [C; 2],
    pub(crate) cond_add_gen_multiples_acc: [C; 2],
    pub(crate) cond_add_vrfout_acc: [C; 2],
    pub(crate) phantom: PhantomData<F>,
}

impl<F: PrimeField, C: Commitment<F>> ColumnsCommited<F, C> for RingCommitments<F, C> {
    fn to_vec(self) -> Vec<C> {
        vec![
            self.bits,
            self.inn_prod_acc,
            self.cond_add_acc[0].clone(),
            self.cond_add_acc[1].clone(),
        ]
    }
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct RingEvaluations<F: PrimeField> {
    pub(crate) pubkey_points: [F; 2],
    pub(crate) ring_selector: F,
    pub(crate) signer_index: F,
    pub(crate) signer_secret_key_bits: F,
    pub(crate) sole_signer_inn_prod_acc: F,
    pub(crate) cond_add_pubkey_acc: [F; 2],
    pub(crate) cond_add_gen_multiples_acc: [F; 2],
    pub(crate) cond_add_vrfout_acc: [F; 2],
}

impl<F: PrimeField> ColumnsEvaluated<F> for RingEvaluations<F> {
    fn to_vec(self) -> Vec<F> {
        vec![
            self.points[0],
            self.points[1],
            self.ring_selector,
            self.bits,
            self.inn_prod_acc,
            self.cond_add_acc[0],
            self.cond_add_acc[1],
        ]
    }
}

// Columns commitment to which the verifier knows (or trusts).
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct FixedColumns<F: PrimeField, G: AffineRepr<BaseField = F>> {
    // Public keys of the ring participants in order,
    // pk_1, ..., pk_n, H, 2H, 4H, ..., 2^sH
    pubkey_points: AffineColumn<F, G>,
    // the powers-of-2 multiples of the prime subgroup generator
    // 1          n                     n+s+1
    power_of_2_multiples_of_gen: AffineColumn<F, G>,
}

// Commitments to the fixed columns (see above).
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize, PartialEq, Eq, Debug)]
pub struct FixedColumnsCommitted<F: PrimeField, C: Commitment<F>> {
    pub points: [C; 2],
    ring_selector: C,
    phantom: PhantomData<F>,
}

impl<F: PrimeField, C: Commitment<F>> FixedColumnsCommitted<F, C> {
    fn as_vec(&self) -> Vec<C> {
        vec![
            self.points[0].clone(),
            self.points[1].clone(),
            self.ring_selector.clone(),
        ]
    }
}

impl<E: Pairing> FixedColumnsCommitted<E::ScalarField, KzgCommitment<E>> {
    pub fn from_ring<P: AffineRepr<BaseField = E::ScalarField>>(
        ring: &Ring<E::ScalarField, E, P>,
    ) -> Self {
        let cx = KzgCommitment(ring.cx);
        let cy = KzgCommitment(ring.cy);
        Self {
            points: [cx, cy],
            ring_selector: KzgCommitment(ring.selector),
            phantom: Default::default(),
        }
    }
}

impl<F: PrimeField, P: AffineRepr<BaseField = F>> FixedColumns<F, P> {
    fn commit<CS: PCS<F>>(&self, ck: &CS::CK) -> FixedColumnsCommitted<F, CS::C> {
        let points = [
            CS::commit(ck, self.points.xs.as_poly()),
            CS::commit(ck, self.points.ys.as_poly()),
        ];
        let ring_selector = CS::commit(ck, self.ring_selector.as_poly());
        FixedColumnsCommitted {
            points,
            ring_selector,
            phantom: Default::default(),
        }
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct ProverKey<F: PrimeField, CS: PCS<F>, P: AffineRepr<BaseField = F>> {
    pub(crate) pcs_ck: CS::CK,
    pub(crate) fixed_columns: FixedColumns<F, P>,
    pub(crate) verifier_key: VerifierKey<F, CS>, // used in the Fiat-Shamir transform
}

#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerifierKey<F: PrimeField, CS: PCS<F>> {
    pub(crate) pcs_raw_vk: <CS::Params as PcsParams>::RVK,
    pub(crate) fixed_columns_committed: FixedColumnsCommitted<F, CS::C>,
    //TODO: domain
}

impl<E: Pairing> VerifierKey<E::ScalarField, KZG<E>> {
    pub fn from_ring_and_kzg_vk<P: AffineRepr<BaseField = E::ScalarField>>(
        ring: &Ring<E::ScalarField, E, P>,
        kzg_vk: RawKzgVerifierKey<E>,
    ) -> Self {
        Self::from_commitment_and_kzg_vk(FixedColumnsCommitted::from_ring(ring), kzg_vk)
    }

    pub fn from_commitment_and_kzg_vk(
        commitment: FixedColumnsCommitted<E::ScalarField, KzgCommitment<E>>,
        kzg_vk: RawKzgVerifierKey<E>,
    ) -> Self {
        Self {
            pcs_raw_vk: kzg_vk,
            fixed_columns_committed: commitment,
        }
    }

    pub fn commitment(&self) -> FixedColumnsCommitted<E::ScalarField, KzgCommitment<E>> {
        self.fixed_columns_committed.clone()
    }
}

pub fn index<F: PrimeField, CS: PCS<F>, P: AffineRepr<BaseField = F>>(
    pcs_params: &CS::Params,
    piop_params: &PiopParams<F, P>,
    keys: &[P],
) -> (ProverKey<F, CS, P>, VerifierKey<F, CS>) {
    let pcs_ck = pcs_params.ck();
    let pcs_raw_vk = pcs_params.raw_vk();
    let fixed_columns = piop_params.fixed_columns(keys);
    let fixed_columns_committed = fixed_columns.commit::<CS>(&pcs_ck);
    let verifier_key = VerifierKey {
        pcs_raw_vk: pcs_raw_vk.clone(),
        fixed_columns_committed: fixed_columns_committed.clone(),
    };
    let prover_key = ProverKey {
        pcs_ck,
        fixed_columns,
        verifier_key,
    };
    let verifier_key = VerifierKey {
        pcs_raw_vk,
        fixed_columns_committed,
    };
    (prover_key, verifier_key)
}
