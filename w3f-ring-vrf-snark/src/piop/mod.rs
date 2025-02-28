use ark_ec::pairing::Pairing;
use ark_ec::twisted_edwards::{Affine, TECurveConfig};
use ark_ec::AffineRepr;
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::marker::PhantomData;
use ark_std::{vec, vec::Vec};
use w3f_pcs::pcs::kzg::commitment::KzgCommitment;
use w3f_pcs::pcs::kzg::params::RawKzgVerifierKey;
use w3f_pcs::pcs::kzg::KZG;
use w3f_pcs::pcs::{Commitment, PcsParams, PCS};

pub(crate) use prover::PiopProver;
pub(crate) use verifier::PiopVerifier;
use w3f_plonk_common::gadgets::ec::AffineColumn;
use w3f_plonk_common::{Column, ColumnsCommited, ColumnsEvaluated, FieldColumn};

use crate::ring_vrf::Ring;
use crate::PiopParams;

pub mod params;
mod prover;
mod verifier;

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct RingCommitments<F: PrimeField, C: Commitment<F>> {
    // doublings_of_g are prepented by the plonk verifier
    pub(crate) sk_bits: C,
    pub(crate) pk_from_sk: [C; 2],
    pub(crate) doublings_of_vrf_in: [C; 2],
    // pub(crate) signer_index: C,
    // pub(crate) ring_selector: C,
    // pub(crate) sole_signer_inn_prod_acc: C,
    // pub(crate) cond_add_pubkey_acc: [C; 2],
    // pub(crate) vrf_out_acc: [C; 2],
    pub(crate) phantom: PhantomData<F>,
}

impl<F: PrimeField, C: Commitment<F>> ColumnsCommited<F, C> for RingCommitments<F, C> {
    fn to_vec(self) -> Vec<C> {
        vec![
            // doublings_of_g are prepented by the plonk verifier
            self.sk_bits,
            self.pk_from_sk[0].clone(),
            self.pk_from_sk[1].clone(),
            // self.signer_index,
            // self.sole_signer_inn_prod_acc,
            // self.cond_add_pubkey_acc[0].clone(),
            // self.cond_add_pubkey_acc[1].clone(),
            // self.vrf_out_acc[0].clone(),
            // self.vrf_out_acc[1].clone(),
        ]
    }
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct RingEvaluations<F: PrimeField> {
    pub(crate) doublings_of_g: [F; 2],
    pub(crate) sk_bits: F,
    pub(crate) pk_from_sk: [F; 2],
    pub(crate) doublings_of_vrf_in: [F; 2],
    // pub(crate) pks: [F; 2],
    // pub(crate) ring_selector: F,
    // pub(crate) signer_index: F,
    // pub(crate) powers_of_in: [F; 2],
    // pub(crate) k_is_one_bit: F,
    // pub(crate) pk_from_k_acc: [F; 2],
    // pub(crate) vrf_out_acc: [F; 2],
}

impl<F: PrimeField> ColumnsEvaluated<F> for RingEvaluations<F> {
    fn to_vec(self) -> Vec<F> {
        vec![
            self.doublings_of_g[0],
            self.doublings_of_g[1],
            self.sk_bits,
            self.pk_from_sk[0],
            self.pk_from_sk[1],
            // self.pks[0],
            // self.pks[1],
            // self.ring_selector,
            // self.signer_index,
            // self.k_is_one_bit,
            // self.pk_from_k_acc[0],
            // self.pk_from_k_acc[1],
            // self.vrf_out_acc[0],
            // self.vrf_out_acc[1],
        ]
    }
}

// Columns commitment to which the verifier knows (or trusts).
// TODO: comments
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct FixedColumns<F: PrimeField, G: AffineRepr<BaseField = F>> {
    // Public keys of the ring participants in order,
    // followed by the powers-of-2 multiples of the second Pedersen base.
    // pk_1, ..., pk_n, H, 2H, 4H, ..., 2^sH
    // 1          n                     n+s+1
    // pubkey_points: AffineColumn<F, G>,
    // The powers-of-2 multiples of the prime subgroup generator.
    doublings_of_g: AffineColumn<F, G>,
    // Binary column that highlights which rows of the table correspond to the ring.
    // 1, 1, ..., 1, 0, 0, ..., 0
    // 1          n
    // ring_selector: FieldColumn<F>,
}

// Commitments to the fixed columns (see above).
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize, PartialEq, Eq, Debug)]
pub struct FixedColumnsCommitted<F: PrimeField, C: Commitment<F>> {
    pub doublings_of_g: [C; 2],
    // pub ring_selector: C,
    pub phantom: PhantomData<F>,
}

impl<F: PrimeField, C: Commitment<F>> FixedColumnsCommitted<F, C> {
    fn as_vec(&self) -> Vec<C> {
        vec![
            self.doublings_of_g[0].clone(),
            self.doublings_of_g[1].clone(),
            // self.ring_selector.clone(),
        ]
    }
}

impl<E: Pairing> FixedColumnsCommitted<E::ScalarField, KzgCommitment<E>> {
    pub fn from_ring<G: TECurveConfig<BaseField = E::ScalarField>>(
        ring: &Ring<E::ScalarField, E, G>,
    ) -> Self {
        let cx = KzgCommitment(ring.cx);
        let cy = KzgCommitment(ring.cy);
        Self {
            doublings_of_g: [cx, cy],
            // ring_selector: KzgCommitment(ring.selector),
            phantom: Default::default(),
        }
    }
}

impl<F: PrimeField, G: AffineRepr<BaseField = F>> FixedColumns<F, G> {
    fn commit<CS: PCS<F>>(&self, ck: &CS::CK) -> FixedColumnsCommitted<F, CS::C> {
        let points = [
            CS::commit(ck, self.doublings_of_g.xs.as_poly()).unwrap(),
            CS::commit(ck, self.doublings_of_g.ys.as_poly()).unwrap(),
        ];
        // let ring_selector = CS::commit(ck, self.ring_selector.as_poly()).unwrap();
        FixedColumnsCommitted {
            doublings_of_g: points,
            // ring_selector,
            phantom: Default::default(),
        }
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct ProverKey<F: PrimeField, CS: PCS<F>, G: AffineRepr<BaseField = F>> {
    pub(crate) pcs_ck: CS::CK,
    pub(crate) fixed_columns: FixedColumns<F, G>,
    pub(crate) verifier_key: VerifierKey<F, CS>, // used in the Fiat-Shamir transform
}

#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerifierKey<F: PrimeField, CS: PCS<F>> {
    pub(crate) pcs_raw_vk: <CS::Params as PcsParams>::RVK,
    pub(crate) fixed_columns_committed: FixedColumnsCommitted<F, CS::C>,
    //TODO: domain
}

impl<E: Pairing> VerifierKey<E::ScalarField, KZG<E>> {
    pub fn from_ring_and_kzg_vk<G: TECurveConfig<BaseField = E::ScalarField>>(
        ring: &Ring<E::ScalarField, E, G>,
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

pub fn index<F: PrimeField, CS: PCS<F>, Curve: TECurveConfig<BaseField = F>>(
    pcs_params: &CS::Params,
    piop_params: &PiopParams<F, Curve>,
    keys: &[Affine<Curve>],
) -> (ProverKey<F, CS, Affine<Curve>>, VerifierKey<F, CS>) {
    let pcs_ck = pcs_params.ck();
    let pcs_raw_vk = pcs_params.raw_vk();
    let fixed_columns = piop_params.fixed_columns(&keys);
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
