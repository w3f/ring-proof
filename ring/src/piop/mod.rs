use std::marker::PhantomData;

use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use fflonk::pcs::Commitment;

use common::{ColumnsCommited, ColumnsEvaluated};
pub(crate) use prover::PiopProver;
pub(crate) use verifier::PiopVerifier;

mod prover;
mod verifier;
pub mod params;

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct RingCommitments<F: PrimeField, C: Commitment<F>> {
    pub(crate) bits: C,
    pub(crate) cond_add_acc: [C; 2],
    pub(crate) inn_prod_acc: C,
    pub(crate) phantom: PhantomData<F>,
}

impl<F: PrimeField, C: Commitment<F>> ColumnsCommited<F, C> for RingCommitments<F, C> {
    fn to_vec(self) -> Vec<C> {
        vec![
            self.bits,
            self.cond_add_acc[0].clone(),
            self.cond_add_acc[1].clone(),
            self.inn_prod_acc,
        ]
    }
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct RingEvaluations<F: PrimeField> {
    pub(crate) points: [F; 2],
    pub(crate) bits: F,
    pub(crate) cond_add_acc: [F; 2],
    pub(crate) inn_prod_acc: F,
}

impl<F: PrimeField> ColumnsEvaluated<F> for RingEvaluations<F> {
    fn to_vec(self) -> Vec<F> {
        vec![
            self.points[0],
            self.points[1],
            self.bits,
            self.cond_add_acc[0],
            self.cond_add_acc[1],
            self.inn_prod_acc,
        ]
    }
}
