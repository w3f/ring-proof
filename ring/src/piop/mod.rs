use std::marker::PhantomData;
use ark_ec::AffineRepr;
use ark_ec::short_weierstrass::{Affine, SWCurveConfig};

use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use fflonk::pcs::{Commitment, PCS};

use common::{ColumnsCommited, ColumnsEvaluated, FieldColumn};
use common::gadgets::sw_cond_add::AffineColumn;
use common::setup::Setup;
pub(crate) use prover::PiopProver;
pub(crate) use verifier::PiopVerifier;
use crate::piop::params::PiopParams;

mod prover;
mod verifier;
pub mod params;
mod lagrange;

pub struct FixedColumns<F: PrimeField, Curve: SWCurveConfig<BaseField=F>> {
    selector: FieldColumn<F>,
    pub points: AffineColumn<F, Affine<Curve>>,
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct FixedColumnsCommitted<F: PrimeField, C: Commitment<F>> {
    selector_committed: C,
    pub points_committed: [C; 2],
    phantom: PhantomData<F>,
}

impl<F: PrimeField, C: Commitment<F>> ColumnsCommited<F, C> for FixedColumnsCommitted<F, C> {
    fn to_vec(self) -> Vec<C> {
        vec![
            self.selector_committed,
            self.points_committed[0].clone(),
            self.points_committed[1].clone(),
        ]
    }
}


pub struct Indexer<F: PrimeField, Curve: SWCurveConfig<BaseField=F>, CS: PCS<F>> {
    pub fixed_columns: FixedColumns<F, Curve>,
    pub fixed_columns_committed: FixedColumnsCommitted<F, CS::C>,
}

impl<F: PrimeField, Curve: SWCurveConfig<BaseField=F>, CS: PCS<F>> Indexer<F, Curve, CS> {
    pub fn new(setup: &Setup<F, CS>, params: &PiopParams<F, Curve>, keys: &[Affine<Curve>]) -> Self {
        let points = Self::points_column(params, keys);
        let keyset_part_selector = params.keyset_part_selector();
        let selector = params.domain.selector(keyset_part_selector);
        let selector_committed = setup.commit_to_column(&selector);
        let points_committed = [setup.commit_to_column(&points.xs), setup.commit_to_column(&points.ys)];
        let fixed_columns = FixedColumns { selector, points };
        let fixed_columns_committed = FixedColumnsCommitted { selector_committed, points_committed, phantom: PhantomData };
        Self {
            fixed_columns,
            fixed_columns_committed,
        }
    }

    fn points_column(params: &PiopParams<F, Curve>, keys: &[Affine<Curve>]) -> AffineColumn<F, Affine<Curve>> {
        assert!(keys.len() <= params.keyset_part_size);
        let padding_len = params.keyset_part_size - keys.len();
        let padding_point = Affine::<Curve>::generator(); //TODO: NUMS
        let padding = vec![padding_point; padding_len];
        let points = [
            keys,
            &padding,
            &params.power_of_2_multiples_of_h(),
        ].concat();
        assert_eq!(points.len(), params.domain.capacity - 1);
        AffineColumn::init(points, &params.domain)
    }
}


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
    pub(crate) selector: F,
    pub(crate) points: [F; 2],
    pub(crate) bits: F,
    pub(crate) cond_add_acc: [F; 2],
    pub(crate) inn_prod_acc: F,
}

impl<F: PrimeField> ColumnsEvaluated<F> for RingEvaluations<F> {
    fn to_vec(self) -> Vec<F> {
        vec![
            self.selector,
            self.points[0],
            self.points[1],
            self.bits,
            self.cond_add_acc[0],
            self.cond_add_acc[1],
            self.inn_prod_acc,
        ]
    }
}
