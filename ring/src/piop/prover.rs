use std::marker::PhantomData;

use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::PrimeField;
use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain};
use ark_poly::univariate::DensePolynomial;
use ark_std::{test_rng, UniformRand};
use fflonk::pcs::Commitment;

use common::Column;
use common::gadgets::booleanity::{BitColumn, Booleanity};
use common::gadgets::fixed_cells::FixedCells;
use common::gadgets::inner_prod::InnerProd;
use common::gadgets::ProverGadget;
use common::gadgets::sw_cond_add::{AffineColumn, CondAdd};
use common::piop::ProverPiop;

use crate::piop::{RingCommitments, RingEvaluations, SelectorColumns};
use crate::piop::params::PiopParams;


// The 'table': columns representing the execution trace of the computation
// and the constraints -- polynomials that vanish on every 2 consecutive rows.
pub struct PiopProver<F: PrimeField, Curve: SWCurveConfig<BaseField=F>> {
    bits: BitColumn<F>,
    points: AffineColumn<F, Affine<Curve>>,
    selectors: SelectorColumns<F>,
    inner_prod: InnerProd<F>,
    cond_add: CondAdd<F, Affine<Curve>>,
    booleanity: Booleanity<F>,
    fixed_cells_acc_x: FixedCells<F>,
    fixed_cells_acc_y: FixedCells<F>,
}


impl<F: PrimeField, Curve: SWCurveConfig<BaseField=F>> PiopProver<F, Curve>
{
    pub fn init(params: &PiopParams<F, Curve>,
                selectors: SelectorColumns<F>,
                points: AffineColumn<F, Affine<Curve>>,
                prover_index_in_keys: usize,
                secret: Curve::ScalarField) -> Self {
        let bits = Self::bits_column(params, prover_index_in_keys, secret);
        let inner_prod = InnerProd::init(selectors.ring_selector.clone(), bits.col.clone(), selectors.l_last.clone());
        let cond_add = CondAdd::init(bits.clone(), points.clone(), selectors.not_last.clone());
        let booleanity = Booleanity::init(bits.clone());
        let fixed_cells_acc_x = FixedCells::init(cond_add.acc.xs.clone(), selectors.l_first.clone(), selectors.l_last.clone());
        let fixed_cells_acc_y = FixedCells::init(cond_add.acc.ys.clone(), selectors.l_first.clone(), selectors.l_last.clone());
        Self {
            bits,
            points,
            selectors,
            inner_prod,
            cond_add,
            booleanity,
            fixed_cells_acc_x,
            fixed_cells_acc_y,
        }
    }

    pub fn keyset_column(params: &PiopParams<F, Curve>, keys: &[Affine<Curve>]) -> AffineColumn<F, Affine<Curve>> {
        assert!(keys.len() <= params.keyset_part_size);
        let padding_len = params.keyset_part_size - keys.len();
        let padding_point = Affine::<Curve>::rand(&mut test_rng()); //TODO: Ask Al
        let padding = vec![padding_point; padding_len];
        let points = [
            keys,
            &padding,
            &params.powers_of_h,
        ].concat();
        assert_eq!(points.len(), params.domain.size());
        AffineColumn::init(points)
    }

    fn bits_column(params: &PiopParams<F, Curve>, index_in_keys: usize, secret: Curve::ScalarField) -> BitColumn<F> {
        let mut keyset_part = vec![false; params.keyset_part_size];
        keyset_part[index_in_keys] = true;
        let scalar_part = params.scalar_part(secret);
        let bits = [
            keyset_part,
            scalar_part
        ].concat();
        BitColumn::init(bits)
    }
}

impl<F, C, Curve> ProverPiop<F, C> for PiopProver<F, Curve>
    where
        F: PrimeField,
        C: Commitment<F>,
        Curve: SWCurveConfig<BaseField=F>,

{
    type Commitments = RingCommitments<F, C>;
    type Evaluations = RingEvaluations<F>;
    type PublicInput = Affine<Curve>;

    fn committed_columns<Fun: Fn(&DensePolynomial<F>) -> C>(&self, commit: Fun) -> Self::Commitments {
        let bits = commit(self.bits.as_poly());
        let cond_add_acc = [
            commit(self.cond_add.acc.xs.as_poly()),
            commit(self.cond_add.acc.ys.as_poly())
        ];
        let inn_prod_acc = commit(self.inner_prod.acc.as_poly());
        Self::Commitments {
            bits,
            cond_add_acc,
            inn_prod_acc,
            phantom: PhantomData,
        }
    }

    fn columns(&self) -> Vec<DensePolynomial<F>> {
        vec![
            self.points.xs.as_poly().clone(),
            self.points.ys.as_poly().clone(),
            self.bits.as_poly().clone(),
            self.cond_add.acc.xs.as_poly().clone(),
            self.cond_add.acc.ys.as_poly().clone(),
            self.inner_prod.acc.as_poly().clone(),
        ]
    }

    fn columns_evaluated(&self, zeta: &F) -> Self::Evaluations {
        let bits = self.bits.evaluate(zeta);
        let cond_add_acc = [
            self.cond_add.acc.xs.evaluate(zeta),
            self.cond_add.acc.ys.evaluate(zeta),
        ];
        let inn_prod_acc = self.inner_prod.acc.evaluate(zeta);
        let points = [
            self.points.xs.evaluate(zeta),
            self.points.ys.evaluate(zeta),
        ];
        Self::Evaluations {
            bits,
            cond_add_acc,
            inn_prod_acc,
            points,
        }
    }

    fn constraints(&self) -> Vec<Evaluations<F>> {
        vec![
            self.inner_prod.constraints(),
            self.cond_add.constraints(),
            self.booleanity.constraints(),
            self.fixed_cells_acc_x.constraints(),
            self.fixed_cells_acc_y.constraints(),
        ].concat()
    }

    fn constraints_lin(&self, zeta: &F) -> Vec<DensePolynomial<F>> {
        vec![
            self.inner_prod.constraints_linearized(zeta),
            self.cond_add.constraints_linearized(zeta),
            self.booleanity.constraints_linearized(zeta),
            self.fixed_cells_acc_x.constraints_linearized(zeta),
            self.fixed_cells_acc_y.constraints_linearized(zeta),
        ].concat()
    }

    fn domain(&self) -> GeneralEvaluationDomain<F> {
        self.inner_prod.domain()
    }
}
