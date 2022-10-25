use std::marker::PhantomData;

use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::PrimeField;
use ark_poly::Evaluations;
use ark_poly::univariate::DensePolynomial;
use ark_std::{test_rng, UniformRand};
use fflonk::pcs::Commitment;

use common::{Column, FieldColumn};
use common::domain::Domain;
use common::gadgets::booleanity::{BitColumn, Booleanity};
use common::gadgets::fixed_cells::FixedCells;
use common::gadgets::inner_prod::InnerProd;
use common::gadgets::ProverGadget;
use common::gadgets::sw_cond_add::{AffineColumn, CondAdd};
use common::piop::ProverPiop;

use crate::piop::{FixedColumns, RingCommitments, RingEvaluations};
use crate::piop::params::PiopParams;

// The 'table': columns representing the execution trace of the computation
// and the constraints -- polynomials that vanish on every 2 consecutive rows.
pub struct PiopProver<F: PrimeField, Curve: SWCurveConfig<BaseField=F>> {
    domain: Domain<F>,
    bits: BitColumn<F>,
    points: AffineColumn<F, Affine<Curve>>,
    selector: FieldColumn<F>,
    inner_prod: InnerProd<F>,
    cond_add: CondAdd<F, Affine<Curve>>,
    booleanity: Booleanity<F>,
    cond_add_acc_x: FixedCells<F>,
    cond_add_acc_y: FixedCells<F>,
    inner_prod_acc: FixedCells<F>,
}


impl<F: PrimeField, Curve: SWCurveConfig<BaseField=F>> PiopProver<F, Curve>
{
    pub fn build(params: &PiopParams<F, Curve>,
                 fixed_columns: &FixedColumns<F, Curve>,
                 prover_index_in_keys: usize,
                 secret: Curve::ScalarField) -> Self {
        let domain = params.domain.clone();
        let bits = params.bits_column(prover_index_in_keys, secret);
        let inner_prod = InnerProd::init(fixed_columns.selector.clone(), bits.col.clone(), &domain);
        let cond_add = CondAdd::init(bits.clone(), fixed_columns.points.clone(), &domain);
        let booleanity = Booleanity::init(bits.clone());
        let cond_add_acc_x = FixedCells::init(cond_add.acc.xs.clone(), &domain);
        let cond_add_acc_y = FixedCells::init(cond_add.acc.ys.clone(), &domain);
        let inner_prod_acc = FixedCells::init(inner_prod.acc.clone(), &domain);
        Self {
            domain,
            bits,
            points: fixed_columns.points.clone(),
            selector: fixed_columns.selector.clone(),
            inner_prod,
            cond_add,
            booleanity,
            cond_add_acc_x,
            cond_add_acc_y,
            inner_prod_acc,
        }
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
    type Instance = Affine<Curve>;

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
            self.selector.as_poly().clone(),
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
        let selector = self.selector.evaluate(zeta);
        Self::Evaluations {
            selector,
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
            self.cond_add_acc_x.constraints(),
            self.cond_add_acc_y.constraints(),
            self.inner_prod_acc.constraints(),
        ].concat()
    }

    fn constraints_lin(&self, zeta: &F) -> Vec<DensePolynomial<F>> {
        vec![
            self.inner_prod.constraints_linearized(zeta),
            self.cond_add.constraints_linearized(zeta),
            self.booleanity.constraints_linearized(zeta),
            self.cond_add_acc_x.constraints_linearized(zeta),
            self.cond_add_acc_y.constraints_linearized(zeta),
            self.inner_prod_acc.constraints_linearized(zeta),
        ].concat()
    }

    fn domain(&self) -> &Domain<F> {
        &self.domain
    }

    fn result(&self) -> Self::Instance {
        self.cond_add.result
    }
}
