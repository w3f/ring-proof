use ark_ec::{AffineRepr, CurveConfig};
use ark_ff::PrimeField;
use ark_poly::Evaluations;
use ark_poly::univariate::DensePolynomial;
use ark_std::marker::PhantomData;
use ark_std::{vec, vec::Vec};
use fflonk::pcs::Commitment;

use common::{Column, FieldColumn};
use common::domain::Domain;
use common::gadgets::booleanity::{BitColumn, Booleanity};
use common::gadgets::fixed_cells::FixedCells;
use common::gadgets::inner_prod::InnerProd;
use common::gadgets::ProverGadget;
use common::gadgets::cond_add::{AffineColumn, CondAdd};
use common::piop::ProverPiop;

use crate::piop::{RingCommitments, RingEvaluations};
use crate::piop::FixedColumns;
use crate::piop::params::PiopParams;

// The 'table': columns representing the execution trace of the computation
// and the constraints -- polynomials that vanish on every 2 consecutive rows.
pub struct PiopProver<F: PrimeField, P: AffineRepr<BaseField=F>, CondAddT: CondAdd<F,P>> {
    domain: Domain<F>,
    // Fixed (public input) columns:
    points: AffineColumn<F, P>,
    ring_selector: FieldColumn<F>,
    // Private input column.
    bits: BitColumn<F>,
    // Gadgets:
    booleanity: Booleanity<F>,
    inner_prod: InnerProd<F>,
    inner_prod_acc: FixedCells<F>,
    cond_add: CondAddT,
    cond_add_acc_x: FixedCells<F>,
    cond_add_acc_y: FixedCells<F>,
}

impl<F: PrimeField, P: AffineRepr<BaseField=F>, CondAddT: CondAdd<F, P>> PiopProver<F, P, CondAddT>
{
    pub fn build(params: &PiopParams<F, P>,
                 fixed_columns: FixedColumns<F, P>,
                 prover_index_in_keys: usize,
                 secret: P::ScalarField) -> Self {
        let domain = params.domain.clone();
        let FixedColumns { points, ring_selector } = fixed_columns;
        let bits = Self::bits_column(&params, prover_index_in_keys, secret);
        let inner_prod = InnerProd::init(ring_selector.clone(), bits.col.clone(), &domain);
        let cond_add = CondAddT::init(bits.clone(), points.clone(), params.seed, &domain);
        let booleanity = Booleanity::init(bits.clone());
        let acc = cond_add.get_acc();
        let cond_add_acc_x = FixedCells::init(acc.xs.clone(), &domain);
        let cond_add_acc_y = FixedCells::init(acc.ys.clone(), &domain);
        let inner_prod_acc = FixedCells::init(inner_prod.acc.clone(), &domain);
        Self {
            domain,
            points,
            ring_selector,
            bits,
            inner_prod_acc,
            cond_add_acc_x,
            cond_add_acc_y,
            booleanity,
            inner_prod,
            cond_add,
        }
    }

    fn bits_column(params: &PiopParams<F, P>, index_in_keys: usize, secret: P::ScalarField) -> BitColumn<F> {
        let mut keyset_part = vec![false; params.keyset_part_size];
        keyset_part[index_in_keys] = true;
        let scalar_part = params.scalar_part(secret);
        let bits = [
            keyset_part,
            scalar_part
        ].concat();
        assert_eq!(bits.len(), params.domain.capacity - 1);
        BitColumn::init(bits, &params.domain)
    }
}

impl<F, C, P, CondAddT> ProverPiop<F, C> for PiopProver<F, P, CondAddT>
    where
        F: PrimeField,
        C: Commitment<F>,
    P: AffineRepr<BaseField=F>,
    CondAddT: CondAdd<F, P> + ProverGadget<F>,

{
    type Commitments = RingCommitments<F, C>;
    type Evaluations = RingEvaluations<F>;
    type Instance = P;

    fn committed_columns<Fun: Fn(&DensePolynomial<F>) -> C>(&self, commit: Fun) -> Self::Commitments {
        let bits = commit(self.bits.as_poly());
        let cond_add_acc = super::ArrayWrap([
            commit(self.cond_add.get_acc().xs.as_poly()),
            commit(self.cond_add.get_acc().ys.as_poly())
        ]);
        let inn_prod_acc = commit(self.inner_prod.acc.as_poly());
        Self::Commitments {
            bits,
            cond_add_acc,
            inn_prod_acc,
            phantom: PhantomData,
        }
    }

    // Should return polynomials in the consistent with
    // Self::Evaluations::to_vec() and Self::Commitments::to_vec().
    fn columns(&self) -> Vec<DensePolynomial<F>> {
        vec![
            self.points.xs.as_poly().clone(),
            self.points.ys.as_poly().clone(),
            self.ring_selector.as_poly().clone(),
            self.bits.as_poly().clone(),
            self.inner_prod.acc.as_poly().clone(),
            self.cond_add.get_acc().xs.as_poly().clone(),
            self.cond_add.get_acc().ys.as_poly().clone(),
        ]
    }

    fn columns_evaluated(&self, zeta: &F) -> Self::Evaluations {
        let points = super::ArrayWrap([
            self.points.xs.evaluate(zeta),
            self.points.ys.evaluate(zeta),
        ]);
        let ring_selector = self.ring_selector.evaluate(zeta);
        let bits = self.bits.evaluate(zeta);
        let inn_prod_acc = self.inner_prod.acc.evaluate(zeta);
        let cond_add_acc = super::ArrayWrap([
            self.cond_add.get_acc().xs.evaluate(zeta),
            self.cond_add.get_acc().ys.evaluate(zeta),
        ]);
        Self::Evaluations {
            points,
            ring_selector,
            bits,
            inn_prod_acc,
            cond_add_acc,
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
        self.cond_add.get_result()
    }
}
