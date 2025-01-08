use ark_ec::AffineRepr;
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::Evaluations;
use ark_std::marker::PhantomData;
use ark_std::{vec, vec::Vec};
use ark_std::cmp::max;
use fflonk::pcs::Commitment;

use common::domain::Domain;
use common::gadgets::booleanity::{BitColumn, Booleanity};
use common::gadgets::cond_add::CondAdd;
use common::gadgets::fixed_cells::FixedCells;
use common::gadgets::inner_prod::InnerProd;
use common::gadgets::ProverGadget;
use common::piop::ProverPiop;
use common::{AffineColumn, Column, FieldColumn};

use crate::piop::params::PiopParams;
use crate::piop::FixedColumns;
use crate::piop::{RingCommitments, RingEvaluations};

// The 'table': columns representing the execution trace of the computation
// and the constraints -- polynomials that vanish on every 2 consecutive rows.
pub struct PiopProver<F: PrimeField, P: AffineRepr<BaseField = F>, CondAddT: CondAdd<F, P>> {
    domain: Domain<F>,
    // Fixed (public input) columns:
    pubkey_points: AffineColumn<F, P>,    // Private input column.
    signer_index: BitColumn<F>,
    // Gadgets:
    booleanity: Booleanity<F>, //this to prove the bit column is actually holding bits 
    //inner_prod: InnerProd<F>, TODO chcek! skalman: we do not have a ring selector anymore. We alreday kan check that cond add result in Pubkey but we don't know which pub key but is one of the pubkeys. So I assume we don't need the inner product test of k.<pubkeys>
    //inner_prod_acc: FixedCells<F>,
 
    signer_secret_key_bits: BitColumn<F>,

    cond_add_pubkey: CondAddT,
    cond_add_pubkey_acc_x: FixedCells<F>,
    cond_add_pubkey_acc_y: FixedCells<F>,

    cond_add_vrfout: CondAddT,
    cond_add_vrfout_acc_x: FixedCells<F>,
    cond_add_vrfout_acc_y: FixedCells<F>,

}

impl<F: PrimeField, P: AffineRepr<BaseField = F>, CondAddT: CondAdd<F, P>>
    PiopProver<F, P, CondAddT>
{
    pub fn build(
        params: &PiopParams<F, P>,
        fixed_columns: FixedColumns<F, P>,
        prover_index_in_keys: usize,
        secret: P::ScalarField,
    ) -> Self {
        let domain = params.domain.clone();
        let FixedColumns {
            pubkey_points,
            power_of_2_multiples_of_gen,
        } = fixed_columns;
        let bits = Self::bits_columns(params, prover_index_in_keys, secret);
        let cond_add_pubkey = CondAddT::init(bits.clone(), pubkey_points.clone(), params.seed, &domain);
        let booleanity = Booleanity::init(bits.clone());
        let pubkey_acc = cond_add_pubkey.get_acc();
        let cond_add_acc_x = FixedCells::init(pubkey_acc.xs.clone(), &domain);
        let cond_add_acc_y = FixedCells::init(pubkey_acc.ys.clone(), &domain);
        Self {
            domain,
            pubkey_points,
            bits,
            cond_add_acc_x,
            cond_add_acc_y,
            booleanity,
            cond_add_pubkey,
        }
    }

    fn bits_columns(
        params: &PiopParams<F, P>,
        index_in_keys: usize,
        secret: P::ScalarField,
    ) -> (BitColumn<F>, BitColumn<F>) {
        let mut signer_selector_bit_vector = vec![false; params.ring_size];
        signer_selector_bit_vector[index_in_keys] = true;
        let secret_key_bit_vector = params.scalar_part(secret);
        assert_eq!(max(signer_selector_bit_vector, secret_key_bit_vector.len()), params.domain.capacity - 1);
        (BitColumn::init(signer_selector_bit_vector, &params.domain),
         BitColumn::init(secret_key_bit_vector, &params.domain))
    }
}

impl<F, C, P, CondAddT> ProverPiop<F, C> for PiopProver<F, P, CondAddT>
where
    F: PrimeField,
    C: Commitment<F>,
    P: AffineRepr<BaseField = F>,
    CondAddT: CondAdd<F, P>,
{
    type Commitments = RingCommitments<F, C>;
    type Evaluations = RingEvaluations<F>;
    type Instance = P;

    fn committed_columns<Fun: Fn(&DensePolynomial<F>) -> C>(
        &self,
        commit: Fun,
    ) -> Self::Commitments {
        let bits = commit(self.bits.as_poly());
        let cond_add_acc = [
            commit(self.cond_add.get_acc().xs.as_poly()),
            commit(self.cond_add.get_acc().ys.as_poly()),
        ];
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
        let points = [self.points.xs.evaluate(zeta), self.points.ys.evaluate(zeta)];
        let ring_selector = self.ring_selector.evaluate(zeta);
        let bits = self.bits.evaluate(zeta);
        let inn_prod_acc = self.inner_prod.acc.evaluate(zeta);
        let cond_add_acc = [
            self.cond_add.get_acc().xs.evaluate(zeta),
            self.cond_add.get_acc().ys.evaluate(zeta),
        ];
        Self::Evaluations {
            points,
            ring_selector,
            bits,
            inn_prod_acc,
            cond_add_acc,
        }
    }

    fn constraints(&self) -> Vec<Evaluations<F>> {
        [
            self.inner_prod.constraints(),
            self.cond_add.constraints(),
            self.booleanity.constraints(),
            self.cond_add_acc_x.constraints(),
            self.cond_add_acc_y.constraints(),
            self.inner_prod_acc.constraints(),
        ]
        .concat()
    }

    fn constraints_lin(&self, zeta: &F) -> Vec<DensePolynomial<F>> {
        [
            self.inner_prod.constraints_linearized(zeta),
            self.cond_add.constraints_linearized(zeta),
            self.booleanity.constraints_linearized(zeta),
            self.cond_add_acc_x.constraints_linearized(zeta),
            self.cond_add_acc_y.constraints_linearized(zeta),
            self.inner_prod_acc.constraints_linearized(zeta),
        ]
        .concat()
    }

    fn domain(&self) -> &Domain<F> {
        &self.domain
    }

    fn result(&self) -> Self::Instance {
        self.cond_add.get_result()
    }
}
