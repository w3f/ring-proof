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
    signer_index: BitColumn<F>,
    ring_selector: BitColumn<F>, //all one vector
    
    pubkey_points: AffineColumn<F, P>,    // Private input column.
    signer_secret_key_bits: BitColumn<F>,
    power_of_2_multiples_of_gen: AffineColumn<F, P>,
    power_of_2_multiples_of_vrf_in: AffineColumn<F, P>,
    // Gadgets:
    booleanity_of_signer_index: Booleanity<F>, //this to prove the bit column is actually holding bits
    booleanity_of_secret_key_bits: Booleanity<F>,

    sole_signer_inner_prod: InnerProd<F>, //This is a binary cond add making sure \sum signer_index[i] = 1
    sole_signer_inner_prod_acc: FixedCells<F>,
    
    cond_add_pubkey: CondAddT,
    cond_add_pubkey_acc_x: FixedCells<F>,
    cond_add_pubkey_acc_y: FixedCells<F>,

    cond_add_gen_multiples: CondAddT,
    cond_add_gen_multiples_acc_x: FixedCells<F>,
    cond_add_gen_multiples_acc_y: FixedCells<F>,

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
        vrf_input: P,
        
    ) -> Self {
        let domain = params.domain.clone();
        let FixedColumns {
            pubkey_points,
            power_of_2_multiples_of_gen,
            ring_selector,
        } = fixed_columns;
        let (signer_index, signer_secret_key_bits, ring_selector) = Self::bits_columns(params, prover_index_in_keys, secret); 
        let sole_signer_inner_prod = InnerProd::init(ring_selector.col.clone(), signer_index.col.clone(), &domain);
        
        let cond_add_pubkey = CondAddT::init(signer_index.clone(), pubkey_points.clone(), params.seed, &domain);
        
        let booleanity_of_signer_index = Booleanity::init(signer_index.clone());
        let booleanity_of_secret_key_bits = Booleanity::init(signer_secret_key_bits.clone());

        let sole_signer_inner_prod_acc = FixedCells::init(sole_signer_inner_prod.acc.clone(), &domain);

        let pubkey_acc = cond_add_pubkey.get_acc();
        let cond_add_pubkey_acc_x = FixedCells::init(pubkey_acc.xs.clone(), &domain);
        let cond_add_pubkey_acc_y = FixedCells::init(pubkey_acc.ys.clone(), &domain);

        //we also need to add the accumulation of generator multiples in respect to to secret key's bits
        let cond_add_gen_multiples = CondAddT::init(signer_secret_key_bits.clone(), power_of_2_multiples_of_gen.clone(), params.seed, &domain);
        let gen_multiples_acc = cond_add_gen_multiples.get_acc();
        let cond_add_gen_multiples_acc_x = FixedCells::init(gen_multiples_acc.xs.clone(), &domain);
        let cond_add_gen_multiples_acc_y = FixedCells::init(gen_multiples_acc.ys.clone(), &domain);

        //And also we need compute VRFout in the snark and verify it has been generated using the same secret key.
	let power_of_2_multiples_of_vrf_in = PiopParams::power_of_2_multiples_of(vrf_input, params.scalar_bitlen);
	let power_of_2_multiples_of_vrf_in = AffineColumn::public_column(power_of_2_multiples_of_vrf_in, &params.domain);

        
        let cond_add_vrfout = CondAddT::init(signer_secret_key_bits.clone(), power_of_2_multiples_of_vrf_in.clone(), params.seed, &domain);
        let vrfout_acc = cond_add_vrfout.get_acc();
        let cond_add_vrfout_acc_x = FixedCells::init(vrfout_acc.xs.clone(), &domain);
        let cond_add_vrfout_acc_y = FixedCells::init(vrfout_acc.ys.clone(), &domain);

        
        Self {
            domain,
            pubkey_points,
            signer_index,
            ring_selector,
            signer_secret_key_bits,            

            power_of_2_multiples_of_gen,
            power_of_2_multiples_of_vrf_in,
            
            booleanity_of_signer_index,
            booleanity_of_secret_key_bits,

            sole_signer_inner_prod,
            sole_signer_inner_prod_acc,
            
            cond_add_pubkey_acc_x,
            cond_add_pubkey_acc_y,
            cond_add_pubkey,

            cond_add_gen_multiples_acc_x,
            cond_add_gen_multiples_acc_y,
            cond_add_gen_multiples,

            cond_add_vrfout_acc_x,
            cond_add_vrfout_acc_y,
            cond_add_vrfout,
        }
    }

    fn bits_columns(
        params: &PiopParams<F, P>,
        index_in_keys: usize,
        secret: P::ScalarField,
    ) -> (BitColumn<F>, BitColumn<F>, BitColumn<F>) {
        //TODO: this is wrong it should be number of secret key bits.
        let mut signer_selector_bit_vector = vec![false; params.padded_keyset_size];
        signer_selector_bit_vector[index_in_keys] = true;

        let select_all_bit_vector = vec![true; params.padded_keyset_size];

        let secret_key_bit_vector = params.scalar_part(secret);
        assert_eq!(max(signer_selector_bit_vector.len(), secret_key_bit_vector.len()), params.domain.capacity - 1);
        (BitColumn::init(signer_selector_bit_vector, &params.domain),
         BitColumn::init(secret_key_bit_vector, &params.domain),
         BitColumn::init(select_all_bit_vector, &params.domain),
        )
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
        let signer_index = commit(self.signer_index.as_poly());
        let signer_secret_key_bits = commit(self.signer_secret_key_bits.as_poly());
        let ring_selector = commit(self.ring_selector.as_poly());

        let cond_add_pubkey_acc = [
            commit(self.cond_add_pubkey.get_acc().xs.as_poly()),
            commit(self.cond_add_pubkey.get_acc().ys.as_poly()),
        ];

        let sole_signer_inn_prod_acc = commit(self.sole_signer_inner_prod.acc.as_poly());
        let cond_add_gen_multiples_acc = [
            commit(self.cond_add_gen_multiples.get_acc().xs.as_poly()),
            commit(self.cond_add_gen_multiples.get_acc().ys.as_poly()),
        ];

        let cond_add_vrfout_acc = [
            commit(self.cond_add_vrfout.get_acc().xs.as_poly()),
            commit(self.cond_add_vrfout.get_acc().ys.as_poly()),
        ];

        Self::Commitments {
            signer_index,
            signer_secret_key_bits,
            ring_selector,
            sole_signer_inn_prod_acc,
            cond_add_pubkey_acc,
            cond_add_gen_multiples_acc,
            cond_add_vrfout_acc,
            phantom: PhantomData,
        }
    }

    // Should return polynomials in the consistent with
    // Self::Evaluations::to_vec() and Self::Commitments::to_vec().
    fn columns(&self) -> Vec<DensePolynomial<F>> {
        vec![
            self.pubkey_points.xs.as_poly().clone(),
            self.pubkey_points.ys.as_poly().clone(),
            self.ring_selector.as_poly().clone(),
            self.signer_index.as_poly().clone(),
            self.signer_secret_key_bits.as_poly().clone(),
            
            self.sole_signer_inner_prod.acc.as_poly().clone(),

            self.cond_add_pubkey.get_acc().xs.as_poly().clone(),
            self.cond_add_pubkey.get_acc().ys.as_poly().clone(),

            self.cond_add_gen_multiples.get_acc().xs.as_poly().clone(),
            self.cond_add_gen_multiples.get_acc().ys.as_poly().clone(),

            self.cond_add_vrfout.get_acc().xs.as_poly().clone(),
            self.cond_add_vrfout.get_acc().ys.as_poly().clone(),
        ]
    }

    fn columns_evaluated(&self, zeta: &F) -> Self::Evaluations {
        let pubkey_points = [self.pubkey_points.xs.evaluate(zeta), self.pubkey_points.ys.evaluate(zeta)];
        let ring_selector = self.ring_selector.evaluate(zeta);
        let signer_index = self.signer_index.evaluate(zeta);
        let signer_secret_key_bits = self.signer_secret_key_bits.evaluate(zeta);

        let sole_signer_inn_prod_acc = self.sole_signer_inner_prod.acc.evaluate(zeta);

        let cond_add_pubkey_acc = [
            self.cond_add_pubkey.get_acc().xs.evaluate(zeta),
            self.cond_add_pubkey.get_acc().ys.evaluate(zeta),
        ];

        let cond_add_gen_multiples_acc = [
            self.cond_add_gen_multiples.get_acc().xs.evaluate(zeta),
            self.cond_add_gen_multiples.get_acc().ys.evaluate(zeta),
        ];

        let cond_add_vrfout_acc = [
            self.cond_add_vrfout.get_acc().xs.evaluate(zeta),
            self.cond_add_vrfout.get_acc().ys.evaluate(zeta),
        ];

        Self::Evaluations {
            pubkey_points,
            ring_selector,
            signer_index,
            signer_secret_key_bits,
            sole_signer_inn_prod_acc,
            cond_add_pubkey_acc,
            cond_add_gen_multiples_acc,
            cond_add_vrfout_acc,
        }
    }

    fn constraints(&self) -> Vec<Evaluations<F>> {
        [
            self.sole_signer_inner_prod.constraints(),
            
            self.cond_add_pubkey.constraints(),
            self.cond_add_gen_multiples.constraints(),
            self.cond_add_vrfout.constraints(),
            
            self.booleanity_of_signer_index.constraints(),
            self.booleanity_of_secret_key_bits.constraints(),
            
            self.cond_add_pubkey_acc_x.constraints(),
            self.cond_add_pubkey_acc_y.constraints(),
            
            self.cond_add_gen_multiples_acc_x.constraints(),
            self.cond_add_gen_multiples_acc_y.constraints(),

            self.cond_add_vrfout_acc_x.constraints(),
            self.cond_add_vrfout_acc_y.constraints(),

            self.sole_signer_inner_prod_acc.constraints(),
        ]
        .concat()
    }

    fn constraints_lin(&self, zeta: &F) -> Vec<DensePolynomial<F>> {
        [
            self.sole_signer_inner_prod.constraints_linearized(zeta),

            self.cond_add_pubkey.constraints_linearized(zeta),
            self.cond_add_gen_multiples.constraints_linearized(zeta),
            self.cond_add_vrfout.constraints_linearized(zeta),

            self.booleanity_of_signer_index.constraints_linearized(zeta),
            self.booleanity_of_secret_key_bits.constraints_linearized(zeta),
            
            self.cond_add_pubkey_acc_x.constraints_linearized(zeta),
            self.cond_add_pubkey_acc_y.constraints_linearized(zeta),

            self.cond_add_gen_multiples_acc_x.constraints_linearized(zeta),
            self.cond_add_gen_multiples_acc_y.constraints_linearized(zeta),

            self.cond_add_vrfout_acc_x.constraints_linearized(zeta),
            self.cond_add_vrfout_acc_y.constraints_linearized(zeta),


            self.sole_signer_inner_prod_acc.constraints_linearized(zeta),
        ]
        .concat()
    }

    fn domain(&self) -> &Domain<F> {
        &self.domain
    }

    fn result(&self) -> Self::Instance {
        self.cond_add_pubkey.get_result()
    }
}
