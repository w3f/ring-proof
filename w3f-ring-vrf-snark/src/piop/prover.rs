use ark_ec::twisted_edwards::{Affine, TECurveConfig};
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::Evaluations;
use ark_std::marker::PhantomData;
use ark_std::{vec, vec::Vec};
use w3f_pcs::pcs::Commitment;

use crate::piop::params::PiopParams;
use crate::piop::FixedColumns;
use crate::piop::{RingCommitments, RingEvaluations};
use w3f_plonk_common::domain::Domain;
use w3f_plonk_common::gadgets::booleanity::{BitColumn, Booleanity};
use w3f_plonk_common::gadgets::ec::AffineColumn;
use w3f_plonk_common::gadgets::ec::CondAdd;
use w3f_plonk_common::gadgets::fixed_cells::FixedCells;
use w3f_plonk_common::gadgets::inner_prod::InnerProd;
use w3f_plonk_common::gadgets::ProverGadget;
use w3f_plonk_common::piop::ProverPiop;

use w3f_plonk_common::Column;

/// The prover's private input is its secret key `sk`.
/// The public inputs are:
/// 1. the ring `R`, a vector of members' public keys,
/// 2. EC VRF input point `vrf_in`,
/// 3. EC VRF output point `vrf_out`.

/// The statement to prove is that:
/// 1. The prover is a member of the ring.
///    That means its public key `sk.G` is an element of `R`, or
///    `sk.G == R[k]` for some `k` witnessed by the prover, and the curve generator `G`.
/// 2. `vrf_out` is the valid prover's EC VRF pre-output for the input point `vrf_in`:
///    `vrf_out == sk.vrf_in`.

/// The `sk` is represented as a binary vector of length `L`. To address #1 we compute the prover's
/// public key `PK = sk.G`. As the point `G` is publicly known...

// 1. k is binary of norm 1
// 2. sk is binary
// 3. PK_1 := <sk, 2^i.G>
// 4. PK_2 := <k, R>
// 5. PK_1 == PK_2
// 6. 2^i.H is well formed
// 7. O := <sk, 2^i.H>
// TODO: selector

// The 'table': columns representing the execution trace of the computation
/// | 1              | 2           | 3           | 4              | 5              | 6                        | 7&8      | 9&10            | 11&12    | 13&14           |
/// | --------       | --------    | --------    | --             | -              | -                        | -        | -               | -        | -               |
/// | $k$            | $pk_x$      | $pk_y$      | $acc_{pk_x}$   | $acc_{pk_y}$   | $sk$                     | $2^iG$x2 | $acc_{sk}$x2    | $2^iH$x2 | $acc_{out}$x2   |
/// | --------       | --------    | --------    | --             | -              | -                        | -        | -               | -        | -               |
/// | signer's index | x of pubkey | y of pubkey | $\sum k_ipk_x$ | $\sum k_ipk_y$ | binary rep of secret key |          | $\sum sk_i2^iG$ |          | $\sum sk_i2^iH$ |
pub struct PiopProver<F: PrimeField, Curve: TECurveConfig<BaseField = F>> {
    domain: Domain<F>,
    // advice columns:
    // ring_selector: BitColumn<F>,
    doublings_of_g: AffineColumn<F, Affine<Curve>>,
    // pubkey_points: AffineColumn<F, Affine<Curve>>, // `R` above

    // private columns
    // signer_index: BitColumn<F>, // `k` above
    sk_bits: BitColumn<F>,

    // Gadgets:
    // booleanity_of_signer_index: Booleanity<F>, //this to prove the bit column is actually holding bits
    sk_bits_bool: Booleanity<F>,

    // sole_signer_inner_prod: InnerProd<F>, //This is a binary cond add making sure \sum signer_index[i] = 1
    // sole_signer_inner_prod_acc: FixedCells<F>,

    // pk_from_k: CondAdd<F, Affine<Curve>>,
    // pk_from_k_x: FixedCells<F>,
    // pk_from_k_y: FixedCells<F>,

    pk_from_sk: CondAdd<F, Affine<Curve>>,
    pk_from_sk_x: FixedCells<F>,
    pk_from_sk_y: FixedCells<F>,

    // powers_of_in: PowersOfTwoMultiplesTE<F, Curve>, // TODO; reconcile with addition
    //
    // vrf_out: CondAdd<F, Affine<Curve>>,
    // vrf_out_x: FixedCells<F>,
    // vrf_out_y: FixedCells<F>,
}


impl<F: PrimeField, Curve: TECurveConfig<BaseField = F>> PiopProver<F, Curve> {
    pub fn build(
        params: &PiopParams<F, Curve>,
        fixed_columns: FixedColumns<F, Affine<Curve>>, // TODO: rename to AdviceColumns
        prover_index_in_keys: usize,
        secret: Curve::ScalarField,
        vrf_input: Affine<Curve>,
    ) -> Self {
        let domain = params.domain.clone();
        // let pubkey_points = fixed_columns.pubkey_points;
        let powers_of_g = fixed_columns.doublings_of_g;
        // let (signer_index, signer_secret_key_bits, ring_selector) =
        //     Self::bits_columns(params, prover_index_in_keys, secret); // TODO: that's ugly
        let mut secret_key_bit_vector = params.scalar_part(secret);
        assert!(secret_key_bit_vector.len() + 1 <= domain.capacity);
        secret_key_bit_vector.resize(domain.capacity - 1, false);
        let signer_secret_key_bits = BitColumn::init(secret_key_bit_vector, &params.domain);

        // C1: `signer_index` column is boolean
        // let signer_index_is_bool = Booleanity::init(signer_index.clone());
        // C2 + С3: <signer_index, ring_selector> == 1
        // let scalar_inner_prod = InnerProd::init(ring_selector.col.clone(), signer_index.col.clone(), &domain);
        // let inner_prod_acc = FixedCells::init(scalar_inner_prod.acc.clone(), &domain);

        // C4: PK_k := <k, R>
        // let pk_from_k = CondAdd::init(
        //     signer_index.clone(),
        //     pubkey_points.clone(),
        //     params.seed,
        //     &domain,
        // );

        // C5: `sk` is boolean
        let booleanity_of_secret_key_bits = Booleanity::init(signer_secret_key_bits.clone());

        // C6: `PK_sk := sk.G`
        // assert_eq!(signer_secret_key_bits.bits.len(), powers_of_g.xs.vals().len(), "PIZDA");
        let pk_from_sk = CondAdd::init(
            signer_secret_key_bits.clone(),
            powers_of_g.clone(),
            params.seed,
            &domain,
        );

        // TODO: actually we want to prove that `PK_sk = PK_k` w/o exposing PK to the verifier
        // C7 + C8: `PK_sk == PK`
        let pk_from_sk_acc = pk_from_sk.acc.clone();
        let pk_from_sk_x = FixedCells::init(pk_from_sk_acc.xs.clone(), &domain);
        let pk_from_sk_y = FixedCells::init(pk_from_sk_acc.ys.clone(), &domain);
        // C9 + C10: `PK_k == PK`
        // let pk_from_k_acc = pk_from_k.acc.clone();
        // let pk_from_k_x = FixedCells::init(pk_from_k_acc.xs.clone(), &domain);
        // let pk_from_k_y = FixedCells::init(pk_from_k_acc.ys.clone(), &domain);

        // C11: 2-adic powers of the `vrf_input`
        // let powers_of_in = PowersOfTwoMultiplesTE::init(vrf_input, &domain);
        //
        // // C12: `snark_out := <sk, powers_of_in>`
        // let vrf_out = CondAdd::init(
        //     signer_secret_key_bits.clone(),
        //     powers_of_in.multiples.clone(),
        //     params.seed,
        //     &domain,
        // );
        // // C13+14: `snark_out == vrf_out`
        // let vrf_out_x = FixedCells::init( vrf_out.acc.xs.clone(), &domain);
        // let vrf_out_y = FixedCells::init( vrf_out.acc.ys.clone(), &domain);

        Self {
            domain,
            // pubkey_points,
            // signer_index,
            // ring_selector,
            sk_bits: signer_secret_key_bits,

            doublings_of_g: powers_of_g,

            // booleanity_of_signer_index: signer_index_is_bool,
            sk_bits_bool: booleanity_of_secret_key_bits,

            // sole_signer_inner_prod: scalar_inner_prod,
            // sole_signer_inner_prod_acc: inner_prod_acc,

            // pk_from_k_x,
            // pk_from_k_y,
            // pk_from_k,

            // powers_of_in,

            pk_from_sk_x,
            pk_from_sk_y,
            pk_from_sk,
            //
            // vrf_out_x,
            // vrf_out_y,
            // vrf_out,
        }
    }

    fn bits_columns(
        params: &PiopParams<F, Curve>,
        signer_index: usize,
        secret: Curve::ScalarField,
    ) -> (BitColumn<F>, BitColumn<F>, BitColumn<F>) {
        //TODO: this is wrong it should be number of secret key bits.
        let mut signer_selector_bit_vector = vec![false; params.keyset_part_size];
        signer_selector_bit_vector[signer_index] = true;

        let select_all_bit_vector = vec![true; params.keyset_part_size];

        let secret_key_bit_vector = params.scalar_part(secret);
        // assert_eq!(max(signer_selector_bit_vector.len(), secret_key_bit_vector.len()), params.domain.capacity - 1);
        (
            BitColumn::init(signer_selector_bit_vector, &params.domain),
            BitColumn::init(secret_key_bit_vector, &params.domain),
            BitColumn::init(select_all_bit_vector, &params.domain),
        )
    }
}

impl<F, C, Curve> ProverPiop<F, C> for PiopProver<F, Curve>
where
    F: PrimeField,
    C: Commitment<F>,
    Curve: TECurveConfig<BaseField = F>,
{
    type Commitments = RingCommitments<F, C>;
    type Evaluations = RingEvaluations<F>;
    type Instance = Affine<Curve>;

    fn committed_columns<Fun: Fn(&DensePolynomial<F>) -> C>(
        &self,
        commit: Fun,
    ) -> Self::Commitments {
        // let signer_index = commit(self.signer_index.as_poly());
        let signer_secret_key_bits = commit(self.sk_bits.as_poly());
        // let ring_selector = commit(self.ring_selector.as_poly());

        // let cond_add_pubkey_acc = [
        //     commit(self.pk_from_k.acc.xs.as_poly()),
        //     commit(self.pk_from_k.acc.ys.as_poly()),
        // ];

        // let sole_signer_inn_prod_acc = commit(self.sole_signer_inner_prod.acc.as_poly());

        let pk_from_sk_acc = [
            commit(self.pk_from_sk.acc.xs.as_poly()),
            commit(self.pk_from_sk.acc.ys.as_poly()),
        ];

        // let vrf_out_acc = [
        //     commit(self.vrf_out.acc.xs.as_poly()),
        //     commit(self.vrf_out.acc.ys.as_poly()),
        // ];

        RingCommitments {
            // signer_index,
            sk_bits: signer_secret_key_bits,
            // ring_selector,
            // sole_signer_inn_prod_acc,
            // cond_add_pubkey_acc,
            pk_from_sk: pk_from_sk_acc,
            // vrf_out_acc,
            phantom: PhantomData,

        }
    }

    // Should return polynomials in the consistent with
    // Self::Evaluations::to_vec() and Self::Commitments::to_vec().
    fn columns(&self) -> Vec<DensePolynomial<F>> {
        vec![
            self.sk_bits.as_poly().clone(),
            self.pk_from_sk.acc.xs.as_poly().clone(),
            self.pk_from_sk.acc.ys.as_poly().clone(),
            // self.pubkey_points.xs.as_poly().clone(),
            // self.pubkey_points.ys.as_poly().clone(),
            // self.ring_selector.as_poly().clone(),
            // self.signer_index.as_poly().clone(),
            // self.sole_signer_inner_prod.acc.as_poly().clone(),
            // self.pk_from_k.acc.xs.as_poly().clone(),
            // self.pk_from_k.acc.ys.as_poly().clone(),
            // self.vrf_out.acc.xs.as_poly().clone(),
            // self.vrf_out.acc.ys.as_poly().clone(),
        ]
    }

    fn columns_evaluated(&self, zeta: &F) -> RingEvaluations<F> {
        let doublings_of_g = [
            self.doublings_of_g.xs.evaluate(zeta),
            self.doublings_of_g.ys.evaluate(zeta),
        ];
        let sk_bits = self.sk_bits.evaluate(zeta);
        let pk_from_sk = [
            self.pk_from_sk.acc.xs.evaluate(zeta),
            self.pk_from_sk.acc.ys.evaluate(zeta),
        ];
        // let pubkey_points = [
        //     self.pubkey_points.xs.evaluate(zeta),
        //     self.pubkey_points.ys.evaluate(zeta),
        // ];
        // let ring_selector = self.ring_selector.evaluate(zeta);
        // let signer_index = self.signer_index.evaluate(zeta);


        // let sole_signer_inn_prod_acc = self.sole_signer_inner_prod.acc.evaluate(zeta);

        // let cond_add_pubkey_acc = [
        //     self.pk_from_k.acc.xs.evaluate(zeta),
        //     self.pk_from_k.acc.ys.evaluate(zeta),
        // ];



        // let cond_add_vrfout_acc = [
        //     self.vrf_out.acc.xs.evaluate(zeta),
        //     self.vrf_out.acc.ys.evaluate(zeta),
        // ];



        // let powers_of_in = [
        //     self.powers_of_in.multiples.xs.evaluate(zeta),
        //     self.powers_of_in.multiples.ys.evaluate(zeta),
        // ];

        // let pk_from_k_acc = [
        //     self.pk_from_k.acc.xs.evaluate(zeta),
        //     self.pk_from_k.acc.ys.evaluate(zeta),
        // ];

        RingEvaluations {
            doublings_of_g,
            sk_bits,
            pk_from_sk,
            // pks: pubkey_points,
            // ring_selector,
            // signer_index,
            // powers_of_in,
            // k_is_one_bit: sole_signer_inn_prod_acc,
            // vrf_out_acc: cond_add_vrfout_acc,
            // pk_from_k_acc,
        }
    }

    fn constraints(&self) -> Vec<Evaluations<F>> {
        [
            self.pk_from_sk.constraints(),
            self.sk_bits_bool.constraints(),
            // self.sole_signer_inner_prod.constraints(),
            // self.pk_from_k.constraints(),
            // self.vrf_out.constraints(),
            // self.booleanity_of_signer_index.constraints(),
            // self.pk_from_k_x.constraints(),
            // self.pk_from_k_y.constraints(),
            // self.pk_from_sk_x.constraints(),
            // self.pk_from_sk_y.constraints(),
            // self.vrf_out_x.constraints(),
            // self.vrf_out_y.constraints(),
            // self.sole_signer_inner_prod_acc.constraints(),
        ]
        .concat()
    }

    fn constraints_lin(&self, zeta: &F) -> Vec<DensePolynomial<F>> {
        [
            self.pk_from_sk.constraints_linearized(zeta),
            self.sk_bits_bool .constraints_linearized(zeta),
            // self.sole_signer_inner_prod.constraints_linearized(zeta),
            // self.pk_from_k.constraints_linearized(zeta),
            // self.vrf_out.constraints_linearized(zeta),
            // self.booleanity_of_signer_index.constraints_linearized(zeta),

            // self.pk_from_k_x.constraints_linearized(zeta),
            // self.pk_from_k_y.constraints_linearized(zeta),
            // self.pk_from_sk_x
            //     .constraints_linearized(zeta),
            // self.pk_from_sk_y
            //     .constraints_linearized(zeta),
            // self.vrf_out_x.constraints_linearized(zeta),
            // self.vrf_out_y.constraints_linearized(zeta),
            // self.sole_signer_inner_prod_acc.constraints_linearized(zeta),
        ]
        .concat()
    }

    fn domain(&self) -> &Domain<F> {
        &self.domain
    }

    fn result(&self) -> Self::Instance {
        self.pk_from_sk.result
    }
}
