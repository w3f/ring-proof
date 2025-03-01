use ark_ec::twisted_edwards::{Affine, TECurveConfig};
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::Evaluations;
use ark_std::{vec, vec::Vec};
use std::rc::Rc;
use w3f_pcs::pcs::Commitment;

use crate::piop::params::PiopParams;
use crate::piop::FixedColumns;
use crate::piop::{RingCommitments, RingEvaluations};
use w3f_plonk_common::domain::Domain;
use w3f_plonk_common::gadgets::booleanity::{BitColumn, Booleanity};
use w3f_plonk_common::gadgets::ec::AffineColumn;
use w3f_plonk_common::gadgets::ec::CondAdd;
use w3f_plonk_common::gadgets::fixed_cells::FixedCells;
use w3f_plonk_common::gadgets::ProverGadget;
use w3f_plonk_common::piop::ProverPiop;

use w3f_plonk_common::gadgets::ec::te_doubling::Doubling;
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
    doublings_of_g: Rc<AffineColumn<F, Affine<Curve>>>,
    sk_bits: Rc<BitColumn<F>>,
    // gadgets
    sk_bits_bool: Booleanity<F>,
    pk_from_sk: CondAdd<F, Affine<Curve>>,
    doublings_of_in_gadget: Doubling<F, Affine<Curve>>,
    out_from_in: CondAdd<F, Affine<Curve>>,
    out_from_in_x: FixedCells<F>,
    out_from_in_y: FixedCells<F>,
}

impl<F: PrimeField, Curve: TECurveConfig<BaseField = F>> PiopProver<F, Curve> {
    pub fn build(
        params: &PiopParams<F, Curve>,
        fixed_columns: FixedColumns<F, Affine<Curve>>, // TODO: rename to AdviceColumns
        _signer_index: usize, //TODO
        sk: Curve::ScalarField,
        vrf_in: Affine<Curve>,
    ) -> Self {
        let domain = params.domain.clone();

        let doublings_of_g = fixed_columns.doublings_of_g;
        let doublings_of_g = Rc::new(doublings_of_g);

        let sk_bits = {
            let mut sk_bits = params.scalar_part(sk); //TODO: return right thing
            assert!(sk_bits.len() <= domain.capacity - 1);
            sk_bits.resize(domain.capacity - 1, false);
            let sk_bits = BitColumn::init(sk_bits, &params.domain);
            Rc::new(sk_bits)
        };
        let sk_bits_bool = Booleanity::init(sk_bits.clone());

        // `PK_sk := sk.G`
        let pk_from_sk = CondAdd::init(
            sk_bits.clone(),
            doublings_of_g.clone(),
            params.seed,
            &domain,
        );

        let doublings_of_in_gadget = Doubling::init(vrf_in, &domain);
        let doublings_of_in = doublings_of_in_gadget.doublings.clone();
        let out_from_in = CondAdd::init(
            sk_bits.clone(),
            doublings_of_in,
            params.seed,
            &domain,
        );
        let out_from_in_x = FixedCells::init(out_from_in.acc.xs.clone(), &domain);
        let out_from_in_y = FixedCells::init(out_from_in.acc.ys.clone(), &domain);

        Self {
            domain,
            doublings_of_g,
            sk_bits,
            sk_bits_bool,
            pk_from_sk,
            doublings_of_in_gadget,
            out_from_in,
            out_from_in_x,
            out_from_in_y,
        }
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
        let sk_bits = commit(self.sk_bits.as_poly());
        let pk_from_sk = [
            commit(self.pk_from_sk.acc.xs.as_poly()),
            commit(self.pk_from_sk.acc.ys.as_poly()),
        ];
        let doublings_of_in = [
            commit(self.doublings_of_in_gadget.doublings.xs.as_poly()),
            commit(self.doublings_of_in_gadget.doublings.ys.as_poly()),
        ];
        let out_from_in = [
            commit(self.out_from_in.acc.xs.as_poly()),
            commit(self.out_from_in.acc.ys.as_poly()),
        ];
        RingCommitments {
            sk_bits,
            pk_from_sk,
            doublings_of_in,
            out_from_in,
            phantom: Default::default(),
        }
    }

    // Should return polynomials in order the consistent with
    // Self::Evaluations::to_vec() and Self::Commitments::to_vec().
    fn columns(&self) -> Vec<DensePolynomial<F>> {
        vec![
            self.sk_bits.as_poly().clone(),
            self.pk_from_sk.acc.xs.as_poly().clone(),
            self.pk_from_sk.acc.ys.as_poly().clone(),
            self.doublings_of_in_gadget.doublings.xs.as_poly().clone(),
            self.doublings_of_in_gadget.doublings.ys.as_poly().clone(),
            self.out_from_in.acc.xs.as_poly().clone(),
            self.out_from_in.acc.ys.as_poly().clone(),
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
        let doublings_of_in = [
            self.doublings_of_in_gadget.doublings.xs.evaluate(zeta),
            self.doublings_of_in_gadget.doublings.ys.evaluate(zeta),
        ];
        let out_from_in = [
            self.out_from_in.acc.xs.evaluate(zeta),
            self.out_from_in.acc.ys.evaluate(zeta),
        ];
        RingEvaluations {
            doublings_of_g,
            sk_bits,
            pk_from_sk,
            doublings_of_in,
            out_from_in,
        }
    }

    fn constraints(&self) -> Vec<Evaluations<F>> {
        [
            self.sk_bits_bool.constraints(),
            self.pk_from_sk.constraints(),
            self.doublings_of_in_gadget.constraints(),
            // self.out_from_in.constraints(),
            // self.out_from_in_x.constraints(),
            // self.out_from_in_y.constraints(),
        ]
        .concat()
    }

    fn constraints_lin(&self, zeta: &F) -> Vec<DensePolynomial<F>> {
        [
            self.sk_bits_bool.constraints_linearized(zeta),
            self.pk_from_sk.constraints_linearized(zeta),
            self.doublings_of_in_gadget.constraints_linearized(zeta),
            // self.out_from_in.constraints_linearized(zeta),
            // self.out_from_in_x.constraints_linearized(zeta),
            // self.out_from_in_y.constraints_linearized(zeta),
        ]
        .concat()
    }

    fn domain(&self) -> &Domain<F> {
        &self.domain
    }

    fn result(&self) -> Self::Instance {
        self.out_from_in.result
    }
}
