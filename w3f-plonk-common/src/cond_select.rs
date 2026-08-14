//! Constant-time handling of secret witness bits (the prover's ring
//! position and the blinding scalar bits). The secret bit is wrapped
//! into a `subtle::Choice` (an optimization barrier) on entry and
//! selection is limb-wise on the raw Montgomery representation, so no
//! field arithmetic -- and none of arkworks' data-dependent conditional
//! reductions -- ever touches the bits. Defense in depth, not a complete
//! countermeasure: the unconditional point additions run variable-time
//! arkworks field arithmetic on secret-derived accumulator values, and
//! the column FFTs and polynomial commitment MSMs downstream still
//! process the witness in variable time.

use ark_ec::short_weierstrass::{Projective as SwProjective, SWCurveConfig};
use ark_ec::twisted_edwards::{Projective as TeProjective, TECurveConfig};
use ark_ff::{BigInt, Field, Fp, FpConfig};
use ark_std::marker::PhantomData;
use subtle::ConditionallySelectable;

pub use subtle::Choice;

/// Conversion into [`Choice`]. The `bool` impl applies subtle's
/// optimization barrier, preventing the compiler from branching on the
/// bit later. (`From<bool> for Choice` does not exist upstream, hence
/// this local trait.)
pub trait IntoChoice {
    fn into_choice(self) -> Choice;
}

impl IntoChoice for Choice {
    fn into_choice(self) -> Choice {
        self
    }
}

impl IntoChoice for bool {
    fn into_choice(self) -> Choice {
        Choice::from(self as u8)
    }
}

/// Lifts a bit to a field element without branching on its value.
pub fn bit_to_field<F: Field + CondSelect>(bit: bool) -> F {
    F::select(bit, &F::one(), &F::zero())
}

/// Constant-time two-way select.
pub trait CondSelect: Sized {
    fn select(bit: impl IntoChoice, if_true: &Self, if_false: &Self) -> Self;
}

impl<P: FpConfig<N>, const N: usize> CondSelect for Fp<P, N> {
    fn select(bit: impl IntoChoice, if_true: &Self, if_false: &Self) -> Self {
        let limbs =
            <[u64; N]>::conditional_select(&if_false.0 .0, &if_true.0 .0, bit.into_choice());
        Fp(BigInt(limbs), PhantomData)
    }
}

impl<C: TECurveConfig> CondSelect for TeProjective<C>
where
    C::BaseField: CondSelect,
{
    fn select(bit: impl IntoChoice, if_true: &Self, if_false: &Self) -> Self {
        let choice = bit.into_choice();
        Self::new_unchecked(
            CondSelect::select(choice, &if_true.x, &if_false.x),
            CondSelect::select(choice, &if_true.y, &if_false.y),
            CondSelect::select(choice, &if_true.t, &if_false.t),
            CondSelect::select(choice, &if_true.z, &if_false.z),
        )
    }
}

impl<C: SWCurveConfig> CondSelect for SwProjective<C>
where
    C::BaseField: CondSelect,
{
    fn select(bit: impl IntoChoice, if_true: &Self, if_false: &Self) -> Self {
        let choice = bit.into_choice();
        Self::new_unchecked(
            CondSelect::select(choice, &if_true.x, &if_false.x),
            CondSelect::select(choice, &if_true.y, &if_false.y),
            CondSelect::select(choice, &if_true.z, &if_false.z),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ed_on_bls12_381_bandersnatch::{EdwardsProjective, Fq};
    use ark_std::{test_rng, UniformRand};

    #[test]
    fn bit_lift_is_exact() {
        assert_eq!(bit_to_field::<Fq>(false), Fq::from(0));
        assert_eq!(bit_to_field::<Fq>(true), Fq::from(1));
    }

    // The selected values feed committed columns, so they are
    // consensus-critical: the select must return the operand bit for bit,
    // not merely an equivalent representation.
    #[test]
    fn select_returns_exact_operand() {
        let rng = &mut test_rng();
        let a = Fq::rand(rng);
        let b = Fq::rand(rng);
        assert_eq!(Fq::select(true, &a, &b), a);
        assert_eq!(Fq::select(false, &a, &b), b);

        let p = EdwardsProjective::rand(rng);
        let q = EdwardsProjective::rand(rng);
        let s = EdwardsProjective::select(true, &p, &q);
        assert_eq!((s.x, s.y, s.t, s.z), (p.x, p.y, p.t, p.z));
        let s = EdwardsProjective::select(false, &p, &q);
        assert_eq!((s.x, s.y, s.t, s.z), (q.x, q.y, q.t, q.z));
    }
}
