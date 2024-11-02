use ark_ec::AffineRepr;
use ark_ff::{FftField, Field};

use crate::domain::Domain;
use crate::gadgets::booleanity::BitColumn;
use crate::{AffineColumn, FieldColumn};

use super::{ProverGadget, VerifierGadget};

/// Affine point with conditional add implementation.
///
/// Currently supported for Arkworks Short Weierstrass and Twisted Edwards affine points.
pub trait AffineCondAdd: AffineRepr
where
    FieldFor<Self>: FftField,
{
    /// Conditional addition operation
    type CondAddT: CondAdd<FieldFor<Self>, Self>;
}

// Conditional affine addition.
//
// If the bit is set for a point, add the point to the acc and store,
// otherwise copy the acc value
pub trait CondAdd<F, P>: ProverGadget<F>
where
    F: FftField,
    P: AffineRepr<BaseField = F>,
{
    type Values: CondAddValues<F>;

    fn init(bitmask: BitColumn<F>, points: AffineColumn<F, P>, seed: P, domain: &Domain<F>)
        -> Self;

    fn evaluate_assignment(&self, z: &F) -> Self::Values;

    fn get_acc(&self) -> AffineColumn<F, P>;

    fn get_result(&self) -> P;
}

pub trait CondAddValues<F>: VerifierGadget<F>
where
    F: Field,
{
    fn acc_coeffs_1(&self) -> (F, F);

    fn acc_coeffs_2(&self) -> (F, F);

    fn init(bitmask: F, points: (F, F), not_last: F, acc: (F, F)) -> Self;
}

type FieldFor<P> = <P as AffineRepr>::BaseField;

pub struct CondAddGen<P>
where
    P: AffineRepr,
    <P as AffineRepr>::BaseField: FftField,
{
    pub(super) bitmask: BitColumn<FieldFor<P>>,
    pub(super) points: AffineColumn<FieldFor<P>, P>,
    pub(super) not_last: FieldColumn<FieldFor<P>>,
    pub acc: AffineColumn<FieldFor<P>, P>,
    pub result: P,
}

pub struct CondAddValuesGen<P: AffineRepr> {
    pub bitmask: FieldFor<P>,
    pub points: (FieldFor<P>, FieldFor<P>),
    pub not_last: FieldFor<P>,
    pub acc: (FieldFor<P>, FieldFor<P>),
}
