use ark_ec::AffineRepr;
use ark_ff::{FftField, Field};

use crate::domain::Domain;
use crate::gadgets::booleanity::BitColumn;
use crate::AffineColumn;

pub trait CondAdd<F, AffinePoint>
where
    F: FftField,
    AffinePoint: AffineRepr<BaseField = F>,
{
    type CondAddValT: CondAddValues<F>;
    fn init(
        bitmask: BitColumn<F>,
        points: AffineColumn<F, AffinePoint>,
        seed: AffinePoint,
        domain: &Domain<F>,
    ) -> Self;

    fn evaluate_assignment(&self, z: &F) -> Self::CondAddValT;
    fn get_acc(&self) -> AffineColumn<F, AffinePoint>;
    fn get_result(&self) -> AffinePoint;
}

pub trait CondAddValues<F>
where
    F: Field,
{
    fn acc_coeffs_1(&self) -> (F, F);
    fn acc_coeffs_2(&self) -> (F, F);

    fn init(bitmask: F, points: (F, F), not_last: F, acc: (F, F)) -> Self;
}
