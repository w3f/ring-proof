use ark_ec::{AffineRepr, CurveConfig};
use ark_ff::{FftField, Field};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use crate::{Column, FieldColumn};
use crate::domain::Domain;
use crate::gadgets::booleanity::BitColumn;

// A vec of affine points from the prime-order subgroup of the curve whose base field enables FFTs,
// and its convenience representation as columns of coordinates over the curve's base field.
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct AffineColumn<F: FftField, P: AffineRepr<BaseField=F>> {
    pub (super) points: Vec<P>,
    pub xs: FieldColumn<F>,
    pub ys: FieldColumn<F>,
}

impl<F: FftField, P: AffineRepr<BaseField=F>> AffineColumn<F, P> {

    fn column(points: Vec<P>, domain: &Domain<F>, hidden: bool) -> Self {
        assert!(points.iter().all(|p| !p.is_zero()));
        let (xs, ys) = points.iter()
            .map(|p| p.xy().unwrap())
            .unzip();
        let xs = domain.column(xs, hidden);
        let ys = domain.column(ys, hidden);
        Self { points, xs, ys }
    }
    pub fn private_column(points: Vec<P>, domain: &Domain<F>) -> Self {
        Self::column(points, domain, true)
    }

    pub fn public_column(points: Vec<P>, domain: &Domain<F>) -> Self {
        Self::column(points, domain, false)
    }

    pub fn evaluate(&self, z: &F) -> (F, F) {
        (self.xs.evaluate(z), self.ys.evaluate(z))
    }
}

pub trait CondAdd<F, Curve, AffinePoint, ContAddVal> where
    F: FftField,
    Curve: CurveConfig,
    AffinePoint: AffineRepr<BaseField=F>,
    ContAddVal: CondAddValues<F>
{
    fn init(bitmask: BitColumn<F>,
                points: AffineColumn<F, AffinePoint>,
                seed: AffinePoint,
                domain: &Domain<F>) -> Self;
    
    fn evaluate_assignment(&self, z: &F) -> ContAddVal;

}

pub trait CondAddValues<F>
    where F: Field
{
    fn acc_coeffs_1(&self) -> (F, F);
    fn acc_coeffs_2(&self) -> (F, F);

}
