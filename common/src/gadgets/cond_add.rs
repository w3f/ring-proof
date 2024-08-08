use ark_ec::{AffineRepr};
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


// Conditional affine addition:
// if the bit is set for a point, add the point to the acc and store,
// otherwise copy the acc value
pub struct CondAdd<F: FftField, P: AffineRepr<BaseField=F>> {
    pub(super)bitmask: BitColumn<F>,
    pub(super)points: AffineColumn<F, P>,
    // The polynomial `X - w^{n-1}` in the Lagrange basis
    pub(super)not_last: FieldColumn<F>,
    // Accumulates the (conditional) rolling sum of the points
    pub acc: AffineColumn<F, P>,
    pub result: P,
}

pub struct CondAddValues<F: Field> {
    pub bitmask: F,
    pub points: (F, F),
    pub not_last: F,
    pub acc: (F, F),
}


