use std::iter;
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::{FftField, Field};
use ark_poly::{Evaluations, GeneralEvaluationDomain};
use ark_poly::univariate::DensePolynomial;
use ark_std::{test_rng, UniformRand, Zero};
use crate::{Column, const_evals, FieldColumn};
use crate::gadgets::booleanity::BitColumn;
use crate::gadgets::{ProverGadget, VerifierGadget};


// A vec of affine points from the prime-order subgroup of the curve whose base field enables FFTs,
// and its convenience representation as columns of coordinates over the curve's base field.
#[derive(Clone)]
pub struct AffineColumn<F: FftField, P: AffineCurve<BaseField=F>> {
    points: Vec<P>,
    pub xs: FieldColumn<F>,
    pub ys: FieldColumn<F>,
}

impl<F: FftField, P: AffineCurve<BaseField=F>> AffineColumn<F, P> {
    pub fn init(points: Vec<P>) -> Self {
        assert!(points.iter().all(|p| !p.is_zero()));
        let (xs, ys) = points.iter()
            .map(|p| p.xy().unwrap())
            .unzip();
        let xs = FieldColumn::init(xs);
        let ys = FieldColumn::init(ys);
        Self { points, xs, ys }
    }

    pub fn evaluate(&self, z: &F) -> (F, F) {
        (self.xs.evaluate(z), self.ys.evaluate(z))
    }
}


// Conditional affine addition:
// if the bit is set for a point, add the point to the acc and store,
// otherwise copy the acc value
pub struct CondAdd<F: FftField, P: AffineCurve<BaseField=F>> {
    bitmask: BitColumn<F>,
    points: AffineColumn<F, P>,
    // The polynomial `X - w^{n-1}` in the Lagrange basis
    not_last: FieldColumn<F>,
    pub acc: AffineColumn<F, P>, // accumulates the (conditional) rolling sum of the points
    pub result: P,
}

pub struct CondAddValues<F: Field> {
    pub bitmask: F,
    pub points: (F, F),
    pub not_last: F,
    pub acc: (F, F),
}


impl<F, Curve> CondAdd<F, Affine<Curve>> where
    F: FftField,
    Curve: SWCurveConfig<BaseField=F>,
{
    // Populates the acc column starting from the supplied initial point (as 0 doesn't have an affine representation).
    // The last point of the input column is ignored, as adding it would made the acc column overflow due the initial point.
    pub fn init(bitmask: BitColumn<F>,
                points: AffineColumn<F, Affine<Curve>>,
                not_last: FieldColumn<F>) -> Self {
        let n = bitmask.size();
        assert_eq!(points.points.len(), n);
        assert_eq!(not_last.evals.evals.len(), n);
        let init = Self::point_in_g1_complement();
        assert!(!init.is_zero());
        let acc = bitmask.bits.iter()
            .zip(points.points.iter())
            .take(n - 1) // last cells are ignored
            .scan(init.clone(), |acc, (&b, point)| {
                if b {
                    *acc += point;
                }
                Some(*acc)
            });
        let acc: Vec<_> = iter::once(init)
            .chain(acc)
            .collect();
        let init_plus_result = acc.last().unwrap();
        let result = init_plus_result.into_projective() - init.into_projective();
        let result = result.into_affine();
        let acc = AffineColumn::init(acc);

        Self { bitmask, points, acc, not_last, result }
    }

    fn evaluate_assignment(&self, z: &F) -> CondAddValues<F> {
        CondAddValues {
            bitmask: self.bitmask.evaluate(z),
            points: self.points.evaluate(z),
            not_last: self.not_last.evaluate(z),
            acc: self.acc.evaluate(z),
        }
    }

    //TODO: find
    pub fn point_in_g1_complement() -> Affine<Curve> {
        Affine::<Curve>::prime_subgroup_generator()
    }
}


impl<F, Curve> ProverGadget<F> for CondAdd<F, Affine<Curve>>
    where
        F: FftField,
        Curve: SWCurveConfig<BaseField=F>,
{
    fn witness_columns(&self) -> Vec<DensePolynomial<F>> {
        vec![self.acc.xs.poly.clone(), self.acc.ys.poly.clone()]
    }

    fn constraints(&self) -> Vec<Evaluations<F>> {
        let domain = self.bitmask.domain_4x();
        let b = &self.bitmask.col.evals_4x;
        let one = &const_evals(F::one(), domain);
        let (x1, y1) = (&self.acc.xs.evals_4x, &self.acc.ys.evals_4x);
        let (x2, y2) = (&self.points.xs.evals_4x, &self.points.ys.evals_4x);
        let (x3, y3) = (&self.acc.xs.shifted_4x(), &self.acc.ys.shifted_4x());

        let mut c1 =
            &(
                b *
                    &(
                        &(
                            &(
                                &(x1 - x2) * &(x1 - x2)
                            ) *
                                &(
                                    &(x1 + x2) + x3
                                )
                        ) -
                            &(
                                &(y2 - y1) * &(y2 - y1)
                            )
                    )
            ) +
                &(
                    &(one - b) * &(y3 - y1)
                );

        let mut c2 =
            &(
                b *
                    &(
                        &(
                            &(x1 - x2) * &(y3 + y1)
                        ) -
                            &(
                                &(y2 - y1) * &(x3 - x1)
                            )
                    )
            ) +
                &(
                    &(one - b) * &(x3 - x1)
                );

        let not_last = &self.not_last.evals_4x;
        c1 *= not_last;
        c2 *= not_last;

        vec![c1, c2]
    }

    fn constraints_linearized(&self, z: &F) -> Vec<DensePolynomial<F>> {
        let vals = self.evaluate_assignment(z);
        let acc_x = self.acc.xs.as_poly();
        let acc_y = self.acc.ys.as_poly();

        let (c_acc_x, c_acc_y) = vals.acc_coeffs_1();
        let c1_lin = acc_x * c_acc_x + acc_y * c_acc_y;

        let (c_acc_x, c_acc_y) = vals.acc_coeffs_2();
        let c2_lin = acc_x * c_acc_x + acc_y * c_acc_y;

        vec![c1_lin, c2_lin]
    }

    fn domain(&self) -> GeneralEvaluationDomain<F> {
        self.bitmask.domain()
    }
}


impl<F: Field> VerifierGadget<F> for CondAddValues<F> {
    fn evaluate_constraints_main(&self) -> Vec<F> {
        let b = self.bitmask;
        let (x1, y1) = self.acc;
        let (x2, y2) = self.points;
        let (x3, y3) = (F::zero(), F::zero());

        let mut c1 =
            b * (
                (x1 - x2) * (x1 - x2) * (x1 + x2 + x3)
                    - (y2 - y1) * (y2 - y1)
            ) + (F::one() - b) * (y3 - y1);

        let mut c2 =
            b * (
                (x1 - x2) * (y3 + y1)
                    - (y2 - y1) * (x3 - x1)
            ) + (F::one() - b) * (x3 - x1);

        c1 *= self.not_last;
        c2 *= self.not_last;

        vec![c1, c2]
    }
}


impl<F: Field> CondAddValues<F> {
    pub fn acc_coeffs_1(&self) -> (F, F) {
        let b = self.bitmask;
        let (x1, _y1) = self.acc;
        let (x2, _y2) = self.points;

        let mut c_acc_x = b * (x1 - x2) * (x1 - x2);
        let mut c_acc_y = F::one() - b;

        c_acc_x *= self.not_last;
        c_acc_y *= self.not_last;

        (c_acc_x, c_acc_y)
    }

    pub fn acc_coeffs_2(&self) -> (F, F) {
        let b = self.bitmask;
        let (x1, y1) = self.acc;
        let (x2, y2) = self.points;

        let mut c_acc_x = b * (y1 - y2) + F::one() - b;
        let mut c_acc_y = b * (x1 - x2);

        c_acc_x *= self.not_last;
        c_acc_y *= self.not_last;

        (c_acc_x, c_acc_y)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{test_rng, UniformRand};
    use ark_ed_on_bls12_381_bandersnatch::{Fq, SWAffine};
    use ark_poly::{GeneralEvaluationDomain, Polynomial};
    use crate::gadgets::tests::test_gadget;
    use crate::not_last;
    use crate::test_helpers::cond_sum;
    use crate::test_helpers::*;
    use ark_poly::EvaluationDomain;


    #[test]
    fn test_sw_cond_add_gadget() {
        let rng = &mut test_rng();

        let log_n = 10;
        let n = 2usize.pow(log_n);
        let domain = GeneralEvaluationDomain::<Fq>::new(n).unwrap();

        let bitmask = random_bitvec(n, 0.5, rng);
        let points = random_vec::<SWAffine, _>(n, rng);
        let init = CondAdd::point_in_g1_complement();
        let expected_res = init + cond_sum(&bitmask[..n - 1], &points[..n - 1]);

        let bitmask_col = BitColumn::init(bitmask);
        let points_col = AffineColumn::init(points);
        let not_last_col = FieldColumn::from_poly(not_last(domain), n);
        let gadget = CondAdd::init(bitmask_col, points_col, not_last_col);
        let res = gadget.acc.points.last().unwrap();
        assert_eq!(res, &expected_res);

        let cs = gadget.constraints();
        let (c1, c2) = (&cs[0], &cs[1]);
        let c1 = c1.interpolate_by_ref();
        let c2 = c2.interpolate_by_ref();
        assert_eq!(c1.degree(), 4 * n - 3);
        assert_eq!(c2.degree(), 3 * n - 2);

        test_gadget(gadget);
    }
}