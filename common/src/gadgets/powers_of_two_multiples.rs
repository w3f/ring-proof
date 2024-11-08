use core::marker::PhantomData;

use ark_ec::{AffineRepr, CurveGroup, Group};
use ark_ff::{FftField, Field, PrimeField};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{Evaluations, GeneralEvaluationDomain};
use ark_std::{vec, vec::Vec};

use crate::domain::Domain;
use crate::gadgets::{ProverGadget, VerifierGadget};
use crate::{const_evals, AffineColumn, Column, FieldColumn};

pub trait PowersOfTwoMultiples<F, AffinePoint>
where
    F: FftField,
    AffinePoint: AffineRepr<BaseField = F>,
{
    type PowersOfTwoMultipleValT: PowersOfTwoMultipleValues<F>;
    fn init(point: AffinePoint, domain: &Domain<F>) -> Self;

    fn evaluate_assignment(&self, z: &F) -> Self::PowersOfTwoMultipleValT;
}

pub trait PowersOfTwoMultipleValues<F>
where
    F: Field,
{
    fn acc_coeffs_1(&self) -> F;
    fn acc_coeffs_2(&self) -> F;

    fn init(points: (F, F), not_last: F, acc: (F, F)) -> Self;
}

use ark_ec::twisted_edwards::{Affine as TEAffine, TECurveConfig};

// Conditional affine addition:
// if the bit is set for a point, add the point to the acc and store,
// otherwise copy the acc value
pub struct PowersOfTwoMultiplesTE<F: FftField, Curve: TECurveConfig<BaseField = F>> {
    // The  Base point we are computing its multiples
    pub(crate) point: TEAffine<Curve>,
    // The polynomial `X - w^{n-1}` in the Lagrange basis
    pub(crate) not_last: FieldColumn<F>,
    // 2^i * point
    pub(crate) multiples: AffineColumn<F, TEAffine<Curve>>,
}

impl<F, Curve> PowersOfTwoMultiples<F, TEAffine<Curve>> for PowersOfTwoMultiplesTE<F, Curve>
where
    F: FftField,
    Curve: TECurveConfig<BaseField = F>,
{
    type PowersOfTwoMultipleValT = PowersOfTwoMultipleValuesTE<F, Curve>;
    fn init(point: TEAffine<Curve>, domain: &Domain<F>) -> Self {
        let scalar_bitlen = <TEAffine<Curve> as AffineRepr>::ScalarField::MODULUS_BIT_SIZE as usize;
        let not_last = domain.not_last_row.clone();
        let mut multiples = Vec::with_capacity(scalar_bitlen);
        let mut point_multiple = point.into_group();
        multiples.push(point_multiple);
        for _ in 1..domain.capacity {
            point_multiple.double_in_place();
            multiples.push(point_multiple);
        }

        let multiples =
            AffineColumn::public_column(CurveGroup::normalize_batch(&multiples), domain);
        Self {
            point,
            not_last,
            multiples,
        }
    }

    fn evaluate_assignment(&self, z: &F) -> PowersOfTwoMultipleValuesTE<F, Curve> {
        PowersOfTwoMultipleValuesTE {
            point: (*self.point.x().unwrap(), *self.point.y().unwrap()),
            not_last: self.not_last.evaluate(z),
            multiples: self.multiples.evaluate(z),
            _curve: PhantomData,
        }
    }
}

impl<F: FftField, Curve: TECurveConfig<BaseField = F>> ProverGadget<F>
    for PowersOfTwoMultiplesTE<F, Curve>
{
    fn witness_columns(&self) -> Vec<DensePolynomial<F>> {
        vec![
            self.multiples.xs.poly.clone(),
            self.multiples.ys.poly.clone(),
        ]
    }

    fn constraints(&self) -> Vec<Evaluations<F>> {
        let domain = self.not_last.domain_4x();
        let two = &const_evals(F::one() + F::one(), domain);
        let te_a_coeff = &const_evals(Curve::COEFF_A, domain);
        let (x1, y1) = (&self.multiples.xs.evals_4x, &self.multiples.ys.evals_4x);
        let (x2, y2) = (
            &self.multiples.xs.shifted_4x(),
            &self.multiples.ys.shifted_4x(),
        );

        //x_2(y_1^2+a x_1^2)-2 x_1 y_1=0
        #[rustfmt::skip]
        let mut c1 =
            &(
                x2 *
                    &(
                        &(y1 * y1)
                             +
                             &(te_a_coeff *
                               &(x1 * x1)
                             )
                     )
                )
                    -
                    &(
                        two *
                            &(x1 * y1)
                    );

        //y_2 (2-y_1^2-a x_1^2)-y_1^2+a x_1^2=0
        #[rustfmt::skip]
        let mut c2 =
                &(
                    y2 *
                        &(
                            two -
                                &(
                                    &(y1 * y1)
                                        +
                                        &( te_a_coeff *
                                           &( x1 * x1)
                                        )
                                )
                        )
                )
                    -
                    &(
                        &(y1 * y1)
                            -
                            &(
                                te_a_coeff *
                                    &( x1 * x1 )
                            )
                    );

        let not_last = &self.not_last.evals_4x;
        c1 *= not_last;
        c2 *= not_last;

        vec![c1, c2]
    }

    fn constraints_linearized(&self, z: &F) -> Vec<DensePolynomial<F>> {
        let vals = self.evaluate_assignment(z);
        let multiples_x = self.multiples.xs.as_poly();
        let multiples_y = self.multiples.ys.as_poly();

        let c1_lin = multiples_x * vals.acc_coeffs_1(); // multiples_y coeff is 0

        let c2_lin = multiples_y * vals.acc_coeffs_2(); //multiples_x coeff is 0

        vec![c1_lin, c2_lin]
    }

    fn domain(&self) -> GeneralEvaluationDomain<F> {
        self.not_last.domain()
    }
}

pub struct PowersOfTwoMultipleValuesTE<F: Field, Curve: TECurveConfig<BaseField = F>> {
    pub point: (F, F),
    pub not_last: F,
    pub multiples: (F, F),
    pub _curve: PhantomData<Curve>,
}

impl<F: Field, Curve: TECurveConfig<BaseField = F>> PowersOfTwoMultipleValues<F>
    for PowersOfTwoMultipleValuesTE<F, Curve>
{
    fn init(point: (F, F), not_last: F, multiples: (F, F)) -> Self {
        PowersOfTwoMultipleValuesTE::<F, Curve> {
            point,
            not_last,
            multiples,
            _curve: PhantomData,
        }
    }

    fn acc_coeffs_1(&self) -> F {
        let (x1, y1) = self.point;

        //y_2 (2-y_1^2-a x_1^2)-y_1^2+a x_1^2=0
        (y1 * y1 - Curve::COEFF_A * x1 * x1) * self.not_last
    }

    fn acc_coeffs_2(&self) -> F {
        let (x1, y1) = self.multiples;

        //y_2 (2-y_1^2-a x_1^2)-y_1^2+a x_1^2=0
        (F::one() + F::one() - y1 * y1 - Curve::COEFF_A * x1 * x1) * self.not_last
    }
}

impl<F: Field, Curve: TECurveConfig<BaseField = F>> VerifierGadget<F>
    for PowersOfTwoMultipleValuesTE<F, Curve>
{
    fn evaluate_constraints_main(&self) -> Vec<F> {
        let (x1, y1) = self.multiples;
        let (x2, y2) = (F::zero(), F::zero());
        let two = F::one() + F::one();

        //x_2(y_1^2+a x_1^2)-2 x_1 y_1=0
        let mut c1 = x2 * (y1 * y1 + Curve::COEFF_A * x1 * x1) - two * x1 * y1;

        //y_2 (2-y_1^2-a x_1^2)-y_1^2+a x_1^2=0
        let mut c2 =
            y2 * (two - y1 * y1 - Curve::COEFF_A * x1 * x1) - y1 * y1 + Curve::COEFF_A * x1 * x1;

        c1 *= self.not_last;
        c2 *= self.not_last;

        vec![c1, c2]
    }
}

#[cfg(test)]
mod tests {
    use ark_ed_on_bls12_381_bandersnatch::{EdwardsAffine, EdwardsConfig, Fq};
    use ark_poly::Polynomial;
    use ark_std::test_rng;
    use ark_std::UniformRand;

    use crate::test_helpers::power_of_two_multiple;

    use super::*;

    fn _test_powers_of_two_multiples_te_gadget(
        hiding: bool,
    ) -> (
        Domain<Fq>,
        PowersOfTwoMultiplesTE<Fq, EdwardsConfig>,
        Vec<Evaluations<Fq>>,
    ) {
        let rng = &mut test_rng();

        let log_n = 8;
        let n = 2usize.pow(log_n);
        let domain = Domain::new(n, hiding);

        let point = EdwardsAffine::rand(rng);

        let scalar_bitlen = domain.capacity; //<EdwardsAffine as AffineRepr>::ScalarField::MODULUS_BIT_SIZE as usize;
        let expected_res = power_of_two_multiple(point, scalar_bitlen);

        let gadget = PowersOfTwoMultiplesTE::init(point, &domain);
        let res = gadget.multiples.points.last().unwrap();
        assert_eq!(res, &expected_res);

        let cs = gadget.constraints();
        let (c1, c2) = (&cs[0], &cs[1]);
        let c1 = c1.interpolate_by_ref();
        let c2 = c2.interpolate_by_ref();

        assert_eq!(c1.degree(), 3 * n - 2);
        assert_eq!(c2.degree(), 3 * n - 2);

        domain.divide_by_vanishing_poly(&c1);
        domain.divide_by_vanishing_poly(&c2);

        (domain, gadget, cs)
    }

    #[test]
    fn test_te_cond_add_gadget() {
        _test_powers_of_two_multiples_te_gadget(false);
        _test_powers_of_two_multiples_te_gadget(true);
    }
}
