use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_std::rand::Rng;
use ark_std::UniformRand;
use std::ops::Mul;


pub struct PiopParams<F: PrimeField, Curve: SWCurveConfig<BaseField=F>> {
    // Columns' length.
    pub domain: GeneralEvaluationDomain<F>,

    // Number of bits used to represent a scalar from jubjub.
    scalar_bitlen: usize,

    // Length of the part of the column representing the public keys.
    pub(crate) keyset_part_size: usize,

    // Length of the part of the column representing the secret scalar.
    // One bit longer than `self.scalar_bitlen` as the affine_addition gadget
    // ignores the last cells of the input columns.
    scalar_part_size: usize,

    // The blinding base, a point from jubjub
    pub h: Affine<Curve>,

    // A vec of `self.scalar_part` jubjub points of the form
    // H, 2H, 4H, ..., 2^(self.scalar_part-1)H
    pub powers_of_h: Vec<Affine<Curve>>,

    pub(crate) init: Affine<Curve>,
}

impl<F: PrimeField, Curve: SWCurveConfig<BaseField=F>> PiopParams<F, Curve> {
    pub fn setup<R: Rng>(domain_size: usize, rng: &mut R) -> Self {
        let domain = GeneralEvaluationDomain::new(domain_size)
            .expect("field is not smooth enough to construct domain");
        let scalar_bitlen = Curve::ScalarField::MODULUS_BIT_SIZE as usize;
        // 1 accounts for the last cells of the point columns that remain unconstrained
        let scalar_part_size = scalar_bitlen + 1;
        let keyset_part_size = domain_size - scalar_part_size;

        let h = Affine::<Curve>::rand(rng);
        let powers_of_h = Self::power_of_2_multiples(scalar_part_size, h.into_projective());
        let powers_of_h = ProjectiveCurve::batch_normalization_into_affine(&powers_of_h);

        let init = Affine::<Curve>::rand(rng);

        Self {
            domain,
            scalar_bitlen,
            keyset_part_size,
            scalar_part_size,
            h,
            powers_of_h,
            init,
        }
    }

    fn power_of_2_multiples(scalar_part_size: usize, mut h: Projective::<Curve>) -> Vec<Projective::<Curve>> {
        let mut multiples = Vec::with_capacity(scalar_part_size);
        multiples.push(h);
        for _ in 1..scalar_part_size {
            h.double_in_place();
            multiples.push(h);
        }
        multiples
    }

    pub fn scalar_part(&self, e: Curve::ScalarField) -> Vec<bool> {
        let bits_with_trailing_zeroes = e.into_bigint().to_bits_le();
        let significant_bits = &bits_with_trailing_zeroes[..self.scalar_bitlen];
        let mut padded_bits = significant_bits.to_vec();
        padded_bits.push(false);
        padded_bits
    }
}

#[cfg(test)]
mod tests {
    use ark_ec::AffineCurve;
    use ark_ed_on_bls12_381_bandersnatch::{BandersnatchParameters, Fq, Fr};
    use ark_std::{test_rng, UniformRand};
    use std::ops::Mul;
    use common::test_helpers::cond_sum;
    use crate::piop::params::PiopParams;

    #[test]
    fn test_powers_of_h() {
        let rng = &mut test_rng();
        let params = PiopParams::<Fq, BandersnatchParameters>::setup(1024, rng);
        let t = Fr::rand(rng);
        let t_bits = params.scalar_part(t);
        let th = cond_sum(&t_bits, &params.powers_of_h);
        assert_eq!(th, params.h.mul(t));
    }
}