use ark_ec::{AffineRepr, CurveGroup, Group};
use ark_ff::{BigInteger, PrimeField};
use ark_std::{vec, vec::Vec};

use common::domain::Domain;
use common::AffineColumn;

use crate::piop::FixedColumns;

#[derive(Clone)]
pub struct PiopParams<F: PrimeField, P: AffineRepr<BaseField = F>> {
    // Domain over which the piop is represented.
    pub(crate) domain: Domain<F>,

    // Number of bits used to represent a jubjub scalar.
    pub(crate) scalar_bitlen: usize,

    // Length of the part of the column representing the public keys (including the padding).
    pub keyset_part_size: usize,

    // The blinding base, a point from jubjub.
    pub(crate) h: P,

    // The point to start the summation from (as zero doesn't have a SW affine representation),
    // should be from the jubjub prime-order subgroup complement.
    pub(crate) seed: P,

    // The point used to pad the actual list of public keys. Should be of an unknown dlog.
    pub(crate) padding_point: P,
}

impl<F: PrimeField, P: AffineRepr<BaseField = F>> PiopParams<F, P> {
    pub fn setup(domain: Domain<F>, h: P, seed: P, padding_point: P) -> Self {
        let scalar_bitlen = P::ScalarField::MODULUS_BIT_SIZE as usize;
        // 1 accounts for the last cells of the points and bits columns that remain unconstrained
        let keyset_part_size = domain.capacity - scalar_bitlen - 1;
        Self {
            domain,
            scalar_bitlen,
            keyset_part_size,
            h,
            seed,
            padding_point,
        }
    }

    pub fn fixed_columns(&self, keys: &[P]) -> FixedColumns<F, P> {
        let ring_selector = self.keyset_part_selector();
        let ring_selector = self.domain.public_column(ring_selector);
        let points = self.points_column(keys);
        FixedColumns {
            points,
            ring_selector,
        }
    }

    pub fn points_column(&self, keys: &[P]) -> AffineColumn<F, P> {
        assert!(keys.len() <= self.keyset_part_size);
        let padding_len = self.keyset_part_size - keys.len();
        let padding = vec![self.padding_point; padding_len];
        let points = [keys, &padding, &self.power_of_2_multiples_of_h()].concat();
        assert_eq!(points.len(), self.domain.capacity - 1);
        AffineColumn::public_column(points, &self.domain)
    }

    pub fn power_of_2_multiples_of_h(&self) -> Vec<P> {
        let mut h = self.h.into_group();
        let mut multiples = Vec::with_capacity(self.scalar_bitlen);
        multiples.push(h);
        for _ in 1..self.scalar_bitlen {
            h.double_in_place();
            multiples.push(h);
        }
        CurveGroup::normalize_batch(&multiples)
    }

    pub fn scalar_part(&self, e: P::ScalarField) -> Vec<bool> {
        let bits_with_trailing_zeroes = e.into_bigint().to_bits_le();
        let significant_bits = &bits_with_trailing_zeroes[..self.scalar_bitlen];
        significant_bits.to_vec()
    }

    pub fn keyset_part_selector(&self) -> Vec<F> {
        [
            vec![F::one(); self.keyset_part_size],
            vec![F::zero(); self.scalar_bitlen],
        ]
        .concat()
    }

    pub fn padding_point(&self) -> P {
        self.padding_point
    }
}

#[cfg(test)]
mod tests {
    use ark_ec::AffineRepr;
    use ark_ed_on_bls12_381_bandersnatch::{EdwardsAffine, Fq, Fr, SWAffine};
    use ark_std::{test_rng, UniformRand};

    use common::domain::Domain;
    use common::test_helpers::cond_sum;

    use crate::piop::params::PiopParams;

    fn _test_powers_of_h<P: AffineRepr<BaseField = Fq, ScalarField = Fr>>() {
        let rng = &mut test_rng();
        let h = P::rand(rng);
        let seed = P::rand(rng);
        let pad = P::rand(rng);
        let domain = Domain::new(1024, false);
        let params = PiopParams::<Fq, P>::setup(domain, h, seed, pad);
        let t = Fr::rand(rng);
        let t_bits = params.scalar_part(t);
        let th = cond_sum(&t_bits, &params.power_of_2_multiples_of_h());
        assert_eq!(th, params.h.mul(t).into());
    }

    #[test]
    fn test_powers_of_h_te() {
        _test_powers_of_h::<EdwardsAffine>();
    }

    #[test]
    fn test_powers_of_h_sw() {
        _test_powers_of_h::<SWAffine>();
    }
}
