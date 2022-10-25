use ark_ec::{AffineRepr, CurveGroup, Group, VariableBaseMSM};
use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::{BigInteger, PrimeField};
use ark_std::rand::Rng;
use ark_std::UniformRand;
use fflonk::pcs::kzg::KZG;
use fflonk::pcs::PCS;
use common::domain::Domain;
use common::gadgets::booleanity::BitColumn;


#[derive(Clone)]
pub struct PiopParams<F: PrimeField, Curve: SWCurveConfig<BaseField=F>> {
    // Domain over which the piop is represented.
    pub(crate) domain: Domain<F>,

    // Length of the part of the column representing the public keys (including the padding).
    pub(crate) keyset_part_size: usize,

    // The blinding base, a point from jubjub.
    pub h: Affine<Curve>,

}

impl<F: PrimeField, Curve: SWCurveConfig<BaseField=F>> PiopParams<F, Curve> {
    // Number of bits used to represent a jubjub scalar.
    const SC_LEN: usize = Curve::ScalarField::MODULUS_BIT_SIZE as usize;

    pub fn setup<R: Rng>(domain: Domain<F>, rng: &mut R) -> Self {
        // 1 accounts for the last cells of the points and bits columns that remain unconstrained
        let max_points = domain.capacity - 1;
        assert!(max_points > Self::SC_LEN);
        let max_keys = max_points - Self::SC_LEN;

        let h = Affine::<Curve>::rand(rng);


        Self {
            domain,
            keyset_part_size: max_keys,
            h,
        }
    }

    pub fn bits_column(&self, index_in_keys: usize, secret: Curve::ScalarField) -> BitColumn<F> {
        let keyset_part = self.ring_part(index_in_keys);
        let scalar_part = Self::scalar_part(secret);
        let bits = [
            keyset_part,
            scalar_part
        ].concat();
        assert_eq!(bits.len(), self.domain.capacity - 1);
        BitColumn::init(bits, &self.domain)
    }

    pub fn ring_part(&self, index_in_keys: usize) -> Vec<bool> {
        let mut keyset_part = vec![false; self.keyset_part_size];
        keyset_part[index_in_keys] = true;
        keyset_part
    }

    pub fn power_of_2_multiples_of_h(&self) -> Vec<Affine::<Curve>> {
        let mut h = self.h.into_group();
        let mut multiples = Vec::with_capacity(Self::SC_LEN);
        multiples.push(h);
        for _ in 1..Self::SC_LEN {
            h.double_in_place();
            multiples.push(h);
        }
        CurveGroup::normalize_batch(&multiples)
    }

    pub fn scalar_part(e: Curve::ScalarField) -> Vec<bool> {
        let bits_with_trailing_zeroes = e.into_bigint().to_bits_le();
        let significant_bits = &bits_with_trailing_zeroes[..Self::SC_LEN];
        significant_bits.to_vec()
    }

    pub fn keyset_part_selector(&self) -> Vec<F> {
        [
            vec![F::one(); self.keyset_part_size],
            vec![F::zero(); Self::SC_LEN]
        ].concat()
    }
}

#[cfg(test)]
mod tests {
    use ark_ed_on_bls12_381_bandersnatch::{BandersnatchParameters, Fq, Fr};
    use ark_std::{test_rng, UniformRand};
    use std::ops::Mul;
    use common::domain::Domain;
    use common::test_helpers::cond_sum;
    use crate::piop::params::PiopParams;

    #[test]
    fn test_powers_of_h() {
        let rng = &mut test_rng();
        let domain = Domain::new(1024, false);
        let params = PiopParams::<Fq, BandersnatchParameters>::setup(domain, rng);
        let t = Fr::rand(rng);
        let t_bits = PiopParams::<Fq, BandersnatchParameters>::scalar_part(t);
        let th = cond_sum(&t_bits, &params.power_of_2_multiples_of_h());
        assert_eq!(th, params.h.mul(t));
    }
}