use ark_ec::twisted_edwards::{Affine, TECurveConfig};
use ark_ec::{AdditiveGroup, AffineRepr, CurveGroup};
use ark_ff::{BigInteger, PrimeField};
use ark_std::{vec, vec::Vec};

use w3f_plonk_common::domain::Domain;
use w3f_plonk_common::gadgets::booleanity::BitColumn;
use w3f_plonk_common::gadgets::ec::AffineColumn;

use crate::piop::FixedColumns;

/// Plonk Interactive Oracle Proofs (PIOP) parameters.
#[derive(Clone)]
pub struct PiopParams<F: PrimeField, Curve: TECurveConfig<BaseField = F>> {
    /// Domain over which the piop is represented.
    pub(crate) domain: Domain<F>,
    /// Number of bits used to represent a jubjub scalar.
    pub(crate) scalar_bitlen: usize,
    /// Length of the part of the column representing the public keys (including the padding).
    pub keyset_part_size: usize,
    /// The generator used to compute public keys, `pk=sk.G`.
    pub(crate) g: Affine<Curve>,
    /// Blinding base point.
    pub(crate) h: Affine<Curve>,
    /// Summation base point.
    pub(crate) seed: Affine<Curve>,
    /// The point used to pad the list of public keys.
    pub(crate) padding: Affine<Curve>,
}

impl<F: PrimeField, Curve: TECurveConfig<BaseField = F>> PiopParams<F, Curve> {
    /// Initialize PIOP parameters.
    ///
    /// - `domain`: polynomials evaluation domain.
    /// - `h`: Blinding base point.
    /// - `seed`: Accumulation base point
    /// - `padding`: The point used to pad the list of public keys.
    ///
    /// All points should be of an unknown discrete log.
    pub fn setup(
        domain: Domain<F>,
        h: Affine<Curve>,
        seed: Affine<Curve>,
        padding: Affine<Curve>,
    ) -> Self {
        let scalar_bitlen = Curve::ScalarField::MODULUS_BIT_SIZE as usize;
        // 1 accounts for the last cells of the points and bits columns that remain unconstrained
        let keyset_part_size = domain.capacity - scalar_bitlen - 1;
        Self {
            domain,
            scalar_bitlen,
            keyset_part_size,
            g: Affine::<Curve>::generator(),
            h,
            seed,
            padding,
        }
    }

    pub fn fixed_columns(&self, pks: Vec<Affine<Curve>>) -> FixedColumns<F, Affine<Curve>> {
        let doublings_of_g = self.doublings_of_g_col();
        let mut pks = pks;
        assert!(pks.len() <= self.domain.capacity - 1);
        pks.resize(self.domain.capacity - 1, self.padding);
        let pks = AffineColumn::public_column(pks, &self.domain); // TODO: here we assume pk.len() > 256
        FixedColumns {
            pks,
            doublings_of_g,
        }
    }

    fn doublings_of_g_col(&self) -> AffineColumn<F, Affine<Curve>> {
        let mut doublings_of_g = self.doublings_of(self.g);
        doublings_of_g.resize(self.domain.capacity - 1, self.g); //TODO: eh, may be zeros
        AffineColumn::public_column(doublings_of_g, &self.domain)
    }

    pub fn pk_index_col(&self, index: usize) -> BitColumn<F> {
        assert!(index <= self.domain.capacity - 1);
        let mut col = vec![false; self.domain.capacity - 1];
        col[index] = true;
        BitColumn::init(col, &self.domain)
    }

    pub fn sk_bits(&self, sk: Curve::ScalarField) -> Vec<bool> {
        let bits_with_trailing_zeroes = sk.into_bigint().to_bits_le();
        let significant_bits = &bits_with_trailing_zeroes[..self.scalar_bitlen];
        significant_bits.to_vec()
    }

    pub fn doublings_of(&self, p: Affine<Curve>) -> Vec<Affine<Curve>> {
        let mut p = p.into_group();
        let mut doublings = Vec::with_capacity(self.domain.capacity);
        doublings.push(p);
        for _ in 1..self.scalar_bitlen {
            p.double_in_place();
            doublings.push(p);
        }
        CurveGroup::normalize_batch(&doublings)
    }
}

// #[cfg(test)]
// mod tests {
//     use ark_ed_on_bls12_381_bandersnatch::{BandersnatchConfig, EdwardsAffine, Fq, Fr};
//     use ark_std::ops::Mul;
//     use ark_std::{test_rng, UniformRand};
//
//     use w3f_plonk_common::domain::Domain;
//     use w3f_plonk_common::test_helpers::cond_sum;
//
//     use crate::piop::params::PiopParams;
//
//     #[test]
//     fn test_powers_of_h() {
//         let rng = &mut test_rng();
//         let h = EdwardsAffine::rand(rng);
//         let seed = EdwardsAffine::rand(rng);
//         let padding = EdwardsAffine::rand(rng);
//         let domain = Domain::new(1024, false);
//
//         let params = PiopParams::<Fq, BandersnatchConfig>::setup(domain, h, seed, padding);
//         let t = Fr::rand(rng);
//         let t_bits = params.scalar_part(t);
//         let th = cond_sum(&t_bits, &params.power_of_2_multiples_of_h());
//         assert_eq!(th, params.h.mul(t));
//     }
// }
