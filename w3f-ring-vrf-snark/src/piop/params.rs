use ark_ec::{AdditiveGroup, AffineRepr, CurveGroup,};
use ark_ff::{BigInteger, PrimeField};
use ark_std::{vec, vec::Vec};

use common::domain::Domain;
use common::AffineColumn;

use crate::piop::FixedColumns;

///
/// | 1              | 2           | 3           | 4              | 5              | 6                        | 7&8      | 9&10            | 11&12    | 13&14           |
/// | --------       | --------    | --------    | --             | -              | -                        | -        | -               | -        | -               |
/// | $k$            | $pk_x$      | $pk_y$      | $acc_{pk_x}$   | $acc_{pk_y}$   | $sk$                     | $2^iG$x2 | $acc_{sk}$x2    | $2^iH$x2 | $acc_{out}$x2   |
/// | --------       | --------    | --------    | --             | -              | -                        | -        | -               | -        | -               |
/// | signer's index | x of pubkey | y of pubkey | $\sum k_ipk_x$ | $\sum k_ipk_y$ | binary rep of secret key |          | $\sum sk_i2^iG$ |          | $\sum sk_i2^iH$ |
///
/// We do not have ring selector part so I assume that I we are going to generate two different polynomials for essentially two different table.
/// TODO: But for now we use a fake selector so we can use the inner product gadget.
/// Later with we will add a bit vector sum gadget and get rid of selector
///
#[derive(Clone)]
pub struct PiopParams<F: PrimeField, P: AffineRepr<BaseField = F>> {
    // Domain over which the piop is represented.
    pub(crate) domain: Domain<F>,

    // Number of bits used to represent a jubjub scalar.
    pub(crate) scalar_bitlen: usize,

    // size of the padded ring
    pub padded_keyset_size: usize,

    // The point to start the summation from (as zero doesn't have a SW affine representation),
    // should be from the jubjub prime-order subgroup complement.
    pub(crate) seed: P,

    // The point used to pad the actual list of public keys. Should be of an unknown dlog.
    pub(crate) padding_point: P,

}

impl<F: PrimeField, P: AffineRepr<BaseField = F>> PiopParams<F, P> {
    pub fn setup(domain: Domain<F>, seed: P, padding_point: P) -> Self {        
        let scalar_bitlen = P::ScalarField::MODULUS_BIT_SIZE as usize;
        //make sure the domain size is big enough to at least deal with VRF SNARK
        //We will later check if the domain is big enough to accomodate the ring
        assert!(domain.capacity > scalar_bitlen);
        // 1 accounts for the last cells of the points and bits columns that remain unconstrained
        let padded_keyset_size = domain.capacity - 1;
        Self {
            domain,
            scalar_bitlen,
            padded_keyset_size,
            seed,
            padding_point,
        }
    }

    pub fn fixed_columns(&self, keys: &[P]) -> FixedColumns<F, P> {
        let ring_selector = self.domain.public_column(self.keyset_part_selector());
        let pubkey_points = self.pubkey_points_column(keys);
        let power_of_2_multiples_of_gen = self.gen_multiples_column();

        FixedColumns {
            pubkey_points,
	        power_of_2_multiples_of_gen,
            ring_selector,
        }
    }

    pub fn pubkey_points_column(&self, keys: &[P]) -> AffineColumn<F, P> {
        assert!(keys.len() <= self.padded_keyset_size);
        let padding_len = self.padded_keyset_size - keys.len();
        let padding = vec![self.padding_point; padding_len];
        let points = [keys, &padding].concat();
        assert!(points.len() < self.domain.capacity); //check if it fits the domain.
        AffineColumn::public_column(points, &self.domain)
    }

    pub fn gen_multiples_column(&self) -> AffineColumn<F, P> {
	    let prime_subgroup_gen = P::generator(); //TODO: should this be fed as an input param of the ring?
	    let power_of_2_multiples_of_gen = Self::power_of_2_multiples_of(prime_subgroup_gen, self.scalar_bitlen);
        //TODO: we might need different domain for different columns
        AffineColumn::public_column(power_of_2_multiples_of_gen, &self.domain)
    }
    
    pub fn power_of_2_multiples_of(base_point: P, scalar_bitlen: usize) -> Vec<P> {
        let mut h = base_point.into_group();
        let mut multiples = Vec::with_capacity(scalar_bitlen);
        multiples.push(h);
        for _ in 1..scalar_bitlen {
            h.double_in_place();
            multiples.push(h);
        }
        CurveGroup::normalize_batch(&multiples)
    }

    pub fn power_of_2_multiples_of_h(&self) {
        
    }
    pub fn keyset_part_selector(&self) -> Vec<F> {
        [
            vec![F::one(); self.padded_keyset_size],
        ]
        .concat()
    }
    
    pub fn scalar_to_bitvec(&self, e: P::ScalarField) -> Vec<bool> {
        let bits_with_trailing_zeroes = e.into_bigint().to_bits_le();
        let significant_bits = &bits_with_trailing_zeroes[..self.scalar_bitlen];
        significant_bits.to_vec()
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
