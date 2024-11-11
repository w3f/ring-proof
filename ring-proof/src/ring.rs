use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ff::PrimeField;
use ark_poly::EvaluationDomain;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::fmt;
use ark_std::iter;
use ark_std::ops::Range;
use ark_std::vec::Vec;
use fflonk::pcs::kzg::urs::URS;
use fflonk::pcs::PcsParams;

use common::domain::ZK_ROWS;

use crate::PiopParams;

const IDLE_ROWS: usize = ZK_ROWS + 1;

/// Commitment to a list of VRF public keys as is used as a public input to the ring proof SNARK verifier.

/// The VRF keys are (inner) curve points that we represent in the affine short Weierstrass coordinates.
/// We commit to the coordinate vectors independently using KZG on the outer curve. To make the commitment
/// updatable we use SRS in the Lagrangian form: `L1, ..., Ln`, where `Li = L_i(t)G`.
/// The commitment to a vector `a1, ..., an` is then `a1L1 + ... + anLn`.

/// We pad the list of keys with a `padding` point with unknown dlog up to a certain size.
/// Additionally, to make the commitment compatible with the snark,
/// we append the power-of-2 powers of the VRF blinding Pedersen base
/// `H, 2H, 4H, ..., 2^(s-1)H`, where `s` is the bitness of the VRF curve scalar field.
/// The last `IDLE_ROWS = 4` elements are set to `(0, 0)`.

/// Thus, the vector of points we commit to coordinatewise is
/// `pk1, ..., pkn, padding, ..., padding, H, 2H, ..., 2^(s-1)H, 0, 0, 0, 0`

// `KzgCurve` -- outer curve, subgroup of a pairing-friendly curve. We instantiate it with bls12-381 G1.
// `VrfCurveConfig` -- inner curve, the curve used by the VRF, in SW form. We instantiate it with Bandersnatch.
// `F` shared scalar field of the outer and the base field of the inner curves.
#[derive(Clone, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Ring<F: PrimeField, KzgCurve: Pairing<ScalarField = F>, P: AffineRepr<BaseField = F>> {
    // KZG commitments to the coordinates of the vector described above
    pub cx: KzgCurve::G1Affine,
    pub cy: KzgCurve::G1Affine,
    // KZG commitment to a bitvector highlighting the part of the vector corresponding to the public keys.
    pub selector: KzgCurve::G1Affine,
    // maximal number of keys the commitment can "store". For domain of size `N` it is `N - (s + IDLE_ROWS)`
    pub max_keys: usize,
    // the number of keys "stored" in this commitment
    pub curr_keys: usize,
    // a parameter
    pub padding_point: P,
}

impl<F: PrimeField, KzgCurve: Pairing<ScalarField = F>, P: AffineRepr<BaseField = F>> fmt::Debug
    for Ring<F, KzgCurve, P>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Ring(curr_keys={}, max_keys{})",
            self.curr_keys, self.max_keys
        )
    }
}

impl<F: PrimeField, KzgCurve: Pairing<ScalarField = F>, P: AffineRepr<BaseField = F>>
    Ring<F, KzgCurve, P>
{
    // Builds the commitment to the vector
    // `padding, ..., padding, H, 2H, ..., 2^(s-1)H, 0, 0, 0, 0`.
    // We compute it as a sum of commitments of 2 vectors:
    // `padding, ..., padding`, and
    // `0, ..., 0, (H - padding), (2H - padding), ..., (2^(s-1)H  - padding), -padding, -padding, -padding, -padding`.
    // The first one is `padding * G`, the second requires an `(IDLE_ROWS + s)`-msm to compute.
    pub fn empty(
        // SNARK parameters
        piop_params: &PiopParams<F, P>,
        // Should return `srs[range]` for `range = (piop_params.keyset_part_size..domain_size)`
        srs: impl Fn(Range<usize>) -> Result<Vec<KzgCurve::G1Affine>, ()>,
        // generator used in the SRS
        g: KzgCurve::G1,
    ) -> Self {
        let padding_point = piop_params.padding_point;
        let (padding_x, padding_y) = padding_point.xy().unwrap(); // panics on inf, never happens
        let c1x = g * padding_x;
        let c1y = g * padding_y;

        let powers_of_h = piop_params.power_of_2_multiples_of_h();
        let (mut xs, mut ys): (Vec<F>, Vec<F>) = powers_of_h
            .iter()
            .map(|p| p.xy().unwrap())
            .map(|(&x, &y)| (x - padding_x, y - padding_y))
            .unzip();
        xs.resize(xs.len() + IDLE_ROWS, -*padding_x);
        ys.resize(ys.len() + IDLE_ROWS, -*padding_y);
        let domain_size = piop_params.domain.domain().size();
        let srs_segment = &srs(piop_params.keyset_part_size..domain_size).unwrap();
        let c2x = KzgCurve::G1::msm(srs_segment, &xs).unwrap();
        let c2y = KzgCurve::G1::msm(srs_segment, &ys).unwrap();

        let selector_inv = srs_segment.iter().sum::<KzgCurve::G1>();
        let selector = g - selector_inv;

        let (cx, cy, selector) = {
            let affine = KzgCurve::G1::normalize_batch(&[c1x + c2x, c1y + c2y, selector]);
            (affine[0], affine[1], affine[2])
        };

        Self {
            cx,
            cy,
            selector,
            max_keys: piop_params.keyset_part_size,
            curr_keys: 0,
            padding_point,
        }
    }

    pub fn append(
        &mut self,
        keys: &[P],
        // Should return `srs[range]` for `range = (self.curr_keys..self.curr_keys + keys.len())`
        srs: impl Fn(Range<usize>) -> Result<Vec<KzgCurve::G1Affine>, ()>,
    ) {
        let new_size = self.curr_keys + keys.len();
        assert!(new_size <= self.max_keys);
        let (padding_x, padding_y) = self.padding_point.xy().unwrap();
        let (xs, ys): (Vec<F>, Vec<F>) = keys
            .iter()
            .map(|p| p.xy().unwrap())
            .map(|(&x, &y)| (x - padding_x, y - padding_y))
            .unzip();
        let srs_segment = &srs(self.curr_keys..self.curr_keys + keys.len()).unwrap();
        let cx_delta = KzgCurve::G1::msm(srs_segment, &xs).unwrap();
        let cy_delta = KzgCurve::G1::msm(srs_segment, &ys).unwrap();

        let (new_cx, new_cy) = {
            let affine = KzgCurve::G1::normalize_batch(&[self.cx + cx_delta, self.cy + cy_delta]);
            (affine[0], affine[1])
        };

        self.cx = new_cx;
        self.cy = new_cy;
        self.curr_keys = new_size;
    }

    // Builds the ring from the keys provided with 2 MSMs of size `keys.len() + scalar_bitlen + 5`.
    // In some cases it may be beneficial to cash the empty ring, as updating it costs 2 MSMs of size `keys.len()`.
    pub fn with_keys(
        // SNARK parameters
        piop_params: &PiopParams<F, P>,
        keys: &[P],
        // full-size Lagrangian srs
        srs: &RingBuilderKey<F, KzgCurve>,
    ) -> Self {
        let padding_point = piop_params.padding_point;
        let (&padding_x, &padding_y) = padding_point.xy().unwrap(); // panics on inf, never happens
        let powers_of_h = piop_params.power_of_2_multiples_of_h();

        // Computes
        // [(pk1 - padding), ..., (pkn - padding),
        //  (H - padding), ..., (2^(s-1)HH - padding),
        //  -padding, -padding, -padding, -padding,
        //  padding].
        let (xs, ys): (Vec<F>, Vec<F>) = keys
            .iter()
            .chain(&powers_of_h)
            .map(|p| p.xy().unwrap())
            .map(|(&x, &y)| (x - padding_x, y - padding_y))
            .chain(iter::repeat((-padding_x, -padding_y)).take(4))
            .chain(iter::once((padding_x, padding_y)))
            .unzip();

        // Composes the corresponding slices of the SRS.
        let bases = [
            &srs.lis_in_g1[..keys.len()],
            &srs.lis_in_g1[piop_params.keyset_part_size..],
            &[srs.g1.into()],
        ]
        .concat();

        let cx = KzgCurve::G1::msm(&bases, &xs).unwrap();
        let cy = KzgCurve::G1::msm(&bases, &ys).unwrap();
        let selector_inv = srs.lis_in_g1[piop_params.keyset_part_size..]
            .iter()
            .sum::<KzgCurve::G1>();
        let selector = srs.g1 - selector_inv;

        let (cx, cy, selector) = {
            let affine = KzgCurve::G1::normalize_batch(&[cx, cy, selector]);
            (affine[0], affine[1], affine[2])
        };

        Self {
            cx,
            cy,
            selector,
            max_keys: piop_params.keyset_part_size,
            curr_keys: keys.len(),
            padding_point,
        }
    }

    pub fn slots_left(&self) -> usize {
        self.max_keys - self.curr_keys
    }

    pub const fn empty_unchecked(
        domain_size: usize,
        cx: KzgCurve::G1Affine,
        cy: KzgCurve::G1Affine,
        selector: KzgCurve::G1Affine,
        padding_point: P,
    ) -> Self {
        let max_keys = domain_size - (P::ScalarField::MODULUS_BIT_SIZE as usize + IDLE_ROWS);
        Self {
            cx,
            cy,
            selector,
            max_keys,
            curr_keys: 0,
            padding_point,
        }
    }
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct RingBuilderKey<F: PrimeField, KzgCurve: Pairing<ScalarField = F>> {
    // Lagrangian SRS
    pub lis_in_g1: Vec<KzgCurve::G1Affine>,
    // generator used in the SRS
    pub g1: KzgCurve::G1,
}

impl<F: PrimeField, KzgCurve: Pairing<ScalarField = F>> RingBuilderKey<F, KzgCurve> {
    pub fn from_srs(srs: &URS<KzgCurve>, domain_size: usize) -> Self {
        let g1 = srs.powers_in_g1[0].into_group();
        let ck = srs.ck_with_lagrangian(domain_size);
        let lis_in_g1 = ck.lagrangian.unwrap().lis_in_g;
        Self { lis_in_g1, g1 }
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr as Bls12_381_Fr, G1Affine};
    use ark_ed_on_bls12_381_bandersnatch::{EdwardsAffine, SWAffine};
    use ark_std::test_rng;
    use fflonk::pcs::kzg::urs::URS;
    use fflonk::pcs::kzg::KZG;
    use fflonk::pcs::PCS;

    use common::domain::Domain;
    use common::test_helpers::random_vec;

    use crate::ring::Ring;
    use crate::PiopParams;

    use super::*;

    type TestRing<P> = Ring<Bls12_381_Fr, Bls12_381, P>;

    fn _test_ring_mgmt<P: AffineRepr<BaseField = Bls12_381_Fr>>() {
        let rng = &mut test_rng();

        let domain_size = 1 << 9;

        let pcs_params = KZG::<Bls12_381>::setup(domain_size - 1, rng);
        let ring_builder_key = RingBuilderKey::from_srs(&pcs_params, domain_size);
        let srs = |range: Range<usize>| Ok(ring_builder_key.lis_in_g1[range].to_vec());

        // piop params
        let h = P::rand(rng);
        let seed = P::rand(rng);
        let pad = P::rand(rng);
        let domain = Domain::new(domain_size, true);
        let piop_params = PiopParams::setup(domain, h, seed, pad);

        let mut ring = TestRing::empty(&piop_params, srs, ring_builder_key.g1);
        let (monimial_cx, monimial_cy) = get_monomial_commitment(&pcs_params, &piop_params, &[]);
        assert_eq!(ring.cx, monimial_cx);
        assert_eq!(ring.cy, monimial_cy);

        let keys = random_vec::<P, _>(ring.max_keys, rng);
        ring.append(&keys, srs);
        let (monimial_cx, monimial_cy) = get_monomial_commitment(&pcs_params, &piop_params, &keys);
        assert_eq!(ring.cx, monimial_cx);
        assert_eq!(ring.cy, monimial_cy);

        let same_ring = TestRing::with_keys(&piop_params, &keys, &ring_builder_key);
        assert_eq!(ring, same_ring);
    }

    #[test]
    fn test_ring_mgmt_sw() {
        _test_ring_mgmt::<SWAffine>();
    }

    #[test]
    fn test_ring_mgmt_te() {
        _test_ring_mgmt::<EdwardsAffine>();
    }

    fn _test_empty_rings<P: AffineRepr<BaseField = Bls12_381_Fr>>() {
        let rng = &mut test_rng();

        let domain_size = 1 << 9;

        let pcs_params = KZG::<Bls12_381>::setup(domain_size - 1, rng);
        let ring_builder_key = RingBuilderKey::from_srs(&pcs_params, domain_size);
        let srs = |range: Range<usize>| Ok(ring_builder_key.lis_in_g1[range].to_vec());

        // piop params
        let h = P::rand(rng);
        let seed = P::rand(rng);
        let pad = P::rand(rng);
        let domain = Domain::new(domain_size, true);
        let piop_params = PiopParams::setup(domain, h, seed, pad);

        let ring = TestRing::<P>::empty(&piop_params, srs, ring_builder_key.g1);
        let same_ring = TestRing::with_keys(&piop_params, &[], &ring_builder_key);
        assert_eq!(ring, same_ring);
    }

    #[test]
    fn test_empty_rings_sw() {
        _test_empty_rings::<SWAffine>();
    }

    #[test]
    fn test_empty_rings_te() {
        _test_empty_rings::<EdwardsAffine>();
    }

    fn get_monomial_commitment<P: AffineRepr<BaseField = Bls12_381_Fr>>(
        pcs_params: &URS<Bls12_381>,
        piop_params: &PiopParams<<P as AffineRepr>::BaseField, P>,
        keys: &[P],
    ) -> (G1Affine, G1Affine) {
        let (_, verifier_key) =
            crate::piop::index::<_, KZG<Bls12_381>, _>(pcs_params, piop_params, keys);
        let [monimial_cx, monimial_cy] = verifier_key.fixed_columns_committed.points;
        (monimial_cx.0, monimial_cy.0)
    }
}
