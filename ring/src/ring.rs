use ark_ec::{AffineRepr, VariableBaseMSM};
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::PrimeField;
use ark_poly::EvaluationDomain;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::fmt;
use ark_std::iter;
use ark_std::ops::{Index, Range};
use ark_std::vec::Vec;
use fflonk::pcs::kzg::urs::URS;
use fflonk::pcs::PcsParams;

use crate::PiopParams;

/// Commitment to a list of VRF public keys as is used as a public input to the ring proof SNARK verifier.

/// The VRF keys are (inner) curve points that we represent in the affine short Weierstrass coordinates.
/// We commit to the coordinate vectors independently using KZG on the outer curve. To make the commitment
/// updatable we use SRS in the Lagrangian form: `L1, ..., Ln`, where `Li = L_i(t)G`.
/// The commitment to a vector `a1, ..., an` is then `a1L1 + ... + anLn`.

/// We pad the list of keys with a `padding` point with unknown dlog up to a certain size.
/// Additionally, to make the commitment compatible with the snark,
/// we append the power-of-2 powers of the VRF blinding Pedersen base
/// `H, 2H, 4H, ..., 2^(s-1)H`, where `s` is the bitness of the VRF curve scalar field.
/// The last `4` elements are set to `(0, 0)`.

/// Thus, the vector of points we commit to coordinatewise is
/// `pk1, ..., pkn, padding, ..., padding, H, 2H, ..., 2^(s-1)H, 0, 0, 0, 0`

// `KzgCurve` -- outer curve, subgroup of a pairing-friendly curve. We instantiate it with bls12-381 G1.
// `VrfCurveConfig` -- inner curve, the curve used by the VRF, in SW form. We instantiate it with Bandersnatch.
// `F` shared scalar field of the outer and the base field of the inner curves.
#[derive(Clone, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Ring<F: PrimeField, KzgCurve: Pairing<ScalarField=F>, VrfCurveConfig: SWCurveConfig<BaseField=F>> {
    // KZG commitments to the coordinates of the vector described above
    pub cx: KzgCurve::G1,
    pub cy: KzgCurve::G1,
    // KZG commitment to a bitvector highlighting the part of the vector corresponding to the public keys.
    pub selector: KzgCurve::G1,
    // maximal number of keys the commitment can "store". For domain of size `N` it is `N-(s+4)`
    max_keys: usize,
    // the number of keys "stored" in this commitment
    pub curr_keys: usize,
    // a parameter
    padding_point: Affine<VrfCurveConfig>,
}

impl<F: PrimeField, KzgCurve: Pairing<ScalarField=F>, VrfCurveConfig: SWCurveConfig<BaseField=F>> fmt::Debug for Ring<F, KzgCurve, VrfCurveConfig> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Ring(curr_keys={}, max_keys{})", self.curr_keys, self.max_keys)
    }
}

impl<F: PrimeField, KzgCurve: Pairing<ScalarField=F>, VrfCurveConfig: SWCurveConfig<BaseField=F>> Ring<F, KzgCurve, VrfCurveConfig> {
    // Builds the commitment to the vector
    // `padding, ..., padding, H, 2H, ..., 2^(s-1)H, 0, 0, 0, 0`.
    // We compute it as a sum of commitments of 2 vectors:
    // `padding, ..., padding`, and
    // `0, ..., 0, (H - padding), (2H - padding), ..., (2^(s-1)H  - padding), -padding, -padding, -padding, -padding`.
    // The first one is `padding * G`, the second requires an `(4+s)`-msm to compute.
    pub fn empty<Srs: Index<Range<usize>, Output=[KzgCurve::G1Affine]>>(
        // SNARK parameters
        piop_params: &PiopParams<F, VrfCurveConfig>,
        // MUST contain `srs[piop_params.keyset_part_size..domain_size]`
        srs: &Srs,
        // generator used in the SRS
        g: KzgCurve::G1,
    ) -> Self {
        let padding_point = piop_params.padding_point;
        let (padding_x, padding_y) = padding_point.xy().unwrap(); // panics on inf, never happens
        let c1x = g * padding_x;
        let c1y = g * padding_y;

        let powers_of_h = piop_params.power_of_2_multiples_of_h();
        let (mut xs, mut ys): (Vec<F>, Vec<F>) = powers_of_h.iter()
            .map(|p| p.xy().unwrap())
            .map(|(&x, &y)| (x - padding_x, y - padding_y))
            .unzip();
        xs.resize(xs.len() + 4, -*padding_x);
        ys.resize(ys.len() + 4, -*padding_y);
        let domain_size = piop_params.domain.domain().size();
        let srs_segment = &srs[piop_params.keyset_part_size..domain_size];
        let c2x = KzgCurve::G1::msm(srs_segment, &xs).unwrap();
        let c2y = KzgCurve::G1::msm(srs_segment, &ys).unwrap();

        let selector_inv = srs_segment.iter().sum::<KzgCurve::G1>();
        let selector =  g - selector_inv;

        Self {
            cx: c1x + c2x,
            cy: c1y + c2y,
            selector,
            max_keys: piop_params.keyset_part_size,
            curr_keys: 0,
            padding_point,
        }
    }

    pub fn append<Srs: Index<Range<usize>, Output=[KzgCurve::G1Affine]>>(
        self,
        keys: &[Affine<VrfCurveConfig>],
        // MUST contain `srs[ring.curr_keys..ring.curr_keys + keys.len()]`
        srs: &Srs,
    ) -> Self {
        let new_size = self.curr_keys + keys.len();
        assert!(new_size <= self.max_keys);
        let (padding_x, padding_y) = self.padding_point.xy().unwrap();
        let (xs, ys): (Vec<F>, Vec<F>) = keys.iter()
            .map(|p| p.xy().unwrap())
            .map(|(&x, &y)| (x - padding_x, y - padding_y))
            .unzip();
        let srs_segment = &srs[self.curr_keys..self.curr_keys + keys.len()];
        let cx_delta = KzgCurve::G1::msm(srs_segment, &xs).unwrap();
        let cy_delta = KzgCurve::G1::msm(srs_segment, &ys).unwrap();
        Self {
            cx: self.cx + cx_delta,
            cy: self.cy + cy_delta,
            curr_keys: new_size,
            ..self
        }
    }

    // Builds the ring from the keys provided with 2 MSMs of size `keys.len() + scalar_bitlen + 5`.
    // In some cases it may be beneficial to cash the empty ring, as updating it costs 2 MSMs of size `keys.len()`.
    pub fn with_keys(
        // SNARK parameters
        piop_params: &PiopParams<F, VrfCurveConfig>,
        keys: &[Affine<VrfCurveConfig>],
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
        let (xs, ys): (Vec<F>, Vec<F>) = keys.iter()
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
        ].concat();

        let cx = KzgCurve::G1::msm(&bases, &xs).unwrap();
        let cy = KzgCurve::G1::msm(&bases, &ys).unwrap();
        let selector_inv = srs.lis_in_g1[piop_params.keyset_part_size..].iter().sum::<KzgCurve::G1>();
        let selector =  srs.g1 - selector_inv;

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
}

pub struct SrsSegment<'a, KzgCurve: Pairing> {
    slice: &'a [KzgCurve::G1Affine],
    offset: usize,
}

impl<'a, KzgCurve: Pairing> SrsSegment<'a, KzgCurve> {
    pub fn shift(slice: &'a [KzgCurve::G1Affine], offset: usize) -> Self {
        Self {
            slice,
            offset,
        }
    }
}

impl <'a, KzgCurve: Pairing> Index<Range<usize>> for SrsSegment<'a, KzgCurve> {
    type Output = [KzgCurve::G1Affine];

    fn index(&self, index: Range<usize>) -> &Self::Output {
        &self.slice[index.start - self.offset..index.end - self.offset]
    }
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct RingBuilderKey<F: PrimeField, KzgCurve: Pairing<ScalarField=F>> {
    // Lagrangian SRS
    pub lis_in_g1: Vec<KzgCurve::G1Affine>,
    // generator used in the SRS
    pub g1: KzgCurve::G1,
}

impl<F: PrimeField, KzgCurve: Pairing<ScalarField=F>> RingBuilderKey<F, KzgCurve> {
    pub fn from_srs(srs: &URS<KzgCurve>, domain_size: usize) -> Self {
        let g1 = srs.powers_in_g1[0].into_group();
        let ck = srs.ck_with_lagrangian(domain_size);
        let lis_in_g1 = ck.lagrangian.unwrap().lis_in_g;
        Self { lis_in_g1, g1 }
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr, G1Affine};
    use ark_ed_on_bls12_381_bandersnatch::{BandersnatchConfig, SWAffine};
    use ark_std::{test_rng, UniformRand, vec};
    use fflonk::pcs::kzg::KZG;
    use fflonk::pcs::kzg::urs::URS;
    use fflonk::pcs::PCS;

    use common::domain::Domain;
    use common::test_helpers::random_vec;

    use crate::PiopParams;
    use crate::ring::Ring;

    use super::*;

    #[test]
    fn test_ring_mgmt() {
        let rng = &mut test_rng();

        let domain_size_log = 9;
        let domain_size = 1 << domain_size_log;

        let pcs_params = KZG::<Bls12_381>::setup(domain_size - 1, rng);
        let ring_builder_key = RingBuilderKey::from_srs(&pcs_params, domain_size);

        // piop params
        let h = SWAffine::rand(rng);
        let seed = SWAffine::rand(rng);
        let domain = Domain::new(domain_size, true);
        let piop_params = PiopParams::setup(domain, h, seed);

        let ring = Ring::<_, Bls12_381, _>::empty(&piop_params, &ring_builder_key.lis_in_g1, ring_builder_key.g1);
        let (monimial_cx, monimial_cy) = get_monomial_commitment(pcs_params.clone(), &piop_params, vec![]);
        assert_eq!(ring.cx, monimial_cx);
        assert_eq!(ring.cy, monimial_cy);

        let srs_segment = &ring_builder_key.lis_in_g1[piop_params.keyset_part_size..domain_size];
        let srs_segment = SrsSegment::<Bls12_381>::shift(srs_segment, piop_params.keyset_part_size);
        let ring2 = Ring::<_, Bls12_381, _>::empty(&piop_params, &srs_segment, ring_builder_key.g1);
        assert_eq!(ring2.cx, ring.cx);
        assert_eq!(ring2.cy, ring.cy);

        let keys = random_vec::<SWAffine, _>(ring.max_keys, rng);
        let ring = ring.append(&keys, &ring_builder_key.lis_in_g1);
        let (monimial_cx, monimial_cy) = get_monomial_commitment(pcs_params, &piop_params, keys.clone());
        assert_eq!(ring.cx, monimial_cx);
        assert_eq!(ring.cy, monimial_cy);

        let srs_segment2 = &ring_builder_key.lis_in_g1[ring2.curr_keys..ring2.curr_keys + keys.len()];
        let srs_segment2 = SrsSegment::<Bls12_381>::shift(srs_segment2, ring2.curr_keys);
        let ring2 = ring2.append(&keys, &srs_segment2);
        assert_eq!(ring2.cx, ring.cx);
        assert_eq!(ring2.cy, ring.cy);

        let same_ring = Ring::<_, Bls12_381, _>::with_keys(&piop_params, &keys, &ring_builder_key);
        assert_eq!(ring, same_ring);
    }

    #[test]
    fn test_empty_rings() {
        let rng = &mut test_rng();

        let domain_size_log = 9;
        let domain_size = 1 << domain_size_log;

        let pcs_params = KZG::<Bls12_381>::setup(domain_size - 1, rng);
        let ring_builder_key = RingBuilderKey::from_srs(&pcs_params, domain_size);

        // piop params
        let h = SWAffine::rand(rng);
        let seed = SWAffine::rand(rng);
        let domain = Domain::new(domain_size, true);
        let piop_params = PiopParams::setup(domain, h, seed);

        let ring = Ring::<_, Bls12_381, _>::empty(&piop_params, &ring_builder_key.lis_in_g1, ring_builder_key.g1);
        let same_ring = Ring::<_, Bls12_381, _>::with_keys(&piop_params, &[], &ring_builder_key);
        assert_eq!(ring, same_ring);
    }

    fn get_monomial_commitment(pcs_params: URS<Bls12_381>, piop_params: &PiopParams<Fr, BandersnatchConfig>, keys: Vec<SWAffine>) -> (G1Affine, G1Affine) {
        let (_, verifier_key) = crate::piop::index::<_, KZG::<Bls12_381>, _>(pcs_params, &piop_params, keys);
        let [monimial_cx, monimial_cy] = verifier_key.fixed_columns_committed.points;
        (monimial_cx.0, monimial_cy.0)
    }
}