use ark_ec::{AffineRepr, VariableBaseMSM};
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::PrimeField;
use ark_std::vec::Vec;

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
pub struct Ring<F: PrimeField, KzgCurve: Pairing<ScalarField=F>, VrfCurveConfig: SWCurveConfig<BaseField=F>> {
    // KZG commitments to the coordinates of the vector described above
    pub cx: KzgCurve::G1,
    pub cy: KzgCurve::G1,
    // maximal number of keys the commitment can "store". For domain of size `N` it is `N-(s+4)`
    max_keys: usize,
    // the number of keys "stored" in this commitment
    curr_keys: usize,
    // a parameter
    padding_point: Affine<VrfCurveConfig>,
}

impl<F: PrimeField, KzgCurve: Pairing<ScalarField=F>, VrfCurveConfig: SWCurveConfig<BaseField=F>> Ring<F, KzgCurve, VrfCurveConfig> {
    // Builds the commitment to the vector
    // `padding, ..., padding, H, 2H, ..., 2^(s-1)H, 0, 0, 0, 0`.
    // We compute it as a sum of commitments of 2 vectors:
    // `padding, ..., padding`, and
    // `0, ..., 0, (H - padding), (2H - padding), ..., (2^(s-1)H  - padding), -padding, -padding, -padding, -padding`.
    // The first one is `padding * G`, the second requires an `(4+s)`-msm to compute.
    pub fn empty(
        // SNARK parameters
        piop_params: &PiopParams<F, VrfCurveConfig>,
        // MUST be set to `SRS[piop_params.keyset_part_size..]`
        srs_segment: &[KzgCurve::G1Affine],
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
        let c2x = KzgCurve::G1::msm(srs_segment, &xs).unwrap();
        let c2y = KzgCurve::G1::msm(srs_segment, &ys).unwrap();

        Self {
            cx: c1x + c2x,
            cy: c1y + c2y,
            max_keys: piop_params.keyset_part_size,
            curr_keys: 0,
            padding_point,
        }
    }

    pub fn append(
        self,
        keys: &[Affine<VrfCurveConfig>],
        // MUST be set to `SRS[ring.curr_keys..ring.curr_keys + keys.len()]`
        srs_segment: &[KzgCurve::G1Affine],
    ) -> Self {
        let new_size = self.curr_keys + keys.len();
        assert!(new_size <= self.max_keys);
        assert_eq!(keys.len(), srs_segment.len());
        let (padding_x, padding_y) = self.padding_point.xy().unwrap();
        let (xs, ys): (Vec<F>, Vec<F>) = keys.iter()
            .map(|p| p.xy().unwrap())
            .map(|(&x, &y)| (x - padding_x, y - padding_y))
            .unzip();
        let cx_delta = KzgCurve::G1::msm(srs_segment, &xs).unwrap();
        let cy_delta = KzgCurve::G1::msm(srs_segment, &ys).unwrap();
        Self {
            cx: self.cx + cx_delta,
            cy: self.cy + cy_delta,
            curr_keys: new_size,
            ..self
        }
    }

    pub fn slots_left(&self) -> usize {
        self.max_keys - self.curr_keys
    }
}


#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr, G1Affine};
    use ark_ec::AffineRepr;
    use ark_ed_on_bls12_381_bandersnatch::{BandersnatchConfig, SWAffine};
    use ark_std::{test_rng, UniformRand, vec};
    use fflonk::pcs::{PCS, PcsParams};
    use fflonk::pcs::kzg::KZG;
    use fflonk::pcs::kzg::urs::URS;

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
        let ck = pcs_params.ck_with_lagrangian(domain_size);
        let g1 = pcs_params.powers_in_g1[0].into_group();
        let lagrangian_srs = ck.lagrangian.unwrap().lis_in_g;

        // piop params
        let h = SWAffine::rand(rng);
        let seed = SWAffine::rand(rng);
        let domain = Domain::new(domain_size, true);
        let piop_params = PiopParams::setup(domain, h, seed);

        let ring = Ring::<_, Bls12_381, _>::empty(&piop_params, &lagrangian_srs[piop_params.keyset_part_size..], g1);
        let (monimial_cx, monimial_cy) = get_monomial_commitment(pcs_params.clone(), &piop_params, vec![]);
        assert_eq!(ring.cx, monimial_cx);
        assert_eq!(ring.cy, monimial_cy);

        let keys = random_vec::<SWAffine, _>(ring.max_keys, rng);
        let srs_segment = &lagrangian_srs[ring.curr_keys..ring.curr_keys + keys.len()];
        let ring = ring.append(&keys, srs_segment);
        let (monimial_cx, monimial_cy) = get_monomial_commitment(pcs_params, &piop_params, keys);
        assert_eq!(ring.cx, monimial_cx);
        assert_eq!(ring.cy, monimial_cy);
    }

    fn get_monomial_commitment(pcs_params: URS<Bls12_381>, piop_params: &PiopParams<Fr, BandersnatchConfig>, keys: Vec<SWAffine>) -> (G1Affine, G1Affine) {
        let (_, verifier_key) = crate::piop::index::<_, KZG::<Bls12_381>, _>(pcs_params, &piop_params, keys);
        let [monimial_cx, monimial_cy] = verifier_key.fixed_columns_committed.points;
        (monimial_cx.0, monimial_cy.0)
    }
}