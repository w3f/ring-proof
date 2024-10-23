use ark_ec::{AffineRepr, CurveGroup, Group};
use ark_ff::PrimeField;
use ark_std::rand::Rng;
use ark_std::vec::Vec;
use ark_std::UniformRand;

pub fn random_bitvec<R: Rng>(n: usize, density: f64, rng: &mut R) -> Vec<bool> {
    (0..n).map(|_| rng.gen_bool(density)).collect()
}

pub fn random_vec<X: UniformRand, R: Rng>(n: usize, rng: &mut R) -> Vec<X> {
    (0..n).map(|_| X::rand(rng)).collect()
}

pub fn cond_sum<P: AffineRepr>(bitmask: &[bool], points: &[P]) -> P {
    assert_eq!(bitmask.len(), points.len());
    bitmask
        .iter()
        .zip(points.iter())
        .map(|(&b, &p)| if b { p } else { P::zero() })
        .sum::<P::Group>()
        .into_affine()
}

pub fn power_of_two_multiple<P>(point: P, power: usize) -> P
where
    P: AffineRepr,
{
    let mut point_multiple = point.into_group();
    for _ in 1..power {
        point_multiple.double_in_place();
    }

    point_multiple.into()
}

pub fn find_random_point<F: PrimeField, P: AffineRepr<BaseField = F>>() -> P {
    let mut x: u8 = 0;
    loop {
        let p = P::from_random_bytes(&[x]);
        if let Some(p) = p {
            let p = p.clear_cofactor();
            if !p.is_zero() {
                return p
            }
        }
        x = x + 1;
    }
}
