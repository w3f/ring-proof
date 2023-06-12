use ark_ec::CurveGroup;
use ark_ec::hashing::map_to_curve_hasher::{MapToCurve, MapToCurveBasedHasher};
use ark_ff::{Field, PrimeField};
use ark_ff::field_hashers::DefaultFieldHasher;
use ark_std::rand::Rng;
use ark_std::UniformRand;
use sha2::Sha256;
use ark_ec::hashing::HashToCurve;
use ark_serialize::CanonicalSerialize;


trait Params {
    type Point: CurveGroup<ScalarField=Self::Scalar>;
    type Scalar: PrimeField;

    fn g(&self) -> Self::Point;

    fn k(&self) -> Self::Point;

    fn hash_to_curve(&self, msg: &[u8]) -> <Self::Point as CurveGroup>::Affine;

    fn hash_to_field(&self,
                     ass: &[u8],
                     inbase: &<Self::Point as CurveGroup>::Affine,
                     compk: &Self::Point,
                     preout: &Self::Point,
                     r: &Self::Point,
                     rm: &Self::Point) -> Self::Scalar;

    fn hash_to_bytes(&self,
                     inbase: &<Self::Point as CurveGroup>::Affine,
                     preout: &Self::Point) -> [u8; 32];
}

struct PedersenBases<G: CurveGroup, M: MapToCurve<G>> {
    g: G,
    k: G,
    to_curve_hasher: MapToCurveBasedHasher<G, DefaultFieldHasher<Sha256, 128>, M>,
}

impl<G: CurveGroup, M: MapToCurve<G>> PedersenBases<G, M> {
    pub fn new<R: Rng>(rng: &mut R) -> Self {
        let to_curve_hasher =
            MapToCurveBasedHasher::<G, DefaultFieldHasher<Sha256, 128>, M>::new(b"ped-vrf-test")
                .unwrap();
        Self {
            g: G::generator(),
            k: G::rand(rng),
            to_curve_hasher,
        }
    }
}

impl<G: CurveGroup, M: MapToCurve<G>> Params for PedersenBases<G, M> {
    type Point = G;
    type Scalar = G::ScalarField;

    fn g(&self) -> Self::Point {
        self.g
    }

    fn k(&self) -> Self::Point {
        self.k
    }

    fn hash_to_curve(&self, msg: &[u8]) -> <Self::Point as CurveGroup>::Affine {
        self.to_curve_hasher.hash(msg).unwrap()
    }

    fn hash_to_field(&self,
                     ass: &[u8],
                     inbase: &<Self::Point as CurveGroup>::Affine,
                     compk: &Self::Point,
                     preout: &Self::Point,
                     r: &Self::Point,
                     rm: &Self::Point) -> Self::Scalar {
        let mut t = merlin::Transcript::new(b"ped-vrf-test-fs");
        t.append_message(b"ass", ass);
        append_serializable(&mut t, b"inbase", inbase);
        append_serializable(&mut t, b"compk", compk);
        append_serializable(&mut t, b"preout", preout);
        append_serializable(&mut t, b"r", r);
        append_serializable(&mut t, b"rm", rm);
        let mut buf = [0u8; 16];
        t.challenge_bytes(b"challenge", &mut buf);
        Self::Scalar::from_random_bytes(&buf).unwrap()
    }

    fn hash_to_bytes(&self,
                     inbase: &<Self::Point as CurveGroup>::Affine,
                     preout: &Self::Point) -> [u8; 32] {
        let mut t = merlin::Transcript::new(b"ped-vrf-test-out");
        append_serializable(&mut t, b"inbase", inbase);
        append_serializable(&mut t, b"preout", preout);
        let mut buf = [0u8; 32];
        t.challenge_bytes(b"out", &mut buf);
        buf
    }
}

fn append_serializable<M: CanonicalSerialize>(t: &mut merlin::Transcript, label: &'static [u8], message: &M) {
    let mut buf = vec![0; message.uncompressed_size()];
    message.serialize_uncompressed(&mut buf).unwrap();
    t.append_message(label, &buf);
}

struct SecretKey<P: Params>(P::Scalar);

impl<P: Params> SecretKey<P> {
    pub fn new<R: Rng>(rng: &mut R) -> Self {
        let sk = P::Scalar::rand(rng);
        SecretKey(sk)
    }

    pub fn commit_key<R: Rng>(&self, params: &P, rng: &mut R) -> SigningKey<P> {
        let sk = self.0;
        let oppk = P::Scalar::rand(rng);
        let compk = params.g() * sk + params.k() * oppk;
        SigningKey { sk, oppk, compk }
    }
}

struct SigningKey<P: Params> {
    /// Long-term signer's secret key.
    sk: P::Scalar,
    /// Blinding factor used to produce the Pedersen commitment to the secret key.
    oppk: P::Scalar,
    compk: P::Point,
}

#[derive(Clone)]
struct Signature<G: CurveGroup> {
    preout: G,
    c: G::ScalarField,
    s1: G::ScalarField,
    s2: G::ScalarField,
}

impl<P: Params> SigningKey<P> {
    /// Evaluates the VRF using `inbase` point as the input, and `ass` as the associated data.
    pub fn sign<R: Rng>(&self, params: &P, input: &[u8], ass: &[u8], rng: &mut R) -> Signature<P::Point> {
        let inbase = params.hash_to_curve(input);
        let preout = inbase * self.sk;
        let (r1, r2) = (P::Scalar::rand(rng), P::Scalar::rand(rng));
        let r = params.g() * r1 + params.k() * r2;
        let rm = inbase * r1;
        let c = params.hash_to_field(ass, &inbase, &self.compk, &preout, &r, &rm);
        let s1 = r1 + c * self.sk;
        let s2 = r2 + c * self.oppk;
        Signature { preout, c, s1, s2 }
    }
}

struct VerifyingKey<P: Params> {
    compk: P::Point,
}

impl<P: Params> VerifyingKey<P> {
    pub fn verify(&self, params: &P, input: &[u8], ass: &[u8], sig: Signature<P::Point>) -> Option<[u8; 32]> {
        let Signature { preout, c, s1, s2 } = sig;
        let inbase = params.hash_to_curve(input);
        let r = params.g() * s1 + params.k() * s2 - self.compk * c;
        let rm = inbase * s1 - preout * c;
        if c == params.hash_to_field(ass, &inbase, &self.compk, &preout, &r, &rm) {
            let out = params.hash_to_bytes(&inbase, &preout);
            Some(out)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use ark_std::test_rng;
    use ark_ec::hashing::curve_maps::swu::SWUMap;
    use crate::vrf::{PedersenBases, SecretKey, VerifyingKey};

    type P = PedersenBases::<
        ark_ed_on_bls12_381_bandersnatch::SWProjective,
        SWUMap<ark_ed_on_bls12_381_bandersnatch::BandersnatchConfig>
    >;

    #[test]
    fn test_vrf() {
        let rng = &mut test_rng();
        let params = P::new(rng);
        let secret_key = SecretKey::new(rng);
        let signing_key = secret_key.commit_key(&params, rng);
        let vk = VerifyingKey { compk: signing_key.compk };
        let sig = signing_key.sign(&params, b"input", b"ass", rng);
        assert!(vk.verify(&params, b"input", b"ass", sig.clone()).is_some());
        assert!(vk.verify(&params, b"not-input", b"ass", sig.clone()).is_none());
        assert!(vk.verify(&params, b"input", b"not-ass", sig).is_none());
    }
}