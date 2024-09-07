pub use common::transcript::PlonkTranscript;

#[cfg(feature = "ark-transcript")]
pub use ark::Transcript as ArkTranscript;

#[cfg(feature = "merlin-transcript")]
pub use mer::Transcript as MerlinTranscript;

#[cfg(test)]
pub use test::Transcript as TestTranscript;

#[cfg(feature = "ark-transcript")]
mod ark {
    use super::*;
    use ark_ff::PrimeField;
    use ark_serialize::CanonicalSerialize;
    use ark_std::rand::RngCore;
    use fflonk::pcs::PCS;
    
    #[derive(Clone)]
    pub struct Transcript(ark_transcript::Transcript);

    impl Transcript {
        pub fn new(label: &'static [u8]) -> Self {
            Self(ark_transcript::Transcript::new_labeled(label))
        }
    }

    impl<F: PrimeField, CS: PCS<F>> PlonkTranscript<F, CS> for Transcript {
        fn _128_bit_point(&mut self, label: &'static [u8]) -> F {
            self.0.challenge(label).read_reduce()
        }

        fn _add_serializable(&mut self, label: &'static [u8], message: &impl CanonicalSerialize) {
            self.0.label(label);
            self.0.append(message);
        }

        fn to_rng(mut self) -> impl RngCore {
            self.0.challenge(b"transcript_rng")
        }
    }
}

#[cfg(feature = "merlin-transcript")]
mod mer {
    use super::*;
    use ark_ff::PrimeField;
    use ark_serialize::CanonicalSerialize;
    use ark_std::rand::{RngCore, SeedableRng};
    use fflonk::pcs::PCS;
    
    #[derive(Clone)]
    pub struct Transcript(merlin::Transcript);

    impl Transcript {
        pub fn new(label: &'static [u8]) -> Self {
            Self(merlin::Transcript::new(label))
        }
    }

    impl<F: PrimeField, CS: PCS<F>> PlonkTranscript<F, CS> for Transcript {
        fn _128_bit_point(&mut self, label: &'static [u8]) -> F {
            let mut buf = [0u8; 16];
            self.0.challenge_bytes(label, &mut buf);
            F::from_random_bytes(&buf).unwrap()
        }

        fn _add_serializable(&mut self, label: &'static [u8], message: &impl CanonicalSerialize) {
            let mut buf = vec![0; message.uncompressed_size()];
            message.serialize_uncompressed(&mut buf).unwrap();
            self.0.append_message(label, &buf);
        }

        fn to_rng(mut self) -> impl RngCore {
            let mut buf = [0u8; 32];
            self.0.challenge_bytes(b"transcript_rng", &mut buf);
            rand_chacha::ChaCha20Rng::from_seed(buf)
        }
    }
}

// Simple transcript used by tests
#[cfg(test)]
mod test {
    use super::*;
    use ark_ff::PrimeField;
    use ark_serialize::CanonicalSerialize;
    use ark_std::rand::{RngCore, SeedableRng};
    use blake2::Digest;
    use fflonk::pcs::PCS;
    
    #[derive(Clone)]
    pub struct Transcript([u8; 32]);

    impl Transcript {
        pub fn new(label: &'static [u8]) -> Self {
            Self(blake2::Blake2b::digest(label).into())
        }

        pub fn raw_append(&mut self, raw: &[u8]) {
            let len = (raw.len() as u32).to_le_bytes();
            self.0 = blake2::Blake2b::digest([&self.0[..], raw, &len[..]].concat()).into();
        }
    }

    impl<F: PrimeField, CS: PCS<F>> PlonkTranscript<F, CS> for Transcript {
        fn _128_bit_point(&mut self, label: &'static [u8]) -> F {
            self.raw_append(label);
            let mut rng = rand_chacha::ChaCha20Rng::from_seed(self.0);
            F::rand(&mut rng)
        }

        fn _add_serializable(&mut self, label: &'static [u8], message: &impl CanonicalSerialize) {
            let mut buf = vec![0; message.uncompressed_size()];
            message.serialize_uncompressed(&mut buf).unwrap();
            self.raw_append(label);
            self.raw_append(&buf);
        }

        fn to_rng(self) -> impl RngCore {
            rand_chacha::ChaCha20Rng::from_seed(self.0)
        }
    }
}
