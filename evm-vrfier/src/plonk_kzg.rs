use alloy::primitives::{FixedBytes, U256};
use ark_ff::fields::PrimeField;
use ark_ff::BigInteger;

alloy::sol!(
    #[sol(rpc)]
    PlonkKzg,
    "contracts/out/PlonkKzg.sol/PlonkKzg.json"
);

/// Encodes a BLS12-381 base field element (381 bits) as specified in
/// [eip-2537](https://eips.ethereum.org/EIPS/eip-2537#fine-points-and-encoding-of-base-elements):
/// > A base field element (Fp) is encoded as 64 bytes
/// > by performing the BigEndian encoding of the corresponding (unsigned) integer.
/// > Due to the size of p, the top 16 bytes are always zeroes.
pub fn bls_base_field_to_bytes(fq: ark_bls12_381::Fq) -> [FixedBytes<32>; 2] {
    let be_bytes = fq.into_bigint().to_bytes_be(); // 48 bytes
    let high_bytes = FixedBytes::left_padding_from(&be_bytes[..16]); // 16 bytes
    let low_bytes = FixedBytes::from_slice(&be_bytes[16..]); // 32 bytes
    [high_bytes, low_bytes]
}

pub fn bls_scalar_field_to_uint256(fr: ark_bls12_381::Fr) -> U256 {
    let be_bytes = fr.into_bigint().to_bytes_be();
    U256::from_be_slice(&be_bytes)
}

pub fn bls_scalar_field_to_bytes32(fr: ark_bls12_381::Fr) -> FixedBytes<32> {
    let be_bytes = fr.into_bigint().to_bytes_be();
    FixedBytes::left_padding_from(&be_bytes)
}

pub fn encode_bls_g1(p: ark_bls12_381::G1Affine) -> BLS::G1Point {
    let [x_a, x_b] = bls_base_field_to_bytes(p.x);
    let [y_a, y_b] = bls_base_field_to_bytes(p.y);
    BLS::G1Point { x_a, x_b, y_a, y_b }
}

pub fn encode_bls_g2(p: ark_bls12_381::G2Affine) -> BLS::G2Point {
    let [x_c0_a, x_c0_b] = bls_base_field_to_bytes(p.x.c0);
    let [x_c1_a, x_c1_b] = bls_base_field_to_bytes(p.x.c1);
    let [y_c0_a, y_c0_b] = bls_base_field_to_bytes(p.y.c0);
    let [y_c1_a, y_c1_b] = bls_base_field_to_bytes(p.y.c1);
    BLS::G2Point {
        x_c0_a,
        x_c0_b,
        x_c1_a,
        x_c1_b,
        y_c0_a,
        y_c0_b,
        y_c1_a,
        y_c1_b,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plonk_kzg::PlonkKzg;
    use ark_bls12_381::{Bls12_381, Fr, G2Affine};
    use ark_ec::pairing::Pairing;
    use ark_ec::{AffineRepr, PrimeGroup};
    use ark_std::rand::Rng;
    use ark_std::{test_rng, UniformRand};
    use w3f_pcs::aggregation::single::aggregate_polys;
    use w3f_pcs::pcs::kzg::params::RawKzgVerifierKey;
    use w3f_pcs::pcs::kzg::urs::URS;
    use w3f_pcs::pcs::kzg::KZG;
    use w3f_pcs::pcs::{PcsParams, PCS};
    use w3f_pcs::DenseUVPolynomial;
    use w3f_pcs::Poly;
    use w3f_pcs::Polynomial;

    struct ArksBatchKzgOpenning<E: Pairing> {
        polys_z1: Vec<E::G1Affine>,
        poly_z2: E::G1Affine,
        z1: E::ScalarField,
        z2: E::ScalarField,
        evals_at_z1: Vec<E::ScalarField>,
        eval_at_z2: E::ScalarField,
        kzg_proof_at_z1: E::G1Affine,
        kzg_proof_at_z2: E::G1Affine,
    }

    struct EthBatchKzgOpenning {
        polys_z1: Vec<BLS::G1Point>,
        poly_z2: BLS::G1Point,
        z1: U256,
        z2: U256,
        evals_at_z1: Vec<U256>,
        eval_at_z2: U256,
        kzg_proof_at_z1: BLS::G1Point,
        kzg_proof_at_z2: BLS::G1Point,
    }

    impl ArksBatchKzgOpenning<Bls12_381> {
        fn encode(self) -> EthBatchKzgOpenning {
            EthBatchKzgOpenning {
                polys_z1: self.polys_z1.into_iter().map(encode_bls_g1).collect(),
                poly_z2: encode_bls_g1(self.poly_z2),
                z1: bls_scalar_field_to_uint256(self.z1),
                z2: bls_scalar_field_to_uint256(self.z2),
                evals_at_z1: self
                    .evals_at_z1
                    .into_iter()
                    .map(bls_scalar_field_to_uint256)
                    .collect(),
                eval_at_z2: bls_scalar_field_to_uint256(self.eval_at_z2),
                kzg_proof_at_z1: encode_bls_g1(self.kzg_proof_at_z1),
                kzg_proof_at_z2: encode_bls_g1(self.kzg_proof_at_z2),
            }
        }
    }

    fn random_opening<E: Pairing, R: Rng>(
        d: usize,
        k: usize,
        rng: &mut R,
    ) -> (
        ArksBatchKzgOpenning<E>,
        Vec<E::ScalarField>,
        RawKzgVerifierKey<E>,
    ) {
        // KZG setup
        let urs = URS::from_trapdoor(
            E::ScalarField::rand(rng),
            d + 1,
            2,
            E::G1::generator(),
            E::G2::generator(),
        );
        let (ck, rvk) = (urs.ck(), urs.raw_vk());

        // Polynomials
        let polys_z1: Vec<Poly<E::ScalarField>> = (0..k)
            .map(|_| Poly::<E::ScalarField>::rand(d, rng))
            .collect();
        let poly_z2 = Poly::<E::ScalarField>::rand(d, rng);

        // Aggregate polynomial
        let nus: Vec<E::ScalarField> = (0..k).map(|_| E::ScalarField::rand(rng)).collect();
        let agg_poly_z1 = aggregate_polys(&polys_z1, &nus);

        // Evaluation points
        let z1 = E::ScalarField::rand(rng);
        let z2 = E::ScalarField::rand(rng);

        // Proofs
        let kzg_proof_at_z1 = KZG::<E>::open(&ck, &agg_poly_z1, z1).unwrap();
        let kzg_proof_at_z2 = KZG::<E>::open(&ck, &poly_z2, z2).unwrap();

        // Evaluations
        let evals_at_z1: Vec<E::ScalarField> = polys_z1.iter().map(|p| p.evaluate(&z1)).collect();
        let eval_at_z2 = poly_z2.evaluate(&z2);

        // Commitments
        let polys_z1: Vec<E::G1Affine> = polys_z1
            .iter()
            .map(|p| KZG::<E>::commit(&ck, p).unwrap().0)
            .collect();
        let poly_z2 = KZG::<E>::commit(&ck, &poly_z2).unwrap().0;

        (
            ArksBatchKzgOpenning {
                polys_z1,
                poly_z2,
                z1,
                z2,
                evals_at_z1,
                eval_at_z2,
                kzg_proof_at_z1,
                kzg_proof_at_z2,
            },
            nus,
            rvk,
        )
    }

    #[tokio::test]
    async fn test_batch_openning() -> Result<(), Box<dyn std::error::Error>> {
        let provider = alloy::providers::builder()
            .with_recommended_fillers()
            .on_anvil_with_wallet_and_config(|anvil| anvil.prague())?;

        let (test_openning, nus, kzg_vk) = random_opening::<Bls12_381, _>(123, 1, &mut test_rng());
        let test_openning = test_openning.encode();

        let plonk_kzg = PlonkKzg::deploy(&provider, encode_bls_g2(kzg_vk.tau_in_g2)).await?;

        let res = plonk_kzg
            .verify_plonk_kzg(
                test_openning.polys_z1,
                test_openning.poly_z2,
                test_openning.z1,
                test_openning.z2,
                test_openning.evals_at_z1,
                test_openning.eval_at_z2,
                test_openning.kzg_proof_at_z1,
                test_openning.kzg_proof_at_z2,
                nus.into_iter().map(bls_scalar_field_to_bytes32).collect(),
                bls_scalar_field_to_uint256(Fr::from(1)),
            )
            .call()
            .await?;
        assert!(res);

        Ok(())
    }

    #[tokio::test]
    async fn test_single_openning() -> Result<(), Box<dyn std::error::Error>> {
        let provider = alloy::providers::builder()
            .with_recommended_fillers()
            .on_anvil_with_wallet_and_config(|anvil| anvil.prague())?;

        let (test_openning, _, kzg_vk) = random_opening::<Bls12_381, _>(123, 0, &mut test_rng());
        let test_openning = test_openning.encode();

        let plonk_kzg = PlonkKzg::deploy(&provider, encode_bls_g2(kzg_vk.tau_in_g2)).await?;

        let res = plonk_kzg
            .verify(
                test_openning.poly_z2,
                test_openning.z2,
                test_openning.eval_at_z2,
                test_openning.kzg_proof_at_z2,
            )
            .call()
            .await?;
        assert!(res);

        Ok(())
    }

    #[tokio::test]
    async fn test_pairing() -> Result<(), Box<dyn std::error::Error>> {
        let provider = alloy::providers::builder()
            .with_recommended_fillers()
            .on_anvil_with_wallet_and_config(|anvil| anvil.prague())?;

        let _tau_g2 = G2Affine::generator();
        let plonk_kzg = PlonkKzg::deploy(&provider, encode_bls_g2(_tau_g2)).await?;

        let res = plonk_kzg
            .pairing2(
                encode_bls_g1(-ark_bls12_381::G1Affine::generator()),
                encode_bls_g2(ark_bls12_381::G2Affine::generator()),
                encode_bls_g1(ark_bls12_381::G1Affine::generator()),
                encode_bls_g2(ark_bls12_381::G2Affine::generator()),
            )
            .call()
            .await?;
        assert!(res);

        Ok(())
    }
}
