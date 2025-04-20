use alloy::primitives::{FixedBytes, U256};
use ark_ff::{BigInteger, PrimeField};

pub mod plonk_kzg;
mod constraints;

/// Encodes a BLS12-381 base field element (381 bits) into 2 bytes32 as specified in
/// [eip-2537](https://eips.ethereum.org/EIPS/eip-2537#fine-points-and-encoding-of-base-elements):
/// > A base field element (Fp) is encoded as 64 bytes
/// > by performing the BigEndian encoding of the corresponding (unsigned) integer.
/// > Due to the size of p, the top 16 bytes are always zeroes.
pub fn fq_to_bytes(fq: ark_bls12_381::Fq) -> [FixedBytes<32>; 2] {
    let be_bytes = fq.into_bigint().to_bytes_be(); // 48 bytes
    let high_bytes = FixedBytes::left_padding_from(&be_bytes[..16]); // 16 bytes
    let low_bytes = FixedBytes::from_slice(&be_bytes[16..]); // 32 bytes
    [high_bytes, low_bytes]
}

pub fn fr_to_bytes(fr: ark_bls12_381::Fr) -> FixedBytes<32> {
    let be_bytes = fr.into_bigint().to_bytes_be();
    FixedBytes::left_padding_from(&be_bytes)
}

pub fn fr_to_uint(fr: ark_bls12_381::Fr) -> U256 {
    let be_bytes = fr.into_bigint().to_bytes_be();
    U256::from_be_slice(&be_bytes)
}

pub fn unit_to_fr(f: U256) -> ark_bls12_381::Fr {
    ark_bls12_381::Fr::from_le_bytes_mod_order(&f.as_le_bytes())
}

#[cfg(test)]
mod tests {
    use alloy::primitives::Uint;

    alloy::sol!(
        #[sol(rpc)]
        Counter,
        "contracts/out/Counter.sol/Counter.json"
    );

    #[tokio::test]
    async fn it_works() -> Result<(), Box<dyn std::error::Error>> {
        let provider = alloy::providers::builder()
            .with_recommended_fillers()
            .on_anvil_with_wallet();

        let counter = Counter::deploy(&provider).await?;

        let number = counter.number().call().await?;
        let _ = counter.increment().send().await?;
        assert_eq!(counter.number().call().await?, number + Uint::from(1));

        Ok(())
    }
}

