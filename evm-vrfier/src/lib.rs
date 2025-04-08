alloy::sol!(
    #[allow(missing_docs)]
    #[sol(rpc)]
    Counter,
    "contracts/out/Counter.sol/Counter.json"
);

#[cfg(test)]
mod tests {
    use alloy::primitives::Uint;

    #[tokio::test]
    async fn it_works() -> Result<(), Box<dyn std::error::Error>> {
        let provider = alloy::providers::builder()
            .with_recommended_fillers()
            .on_anvil_with_wallet();

        let counter = crate::Counter::deploy(&provider).await?;

        let number = counter.number().call().await?._0;
        let _ = counter.increment().send().await?;
        assert_eq!(counter.number().call().await?._0, number + Uint::from(1));

        Ok(())
    }
}
