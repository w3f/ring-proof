#[cfg(test)]
mod tests {
    use alloy::primitives::U256;
    use ark_bls12_381::Fr;
    use ark_ec::twisted_edwards::TECurveConfig;
    use ark_ed_on_bls12_381_bandersnatch::{BandersnatchConfig, Fq};
    use ark_ff::{Field, One};
    use ark_std::{test_rng, UniformRand};
    use w3f_plonk_common::domain::Domain;
    use w3f_plonk_common::gadgets::ec::te_cond_add::cond_te_addition;
    use crate::plonk_kzg::bls_scalar_field_to_uint256;

    alloy::sol!(
        #[sol(rpc)]
        Constraints,
        "contracts/out/Constraints.sol/Constraints.json"
    );

    #[tokio::test]
    async fn constraints() -> Result<(), Box<dyn std::error::Error>> {
        let provider = alloy::providers::builder()
            .with_recommended_fillers()
            .on_anvil_with_wallet();
        let constraints = Constraints::deploy(&provider).await?;

        let rng = &mut test_rng();

        let z = Fr::rand(rng);
        let z_n_rust = z.pow([256]);
        let z_n_sol = constraints.mod_exp(
            bls_scalar_field_to_uint256(z),
            U256::from(256),
        ).call().await?;
        assert_eq!(bls_scalar_field_to_uint256(z_n_rust), z_n_sol);

        let te_coeff_a = BandersnatchConfig::COEFF_A;
        let b = Fr::rand(rng);
        let x1 = Fr::rand(rng);
        let y1 = Fr::rand(rng);
        let x2 = Fr::rand(rng);
        let y2 = Fr::rand(rng);
        let x3 = Fr::rand(rng);
        let y3 = Fr::rand(rng);
        let cx_rust= cond_te_addition(
            te_coeff_a,
            &b,
            &x1,
            &y1,
            &x2,
            &y2,
            x3,
            y3,
            &Fr::one(),
            Fr::one(),
        )[0];
        let cx_sol = constraints.cond_te_addition(
            bls_scalar_field_to_uint256(b),
            bls_scalar_field_to_uint256(x1),
            bls_scalar_field_to_uint256(y1),
            bls_scalar_field_to_uint256(x2),
            bls_scalar_field_to_uint256(y2),
            bls_scalar_field_to_uint256(x3),
            bls_scalar_field_to_uint256(y3),
            bls_scalar_field_to_uint256(Fr::one()),
        ).call().await?;
        assert_eq!(bls_scalar_field_to_uint256(cx_rust), cx_sol);

        Ok(())
    }

    #[tokio::test]
    async fn domain_at_zeta() -> Result<(), Box<dyn std::error::Error>> {
        let provider = alloy::providers::builder()
            .with_recommended_fillers()
            .on_anvil_with_wallet();

        let constraints = Constraints::deploy(&provider).await?;

        let rng = &mut test_rng();

        let domain = Domain::new(256, true);
        let z = Fq::rand(rng);
        let domain_at_z = domain.evaluate(z);

        let w = domain.omega();
        let w_inv = domain.omega_inv();
        let w_inv_2 = w_inv * w_inv;
        let w_inv_3 = w_inv_2 * w_inv;
        let w_inv_4 = w_inv_3 * w_inv;
        println!("\tuint256 constant w = {};", bls_scalar_field_to_uint256(w));
        println!("\tuint256 constant w_inv = {};", bls_scalar_field_to_uint256(w_inv));
        println!("\tuint256 constant w_inv_2 = {};", bls_scalar_field_to_uint256(w_inv_2));
        println!("\tuint256 constant w_inv_3 = {};", bls_scalar_field_to_uint256(w_inv_3));
        println!("\tuint256 constant w_inv_4 = {};", bls_scalar_field_to_uint256(w_inv_4));

        let z_inv = constraints.inv(
            bls_scalar_field_to_uint256(z),
        ).call().await?;
        assert_eq!(z_inv, bls_scalar_field_to_uint256(z.inverse().unwrap()));

        let v_inv_hiding_at = constraints.v_inv_hiding_at(
            bls_scalar_field_to_uint256(z),
        ).call().await?;
        assert_eq!(v_inv_hiding_at, bls_scalar_field_to_uint256(domain_at_z.vanishing_polynomial_inv));

        let not_last_row = constraints.not_last_row(
            bls_scalar_field_to_uint256(z),
        ).call().await?;
        assert_eq!(not_last_row, bls_scalar_field_to_uint256(domain_at_z.not_last_row));

        Ok(())
    }
}
