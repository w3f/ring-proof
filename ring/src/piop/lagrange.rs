// use ark_ec::pairing::Pairing;
// use ark_ec::VariableBaseMSM;
// use ark_poly::Evaluations;
// use ark_std::iterable::Iterable;
//
//
// #[cfg(test)]
// mod tests {
//     use ark_bls12_381::{Bls12_381, Fr, G1Affine};
//     use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
//     use ark_ec::mnt6::G1Projective;
//     use ark_poly::{Evaluations, GeneralEvaluationDomain};
//     use ark_std::test_rng;
//     use fflonk::pcs::kzg::KZG;
//     use fflonk::pcs::kzg::urs::URS;
//     use fflonk::pcs::{PCS, PcsParams};
//     use ark_poly::EvaluationDomain;
//     use ark_ff::UniformRand;
//     use ark_poly::univariate::DensePolynomial;
//     use ark_poly::DenseUVPolynomial;
//
//
//     #[test]
//     fn lag() {
//         let rng = &mut test_rng();
//         let n = 2usize.pow(10);
//         let urs = URS::<Bls12_381>::generate(n, 0, &mut test_rng());
//         let tau = Fr::rand(rng);
//         let monomials = urs.powers_in_g1.clone();
//         assert_eq!(monomials.len(), n);
//         let g1 = monomials[0];
//         assert_eq!(g1 * tau, monomials[1]);
//         let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
//         let evals = domain.evaluate_all_lagrange_coefficients(tau);
//         let base1: Vec<_> = evals.iter().map(|x| g1 * x).collect();
//         let proj: Vec<_> = monomials.iter().map(|p| p.into_group()).collect();
//         let base2 = domain.ifft(&proj);
//         assert_eq!(base1, base2);
//
//         let poly = DensePolynomial::rand(n - 1, rng);
//         let lag = poly.evaluate_over_domain_by_ref(domain);
//
//         let mono_c = KZG::<Bls12_381>::commit(&urs.ck(), &poly).0;
//
//         let base = ark_bls12_381::G1Projective::normalize_batch(&base1);
//         let lag_c: ark_bls12_381::G1Projective = VariableBaseMSM::msm(
//             &base,
//             &lag.evals,
//         );
//
//         assert_eq!(mono_c.into_group(), lag_c);
//     }
// }
