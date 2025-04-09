pragma solidity ^0.8.24;

import "./SoladyBls.sol";
import "./BlsGenerators.sol";

contract PlonkKzg {
    // The trapdoor `tau` in G2, part of the standard KZG verification key.
    BLS.G2Point tau_g2;

    constructor(BLS.G2Point memory tau_g2_) {
        tau_g2 = tau_g2_;
    }

    function verify_plonk_kzg(
        BLS.G1Point[] memory polys_z1,
        BLS.G1Point memory poly_z2,
        uint256 z1,
        uint256 z2,
        uint256[] memory evals_at_z1,
        uint256 eval_at_z2,
        BLS.G1Point memory proof_z1,
        BLS.G1Point memory proof_z2,
        bytes32[] memory nus,
        uint256 r
    ) public view returns (bool) {
        assert(evals_at_z1.length == polys_z1.length);
        assert(nus.length == polys_z1.length);

        uint256 n_bases = polys_z1.length + 4;

        BLS.G1Point[] memory msm_bases = new BLS.G1Point[](n_bases);
        bytes32[] memory msm_scalars = new bytes32[](n_bases);

        uint256 i;
        for (i = 0; i < polys_z1.length; i++) {
            msm_bases[i] = polys_z1[i];
        }

        for (i = 0; i < polys_z1.length; i++) {
            msm_scalars[i] = nus[i];
        }

        uint256 agg_at_z = 0;
        for (i = 0; i < polys_z1.length; i++) {
            agg_at_z = BlsGenerators.add_fr(agg_at_z, BlsGenerators.mul_fr(uint256(nus[i]), evals_at_z1[i]));
        }
        msm_bases[i] = BlsGenerators.G1();
        msm_scalars[i] = bytes32(BlsGenerators.q - BlsGenerators.add_fr(agg_at_z, BlsGenerators.mul_fr(r, eval_at_z2)));

        msm_bases[++i] = proof_z1;
        msm_scalars[i] = bytes32(z1);

        msm_bases[++i] = poly_z2;
        msm_scalars[i] = bytes32(r);

        msm_bases[++i] = proof_z2;
        msm_scalars[i] = bytes32(BlsGenerators.mul_fr(r, z2));

        BLS.G1Point memory agg_acc = BLS.msm(msm_bases, msm_scalars);
        BLS.G1Point memory acc_proof = BLS.add(proof_z1, BlsGenerators.g1_mul(proof_z2, bytes32(r)));
        return verify(agg_acc, acc_proof);
    }

    function verify(BLS.G1Point memory acc, BLS.G1Point memory acc_proof) public view returns (bool) {
        BLS.G1Point[] memory g1_points = new BLS.G1Point[](2);
        BLS.G2Point[] memory g2_points = new BLS.G2Point[](2);
        g1_points[0] = acc;
        g2_points[0] = BlsGenerators.G2_NEG();
        g1_points[1] = acc_proof;
        g2_points[1] = tau_g2;
        return BLS.pairing(g1_points, g2_points);
    }
}
