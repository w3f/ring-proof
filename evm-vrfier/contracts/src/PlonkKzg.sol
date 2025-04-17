pragma solidity ^0.8.24;

import "./BlsGenerators.sol";

library PlonkKzg {
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
        uint256 r,
        BLS.G2Point memory tau_g2
    ) internal view returns (bool) {
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
        return verify_acc(agg_acc, acc_proof, tau_g2);
    }

    function verify(BLS.G1Point memory c, uint256 z, uint256 v, BLS.G1Point memory proof, BLS.G2Point memory tau_g2)
        internal
        view
        returns (bool)
    {
        bytes32[] memory msm_scalars = new bytes32[](2);
        BLS.G1Point[] memory msm_bases = new BLS.G1Point[](2);
        msm_scalars[0] = bytes32(z);
        msm_bases[0] = proof;
        msm_scalars[1] = bytes32(BlsGenerators.q - v);
        msm_bases[1] = BlsGenerators.G1();
        BLS.G1Point memory acc = BLS.msm(msm_bases, msm_scalars);
        acc = BLS.add(acc, c);
        return verify_acc(acc, proof, tau_g2);
    }

    function verify_acc(BLS.G1Point memory acc, BLS.G1Point memory acc_proof, BLS.G2Point memory tau_g2)
        internal
        view
        returns (bool)
    {
        return pairing2(acc, BlsGenerators.G2_NEG(), acc_proof, tau_g2);
    }

    function pairing2(
        BLS.G1Point memory g1_1,
        BLS.G2Point memory g2_1,
        BLS.G1Point memory g1_2,
        BLS.G2Point memory g2_2
    ) internal view returns (bool result) {
        BLS.G1Point[] memory g1_points = new BLS.G1Point[](2);
        BLS.G2Point[] memory g2_points = new BLS.G2Point[](2);
        g1_points[0] = g1_1;
        g2_points[0] = g2_1;
        g1_points[1] = g1_2;
        g2_points[1] = g2_2;
        return BLS.pairing(g1_points, g2_points);
    }
}
