pragma solidity ^0.8.24;

import {Test, console} from "forge-std/Test.sol";
import "../src/SoladyBls.sol";
import "../src/BlsGenerators.sol";
import "../src/PlonkKzg.sol";

contract PlonkKzgTest is Test {
    PlonkKzg kzg;
    BLS.G1Point[] srs_g1;
    BLS.G2Point tau_g2;

    function setUp() public {
        uint256 n = 3;
        bytes32 tau = bytes32(uint256(123));
        srs_g1.push(BlsGenerators.G1());
        for (uint256 i = 1; i < n; i++) {
            srs_g1.push(BlsGenerators.g1_mul(srs_g1[i - 1], tau));
        }
        tau_g2 = BlsGenerators.g2_mul(BlsGenerators.G2(), tau);
        kzg = new PlonkKzg(tau_g2);
    }

    function commit(bytes32[] memory coeffs) internal view returns (BLS.G1Point memory c) {
        assert(coeffs.length == srs_g1.length);
        c = BLS.msm(srs_g1, coeffs);
    }

    function div_by_x_minus_z(bytes32[] memory coeffs, uint256 z) internal pure returns (bytes32[] memory res_coeffs) {
        uint256 l = coeffs.length;
        res_coeffs = new bytes32[](l);
        res_coeffs[l - 1] = 0;
        res_coeffs[l - 2] = bytes32(coeffs[l - 1]);
        for (uint256 j = l - 2; j > 0; j--) {
            res_coeffs[j - 1] =
                bytes32(BlsGenerators.add_fr(BlsGenerators.mul_fr(uint256(res_coeffs[j]), z), uint256(coeffs[j])));
        }
    }

    function prove(bytes32[] memory coeffs, uint256 z) internal view returns (BLS.G1Point memory proof) {
        bytes32[] memory q = div_by_x_minus_z(coeffs, z);
        proof = commit(q);
    }

    function eval(bytes32[] memory coeffs, uint256 z) internal pure returns (uint256 v) {
        uint256 l = coeffs.length;
        v = uint256(coeffs[l - 1]);
        for (uint256 j = l - 1; j > 0; j--) {
            v = BlsGenerators.add_fr(BlsGenerators.mul_fr(v, z), uint256(coeffs[j - 1]));
        }
    }

    function verify_single(BLS.G1Point memory c, uint256 z, uint256 v, BLS.G1Point memory proof)
        internal
        view
        returns (bool)
    {
        return kzg.verify(c, z, v, proof);
    }

    function test_div() public pure {
        // 3x^2 + 2x + 1 = [1, 2, 3]
        bytes32[] memory poly = new bytes32[](3);
        for (uint256 i = 0; i < poly.length; i++) {
            poly[i] = bytes32(i + 1);
        }
        bytes32[] memory res = div_by_x_minus_z(poly, 1);
        // 3x + 5 = (3x^2 + 2x + 1) // (x - 1) = [5, 3, 0]
        assertEq(res.length, 3);
        assertEq(res[0], bytes32(uint256(5)));
        assertEq(res[1], bytes32(uint256(3)));
        assertEq(res[2], bytes32(uint256(0)));
    }

    function test_kzg() public view {
        bytes32[] memory poly = new bytes32[](3);
        for (uint256 i = 0; i < poly.length; i++) {
            poly[i] = bytes32(i + 10);
        }
        BLS.G1Point memory c = commit(poly);
        uint256 z = 777;
        BLS.G1Point memory proof = prove(poly, z);
        uint256 v = eval(poly, z);
        assert(kzg.verify(c, z, v, proof));
    }

    function test_plonk_kzg() public view {
        bytes32[] memory poly = new bytes32[](3);
        for (uint256 i = 0; i < poly.length; i++) {
            poly[i] = bytes32(i + 10);
        }
        BLS.G1Point memory c = commit(poly);
        uint256 z1 = 666;
        uint256 z2 = 777;
        BLS.G1Point memory proof_z1 = prove(poly, z1);
        BLS.G1Point memory proof_z2 = prove(poly, z2);
        uint256 eval_at_z1 = eval(poly, z1);
        uint256 eval_at_z2 = eval(poly, z2);

        assert(kzg.verify(c, z1, eval_at_z1, proof_z1));
        assert(kzg.verify(c, z2, eval_at_z2, proof_z2));

        BLS.G1Point[] memory polys_z1 = new BLS.G1Point[](1);
        uint256[] memory evals_at_z1 = new uint256[](1);
        polys_z1[0] = c;
        evals_at_z1[0] = eval_at_z1;

        BLS.G1Point memory poly_z2 = c;

        bytes32 nu = bytes32(uint256(123));
        uint256 r = 345;
        bytes32[] memory nus = new bytes32[](1);
        nus[0] = nu;
        proof_z1 = BlsGenerators.g1_mul(proof_z1, nu);

        assert(kzg.verify_plonk_kzg(polys_z1, poly_z2, z1, z2, evals_at_z1, eval_at_z2, proof_z1, proof_z2, nus, r));
    }
}
