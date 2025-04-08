pragma solidity ^0.8.24;

import "./SoladyBls.sol";

contract SoladyBlsTest {
    function G1() internal pure returns (BLS.G1Point memory) {
        return BLS.G1Point(
            bytes32(uint256(31827880280837800241567138048534752271)),
            bytes32(uint256(88385725958748408079899006800036250932223001591707578097800747617502997169851)),
            bytes32(uint256(11568204302792691131076548377920244452)),
            bytes32(uint256(114417265404584670498511149331300188430316142484413708742216858159411894806497))
        );
    }

    function G2_NEG() internal pure returns (BLS.G2Point memory) {
        return BLS.G2Point(
            bytes32(uint256(3045985886519456750490515843806728273)),
            bytes32(uint256(89961632905173714226157479458612185649920463576279427516307505038263245192632)),
            bytes32(uint256(26419286191256893424348605754143887205)),
            bytes32(uint256(40446337346877272185227670183527379362551741423616556919902061939448715946878)),
            bytes32(uint256(17421388336814597573762763446246275004)),
            bytes32(uint256(82535940630695547964844822885348920226556672312706312698172214783216175252138)),
            bytes32(uint256(26554973746327433462396120515077546301)),
            bytes32(uint256(69304817850384178235384652711014277219752988873539414788182467642510429663469))
        );
    }

    uint256 constant q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001;

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
            agg_at_z = addmod(agg_at_z, mulmod(uint256(nus[i]), evals_at_z1[i], q), q);
        }
        msm_bases[i] = G1();
        msm_scalars[i] = bytes32(q - addmod(agg_at_z, mulmod(r, eval_at_z2, q), q));

        msm_bases[++i] = proof_z1;
        msm_scalars[i] = bytes32(z1);

        msm_bases[++i] = poly_z2;
        msm_scalars[i] = bytes32(r);

        msm_bases[++i] = proof_z2;
        msm_scalars[i] = bytes32(mulmod(r, z2, q));

        BLS.G1Point memory agg_acc = BLS.msm(msm_bases, msm_scalars);
        BLS.G1Point memory acc_proof = BLS.add(proof_z1, g1_mul(proof_z2, bytes32(r)));
        return verify(agg_acc, acc_proof);
    }

    function verify(BLS.G1Point memory acc, BLS.G1Point memory acc_proof) public view returns (bool) {
        BLS.G1Point[] memory g1_points = new BLS.G1Point[](2);
        BLS.G2Point[] memory g2_points = new BLS.G2Point[](2);
        g1_points[0] = acc;
        g2_points[0] = G2_NEG();
        g1_points[1] = acc_proof;
        g2_points[1] = tau_g2;
        return BLS.pairing(g1_points, g2_points);
    }

    function g1_mul(BLS.G1Point memory point, bytes32 scalar) private view returns (BLS.G1Point memory) {
        BLS.G1Point[] memory points = new BLS.G1Point[](1);
        bytes32[] memory scalars = new bytes32[](1);
        points[0] = point;
        scalars[0] = scalar;
        return BLS.msm(points, scalars);
    }
}
