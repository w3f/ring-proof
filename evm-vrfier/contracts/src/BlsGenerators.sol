pragma solidity ^0.8.24;

import "./SoladyBls.sol";

library BlsGenerators {
    uint256 constant q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001;

    function add_fr(uint256 a, uint256 b) internal pure returns (uint256 c) {
        c = addmod(a, b, q);
    }

    function mul_fr(uint256 a, uint256 b) internal pure returns (uint256 c) {
        c = mulmod(a, b, q);
    }

    function G1() internal pure returns (BLS.G1Point memory) {
        return BLS.G1Point(
            bytes32(uint256(31827880280837800241567138048534752271)),
            bytes32(uint256(88385725958748408079899006800036250932223001591707578097800747617502997169851)),
            bytes32(uint256(11568204302792691131076548377920244452)),
            bytes32(uint256(114417265404584670498511149331300188430316142484413708742216858159411894806497))
        );
    }

    function G2() internal pure returns (BLS.G2Point memory) {
        return BLS.G2Point(
            bytes32(uint256(3045985886519456750490515843806728273)),
            bytes32(uint256(89961632905173714226157479458612185649920463576279427516307505038263245192632)),
            bytes32(uint256(26419286191256893424348605754143887205)),
            bytes32(uint256(40446337346877272185227670183527379362551741423616556919902061939448715946878)),
            bytes32(uint256(17144095208600308495026432580569150746)),
            bytes32(uint256(78698209480990513415779284580404715789311803115477290401294577488850054555649)),
            bytes32(uint256(8010509799087472606393075511737879449)),
            bytes32(uint256(91929332261301883145239454754739358796115486554644188311284324629555800144318))
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

    function g1_mul(BLS.G1Point memory point, bytes32 scalar) internal view returns (BLS.G1Point memory) {
        BLS.G1Point[] memory points = new BLS.G1Point[](1);
        bytes32[] memory scalars = new bytes32[](1);
        points[0] = point;
        scalars[0] = scalar;
        return BLS.msm(points, scalars);
    }

    function g2_mul(BLS.G2Point memory point, bytes32 scalar) internal view returns (BLS.G2Point memory) {
        BLS.G2Point[] memory points = new BLS.G2Point[](1);
        bytes32[] memory scalars = new bytes32[](1);
        points[0] = point;
        scalars[0] = scalar;
        return BLS.msm(points, scalars);
    }
}
