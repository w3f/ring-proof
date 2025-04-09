pragma solidity ^0.8.24;

import {Test, console} from "forge-std/Test.sol";
import "../src/SoladyBls.sol";
import "../src/BlsGenerators.sol";

contract SoladyBlsTest is Test {
    function test_pairing() public view {
        BLS.G1Point[] memory g1_points = new BLS.G1Point[](2);
        BLS.G2Point[] memory g2_points = new BLS.G2Point[](2);
        g1_points[0] = BlsGenerators.G1();
        g2_points[0] = BlsGenerators.G2_NEG();
        g1_points[1] = BlsGenerators.G1();
        g2_points[1] = BlsGenerators.G2();
        assertEq(BLS.pairing(g1_points, g2_points), true);
    }
}
