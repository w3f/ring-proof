pragma solidity ^0.8.24;

library Constraints {
    uint256 constant r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001;
    uint256 constant te_coeff_a = r - 5;

    function add(uint256 a, uint256 b) internal pure returns (uint256 c) {
        c = addmod(a, b, r);
    }

    function mul(uint256 a, uint256 b) internal pure returns (uint256 c) {
        c = mulmod(a, b, r);
    }

    function cond_te_addition(
        uint256 b,
        uint256 x1,
        uint256 y1,
        uint256 x2,
        uint256 y2,
        uint256 x3,
        uint256 y3,
        uint256 not_last
    ) public pure returns (uint256 cx) {
        /// `cx = {[(a.x1.x2 + y1.y2).x3 - x1.y1 - x2.y2].b + (x3 - x1).(1 - b)}.not_last`
        uint256 x1y1 = mul(x1, y1);
        uint256 x2y2 = mul(x2, y2);
        uint256 lx = mul(add(mul(add(mul(te_coeff_a, mul(x1, x2)), mul(y1, y2)), x3), r - add(x1y1, x2y2)), b);
        cx = add(lx, mul(add(x3, r - x1), add(1, r - b)));
    }

//    function not_last_row(unit256 domain_size) internal pure returns uint256 {
//
//    }
}
