pragma solidity ^0.8.24;

import {BLS, PlonkKzg} from "src/PlonkKzg.sol";
import {Constraints} from "src/Constraints.sol";

contract Plonk {
    // The trapdoor `tau` in G2, part of the standard KZG verification key.
    BLS.G2Point tau_g2;

    constructor(BLS.G2Point memory tau_g2_) {
        tau_g2 = tau_g2_;
    }

    function verify_proof(
        BLS.G1Point[] memory columns,
        BLS.G1Point memory quotient,
        uint256 z1,
        uint256 z2,
        uint256[] memory columns_at_z1,
        uint256 column_at_z2,
        BLS.G1Point memory kzg_proof_at_z1,
        BLS.G1Point memory kzg_proof_at_z2,
        bytes32[] memory nus
    ) public view returns (bool) {
        uint256 k = columns.length;
        require(columns_at_z1.length == k);

        BLS.G1Point[] memory polys_z1 = new BLS.G1Point[](k + 1);
        for (uint256 i = 0; i < k; i++) {
            polys_z1[i] = columns[i];
        }
        polys_z1[k] = quotient;

        uint256[] memory evals_at_z1 = new uint256[](k + 1);
        for (uint256 i = 0; i < k; i++) {
            evals_at_z1[i] = columns_at_z1[i];
        }
        evals_at_z1[k] = compute_quotient(columns_at_z1, column_at_z2, z1);

        return PlonkKzg.verify_plonk_kzg(
            polys_z1, polys_z1[0], z1, z2, evals_at_z1, column_at_z2, kzg_proof_at_z1, kzg_proof_at_z2, nus, 123, tau_g2
        );
    }

    function compute_quotient(uint256[] memory columns_at_z1, uint256 column_at_z2, uint256 z)
        internal
        view
        returns (uint256)
    {
        uint256 not_last = Constraints.not_last_row(z);
        uint256 c = Constraints.cond_te_addition(
            columns_at_z1[4],
            columns_at_z1[0],
            columns_at_z1[1],
            columns_at_z1[2],
            columns_at_z1[3],
            column_at_z2,
            column_at_z2,
            not_last
        );
        return Constraints.quotient_at(c, z);
    }
}
