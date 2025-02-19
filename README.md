### The code has not been audited. Use on your own risk.

## Contents

* [`w3f-plonk-common`](w3f-plonk-common) provides infrastructure for creating plonkish proofs.
* [`w3f-ring-proof`](w3f-ring-proof) for a vector commitment to a list of public keys, and a Pedersen commitment to one of the secret keys,
  implements a zk proof of knowledge of the blinding factor for the Pedersen commitment, and the position of the
  corresponding public key in the list.
