# Ring VRF In SNARK

We are proving that the prover only knows one sk such that 

$sk G = pk_i $ 

for an unknown i among $\{ pk_1, \dot, pk_n\}$ in the ring and that the public VRF out put $VRF_{out}$ has been generated by the same $sk$ as 

$sk H = VRF_{out}$

So it is sort of proving DLEQ in SNARK in addition to proving ring membership.