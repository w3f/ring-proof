# Ring Proof

for a vector commitment to a list of public keys, and a Pedersen commitment to one of the secret keys,
  implements a zk proof of knowledge of the blinding factor for the Pedersen commitment, and the position of the
  corresponding public key in the list.

### Preliminaries 

### Plonk Prover


### Plonk Verifier

#### Plonk.Verify\
**Inputs**:\
  - $Piop$: an object of Piop type\
  - $Proof$: a proof tuple as defined in ???\
  - $Challenges: ([\alpha_1,...,\alpha_n, \zeta, [\nu_1,..,nu_n])$ A Plonk Verifier challenge defined in ???\
  - $H$: R Random oracle\
**Output**:\
  - A boolean value indicating if the $Proof$ represents a correct proof\

*****  
$C \leftarrow EvaluateConsttrain(Poip)$  
$E \leftarrow \sum_i^n \alpha[i] * C[i]$  
$D \leftarrow DomainEvaluated(Piop)$  
$q_zeta \leftarrow \frac{D}{\Omega(C + Proof.x)}$  
**return** $BatchVerify()$  
     
*****    
### Challenge  

**Definition**: *Plonk verifier challange* is defined as triple:
$$([\alpha_1,...,\alpha_n, \zeta, [\nu_1,..,nu_n])$$
where \alpha_i, zeta and nu_i are all elements of Plonk Scalar Field.



### Ring Prover



### Ring Verifier



