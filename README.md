SECP256K1-CL
========

SECP256K1-CL is a fork of sipa's (Pieter Wuille) optimized 
ECDSA library for Bitcoin.

The original SECP256K1 is the fastest crypto library working
on Bitcoin's curve. On my computer (i7 3770K), it is able to verify
a signature in 0.08 ms.

It achieves this considerable speed by a combination of mathematical
optimizations and good old programming tricks.

For example, the verification involves calculating hG + xQ where G 
is the Generator point and Q is the public key. It uses the windowed
non adjacent form (WNAF) to reduce the number of multiplication required. 
For G a large array of precomputed multiple are used since G is a constant.
hG is calculated in cartesian coordinates and xQ in Jacobian coordinates because
cartesian coordinates are faster to add but slower to compute multiples, and
Jacobian coordinates are the opposite. 

For this OpenCL version, I have done the following modifications. First of all, 
verifications are batched. The latency of the GPU is very bad. It excels at doing
calculations in parallel so we need to give them in a big batch.

I split the verifications into 3 phases. Phase 1 performs the setup of h and x and finds 
their WNAF. Phase 2 computes the array of multiples of Q and calculates hG + xQ.
Phase 3 checks that the x-coord of the resulting point matches r.
Only phase 2 is run on the GPU because it is by far the most expensive operation.

Also splitting the task into CPU/GPU reduces the size of the kernels. On the GPU, we have
the group multiplication code and the field operations. The multiple precision arithmetic
library (GMP) is not required. All the math operators are run as auxilliary functions
and use local storage. However, openCL does not have a native 128 bit integer type so I had to
use the smaller size (FIELD10x26).

One big problem remains though. The original code optimizes for a single execution at a 
time so the size of the data structures are mostly irrelevant. 
Unfortunately, once we create a large batch they matter. The WNAF is sparse. It is 
guaranteed to have non zero values separated by at least the window size. For G, it's
14 and for Q it's 5. So we have big arrays that have many zeroes and they have to 
be transferred to the GPU. I believe memory bandwidth could be an issue here.

Another point of concern is the trade off between precalculation and memory usage. It
is possible that for a GPU, it is worth using smaller window sizes. This should be
profiled.

Some results: On my Nvidia GTX 560Ti, 262144 samples took
* prepare -> 2518 (ms)
* compute -> 15201
* check   -> 35
* total   -> 17755

So about 0.067733 ms per sample.

TODO:
* profile
* look at using the endormorphism code