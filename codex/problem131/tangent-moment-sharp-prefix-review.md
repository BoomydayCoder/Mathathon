# Independent digit-arithmetic audit of the sharp tangent-moment prefix

Status: **accepted; no correction required** for tangent-moment-sharp-prefix.md. This review independently checks the finite-field operator and the exact digit threshold; it makes no claim about determinant denominators or irrationality.

Let s>=3 be odd, A_e=T_{2e+1} binom(s+2e,s-1), and Q=max_{p|s}p^{v_p(s)}. The claimed first index of a gcd-one prefix is exactly (Q-1)/2.

Fix p|s, h=v_p(s), P=p^h, and s=Pu with p not dividing u. Put r=2e+1. A slightly stronger carry statement is available for every 1<=r<=P:

v_p binom(s+r-1,s-1)=h-v_p(r).

For r<P, the low h digits of s-1 are all p-1. In the addition (s-1)+r, the first nonzero digit of r, at j=v_p(r), initiates a carry. Every position j,...,h-1 then carries. At position h the carry stops: r has digit zero there, while (s-1)_h=(u mod p)-1<=p-2. There are no further carries. Kummer's count is therefore h-j. Equivalently, Lucas alone gives the lower-index divisibility needed in the source proof. At r=P there is no carry, and Lucas additionally gives

binom(s+P-1,P)=u modulo p.

Thus all e<(P-1)/2 have binomial factor divisible by p. At the threshold, that factor is a unit; the tangent factor must still be checked rather than assumed.

For odd p, take i^2=-1 in F_{p^2} and delta=(1+y^2)d/dy on the rational-function field. The coordinate v=(y-i)/(y+i) has delta(v)=2iv. Hence delta=2i E with E=v d/dv. In characteristic p, E^p is a derivation: the iterated Leibniz middle terms vanish modulo p. It agrees with E on v and on the constant field, so E^p=E on the entire rational-function field. Consequently

delta^p=(2i)^{p-1}delta=(-1)^{(p-1)/2}delta.

The scalar lies in F_p. Iterating gives delta^{p^h}=(-1)^{h(p-1)/2}delta. The tangent differential equation y'=1+y^2 identifies its derivative numbers by T_r=(delta^r y)(0) for odd r, so

T_P=(-1)^{h(p-1)/2} modulo p.

This is nonzero. Evaluation is performed on the polynomial delta^P y, so the rational Möbius coordinate introduces no pole-evaluation issue. The proof specifically requires odd p; no analogous assertion at p=2 is being silently used.

Therefore each prime p|s first ceases to divide the running prefix exactly at (P-1)/2. Any prefix gcd divides A_0=s, so no other prime can occur. Since s is odd, 2 is excluded at the outset. Before (Q-1)/2 the prime attaining Q still divides every prefix entry; at that index every prime divisor of s has an explicit unit witness. The first gcd-one index is exactly (Q-1)/2. Finally Q=s iff s is an odd prime power, proving the claimed sharpness criterion for the older (s-1)/2 bound.

## Independent computations

check_tangent_sharp_prefix_review.py does not import the source checker. It generates tangent derivatives from the nonlinear ODE's binomial convolution rather than repeated delta, tests all odd s through401, verifies the *exact* binomial valuations by both integer division and digit carries through each local threshold, and independently checks delta^p on polynomial monomials and delta^{p^h}y at prime powers through401. It also checks the finite-field tangent unit and binomial witness residues separately. Finite tests supplement the proof above.
