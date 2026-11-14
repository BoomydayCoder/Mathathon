# A standard-library-only certificate for the equilibrium energy

Status: exact interval computation completed; independent primitive and assembly
review is being recorded in rational-energy-review.md. This is a third numerical
implementation, using integer arithmetic rather than Arb. It shares the already
proved endpoint enclosure and virial identity with the earlier route. It is not
an independent mathematical reconstruction of those inputs.

The script is computations/certify_equilibrium_energy_rational.py. Run it with
ordinary Python3; it imports only the standard library and the local exact
interval primitives in certify_equilibrium_endpoints.py. No numerical package,
floating-point optimizer, numerical integration tolerance, or value of pi enters
any certificate inequality. Wall-clock timing is printed for convenience only.

## Certified statement

For the actual endpoint pair in
\[
a\in[0.00007119508,0.00007119517],\qquad
b\in[0.75394957,0.75394967],
\]
the equilibrium objective satisfies the outward-rounded enclosure
\[
\boxed{-1.385510140104\le J\le-1.385509421992
<-\frac{2771}{2000}<-\frac{69}{50}.}
\]
The corresponding integral is
\[
\int V\,d\rho\in[1.915868444229,1.915869572215].
\]
The JSON result records the exact dyadic integer endpoints and denominator2^120,
so the displayed decimals are not the stored certificate precision.

## Why the elementary numerical algorithm is rigorous

All interval arithmetic rounds outward at120 binary places. Square roots use
integer isqrt with an explicit one-integer correction for the upper endpoint.
For log x, first write x=2^k u with1<=u<=2. Then
\[
\log u=2\sum_{j=0}^{N-1}\frac{z^{2j+1}}{2j+1}+R_N,
\quad z=\frac{u-1}{u+1},\quad
0\le R_N\le\frac{2z^{2N+1}}{(2N+1)(1-z^2)}.
\]
The same formula encloses log2; interval multiplication handles negative k.
N=42 is fixed, and monotonicity encloses the logarithm of an entire input interval.

Fourth-order Taylor jets propagate interval bounds on derivatives over each
whole quadrature panel. For logarithms their coefficients follow from
\((\log f)'f=f'\), so no numerical differentiation is performed. A panel of width h
has Simpson remainder at most
\[
\frac{h^5}{2880}\sup|f^{(4)}|
=\frac{h^5}{120}\sup|[z^4]f(x+z)|.
\]
The code uses384 panels on[0,3/40],384 on[3/40,1], and192 on the transformed tail.

The tail logarithm has a removable singularity. Instead of dividing by a tiny
square and accepting cancellation, write
\[
-\frac{\log(1-z)}z=\sum_{k=0}^{14}\frac{z^k}{k+1}+R(z),
\quad0\le R(z)\le\frac{z^{15}}{16(1-z)}.
\]
The tail argument has0<=z<=(b-a)/4<1, and its extra prefactor is bounded by the
same number. The total omitted integral is therefore at most
\[
(1+6\alpha)\frac{((b-a)/4)^{16}}{16(1-(b-a)/4)}.
\]
Simpson bounds are applied to the finite polynomial expression; this uniform
remainder is added afterwards. Thus a derivative bound on an infinite series is
unnecessary, and the endpoint at0 is included without division by zero.

The remaining formulas are precisely the already proved virial/potential
reduction in fauzan-analytic-simplification.md. They consist of the three bounded
integrals above and logarithms/square roots of the endpoint variables. Outward
interval evaluation over the entire endpoint rectangle contains the value at
its actual root. The formulas need not describe an equilibrium measure at
other points in the rectangle.

## A quantitative improvement with the simpler arithmetic normalizer

Retain the new plain-rank normalizer from fauzan-simplified-proof-audit.md, whose
arithmetic exponent is at most2203747/1600000. The stronger certificate gives
\[
\frac{2771}{2000}-\frac{2203747}{1600000}
=\frac{13053}{1600000}.
\]
Therefore, with K=40n and the same degree37n integer polynomials,
\[
\limsup_{n\to\infty}\frac1{n^2}\log\widetilde Q_{40n}(\zeta(5))
\le-\frac{13053}{1000}<-13.
\]
In particular the eventual bound exp(-13n²) follows. The earlier exp(-4n²)
bound remains correct; either suffices for the rationality contradiction.
This improvement does not provide an explicit threshold n or a new irrationality
measure. All arithmetic and determinant-transfer dependencies remain as audited.

The point of this implementation is reproducibility and transparency. The
mathematical proof still uses the determinant construction, local prime estimates,
positive moment formula, potential comparison and prime number theorem. An
integer-only certificate does not make that entire argument elementary.
