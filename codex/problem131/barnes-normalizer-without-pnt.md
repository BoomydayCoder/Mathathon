# A superfactorial majorant removes prime-number-theorem averaging

Status: exact rational finite-profile certificate passed; the complete
normalizer and elementary growth proof below passed independent audit in
[`pnt-free-normalizer-independent-review.md`](pnt-free-normalizer-independent-review.md).
This replaces the prime-number-theorem step in Fauzan's arithmetic estimate,
conditional on the source's local valuation lemmas. It does not independently
verify the entire preprint. The certificate is elementary but sizeable; it is
not yet a short hand proof. No historical novelty is asserted.

## 1. The outcome and dependencies

There is an explicitly specified rational superfactorial product, with an
elementary small-prime repair, whose arithmetic cost is less than 1379/1000=1.379
in the K² exponent. Combined with the separately certified real bound −69/50,
this leaves a strict negative margin. No asymptotic estimate for the number
or distribution of primes is used in this construction.

The input is the original source's **refined** outer valuation envelope and
inner allocation. The recent outer-kernel improvement is not required.
`check_sparse_barnes_certificate.py` checks the rational breakpoints, the exact
quadratic cancellation, a positive periodic-extension slope, and rational
logarithm bounds. An independent check here passed all 561 necessary
breakpoints. The exploratory LP script is not part of the certificate. The
earlier 100-frequency certificate remains available as a separate fallback.

## 2. Why superfactorials fit the local floor functions

Define

\[
 J(v)=\lfloor v\rfloor v-\frac{\lfloor v\rfloor(\lfloor v\rfloor+1)}2,
 \qquad B(N)=\prod_{j=0}^{N-1}j!.
\]

The exact valuation identity is

\[
 v_p(B(N))=\sum_{a\ge1}p^a J(N/p^a). \tag{2.1}
\]

Indeed, count the contributions of each p^a to all factorials j!, j<N.
The coefficient is ∑_{k≥1}(N−kp^a)₊. The convention B(N), rather than the
product through N!, avoids an extra floor term.

For α=3/40, λ=37/40 and H=23/20, direct expansion of the source inner
allocation yields the exact identity

\[
\begin{aligned}
\Gamma(x)={}&-\frac{197}{40}x+\frac12J(23x/10)-\frac12J(2x)
 +\frac12J(3x/10)\\
 &-3J(37x/40)+3J(43x/40)-9J(3x/20). \tag{2.2}
\end{aligned}
\]

One useful derivation uses B₂(v)=v²−v+1/6 and

\[
 \int_0^{1/2} A(z)B(z)\,dz
 =\frac{B_2(\{(1-\alpha)x\})-B_2(\{(1+\alpha)x\})}{2}.
\]

The variance is [B₂(0)−B₂({2αx})]/2. The remaining allocation term is
[−B₂({2Hx})+B₂({2x})−B₂({2(H−1)x})+B₂(0)]/4.
Finally use
B₂({v})=v²−v+1/6−2J(v). Thus the local floor structure is already closely
matched to exact superfactorial valuations; the difficult part is the
pointwise majorization across the outer/inner transition.

## 3. The exact finite-profile certificate

Use the 25 nonzero frequency/coefficient pairs in
`codex/computations/sparse_barnes_certificate_100.json`. The frequencies are
multiples of 1/40 and at most 9/8. All coefficients are hundredths except the
coefficient −457/40000 at frequency one. (The JSON also includes one zero
coefficient, which can be discarded.) Put

\[
 H_*(x)=\sum_j a_j J(c_jx),\qquad
 G_2=10\alpha\lambda=111/160.
\]

The script proves

\[
 \frac12\sum_j a_jc_j^2=-G_2,\qquad
 \beta=-\frac12\sum_j a_jc_j-5\lambda
       =\frac{9377}{80000}>0. \tag{3.1}
\]

Define the required raw determinant cost profile by

\[
 T(x)=\begin{cases}
 0,&0\le x\le1,\\
 x[R_0(1/x)-d(1/x)],&1<x<3,\\
 -\Gamma(x),&x>3.
 \end{cases}
\]

At x=3 use the larger one-sided value. R₀,d are exactly the source's outer
functions, reproduced using rational arithmetic in the script. The certificate
checks H_*(x)≥T(x) at every breakpoint on [0,43]. Both sides are affine
between these breakpoints, so this proves the inequality throughout that
interval. All slopes c_j have denominator 40, and (2.2) has the same period
in its fractional parts. The quadratic equality (3.1) implies

\[
 [H_*+\Gamma](x+40)-[H_*+\Gamma](x)=40\beta>0.
\]

Thus checking [3,43] proves H_*≥−Γ for every x≥3. The finite certificate is
a genuine all-x inequality, not a numerical mesh test.

The coefficient common denominator is

\[
 q=40000,\qquad q a_j\in\mathbb Z.
\]

Raising the determinant to this fixed power preserves the degree-O(K)
irrationality criterion.

## 4. Explicit normalizer, including all repairs

Let K run through multiples of 40, N=3K/40, h=37K/40, and let Δ_K be the
source's raw positive determinant. Write S_K for its original scalar so that
F_K=S_KΔ_K. Define

\[
 T_K=\prod_j B(c_jK)^{q a_j}.
\]

Negative integer exponents are permitted: T_K is a positive rational number.
Take M=10⁷, H_K=⌊K/M⌋, C=148 and D_n=lcm(1,…,n). Set

\[
 \boxed{P_K(X)=
 \left[D_{H_K}^{3CK}\,((5K)!)^{\lceil\sqrt K\rceil}\right]^q
                   T_K\,\Delta_K(X)^q.} \tag{4.1}
\]

This is explicit and is eventually a polynomial in Z[X]. The proof follows.
Its degree is at most qh and its value at ζ(5) is positive. The extra factorial
repair is deliberately excessive but contributes only o(K²) to the logarithm.

### Large primes p>K/M

For sufficiently large K, p²>5K. Formula (2.1) then gives the exact identity
v_p(T_K)/q=pH_*(K/p). The source's elementary counting estimates, before
any appeal to the prime number theorem, give uniformly

\[
 v_p^G(\Delta_K)\ge-pT(K/p)-O_M(1).
\]

For p>K its bound is simply integrality. The pointwise certificate proves
v_p^G(T_KΔ_K^q)≥−q O_M(1). Every relevant prime is at most 5K, so the
factor ((5K)!)^{q⌈√K⌉} eventually clears these bounded deficits. Its logarithm
is O(qK^{3/2}log K)=o(qK²). This avoids an unspecified repair exponent:
⌈√K⌉ eventually exceeds the fixed uniform constant for the chosen M.

### Small primes p≤H_K, including every higher prime power

Put L_p=⌊log_p(5K)⌋. An elementary floor estimate from (2.1), valid for
0≤n≤5K and p≤5K, is

\[
 \left|v_p(B(n))-\frac{n^2}{2(p-1)}\right|\le2nL_p.
\]

To verify it, use J(v)=v²/2−v/2+{v}(1−{v})/2 for p^a≤n, and bound the
remaining geometric tail. The elementary factorial digit-sum formula also
gives |v_p(n!)−n/(p−1)|≤2L_p. Applying it to the exact S_K formula gives

\[
 \left|v_p(S_K)+\frac{G_2K^2}{p-1}\right|\le36KL_p.
\]

For the last bound, the factorial-error coefficients sum to at most 32hL_p;
the endpoint correction 2h/(p−1) and the power of four add at most 4hL_p.
Let W₁=∑|a_j|c_j. Quadratic cancellation (3.1) now gives

\[
 v_p(T_K)/q-v_p(S_K)\ge-(2W_1+36)KL_p.
\]

The source's all-prime bound v_p^G(F_K)≥−6hL_p−h v_p(24) is at least
−9KL_p. The exact coefficient list has W₁=2053577/40000 and satisfies
C=⌈2W₁+45⌉=148. Therefore

\[
 v_p^G(T_K\Delta_K^q)\ge-qCKL_p.
\]

For K≥200M², one has 5K≤H_K². If p≤H_K and r=⌊log_pH_K⌋≥1, then
L_p≤2r+1≤3r. Hence D_{H_K}^{3qCK} clears every small-prime deficit. This
argument includes all higher prime powers, without counting primes or using
a prime-density estimate.

Primes above 5K divide none of the factorial factors and the source already
proves Δ_K integral there. This completes eventual coefficientwise integrality.
Gauss valuation is multiplicative, so raising Δ_K to q creates no ambiguity.

## 5. Growth from Stirling and an elementary lcm bound

Summing Stirling's formula gives

\[
 \log B(n)=\frac{n^2}{2}\log n-\frac{3n^2}{4}+O(n\log n).
\]

Thus

\[
 \frac1q\log(T_K/S_K^q)=A K^2+O(K\log K),
\]
\[
 A=\frac12\sum_j a_jc_j^2\log c_j+\frac32G_2-C_*,
\]
\[
 C_*=-2\lambda+12\alpha\lambda(1-\log\alpha)
                  +3\lambda^2-2\lambda^2\log(2\lambda).
\]

The K²log K terms cancel exactly by (3.1). The certificate bounds each logarithm
using the positive atanh series after multiplication by a power of two, with
an explicit rational tail. The diagnostic decimal for its upper endpoint is
approximately 1.3785203557; the combined bound below is certified rationally.

Only an elementary lcm bound is needed for the repair:

\[
 \log D_n\le n\log4+O((\log n)^2).
\]

For completeness, every prime power in (⌊n/2⌋,n] contributes to
binom(n,⌊n/2⌋), apart from the possible single value (n+1)/2. Hence
ψ(n)≤ψ(⌊n/2⌋)+n log2+log n. Iteration proves the displayed estimate, with
no prime number theorem. Applying it to (4.1) gives the repair cost
3C log4/M in the quadratic exponent. The exact certificate proves

\[
 \boxed{A+\frac{3C\log4}{M}<\frac{1379}{1000}=1.379<1.38.} \tag{5.1}
\]

The diagnostic sum is approximately 1.3785819071. All other repairs and floor
errors are o(K²). Combining (5.1) with the separate real estimate for F_K
gives a negative exponent for P_K(ζ(5)), subject to its already stated analytic
and local-arithmetic dependencies. Positive integral polynomial values then
contradict rationality by the usual degree-O(K) criterion.

## 6. What has and has not been simplified

Prime-density averaging and the prime number theorem are absent. All arithmetic
growth comes from exact factorial valuations, Stirling, and the elementary
central-binomial lcm bound. The local p-adic determinant lemmas remain essential.
The price is a sizeable rational certificate and a very large fixed power q.
This is an elementary alternative with a reproducible finite verification,
not yet a claim that the whole proof is shorter or easier to read.

The sparse certificate was developed by the parent and independently checked
here. The earlier 100-frequency coefficient list in
`certify_barnes_majorant.py` also passes exact verification and gives the
stronger bound 1.372, with the larger common power 800000000 and M=10⁶.
It is retained as a fallback rather than silently discarded. The exploratory
LP additionally found that the same finite ansatz
with the plain outer bound gives about 1.399 and does not suffice; no failure
claim is made for all possible plain-bound majorants.
