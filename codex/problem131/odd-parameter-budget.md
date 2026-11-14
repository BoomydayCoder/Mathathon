# Changing the root multiplicity: a general budget and its compensating cost

Status: exact algebraic formulas and a proved monotonicity observation, plus a
bounded exploratory numerical experiment. General local p-adic estimates below
are extensions of the source's proof architecture and still require independent
review. Numerical equilibrium energies are **not interval certificates**. No
ζ(7) irrationality conclusion follows.

## 1. Parameters and normalization

Take an even numerator power r=2m, N/K→α, h/K→λ and form

\[
 G_{ij}(X)=\mu_{s,X}\left(\frac{D_N(t)^r t^{i+j}}{D_K(t)}\right),
 \quad 0\le i,j<h,
\]
\[
 S_K=\frac{(K!)^{2h}4^{h-1}}
 {(N!)^{2rh}\prod_{i=1}^{h-1}((2i)!)^2}.
\]

The integer-valued row basis uses (D_N/(N!)²)^m times the source's q_i.
The scalar normalization cancels the K²log K terms for arbitrary fixed
α,λ,r. Its remaining quadratic constant is

\[
 C_{r,\alpha,\lambda}=-2\lambda+2r\alpha\lambda(1-\log\alpha)
       +3\lambda^2-2\lambda^2\log(2\lambda).
\]

The positive real determinant has external field
V(x)=2π√x+∫₀¹log(x+u²)du−r∫₀^αlog(x+u²)du. Its equilibrium exists for
r>3 by [`general-equilibrium-family.md`](general-equilibrium-family.md).
The leading energy is independent of fixed odd s, but not of r,α,λ.

The full-pole outer arithmetic basis below assumes λ=1−α, so h=K−N.
The inner algebra does not need this relation. If λ changes independently,
one cannot silently reuse the source's outer basis of K−N rows. A subsequent
restriction/minor lemma now supplies a sufficient bound; see
[`outer-subspace-minors.md`](outer-subspace-minors.md). Positivity and a degree
bound at most h still hold. Exact degree h is unnecessary for the basic polynomial
irrationality criterion if one already has nonzero values and degree ≤h.

## 2. General small/inner envelope

Write H=λ+mα, b(z)=mℓ(αx,z), T=⌊2Hx⌋, ξ=Hx−T/2,
q=⌊2x⌋ and n₊=(2x−q)/2. The same allocation gives

\[
 \Gamma_{s,r}=\int_0^{1/2}(T-b)(T+b-\ell-s)\,dz
       +\xi(2T-q-s)+(\xi-n_+)_+.
\]

The scalar factorial contribution is
N_r(x)=2λx⌊x⌋−2rλx⌊αx⌋−2J(λx), with J(u)=⌊2u⌋u−⌊2u⌋(⌊2u⌋+1)/4.
Expanding shows that their quadratic terms are opposite:
Γ has coefficient 2λ(λ+rα−1), while N has its negative. Consequently

\[
 R_{s,r}(x)=xF_{s,r}(x)+Q_{s,r}(x),\qquad
 F_{s,r}=\lambda[(s-1)+2\{x\}-2r\{\alpha x\}],
\]
\[
 Q_{s,r}=\frac{\tau(\tau-\sigma)-(\tau-\sigma)_+
                    +\eta(1-\eta)}2
          +m^2\int B^2-m\int AB,
\]

where A=ℓ(x,z)−2x, B=ℓ(αx,z)−2αx and τ,σ,η are as before. Thus

\[
 -\frac{m+1}{8}\le Q_{s,r}\le\frac{m^2+m+1}{8}.
\]

For rational α, the periodic average of F is (s−r)λ. This suggests the natural
choice **r=s+1**: for ζ(7), change r from 6 to 8 so that the mean remains −λ.
This is a concrete compensating mechanism, rather than a change of notation.
However its real cost must also be included, as Section 4 proves.

Assume α<1/3, rα≤1 and the inner allocation is valid for x≥3. Exact fractional
part integration and γ>H₃−log3−1/6 give

\[
 \int_3^\infty F_{s,r}(x)x^{-2}\,dx
 \le\lambda\left[\frac{s}{3}+\frac{4r\alpha}{3}
                         +2r\alpha\log\alpha\right].
\]

For r=s+1 and the experiment below, F≥−(2r−s+1)λ. The generalized small-prime
cost is (s+1)λ/M; together with this tail bound it leaves (2r+2)λ/M.
The resulting sufficient arithmetic budget is

\[
 A_{s,r,\alpha}(M)\le I^{\rm plain}_{\rm out}(s,r,\alpha)
 +\frac{m^2+m+1}{144}
 +\lambda\left[\frac{s}{3}+\frac{4r\alpha}{3}
                         +2r\alpha\log\alpha\right]
 +\frac{(2r+2)\lambda}{M}. \tag{B}
\]

Here A denotes the bound derived from the local envelopes, not a minimal
normalizer. Formula (B) is conditional on their validity for the chosen
parameters. In particular nonnegative allocation dimensions, degree cutoffs,
and prime ranges must be checked. For much smaller λ an inner boundary at 3
may no longer suffice; it must then be moved and the outer range rederived.

## 3. Exact outer integration for the full-pole basis

For α<1/6 and λ=1−α, put x=1/y, z∈(0,1/2),
ℓ=⌊x−z⌋+⌊x+z⌋+1 and δ=1_{z<αx}. The ordinary residue cost per class is

\[
 c_{s,r}(\ell,\delta)=\sum_{0\le i<\ell-2}
               (\ell+s-1-r\delta-2i)_+.
\]

The plain polynomial correction has scaled rank
(1+(r−2)α−(s−1)y/2)₊. Thus the outer envelope is

\[
 I^{\rm plain}_{\rm out}
 =\int_{1/3}^{2\lambda}\left[
 y\int_0^{1/2}c_{s,r}(\ell,\delta)\,dz
 +(1+(r-2)\alpha-(s-1)y/2)_+
 -2\lambda\lfloor1/y\rfloor
 +\sum_{j\ge1}(2\lambda-jy)_+\right]dy. \tag{O}
\]

The z intervals change at z={x},1−{x},αx and 1/2. Their ordering changes only
when x is an integer or half-integer or x=j/(1±α). After multiplying lengths
by y, the integrand is piecewise affine in y. Along with the rank and scalar
hinge zeros, these rational breakpoints give an exact rational integral when
α is rational. The accompanying script performs this integration and checks
affineness at two additional interior points in every interval. At s=5,r=6,
α=3/40 it reproduces the separately reviewed value 129101/96000.

## 4. A rigorous compensating cost: increasing r worsens the real exponent

For fixed α and λ, let

\[
 E(r)=\sup_{\rho([0,\infty))=\lambda}
        \left[I(\rho)-\int V_r\,d\rho\right]+C_{r,\alpha,\lambda}.
\]

Combining the r-dependent field and scalar-normalization terms gives exactly

\[
 r\int\rho(dx)\int_0^\alpha\log(1+x/u^2)\,du.
\]

This quantity is strictly increasing in r for every nontrivial positive
equilibrium measure. Evaluating the larger-r variational problem at the
smaller-r optimizer proves **E(r) is strictly increasing**. Where the envelope
is differentiable, the derivative is

\[
 \frac{dE}{dr}=\int\rho_r(dx)\int_0^\alpha\log(1+x/u^2)\,du>0.
\]

There is also a finite-dimensional explanation. Increasing r multiplies the
normalized Gram weight by a power of
D_N(t)/(N!)²=∏_{j=1}^N(1+t/j²)>1 for t>0. The Gram matrix increases in
positive definite order, so its determinant increases. The factorial
normalization therefore does not make the extra positive roots analytically
free. This monotonicity is proved independently of the numerical experiment.

The large apparent αlogα terms cancel between the arithmetic gain and real
cost. For example the derivative of E plus the gamma-based inner bound in
(B), at fixed α,λ, is

\[
 \int_0^\alpha U^{\rho_r}(-u^2)\,du+\frac{10}{3}\alpha\lambda,
\]

before including outer costs and the variance bound. Thus reducing the mean
inner arithmetic cost alone is not a reliable optimization target.

## 5. Bounded exploratory experiment

The script `python codex/computations/explore_odd_parameters.py` checks the
known s=5 case, then nine s=7 choices: r=6,8,10 and α=1/20,3/40,1/10.
The outer integral is exact rational arithmetic. Endpoint roots and equilibrium
energies use ordinary SciPy floating-point quadrature and are not certificates.
The arithmetic column is the right side of (B) with the cutoff term omitted;
every finite M adds a positive correction. All displayed totals are positive.

| s | r | α | Arithmetic budget at M=∞ | Real energy, numerical | Sum |
|---|---|---|---:|---:|---:|
| 5 | 6 | 3/40 | 1.375349 | −1.385510 | −0.010161 |
| 7 | 6 | 1/20 | 2.895794 | −1.880434 | 1.015360 |
| 7 | 6 | 3/40 | 2.339151 | −1.385510 | 0.953641 |
| 7 | 6 | 1/10 | 1.902486 | −0.971358 | 0.931128 |
| 7 | 8 | 1/20 | 2.426952 | −1.406925 | 1.020027 |
| 7 | 8 | 3/40 | 1.745439 | −0.761938 | 0.983500 |
| 7 | 8 | 1/10 | 1.224944 | −0.217239 | 1.007705 |
| 7 | 10 | 1/20 | 2.016374 | −0.912542 | 1.103832 |
| 7 | 10 | 3/40 | 1.234209 | −0.100072 | 1.134137 |
| 7 | 10 | 1/10 | 0.655458 | 0.594512 | 1.249970 |

This finite test rules out none of the untested parameters or sharper estimates.
Its useful result is that the motivated r=s+1 adjustment does not immediately
repair the budget: at α=3/40 its arithmetic gain is slightly smaller than the
real loss. There is no numerical candidate here warranting an interval
certificate or an irrationality claim.

## 6. Strongest next gate

Two changes would actually escape the identified constraint:

1. Take λ<1−α and test a genuinely different local estimate. The subsequent
   `outer-subspace-minors.md` proves a sufficient bound using the worst-weight
   minors, but its bounded experiment still produces no successful budget.
   Its Hermite-interpolation lemma also proves that a worst-weight minor is
   always a unit here: the embedding pattern alone cannot improve that bound.
2. Improve ordinary pole-class valuations before determinant expansion, using
   cancellation not captured by the independent row weights. Increasing r only
   buys arithmetic improvement by paying the proved positive real cost.

A new scalar normalization without stronger valuation information cannot create
a gain: multiplication by a fixed rational scalar changes arithmetic and real
logarithmic costs by opposite amounts. A different basis can help only insofar
as it exposes arithmetic cancellation previously missed by the bounds.
