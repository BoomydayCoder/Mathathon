# Classifying the local annihilator and lifting integer-parity kernels

Status: proved abstract parity lemma and polynomial annihilator identities;
complete finite floor-pattern rank check on 3<x<10 for s=7,r=6. The global
finite-K application remains in the separately reviewed proposal
[`inner-leading-block-cancellation.md`](inner-leading-block-cancellation.md).
This note does not assert a new integral normalizer.

## 1. A standalone parity-weight lemma

Let π²=p in a ramified quadratic extension of Q_p. Suppose a symmetric matrix
C has entries indexed by two sets P,Q, with C_PP and C_QQ integral over Q_p
and C_PQ contained in πZ_p. The same statement allows rational functions in
X with Gauss valuation. Suppose a d-dimensional subspace of the reduction
of C_PP modulo p lies in its kernel, and its basis lifts to a primitive
rank-d submodule over Z_p, independent of X when X is present.

Then

\[
 v_p(\det C)\ge d.
\]

**Proof.** Extend the lifted kernel vectors to a unimodular basis of the P
coordinates and perform the same change on rows and columns. A kernel row
has p-integral entries divisible by p against all P columns. Against Q its
entries are divisible by π. In particular kernel–kernel entries are divisible
by p. Give each of the d kernel indices an extra half-weight 1/2 and the
remaining indices weight zero. Every entry has valuation at least the sum of
its new row weights, so every determinant term has valuation at least d.
No nonsingularity assumption on either complementary block is needed.

This is stronger than merely noticing a d-dimensional kernel modulo π, which
would only give d/2 by a generic rank argument. The gain comes from knowing
that all selected vectors belong to the **same weight-parity sector**.

For an original matrix B satisfying v(B_ij)≥w_i+w_j with half-integer weights,
set C_ij=p^{-w_i-w_j}B_ij. Rows with integral w_i form P and those with
nonintegral half-integer w_i form Q. Same-parity entries have integral valuation;
mixed-parity entries have valuation at least 1/2. The lemma upgrades
v(det B)≥2∑w_i to v(det B)≥2∑w_i+d whenever the stated kernel is established.

## 2. General annihilators from the two zero pole values

Let I_N⊂I_K be the near-pole intervals and put
R(u)=U_N(u)^r/U_K(u). After cancellation write

\[
 R(u)=\frac{U_N(u)^{r-1}}{\prod_{j\in I_K\setminus I_N}(u-j)}.
\]

Let Z={0,−1}∩(I_K\I_N), e=|Z| and

\[
 q_0(u)=\prod_{j\in(I_K\setminus I_N)\setminus Z}(u-j),
 \quad P=\deg q_0=|I_K|-|I_N|-e.
\]

Multiplication by q₀ kills every pole whose τ_{s,0} value might be nonzero:

\[
 q_0(u)R(u)=\frac{U_N(u)^{r-1}}{\prod_{j\in Z}(u-j)}.
\]

The remaining poles have zero functional values. Its polynomial quotient
after multiplication by u^{k+j} has degree at most
k+j+(r−1)|I_N|−e. Therefore

\[
 \tau_{s,0}(u^{k+j}q_0R)=0
 \quad\text{if}\quad k+j\le s-3-(r-1)|I_N|+e. \tag{1}
\]

For an L×(L+1) bordered block, testing column degrees 0,…,L, this supplies
the following number of independent explicit polynomial kernel vectors:

\[
 \boxed{d_{\rm explicit}=
 \max\{0,\,1+\min(L-1-P,\ s-3-(r-1)|I_N|+e-L)\}.} \tag{2}
\]

This is a lower bound on kernel dimension, not a general equality theorem.
An unrelated exact cancellation could in principle create more kernel vectors.

For I_N empty, e=2 and P=|I_K|−2, so (2) becomes
max(0,1+min(L−|I_K|+1,s−1−L)). At s=7, |I_K|=6,L=6 it equals one.
If I_N={0}, then e=1 and the degree threshold in (1) is s−r−1. For s=7,r=6
it is zero, too small for the relevant bordered blocks. Thus the actual
extra-root classes do not inherit the same useful annihilator.

## 3. Exact finite classification, including possible extra rows

Fix α=3/40, λ=37/40, s=7,r=6, H=23/20 and 3<x<10. For each class z∈(0,1/2),
form the actual intervals I_K,I_N and row sizes
L=⌊2Hx⌋−3|I_N|+ε with ε∈{0,1}. Including both ε values is conservative:
some sizes need not be selected by a particular global allocation.

Every floor pattern is constant between the following rational x breakpoints:

- integer and half-integer x;
- half-integer αx;
- integer (1−α)x or (1+α)x, where the two class thresholds exchange order;
- integer 2Hx, where the base row allocation changes.

Inside each x interval, split z at {x},1−{x},{αx},1−{αx}, keeping only points
inside (0,1/2). Rational midpoint sampling now enumerates **all** open-cell
patterns, rather than sampling a numerical grid. Boundary classes affect only
O(1) counts per prime and do not contribute to the leading prime-density term.

There are 45 x intervals and 110 distinct bordered matrices. The script
`python codex/computations/classify_inner_border_blocks.py` constructs their
moments exactly over Q. For 109 matrices, elimination modulo 1,000,003 gives
full row rank. Their denominators are units at this prime, so a nonzero modular
minor certifies a nonzero rational minor. The remaining matrix is checked by
exact rational rank and is precisely

\[
 I_K=\{-3,-2,-1,0,1,2\},\quad I_N=\varnothing,\quad L=6,
 \qquad\operatorname{rank}=5.
\]

Its one-dimensional kernel is exactly the explicit space from (2). Thus there
are no other rational bordered rank defects in this bounded region for these
parameters. This does not classify regions x≥10, other parameters, higher
p-adic lifts, or global cancellations not detected by a local bordered block.

## 4. Restoring the unit factors in the actual CRT block

The actual leading self matrix is a common unit times
diag(d_i) B diag(d_i), where d_i=(−2a)^i times any fixed unit normalization
of the rows. If q(u)=∑q_i u^i is the polynomial null vector, its actual
normalized row coefficients are q_i/d_i. All d_i are p-units because
0<a<p/2. These adjustments preserve the nullspace and primitivity.

For q=(u+3)(u+2)(u−1)(u−2), the coefficient at degree four is one, so the
adjusted degree-four coefficient is a unit. One may replace that row with
the kernel vector by a unimodular operation. Distinct selected classes use
disjoint rows, so these operations can be performed simultaneously. The
coefficients are independent of X. This verifies the algebraic hypotheses
of the parity-weight lemma once the global leading couplings vanish as
described in the earlier note.

## 5. Interpretation

The exact higher-odd cancellation exists, but its first instance is narrow.
The small candidate saving 1/25200 does not suggest that a large missing gain
is hidden throughout the tested interval. The exhaustive finite classification
instead gives a clear next boundary: larger x, another allocation geometry,
or higher-order p-adic information is required for a substantially stronger
effect. No such improvement is claimed here.
