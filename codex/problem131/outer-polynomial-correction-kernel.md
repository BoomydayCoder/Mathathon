# A kernel in the actual outer polynomial correction

Status: independently derived finite-field kernel and valuation-lifting proof;
adversarial review passed in `outer-correction-adversarial-review.md`. An explicit
Hermite-ideal proof is now in `outer-correction-Hermite-ideal.md`, with independent
review. This improves a local arithmetic estimate; the other prime ranges and
real determinant estimate remain separate dependencies.

**Attribution of the research lead.** The parent read, without editing, Claude's
published `SUMMARY` on `origin/research-claudefc6da97` and relayed its suggestion
of possible additional outer-prime savings around 2.47<K/p<3. This note's
derivation began from Fauzan's original local functional and our parity-kernel
method. Claude's proof drafts were not read. The saving derived here is
27/3200≈0.0084375, not the approximately 0.02 suggested in that summary.
No historical novelty is asserted.

## 1. Precise proposed finite-prime improvement

Use K∈40Z_{>0}, N=3K/40, h=K−N and the source's s=5,r=6 determinant Δ_K.
Let p satisfy the source's outer hypotheses (4.9), with 2p<K<3p, and suppose

\[
 d:=K+3N-3p+1>0.
\]

Let γ_p^{out} be the source's original bound (including min(rank,zero weights)).
The proposed strengthened estimate is

\[
 \boxed{v_p^G(\Delta_K)\ge\gamma_p^{out}+d.} \tag{1}
\]

The finite +1 is retained. Its limiting positive range begins at
K/p>120/49. The proof uses the structure of the actual polynomial correction,
not merely its rank, and its kernel is independent of the indeterminate X.

## 2. A finite-field formula for the polynomial correction

For a polynomial Q(t) with p-integral coefficients, von Staudt–Clausen gives

\[
 p\mu(Q)\equiv\frac1{24}\sum_{a\in\mathbb F_p}
       \left.\frac{d^3}{dv^3}\big[v^5Q(-v^2)\big]\right|_{v=a}
            \pmod p. \tag{2}
\]

Indeed pB_k≡∑_{a∈F_p}a^k mod p for every k≥0, including k=0 and k=1.
Apply this to τ(P)=L(P''')/24 and μ(Q)=τ(v⁵Q(−v²)). All denominators 24 are
units since p≥7. Formula (2) localizes the otherwise global Bernoulli correction
to third derivatives at finite-field nodes.

Write D=D_tail=D_K/D_N and W=D_N⁵. The source's CRT row polynomials are
g_{a,i}=P_a q_{a,i}, P_a=D/Q_a. Modulo p,

\[
 Q_a=(t+a^2)^{\ell-\delta},\qquad q_{a,i}=(t+a^2)^i,
\]

where δ=1 for a removed class 1≤a≤N and zero otherwise. The latter identity
also holds for the special rows containing the bad-pole factor E_a. In a
removed class, W has order five at −a² and ℓ−δ≤5. Thus, for every selected
row g_c defined below,

\[
 H_c:=Wg_c/D\quad\text{is a polynomial modulo }p. \tag{3}
\]

Polynomial division commutes with reduction because D is monic. Consequently
the quotient polynomial of Wg_cg_j/D reduces exactly to H_cg_j, with no
remainder. Formula (2) therefore computes the leading correction entry for
these rows. This avoids assigning finite-field values to unresolved poles.

## 3. Full weight chart and selected rows

In this band ordinary class sizes are ℓ=4,5,6. The source's weights are:

| Class | Row weights, in local order |
|---|---|
| unremoved, ℓ=4 | −4, −3, 0, 0 |
| unremoved, ℓ=5 | −9/2, −7/2, −5/2, 0, 0 |
| unremoved, ℓ=6 | −5, −4, −3, −2, 0, 0 |
| removed, ℓ=5 | −3/2, −1/2, 0, 0 |
| removed, ℓ=6 | −2, −1, 0, 0, 0 |
| zero class p,2p | −2, 0 |

No removed ℓ=4 class occurs in the stated positive-d range. Let Z be all
zero-weight rows. Let C⊂Z consist of the removed-class zero rows: local orders
2,3 for ℓ=5, and 2,3,4 for ℓ=6. Let T be the old weight −1/2 rows; these are
exactly local order one in removed ℓ=5 classes.

The source decomposes its matrix as A+p⁻¹L, with A satisfying the listed
weights and L integral. Its reduction L mod p is precisely (2) applied to the
polynomial quotient, after the CRT change of basis.

## 4. Which correction entries can survive?

Let c∈C. In a class other than the origins of g_c and g_j, both row factors
have full local order, so the reduced rational product H_cg_j has order at
least four. At the origin of g_c, the orders are also at least four whenever
j∈Z or j∈T. For example in the only potentially delicate case, a removed
ℓ=5 class with orders 2 and 1, the order is
6+2+1−5=4. In a removed ℓ=6 class, two zero rows of order two give
6+2+2−6=4. Thus these nodes contribute nothing to the third derivative.

At an **unremoved** origin of a zero row g_j, the full factor in g_c cancels
the local denominator, leaving exactly the order of g_j. Consequently only
the following target zero columns can contribute:

- both zero rows, of orders 2 and 3, in an unremoved ℓ=4 class;
- the first zero row, of order 3, in an unremoved ℓ=5 class;
- none in an unremoved ℓ=6 class.

The zero node contributes nothing because v⁵ times a polynomial has zero
third derivative at v=0. It follows that

\[
 L_{C,C}\equiv0,\qquad L_{C,T}\equiv0\pmod p,
\]

and L_{C,Z} is supported on at most the explicitly listed target columns.
This includes the old −1/2 sector, which cannot be ignored in the valuation
lift. The vanishing uses the actual third-derivative functional; an arbitrary
matrix with the same row weights need not satisfy it.

## 5. Exact class counts

Put u=(K+N−3p+1)₊. There are u removed ℓ=6 classes and N−u removed ℓ=5
classes, hence

\[
 |C|=2N+u.
\]

The number of potentially nonzero target columns from Section 4 is

\[
 b=(3p-K-N-1)_+.
\]

For K>(5p)/2, this counts unremoved ℓ=5 classes. For K<(5p)/2, it equals
the number of unremoved ℓ=5 classes plus twice the number of ℓ=4 classes:
(K−2p−N)+2((p−1)/2−(K−2p))=3p−K−N−1. These expressions also show directly
that removed ℓ=4 classes do not arise here. Therefore

\[
 \dim\ker\big(L_{C,Z}^{\mathsf T}\big)
       \ge |C|-b=K+3N-3p+1=d. \tag{4}
\]

This notation means row combinations supported on C annihilating all columns
in Z. They automatically annihilate T as well. The kernel is over F_p and
independent of X; lift a basis to primitive p-integral vectors.

## 6. Recover a full p-power per kernel vector

In this range the source's rank bound r_p=K+4N−2p+2 exceeds the number
z=|Z|=p+u: if u=0 the difference is d+N+1>0; if u>0 it is 3N+1>0.
Thus its original exponent is γ_p^{out}=2∑w_i−z.

Implement this plain correction loss entrywise by replacing every old zero
weight with −1/2, leaving negative weights unchanged. Call the new weights
w_i'. Every entry of A+p⁻¹L has valuation at least w_i'+w_j', since all new
weights are at most −1/2 and L is integral. Its baseline determinant bound
is exactly γ_p^{out}.

After scaling out the new weights, all C rows belong to the same half-integer
parity sector. Against any old zero row, the leading scaled entry is L_{C,Z}
mod p; the A part is divisible by p. Against an old −1/2 row it is L_{C,T}
mod p, again with the A part divisible by p. Against an old negative weight
−k−1/2 with k≥1, both contributions are automatically divisible by p.
Against integer-weight columns the scaled entries have valuation at least
1/2 automatically.

The d kernel combinations from (4) therefore satisfy the standalone parity
lemma in [`inner-annihilator-classification.md`](inner-annihilator-classification.md):
they gain half a row weight each, hence **one full power of p each** in the
determinant valuation. This proves the proposed bound (1), subject to the
independent audit of the source local decomposition and the finite-field
identifications used above. No assumption of invertibility of a complementary
block is needed, and the argument applies coefficientwise in X.

## 7. Reproducible finite checks and asymptotic consequence

`python codex/computations/check_outer_correction_kernel.py` constructs the
actual reduced CRT polynomials, computes the polynomial correction by exact
finite-field arithmetic, checks every forbidden target column including T,
and computes its rank. The seven independently selected examples give:

| K | p | selected rows | allowed target columns | actual correction rank | kernel dimension d |
|---:|---:|---:|---:|---:|---:|
| 80 | 29 | 12 | 0 | 0 | 12 |
| 80 | 31 | 12 | 6 | 6 | 6 |
| 120 | 43 | 19 | 0 | 0 | 19 |
| 120 | 47 | 18 | 11 | 11 | 7 |
| 160 | 61 | 24 | 10 | 10 | 14 |
| 200 | 73 | 30 | 3 | 3 | 27 |
| 240 | 97 | 36 | 32 | 32 | 4 |

These are finite correction-matrix checks, not computations of the whole
determinant or a replacement for the all-prime argument.

If the local proposition is accepted, the same prime-number-theorem conversion
used by the source saves

\[
 \int_{120/49}^{3}\frac{(49/40)x-3}{x^3}\,dx
 =\frac{(49/40-1)^2}{6}=\boxed{\frac{27}{3200}}.
\]

Equivalently integrate [(1+3α)−3y]₊ over 1/3<y<(1+3α)/3.
The O(1) finite-prime rounding contributes o(K²). No saving outside this band
or beyond this explicitly proved kernel is included. The original source
integrality result, other local prime ranges, and real determinant estimate
remain separate dependencies.
