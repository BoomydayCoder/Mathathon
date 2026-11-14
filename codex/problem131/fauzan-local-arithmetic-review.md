# Independent audit of Fauzan's local arithmetic, with a simpler rank bound

**Status.** The proofs of Propositions 4.1, 4.3 and 5.1 have been independently reconstructed below from the primary source `../sources/fauzan-zeta5.pdf`, §§3–5. No defect was found in those three propositions. The local functional argument needed for the large primes was also checked. A subsequent reconstruction of Lemma 3.3, appended in §8, checks the small-prime estimate as well. Standard Bernoulli and integer-valued-polynomial identities are stated explicitly below. This is not a certification of the entire preprint or a historical-novelty claim.

There is a useful simplification: replace the refined correction `min(rank, number of zero weights)` by the rank alone. The extra asymptotic arithmetic cost is exactly **9/640**. Together with the separate table-free inner-integral estimate, the arithmetic bound becomes **2203747/1600000 = 1.377341875**, strictly below the team's independently certified analytic target **69/50 = 1.38**. This removes the zero-weight enumeration and the refined rank lemma from the proposed proof. The analytic target and table-free estimate are dependencies on the separate notes, not results established here.

## 1. Objects and motivation

Write

\[
 K=40n,\quad N=3n,\quad h=37n,\quad
 \alpha=3/40,\quad\lambda=37/40,\quad H=1+2\alpha=23/20.
\]

Let \(D_A(t)=\prod_{j=1}^A(t+j^2)\). The matrix is

\[
 G_K(X)_{ij}=\mu_X\left(\frac{D_N(t)^6t^{i+j}}{D_K(t)}\right),
 \qquad 0\le i,j<h,
 \quad\Delta_K=\det G_K.
\]

Its moments are

\[
 \mu_X(t^e)=\frac{(-1)^e B_{2e+2}(2e+3)(2e+4)(2e+5)}{24},
\]
\[
 \mu_X((t+j^2)^{-1})
 =j^4(X-H_j^{(5)})-\frac14+\frac1{2j}.
\]

The scalar normalization is

\[
 S_K=\frac{(K!)^{2h}4^{h-1}}
 {(N!)^{12h}\prod_{i=1}^{h-1}((2i)!)^2},\qquad F_K=S_K\Delta_K.
\]

Each prime presents a different collision pattern among the nodes \(-j^2\). The local bases exploit these collisions. A row assigned to one residue class vanishes strongly at all other classes, making its contribution uniformly divisible. The determinant then adds the assigned row divisibilities. This is a local interpolation construction, rather than a direct estimate of the original monomial matrix.

The inner range has many colliding poles and uses balanced allocation. The outer range has only a few poles per class and uses direct cancellation, with the nonintegral polynomial part recorded as a low-rank matrix. Neither mechanism proves smallness over the real numbers; it independently controls the arithmetic cost of making a polynomial integral.

## 2. Local functional input checked for the large-prime arguments

The pullback identity is

\[
 \mu_X(R(t))=\tau_X[x^5R(-x^2)],
\]

where \(\tau(P)=L(P''')/24\), \(L(x^k)=B_k\), and

\[
 \tau_X((x-r)^{-1})=H_{d(r)}^{(5)}-X,
 \quad d(r)=r\ (r\ge0),\quad d(r)=-r-1\ (r<0).
\]

For \(p\ge7\), von Staudt–Clausen gives valuation at least \(-1\) for the polynomial moments of \(\tau\), and integrality through degree \(p+1\). Division by a monic polynomial with integer roots separated by units preserves integral coefficients; its partial-fraction denominators are units. Thus the source's Lemma 3.1 is valid for a leading numerator of degree at most \(p+1\), simple integer poles separated by units, and harmonic indices below \(p\). Higher terms in a convergent expansion \(\sum_{j\ge0}p^jU_j/T\) absorb the possible one-factor \(p\) loss of \(\tau\). The operator bound is independent of the polynomial degrees, so varying degrees in the expansion cause no convergence problem.

The distribution identity is

\[
 \tau_X(g)=p^{-4}\sum_{c=0}^{p-1}
 \tau_Y^{\rm ext}[g(c+px)],\qquad
 Y=p^5X+C_p,\quad C_p\in p^5\mathbf Z_p.
\]

On polynomials, the factor \(p^4\) comes from three derivatives and Bernoulli multiplication. On the simple pole \(1/x\), the identity is exactly the definition of \(C_p\); reflection and translation give the other integer poles. For the asserted divisibility of \(C_p\), its zeroth expansion term gains a factor \(p\) because \(\sum_{a\in\mathbf F_p^\times}a^{-4}=0\) for \(p\ge7\); higher terms have their required powers of \(p\), even allowing the one-factor Bernoulli denominator. In particular, substituting \(Y=p^5X+C_p\) preserves all coefficientwise integrality statements. The argument controls the **Gauss valuation**, not merely values at integer \(X\).

## 3. Proposition 4.1: balanced allocation for the inner primes

Assume

\[
 M\ge40,\quad K\ge200M^2,\quad K/M<p\le K/3.
\]

Then \(p>200M\) and \(p^2>5K\). Put

\[
 \ell_A(a)=\#\{1\le j\le A:j\equiv\pm a\pmod p\},\quad
 m_A=\lfloor A/p\rfloor.
\]

If row polynomials vanish to orders \(\nu_i(a),\nu_j(a)\) at \(-a^2\), the source valuations after the distribution formula are

\[
 \nu_i(a)+\nu_j(a)+6\ell_N(a)-\ell_K(a)-4\quad(a\ne0),
\]
\[
 2\nu_i(0)+2\nu_j(0)+12m_N-2m_K+1\quad(a=0).
\]

**Why these exponents are correct.** At a nonzero class, \(-x^2+a^2\) has a simple zero, each near denominator factor contributes one \(p\), and the functional distribution loses four. At zero, the row zeros and denominator zeros have doubled order under \(t=-x^2\), while the factor \(x^5\) contributes five powers; five minus four explains the final \(+1\). Confusing the zero class with an ordinary one would invalidate the construction.

Set \(L_0=4M+10\), \(m=(p-1)/2\), \(b_a=3\ell_N(a)\), and choose integers \(T,E\) by

\[
 mT+E=h-L_0+3(N-m_N),\qquad 0\le E<m.
\]

Give the first \(E\) classes, ordered by decreasing \(\ell_K\), an extra \(\epsilon_a=1\); set

\[
 L_a=T-b_a+\epsilon_a,\qquad Z_a=L_a+b_a=T+\epsilon_a.
\]

All dimensions are nonnegative and sum to \(h\). To check the constants, for \(x=K/p\),

\[
 T=\left\lfloor\frac{2HK-2L_0-6m_N}{p-1}\right\rfloor,
 \quad 2Hx-21/20<T<2Hx,
 \quad b_a\le6\alpha x+3.
\]

The first error follows from \((2L_0+6m_N)/(p-1)<1/20\); the upper bound uses \(2Hx<2L_0+6m_N\). Consequently \(T-b_a>2\lambda x-81/20\ge3/2\).

The row polynomials

\[
 E_{a,i}(t)=\prod_{c\ne a}(t+c^2)^{L_c}(t+a^2)^i,
 \qquad 0\le i<L_a,
\]

form a \(\mathbf Z_p\)-unimodular basis: reduction into the pairwise coprime powers \((t+a^2)^{L_a}\) gives triangular local blocks with unit diagonal. This includes the block at zero.

The assigned weights are

\[
 w_{a,i}=i+b_a-(\ell_K(a)+4)/2\quad(a>0),
\]
\[
 w_{0,i}=\min\left(2i+6m_N-m_K+\tfrac12,
 \min_{c>0}\{Z_c-(\ell_K(c)+4)/2\}\right).
\]

For a row belonging to another ordinary block, the comparison reduces to

\[
 Z_c-Z_a+1+(\ell_K(a)-\ell_K(c))/2\ge0.
\]

The counts \(\ell_K\) take two consecutive values. The only potentially negative case has an extra row at \(a\) but not at \(c\); decreasing-order allocation then ensures \(\ell_K(a)\ge\ell_K(c)\). Thus the ordering of extras is a necessary part of this proof, not decoration.

At the zero source, every ordinary row has half-weight at least

\[
 2L_0+6m_N-m_K+\tfrac12\ge7M+\tfrac{41}2,
\]

larger than its assigned ordinary weight. The definition of \(w_{0,i}\) handles all comparisons involving zero-block rows.

The local analytic hypotheses also hold, not merely the formal valuations. Near poles have absolute size at most \(M+1\), distinct differences are units since \(2M+2<p\), and harmonic indices are below \(p\). Leading numerator degrees are bounded by

\[
 2(T+1)<(23/5)M+2<p+1
\]

at ordinary sources, and

\[
 5+4L_0+12m_N\le(169/10)M+45<p+1
\]

at zero. The far factors expand in integral powers of \(pz\). Lemma 3.1 therefore applies to the full local summands.

Every transformed matrix entry has Gauss valuation at least the sum of its row weights. Every determinant term has valuation at least \(2\sum w_{a,i}\), proving Proposition 4.1.

## 4. Proposition 4.3: outer primes and actual rank

Assume the source's conditions

\[
 p\ge7,\quad p\le K<3p,\quad p^2>2K,
 \quad2N<p,\quad5N\le2p-2.
\]

### Polynomial correction

Polynomial moments have valuation at least \(-1\), and are integral for \(e<2p-3\). Indeed, potential Bernoulli denominators have \(p-1\mid2e+2\). The first, second and third positive multiples are canceled, respectively, by \(2e+3\), \(2e+4\), \(2e+5\). The fourth multiple first appears at \(e=2p-3\).

Subtracting \(c_e/p\) from each nonintegral polynomial moment, with \(c_e\in\mathbf Z_p\), produces

\[
 G_K=A+p^{-1}L,\qquad
 \operatorname{rank}_{\mathbf Q_p}L\le r_p=\max(0,K+4N-2p+2).
\]

This is an actual rank statement over \(\mathbf Q_p\), which is sufficient for the direct multilinearity proof below. Actual rank is not necessary: a rank bound modulo \(p\) also suffices after lifting a unit pivot block and absorbing an integral remainder into \(A\); see [refined-rank-two-bounds.md](refined-rank-two-bounds.md), Section 3. The quotient of \(D_N^5t^{i+j}\) by \(D_K/D_N\) has degree at most \(i+j+6N-K\), so

\[
 L_{ij}=0\quad\text{if}\quad i+j<K-6N+2p-3.
\]

The first \(h-r_p\) rows and columns vanish identically, so all larger minors vanish exactly. This directly justifies the actual-rank proof; the modular-rank version first replaces \(L\) by a suitable lift of small actual rank. Unimodular basis changes preserve both integrality of \(L\) and its rank.

### Residues and the choice of rows

For a nonzero square class let \(\ell=\ell_K(a)\), \(\delta=1_{a\le N}\), and let \(Q_a\) be its remaining denominator factors. Write \(P_a=(D_K/D_N)/Q_a\). There are \(\ell-\delta\) remaining nodes; exactly \(\ell-2\) original nodes have index \(j>p\). Let \(E_a\) be their product.

Use rows \(P_aq_{a,i}\), where

\[
 q_{a,i}=(t+a^2)^i\quad(i<\ell-2),\qquad
 q_{a,i}=E_a(t)(t+a^2)^{i-(\ell-2)}\quad(i\ge\ell-2).
\]

They have successive monic degrees in each local block, so the resulting basis is unimodular. Cross-class entries are polynomials and thus integral for the corrected functional.

If either index is at least \(\ell-2\), all poles with \(j>p\) cancel. The two possible remaining moments, at \(a\) and \(p-a\), are integral and congruent modulo \(p\). The required identity is

\[
 H_{p-a}^{(5)}\equiv H_{a-1}^{(5)}\pmod p,
\]

which follows from \(\sum_{u\in\mathbf F_p^\times}u^{-5}=0\) for \(p\ge7\). For the constant terms their difference reduces to \(-a^4/a^5+1/a=0\); the \(X\)-coefficients are congruent as well. The node difference has valuation exactly one. Their divided difference is integral, and integral polynomial numerators preserve this congruence.

For the uncanceled low rows the valuation is at least

\[
 i+j+6\delta-\ell-4.
\]

The accounting is: row zeros contribute \(i+j\), the numerator \(D_N^5\) contributes \(5\delta\), the node derivative loses \(\ell-\delta-1\), and the harmonic value loses at most five. The node differences have valuation exactly one because \(p^2>2K\) and nodes in a fixed square class cannot be simultaneously congruent with both signs. Thus valid nonpositive weights are

\[
 w_{a,i}=\min(0,i+3\delta-(\ell+4)/2)\quad(i<\ell-2),
 \qquad w_{a,i}=0\quad(i\ge\ell-2).
\]

The zero class with one pole has weight \(-1/2\). With two poles \(p,2p\), the node difference has valuation two; in the rows \(1,t+p^2\) the entry bounds are \(-3,-1,0\). Weights \(-2,0\) are valid and nonpositive.

### The refined lemma is valid, but dispensable

**Subsequent simplification.** The original refined bound has a shorter
proof in [refined-rank-two-bounds.md](refined-rank-two-bounds.md):
multilinearity gives loss at most \(r\); independently, replacing every
zero weight by \(-1/2\) makes the correction integral after scaling and
gives loss at most \(z\). Taking the better bound gives \(\min(r,z)\)
without sorting minors. The older argument and plain-rank alternative
below remain valid; dropping the refinement is no longer needed merely
to shorten the lemma's proof.

For complementary row/column sets \(I,J\), a term using a \(k\)-minor of \(L\) has valuation at least

\[
 2\sum_iw_i-k-\sum_{i\in I}w_i-\sum_{j\in J}w_j.
\]

Terms with \(k>r\) vanish. If there are \(z\) zero weights, sorting the weights shows the maximal loss is \(\min(r,z)\): adding a zero weight increases the loss by one, while adding a weight at most \(-1/2\) gives a nonpositive increment. This proves the original Lemma 4.2, including coefficientwise Gauss valuation.

For a simpler proof, just use \(w_i\le0\) in the displayed inequality:

\[
 \boxed{\ v^G\det(A+p^{-1}L)\ge2\sum_iw_i-r.\ }
\]

This needs no count of zero weights and no optimization of the minors.

For completeness, the original outer aggregation was independently tested for all admissible primes and \(K\in40\mathbf Z\), \(40\le K\le1600\): **3335 exact cases**. The direct class-weight sums and zero-weight counts agree with (4.14). This finite check supplements the proof; it does not replace it.

For \(p>K\), each square class has at most two nodes, all pole moments are integral, the same divided-difference argument applies, and the polynomial correction vanishes because \(K+4N+1<2p\). Hence \(\Delta_K\) is integral at such primes.

## 5. Proposition 5.1: coefficientwise normalization

The small-prime input, reconstructed in §8 from Lemma 3.3, is

\[
 v_p^G(F_K)\ge-6h\lfloor\log_p(5K)\rfloor-hv_p(24).
\]

For larger primes use \(v_p(S_K)+\gamma_p\), with either the original or relaxed outer exponent. Define \(L_p\) by these piecewise lower bounds and

\[
 m_{K,M}=\prod_{p\le2h}p^{-L_p}.
\]

**The signs must be retained:** a positive lower valuation permits division by a power of \(p\). Replacing \(L_p\) by its negative part would change the asymptotic arithmetic budget.

At \(p\le2h\), every coefficient of \(mF_K\) has nonnegative valuation. At \(p>2h>K\), all factorials occurring in \(S_K\) are units and \(\Delta_K\) is integral by the preceding argument. Thus all rational coefficients are integral at every prime, so \(mF_K\in\mathbf Z[X]\). This proves Proposition 5.1 from its explicit local inputs. No inference from merely integer values of a polynomial is being made.

## 6. Exact budget for dropping the refined rank correction

Put \(y=p/K\). The original saving is

\[
 d(y)=\left(1+4\alpha-3y-(1+\alpha-3y)_+\right)_+
 \mathbf1_{(1/3,1/2)}(y).
\]

For \(\alpha=3/40\), it equals \(9/40\) on \((1/3,43/120)\), then \(13/10-3y\) on \((43/120,13/30)\), and zero elsewhere. Therefore

\[
 \int d(y)\,dy
 =\frac9{40}\left(\frac{43}{120}-\frac13\right)
 +\frac12\frac9{40}\left(\frac{13}{30}-\frac{43}{120}\right)
 =\frac9{640}.
\]

The relaxed outer integral is consequently

\[
 I_{\rm out}^{\rm relaxed}=\frac{127751}{96000}+\frac9{640}
 =\frac{129101}{96000}.
\]

Using the original source's exact finite inner integral and tail bounds, this already gives a bound of approximately \(1.3715800064<1.38\) at \(M=160\). A cleaner combined option uses the independently derived table-free inner bound from `fauzan-arithmetic-simplification.md`:

\[
 \frac{127751}{96000}+\frac{13}{144}
 +\frac{37}{40}\left(\frac53+\frac35-
 \frac9{10}\frac{33559}{12960}\right)
 +\frac{14(37/40)}M+\frac9{640}.
\]

At \(M=10000\), this is exactly

\[
 \boxed{\frac{2203747}{1600000}=1.377341875
 <\frac{69}{50},\qquad
 \text{gap}=\frac{4253}{1600000}.\ }
\]

As \(K=40n\), the combined exponent has margin \((4253/1000)n^2\). If the analytic and arithmetic limsup estimates are established as stated, one can absorb lower-order errors and eventually obtain \(0<|mF_K(\zeta(5))|<e^{-4n^2}\). This remains conditional on the complete analytic and small-prime inputs, as well as the polynomial nonvanishing argument; it is not an unconditional conclusion of this local review.

**Source-reading correction.** PDF text extraction displaces the scalable parentheses in §5.2. The actual PDF page 15 correctly reads
\(-\gamma_p^{\rm out}=K(R_0(p/K)-d(p/K))+O(1)\), with \(K\) multiplying both terms. The rendered page was inspected; there is no missing-factor typo in the source.

## 7. What remains genuinely nontrivial

1. The small-prime bound (3.12), including the integer-valued-polynomial argument in Lemma 3.3, remains necessary; §8 below supplies its independent reconstruction. It uses the Bernoulli generating identity and properties of the binomial polynomial basis.
2. The Bernoulli distribution identity and local functional extension are structurally indispensable; the checks above explain their constants and applicability, but a self-contained final proof should reproduce their complete definitions.
3. The balanced CRT allocation for the inner primes remains necessary with this budget. Removing the refined outer rank saving does not remove that construction.
4. The transition from prime sums to integrals uses the prime number theorem. No unproved assumption about exceptional-prime density is present, but the proof is not elementary in the sense of avoiding analytic prime asymptotics.
5. The table-free integral and analytic determinant estimate are independently documented elsewhere. Their combination must be reviewed together with positivity/nonvanishing and the final rationality contradiction before claiming a full theorem.

The useful conceptual separation is: **local interpolation supplies divisibility; a low-rank correction measures the failure of polynomial moments to be integral; global prime averaging prices those local corrections.** The stronger real estimate permits the rank price to be paid without optimizing zero-weight directions.

Reproduce the finite and rational checks with `python codex/computations/check_fauzan_local_arithmetic_review.py`.


## 8. Follow-up: the small-prime estimate also checks out

This section removes the earlier local-review dependency on Lemma 3.3 by giving the argument in detail. It is valid for **every prime**, including 2 and 3; their extra denominator is recorded by \(v_p(24)\).

### Integer-valued polynomials and the Bernoulli functional

If \(A(\mathbf Z_p)\subset\mathbf Z_p\) and \(\deg A\le d\), its binomial expansion
\(A(x)=\sum_{j=0}^d a_j\binom{x}{j}\) has integral coefficients: \(a_j=\Delta^jA(0)\).

For \(u\in\mathbf Z_p\),
\[
 \binom{u}{j}=\frac uj\binom{u-1}{j-1},\qquad
 v_p\binom uj\ge v_p(u)-\lfloor\log_p j\rfloor.
\]
Vandermonde's identity therefore gives, with \(L_d=\lfloor\log_p\max(1,d)\rfloor\),
\[
 v_p(A(x)-A(y))\ge v_p(x-y)-L_d.
\]

The coefficient of \(z^r\) in \(\binom zj\) is a sum of products of \(r\) reciprocals of integers at most \(j\), because
\[
 \binom zj=\frac{(-1)^{j-1}z}{j}
 \prod_{m=1}^{j-1}\left(1-\frac zm\right).
\]
Another application of Vandermonde proves that the \(r\)-th derivative of an integer-valued polynomial has values of valuation at least \(-rL_d\); multiplication by \(r!\) only improves this estimate.

Finally,
\[
 L((1+z)^x)=\frac{\log(1+z)}z,
 \qquad L\binom xj=\frac{(-1)^j}{j+1}.
\]
Thus, if \(P(\mathbf Z_p)\subset p^{-A}\mathbf Z_p\), \(\deg P\le d\), then
\[
 v_p\tau(P)\ge-A-4\lfloor\log_p(d+1)\rfloor-v_p(24).
\]
Three logarithmic losses come from differentiating three times; the fourth comes from applying \(L\) to the binomial basis. There is no appeal here to a bound on the monomial coefficients of an integer-valued polynomial.

### Root balls: proving Lemma 3.3

Let
\[
 g(x)=\frac{(K!)^2A(x)}{\prod_{-K\le r\le K,\ r\ne0}(x-r)},
 \qquad A(\mathbf Z_p)\subset\mathbf Z_p,
\]
and write \(g=P+\sum_{r\ne0}c_r/(x-r)\). Put
\[
 L_0=\lfloor\log_p(2K)\rfloor,
 \qquad L=\lfloor\log_p\max(2K,d+1)\rfloor.
\]
The residue formula is
\[
 c_r=\pm\frac{(K!)^2rA(r)}{(K+r)!(K-r)!},\qquad v_p(c_r)\ge-L_0.
\]
Indeed, at each prime-power level the loss between the two factorial floor sums is at most one, and there are at most \(L_0\) relevant levels.

Outside the balls \(r+p^{L_0+1}\mathbf Z_p\), each distance to a root has valuation at most \(L_0\). At the level \(p^a\), the number of roots in any residue class is at most
\(2\lfloor K/p^a\rfloor+2\): count in the full interval \([-K,K]\), then remove zero if necessary. The numerator contributes \(2\lfloor K/p^a\rfloor\). Therefore
\[
 v_p(g(x))\ge-2L_0.
\]
Each pole term has the same lower bound, so \(v_p(P(x))\ge-2L_0\) there.

The root balls are disjoint, since distinct integer roots have difference of absolute size at most \(2K\). In the ball around \(r\), define
\[
 F_r(x)=\frac{(K!)^2}{\prod_{s\ne r,\ s\ne0}(x-s)}.
\]
The residue estimate without \(A(r)\) gives \(v_p(F_r(r))\ge-L_0\). If \(\delta=x-r\), then
\[
 \frac{F_r(x)}{F_r(r)}=
 \prod_{s\ne r,\ s\ne0}
 \left(1+\frac\delta{r-s}\right)^{-1}.
\]
Every factor is a unit congruent to 1 modulo a quantity of valuation at least \(v_p(\delta)-L_0\). Hence
\[
 v_p(F_r(x))\ge-L_0,\qquad
 v_p(F_r(x)-F_r(r))\ge v_p(\delta)-2L_0.
\]
Subtract the principal part at \(r\):
\[
 \frac{A(x)F_r(x)-A(r)F_r(r)}{x-r}.
\]
The integer-valued-polynomial Lipschitz estimate bounds its valuation below by
\(-L_0-\max(L_0,L_d)\). Every other pole term has valuation at least \(-2L_0\), because its denominator distance equals that from its center to \(r\). Therefore
\[
 P(\mathbf Z_p)\subset p^{-L_0-\max(L_0,L_d)}\mathbf Z_p
 \subset p^{-2L}\mathbf Z_p.
\]
Continuity includes the centers of the balls. Since \(\deg P\le d\), the preceding Bernoulli-functional estimate gives
\[
 v_p\tau(P)\ge-6L-v_p(24).
\]
For each pole value, \(d(r)\le K\), so
\[
 v_pH_{d(r)}^{(5)}\ge-5\lfloor\log_pK\rfloor.
\]
Together with \(v_p(c_r)\ge-L_0\), this gives the same \(-6L\) bound for each coefficient of every pole contribution. It follows that
\[
 \boxed{v_p^G\tau_X(g)\ge-6\lfloor\log_p\max(2K,d+1)\rfloor-v_p(24).}
\]
This is precisely Lemma 3.3.

### Applying it to the determinant

The source uses
\[
 q_0(t)=1,\quad q_i(t)=\frac{(-1)^i2tD_{i-1}(t)}{(2i)!},\quad
 f_i(t)=\left(\frac{D_N(t)}{(N!)^2}\right)^3q_i(t).
\]
The identities
\[
 q_i(-x^2)=\binom{x+i}{2i}+\binom{x+i-1}{2i},
 \qquad
 \frac{D_N(-x^2)}{(N!)^2}=\binom{N-x}{N}\binom{N+x}{N}
\]
show that each \(f_i(-x^2)\) is integer-valued on \(\mathbf Z_p\).
The matrix of \((K!)^2\mu_X(f_if_j/D_K)\) has determinant \(F_K\): the squared triangular basis-change determinant contributes exactly the factors \(4^{h-1}\), \((N!)^{-12h}\), and \(\prod((2i)!)^{-2}\) in \(S_K\).

Its pullback numerator is \(x^5f_i(-x^2)f_j(-x^2)\), an integer-valued polynomial of degree at most
\[
 d=12N+4h+1=\frac{23}{5}K+1,
 \qquad\max(2K,d+1)\le5K.
\]
Lemma 3.3 bounds each matrix entry below by
\(-6\lfloor\log_p(5K)\rfloor-v_p(24)\). Summing along determinant terms proves the claimed bound
\[
 v_p^G(F_K)\ge-6h\lfloor\log_p(5K)\rfloor-hv_p(24).
\]

**Updated audit boundary.** The local normalization, including the small-prime estimate, is now independently reconstructed. The separate real determinant bound, positivity/nonvanishing, and full assembly of the rationality contradiction still require their own audits; nothing in this update changes that distinction.
