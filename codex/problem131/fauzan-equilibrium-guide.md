# Simplifying the analytic and arithmetic estimates in Fauzan’s construction

This note explains a concrete simplification of the analytic part of Fauzan’s 2026 \(\zeta(5)\) construction. An explicitly constructed equilibrium measure replaces the preprint’s mixture of sixteen arcsine measures and its many potential comparisons. A one-sign-change argument proves the global inequality; rigorous interval computations certify the two endpoint equations and the remaining energy bound. A separate elementary estimate removes the source’s arithmetic integration tables, conditional on its local valuation envelope. The analytic construction also extends to every fixed odd index, with a clearly identified arithmetic gap.

The construction and certificates were developed and independently checked within this session. Historical novelty has not been established. This is **not an independent audit of the complete Fauzan preprint**. The moment and local arithmetic arguments now have separate independent reconstructions, including the small-prime estimate. The combined final proof is still undergoing assembly review; this guide does not promote the entire preprint to an independently audited theorem. The source is [Fauzan’s preprint, version 1](https://doi.org/10.5281/zenodo.22826419), Sections 6–7 and Appendix A; detailed derivations are in [fauzan-analytic-simplification.md](fauzan-analytic-simplification.md).

## 1. Start from the quantity that actually needs to be bounded

Keep the source parameters

\[
\alpha=\frac3{40},\qquad \lambda=\frac{37}{40},
\]

and the external field

\[
V(x)=2\pi\sqrt{x}
 +\int_0^1\log(x+u^2)\,du
 -6\int_0^\alpha\log(x+u^2)\,du,
\qquad x\ge0.
\tag{1}
\]

For a positive measure \(\rho\) of mass \(\lambda\), write

\[
U_\rho(x)=\int\log|x-t|\,d\rho(t),\qquad
I(\rho)=\iint\log|x-t|\,d\rho(x)d\rho(t).
\]

The relevant analytic bound is

\[
\mathcal J(\rho)=\lambda M-I(\rho)+C_*,
\qquad
M=\sup_{x\ge0}\{2U_\rho(x)-V(x)\},
\tag{2}
\]

where

\[
C_*=-2\lambda+12\alpha\lambda(1-\log\alpha)
     +3\lambda^2-2\lambda^2\log(2\lambda).
\]

Why this functional? If \(\sigma\) is another measure of mass \(\lambda\), the zero-mass logarithmic-energy inequality gives

\[
I(\sigma)-\int V\,d\sigma
\le\int(2U_\rho-V)\,d\sigma-I(\rho)
\le\lambda M-I(\rho).
\tag{3}
\]

Indeed, apply \(I(\delta)\le0\) to \(\delta=\sigma-\rho\). For finite-energy compactly supported measures, the inequality follows by regularizing the Fourier representation
\(I(\delta)=-\int_0^\infty|\widehat\delta(t)|^2\,dt/t\); the zero total mass removes the otherwise divergent common constant. The usual limiting argument covers the measures used here.

An equilibrium measure makes \(2U_\rho-V\) constant on its support and no larger elsewhere. It therefore gives the best possible comparison of this form. Instead of optimizing a large table of trial weights, we can try to construct that measure directly.

## 2. The field hides a particularly simple signed charge distribution

Differentiating (1) gives

\[
V'(x)=\int_0^\infty\frac{q(u)}{x+u^2}\,du,
\quad
q(u)=
\begin{cases}
-3,&0<u<\alpha,\\
3,&\alpha<u<1,\\
2,&u>1.
\end{cases}
\tag{4}
\]

The identity \(2\int_0^\infty(x+u^2)^{-1}\,du=\pi/\sqrt x\) accounts for the square-root term. Each denominator \(x+u^2\) is the derivative of a logarithmic potential from a charge at \(-u^2\).

The crucial feature is the **single change from negative to positive charge at \(u=\alpha\)**. The later change from \(3\) to \(2\) preserves the sign. This will control the potential on the whole positive axis; no grid of sample points will be needed to establish its global shape.

The sign pattern also suggests why the support begins at a positive left endpoint. Near zero the negative portion dominates \(V'\); far away the positive square-root contribution dominates. This is motivation for the construction, not an assumption replacing its proof.

## 3. Sweep each charge onto an interval

Choose provisional endpoints \(0<a<b\), and let \(\mu\) be the unit arcsine measure on \([a,b]\):

\[
d\mu(t)=\frac{dt}{\pi\sqrt{(t-a)(b-t)}}.
\]

For each \(u\ge0\), set

\[
C_u=\sqrt{(u^2+a)(u^2+b)},\qquad
 dB_u(t)=\frac{C_u}{t+u^2}\,d\mu(t).
\tag{5}
\]

The elementary arcsine transform
\(\int(t+u^2)^{-1}d\mu(t)=C_u^{-1}\) shows that \(B_u\) has mass one. It is the measure obtained by sweeping the charge at \(-u^2\) onto the interval, often called *balayage*. Formula (5), rather than the terminology, is what the proof uses.

Define

\[
\rho=\lambda\mu+rac12\int_0^\infty q(u)(B_u-\mu)\,du.
\tag{6}
\]

Each difference \(B_u-\mu\) has mass zero, so \(\rho\) has mass \(\lambda\). The integral converges: uniformly on the interval,
\(C_u/(t+u^2)-1=O(u^{-2})\) at infinity, and \(a>0\) controls the origin.

At this stage \(\rho\) could be signed. Its density is

\[
d\rho(x)=\frac{F(x)}{\pi\sqrt{(x-a)(b-x)}}\,dx,
\quad
F(x)=\lambda+rac12\int_0^\infty q(u)
 \left(\frac{C_u}{x+u^2}-1\right)du.
\tag{7}
\]

The remaining choice of endpoints will make this density positive.

## 4. Two endpoint equations turn the signed construction positive

Require

\[
I_0(a,b):=\int_0^\infty\frac{q(u)}{C_u}\,du=0,
\tag{8}
\]

and

\[
H(a,b):=\lambda+rac12\int_0^\infty q(u)
 \left(\sqrt{\frac{u^2+a}{u^2+b}}-1\right)du=0.
\tag{9}
\]

These are precisely the conditions \(F(a)=F(b)=0\): we have \(H=F(b)\) and
\(F(a)-F(b)=(b-a)I_0/2\).

For \(x>0\), define

\[
K(x)=\int_0^\infty\frac{q(u)}{(x+u^2)C_u}\,du.
\]

Because \(I_0=0\), we may subtract the constant \(1/(x+\alpha^2)\) inside the integral:

\[
K(x)=\int_0^\infty\frac{q(u)}{C_u}
 \left(\frac1{x+u^2}-\frac1{x+\alpha^2}\right)du<0.
\tag{10}
\]

Before \(\alpha\), the charge is negative and the bracket positive. After \(\alpha\), the charge is positive and the bracket negative. Both regions contribute strictly negatively.

The endpoint conditions and an exact divided-difference calculation give

\[
F(x)=\frac12(x-a)(x-b)K(x).
\]

Consequently

\[
\boxed{
 d\rho(x)=-\frac{\sqrt{(x-a)(b-x)}}{2\pi}K(x)\,dx>0,
 \qquad a<x<b.
}
\tag{11}
\]

The density is bounded and vanishes like a square root at the endpoints. Positivity and finite logarithmic energy are now consequences of one sign-change identity, rather than numerical observations.

## 5. The same sign proves the global potential inequality

The arcsine Stieltjes transform is

\[
G_\mu(z)=\int\frac{d\mu(t)}{z-t}
 =\frac1{\sqrt{(z-a)(z-b)}},
\]

where the square root behaves as \(z\) at infinity. Partial fractions in (5) yield

\[
G_{B_u}(z)=\frac{1+C_uG_\mu(z)}{z+u^2}.
\]

Outside the interval, these identities imply

\[
\frac{d}{dx}(2U_\rho(x)-V(x))
 =\sqrt{(x-a)(x-b)}\,K(x).
\tag{12}
\]

Inside the interval, the real part of \(G_\mu\) is zero, so the potential difference is constant there. On \((0,a)\), the branch of the square root in (12) is negative; on \((b,\infty)\), it is positive. Since \(K<0\), the potential difference increases up to the support and decreases after it. By continuity,

\[
2U_\rho(x)-V(x)\le M\quad(x\ge0),
\]

with equality on \([a,b]\). This is the required global variational inequality. Equation (3) proves maximality of the measure; strict zero-mass energy negativity gives uniqueness among the relevant finite-energy competitors.

Thus a single interval really is enough. We have constructed it and verified the global condition, rather than assuming in advance that an equilibrium measure has one interval of support.

## 6. Endpoint existence is certified, not inferred from a root finder

The endpoint equations have exactly one solution in the rational box

\[
0.00007119508<a<0.00007119517,
\qquad
0.75394957<b<0.75394967.
\tag{13}
\]

Every decimal endpoint in (13) is an exact terminating rational number. The [endpoint certificate](equilibrium-endpoint-certificate.md) proves opposite signs on the **entire four faces**, not just at their corners. Its signs are

| Face | Certified sign |
|---|---|
| \(a=a_{\min}\), all allowed \(b\) | \(I_0<0\) |
| \(a=a_{\max}\), all allowed \(b\) | \(I_0>0\) |
| \(b=b_{\min}\), all allowed \(a\) | \(H>0\) |
| \(b=b_{\max}\), all allowed \(a\) | \(H<0\) |

On the curve \(I_0=0\), differentiation gives

\[
\partial_a I_0=-K(a)/2>0,\qquad
\partial_b I_0=-K(b)/2>0,
\]

and

\[
\frac{d}{db}H(a(b),b)=\frac{b-a}{4}K(b)<0.
\]

The sign-change argument also proves uniqueness of the \(a\)-root at fixed \(b\). Therefore the four face signs and ordinary one-dimensional intermediate-value arguments establish a unique pair in the box.

The [standard-library Python certificate](../computations/certify_equilibrium_endpoints.py) uses outward-rounded integer interval arithmetic, integer square roots, and Simpson remainders bounded with interval Taylor coefficients. Approximate roots were used only to select the box; they are not part of the proof.

## 7. Compute the energy without integrating a logarithmic singularity

Let \(\operatorname{cap}=(b-a)/4\),
\(V_0=-2+12\alpha(1-\log\alpha)\), and

\[
\kappa(u)=\log\frac{4u^2}
 {(\sqrt{u^2+a}+\sqrt{u^2+b})^2}.
\]

The plateau constant has the closed formula

\[
M=2\lambda\log\operatorname{cap}-V_0
  +\kappa(1)-6\alpha\kappa(\alpha)-4\lambda.
\tag{14}
\]

It follows from the explicit arcsine and swept-charge potentials, followed by integration by parts across the jumps of \(q\). The endpoint equations eliminate the remaining integral. This avoids numerical integration of \(\log(b-x)\) near an endpoint.

Since the potential is constant on the support,

\[
2I(\rho)-\int V\,d\rho=\lambda M,
\qquad
\mathcal J(\rho)=\frac{\lambda M-\int V\,d\rho}{2}+C_*.
\tag{15}
\]

Two different integral representations were used to certify the remaining moment.

The first uses (11), an angular substitution that removes square-root endpoints, and Carlson \(R_F,R_J\) expressions for \(K(x)\). The second uses the virial identity

\[
\int xV'(x)\,d\rho(x)=\lambda^2,
\]

obtained by symmetrizing the equilibrium force equation. Direct differentiation of (1) gives

\[
V(x)-2xV'(x)=\log(x+1)-6\alpha\log(x+\alpha^2)-2+12\alpha.
\]

Therefore

\[
\int V\,d\rho=2\lambda^2+U_\rho(-1)
 -6\alpha U_\rho(-\alpha^2)+\lambda(-2+12\alpha).
\tag{16}
\]

The off-support potentials reduce to smooth one-dimensional integrals. For \(s>0\), put

\[
D_s=\frac{(\sqrt{s+a}+\sqrt{s+b})^2}{4},
\qquad
\tau_s=\frac{\sqrt{s+b}-\sqrt{s+a}}
 {\sqrt{s+b}+\sqrt{s+a}}.
\]

Then the arcsine/Poisson-kernel identity gives

\[
U_\rho(-s)=\lambda\log D_s+
 \int_0^\infty q(u)\log(1-\tau_s\tau_{u^2})\,du.
\tag{17}
\]

This second route avoids the Carlson elliptic-integral reduction. It uses the same certified endpoint box and the same Arb arithmetic library, so the checks are independent integral representations, not wholly independent software stacks.

## 8. The rigorous numerical conclusion and what it improves

Both validated calculations prove the rational target

\[
\boxed{\mathcal J(\rho)<-\frac{69}{50}=-1.38.}
\tag{18}
\]

The saved enclosures are approximately

| Representation | Enclosure for \(\mathcal J\) |
|---|---|
| Carlson density integral | \(-1.3855\pm0.0000964\) |
| Virial/off-support potentials | \(-1.3855\pm0.0000295\) |

These are interval certificates, not agreement of sampled numerical approximations. Reproduce them with [the Carlson script](../computations/certify_equilibrium_energy.py) and [the virial script](../computations/certify_equilibrium_energy_elementary.py), using `python-flint`. Their saved outputs are [here](../computations/certify_equilibrium_energy.json) and [here](../computations/certify_equilibrium_energy_elementary.json). The computations use ball-arithmetic analytic quadrature and explicitly assert (18); they are not Lean formalizations.

The preprint’s exact arithmetic constant at cutoff \(200\) is

\[
A_{200}=
\frac{127125602969131786927559}
 {94195881588024216000000}.
\]

Exact rational arithmetic gives

\[
1600\left(\frac{69}{50}-A_{200}\right)
 =\frac{2864713622341631152441}{58872425992515135000}
 >48.
\tag{19}
\]

**Conditional consequence.** Retain the preprint’s integer-polynomial construction, its moment representation, and its separate arithmetic-normalization bound. Replacing its analytic comparison by the measure above gives

\[
\limsup_{n\to\infty}\frac1{n^2}
 \log Q_{40n,200}(\zeta(5))
 \le1600\left(A_{200}-\frac{69}{50}\right)<-48.
\]

Consequently, for all sufficiently large \(n\),

\[
0<Q_{40n,200}(\zeta(5))<e^{-48n^2},
\tag{20}
\]

conditionally on those retained arguments. The source’s stated bound is \(e^{-(139/5)n^2}\), whose exponent constant is \(27.8\).

This comparison concerns the limiting quadratic coefficient. It does not assert that the source’s particular finite remainder constants remain unchanged, nor provide an explicit starting index. The fixed bounded equilibrium density and the growth of \(V\) at infinity permit the same regularization-and-tail strategy with lower-order errors; the separate arithmetic and moment proofs still carry their own obligations.

The central simplification is conceptual: a signed field with one sign change determines a positive equilibrium density and its global potential shape. A large trial-measure table is replaced by two endpoint equations, a sign argument, and two reproducible scalar energy certificates. It is a local improvement to one analytic component, not a new fully independent elementary proof of \(\zeta(5)\)’s irrationality.

## 9. The same leading analytic budget works for every fixed odd index

The equilibrium construction has a useful extension that is stronger than a numerical observation but weaker than a general irrationality theorem.

Fix one odd integer \(s\ge3\). Put

\[
w_s(y)=\frac{2y^s}{(s-1)!}
 \frac{d^{s-1}}{dy^{s-1}}\frac1{e^{2\pi y}-1}.
\]

This weight is positive. On rational functions with polynomial parts and only simple poles at negative integer squares, its moment functional is affine in a single constant:

\[
\begin{aligned}
\mu_{s,X}(t^e)
 &=\frac{(-1)^eB_{2e+2}(2e+s)!}{(s-1)!(2e+2)!},\\
\mu_{s,X}\!\left(\frac1{t+j^2}\right)
 &=j^{s-1}\bigl(X-H_j^{(s)}\bigr)
   -\frac1{s-1}+\frac1{2j}.
\end{aligned}
\]

At \(X=\zeta(s)\), it is the positive integral \(\int_0^\infty R(y^2)w_s(y)\,dy\). The [moment theorem](odd-moment-generalization.md) proves the identity using Hermite’s integral formula.

Set \(K=40n\), \(N=3n\), \(h=37n\), and

\[
\mathscr D_m(t)=\prod_{j=1}^m(t+j^2),
\qquad
\Delta_{K,s}(X)=\det\left[
 \mu_{s,X}\!\left(\frac{\mathscr D_N(t)^6t^{i+j}}
                         {\mathscr D_K(t)}\right)
                         \right]_{0\le i,j<h}.
\]

Use the same rational normalizer as the source:

\[
S_K=\frac{(K!)^{2h}4^{h-1}}
 {(N!)^{12h}\prod_{i=1}^{h-1}((2i)!)^2}.
\]

**Fixed-odd analytic theorem.** The determinant has degree exactly \(h\), its value at \(\zeta(s)\) is positive, and

\[
\boxed{\quad
\limsup_{\substack{K\to\infty\\40\mid K}}
 K^{-2}\log\bigl[S_K\Delta_{K,s}(\zeta(s))\bigr]
 \le-\frac{69}{50}.
\quad}
\tag{21}
\]

The proof and [independent review](common-odd-zeta-budget-review.md) are in [the common-budget note](common-odd-zeta-analytic-budget.md). The reason is concrete: every fixed \(s\) changes only the polynomial factor in

\[
w_s(y)\le C_s(1+y)^s e^{-2\pi y}.
\]

After scaling \(y=K\sqrt t\), those changes contribute \(O_s(K\log K)\), while the common external field controls the \(K^2\) term. The proof also bounds the Riemann-sum error uniformly at \(t=0\), retains an integrable tail on the unbounded domain, and checks cancellation of every \(K^2\log K\) term. The conclusion is for each **fixed** \(s\); it is not uniform when \(s\) grows with \(K\).

The remaining arithmetic target can now be stated without ambiguity. One would need positive rational multipliers \(\nu_{K,s}\) such that

\[
\nu_{K,s}S_K\Delta_{K,s}(X)\in\mathbb Z[X],
\qquad
\limsup K^{-2}\log\nu_{K,s}<\frac{69}{50}.
\tag{22}
\]

Then positivity and the degree \(h=37K/40\) would give the polynomial integer contradiction. Equation (22) is a **sufficient condition to be proved**, not a result established for all odd \(s\). A common analytic field does not make the prime-adic arithmetic independent of \(s\).

## 10. An exact index-shift identity explains why the arithmetic does not transfer for free

There is an analytic operator shifting \(s\) to \(s+2\). With \(\theta=t\,d/dt\),

\[
\mu_{s+2}(R)=\mu_s(\mathcal L_sR),
\qquad
\mathcal L_s=
 \frac{(2\theta+s+1)(2\theta+s+2)}{s(s+1)}.
\tag{23}
\]

Two integrations by parts prove the identity. On polynomials it is particularly harmless:

\[
\mathcal L_s(t^e)=
 \frac{(2e+s+1)(2e+s+2)}{s(s+1)}t^e.
\]

On a simple pole, however,

\[
\mathcal L_s\!\left(\frac1{t+a}\right)
 =\frac{s-1}{s+1}\frac1{t+a}
 +\frac{2(2s-3)a}{s(s+1)}\frac1{(t+a)^2}
 +\frac{8a^2}{s(s+1)}\frac1{(t+a)^3}.
\tag{24}
\]

The repeated poles leave the domain of the original single-zeta arithmetic lemma. Their moment values contain additional zeta indices. The cancellation in the final analytic identity is valid, but it cannot be obtained by pretending that the old simple-pole coefficient bounds still apply unchanged.

There is a precise obstruction to avoiding this issue universally. Any finite-order differential operator over \(\mathbb Q(t)\) that preserves the space of rational functions with arbitrary polynomial part and only simple poles among \(-1,-4,-9,\ldots\) must be multiplication by a polynomial. The highest derivative would otherwise create a higher-order pole at infinitely many locations, forcing its rational coefficient to have infinitely many zeros. No such polynomial multiplier realizes the index shift: comparing large-degree moments would force a linear coefficient \(4\pi^2/[s(s+1)]\), which is not rational.

The arithmetic loss is also visible without an asymptotic guess. Write the rational constant in the simple-pole functional as

\[
h_s(j)=-j^{s-1}H_j^{(s)}-\frac1{s-1}+\frac1{2j}.
\]

For every prime \(p>s+1\),

\[
v_p(h_s(p+1))=-s,
\qquad
v_p(h_{s+2}(p+1))=-(s+2).
\tag{25}
\]

Thus increasing the index can genuinely cost two additional prime powers. These statements are proved in [odd-functional-shift.md](odd-functional-shift.md).

The scope matters. This rules out a cost-free transfer by one universal simple-pole-preserving rational differential operator. It does **not** rule out a \(K\)-dependent operator, determinant-specific cancellations, a different basis, or a fresh local arithmetic argument. Those are the appropriate next targets after the common analytic budget has been established.


## 11. A table-free arithmetic budget for the source's zeta-five construction

There is a second simplification, separate from equilibrium theory. It reduces the integration of the source's arithmetic envelope to variance bounds and fractional-part identities. Its precise status matters: the envelope and the passage from local prime valuations to its integral are separate from this integration calculation. The [local arithmetic reconstruction](fauzan-local-arithmetic-review.md) now checks the local inputs, including small primes and coefficientwise integrality, while [a separate review](fauzan-arithmetic-simplification-review.md) checks the envelope integration. The full assembly remains under review. Neither component extends the arithmetic normalization to every odd index.

For the source normalizer \(m_{K,M}\), take the fixed cutoff \(M=1000\). The [arithmetic note](fauzan-arithmetic-simplification.md) derives

\[
\boxed{\quad
\limsup_{\substack{K\to\infty\\40\mid K}}
 K^{-2}\log m_{K,1000}
 \le \frac{439979}{320000}
 =1.374934375<\frac{11}{8}.
\quad}
\tag{26}
\]

The key design choice is to cancel the large terms before estimating anything. In the source's inner envelope, the quadratic terms cancel exactly, leaving

\[
R(x)=xF(x)+Q(x),\qquad
F(x)=4\lambda+2\lambda\{x\}-12\lambda\{\alpha x\},
\qquad Q(x)\le\frac{13}{8}.
\]

The bound on \(Q\) follows from the variance of a function taking only two adjacent integer values, followed by Cauchy–Schwarz. This is where a long subdivision into floor-function cases can be avoided: one bounds the fluctuation collectively rather than integrating every case exactly.

The remaining fractional parts have a closed integral. For an integer \(n\ge1\),

\[
\int_n^\infty\frac{\{x\}}{x^2}\,dx
 =H_n-\log n-\gamma.
\]

The convexity of \(1/x\) gives \(\gamma>5/3-\log3\). Combining this with the source's outer envelope, integrated using positive-part primitives, yields

\[
\limsup K^{-2}\log m_{K,M}
 \le \frac{127751}{96000}+\frac{13}{144}
 +\lambda\left(\frac53+8\alpha+12\alpha\log\alpha\right)
 +\frac{14\lambda}{M}.
\tag{27}
\]

Here \(\alpha=3/40\), \(\lambda=37/40\). Even the logarithm needs no decimal certificate: the positive series for \(\operatorname{atanh}\), truncated after five terms in total, proves

\[
\log\alpha
 =-2\operatorname{atanh}(1/4)-6\operatorname{atanh}(1/3)
 <-\frac{33559}{12960}.
\]

Substituting this into (27) proves (26) by exact fractions. The [standard-library computation](../computations/check_table_free_arithmetic.py) verifies the final rational arithmetic; the accompanying note proves the integral estimates. This replaces the source's 143-piece inner integration and 11-piece outer table by a weaker but sufficient bound.

The analytic and arithmetic constants then leave the exact margin

\[
\frac{69}{50}-\frac{439979}{320000}
 =\frac{1621}{320000}>\frac1{200}.
\tag{28}
\]

Thus, retaining the source's moment identities, local arithmetic envelope, prime-sum conversion, and integrality proposition, these two replacements give a negative quadratic exponent for its normalized zeta-five determinant. The tradeoff for the simpler arithmetic estimate is a larger fixed cutoff: the source requirement \(K\ge200M^2\) becomes \(K\ge200{,}000{,}000\). This is acceptable for an asymptotic argument but not a useful small-approximant computation.

The result is a substantial simplification of two quantitative components, with their dependencies exposed. It is not a complete independent audit of the source proof, and it does not remove the arithmetic gate in Section 9 for higher odd zeta values.


## 12. Spend a little analytic margin to simplify the local arithmetic too

The strongest simplification is not the smallest numerical constant. The extra analytic margin lets us remove one of the source's arithmetic optimizations as well as its integration tables.

In the outer-prime argument the matrix is written \(A+p^{-1}L\), with row weights \(w_i\le0\), entry bounds \(v_p(A_{ij})\ge w_i+w_j\), integral \(L\), and \(r=\operatorname{rank}L\). Complementary-minor expansion immediately gives

\[
v_p^{\mathrm G}\det(A+p^{-1}L)\ge2\sum_iw_i-r.
\tag{29}
\]

Here \(v_p^{\mathrm G}\) is the minimum valuation of the polynomial coefficients. A term using a \(k\)-minor of \(L\) vanishes if \(k>r\); otherwise the factor \(p^{-k}\) costs at most \(r\), and deleting nonpositive row weights cannot worsen the lower bound. This proof avoids counting zero-weight rows and avoids the source's finer rank correction.

The [local arithmetic review](fauzan-local-arithmetic-review.md) reconstructs the actual bases and local functional estimates needed to apply (29). The simplification increases the asymptotic arithmetic cost by exactly \(9/640\). Choose the fixed cutoff \(M=10000\), instead of \(1000\), to reduce the tail cost. The table-free calculation then gives, for the resulting positive rational normalizer \(\widetilde m_{K,10000}\),

\[
\widetilde m_{K,10000}S_K\Delta_{K,5}(X)\in\mathbb Z[X],
\qquad
\limsup K^{-2}\log\widetilde m_{K,10000}
 \le\frac{2203747}{1600000}=1.377341875.
\tag{30}
\]

The normalizer retains signed prime exponents: positive local divisibility permits division by a prime power. Taking only denominator-clearing exponents would discard useful arithmetic savings and would require a new budget calculation.

The exact remaining margin is

\[
\frac{69}{50}-\frac{2203747}{1600000}
 =\frac{4253}{1600000}>\frac1{400}.
\tag{31}
\]

For \(K=40n\), this corresponds to a leading negative exponent at most \(-4.253n^2\). Once the analytic, normalization, and positivity statements are assembled with their stated hypotheses, lower-order terms can be absorbed to give an eventual bound \(0<\widetilde m S_K\Delta_{K,5}(\zeta(5))<e^{-4n^2}\). The degree is \(37n\), so this decay is strong enough to beat the factor \(q^{37n}\) arising from a hypothetical rational value with denominator \(q\).

That last paragraph describes the assembly implication, which is being independently reviewed; it is not a declaration that every component of the source preprint has been certified. The source threshold is now \(K\ge200M^2=20{,}000{,}000{,}000\), and no explicit practical starting index for the final decay is supplied. The quantitative goal is eventual decay; the conceptual gain is that stronger real-variable control can pay for simpler local arithmetic.
