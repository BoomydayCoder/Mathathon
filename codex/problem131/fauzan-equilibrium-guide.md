# Discovering an equilibrium measure for Fauzan’s analytic estimate

This note explains a concrete simplification of the analytic part of Fauzan’s 2026 \(\zeta(5)\) construction. An explicitly constructed equilibrium measure replaces the preprint’s mixture of sixteen arcsine measures and its many potential comparisons. A one-sign-change argument proves the global inequality; rigorous interval computations certify the two endpoint equations and the remaining energy bound.

The construction and certificates were developed and independently checked within this session. Historical novelty has not been established. This is **not an independent audit of the complete Fauzan preprint**. The final consequence for its integer polynomials remains conditional on the separate arithmetic and moment arguments. The source is [Fauzan’s preprint, version 1](https://doi.org/10.5281/zenodo.22826419), Sections 6–7 and Appendix A; detailed derivations are in [fauzan-analytic-simplification.md](fauzan-analytic-simplification.md).

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
