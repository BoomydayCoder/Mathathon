# Independent review of the modular and elliptic Apéry construction

Status: reviewed the proof and reran `check_apery_modular.py`; all checks pass.
No incorrect formula was found. The independent additional script
`check_apery_modular_review.py` checks the eta cusp table, the Fricke permutation,
the moment substitutions below, and the elliptic discriminant exactly. It also
passes. The derivation in Section 4 below removes the need to import Beukers's
general Eichler-integral transformation theorem for this particular example.
The global branching/continuation statement was initially imported from Beukers;
the subsequent direct proof is independently audited in Section 7 below.
No historical novelty is claimed.

## 1. Eta quotients and the finite certificate

For level \(N=6\), exponent vector \((r_1,r_2,r_3,r_6)\), and a cusp indexed
by denominator \(c\mid6\), the order in its local parameter is
\[
\operatorname{ord}_{c}(\prod_{d\mid6}\eta(d\tau)^{r_d})
=\frac6{24\gcd(c^2,6)}\sum_{d\mid6}\frac{r_d\gcd(c,d)^2}{d}.
\]
This gives precisely
\[
\operatorname{ord}(t)=(1,-1,-1,1),\qquad
\operatorname{ord}(E)=(0,1,1,0)
\]
at \(c=1,2,3,6\). The two congruence conditions modulo 24 and the square
character products are correct. Thus the usual eta-quotient theorem gives a
weight-zero modular function \(t\) and weight-two modular form \(E\), both
with trivial character on \(\Gamma_0(6)\). The eta function never vanishes in
the upper half-plane, so \(E\) has no interior zeros and \(t\) has no interior
zeros or poles.

Let \(D=(2\pi i)^{-1}d/d\tau\), \(h=D\log t\), and \(G=D\log E\).
The logarithmic derivative of a nonzero modular function is a weight-two
modular form. At each cusp its transformed function is a local power of the
cusp parameter times a holomorphic unit, so its logarithmic derivative is
bounded. This proves cusp holomorphy of \(h\), even at poles of \(t\).

For any nonvanishing weight-two form \(f\), define
\[
\mathcal Q(f)=2D(D\log f)-(D\log f)^2.
\]
Differentiating the weight-two transformation formula proves
\(\mathcal Q(f)|_4\sigma=\mathcal Q(f|_2\sigma)\). This covariance also
holds for cusp scaling matrices, not just matrices in the original group.
Since \(E|_2\sigma\) is a local power of the cusp parameter times a
holomorphic unit, \(\mathcal Q(E)\) is bounded at every cusp. Thus
\(2DG-G^2\) is holomorphic of weight four. This step is essential: the
quasimodularity of \(G\) alone would not supply cusp holomorphy.

The cusp orders of \(E^2t^j\), \(j=0,1,2\), are nonnegative. Both differences
in the proposed identities are therefore holomorphic weight-four modular forms.
The index is 12, so the valence bound gives order at infinity at most
\(4\cdot12/12=4\) for a nonzero form. Vanishing of the five coefficients
\(q^0,\ldots,q^4\) suffices; testing through \(q^{24}\) is redundant.
The main script correctly verifies these five coefficients using exact rational
arithmetic. The extra Eisenstein weight-two identity similarly needs only
coefficients through \(q^2\), with its anomalous transformation cancelled.

These proofs rely on the standard eta-quotient theorem and valence bound.
They do not rely on a conjecture inferred from matching many coefficients.

## 2. Hypergeometric pullback and elliptic interpretation

The script checks both transformed coefficients of the normalized second-order
ODE as rational-function identities. The first substitution uses
\(u(T(x))=\sqrt{1-x}\,y(x)\); the second uses
\(y(x)=H(Z(x))/(1-3x)\). The signs and gauge factors are correct. Although
\(Z'(0)=0\), the substituted hypergeometric series is analytic at zero and
has value one. The resulting equation's recurrence has nonzero leading
coefficient \((n+1)^2\), so the normalized analytic solution is unique.
There is no invalid inverse substitution at this ramified point.

The constant-term integral gives
\(\sum (3n)!w^{3n}/(n!)^3=H(27w^3)\). Integrating the inner variable yields
\[
\frac1{2\pi i}\oint\frac{dX}
{\sqrt{X^2(1-wX)^2-4w^2X}}.
\]
For small \(w\), the branch is fixed by continuation from its value at zero.
The quartic discriminant is exactly
\(-256w^9(27w^3-1)\). Away from the two indicated degenerations its smooth
compactification has genus one, and \(dX/\sqrt{\text{quartic}}\) is holomorphic.
The claimed elliptic-period interpretation is therefore explicit and valid.

The explanation of the quadratic substitution also checks: the finite nonzero
singularities of the second-order Apéry equation have local exponents
\(0,1/2\); the discriminant square becomes rational after \(t=T(x)\).
The verified differential identities make the removal of those singularities
precise, rather than relying on exponent-counting alone.

## 3. The Fricke identities can be derived directly

The eta transformation \(\eta(-1/z)=\sqrt{-iz}\,\eta(z)\), with its standard
holomorphic branch, exchanges the factors at \(d\) and \(6/d\). For each
paired product the transformation factor is \(-i\sqrt6\,\tau\). Hence
\[
t\left(-\frac1{6\tau}\right)=t(\tau),\qquad
E\left(-\frac1{6\tau}\right)=-6\tau^2E(\tau).
\]
The sign is fixed on the positive imaginary axis and then by analytic
continuation; it is not an arbitrary choice of square-root branches.

Write \(40F=\sum_{d\mid6}a_dE_4(d\tau)\), with
\((a_1,a_2,a_3,a_6)=(1,-28,63,-36)\). The full-level weight-four transformation
of \(E_4\) gives
\[
40F\left(-\frac1{6\tau}\right)
=\tau^4\sum_{d\mid6}a_d(6/d)^4E_4((6/d)\tau).
\]
The four exact identities \(a_{6/d}d^4=-36a_d\) prove
\[
F\left(-\frac1{6\tau}\right)=-36\tau^4F(\tau).
\]
This independently confirms the corrected argument \(6\tau\) in the source's
last Eisenstein term. Using \(36\tau\) would not produce this coefficient
permutation or a level-six form.

## 4. The triple-primitive transformation without a general Eichler theorem

Put \(\mathcal F(u)=F(iu)\) for \(u>0\), and
\(C=(2\pi)^3/2\). Because \(F\) has zero constant term at infinity,
\(\mathcal F(u)=O(e^{-2\pi u})\) as \(u\to\infty\). The Fricke identity gives
\[
\mathcal F(1/(6u))=-36u^4\mathcal F(u),
\]
so \(\mathcal F(u)=O(u^{-4}e^{-\pi/(3u)})\) as \(u\downarrow0\).
All moments below converge absolutely, including with any fixed real power.

The triple primitive has the elementary integral representation
\[
f(it)=C\int_t^\infty (u-t)^2\mathcal F(u)\,du.                 \tag{A}
\]
It follows either by integrating the Fourier series three times, uniformly on
\(u\ge t>0\), or by differentiation with the boundary values at infinity.
Let \(M_j=\int_0^\infty u^j\mathcal F(u)\,du\). The substitution
\(u=1/(6v)\) yields
\[
M_j=-6^{1-j}M_{2-j};\qquad M_1=0,\quad M_0=-6M_2.           \tag{B}
\]
For \(\Re z>4\), termwise Mellin integration gives
\[
\int_0^\infty u^{z-1}\mathcal F(u)\,du
=\frac{\Gamma(z)}{(2\pi)^z}L(F,z).
\]
The left side is entire in \(z\) by the two exponential endpoint bounds.
Analytic continuation and the explicit Dirichlet factorization therefore give
\[
CM_2=L(F,3)=\zeta(3).                                      \tag{C}
\]
Here the last equality uses the usual continued value \(\zeta(0)=-1/2\).
It is not a termwise evaluation of a non-absolutely-convergent Dirichlet series.

Apply \(v=1/(6u)\) to (A) at \(1/(6t)\). Direct simplification gives
\[
6t^2f(i/(6t))=-C\int_0^t(u-t)^2\mathcal F(u)\,du.
\]
Subtracting this from (A), and using (B)–(C), proves
\[
f(it)-6t^2f(i/(6t))
=C(M_2-2tM_1+t^2M_0)=\zeta(3)(1-6t^2).
\]
Both sides extend holomorphically to the upper half-plane, and the identity
theorem now gives exactly
\[
6\tau^2\left[f\left(-\frac1{6\tau}\right)-\zeta(3)\right]
=-[f(\tau)-\zeta(3)].                                      \tag{D}
\]
Thus the required special case of Beukers's Proposition 1.2 follows from
Fricke symmetry, a quadratic integral kernel, and three moments. The vanishing
of the middle moment is the concrete reason no extra period constant appears.
Multiplication by \(E\) proves the invariance of \(E(f-\zeta(3))\).

## 5. Initial dependency assessment, and the nonpolynomiality check

At the time of the first audit, the identities above did not by themselves
classify every branch of the inverse
modular parameter. The assertion that the relevant nonzero branching values
are \((\sqrt2-1)^4\) and \((\sqrt2+1)^4\), and that invariance removes the
first branch on the chosen continuation domain, still uses Beukers's
Proposition 2.1 and its continuation argument. That dependence is properly
identified in the main note. Local removal at the interior fixed point is
consistent with the analyticity of \(E\) and \(f\); the global radius requires
the covering description.

The added nonpolynomiality argument is sound. Uniformly in the real part,
\(|F_n|/n^3\) is bounded, so \(f(\tau)=O(1/\Im\tau)\) near the real axis.
If \(E(f-\zeta(3))\) were a polynomial in the modular function \(t\), then
\(g=f-\zeta(3)\) would transform with weight \(-2\). At any cusp represented
by \(\sigma\), the transformed function \(g|_{-2}\sigma\) has at most
polynomial growth as the imaginary part tends to infinity: the non-infinity
cusps give the explicit bound \(O(y^3)\). It is holomorphic and periodic there.
A negative Fourier power would grow exponentially; equivalently, estimating
its Fourier coefficient at height \(y\) and letting \(y\to\infty\) forces
that coefficient to vanish. Thus \(g\) would be holomorphic at every cusp,
contradicting the absence of nonzero holomorphic modular forms of negative
weight and \(g(i\infty)=-\zeta(3)\ne0\).

At that stage the modular irrationality proof was complete **with the explicitly
imported branching theorem**, the standard modular-form input, and the stated
elementary lcm estimate. Section 7 now audits the replacement of that imported
branching theorem by a direct proof. Its core recurrence, elliptic interpretation, arithmetic
integrality and special triple-primitive transformation are proved directly.

## 6. A useful presentation for the team

There are three distinct mechanisms, each solving a concrete difficulty:

1. Squaring an elliptic period explains why a third-order differential equation
   and cubic recurrence coefficients occur.
2. Integer eta products and an inverse series with linear coefficient one give
   integrality; triple integration gives the denominator \(L_n^3\).
3. Fricke symmetry makes one middle moment vanish and identifies the remaining
   constant as \(\zeta(3)\), cancelling the nearest branch singularity.

The third point is more concrete when presented through (A)–(D) than as an
unexplained invocation of an Eichler integral. It also pinpoints the obstacle
to generalization: a higher primitive has more moments, and one involution
need not remove all the additional independent constants.

## 7. Addendum: the direct covering proof is valid

The main note subsequently replaced its imported branching theorem by an
explicit modular-curve argument. I independently checked that argument; no
error was found. This supersedes the imported-branching dependency recorded in
Section 5, but still uses standard compact modular curves and covering theory.

The eta quotient
\(x=\eta_1^4\eta_6^8/(\eta_2^8\eta_3^4)\) has cusp orders
\((0,-1,0,1)\), trivial character, and no poles in the upper half-plane. Its
single simple pole gives degree one on the compact curve, so it is a global
coordinate. In particular the four cusps have distinct \(x\)-values. The
rational shape of \(t\) follows from its divisor, not from a finite series
conjecture: it has simple poles at the cusp \(x=\infty\) and one finite cusp,
and a simple zero at \(x=0\). The first three coefficients then determine
\(T(x)=x(1-9x)/(1-x)\) exactly. Its four cusp preimages are
\(x=0,1/9,1,\infty\), mapping to \(t=0,\infty\).

The torsion-free assertion for \(\Gamma_0(6)/\{\pm1\}\) is correct: an
order-two stabilizer would require a square root of \(-1\) modulo 3, and an
order-three stabilizer a root of \(X^2+X+1\) modulo 2. Thus the upper-half-plane
covering of the noncompact modular curve is unramified. Removing the critical
points of \(T\) gives a genuine covering over the sphere punctured at
\(0,\infty,\rho,\rho^{-1}\). There are no additional finite obstructions to
continuation hidden in a cusp or an apparent differential-equation singularity.
This covering description, rather than a differential-equation assertion about
apparent singularities, is sufficient for the proof.

The two critical points of \(T\) are
\(1\mp2\sqrt2/3\), with respective values \(\rho\) and \(\rho^{-1}\).
Both are simple critical points. The nontrivial degree-two deck transformation
in this coordinate is explicitly
\[
J(x)=\frac{x-1/9}{x-1},\qquad T(J(x))=T(x),\quad J(J(x))=x.
\]
It agrees with Fricke: the latter preserves \(t\) and is nontrivial on the
compact curve, since it exchanges the cusps zero and infinity.

The bound on \(x(i/\sqrt6)\) checks. Discarding numerator factors and using
\(-\log(1-z)\le z/(1-z)\) gives precisely the displayed exponential bound.
Its exponent at \(q=1/8\) is strictly less than \(1/4\), by exact rational
arithmetic. Thus the fixed point reached from the imaginary axis is the smaller
critical point. Along the intervening imaginary-axis segment, \(x\) is real
and its derivative cannot vanish, so this identifies the lift from the chosen
small positive base value, not merely an unrelated critical orbit.

For completeness, the monodromy generators have the needed specific lifts.
Take a base value \(0<t_*<\rho\), choose the inverse near \(q=0\), and choose
simple based generators in the disk \(|t|<\rho^{-1}\) punctured at
\(0,\rho\):

- The loop about zero continues the logarithm \(\tau=(2\pi i)^{-1}\log q\)
  by \(\tau\mapsto\tau+1\).
- Reach a neighborhood of \(\rho\) along the positive interval. Its lift
  approaches \(i/\sqrt6\), and the local interchange of the two sheets is
  exactly Fricke. Returning along the same interval therefore applies that
  involution to the original inverse germ.

Both operations preserve \(W\) as identities of holomorphic functions on the
upper half-plane. They consequently preserve every continuation obtained by
words in the two generators. This proves trivial monodromy of \(W\), without
assuming it is invariant under all of \(\Gamma_0(6)\).

The extension at zero follows from the original convergent series. At the
Fricke fixed point use \(v=(\tau-\tau_0)/(\tau+\tau_0)\), so Fricke is
\(v\mapsto-v\). Both \(W\) and \(t-\rho\) are even analytic functions of
\(v\); the latter has a nonzero quadratic coefficient. Inverting its analytic
series in \(v^2\) proves that \(W\) is analytic in \(t\) at \(\rho\).
Hence the enlarged disk is proved directly, with no additional analytic
continuation theorem imported from the primary paper.

## 8. Addendum: the restricted weight-six obstruction is valid

I independently checked all signs and normalizations in Section 9 of the main
note. In the four-dimensional span of \(E_6(d\tau)\), \(d\mid6\), Fricke
compatibility with the prefactor \(E^2\) requires
\[
a_6=216a_1,\qquad a_3=\frac{27}{8}a_2.
\]
The zero constant term then imposes
\(217a_1+(35/8)a_2=0\). The four Eisenstein functions are independent by
successively examining their first new coefficients. Thus the displayed vector
\((5,-248,-837,1080)\) is unique up to scale within this specified space.
No claim of uniqueness among cusp forms or other modular constructions follows.

On the positive imaginary axis its Fricke law has a **minus** sign:
\[
F_6(i/(6u))=-216u^6F_6(iu),
\]
since \(i^6=-1\). The absolutely convergent moments therefore satisfy
\[
M_j=-6^{2-j}M_{4-j}.
\]
The Dirichlet factors at 4 and 5 are respectively \(-20\) and \(-109/18\).
Using the Eisenstein coefficient \(-504\) and the continued zeta values gives
\(L(F_6,4)=-840\zeta(4)\) and \(L(F_6,5)=-1526\zeta(5)\).
The Mellin normalizations give exactly \(M_3=-7/2\), \(M_1=21\) and
\(CM_4=-1526\zeta(5)\), where \(C=(2\pi)^5/24\).

The fifth-primitive substitution has factor \(36v^4\), and its defect is
\[
A(1-36v^4)+14C(v-6v^3),\qquad A=-1526\zeta(5)-c.
\]
These two polynomials are independent, so no scalar \(c\) makes the defect
identically zero. If the resulting function were analytic in \(t\) at the
first branch value, it would be invariant under the local Fricke involution.
Its holomorphic defect would then vanish on a neighborhood and hence
identically. Thus the polynomial obstruction really excludes local branch
removal, rather than merely global modular invariance.

A sharper local check makes the failure concrete. At \(v_0=1/\sqrt6\), the
unique choice cancelling the first derivative is
\(A=-7C/\sqrt6\). The defect and its first two derivatives vanish there,
but its third derivative equals \(504C\ne0\). The best scalar choice leaves
a cubic local defect. This exact additional calculation is included in the
independent script.

This is an obstruction to this specific scalar-correction construction. It is
not an irrationality obstruction for \(\zeta(5)\), nor a statement that another
modular form, several constants, or stronger denominator cancellation cannot
work.
