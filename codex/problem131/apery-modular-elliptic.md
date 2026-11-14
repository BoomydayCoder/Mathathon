# Apéry's recurrence from an elliptic period and two modular identities

Status: the eta-product identity and the algebraic hypergeometric pullback below are proved, with exact reproducible checks. The concluding modular irrationality proof is reconstructed from the primary source using an explicit covering argument and a direct three-integration calculation, with its standard modular-form inputs identified. It is an advanced explanation of classical structure, not a claim of historical novelty or a new proof of an open odd-zeta case.

Primary source: F. Beukers, *Irrationality proofs using modular forms*, Astérisque 147–148 (1987), 271–283, especially Propositions 1.1, 1.2, 2.1 and Theorem 1. [Original PDF](https://www.numdam.org/item/AST_1987__147-148__271_0.pdf); downloaded as `../sources/beukers-modular-1987.pdf` with text extraction beside it. The earlier author survey `../sources/beukers-caen.txt`, pp. 2–3, motivates the quadratic change of variable. A visual check of the original printed p. 275 confirms a small typographical error: its definition of F prints E₄(36τ), but E₄(6τ) is required by the immediately following Dirichlet series, the claimed level, and the Fricke transformation law. We use the corrected E₄(6τ) throughout and independently verify the relevant identities; this is not merely an OCR inference.

Reproduce the finite modular checks and exact rational differential substitutions with `python codex/computations/check_apery_modular.py`. The calculation is supplemented by the all-orders arguments here; coefficient matching alone would not establish these identities.

## 1. What the explanation should accomplish

The cubic coefficients in Apéry's recurrence look unnatural when one starts with integers. Two changes of viewpoint make them expected:

* A period of an elliptic curve satisfies a second-order differential equation. The square of that period satisfies a third-order equation. Third derivatives become cubic polynomials in the index when one takes coefficients.
* A modular parameter makes the arithmetic transparent: products with integer coefficients and triple integration produce integers and least-common-multiple denominators automatically. A modular symmetry then cancels the nearest analytic singularity, producing small linear forms.

Neither statement by itself proves irrationality. Below we exhibit the actual period, prove the parametrization, and specify the additional analytic input that completes the proof.

## 2. An exact elliptic hypergeometric identity

Let A(t)=Σ Aₙtⁿ, where A₀=1,A₁=5 and

(n+1)³Aₙ₊₁−(34n³+51n²+27n+5)Aₙ+n³Aₙ₋₁=0.

Define

T(x)=x(1−9x)/(1−x),   Z(x)=27x²(1−x)/(1−3x)³,

and H(z)={}_2F_1(1/3,2/3;1;z). Then, as analytic germs at zero (equivalently as formal series),

**A(T(x))=(1−x)/(1−3x)² · H(Z(x))².**

This is a concrete theorem, not just the observation that the Apéry equation resembles a period equation.

Proof. The symmetric-square computation in [symmetric-square.md](symmetric-square.md) proves A=u², where u(0)=1 and

(t³−34t²+t)u''+(2t²−51t+1)u'+(t−10)u/4=0.       (1)

Make the substitution u(T(x))=√(1−x)y(x). Direct differentiation transforms (1) into

x(x−1)(9x−1)y''+(27x²−20x+1)y'+(9x−3)y=0.       (2)

The exact constant −3 fills the unspecified constant in the author's brief survey. Equation (2) determines a unique power series with y(0)=1; its first coefficients are 1,3,15,93,639. More explicitly its recurrence is

(n+1)²bₙ₊₁−(10n²+10n+3)bₙ+9n²bₙ₋₁=0.

The hypergeometric differential equation is

z(1−z)H''+(1−2z)H'−(2/9)H=0.

Substitute z=Z(x) and y=H(Z(x))/(1−3x). Rational simplification gives exactly (2), and y(0)=1. Uniqueness proves the identity. The script checks both substitutions as identities of rational functions, rather than just checking series to finite order.

The quadratic substitution is motivated by the local exponents of (1). Its two nonzero finite singularities have exponents 0 and 1/2. A degree-two map ramified there doubles their exponent difference to 1 and removes those singularities after the appropriate change of dependent variable. Indeed

(1−34T(x)+T(x)²)=[(1−18x+9x²)/(1−x)]²,

so precisely the troublesome square root becomes rational. The inverse branch is

x=(1+t−√(1−34t+t²))/18.

This explains why that quadratic irrational, rather than a randomly chosen substitution, appears.

## 3. Why H is an actual elliptic period

For sufficiently small w, expand a uniformly convergent geometric series on |X|=|Y|=1:

P(w)=(2πi)⁻²∮∮ [1−w(X+Y+1/(XY))]⁻¹ dX/X dY/Y.

A constant term can occur only when the three selected monomials have equal multiplicities. Thus

P(w)=Σₙ≥₀ (3n)!/(n!)³ · w³ⁿ=H(27w³).

Integrating first in Y gives

P(w)=(2πi)⁻¹∮ dX/√[X²(1−wX)²−4w²X],

with the square-root branch chosen by continuation from small w and the lifted contour. The quartic under the square root has discriminant

−256w⁹(27w³−1).

For w≠0 and 27w³≠1 it has distinct roots, so its smooth double cover is an elliptic curve, and dX divided by that square root is a holomorphic differential. This proves the elliptic-period interpretation directly. At w=0 the curve degenerates, while its normalized period tends to 1.

Combining this with Section 2 realizes the Apéry generating function as a rational factor times the square of an explicitly specified elliptic period. The “cubic” recurrence is therefore compatible with a two-dimensional period space: symmetric products of two solutions form a three-dimensional space. This supplies a structural motivation for the recurrence, though integrality and identification of ζ(3) still need their own arguments.

## 4. The eta-product parametrization, proved with a small finite certificate

Put q=e^{2πiτ}, D=q d/dq, and η(τ)=q^{1/24}∏ₙ≥₁(1−qⁿ). Write η_d=η(dτ). Define

t=(η₁η₆/(η₂η₃))¹²,

E=(η₂η₃)⁷/(η₁η₆)⁵.

Then **A(t(q))=E(q)**. We give an all-orders proof using only standard basic modular-form facts and finite rational checks.

The advanced inputs are the eta transformation law and its usual eta-quotient criterion, together with the modular-form valence bound. They are independent mathematical theorems; omitting them would turn the finite computations into evidence rather than proof.

The group Γ₀(6) has index 12 in SL₂(Z). Its four cusps may be indexed by denominators c=1,2,3,6, with widths 6,3,2,1. The eta-quotient criterion shows t is a modular function with trivial character and E is a modular form of weight 2 with trivial character. The congruence checks for the exponent vectors (12,−12,−12,12) and (−5,7,7,−5) are respectively (24,24) and (0,0), both zero modulo 24. The character products are 1 and 6², respectively.

The cusp orders from the eta formula are

| c | 1 | 2 | 3 | 6 |
|---|---|---|---|---|
| ord(t) | 1 | −1 | −1 | 1 |
| ord(E) | 0 | 1 | 1 | 0 |

There are no zeros or poles of either eta quotient in the upper half-plane. Beukers uses Γ₁(6); for these even weights this makes no difference, since Γ₀(6)=±Γ₁(6).

Set h=Dlog t and G=Dlog E. The two crucial identities are

h²=E²(1−34t+t²),                              (3)

2DG−G²=t(10−t)E².                            (4)

Here is why a finite check proves them. The logarithmic derivative h is a holomorphic weight-2 modular form, including at every cusp. Although G itself is quasimodular, its transformation is

G(γτ)=(cτ+d)²G(τ)+2c(cτ+d)/(2πi).

Differentiating shows that 2DG−G² is a genuine weight-4 modular form. It is holomorphic in the upper half-plane because E has no zeros there. At a cusp, apply the same cancellation to the logarithmic derivative of E|₂σ; its local form is q_c^m times a nonvanishing analytic series, so the logarithmic derivative and its derivative are bounded. Thus this weight-4 expression is holomorphic there too.

The cusp table shows E²t and E²t² are holomorphic at every cusp: E's zeros cancel the poles of t. Therefore both differences in (3) and (4) are holomorphic weight-4 forms on Γ₀(6).

A nonzero such form has vanishing order at infinity at most 4, by the valence bound k[SL₂(Z):Γ₀(6)]/12. Consequently vanishing of coefficients q⁰ through q⁴ proves it is zero. One can also obtain this bound by taking the product over the 12 cosets to form a full-level weight-48 form and applying the full-level valence bound.

The exact eta expansions start

t=q−12q²+66q³−220q⁴+495q⁵+⋯,

E=1+5q+13q²+23q³+29q⁴+30q⁵+⋯.

The script performs the indicated finite checks using rational arithmetic. Both differences vanish through q⁴; the additional check through q²⁴ is redundant. The holomorphy and vanishing bound, not the extra terms, make the proof complete.

## 5. How the modular identities force Apéry's recurrence

Set d=1−34t+t² and θ=t d/dt=D/h. Differentiating (3) gives

Dh/h=G+t(t−17)h/d.

Let u=√E, with branch u(0)=1. Since Du/u=G/2, direct differentiation gives

θ²u/u=DG/(2h²)+G²/(4h²)−G Dh/(2h³).

Substituting the preceding derivative identity and (4) yields

d θ²u+t(t−17)θu+t(t−10)u/4=0,

which is precisely (1). Its normalized power-series solution is unique. Taking the symmetric square therefore proves E(q)=A(t(q)) and, on taking coefficients, proves the displayed cubic recurrence for the coefficients obtained from the modular products.

This explains the coefficients 34 and 10 in a compact way: they occur in two identities inside the small space of weight-4 modular forms. The coefficient 5 and the cubic recurrence then follow by differentiation and symmetric squaring. The modular setting does not make their discovery automatic, but reduces their verification to a short, justified finite calculation.

A further useful identity is

24E=−5E₂(τ)+2E₂(2τ)−3E₂(3τ)+30E₂(6τ).

The anomalous transformations cancel since −5+2/2−3/3+30/6=0. Both sides are holomorphic weight-2 forms; their coefficients through q² agree, which proves equality by the weight-2 vanishing bound 2. This connects our eta product exactly to Beukers's E in Theorem 1.

## 6. Arithmetic becomes easy in q, rather than t

Both products have integer coefficients, and t=q+O(q²). The inverse series q(t)=t+O(t²) also has integer coefficients: solve successively for each coefficient, whose coefficient in the inversion equation is 1. Thus E(q(t)) has integer coefficients. The modular parametrization alone gives Aₙ∈Z, without the binomial-sum expression.

For the companion linear forms define the holomorphic weight-4 form

F=[E₄(τ)−36E₄(6τ)−28E₄(2τ)+63E₄(3τ)]/40,

where E₄=1+240Σσ₃(n)qⁿ. Its constant term is zero and every other coefficient Fₙ is an integer. Put

f(q)=Σₙ≥₁Fₙqⁿ/n³,   so D³f=F.

If Lₙ=lcm(1,…,n), then the coefficient of tⁿ in E(q(t))f(q(t)) has denominator dividing Lₙ³. Indeed only indices j≤n enter, and each denominator j³ divides Lₙ³; composition and multiplication use integer coefficients. This is a much more natural origin for Lₙ³ than estimating denominators from the cubic recurrence, which would only give the much worse (n!)³ bound directly.

The first terms are

E(q(t))=1+5t+73t²+1445t³+⋯,

E(q(t))f(q(t))=6t+(351/4)t²+⋯.

The three integrations explain both the denominator exponent 3 and the eventual occurrence of ζ(3).

## 7. The analytic step: cancel the nearest singularity

The background facts used here are that the compact modular curve X₀(6) adds its four cusps to the quotient of the upper half-plane, and that the quotient is unramified away from elliptic points and cusps. We spell out the relevant covering and continuation argument, following the mechanism of Beukers's Proposition 2.1, rather than importing its picture as an unexplained assertion.

The Dirichlet series of F factors as

L(F,s)=6(1−6^{2−s}−7·2^{2−s}+7·3^{2−s})ζ(s)ζ(s−3).

At s=3 this is ζ(3), since the elementary factor is −1/3 and ζ(0)=−1/2. This evaluation uses analytic continuation of the Dirichlet series; it is not an assertion of absolute convergence of its defining sum at s=3.

The Fricke transformation τ↦−1/(6τ) satisfies

t(−1/(6τ))=t(τ),   E(−1/(6τ))=−6τ²E(τ),

F(−1/(6τ))=−36τ⁴F(τ).

The following special case of Beukers's Proposition 1.2 has a direct proof by three integrations, provided in [apery-modular-review.md](apery-modular-review.md), Section 4; no general Eichler-integral transformation theorem is necessary. It gives

6τ²[f(−1/(6τ))−ζ(3)]=−[f(τ)−ζ(3)].

Multiplying by E shows that W=E(f−ζ(3)) is invariant under this involution. It is also invariant under τ↦τ+1, from its q-series. We do **not** assume that W is invariant under all of Γ₀(6).

Here is an explicit proof that these two invariances give the required enlarged disk. Define

x(τ)=η₁⁴η₆⁸/(η₂⁸η₃⁴).

The eta criterion makes x a modular function whose cusp orders, in the order c=1,2,3,6, are (0,−1,0,1). It has no interior poles. A meromorphic function with exactly one simple pole has degree 1 as a map from the compact curve to the sphere; thus x is a global coordinate on X₀(6). The divisor of t then forces t, as a rational function of x, to have the shape x(1−a x)/(1−b x): its poles are the cusp x=∞ and one other simple finite pole, and one of its zeros is x=0. The leading coefficient is fixed by t∼x∼q. The expansions

x=q−4q²+10q³+⋯,   t=q−12q²+66q³+⋯

give b−a=−8 and 74−8b=66, hence a=9,b=1. Therefore t=T(x)=x(1−9x)/(1−x), exactly the quadratic substitution already used in the period calculation.

The quotient Γ₀(6)/{±1} is torsion free: an order-2 elliptic point would require a square root of −1 modulo 3, while an order-3 elliptic point would require a root of X²+X+1 modulo 2. Neither exists. Consequently the inverse τ(x) has no interior branching; its only excluded points are the four cusps. Under T those cusps map to 0 and ∞. The two remaining finite branch values are the critical values of T:

ρ=(√2−1)⁴,   ρ⁻¹=(√2+1)⁴.

The first is the value at the Fricke fixed point τ₀=i/√6. To select the sign without numerical approximation, note that x(τ₀)>0 and x(τ₀)<1. Indeed q=e^{−2π/√6}<1/8, and the product gives

x≤q exp[8q²/(1−q²)²+4q³/(1−q³)²]<q e^{1/4}<1/4.

Fricke invariance of t forces t'(τ₀)=0; x'(τ₀)≠0 by the unramified quotient. Thus x(τ₀) is the smaller critical point 1−2√2/3, and T(x(τ₀))=ρ.

On the disk |t|<ρ⁻¹ with 0 and ρ removed, the inverse continues along every path by this covering description. Loops around the two punctures generate its fundamental group. Starting from small positive t, a loop around 0 acts by the cusp translation τ↦τ+1, and a loop around ρ acts by the Fricke involution, with the lift reached along the imaginary axis. Both preserve W. The monodromy theorem therefore makes W single valued on the punctured disk. It extends over 0 by its convergent q-series. Near τ₀, the coordinate v=(τ−τ₀)/(τ+τ₀) changes sign under Fricke; both t−ρ and W are even analytic functions of v, and t−ρ has a nonzero quadratic leading term. Hence W is analytic in t at ρ as well.

This proves W is analytic throughout |t|<ρ⁻¹. Nonpolynomiality can be justified carefully as follows. The coefficients Fₙ/n³ are bounded, since σ₃(n)/n³≤ζ(3), so f(τ)=O(1/Im τ) as Im τ decreases to zero. If W were a polynomial in t, modularity of t and E would force g=f−ζ(3) to transform with weight −2 on Γ₀(6). At every cusp, g|₋₂σ then has at most polynomial growth in Im τ (the preceding bound gives O((Im τ)³)). A holomorphic periodic function with such growth has no negative Fourier powers, so it is holomorphic at the cusp. Thus g would be a holomorphic modular form of negative weight, hence zero by the valence bound. This contradicts g(i∞)=−ζ(3). Therefore W is not a polynomial. Its coefficients Bₙ−ζ(3)Aₙ consequently satisfy exponential bounds with any rate greater than ρ, and infinitely many are nonzero.

If ζ(3)=a/b were rational, bLₙ³(Bₙ−ζ(3)Aₙ) would be integers tending to zero along those nonzero coefficients. This is impossible using the elementary estimate log Lₙ≤c n+o(n), c≈1.10555043, proved in the companion Zudilin exposition: 3c+log ρ<0. The simpler bound Lₙ≤4ⁿ would not suffice. Thus the modular argument can avoid the prime number theorem at its final denominator-growth step, although Beukers's original presentation uses Lₙ=exp(n+o(n)).

The exact reason the error is small is now visible: a special constant, determined by a modular period integral, cancels the nearest singularity. This is the complex-analytic counterpart of constructing an auxiliary integral with a small positive kernel.

## 8. Scope and pedagogical recommendation

For a first complete proof of ζ(3) irrationality, the positive Beukers integral remains less demanding. The modular route substantially simplifies the conceptual explanation of the recurrence and denominator bounds, at the cost of modular transformations, the valence bound, and elementary covering theory for the modular curve. The particular triple-primitive transformation is proved directly in the independent review.

For the team's project, a manageable advanced module is Sections 2–5: prove the explicit hypergeometric identity by rational substitution, realize its hypergeometric function as an elliptic period, and prove the two weight-4 identities using the small vanishing bound. This produces concrete theorems and reproducible calculations, rather than only discussing a possible geometric connection.

For higher odd zeta values, merely replacing the triple primitive by a higher primitive changes the transformation's period polynomial: several critical L-values generally appear, and the simple scalar cancellation need not survive. The method therefore suggests a research question—find modular data with unusually sparse period polynomials and adequate denominator/decay balance—but does not prove ζ(5), or all odd zeta values, irrational.

## 9. A concrete obstruction to the direct ζ(5) analogue

The warning about extra moments can be made exact. Keep the same parameter t(τ) and try the weight-four prefactor E(τ)², so that its Fricke factor is 36τ⁴. The compatible weight-six Eisenstein input, with zero constant term at infinity, is unique up to scale in the span of E₆(dτ), d=1,2,3,6:

F₆=5E₆(τ)−248E₆(2τ)−837E₆(3τ)+1080E₆(6τ).

Here E₆=1−504Σσ₅(n)qⁿ. Its coefficients add to zero, and the exact Fricke permutation gives F₆(−1/(6τ))=216τ⁶F₆(τ). The four E₆(dτ) are linearly independent: inspect successively their q, q², q³ and q⁶ coefficients. The compatible Fricke eigenspace has dimension two, with a₆=216a₁ and a₃=(27/8)a₂. Zero constant term imposes the nonzero linear equation 217a₁+(35/8)a₂=0, leaving exactly one dimension and forcing the displayed vector. Thus the nonzero fourth critical value cannot be removed by retuning coefficients within this span. This is not a uniqueness assertion among all modular forms or all irrationality constructions.

Let f₆=Σ(F₆)ₙqⁿ/n⁵ and C=(2π)⁵/24. The same absolutely convergent integral argument gives, for v>0,

f₆(iv)=C∫ᵥ^∞(u−v)⁴F₆(iu)du.

Define Mⱼ=∫₀^∞uʲF₆(iu)du. Fricke substitution yields

Mⱼ=−6^{2−j}M₄₋ⱼ,   hence M₂=0, M₁=−6M₃, M₀=−36M₄.

The Dirichlet factorization is

L(F₆,s)=−504(5−248·2^{−s}−837·3^{−s}+1080·6^{−s})ζ(s)ζ(s−5).

Its two relevant values are

L(F₆,5)=−1526ζ(5),   L(F₆,4)=−840ζ(4)=−(28/3)π⁴.

Mellin integration therefore gives CM₄=−1526ζ(5), M₃=−7/2 and M₁=21. Substituting u↦1/(6u) into the fifth primitive and combining the two integration intervals proves the exact transformation defect

f₆(iv)−36v⁴f₆(i/(6v))
 =−1526ζ(5)(1−36v⁴)+14C v(1−6v²).                 (5)

The second term is nonzero. Subtracting the ζ(5) constant therefore does not produce Fricke invariance of E²(f₆+1526ζ(5)). In fact **no scalar constant c makes E²(f₆−c) Fricke invariant**: its defect is

(−1526ζ(5)−c)(1−36v⁴)+14C(v−6v³),

and these two polynomials are linearly independent. Since analyticity in the local coordinate t at the Fricke fixed point requires invariance under the local involution, no choice of c removes that branch completely. A suitable c can cancel its first odd Taylor coefficient, but cannot cancel all the odd coefficients; the radius cannot jump past the first branch by this scalar correction.

This is a proved obstruction for the most direct weight-six Eisenstein analogue. It isolates the precise extra quantity: the fourth critical L-value, here an explicit nonzero multiple of π⁴. It does not rule out adding other modular forms, using several linear forms or constants, choosing another modular group, or changing the approximation construction. Even a hypothetical scalar cancellation with the same radius and fifth-power lcm denominator would not by itself give the ζ(3) balance: the available inverse radius ρ⁻¹≈33.97 is smaller than e⁵≈148.41. That last comparison only diagnoses insufficiency of this naive estimate; it is not a lower bound excluding stronger arithmetic cancellation.

## 10. A finite-field consequence with an elementary proof

The same quadratic substitution T=x(1−9x)/(1−x) has involution M(x)=(x−1/9)/(x−1). This symmetry also determines the square factor of the Apéry digit polynomial modulo every prime p≥5:

Σₙ₌₀ᵖ⁻¹Aₙtⁿ = (1−34t+t²)^{(1−(−6/p))/2} Qₚ(t)².

A complete elementary proof is given in [apery-cm-character.md](apery-cm-character.md). It uses the second-order equation, a central-binomial sequence, Lucas congruences and two special values; the CM interpretation is optional. The character (−6/p) arises as the product of the involution's two signs, (2/p) and (−3/p), at its fixed points. Geometrically, those fixed points correspond to elliptic curves with CM discriminant −24. This is a concrete connection among the same recurrence, its elliptic reduction and prime arithmetic, rather than a conjectural analogy.
