# Higher Rodrigues polynomials: a Laurent geometry and a sharp oscillating ζ(5) asymptotic

**Status:** exact construction, critical-point classification and sharp ζ(5)
contour asymptotic proved below. The contour asymptotic has passed independent
review in [higher-rodrigues-asymptotic-independent-review.md](higher-rodrigues-asymptotic-independent-review.md),
including direct Barnes-integral checks independent of the coefficient formula. Symbolic identities and high-precision corroboration are in
`../computations/check_higher_rodrigues_laurent.py`. No historical novelty is
claimed. This is not an irrationality proof: the available denominator bound
remains too expensive.

The starting polynomials, Mellin transforms and rational forms are established
in [multiplicity-lifted-zeta-family.md](multiplicity-lifted-zeta-family.md). The
new point is that their positive integer coefficient and their oscillating
small form can be understood through the same explicit Laurent polynomial.
For ζ(3) the relevant small critical value is real. Here the leading small
values are a complex conjugate pair, and a contour argument actually selects
them. Their presence alone would not have proved this selection.

**Primary-source context.** The underlying rational forms, finite Fourier/Barnes
decomposition and positive central integral are exact specializations of
Nesterenko’s published general constructions; the logarithmic polynomials
are confluent Jacobi–Piñeiro polynomials. See the explicit parameter dictionaries
in [rodrigues-primary-source-comparison.md](rodrigues-primary-source-comparison.md).
The specific geometric analysis and rate certificates here were independently
developed, but their historical priority remains unresolved.

## 1. An exact constant-term construction for every multiplicity

For a positive integer m put d=2m+1, introduce x₁,…,x₂ₘ,z, and set

\[
 C=\prod_{i=1}^{2m}(1+x_i),\qquad P=\prod_{i=1}^{2m}x_i,
 \qquad \boxed{\Lambda_m=(1+z)^mC\left(1+\frac C{zP}\right)^m.} \tag{1}
\]

This is a Laurent polynomial with positive integer coefficients. Expanding
the two factors that involve z gives

\[
 \operatorname{CT}_z\Lambda_m^n
 =\sum_{k=0}^{mn}\binom{mn}{k}^2C^{n+k}P^{-k}.
\]

Each remaining variable extraction is independent, so

\[
 \boxed{a_{n,m}:=\operatorname{CT}\Lambda_m^n
 =\sum_{k=0}^{mn}\binom{mn}{k}^2\binom{n+k}{k}^{2m}.} \tag{2}
\]

The coefficient of ζ(2m+1) in the previous pure form is exactly 2m aₙ,ₘ.
This explains the squared binomial sum without reverse-engineering a recurrence.
For m=1, (1) is precisely the previously studied Apéry Laurent polynomial.

## 2. Classify every critical point with nonzero critical value

All statements in this section concern the complex algebraic torus: every
coordinate is nonzero. A nonzero value further requires xᵢ≠−1, z≠−1, and
z+Q≠0, where Q=C/P. Logarithmic differentiation is therefore legitimate.
The z equation is

\[
 \frac z{1+z}=\frac Q{z+Q},\qquad\text{equivalently }z^2=Q.
\]

The xᵢ equation is

\[
 \frac{x_i}{1+x_i}-\frac{mQ}{(z+Q)(1+x_i)}=0.
\]

Consequently all xᵢ equal the same x, and

\[
 x=\frac{mz}{1+z},\qquad z=\frac x{m-x},\qquad
 \boxed{(m-x)(1+x)^m=\varepsilon x^{m+1},\quad\varepsilon\in\{1,-1\}.} \tag{3}
\]

Conversely every root of either equation (3) gives a nonzero torus critical
point by these formulas. None can equal 0, −1 or m, as direct substitution
shows; z=−1 would force m=0. The two signs cannot share a root. Thus this is a
complete classification of **nonzero** critical values, without a claim about
the possibly singular zero level. Their values are

\[
 \boxed{\Lambda_m=\left(\frac{m(1+x)}{m-x}\right)^{2m}.} \tag{4}
\]

There are exactly **2m nonzero critical points**, all nondegenerate. Indeed
let f(x)=(m−x)(1+x)^m/x^(m+1). Its logarithmic derivative is

\[
 \frac{f'(x)}{f(x)}=-\frac{m(m+1)}{x(1+x)(m-x)}.
\]

Thus all roots of f=±1 are simple. The plus equation has degree m+1.
The minus equation has degree m−1: both its leading degree-m+1 term and its
degree-m term cancel, and its next coefficient is m(m+1)/2. For m=1 it
is a nonzero constant and contributes no roots. The determinant formula (8)
holds algebraically at every complex critical point, and its factors are all
nonzero there. This proves nondegeneracy, not just simplicity of the reduced
equation. Critical values from different points need not be distinct for
arbitrary m; no such assertion is required here.

There is also an exact product law:

\[
 \boxed{\prod_{\text{plus branch}}\Lambda_m=1,\qquad
 \prod_{\text{minus branch}}\Lambda_m=m^{-4m}.} \tag{4a}
\]

Products count critical points, with the empty m=1 minus product equal to 1.
For the plus polynomial P₊ of degree m+1, evaluate P₊(−1)=(−1)^m and
P₊(m)=−m^(m+1). The product of (1+x)/(m−x) over its roots is therefore
m^(−m−1). Multiplying each factor by m and raising the product to 2m
proves the first identity. For P₋, of degree m−1, the evaluations are
P₋(−1)=(−1)^(m+1), P₋(m)=m^(m+1), giving the same unscaled ratio product
but two fewer factors of m. This proves the second identity. In particular
the reciprocal pair at m=1 and the real-times-conjugate-pair relation at m=2
are two instances of the same product law.

For m=1 the minus equation has no roots and the plus equation gives
x=±1/√2, recovering 17±12√2. For m=2 the two equations become

\[
 2x^3-3x-2=0\quad(\varepsilon=1),\qquad x=-2/3\quad(\varepsilon=-1).
\]

The latter gives z=−1/4 and critical value 1/256. Eliminating x from the
cubic and L(2−x)⁴=16(1+x)⁴ gives, exactly,

\[
 \operatorname{Res}_x(2x^3-3x-2,L(2-x)^4-16(1+x)^4)
 =16\bigl(256L^3-2036001L^2-17376L-256\bigr). \tag{5}
\]

The three values are one positive real number λ and a conjugate pair η,η̄:

\[
 \lambda=7953.13744063360\ldots,\qquad
 \eta=-0.00426719179966382+0.0103695523609027i\ldots.
\]

These decimals only identify the algebraic roots. The cubic's product of roots
is 1, hence

\[
 \boxed{|\eta|=\lambda^{-1/2}=0.011213\ldots.} \tag{6}
\]

The three values are distinct and the conjugate pair is nonreal: the value
cubic has negative discriminant (also exactly checked by the script). No
coincidence with 1/256 occurs, by substitution in (5).

## 3. Positive growth and the inverse real period

On positive coordinates let uᵢ=log xᵢ and w=log z. Then

\[
 \log\Lambda_m=
 m\log(1+e^w)+\sum_i\log(1+e^{u_i})+
 m\log\left(1+\exp\{\sum_i\log(1+e^{-u_i})-w\}\right).
\]

The first two terms have positive definite diagonal Hessian, and the last is
convex (a convex increasing softplus composed with a convex function).
Therefore log Λₘ is strictly convex. It is also coercive. Writing
S₊=Σmax(uᵢ,0) and S₋=Σmax(−uᵢ,0), it is at least

\[
 S_++m\max(w,0)+m\max(S_--w,0)
 \ge\min(1,m/2)(S_++S_-+|w|).
\]

There is exactly one minimum. Its x is the unique positive root in (0,m) of
the plus equation (3): the positive ratio
(m−x)(1+x)^m/x^(m+1) decreases from infinity to zero. Denote the minimum by λₘ.
Coefficient positivity, or the torus integral on the minimizing radii, gives

\[
 0<a_{n,m}\le\lambda_m^n. \tag{7}
\]

The positive real inverse period

\[
 I_{n,m}=\int_{\mathbb R^d}\Lambda_m(e^{u_1},\ldots,e^{u_{2m}},e^w)^{-n-1}
 \,du_1\cdots du_{2m}\,dw
\]

converges by the coercive bound, and 0<Iₙ,ₘ≤I₀,ₘ λₘ^(−n). It is not being
identified with the pure ζ(2m+1) form. For m=2 that identification would be
false: the pure form changes sign and has a much slower sharp decay.

Here the Gaussian constants are also explicit. At the minimum put x=xₘ and
let Σ be the Hessian of log Λ in the log coordinates. Its x-coordinate block,
cross entries and last entry are

\[
 \Sigma_{xx}=\frac{x}{1+x}I_{2m}+
 \frac{x(m-x)}{m(1+x)^2}{\bf1}{\bf1}^{T},\quad
 \Sigma_{xw}=\frac{x(m-x)}{m(1+x)}{\bf1},\quad
 \Sigma_{ww}=\frac{2x(m-x)}m.
\]

The block determinant is

\[
 \boxed{\det\Sigma=\frac{2(m+1)(m-x)}m
                 \left(\frac{x}{1+x}\right)^{2m+1}>0.} \tag{8}
\]

Ordinary Laplace expansion and the torus Gaussian expansion give

\[
 a_{n,m}\sim\frac{\lambda_m^n}{(2\pi n)^{d/2}\sqrt{\det\Sigma}},\qquad
 I_{n,m}\sim\frac{(2\pi/n)^{d/2}}{\lambda_m^{n+1}\sqrt{\det\Sigma}}.
 \tag{9}
\]

For completeness the compact saddle is unique: the support of Λ contains 0
and every positive coordinate unit vector; equality in the triangle inequality
forces all angles to vanish. The real saddle is unique by strict convexity,
and coercivity controls the tails. Thus no unproved saddle-selection
assumption is concealed in (9). In particular their product is asymptotic to
λₘ^(−1)(det Σ)^(−1)n^(−d).

## 4. Exact Mellin–Barnes decomposition of the pure ζ(5) form

Now set m=2 and use the previously proved normalization

\[
 \mathscr L_n=4a_{n,2}\zeta(5)-B_{n,2}
 =-\frac1{12\pi i}\int_{n-i\infty}^{n+i\infty}
   (\pi\cot\pi z)'''F_n(z)^2\,dz,
\]

where

\[
 F_n(z)=\binom{2n}n\frac{\Gamma(z)^3}
                {\Gamma(z-n)^2\Gamma(z+2n+1)}.
\]

Gamma reflection, including its sign which disappears after squaring, gives

\[
 F_n(z)^2=\binom{2n}n^2\frac{\sin^4\pi z}{\pi^4}H_n(z),\qquad
 H_n(z)=\frac{\Gamma(z)^6\Gamma(n+1-z)^4}{\Gamma(z+2n+1)^2}.
\]

The elementary identity
(πcot πz)''' sin⁴(πz)/π⁴=−4−2cos(2πz) therefore gives

\[
 \boxed{\mathscr L_n=\frac{\binom{2n}n^2}{6\pi i}
       \int_{c-i\infty}^{c+i\infty}(2+\cos2\pi z)H_n(z)\,dz,
       \quad0<c<n+1.} \tag{10}
\]

Moving the original line to any such c crosses no poles: the reflected
integrand is holomorphic in this strip and has exponential decay at imaginary
infinity, even after multiplication by cos. Thus (10) is an equality, not
merely a formal contour.

Define for j=0,±1

\[
 T_{n,j}=\frac{\binom{2n}n^2}{2\pi i}
       \int_{c-i\infty}^{c+i\infty}H_n(z)e^{2\pi i jz}\,dz.
\]

The integrals are independent of c within the strip and Tₙ,₊₁ and Tₙ,₋₁
are complex conjugates. Equation (10) becomes

\[
 \boxed{\mathscr L_n=\frac23T_{n,0}+\frac13\Re T_{n,-1}.} \tag{11}
\]

This decomposition will select the three small Laurent values exactly.

## 5. A global maximum calculation supplies the missing contour justification

All logarithms below are their principal values in the strip 0<Re t<1. Put

\[
 \psi(t)=6t\log t+4(1-t)\log(1-t)-2(t+2)\log(t+2),\quad
 f(t)=\frac{(1-t)^2}{t^3(t+2)}.
\]

On a fixed vertical line t=a+iy, 0<a<1, uniform Stirling gives

\[
 \binom{2n}n^2H_n(nt)
 =16\pi^3 n^{-3}f(t)
   e^{n(\log16+\psi(t))}\bigl(1+O(n^{-1})\bigr). \tag{12}
\]

The relative error is uniform over the entire line: each gamma argument has
positive real part at least a fixed positive multiple of n, and stays in the
closed right half-plane, where the standard complex Stirling remainder is
O(1/|argument|). The modulus of f(a+iy) is bounded and O(|y|⁻²) at infinity.

The decisive global identity is

\[
 \frac d{dy}\Im\psi'(a+iy)
 =\frac{12\bigl((3a+1)y^2+a(1-a)(a+2)\bigr)}
 {(y^2+a^2)(y^2+(1-a)^2)(y^2+(a+2)^2)}>0. \tag{13}
\]

Moreover Im ψ′ tends to ±4π at the two ends. Thus for each j=0,±1,
the real part of ψ(a+iy)+2πij(a+iy) has precisely one global maximum,
where Im ψ′=−2πj. It decreases at least linearly sufficiently far into either
tail (its limiting slopes are nonzero). Consequently there is a uniform
strict exponential gap outside every fixed neighborhood of its maximum.

For j=0 choose a=2/3. The complete complex derivative vanishes there, since
ψ′(2/3)=0. For j=−1 let t★ be the upper-half-plane root

\[
 2t^3-3t+2=0,\qquad
 t_\star=0.73784325889786+0.365017840801028i\ldots,
\]

and choose a=Re t★. Then ψ′(t★)=2πi. To check the branch without trusting
a numerical logarithm, exponentiation gives e^(ψ′/2)=−1 from the cubic.
For 0<a<1,y>0, (13) gives 0<Im ψ′<4π, while its real part is zero by the
modulus of that same algebraic equality. Hence Im ψ′ must be 2π.
The lower root similarly selects j=+1.

These choices lie strictly inside the strip: if r is the positive root of
2r³−3r−2, then 1<r<2, the other roots are −r/2±ib, and t★ has real part
r/2 in (1/2,1). Identity (13) also implies Re ψ″(t★)>0, so the Gaussian
has a unique square-root branch with positive real part.

Expanding locally on the vertical line, and using (12), (13) and the tail gap,
now gives ordinary Gaussian steepest descent with relative error O(1/n)
for each individual T. Odd local terms integrate to zero; bounded higher
derivatives and Gaussian domination control the next term. This is a proof
of saddle selection, not an inference from a list of algebraic critical values.

## 6. The sharp oscillating asymptotic and its geometric meaning

For a saddle t define

\[
 \mathcal C(t)=8\pi^2\sqrt{2\pi}\,
          \frac{f(t)}{\sqrt{\psi''(t)}},\qquad
 \psi''(t)=\frac6t+\frac4{1-t}-\frac2{t+2}, \tag{14}
\]

using the square root with positive real part. The preceding Gaussian
calculation, including dz=in dy, gives

\[
 T_{n,0}=\frac{\pi^{5/2}}{2\sqrt2}\,n^{-5/2}256^{-n}
              (1+O(n^{-1})),
\]

\[
 T_{n,-1}=\mathcal C(t_\star)n^{-5/2}\eta^n(1+O(n^{-1})),\qquad
 \eta=16\left(\frac{1-t_\star}{t_\star+2}\right)^4. \tag{15}
\]

Indeed at a full saddle the terms multiplying t cancel in
log16+ψ(t)+2πijt, leaving log16+4log(1−t)−4log(t+2).
For t=2/3 this is log(1/256). For t=t★ and x=−t★, it is exactly (4).
This proves that the same Laurent critical values control these contour
components, with all normalizations accounted for.

Combining (11) and (15) gives the sharp form

\[
 \boxed{\mathscr L_n=
 \frac{n^{-5/2}}3\Re\bigl(\mathcal C(t_\star)\eta^n\bigr)
       +O\bigl(n^{-7/2}|\eta|^n\bigr).} \tag{16}
\]

The real 256⁻ⁿ component is exponentially smaller and is absorbed in this
error. This comparison does not depend on the decimals: the positive cubic
root r lies in (1,3/2), so λ=[2(1+r)/(2−r)]⁴<10⁴<256² and therefore
|η|=λ^(−1/2)>1/256. If it is of interest separately, its contribution to (11) is
π^(5/2)/(3√2) n^(−5/2)256^(−n)(1+O(1/n)).

In particular

\[
 \mathscr L_n=O(n^{-5/2}\lambda^{-n/2}),\qquad
 \limsup_{n\to\infty}|\mathscr L_n|^{1/n}=\lambda^{-1/2}. \tag{17}
\]

There are infinitely many positive and infinitely many negative terms. Here
one need not prove that arg η/π is irrational. If its rotation is irrational,
the phases are dense on the circle. If rational, their finite orbit has at
least three points (η is not real); the real parts of C times that orbit sum
to zero and are not all zero, so both signs occur. Along corresponding
subsequences the main term stays a fixed fraction of |Cηⁿ| and dominates
the error in (16). This also proves the limsup assertion.

The independent review gives an optional algebraic strengthening: the positive
and negative terms each have natural density one half, and zeros have density
zero. Its Section8 proves arg η/π irrational using the S₃ splitting field of
2x³−3x−2 (the shifted polynomial is Eisenstein at3 and the only quadratic
subfield is Q(√−6)), then applies irrational-rotation equidistribution to (16).
The simpler infinite-sign argument above does not require this extra machinery.

The difference from the Apéry case is now precise. The integer coefficient
grows like λⁿ, the positive inverse real period decays like λ⁻ⁿ, but the pure
ζ(5) contour selects a conjugate pair of modulus λ^(−1/2). The coefficient
and small-form geometric construction still match; the relevant integration
cycles do not select reciprocal real values as they did for ζ(3).

The new sharp analytic logarithmic rate is (log λ)/2≈4.4907, improving the
previous proved log(729/16)≈3.819. The available sufficient normalization
D₂ₙ⁵ still has rate 10. Nothing here proves a lower bound of rate 10 for the
actual primitive denominators, nor excludes a better arithmetic construction.

## 7. What remains open in this construction

* Extend the exact geometric recurrence construction; the uniform contour
  asymptotic (10)–(16) has now passed independent review.
* Derive an exact recurrence for these periods and identify the boundary term
  appropriate to the pure form. Critical-value counts alone do not determine
  the recurrence order or prove that two sequences obey the same recurrence.
* Find substantially cheaper rational denominators or a different cycle/design;
  the new sharp decay alone does not close the irrationality budget.
* For m>2, analyze the corresponding higher cotangent derivative and reflection
  factors. The general Laurent classification does not by itself prove any
  higher-odd-zeta asymptotic or arithmetic conclusion.
