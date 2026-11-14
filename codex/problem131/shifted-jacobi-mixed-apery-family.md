# A shifted-Jacobi family: the same mixed telescope, kernel, and reciprocal rates

Status: independent development of the root agent's proposed shifted rational
kernel. General rational identities, diagonal elimination, finite arithmetic
samples, and differentiated-series samples are checked exactly/numerically by
`check_shifted_jacobi_mixed_recurrence.py`. The analytic and denominator proofs
are given below; independent review is requested. No historical novelty or
irrationality assertion for general Hurwitz zeta values is made.

The useful generalization is concrete: the auxiliary Legendre polynomials can
be replaced by shifted Jacobi polynomials, while the Beukers kernel itself
stays unchanged. The shift enters as a positive prefactor and modifies the
finite recurrence coefficients, but not their two exponential scales.

## 1. A single shifted rational kernel

Fix a real a≥0 and put
\[
 P_n(t)=\frac{\prod_{j=1}^n(t-j)}{\prod_{j=0}^n(t+a+j)},
 \qquad F_n(a)=-\sum_{k=1}^\infty(P_n(t)^2)'\big|_{t=k}.
\]
Empty products are1. These derivative series converge absolutely because
Pₙ(t)=t⁻¹+O(t⁻²) at infinity. Introduce
\[
 S_{nm}(t)=(-1)^{n+m}P_n(t)P_m(t),\qquad
 H_{nm}=-\sum_{k=1}^\infty S_{nm}'(k).
\]
Thus Hₙₙ=Fₙ. Let Δ denote forward difference in t.

There are two complementary linear functionals on these rational functions:

- the analytic derivative sum \(\mathscr D(R)=-\sum_{k\ge1}R'(k)\);
- the algebraic sum \(\mathscr C(R)\) of the coefficients of all double poles
  in R's partial-fraction expansion (simple poles and polynomial parts ignored).

Translation preserves the latter sum, so \(\mathscr C(\Delta R)=0\).
For the former,
\(\mathscr D(\Delta R)=R'(1)\) when R′(∞)=0. This explicitly separates the
formal residue calculation from the summation boundary calculation. It will
prove the same recurrence for both the growing coefficient and small error,
without assuming that the target Hurwitz zeta value is irrational.

## 2. The adjacent and mixed identities, with clean coefficients

For n≥1 write
\[
 K_n=6+\frac{2a^2}{n(n+a)},\qquad
 T_n(t)=\frac{(2n+a)^2}{n^2(n+a)^2}(t+n+a)^2P_n(t)^2.
\]
Then the proposed adjacent identity is exactly
\[
 \boxed{P_n^2+P_{n-1}^2+K_nP_nP_{n-1}=\Delta T_n.}      \tag{1}
\]
A two-line derivation explains the extra a² term. Divide by Pₙ₋₁² and use
Pₙ/Pₙ₋₁=(t−n)/(t+n+a). After multiplying by (t+n+a)², set
X=t(t+a), N=n(n+a). The left side becomes
\[
 (K_n+2)X+a^2+(2-K_n)N,
\]
whereas the right side becomes
\[
 \frac{(2n+a)^2}{N^2}\{X^2-(X-N)^2\}.
\]
Since (2n+a)²=4N+a², the expressions agree.

For the mixed telescope, define
\[
 \lambda_j=j(j+a+1),\quad d_n=2n+a+1,
\]
\[
 u_n=\frac{2(n+1)^2(n+a+1)^2}{2n+a+2},\qquad
 w_n=\frac{2n^2(n+a)^2}{2n+a},
\]
\[
 \boxed{v_{nm}=2d_n(\lambda_n+\lambda_m+a+1)-u_n-w_n.}   \tag{2}
\]
Then, with Qₙₘ=(t−m−1)(t+m+a)Sₙₘ,
\[
 \boxed{u_nS_{n+1,m}-v_{nm}S_{nm}+w_nS_{n-1,m}
                       =-2d_n\Delta Q_{nm}.}            \tag{3}
\]
This holds as a rational identity for n≥1, m≥0.

To derive it, the two n-shift ratios are
\[
 S_{n+1,m}/S_{nm}=-(t-n-1)/(t+n+a+1),
\quad S_{n-1,m}/S_{nm}=-(t+n+a)/(t-n).
\]
The t-shift ratio is
\[
 S_{nm}(t+1)/S_{nm}(t)=
 \frac{t^2(t+a)^2}
 {(t-n)(t+n+a+1)(t-m)(t+m+a+1)}.
\]
Multiplication by (t−m)(t+m+a+1) cancels the m-dependent denominator.
After putting both sides over (t−n)(t+n+a+1), the only remaining numerator is
\[
 t^2(t+a)^2-(t-m-1)(t+m+a)(t-n)(t+n+a+1).
\]
It is quadratic. Its t² coefficient is λₙ+λₘ+a+1, giving (2); its other
coefficients agree by direct expansion. The n±1 coefficients can also be
read off at the two denominator roots, rather than guessed.

### Boundary audit

For n≥1, Tₙ has a double zero at t=1, so Tₙ′(1)=0. For (3), if m≥1 the
product PₙPₘ has a double zero at1. If m=0 it has one zero and the extra
factor t−m−1 supplies the second. Thus Qₙₘ′(1)=0 in both cases. Tₙ and Qₙₘ
have finite limits at infinity and derivatives O(t⁻²). Differentiating
(1)/(3) and summing therefore kills their boundary terms legitimately.
Consequently both H and the analogous double-pole coefficient family satisfy
the mixed recurrence and
\[
 H_{n,n-1}=\frac{F_n+F_{n-1}}{K_n}.                       \tag{4}
\]
The sign in (4) comes from Sₙ,ₙ₋₁=−PₙPₙ₋₁.

## 3. The diagonal three-term recurrence and the a=0 specialization

Set
\[
 h_n=3n^2+3an+a^2,
\]
\[
 \alpha_n=\frac{(n+1)^3(n+a+1)^3}{(2n+a+2)h_{n+1}},
 \qquad
 \gamma_n=\frac{n^3(n+a)^3}{(2n+a)h_n},
\]
\[
 \boxed{\beta_n=v_{nn}-\alpha_n-\gamma_n.}               \tag{5}
\]
Substitute (4) and its n+1 version into the mixed recurrence at m=n. This gives
\[
 \boxed{\alpha_nF_{n+1}-\beta_nF_n+\gamma_nF_{n-1}=0,
 \qquad n\ge1.}                                        \tag{6}
\]
Here αₙ=uₙ/Kₙ₊₁ and γₙ=wₙ/Kₙ. All denominators displayed are nonzero for
n≥1, a≥0. This factored coefficient form is cleaner than a large expanded
polynomial numerator for βₙ.

For a=0, multiplying by6 gives exactly
\[
 (n+1)^3F_{n+1}-(2n+1)(17n^2+17n+5)F_n+n^3F_{n-1}=0.
\]
For every fixed a, the leading coefficients remain
\[
 \alpha_n/n^3\to1/6,\quad \beta_n/n^3\to17/3,
 \quad\gamma_n/n^3\to1/6.
\]
The characteristic polynomial is therefore X²−34X+1, independently of the
shift. This observation alone does not prove asymptotics; the common positive
kernel below does identify the two scales.

## 4. Shifted Jacobi polynomials give the same positive Beukers kernel

Define
\[
 L_n^{(a)}(x)=\frac{x^{-a}}{n!}
                \frac{d^n}{dx^n}[x^{n+a}(1-x)^n]
 =\sum_{j=0}^n(-1)^j\binom nj\binom{a+n+j}n x^j.        \tag{7}
\]
This is the shifted Jacobi Rodrigues polynomial (the formula fixes all signs
and conventions). Integrating by parts n times gives, for t+a>0,
\[
 \int_0^1x^{t+a-1}L_n^{(a)}(x)\,dx=(-1)^nP_n(t).        \tag{8}
\]
At zero, each boundary product vanishes as xᵗ⁺ᵃ; at1 Rodrigues's factor has
an n-fold zero. Formula (8) is also a direct beta-integral identity.

Let K(v)=−log(v)/(1−v). Squaring the moment, differentiating, and summing a
geometric series proves
\[
 F_n(a)=\int_0^1\!\int_0^1
 (xy)^aL_n^{(a)}(x)L_n^{(a)}(y)K(xy)\,dx\,dy.            \tag{9}
\]
All interchanges are absolutely dominated: the Jacobi polynomials are bounded
on [0,1], (xy)ᵃ≤1, and the K integral is finite.

Insert the elementary Stieltjes identity
\[
 K(xy)=\int_0^\infty\frac{du}{(1+u)(xy+u)}.
\]
Integrate by parts n times in x, then substitute u=yv before doing the y
integrations by parts. The powers of y cancel exactly as in the Legendre
case, leaving
\[
 F_n(a)=\int_{(0,1)^2}\!\int_0^\infty
 \frac{(xy)^a[x(1-x)y(1-y)v]^n}
 {(1+yv)^{n+1}(x+v)^{n+1}}\,dv\,dx\,dy.
\]
Finally set v=(1−c)/(yc). The resulting identity is
\[
 \boxed{F_n(a)=\int_{(0,1)^3}(xy)^a
 \frac{[x(1-x)y(1-y)c(1-c)]^n}
      {(1-c+xyc)^{n+1}}\,dx\,dy\,dc>0.}                \tag{10}
\]
The endpoint and absolute-convergence argument is the same as the elementary
all-q proof in `deformed-laurent-period-reduction.md`: the extra powers xᵃ,yᵃ
are bounded by1; Rodrigues derivatives through order n−1 still vanish at0
like positive powers, and at1 by the factor (1−x)ⁿ. After the first integration
by parts the absolute integrand is bounded by a constant multiple of the
initial integrable Stieltjes kernel. The substitution preserves this bound;
the second pointwise integration by parts gives a nonnegative integral bounded
by the preceding absolute integral. Thus no conditional Fubini operation is
hidden in (10).

Let ρ=(√2−1)⁴ and λ=ρ⁻¹=17+12√2. The already established sharp bound for the
unchanged Beukers kernel gives
\[
 0<F_n(a)\le2\zeta(3,a+1)\rho^n.                        \tag{11}
\]
Thus the same decay follows for all a≥0 directly from the auxiliary-integral
construction, rather than from recurrence asymptotics alone.

## 5. Residues identify the growing solution and the target

At t=−a−j the double-pole coefficient of Pₙ² is
\[
 A_{n,j}(a)=\binom nj^2\binom{a+n+j}n^2.
\]
Therefore
\[
 \boxed{A_n(a)=\mathscr C(P_n^2)
     =\sum_{j=0}^n\binom nj^2\binom{a+n+j}n^2.}          \tag{12}
\]
The simple-pole coefficient is
\[
 B_{n,j}(a)=2A_{n,j}(a)
 \left[H_j-H_{n-j}-\sum_{h=1}^n\frac1{a+j+h}\right].
\]
Since Pₙ²=O(t⁻²), one has ΣⱼBₙ,ⱼ=0. Summing the differentiated partial
fractions consequently cancels the Hurwitz ζ2 term and gives
\[
 \boxed{F_n(a)=2A_n(a)\zeta(3,a+1)-V_n(a),}             \tag{13}
\]
\[
 V_n=\sum_{j=0}^n\left[
 2A_{n,j}\sum_{h=1}^j\frac1{(a+h)^3}
 +B_{n,j}\sum_{h=1}^j\frac1{(a+h)^2}\right].
\]
For rational a, both Aₙ and Vₙ are rational. The formal double-pole functional
in Section1 proves recurrence (6) for Aₙ independently of any arithmetic
property of ζ(3,a+1). Subtracting then proves it for Vₙ as well.

Initial values, including every normalization, are
\[
 A_0=1,\quad V_0=0,\quad F_0=2\zeta(3,a+1),
\]
\[
 A_1=2a^2+6a+5,\qquad
 V_1=\frac{2(a+2)(a^2+3a+3)}{(a+1)^3}.
\]
At a=0 these are the familiar A₁=5 and V₁=12.

### Reciprocal scales and constants

The exact common kernel gives the stronger limits
\[
 \frac{F_n(a)}{F_n(0)}\longrightarrow(\sqrt2-1)^{2a},
 \qquad
 \frac{A_n(a)}{A_n(0)}\longrightarrow(1+\sqrt2)^{2a}.     \tag{14}
\]
For the first, the positive integrals concentrate at the unique Beukers
maximum x=y=√2−1, c=1/√2, already proved in the reciprocal-Laurent note.
Multiplying its integrand by the bounded continuous weight (xy)ᵃ multiplies
the leading Laplace constant by its value there. The complement of every
fixed neighborhood is exponentially negligible; the original n=0 integral
is an integrable dominating factor.

For the second, the unshifted nonnegative summands in (12) concentrate at
j/n=1/√2. Their consecutive-term ratio tends to
((1−θ²)/θ²)², giving the unique maximum; Stirling's formula gives an
exponential gap off any neighborhood of it. In such a neighborhood,
\[
 \log\frac{\binom{a+n+j}n}{\binom{n+j}n}
 =\sum_{h=1}^n\log\left(1+\frac a{j+h}\right)
 \longrightarrow a\log\frac{1+\theta}{\theta}
\]
uniformly. Everywhere this ratio is at most exp(aHₙ)≤eᵃnᵃ, so the polynomial
factor cannot upset the exponential tail gap. This proves the second limit.

In particular the exponential scales are λ and ρ for every fixed shift. Using
the already proved base asymptotics in `apery-laurent-reciprocal.md`, their
product has exactly the same leading constant:
\[
 \boxed{A_n(a)F_n(a)\sim\frac1{2\sqrt2\,n^3}.}          \tag{15}
\]
The shift changes the two leading constants by reciprocal factors. It is a
concrete deformation preserving the large/small-solution balance.

## 6. Exact denominator statements: integer and rational shifts differ

### Integer shift a=ℓ≥0

Every coefficient in (7), and every Aₙ(ℓ), is an integer. Set d_M=lcm(1,…,M),
d₀=1. Then
\[
 \boxed{d_{n+\ell}^{\,3}V_n(\ell)\in\mathbb Z.}         \tag{16}
\]
This is a proved sufficient denominator bound; least denominators are not
claimed. A direct polynomial-moment proof avoids a spurious index2n.
For i,j∈{0,…,n}, let
\(J_{ij}(a)=\iint x^{a+i}y^{a+j}K(xy)dxdy\). Then
\[
 J_{ii}=2\zeta(3,a+i+1),
\]
\[
 J_{ij}=\frac{\zeta(2,a+i+1)-\zeta(2,a+j+1)}{j-i}
       =\frac{\sum_{h=i+1}^j(a+h)^{-2}}{j-i}\quad(i<j).
\]
The off-diagonal identity follows by differentiating the elementary telescoping
partial fraction in the two integer-spaced shifts; it can also be obtained by
summing the positive geometric kernel. Diagonal subtraction from ζ(3,a+1)
uses finite cubic harmonic sums. For ℓ integer all harmonic denominators are
≤ℓ+n; the factor |j−i|≤n is also a divisor of dₙ₊ℓ. Thus dₙ₊ℓ³ clears every
rational correction in (9), proving (16).

This reproduces the classical denominator scale for every fixed integer shift.
Its target is ζ(3,ℓ+1)=ζ3−Hℓ⁽³⁾, so it supplies shifted auxiliary constructions
for the same irrational number, not a new independent Hurwitz irrationality
result.

### Rational shift a=p/q≥0 in lowest terms

Let
\[
 D_{q,n}=\prod_{\pi\mid q}\pi^{n v_\pi(q)+v_\pi(n!)},
\quad L_{p,q,n}=\operatorname{lcm}_{1\le h\le n}(p+qh),
\]
with empty lcm/product1. Then D_q,n is the exact least common denominator of
the Jacobi polynomial's coefficients. Indeed for π|q every numerator factor
p+q(i+h) is a π-adic unit, while the denominator of binom(a+n+i,n) is qⁿn!.
For π∤q, each block of πᵏ consecutive h supplies at least floor(n/πᵏ)
factors divisible by πᵏ, canceling vπ(n!). The coefficient at i=0 proves
sharpness at the primes dividing q.

A proved sufficient rational-moment normalizer is
\[
 \boxed{C_{p,q,n}=D_{q,n}^2 L_{p,q,n}^2
                       \operatorname{lcm}(L_{p,q,n},d_n).} \tag{17}
\]
It clears both Aₙ(a) and Vₙ(a). To see this, first clear the two polynomial
coefficients in each term of (9) by D_q,n². A diagonal correction has denominator
dividing L³; an off-diagonal correction has denominator dividing dₙL². Their
lcm divides the remaining factor in (17). No minimality is claimed for C.
For q=1 the familiar integer-shift bound (16) is a simpler sufficient bound;
(17) can be smaller and also divides d_(n+ℓ)^3.

This identifies an actual extra arithmetic cost for noninteger shifts, rather
than inferring useful integer forms from positive decay. For a proposed
irrationality application one would still have to show that the relevant
normalizer times Fₙ tends to0. That has not been proved here for arbitrary
p/q, and no blanket irrationality claim for ζ(3,a+1) follows.

## Reproduction and next steps

```sh
python codex/computations/check_shifted_jacobi_mixed_recurrence.py
python codex/computations/check_shifted_jacobi_mixed_recurrence.py --numeric
```

The checker proves the general rational identities symbolically, checks their
exact a=0 specialization, and derives A,V independently from finite Jacobi
polynomial moments for a=0,1,3,1/2,2/3 and n≤7. Both coefficient sequences
satisfy the displayed diagonal recurrence. It checks the denominator bounds
on these samples and the parameter-independent leading recurrence coefficients.
Optional differentiated-rational-series tests independently match the initial
linear forms at40-digit precision. The integral positivity and asymptotic
arguments are proofs in the note, not claims of script verification.

Useful next steps are a short direct constant-term model for the integer-shift
Aₙ(ℓ), and a sharper arithmetic analysis of selected small noninteger shifts.
The Jacobi deformation already supplies a motivated family with exact mixed
closure, positive errors, controlled integer-shift denominators, and reciprocal
asymptotic constants. Historical identification with existing shifted Apéry
families remains unresolved.
