# One Laurent polynomial contains both the Apéry integers and the small integral

Status: direct identities and proofs below, independently reviewed in
[apery-laurent-independent-review.md](apery-laurent-independent-review.md).
Both the main exact script and the independent covariance/measure checks pass.
The exact Laurent model was proposed by the parent agent and derived here from
the binomial sum. The construction is presented as an explanatory connection,
not a claim of historical novelty. No existing note in `codex/problem131/`
contained this three-variable model when this investigation began.

Define
\[
\boxed{\quad
\Lambda(x,y,z)=
\frac{(1+x)(1+y)(1+z)\big((1+x)(1+y)+xyz\big)}{xyz}.
\quad}                                                     \tag{1}
\]
Let \(\operatorname{CT}\) denote the constant term in all three variables.
The central bridge is
\[
\boxed{\quad
A_n=\operatorname{CT}\Lambda^n,\qquad
I_n=\int_{(0,\infty)^3}\Lambda^{-n-1}\,
          \frac{dx\,dy\,dz}{xyz}.
\quad}                                                     \tag{2}
\]
Here
\(A_n=\sum_{j=0}^n\binom nj^2\binom{n+j}j^2\), and \(I_n\) is the
positive Beukers integral, normalized by \(I_0=2\zeta(3)\). The two objects
are powers and inverse powers of **the same positive Laurent polynomial**.
This gives a concrete reason that their exponential scales are reciprocal.

## 1. Deriving the model instead of guessing it

Represent three of the four binomial factors in a summand by coefficient
extraction:
\[
\binom nj=\operatorname{CT}_z(1+z)^n z^{-j},\qquad
\binom{n+j}j=\operatorname{CT}_x(1+x)^{n+j}x^{-j}.
\]
Use a second variable \(y\) for the other \(\binom{n+j}j\), leaving one
\(\binom nj\) outside. With \(A=(1+x)(1+y)\), the binomial theorem gives
\[
\begin{aligned}
A_n
&=\operatorname{CT}_{x,y,z}
 A^n(1+z)^n\sum_{j=0}^n\binom nj\left(\frac A{xyz}\right)^j\\
&=\operatorname{CT}_{x,y,z}
 \left[A(1+z)\left(1+\frac A{xyz}\right)\right]^n
 =\operatorname{CT}\Lambda^n.
\end{aligned}
\]
All operations concern finite Laurent polynomials, so no convergence issue is
involved. This is a reusable design procedure: leave one binomial coefficient
for the binomial theorem and realize the others by independent coefficient
extractions.

The expanded polynomial has 22 distinct monomials, all with nonnegative integer
coefficients, total coefficient sum \(\Lambda(1,1,1)=40\), and constant
coefficient 5. Every exponent vector is in \(\{-1,0,1\}^3\), and the support
contains zero and all six signed coordinate vectors. Thus \(A_n\) counts
closed walks of length \(n\) in \(\mathbb Z^3\) using 40 colored steps.
Integrality, positivity, \(A_0=1\), and \(A_1=5\) become immediate.

## 2. The exact reciprocal substitution, including its measure

Use distinct letters \(a,b,c\in(0,1)\) for the Beukers cube variables. Set
\[
a=\frac{x}{1+x},\qquad b=\frac{y}{1+y},\qquad c=\frac{z}{1+z}.
\tag{3}
\]
This is a bijection from the positive octant onto the open cube. Write
\(A=(1+x)(1+y)\) and \(B=xy\). Then
\[
Q=1-c+abc=\frac{A+Bz}{A(1+z)},
\]
\[
a(1-a)b(1-b)c(1-c)=\frac{Bz}{A^2(1+z)^2}.
\]
Consequently the positive kernel satisfies the exact identity
\[
\boxed{\quad
h(a,b,c):=\frac{a(1-a)b(1-b)c(1-c)}{1-c+abc}
=\frac1{\Lambda(x,y,z)}.
\quad}                                                     \tag{4}
\]
The Jacobian is \(da\,db\,dc=dx\,dy\,dz/[A^2(1+z)^2]\). Hence
\[
\frac{da\,db\,dc}{Q}
=\frac{dx\,dy\,dz}{A(1+z)(A+Bz)}
=\frac1\Lambda\frac{dx\,dy\,dz}{xyz}.
\tag{5}
\]
Applying (4)–(5) to
\[
I_n=\int_{(0,1)^3}\frac{h(a,b,c)^n}{Q}\,da\,db\,dc
\]
proves the second identity in (2), with the exponent **−n−1** and no missing
normalizing factor. Positivity makes the change of variables legitimate by
monotone convergence; finiteness also follows from the bound in Section 3.

The separate arithmetic calculation in `beukers-explained.md` identifies
\(I_n=2A_n\zeta(3)-v_n\), with \(v_n\in\mathbb Q\) and
\(\operatorname{lcm}(1,\ldots,n)^3v_n\in\mathbb Z\). That arithmetic
identification is not being inferred merely from the reciprocal geometry.

## 3. One elementary minimization proves both exponential bounds

For \(x,y,z>0\), write \(A=(1+x)(1+y)\), \(B=xy\),
\(q=\sqrt A\), and \(r=\sqrt B\). Expanding (1) as a quadratic reciprocal
expression in \(z\) gives
\[
\Lambda=A z+\frac{A^2}{Bz}+A+\frac{A^2}{B}
\ge A\left(1+\sqrt{A/B}\right)^2.
\]
Equality holds precisely when \(z=\sqrt{A/B}\). Since
\(x+y\ge2\sqrt{xy}\), one has \(q\ge1+r\). Therefore
\[
\sqrt\Lambda\ge q+\frac{q^2}{r}
\ge(1+r)+\frac{(1+r)^2}{r}
=3+2r+\frac1r\ge3+2\sqrt2.
\]
All equality conditions can hold simultaneously, uniquely at
\[
x=y=\frac1{\sqrt2},\qquad z=1+\sqrt2.
\tag{6}
\]
Thus
\[
\min_{(0,\infty)^3}\Lambda
=(3+2\sqrt2)^2=17+12\sqrt2=:\lambda,
\qquad \rho:=\lambda^{-1}=(\sqrt2-1)^4.
\tag{7}
\]
The substitution (3) sends (6) to the Beukers maximum
\(a=b=\sqrt2-1,c=1/\sqrt2\). This proves
\(\max h=\rho\) without a separate three-variable optimization.

Because the coefficients of \(\Lambda^n\) are nonnegative, its constant term
is bounded by its value at any positive point: each other monomial contributes
a positive number, while the constant term contributes \(A_n\). At (6),
\[
A_n\le\lambda^n.
\]
The real integral obeys
\[
0<I_n\le\lambda^{-n}I_0=2\zeta(3)\rho^n.
\tag{8}
\]
Thus the growth and decay bounds are produced by the same minimization. The
last equality uses the elementary geometric-series evaluation of \(I_0\).
Together with the separate arithmetic denominator estimate, (8) is sufficient
for the usual Beukers irrationality argument.

## 4. A two-line geometric reason for Lucas congruences

For every prime \(p\), every \(m\ge0\), and \(0\le r<p\), Frobenius gives
\[
\Lambda(x,y,z)^{mp+r}
\equiv\Lambda(x,y,z)^r\Lambda(x^p,y^p,z^p)^m\pmod p.
\]
Every exponent of \(\Lambda^r\) lies in \([-r,r]^3\). The only vector in
this box whose coordinates are all divisible by \(p\) is zero. Consequently
constant-term extraction separates exactly, giving
\[
\boxed{A_{mp+r}\equiv A_m A_r\pmod p.}                       \tag{9}
\]
Iteration gives the usual digit product. This proof uses the small exponent
box directly, rather than a summation involving carries. For example \(A_n\)
is always odd because \(A_0\equiv A_1\equiv1\pmod2\).

The support contains all signed coordinate vectors, so zero lies in the
interior of its Newton polytope. Since the polytope lies in \([-1,1]^3\),
zero is its unique interior lattice point. This places the model within the
hypotheses of Samol–van Straten's general Dwork congruence theorem. Their
primary paper, §3–4, was retrieved and inspected; (9) is proved directly here,
while their stronger prime-power congruences require their additional theorem.

## 5. The two nonzero critical values explain the discriminant 34

Here allow complex \(x,y,z\ne0\), but exclude critical points with
\(\Lambda=0\). Regard \(\Lambda=F(A,B,z)\), where
\(A=(1+x)(1+y)\), \(B=xy\). Then
\[
x\partial_x\Lambda-y\partial_y\Lambda=(x-y)F_A.
\]
If a critical point had \(x\ne y\), it would have \(F_A=0\), and then
\(x\partial_x\Lambda=BF_B\). But
\[
BF_B=-\frac{A^2(1+z)}{Bz}\ne0
\]
when \(\Lambda\ne0\). This contradiction proves \(x=y=u\).

The remaining \(z\)-equation is
\(z^2=(1+u)^2/u^2\). On its two branches, substitution gives
\[
\Lambda\left(u,u,\frac{1+u}{u}\right)
=\left(2u+3+\frac1u\right)^2,
\]
\[
\Lambda\left(u,u,-\frac{1+u}{u}\right)
=\left(1+\frac1u\right)^2.
\]
Since \(\partial_z\Lambda=0\) on these branches, differentiation with
respect to \(u\) tests the remaining critical equation. The second branch
has no allowable nonzero critical point. On the first, nonzero \(\Lambda\)
forces \(u^2=1/2\). Therefore the only nonzero critical points on the torus are
\[
(1/\sqrt2,1/\sqrt2,1+\sqrt2),\qquad
(-1/\sqrt2,-1/\sqrt2,1-\sqrt2),
\]
with values \(\lambda\) and \(\rho\), respectively.

The hypersurface \(1-t\Lambda=0\) can consequently develop an interior
critical point at exactly \(t=\rho,\lambda\). Their polynomial is
\[
(1-\lambda t)(1-\rho t)=1-34t+t^2.                          \tag{10}
\]
This explains the quadratic discriminant appearing in the Apéry differential
equation directly from the Laurent model. The constant 5 is already the
constant coefficient of \(\Lambda\). These observations do not, by themselves,
prove every coefficient or the order of the Picard–Fuchs equation, nor exclude
singular behavior at a toric boundary. The full recurrence is proved separately
by the finite-balance kernel argument. Here (10) supplies a precise geometric
motivation for its singular values, without pretending that singularity data
alone determine an arbitrary differential equation.

## 6. Sharp asymptotics from the same Hessian on two integration domains

This section is optional for irrationality, but makes the reciprocal behavior
particularly explicit. Put \(r_*=(1/\sqrt2,1/\sqrt2,1+\sqrt2)\) and write
\(\Lambda=\sum_{e\in\mathbb Z^3}c_e X^e\). Define a probability distribution
on its finite support by
\[
p_e=\frac{c_e r_*^e}{\lambda}.
\]
The critical equations say that its mean is zero. Its covariance matrix is the
Hessian of \(\log\Lambda\) in logarithmic coordinates at \(r_*\):
\[
\Sigma=\begin{pmatrix}
-8+6\sqrt2&-7+5\sqrt2&-2+3\sqrt2/2\\
-7+5\sqrt2&-8+6\sqrt2&-2+3\sqrt2/2\\
-2+3\sqrt2/2&-2+3\sqrt2/2&-1+\sqrt2
\end{pmatrix},
\qquad
\det\Sigma=34\sqrt2-48=\frac{2\sqrt2}{\lambda}.
\tag{11}
\]
It is positive definite: all three coordinate directions occur with positive
probability, as does zero, so no nonzero linear functional of the steps is
constant. The script verifies the matrix and determinant exactly.

For the torus integral, set
\(\varphi(\theta)=\Lambda(r_*e^{i\theta})/\lambda\). Fourier extraction gives
\[
A_n=\frac{\lambda^n}{(2\pi)^3}
          \int_{[-\pi,\pi]^3}\varphi(\theta)^n\,d\theta.
\]
The zero step and the three positive coordinate steps show that
\(|\varphi(\theta)|=1\) only at zero modulo \(2\pi\). Near zero,
\[
\log\varphi(\theta)=-\tfrac12\theta^T\Sigma\theta+O(|\theta|^3).
\]
Positive definiteness gives \(|\varphi(\theta)|\le e^{-c|\theta|^2}\) on a
small fixed neighborhood, and compactness gives a strict bound below one on
its complement. Rescaling \(\theta=v/\sqrt n\) and applying dominated
convergence therefore proves the Gaussian integral formula, including its
constant, without assuming a probabilistic local limit theorem:
\[
\boxed{\quad
A_n\sim\frac{(1+\sqrt2)^2}{2^{9/4}\pi^{3/2}}
             \frac{\lambda^n}{n^{3/2}}.
\quad}                                                     \tag{12}
\]

For the real integral, put \(X=r_*e^v\). The logarithmic measure in (2) is
exactly \(d^3v\), and
\[
\log\Lambda(r_*e^v)=\log\lambda+	frac12v^T\Sigma v+O(|v|^3).
\]
The minimum is unique. The signed coordinate monomials make \(\Lambda(r_*e^v)\)
coercive as \(|v|\to\infty\), so outside a fixed neighborhood it is at least
\(\lambda+\delta\). The corresponding integral is bounded by
\((\lambda+\delta)^{-n}I_0\). The same local Gaussian calculation now gives
\[
\boxed{\quad
I_n\sim\frac{2^{3/4}\pi^{3/2}}{(1+\sqrt2)^2}
             \frac{\rho^n}{n^{3/2}}.
\quad}                                                     \tag{13}
\]
Here the factor \(\lambda^{-1}\) from the exponent \(-n-1\) is essential.
Multiplying the constants yields the particularly simple relation
\[
\boxed{\quad A_n I_n\sim\frac1{2\sqrt2\,n^3}.\quad}         \tag{14}
\]
Thus the power \(n^{-3}\) in the product of the large denominator and small
error has a direct explanation: each of the two three-dimensional saddle
integrals contributes \(n^{-3/2}\).

## 7. The period viewpoint and the remaining arithmetic

As an analytic germ at zero,
\[
\sum_{n\ge0}A_nt^n=
\frac1{(2\pi i)^3}\int_{|X_j|=(r_*)_j}
 \frac1{1-t\Lambda(X)}\,\frac{dx\,dy\,dz}{xyz}.
\]
This is the standard torus period of a Laurent hypersurface family. By contrast,
\(I_n\) integrates inverse powers over the positive real domain. One should
not silently identify those integration domains or conclude their arithmetic
from the shared function. The Legendre integration-by-parts argument supplies
the rational linear-form identity and its lcm denominator bound.

The useful design picture is concrete: coefficient extraction gives integers;
the positive reciprocal kernel gives a small nonzero quantity; a rational
change of variables shows why the two rates fit together. The difficult
arithmetic bridge remains explicit, but the choice of analytic kernel is much
less mysterious once its reciprocal is recognized as the Laurent model of the
integer sequence.

Primary references and status:

- W. Zudilin, *An elementary proof of Apéry's theorem*, arXiv:math/0202159,
  supplies the classical Apéry binomial sum and recurrence used elsewhere here.
- F. Beukers, *A note on the irrationality of ζ(2) and ζ(3)* (1979), supplies the
  auxiliary-integral construction; retrieval and independent derivation are
  documented in `beukers-explained.md`.
- K. Samol and D. van Straten, *Dwork Congruences and reflexive Polytopes*,
  [arXiv:0911.0797](https://arxiv.org/abs/0911.0797), §3–4, gives the primary
  constant-term/period framework and the general congruence theorem. Retrieved
  files: `codex/sources/samol-vanstraten-09110797.pdf` and `.txt`.

The last paper confirms the general framework, not historical priority for this
particular displayed polynomial or reciprocal substitution. That provenance
has not been resolved here. All specific identities needed above are proved
in the note and checked in `check_apery_laurent_reciprocal.py`.
