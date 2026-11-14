# One Bernoulli functional, two geometries

Status: exact algebraic identities and local distribution mechanism derived below
from the primary preprint's Section 3, generalized to every fixed odd index.
The real moment representation was independently proved in odd-moment-generalization.md.
This is a conceptual bridge, not a proof of the required determinant valuation bounds.
The p-adic interpretation is standard; historical novelty is not asserted.

## 1. A generating function replaces a list of Bernoulli identities

Let \(s\ge3\) be odd and define a rational linear functional on polynomials by
\[
L(x^k)=B_k,\qquad
\tau_s(g)=\frac{L(g^{(s-2)})}{(s-1)!}.
\]
Here \(B_1=-1/2\). Its formal exponential generating function is
\[
\tau_s(e^{tx})=\frac{t^{s-1}}{(s-1)!(e^t-1)}.
\tag{1}
\]
All expressions are formal series; no analytic convergence assertion is needed.
Three consequences follow simply by multiplication of this generating function:
\[
\tau_s(g(x+1)-g(x))=\frac{g^{(s-1)}(0)}{(s-1)!},
\tag{2}
\]
\[
\tau_s(g(-1-x))=-\tau_s(g(x)),
\tag{3}
\]
\[
\tau_s(g)=p^{-(s-1)}\sum_{a=0}^{p-1}\tau_s(g(a+px)).
\tag{4}
\]
For (4), the geometric sum \(\sum_a e^{at}=(e^{pt}-1)/(e^t-1)\)
cancels the denominator after substitution. Comparing powers of \(t\) proves the
identities for every polynomial. The sign in (3) uses that \(s-1\) is even.

Equation (2) says this functional turns a finite difference into a derivative at
one boundary point. This is exactly the sort of boundary cancellation used in
Zudilin's telescoping recurrence, although the functions and resulting estimates
differ. It gives a reason to introduce Bernoulli numbers: their translation law
is the desired mechanism, rather than an unrelated special-number trick.

## 2. The harmonic pole values are forced by that translation law

Adjoin rational functions with simple integer poles. Set
\(\tau_{s,X}(1/x)=-X\), and require (2) at functions regular at zero and (3).
These rules force
\[
\tau_{s,X}\!\left(\frac1{x-r}\right)=H_{d(r)}^{(s)}-X,
\qquad
 d(r)=\begin{cases}r&r\ge0,\\-r-1&r<0.\end{cases}
\tag{5}
\]
For positive \(r\), apply (2) to \(g=1/(x-r)\): its right side is
\(-r^{-s}\), and its left side is the difference between the pole values at
\(r-1\) and \(r\). Starting with the pole at zero yields (5). Reflection gives
the negative integers. Polynomial division makes this a uniquely defined rational
functional on the stated domain. Thus the harmonic sums are consequences of the
boundary rule and the single free parameter \(X\).

## 3. Why negative squares and the factor x^s appear together

Define
\[
\mu_{s,X}(R)=\tau_{s,X}\bigl(x^s R(-x^2)\bigr).
\tag{6}
\]
This transforms simple poles at negative integer squares into simple integer poles.
For monomials,
\[
\mu_{s,X}(t^e)=
\frac{(-1)^eB_{2e+2}(2e+s)!}{(s-1)!(2e+2)!}.
\tag{7}
\]
For a simple pole, polynomial division gives
\[
\frac{x^s}{j^2-x^2}=
-\sum_{\ell=0}^{(s-3)/2}j^{2\ell}x^{s-2-2\ell}
-\frac{j^{s-1}}2\left(\frac1{x-j}+\frac1{x+j}\right).
\]
Only the highest polynomial term survives \(s-2\) derivatives. Using (5),
\[
\boxed{\quad
\mu_{s,X}\!\left(\frac1{t+j^2}\right)
=j^{s-1}(X-H_j^{(s)})-\frac1{s-1}+\frac1{2j}.
\quad}
\tag{8}
\]
This explains every correction term in the pole formula. In particular the
\(1/(2j)\) is the mismatch between \(H_j^{(s)}\) and \(H_{j-1}^{(s)}\),
not an independently guessed adjustment.

At \(X=\zeta(s)\), the same rational data have a positive real interpretation:
\[
\mu_{s,\zeta(s)}(R)=\int_0^\infty R(y^2)w_s(y)\,dy,
\qquad
w_s(y)=\frac{2y^s}{(s-1)!}
\left(\frac{d}{dy}\right)^{s-1}\frac1{e^{2\pi y}-1}>0.
\tag{9}
\]
Formula (9) uses Euler's even-zeta formula and Hermite's integral formula;
its full convergence and integration-by-parts proof is in the moment note.
Positivity follows immediately from the positive exponential series for the even
order derivative. Equations (7)--(8), rather than an analogy, identify the real
and arithmetic functionals exactly.

## 4. A p-adic integral interpretation explains residue-class subdivision

For a polynomial \(g\) with rational coefficients and any prime \(p\),
\[
L(g)=\lim_{m\to\infty}\frac1{p^m}
\sum_{a=0}^{p^m-1}g(a)\qquad\text{in }\mathbb Q_p.
\tag{10}
\]
For a monomial this follows from the Bernoulli polynomial formula for the power
sum: divide by \(p^m\); the constant term is \(B_k\), and all remaining terms
contain positive powers of \(p^m\). Rational coefficients do not affect the limit.
This polynomial version of the Volkenborn integral is sufficient here.
Consequently \(\tau_s\) is a normalized p-adic average of an \((s-2)\)-nd
derivative. Splitting the finite average into residue classes explains (4);
the chain rule accounts for its exponent \(s-1\).

This is an unexpected but concrete link: the positive real weight in (9) and
the p-adic averaging operation in (10) encode the same rational moments.
The first controls the size and sign of a determinant. The second controls its
denominators. Neither viewpoint replaces the other.

## 5. What changes when poles are present

For \(p>s\), the polynomial functional extends to restricted p-adic power series:
its monomial values have valuation at least \(-1\), by von Staudt--Clausen and
because \((s-1)!\) is a p-unit. For a pole \(u\) with \(v_p(u)<0\), its
analytic expansion therefore has value
\[
\tau_s^{\rm an}\!\left(\frac1{x-u}\right)=
-\frac1{s-1}\sum_{k\ge0}
\binom{k+s-2}{s-2}B_k\,u^{-k-s+1}.
\tag{11}
\]
The series converges p-adically, since the valuation of its terms tends to infinity.
Let
\[
C_{p,s}=\sum_{a=1}^{p-1}\tau_s^{\rm an}\!\left(\frac1{x+a/p}\right),
\qquad Y=p^sX+C_{p,s}.
\tag{12}
\]
Then \(C_{p,s}\in p^s\mathbb Z_p\). Indeed the \(k=0\) contribution contains
\(p^{s-1}\sum_{a=1}^{p-1}a^{-(s-1)}\), whose sum is zero modulo \(p\),
since \(0<s-1<p-1\). The \(k=1\) terms have valuation at least \(s\);
for \(k\ge2\), the one possible Bernoulli denominator gives valuation at least
\(k+s-2\ge s\). The fixed denominator \(s-1\) is a p-unit.

With integer poles evaluated by (5) at parameter \(Y\), and far poles by (11),
the distribution identity extends to the rational domain:
\[
\tau_{s,X}(g)=p^{-(s-1)}\sum_{a=0}^{p-1}
\tau_{s,Y}^{\rm ext}(g(a+px)).
\tag{13}
\]
To prove it, it suffices to check polynomials and simple integer poles.
Polynomials obey (4). At \(g=1/x\), the right side is
\(p^{-s}(-Y+C_{p,s})=-X\). Translation and reflection then determine every
integer pole, just as in Section 2. More explicitly, put
\(D(g)=p^{-(s-1)}\sum_a\tau_{s,Y}^{\rm ext}(g(a+px))\). Its
translation difference telescopes to
\(p^{-(s-1)}\tau_{s,Y}^{\rm ext}(g(p(x+1))-g(px))\); the chain rule
cancels the factor \(p^{-(s-1)}\). Reflection follows by reindexing
\(a\mapsto p-1-a\). Hence \(D\) obeys exactly the same boundary and
reflection laws as \(\tau_{s,X}\), making the uniqueness argument applicable.
Translation of far-pole terms respects these identities because the functional
is bounded, with norm at most \(p\), on restricted power series equipped with
the Gauss norm. Translation and reflection are continuous there, so polynomial
approximation proves the required identities. This does not assert boundedness
on all continuous p-adic functions. Partial fractions complete the argument.
The full continuity and telescoping proof, with exact tests for
\(s=3,5,7,9\), is in
[bernoulli-dictionary-review.md](bernoulli-dictionary-review.md). This is the general-index
version of the preprint's distribution lemma, not a determinant estimate.

The powers \(p^{-(s-1)}\) and \(p^sX\) are forced by the boundary rule and
pole rescaling. They explain why changing the odd index affects local arithmetic,
even though the leading real energy estimate can remain unchanged.

## 6. Remaining construction problem

After (13), a useful determinant bound requires a basis whose zeros are distributed
among the residue classes. Local basis dimensions, ramification of \(x\mapsto-x^2\),
polynomial-part degree restrictions, and coefficientwise valuations still have to
be proved. The source's s=5 estimates cannot be transferred merely by writing
(13) for general s. The advantage of the dictionary is that it isolates precisely
which mechanism is universal and which estimates must be rebuilt.

Primary source: Fauzan, supplied September2026 preprint, Section3, saved in
codex/sources/fauzan-zeta5.pdf. Standard Bernoulli and Hermite identities are
referenced there and in odd-moment-generalization.md. No claim about the whole
preprint follows from the identities in this note.
