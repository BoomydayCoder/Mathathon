# Independent review of the Bernoulli local/global dictionary

Status: the algebraic statements (1)–(13) in
[bernoulli-local-global-dictionary.md](bernoulli-local-global-dictionary.md)
were independently derived and checked. No incorrect formula was found. The
rational distribution argument needs the continuity and telescoping details
below; these provide a proof, not just a numerical check. This review does not
establish the determinant valuation estimates required for irrationality.

## 1. Exact symbolic and rational checks

Run `python codex/computations/check_bernoulli_dictionary.py` from the repository
root. The script uses rational arithmetic throughout; it explicitly uses
Bernoulli numbers with \(B_1=-1/2\).

For \(s=3,5,7,9\), it verifies:

- 60 monomials against the formal exponential generating function, translation,
  reflection, and distribution with subdivision factors 2, 3 and 5;
- the symbolic partial-fraction decomposition of \(x^s/(j^2-x^2)\), followed by
  16 exact pullback pole evaluations at \(j=1,2,3,5\);
- 152 instances of rational distribution modulo \(p^{10}\), for primes
  \(s<p<20\), including positive and negative pole indices crossing several
  residue-class blocks;
- the far-pole translation and reflection laws, and the valuation of
  \(C_{p,s}\), in those prime ranges.

All checks pass. Finite checks supplement the following proof, which applies to
every fixed odd \(s\ge3\).

## 2. The precise analytic extension

Let \(\mathbb Q_p\langle x\rangle\) denote the restricted power series
\(f(x)=\sum_{k\ge0}c_kx^k\) with \(c_k\to0\) p-adically. Give this algebra
the Gauss norm \(\|f\|=\sup_k|c_k|_p\). For \(p>s\),
\[
\tau_s(x^k)=
\begin{cases}
0,&k<s-2,\\
\displaystyle\frac1{s-1}\binom{k}{s-2}B_{k-s+2},&k\ge s-2.
\end{cases}
\]
The von Staudt–Clausen denominator theorem gives
\(v_p(\tau_s(x^k))\ge-1\). Therefore
\[
\tau_s(f)=\sum_{k\ge0}c_k\tau_s(x^k)
\]
converges and satisfies \(|\tau_s(f)|_p\le p\|f\|\). This is a bounded
functional on the restricted power-series algebra.

Translation \(f(x)\mapsto f(x+1)\) and reflection
\(f(x)\mapsto f(-1-x)\) have Gauss norm at most one: their coefficient formulas
use integer binomial coefficients, and the tails are uniformly bounded by the
original coefficient tails. They consequently extend continuously from
polynomials. Evaluation of any fixed derivative at zero is also continuous.
Thus the translation and reflection identities pass from polynomials to all
restricted power series by polynomial approximation.

This argument concerns analytic functions on the closed p-adic unit disc. It does
**not** claim that the Volkenborn functional is bounded on all continuous
functions on \(\mathbb Z_p\) with their uniform norm.

If \(v_p(u)<0\), then
\(1/(x-u)=-\sum_{j\ge0}x^ju^{-j-1}\) lies in this algebra. Substitution gives
\[
S_s(u):=\tau_s^{\rm an}\left(\frac1{x-u}\right)
=-\frac1{s-1}\sum_{k\ge0}
\binom{k+s-2}{s-2}B_k u^{-k-s+1}.
\]
In particular, continuity proves the exact identities
\[
S_s(u-1)-S_s(u)=-u^{-s},\qquad S_s(-1-u)=S_s(u).
\]
Both translated poles remain outside the unit disc.

For the script's poles \(v_p(u)=-1\), omitting all terms \(k\ge K\) introduces
an error of valuation at least \(K+s-2\). Rational distribution subsequently
divides by \(p^s\), so its final error has valuation at least \(K-2\).
Taking \(K=12\) therefore rigorously supports checks modulo \(p^{10}\).

## 3. A complete proof of rational distribution

Let \(\mathcal R\) be the space of rational functions with rational coefficients,
polynomial part arbitrary, and simple poles only at integers. The source
functional \(\tau_{s,X}\) acts on this space using the harmonic pole prescription.
On the target side allow both integer poles, with parameter \(Y\), and far poles,
with the analytic values above. All evaluations may be regarded as taking values
in \(\mathbb Q_p[X]\); in particular \(C_{p,s}\) is a p-adic constant, not an
asserted rational or real constant.

Set
\[
C_{p,s}=\sum_{a=1}^{p-1}S_s(-a/p),\quad
Y=p^sX+C_{p,s},\quad
D(g)=p^{-(s-1)}\sum_{a=0}^{p-1}\tau_{s,Y}^{\rm ext}(g(a+px)).
\]
For each integer pole of \(g\), precisely one of its transformed poles is an
integer and the others have valuation \(-1\). Hence every summand is in the
stated extended domain.

First, \(D\) agrees with \(\tau_s\) on polynomials by polynomial distribution.
At the pole zero, simple-pole rescaling supplies the additional factor \(p^{-1}\):
\[
D(1/x)=p^{-s}(-Y+C_{p,s})=-X.
\]
This explains why \(Y\) uses \(p^s\), although polynomial subdivision uses
\(p^{s-1}\).

For \(g\in\mathcal R\) regular at zero, the finite sum telescopes:
\[
D(g(x+1)-g(x))
=p^{-(s-1)}\tau_{s,Y}^{\rm ext}\bigl(g(p(x+1))-g(px)\bigr).
\]
Apply the extended translation identity to \(h(x)=g(px)\), which is regular
at zero. The chain rule gives \(h^{(s-1)}(0)=p^{s-1}g^{(s-1)}(0)\), so
\[
D(g(x+1)-g(x))=\frac{g^{(s-1)}(0)}{(s-1)!}.
\]
For reflection, reindex \(a\mapsto p-1-a\) and use
\[
-1-a-px=(p-1-a)+p(-1-x).
\]
The extended reflection law yields \(D(g(-1-x))=-D(g(x))\).

The boundary law, reflection and the single value \(D(1/x)=-X\) force
\[
D\left(\frac1{x-r}\right)=H_{d(r)}^{(s)}-X
\]
for every integer \(r\), exactly as in the original construction. Partial
fractions then prove \(D(g)=\tau_{s,X}(g)\) on all of \(\mathcal R\). There
is no tacit interchange of a divergent real series here.

Finally, \(C_{p,s}\in p^s\mathbb Z_p\): the \(k=0\) coefficient is a p-unit
multiple of \(p^{s-1}\sum_{a=1}^{p-1}a^{-(s-1)}\), whose residue is zero since
\(0<s-1<p-1\). The \(k=1\) terms already have valuation at least \(s\), and
for \(k\ge2\) the Bernoulli denominator bound gives valuation at least
\(k+s-2\ge s\).

## 4. What this explains conceptually

The useful starting point is a desired boundary operation:
**a discrete difference should become one derivative at one endpoint.** Its
exponential generating function forces the Bernoulli functional. Extending the
same rule across simple integer poles forces harmonic sums, leaving precisely
one free constant \(X\).

The substitution \(t=-x^2\) pairs positive and negative integer poles. The
factor \(x^s\) shifts the derivative index so that the resulting polynomial
moments use the even Bernoulli numbers. Polynomial division then forces the
correction terms \(-1/(s-1)+1/(2j)\); they are consequences of the design.

At \(X=\zeta(s)\), these exact rational moments also come from a positive real
integral. In the p-adic geometry, residue-class subdivision of the same
functional exposes denominator structure. Thus positivity and divisibility
are two realizations of identical rational data, connected by exact formulas.
The next difficult step is selecting a determinant basis that turns those local
identities into sufficiently strong valuation bounds. The dictionary itself
supplies no such bound and no irrationality proof for a new odd zeta value.
