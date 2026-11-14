# Every fixed coprime depletion has the same new-prime obstruction

Status: proved extension of the exact \((2,3)\) local theorem. The
positive density and integral-moment construction are those in
cyclotomic-depletion-connection.md. The new contribution here is an
all-coprime-pair sign bound for the harmonic step values, an exact
coefficientwise determinant denominator theorem with fixed exceptions,
and its prime-sum consequence. No historical novelty is asserted.

Fix odd \(s\ge3\), coprime integers \(a,b\ge2\), and write
\[
m=s-1,\qquad q=ab\ge6,\qquad C=(a^m-1)(b^m-1).
\]
All of \(s,a,b\) remain fixed as \(K\) grows.

## 1. Normalized moments and poles

The positive density
\[
\sigma_{s;a,b}(y)
=w_s(y)-a^{-1}w_s(ay)-b^{-1}w_s(by)+q^{-1}w_s(qy)
\]
has the original real exponential decay \(2\pi\). Its polynomial
moments belong to \(\mathbb Z[1/q]\). For example, if
\[
c_k=\left.\frac{d^k}{dz^k}
\log\frac{(e^{qz}-1)(e^z-1)}
 {(e^{az}-1)(e^{bz}-1)}\right|_{z=0},
\]
then \(c_k\in\mathbb Z\) and
\[
\Sigma_{s;a,b}(t^e)
=\frac{(-1)^e c_{2e+2}}{q^{2e+2}}
    \binom{2e+s}{s-1}.
\]
The integrality follows from the integer cyclotomic quotient with value
one at one, as in the cited construction.

The pole value is \(Cr^mX+c_r\), with
\[
c_r=r^m[-H_r^{(s)}+a^mH_{ar}^{(s)}
                   +b^mH_{br}^{(s)}-q^mH_{qr}^{(s)}]
       +\frac{(a-1)(b-1)}{2qr}.                          \tag{1}
\]

## 2. A uniform positive step value

For \(0<x<1/2\), define
\[
R_s(x)=H_{\lfloor x\rfloor}^{(s)}
 -a^mH_{\lfloor ax\rfloor}^{(s)}
 -b^mH_{\lfloor bx\rfloor}^{(s)}
 +q^mH_{\lfloor qx\rfloor}^{(s)}.
\]
Here the first term is zero. If \(qx<1\), all terms are zero.
If \(qx\ge1\), then
\[
\frac{R_s(x)}{q^m}
\ge1-(a^{-m}+b^{-m})\zeta(s)
>\boxed{\frac{79}{144}}.                                \tag{2}
\]
Indeed, coprimality excludes \(a=b=2\), so
\(a^{-m}+b^{-m}\le1/4+1/9=13/36\).
Also \(\zeta(s)\le\zeta(3)<5/4\), the latter following from
\(\sum_{n\ge2}n^{-3}<1/8+\int_2^\infty x^{-3}\,dx=1/4\).
This proves (2), including its sign, for every allowed pair and every
odd \(s\ge3\).

Since \(a,b\mid q\), all breakpoints are multiples of \(1/q\).
The nonzero values are among the explicitly finite positive rational
numbers
\[
R_k=q^mH_k^{(s)}
      -a^mH_{\lfloor k/b\rfloor}^{(s)}
      -b^mH_{\lfloor k/a\rfloor}^{(s)},\qquad
1\le k\le\lfloor(q-1)/2\rfloor.                          \tag{3}
\]

Define the fixed exceptional prime set \(\mathcal E_{s,a,b}\) to contain
the prime divisors of \(2qC\) and the numerators and denominators in
lowest terms of all \(R_k\). It is finite because every \(R_k>0\).
This explicit set is sufficient, not necessarily minimal. Including
denominators avoids any hidden assumption that \(p>q\) in the finite
theorem.

## 3. Exact denominator theorem beyond node collisions

Let \(0\le N<K\), \(h=K-N\), \(D_j(t)=\prod_{r=1}^j(t+r^2)\), and
\[
\Delta(X)=\det_{0\le i,j<h}
\Sigma_{s;a,b,X}\!\left(D_N(t)^6t^{i+j}/D_K(t)\right).
\]
For every prime
\[
2K<p\le qK,\qquad p\notin\mathcal E_{s,a,b},
\]
put \(b_p=K-\max(N,\lfloor p/q\rfloor)\). Then
\[
\boxed{v_p^G(\Delta)=-sb_p,\qquad
v_p([X^{h-b_p}]\Delta)=-sb_p.}                           \tag{4}
\]

To prove this, separate the multiples of \(p\) in each harmonic sum:
\[
H_M^{(s)}=p^{-s}H_{\lfloor M/p\rfloor}^{(s)}
                   +\sum_{\substack{1\le n\le M\\p\nmid n}}n^{-s}.
\]
The second summand is \(p\)-integral. This equality is exact even when
some multiples of \(p^2\) occur; they remain inside the first summand.
For \(r\le K<p/2\), equation (1) therefore gives
\[
c_r=-r^mp^{-s}R_s(r/p)+(\text{\(p\)-integral remainder}).
\]
If \(qr<p\), it is integral. Otherwise \(R_s(r/p)\) is one of the
values (3), and is a unit by the fixed exclusions. Its valuation is
therefore exactly \(-s\). There are exactly \(b_p\) such surviving
poles. All pole \(X\)-coefficients \(Cr^m\) are units.

The rest is a full-matrix argument. The polynomial quotient has integral
coefficients, and its moments are \(p\)-integral since \(p\nmid q\).
Every pole residue factor
\(\rho_r=D_N(-r^2)^6/D_K'(-r^2)\) is a unit, and the square
Vandermonde \(V_{i,r}=(-r^2)^i\) is unimodular, because \(p>2K\).
Conjugate the moment matrix to
\[
P+\operatorname{diag}(\rho_r(Cr^mX+c_r)),
\qquad P\text{ constant and integral}.
\]
Multiply the \(b_p\) bad rows by \(p^s\). Modulo \(p\), they have only
nonzero diagonal constant entries. Expanding along them leaves a good
block whose \(X^{h-b_p}\) coefficient is the product of unit diagonal
\(X\)-coefficients. This proves both equalities (4), including the
case \(b_p=h\).

## 4. Monic normalization and the exact band cost

The leading coefficient is the nonzero integer
\[
\ell_K=(-1)^{h(h-1)/2}C^h
\prod_{r=N+1}^K r^mD_N(-r^2)^5.
\]
It is a unit at every prime in (4). Let \(F_K=\Delta/\ell_K\), and
let \(M_K\) be the least positive integer clearing all coefficients of
this monic polynomial. Then \(M_KF_K\) is primitive and
\[
\boxed{v_p(M_K)=s[K-\max(N,\lfloor p/q\rfloor)]}
\]
at those primes. Since \(M_K\) is an integer, other primes cannot
offset this contribution by negative valuations.

Suppose \(N/K\to\alpha\in[0,1)\). The prime number theorem, with the
finite exceptional set eventually outside the band, gives
\[
\liminf_{K\to\infty}\frac{\log M_K}{K^2}
\ge sB_q(\alpha),
\qquad
B_q(\alpha)=\int_2^q[1-\max(\alpha,u/q)]\,du.
\]
Computing the elementary integral yields
\[
\boxed{
B_q(\alpha)=
\begin{cases}
(q-2)^2/(2q),&0\le\alpha\le2/q,\\
(1-\alpha)\,[q(1+\alpha)/2-2],&2/q<\alpha<1.
\end{cases}}                                           \tag{5}
\]
The formulas agree at the breakpoint. The floor error contributes only
\(O(K)\); replacing \(N/K\) by its limit contributes \(o(K^2)\).
These asymptotic statements use the prime number theorem, whereas (4)
is an exact finite theorem.

For every \(q\ge6\),
\[
\boxed{B_q(\alpha)\ge B_6(\alpha).}                      \tag{6}
\]
Indeed, on \(2\le u\le6\), \(u/q\le u/6\), so the integrand for \(q\)
is at least the one for six. The additional interval \(6\le u\le q\)
has nonnegative integrand. In particular, at \(\alpha=3/40\),
every fixed coprime pair incurs band cost at least \(4s/3\).

## 5. Why no new real calculation is needed

The ratios of each normalized density to \(w_s\) tend to positive
constants at zero and to one at infinity; they are strictly positive
and continuous in between. Thus any two fixed coprime choices have
pointwise density ratios bounded above and below by positive constants.
Their raw Gram determinants, formed with the same \(D_N,D_K\) and
dimension \(h\), differ by factors between \(c^h\) and \(C_0^h\).
Their leading coefficients differ only by
\((C_{s;a,b}/C_{s;2,3})^h\). Therefore their **monic** values satisfy
\[
\log|F_{K;a,b}(\zeta(s))|
-\log|F_{K;2,3}(\zeta(s))|=O_{s,a,b}(h).
\]
For \(h=O(K)\) this has no effect on the \(K^2\) rate.

Consequently the separately audited monic Jensen lower bound for
\(K=40n,N=3n\), together with (6), proves the same primitive-scalar
obstruction for every fixed coprime pair \(a,b\ge2\).
This does not compare their primitive polynomials term by term, or
assert that their complete denominator lcms are ordered; it compares
the proved real rate and the sufficient band lower bounds.

The restriction that \(a,b,s\) be fixed is essential to this conclusion:
the density constants, finite exceptional primes, and asymptotic errors
have not been controlled uniformly for growing depletion parameters.

## Exact verification

The separate script check_general_coprime_depletion_band.py generates
polynomial moments from the integer cyclotomic quotient and pole constants
from the divisibility-weight sum. It checks 492 rational step values
against the uniform \(79/144\) bound, 84 exact piecewise band integrals,
and 41 actual determinant-prime cases for coprime pairs
\((2,5),(3,4),(3,5),(4,5),(4,7),(5,6)\).
All checks pass, including the specified coefficient and primitive monic
denominator exponents. These computations supplement the all-parameter
proof; they are not its basis.
