# Coprime depletion: integral polynomial moments, but exact new pole costs

Status: exact identities and local determinant valuation theorem proved
below; independently generated finite checks pass in
../computations/check_coprime_depletion_poles.py. Independent adversarial
review is requested. No general failure theorem or irrationality claim.

This develops the parent's coprime-depletion proposal independently at
the pole and determinant level. The basic positivity/integrality theorem
is in the parent's depletion notes; the cyclotomic interpretation is in
cyclotomic-depletion-connection.md.

## 1. Normalize the real scale before comparing arithmetic

Fix odd \(s\ge3\), write \(m=s-1\), and let \(\mu_s\) be the original
positive Hermite functional with density \(w_s\). For coprime \(a,b>1\),
the positive depleted density is

\[
W(y)=ab\,w_s(y/(ab))-a\,w_s(y/a)-b\,w_s(y/b)+w_s(y).
\]

Its frequencies are exactly those divisible by neither \(a\) nor \(b\).
Its \(t^e\) moment is
\(((ab)^{2e+2}-a^{2e+2}-b^{2e+2}+1)\mu_s(t^e)\).
Equivalently, it is the product
\((a^{2e+2}-1)(b^{2e+2}-1)\mu_s(t^e)\).
The two single-depletion integrality results put this number in both
\(\mathbb Z[1/a]\) and \(\mathbb Z[1/b]\), hence in \(\mathbb Z\).

Restore the original exponential scale by defining

\[
\sigma(y)=\frac1{ab}W(ab\,y)
=w_s(y)-a^{-1}w_s(ay)-b^{-1}w_s(by)+(ab)^{-1}w_s(ab\,y)>0.
\]

Write \(\Sigma_s\) for integration against this density. Then

\[
\Sigma_s(t^e)
=(1-a^{-2e-2})(1-b^{-2e-2})\mu_s(t^e)
\in\mathbb Z[1/(ab)].
\]

At zero the density ratio \(\sigma(y)/w_s(y)\) tends to
\((1-1/a)(1-1/b)>0\), and at infinity it tends to one. It is positive
and continuous in between. Consequently there exist constants
\(0<c_{s,a,b}\le C_{s,a,b}<\infty\), independent of \(K,N,h\), such that

\[
c_{s,a,b}w_s(y)\le \sigma(y)\le C_{s,a,b}w_s(y).
\]

**Exact Gram comparison.** For any common positive rational weight and
any \(h\)-dimensional polynomial basis, their positive Gram determinants
satisfy

\[
c_{s,a,b}^{\,h}\Delta_\mu
\le\Delta_\Sigma\le C_{s,a,b}^{\,h}\Delta_\mu.
\tag{1}
\]

Proof: the pointwise inequalities imply the corresponding positive
definite matrix inequalities. Conjugate by the square root of the
original Gram matrix and multiply the resulting eigenvalue bounds.

In particular, when \(h=O(K)\), both determinants have the same leading
\(K^2\) logarithmic asymptotics. This is a theorem, requiring no new
equilibrium calculation. It does not compare their arithmetic normalizers.

## 2. Every integer pole still isolates one zeta value

Let \(H_j^{(s)}=\sum_{k=1}^j k^{-s}\) and set
\(C_s=(a^m-1)(b^m-1)\). Change of variables gives

\[
\Sigma_s((t+r^2)^{-1})
=\mu_s((t+r^2)^{-1})-\mu_s((t+a^2r^2)^{-1})
 -\mu_s((t+b^2r^2)^{-1})+\mu_s((t+a^2b^2r^2)^{-1}).
\]

Hermite's pole identity therefore gives exactly

\[
\boxed{
\begin{aligned}
\Sigma_s((t+r^2)^{-1})={}&C_s r^m\zeta(s)\\
&+r^m[-H_r^{(s)}+a^mH_{ar}^{(s)}
                  +b^mH_{br}^{(s)}-(ab)^mH_{abr}^{(s)}]\\
&+\frac{(a-1)(b-1)}{2ab\,r}.
\end{aligned}}
\tag{2}
\]

All four constants \(-1/(s-1)\) cancel. At \(a=2,b=3\), the final
constant is \(1/(6r)\), and \(C_s=(2^m-1)(3^m-1)\).
This agrees with evaluating the unscaled depleted functional at the pole
\((t+(6r)^2)^{-1}\).

The unscaled determinant with nodes \(6,12,\ldots,6K\), numerator power
six, and \(h=K-N\), equals \(36^{5Nh}\Delta_\Sigma\). Its primitive
integer polynomial is therefore the same up to sign as that of
\(\Delta_\Sigma\). Scaling alone supplies no arithmetic improvement.

## 3. Exact small-prime costs return after normalization

From now on \(a=2,b=3\). Put \(q=e+1\) and
\(B=\binom{2e+s}{s-2}\). Von Staudt–Clausen and elementary lifting of
exponents give

\[
\boxed{
\begin{aligned}
v_2(\Sigma_s(t^e))&=-2e+v_2(e+1)-v_2(s-1)+v_2(B),\\
v_3(\Sigma_s(t^e))&=-2e-2+v_3(e+1)-v_3(s-1)+v_3(B).
\end{aligned}}
\tag{3}
\]

Indeed, \(v_2(B_{2q})=v_3(B_{2q})=-1\),
\(v_2(9^q-1)=3+v_2(q)\), and
\(v_3(4^q-1)=1+v_3(q)\); the other two factors are units.
Both valuations are \(-2e+O_s(\log(e+2))\). These are individual-moment
identities, not determinant denominator lower bounds.

There is also an exact pole identity at two. If \(c_r\) is the rational
constant in (2), then

\[
\boxed{v_2(c_r)=m(1+v_2(r))-s\lfloor\log_2(6r)\rfloor.}
\tag{4}
\]

To prove it, rewrite the harmonic part as

\[
-(6r)^m\sum_{n=1}^{6r}
\frac{(1-2\,1_{2\mid n})(1-3\,1_{3\mid n})}{n^s}.
\]

The unique index having maximal 2-adic order is
\(n=2^{\lfloor\log_2(6r)\rfloor}\), whose numerator is \(-1\).
Other terms have strictly larger valuation. This leading term also has
strictly smaller valuation than \(1/(6r)\), since
\(\lfloor\log_2(6r)\rfloor\ge v_2(r)+2\).

## 4. An exact determinant obstruction at primes beyond all node collisions

Here is a stronger statement than the existence of bad individual poles.
Let

\[
D_j(t)=\prod_{r=1}^j(t+r^2),\qquad h=K-N>0,
\]

and define the affine functional \(\Sigma_{s,X}\) by replacing \(\zeta(s)\)
with \(X\) in (2). Set

\[
\Delta_{s,K,N}(X)
=\det_{0\le i,j<h}
\Sigma_{s,X}\!\left(\frac{D_N(t)^6t^{i+j}}{D_K(t)}\right).
\tag{5}
\]

The Gauss valuation \(v_p^G\) of a rational polynomial means the minimum
valuation of its coefficients.

**Theorem.** Suppose \(p>\max(3,2K)\) is prime, \(p\le6K\), and
\[
p\nmid(2^m-1)(3^m-1)(2^s-1).
\]
Then, with \(b=K-\max(N,\lfloor p/6\rfloor)\),

\[
\boxed{v_p^G(\Delta_{s,K,N})=-sb.}
\tag{6}
\]

Moreover the coefficient of \(X^{h-b}\) itself has valuation exactly
\(-sb\). Thus this is a coefficientwise necessary denominator cost,
not an upper-bound artifact or a numerical observation.

**Proof of the pole valuations.** For \(N<r\le K<p/2\), all harmonic
terms in \(H_r\) and \(H_{2r}\) are \(p\)-integral. Because \(6r<3p\),
the only possible nonintegral summands in the other tails are at \(p,2p\).
Modulo integral terms, the constant \(c_r\) is \(r^m p^{-s}\) times

\[
\begin{cases}
0,&r<p/6,\\
-6^m,&p/6<r<p/3,\\
3^m-6^m(1+2^{-s})=-3^m(2^s-1)/2,&p/3<r<p/2.
\end{cases}
\tag{7}
\]

There are no endpoint equalities since \(p>3\). Thus exactly \(b\) poles
have constant valuation \(-s\), the others have integral constants.
Every \(X\)-coefficient \(C_s r^m\) is a unit under the stated exclusion.

**Proof for the determinant.** Divide each rational integrand in (5)
into a polynomial and its surviving simple poles \(r=N+1,\ldots,K\).
All polynomial coefficients and their moments are \(p\)-integral:
the division is monic over \(\mathbb Z[t]\), and polynomial moment
denominators are supported only at two and three.

The residue at the pole \(-r^2\) is
\[
\rho_r(-r^2)^{i+j},\qquad
\rho_r=\frac{D_N(-r^2)^6}{D_K'(-r^2)}.
\]
Each \(\rho_r\) is a \(p\)-adic unit because every node difference
\((r-j)(r+j)\) has factors of absolute value less than \(p\). Let
\(V_{i,r}=(-r^2)^i\). This square \(h\times h\) Vandermonde matrix is
unimodular over \(\mathbb Z_p\), and the moment matrix has the form

\[
P+V\,\operatorname{diag}(\rho_r(C_s r^mX+c_r))\,V^T,
\]

with \(P\) a constant integral matrix. Conjugating by \(V^{-1}\) and
\(V^{-T}\) changes the determinant by a unit and leaves an integral
constant matrix plus the displayed affine diagonal.

Multiply each of the \(b\) bad rows by \(p^s\). Modulo \(p\), these
rows have only one nonzero entry, the corresponding diagonal constant.
Expansion along them leaves the complementary good block. Its
\(X^{h-b}\) coefficient is the product of its diagonal \(X\)-coefficients,
all units. Therefore the scaled determinant is integral and nonzero
modulo \(p\), including the specified coefficient. This proves (6).

For \(s=5\), the only fixed odd exceptions are \(5,31\); for \(s=7\)
they are \(7,13,127\). For sufficiently large \(K\), none lies in the band.

If \(N/K<1/3\) eventually, then \(b=K-\lfloor p/6\rfloor\) throughout
this band. Using the prime number theorem merely to summarize the size
of this **part** of the denominator gives

\[
\sum_{2K<p\le6K} sb\log p
=\left(\frac{4s}{3}+o(1)\right)K^2.
\tag{8}
\]

The exact finite theorem (6) does not use the prime number theorem.
Equation (8) is not a lower bound on the total magnitude of a primitive
rational normalizing scalar: common coefficient factors at other primes
can contribute in the opposite direction. Nor does it include any
additional factorial scalar chosen for the analytic normalization.

### 4a. Monic normalization turns the band cost into a global lower bound

The parent observed that the preceding caveat can be removed by normalizing
to a monic polynomial first. Here is the exact calculation.

The leading coefficient of (5) is

\[
\boxed{\ell_K=(-1)^{h(h-1)/2}C_s^h
\prod_{r=N+1}^K r^m D_N(-r^2)^5\in\mathbb Z\setminus\{0\}.}
\tag{9}
\]

Indeed, the coefficient matrix of \(X\) is
\(V\operatorname{diag}(\rho_r C_s r^m)V^T\).
Writing \(D_K'=D_ND_{\rm tail}'\) at the surviving nodes leaves power
five of \(D_N\). The product of the derivatives of the monic tail
polynomial at its \(h\) roots is
\((-1)^{h(h-1)/2}\det(V)^2\), canceling the Vandermonde square.
This also proves that (5) has degree exactly \(h\).

Set \(F_K=\Delta_{s,K,N}/\ell_K\), a monic rational polynomial, and let
\(M_K\) be the least common multiple of its coefficient denominators.
Then \(P_K=M_KF_K\) is primitive in \(\mathbb Z[X]\).
To check primitivity, if a prime divided every coefficient of \(P_K\),
then \(M_K/p\) would also clear all denominators, contradicting minimality.
Conversely every primitive integer representative of \(F_K\) is
\(\pm P_K\). In particular \(M_K\) is a positive **integer**; there are
no negative prime valuations available to offset the band cost.

For primes in (6), \(\ell_K\) is a unit: the factors \(r\), \(r-j\),
and \(r+j\) are all less than \(p\), and \(p\nmid C_s\).
Consequently

\[
v_p(M_K)=s[K-\max(N,\lfloor p/6\rfloor)].
\tag{10}
\]

For every asymptotic sequence \(N/K\to\alpha<1/3\), (8) therefore gives
the actual global lower bound

\[
\boxed{\liminf_{K\to\infty}\frac{\log M_K}{K^2}\ge\frac{4s}{3}.}
\tag{11}
\]

This corollary uses the prime number theorem. Indeed, replacing
\(\lfloor p/6\rfloor\) by \(p/6\) costs \(O(K)\) by the theorem, and
partial summation gives
\(\int_2^6(1-u/6)\,du=4/3\).
The fixed excluded primes are eventually outside the band.

Combining (11) with a lower bound on \(\log|F_K(\zeta(s))|\), rather
than on the unnormalized determinant, proves a failure theorem for this
particular family. The sign must be retained: \(\ell_K\) can be negative
although the raw integral determinant is positive. Absolute values make
the comparison unambiguous. The parent has now proved the requisite
Jensen lower bound in coprime-depletion-primitive-obstruction.md, and it
passed the independent review coprime-primitive-obstruction-independent-review.md.
Its conclusion is that at \(K=40n,N=3n\), for every fixed odd \(s\ge3\),
the positive primitive value eventually exceeds \(\exp(K^2/3)\).
That analytic input is not needed for (9)--(11) themselves.

## 5. Reflected-node divided differences also lose congruence

For \(p\ge7\), compare the poles \(r=1\) and \(r'=p-1\). The first is
\(p\)-integral. The second constant's leading \(p^{-s}\) coefficient is
\((p-1)^m E_s\), where

\[
E_s=3^m(1+2^{-s})+2^m
       -6^m\sum_{k=1}^{5}k^{-s}<0.
\]

For every prime not dividing the numerator or denominator of this fixed
rational number, the second constant has valuation exactly \(-s\).
The squared-node difference is \(p(p-2)\), so the constant coefficient
of the divided difference has valuation exactly \(-s-1\).

The \(X\)-coefficient does remain integral: \(m\) is even, so
\((p-1)^m-1\) is divisible by \(p\). Unlike the all-integer Euler model,
this failure comes from rational harmonic tails, not a sign mismatch.

## 6. What has actually improved, and the remaining design problem

All odd-prime polynomial corrections disappear completely. This is a
real algebraic simplification, and the normalized positive determinant
has exactly the same leading analytic decay as before by (1).

However, the harmonic tails create a large explicit new prime band, and
(6) shows that its cost survives the determinant, even when no two pole
nodes collide modulo the prime. A successful design must compensate for
this exact cost elsewhere or use a different pole basis/family.
For the specific ratio \(N/K=3/40\), the accompanying monic/Jensen
argument now rules out every rational scalar normalization: the actual
primitive values grow. Other ratios and constructions remain open.
