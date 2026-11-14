# Catalan numbers clear the logarithmic endpoint with one lcm

**Status: complete elementary arithmetic proof and irrationality assembly,
accepted in [independent review](log2-catalan-independent-review.md).** The positive-measure and isometry inputs are
the independently reviewed
[logarithm–zeta deformation](logarithm-zeta3-deformation-independent-review.md).
The Wilson family and logarithmic irrationality are classical. This is a
reconstructed explanation of the arithmetic at its endpoint, with no
historical novelty claim.

The key discovery is an exact finite-sum factor: the associated numerator
contains a Catalan number times a beta integral of an integer polynomial.
One lcm clears that integral. At the coalesced zeta-three member, the same
finite-quotient calculation instead gives the familiar cubic-denominator
terms.

## 1. The endpoint and the theorem to prove

Fix \(n\ge1\), and write
\[
 J=2n+1,\qquad \ell_J=\operatorname{lcm}(1,\ldots,J),\qquad
 C_{n,k}=\binom nk\binom{n+k}k.
\]
At \(\tau=1/2\), use the positive probability measure
\[
 d\nu(y^2)=\frac{16y}{\sinh(2\pi y)}\,dy,
\]
and set
\[
 D_k(t)=\prod_{r=0}^{k-1}\bigl(t+(r+1/2)^2\bigr),\qquad
 \delta_k=(1/2)_k(3/2)_k
 =\frac{(2k)!(2k+1)!}{16^k(k!)^2}.
\]
The reviewed Hilbert Gram identity is
\[
 \nu(D_jD_k)=\frac{\delta_j\delta_k}{j+k+1}.
\]
Consequently
\[
 \Phi_n(t)=\sum_{k=0}^n(-1)^kC_{n,k}\frac{D_k(t)}{\delta_k}
\]
is orthogonal to all lower degrees and has squared norm \(1/(2n+1)\).
At \(s_n=-(n+1/2)^2\), its positive value is
\[
 A_n=\Phi_n(s_n)=\sum_{k=0}^n
 \frac{16^k}{2k+1}\binom{n+k}{2k}^{\!2}.
 \tag{1}
\]

Define the rational associated numerator and finite tail by
\[
 \mathcal N_n=\nu\!\left(\frac{A_n-\Phi_n(t)}{t-s_n}\right),\qquad
 T_n=\sum_{r=1}^n\frac1{r(r^2-1/4)}.
\]
The quotient in \(\mathcal N_n\) is a polynomial. Put
\[
 b_n=A_nT_n+\frac{\mathcal N_n}{2},\qquad
 B_n=\frac{b_n+4A_n}{8}.
 \tag{2}
\]
The endpoint target is
\[
 Z_{1/2}=\sum_{r\ge1}\frac1{r(r^2-1/4)}=8\log2-4.
\]
The positive Padé error therefore gives \(B_n/A_n<\log2\).

**Arithmetic theorem.** For every \(n\ge1\),
\[
 \boxed{JA_n\in\mathbb Z,\qquad J\ell_JB_n\in\mathbb Z.}
 \tag{3}
\]
The restriction \(n\ge1\) matters: at \(n=0\), \(A_0=1,B_0=1/2\), so the
second assertion with \(\ell_1=1\) would be false.

## 2. The binomial evaluation clears \(A_n\)

For \(0\le k\le n\),
\[
 \frac{J}{2k+1}\binom{n+k}{2k}
 =\binom{n+k}{2k}+2\binom{n+k}{2k+1}\in\mathbb Z.
 \tag{4}
\]
Use the convention that the last binomial is zero when \(k=n\).
Multiplying (4) by the second binomial factor and \(16^k\) proves
\(JA_n\in\mathbb Z\) term by term.

It remains to clear the associated numerator; this is where the Catalan
factor enters.

## 3. A finite Newton quotient reveals the Catalan factor

Write
\[
 w_{n,k}=\frac1{2k+1}\binom{n+k}{2k}^{\!2},
 \qquad A_n=\sum_{k=0}^n16^kw_{n,k}.
\]
The exact divided-product identity
\[
 \frac{D_k(t)-D_k(s)}{t-s}
 =\sum_{j=0}^{k-1}D_j(t)
 \prod_{r=j+1}^{k-1}\bigl(s+(r+1/2)^2\bigr)
 \tag{5}
\]
is a telescoping finite product. At \(s=s_n\), the last factors are
\((r-n)(r+n+1)\), and the Gram identity gives
\(\nu(D_j)=\delta_j/(j+1)\).
Inserting (5) into the definition of \(\mathcal N_n\) and setting \(m=j+1\)
therefore gives
\[
 \mathcal N_n=
 \sum_{k=1}^n16^kw_{n,k}
 \sum_{m=1}^k(-1)^{m-1}
 \frac{\delta_{m-1}(n-m)!}{m(n+m)!}.
 \tag{6}
\]
For clarity, the product tail used here is
\[
 \prod_{r=m}^{k-1}(r-n)(r+n+1)
 =(-1)^{k-m}\frac{(n-m)!}{(n-k)!}
                    \frac{(n+k)!}{(n+m)!}.
\]
Its sign combined with the minus sign in \(A_n-\Phi_n\) gives
\((-1)^{m-1}\) in (6).

Now define
\[
 H_{n,m}=
 \frac{(2m-2)!(2m-1)!(n-m)!}
 {m((m-1)!)^2(n+m)!}.
\]
The factorial identity for \(\delta_{m-1}\) turns (6) into
\[
 \boxed{
 \frac{\mathcal N_n}{16}
 =\sum_{k=1}^n\sum_{m=1}^k
 (-1)^{m-1}16^{k-m}w_{n,k}H_{n,m}.
 }
 \tag{7}
\]

The helpful factorization is
\[
 \boxed{
 H_{n,m}=\operatorname{Cat}_{m-1}
 \int_0^1x^{n-m}(1-x)^{2m-1}\,dx,\qquad
 \operatorname{Cat}_{m-1}=\frac1m\binom{2m-2}{m-1}.
 }
 \tag{8}
\]
The beta integral proves the equality directly. The Catalan factor is an
integer, since
\[
 \operatorname{Cat}_{m-1}
 =\binom{2m-2}{m-1}-\binom{2m-2}{m-2}.
\]
The integrand in (8) is an integer polynomial of degree \(n+m-1\).
Integrating its monomials introduces only denominators at most
\(n+m\le2n<J\). Therefore
\[
 \ell_JH_{n,m}\in\mathbb Z.
 \tag{9}
\]
Equation (4) already gives \(Jw_{n,k}\in\mathbb Z\), and
\(16^{k-m}\) is an integer because \(m\le k\).
Thus every summand in (7) is individually cleared:
\[
 \boxed{J\ell_J\,\frac{\mathcal N_n}{16}\in\mathbb Z.}
 \tag{10}
\]
No cancellation between summands and no experimental divisibility pattern
is needed.

## 4. The harmonic contribution costs the same single lcm

Elementary partial fractions give
\[
 \frac1{r(r^2-1/4)}
 =\frac4{2r-1}+\frac4{2r+1}-\frac4r.
\]
Writing \(H_j=\sum_{r=1}^j1/r\), summation gives
\[
 \frac{T_n+4}{8}
 =H_J-H_n-\frac1{2J}.
 \tag{11}
\]
For \(n\ge1\), \(J\ge3\) is odd, and \(\ell_J\) is divisible by both
\(2\) and \(J\). Hence \(2J\mid\ell_J\), so \(\ell_J\) clears the entire
right side of (11).

From (2),
\[
 B_n=A_n\left(H_J-H_n-\frac1{2J}\right)
                  +\frac{\mathcal N_n}{16}.
\]
Combine \(JA_n\in\mathbb Z\), (10), and (11) to prove (3).
For example,
\[
 A_1=\frac{19}{3},\qquad \mathcal N_1=\frac83,\qquad
 B_1=\frac{79}{18},\qquad \frac{B_1}{A_1}=\frac{79}{114}.
\]
Here \(J\ell_J=3\cdot6=18\) clears both \(A_1\) and \(B_1\).

## 5. A complete irrationality contradiction

The reviewed positive error identity and norm give
\[
 0<A_nZ_{1/2}-b_n\le\frac{2}{J^3A_n}.
\]
Since \(A_nZ_{1/2}-b_n=8(A_n\log2-B_n)\), this yields
\[
 0<J\ell_J(A_n\log2-B_n)
 \le\frac{\ell_J}{4J^2A_n}.
 \tag{12}
\]

The endpoint factor in each summand of the deformed denominator is at
least one. Thus \(A_n\) is at least the ordinary Apéry binomial sum.
The [positive Delannoy generating-function argument](resolvent-legendre-pade-proof.md)
already proves the elementary bound
\[
 A_n\ge
 \frac{\lambda^n}{(n+1)(2n+1)^2},\qquad
 \lambda=17+12\sqrt2.
 \tag{13}
\]
For completeness, the unsquared binomial sum has generating function
\((1-6z+z^2)^{-1/2}\); its positive factorization gives a lower bound
\((3+2\sqrt2)^n/(2n+1)\). Cauchy–Schwarz then gives (13).

A short lcm estimate is sufficient here:
\[
 \ell_M\le4^{M-1}\qquad(M\ge1).
 \tag{14}
\]
To prove it, put \(r=\lceil M/2\rceil\).
Every prime power in \((r,M]\) contributes its prime to
\(\binom{M}{\lfloor M/2\rfloor}\), by the factorial valuation formula.
There is at most one such power of a given prime. Hence
\[
 \ell_M\mid\ell_r\binom{M}{\lfloor M/2\rfloor}.
\]
Also \(\binom{M}{k}\le2^{M-1}\) for \(M\ge2\), by Pascal's identity
and the sum of the previous row. Induction, starting with \(\ell_1=1\),
gives
\[
 \ell_M\le4^{r-1}2^{M-1}\le4^{M-1}.
\]
The last inequality holds for both parities of \(M\).
This proof uses neither a prime asymptotic nor a decimal estimate.

Since \(J=2n+1\), equations (12)–(14) imply
\[
 \boxed{
 0<J\ell_J(A_n\log2-B_n)
 \le\frac{n+1}{4}\left(\frac{16}{\lambda}\right)^n
 \longrightarrow0.
 }
 \tag{15}
\]
Here \(\lambda>17>16\) gives the strict advantage.
If \(\log2=p/q\) were rational, multiplying (15) by \(q\) would give
a positive integer tending to zero: \(J\ell_JA_n\) and \(J\ell_JB_n\)
are integers by (3). This is impossible. Therefore
\[
 \boxed{\log2\ \text{is irrational}.}
\]

## 6. What the connection explains

The Catalan number in (8) is exactly the number of rooted plane full binary
trees with \(m-1\) internal vertices, the \(q=2\) case of the team's
Fuss–Catalan strand. No tree valuation theorem is needed for this proof:
its elementary integrality converts the associated numerator into the
integral of an integer polynomial. The connection is an exact factorial
factorization, not an analogy inferred from numerical values.

At the coalesced parameter \(\tau=0\), the same finite-quotient coefficient
is instead
\[
 \frac{\delta_{m-1}(n-m)!}{m(n+m)!}
 =\frac1{m^3C_{n,m}},
 \qquad \delta_{m-1}=((m-1)!)^2.
\]
That is the coefficient treated by the cubic-lcm arithmetic in the
zeta-three proof. At \(\tau=1/2\), (7)–(9) reveal the single-lcm beta
integral. These are two separately proved arithmetic mechanisms within the
same positive family.

Analytically, the two separated logarithmic arguments coalesce into the
third-power zeta sum. Arithmetically, the endpoint and coalesced member
need their own denominator proofs. Passing to a limit of irrational
members would not establish irrationality of the limit.

## Reproduction and scope

[check_log2_catalan_arithmetic.py](../computations/check_log2_catalan_arithmetic.py)
checks the finite quotient against independently reconstructed endpoint
moments, every termwise clearing claim through \(n=30\), the first
approximants, and the elementary lcm recursion.
These checks supplement the all-index proof.
No arithmetic theorem for arbitrary intermediate deformation parameters,
and no historical priority claim for this logarithmic proof, is asserted.
