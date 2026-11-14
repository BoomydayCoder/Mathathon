# Independent audit: upper-prime cancellation and the fixed-grid normalizer

**Status: accepted under the stated properness and monic degree-drop
hypotheses.** This reviews
[rank-one-prime-band-cancellation.md](rank-one-prime-band-cancellation.md)
and Section 4 of
[balanced-rank-one-trial-experiment.md](balanced-rank-one-trial-experiment.md).
The refined normalizer is valid even when a band prime divides the
coefficient denominator. The finite numerical error data remain diagnostics;
no asymptotic arithmetic obstruction or irrationality result follows.

## 1. Why the paired pole values agree

Let \(s=2m+1\ge3\), \(K\ge m+1\), \(F\in\mathbb Z[t]\), and
\(\deg F<K\). Write
\[
 D_K(t)=\prod_{j=1}^K(t+j^2),\qquad
 \mu_s(F/D_K)=a_F\zeta(s)+b_F.
\]
Fix a prime \(p>K\). Then \(p\) is odd and \(p>m\), so \(2m\) is a unit
modulo \(p\). All harmonic denominators in
\[
 \alpha_j=j^{2m},\qquad
 \beta_j=-j^{2m}H_j^{(2m+1)}-\frac1{2m}+\frac1{2j}
\]
are also \(p\)-units for \(1\le j\le K\). This remains true when \(p\)
divides \(s\); no exclusion of that case is needed.

For \(K<p\le2K\), the only distinct colliding nodes \(-j^2,-k^2\)
modulo \(p\) have \(j+k=p\). No node has two distinct partners.
Oddness of \(s\) gives \(H_{p-1}^{(s)}\equiv0\pmod p\), simply by
pairing nonzero residues. For \(k=p-j\),
\[
 H_k^{(s)}\equiv H_{j-1}^{(s)}
 =H_j^{(s)}-j^{-s}\pmod p.
\]
The parity of \(2m\) then gives \(\alpha_k\equiv\alpha_j\).
In \(\beta_k\), the added \(j^{-1}\) from the harmonic correction and
the reflected \(-1/(2j)\) combine to \(+1/(2j)\). Thus
\[
 \beta_k\equiv\beta_j\pmod p.
\]
The half-pole correction is indispensable to this equality.

## 2. Why a pair cancels its entire possible denominator

Properness gives a pure partial-fraction expansion, with no polynomial part:
\[
 \frac{F(t)}{D_K(t)}
 =\sum_{j=1}^K\frac{c_j}{t+j^2},\qquad
 c_j=\frac{F(-j^2)}{D_K'(-j^2)}.
\]
The derivative formula
\[
 D_K'(-j^2)=
 \frac{(-1)^{j-1}(K-j)!(K+j)!}{2j^2}
\]
is correct, including its sign. Since \(K+j<2p\), its \(p\)-valuation
is at most one. It is one exactly when \(j\ge p-K\), the colliding range.

For a pair \(j,k=p-j\), put \(z_j=-j^2\). Factoring its one vanishing
difference gives
\[
 D_K'(z_j)=(z_j-z_k)B_j,\qquad
 D_K'(z_k)=-(z_j-z_k)B_k.
\]
Here \(B_j,B_k\) are \(p\)-units and are congruent modulo \(p\).
Also \(F(z_j)\equiv F(z_k)\), and
\[
 z_j-z_k=p(p-2j),\qquad v_p(z_j-z_k)=1.
\]
Thus \(pc_j\equiv-pc_k\pmod p\). Combined with the equal pole-value
residues, this proves that both \(c_j\alpha_j+c_k\alpha_k\) and
\(c_j\beta_j+c_k\beta_k\) are \(p\)-integral.
There is no unaddressed higher power of \(p\): each summand had valuation
at least \(-1\). Noncolliding terms are integral individually.
For \(p>2K\), every residue denominator is already a \(p\)-unit.
This proves \(a_F,b_F\in\mathbb Z_{(p)}\) for every \(p>K\).

The optional stronger statement \(a_F\in\mathbb Z\) is also correct:
\[
 a_F=(-1)^m[t^{K-1}]
 \operatorname{rem}_{D_K}(t^mF(t)).
\]
The remainder has integer coefficients because \(D_K\) is monic.

## 3. The exact fixed divisor and the original clearing integer

Let \(p_{\rm trial}=P/T\) be rational monic of degree \(d\ge0\), with
\(T\) its least positive coefficient denominator. Then \(P\) is primitive
integer, with leading coefficient \(T\). Set
\[
 K=2d+m+1,\quad
 \ell_K=\operatorname{lcm}(1,\ldots,K),\quad
 G=\gcd_{0\le j\le2d}P(-j^2)>0.
\]
The polynomial \(f(j)=P(-j^2)\) has degree \(2d\). Its Newton coefficients
\(\Delta^r f(0)\) have the same gcd as these initial values, because the
two bases are related by inverse integral triangular matrices.
Generalized binomial coefficients are integers at every integer argument,
including negative arguments, so \(G\) is exactly the fixed divisor on
the full integer grid.

The top Newton coefficient is \((-1)^dT(2d)!\), giving
\[
 G^2\mid T^2((2d)!)^2\mid T^2(2K)!.
\]
The second divisibility follows from the integer central binomial
coefficient and \(2K\ge4d\). Hence
\[
 C_{\rm grid}=\frac{T^2m(2K)!\ell_K^{2m+1}}{G^2}
\]
is a positive integer. This includes \(d=0\), where \(P=T=G=1\).

For the clearing property, insert \(P(-j^2)=Gu_j\) into the residue
formula. Each \((K-j)!(K+j)!\) divides \((2K)!\).
The harmonic denominator divides \(\ell_K^{2m+1}\).
The residue's factor \(2j^2\) cancels the two in \(1/(2m)\), leaving
only \(m\), and cancels the full denominator \(2j\) of the other correction.
Thus every rational pole contribution becomes integral after multiplication
by \(C_{\rm grid}\).

The zeta coefficient is exactly \((-1)^m\): the polynomial
\(t^mp_{\rm trial}(t)^2\) is monic of degree \(K-1\).
Therefore the corresponding rational root \(q_{\rm trial}\) needs no
additional affine-slope division, and \(C_{\rm grid}q_{\rm trial}\) is
integral.

## 4. Removing the band, including primes dividing \(T\)

The stronger fixed-divisor statement \(G\mid(2d)!\) uses primitivity, and
is valid. Indeed the Newton expansion shows that \(G\) divides every
coefficient of \((2d)!f(j)\). The coefficient content of \(f\) is one,
so the gcd of those coefficients is exactly \((2d)!\).
Consequently no prime above \(2d\), and therefore none above \(K\), divides
\(G\).

For each prime \(K<p\le2K\), the factors \(m,\ell_K,G\) are \(p\)-units,
and \(v_p((2K)!)=1\). Thus
\[
 v_p(C_{\rm grid})=1+2v_p(T).
\]
Let \(\Pi_K=\prod_{K<p\le2K}p\). Division by this squarefree product leaves
a positive integer
\[
 \boxed{C_{\rm band}=C_{\rm grid}/\Pi_K}.
\]
Apply the proper-numerator theorem to \(P^2\), whose degree is
\(2d<K\), and then divide by \(T^2\). It gives
\[
 v_p(q_{\rm trial})\ge-2v_p(T)\qquad(p>K).
\]
In the removed band, \(C_{\rm band}\) retains exactly \(2v_p(T)\)
factors of \(p\), enough to clear this denominator even when \(p\mid T\).
At every other prime the original clearing valuation is unchanged.
Therefore \(C_{\rm band}q_{\rm trial}\in\mathbb Z\).

The support conclusion is equally precise: denominator primes above \(K\)
can only be inherited from \(T\), with exponent at most \(2v_p(T)\).
It does not prove that the special balanced polynomials have
\(K\)-smooth coefficient denominators at all indices.

## 5. Two necessary scope boundaries

The improper example is correct:
\[
 \frac{t^4}{D_2(t)}
 =t^2-5t+21+\frac{-85t-84}{D_2(t)},\qquad
 \mu_3(t^4/D_2)=\frac{4919}{12}-341\zeta(3).
\]
The proper part contributes \(405-341\zeta(3)\), while the polynomial
quotient contributes \(59/12\), restoring the forbidden denominator prime
\(3>K=2\). The theorem excludes precisely this polynomial contribution.

There is a second boundary even for a proper numerator: \(p\)-integrality
of both affine coefficients does not automatically imply \(p\)-integrality
of their root. For example
\[
 \mu_3\!\left(\frac{t+2}{D_2(t)}\right)
 =3\zeta(3)-\frac72.
\]
Its affine coefficients are \(3\)-integral, but the root is \(7/6\).
The monic degree-drop condition, making the normalized slope \(\pm1\),
is therefore essential to the approximant support statement.

## 6. Reproducibility and finite experimental claims

The original exact checker passed 1,742 proper-monomial/prime tests,
5,214 colliding pairs, the improper counterexample, and all eighteen saved
normalizers. Only its report write was suppressed; the expensive
high-precision experiment generator was not rerun.

A separate [direct-node-product checker](../computations/check_prime_band_independent.py)
tests 39 arbitrary rational monic trials, including coefficient denominators
with primes in the removed band and above \(2K\), and squares of such primes.
It finds 27 cases where denominator primes are actually inherited.
It verifies both normalizers and the support bound without using the
factorial residue formula. It also checks both counterexamples above and
reconstructs the first balanced value \(84127/69984\), whose exact deficit
from the optimized \(577/480\) is \(1/174960\). All checks passed.

The reported errors and primitive-error growth in the eighteen-case
experiment are explicitly high-precision diagnostics, not interval-certified
inequalities or asymptotics. The exact rational variational deficits and
normalizer divisibilities have a different, exact status. The experiment
does not prove eventual failure of the balanced sequence; nor does removing
the primorial band prove that this sufficient normalizer is minimal or
effective for an irrationality contradiction.
