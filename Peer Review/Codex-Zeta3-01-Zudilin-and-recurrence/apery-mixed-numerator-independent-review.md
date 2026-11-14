# Independent review: two mixed telescopes and the Apéry recurrence

**Status: accepted after independent mathematical review and exact checks.**
This reviews [the mixed-numerator proof](apery-mixed-numerator-balance.md),
including its optional Bernoulli-functional identification. No historical
novelty conclusion is implied. The result does not identify a varying-measure
Hankel Schur root with the associated numerator quotient.

## 1. Rational identities and all boundary cases

Write

\[
 P_n(t)=\frac{\prod_{j=1}^n(t-j)}{\prod_{j=0}^n(t+j)},\quad
 S_{n,m}=(-1)^{n+m}P_nP_m,\quad
 Q_{n,m}=(t-m-1)(t+m)S_{n,m}.
\]

The first certificate is valid for every pair of nonnegative integers
\(n,m\), omitting the \(n^3\) term at \(n=0\). Independently multiplying by
\((t-n)(t+n+1)/S_{n,m}\), with \(N=n(n+1), M=m(m+1)\), reduces its right
side to

\[
 -2(2n+1)\bigl[(N+M+1)t^2+(M-N)t-NM\bigr].
\]

Direct expansion of the three left terms gives the same quadratic. Ratios
at zeros are used only to prove a rational identity; no evaluation at a
removable singularity is required. The identity is regular at every positive
integer in its original product form, including \(t=n\).

Since \(S=O(t^{-2})\), its derivative is \(O(t^{-3})\). Also
\(Q'=O(t^{-2})\). Differentiate the finite identity, sum from \(1\) to
\(T\), and then let \(T\to\infty\); this justifies telescoping without any
interchange of a parameter derivative and an infinite series.

For \(n\ge1,m\ge1\), the double zero of \(P_nP_m\) gives \(Q'(1)=0\).
For \(n\ge1,m=0\), the factor \(t-1\) in \(Q\) supplies the missing zero.
For \(n=m=0\), \(Q=(t-1)/t\), so \(Q'(1)=1\). For \(n=0,m\ge1\),

\[
 P_m'(1)=\frac{(-1)^{m-1}}{m(m+1)},\qquad Q_{0,m}'(1)=1.
\]

Thus with \(F_{n,m}=-\tfrac12\sum_{t\ge1}S'_{n,m}(t)\), the homogeneous
mixed recurrence holds only for \(n\ge1\), and the exceptional equation is

\[
 F_{1,m}-[2m(m+1)+1]F_{0,m}=-1.
\]

The sign is essential: summing the difference of derivatives gives
\(-Q'(1)\), and multiplication by \(-1/2\) gives the displayed source.

For the second certificate, let \(r=P_n^2/(t-n)^2\). Its left side is
\(4(2t^2-n^2)r\), while its right side is
\(4[t^4-(t^2-n^2)^2]r/n^2\), the same expression. This holds for
\(n\ge1\). The exceptional smallest case presents no boundary problem:
\(T_1=4(t-1)^2/t^2\), so \(T_1'(1)=0\). For every \(n\ge1\),
\(T_n'(1)=0\) and \(T_n'(t)\to0\). Because
\(S_{n,n-1}=-P_nP_{n-1}\), the resulting balance has sign

\[
 F_{n,n}+F_{n-1,n-1}=6F_{n,n-1}.
\]

## 2. Diagonal closure and normalization

Symmetry gives \(F_{n+1,n}=F_{n,n+1}\). The mixed recurrence at \(m=n\)
therefore closes directly after using the adjacent balance twice. Its middle
coefficient is

\[
 6(2n+1)[3n(n+1)+1]-(n+1)^3-n^3
 =(2n+1)[17n(n+1)+5].
\]

The diagonal initial values follow from the source, not a homogeneous
recurrence incorrectly extended to zero:

\[
 F_{0,0}=\zeta(3),\quad F_{0,1}=\zeta(3)-1,\quad
 F_{1,1}=5\zeta(3)-6.
\]

The discovery guide and Beukers integral use \(2F_{n,n}\); the recurrence is
unchanged, but their initial values are \(2\zeta(3)\) and
\(10\zeta(3)-12\).

The separate denominator recurrence is logically sound. The residue of
\(P_n\) at \(-k\) is \((-1)^{n+k}\phi_k(n)\). Consequently the sum of
double-pole coefficients of \(S_{n,m}\) is
\(\Pi_n(m)=\sum_k\phi_k(n)\phi_k(m)\). Translation preserves the sum of
double-pole coefficients, even when the translated rational function has a
nonzero constant polynomial part. Applying this operation to the second
certificate proves the denominator balance. Applying it to the first also
proves the mixed polynomial recurrence at every integer \(m\ge0\), and
hence at every \(y\), since both sides are polynomials in \(y\).
This provides an alternative to importing the Newton-basis derivation.
Subtracting the independently established denominator recurrence from the
linear-form recurrence proves the numerator recurrence. No comparison of
coefficients in an assumed-irrational constant is made.

## 3. Optional associated-numerator identification

Both Bernoulli identities follow coefficientwise from
\(\mathcal B(e^{zx})=z/(e^z-1)\), with \(B_1=-1/2\). The reflection sign
for \(\mathcal L(f)=\mathcal B(f')\) is negative, as required. The even
polynomial \(f_k\) in the proof satisfies its stated finite-difference
identity, giving \(\mathcal L(x\phi_k)=(-1)^k/(k+1)\), including \(k=0\).
The shifted Legendre integral then gives \(\mathcal L(x\Pi_n)=0\) for
\(n\ge1\), but equals \(1\) for \(n=0\).

Taking a divided difference of \(x(x+1)\Pi_n(x)\) produces exactly the
extra term \((x+y+1)\Pi_n(x)\). Its functional vanishes for \(n\ge1\),
so the associated numerator recurrence and initial values are correct.
The fixed-parameter-zero mixed recurrence gives
\(F_{0,m}=\zeta(3)-H_m^{(3)}\); the source determines the next initial
value. Uniqueness then proves

\[
 F_{n,m}=\Pi_n(m)\zeta(3)-C_n(m)
\]

for all nonnegative \(n,m\). Together with the two separately proved
balances, this proves the proposed mixed-numerator theorem and quotient.
The quotient denominator is strictly positive: at the two adjacent integer
arguments every Newton summand is nondecreasing and its \(k=1\) summand
strictly increases, also when \(n=1\).

## 4. Exact computational checks and limits

The original checker was rerun successfully. It verifies symbolic identities
in indeterminate \(n,m,t\), Bernoulli functionals through \(n=8\), and
independent partial-fraction pole jets for 36 mixed pairs.

A [second checker](https://github.com/BoomydayCoder/Mathathon/blob/03cc75d694220d0fba888c637e573a29799e6d5b/codex/computations/check_apery_mixed_independent.py) formed the unreduced products directly, rather than using the
cleared certificate polynomial: 25 pairs \(0\le n,m\le4\), including both
asymmetric boundaries, and seven adjacent pairs \(1\le n\le7\). Every
rational identity and each lower and upper derivative endpoint passed.
These checks corroborate the general algebra and endpoint proofs above;
finite examples alone would not establish the theorem.

**Conclusion.** The two certificates provide a complete short recurrence
transfer, including its source term. The associated Bernoulli-functional
extension is also proved. No unresolved mathematical gap was found in the
reviewed statements; historical priority and the separate Hankel-root
identification remain outside this conclusion.
