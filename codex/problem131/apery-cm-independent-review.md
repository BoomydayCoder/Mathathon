# Independent review of the Apéry character rule and optional CM interpretation

**Status: the elementary proof in `apery-cm-character.md`, §§1–5, checks out.** It proves, for every prime \(p\ge5\),
\[
 P_p(t)=\Delta(t)^{(1-(\frac{-6}{p}))/2}Q_p(t)^2,
 \quad\Delta(t)=1-34t+t^2,
\]
using the already independently proved square-factor theorem, Lucas congruences, and the exact differential pullback. No complex multiplication or good-reduction theorem is needed for this main result.

The optional geometric interpretation also checks out, conditional on its explicitly stated standard modular, Cartier, and good-reduction endomorphism facts. Particular care is needed with these hypotheses; details are recorded below. No historical-novelty conclusion is made.

## 1. Finite sums, the differential equation, and the bridge

For
\[
 b_n=\sum_k\binom nk^2\binom{2k}k,
\]
the stated telescoping certificate is a rational identity away from its intermediate quotient poles. Its original finite-sum form also holds at those boundary positions. An independent check evaluated all \(0\le k\le n+1\), \(1\le n\le30\), without division by \(f(n,k)\): **525 exact identities**, including the potentially problematic \(k=n,n+1\) cases. For example, at \(k=n+1\),
\(G(n,n+1)=-(n+1)^2f(n+1,n+1)\), while \(G(n,n+2)=0\), giving the correct last summand. The proof does not require a meaning for a divided certificate at a zero of the summand.

The recurrence correctly yields
\[
 x(x-1)(9x-1)B''+(27x^2-20x+1)B'+(9x-3)B=0.
\tag{1}
\]
The auxiliary Lucas argument is valid even when doubling a low digit produces a carry: when \(2s\ge p\), the low digit of \(2k\) is \(2s-p<s\), so \(\binom{2k}k\equiv0\). When there is no carry, the expected product of the two digit factors holds. Combined with the binomial factors \(\binom nk^2\), this proves the asserted sequence Lucas identity.

The bridge is obtained by dividing the exact integer-series identity
\(A(T(x))=(1-x)B(x)^2\) by its Frobenius-substituted version. The factors being divided have constant term one. In characteristic \(p\),
\[
 \frac{A(T(x))}{A(T(x^p))}=P_p(T(x)),\qquad
 \frac{1-x}{1-x^p}=(1-x)^{1-p},
\]
so
\[
 (1-x)^{p-1}P_p(T(x))=B_p(x)^2.
\tag{2}
\]
The series calculation is legitimate as a formal identity near zero; both sides are rational functions, so it is an identity in \(\mathbf F_p(x)\). No analytic continuation is being assumed in characteristic \(p\).

The earlier proof of the algebraic pullback in `apery-modular-elliptic.md`, §2, uses exact rational differential substitutions and uniqueness in characteristic zero. The name of that note does not introduce a modular-theory dependency into this algebraic step.

## 2. The two special values are purely algebraic

The identity
\[
 \operatorname{CT}\bigl((1+X+Y)(1+X^{-1}+Y^{-1})\bigr)^n=b_n
\]
follows by matching the three multinomial counts; summing over two of them uses \(\sum_j\binom kj^2=\binom{2k}k\).

For \(m=(p-1)/2\), the coefficient identity used to eliminate \(Y\) follows term by term:
\[
 \operatorname{CT}_Y(a-bY-c/Y)^{p-1}
 =\sum_{j=0}^m\binom{p-1}{2j}\binom{2j}j
 a^{p-1-2j}(bc)^j
 \equiv(a^2-4bc)^m.
\]
Here \(\binom{p-1}{2j}\equiv1\) and
\(\binom mj(-4)^j\equiv\binom{2j}j\). Multiplying by \(X^{p-1}\) to replace the Laurent discriminant by its quartic therefore gives the stated coefficient formula for \(B_p\).

At \(x=1\), the quartic is \((X^2-1)^2\), and its coefficient calculation gives exactly \(B_p(1)=1\). The degeneration causes no difficulty: this is a polynomial identity, not an application of smooth-curve Hasse theory at a singular curve.

Also
\[
 b_{p-1}\equiv\sum_{k=0}^{(p-1)/2}\binom{2k}k
 =(-3)^{(p-1)/2}=\left(\frac{-3}{p}\right).
\]
The leading coefficient is nonzero for \(p\ge5\), so \(\deg B_p=p-1\). This degree bound is essential later.

## 3. The involution proof has no characteristic-p gap

Let
\[
 M(x)=\frac{x-1/9}{x-1},\quad h(x)=\frac1{1-x}.
\]
An independent rational calculation confirms \(M\circ M=x\), \(T\circ M=T\), and the two normalized differential-operator identities displayed in the main note. Thus \(y\mapsto h\,y\circ M\) preserves (1). All constant denominators used here are products of powers of 2 and 3; after clearing the rational-function denominators, the identity specializes for every \(p\ge5\).

The claim that \(B_p\) solves (1) is justified by Lucas factorization, rather than naively discarding the higher coefficients of a characteristic-zero solution. If \(B=B_pB(x^p)\), all derivatives of the second factor vanish. Since it is invertible, applying the differential operator gives the equation for \(B_p\).

Now
\[
 R=(1-x)^{p-1}B_p(M)
\]
is a polynomial of degree at most \(p-1\). As a rational function it equals the transformed solution times \((1-x)^p\), whose derivative is zero. Hence it solves the same equation.

The uniqueness statement is exactly strong enough and does not overreach: **among polynomials of degree below \(p\)**, the coefficients are determined by the constant term. At indices \(0\le n\le p-2\), the coefficient \((n+1)^2\) in the recurrence is a unit. Arbitrary formal power series would not enjoy this uniqueness, because characteristic-\(p\) constants can contain \(x^p\); the proof correctly avoids that claim.

Therefore \(R=cB_p\). Evaluating the polynomial expression at \(x=1\), only its term of degree \(p-1\) before the Möbius substitution survives. Its value is
\[
 b_{p-1}(8/9)^{p-1}=b_{p-1},
\]
by Fermat's theorem. With \(B_p(1)=1\), this proves the exact functional equation with multiplier \((\frac{-3}{p})\).

## 4. Fixed-point multiplicity and the character

The fixed-point polynomial is \(D_x=9x^2-18x+1\), with discriminant \(288\ne0\) modulo \(p\ge5\). Its values at the singular points \(0,1,1/9\) are \(1,-8,-8/9\), respectively, so either fixed point is ordinary for (1).

If \(B_p\) vanishes there with multiplicity \(r\), then \(1\le r\le p-1\). The lowest coefficient in (1) is its nonzero leading differential coefficient times \(r(r-1)\) and the leading coefficient of the local polynomial. Thus \(r=1\). This is a direct second-order ODE argument; no Hasse invariant or CM property is involved.

At a fixed point \(x_0\),
\[
 (1-x_0)^2=8/9,\qquad M'(x_0)=-1,
\]
so the multiplier \((1-x_0)^{p-1}\) equals \((\frac2p)\), even when \(x_0\notin\mathbf F_p\). The equality is obtained by raising its square, which lies in \(\mathbf F_p\), to \((p-1)/2\); no unjustified assumption that the fixed point is rational is made.

If \(B_p(x_0)\ne0\), evaluating the functional equation gives
\((\frac2p)=(\frac{-3}p)\). If it vanishes, its derivative is nonzero, and differentiating the functional equation gives the extra factor \(M'(x_0)=-1\), hence
\(-(\frac2p)=(\frac{-3}p)\).
These exhaustive cases prove the claimed equivalence with \((\frac{-6}p)\).

The two critical points map to the two roots of \(\Delta\). One can check distinctness directly: if both mapped to a common value, the degree-two equation \(T(x)=t_0\) would have two distinct double roots, impossible. Equivalently, solve the explicit quadratic or use \(\Delta(T)=D_x^2/(1-x)^2\) together with the nonzero derivative of \(D_x\). The bridge (2) and the established square-factor theorem finish the classification.

**There are no exceptions at 5 or 7.** All excluded constants are supported at 2 and 3, the polynomial-solution uniqueness remains valid, and both fixed points are ordinary. Direct checks give \((\frac{-6}5)=(\frac{-6}7)=1\). For \(p=5\), the remainder of \(B_p\) modulo \(D_x\) is \(x+2\), coprime to \(D_x\); for \(p=7\) it is the nonzero constant 2.

## 5. Optional Hasse and modular interpretation: audited details

The quartic coefficient formula has the stated Hasse interpretation at smooth parameters. In the differential \(dX/V=\Phi_x^{(p-1)/2}dX/V^p\), the polynomial numerator has degree \(2p-2\). The only exponent of the form \(jp-1\) in this range is \(p-1\). Thus the Cartier calculation leaves precisely the coefficient claimed. The differential is holomorphic on the smooth projective genus-one quartic, including infinity. The equivalence of Cartier zero and supersingularity is an external standard elliptic-curve theorem, explicitly acknowledged in the note.

The binary-quartic invariant normalization gives the correct \(j\)-invariant. The script's factor 6912 is consistent with \(\operatorname{disc}\Phi=(4I^3-J^2)/27\), hence \(j=256I^3/\operatorname{disc}\Phi\). Its constant denominators introduce no primes beyond 2 and 3.

The finite modular-function certificate is sufficient for a proof:

- The coordinate \(x\) is a degree-one meromorphic function on \(X_0(6)\), as proved from its single simple pole in the companion note.
- Its cusp values in denominator order \(1,2,3,6\) are \(1/9,\infty,1,0\). The rational \(j(x)\) has pole orders at most \(1,2,3,6\) there and no interior poles.
- The function \(j(6\tau)\) is invariant under \(\Gamma_0(6)\): the corresponding full-level matrix is \(\left(\begin{smallmatrix}a&6b\\c/6&d\end{smallmatrix}\right)\). Its cusp orders are also \(1,2,3,6\). At a cusp of denominator \(c\mid6\), the width is \(6/c\), and the scaling yields order \(c\).
- The exact comparison through powers \(q^{-6},\ldots,q^6\) makes the difference vanish to order at least 7 at infinity. Its total remaining pole degree is at most \(1+2+3=6\). A nonzero meromorphic function on a compact curve cannot have that divisor. No Sturm bound or unproved inference from series agreement is being substituted here.
- The code computes the first 13 coefficients of \(q^6j(x(q))\). The eta-series precision 24 is more than sufficient to recover the required coefficients after this normalization.

At \(\tau_0=i/\sqrt6\), the previously established Fricke and coordinate argument selects \(x_-=1-2\sqrt2/3\). The identity \(j(x(\tau))=j(6\tau)\) then identifies the corresponding complex curve with the lattice \(\mathbf Z+\mathbf Z\sqrt{-6}\). Its endomorphism ring is precisely \(\mathbf Z[\sqrt{-6}]\), the maximal order of discriminant \(-24\). Algebraic conjugation transports its endomorphisms to the other curve; it is not necessary to assert that the two complex curves are isomorphic.

## 6. Optional good reduction and CM criterion: audited assumptions

At every place above \(p\ge5\), the values \(x_\pm=(3\pm2\sqrt2)/3\) are integral. Their reductions avoid \(0,1,1/9\), by the fixed-point polynomial calculation. The quartic has unit leading coefficient \(x^2\) and unit discriminant
\(256x^6(x-1)^3(9x-1)\). Its smooth proper double-cover model therefore gives good reduction. The point \((0,x)\) is a smooth integral section because \(2x\) is a unit. With this section the model is an elliptic scheme. Passing to a finite field extension where all endomorphisms are defined preserves good reduction.

The two critical **j-invariants can coincide modulo 13 or 19**: the discriminant of their class polynomial contains these primes. That does not signal bad reduction of either curve, whose quartic discriminants remain units. The argument correctly uses the quartic discriminant, not the class-polynomial discriminant, to establish good reduction.

The specialization theorem must preserve degree, composition and addition of endomorphisms; the standard good-reduction homomorphism does. This ensures that a degree-\(p^h\) CM endomorphism remains a nonconstant isogeny of exactly that degree after reduction.

For inert primes, the argument on ordinary geometric \(p\)-torsion is valid: this group is cyclic of order \(p\), so every endomorphism acts by a scalar in \(\mathbf F_p\). The specialized relation \(\alpha^2=[-6]\) would give a square root of \(-6\), a contradiction.

For split primes, choosing \((\beta)=\mathfrak p^h\) using the finite ideal class group gives an integral \(\beta\), norm \(p^h\), and trace not divisible by \(p\). The trace assertion follows by reduction modulo \(\mathfrak p\): \(\beta\) vanishes and \(\bar\beta\) is a unit. In a supersingular curve an isogeny of \(p\)-power degree has separable degree one; otherwise its geometric kernel would contain a nontrivial \(p\)-power-order point, hence a nontrivial \(p\)-torsion point. Thus both specialized \(\beta\) and \(\bar\beta\) are purely inseparable of degree greater than one and have zero tangent maps. Addition of endomorphisms adds tangent maps, contradicting the nonzero differential of multiplication by their trace.

This checks the direct CM reduction argument. It uses the stated standard specialization and ordinary/supersingular facts, but does not require a classification of quaternionic endomorphism algebras or an unstated separability hypothesis. It is now optional, since the involution proof establishes the character rule independently.

## 7. Reproducibility and remaining scope

The main script `check_apery_cm_character.py` was rerun successfully, including its exact rational identities, modular finite certificate, Hasse identities, finite-field involution examples, and prime-pattern checks. The independent script `check_apery_cm_independent_review.py` separately checks the finite certificate at its boundaries, the rational involution substitution, and the small-prime coprimality cases.

No proof gap was found in the revised main argument. The standard facts in the optional interpretation remain named dependencies, not newly established foundational theorems. Literature novelty has not been audited; the result should be presented as a proved reconstruction and explanatory connection, not a new discovery in the historical sense.
