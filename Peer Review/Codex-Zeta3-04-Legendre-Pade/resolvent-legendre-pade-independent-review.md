# Independent full-proof audit of the resolvent–Legendre Padé route

**Verdict: accepted as a complete proof of the irrationality of
\(\zeta(3)\).** I found no mathematical gap in
[the proof draft](resolvent-legendre-pade-proof.md).
The conclusion uses the explicitly proved denominator and growth
estimates, not just positivity or formal Padé approximation.
The underlying Wilson–Padé proof is classical and is correctly credited
to Prévost (1996). This is not a novelty claim.

The new independent exact checker is
[check_resolvent_pade_independent.py](computations/check_resolvent_pade_independent.py).
It uses only standard-library rational arithmetic. It does not import
the author's checker, symbolic algebra, a Wilson orthogonality theorem,
an Apéry recurrence, or numerical zeta values.

## 1. The positive measure and its Hurwitz resolvent

Under \(t=y^2\) the measure becomes
\[
 d\nu(y^2)=-\pi y(\operatorname{sech}^2\pi y)'\,dy.
\]
It is positive, has all moments, and integration by parts gives mass
\(\pi\int_0^\infty\operatorname{sech}^2(\pi y)\,dy=1\).
Neither the mass nor the moment existence imports an irrationality
result.

For \(u>0\), the stated rectangular contour for
\(e^{iuz}/\cosh^2(\pi z)\) has only the double pole \(z=i/2\).
Writing \(z=i/2+w\) gives
\(\cosh^2(\pi z)=-\pi^2w^2(1+O(w^2))\), so the residue is
\(-iu e^{-u/2}/\pi^2\). Hence the contour integral is
\(2u e^{-u/2}/\pi\). The top edge reverses orientation and multiplies
the bottom integral by \(e^{-u}\); vertical edges vanish. This yields
the stated half-line transform
\[
 H(u)=\frac{u}{2\pi\sinh(u/2)}.
\]
Continuity supplies \(H(0)=1/\pi\). Differentiation under its real
integral is justified by exponential decay.

Integration by parts against the measure gives
\(\int\cos(u\sqrt t)\,d\nu=\pi(H+uH')\), with a plus sign before
\(uH'\). The Laplace representation of \(1/(a^2+y^2)\) has an
absolute double-integral bound \(a^{-2}\) for every fixed \(a>0\).
Thus Fubini is valid even though the cosine is not positive.
After the second integration by parts, both boundary terms vanish:
\(uH(u)\to0\) at zero, and the exponential controls infinity.
The remaining integrand
\[
 \frac{u^2e^{-(a+1/2)u}}{1-e^{-u}}
\]
is positive and integrable. Tonelli and
\(\int_0^\infty u^2e^{-bu}\,du=2/b^3\) prove the exact factor
\(2\zeta(3,a+1/2)\). There is no missing factor of \(\pi\), \(a\),
or \(2\).

## 2. Bernoulli functional and the Hilbert Gram isometry

The expansion of \(u/(2\sinh(u/2))\) and the cosine transform give
\[
 \nu(t^j)=(-1)^j(2j+1)B_{2j}(1/2).
\]
This is a consequence of the Fourier formula; it does not assume an
even-zeta evaluation. For
\(\mathcal L(f)=\mathcal B(f')\), the pullback identity follows first
on monomials and then by linearity. Both identities
\[
 \mathcal L(f(x+1)-f(x))=f''(0),\qquad
 \mathcal L(f(-x-1))=-\mathcal L(f(x))
\]
hold with the convention \(B_1=-1/2\). In particular the derivative
on the right of the first formula is the **second** derivative.

I checked the finite-product certificate in full, including zero
indices. The roots of \(\phi_j\phi_k\) are the integer intervals
\([-j,j-1]\) and \([-k,k-1]\), with overlap counted twice.
At each such integer, both \(E_{jk}(x)\) and \(E_{jk}(x+1)\)
have at least the required multiplicity. At the endpoint transitions
\(0,-1\), the factor \(x^2\) in \(E_{jk}\) supplies the necessary
multiplicity. Evenness supplies the additional root \(-1/2\);
none of the integer roots is \(-1/2\). Their total multiplicity is
the degree \(2j+2k+1\) of the difference. Its leading coefficient
is \(2(j+k+1)\), proving the displayed polynomial identity.

At zero,
\[
 E_{jk}''(0)=2(-1)^{j+k}(j!k!)^2.
\]
Applying \(\mathcal L\) to the difference and using the pullback
therefore gives
\(\nu(Q_jQ_k)=1/(j+k+1)\), with both signs and factorial factors
canceling. In particular \(j=k=0\) agrees with mass one.

The \(Q_k\) have exact degree \(k\) and nonzero leading coefficients,
so the map \(u^k\mapsto Q_k\) maps every lower-degree polynomial space
bijectively to the corresponding lower-degree space. Its Gram identity
is an isometry of inner products; multiplicativity is neither true nor
used. Rodrigues' formula and integration by parts give the shifted
Legendre norm \(1/(2n+1)\), and the isometry transfers its full
orthogonality and norm to \(\Phi_n\). No Wilson theorem is required.

The special evaluation
\[
 Q_k(-(n+1/2)^2)=(-1)^k\binom nk\binom{n+k}k
\]
then gives \(A_n=\sum_k C_{n,k}^2\). This derives the squared binomial
denominator independently of the Apéry recurrence.

## 3. Exact positive error, including \(n=0\)

Because \(\Phi_n(-a^2)=A_n\), both divided differences used in the
proof are polynomials of degree at most \(n-1\). Orthogonality yields
\[
 \int\frac{\Phi_n^2}{t+a^2}\,d\nu
 = A_n\int\frac{\Phi_n}{t+a^2}\,d\nu.
\]
The resolvent is \(2(\zeta(3)-H_n^{(3)})\) when \(a=n+1/2\).
Combining it with the definition of \(N_n\) gives exactly
\[
 A_n\zeta(3)-b_n
 =\frac1{2A_n}\int\frac{\Phi_n(t)^2}{t+(n+1/2)^2}\,d\nu(t).
\]
The integrand is positive off finitely many zeros, so the error is
strictly positive. Bounding the denominator below by \((n+1/2)^2\)
and using the exact norm gives \(2/((2n+1)^3A_n)\).

For \(n=0\), the divided differences are zero, \(A_0=1\), \(N_0=b_0=0\),
and the same identity reduces to the resolvent at \(a=1/2\). There is
no implicit division by \(n\) or use of a negative-degree polynomial.

The off-diagonal associated-polynomial identity also checks. After
pullback, its quotient splits into half the sum with denominators
\(x-m\) and \(x+m+1\). Reflection maps one quotient to the negative
of the other; the antisymmetry of \(\mathcal L\) makes their functional
values equal. This proof covers all nonnegative integers \(m,n\),
including \(m<n\), and does not infer equality from initial values.

## 4. The exact numerator and termwise denominator parity

The product divided-difference identity for \(D_k\) is a standard
telescoping product and is valid as written. Substitution of
\(\nu(D_j)=(j!)^2/(j+1)\), followed by \(m=j+1\), gives the stated
double sum for \(N_n\). I checked the signs and the factors \((k!)^4\)
before simplification; the result is
\[
 N_n=\sum_{k=1}^n\sum_{m=1}^k
       \frac{(-1)^{m-1}C_{n,k}^2}{m^3C_{n,m}}.
\]

Here is a separate check of the valuation argument, including primes
outside the lcm. Put \(v=\lfloor\log_p n\rfloor\), \(e=v_p(m)\),
and write \(c_\ell(s)\) for the displayed factorial-floor summand.
Writing the residues of \(n,s\) modulo \(p^\ell\) shows directly that
\(c_\ell(s)\in\{0,1,2\}\). When \(\ell\le e\), \(s=m\) is a
multiple of \(p^\ell\), so \(c_\ell(m)=0\). When \(\ell>v\),
\[
 c_\ell(s)=\lfloor(n+s)/p^\ell\rfloor
\]
is nondecreasing in \(s\). It follows that
\[
 v_p(C_{n,m})\le v_p(C_{n,k})+2(v-e),
\]
and the claimed rational summand multiplied by \(d_n^3\) has valuation
at least \(v_p(C_{n,k})+v-e\).

For **\(p>n\)**, \(v=e=0\); every relevant level lies in the monotone
range. The same proof thus shows integrality at these primes too.
Restricting the argument only to primes dividing \(d_n\) would have
been insufficient, but the draft does not make that restriction.

For \(p=2\), if \(e<v\) the lower bound is at least one. If \(e=v\),
then \(m=2^v\) since \(m\le n<2^{v+1}\). Consequently \(n+k\ge2^{v+1}\),
and the level \(v+1\) contributes one to \(v_2(C_{n,k})\).
This includes \(n=1\). Therefore each unsigned summand becomes an
**even** integer, proving \(d_n^3N_n/2\in\mathbb Z\), not merely a
bound with an extra denominator two. Adding the harmonic part proves
\(d_n^3b_n\in\mathbb Z\).

## 5. Growth and the final integer contradiction

The positive Delannoy generating series follows formally from
\[
 \sum_{n\ge k}\binom{n+k}{2k}z^n=\frac{z^k}{(1-z)^{2k+1}}
\]
and the central-binomial generating function. In its factored form
both series have positive coefficients. Keeping only the constant
coefficient of the smaller-root factor proves the required **lower**
bound for \(S_n\), not just an upper growth estimate. The central
binomial bound and Cauchy–Schwarz then imply
\[
 A_n\ge\frac{(17+12\sqrt2)^n}{(n+1)(2n+1)^2}.
\]
This makes the positive error less than \(2(17+12\sqrt2)^{-n}\).

In the lcm argument the floor profile is nonnegative everywhere and
equals one on \(1\le y<3\). Thus it majorizes the exact prime-power
indicator in \(\log d_{\lfloor x\rfloor}-\log d_{\lfloor x/3\rfloor}\).
The endpoints \(y=1,3\) have the correct inclusion conventions.
The factorial integral comparison gives error \(O(\log(x+2))\)
at each scale; the geometric iteration gives \(O(\log^2(n+2))\).
No prime asymptotic, Stirling expansion, or Chebyshev estimate is hidden.

The final strict comparison is rational after taking fourth powers:
\[
 e^{12c}=64\cdot3^9
 <(169/5)^4<(17+12\sqrt2)^4.
\]
Hence \(3c<\log(17+12\sqrt2)\), and the sublinear logarithmic error
does not affect the exponential decay.
If \(\zeta(3)=r/q\) with \(q>0\), then
\(qd_n^3(A_n\zeta(3)-b_n)\) is a positive integer tending to zero.
Every part needed for this contradiction has been proved in the note.

## 6. Source check, reproducibility and scope

I visually inspected the saved scans of Prévost's printed pp.230–231.
Equations (57)–(58) identify the same denominator and associated
numerator; pp.230–231 explicitly finish a Wilson–Padé irrationality
proof. The scan's \(a_n,b_n\) names are reversed relative to this note.
The saved PDF hash matches its provenance:
e144782b86aef226b603078991d4f93350b10b8a0ba13f233b41caa372ca962c.
This verifies the source credit without relying on OCR.

The independent checker passes:

- 62 Bernoulli symmetry identities and 16 real-functional pullbacks;
- 81 complete difference and Hilbert-Gram identities, including zero
  indices;
- 81 off-diagonal associated-numerator and orthogonality cancellations;
- 37,820 exact even-integer summands and 6,843 prime cutoff/parity edges;
- 200 exact factorial/lcm divisibilities, the period-six floor profile,
  the Delannoy generating series, and the strict exponential comparison.

The author's checker also passed, including the five Hurwitz-resolvent
and six positive-error quadratures. Its recurrence comparison is only a
diagnostic and is not a proof input.

The result is a complete elementary reconstruction of a classical Padé
route, with one short Fourier contour calculation. It does not establish
irrationality of other individual odd zeta values or prove historical
novelty of the isometry or auxiliary estimates.
