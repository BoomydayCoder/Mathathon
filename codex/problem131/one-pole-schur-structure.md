# The one-pole zeta(3) endpoints are Wilson-polynomial identities

**Status: proved**, using the classical Wilson-polynomial orthogonality formula. The proposed identities are not just a small-matrix pattern. They were checked exactly through matrix size 15, and the proof below works for every size. The result identifies this one-pole construction with explicit harmonic-sum tail bounds; it does not prove irrationality. Historical novelty is not claimed.

The stronger conclusion is conceptual: the Schur-complement minimizers are normalized Wilson polynomials. Their evaluations and squared norms reduce to factorials, so the successive improvements telescope to harmonic sums.

## 1. Exact theorem

Use the positive functional for \(\zeta(3)\), whose polynomial moments are
\[
 M_e=\mu_3(t^e)=\frac{(-1)^e(2e+3)B_{2e+2}}2.
\]
For a positive integer \(m\), form the two size-\(m\) pencils
\[
 G_m(X)_{ij}=\mu_{3,X}\left(\frac{t^{i+j}}{t+1}\right),
 \qquad
 G_m^{(1)}(X)_{ij}=\mu_{3,X}\left(\frac{t^{i+j+1}}{t+1}\right).
\]
Their common positive-definite interval is exactly \((L_m,U_m)\), where
\[
 \boxed{L_m=H_m^{(3)}+\frac{m-1}{2m^3},\qquad
 U_m=H_m^{(3)}+\frac1{2m(m+1)}.}
\tag{1}
\]
In particular,
\[
 L_m<\zeta(3)<U_m,\qquad
 U_m-L_m=\frac1{2m^3(m+1)}.
\tag{2}
\]
The endpoints are strictly increasing and decreasing, respectively, and converge to \(\zeta(3)\).

The same elementary bounds can of course be proved directly by telescoping rational estimates for the zeta tail. What requires explanation here is why these precise expressions emerge from the Bernoulli Hankel matrices and the optimized Schur complements.

## 2. The Schur minimizer is an orthogonal polynomial

Let \(\rho\) be the positive measure representing \(\mu_3\) on \([0,\infty)\), and let \(\nu=t\rho\). Put \(n=m-1\). For either measure \(\omega\in\{\rho,\nu\}\), let \(P_n^\omega\) be the monic orthogonal polynomial of degree \(n\), with squared norm \(h_n^\omega\). Its zeros are positive, so \(P_n^\omega(-1)\ne0\).

The lower Schur error is
\[
 \zeta(3)-L_m=
 \min_{\deg p<n}\int(t+1)
 \left(p(t)-\frac1{t+1}\right)^2d\rho(t).
\tag{3}
\]
The upper error is the same minimization with \(d\rho\) replaced by \(d\nu=t\,d\rho\). Define the residual \(q(t)=1-(t+1)p(t)\). It has degree at most \(n\) and satisfies \(q(-1)=1\). Varying \(p\) through polynomials of degree less than \(n\) shows the minimizer satisfies
\[
 \int q(t)r(t)\,d\omega(t)=0\qquad(\deg r<n).
\]
Consequently
\[
 q_n^\omega(t)=\frac{P_n^\omega(t)}{P_n^\omega(-1)}.
\tag{4}
\]
Orthogonality also gives
\[
 \int\frac{(q_n^\omega)^2}{t+1}\,d\omega
 =\int\frac{q_n^\omega}{t+1}\,d\omega,
\tag{5}
\]
because \((q_n^\omega-1)/(t+1)\) has degree less than \(n\).

Put
\[
 R_n^\omega=\int\frac{1-q_n^\omega(t)}{t+1}\,d\omega(t).
\]
Since \(\mu_{3,X}((t+1)^{-1})=X-1\) and \(M_0=1/4\),
\[
 L_m=1+R_{m-1}^{\rho},\qquad
 U_m=\frac54-R_{m-1}^{\nu}.
\tag{6}
\]
These expressions contain only polynomial moments, because the numerator in \(R_n\) vanishes at \(-1\).

## 3. The useful general increment identity

For any positive measure with all moments and existing monic orthogonal
polynomials \(P_n,P_{n+1}\), provided
\(P_n(-1)P_{n+1}(-1)\ne0\),
\[
 \boxed{R_{n+1}-R_n=
 -\frac{h_n}{P_n(-1)P_{n+1}(-1)}.}
\tag{7}
\]
The right side is positive for a measure supported on \([0,\infty)\).

To verify the constant and sign, define
\[
 Q_n(z)=\int\frac{P_n(z)-P_n(t)}{z-t}\,d\omega(t).
\]
The monic three-term recurrence, or the Christoffel–Darboux identity integrated in \(t\), gives
\[
 Q_{n+1}(z)P_n(z)-Q_n(z)P_{n+1}(z)=h_n.
\]
At \(z=-1\), \(R_n=-Q_n(-1)/P_n(-1)\), which proves (7). This is the small calculation that turns orthogonal-polynomial norms into improvements of rational bounds.

## 4. Why Wilson polynomials appear

Writing \(t=y^2\), the \(\zeta(3)\) density is
\[
 d\rho= w_3(y)\,dy,\qquad
 w_3(y)=y^3\frac{d^2}{dy^2}\frac1{e^{2\pi y}-1}
 =\pi^2y^3\coth(\pi y)\operatorname{csch}^2(\pi y).
\tag{8}
\]
The gamma identities
\[
 |\Gamma(1+iy)|^2=\frac{\pi y}{\sinh\pi y},\quad
 |\Gamma(iy)|^2=\frac\pi{y\sinh\pi y},\quad
 |\Gamma(2iy)|^2=\frac\pi{2y\sinh2\pi y}
\]
therefore imply
\[
 w_3(y)=\frac1{4\pi}
 \left|\frac{\Gamma(1+iy)^3\Gamma(iy)}{\Gamma(2iy)}\right|^2.
\tag{9}
\]
Multiplying by \(t=y^2\) replaces \(\Gamma(iy)\) by \(\Gamma(1+iy)\). Thus the two measures are Wilson weights with parameters
\[
 \rho:\ (a,b,c,d)=(1,1,1,0),\qquad
 \nu:\ (a,b,c,d)=(1,1,1,1),
\]
with half the conventional Wilson normalization \(1/(2\pi)\).

This is a concrete link between the Bose-type factor \((e^{2\pi y}-1)^{-1}\), Bernoulli moments, gamma products, and hypergeometric orthogonal polynomials. It is the reason factorial norm formulas are available.

### Precise classical formulas used

Let \(S=a+b+c+d\). Wilson polynomials are
\[
 W_n(y^2;a,b,c,d)=(a+b)_n(a+c)_n(a+d)_n
 {}_4F_3\!\left(\begin{matrix}-n,n+S-1,a+iy,a-iy\\
 a+b,a+c,a+d\end{matrix};1\right).
\tag{10}
\]
Their leading coefficient in \(t=y^2\) is \((-1)^n(n+S-1)_n\). Their norm for the weight with factor \(1/(2\pi)\) is
\[
 \frac{n!\prod_{u<v}\Gamma(n+u+v)}
 {(2n+S-1)\Gamma(n+S-1)},
 \qquad \{u,v\}\subset\{a,b,c,d\}.
\tag{11}
\]
They are mutually orthogonal for positive parameters.

For \(d=0\), take the limit \(d\downarrow0\). No additional mass appears at the boundary: near \(y=0\), the potentially singular ratio is bounded by a constant times \(y^2/(d^2+y^2)\), and away from zero the gamma-product weight has a common exponential majorant. Fixed polynomial moments therefore converge by dominated convergence. The polynomial coefficients and the norm formula also converge, so (10)–(11) apply to \((1,1,1,0)\).

The formulas were checked against the retrieved Koekoek–Swarttouw report, §1.1, equations (1.1.1)–(1.1.2); its norm expression is algebraically equivalent to (11). See the references below for the original Wilson paper and the available source. The proof here uses this standard orthogonality theorem; it does not pretend that theorem was derived from scratch.

## 5. Factorial evaluations and the telescoping calculation

At \(t=-1\), one of \(1+iy,1-iy\) vanishes, so the terminating hypergeometric series in (10) equals one. After dividing by the leading coefficient, the monic evaluations are
\[
 P_n^\rho(-1)=(-1)^n\frac{(n+1)!^3n!}{(2n+1)!},
\tag{12}
\]
\[
 P_n^\nu(-1)=(-1)^n\frac{(n+1)!^3(n+2)!}{(2n+2)!}.
\tag{13}
\]
Taking half of (11) and dividing by the square of the leading coefficient gives
\[
 h_n^\rho=\frac{(n+1)!^4n!^4}{4(n+1)(2n+1)!^2},
\tag{14}
\]
\[
 h_n^\nu=\frac{n!(n+1)!^6(n+2)!}{2(2n+3)(2n+2)!^2}.
\tag{15}
\]
For a normalization check, \(h_0^\rho=1/4\) and \(h_0^\nu=1/12\), exactly \(M_0\) and \(M_1\).

Substitution into (7), followed by \(n=m-1\), yields
\[
 \boxed{L_{m+1}-L_m=\frac{2m+1}{2m^3(m+1)^3},}
\tag{16}
\]
\[
 \boxed{U_m-U_{m+1}=\frac1{m(m+1)^3(m+2)}.}
\tag{17}
\]
Finally, the proposed formulas (1) have precisely these increments:
\[
 \frac1{(m+1)^3}+\frac{m}{2(m+1)^3}-\frac{m-1}{2m^3}
 =\frac{2m+1}{2m^3(m+1)^3},
\]
\[
 \frac1{2m(m+1)}-\frac1{(m+1)^3}
 -\frac1{2(m+1)(m+2)}
 =\frac1{m(m+1)^3(m+2)}.
\]
Both start at \(L_1=1,U_1=5/4\). Induction proves (1) for all \(m\).

## 6. The minimizers themselves are explicit

For \(n=m-1\), let \(D_k(t)=\prod_{j=1}^k(t+j^2)\), including \(D_0=1\). Normalizing (10) at \(-1\) gives
\[
 q_n^\rho(t)=\sum_{k=0}^n
 \frac{(-n)_k(n+2)_k}{(k+1)!^2k!^2}D_k(t),
\tag{18}
\]
\[
 q_n^\nu(t)=\sum_{k=0}^n
 \frac{(-n)_k(n+3)_k}{(k+1)!^3k!}D_k(t).
\tag{19}
\]
The corresponding best approximation is
\(p_n=(1-q_n)/(t+1)\). The divisibility is exact because every \(D_k\), \(k\ge1\), contains \(t+1\). These finite formulas make the Schur optimization reproducible without inverting a matrix.

For example, \(q_1^\rho=(1-3t)/4\), so the optimal constant \(p_1=3/4\), and the lower endpoint is \(1+3M_0/4=19/16\). The denominator \((k+1)!^2k!^2\) in (18) is explained by Wilson's hypergeometric parameters rather than guessed from numerical matrix inverses.

## 7. What this says about zeta(5), and what it does not

The exact polynomial moments satisfy
\[
 \mu_5(P)=\frac1{12}\mu_3\bigl((2\Theta+4)(2\Theta+5)P\bigr),
 \qquad\Theta=t\frac d{dt}.
\tag{20}
\]
This follows on each monomial from
\(M_e^{(5)}/M_e^{(3)}=(2e+4)(2e+5)/12\).
It is a **differential transformation of the functional**, not multiplication of its measure by a polynomial. In particular, it does not automatically preserve the Wilson family.

The densities make the obstacle explicit:
\[
 w_5(y)=\frac{\pi^2}{3}y^2
 \bigl(1+3\operatorname{csch}^2(\pi y)\bigr)w_3(y).
\tag{21}
\]
This identity follows by differentiating the Bose factor four times:
\(f''''=4\pi^4\coth(\pi y)\operatorname{csch}^2(\pi y)
 (1+3\operatorname{csch}^2(\pi y))\).
The extra hyperbolic factor prevents interpreting \(w_5\) as a mere fixed polynomial Christoffel transform of the same Wilson weight. A more general orthogonal-polynomial or differential-functional construction may still be possible; no closed endpoint formula for \(\zeta(5)\) is established here.

For \(\zeta(3)\), (2) shows that this one-pole sequence captures a polynomial error of order at most \(m^{-4}\). It does not reproduce Apéry's exponential mechanism. The many-pole local-arithmetic construction addresses a substantially stronger issue. The present identity is useful as a transparent control example: a large-looking determinant can collapse to familiar optimized tail bounds.

## 8. References and reproducibility

- James A. Wilson, *Some Hypergeometric Orthogonal Polynomials*, SIAM Journal on Mathematical Analysis **11** (1980), 690–701, [doi:10.1137/0511064](https://doi.org/10.1137/0511064). Metadata independently checked through Crossref. The publisher PDF returned HTTP 403 in this session, so no claim of having inspected its full text is made.
- Roelof Koekoek and René F. Swarttouw, *The Askey-scheme of hypergeometric orthogonal polynomials and its q-analogue*, report dated February 20, 1996, [arXiv:math/9602214](https://arxiv.org/abs/math/9602214), §1.1, equations (1.1.1)–(1.1.2). Retrieved and read directly; local files are `../sources/koekoek-swarttouw-9602214.pdf` and `.txt`.
- `localizing-zeta-enclosures.md` proves the Schur-complement formulas and general fixed-denominator convergence theorem used as context here.

Run `python codex/computations/check_one_pole_schur_structure.py`. It checks the Schur endpoints exactly through \(m=15\), the explicit Wilson residuals and norms through degree 6, and the exact increment identities. These computations supplement the all-degree proof above.

## 9. A rigorous obstruction using the actual reduced denominators

The explicit formulas permit a stronger conclusion than observing slow convergence or a large sufficient clearing factor. For this family, the **least denominators themselves** grow exponentially, while the errors are only polynomially small.

Write the endpoints in lowest terms and let
\[
 d_m^- = \operatorname{den}(L_m),\qquad
 d_m^+ = \operatorname{den}(U_m).
\]
Here \(\operatorname{den}\) is the positive denominator after all cancellation.

**Theorem.** For every prime \(p\ge3\) with \(m/2<p<m\),
\[
 \boxed{v_p(L_m)=v_p(U_m)=-3.}
\tag{22}
\]
Consequently,
\[
 \log d_m^\pm\ge\left(\frac32+o(1)\right)m.
\tag{23}
\]
Moreover,
\[
 \zeta(3)-L_m\sim\frac1{4m^4},\qquad
 U_m-\zeta(3)\sim\frac1{4m^4},
\tag{24}
\]
and therefore
\[
 \boxed{d_m^-(\zeta(3)-L_m)\longrightarrow+\infty,
 \qquad d_m^+(U_m-\zeta(3))\longrightarrow+\infty.}
\tag{25}
\]

**Proof of the exact prime valuations.** In \(H_m^{(3)}\), the unique denominator divisible by \(p\) is the term \(1/p^3\). Every other summand is \(p\)-integral. Thus
\(v_p(H_m^{(3)})=-3\).

For \(L_m\), the correction \((m-1)/(2m^3)\) is \(p\)-integral: \(p\ge3\), and \(p\mid m\) would imply \(m\ge2p\) or \(m=p\), both excluded. Adding an integral number cannot cancel the valuation \(-3\).

For \(U_m\), the correction is \(1/[2m(m+1)]\). If \(p\nmid m+1\), it is integral. The only other possibility is \(m+1=2p\); in that case its valuation is \(-1\), still strictly greater than \(-3\). Thus the valuation of the sum is again exactly \(-3\). In particular **no exceptional prime at \((m+1)/2\) has to be discarded**.

It follows that each least denominator is divisible by
\[
 \prod_{m/2<p<m,\ p\ge3}p^3.
\]
Writing \(\vartheta(x)=\sum_{p\le x}\log p\), the prime number theorem gives
\[
 \log d_m^\pm\ge3(\vartheta(m)-\vartheta(m/2))-O(\log m)
 =\left(\frac32+o(1)\right)m.
\]
The endpoint omission \(p=m\), when relevant, contributes at most \(O(\log m)\); it does not affect the exponent.

**Proof of the error asymptotics without Euler–Maclaurin.** The exact increments (16)–(17), together with the explicit limits in (1), give
\[
 \zeta(3)-L_m
 =\sum_{k=m}^{\infty}\frac{2k+1}{2k^3(k+1)^3},
\tag{26}
\]
\[
 U_m-\zeta(3)
 =\sum_{k=m}^{\infty}\frac1{k(k+1)^3(k+2)}.
\tag{27}
\]
Each summand is asymptotic to \(k^{-5}\). More explicitly,
\[
 \frac1{(k+1)^5}\le\frac{2k+1}{2k^3(k+1)^3}\le\frac1{k^5},
\]
\[
 \frac1{(k+2)^5}\le\frac1{k(k+1)^3(k+2)}\le\frac1{k^5}.
\]
Comparison with the integral of \(x^{-5}\) proves (24). It also gives fully explicit lower estimates
\[
 \zeta(3)-L_m\ge\frac1{4(m+1)^4},\qquad
 U_m-\zeta(3)\ge\frac1{4(m+2)^4}.
\tag{28}
\]
Combining these estimates with (23) proves divergence (25).

**Exact scope of the obstruction.** Every integer multiplier that clears the rational endpoint is a multiple of its least denominator. Thus no choice of a smaller, more cleverly canceled integer normalizer can turn either of these endpoint errors into a sequence tending to zero. Under a hypothetical rationality assumption \(q\zeta(3)\in\mathbb Z\), multiplying by the fixed integer \(q\) only increases the displayed positive quantities. Hence the direct nonzero-integer-below-one argument fails for both optimal one-pole endpoint sequences, and for every subsequence tending to infinity.

This is an obstruction to this specific family of approximants, not a proof that every possible fixed-denominator construction fails. It explains why simply increasing the matrix dimension here cannot replace Apéry's simultaneous analytic and arithmetic design. Moving poles and changing the rational-function family are substantive changes to the mechanism, not cosmetic improvements to this sequence.
