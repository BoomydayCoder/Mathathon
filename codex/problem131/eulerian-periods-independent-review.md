# Independent review of the all-m Eulerian period decomposition

**Verdict.** Sections 1–3 of `eulerian-rodrigues-period-decomposition.md` are correct: the exact Eulerian coefficients, gamma reflection, contour normalization, global vertical maximum, conditional component asymptotic, and unconditional central-component theorem all pass independent review. The statements hold for each fixed positive integer \(m\); no uniform estimate as \(m\) grows is asserted. The main-form asymptotic for \(m=2\) has a separate reviewer. Any subsequent finite root certificates for \(m=3,4,5\) are outside this review's scope unless separately identified below.

There is also a proved strengthening: the central component is **strictly positive at every index**, not merely eventually positive. Section 5 below supplies a positive Mellin-convolution representation.

## 1. Eulerian coefficients and all contour signs

Let \(q=2m-1\) and \(w=e^{2\pi iz}\). In the upper half-plane,
\[
 \pi\cot\pi z=-i\pi(1+2\sum_{\ell\ge1}w^\ell),\qquad
 \sin^{2m}\pi z=\frac{(-1)^m(1-w)^{2m}}{2^{2m}w^m}.
\]
Differentiating \(q\) times multiplies each summand by \((2\pi i\ell)^q\). Combining all factors gives exactly
\[
 -\pi^{-2m}\sin^{2m}\pi z\,(\pi\cot\pi z)^{(q)}
 =w^{1-m}A_{2m-1}(w).
\]
The minus sign and the absence of an extra factorial on this identity are correct. The factorial enters only from the original higher-log functional. Eulerian symmetry gives conjugate component pairs, and the total Eulerian sum is \((2m-1)!\), so the displayed normalized weights are positive and sum to one.

The rational squared Mellin transform has zeros of order \(2m\) at \(1,\ldots,n\). Thus its derivative of order \(q\) vanishes at every early integer. A clockwise right closure of an upward cotangent line gives the minus sign in the residue sum. The original factor \(-1/q!\) converts this to \(+1/(2\pi iq!)\int K R^{(q)}\); the \(q\) integrations by parts, with \(q\) odd, then yield the displayed negative integral of \(K^{(q)}R\).

An intermediate boundary product is \(K^{(r)}R^{(q-1-r)}\). At an early integer its pole order is at most \(r+1\), while the second factor has a zero of order at least \(r+2\). It is therefore regular and vanishes. The final product \(K^{(q)}R\) is regular too. Hence moving through an early integer introduces no missing residue or endpoint contribution.

Gamma reflection gives
\[
 R=\kappa^2\sin^{2m}(\pi z)H/\pi^{2m};
\]
the sign from the integer shift disappears in the even power. The Fourier components are holomorphic on \(0<\Re z<n+1\). At imaginary infinity, \(H\) has exponential factor \(e^{-2m\pi|\Im z|}\); the largest permitted Fourier tilt costs only \(e^{2\pi(m-1)|\Im z|}\). A strictly decaying exponential remains for every component. This justifies splitting and shifting the individual integrals, rather than relying on cancellation between nonconvergent pieces.

In particular the exact decomposition is
\[
 \mathscr L_{n,m}=\sum_{j=1-m}^{m-1}
 \frac{E(2m-1,j+m-1)}{(2m-1)!}\,T_{n,m,j}.
\]
Conjugation of the parametrized upward-line integral, followed by \(y\mapsto-y\), gives \(\overline{T_j}=T_{-j}\); it introduces no orientation minus sign.

## 2. Uniform Stirling and the global maximum

Independently collecting the Stirling powers gives
\[
 \kappa^2\sim m(2\pi n)^{1-m}m^{2mn},\qquad
 H(nt)\sim(2\pi)^{2m}n^{-2}f_m(t)e^{n\psi_m(t)}.
\]
Their product is precisely formula (4) of the main note. For a fixed vertical line \(t=a+iy\), \(0<a<1\), each gamma argument has real part at least a positive constant times \(n\), and lies in the right half-plane. The usual sectorial Stirling remainder is uniformly \(O(1/n)\) there, including unbounded \(y\). The same conclusion follows by using \(\Gamma(w+1)\)'s shifted Stirling form, so no nonuniform replacement of \(w+1\) by \(w\) is being made.

Direct rational simplification gives
\[
 \Re\psi_m''(a+iy)=
 \frac{2m(m+1)[(3a+m-1)y^2+a(1-a)(a+m)]}
 {(y^2+a^2)(y^2+(1-a)^2)(y^2+(a+m)^2)}>0.
\]
Consequently \(\Im\psi_m'(a+iy)\) is strictly increasing from \(-2m\pi\) to \(2m\pi\). The derivative with respect to \(y\) of the tilted real action is
\(-\Im\psi_m'(a+iy)-2\pi j\). It has exactly one zero for \(|j|\le m-1\), is positive before that zero and negative after it, and has nonzero limiting slopes. This proves the claimed global maximum and the tail gap. It does **not** by itself say that the maximum solves the full complex saddle equation; the main theorem correctly imposes that additional condition.

A useful further consequence is that \(\psi_m'\) is injective on the strip. For distinct \(t_1,t_2\), convexity of the strip gives
\[
 \frac{\psi_m'(t_1)-\psi_m'(t_2)}{t_1-t_2}
 =\int_0^1\psi_m''(t_2+s(t_1-t_2))ds,
\]
whose real part is positive. Thus a given Fourier component has **at most one** full saddle in this strip. This does not prove its existence.

## 3. The conditional Gaussian theorem has the stated domain

Suppose a fixed \(t_j\) lies strictly in the strip and obeys \(\psi_m'(t_j)+2\pi ij=0\) with the specified logarithms. Choosing the line through it makes the unique real-action maximum a full saddle. Along this line the quadratic term is \(-\psi_m''(t_j)(y-y_j)^2/2\), with positive real part of its coefficient. Therefore the Gaussian integral uses the unique square root having positive real part.

The prefactor is analytic near the saddle, nonzero there, bounded on the line, and \(O(|y|^{-2})\) at infinity. The uniform Stirling error and the global tail gap justify localization. On a symmetric local interval, odd expansion terms integrate to zero against the even complex Gaussian; the remaining error is relative \(O(1/n)\). The Gaussian factor \(\sqrt{2\pi/(n\psi_m'')}\), the contour factor \(n/(2\pi)\), and the Stirling prefactor combine to
\[
 \mathcal C_m(t)=m(2\pi)^{m+1/2}\frac{f_m(t)}{\sqrt{\psi_m''(t)}}.
\]
This checks both the power \(n^{-m-1/2}\) and the constant.

The action identity
\[
 \psi_m(t)-t\psi_m'(t)=2m[\log(1-t)-\log(t+m)]
\]
proves the stated exponential value exactly after applying the full saddle equation. Exponentiating half that equation gives the Laurent critical-point equation with sign \((-1)^{j+m+1}\). This implication is valid; the converse would require the logarithm branch and strip membership, exactly as the note warns.

## 4. The central theorem is unconditional for every fixed m

On the real interval \((0,1)\), the ratio
\(t^{m+1}/[(1-t)^m(t+m)]\) increases strictly from zero to infinity. It has one point \(t_0\) where the ratio is 1. All logarithms are real there, so it solves the full \(j=0\) saddle equation, not merely its exponentiation. Its \(f_m\) and \(\psi_m''\) are positive. Also
\[
 0<\eta_0=\left[\frac{m(1-t_0)}{t_0+m}\right]^{2m}<1.
\]
The conditional theorem therefore yields an unconditional positive-leading-coefficient exponentially decaying central component for every fixed positive integer \(m\).

Two exact normalization checks are:

- For \(m=1\), \(t_0=1/\sqrt2\), \(\eta_0=17-12\sqrt2\), and
  \(\mathcal C_1=2^{3/4}\pi^{3/2}/(1+\sqrt2)^2\). This is the established Beukers/Apéry integral constant for the normalization \(2A_n\zeta(3)-2B_n\).
- For \(m=2\), \(t_0=2/3\), \(\eta_0=1/256\), and
  \(\mathcal C_2=\pi^{5/2}/(2\sqrt2)\). Its Eulerian weight is \(4/6\), matching the separately derived central contribution in the ζ(5) decomposition.

Neither result bounds the other components by the central one. Positive coefficients in a sum of complex components do not give positivity or a central-dominance estimate for the full form.

## 5. Exact positive convolution for the central component

A direct beta-factor identity strengthens the central conclusion:
\[
 \boxed{\kappa_{n,m}^2H_{n,m}(z)
 =B(z,n+1-z)^{2m}B(z,mn+1)^2.}
\]
This follows by canceling the factorials and gamma factors. The two beta factors are Mellin transforms of
\[
 f_n(u)=(1+u)^{-n-1},\quad u>0;\qquad
 g_{mn}(x)=(1-x)^{mn}\mathbf1_{(0,1)}(x).
\]
Choose any \(0<c<n+1\). After taking logarithms of the variables and tilting by \(e^{cv}\), both are positive integrable bounded functions. Their repeated additive convolution is continuous and strictly positive at zero. Fourier inversion, justified also by the exponential decay of the product gamma transform, identifies that value with \(T_{n,m,0}\). Therefore
\[
 \boxed{T_{n,m,0}>0\quad\text{for every }n\ge1,m\ge1.}
\]

An explicit positive integral is obtained by putting \(V=\prod_{i=1}^{2m-1}u_i\):
\[
 T_{n,m,0}=\int_{(0,\infty)^{2m-1}}\int_0^1\int_0^1
 [(1-x)(1-y)]^{mn}
 \left(\frac{Vxy}{1+Vxy}\right)^{n+1}
 \prod_{i=1}^{2m-1}(1+u_i)^{-n-1}
 \frac{dx}{x}\frac{dy}{y}\prod_i\frac{du_i}{u_i}.
\]
The Mellin-convolution argument proves its finiteness as well as its value. This identifies a positive central period; it does not identify the whole pure higher-zeta form with that period.

## 6. Verification and limits of this review

`python codex/computations/check_eulerian_periods_independent.py` independently generates the cotangent derivative through the polynomial recurrence in \(\cot\pi z\), and the Eulerian coefficients through their inclusion–exclusion formula. Six identities pass, along with the general rational curvature identity, the exact \(m=1,2\) central constants, and the beta-convolution factorization. The original `check_eulerian_rodrigues_periods.py` also passes.

The numerical candidate table in the initial main note is diagnostic only. This review does not use it to establish root locations, logarithm branches, dominant components, or a full asymptotic for \(m>2\). In particular a candidate outside the strip cannot be substituted into the conditional theorem. Crossing \(z=n+1\) would encounter actual poles of the individual reflected components, and requires additional residue analysis. Any new finite certificates or higher-component claims require their own stated verification.
