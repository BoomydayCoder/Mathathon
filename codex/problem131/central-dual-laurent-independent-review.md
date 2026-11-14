# Independent review of the positive central period and its dual Laurent model

**Verdict: the new identities are correct.** This review checks `central-period-dual-laurent.md` independently of its own checker: the coordinate change and weight, constant term, unique positive minimum, index-zero contour closure and lower-zeta coefficients all pass. It also gives a simpler proof of the constant-term exponential rate and its leading asymptotic constant. No arithmetic independence or general impossibility theorem is asserted.

## 1. The Jacobian produces exactly the stated weight

For \(x=a/(1+a)\),
\[
 \frac{dx}{x}=\frac{da}{a(1+a)},\qquad 1-x=(1+a)^{-1},
\]
and similarly for \(y,b\). Put \(U=\prod u_i\) and
\(D=(1+a)(1+b)+Uab\). The transformed integrand, relative to logarithmic measure in all variables, is
\[
 (1+a)^{-mn-1}(1+b)^{-mn-1}
 \left(\frac{Uab}{D}\right)^{n+1}
 \prod_i(1+u_i)^{-n-1}.
\]
Thus it equals
\[
 \widetilde\Lambda_m^{-n-1}(1+a)^{m-1}(1+b)^{m-1},
\quad
 \widetilde\Lambda_m=(1+a)^m(1+b)^m\prod_i(1+u_i)\frac{D}{Uab}.
\]
The nonconstant weight for \(m>1\) is necessary. Dropping it would change the integral.

Convergence also holds at \(n=0\), although the earlier convolution theorem was initially stated for \(n\ge1\). In logarithmic variables and with any tilt \(0<c<n+1\), the factors
\[
 e^{cv}(1+e^v)^{-n-1},\qquad
 e^{cv}(1-e^v)^{mn}\mathbf1_{v<0}
\]
are bounded and integrable, including when \(n=0\). Repeated convolution is bounded, continuous and strictly positive at zero. The product of their Mellin transforms decays exponentially on the vertical line, so Fourier inversion applies. This proves finiteness and positivity of the starting integral, hence of its exact coordinate transform.

## 2. Extracting the constant term

Split the last factor as
\(D/(Uab)=1+(1+a)(1+b)/(Uab)\).
Selecting its second summand \(k\) times contributes \(\binom nk\). The \(u_i\)-factor is \((1+u_i)^nu_i^{-k}\), whose constant term is \(\binom nk\), independently for each of the \(2m-1\) variables. Each of \(a,b\) contributes
\((1+a)^{mn+k}a^{-k}\), with constant term \(\binom{mn+k}{k}\). Therefore
\[
 \boxed{b_{n,m}=\operatorname{CT}\widetilde\Lambda_m^n
 =\sum_{k=0}^n\binom nk^{2m}\binom{mn+k}{k}^{2}.}
\]
This is genuinely the transposed sum, rather than the original Rodrigues coefficient. At \(m=1\) both models agree after relabeling variables; at \(m=2\), the first values are \(1,10,626\), while the original coefficient at \(n=1\) is 146.

## 3. Unique positive minimum and the Newton-polytope argument

Let \(Q=Uab/((1+a)(1+b))\). Differentiating \(\log\widetilde\Lambda_m\) in logarithmic coordinates gives
\[
 \frac{u_i}{1+u_i}-\frac1{1+Q}=0,
 \qquad
 \frac{ma}{1+a}-\frac1{(1+Q)(1+a)}=0,
\]
with the analogous equation for \(b\). Thus all \(u_i=s=1/Q\), and
\(a=b=s/(m(1+s))\). Substituting these into the definition of \(Q\), taking only the positive square root, gives
\[
 s^{m+1}=(m+1)s+m.
\]
The polynomial is negative on \([0,1]\), strictly increasing on \((1,\infty)\), and tends to infinity, so there is exactly one positive solution, lying above 1. With \(t=s/(1+s)\), this becomes
\[
 t^{m+1}=(1-t)^m(t+m).
\]
Evaluation of the Laurent polynomial gives
\[
 \widetilde\Lambda_m=\left(\frac{m+t}{m(1-t)}\right)^{2m}=\eta_0^{-1}.
\]

The global argument is valid. A positive Laurent polynomial in logarithmic coordinates is a finite exponential sum; its log-Hessian is the covariance matrix of its exponent vectors under strictly positive weights. The first summand of the split polynomial has exponent 0 and every positive coordinate unit vector, so the support has full affine span. The covariance is positive definite everywhere. Hence the stationary point is the unique global minimum.

Furthermore, its zero gradient expresses 0 as a convex combination with **all** weights positive. A supporting hyperplane through 0 would then contain the entire support, contradicting full dimension. Thus 0 is in the interior of the Newton polytope, and the log-Laurent polynomial is coercive. There is no unproved boundary minimum assumption.

For \(m=2\), the critical point is \(a=b=1/3\), all three \(u_i=2\), and the value is 256. An exact independent expansion gives log-Hessian leading principal minors
\[
 5/8,\quad3/8,\quad5/36,\quad7/162,\quad1/81,
\]
consistent with strict positivity.

## 4. Elementary growth proof and a quantitative refinement

A torus local-limit argument would justify the root rate, but the positive binomial sum is enough. For fixed \(m\), uniform elementary factorial bounds give
\[
 \log\left(\binom nk^{2m}\binom{mn+k}{k}^{2}\right)
 =nF_m(k/n)+O_m(\log(n+1)),
\]
including the endpoints, where \(0\log0=0\) and
\[
 F_m(t)=2m[-t\log t-(1-t)\log(1-t)]
 +2[(m+t)\log(m+t)-m\log m-t\log t].
\]
Its derivative is
\[
 F_m'(t)=2\log\frac{(1-t)^m(m+t)}{t^{m+1}},
\]
and its second derivative is strictly negative. Hence its unique maximum is the same \(t_0\), and the saddle equation simplifies its value to
\[
 F_m(t_0)=2m\log\frac{m+t_0}{m(1-t_0)}=\log\eta_0^{-1}.
\]
An upper bound by \(n+1\) times the largest summand and a lower bound by the summand nearest \(nt_0\) establish
\[
 \boxed{\lim_{n\to\infty}b_{n,m}^{1/n}=\eta_0^{-1}.}
\]

There is also a standard, fully justified real Laplace refinement. Let
\[
 D_m=\frac{2m+2}{t_0}+\frac{2m}{1-t_0}-\frac2{m+t_0}
 =-F_m''(t_0)=\psi_m''(t_0)>0.
\]
Uniform Stirling on a compact neighborhood of \(t_0\) gives a summand prefactor
\[
 (2\pi n)^{-m-1}\frac{m+t}{m\,t^{m+1}(1-t)^m}.
\]
Outside that neighborhood the entropy has a strict gap. Inside it, Taylor expansion and the Gaussian Riemann sum, whose mesh is \(n^{-1/2}\), give
\[
 \boxed{b_{n,m}\sim \mathcal B_m n^{-m-1/2}\eta_0^{-n},\qquad
 \mathcal B_m=(2\pi)^{-m-1/2}
 \frac{m+t_0}{m\,t_0^{m+1}(1-t_0)^m\sqrt{D_m}}.}
\]
Gaussian domination follows from strict negative curvature on a sufficiently small fixed neighborhood; the entropy gap handles all remaining summands. This proves the asymptotic for every fixed \(m\), with no uniformity in \(m\) claimed.

Combining with the previously independently checked central-period constant \(\mathcal C_m\) gives
\[
 \boxed{\mathcal B_m\mathcal C_m
 =\frac1{t_0^{2m+2}D_m},\qquad
 b_{n,m}T_{n,m,0}\sim
 \frac{n^{-2m-1}}{t_0^{2m+2}D_m}.}
\]
For \(m=2\), \(D_2=81/4\),
\(\mathcal B_2=9/(2\pi)^{5/2}\), and the product constant is \(9/16\). This is a useful exact relation between the dual integer period and the positive central period; it makes no denominator or irrationality claim.

## 5. Initial period: convergence, contour orientation and residues

At \(n=0\), the reflection identity gives exactly
\[
 T_{0,m,0}=\frac1{2\pi i}\int_{c-i\infty}^{c+i\infty}
 \frac{\pi^{2m}}{z^2\sin^{2m}(\pi z)}\,dz,\qquad0<c<1.
\]
For the right edge \(R=N+1/2\),
\[
 |\sin\pi(R+iy)|=\cosh\pi y,
\]
so its absolute integral is bounded by a constant depending on \(m\) times \(R^{-2}\). For fixed \(R\), the two horizontal edges vanish exponentially as their heights tend to infinity. First taking those heights to infinity and then \(N\to\infty\) is a rigorous closure procedure. The upward left edge, closed to the right, has clockwise orientation, so the result is minus the sum of residues at positive integers. The pole at zero stays outside.

With \(h=z-k\),
\[
 \frac{\pi^{2m}}{\sin^{2m}(\pi z)}
 =h^{-2m}\sum_{\ell\ge0}c_{m,\ell}\pi^{2\ell}h^{2\ell},
 \qquad
 (k+h)^{-2}=\sum_{r\ge0}(-1)^r(r+1)k^{-r-2}h^r.
\]
For a residue, \(r=2m-2\ell-1\), which is odd. This yields precisely
\[
 -\sum_{\ell=0}^{m-1}(2m-2\ell)c_{m,\ell}\pi^{2\ell}
 k^{-(2m-2\ell+1)}.
\]
All these series over \(k\) converge absolutely, since their smallest exponent is 3. Thus
\[
 T_{0,m,0}=\sum_{\ell=0}^{m-1}(2m-2\ell)c_{m,\ell}\pi^{2\ell}
 \zeta(2m-2\ell+1).
\]
The checked instances are
\[
 T_{0,2,0}=4\zeta(5)+\frac{4\pi^2}{3}\zeta(3),
\]
\[
 T_{0,3,0}=6\zeta(7)+4\pi^2\zeta(5)
 +\frac{16\pi^4}{15}\zeta(3).
\]
These exhibit the lower periods that the noncentral Fourier components must cancel to recover the original initial value \(2m\zeta(2m+1)\). In particular the positive central period is not that original initial value. This does not prove arithmetic independence of the displayed constants or rule out other positive constructions.

## 6. Independent checks

`codex/computations/check_central_dual_independent_review.py` builds the Laurent polynomial as an explicit multivariate coefficient dictionary and multiplies it directly. Twelve constant-term cases pass without using the binomial extraction formula to generate the polynomial. It also verifies the exact \(m=2\) stationary gradient and covariance, independently in all five variables.

A separate exact power-series inversion checks the initial residue coefficients. Four numerical Barnes quadratures at 55-digit precision, for \(m=1,2,3,4\), agree with the zeta expressions to better than 45 relative decimal digits. These quadratures are independent numerical diagnostics; the contour proof above establishes the identities. The author's original checker was also rerun successfully.
