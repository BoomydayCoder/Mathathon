# Independent review: nested Gaussian–Radau bounds for odd zeta values

**Verdict:** accepted, with no correction required. This review independently
checks the proof in [nested-higher-zeta-gauss-radau.md](nested-higher-zeta-gauss-radau.md),
including its normalizations, endpoint, strictness, derivative counts, and
fixed-order convergence rate. It uses the already reviewed positive measure
in [higher-order-positive-stieltjes-zeta.md](higher-order-positive-stieltjes-zeta.md).
It is a mathematical review with exact finite corroboration, not formal
verification or a claim of historical novelty. There is no new denominator
bound or irrationality conclusion here.

## 1. Hypotheses and normalization

Fix an integer \(m\ge1\), and write \(s=2m+1\). The supplied measure
\(\mu(dt)=\sigma_m(t)\,dt\) has positive density on \((0,\infty)\),
all moments finite, mass \(M_0=1/m\), and
\[
 \int_0^\infty\frac{d\mu(t)}{(z+t)^m}
 =2\zeta(2m+1,\sqrt z+1/2).
\]
Its exact moment formula is
\[
 M_k=\frac{(-1)^k B_{2k}(1/2)}{m(2m-1)!!}
       \prod_{j=0}^{m-1}(2k+2j+1).
\]
In particular all Gaussian and Radau rational functions have rational
coefficients. The usual monic orthogonal polynomials exist in every degree:
positive density makes every nonzero polynomial have strictly positive
squared norm. Their zeros are positive and simple. The same facts hold
for the measure \(t\,d\mu(t)\).

Let \(D_r=(-1)^r\partial_z^r/r!\). Then
\(D_{m-1}(z+t)^{-1}=(z+t)^{-m}\), with no missing factorial. At
\(z_n=(n+1/2)^2\), the target integral is
\(2(\zeta(s)-H_n^{(s)})\). These two facts check the normalization
of both proposed rational bounds.

## 2. The Wronskian constant really is the monic norm

For a positive measure with the stated moments, let \(p_n\) be monic,
\(h_n=\int p_n^2\,d\mu\), \(Q_n(z)=(-1)^np_n(-z)\), and
\[
 P_n(z)=\int\frac{Q_n(z)-Q_n(-t)}{z+t}\,d\mu(t).
\]
The three-term recurrence has coefficient \(\beta_n=h_n/h_{n-1}\).
Applying its difference quotient gives the same recurrence for \(P_n\)
when \(n\ge1\); the additional term is
\(+\int Q_n(-t)\,d\mu(t)=0\). At \(n=0\), it is instead \(h_0\).
Thus \(P_0=0\), \(P_1=h_0\), and the cross-difference obeys
\[
 W_n:=P_{n+1}Q_n-P_nQ_{n+1}=\beta_nW_{n-1},
 \qquad W_0=h_0.
\]
Consequently \(W_n=h_n\) in every degree, and
\[
 R_{n+1}-R_n=\frac{h_n}{Q_nQ_{n+1}}.
\]
This verifies the sign, the mass \(1/m\), and the potentially delicate
initial case simultaneously. Replacing \(h_n\) by a normalized norm or
by the mass would be incorrect; the manuscript does neither.

## 3. Complete monotonicity and the moving point

The exact Gaussian error is
\[
 E_n(z)=\int\frac{p_n(t)^2}{(z+t)Q_n(z)^2}\,d\mu(t).
\]
Every integrand is a nonnegative coefficient times \(2n+1\) reciprocal
positive linear factors. Every alternating derivative is nonnegative,
and its integral is strictly positive: only finitely many zeros of
\(p_n\) are excluded from a measure with positive density. Differentiation
under the integral is justified locally on \(z\ge\varepsilon>0\):
all differentiated factors are bounded by constants depending on
\(\varepsilon,n,r\), while \(\int p_n^2\,d\mu<\infty\).

The same argument applies to the Wronskian quotient. Hence
\(D_{m-1}E_n\) decreases strictly with \(n\). Its derivative with
respect to \(z\) is \(-mD_mE_n<0\), so it also decreases strictly
with \(z\). Applying these two inequalities consecutively at
\(z_{n+1}>z_n\) proves the strict increase of the moving lower bound.
It does not rely on an unjustified comparison of rational functions
evaluated at different points.

## 4. The zero-node rule is genuinely positive Radau quadrature

Write tildes for the measure \(t\,d\mu\). Since
\(\widetilde S=M_0-zS\),
\[
 R_n^+-S=\frac{\widetilde E_n}{z},\qquad
 R_n^+=\frac{M_0-\widetilde R_n}{z}.
\]
If its Gaussian nodes and weights are \(x_i,\lambda_i>0\), then
\[
 R_n^+(z)=\frac{w_0}{z}+\sum_i\frac{\lambda_i/x_i}{z+x_i},
 \qquad
 w_0=M_0-\sum_i\frac{\lambda_i}{x_i}
     =\frac{\int\widetilde p_n(t)^2\,d\mu(t)}
             {\widetilde Q_n(0)^2}>0.
\]
For the last equality, take \(z\downarrow0\) in the exact error.
The factor \(t/(z+t)\le1\) and the lower bound
\(\widetilde Q_n(z)\ge\widetilde Q_n(0)>0\) provide an integrable
majorant \(\widetilde p_n^2/\widetilde Q_n(0)^2\). The original
measure has no atom at zero, so the pointwise limit is the stated
integral. This also proves strict positivity, not merely nonnegativity.

The rule integrates every polynomial through degree \(2n\) exactly:
write it as \(f(0)+t g(t)\), with \(\deg g\le2n-1\), and use
Gaussian exactness for \(t\,d\mu\). Thus node zero and the positive
weights are justified independently of the resolvent notation.

Moreover
\[
 R_n^+-R_{n+1}^+
 =\frac{\widetilde h_n}{z\widetilde Q_n\widetilde Q_{n+1}}.
\]
This quotient and \(\widetilde E_n/z\) are strictly completely
monotone. The same two-step argument proves strict decrease of the
moving upper bounds. The conventions at \(n=0\) give
\(L_{0,m}=0\) and \(U_{0,m}=4^m/(2m)\); all strict nesting
arguments remain valid there.

## 5. Both rational radii and their common rate

For a product of \(N\) reciprocal linear factors, its divided
alternating derivative of order \(r\) has one term for each weak
composition of \(r\) into \(N\) parts. There are
\(\binom{N+r-1}{r}\) terms, with coefficient one. This explains,
and checks separately, the two combinatorial factors:
\[
\begin{aligned}
0<\zeta(s)-L_{n,m}
&\le\frac12\binom{2n+m-1}{m-1}
 \frac{h_n}{z_n^mQ_n(z_n)^2}
\le\frac12\binom{2n+m-1}{m-1}
 \frac{M_{2n}}{z_n^{2n+m}},\\
0<U_{n,m}-\zeta(s)
&\le\frac12\binom{2n+m}{m-1}
 \frac{\widetilde h_n}{z_n^{m+1}\widetilde Q_n(z_n)^2}
\le\frac12\binom{2n+m}{m-1}
 \frac{M_{2n+1}}{z_n^{2n+m+1}}.
\end{aligned}
\]
The upper error has one more reciprocal factor, namely \(1/z\).
The final inequalities use \(Q_n(z)\ge z^n\) and monic norm
minimality with trial \(t^n\): \(h_n\le M_{2n}\) and
\(\widetilde h_n\le M_{2n+1}\).

For fixed \(m\), the moment formula gives
\(M_{2n}=O_m(n^m(4n)!/(2\pi)^{4n})\). The lower radius has
polynomial exponent \((m-1)+m+1/2-2m=-1/2\) after Stirling's
estimate. The upper one uses \((4n+2)!\), adding two powers of
\(n\), and its extra \(z_n\) removes precisely those powers.
Its exponent is likewise
\((m-1)+m+2+1/2-(2m+2)=-1/2\). Both radii, and hence their
enclosure width, are
\[
 O_m\!\left(n^{-1/2}\left(\frac2{\pi e}\right)^{4n}\right).
\]
This is a fixed-\(m\) statement, not a uniform estimate for growing
zeta order. The explicit radii prove convergence without any separate
moment-determinacy theorem. Rational width supplies a terminating
approximation algorithm for each fixed \(m\), but supplies no bound
on the denominator after reduction.

## 6. Independent exact corroboration

[check_nested_radau_independent.py](../computations/check_nested_radau_independent.py)
uses only Python's standard library. It reconstructs Bernoulli numbers,
moments, monic polynomials, associated numerators, and norms independently.
Higher derivatives are obtained by exact Taylor-series division, rather
than the author's rational-function differentiation recurrence. It checks:

- 56 Wronskians for both measures, \(m=1,\ldots,4\), \(n=0,\ldots,6\);
- 28 rational pairs, including degree zero, strict nesting, endpoint
  weight identities, and both refined and crude radii;
- four actual enclosures and both individual error bounds against
  independent rational zeta intervals obtained from a 50-term positive
  series and integral bounds for its tail.

All checks pass. The author's checker was also rerun successfully:
80 Wronskians, 45 pairs, and five rational series enclosures. These
finite computations corroborate the proof; the all-degree conclusions
come from the preceding norm, factorization, and domination arguments.
