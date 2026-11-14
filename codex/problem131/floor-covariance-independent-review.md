# Independent review of the floor-covariance derivation

Reviewed source: floor-covariance-barnes-design.md and its exact checker.
Verdict: the centered-indicator formula, covariance identity, allocation
expansion, six-\(J\) profile, and superfactorial valuation formula are
correct. The derivation below checks the symbolic identities independently
of the certificate computation. Use \(x\ge0\), as in the allocation
application; \(0<\alpha<1\).

## 1. Indicator and covariance

Write \(u=k+r\), \(0\le r<1\). For \(0<z<1/2\), away from endpoints,
if \(r\le1/2\) the sum
\(\ell_u(z)=\lfloor u-z\rfloor+\lfloor u+z\rfloor+1\) is \(2k+1\)
when \(z<r\), and \(2k\) when \(z>r\). If \(r>1/2\), it is \(2k+1\)
when \(z<1-r\), and \(2k+2\) when \(z>1-r\). Thus, with
\(\rho=\min(r,1-r)\) and the indicated sign \(\epsilon\),

\[
\ell_u(z)-2u=\epsilon(\mathbf1_{z<\rho}-2\rho).
\]

At \(r=0\) or \(r=1/2\) the centered function is zero almost everywhere,
as required. Its integral is zero. Multiplication gives

\[
\int_0^{1/2}(\ell_u-2u)(\ell_v-2v)\,dz
=\epsilon_u\epsilon_v(\min(\rho_u,\rho_v)-2\rho_u\rho_v).
\]

For \(a,b\in[0,1/2]\), assume first \(a\ge b\). Direct subtraction gives

\[
B_2(a-b)-B_2(a+b)=2b-4ab.
\]

The other ordering is symmetric. Replacing either original fractional part
by its reflection \(r\mapsto1-r\) changes the sign of
\(B_2(\{u-v\})-B_2(\{u+v\})\), by \(B_2(1-t)=B_2(t)\), matching
the change of \(\epsilon\). This proves the claimed covariance formula in
every case, with endpoints understood by periodic continuity. The diagonal
variance formula follows by setting \(u=v\). The Brownian-bridge comparison
uses the length-\(1/2\) covariance \(\min(a,b)-ab/(1/2)\); the factor two is
correct, and no probability result is needed.

## 2. Independent expansion of the allocation

Let \(H=1+2\alpha\), \(\lambda=1-\alpha\),
\(T=2Hx-2\xi\), \(q=2x-2\eta\), with \(\xi,\eta\in[0,1/2)\).
Write \(C=\operatorname{Cov}(x,\alpha x)\) and
\(V=\operatorname{Var}(\alpha x)\). Expanding the original integrand first,
before making any Bernoulli substitutions, gives

\[
\begin{aligned}
\int_0^{1/2}(T-3\ell_{\alpha x})
 (T+3\ell_{\alpha x}-\ell_x-5)\,dz
={}&\tfrac12T^2-Tx-\tfrac52T\\
 &+6\alpha(1-3\alpha)x^2+15\alpha x+3C-9V.
\end{aligned}
\]

Adding \(\xi(2T-q-5)+(\xi-\eta)_+\), then substituting \(T,q\),
cancels every \(x\xi\) and linear-\(\xi\) term and yields exactly

\[
\Gamma_\alpha(x)=10\alpha(1-\alpha)x^2-5(1-\alpha)x
 +3C-9V-2\xi^2+2\xi\eta+(\xi-\eta)_+.
\]

For the remaining term, use \(2\xi=\{2Hx\}\), \(2\eta=\{2x\}\), and
\(\{4\alpha x\}=\{2\xi-2\eta\}\). If \(\xi\ge\eta\) the last
fractional part is \(2(\xi-\eta)\); otherwise it is \(1+2(\xi-\eta)\).
Expanding \(B_2\) in these two cases proves the displayed identity for \(E\)
in the reviewed note.

Now substitute the covariance and variance formulas, followed by
\(B_2(\{v\})=v^2-v+1/6-2J(v)\). The coefficient of \(x^2\) is

\[
10\alpha-10\alpha^2-6\alpha+18\alpha^2
 -(1+2\alpha)^2+1-4\alpha^2=0.
\]

The coefficient of \(x\) is \(-5+\alpha\), and the constants cancel.
The six remaining terms are therefore precisely

\[
\boxed{\begin{aligned}
\Gamma_\alpha(x)={}&-(5-\alpha)x
 +\tfrac12J(2(1+2\alpha)x)-\tfrac12J(2x)
 +\tfrac12J(4\alpha x)\\
 &-3J((1-\alpha)x)+3J((1+\alpha)x)-9J(2\alpha x).
\end{aligned}}
\]

This proves the profile for arbitrary \(0<\alpha<1\) and \(x\ge0\),
not just the rational inputs checked by the program. It does not extend
the validity range of the underlying source row allocations.

## 3. Why the profile belongs to superfactorials

For \(B(N)=\prod_{j=0}^{N-1}j!\), count the multiples of each \(p^a\)
in each factorial:

\[
\begin{aligned}
v_p(B(N))
 &=\sum_{a\ge1}\sum_{j=0}^{N-1}\Big\lfloor\frac j{p^a}\Big\rfloor\\
 &=\sum_{a\ge1}\sum_{k\ge1}(N-kp^a)_+
 =\sum_{a\ge1}p^aJ(N/p^a).
\end{aligned}
\]

If \(p^a\mid N\), the final term with \(k=N/p^a\) is zero, so no
endpoint correction is missing. \(N=0\) is also valid by empty products.
Thus the floor profile has exactly the prime-power shape of a finite
superfactorial product. The small-prime repair and majorant certificate
remain separate requirements; this identity alone does not prove an
integral global normalization or its asymptotic size.

## Verification

Independently reran the standard-library checker successfully:
1,334 covariance integrals, 5,995 allocation integrals, and 400
superfactorial valuations. The symbolic arguments above, rather than
these finite tests, establish the all-parameter identities.
