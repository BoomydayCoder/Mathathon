# Strictly nested rational Gaussian–Radau bounds for every odd zeta value

Status: full proof developed from the root agent's strengthening of the
audited higher-order positive-kernel construction. Exact rational checks
are supplied, and the [independent audit](nested-higher-zeta-radau-independent-review.md)
has accepted the full proof. No historical
novelty or irrationality theorem is asserted.

Fix an integer \(m\ge1\), put \(s=2m+1\), and let
\(\mu(dt)=\sigma_m(t)\,dt\) be the positive measure proved in
[higher-order-positive-stieltjes-zeta.md](higher-order-positive-stieltjes-zeta.md).
Its moments \(M_k\) are rational, its mass is \(M_0=1/m\), and
\[
 T_m(z):=\int_0^\infty(z+t)^{-m}\,d\mu(t)
             =2\zeta(2m+1,\sqrt z+1/2),\qquad z>0.
 \tag{1}
\]

Let \(R_n\) be the ordinary \(n\)-node Gaussian resolvent for \(\mu\),
and let \(\widetilde R_n\) be the corresponding Gaussian resolvent
for the positive measure \(t\,d\mu(t)\). Define
\[
 R_n^+(z)=\frac{M_0-\widetilde R_n(z)}z,\qquad
 D_rF(z)=\frac{(-1)^r}{r!}F^{(r)}(z).
 \tag{2}
\]
For \(n=0\), use \(R_0=\widetilde R_0=0\).
With \(a_n=n+1/2\), set
\[
 L_{n,m}=H_n^{(s)}+\frac12D_{m-1}R_n(a_n^2),\qquad
 U_{n,m}=H_n^{(s)}+\frac12D_{m-1}R_n^+(a_n^2).
 \tag{3}
\]

**Theorem.** These are rational numbers, and for every \(n\ge0\),
\[
 \boxed{
 L_{n,m}<L_{n+1,m}<\zeta(2m+1)
       <U_{n+1,m}<U_{n,m}.}
 \tag{4}
\]
Both sequences converge to the target. For fixed \(m\), their enclosure
width is
\[
 O_m\!\left(n^{-1/2}\left(\frac2{\pi e}\right)^{4n}\right).
 \tag{5}
\]
The theorem concerns approximation and nesting. It does not bound the
reduced denominators sufficiently for an irrationality argument.

## 1. The exact associated-polynomial Wronskian

This subsection works for any positive measure on \((0,\infty)\) with
infinite support and all moments finite. Let \(p_n\) be its monic
orthogonal polynomials, \(h_n=\int p_n^2\,d\mu>0\), and
\[
 Q_n(z)=(-1)^np_n(-z),\qquad
 P_n(z)=\int\frac{Q_n(z)-Q_n(-t)}{z+t}\,d\mu(t).
 \tag{6}
\]
Then \(R_n=P_n/Q_n\), \(P_0=0,Q_0=1,P_1=h_0\).

Orthogonality gives the monic recurrence
\[
 p_{n+1}(t)=(t-\alpha_n)p_n(t)-\beta_n p_{n-1}(t),
 \qquad \beta_n=h_n/h_{n-1}\quad(n\ge1).
 \tag{7}
\]
For completeness, expand \(tp_n-p_{n+1}\) in the lower orthogonal
basis. Its inner products with \(p_j\), \(j\le n-2\), vanish because
\(\langle tp_n,p_j\rangle=\langle p_n,tp_j\rangle=0\).
Its \(p_{n-1}\) coefficient is \(h_n/h_{n-1}\), since
\(tp_{n-1}=p_n+\) lower-degree terms. Its remaining coefficient is
\(\alpha_n=\int tp_n^2\,d\mu/h_n\). This proves (7) without
importing a continued-fraction theorem.

Consequently
\[
 Q_{n+1}=(z+\alpha_n)Q_n-\beta_nQ_{n-1}.
 \]
The same recurrence holds for \(P_n\) when \(n\ge1\):
inserting (6) adds the term \(\int Q_n(-t)\,d\mu=0\).
For \(n=0\), that extra term is exactly \(h_0\), explaining
the initial numerator \(P_1=h_0\).

The cross-difference therefore satisfies
\[
 P_{n+1}Q_n-P_nQ_{n+1}
   =\beta_n(P_nQ_{n-1}-P_{n-1}Q_n).
 \]
Its initial value at \(n=0\) is \(h_0\), so telescoping gives
\[
 \boxed{
 R_{n+1}(z)-R_n(z)=\frac{h_n}{Q_n(z)Q_{n+1}(z)}.}
 \tag{8}
\]
This constant is \(h_n\), not the measure's mass or a normalized
orthogonal norm. In particular the formula applies unchanged when
the mass is \(1/m\).

Every \(Q_n(z)\) is a product of \(n\) positive factors \(z+x_i\)
on \(z>0\). Thus the right side of (8) is strictly completely
monotone: every alternating derivative is positive. It follows that
\(D_rR_{n+1}(z)>D_rR_n(z)\) for every integer \(r\ge0\).

## 2. Why moving the evaluation point preserves strict monotonicity

The exact ordinary Gaussian error is
\[
 E_n(z)=S(z)-R_n(z)
  =\int\frac{p_n(t)^2}
   {(z+t)Q_n(z)^2}\,d\mu(t),\qquad
 S(z)=\int(z+t)^{-1}\,d\mu(t).
 \tag{9}
\]
Its integrand is a nonnegative coefficient times a product of
positive reciprocal linear factors. The density is positive on
intervals, so every alternating derivative of \(E_n\) is strictly
positive. Differentiation under the integral follows from local
positive lower bounds on \(z\) and the polynomial moments.

In particular
\[
 \mathcal E_{n,m}(z):=D_{m-1}E_n(z)
      =T_m(z)-D_{m-1}R_n(z)
 \]
is positive and strictly decreasing in \(z\).
Equation (8) also gives
\(\mathcal E_{n+1,m}(z)<\mathcal E_{n,m}(z)\) for fixed \(z>0\).
Since \(a_{n+1}^2>a_n^2\),
\[
 0<\mathcal E_{n+1,m}(a_{n+1}^2)
       <\mathcal E_{n,m}(a_n^2).
 \]
But (1) and (3) identify this error as
\(2(\zeta(s)-L_{n,m})\). Hence the lower bounds increase strictly,
despite both the quadrature rule and its evaluation point changing.

## 3. A positive Gauss–Radau rule produces upper bounds

Put \(\widetilde\mu(dt)=t\,d\mu(t)\), and denote its ordinary
resolvent and Gaussian error by \(\widetilde S,\widetilde E_n\).
The elementary identity
\[
 \widetilde S(z)=M_0-zS(z)
 \]
gives
\[
 \boxed{R_n^+(z)-S(z)=\widetilde E_n(z)/z>0.}
 \tag{10}
\]
The right side is strictly completely monotone, so its alternating
derivative of order \(m-1\) is positive.

There is also an explicit positive quadrature interpretation. If
\(\widetilde x_i,\widetilde\lambda_i\) are the Gaussian nodes and
weights of \(t\,d\mu\), then
\[
 R_n^+(z)=\frac{w_0}{z}
       +\sum_{i=1}^n\frac{\widetilde\lambda_i/\widetilde x_i}
                              {z+\widetilde x_i},
 \quad
 w_0=M_0-\sum_i\frac{\widetilde\lambda_i}{\widetilde x_i}>0.
 \tag{11}
\]
The positivity of \(w_0\) is exact:
\[
 w_0=\widetilde E_n(0)
 =\frac1{\widetilde Q_n(0)^2}
       \int\widetilde p_n(t)^2\,d\mu(t)>0.
 \]
The limit at zero is legitimate because the factor \(t\) in
\(\widetilde\mu\) cancels the resolvent denominator and all remaining
polynomial moments exist. For \(n=0\), use \(w_0=M_0\).
The rule is exact through degree \(2n\): write a test polynomial
as its constant value plus \(t\) times a polynomial of degree
at most \(2n-1\), then use the Gaussian rule for \(t\,d\mu\).
Thus (11) is precisely Gauss–Radau quadrature with prescribed node zero.

Applying (8) to \(t\,d\mu\) gives
\[
 R_n^+(z)-R_{n+1}^+(z)
       =\frac{\widetilde h_n}
              {z\widetilde Q_n(z)\widetilde Q_{n+1}(z)},
 \tag{12}
\]
again strictly completely monotone. Together with (10), this proves
that the upper error decreases strictly in both \(n\) and \(z\).
Its value at \(a_n^2\) equals \(2(U_{n,m}-\zeta(s))\).
This proves the upper half of (4). All rationality statements follow
from rational moments and finite rational polynomial linear algebra.

## 4. Explicit rational radii and convergence

For \(n\ge1\), the lower-error bound already proved is
\[
 0<\zeta(s)-L_{n,m}
 \le\frac12\binom{2n+m-1}{m-1}
           \frac{h_n}{z_n^mQ_n(z_n)^2}
 \le\frac12\binom{2n+m-1}{m-1}
           \frac{M_{2n}}{z_n^{2n+m}},
 \quad z_n=(n+1/2)^2.
 \tag{13}
\]
For the upper error, equation (10) has \(2n+2\) reciprocal linear
factors, including the extra \(1/z\). Its alternating divided
derivative of order \(m-1\) has
\(\binom{2n+m}{m-1}\) terms. Bounding the added factors by \(1/z\)
and the integral's original resolvent by \(1/z\) gives
\[
 \boxed{
 0<U_{n,m}-\zeta(s)
 \le\frac12\binom{2n+m}{m-1}
       \frac{\widetilde h_n}
            {z_n^{m+1}\widetilde Q_n(z_n)^2}
 \le\frac12\binom{2n+m}{m-1}
       \frac{M_{2n+1}}{z_n^{2n+m+1}}.}
 \tag{14}
\]
The last step uses the monic trial \(t^n\) for the measure \(t\,d\mu\)
and its positive Gaussian nodes. Every displayed bound is rational.

The lower crude bound is
\(O_m(n^{-1/2}(2/(\pi e))^{4n})\) by the previously proved
Bernoulli moment formula and Stirling estimate.
The upper moment instead contains \(B_{4n+2}(1/2)\) and
\((4n+2)!\); relative to the lower calculation this adds two powers
of \(n\), exactly canceled by the additional \(z_n\) in the
denominator. The remaining polynomial exponent is again
\((m-1)+m+2+1/2-(2m+2)=-1/2\).
Thus (14) has the same rate, proving convergence and (5).

The finite difference \(U_{n,m}-L_{n,m}\) is itself an exact rational
certificate of enclosure width. One can therefore compute until this
width meets any prescribed rational tolerance without evaluating a
zeta constant numerically. The convergence theorem guarantees termination
for fixed \(m\). This effectiveness concerns real approximation, not
the integer linear forms needed for an irrationality proof.

## Reproduction

[check_nested_higher_zeta_bounds.py](../computations/check_nested_higher_zeta_bounds.py)
checks the exact Wronskian constants for both measures, positive Radau
endpoint weights, rational nesting, and the two rational radii.
Finite checks corroborate the proof and do not replace the all-degree
Wronskian and complete-monotonicity arguments.
