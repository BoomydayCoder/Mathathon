# Pearson reciprocity and compatibility close the deformed Apéry recurrence

Status: full proof developed from the root agent's polynomial Pearson
proposal. The reciprocal kernel calculation and compatibility argument
below prove the associated-numerator recurrence for every index, not only
finite examples. Exact checks are supplied; a fresh independent audit is
requested. The Wilson family and deformed Apéry recurrences are classical;
no historical novelty or general irrationality theorem is asserted.

We use the already reviewed positive deformation for
\(0\le\sigma=\tau^2\le1/4\), its polynomial moments, and its resolvent.
Define
\[
 v(x)=x(x^2-\sigma),\quad u_n=v(n+1),\quad
 \pi_n=2n+1,\quad \lambda_n=n(n+1),
\]
\[
 \Pi_n(x)=\Phi_n(-(x+1/2)^2),\quad
 K_{nm}=\Pi_n(m)=K_{mn},\quad
 Z_\sigma=\sum_{r\ge1}\frac1{r(r^2-\sigma)},\quad
 T_m=\sum_{r=1}^m\frac1{r(r^2-\sigma)}.
\]
The finite-sum symmetry of \(K\) is an input already proved directly.
The integral remainders are
\[
 F_{nm}=\frac12\int\frac{\Phi_n(t)}{t+(m+1/2)^2}\,d\nu_\sigma(t)
       =K_{nm}(Z_\sigma-T_m)-\frac12N_{nm},
\]
\[
 N_{nm}=\nu_\sigma\!\left(
       \frac{K_{nm}-\Phi_n(t)}{t+(m+1/2)^2}\right).
 \tag{1}
\]

**Theorem.** For every pair \(n,m\ge0\),
\[
 \boxed{F_{nm}=F_{mn}.}
 \tag{2}
\]
For every \(n\ge1\),
\[
 \boxed{F_{nn}+F_{n-1,n-1}
       =\frac{2(3n^2-\sigma)}{n^2-\sigma}F_{n,n-1}.}
 \tag{3}
\]
Consequently the diagonal associated numerator satisfies the same
homogeneous recurrence as the denominator for every \(n\ge1\).
The source at \(n=0\) is recorded explicitly below.

## 1. A polynomial Pearson functional

Let \(\mathcal B(x^j)=B_j\), with \(B_1=-1/2\).
For a polynomial \(R\), choose any primitive \(H'=R\) and define
\[
 \mathcal L_\sigma R
 =\frac1{\tau^2}\mathcal B[
        H(x+\tau)+H(x-\tau)-2H(x)].
 \tag{4}
\]
The additive constant in \(H\) cancels, and the expression is a
polynomial in \(\sigma\). At \(\sigma=0\), it is
\(\mathcal B(R')\). The reviewed moment formula gives
\[
 \nu_\sigma(g(t))
 =\mathcal L_\sigma[(x+1/2)g(-(x+1/2)^2)].
 \tag{5}
\]
The Bernoulli translation and reflection identities imply
\[
 \mathcal L_\sigma1=0,\quad\mathcal L_\sigma x=1,\quad
 \mathcal L_\sigma R(-x-1)=-\mathcal L_\sigma R(x),
\]
\[
 \mathcal L_\sigma\Delta f
    =\frac{f(\tau)+f(-\tau)-2f(0)}{\tau^2},
 \qquad \Delta f(x)=f(x+1)-f(x).
 \tag{6}
\]
The limiting identity at \(\tau=0\) is \(f''(0)\).
Since \(v(0)=v(\tau)=v(-\tau)=0\), this proves the useful Pearson rule
\[
 \boxed{\mathcal L_\sigma\Delta[v(x)h(x)]=0}
 \tag{7}
\]
for every polynomial \(h\), including at \(\sigma=0\).

Define the polynomial difference operator
\[
 LR(x)=v(x+1)[R(x+1)-R(x)]
             +v(x)[R(x-1)-R(x)].
 \tag{8}
\]
The dual polynomial eigen-equation is
\[
 L\Pi_n(x)=2(2x+1)\lambda_n\Pi_n(x).
 \tag{9}
\]
Here is a short independent verification. Put
\(\phi_k(x)=\binom xk\binom{x+k}k\) and \(D=L/[2(2x+1)]\).
Reflection makes \(D\) preserve polynomials in \(x(x+1)\).
Root multiplicities and the leading term give
\[
 D\phi_k=k(k+1)\phi_k+(k^2-\sigma)\phi_{k-1}.
 \]
For \(k\ge1\), the lower coefficient is fixed by evaluating at
\(x=k-1\); all smaller nonnegative integer roots remain zeros.
The leading coefficient follows by expanding the degree-\(2k\) term.
For \(k=0\), \(D1=0\).
In the expansion
\[
 \Pi_n(x)=\sum_{k=0}^n a_{nk}\phi_k(x),\qquad
 a_{nk}=\frac{(n+k)!}{(n-k)!\prod_{r=1}^k(r^2-\sigma)},
\]
the coefficient ratio obeys
\[
 ((k+1)^2-\sigma)a_{n,k+1}
      =(\lambda_n-k(k+1))a_{nk}.
 \]
Thus the coefficients give (9). Equivalently, (9) is the known
dual grid relation from the finite symmetric sum.

## 2. A universal second-kind commutator

For an arbitrary polynomial \(R\), set
\[
 \Theta_R(y)=\mathcal L_{\sigma,x}
                    \frac{R(x)-R(y)}{x-y}.
 \tag{10}
\]
Then the following is a polynomial identity:
\[
 \boxed{
 L_y\Theta_R(y)-\Theta_{LR}(y)
 =2[R(y-1)-R(y+1)]+\mathcal L_\sigma R.}
 \tag{11}
\]

We prove it with formal Laurent series at \(y=\infty\). These are
only bookkeeping: each coefficient to which \(\mathcal L_\sigma\)
is applied is a polynomial in \(x\), and the conclusion is the
finite polynomial identity (11). No convergence or extension of
the functional to actual rational integrals is assumed.

Write
\[
 S(y)=\mathcal L_{\sigma,x}\frac1{x-y}.
 \]
Apply (7) coefficientwise to \(v(x)/(x-y)\). Polynomial division gives
\[
 0=2+v(y)[S(y-1)-S(y)],
 \qquad
 v(y+1)[S(y+1)-S(y)]=2.
 \tag{12}
\]
The constant two comes from
\(\mathcal L_\sigma(2x+y+1)=2\).

The same Pearson rule shifts terms of \(LR\) across this formal kernel.
Its adjoint kernel is
\[
 L_x^*\frac1{x-y}
 =\frac{v(x)}{x-y-1}+\frac{v(x+1)}{x-y+1}
       -\frac{v(x)+v(x+1)}{x-y}.
 \]
A direct cubic divided-difference calculation gives the exact identity
\[
 \boxed{L_x^*(x-y)^{-1}=L_y(x-y)^{-1}-1.}
 \tag{13}
\]
Indeed subtract the right-hand difference operator; the three remaining
divided differences of \(v\) reduce to \(-1\). The coefficient
\(-\sigma x\) cancels.

Since
\(\mathcal L[R(x)/(x-y)]=\Theta_R(y)+R(y)S(y)\),
equation (13) gives
\[
 \Theta_{LR}+(LR)S=L_y(\Theta_R+RS)-\mathcal LR.
 \]
Expanding the two differences of \(RS\) and using (12) proves (11).

## 3. Apply the commutator to the associated numerator

The polynomials \(\Pi_n\) obey \(\Pi_n(-x-1)=\Pi_n(x)\), so
\(\mathcal L_\sigma\Pi_n=0\).
Equation (9), together with
\[
 \frac{(2x+1)R(x)-(2y+1)R(y)}{x-y}
       =(2y+1)\frac{R(x)-R(y)}{x-y}+2R(x),
 \]
implies
\[
 \Theta_{L\Pi_n}(y)=2\lambda_n(2y+1)\Theta_{\Pi_n}(y).
 \]
Therefore, writing \(\Theta_n=\Theta_{\Pi_n}\),
\[
 \boxed{
 L_y\Theta_n-2(2y+1)\lambda_n\Theta_n
       =2[\Pi_n(y-1)-\Pi_n(y+1)].}
 \tag{14}
\]

The polynomial \(N_{nm}\) in (1) is exactly \(\Theta_n(m)\).
To see this, use (5) and split its pullback into
\[
 \frac12\mathcal L_x\left[
   \frac{\Pi_n(x)-\Pi_n(m)}{x-m}
  +\frac{\Pi_n(x)-\Pi_n(m)}{x+m+1}\right].
 \]
Reflection makes the two functional values equal. This proves the
identity at every index, not only through initial-value comparison.

Let \(f_m=Z_\sigma-T_m\). For \(m\ge1\),
\[
 v(m+1)(f_{m+1}-f_m)=-1,\qquad
 v(m)(f_{m-1}-f_m)=1.
 \]
Expanding \(L_m(\Pi_n(m)f_m)\), then using (14), gives
\[
 \boxed{
 L_mF_{nm}=2\pi_m\lambda_nF_{nm}-\delta_{m0}.}
 \tag{15}
\]
At \(m=0\) no negative-index value is needed, since \(v(0)=0\).
The first product difference contributes \(-\Pi_n(1)\).
Equation (14) contributes
\(-[\Pi_n(-1)-\Pi_n(1)]\).
Because \(\Pi_n(-1)=\Pi_n(0)=1\), their sum is exactly \(-1\).
This proves the boundary source in (15).

The original polynomial three-term relation, integrated against the
resolvent, likewise gives
\[
 \boxed{
 L_nF_{nm}=2\pi_n\lambda_mF_{nm}-\delta_{n0}.}
 \tag{16}
\]
Its source is \(-\pi_n\nu(\Phi_n)=-\delta_{n0}\).
It can also be derived by using finite-sum symmetry in (9) and
interpolating the resulting polynomial identity in the evaluation point.

At \(m=0\), equation (16) has constant flux:
\(v(n+1)(F_{n+1,0}-F_{n0})=-1\).
Together with \(F_{00}=Z_\sigma\), this proves
\(F_{n0}=Z_\sigma-T_n=F_{0n}\).
Now (15) says that the transposed array satisfies the same \(n\)-grid
recurrence (16), with the same boundary values and source.
Since \(v(n+1)>0\), forward recurrence gives uniqueness.
This proves the full symmetry (2).

## 4. Compatibility forces the adjacent balance

This step uses both grid equations, not symmetry alone.
The difference operators in \(n\) and \(m\) commute.
Apply \(L_m\) to (16) and \(L_n\) to (15), then subtract.
The sources create no extra terms: \(L1=0\) and
\(\lambda_j\delta_{j0}=0\). It follows that
\[
 \pi_n[L_m,\lambda_m]F
       =\pi_m[L_n,\lambda_n]F.
 \]
The commutator with the quadratic grid is explicit:
\[
 [L_n,\lambda_n]F
 =2(n+1)v(n+1)F_{n+1,m}-2nv(n)F_{n-1,m}.
 \]
Hence
\[
 \pi_n[(m+1)v(m+1)F_{n,m+1}-mv(m)F_{n,m-1}]
 =\pi_m[(n+1)v(n+1)F_{n+1,m}-nv(n)F_{n-1,m}].
 \tag{17}
\]

Set \(m=n-1\). Write \(d=F_{nn}\), \(d_-=F_{n-1,n-1}\),
\(e=F_{n,n-1}\), and \(w=v(n)=n(n^2-\sigma)\).
The two grid equations express the two distance-two entries as
\[
 v(n+1)F_{n+1,n-1}
   =\pi_n(3n^2-n+1-\sigma)e-wd_-,
\]
\[
 v(n-1)F_{n,n-2}
   =\pi_{n-1}(3n^2+n+1-\sigma)e-wd-\delta_{n1}.
 \]
For \(n=1\), the second negative-index term is omitted; its coefficient
vanishes. Substitute these identities into (17). The boundary source is
multiplied by \(n-1\), so contributes zero even when \(n=1\).
After cancellation,
\[
 w(d+d_-)
 =[(n+1)(3n^2-n+1-\sigma)
       +(n-1)(3n^2+n+1-\sigma)]e
 =2n(3n^2-\sigma)e.
 \]
This is (3). The same compatibility proof applies to the finite kernel
\(K\), whose two grid equations have no sources. Thus it also gives
the denominator adjacent balance without a Gosper calculation; the
existing finite telescope remains an independent certificate.

## 5. The complete diagonal recurrence and its initial source

Let \(d_n=F_{nn}=A_nZ_\sigma-b_n\), and set
\[
 \alpha_n=\frac{(n+1)((n+1)^2-\sigma)^2}
                   {2(3(n+1)^2-\sigma)},\quad
 \gamma_n=\frac{n(n^2-\sigma)^2}{2(3n^2-\sigma)},
\]
\[
 \beta_n=(2n+1)[3n(n+1)+1-\sigma]-\alpha_n-\gamma_n.
 \]
Use \(\gamma_0=0\), including when \(\sigma=0\).
Symmetry and (3), substituted into (16) on the diagonal, give
\[
 \boxed{\alpha_nd_{n+1}-\beta_nd_n+\gamma_nd_{n-1}
                 =-\delta_{n0}.}
 \tag{18}
\]
The same calculation for \(K\) gives the homogeneous recurrence
for \(A_n\), including at zero. Consequently
\[
 \boxed{
 \alpha_nb_{n+1}-\beta_nb_n+\gamma_nb_{n-1}=0
 \quad(n\ge1),}
 \tag{19}
\]
while the initial source is
\[
 b_0=0,\qquad
 \alpha_0b_1=1,\qquad
 b_1=\frac{2(3-\sigma)}{(1-\sigma)^2}.
 \]
At \(\sigma=0\), multiplying the recurrence by six gives the
ordinary Apéry coefficients and \(b_1=6\).
There is no homogeneous numerator equation at \(n=0\).

All associated numerator identities are rational-function identities
in \(\sigma\): they hold on the positive interval, and hence identically
where their displayed finite denominators are defined. This extension
does not assert existence of a positive measure at other parameters.

The mechanism is conceptual: the cubic Pearson zeros produce a
second-kind commutator; it yields the missing grid equation; compatibility
of the two quadratic-grid equations forces the adjacent coefficient.
No recurrence was fitted to numerical values.

## Reproduction and independent-review boundary

[check_deformed_pearson_symmetry.py](../computations/check_deformed_pearson_symmetry.py)
uses exact rational arithmetic at six parameters, including the coalesced
member, logarithmic endpoint, and two parameters outside the positive range.
It checks 114 universal polynomial commutators (including polynomials
without reflection symmetry), 90 Pearson identities, the formal Laurent
normalization, 48 adjoint-kernel identities, 294 independently constructed
associated numerators, 432 sourced grid equations, 36 adjacent balances,
and 36 diagonal recurrences including their index-zero sources.
All pass. The out-of-range checks concern rational algebra only.

The all-index proof is the Pearson calculation and compatibility
argument above, not extrapolation from those computations. A fresh
independent audit should particularly check the formal-kernel shifts,
the constant \(-1\) in (13), and the source terms in sections 3–5.
