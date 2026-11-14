# Deriving the Apéry operator from the Laurent geometry, including its boundary

Status: explicit rational certificates and analytic boundary argument have
[passed independent review](geometric-divergence-independent-review.md).
The operator coefficients were **unknowns** in a finite exact linear system;
they were not inserted from the known recurrence. The full derivation and
certificate data are reproducible. This establishes the previously missing
geometric exact-form reduction, including the real-chain boundary.
It is not a historical novelty claim. The certificate is sizeable; this is an
explanation of the mechanism plus a computational algebra appendix, not yet a
short hand derivation of all coefficients.

## 1. The conceptual steps

1. The Laurent model is quadratic in one variable. Eliminating that variable
   turns its torus period into a double-cover integral with a square root.
2. Differentiate this algebraic integrand three times with respect to the
   parameter, and seek a linear combination that is a divergence in the other
   two variables. A finite rational coefficient calculation determines the
   Apéry operator, including 34, 51, 27, and 5.
3. A simple identity for a quadratic denominator lifts the double-cover
   divergence to an actual rational three-variable divergence. This matters:
   the real chain gives the logarithmic Legendre Q solution, not merely the
   algebraic P branch.
4. Closed torus cycles kill the divergence. On the positive noncompact chain,
   a boundary survives. A second exact rational reduction makes its integral
   elementary: two infinite edges contribute 6 each. This produces the
   inhomogeneous term, rather than inferring it from initial data.

All variables and normalizations agree with
`apery-laurent-reciprocal.md` and `laurent-legendre-relative-periods.md`.

## 2. Quadratic elimination and an operator with unknown coefficients

Write
\[
 s=xy,\qquad p=(1+x)(1+y),\qquad r=p-s=1+x+y,
\]
\[
 \Lambda=pz+p+p^2/s+p^2/(sz).
\]
For the torus generating function the rational differential form is
\[
 \frac{dx\,dy\,dz}{f},\qquad
 f=az^2+bz+c,
\]
\[
 a=-tsp,\qquad b=s-tp(s+p),\qquad c=-tp^2.
\]
Its quadratic discriminant is
\[
 D=b^2-4ac=s^2-2tsp(s+p)+t^2p^2r^2.                 \tag{1}
\]
Taking the z residue gives, with the branch fixed at t=0,
\[
 \operatorname{CT}_z\frac1{1-t\Lambda}=\frac{s}{\sqrt D}.
\]
Consequently the remaining differential form is \(dx\,dy/\sqrt D\).

Let \(\theta=t\partial_t\). Define explicit polynomials
\[
 N_0=1,\qquad
 N_{k+1}=t\left(D\partial_tN_k-(k+\tfrac12)N_k\partial_tD\right).
                                                               \tag{2}
\]
Then \(\theta^kD^{-1/2}=N_kD^{-k-1/2}\).
Seek an operator
\[
 \mathcal M=\theta^3+
 \sum_{h=1}^2\sum_{k=0}^3 c_{hk}t^h\theta^k                 \tag{3}
\]
and a polynomial
\[
 U=\sum_{k=1}^5\sum_{i=1}^9\sum_{j=0}^8u_{ijk}x^iy^jt^k,
 \qquad V(x,y,t)=U(y,x,t),                                  \tag{4}
\]
such that
\[
 \mathcal M(D^{-1/2})=
 \partial_x(U/D^{5/2})+\partial_y(V/D^{5/2}).                 \tag{5}
\]
The pole orders motivate this ansatz: a third parameter derivative has pole
order \(7/2\), and differentiating a primitive with pole order \(5/2\)
produces that order. The x/y symmetry halves the search. The finite degree
bounds in (4) are a successful search space, not a theorem that every possible
certificate must have those bounds.

Multiplying (5) by \(2D^{7/2}\) yields a polynomial coefficient system,
linear in all 405 coefficients of U and all eight coefficients of M:
\[
 2D(U_x+V_y)-5(UD_x+VD_y)
 =2\left(N_3+\sum c_{hk}t^hN_kD^{3-k}\right).              \tag{6}
\]
There are 1,180 coefficient equations. Exact rational row reduction gives rank
383, leaving 30 choices of primitive but uniquely fixing the eight operator
coefficients **within this ansatz**:
\[
 (c_{13},c_{12},c_{11},c_{10})=(-34,-51,-27,-5),
\]
\[
 (c_{23},c_{22},c_{21},c_{20})=(1,3,3,1).
\]
Thus the geometrically found operator is
\[
 \boxed{\mathcal L=\theta^3
 -t(2\theta+1)(17\theta^2+17\theta+5)+t^2(\theta+1)^3.}    \tag{7}
\]
No global uniqueness/minimal-order result is asserted from this restricted
linear system. Existence of the displayed exact identity is all the proof
needs. One choice of U has 231 nonzero rational coefficients, denominators
dividing 12, in `geometric_certificate_U.json`. The identity is checked by
expansion, so the proof does not depend on trusting a row-reduction algorithm.

For a sufficiently small parameter, the square-root branch and the primitives
in (5) are single valued on the integration torus. Integration of (5) therefore
proves \(\mathcal L\sum A_nt^n=0\). Extracting coefficients gives the full
Apéry recurrence. This derives its coefficients from the geometry, not from
an already established binomial-sum recurrence.

## 3. Lifting back to a rational three-variable divergence

For each \(v\in\{x,y,t\}\), put
\[
 h_v=-\frac{D_v}{2D},\qquad
 R_v=\left(\frac{a_v}{a}-\frac{D_v}{2D}\right)z
       +\frac{b_v-bD_v/(2D)}{2a},\qquad q_v=R_v/f.
\]
Direct differentiation, using \(D=b^2-4ac\), gives the elementary identity
\[
 \boxed{\partial_v(f^{-1})=h_vf^{-1}+\partial_zq_v.}       \tag{8}
\]
This quadratic-resolvent identity is the reason the algebraic residue
calculation can be lifted to a rational differential form.

Let \(S_k=N_k/D^k\), and define rational functions
\[
 T_0=0,\qquad T_{k+1}=\theta T_k+tS_kq_t.
\]
An induction with (8) proves
\[
 \theta^k(f^{-1})=S_kf^{-1}+\partial_zT_k.
\]
Write \(\mathcal L=\sum_{k=0}^3\ell_k(t)\theta^k\), where
\[
 (\ell_3,\ell_2,\ell_1,\ell_0)
 =(1-34t+t^2,-51t+3t^2,-27t+3t^2,-5t+t^2).
\]
Put
\[
 T_{\mathcal L}=\sum\ell_kT_k,\quad
 B=T_{\mathcal L}-\frac U{D^2}q_x-\frac V{D^2}q_y.
\]
Equations (6) and (8) now prove the explicit rational certificate
\[
 \boxed{\mathcal L(f^{-1})=
 \partial_x\!\left(\frac U{D^2f}\right)
 +\partial_y\!\left(\frac V{D^2f}\right)+\partial_zB.}     \tag{9}
\]
The formulas define B explicitly without expanding a large rational function.
Every q_v is O(1/z) at infinity, and parameter derivatives preserve that bound,
so \(B(\infty)=0\) for fixed positive x,y,t away from its poles.

## 4. The exact z=0 boundary

Use the same formulas with a real parameter t greater than
\(\lambda=17+12\sqrt2\). Then \(f<0\) on the positive octant, with no pole.
The boundary value \(W_0=B|_{z=0}\) is rational. In particular,
\[
 q_t(0)=\frac{2s}{tD},\qquad
 q_x(0)=\frac{B_x}{xpD},\qquad
 B_x=b(p+3xp_x)-2xp\,b_x,                                 \tag{10}
\]
and \(q_y(0)\) is the x/y transpose. Here the symbol B_x in (10) denotes a
polynomial, not the derivative of the three-variable primitive B.
For clarity the verification script uses precisely this polynomial convention.

The boundary values of the T primitives are particularly short:
\[
 T_1(0)=\frac{2s}D,\qquad T_2(0)=-\frac{3stD_t}{D^2},
\]
\[
 T_3(0)=-\frac{4s(tD_t+t^2D_{tt})}{D^2}
          +\frac{15st^2D_t^2}{2D^3}.                       \tag{11}
\]
Thus \(W_0=P/(pD^3)\), where P is explicitly obtained from (10)–(11) and U.
A second exact rational coefficient calculation gives a polynomial E with
188 nonzero coefficients such that
\[
 \boxed{W_0=
 \partial_x\!\left(\frac E{pD^2}\right)
 +\partial_y\!\left(\frac {E(y,x,t)}{pD^2}\right).}         \tag{12}
\]
Its full coefficients are in `geometric_boundary_E.json`; the polynomial
identity obtained by clearing denominators in (12) is checked exactly.

The support bounds needed for the improper limits are
\[
\begin{array}{c|ccc}
 &\deg_x&\deg_y&\text{total degree in }x,y\\\hline
 U&9&8&14\\
 E&8&7&13\\
 P&11&11&17
\end{array}
\]
and both U and E are divisible by x. The highest homogeneous pieces in x,y
are exceptionally simple:
\[
 D_6=t x^2y^2\big(t(x+y)^2-4xy\big),
\]
\[
 E_{13}=12t^3x^6y^5\big(t(x+y)^2-2(x^2+y^2)\big).          \tag{13}
\]
These support facts and the leading pieces are also checked from the saved
coefficient tables, not inferred from numerical fitting.

## 5. The surviving corner gives 12 by one elementary antiderivative

Integrate (12) over the square \([\varepsilon,R]^2\). Its finite-coordinate
boundaries tend to zero. On the infinite edge x=R, set y=Rv. Equations (13)
give the limiting integrated flux
\[
 \int_0^1
 \frac{12t\big(t(1+v)^2-2(1+v^2)\big)}
      {\big(t(1+v)^2-4v\big)^2}\,dv.
\]
Its antiderivative is
\[
 \frac{6t(v^2-1)}{t(1+v)^2-4v},
\]
so this edge contributes 6. The transposed edge contributes another 6. Hence
\[
 \boxed{\int_0^\infty\!\int_0^\infty W_0\,dx\,dy=12.}    \tag{14}
\]
The simultaneous x,y→infinity corner is essential. Taking a pointwise boundary
limit before integrating the other variable would miss it.

Here are quantitative justifications of all improper limits. Let
\(a_*=p+p^2/s\), \(b_*=2\sqrt{p^3/s}\), so
\(\Lambda=a_*+b_*\cosh v\) after translating log z. The positive minimum theorem bounds \(a_*+b_*\), not
\(a_*-b_*\); a different elementary estimate is needed here. Put
\(A=\sqrt p\) and \(v_0=\sqrt s\). Since \(x+y\ge2\sqrt{xy}\),
\(A\ge1+v_0\). Hence
\[
 \sqrt{p^2/s}-\sqrt p=A^2/v_0-A
 =A(A-v_0)/v_0\ge1+1/v_0>1.
\]
Therefore \(a_*-b_*=(\sqrt{p^2/s}-\sqrt p)^2>1\), and
\(a_*^2-b_*^2=(a_*-b_*)(a_*+b_*)\ge a_*\).
Using \(s^2(a_*^2-b_*^2)=p^2r^2\), formula (1) gives, for fixed t>λ>2,
\[
 D=t^2s^2(a_*^2-b_*^2)-2ts^2a_*+s^2
 \ge t^2\left(1-\frac2t\right)p^2r^2
 =:c_t p^2r^2>0.                                         \tag{15}
\]
This analytic inequality is proved above; it is not among the polynomial
identities checked by the verification script.
Consequently \(|E/(pD^2)|\le C_t|E|/(p^5r^4)\). The bound
\(\deg_y E\le7\), together with x divisibility, makes the x=ε edge integral
O(ε). On x=R, every monomial of total degree at most 12 contributes
O(log R/R): split y into (0,1) and (1,R), and bound its integral by a constant
times
\(R^{i-9}\int_1^R y^{j-5}dy\).
Only (13) remains. Its rescaled integrand is bounded on 0≤v≤1, so dominated
convergence proves the stated corner limit.

Similarly, \(|W_0|\le C_t|P|/(p^7r^6)\). This is integrable: on x≥y≥1,
a monomial is bounded by \(x^{i-13}y^{j-7}\); its y integral followed by its
x integral converges from i≤11 and i+j≤17. Regions where one variable is
bounded are easier, and the origin has no pole. Thus square exhaustion gives
the actual absolutely convergent integral in (14).

## 6. The real-chain operator and its inhomogeneous term

Define the convergent real resolvent
\[
 H(t)=\int_{(0,\infty)^3}\frac{dx\,dy\,dz}{f},\qquad t>\lambda.
\]
Its relation to the reciprocal generating function is
\[
 F(u)=\int\frac1{\Lambda-u}\,\frac{dx\,dy\,dz}{xyz}
     =\sum_{n\ge0} I_nu^n,\qquad H(t)=-t^{-1}F(1/t).       \tag{16}
\]
Parameter derivatives of H are dominated by a constant times
\(\Lambda^{-1}dx\,dy\,dz/(xyz)\), which is integrable. Thus differentiation
under the integral is justified.

First integrate (9) in z at fixed x,y, where the rational primitives have the
limits already specified. The z flux is −W₀. The x/y fluxes vanish on
exhausting squares. For the latter assertion, the elementary positive
z integral at u=0 is
\[
 \int_0^\infty\frac{dz}{xyz\Lambda}
   =\frac{\log(p/s)}{pr},\qquad
 \log(p/s)=\log(1+1/x)+\log(1+1/y).
\]
For the fixed parameter in (16) the absolute z integral is at most a constant
multiple of this expression. With (15) the x flux is therefore bounded by
\[
 C_t\frac{|U|\log(p/s)}{p^5r^5}.
\]
On x=ε its integrated bound is O(ε log(1/ε)); on x=R it is O(log R/R),
using i≤9, j≤8, and i+j≤14. For y≥1 the extra logarithm is bounded by
\(1/R+1/y\), which supplies the necessary decay; for y<1 its logarithmic
singularity is integrable. The y flux is identical by symmetry. This also
avoids an unjustified exchange of all noncompact boundary limits at once:
integrate z first, then exhaust the remaining square.

Equations (9) and (14) consequently give
\[
 \boxed{\mathcal L_t H(t)=-12.}                            \tag{17}
\]
Finally the elementary reciprocal operator identity
\[
 \mathcal L_u[u^{-1}h(1/u)]
       =-u(\mathcal L_t h)(1/u)
\]
and (16) imply
\[
 \boxed{\mathcal L_uF(u)=-12u.}                           \tag{18}
\]
Initially this holds on a positive interval near zero. Since F is an analytic
germ at zero (indeed the positive kernel gives convergence for |u|<λ), the
identity holds as a generating-function identity there. Its coefficients give
the homogeneous recurrence for I_n at n≥1 and the exceptional n=0 boundary
relation. The constant −12 has been computed from the exact-form boundary,
without using the already known recurrence or evaluating ζ(3).

## 7. Reproduction and limits of the result

From the repository root:

```sh
python codex/computations/derive_geometric_operator.py
python codex/computations/check_geometric_certificate.py
```

The first command uses `python-flint` and SymPy to rediscover all eight operator
coefficients and a primitive. The second uses SymPy only and the saved JSON
tables to verify (6), (8), (10)–(13), the support bounds, and the antiderivative
in Section 5. `search_geometric_boundary.py 1 --exact` reproduces the second
primitive; the exploratory modular searches are not required for verification.

The substantial remaining exposition problem is to compress the 231-term U
and 188-term E into a few natural vector fields or a short elliptic-family
calculation. The present derivation gives exact checkable geometric closure
and a small boundary calculation, but it does not yet make the dimension-three
closure itself as transparent as the one-variable Legendre reduction. The
successful search bounds should not be presented as a conceptual explanation
of that dimension. No literature-priority conclusion follows from these
independently derived certificates.
