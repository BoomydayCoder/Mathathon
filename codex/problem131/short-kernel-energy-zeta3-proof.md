# A complete elementary ζ(3) proof through one logarithmic kernel and positive energy

Status: complete proof, including the shortened arithmetic and ordinary-coefficient presentation, accepted in [the independent review](kernel-energy-independent-review.md). Exact checks are in [the kernel-grid checker](../computations/check_kernel_commuting_grids.py) and [the bridge checker](../computations/check_log_kernel_legendre_bridge.py). The Apéry integers, Legendre polynomials, and logarithmic kernel are classical. No historical novelty is asserted.

One logarithmic kernel supplies both the arithmetic and a positive differential equation. Legendre orthogonality then determines the normalization and makes the error small. The proof uses no recurrence asymptotics, prime number theorem, Bernoulli numbers, hyperbolic measure, spectral theorem, or auxiliary polynomial transform. The final section explicitly recovers Zudilin's rational interpolation product.

## 1. Expanding the kernel gives the arithmetic immediately

For n>=1 put
\[
 L_n(x)=\frac1{n!}\frac{d^n}{dx^n}[x^n(1-x)^n]
       =\sum_{j=0}^n c_jx^j,\qquad
 c_j=(-1)^jC_{nj},\qquad C_{nj}=\binom nj\binom{n+j}j.
\]
Rodrigues integration by parts proves
\[
 \int_0^1L_n(x)p(x)\,dx=0\quad(\deg p<n),\qquad
 \int_0^1L_n(x)^2dx=\frac1{2n+1}.                         \tag{1}
\]
For the norm, the nth derivative of L_n is (-1)^n n! binom(2n,n); the beta integral integral x^n(1-x)^n dx=(n!)^2/(2n+1)! gives the displayed value.

Set
\[
 K(v)=-\frac{\log v}{1-v},\qquad
 I_n=\iint_{(0,1)^2}L_n(x)L_n(y)K(xy)\,dx\,dy.
\]
The removable value K(1)=1 is understood. Write H_j^(r)=sum_(k=1)^j k^(-r), with H_0^(r)=0. The geometric series for K gives the elementary monomial integrals
\[
 \boxed{\iint x^iy^jK(xy)\,dx\,dy=
 \begin{cases}
  2(\zeta(3)-H_i^{(3)}),&i=j,\\[2mm]
  (H_j^{(2)}-H_i^{(2)})/(j-i),&i\ne j.
 \end{cases}}                                               \tag{2}
\]
Indeed the summand is
\[
 \frac1{(r+i+1)^2(r+j+1)}
 +\frac1{(r+i+1)(r+j+1)^2}.
\]
For i!=j this equals
[1/(r+i+1)^2-1/(r+j+1)^2]/(j-i), and its sum telescopes to the displayed harmonic difference. For i=j it is twice a cubic reciprocal. All interchanges are justified by positivity for monomials; bounded polynomial factors are then dominated by the integrable kernel, whose double integral is 2zeta(3).

Expand L_n twice in (2). Pairing the two off-diagonal harmonic terms by interchanging i,j yields
\[
 \boxed{I_n=2(A_n\zeta(3)-B_n),\qquad
 A_n=\sum_{j=0}^n c_j^2=\sum_{j=0}^n C_{nj}^2
       \in\mathbb Z_{>0},}                                \tag{3}
\]
\[
 B_n=\sum_{i=0}^n c_i^2H_i^{(3)}
      +\sum_{i=0}^n c_iH_i^{(2)}
                          \sum_{\substack{0\le j\le n\\j\ne i}}
                                      \frac{c_j}{j-i}.     \tag{4}
\]
Let d_n=lcm(1,...,n). Every nonzero difference j-i has absolute value at most n, so d_n times the inner sum is integral. The harmonic sums need at most the remaining two powers. Hence
\[
 \boxed{d_n^3B_n\in\mathbb Z.}                            \tag{5}
\]
This proves the exact factor of two and the full arithmetic statement. Positivity of I_n remains to be proved; the signs of L_n do not supply it.

## 2. The kernel gives a positive differential equation

Define
\[
 k_n(x)=\int_0^1L_n(y)K(xy)\,dy,\qquad
 H=-\frac d{dx}\left[x^2(1-x)\frac d{dx}\right]+x,\qquad
 J=-\frac d{dy}\left[y(1-y)\frac d{dy}\right].
\]
Thus I_n=integral L_n(x)k_n(x)dx. Differentiating
v[(v-1)K(v)]'=1 gives
\[
 (v^2-v)K''+(3v-1)K'+K=0.
\]
Substitute v=xy to obtain
\[
 \boxed{(H_x+J_y)K(xy)=0.}                               \tag{6}
\]
The Legendre equation is JL_n=n(n+1)L_n. Its Green identity against K has one boundary contribution:
\[
 \begin{aligned}
 \int_0^1L_nJ_yK(xy)\,dy
 ={}&n(n+1)k_n(x)\\
 &+[-L_n y(1-y)\partial_yK(xy)+y(1-y)L_n'K(xy)]_0^1\\
 ={}&n(n+1)k_n(x)-1.
 \end{aligned}
\]
At y=1 the bracket vanishes. At y=0, y partial_y K(xy) tends to -1, L_n(0)=1, and yK(xy) tends to zero. Thus
\[
 \boxed{(H+n(n+1))k_n=1.}                               \tag{7}
\]
This is an equation for an explicit integral, not an assertion about an abstract inverse operator.

For the later integrations by parts, direct differentiation gives
\[
 k_n(x)=O_n(1+|\log x|),\qquad
 k_n'(x)=O_n(x^{-1}+|\log x|)\quad(x\downarrow0).
\]
At x=1 both are bounded, since K and its derivatives are regular at 1. For x<=1/2 the bounds follow directly from
K(xy)=-(log x+log y)/(1-xy). Local differentiation under the integral is dominated by an integrable multiple of 1+|log y|. Consequently k_n and its polynomial modifications lie in L^2, their weighted derivative energies are finite, and x^2(1-x)f'f vanishes at both endpoints.

## 3. A coefficient identity forces the right normalization

We want to modify A k_n by a polynomial r of degree below n so that its forcing term becomes L_n:
\[
 (H+n(n+1))(r+A k_n)=L_n.
\]
The needed polynomial division is particularly simple, because
\[
 Hx^j=(j+1)^2x^{j+1}-j(j+1)x^j.                          \tag{8}
\]
Starting with the coefficient of x^n and working down to x, the nonzero leading coefficients n^2,(n-1)^2,...,1 determine a unique polynomial r_n of degree below n and a unique constant A such that
\[
 L_n=(H+n(n+1))r_n+A.
\]

To identify the constant without solving for r_n, define the linear functional
\[
 E_n(x^j)=(-1)^jC_{nj}\quad(0\le j\le n).
\]
The elementary ratio
\[
 (j+1)^2C_{n,j+1}=(n-j)(n+j+1)C_{nj}
\]
and (8) prove
\[
 E_n((H+n(n+1))x^j)=0\quad(j<n),\qquad E_n(1)=1.         \tag{9}
\]
Applying E_n to the division identity forces
\[
 \boxed{A=E_n(L_n)=\sum_j C_{nj}^2=A_n,\qquad
 L_n=(H+n(n+1))r_n+A_n.}                                \tag{10}
\]
The same integer as in the arithmetic is therefore forced by the differential equation and the Legendre coefficients.

Put g_n=r_n+A_n k_n. Equations (7) and (10) give
(H+n(n+1))g_n=L_n. Orthogonality (1) gives
< L_n,g_n >=A_n I_n. Integration by parts, justified by the endpoint bounds above, now supplies the positive energy:
\[
 \boxed{A_nI_n=
 \int_0^1x^2(1-x)g_n'(x)^2dx
 +\int_0^1[x+n(n+1)]g_n(x)^2dx>0.}                       \tag{11}
\]
The inequality is strict because g_n solves an equation with nonzero right side L_n. Cauchy–Schwarz yields
\[
 n(n+1)\|g_n\|_2^2\le A_nI_n
 =\langle L_n,g_n\rangle\le\|L_n\|_2\|g_n\|_2.
\]
Using (1) proves
\[
 \boxed{0<I_n\le\frac1{n(n+1)(2n+1)A_n}.}                \tag{12}
\]
Thus the error is positive and at most a polynomial factor times the reciprocal of the integer A_n.

## 4. Two elementary growth bounds complete the contradiction


Set S_n=sum_(k=0)^n C_(nk). Since C_(nk)=binom(2k,k)binom(n+k,2k), ordinary formal power-series summation gives

\[
 \sum_{n>=0}S_nz^n
 =\frac1{\sqrt{1-6z+z^2}}
 =\frac1{\sqrt{1-\alpha z}\sqrt{1-\alpha^{-1}z}},
 \qquad \alpha=3+2\sqrt2.
\]

Both final factors have positive coefficients, so

\[
 S_n\ge\alpha^n\frac{\binom{2n}n}{4^n}
       \ge\frac{\alpha^n}{2n+1}.
\]

The last inequality uses the largest of the 2n+1 binomial coefficients summing to 4^n. Cauchy–Schwarz and (12) now imply, with Lambda=alpha^2=17+12sqrt2,

\[
 A_n\ge\frac{\Lambda^n}{(n+1)(2n+1)^2},\qquad
 \boxed{0<I_n\le\frac{2n+1}{n}\Lambda^{-n}\le3\Lambda^{-n}.}    \tag{13}
\]

Here is a complete elementary estimate for the lcm. Put psi(x)=log d_floor(x), with d_0=d_1=1, and

\[
 V(x)=\log\lfloor x\rfloor!-
       \log\lfloor x/2\rfloor!-2\log\lfloor x/3\rfloor!
       +\log\lfloor x/6\rfloor!.
\]

Counting prime powers in factorials gives

\[
 V(x)=\sum_{p\text{ prime},a>=1}h(x/p^a)\log p,
 \quad h(u)=\lfloor u\rfloor-\lfloor u/2\rfloor
            -2\lfloor u/3\rfloor+\lfloor u/6\rfloor.
\]

The function h has period six and values 0,1,1,0,0,1 on its successive unit intervals. It is nonnegative everywhere and equals one on [1,3), so

\[
 \psi(x)-\psi(x/3)\le V(x).
\]

The elementary integral estimate log(r!)=r log r-r+O(log(r+2)) gives

\[
 V(x)=\left(\frac{\log2}3+\frac{\log3}2\right)x+O(\log(x+2)).
\]

Sum at x,x/3,x/9,..., stopping once x/3^j<1. This proves

\[
 \boxed{\log d_n\le cn+O(\log^2(n+2)),\qquad
 c=\frac{\log2}2+\frac{3\log3}4.}                       \tag{14}
\]

The strict inequality 3c<log Lambda is exact. Since sqrt2>7/5, Lambda>169/5, while

\[
 e^{12c}=64\cdot3^9=1259712
          <\frac{169^4}{5^4}=\frac{815730721}{625}.
\]

Equations (13)–(14) therefore show d_n^3 I_n tends to zero, with every term positive.

If zeta(3)=r/q for positive integers r,q, (3)–(5) make

\[
 \boxed{q d_n^3 I_n=2r d_n^3A_n-2q d_n^3B_n}
\]

a positive integer for every n>=1. It tends to zero, a contradiction. Thus **zeta(3) is irrational**.

## Exact connection to Zudilin's rational interpolation and Beukers' integral

The rational product can be recovered directly, so it remains an explicit part of the motivation. Rodrigues integration, initially for s>n, gives
\[
 \int_0^1x^{s-1}L_n(x)dx=(-1)^nP_n(s),\qquad
 P_n(s)=\frac{\prod_{j=1}^n(s-j)}{\prod_{j=0}^n(s+j)}.
\]
Both sides are rational in s, so the identity holds for every real s>0. Expanding K(xy) geometrically and differentiating the two Mellin factors gives exactly
\[
 I_n=-\sum_{s=1}^{\infty}\frac d{ds}[P_n(s)^2].
\]
The poles are at consecutive integers and the first n positive integer zeros express Legendre orthogonality. The residues of P_n are
(-1)^(n-j)C_(nj); their products give precisely the arithmetic in (4). Thus the direct kernel computation and the rational-square partial fractions are identical arithmetic, expressed in different coordinates.

The Stieltjes identity
K(v)=integral_0^infinity du/[(1+u)(v+u)] and two Rodrigues integrations by parts give the same I_n as Beukers' positive cube integral. The complete changes of variables and convergence justification are in [the hybrid proof](short-elementary-hybrid-proof.md).

The [full kernel bridge](log-kernel-legendre-resolvent-bridge.md) recovers the Newton transform, positive Wilson measure, and Hilbert Gram identity. The [reciprocal Legendre design lemma](reciprocal-legendre-kernel-design.md) explains why H and K arise together: norm-preserving inversion gives H, and the product-kernel equation with regularity and unit boundary flux uniquely determines K. None of these comparisons is required for the contradiction above.
