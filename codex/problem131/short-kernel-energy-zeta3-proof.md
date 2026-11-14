# A complete elementary ζ(3) proof through one logarithmic kernel and positive energy

Status: complete proof accepted in the [independent full-proof audit](kernel-energy-independent-review.md); its new kernel/energy bridge has exact checks in [check_log_kernel_legendre_bridge.py](../computations/check_log_kernel_legendre_bridge.py). The independent audit is complete. The rational functions, Apéry integers, Legendre polynomials, and logarithmic kernel are classical. No historical novelty is asserted.

The proof uses no recurrence asymptotics, prime number theorem, Bernoulli numbers, hyperbolic measure, spectral theorem, or multiple integral transformation. One rational sum exposes the denominators. Its logarithmic integral kernel supplies an elementary positive differential equation. Polynomial division then turns Legendre orthogonality into a small positive error.

## 1. The rational square supplies the arithmetic

For n>=1 put

\[
 L_n(x)=\frac1{n!}\frac{d^n}{dx^n}[x^n(1-x)^n]
       =\sum_{k=0}^n(-1)^k C_{nk}x^k,
 \qquad C_{nk}=\binom nk\binom{n+k}k.
\]

Rodrigues integration by parts proves

\[
 \int_0^1L_n(x)p(x)\,dx=0\quad(\deg p<n),\qquad
 \int_0^1L_n(x)^2dx=\frac1{2n+1}.                         \tag{1}
\]

For the norm, the nth derivative of L_n is (-1)^n n! binom(2n,n); the beta integral integral x^n(1-x)^n dx=(n!)^2/(2n+1)! gives the displayed value.

The same integration by parts, initially for s>n, gives

\[
 \int_0^1x^{s-1}L_n(x)dx=(-1)^nP_n(s),\qquad
 P_n(s)=\frac{\prod_{j=1}^n(s-j)}{\prod_{j=0}^n(s+j)}.     \tag{2}
\]

Both sides of (2) are rational functions of s, so the identity holds for all real s>0. Define the absolutely convergent sum

\[
 I_n=-\sum_{s=1}^{\infty}\frac d{ds}[P_n(s)^2].            \tag{3}
\]

The residue of P_n at -k is c_k=(-1)^(n-k)C_(nk). Hence

\[
 P_n(s)^2=\sum_{k=0}^n
       \left(\frac{a_k}{(s+k)^2}+\frac{b_k}{s+k}\right),
 \quad a_k=c_k^2,\quad
 \frac{b_k}2=c_k\sum_{j\ne k}\frac{c_j}{j-k}.             \tag{4}
\]

Since P_n(s)^2=O(s^-2), its coefficient of s^-1 is zero: sum b_k=0. Write d_n=lcm(1,...,n) and H_k^(r)=sum_(j=1)^k j^-r, with H_0^(r)=0. Differentiating (4) and summing gives the exact normalization

\[
 \boxed{I_n=2(A_n\zeta(3)-B_n),\qquad
 A_n=\sum_{k=0}^nC_{nk}^2\in\mathbb Z_{>0},}              \tag{5}
\]

\[
 B_n=\sum_{k=0}^n\left(a_kH_k^{(3)}+\frac{b_k}2H_k^{(2)}\right),
 \qquad \boxed{d_n^3B_n\in\mathbb Z.}                    \tag{6}
\]

Indeed d_n b_k/2 is integral, because every nonzero difference j-k has absolute value at most n. The harmonic denominators in (6) use at most the remaining two powers of d_n. The zeta(2) coefficient vanishes by sum b_k=0. No arithmetic estimate on an orthogonal polynomial or its moments is being assumed.

## 2. The logarithmic kernel gives a positive differential equation

Put

\[
 K(v)=-\frac{\log v}{1-v},\qquad
 k_n(x)=\int_0^1L_n(y)K(xy)\,dy.
\]

The removable value K(1)=1 is understood. The geometric expansion of K and (2) give

\[
 \boxed{I_n=\iint_{(0,1)^2}L_n(x)L_n(y)K(xy)\,dx\,dy
             =\int_0^1L_n(x)k_n(x)dx.}                  \tag{7}
\]

The logarithm differentiates the Mellin product in (2). Interchanges are absolutely justified: the polynomials are bounded, K is positive, and Tonelli gives integral integral K(xy) dxdy=2sum_(s>=1)s^-3<infinity.

Define two ordinary differential operators

\[
 H=-\frac d{dx}\left[x^2(1-x)\frac d{dx}\right]+x,
 \qquad J=-\frac d{dy}\left[y(1-y)\frac d{dy}\right].
\]

Differentiating v[(v-1)K(v)]'=1 yields

\[
 (v^2-v)K''+(3v-1)K'+K=0.
\]

Substitute v=xy to obtain

\[
 \boxed{(H_x+J_y)K(xy)=0.}                               \tag{8}
\]

The Legendre equation is JL_n=n(n+1)L_n. Its Green identity against K has one boundary contribution:

\[
 \int_0^1L_nJ_yK(xy)\,dy
 =n(n+1)k_n(x)
  +[-L_n y(1-y)\partial_yK(xy)+y(1-y)L_n'K(xy)]_0^1
 =n(n+1)k_n(x)-1.                                       \tag{9}
\]

At y=1 the bracket vanishes. At y=0, y partial_y K(xy) tends to -1, L_n(0)=1, and yK(xy) tends to zero. Equations (8)–(9) prove the concrete inhomogeneous equation

\[
 \boxed{(H+n(n+1))k_n=1.}                               \tag{10}
\]

No operator closure or abstract inverse is needed.

For clarity about later boundary terms, direct differentiation of the defining integral gives

\[
 k_n(x)=O_n(1+|\log x|),\qquad
 k_n'(x)=O_n(x^{-1}+|\log x|)\quad(x\downarrow0).
\]

At x=1 both are bounded, since K and its derivatives are regular at 1. For x<=1/2 these bounds follow immediately from K(xy)=-(log x+log y)/(1-xy). Local differentiation under the integral is dominated by an integrable multiple of 1+|log y|. Consequently k_n and the polynomial modifications below lie in L^2, their weighted derivative energies are finite, and x^2(1-x)f'f tends to zero at both endpoints.

## 3. Polynomial division forces the Apéry coefficient and positivity

The monomial action of H is

\[
 Hx^j=(j+1)^2x^{j+1}-j(j+1)x^j.                          \tag{11}
\]

Introduce a degree-preserving polynomial map U by

\[
 U(x^j)=Q_j(t),\qquad
 Q_j(t)=\frac{\prod_{r=0}^{j-1}(t+r(r+1))}{(j!)^2}.
\]

Equation (11) says exactly UH=tU. Since U1=1, it also gives U^(-1)p(t)=p(H)1. At the negative integer-grid point t=-n(n+1),

\[
 Q_j(-n(n+1))=\prod_{r=0}^{j-1}\frac{(r-n)(r+n+1)}{(r+1)^2}
                 =(-1)^jC_{nj}\quad(0<=j<=n).
\]

Therefore, with Phi_n=UL_n,

\[
 \Phi_n(-n(n+1))=\sum_{j=0}^nC_{nj}^2=A_n.             \tag{12}
\]

Ordinary polynomial division now supplies a polynomial r_n of degree at most n-1 such that

\[
 L_n=(H+n(n+1))r_n+A_n.
\]

Explicitly, apply U^(-1) to (Phi_n(t)-A_n)/(t+n(n+1)). Define

\[
 g_n=r_n+A_nk_n.
\]

Equation (10) gives (H+n(n+1))g_n=L_n. By (1) and (7),

\[
 \langle L_n,g_n\rangle=A_n I_n.
\]

Ordinary integration by parts is legitimate by the endpoint bounds in Section 2. It yields the positive energy identity

\[
 \boxed{A_nI_n=
 \int_0^1x^2(1-x)g_n'(x)^2dx
 +\int_0^1[x+n(n+1)]g_n(x)^2dx>0.}                     \tag{13}
\]

The inequality is strict because g_n solves an equation with nonzero right side L_n. Since n>=1, (13) and Cauchy–Schwarz give

\[
 n(n+1)\|g_n\|_2^2
 \le A_nI_n=\langle L_n,g_n\rangle
 \le\|L_n\|_2\|g_n\|_2.
\]

Together with (1), this proves

\[
 \boxed{0<I_n\le\frac1{n(n+1)(2n+1)A_n}.}             \tag{14}
\]

The same A_n appears in the arithmetic and in the positive estimate because its binomial factors arise twice: as coefficients of L_n and as the Newton evaluations (12). No half-square shift, Bernoulli normalization, or measure identification enters this proof.

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

The last inequality uses the largest of the 2n+1 binomial coefficients summing to 4^n. Cauchy–Schwarz and (14) now imply, with Lambda=alpha^2=17+12sqrt2,

\[
 A_n\ge\frac{\Lambda^n}{(n+1)(2n+1)^2},\qquad
 \boxed{0<I_n\le\frac{2n+1}{n}\Lambda^{-n}\le3\Lambda^{-n}.}    \tag{15}
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
 c=\frac{\log2}2+\frac{3\log3}4.}                       \tag{16}
\]

The strict inequality 3c<log Lambda is exact. Since sqrt2>7/5, Lambda>169/5, while

\[
 e^{12c}=64\cdot3^9=1259712
          <\frac{169^4}{5^4}=\frac{815730721}{625}.
\]

Equations (15)–(16) therefore show d_n^3 I_n tends to zero, with every term positive.

If zeta(3)=r/q for positive integers r,q, (5)–(6) make

\[
 \boxed{q d_n^3 I_n=2r d_n^3A_n-2q d_n^3B_n}
\]

a positive integer for every n>=1. It tends to zero, a contradiction. Thus **zeta(3) is irrational**.

## Exact relationship with the other classical representations

The full [log-kernel bridge](log-kernel-legendre-resolvent-bridge.md) proves that shifting H to T=H-1/4 converts the integer grid j(j+1) to (j+1/2)^2 and recovers the positive Wilson/Padé measure. It derives its Hilbert Gram matrix without Bernoulli moments. The same I_n is the Beukers positive cube integral by the two real integrations by parts already proved in [short-elementary-hybrid-proof.md](short-elementary-hybrid-proof.md).

The proof above does not require either comparison. It makes the rational-square arithmetic and the positive error estimate meet directly at the elementary kernel equation (8).
