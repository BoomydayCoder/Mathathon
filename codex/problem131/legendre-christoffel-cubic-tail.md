# The cubic tail is a finite orthogonal projection

Status: accepted in [the independent review](legendre-christoffel-independent-review.md). This is a shorter, purely polynomial explanation of the normalization in the [log-kernel energy proof](log-kernel-legendre-resolvent-bridge.md). It is equivalent to the finite Jacobi deletion identity in [the boundary-kernel supplement](jacobi-boundary-kernel-spectrum.md), but requires no Jacobi formulas, spectral expansion, completeness theorem, or Legendre derivative ladder. The ingredients are classical; no historical novelty is asserted.

The cubic denominator has a simple source: deleting the kth orthogonal mode uses coefficient 1/k, while that mode has squared norm 2/k. Pythagoras therefore removes exactly 2/k³.

## 1. Adjacent Legendre polynomials provide the weighted modes

Let
\[
 L_n(v)=\frac1{n!}\frac{d^n}{dv^n}[v^n(1-v)^n],\qquad
 \kappa_n=(-1)^n\binom{2n}{n}.
\]
Their standard elementary properties are
\[
 L_n(0)=1,\quad L_n(1)=(-1)^n,\quad
 [v^n]L_n=\kappa_n,\quad
 \int_0^1 L_n^2\,dv=\frac1{2n+1},
\]
and L_n is orthogonal to all polynomials of degree below n. Rodrigues integration by parts proves the orthogonality and norm.

For k>=1 define
\[
 Q_k(v)=\frac{L_k(v)+L_{k-1}(v)}{1-v},\qquad
 \langle f,g\rangle_w=\int_0^1(1-v)f(v)g(v)\,dv.
\]
The numerator vanishes at one, so Q_k is a polynomial of degree k-1 with leading coefficient -kappa_k. For deg q<k-1,
\[
 \langle Q_k,q\rangle_w=\int_0^1(L_k+L_{k-1})q\,dv=0.
\]
Consequently the Q_k are mutually orthogonal. Their exact norms need only a leading-coefficient comparison:
\[
 \begin{aligned}
 \|Q_k\|_w^2
 &=\int_0^1 Q_k(L_k+L_{k-1})\,dv\\
 &=\frac{-\kappa_k}{\kappa_{k-1}}
                \int_0^1L_{k-1}^2\,dv\\
 &=\frac{2(2k-1)}k\frac1{2k-1}
 =\boxed{\frac2k}.                                      \tag{1}
 \end{aligned}
\]
This adjacent-polynomial construction is the elementary Christoffel modification from weight 1 to weight 1-v.

## 2. A Volterra integral forces the deletion coefficient

Set
\[
 K_0(v)=-\frac{\log v}{1-v},\qquad
 G_n(v)=\frac1{1-v}\int_v^1\frac{L_n(x)}x\,dx.
\]
In particular G_0=K_0. If L_n(x)=1+sum_(j=1)^n c_(nj)x^j, then
\[
 G_n-K_0=\sum_{j=1}^n\frac{c_{nj}}j
                       (1+v+\cdots+v^{j-1}).            \tag{2}
\]
Thus G_n-K_0 has degree at most n-1 and leading coefficient kappa_n/n. Both functions are in L²_w: at zero they have logarithmic growth, and at one their apparent quotients have finite limits.

For every polynomial q with deg q<n, absolute Fubini gives
\[
 \langle G_n,q\rangle_w
 =\int_0^1L_n(x)\left[\frac1x\int_0^xq(v)\,dv\right]dx=0. \tag{3}
\]
The bracket is a polynomial of the same degree as q. Absolute convergence follows, for example, by bounding q and L_n on the compact interval: the inner integral of |q| is O(x).

The difference G_n-G_(n-1) is therefore a polynomial of degree n-1 orthogonal to all lower degrees for weight 1-v. It is a multiple of Q_n. Equation (2) gives its leading coefficient kappa_n/n, whereas Q_n has leading coefficient -kappa_n. Hence the multiple is forced:
\[
 \boxed{G_n-G_{n-1}=-\frac1nQ_n,\qquad
        K_0=G_n+\sum_{k=1}^n\frac1kQ_k.}                \tag{4}
\]
This is a finite identity. It does not use an infinite orthogonal expansion or a limiting argument.

## 3. Pythagoras gives the zeta-three tail

The double-kernel boundary row in the energy proof is
\[
 I_{n0}=\int_0^1\int_0^1L_n(x)K_0(xy)\,dx\,dy.
\]
Put v=xy and reverse the two integrations; bounded L_n and the integrable logarithmic kernel justify this absolutely. The result is
\[
 I_{n0}=\int_0^1K_0(v)\int_v^1\frac{L_n(x)}x\,dx\,dv
       =\langle K_0,G_n\rangle_w.
\]
By (3)-(4), G_n is orthogonal to its finite polynomial correction. Therefore
\[
 \begin{aligned}
 I_{n0}=\|G_n\|_w^2
 &=\|K_0\|_w^2-\sum_{k=1}^n\frac{\|Q_k\|_w^2}{k^2}\\
 &=\|K_0\|_w^2-2\sum_{k=1}^n\frac1{k^3}.              \tag{5}
 \end{aligned}
\]
Finally, a nonnegative geometric series and Tonelli give
\[
 \|K_0\|_w^2=\int_0^1\frac{\log^2v}{1-v}\,dv
            =\sum_{j=0}^\infty\frac2{(j+1)^3}=2\zeta(3).
\]
Thus
\[
 \boxed{\displaystyle
 I_{n0}=\|G_n\|_w^2
       =2\left(\zeta(3)-\sum_{k=1}^n\frac1{k^3}\right).} \tag{6}
\]
The norm is strictly positive because G_n retains the logarithmic singularity at zero. The cubic decrement is exactly the squared norm of one deleted mode, not an independently imposed normalization.

For comparison with the Jacobi notation only,
\[
 Q_k(v)=2(-1)^{k-1}P_{k-1}^{(1,0)}(2v-1).
\]
Nothing in the finite proof requires this identification. The kernel energy, integer coefficients, and lcm comparison remain the distinct ingredients that prove irrationality of zeta(3); equation (6) supplies their exact boundary normalization.

Exact checks: [check_christoffel_legendre_deletion.py](../computations/check_christoffel_legendre_deletion.py).
