# The missing Wilson boundary term is explicit

Status: all-parameter boundary identity and rational Newton-moment formula
accepted in the [independent audit](rank-one-boundary-independent-review.md). These characterize the original positive rational modification,
not the formal negative-parameter Wilson orthogonality functional. They do not
yet prove its adjacent-Apéry Schur-root identity. No novelty claim is made.

This note was developed during the 2026-11-14 continued Schur investigation.
Its finite checks are in `../computations/check_rank_one_boundary_pearson.py`.

## 1. A rank-two defect, rather than false orthogonality

Let M be a positive integer, D_M(t)=product_{j=1}^M(t+j²), and define the
functional on even polynomials in x by

\[
 \ell_M(p)=\int_0^\infty p(iy)\frac{d\mu_3(y^2)}{D_M(y^2)},\qquad
 d\mu_3(y^2)=\pi^2y^3\coth(\pi y)\operatorname{csch}^2(\pi y)\,dy.
\]

Put

\[
 A_M(x)=\frac{(x+1)^3(x-M)}{2x(2x+1)},\qquad
 (\mathcal D_Mp)(x)=A_M(x)[p(x+1)-p(x)]
                   +A_M(-x)[p(x-1)-p(x)],
 \qquad c_M=\frac{M}{2(M!)^2}.
\]

The apparent poles at x=0 and x=±1/2 cancel for every even polynomial, so
D_M preserves even polynomials. (Here the calligraphic D is the operator;
the ordinary D_M is the product polynomial.) For every pair of even
polynomials p,q,

\[
 \boxed{\ell_M(p\mathcal D_Mq-q\mathcal D_Mp)
       =-c_M[p(0)q(1)-p(1)q(0)].}                 \tag{1}
\]

In particular,

\[
 \boxed{\ell_M(\mathcal D_Mq)=c_M[q(0)-q(1)].}     \tag{2}
\]

Thus the original continuous measure has a precise boundary defect. It is
not legitimate to replace the right side of (1) by zero and invoke ordinary
Wilson orthogonality at parameters (1,1,1,-M).

Here is a contour proof with the poles accounted for. Set

\[
 w_M(x)=\frac{\pi^2x^3\cot(\pi x)\csc^2(\pi x)}{D_M(-x^2)}.
\]

It is even, is regular at zero with w_M(0)=1/[pi(M!)²], and satisfies

\[
 w_M(x+1)A_M(-x-1)=w_M(x)A_M(x).
\]

For K(x)=p(x)q(x+1)-q(x)p(x+1), put F=w_M A_M K. The preceding
identity and evenness give F(-x-1)=-F(x). The pole of A_M at -1/2 is
removable because K(-1/2)=0. In the closed strip -1<=Re x<=0 the only
remaining poles of F are simple endpoint poles at 0 and -1. At -1 the
factor (x+1)^3 cancels the order-three trigonometric pole, leaving the
single pole from D_M(-x²). Their residues are equal, and the residue at
zero is

\[
 r=-\frac{M}{2\pi(M!)^2}[p(0)q(1)-p(1)q(0)].
\]

Integrate first between vertical lines Re x=epsilon and Re x=-1-epsilon,
0<epsilon<1/4. Horizontal integrals vanish as their heights tend to infinity,
because the trigonometric factor decreases exponentially and all other
factors grow at most polynomially. Reflection gives opposite vertical
integrals, so the residue theorem gives integral_{epsilon+iR} F dx=2pi i r.
As epsilon decreases to zero, the usual explicitly integrable simple pole
r/x contributes i pi r, leaving

\[
 \operatorname{PV}\int_{i\mathbb R}F(x)\,dx=\pi i r.
\]

One can verify that limiting statement by subtracting r/[x(1-x²)] near
the line, or directly integrating r/x on a finite symmetric segment; the
remaining analytic function converges ordinarily. On the imaginary axis,
w_M(iy) equals the positive density divided by D_M. Evenness reduces the
two shifted terms in the left side of (1) to
(1/i) PV integral F dx. This equals pi r, proving (1), with its sign and
normalization. Taking p=1 proves (2).

## 2. Triangular action and an explicit rational Newton moment

Let

\[
 \phi_k(x)=\prod_{j=1}^k(x^2-j^2),\qquad
 Q_k(t)=\frac{D_k(t)}{k!(k+1)!},\qquad
 a_k^{(M)}=\int\frac{Q_k(t)}{D_M(t)}\,d\mu_3(t).
\]

An elementary polynomial identity is

\[
 \boxed{\mathcal D_M\phi_k
  =k(k+2-M)\phi_k+k(k+1)^2(k-M)\phi_{k-1}.}       \tag{3}
\]

The operator preserves degree and is triangular in this Newton basis:
all roots ±1,...,±(k-1) of phi_k remain roots of its image. Its leading
coefficient is k(k+2-M); the remaining coefficient is obtained by evaluating
at x=k, where only A_M(k)phi_k(k+1) survives. This also proves (3) for k=1.
Since phi_k(1)=0 and phi_k(0)=(-1)^k(k!)², equation (2) gives

\[
 \boxed{k(k+2-M)a_k^{(M)}
       +(k+1)(M-k)a_{k-1}^{(M)}=\frac{c_M}{k+1}.} \tag{4}
\]

At k=M-2 the coefficient of a_k vanishes. For M>=3 this determines

\[
 a_{M-3}^{(M)}=\frac{c_M}{2(M-1)^2}.
\]

Backward substitution then determines every a_k with 0<=k<=M-3,
without any zeta value or unknown mass. The resulting exact formula is

\[
 \boxed{
 a_k^{(M)}=
 c_M\frac{k+1}{(M-k-2)(M-k-1)}
 \sum_{r=k+1}^{M-2}\frac{M-r-1}{r(r+1)^2}
 }                                                        \tag{5}
\]

or, equivalently,

\[
 \boxed{
 a_k^{(M)}=c_M
 \frac{M-k-2-M(k+1)(H_{M-1}^{(2)}-H_{k+1}^{(2)})}
 {(M-k-2)(M-k-1)}.}                              \tag{6}
\]

To check the sum, write the backward coefficient in (4) as
b_r=r(M-r-2)/[(r+1)(M-r)]. The product from r=k+1 through r=s-1 is
(k+1)(M-s-1)(M-s)/[s(M-k-2)(M-k-1)]. Multiplication by the source
c_M/[(s+1)²(M-s)] gives (5). Finally

\[
 \frac{M-r-1}{r(r+1)^2}
  =(M-1)(1/r-1/(r+1))-M/(r+1)^2
\]

gives (6). Formula (5) also makes positivity transparent.

For M=2h, these Q_k span every degree through 2h-3. Therefore (5) gives
all lower and cross moments needed to determine the original degree-h-1
Schur minimizer, by a single finite reciprocal-square sum. The next moment
at degree 2h-2 contains the zeta parameter and is not determined by (4):
this is exactly the original rank-one degree drop. This is an all-h
characterization of the original problem, not a replacement measure.

## 3. A further exact characterization of the minimizer

Let p(x)=Q_h(x²) be the original monic degree-h-1 polynomial in x², and
write n=h-1 and M=2h. The leading operator coefficient is lambda_n=-n².
For any even q of degree at most n-1 in x², orthogonality and (1) give

\[
 \boxed{\ell_M((\mathcal D_Mp+n^2p)q)
       =c_M[p(0)q(1)-p(1)q(0)].}                  \tag{7}
\]

The polynomial on the left has degree at most n-1 in x². Thus it is the
explicit combination c_M[p(0)K_{n-1}(x²,1)-p(1)K_{n-1}(x²,0)] of the
reproducing kernels of the lower polynomial space (using its positive
measure on negative x²). This is the correct inhomogeneous equation for
the minimizer. The ordinary negative-parameter Wilson polynomial solves
the homogeneous equation and generally has different coefficients.

## 4. Finite patterns still awaiting proof

The exact computations separately test the following candidate endpoint
values for the original monic Q_h:

\[
 Q_h(0)\stackrel?=\frac{12(2h)!}{h^2(5A_h-A_{h-1})},\qquad
 \frac{Q_h(1)}{Q_h(0)}\stackrel?=\frac{h^2(h+1)}2.    \tag{8}
\]

They hold in the tested range. Their product would give the previously
observed leading numerator coefficient
72(h+1)/[h(5A_h-A_{h-1})²], but (1)-(7) do not yet prove (8).

If R_h=Q_h/Q_h(0) is expanded in the Newton basis
product_{j=1}^k(z-j²), the first two tested coefficients are

\[
 a_0\stackrel?=h^2(h+1)/2,
 \qquad
 a_1\stackrel?=
 \frac{h^2(h+1)(h-1)(h^3+2h^2+3h-2)}{16(2h-1)}.
\]

The second formula concerns h>=2. It is recorded only as a promising
dual-grid target; it was reconstructed from exact values and has no proof
here. The original all-h adjacent-Apéry quotient remains open.
