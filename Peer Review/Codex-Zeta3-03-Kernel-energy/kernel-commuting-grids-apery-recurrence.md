# Two commuting Legendre grids force Apéry’s recurrence

Status: complete proof accepted in the [independent audit](kernel-grids-independent-review.md), with exact checks. This derives the recurrence directly from the logarithmic kernel and Legendre identities. It uses no Bernoulli functional and no rational-summation telescope. The objects and recurrence are classical; no historical novelty is asserted. The compatibility argument specializes the reviewed generic argument in [deformed-associated-pearson-symmetry.md](https://github.com/BoomydayCoder/Mathathon/blob/03cc75d694220d0fba888c637e573a29799e6d5b/codex/problem131/deformed-associated-pearson-symmetry.md), but here both grid equations follow immediately from one symmetric real kernel.

Put

\[
 L_n(x)=\sum_{j=0}^n(-1)^jC_{nj}x^j,
 \quad C_{nj}=\binom nj\binom{n+j}j,
 \quad \lambda_n=n(n+1),\quad\pi_n=2n+1,
\]

\[
 K(v)=-\frac{\log v}{1-v},\qquad
 I_{nm}=\iint_{(0,1)^2}L_n(x)L_m(y)K(xy)\,dx\,dy,
 \qquad I_n=I_{nn}.
\]

All these integrals are absolutely convergent, and I_(nm)=I_(mn). The [standalone kernel-energy proof](short-kernel-energy-zeta3-proof.md) identifies

\[
 I_n=-\sum_{s>=1}(P_n(s)^2)'=2(A_n\zeta(3)-B_n),
 \qquad A_n=\sum_jC_{nj}^2,
\]

with P_n the classical rational product, A_0=1, and B_0=0.

We prove the full sourced recurrence

\[
 \boxed{(n+1)^3I_{n+1}-(2n+1)(17n^2+17n+5)I_n+n^3I_{n-1}
        =-12\delta_{n0}.}                                      \tag{1}
\]

The previous-index term is omitted at n=0. The denominator A_n satisfies the corresponding homogeneous equation including n=0; the numerator B_n satisfies it for n>=1, with B_1=6.

## 1. The differential kernel gives the sourced grid

Let

\[
 J=-D[x(1-x)D],\qquad
 H=-D[x^2(1-x)D]+x.
\]

Ordinary monomial calculation gives

\[
 \boxed{H=\tfrac12(xJ+Jx)+\tfrac12.}                         \tag{2}
\]

The Legendre relations JL_n=lambda_n L_n and

\[
 xL_n=\frac12L_n-\frac{n+1}{2\pi_n}L_{n+1}
                    -\frac n{2\pi_n}L_{n-1}
\]

therefore give the cubic tridiagonal action

\[
 \boxed{2\pi_nHL_n=-(n+1)^3L_{n+1}
                 +\pi_n(\lambda_n+1)L_n-n^3L_{n-1}.}       \tag{3}
\]

The cube is a linear Legendre coefficient multiplied by the quadratic mean (lambda_n+lambda_(n+1))/2=(n+1)^2.

Set k_m(x)=integral L_m(y)K(xy)dy. The elementary kernel equation (H_x+J_y)K(xy)=0 and its unit endpoint flux prove

\[
 (H+\lambda_m)k_m=1.
\]

Explicitly, the Green boundary bracket at y=0 is +1 because y partial_y K(xy) tends to -1, and the bracket at y=1 is zero. These details and all endpoint estimates are proved in the standalone note. H can also be moved between a polynomial and k_m by ordinary integration by parts.

Integrate (3) against k_m and use integral L_n=delta_(n0). The result is

\[
 \boxed{(n+1)^3I_{n+1,m}
  -\pi_n(\lambda_n+2\lambda_m+1)I_{nm}
  +n^3I_{n-1,m}=-2\delta_{n0}.}                            \tag{4}
\]

Symmetry supplies the same equation in m. The source is retained; it is not a homogeneous grid at zero.

## 2. The denominator has the same two grids without a source

Define

\[
 Q_j(t)=\frac{\prod_{r=0}^{j-1}(t+r(r+1))}{(j!)^2},
 \quad U(x^j)=Q_j(t),\quad\Phi_n=UL_n.
\]

Monomial action gives UH=tU. Since Q_j(-lambda_m)=(-1)^j C_(mj) (zero when j>m),

\[
 K_{nm}:=\Phi_n(-\lambda_m)=\sum_jC_{nj}C_{mj}=K_{mn}.
\]

Apply U to (3) and evaluate at -lambda_m. This proves exactly (4) with K in place of I and with zero right side. Its diagonal is K_(nn)=A_n.

Thus the following compatibility argument applies at once to I and K. No separation of coefficients of an unknown irrational constant is used.

## 3. Commuting grid equations force the number six

Let X be either I or K, and let epsilon be respectively 2 or 0. Define the difference operator in the first coordinate by

\[
 (\mathcal D_nX)_{nm}
 =(n+1)^3(X_{n+1,m}-X_{nm})
       +n^3(X_{n-1,m}-X_{nm}).
\]

At n=0 the previous term is absent. Because
(n+1)^3+n^3=pi_n(lambda_n+1), the two grid equations become

\[
 \mathcal D_nX=2\pi_n\lambda_mX-\epsilon\delta_{n0},
 \qquad
 \mathcal D_mX=2\pi_m\lambda_nX-\epsilon\delta_{m0}.       \tag{5}
\]

The two difference operators commute. Apply D_m to the first equation and D_n to the second, then subtract. All source terms vanish: D1=0 and lambda_j delta_(j0)=0. The product terms cancel, leaving

\[
 \pi_n[\mathcal D_m,\lambda_m]X
       =\pi_m[\mathcal D_n,\lambda_n]X.
\]

The commutator is explicit, because consecutive eigenvalues differ by 2(n+1):

\[
 [\mathcal D_n,\lambda_n]X
       =2(n+1)^4X_{n+1,m}-2n^4X_{n-1,m}.
\]

Consequently

\[
 \pi_n[(m+1)^4X_{n,m+1}-m^4X_{n,m-1}]
 =\pi_m[(n+1)^4X_{n+1,m}-n^4X_{n-1,m}].                 \tag{6}
\]

Set m=n-1, n>=1. Abbreviate

\[
 d=X_{nn},\quad d_-=X_{n-1,n-1},\quad e=X_{n,n-1},
 \quad a=3n^2-n+1,\quad b=3n^2+n+1.
\]

The two grid equations express the distance-two entries as

\[
 (n+1)^3X_{n+1,n-1}=\pi_n a e-n^3d_-,
\]

\[
 (n-1)^3X_{n,n-2}=\pi_{n-1}b e-n^3d-\epsilon\delta_{n1}.
\]

In the second line the negative-index term is omitted when n=1. In (6), these equations occur multiplied respectively by n+1 and n-1. Thus the exceptional source is killed by (n-1)delta_(n1)=0. Substitution and cancellation give

\[
 n^3(d+d_-)=[(n+1)a+(n-1)b]e=6n^3e.
\]

We have proved the adjacent identity at every index, including n=1:

\[
 \boxed{X_{nn}+X_{n-1,n-1}=6X_{n,n-1}.}                  \tag{7}
\]

The number six is forced by compatibility of two quadratic grids and their cubic coefficients. It is not inferred from symmetry alone or fitted from initial values.

## 4. Closing the diagonal explains seventeen and five

Set m=n in the original grid. Symmetry and (7) give

\[
 (n+1)^3\frac{d_{n+1}+d_n}{6}
 -\pi_n(3\lambda_n+1)d_n
 +n^3\frac{d_n+d_{n-1}}6=-\epsilon\delta_{n0}.
\]

This is valid at n=0 too: the last term has zero coefficient, and (7) at index 1 supplies the forward term. Since
(n+1)^3+n^3=pi_n(lambda_n+1), the middle coefficient becomes

\[
 6\pi_n(3\lambda_n+1)-\pi_n(\lambda_n+1)
 =\pi_n(17\lambda_n+5).
\]

Thus **17=6·3−1 and 5=6·1−1**. We obtain

\[
 (n+1)^3d_{n+1}-\pi_n(17\lambda_n+5)d_n+n^3d_{n-1}
       =-6\epsilon\delta_{n0}.
\]

For d_n=I_n this is (1); for d_n=A_n it is homogeneous even at zero. Since I_0=2zeta(3), equation (1) gives

\[
 I_1-5I_0=-12,\qquad I_1=10\zeta(3)-12.
\]

Using I_n=2(A_n zeta(3)-B_n) and the separately proved A recurrence now gives

\[
 (n+1)^3B_{n+1}-\pi_n(17\lambda_n+5)B_n+n^3B_{n-1}
       =6\delta_{n0}.
\]

Therefore B_0=0, B_1=6, and both A and B satisfy the usual homogeneous recurrence for n>=1. No rationality or irrationality assumption on zeta(3) was used to derive these recurrences.

The same logarithmic kernel now supplies the rational-square linear form, positive energy estimate, and recurrence. For reproduction, run `python codex/computations/check_kernel_commuting_grids.py`.
