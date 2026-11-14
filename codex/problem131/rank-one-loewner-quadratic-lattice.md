# The original Schur root is a Loewner singularity and a quadratic-lattice Riccati problem

Status: constructive bounded investigation. The determinant, interpolation,
and Riccati identities below are proved for every h. Exact checks through h=6
pass. The desired original-root/adjacent-Apéry identification remains open;
neither these reformulations nor the finite matches prove it. Reviewed
moving-pole Wilson results have not been substituted for the original D₂ₕ
measure. No historical novelty claim is made.

The root agent proposed the Loewner route. The exact determinant factors,
interpolation characterization, and harmonic-free Riccati identity were derived
independently here. They provide a concrete new algebraic target for the gap in
`rank-one-apery-bridge-draft.md`.

## 1. Exact Loewner matrix and its determinant factor

Use the established affine ζ3 functional and write
\[
 f_k(X)=\mu_{3,X}((t+k^2)^{-1})=k^2X+\beta_k,
\quad
 \beta_k=-k^2H_k^{(3)}-\frac12+\frac1{2k}.
\]
For i=1,…,h and b=h+1,…,2h set
\[
 K_h(X)_{i,b}=\frac{f_i(X)-f_b(X)}{b^2-i^2}
             =B_{i,b}-X,
\quad B_{i,b}=\frac{\beta_i-\beta_b}{b^2-i^2}.            \tag{1}
\]
The sign is important: the X coefficient is −1, not +1. These are the cross
Cauchy Gram entries because
\[
 \frac1{(t+i^2)(t+b^2)}=
 \frac{(t+i^2)^{-1}-(t+b^2)^{-1}}{b^2-i^2}.
\]
Let Gₕ(X) be the original h×h monomial Gram matrix
\(\mu_{3,X}(t^{r+s}/D_{2h})\), r,s=0,…,h−1. Define
\[
 \Delta_L=\prod_{1\le i<j\le h}(j^2-i^2),\quad
 \Delta_R=\prod_{h+1\le i<j\le2h}(j^2-i^2).
\]
Then
\[
 \boxed{\det K_h(X)=\Delta_L\Delta_R\det G_h(X).}        \tag{2}
\]

Proof: put A(t)=Dₕ(t), C(t)=D₂ₕ(t)/Dₕ(t). The two polynomial bases
A(t)/(t+i²) and C(t)/(t+b²) have degree≤h−1, and their cross Gram matrix
for μ₃,ₓ/D₂ₕ is exactly (1). Their coefficient-matrix determinants are
Δ_L and Δ_R. To verify their signs, evaluate the first basis at −j²:
the evaluation matrix is diagonal with entries A′(−j²). Dividing its
determinant by the Vandermonde on the decreasing roots −j² gives the positive
Vandermonde Δ_L. The right basis is identical in form. This proves (2).

The original construction proves Gₕ(X)=Gₕ(0)−Xeeᵀ, where e is the final
coordinate vector, and Gₕ(0) is positive definite. Indeed Gₕ(0)=Gₕ(ζ3)+ζ3eeᵀ.
Consequently B is invertible and, with 1 the all-one h-vector,
\[
 \boxed{q_{h,3}=\frac1{\mathbf1^TB^{-1}\mathbf1}.}        \tag{3}
\]
This follows either from the matrix determinant lemma applied to
K=B−X11ᵀ, or by transforming the original rank-one formula using the two
monic polynomial bases. It is an exact all-h identity.

## 2. A unique rational interpolation problem

Let pₕ(t) be the original monic degree-h−1 minimizer for μ₃/D₂ₕ and put
\[
 Q_h(z)=(-1)^{h-1}p_h(-z).
\]
All roots of Qₕ are negative, so Qₕ(j²)≠0 for every positive integer j.
At X=q=qₕ,₃ there is a polynomial Nₕ of degree≤h−1 such that
\[
 \boxed{\frac{N_h(j^2)}{Q_h(j^2)}=f_j(q),
                   \qquad j=1,\ldots,2h.}               \tag{4}
\]
The denominator is monic, and this pair is unique.

Here is the direct nullvector construction, including all normalizations.
Let C_R(z)=∏_{b=h+1}^{2h}(z−b²). The nullvector of Kₕ(q), normalized to
have coordinate sum1, is
\[
 r_b=\frac{Q_h(b^2)}{C_R'(b^2)}.
\]
Lagrange interpolation gives
\[
 Q_h(z)=C_R(z)\sum_b\frac{r_b}{z-b^2},\quad
 N_h(z)=C_R(z)\sum_b\frac{r_bf_b(q)}{z-b^2}.             \tag{5}
\]
These match the right nodes automatically. The left matching equations are
exactly Kₕ(q)r=0, proving (4). Conversely, any such pair gives a nonzero
Loewner nullvector. A denominator of degree<h−1 would give a nullpolynomial
in the fixed positive-definite leading Gram block, which is impossible.
The same observation excludes common numerator/denominator factors.
Uniqueness also follows by cross multiplication: the difference of two
candidate rational functions has numerator degree≤2h−2 and vanishes at2h
points. Monicity then fixes their common scale.

There is an equivalent entirely rational-data version. Put
\(P_h=N_h-qzQ_h\). It is the type(h,h−1) rational interpolant of βⱼ at
z=j², and
\[
 q=-[z^h]P_h.
\]
Thus the unknown ζ variable can be removed before solving any matrix problem.
This is the same unique interpolation problem, not an additional approximation
family.

## 3. Consecutive data have an unexpectedly short difference

The harmonic update Hⱼ₊₁³−Hⱼ³=1/(j+1)³ gives, for every real X,
\[
 \boxed{\frac{f_{j+1}(X)}{(j+1)^2}-\frac{f_j(X)}{j^2}
          =-\frac{2j+1}{2j^3(j+1)^3}.}                 \tag{6}
\]
Indeed the left side is
\(\tfrac12[j^{-2}-(j+1)^{-2}-j^{-3}-(j+1)^{-3}]\).
The X term cancels. This first-order rational difference is the relevant
extra structure beyond a generic Loewner matrix.

Substituting the interpolant (4) in (6) gives zero at j=1,…,2h−1 of
\[
 \begin{aligned}
 \mathcal R_h(j)={}&2j^3(j+1)N_h((j+1)^2)Q_h(j^2)\\
 &-2j(j+1)^3N_h(j^2)Q_h((j+1)^2)\\
 &+(2j+1)Q_h(j^2)Q_h((j+1)^2).
 \end{aligned}                                        \tag{7}
\]
The expression is antisymmetric under j↦−j−1. Hence it is (2j+1) times a
polynomial in the quadratic-lattice coordinate j(j+1). Its degree in that
coordinate is≤2h−1: the leading j⁴ʰ terms in the first two summands cancel.
The2h−1 known roots therefore determine it up to one scalar.

Let cₕ=[z^(h−1)]Nₕ, temporarily allowing cₕ=0. Expanding the first two terms at
infinity gives leading coefficient −4cₕ at j⁴ʰ⁻¹. Consequently
\[
 \boxed{\mathcal R_h(j)=
 -2c_h(2j+1)\prod_{r=1}^{2h-1}
                       [j(j+1)-r(r+1)].}              \tag{8}
\]
Evaluating at j=0 fixes the scalar independently:
\[
 \boxed{c_h=\frac{Q_h(0)Q_h(1)}{2(2h-1)!(2h)!}>0.}     \tag{9}
\]
Thus Nₕ actually has degree h−1, including h=1, and the Riccati remainder is
explicitly normalized. Equations (8)–(9) are all-h algebraic identities, not
conjectures inferred from the Apéry matches.

### A harmonic-free converse

Suppose monic Q of degree h−1 and N of degree≤h−1 satisfy (8), with
c=[z^(h−1)]N and
Q(j²)≠0 at j=1,…,2h. Define q=1+N(1)/Q(1). At the integer roots of the
right side, divide (7) by2j³(j+1)³Q(j²)Q((j+1)²). This gives precisely (6)
for N(j²)/Q(j²), with the initial value f₁(q)=q−1. Induction then gives (4)
at every node. Therefore this pair must be the original Schur pair and q must
be qₕ,₃. This characterizes the original root by a finite polynomial Riccati
problem with no harmonic sums and no zeta constants.

It is not yet a solution of that polynomial problem in a form that exposes
Apéry's recurrence. It does make the remaining task more explicit than an
unidentified equality of two optimized norms.

## 4. Examples and a new finite pattern, clearly distinguished

The first interpolation pairs are
\[
 Q_1=1,\quad N_1=1/4,
\]
\[
 Q_2=z+1/5,\qquad N_2=(z+290)/1200,
\]
\[
 Q_3=z^2+(191/149)z+20/149,
\quad N_3=(z^2+133171z+126948)/532824.
\]
They recover q by 1+N(1)/Q(1). The constants begin
\[
 c_1=1/4,\ c_2=1/1200,\ c_3=1/532824,\
 c_4=1/297243040.
\]

Besides the previously observed adjacent-Apéry q pattern, exact computations
through h=6 give another **conjectural all-h pattern**:
\[
 \boxed{c_h\stackrel?=
 \frac{72(h+1)}{h(5A_h-A_{h-1})^2}.}                   \tag{10}
\]
This is not proved by (9): the original polynomial values Qₕ(0),Qₕ(1) are
still unknown in closed form. Nor does (10) alone establish the q formula.
It may give a useful determinant/Christoffel target, since those two polynomial
values can be expressed as neighboring moment-determinant ratios.

## 5. What this bounded route did and did not resolve

The exact Loewner transformation did not collapse to an already known
tridiagonal determinant merely by applying the harmonic update. The quadratic
lattice instead produces the rational Riccati remainder (8). A proof of the
adjacent-Apéry identification could now proceed by constructing Nₕ,Qₕ solving
(8) from a finite hypergeometric/Christoffel complement, or by deriving an h
recurrence for its normalized solution. Neither step was completed here.

In particular the original D₂ₕ minimizer is not the ordinary positive Wilson
polynomial, and formally assigning a negative Wilson parameter would require
additional discrete masses/contour terms. No such continuation has been used.
The proved moving-pole identity remains a separate theorem.

The new half-parameter Wilson/Stieltjes connection communicated by the root
may help with a genuine Christoffel/Geronimus duality. Its meromorphic density
shift alone does not justify moving an integration contour through the integer
poles. That route has not been asserted as a proof here.

`check_rank_one_loewner_route.py` verifies all determinant factors, the exact
inverse formula, the original-minimizer interpolation, the barycentric
nullvector, and the Riccati identity through h=6. It separately flags the
finite Apéry and c patterns as such. Every tested identity passes. The
original all-h equality

    q_(h,3) = (5B_h−B_(h−1))/(5A_h−A_(h−1))

remains an explicitly unresolved bridge in this checkpoint.
