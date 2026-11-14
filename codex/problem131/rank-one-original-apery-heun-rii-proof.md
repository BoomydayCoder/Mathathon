# The original rank-one Schur root equals the adjacent-Apéry quotient

Status: complete proof, accepted in both `original-schur-apery-root-review.md`
and `rank-one-heun-rii-second-independent-review.md`. The argument below
concerns the ORIGINAL positive measure mu_3/D_(2h). It uses the proved boundary
Pearson identity and a varying-measure polynomial recurrence. It never invokes
ordinary orthogonality for Wilson parameters (1,1,1,-2h). No novelty claim is
made. The earlier moving-pole identity is not used to identify the root.

The theorem is

\[
 \boxed{q_{h,3}=\frac{5B_h-B_{h-1}}{5A_h-A_{h-1}}\qquad(h\ge1),}       \tag{T}
\]

where A_0=1,A_1=5,B_0=0,B_1=6 and both sequences obey

\[
 (r+1)^3u_{r+1}=(2r+1)(17r^2+17r+5)u_r-r^3u_{r-1}.
\]

The exact symbolic and finite positive-Gram checks are in
`../computations/check_rank_one_heun_rii_bridge.py`. They supplement the proof.

## 1. Original polynomials and the accepted inputs

Let p_h(t) be the monic degree-h-1 orthogonal polynomial for dmu_3(t)/D_(2h)(t),
and put Q_h(z)=(-1)^(h-1)p_h(-z). Thus Q_h is monic with strictly negative
zeros. Define

\[
 R_h(z)=Q_h(z)/Q_h(0),\qquad L_h=[z^{h-1}]R_h=1/Q_h(0)>0.
\]

The positive Schur construction gives the original rational q_h=q_(h,3) and
its error q_h-zeta(3)>0. We use the all-h rational interpolation consequence
proved in `rank-one-loewner-quadratic-lattice.md`: there is a polynomial S_h
of degree at most h such that

\[
 S_h(j^2)=\beta_jR_h(j^2)\quad(1\le j\le2h),\qquad
 [z^h]S_h=-q_h L_h,
\]
\[
 \beta_j=-j^2H_j^{(3)}-\frac12+\frac1{2j}.           \tag{1}
\]

For clarity, the construction there first gives N_h/Q_h interpolating
j²q_h+beta_j with deg N_h<=h-1; then S_h=(N_h-q_h zQ_h)/Q_h(0).
Its existence follows from the nullvector of the original cross-Cauchy Gram
matrix, not from the conjectured Apéry formula.

The other input is the elementary finite-difference boundary identity in
`rank-one-real-boundary-proof.md`. Only the degree range stated there is
needed below; the independent all-degree contour proof in
`rank-one-wilson-boundary-pearson.md` is optional. For even polynomials in x, set

\[
 \ell_M(f)=\int_0^\infty f(iy)\frac{d\mu_3(y^2)}{D_M(y^2)},\quad
 A_M(x)=\frac{(x+1)^3(x-M)}{2x(2x+1)},
\]
\[
 \mathscr D_Mf=A_M(x)[f(x+1)-f(x)]
                 +A_M(-x)[f(x-1)-f(x)],\qquad
 c_M=M/[2(M!)^2].
\]

It says

\[
 \ell_M(f\mathscr D_Mg-g\mathscr D_Mf)
       =-c_M[f(0)g(1)-f(1)g(0)].                  \tag{2}
\]

Every apparent operator pole cancels on even polynomials. The finite proof
of (2) holds whenever deg_(x²) f+deg_(x²) g<=M-1. It uses the original
rational pole values, detailed balance on the finite integer grid, and a
vanishing finite difference; it has been independently reviewed. The inner product ell_M(fg) on
even real polynomials is positive definite because f(iy) is real and the
original measure is positive.

## 2. A symmetric Heun operator supplies two endpoint values

Let n=h-1, M=2h, and let X denote multiplication by x². Put

\[
 \kappa_h=1-2(h-1)^2,\qquad
 \mathscr H_h=X\mathscr D_M+\mathscr D_MX-\mathscr D_M-\kappa_hX.
                                                               \tag{3}
\]

**This operator is symmetric for the original inner product.** Indeed applying
(2) to (Xf,g) and to (f,Xg), then adding, gives exactly the same boundary
defect for X D + D X as for D: the endpoint values of X are 0 and 1.
Subtracting D cancels that defect, while multiplication by X is symmetric.
For f,g of degree <=n in x², the pairs (Xf,g),(f,Xg) have degree sums
at most 2n+1=M-1, exactly within the real boundary lemma's range. Hence
ell_M(f H_h g)=ell_M(g H_h f) on the entire space used below.

The leading coefficient of D_M on x^(2k) is
lambda_k=k(k-2n). Consequently the coefficient raising degree k to k+1
in x² under H_h is

\[
 \lambda_k+\lambda_{k+1}-\kappa_h
                  =2(k-n)(k-n+1).                \tag{4}
\]

It vanishes BOTH at k=n-1 and k=n. Thus H_h preserves the polynomial
spaces of degrees <=n-1 and <=n. Symmetry and the defining orthogonality
of Q_h now imply

\[
 \boxed{\mathscr H_h[Q_h(x^2)]=\eta_h Q_h(x^2),\qquad
             \eta_h=1-2h-2h^3(h+1).}             \tag{5}
\]

Here is why the eigenvalue is fixed before knowing any coefficients of Q_h.
The coefficient in degree n of H_h x^(2n) is eta_h by the direct computation
below. A lower term of degree n-1 cannot contribute in degree n because
of the second zero in (4). Thus H_h Q_h-eta_h Q_h has degree <n. It is
orthogonal to every degree <n by symmetry and preservation of the lower
space, so positivity forces it to be zero. At h=1 the lower space is empty
and direct application to the constant polynomial gives eta_1=-5.

For the computation and the next endpoint, use

\[
 \phi_k(x)=\prod_{j=1}^k(x^2-j^2),\quad
 \lambda_k=k(k+2-M),\quad \mu_k=k(k+1)^2(k-M).
\]

The already proved Newton action is D_M phi_k=lambda_k phi_k+mu_k phi_(k-1),
while X phi_k=phi_(k+1)+(k+1)²phi_k. Therefore H_h is tridiagonal:

\[
 H_h\phi_k=r_k\phi_{k+1}+b_k\phi_k+d_k\phi_{k-1},
\]
\[
 r_k=\lambda_k+\lambda_{k+1}-\kappa_h,\quad
 b_k=\lambda_k[2(k+1)^2-1]+\mu_k+\mu_{k+1}-\kappa_h(k+1)^2,
\quad d_k=2k(k+1)\mu_k.                            \tag{6}
\]

In particular b_n=eta_h, proving the coefficient calculation in (5).

Evaluation at zero is especially simple:

\[
 (H_h f)(0)=-2M f(1)+(1-M)f(0).
\]

One may obtain it by expanding A_M at zero, or by applying D_M to
(x²-1)f. Equation (5) therefore gives

\[
 \boxed{R_h(1)=y_h:=\frac{h^2(h+1)}2.}             \tag{7}
\]

For a second value, evaluate the SAME eigen-equation at x=1. The terms
X D-D cancel there, and A_M(-1)=0, so

\[
 \eta_h R_h(1)=A_M(1)[4R_h(4)-R_h(1)]-\kappa_h R_h(1),
 \qquad A_M(1)=4(1-M)/3.
\]

This gives, including h=1,

\[
 \boxed{R_h(4)=z_h:=
 \frac{h^2(h+1)^2(3h^3+3h-2)}{16(2h-1)}.}          \tag{8}
\]

Both endpoint formulas are now all-h theorems, not fitted patterns. This
also independently explains the earlier Newton-coefficient pattern:
if R_h(x²)=sum e_k phi_k(x), then e_1/e_0=[z_h/y_h-1]/3,
which is (h-1)(h³+2h²+3h-2)/[8(2h-1)].

## 3. The three-term recurrence exists for the original varying measures

For every h>=2 there are real coefficients a_h,b_h,c_h such that

\[
 \boxed{R_{h+1}(z)=(a_hz+b_h)R_h(z)
       -c_hF_h(z)R_{h-1}(z),\quad
 F_h(z)=[z-(2h-1)^2][z-(2h)^2].}                  \tag{9}
\]

This standard-looking RII relation needs a proof because the measures vary.
Choose the linear polynomial a_hz+b_h to agree with R_(h+1)/R_h at the
TWO OLD nodes (2h-1)²,(2h)². These are valid values because R_h has no
positive zeros. Then

\[
 W(z)=\frac{R_{h+1}(z)-(a_hz+b_h)R_h(z)}{F_h(z)}
\]

is a polynomial of degree <=h-2. For k<=h-3,

\[
 \ell_{2h-2}(W(x^2)x^{2k})
 =\ell_{2h}([R_{h+1}(x^2)-(a_hx^2+b_h)R_h(x^2)]x^{2k})=0.
\]

The R_h term vanishes by its original orthogonality since k+1<=h-2.
For R_(h+1), convert ell_(2h) to ell_(2h+2) by multiplying by the TWO NEW
factors [(2h+1)²-x²][(2h+2)²-x²]; its test degree is k+2<=h-1, so it
vanishes by the next original orthogonality. Thus W is proportional to
R_(h-1), proving (9). For h=2 the orthogonality conditions are empty and
W is simply a constant, so the argument still applies. No h=1 recurrence
or undefined R_0 is asserted.

The three values z=0,1,4 now determine a_h,b_h,c_h. Write

\[
 d(h)=12h^2-15h+5,\qquad
 P(h)=204h^5+51h^4-134h^3-43h^2+30h+12.
\]

Solving the three linear equations obtained from (7)-(9) gives

\[
 \boxed{\begin{aligned}
 a_h&=\frac{P(h)}{h^2(h+1)^2(2h+1)d(h)},\\
 b_h&=\frac{48h^4-12h^3-34h^2+3h+7}
                {(h+1)(2h+1)d(h)},\\
 c_h&=\frac{(h-1)d(h+1)}
       {4h^2(h+1)(2h-1)(2h+1)d(h)}.
 \end{aligned}}                                                    \tag{10}
\]

For direct checking, those equations are

\[
 1=b_h-c_hF_h(0),\quad
 y_{h+1}=(a_h+b_h)y_h-c_hF_h(1)y_{h-1},
\]
\[
 z_{h+1}=(4a_h+b_h)z_h-c_hF_h(4)z_{h-1}.
\]

With unknowns ordered (a_h,b_h,c_h), their coefficient determinant is
3h^6(h+1)^4 d(h)/8, strictly positive for h>=2. Substitution or exact
elimination yields (10). The polynomial d(h) is strictly positive for real h
(its discriminant is -15), so no displayed denominator vanishes in our range.
This computes the original recurrence using only proven endpoint values.

## 4. The interpolation numerator has exactly the same recurrence

With S_h defined in (1), put

\[
 T(z)=S_{h+1}(z)-(a_hz+b_h)S_h(z)+c_hF_h(z)S_{h-1}(z).
\]

For the first 2h-2 old nodes, all three interpolation identities and (9)
give T(j²)=0. At the remaining two old nodes (2h-1)²,(2h)², F_h vanishes,
and the two remaining interpolation identities and (9) again give zero.
Thus T has 2h distinct zeros and degree at most h+1. For EVERY h>=2,
h+1<2h, so T=0. In particular the borderline h=2 case has four zeros
and degree at most three. This proves

\[
 S_{h+1}=(a_hz+b_h)S_h-c_hF_hS_{h-1}.              \tag{11}
\]

Let Y_h=q_hL_h. Comparing leading coefficients in (9) and (11) gives

\[
 \boxed{V_{h+1}=a_hV_h-c_hV_{h-1}\quad(h\ge2),
             \qquad V=L\ \hbox{or}\ Y.}           \tag{12}
\]

The minus sign in the leading coefficient of S_h cancels on both sides.
There is no forcing term and no missing recurrence at h=1.

## 5. Identify the scalar recurrence, including both initial values

Set s_h=h²/[12(2h)!]. If V_h=s_h C_h, equation (12) is

\[
 C_{h+1}=\alpha_hC_h-\gamma_hC_{h-1},\quad
 \alpha_h=\frac{2P(h)}{(h+1)^3d(h)},\qquad
 \gamma_h=\frac{(h-1)^3d(h+1)}{h(h+1)^2d(h)}.       \tag{13}
\]

For ANY solution u of the Apéry recurrence, the adjacent combination
C_h=5u_h-u_(h-1) obeys (13). This is a two-step algebraic elimination:
write u_(h+1) using the recurrence at h and u_(h-2) using the recurrence
at h-1, and compare the coefficients of u_h,u_(h-1). The determinant
that occurs is

\[
 25(h-1)^3-5(2h-1)[17(h-1)^2+17(h-1)+5]+h^3
                       =-12h\,d(h),
\]

so the elimination is valid for h>=2. It gives exactly the alpha and gamma
in (13). This symbolic identity is also checked independently in the script.

The initial original polynomials and roots are

\[
 R_1=1,\quad R_2=1+5z,\qquad
 q_1=5/4,\quad q_2=577/480.
\]

Hence L_1=1,L_2=5,Y_1=5/4,Y_2=577/96. On the Apéry side,

\[
 5A_1-A_0=24,\quad5A_2-A_1=360,
\]
\[
 5B_1-B_0=30,\quad5B_2-B_1=1731/4,
 \qquad s_1=1/24,\ s_2=1/72.
\]

Thus both initial pairs agree. Uniqueness in the scalar recurrence (12)
proves for all h>=1

\[
 \boxed{L_h=s_h(5A_h-A_{h-1}),\qquad
 Y_h=s_h(5B_h-B_{h-1}).}                          \tag{14}
\]

Since L_h>0, their ratio proves theorem (T).

## 6. Consequences and the precise former gap

In addition to the original quotient identity, the proof establishes

\[
 \boxed{Q_h(0)=\frac{12(2h)!}{h^2(5A_h-A_{h-1})},\qquad
        Q_h(1)/Q_h(0)=h^2(h+1)/2.}
\]

Together with the already proved Riccati normalization in the Loewner note,
this proves its former conjecture for the leading coefficient of N_h:

\[
 [z^{h-1}]N_h=\frac{72(h+1)}{h(5A_h-A_{h-1})^2}.
\]

The previously established moving-pole Wilson family has the SAME quotient.
Therefore equality of its optimized norm with the original D_(2h) norm
follows as a corollary of (T), rather than as an assumed norm duality.

The proof mechanism is: retain the actual boundary defect, cancel it with
a symmetric Heun operator, exploit its two consecutive degree cutoffs, and
use the resulting endpoint values to identify the varying-measure recurrence.
The rational interpolation numerator then inherits that recurrence by an
exact polynomial root count. Finite matching alone is not used at any stage.

Independent review should particularly verify: symmetry in (3), the double
zero in (4), positivity in the eigenpolynomial step, the OLD/NEW denominator
factors in the recurrence proof, h=1 and h=2 boundaries, the numerator degree
count, and the factorial rescaling in (13).
