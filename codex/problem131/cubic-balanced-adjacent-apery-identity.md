# The cubic balanced diagonal is exactly an adjacent Apéry remainder

Status: accepted after a corrected Newton-conversion sign typo; see [the independent review](cubic-balanced-apery-independent-review.md). This note concerns only the cubic measure mu_3 and its ordinary orthogonal polynomials (the m=1 balanced construction). It makes no statement about analogous higher-zeta numerators and no historical novelty claim.

Write epsilon_n=A_n zeta(3)-B_n, with the accepted Apéry normalization A_0=1, A_1=5, B_0=0, B_1=6. Let p_n=p_(1,n) be the monic degree-n polynomial orthogonal for the original positive measure mu_3. Then, for every n>=0,
\[
 \boxed{\displaystyle
 \mu_3\!\left[\frac{p_n(t)/p_n(0)}{t+(n+1)^2}\right]
 =\frac{(n+1)^2}{2}I_{n+1,n}
 =\frac{(n+1)^2}{6}(\epsilon_{n+1}+\epsilon_n).}       \tag{1}
\]
Here I_nm is the accepted logarithmic Legendre-kernel integral. In particular this identifies the **full rational numerator**, not only its zeta coefficient. It explains why the cubic balanced diagonal recovers the Apéry exponential.

The inputs about the original measure are its already proved real resolvent and mass:
\[
 F(a):=\mu_3[(t+a^2)^{-1}]
 =a^2(\zeta(3)-H_a^{(3)})-\frac12+\frac1{2a},\quad a=1,2,\ldots,
 \qquad \mu_3[1]=\frac14.                              \tag{2}
\]
The [real kernel model](cubic-wilson-real-kernel-model.md) now proves the normalized polynomial formula (3) and functional identity (6) using only (2), the independently established logarithmic kernel, a maximum principle, and real moment comparison. It therefore supplies a fully real route to this theorem. The Rodrigues calculation in Section 1 is retained as an optional independent algebraic comparison with the earlier construction; it is no longer needed as an orthogonality input. The real-model proof refers back only to the kernel calculations in Section 2, which are independent of Section 1, so there is no circular dependency.

## 1. A short algebraic extraction of the ordinary Wilson polynomial

Define
\[
 Q_j(t)=\frac{\prod_{k=0}^{j-1}(t+k^2)}{(j!)^2},\qquad U(x^j)=Q_j(t),
\]
with Q_0=1. Put
\[
 S_n(x)=\frac{L_n(x)+L_{n+1}(x)}{2(1-x)}
       =\sum_{j=0}^n(-1)^j\binom nj\binom{n+j+1}j x^j.
\]
Thus S_n is half the Christoffel mode Q_(n+1) from the separate Legendre notation. To avoid confusing those two uses, below the Christoffel mode will be denoted C_a(x)=(L_a+L_(a-1))/(1-x), while Q_j always denotes the Newton product in t.

The normalized ordinary Wilson polynomial is exactly
\[
 \boxed{p_n(t)/p_n(0)=U S_n(t).}                         \tag{3}
\]
Here is a direct coefficient derivation from the accepted Rodrigues operators, rather than an unproved hypergeometric identification. For s>=0 set c=s/2 and
\[
 q_{s,j}(x)=(c+x)_j(c-x)_j.
\]
The m=1 raising operator B_s from that construction obeys
\[
 \boxed{B_s q_{s+1,j}=(s+j+1)^3q_{s,j}
                      -(2s+j+3)q_{s,j+1}.}             \tag{4}
\]
Indeed its two shifted terms factor q_(s,j); their difference is
q_(s,j)[(2s+j+3)x²+(c+1)²(2s+3j+1)]. Replacing x² by (c+j)²-q_(s,j+1)/q_(s,j) gives (4), because the remaining coefficient is (s+j+1)³.

Induction using (4), or direct comparison of its adjacent coefficients, yields
\[
 B_s\cdots B_{s+n-1}1
 =(s+1)_n^3\sum_{j=0}^n
 \frac{(-n)_j(n+2s+2)_j}{(s+1)_j^3j!}\,q_{s,j}.         \tag{5}
\]
For clarity, if f_(s,n,j) denotes the displayed coefficient, the induction is the finite identity
\[
 f_{s,n,j}=(s+j+1)^3f_{s+1,n-1,j}
                    -(2s+j+2)f_{s+1,n-1,j-1},
\]
with out-of-range coefficients zero; canceling adjacent Pochhammer factors proves it. At s=0, q_(0,j)(x)=(j!)² Q_j(-x²), and the value at x=0 is (n!)³. Dividing (5) by that value gives precisely (3).

Since the S_n are orthogonal for weight 1-x and form a polynomial basis, (3), ordinary orthogonality for mu_3, and the mass in (2) force the useful functional identity
\[
 \boxed{\mu_3[Uf]=\frac12\int_0^1(1-x)f(x)\,dx
        \quad\hbox{for every polynomial }f.}            \tag{6}
\]
Both sides equal 1/4 on S_0=1 and vanish on S_n for n>=1. In particular mu_3[Q_j]=1/[2(j+1)(j+2)]. This uses no Bernoulli summation.

Write C_(nj)=binom(n,j)binom(n+j,j). At the diagonal a=n+1, evaluating (3) already gives the denominator identity term by term:
\[
 \frac{p_n(-a^2)}{p_n(0)}
 =\sum_{j=0}^n
 \binom nj\binom{n+j+1}j\binom{n+1}j\binom{n+j}j
 =\sum_{j=0}^n C_{nj}C_{n+1,j}
 =\frac{A_n+A_{n+1}}6.                                 \tag{7}
\]
The last equality is the accepted adjacent balance. The rest of the proof identifies the numerator without extrapolating from (7).

## 2. The correct ordinary-coefficient operator and its explicit adjoint solution

Use the operators and functions from the accepted kernel proof:
\[
 H=-D[x^2(1-x)D]+x,\qquad
 k_m(x)=\int_0^1L_m(y)K(xy)dy,\qquad K(v)=-\frac{\log v}{1-v}.
\]
They satisfy (H+m(m+1))k_m=1. Define
\[
 T_0=H+xD,\qquad T_0^*=H-xD-1.
\]
The star denotes the formal adjoint for ordinary Lebesgue integration. The monomial formula is
\[
 T_0x^j=(j+1)^2x^{j+1}-j^2x^j,
 \qquad UT_0=tU.                                      \tag{8}
\]
For a>=1 put P=L_a+L_(a-1), D_a=L_a-L_(a-1), and
\[
 \widetilde k_a=\frac{k_a+k_{a-1}}2,
 \qquad h_a=\frac{a^2\widetilde k_a-1}2.
\]
The elementary Legendre ladder xP'=aD_a gives, by integration by parts in y,
\[
 (xD+1)(k_a+k_{a-1})=-a(k_a-k_{a-1}).
\]
Combining this with the two equations for k gives
\[
 \boxed{(T_0^*+a^2)\widetilde k_a=1,\qquad
        (T_0^*+a^2)h_a=(1-x)/2.}                       \tag{9}
\]
The finite Christoffel projection has coefficient 1/a on C_a and squared norm 2/a. Thus
\[
 \widetilde k_a(1)=\frac12\int_0^1(-\log y)C_a(y)dy=\frac1{a^2},
 \qquad h_a(1)=0.                                     \tag{10}
\]
The accepted boundary-row tail identity and (2) also give
\[
 \int_0^1 h_a(x)dx
 =\frac{a^2}{4}(I_{a0}+I_{a-1,0})-\frac12=F(a).        \tag{11}
\]
All these are explicit functions, not abstract inverse operators. At zero h_a has at most logarithmic growth and derivative O(1/x+|log x|); at one h_a and its derivative are bounded.

## 3. Polynomial division matches the actual original associated numerator

For fixed a and any polynomial f, the degree-raising formula (8) gives unique r and constant A with
\[
 f=(T_0+a^2)r+A.
\]
Applying U and dividing by t+a² gives
\[
 \mu_3\!\left[\frac{Uf}{t+a^2}\right]
 =\mu_3[Ur]+AF(a)
 =\frac12\int_0^1(1-x)r(x)dx+A\int_0^1h_a(x)dx.
\]
The Green identity has boundary bracket
\[
 [x^2(1-x)(r h_a'-r'h_a)+x r h_a]_0^1=0.
\]
It vanishes at one by (10), and at zero by the logarithmic bounds. Equation (9) therefore proves the exact transform
\[
 \boxed{\displaystyle
 \mu_3\!\left[\frac{Uf}{t+a^2}\right]=\int_0^1f(x)h_a(x)dx.}   \tag{12}
\]
In particular it applies to f=S_n and the original normalized p_n. This is why the rational numerator, including its endpoint correction, agrees.

## 4. At a=n+1 the correction cancels and leaves the mixed kernel

Now set a=n+1, so S_n=C_a/2=P/[2(1-x)]. Write
\[
 F=\iint P(x)P(y)K(xy)dxdy,\quad
 G=\iint D_a(x)D_a(y)K(xy)dxdy,\quad
 J=\iint C_a(x)P(y)K(xy)dxdy.
\]
Equation (12) gives
\[
 \mu_3\!\left[\frac{p_n/p_n(0)}{t+a^2}\right]
 =\frac{a^2}{8}J-\frac14\int_0^1C_a(x)dx.              \tag{13}
\]
The two Legendre ladders are
\[
 xP'=aD_a,\qquad (1-x)D_a'=-aP,
\]
with P(1)=0, D_a(0)=0, and D_a(1)=2(-1)^a. First integrate the identity D_a(y)=yP'(y)/a by parts in y in G. Since y partial_y K=x partial_x K, a second integration by parts in x gives
\[
 \begin{aligned}
 G&=-\frac1a\iint D_a(x)P(y)[K(xy)+x\partial_xK(xy)]dxdy\\
 &=-\frac{2(-1)^a}{a}\int_0^1P(y)K(y)dy
       +\frac1a\iint xD_a'(x)P(y)K(xy)dxdy\\
 &=\frac{4(-1)^{a-1}}{a^3}-(J-F).                    \tag{14}
 \end{aligned}
\]
The boundary integral is 2/a² by (10). All lower boundary terms vanish: yK(xy)->0 and xD_a(x)K(xy)->0. Bounded polynomials and the integrable logarithm justify the rearrangements, or they may first be done on truncated intervals.

Also, integrating D_a'=-aC_a gives
\[
 \int_0^1C_a(x)dx=\frac{2(-1)^{a-1}}a.
\]
Substituting (14) in (13) cancels its endpoint correction exactly. Finally symmetry of K gives F-G=4I_(a,a-1). Hence
\[
 \mu_3\!\left[\frac{p_n/p_n(0)}{t+a^2}\right]
 =\frac{a^2}{8}(F-G)=\frac{a^2}{2}I_{a,a-1}.
\]
The accepted adjacent identity I_aa+I_(a-1,a-1)=6I_(a,a-1), together with I_jj=2 epsilon_j, proves (1).

This includes n=0/a=1: both sides are zeta(3)-1. For every n the form is positive, since both diagonal errors in (1) are positive; the restriction to even n in the general moving-pole positivity argument is unnecessary in this cubic specialization. Its classical Apéry arithmetic and scale now follow from an exact identity, rather than from the observed common saddle.

Exact finite and symbolic checks are in [check_cubic_balanced_apery_bridge.py](../computations/check_cubic_balanced_apery_bridge.py).
