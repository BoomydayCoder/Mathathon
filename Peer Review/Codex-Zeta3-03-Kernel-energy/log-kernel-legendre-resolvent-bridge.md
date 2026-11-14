# A logarithmic kernel joins the rational square, Newton transport, and positive energy

Status: complete proof, accepted in [the independent audit](kernel-energy-independent-review.md), with exact checks in the accompanying script. This is a reconstruction connecting classical Apéry, Beukers, Legendre, and Wilson objects. No historical novelty is asserted. Its main contribution to the present exposition is that **the positive Padé estimate and the Hilbert isometry can both be obtained without Bernoulli polynomials or an imported spectral theorem**.

The positive logarithmic kernel already present in the rational-square proof satisfies a short differential equation. Its endpoint flux is exactly one. Those two facts supply the resolvent that the Newton construction needs, including its normalization. No self-adjoint operator closure is assumed: all functions used below are explicit integrals plus polynomials, and all energy estimates are ordinary integrations by parts.

## 1. Objects and the exact bridge

Write

\[
 L_n(x)=\frac1{n!}\frac{d^n}{dx^n}[x^n(1-x)^n]
       =\sum_{j=0}^n(-1)^j C_{nj}x^j,
 \qquad C_{nj}=\binom nj\binom{n+j}j,
\]

\[
 K(v)=-\frac{\log v}{1-v},\qquad
 k_m(x)=\int_0^1L_m(y)K(xy)\,dy,
 \qquad a_m=m+\tfrac12.
\]

The value at v=1 is defined continuously. Set

\[
 T=-\frac d{dx}\left[x^2(1-x)\frac d{dx}\right]+x-\tfrac14,
 \qquad J=-\frac d{dy}\left[y(1-y)\frac d{dy}\right].
\]

The central identities are

\[
 \boxed{(T+a_m^2)k_m=1,\qquad
        \int_0^1 k_m(x)\,dx=2\zeta(3,m+1).}                 \tag{1}
\]

Define the Newton transform on polynomials by

\[
 U(x^j)=Q_j(t),\qquad
 Q_j(t)=\frac{\prod_{r=0}^{j-1}(t+(r+1/2)^2)}{(j!)^2},
 \qquad \Phi_n=UL_n,
\]

and put

\[
 A_n=\Phi_n(-a_n^2)=\sum_{j=0}^n C_{nj}^2.
\]

Let

\[
 r_n(x)=U^{-1}\left(\frac{\Phi_n(t)-A_n}{t+a_n^2}\right)(x),
 \qquad g_n=r_n+A_n k_n.
\]

Here r_n is a polynomial of degree at most n-1, with r_0=0. Then

\[
 \boxed{(T+a_n^2)g_n=L_n,\qquad
 A_n I_n=\langle L_n,g_n\rangle
        =\langle(T+a_n^2)g_n,g_n\rangle>0,}                \tag{2}
\]

where

\[
 I_n=\iint_{(0,1)^2}L_n(x)L_n(y)K(xy)\,dx\,dy.
\]

The same elementary energy calculation proves

\[
 \boxed{0<I_n\le\frac{1}{(2n+1)a_n^2A_n}
             =\frac4{(2n+1)^3A_n}.}                       \tag{3}
\]

This is exactly the positive Wilson/Padé estimate for the classical rational-square linear form. Sections 2–5 prove it directly from the logarithmic kernel. Section 6 identifies the continuous measure and its Hilbert Gram matrix without Bernoulli moments. Section 7 records the Mellin and Beukers identities, so the connection does not rest on matching initial values or on an assumed recurrence.

## 2. The kernel equation and its boundary flux

A direct differentiation of (v-1)K(v)=log v gives

\[
 (v^2-v)K''(v)+(3v-1)K'(v)+K(v)=0.
\]

Indeed v[(v-1)K(v)]'=1; differentiate once more. Applying the two operators at v=xy therefore gives the pointwise identity

\[
 \boxed{(T_x+J_y+\tfrac14)K(xy)=0.}                       \tag{4}
\]

For fixed 0<x<1, integrate (4) against L_m(y). Its ordinary Legendre equation is JL_m=m(m+1)L_m. The Green boundary term is

\[
 \int_0^1 L_mJ_yK(xy)\,dy
 -\int_0^1K(xy)JL_m\,dy
 =[-L_m y(1-y)\partial_yK(xy)
       +y(1-y)L_m'K(xy)]_0^1=-1.                         \tag{5}
\]

At y=1 both contributions vanish. At y=0, K(xy)=-log(xy)+O_x(y|log y|), y partial_y K(xy) tends to -1, and L_m(0)=1. Thus the bracket tends to +1 at the lower endpoint. Substituting (5) in (4) proves the first part of (1), with source +1 and shift m(m+1)+1/4=a_m^2.

All differentiations under the y integral are valid locally in 0<x<1. For example, x derivatives of K(xy) have at worst a constant times 1+|log y| on compact x intervals. The logarithm is integrable at zero.

For later integrations by parts the following simple endpoint bounds suffice:

\[
 k_m(x)=O_m(1+|\log x|),\qquad
 k_m'(x)=O_m(x^{-1}+|\log x|)\quad(x\downarrow0),          \tag{6}
\]

while k_m and k_m' stay bounded as x approaches 1. At the latter endpoint the apparent singularity of K at 1 is removable, together with its derivatives. Bounds (6) follow directly by writing K(xy)=-(log x+log y)/(1-xy), first for x<=1/2. They ensure that all functions below lie in L^2, that their weighted derivative energies are finite, and that the boundary terms x^2(1-x)f'g and x(1-x)fg vanish whenever f,g are a polynomial plus a fixed k_m.

## 3. Two Legendre ladder identities fix the zeta-tail normalization

This calculation avoids Bernoulli moments and rational partial fractions. Put

\[
 F(y)=\int_0^1K(xy)\,dx=\frac1y\int_0^yK(v)\,dv,
 \qquad s_m=\int_0^1L_m(y)F(y)\,dy=\int_0^1k_m(x)\,dx.
\]

Thus (yF(y))'=K(y), and s_0=2zeta(3), by integrating the positive geometric series for K(xy).

For every m>=1 the elementary Legendre ladder identities are

\[
 y(L_m'+L_{m-1}')=m(L_m-L_{m-1}),
\]

\[
 (1-y)(L_m'-L_{m-1}')=-m(L_m+L_{m-1}).                 \tag{7}
\]

They follow from the usual Legendre derivative relations, or directly by comparing the displayed finite coefficients. Two integrations by parts give

\[
\begin{aligned}
 s_m-s_{m-1}
 &=-\frac1m\int_0^1(L_m+L_{m-1})K(y)\,dy\\
 &=-\frac1{m^2}\int_0^1(L_m'-L_{m-1}')\log y\,dy\\
 &=\frac1{m^2}\int_0^1\frac{L_m-L_{m-1}}y\,dy\\
 &=\frac1{m^3}[L_m+L_{m-1}]_0^1=-\frac2{m^3}.
\end{aligned}                                                       \tag{8}
\]

In the first step the endpoint at y=1 vanishes because L_m(1)+L_{m-1}(1)=0, and yF(y) tends to zero at y=0. In the logarithmic integration by parts, L_m-L_{m-1}=O(y) at zero. Every improper integral is absolutely convergent. Telescoping (8) proves s_m=2(zeta(3)-H_m^(3)), the second part of (1).

This identifies both the half-square sampling point and the cubic tail: the quadratic Legendre eigenvalue supplies a_m^2, and the two first-order ladder identities supply three factors of m in the tail difference.

## 4. Newton division explains the Apéry integer

The monomial action is

\[
 Tx^j=(j+1)^2x^{j+1}-(j+1/2)^2x^j.
\]

The Newton products consequently satisfy

\[
 \boxed{UT=tU,\qquad U^{-1}p(t)=p(T)1.}                   \tag{9}
\]

The second equality follows from the first and U1=1. These are finite polynomial identities. No inner-product or measure assertion is needed to obtain them.

Evaluation at -a_n^2 gives

\[
 Q_j(-a_n^2)=\prod_{r=0}^{j-1}\frac{(r-n)(r+n+1)}{(r+1)^2}
             =(-1)^jC_{nj} \quad(0<=j<=n).
\]

Thus A_n=sum C_{nj}^2. Its two binomial factors are the coefficients of L_n and the evaluation values of U(x^j). The integer is the unique remainder in the polynomial operator division

\[
 \boxed{L_n=(T+a_n^2)r_n+A_n.}                           \tag{10}
\]

Equivalently, define the finite evaluation functional E_n(p)=(Up)(-a_n^2). It obeys E_n(1)=1 and E_n((T+a_n^2)p)=0. Hence applying it to (10) forces the remainder A_n. This states exactly why the normalization is chosen; it is not inferred from a numerical approximation.

Combining (10) with (1) gives (T+a_n^2)g_n=L_n. Orthogonality gives <L_n,r_n>=0, so <L_n,g_n>=A_n<L_n,k_n>=A_n I_n.

## 5. Positive energy replaces the positive spectral measure

For every polynomial plus a fixed k_m, integration by parts using (6) proves

\[
 \begin{aligned}
 \langle Tf,f\rangle
 &=\int_0^1x^2(1-x)f'^2dx+\int_0^1(x-1/4)f^2dx\\
 &=\frac14\int_0^1(1-x)(2xf'+f)^2dx
                       +\frac14\int_0^1xf^2dx\ge0.
 \end{aligned}                                                   \tag{11}
\]

The final integral is strictly positive for every nonzero real f. Therefore

\[
 E_n:=\langle L_n,g_n\rangle
      =\langle(T+a_n^2)g_n,g_n\rangle
      \ge a_n^2\|g_n\|_2^2>0.
\]

The function g_n is nonzero because it solves an equation with nonzero right side L_n. Cauchy–Schwarz now gives

\[
 a_n^2\|g_n\|_2^2\le E_n\le\|L_n\|_2\|g_n\|_2,
 \qquad
 E_n\le\frac{\|L_n\|_2^2}{a_n^2}
       =\frac1{(2n+1)a_n^2}.
\]

Divide by A_n>0 to obtain (2)–(3). This uses neither an abstract inverse operator nor its spectral theorem.

For positivity alone an even simpler form is available:

\[
 A_nI_n=\int_0^1x^2(1-x)g_n'^2dx
                 +\int_0^1[x+n(n+1)]g_n^2dx>0.
\]

The critical-quarter square completion is used only to obtain the sharp convenient norm bound (3).

## 6. Recovering the continuous measure and the Hilbert isometry without Bernoulli polynomials

This section is optional for the irrationality proof but makes the connection to the earlier positive-measure argument exact.

Use the independently established real resolvent identity from [half-shift-elementary-resolvent-proof.md](../Codex-Zeta3-04-Legendre-Pade/half-shift-elementary-resolvent-proof.md): the positive measure

\[
 d\nu(t)=\pi^2\operatorname{sech}^2(\pi\sqrt t)
                    \tanh(\pi\sqrt t)\,dt
\]

has all polynomial moments and

\[
 \int\frac{d\nu(t)}{t+a_m^2}=2\zeta(3,m+1)=\langle1,k_m\rangle.
                                                               \tag{12}
\]

Only that resolvent calculation is imported here; its separate Bernoulli moment evaluation is unnecessary.

Write mu_j=<T^j1,1>. The monomial formula makes every mu_j rational. Repeatedly moving T onto a polynomial in the equation (T+a_m^2)k_m=1 gives, for each fixed d>=0,

\[
 \langle1,k_m\rangle
 =\sum_{j=0}^{d}\frac{(-1)^j\mu_j}{a_m^{2j+2}}
  +\frac{(-1)^{d+1}}{a_m^{2d+2}}
                         \langle T^{d+1}1,k_m\rangle.             \tag{13}
\]

Equation (11) applied to k_m and Cauchy–Schwarz give ||k_m||_2<=a_m^(-2). Thus the remainder in (13) is O_d(a_m^(-2d-4)). All integrations by parts are justified by (6).

The ordinary finite geometric identity gives the same asymptotic expansion for the left side of (12), with coefficients M_j=integral t^j dnu and remainder O_d(a_m^(-2d-4)), since M_(d+1) is finite. Comparing successively along a_m tending to infinity proves

\[
 \boxed{\nu(t^j)=\langle T^j1,1\rangle \quad(j>=0).}              \tag{14}
\]

Now T is symmetric on polynomials. For any spectral polynomials p,q, (9) and (14) give

\[
 \nu(pq)=\langle p(T)q(T)1,1\rangle
         =\langle p(T)1,q(T)1\rangle.
\]

Taking p=Q_j and q=Q_k proves immediately

\[
 \boxed{\nu(Q_jQ_k)=\langle x^j,x^k\rangle=\frac1{j+k+1}.}       \tag{15}
\]

This is the Hilbert Gram identity, obtained without a Bernoulli functional or a special-function orthogonality formula. Rationality of every moment also follows from (14).

Finally divide Phi_n(t) by t+a_m^2 and move T across k_m as in (13). Equations (12) and (14) yield the exact all-index identity

\[
 \boxed{
 \iint L_n(x)L_m(y)K(xy)\,dx\,dy
       =\int_0^\infty\frac{\Phi_n(t)}{t+a_m^2}\,d\nu(t).}         \tag{16}
\]

On the diagonal, orthogonality in (15) turns its right side into
A_n^(-1) integral Phi_n(t)^2/(t+a_n^2) dnu(t), exactly the previous positive Padé formula. Thus (2) is its explicit real differential-equation counterpart, not a new unrelated approximation.

## 7. The same kernel is the rational Mellin square and the Beukers cube

Let

\[
 P_n(s)=\frac{\prod_{j=1}^n(s-j)}{\prod_{j=0}^n(s+j)}.
\]

Rodrigues integration gives integral x^(s-1)L_n(x)dx=(-1)^n P_n(s). Expanding K(xy)=sum_(r>=0) -(log x+log y)(xy)^r and integrating termwise gives

\[
 \boxed{I_n=-\sum_{s=1}^{\infty}\frac d{ds}[P_n(s)^2].}          \tag{17}
\]

Absolute convergence follows by bounding the polynomial factors and the integrable positive logarithmic kernel, whose double integral is 2zeta(3).

The real Stieltjes identity

\[
 K(v)=\int_0^\infty\frac{du}{(1+u)(v+u)}
\]

and two Rodrigues integrations by parts give

\[
 \boxed{I_n=\iiint_{(0,1)^3}
 \frac{[x(1-x)y(1-y)z(1-z)]^n}{(1-z+xyz)^{n+1}}\,dx\,dy\,dz.}  \tag{18}
\]

The exact changes of variables and convergence details are already proved in [short-elementary-hybrid-proof.md](../Codex-Zeta3-02-Beukers-and-hybrid/short-elementary-hybrid-proof.md), Section 3. Identities (16)–(18) directly identify the rational, Wilson/Padé, differential-energy, and Beukers descriptions.

The usual partial fractions in (17) give I_n=2(A_n zeta(3)-B_n), with d_n^3 B_n integral. Bound (3), the Delannoy lower bound for A_n, and the elementary period-six lcm estimate then finish irrationality exactly as in the existing resolvent guide. Alternatively use (18) and its AM–GM decay bound. No new arithmetic assertion is needed.

## 8. What is simplified, and what remains classical

The kernel K is already forced by differentiating the rational Mellin square. Its second-order equation produces the differential operator, the unit endpoint source, and the half-square parameter. Monomial action produces the Newton products and the integer A_n. Positive energy supplies the Padé error estimate. The explicit continuous measure, if desired, is identified through its real resolvent; its full Hilbert Gram identity follows by polynomial symmetry and asymptotic moment comparison.

This removes the Bernoulli functional and moment-normalization computation from this route. It does not simplify the classical cubic-lcm arithmetic or claim a new irrationality result. The Sturm–Liouville family and Wilson polynomials are classical, as documented in [newton-hilbert-selfadjoint-classification.md](https://github.com/BoomydayCoder/Mathathon/blob/03cc75d694220d0fba888c637e573a29799e6d5b/codex/problem131/newton-hilbert-selfadjoint-classification.md). The particular bridge above has not received a separate priority search.

Reproduce the finite symbolic checks with `python codex/computations/check_log_kernel_legendre_bridge.py`. The all-index claims rest on the displayed proofs, not the finite checks.
