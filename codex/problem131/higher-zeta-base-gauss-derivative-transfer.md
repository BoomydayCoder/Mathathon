# Reuse the base zeta(3) Gaussian system by differentiation

Status: exact identities and an absolute error theorem proved below from the
root agent's proposed transfer. This is a companion construction to
[higher-order-positive-stieltjes-zeta.md](higher-order-positive-stieltjes-zeta.md).
No sign of the differentiated approximation error is assumed, and no useful
reduced-denominator growth theorem or irrationality result is asserted.

The advantage is concrete: every higher target uses the same base orthogonal
polynomials and their known rational Gaussian functions. The cost is a higher
function-level pole order, and positivity of the original remainder does not
by itself fix the sign of its a-derivatives.

## 1. Exact operator and the sign-changing base kernel

Let nu=sigma1 be the positive probability measure with transform
S(a²)=2 zeta(3,a+1/2). The higher-density recurrence gives, for admissible
smooth rational f,

\[
 \int f\sigma_m
 =\frac1{m(2m-1)!!}\int
        \prod_{j=1}^{m-1}(2\Theta+2j+1)f\,d\nu,
 \qquad\Theta=t\partial_t.
\]

Consequently the pure higher transform can also be written over nu as

\[
 2\zeta(2m+1,a+1/2)=\int K_m(a^2,t)\,d\nu(t),
\]

where, with C_m=2^(m-1)/[m(2m-1)!!],

\[
 \boxed{K_m(z,t)=C_m t^{-1/2}\partial_t^{m-1}
           \left[\frac{t^{m-1/2}}{(z+t)^m}\right]
 =\frac2{(2m)!}\partial_a^{2m-2}(a^2+t)^{-1},\quad z=a^2.}
 \tag{1}
\]

The half powers in the first expression cancel after differentiation; it is
an ordinary rational kernel. They are essential in the formula: replacing
t^(m-1/2) by t^(m-1) would give the wrong kernel. The operator identity
 t^(-1/2)partial_t^(m-1)t^(m-1/2)
 =product_(j=1)^(m-1)(Theta+j+1/2)
proves the first expression. The second follows either by the kernel raising
recurrence or by direct induction in m.

For example

\[
 K_2(z,t)=\frac{3z-t}{6(z+t)^3},
\]

which is negative for t>3z. Thus using nu instead of sigma_m reintroduces a
signed kernel, although its exact integral is the positive pure target.

## 2. A derivative error bound with no sign assumption

Let p_N be the monic orthogonal polynomial for nu, let R_N(z) be its ordinary
Gaussian rational function, and write

\[
 E_N(z)=S(z)-R_N(z)
 =\int\frac{p_N(t)^2}{(z+t)\prod_{i=1}^N(z+x_i)^2}\,d\nu(t)>0.
\]

For every integer q>=0 and a>0,

\[
 \boxed{\frac{|\partial_a^q E_N(a^2)|}{q!}
 \le\binom{4N+q+1}{q}\,a^{-q}E_N(a^2).}
 \tag{2}
\]

This is an elementary finite-product bound. Factor each positive quadratic
a²+b² as (a+ib)(a-ib). The integrand has 2N+1 quadratic factors, hence
4N+2 linear factors, counted with multiplicity. The q-th divided derivative
is a sum over all weak compositions of q among those factors. Each added
reciprocal factor has modulus at most a^(-1), because |a+ib|>=a. The number
of compositions is binom(4N+q+1,q). The product of the un-differentiated
absolute values is the original positive denominator. Taking absolute values
and integrating proves (2). Only complex-number factorization and the triangle
inequality are used; no complex-analytic estimate or contour argument is
required. The same bound locally dominates the derivatives and justifies
interchanging them with the integral.

The Hurwitz series gives

\[
 \partial_a^{2m-2}S(a^2)=(2m)!\zeta(2m+1,a+1/2).
\]

Thus the rational function

\[
 A_{N,m}(a)=\frac1{(2m)!}\partial_a^{2m-2}R_N(a^2)
\]

satisfies

\[
 \boxed{|\zeta(2m+1,a+1/2)-A_{N,m}(a)|
 \le\frac{\binom{4N+2m-1}{2m-2}}{(2m)(2m-1)}
       a^{-2m+2}E_N(a^2).}
 \tag{3}
\]

A priori this is an absolute bound, not a one-sided approximation theorem.

## 3. Rational convergence bounds at the growing sample point

Take N=n>=1 and a=n+1/2. Put

\[
 Q_{n,m}=H_n^{(2m+1)}+A_{n,m}(n+1/2)\in\mathbb Q.
\]

Let M_(1,2n)=(4n+1)|B_(4n)(1/2)| be the base measure's moment. The monic
trial t^n and positivity of all roots of p_n give
E_n(a²)<=M_(1,2n)/a^(4n+2). Therefore

\[
 \boxed{|\zeta(2m+1)-Q_{n,m}|\le B_{n,m},\qquad
 B_{n,m}=\frac{\binom{4n+2m-1}{2m-2}M_{1,2n}}
               {(2m)(2m-1)(n+1/2)^{4n+2m}}\in\mathbb Q_{>0}.}
 \tag{4}
\]

The pair Q_(n,m)-B_(n,m), Q_(n,m)+B_(n,m) is a rigorous rational enclosure.
For each fixed m, the same Stirling calculation as in the positive sigma_m
construction gives

\[
 B_{n,m}=O_m\left(n^{-1/2}(2/(\pi e))^{4n}\right).
\]

In fact the crude radius (4) equals exactly the crude upper-error bound (12)
in the positive sigma_m construction. Indeed the product formula for its
moments gives the finite identity

\[
 \frac{\binom{4n+2m-1}{2m-2}M_{1,2n}}{(2m)(2m-1)}
 =\frac12\binom{2n+m-1}{m-1}M_{m,2n}.
\]

The two approximants themselves are generally different. One keeps a positive
higher-order kernel and changes the Gaussian measure; the other keeps the
base Gaussian system and differentiates its rational function.

## 4. Exact function-level arithmetic cost

Write R_N=P_N/Q_N in reduced form, with Q_N(z)=product_i(z+x_i). Applying
an even derivative of order 2m-2 in a to R_N(a²) produces a rational function
of z=a² whose reduced denominator is exactly Q_N(z)^(2m-1). At each root
z=-x_i, the highest derivative term has nonzero coefficient proportional to
(4z)^(m-1) and the positive Gaussian weight. Since x_i>0, this coefficient
cannot vanish. All other derivative terms have lower pole orders. Thus each
pole has exact order 2m-1.

This larger pole order is a concrete arithmetic cost of reusing the base
system. It is not a theorem about reduced denominators after rational
specialization: numerical cancellation and denominators of polynomial
coefficients must still be analyzed. An explicit Apéry/Wilson description of
the base R_N may make that investigation feasible, but (1)--(4) alone do not
prove or disprove the needed irrationality balance.


## Reproduction

[check_base_gauss_derivative_transfer.py](../computations/check_base_gauss_derivative_transfer.py)
checks the two exact kernel formulas through m=5, 64 identities equating the
crude error bounds, and 12 rational Gaussian derivative systems with exact
reduced denominator Q^(2m-1). Numerical zeta comparisons corroborate their
rational enclosures. It also checks the generic Beta-density identity from
the companion note. All tests pass. The derivative-error theorem itself is
the finite-product argument in section 2, not an inference from numerical
signs. A separate independent proof audit is in progress.
