# Independent audit of the base-Gaussian derivative transfer

**Verdict: accepted.** The identities, absolute error theorem, rational
enclosures and exact function-pole order in
[the transfer note](higher-zeta-base-gauss-derivative-transfer.md)
are correct. The sign of the transferred error is not assumed.
No irrationality or useful reduced numerical denominator theorem
follows from these results alone.

The independent checks are included in
[check_higher_kernel_experiment_independent.py](../computations/check_higher_kernel_experiment_independent.py).

## 1. Kernel identities and the operator normalization

Integrating the density recurrence against an admissible test function
transfers the derivative to that function and gives the multiplier
\(2\Theta+2j+1\). Iteration from the base measure proves the
displayed Euler-operator identity. For the rational kernels in this
note, the endpoint terms vanish and all integrals are absolutely
convergent. Thus no undefined domain of a formal operator is needed.

The half-power identity
\[
 t^{-1/2}\partial_t^{m-1}\bigl(t^{m-1/2}f(t)\bigr)
 =\prod_{j=1}^{m-1}(\Theta+j+1/2)f(t)
 \]
has the right shifts. Multiplying by \(C_m\) gives exactly the
Euler factors divided by \(m(2m-1)!!\).
The identity with
\(2\partial_a^{2m-2}(a^2+t)^{-1}/(2m)!\) follows by the proved
raising recurrence, or by elementary differentiation. In particular
\[
 K_2(z,t)=\frac{3z-t}{6(z+t)^3}
 \]
has a genuinely negative tail. Replacing its measure by the base
measure does not preserve pointwise positivity of the higher kernel.

The factor \(2/(2m)!\) in the integrated kernel and the factor
\(1/(2m)!\) in the zeta approximation are consistent: the base
resolvent is \(2\zeta(3,a+1/2)\), and its \(2m-2\) derivative is
\((2m)!\zeta(2m+1,a+1/2)\).

## 2. The absolute derivative bound

The ordinary Gaussian error integrand has \(2N+1\) positive quadratic
denominators in \(a\), counted with multiplicity. Splitting each as
\((a+ib)(a-ib)\) gives exactly \(4N+2\) reciprocal linear factors.
The divided derivative of order \(q\) is a sum indexed by weak
compositions of \(q\) among those factors, so its term count is
\(\binom{4N+q+1}{q}\).

Each additional reciprocal factor has modulus at most \(1/a\).
The modulus of the original product is its original positive
quadratic denominator. The triangle inequality therefore proves
the stated bound after integrating. The same estimate with a
fixed positive lower bound on \(a\) supplies local domination for
interchanging derivatives and the integral. Complex numbers are
used only to factor a real quadratic; no complex-analytic theorem
is invoked.

Setting \(q=2m-2\) and dividing by \((2m)!\) gives the factor
\[
 \frac{q!}{(2m)!}=\frac1{(2m)(2m-1)}.
 \]
Thus every normalization in the absolute error bound (3) agrees.
Positivity of the original error does not determine the sign of its
even derivatives, and the note correctly asserts a two-sided enclosure.

## 3. Crude radii, rate and function denominator

The monic trial \(t^n\) and positive Gaussian roots give
\[
 E_n(a^2)\le M_{1,2n}/a^{4n+2}.
 \]
Combining this with the derivative bound gives the rational radius
in equation (4). Its polynomial exponent is
\((2m-2)+1+1/2-2m=-1/2\), and its exponential factor is
\((2/(\pi e))^{4n}\). The displayed identity equating it to the
positive-\(\sigma_m\) crude radius follows by factorial cancellation;
the independent checker verifies a separate 144-case grid as well.
Equality of these upper bounds does not identify the approximants.

An even \(a\)-derivative of order \(2m-2\) is a rational function
of \(z=a^2\). Its highest \(z\)-derivative term is
\[
 (4z)^{m-1}D_z^{2m-2}R_N(z).
 \]
At \(z=-x_i\), its leading pole coefficient is a nonzero multiple
of \((-4x_i)^{m-1}\lambda_i\). Positive nodes and positive weights
make this coefficient nonzero; all other terms have lower pole
order. Consequently the reduced function denominator is exactly
\(Q_N(z)^{2m-1}\). Rational specialization can still introduce
cancellation and does not inherit this exact denominator assertion.

## 4. Additional sharper radius from the already proved Apéry isometry

The base system has more structure than the monic-trial estimate uses.
At \(a=n+1/2\), the independently proved identities are
\[
 \nu(\Phi_n^2)=\frac1{2n+1},\qquad
 \Phi_n(-a^2)=A_n=\sum_k\binom nk^2\binom{n+k}k^2.
 \]
The normalization of the orthogonal polynomial cancels from the
ordinary error identity. Hence
\[
 E_n(a^2)\le\frac1{(2n+1)a^2A_n^2}
          =\frac4{(2n+1)^3A_n^2}.
 \]
Substitution into the audited derivative theorem proves the sharper
explicit rational enclosure radius
\[
 \boxed{
 |\zeta(2m+1)-Q_{n,m}|
 \le
 \frac{\binom{4n+2m-1}{2m-2}}{(2m)(2m-1)}
 (n+1/2)^{-2m+2}
 \frac4{(2n+1)^3A_n^2}.}
 \tag{A}
\]
This is sign-free and uses no Apéry recurrence. For fixed \(m\),
the binomial times the power of \(n+1/2\) is bounded, so (A) is
\(O_m(n^{-3}A_n^{-2})\).
The elementary Delannoy lower bound from the complete Padé proof,
\[
 A_n\ge\frac{\lambda^n}{(n+1)(2n+1)^2},
 \qquad \lambda=17+12\sqrt2,
 \]
then gives
\[
 |\zeta(2m+1)-Q_{n,m}|=O_m(n^3\lambda^{-2n}).
 \]
This improves the elementary exponential upper rate for the transfer
family. It does not improve its proved function-pole cost or settle
reduced numerical denominator growth. No historical novelty is
claimed for this corollary.

## 5. Reproduction

The independent Fraction-only checker verifies 30 exact transferred
function denominators by coprimality, 144 identities of crude radii,
and 72 half-power versus even-\(a\)-derivative kernel identities.
Its even derivatives are constructed with the real operator
\(D_a^2=2D_z+4zD_z^2\), independently of the author's symbolic
differentiation in \(a\).
The author's checker also passed all exact identities and its twelve
numerical enclosure corroborations.

The extra rate in (A) is proved from already audited norm and evaluation
identities, not inferred from finite computations. The analytic
transfer is established; its arithmetic effectiveness remains open.
