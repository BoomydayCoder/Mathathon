# Independent audit: positive higher-order Stieltjes transforms for odd zeta

**Verdict: accepted.** The all-order positive density, its transform,
rational moments, differentiated Gaussian error, explicit convergence
rate and function-denominator statement in
[the construction](higher-order-positive-stieltjes-zeta.md) are correct.
No correction is required. The result gives convergent rational lower
bounds to each fixed odd zeta value; it supplies no irrationality theorem
or sufficient bound for their reduced numerical denominators.

Independent checks are in
[check_higher_order_independent.py](../computations/check_higher_order_independent.py).
They use only standard-library exact arithmetic, including rational
certificates for several actual lower bounds and errors.

## 1. Complete monotonicity and density normalization

Differentiating \(\log\cosh(\pi\sqrt t)\) and using the proved real
tanh partial fractions gives
\[
 (\log\cosh(\pi\sqrt t))'
 =\sum_{n\ge0}(t+(n+1/2)^2)^{-1}.
\]
The apparent square root at zero is removable. The series and each
of its derivative series converge locally uniformly for \(t\ge0\).
Integration from zero proves the displayed product for \(g(t)\);
the exponent is \(-2\), not \(-1\).

The definitions of \(H_r\) give exactly
\(g'=-H_1g\) and \(H_r'=-H_{r+1}\).
The Bell-polynomial recurrence therefore has positive signs in both
terms. It preserves nonnegative integer coefficients and the monomial
\(H_1^m\) with coefficient one. Because \(g,H_1>0\), it proves
\((-1)^mg^{(m)}>0\) everywhere on \(t\ge0\). No representation
theorem for completely monotone functions is being assumed.

At zero, \(g\) is analytic, so
\(\sigma_m(t)=O(t^{m-1/2})\). At infinity each \(t\)-derivative
of \(g\) introduces at most the expected \(t^{-1/2}\) leading factor,
with lower inverse powers accompanying it. Thus
\[
 \sigma_m(t)=O_m(t^{(m-1)/2}e^{-2\pi\sqrt t}).
\]
This also supplies the analogous bounds for its fixed derivatives.
All polynomial moments and the boundary estimates used below follow.

For \(m=1\), the coefficient \(C_1=1\) gives exactly
\(\sigma_1=\pi^2\tanh(\pi\sqrt t)\operatorname{sech}^2(\pi\sqrt t)\),
the previously proved probability density in the \(t\) variable.
There is no extra \(2\sqrt t\) Jacobian in this formula; that factor
appears only when subsequently changing back to \(y=\sqrt t\).

## 2. Raising the kernel and the exact moments

The ratio \(C_{m+1}/C_m=2m/((m+1)(2m+1))\) is correct.
Differentiating \(t^{m-1/2}g^{(m)}\) cancels its non-derivative term
in \((2m-1)\sigma_m-2t\sigma_m'\), leaving precisely the next
density with the stated prefactor.

After integration by parts, the new kernel is
\[
 \frac{m}{(m+1)(2m+1)}
 \left[-(a^2+t)^{-m-1}
       +2(m+1)a^2(a^2+t)^{-m-2}\right].
\]
Direct differentiation of \((a^2+t)^{-m}\) twice in \(a\)
gives the same expression after division by
\((2m+1)(2m+2)\). The boundary term vanishes at zero and infinity
by the density bounds. Differentiation under the integral is locally
dominated for \(a>0\).

The Hurwitz series has the identical raising rule. Induction from
the independently proved \(m=1\) resolvent therefore proves
\[
 \int_0^\infty\frac{\sigma_m(t)}{(a^2+t)^m}\,dt
       =2\zeta(2m+1,a+1/2)
\]
for every integer \(m\ge1\). The leading large-\(a\) normalization
also agrees: both sides have leading coefficient \(1/m\) times
\(a^{-2m}\).

The moment integration by parts gives the multiplier
\(m(2m+2k+1)/((m+1)(2m+1))\). Its iteration proves both
closed forms of equation (6). The factor \(1/m\) in the total mass
is essential and has been retained in all subsequent formulas.
The displayed elementary check of \(\sigma_2>0\) is also correct:
the bracket vanishes at zero and its derivative is
\(2\tanh^2u+6u\operatorname{sech}^2u\tanh u>0\).

## 3. Gaussian quadrature and the alternating derivative of its error

The positive density gives positive definite finite Gram forms.
Orthogonality then gives a unique rational monic polynomial \(p_N\).
The standard sign-change proof places all its roots simply in
\((0,\infty)\). Polynomial division by \(p_N\) and Lagrange
interpolation prove exactness through degree \(2N-1\), and squared
Lagrange polynomials prove that every quadrature weight is positive.
These are finite-dimensional arguments; no moment-determinacy theorem
or spectral realization is required.

The ordinary Gaussian function is rational even though its nodes and
weights need not be rational separately. Its numerator is the integral
of a polynomial difference quotient, hence uses only rational moments
and rational polynomial coefficients.

In this section the measure \(\sigma_m(t)\,dt\) is held fixed while
the kernel is differentiated. In particular \(S_m(z)\), its ordinary
resolvent, is not being identified with the earlier signed transform
\(G_{2m+1}\). The pure zeta value arises after \(m-1\) kernel derivatives.

Polynomial division and orthogonality prove the exact ordinary error
\[
 S_m(z)-R_N(z)
 =\int_0^\infty
   \frac{p_N(t)^2\sigma_m(t)}
        {(z+t)\prod_{i=1}^N(z+x_i)^2}\,dt.
\]
For fixed \(t\), this is a positive constant times \(2N+1\) positive
reciprocal linear factors, counted with multiplicity. Its alternating
divided derivative of order \(m-1\) is therefore strictly positive.
Differentiation under the integral is justified on every compact
positive \(z\) interval by its polynomial moment bounds.
This proves the strict Gaussian lower bound for the order-\(m\) target.

The number of terms in the divided product rule is the number of weak
compositions of \(m-1\) among \(2N+1\) slots:
\[
 \binom{2N+m-1}{m-1}.
\]
Every extra reciprocal factor is at most \(1/z\), and the original
factor \(1/(z+t)\) is also at most \(1/z\). Consequently the bound
in equation (10), including \(z^{-m}\) and its binomial coefficient,
is correct. Counting repeated node factors as separate slots gives
the same multiplicities as differentiating their squares directly.

## 4. Rational lower bounds and the rate

At \(a=n+1/2\), the Hurwitz tail is
\(\zeta(2m+1)-H_n^{(2m+1)}\). Division of the quadrature error by two
therefore gives exactly the claimed lower approximation \(L_{n,m}\).
The original version did not assert monotonicity in \(n\).
The later strengthening in
[nested-higher-zeta-gauss-radau.md](nested-higher-zeta-gauss-radau.md)
now proves strict monotonicity from the exact associated-polynomial
Wronskian and complete monotonicity of the error; see the addendum below.

The orthogonal monic polynomial minimizes the squared norm among monic
degree-\(n\) polynomials: the difference from the trial \(t^n\) has
lower degree and the cross term vanishes. Thus \(H_n\le M_{m,2n}\).
Positive roots give
\(|p_n(-z)|=\prod_i(z+x_i)\ge z^n\).
These facts prove the completely explicit rational bound (12).

For fixed \(m\), the Bernoulli estimate supplies a constant times
\((4n)!/(2\pi)^{4n}\), multiplied by a polynomial of degree \(m\).
The binomial factor has degree \(m-1\). Stirling's square-root factor
contributes \(n^{1/2}\), and the remaining denominator contributes
\(n^{2m}\). Hence the total polynomial exponent is
\[
 (m-1)+m+\tfrac12-2m=-\tfrac12.
\]
The exponential factor is exactly \((2/(\pi e))^{4n}\).
Replacing \(n+1/2\) by \(n\) in an upper bound is valid and changes
no exponent; retaining it gives only a further bounded constant factor.
The even zeta factor is uniformly bounded for \(4n\ge4\).
Thus the stated \(O_m(n^{-1/2}(2/(\pi e))^{4n})\) bound is correct.
The estimate \((2/(\pi e))^4<1/256\) follows from the stated strict
elementary inequalities. This rate is for each fixed \(m\), not a
uniform estimate as \(m\) increases with \(n\).

## 5. The denominator statement and an explanatory density identity

At every Gaussian node, \(R_{N,m}\) has a pole of exact order \(m\)
with nonzero leading coefficient \(\lambda_i>0\).
The other nodes are distinct, so none can cancel that pole.
Consequently its reduced function denominator over \(\mathbb Q[z]\)
is exactly \(Q_N(z)^m\), after monic normalization.
This does not assert anything equally exact about the denominator
of its value at a rational point. Numerical cancellations and
coefficient denominators still have to be analyzed separately.

There is a useful direct connection with the earlier signed density.
Let \(r=m-1\) and write that density in the \(t\) variable as
\(\rho_{2m+1}(t)=v_{2m+1}(\sqrt t)/(2\sqrt t)\).
Then
\[
 \boxed{\rho_{2m+1}(t)
       =\frac{m(-1)^r}{r!}\,t^r\sigma_m^{(r)}(t).}
 \tag{A}
\]
Indeed
\[
 \partial_t^r\!\left(\frac{t^r}{z+t}\right)
      =\frac{r!z^r}{(z+t)^{r+1}}.
\]
Integrating by parts \(r\) times proves that the right side of (A)
has ordinary resolvent
\[
 m z^r\int_0^\infty\frac{\sigma_m(t)}{(z+t)^m}\,dt
 =2m z^{m-1}\zeta(2m+1,\sqrt z+1/2)
 =G_{2m+1}(z).
\]
The endpoint terms vanish: at zero the powers remaining from
\(\sigma_m^{(j)}\) and \(t^r/(z+t)\) have positive total exponent;
at infinity the density derivatives decay exponentially in \(\sqrt t\).
Uniqueness of the finite signed Stieltjes transform identifies the
density with the previously proved \(\rho_{2m+1}\).
Equivalently, their moments satisfy
\[
 M^{\rm signed}_{2m+1,k}
       =m\binom{k+m-1}{k}M_{m,k},
\]
in agreement with both exact moment formulas.
This explains why a positive higher-order representation can coexist
with a signed ordinary one: lowering the kernel order differentiates
the density and destroys positivity.

## 6. Independent checks and conclusion

The independent checker passes:

- 180 rational moment identities and 12 density normalizations;
- the positive Bell-polynomial recurrence through order 12;
- 30 exact Gaussian systems, with rational Sturm counts verifying
  all positive roots and exact coprimality verifying the \(Q_N^m\)
  denominator;
- all refined-versus-crude error comparisons and the two stated
  initial rational approximations;
- six strict lower-bound and error certificates obtained from
  rational upper and lower bounds on the zeta series itself.

The author's checker also passed, including its six direct positive
resolvent quadratures. Finite checks do not replace the all-order
complete-monotonicity, integration-by-parts and Gaussian-error proofs.

The pure-all-odd analytic construction is established. The missing
irrationality step is arithmetic control of the reduced numerical
denominators; neither the function denominator nor the exponential
error bound supplies that control by itself.

## Addendum: strict nesting and Gauss–Radau upper bounds

The root agent proposed the stronger Wronskian/complete-monotonicity
argument after the initial audit. I checked the constants and supplied
the proof in the linked supplement. For monic polynomials and their
associated numerators, the identity is exactly
\[
 P_{n+1}Q_n-P_nQ_{n+1}=h_n.
\]
The initial constant is \(h_0=M_0\), so it applies without changing
normalization to the mass-\(1/m\) measure and to \(t\sigma_m(t)\,dt\).
The Gaussian increments and their alternating derivatives are strictly
positive. The error is strictly decreasing in z, so the growing
half-integer sample preserves strict improvement.

For the upper rule, \(R_n^+=(M_0-\widetilde R_n)/z\) has error
\(\widetilde E_n/z\) and positive endpoint weight
\[
 M_0-\widetilde R_n(0)
   =\frac{\int\widetilde p_n^2\,d\mu}{\widetilde Q_n(0)^2}>0.
\]
This is actual positive Gauss–Radau quadrature, not a formal signed
combination. Its derivative-error product has \(2n+2\) factors,
giving the upper-bound coefficient \(\binom{2n+m}{m-1}\).
The new upper radius has the same fixed-m rate as the lower one.

The independent exact checker
[check_nested_higher_zeta_bounds.py](../computations/check_nested_higher_zeta_bounds.py)
passes 80 Wronskian identities for both measures and 45 rational
lower/upper pairs, with positive endpoint weights, strict nesting and
both rational error radii. Five actual enclosures and their refined
errors are independently certified using rational zeta-series bounds.
The new supplement subsequently passed a
[further independent review](nested-higher-zeta-radau-independent-review.md);
the previous core audit remains unchanged.
