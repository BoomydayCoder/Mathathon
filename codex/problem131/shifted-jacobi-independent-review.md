# Independent audit of the shifted-Jacobi mixed Apéry family

Reviewed shifted-jacobi-mixed-apery-family.md and its checker.
**Accepted.** No gap was found in the rational identities, diagonal
closure, positive integral, residue formulas, denominator theorems, or
reciprocal asymptotic constants.

One denominator-comparison wording issue was corrected during review:
at \(q=1\), the normalizer (17) divides the familiar normalizer in (16),
so (16) is simpler, not necessarily sharper. Neither displayed
integrality theorem needed correction.

## 1. Mixed coefficients and diagonal elimination

The adjacent identity is valid after substituting
\(P_n/P_{n-1}=(t-n)/(t+n+a)\).
With \(X=t(t+a)\) and \(N=n(n+a)\), its constant terms agree because
\((2n+a)^2=4N+a^2\). This checks the extra \(a^2\) term in \(K_n\).

For the mixed identity, the two quadratic factors in the undifferenced
term have linear coefficients \(a-1\) and \(a+1\). Their product cancels
both the quartic and cubic terms of \(t^2(t+a)^2\). The remaining
quadratic coefficient is
\[
1+n(n+a+1)+(m+1)(m+a)
=\lambda_n+\lambda_m+a+1.
\]
This gives exactly the stated \(v_{nm}\), with the correct minus signs
from the signed family \(S_{nm}\).

The factors \(u_n/K_{n+1}\) and \(w_n/K_n\) simplify to the displayed
\(\alpha_n,\gamma_n\). Substituting both adjacent values in the
\(m=n\) mixed recurrence gives
\(\beta_n=v_{nn}-\alpha_n-\gamma_n\). All denominators are positive
for \(n\ge1,a\ge0\), including \(a=0\).

At zero shift, multiplication by six recovers the Apéry recurrence
with its exact coefficients. The characteristic limits are correct,
but the source appropriately does not use them alone to infer
asymptotics.

## 2. Boundary functionals, including a useful n=0 certificate

For rational functions with derivative tending to zero,
\[
\mathscr D(\Delta R)=R'(1).
\]
For \(n\ge1\), \(T_n\) has a double zero at one. In the mixed identity,
\(P_nP_m\) has a double zero when \(m\ge1\); at \(m=0\), the extra
factor \(t-1\) supplies the second zero. Their derivatives decay at
infinity. Thus summation legitimately kills the boundary terms
throughout the stated recurrence range.

The formal double-pole functional is translation invariant even when
shifted poles overlap: the coefficients add, and simple poles do not
turn into double poles. It therefore proves the same recurrence for
the growing coefficient without any irrationality assumption.

The note correctly states the diagonal recurrence only for \(n\ge1\).
At \(n=m=0\), a nonzero boundary explains the initial data. Define
\[
u_0=\frac{2(a+1)^2}{a+2},\qquad
v_{00}=\frac{2(a+1)^3}{a+2},\qquad
Q_{00}(t)=\frac{t-1}{t+a}.
\]
Then directly
\[
u_0S_{10}-v_{00}S_{00}=-2(a+1)\Delta Q_{00},
\qquad Q_{00}'(1)=\frac1{a+1}.
\]
Consequently the derivative-sum boundary is exactly \(-2\), not zero.
Using the adjacent identity at one, and putting
\[
\alpha_0=\frac{(a+1)^3}{(a+2)(a^2+3a+3)},
\]
gives
\[
F_1-(2a^2+6a+5)F_0=-\frac2{\alpha_0}=-V_1.
\]
This verifies the initial normalizations independently and specializes
to \(F_1-5F_0=-12\). It is an inhomogeneous boundary relation, not a
homogeneous extension using an undefined \(\gamma_0\).

## 3. Rodrigues and the positive cube integral

In the Mellin integration by parts, every boundary product at zero
has order \(x^{t+a}\), so the condition \(t+a>0\) is sufficient.
The beta integral gives precisely \((-1)^nP_n(t)\), with the
displayed normalization of the Jacobi polynomial.

There is no missing shift power in the two successive integrations
by parts. After the first one, the factor is
\[
\frac{x^{n+a}(1-x)^n\,y^{n+a}L_n^{(a)}(y)}
{(1+u)(xy+u)^{n+1}}.
\]
Under \(u=yv\), the Jacobian and denominator remove \(y^{n+1}\),
leaving \(y^aL_n^{(a)}(y)\) for the second Rodrigues step. That step
produces \(y^{n+a}(1-y)^nv^n\). Finally,
\(v=(1-c)/(yc)\) cancels its remaining Jacobian powers and yields
exactly the cube weight \((xy)^a\) in (10).

The absolute-integrability justification also works for noninteger
\(a\ge0\). After the first step, the extra ratio
\((xy/(xy+u))^n(1-x)^n\) is at most one, leaving a bounded multiple
of the original integrable Stieltjes kernel. The second integration
by parts is pointwise justified; its nonnegative output is bounded
after integration by the preceding absolute integral.
Thus no conditional Fubini step is hidden here.

At \(n=0\), the cube amplitude integrates to
\(2\zeta(3,a+1)\), so both positivity and the sharp-kernel upper bound
have the correct prefactor.

## 4. Residues and arithmetic

The double residues are the squares of the Jacobi coefficients.
Differentiating the regular part at a pole gives the simple residue
with sign
\[
2A_{n,j}\left(H_j-H_{n-j}
 -\sum_{h=1}^n(a+j+h)^{-1}\right).
\]
Their sum vanishes because \(P_n^2\) has no \(t^{-1}\) term at infinity.
Summing derivatives therefore cancels the common Hurwitz zeta-two
term and yields exactly the stated \(2A_n\zeta(3,a+1)-V_n\).

For integer shifts, the off-diagonal polynomial moment is
\[
\frac{\zeta(2,a+i+1)-\zeta(2,a+j+1)}{j-i}.
\]
Its harmonic denominators are at most \(a+n\), and its extra difference
factor is at most \(n\). Diagonal corrections use cubic harmonic
terms. Integral Jacobi coefficients therefore give the asserted
\(d_{a+n}^3\) normalizer without a spurious \(2n\) index.

For \(a=p/q\) in lowest terms, every numerator factor of
\(\binom{a+n+i}{n}\) is a unit at primes dividing \(q\).
At other primes, each block of \(\pi^k\) consecutive progression
indices contains a multiple of \(\pi^k\), canceling the factorial
denominator. The coefficient at \(i=0\) has no binomial multiplier
that could cancel the primes dividing \(q\), so it attains all of
the claimed denominator exponents. Thus \(D_{q,n}\) is genuinely
the least common denominator of the polynomial coefficients.

After clearing two such coefficients, diagonal corrections divide
\(L^3\), and off-diagonal corrections divide \(d_nL^2\).
Their common sufficient normalizer is
\(D_{q,n}^2L^2\operatorname{lcm}(L,d_n)\), as stated.
It clears \(A_n\) too. No least-denominator assertion for \(V_n\)
or the linear form is implied.

At \(q=1\), this normalizer can be smaller than \(d_{a+n}^3\);
for \(a=3,n=1\), they are respectively 64 and 1728.
The source now accurately calls the latter the simpler/familiar
integer-shift bound.

## 5. Concentration and reciprocal constants

The Beukers kernel extends continuously by zero at the cube boundary,
including the apparent singular corners. Its unique interior maximum
therefore has a strict gap off every neighborhood. The amplitude
\((1-c+xyc)^{-1}\) is integrable, while the extra factor \((xy)^a\)
is bounded and continuous. Standard localization, or its elementary
upper/lower neighborhood proof, gives
\[
\frac{F_n(a)}{F_n(0)}\to(\sqrt2-1)^{2a}.
\]

For the coefficient sum, an independent explicit exponent is useful:
\[
\log\!\left[\binom nj^2\binom{n+j}n^2\right]
=n f(j/n)+O(\log n)
\]
uniformly, with the endpoint convention \(0\log0=0\), where
\[
f(\theta)=2[(1+\theta)\log(1+\theta)
 -(1-\theta)\log(1-\theta)-2\theta\log\theta].
\]
Its derivative is \(2\log((1-\theta^2)/\theta^2)\); its second derivative
is \(-4/[\theta(1-\theta^2)]<0\).
It has the unique maximum at \(\theta=1/\sqrt2\), giving the required
exponential concentration including a gap at both endpoints.

The shifted/unshifted binomial ratio is the exact finite product
\(\prod_{h=1}^n(1+a/(j+h))\). Near that maximum, its logarithm tends
uniformly to \(a\log((1+\theta)/\theta)\); away from it the global
bound \(e^a n^a\) cannot overcome the exponential gap.
This proves the second ratio limit without needing a separate
unproved Gamma asymptotic.

The two ratio constants multiply to one. The base constants in
apery-laurent-reciprocal.md multiply to \(1/(2\sqrt2)\), so the
product asymptotic \(A_n(a)F_n(a)\sim1/(2\sqrt2\,n^3)\) follows with
exactly the stated normalization.

## 6. Independent checks and scope

The new check_shifted_jacobi_independent.py differentiates Rodrigues
expressions directly and computes actual partial-fraction residues,
independently of the source checker's finite-moment construction.
It passes:

- 40 Rodrigues, Mellin, residue, and exact coefficient-LCD cases;
- 48 recurrence checks for independently computed \(A,V\);
- the symbolic nonzero \(n=0\) boundary and all initial constants;
- nine actual adjacent identities and 27 actual mixed identities,
  including \(m=0\) and zero shift.

The shifts include \(0,1,3,1/2,1/4,2/3,7/6,7/10\), so the denominator
tests cover prime powers and multiple prime divisors of \(q\).
The source checker and its optional differentiated-series tests also pass.

For integer shifts the target differs from \(\zeta(3)\) by a rational
number. For general rational shifts, the sufficient normalizer still
needs a decay comparison before any irrationality conclusion.
The note correctly makes no blanket claim for all Hurwitz values.
