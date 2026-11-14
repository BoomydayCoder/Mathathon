# Independent review: deformed Laurent inverse periods

Status: independently checked proof; no substantive error found in
[deformed-laurent-period-reduction.md](deformed-laurent-period-reduction.md).
The all-n reduction, its differential ladder, and the unique positive
coefficient-cancellation parameter are proved. This review adds an explicit
real continuation to every q>0. No new irrationality or historical novelty
claim follows.

## 1. Mellin calculation, signs, and convergence

For the positive cube integral I_n(q), integration in q gives
B(s,n+1-s)(abc)^(-s)(1-c)^(s-n-1). The remaining a,b factors each integrate
to B(n+1-s,n+1), and the c factor to B(n+1-s,s). Thus, initially by Tonelli
for 0<s<n+1 and then by absolute convergence throughout that complex strip,

\[
 M_n(s)=\frac{\Gamma(s)^2\Gamma(n+1-s)^4}
                 {\Gamma(2n+2-s)^2}
       =\pi^2\csc^2(\pi s)r_n(-s),
\quad
 r_n(z)=\frac{\prod_{j=1}^n(z+j)^2}
                {\prod_{j=n+1}^{2n+1}(z+j)^2}.
\]

The line Re s=1/2 is valid even for n=0. Its integrand has exponential
vertical decay and rational factor O(|s|^-2), so Mellin inversion is valid.
Closing to the left for 0<q<1 has positive orientation. The expansion of
pi² csc²(pi s) at s=-k has double coefficient 1 and simple coefficient 0.
Differentiating q^(-s)r_n(-s) therefore gives exactly

\[
 I_n(q)=\sum_{k\ge0}q^k[-r_n(k)\log q-r_n'(k)].
\]

There is no missing factor of two or change of sign. On left half-integer
lines the rational function is uniformly O((N+1+|Im s|)^-2), giving a
vanishing new vertical integral; horizontal segments vanish exponentially.
The q=1 passage is also valid: r_n(k)=O(k^-2), r_n'(k)=O(k^-3), and
|log q| sum q^k |r_n(k)| tends to zero as q increases to 1. The series is
not asserted to converge for q>1 or to be termwise positive.

The original integral is finite for every q>0. Its denominator is at least
min(1,q)(1-c+abc), and the corresponding numerator-to-denominator ratio
for n>0 is bounded on the cube; the n=0 comparison integral is finite.
The shift r_n(k)=R_n(k+n+1), with the double zeros of R_n at 1,...,n,
checks the exact Zudilin-to-Beukers identity without recurrence uniqueness.

## 2. Partial fractions and cancellation at q=1

Independent double-pole extraction gives

\[
 a_i=\binom ni^2\binom{n+i}i^2,\quad
 b_i=2a_i(2H_i-H_{n+i}-H_{n-i}),\quad j_i=n+1+i.
\]

The finite-sum corrections in equation (1) of the main note have the correct
signs: the double pole contributes F-2P_3+log(q)P_2, while the simple pole
contributes G-P_2+log(q)P_1. In particular the F coefficient is
q^(-2n-1)A_n(q), with the reversal of the q exponent as stated.

Because r_n=O(z^-2), sum b_i=0. Hence the limit is

\[
 I_n(1)=2A_n(1)\zeta(3)
 -\sum_i[2a_iH_{n+i}^{(3)}+b_iH_{n+i}^{(2)}].
\]

The n=0 and n=1 normalizations are respectively 2 zeta(3) and
10 zeta(3)-12. The explicit n=1 formula, including its 8 log(q)-12
correction, checks exactly. The product-variable reduction to a positive
one-dimensional integral also checks: its a,b density is
u[(1+u)(-log u)-2(1-u)].

There is a useful exact design statement. For n>=1 put
B_n(q)=sum_i b_i q^(n-i). The quantities

\[
 \delta_i=2H_i-H_{n+i}-H_{n-i}
\]

increase strictly, since

\[
 \delta_{i+1}-\delta_i=\frac2{i+1}-\frac1{n+i+1}+\frac1{n-i}>0.
\]

They begin negative and end positive. Thus the nonzero coefficients of
B_n have exactly one sign change. Descartes's rule and B_n(1)=0 prove that
1 is its unique positive root, with multiplicity one. B_n is positive on
(0,1), negative on (1,infinity). This proof was independently derived in
this review and incorporated in the main note. It concerns the coefficient
of G in this explicit reduction; it does not establish independence of the
four functions or exclude special identities at other parameters.

## 3. The differential ladder is an identity, not a fitted recurrence

For any polynomial P and constant r=r(k), direct differentiation gives

\[
 P(D)\{q^k[-r\log q-r']\}
 =q^k[-P(k)r\log q-(P(k)r)'],\qquad D=q\,d/dq.
\]

The rational identity

\[
 \frac{r_n(k)}{r_{n-1}(k)}
 =\frac{(k+n)^4}{(k+2n)^2(k+2n+1)^2}
\]

therefore proves the stated ladder term by term on 0<q<1, independently of
the Mellin-integral differentiation proof. Absolute convergence locally in
this interval permits all four derivatives. This also checks the derivative
of log q, a term that would be lost if q^k alone were treated as an
eigenfunction. The closed differential rules are DF=G, DG=q log(q)/(q-1),
D log(q)=1. The four-dimensional span is an upper bound, not a minimality
claim. The ladder by itself requires suitable initial/boundary data.

The q=0 endpoint is excluded. The exact leading divergence is
-r_n(0)log q-r_n'(0)+O(q|log q|), with
r_n(0)=(n!)^4/((2n+1)!)². For n=1, r_1'(0)=1/108.

## 4. A real formula on the entire positive axis

No branch choice is needed if the two special functions are continued by

\[
 G(q)=\zeta(2)+\int_1^q\frac{\log u}{u-1}\,du,
\qquad
 F(q)=2\zeta(3)+\int_1^q\frac{G(u)}u\,du,
\qquad q>0.
\]

The quotient log(u)/(u-1) has removable value 1 at u=1. Both functions
are real analytic for q>0 and agree with the main note's polylogarithms
on (0,1), by the derivative identities and their limits at 1.

The positive integral I_n is real analytic on the same connected interval:
for a disk |q-q0|<q0 centered at any q0>0, expand its denominator using
abc/(1-c+q0 abc)<=1/q0. On smaller disks this is dominated by a constant
multiple of the integrable q0 integrand. Hence the formula (1), and the
differential ladder, extend to every real q>0 by analytic uniqueness.
This extends the closed formula, not the q-power series.

Equivalently, boundary values of the principal polylogarithms at q>1
give the same real F and G: the jumps of Li_1, Li_2, Li_3 are proportional
to 1, log q, (log q)^2/2, so they cancel in both combinations. The integral
definition above avoids reliance on a particular sign convention for those
jumps.

## 5. Reproducible independent checks and scope

[check_deformed_period_independent.py](../computations/check_deformed_period_independent.py)
extracts residues by differentiating pole-cancelled rational functions,
checks the actual differential-module identity through n=5, and checks
normalizations and coefficient signs. It also compares the n=0,1 formulas
with direct positive one-dimensional integrals at q=0.4,1.4,2, including
cancellation of the principal-polylogarithm imaginary parts. These finite
checks support the proof; they do not replace its all-n argument.

The main author's checker additionally compares q=1 rational terms against
the Apéry recurrence and checks residue series at four q<1 parameters.
The direct integral check is independent of that recurrence.

The conceptual gain is concrete: Mellin transformation explains why the
same integer binomial coefficients appear in the growing constant-term
period and the decaying inverse period, and why q=1 removes the additional
simple-pole component. This is a proof-design explanation, not a proof of
irrationality at new q values, nor a claim of historical novelty.


## 6. Additional independent audit: the elementary replacement

The rewritten main section 2, equations (4a)--(4c), has now been independently
checked in full. No error was found. It proves the positive integral identity
for every q>0 directly, without complex integration or analytic continuation.
The Mellin proof remains a valid optional second proof.

Write f_n(x)=x^n(1-x)^n and L_n=f_n^(n)/n!. For fixed y,t>0,

\[
 \partial_x^n(qxy+t)^{-1}
 =(-1)^n n!(qy)^n(qxy+t)^{-n-1}.
\]

The n integration-by-parts signs cancel this sign, and Rodrigues's n!
cancels the derivative's factorial. Endpoint terms vanish because all
derivatives of f_n of order less than n vanish at both endpoints. Thus the
first transformed integrand has exactly the factor q^n y^n claimed.
The substitution t=yu contributes dt=y du and cancels y^n against the
power y^(n+1) from qxy+t. It is essential to make this substitution before
the second Rodrigues step.

For fixed x,u>0 the second derivative is

\[
 \partial_y^n(1+yu)^{-1}=(-1)^n n!u^n(1+yu)^{-n-1}.
\]

The same cancellation gives the nonnegative integrand (4b). Finally,
substituting u=(1-c)/(yc), with absolute Jacobian 1/(yc²), gives precisely

\[
 \frac{u^n\,|du|}{(1+yu)^{n+1}(qx+u)^{n+1}}
 =\frac{[c(1-c)]^n\,dc}{(1-c+qxyc)^{n+1}}.
\]

Every y factor cancels. The remaining x,y factors are exactly the numerator
of I_n; the outside factor is q^n, not q^(n+1) or its reciprocal.

The domination argument is sufficient on the full positive q axis.
For completeness, J(q)=integral K(qxy) is finite: on qxy<=1/2,
K(qxy)<=2(|log q|+|log x|+|log y|), and on the remaining region K is
bounded because its apparent singularity at 1 is removable. All logarithms
are integrable on the unit square. Boundedness of L_n therefore justifies
the initial signed triple integral absolutely. After the first integration
by parts,

\[
 \frac{q^nx^ny^n}{(qxy+t)^n}\le1,
\]

so its absolute integrand is bounded by C_n/[(1+t)(qxy+t)]. This has
finite integral C_n J(q). The change t=yu preserves that property. The
second integration by parts is performed pointwise for x,u>0; its new
y-integral is nonnegative and equals the previous signed y-integral,
so it is bounded by the previous absolute y-integral. Tonelli then proves
finiteness and legitimizes the remaining rearrangement. Endpoint x=0 or
y=0 is a null set, and no uniform boundary assertion there is needed.
For n=0 there are no integrations by parts and the same substitutions and
integrability argument still apply.

The polynomial moment formula in (4c) also checks, including its (-1)^n
sign and its shifted denominator. Polynomial expansion makes the integral
a rational function of z; equality first proved for z>n-1 therefore extends
to all z>-1. Differentiation under the integral is dominated locally there.
The square M_n(k)^2=R_n(k+1) has double zeros at k=0,...,n-1. Consequently
both its value and derivative disappear in the geometric expansion; shifting
k=n+ell then cancels exactly the outside q^n. Absolute termwise integration
is justified by C_n² times the positive kernel series, whose integral is
J(q). This proves equation (5) elementarily for 0<q<1. The geometric series
is still not used at q>1.

The rewritten one-crossing proof in section 4 is valid and avoids Descartes:
choose c strictly between negative and positive coefficient indices, or at
the sole possible zero coefficient. For q>1, each nonzero term satisfies
b_i q^(-i)<b_i q^(-c); summing gives a negative value because sum b_i=0.
The inequalities reverse for q<1. At q=1 the derivative is
-sum(i-c)b_i<0, proving simplicity as well as uniqueness.

The independent checker now additionally verifies both differentiated
kernels with their signs and factorials, both substitution Jacobians through
n=5, and the signed Legendre-kernel identity by direct product-variable
convolution for n=0,1,2 at q=0.4,1.4,2. These checks all pass. The exact
all-q identity itself is the preceding real-variable proof.
