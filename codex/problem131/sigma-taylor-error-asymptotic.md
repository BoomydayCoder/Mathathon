# A sinc-squared asymptotic for every fixed odd-zeta parameter coefficient

Status: proof developed from the root agent's uniform-parameter proposal.
The proof below depends on the all-index associated-numerator recurrence in
[deformed-associated-pearson-symmetry.md](deformed-associated-pearson-symmetry.md),
which has a full written proof accepted in its independent audit.
The present analytic argument is accepted in [sigma-taylor-independent-review.md](sigma-taylor-independent-review.md). No historical
novelty or higher-odd irrationality claim is made.

The denominator family and recurrence have an exact Askey–Wilson (1984)
antecedent, as documented in
[koecher-askey-wilson-primary-dictionary.md](koecher-askey-wilson-primary-dictionary.md).
Thus this result is an asymptotic analysis of that classical deformation;
priority of this particular coefficient formulation has not been determined.

## 1. Statement with exact normalization

Write

\[
 C_{nk}=\binom nk\binom{n+k}k,\qquad
 w_k(\sigma)=\prod_{j=1}^k(1-\sigma/j^2)^{-1},\qquad
 A_n(\sigma)=\sum_{k=0}^n C_{nk}^2w_k(\sigma).
\]

Let b_n be the associated rational numerator of the positive deformation,
and put

\[
 q_n(\sigma)=b_n(\sigma)/A_n(\sigma),\quad
 Z(\sigma)=\sum_{h=1}^\infty\frac1{h(h^2-\sigma)}
          =\sum_{r=0}^\infty\zeta(2r+3)\sigma^r,\quad
 E_n(\sigma)=Z(\sigma)-q_n(\sigma).
\]

The initial values are A_0=1, b_0=0,
A_1=(5-sigma)/(1-sigma), and b_1=2(3-sigma)/(1-sigma)².
All Taylor coefficients of q_n at zero are rational.

Then uniformly on the closed disk |sigma|<=1/4,

\[
 \boxed{\frac{E_n(\sigma)}{E_n(0)}\longrightarrow
 \left(\frac{\sin(\pi\sqrt\sigma)}{\pi\sqrt\sigma}\right)^2.} \tag{1}
\]

The expression is interpreted by its entire power series, so there is no
branch ambiguity at zero. In particular, for every fixed integer r>=0,

\[
 \boxed{\zeta(2r+3)-[\sigma^r]q_n(\sigma)
 \sim \frac{(-1)^r2^{2r+1}\pi^{2r}}{(2r+2)!}\,E_n(0).}       \tag{2}
\]

Consequently these rational approximants converge, and their errors have
eventual sign (-1)^r. Thus the first parameter coefficient eventually lies
ABOVE zeta(5), the second BELOW zeta(7), and the signs continue alternately.
The theorem concerns each fixed r; it does not assert a threshold uniform in r,
a sign at every finite n, or control of reduced rational denominators.

## 2. Exact error series from the recurrence

The established common recurrence, for n>=1, is

\[
 \alpha_n y_{n+1}-\beta_n y_n+\alpha_{n-1}y_{n-1}=0,
\qquad
 \alpha_n(\sigma)=
 \frac{(n+1)((n+1)^2-\sigma)^2}{2(3(n+1)^2-\sigma)}.
\]

It holds for y=A and y=b. The initial b_1 is exactly 1/alpha_0.
For W_n=b_n A_(n-1)-b_(n-1)A_n the recurrence gives
alpha_n W_(n+1)=alpha_(n-1)W_n. Therefore

\[
 W_n=1/\alpha_{n-1},\qquad
 q_k-q_{k-1}=h_k(\sigma):=
 \frac1{\alpha_{k-1}(\sigma)A_k(\sigma)A_{k-1}(\sigma)}.       \tag{3}
\]

For real 0<=sigma<=1/4 the previously proved positive Padé error gives

\[
 0<E_n(\sigma)\le
 \frac2{(2n+1)^3A_n(\sigma)^2}\longrightarrow0.
\]

Consequently telescoping (3) yields the exact tail formula on that interval:

\[
 E_n(\sigma)=\sum_{k=n+1}^\infty h_k(\sigma).                 \tag{4}
\]

We will justify its analytic extension, rather than assuming that a positive
real error can be differentiated with its sign preserved.

## 3. A common disk with no denominator zeros

Define the uniformly convergent product on |sigma|<1

\[
 G(\sigma)=\prod_{j=1}^\infty(1-\sigma/j^2)^{-1}
 =\frac{\pi\sqrt\sigma}{\sin(\pi\sqrt\sigma)}.
\]

The last identity is Euler's sine product. Every w_k has nonnegative Taylor
coefficients. On |sigma|<=1/4, coefficientwise triangle bounds give

\[
 |w_k(\sigma)-1|\le w_k(1/4)-1
 \le G(1/4)-1=\pi/2-1=:c_0<1.
\]

The last strict inequality requires only pi<4. Since the coefficients C_nk²
are positive,

\[
 \boxed{\left|\frac{A_n(\sigma)}{A_n(0)}-1\right|\le c_0<1}
 \quad\hbox{for EVERY }n.                              \tag{5}
\]

Thus all denominators are nonzero on the same closed disk and
(1-c_0)A_n(0)<=|A_n(sigma)|<=(1+c_0)A_n(0).
The moment formula defining b_n, or its rational recurrence, has no singularity
on this disk. Hence q_n and E_n are holomorphic on its interior and continuous
on it.

For convergence of (4), note the elementary bound

\[
 A_n(0)\ge\binom{2n}n^2\ge\frac{16^n}{(2n+1)^2}.             \tag{6}
\]

The second inequality follows because the largest coefficient in
(1+x)^(2n) is at least the average of its 2n+1 coefficients.
Also alpha_(k-1)(0)=k³/6, and uniformly on our disk

\[
 \frac{\alpha_{k-1}(0)}{\alpha_{k-1}(\sigma)}
 =\frac{1-\sigma/(3k^2)}{(1-\sigma/k^2)^2}                  \tag{7}
\]

is bounded by an absolute constant. Equations (5)–(7) bound |h_k(sigma)|
by a constant times h_k(0), whose sum converges by (6).
The tail series therefore converges uniformly on the closed disk and is
holomorphic inside. It agrees with E_n on a real interval, so the identity
theorem proves (4) throughout the disk, with its boundary values obtained
by continuity. There is no exchange of divergent series or unproved
nonvanishing assumption.

## 4. Uniform denominator asymptotics need only finite-head escape

The product w_k converges uniformly to G on |sigma|<=1/4. We claim

\[
 \frac{A_n(\sigma)}{A_n(0)}\longrightarrow G(\sigma)
 \quad\hbox{uniformly on that disk}.                         \tag{8}
\]

The quotient is a convex combination of w_k with positive weights
C_nk²/A_n(0). Fix an integer M. For k<=M, C_nk² is a polynomial in n
of degree 4k, while A_n(0) has the exponential lower bound (6). Hence

\[
 \frac{\sum_{k=0}^M C_{nk}^2}{A_n(0)}\longrightarrow0.
\]

The contribution of this fixed head to the difference in (8) tends uniformly
to zero, since all w_k and G are uniformly bounded. The remaining contribution
is bounded by sup_(k>M, |sigma|<=1/4)|w_k(sigma)-G(sigma)|, which tends to zero
as M tends to infinity. Choosing M first and then n proves (8).
No saddle-point estimate or detailed Apéry asymptotic is needed.

## 5. The tail retains the same analytic multiplier

Equations (7) and (8), together with the nonvanishing bounds, imply

\[
 \frac{h_k(\sigma)}{h_k(0)}\longrightarrow G(\sigma)^{-2}
 \quad\hbox{uniformly on }|\sigma|\le1/4.
\]

Because h_k(0)>0, the normalized tail is again a convex combination:

\[
 \frac{E_n(\sigma)}{E_n(0)}
 =\sum_{k>n}\frac{h_k(0)}{\sum_{j>n}h_j(0)}
                           \frac{h_k(\sigma)}{h_k(0)}.
\]

Its uniform difference from G^(-2) is at most
sup_(k>n, |sigma|<=1/4)|h_k(sigma)/h_k(0)-G(sigma)^(-2)|,
which tends to zero. This proves (1).

Apply Cauchy's coefficient formula on any fixed circle of radius between
zero and 1/4. Uniform convergence there justifies coefficient extraction.
Finally

\[
 G(\sigma)^{-2}
 =\frac{1-\cos(2\pi\sqrt\sigma)}{2\pi^2\sigma}
 =\sum_{r\ge0}\frac{(-1)^r2^{2r+1}\pi^{2r}}{(2r+2)!}\sigma^r.
\]

Every displayed coefficient is nonzero, proving the asymptotic equivalence
and eventual sign in (2). The first three multipliers are
1, -pi²/3, and 2pi⁴/45. In particular the pi² factor is the first
coefficient of the reciprocal squared Euler product, not a numerical fit.

## 6. Arithmetic and verification scope

The construction isolates a single odd zeta value at each coefficient order,
and the analytic errors are asymptotically nonzero multiples of the classical
Apéry rational-approximation error. This does not settle their integer-scaled
errors. The local denominator-power analysis in
[sigma-taylor-zeta-experiment.md](sigma-taylor-zeta-experiment.md) remains a
separate arithmetic issue. No irrationality conclusion for a new odd zeta
value follows from this analytic theorem.

[check_sigma_taylor_asymptotic.py](../computations/check_sigma_taylor_asymptotic.py)
checks exact finite Taylor arithmetic, the denominator recurrence against the
independent finite sum, initial normalization, and the Casoratian. High-precision
coefficient-error and complex-parameter tests corroborate the limiting
multipliers. Uniform convergence and the all-index conclusions are proved
above, not inferred from these numerical tests.
