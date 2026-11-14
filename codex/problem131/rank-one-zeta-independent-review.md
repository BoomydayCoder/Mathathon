# Independent review: rank-one subspace approximants

Status: independently reviewed; no substantive error found in the proved
parts of [rank-one-zeta-subspace-approximants.md](rank-one-zeta-subspace-approximants.md).
The rank, Schur-root construction, strict monotone convergence, factorial
bound, trial asymptotic and multipoint-Padé remainder are valid. The final
adjacent-Apéry formula remains a conjecture; this review does not prove it
or any higher-odd-zeta irrationality statement.

## 1. The degree drop and exact Schur normalization

For s=2m+1, the X coefficient at pole -j² is j^(2m)=(-1)^m(-j²)^m.
Multiplying the pole residue of D_N^r t^e/D_K gives exactly the displayed
Lagrange sum. The numerator has degree m+e+rN and is monic. Therefore
coefficients vanish below e0=K-1-rN-m and equal (-1)^m at e0.

If h-1<=e0<=2h-2, put q=2h-e0-1 and a=h-q. The largest index sum
in any of the first a rows is (a-1)+(h-1)=e0-1, so these rows and columns
of the coefficient matrix vanish. In the remaining block its anti-diagonal
occurs at local index sum q-1. Entries before it vanish, and every entry
on it is (-1)^m. Reversing columns makes the block triangular and gives
determinant (-1)^(mq+q(q-1)/2). Consequently the rank is q.

Only the bottom q rows can supply X factors to the determinant. At degree
q their nonzero columns must also be in the bottom block. The top block
therefore contributes its determinant independently, proving the stated
leading coefficient with no omitted cross terms. That top block is a
positive principal Gram block and has nonzero determinant.

For N=0 and K=2h+m-1, q=1 and the sole X coefficient is exactly (-1)^m,
not an additional product of pole differences. Completing the Schur square
therefore yields

\[
 E_{h,s}=(-1)^m(\zeta(s)-q_{h,s})
 =\min_{p\text{ monic},\,\deg p=h-1}
    \int\frac{p(y^2)^2}{D_K(y^2)}w_s(y)\,dy>0.
\]

The normal equations use moments of index at most 2h-3, where the X
coefficient vanishes. Thus their coefficients and the minimizer are rational.
For h=1 the empty matrix convention gives p=1, q=(-1)^(m+1)c, and E equal
to the total positive mass. The side of approximation is as stated.

## 2. Strict monotonicity and convergence

At the next h, K increases by two. The monic trial polynomial
(t+(K+1)²)p_h(t) has its norm multiplied pointwise by
(t+(K+1)²)/(t+(K+2)²)<1. Since the old error density is positive off the
finitely many polynomial zeros and has nonzero mass, its integral decreases
strictly. Minimizing can only decrease it further. This proves strict
monotonicity despite the changing measures.

The trial polynomial D_{h-1} gives

\[
 R_h(0)=\frac{((h-1)!)^4}{((2h+m-1)!)^2},\qquad
 \frac{R_h(y)}{R_h(0)}
 =\frac{\prod_{j<h}(1+y^2/j^2)}
        {\prod_{h\le j\le2h+m-1}(1+y^2/j^2)}
 \le\frac{\sinh\pi y}{\pi y}.
\]

This also covers h=1 (the numerator product is empty). The dominating
function times w_s is integrable: near zero both factors have finite limits,
and at infinity it is O_s(y^(s-1)e^(-pi y)). Thus the explicit factorial
bound is valid for every h>=1, and its right side tends to zero. Together
with strict monotonicity this proves the claimed one-sided convergence.

## 3. Constants J_s and the trial asymptotic

With w_s=c_s y^s sum_{r>=1}r^(s-1)e^(-2pi r y), where
c_s=2(2pi)^(s-1)/(s-1)!, integration of each positive exponential difference
in w_s sinh(pi y)/(pi y) gives

\[
 J_s=\frac{2^{s-1}}{\pi^2}
 \sum_{r\ge1}\frac{r^{s-1}-(r-1)^{s-1}}{(2r-1)^s}.
\]

The endpoint term (r-1)^(s-1) at r=1 is zero. Writing v=2r-1 and expanding
((v+1)/2)^(s-1)-((v-1)/2)^(s-1) leaves precisely the odd binomial terms.
Since sum over positive odd v of v^(-j-1) equals
(1-2^(-j-1)) zeta(j+1), this proves equation (4), including its factor
2/pi². It checks J3=1/2, J5=1+pi²/12, and
J7=3/2+5pi²/12+pi⁴/80.

For fixed y the finite ratio tends to the infinite Euler product: the
omitted denominator logarithm is bounded by y² sum_{j>=h}j^(-2), which
tends to zero. Dominated convergence therefore shows that the trial norm
divided by R_h(0) tends to J_s. The factorial calculation is

\[
 \frac{\Gamma(h)^4}{\Gamma(2h+m)^2}
 \sim 2^{2-2m}\pi\,h^{-2m-1}16^{-h}
 =2^{3-s}\pi\,h^{-s}16^{-h}.
\]

Hence C_s=2^(3-s)pi J_s is correct; in particular C3=pi/2. This is the
asymptotic of the specified trial norm only. The optimized error inherits
an upper bound, not an asymptotic equivalence or a matching lower rate.

## 4. Rational multipoint-Padé data and the exact remainder

Write n=h-1, rho0=mu_s/D_(m+1), and Omega_n=D_(2n+m+1)/D_(m+1).
Then rho0/Omega_n is exactly the Schur measure, so no weight has changed
in the interpretation. At z=-j² with j>m+1, S(z) is minus the functional
of 1/[D_(m+1)(t)(t+j²)]. This denominator has m+2 distinct poles and its
zeta coefficient is the Lagrange sum of a degree-m polynomial, strictly
below m+1. Thus all interpolation values are rational.

The displayed P_n integral is polynomial in z because its numerator
vanishes at z=t. Polynomial division shows that its coefficients of degrees
n,...,2n-1 involve integrals of p_n(t)/Omega_n(t) times polynomials of
degree at most n-1. These vanish by orthogonality. Thus deg P_n<=n-1.
Direct subtraction from p_n(z)S(z) first gives

\[
 S(z)-P_n(z)/p_n(z)
 =\frac{\Omega_n(z)}{p_n(z)}
   \int\frac{p_n(t)}{\Omega_n(t)(z-t)}\,d\rho_0(t).
\]

Multiplying the integral by p_n(z) and replacing p_n(z) by p_n(t) changes
it by an orthogonality integral: (p_n(z)-p_n(t))/(z-t) is polynomial in t
of degree n-1. This proves the squared-denominator remainder (7) with its
stated sign for the convention S(z)=integral (z-t)^(-1)d rho0.

The orthogonal polynomial has simple positive roots by the usual positive
measure theorem, so it does not vanish at the negative interpolation nodes.
The remainder vanishes there. Rationality of p_n was already proved; any
n of the rational samples then determine the n coefficients of P_n by
ordinary rational interpolation. This also establishes rationality without
assuming the integral formula's individual coefficients are rational.

As z tends to negative infinity, Omega_n(z)/p_n(z)² tends to 1 and
0<z/(z-t)<=1 for t>=0. Dominated convergence gives

\[
 \lim_{z\to-\infty}z\left(S(z)-P_n(z)/p_n(z)\right)=E_{h,s}.
\]

The dominating integrand is the finite Schur norm. The 1/z coefficient of
S is its total mass (-1)^m zeta(s)+c0. Therefore the leading coefficient
of P_n equals c0+(-1)^m q_(h,s), as stated. For n=0, P0=0 and the remainder
is S itself; the empty interpolation problem and mass identity are valid.
If coefficient extraction at degree -1 is used in that edge case, it means
zero for the zero polynomial, not an extra Padé coefficient.

## 5. Exact checks and unresolved arithmetic

[check_rank_one_independent.py](../computations/check_rank_one_independent.py)
builds pole residues directly from products D_K'(-j²), rather than the
closed factorial formula used by the original checker. It verifies 20
Schur systems for s=3,5,7,9 and h=1,...,5, their rational orthogonality and
norm identities, and 40 exact Padé remainders at additional rational points.
It constructs P_n from rational interpolation values and checks all 2n nodes
independently. It also checks 54 general anti-diagonal rank/leading-block
cases and the exact J3,J5,J7 constants. All pass.

The original exact Schur-elimination script was inspected. Its factorial
residue equals 1/D_K'(-j²), its affine pole constant has the correct endpoint
terms, and its symmetric elimination computes the Schur complement. A
read-only rerun recomputed all 48 archived exact roots and matched them.
The 180-digit error diagnostics are not used as proof of sign, convergence,
or a denominator asymptotic; those claims are kept separate in the note.

No all-h identity with the Apéry recurrence has been established in this
review. In particular, the adjacent-Apéry formula in section 7 remains
explicitly conjectural despite exact finite agreement. For higher odd s,
the available exponential error upper bound alone gives no usable control
of the reduced rational denominators. The earlier full-pole scalar
obstructions do not automatically apply to this smaller subspace, and no
irrationality conclusion is drawn either way.
