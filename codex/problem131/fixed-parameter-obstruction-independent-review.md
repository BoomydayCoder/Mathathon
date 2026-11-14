# Independent audit of the fixed-parameter odd-zeta obstruction

Reviewed fixed-parameter-odd-zeta-obstruction.md and its checker.
**Accepted: no mathematical gap found.** The result proves growth of
actual primitive integer polynomial values, not merely failure of a
sufficient normalization bound.

The quantifiers are correctly restricted: numerator power \(r\ge1\)
and rational aspect ratio \(0\le\alpha<1\) are fixed, then each odd
\(s\) above the explicit threshold is fixed before \(K\to\infty\).
The argument does not rule out varying \(r\) or \(\alpha\) with \(s\),
or choosing a different subspace, functional, or construction.

## 1. Leading coefficient and rational basis changes

Canceling \(D_N\) from the denominator leaves \(W=D_N^{r-1}\) and
\(D=D_K/D_N\). The surviving poles are exactly \(N<j\le K\).
Vandermonde cancellation gives the leading coefficient with power
\(r-1\), including its sign, as stated. It is a nonzero integer for
every \(r\ge1\), so the degree is exactly \(h\).

The rational polynomials \(D/(t+j^2)\) form a basis even when their
change-of-basis determinant is not a \(p\)-adic unit. The determinant
and its leading coefficient acquire the same square factor, which
cancels in the monic ratio. Thus no unjustified unimodularity assumption
enters the local argument.

The denominator lcm of a monic rational polynomial is a positive
integer, and multiplication by that lcm produces a primitive integer
polynomial. This justifies summing necessary costs at different primes:
numerator content at other primes cannot offset them.

The positivity assumption is also valid for odd numerator powers.
\(D_N(y^2)\) is strictly positive on the integration domain, so all
powers \(r\ge1\) give a positive Gram determinant.

## 2. Polynomial quotient degrees and the Bernoulli cutoff

In the interpolation basis, off-diagonal integrands are the integer
polynomials \(WD/[(t+j^2)(t+k^2)]\). The diagonal terms have one simple
pole and an integer polynomial quotient. Both polynomial degrees are
at most
\[
E=(r-1)N+h-2=K+(r-2)N-2.
\]

For \(r=1\), this is \(h-2\); if \(h=1\), \(E=-1\) and the polynomial
part is absent. The proof does not require a nonnegative degree in that
case.

The Bernoulli cancellation input is used with both of its hypotheses:
\(p>s-1\) and \(p>1+(2E+2)/(s-1)\). To see the second condition
directly, if \(p-1\mid 2e+2\), write \(2e+2=a(p-1)\).
The strict cutoff gives \(1\le a\le s-2<p\), and the numerator product
\((2e+3)\cdots(2e+s)\) contains the multiple \(ap\).
It cancels the unique Bernoulli denominator factor \(p\).
The factorial \((s-1)!\) is a unit. If \(p-1\) does not divide
\(2e+2\), the Bernoulli factor is already integral.

Here \(\kappa=1+(r-2)\alpha>0\), including \(r=1\), where
\(\kappa=1-\alpha\). The threshold's strict condition
\((s-1)\beta>2\kappa\) supplies the needed cutoff throughout the band
for all sufficiently large \(K\). The constant terms in the inequality
do not cause an endpoint failure.

## 3. Two removed nodes and the zero class

The prime band has \(N<p\) and \(K<2p\). An ordinary nonzero squared
residue class contains at most four original nodes and at most two
removed nodes. Since \(p^2>2K\), each ordinary same-class square
difference has valuation exactly one: only one of \(j-u,j+u\) is
divisible by \(p\), and neither can be divisible by \(p^2\).

Consequently
\[
v_p(C_j)=(r-1)\delta+(\ell-\delta-1)
=\ell+(r-2)\delta-1.
\]
For \(r\ge2\), the maximum is at most \(2r-1\); for \(r=1\) it is
at most three. Thus the combined bound \(\max(3,2r-1)\) is valid.
Ignoring the possibility of two removed nodes would be incorrect at
large aspect ratios; the proof explicitly retains it.

There is exactly one zero-class node \(j=p\), and it survives because
\(p>N\). All its node-difference factors are units, so
\(v_p(C_p)=s-1\). Hence \(s>\max(3,2r-1)\) makes every polynomial-part
row of \(p^sC^{-1}B\) vanish modulo \(p\), including the zero row.

The harmonic contribution reduces respectively to \(0,-1/2,-1\) for
\(j<p,j=p,j>p\). At \(j=p\), the two surviving terms are \(-1\) from
\(H_p^{(s)}\) and \(+1/2\) from \(1/(2j)\). The term
\(-1/(s-1)\) vanishes after row scaling. This verifies the sign and
the exceptional middle factor in the reduced polynomial.

## 4. The coefficient valuation is exact

The whole scaled monic polynomial is integral and reduces to
\[
Z^{p-N-1}(Z-1/2)(Z-1)^{K-p}.
\]
For \(b=K-p+1\), its coefficient of \(Z^{h-b}\) is the nonzero residue
\((-1)^b/2\). Since that coefficient before reduction is
\(p^{sb}[X^{h-b}]F_K\), its valuation is exactly \(-sb\).
This implies a necessary denominator exponent at least \(sb\);
the note correctly does not call it the exact denominator of the
entire polynomial.

This remains valid when \(b=h\), when \(p=K\), and when \(h\ge p\).
No inversion of characteristic indices modulo \(p\) is part of the
proof. PNT partial summation then gives exactly
\(s\int_\beta^1(1-u)\,du=s(1-\beta)^2/2\).

## 5. The real lower estimate and explicit threshold

The Jensen estimate is in the required lower-bound direction.
The log-distance singularities and the endpoint logarithm are
integrable. The total \(K\)-power is
\[
2(r-1)Nh+(s-1)h,
\]
matching the leading coefficient exactly. Remaining \(s\)-dependent
terms are lower order for fixed \(s\), so the limiting real lower
constant is independent of that fixed index.

The assertions \(I_\alpha\le0\) and \(I_0=0\) are immediate from
\(0<x^2-u^2\le1\) on \(0\le u\le\alpha\le x\le1\), apart from an
integrable boundary singularity. Thus its contribution can be dropped
in a lower bound for every \(r\ge1\).

With \(b=1/2\),
\[
J_\alpha\ge-2\alpha(\log2+1)\ge-4\alpha,\quad
J_1<1/4,\quad\pi/2<2.
\]
The last two inequalities are strict, so \(L_{r,\alpha}>-C\) remains
strict even at \(\alpha=0\). At \(r=1\), the leading-coefficient
integral disappears without a singular limiting argument.

The stated formula for \(S\) is indeed the smallest odd integer
strictly greater than \(T\). It makes the local row bound, eventual
Bernoulli cutoff, and final positive margin hold simultaneously.
Choosing \(c\) below half the displayed rational margin is safe.

## 6. Independent exact matrix tests

The new checker check_fixed_parameter_independent.py imports no source
helper. It constructs every integer quotient and residue, evaluates
the full matrix in exact fractions, and computes its entire rational
characteristic polynomial. It then verifies the complete reduction
and the claimed exact coefficient valuation, not just the leading
factors or separate moment bounds.

All eight cases pass:

| \(r,K,N,s,p\) | Edge or feature |
|---|---|
| \(1,5,4,5,5\) | \(h=1,E=-1\), no polynomial part |
| \(1,7,0,5,5\) | \(\alpha=0,h>p\) |
| \(1,7,3,5,5\) | power one, two removed nodes |
| \(2,7,3,5,5\) | power two, two removed nodes |
| \(3,11,6,7,7\) | ordinary row cost attains \(2r-1\) |
| \(4,17,10,9,11\) | larger numerator and high aspect ratio |
| \(6,20,12,13,13\) | cost eleven, \(s=p\), full eight-dimensional polynomial |
| \(1,20,18,5,19\) | \(\alpha=9/10\), coefficient with \(b=h\) |

Thirty independent exact threshold checks also pass, including
\(\alpha=99/100\). The original checker passes its nine local cases
and five threshold examples.

## 7. Scope

The conclusion covers all nonzero rational scalar normalizations and
all rational basis changes within the same full polynomial space.
The basis may depend on \(K,s\); its coefficients are rational and
independent of the indeterminate \(X\).

It does not cover an \(X\)-dependent transformation, a different
subspace dimension, a changed functional, or numerator/aspect-ratio
choices depending on \(s\). The PNT input and the order of limits are
explicit. Within that scope, the claimed actual primitive-growth
obstruction follows.
