# Independent review: balanced diagonal arithmetic and decay

Verdict: accepted. Root independently reviewed
[the numerator prime theorem](balanced-diagonal-numerator-primes.md) and
[the positive-integral rate](balanced-diagonal-positive-rate.md).
The resulting limitation concerns the specified lcm multipliers, not primitive
denominators and not arbitrary higher-zeta approximants.

## Arithmetic

At a=n+1 the Andrews sum has integer Z_k factors and its final factor becomes
binom(n+1,k), proving R(a) integral for every m,n. The coefficient theorem gives
R in Z[1/n!][x]; monic division by x-a keeps this property. At p>n+1, a possible
Bernoulli denominator in the r-th derivative, r=2m-1, has k=ell(p-1), with
1<=ell<=r because k<=2mn and p-1>=n+1. The derivative product then contains
k+ell=ell p. This cancels the possible simple denominator. The k=0,1 cases
are separately harmless at these odd primes. At p=n+1 only a single possible
Bernoulli denominator remains. The harmonic term loses at most one p after
the a^(2m) coefficient; fixed factors are covered by (2m)!, including p=2
when n=1. Thus both claimed primewise conclusions are valid.

At p<=n, the earlier Int(Z) theorem clears b using the p-part of d_L with
L=2mn+2m-1; L>=n+1. The stated Lambda retains precisely these parts. At
p>n, the new theorem and the extra n+1 factor suffice. This proves U b integral
prime by prime. If Lambda/d_n has a prime factor p, a power p^e with e>=2
lies in (n,L], so p<=sqrt(L). Its total logarithmic contribution is at most
sqrt(L) log L. Hence the fixed-m asymptotic log U=(4m-1)log d_n+o(n) is
correct. This does not prove the stronger observed d_n^(2m+1) clearing.
The derivative cancellation is the same mechanism already present in the
earlier fixed-odd-index polynomial-moment threshold; it is not asserted new.

The independent checker reconstructs b directly using ordinary rational
moments and polynomial division, rather than the author's derivative formula.
It verifies the diagonal integer coefficient, prime support and complete U.

## The full positive Laplace estimate

The numerator degree in t is exactly (m+1)n; the denominator degree is n+1.
The shifts r+1 through3r+1 and the factorial normalization in equation(1)
are consistent with the previously accepted moving-pole identity. On a compact
positive v interval, the Riemann sums and the hyperbolic density give the stated
Phi, with all log n terms canceled. The density's bounds follow from its finite
positive csch decomposition: it is bounded at zero and has polynomial times
exp(-2pi y) behavior at infinity.

The small-v upper bound uses numerator monotonicity and D(v)>=D(0); its rate
Psi(delta) tends to Phi(0). Since Phi'(0+)=2m pi, this endpoint is strictly
subdominant. The large-v bound has exactly the polynomial prefactor
n^(2m)v^(2m-1) and exponent H_m in the note. For large fixed M, H'_m<=-pi
and H_m(M) tends to minus infinity, which controls the integral of the entire
tail, not just each fixed v. These bounds together with compact convergence
and a neighborhood of the saddle prove the claimed nth-root limit.

The displayed second derivative is strictly negative for m>=1. Hence the
saddle is unique. Integrating the logarithm explicitly cancels the angle and
constant terms and gives rho exactly. For m=1 the sole positive algebraic
root has the required angle pi and yields17-12sqrt(2). For m=2 the unique
root in(2,3) has angle between2pi/3 and17pi/12, so the algebraic root is on
the required branch. The coarse interval already proves rho_2>9/128.

The elementary binomial divisibility bound gives liminf d_n^(1/n)>=2.
Thus the specified fifth-power lcm multiplier has growth rate at least
32rho_2>9/4. At the fifth-value saddle, increasing m increases Phi by more
than2-log2; this proves the corresponding statement for every fixed m>=2.
All these claims concern positive even-n forms. Additional cancellation in
actual reduced denominators remains open. No irrationality impossibility is
inferred from the failure of these sufficient multipliers.

## Fixed aspect ratios, subsequently reviewed

The added section5 is accepted. For a_n/n->c>1/2 the denominator endpoints
stay away from zero, so the same uniform estimates cover their o(1) changes.
The second derivative calculation has strictly negative sign; the displayed
positive numerator verifies the comparison for every endpoint a>=0. At v=1/2
the saddle derivative is positive for m>=2, using the maximal angle difference
atan2. The envelope derivative is therefore between -log6 and0. The elementary
lcm cost is constant to the left of c=1 and has derivative at least5log2 to
the right. The resulting bound is minimized at c=1. This proves the claimed
limitation for every fixed ratio, independent of whether that multiplier
actually clears the associated numerator. The primitive-arithmetic caveat is
essential and correctly retained.

## Uniform all-pole comparison, subsequently reviewed

The added section6 is accepted. Positivity and factorwise comparison give
L(a)>=L(n+1)[K(n+1)/K(a)]² for a>=n+1, while L is decreasing throughout
a>n/2. Midpoint convexity of1/x bounds the logarithmic derivative ofK by
the integral fromr+1/2 to3r+3/2, exactlylog3. Therefore the loss is at most
9^(a-n-1). The elementary lower bound for the lcm compensates uniformly:
for A>=n+1>=5, the ratio of successive lower bounds is at least
(5/3)^5/9=3125/2187>1. The cases A=n,n+1 are checked directly. Hence the
infimum over ALL integer a>n/2 obeys the stated lower bound, with no
aspect-ratio assumption. Combining with the independently accepted diagonal
rate proves the strict9/4 threshold. The proof retains its essential scope:
only the specified lcm multiplier or a larger multiplier, not a lower bound
on actual primitive denominators.
