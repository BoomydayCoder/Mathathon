# Independent review: consecutive-pole decay and specified arithmetic cost

Verdict: accepted. Root reviewed balanced-consecutive-pole-rate.md, including
its all-m comparison. The claim concerns the proved nonprimitive normalizer,
not a lower bound for actual reduced denominators.

The numerator degree is exactly(m+1)n and denominator degree is
B=(m+1)n+m+1. On compact positive scaled intervals the log n terms cancel
and the Riemann sums give Psi exactly. At zero the denominator is(B!)²,
so Stirling controls the endpoint directly despite its logarithmic integral.
At infinity the excess degree, Jacobian and y^(2m+1) density prefactor cancel
all powers of n and leave v^(-1) times the displayed exponential. Its derivative
is less than-2pi, making the integrated tail uniformly negligible. The positive
saddle is strict, unique and larger than the zero-endpoint exponent. These
facts prove the complete nth-root limit, not merely a stationary-point guess.

The second derivative and the angle equation have the displayed factors.
Evaluation at the saddle cancels both linear and arctangent terms, leaving
the exact ratio raised to m+1. At m=2 the imaginary polynomial is
2(35u²-42u+3); the angle range for the larger root in(1,3) selects pi.
The radical simplification of the rate and coarse bound(2/37)³ agree.
The binomial/factorial nth-root cost is84375/2, whose product already exceeds1
without any lcm contribution.

For general m, the nth-root arithmetic factor is
2(2m+1)^(2m+1)(m+1)^(m+1)/m^m. At v=0 its log plus Psi equals the negative
symmetric second difference of t log t. Bounding the second derivative by1/m
on[m,m+1] gives Delta<=1/(4m). Meanwhile Psi''>=-8(m+1). Evaluation at
v=1/(4m) yields pi/2-(pi+2)/(4m)-1/(4m²), increasing for m>=2 and bounded
below by3pi/8-5/16>0. Thus the stated sufficient normalizer grows against the
positive remainder for every fixed m>=2. The m=1 case is not included in this
last claim. No conclusion about smaller primitive clearing follows.
