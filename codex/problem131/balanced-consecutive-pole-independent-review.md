# Independent review: consecutive-pole Rodrigues transport

Verdict: accepted. Root reviewed balanced-consecutive-pole-transport.md.
The accepted balanced Rodrigues/adjoint theorem is the analytic input.
No successful irrationality inequality is claimed.

Putting the two half-shifted rational kernels over their union denominator
gives the displayed two quadratic numerators. Their difference is exactly
2Kx, so Dcal g_(a,K)=K g_(a-1/2,K+1), including the sign and parameter shift.
Iteration gives (K)_n. For a>n/2 the closest shifted pole remains beyond
the required half-width strip at every adjoint step. The prior gamma-weight
estimates apply unchanged and the even-n periodic multipliers return to their
initial values. This justifies replacing the single-pole n! factor by (K)_n
in the previously proved transport. The resulting integrand is positive.

When K=mn+m+1, numerator degree is K-m-1. Expansion at infinity has leading
term t^(-m-1); partial fractions then give sum c_j j^(2m)=(-1)^m. Inserting
the exact resolvent makes the form (-1)^m(zeta-q), with all lower target
coefficients absent by degree. The transport prefactor is exactly1. At the
minimal a=n/2+1 the denominator interval is precisely1 through
(m+1)n+m+1. This verifies all signs, endpoints and normalizations.

For clearing, the missing initial factors contribute
((a-1)!)² binom(j-1,a-1)binom(j+a-1,a-1), up to the already tracked sign.
The full-grid derivative gives2j² binom(2b,b-j)/(2b)!, again with sign only.
Substitution of p(-j²)=(-1)^(mn)(n!)^(2m)R(j)/C leaves M c_j equal to
2m j² d_b^(2m+1) times an integer, by the accepted grid theorem H R(j) inZ.
This clears all three rational parts, including the 1/(2j) correction.
The factorial ratio is itself an integer: its denominator arguments total
2mn+2a-2, leaving2m+2 factorial units in the numerator multinomial.
Its log grows only linearly in n when a=O(n), with fixed m; the binomial
and elementary lcm estimates have the same property. No comparison of that
cost with the actual positive remainder is inferred.

The algebraic target-isolation and normalizer arguments hold for odd n;
positivity is explicitly restricted to even n. The author's checker passed
90 exact iterated identities,48 coefficient/normalizer systems and4 positive
quadratures. The written proof, not those finite tests, supplies the all-index
assertion.
