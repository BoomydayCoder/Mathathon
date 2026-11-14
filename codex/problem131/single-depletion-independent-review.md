# Independent review of the single-depletion scalar obstruction

Status: independently reviewed; no substantive error found in
[single-depletion-monic-obstruction.md](single-depletion-monic-obstruction.md).
The conclusion is an actual primitive-value lower bound for the specified
family, not a comparison with an unnecessarily expensive sufficient
normalizer. The prime number theorem and the previously certified canonical
monic Jensen lower bound are explicit inputs.

Throughout, L>=2 and odd s>=7 are fixed, K=40n, N=3n, h=K-N, and the
numerator power remains six. The result is not uniform when L or s varies
with K; it makes no assertion about s=3,5 or other families.

## 1. Positivity, moments, and the real comparison

Using w_s(y)=c_s y^s sum_{r>=1}r^(s-1)exp(-2 pi r y), subtraction of
L^(-1)w_s(Ly) removes exactly the frequencies divisible by L. Every
remaining term is positive and frequency 1 remains. Since w_s(y) tends
to 1/pi at zero, the ratio sigma_L/w_s tends to 1-1/L there; its limit
at infinity is 1. Positivity and continuity on intermediate compact
intervals give two-sided bounds by fixed positive constants.

Consequently any two corresponding h-dimensional positive Gram
determinants differ logarithmically by O_{s,L}(h). Their monic leading
coefficients differ by C^h, C=1-L^(s-1), so their absolute monic values
also differ logarithmically by O_{s,L}(h). This justifies transferring
the canonical Jensen constant without reversing an upper bound.

The moment formula checks exactly. For d=2e+2, the cumulant of the
uniform distribution on 0,...,L-1 is
kappa_d=(L^d-1)B_d/d. Substitution into the displayed cumulant expression
reproduces (1-L^(-2e-2)) mu_s(t^e). Alternatively, the cumulant recurrence
shows directly that its denominators are supported on primes dividing L.
No p>s condition is needed; the binomial factor is an integer.

The pole formula follows by applying the original Hermite pole identity
at j and Lj and subtracting. The 1/(2j) endpoint terms leave
(1-L^(-1))/(2j). This is essential at the node j=p.

## 2. Finite exceptions and exact residues

For an ordinary node, scaling its pole constant by p^s/(Cj^(s-1)) leaves
exactly the residue in equation (3). At j=p, the endpoint correction
survives. The calculation is

\[
 \frac{p^sc_p}{Cp^{s-1}}
 \equiv\frac{L^{s-1}H_L^{(s)}-1+(1-L^{-1})/2}{C}
 =\frac{L^{s-1}H_L^{(s)}-(1+L^{-1})/2}{C}\pmod p.
\]

Thus the corrected zero-node constant is right. The quantities u_k and
u_* used to exclude numerator and denominator primes are all nonzero
positive rationals, so these are genuinely finite fixed exclusions.
The residues vanish exactly for Lj<p and are otherwise nonzero.

A useful clarification closes a possible small-prime ambiguity: the
listed exceptional set already excludes every prime p<=2L. If p divides
L this is immediate. Otherwise p<=2L-1 and u_p has denominator divisible
by p: H_p has exactly one denominator divisible by p, while
H_floor(p/L) has index at most 1. Since p does not divide L, cancellation
is impossible. Therefore admissible primes automatically satisfy p>2L.
Together with K<2p this implies LK<p². There is no unaccounted p² term
inside H_{Lj}, and the short harmonic sums in the residue formulas have
unit denominators. The note's stated finite hypotheses are sufficient.

## 3. Whole-matrix integrality and the exceptional rows

The integer basis f_j=D/(t+j²) yields
F_K(X)=det(XI+Ccal^(-1)B). Its change-of-basis determinant cancels exactly
between the original determinant and its leading coefficient; it need
not be a local unit. Off-diagonal polynomial quotients and diagonal
polynomial remainders have integer coefficients. Applying the moment
integrality therefore proves integrality of the entire polynomial part,
not just its diagonal entries.

For ordinary squared residue classes the valuation is ell+4 delta-1.
Indeed ell-delta-1 surviving partners contribute one each and the removed
node contributes five. The inequalities p²>2K and 2N<p make these counts
exact and give ell<=4 and delta<=1. At j=p all the node-difference factors
are units and the valuation is s-1. Hence p^s Ccal^(-1)B is integral.
For odd s>=9 every polynomial-part row vanishes after reduction. For
s=7 this remains true except where ell=4 and delta=1.

Those exceptional classes are exactly the removed a satisfying
2p-K<=a<=N; their surviving nodes are p-a,p+a,2p-a. The classes are
disjoint, their number is d=(K+N-2p+1)_+, and all 3d surviving nodes
exceed p/2, hence p/L. They all belong to the counted nonzero pole
residues. This last fact is necessary for the loss bound.

Order nonexceptional rows first. They are literally diagonal, including
zero entries in exceptional columns. The full reduced characteristic
matrix is block triangular even though exceptional rows can couple to
every outside column. Its determinant is the outside diagonal product
times a monic polynomial of degree 3d. If the latter has vanishing order
t0 at zero, then 0<=t0<=3d and the full determinant has vanishing order
h-b_p+t0. Its first nonzero coefficient has descending index
r=b_p-t0>=b_p-3d. Undoing the p^s scaling gives exact valuation -7r.
This proves the actual coefficient-denominator bound without assuming
that the residual block is nonsingular or its constant term nonzero.

For clean rows the same argument gives the particular coefficient
valuation -s b_p exactly. No statement about which coefficient has the
largest denominator is needed.

## 4. Prime sum and primitive normalization

Let alpha=3/40. The clean prime cost for L=2 is

\[
 s\int_{1/2}^{2}(1-u/2)\,du=9s/16.
\]

For L>=2 the corresponding integrand on [1/2,2] increases, and the
additional interval contributes nonnegatively. The s=7 loss is at most

\[
 3s\int_{1/2}^{(1+\alpha)/2}(1+\alpha-2u)\,du
 =3s\alpha^2/4.
\]

It is legitimate to subtract this full loss even if b_p-3d is negative
at an individual prime: replacing that negative lower bound by zero
only improves the result. Floors and the +1 contribute O(K), by
sum_{p<=LK}log p=O_L(K). Fixed exceptions eventually lie below the band.
The prime number theorem and partial summation therefore give the claimed
liminf cost 25011/6400 for s=7; s>=9 has a larger lower bound.

The independently certified Jensen constant is greater than -1789/500.
The exact margin is

\[
 \frac{25011}{6400}-\frac{1789}{500}
 =\frac{10559}{32000}>\frac14.
\]

The density comparison in section 1 transfers that constant to this
normalized model. Thus the logarithm of the actual primitive value is
eventually greater than K²/4.

To check the final arithmetic implication: F_K is monic, so its least
positive coefficient-clearing integer M_K is a positive integer and
M_K F_K is primitive. Otherwise a common prime divisor would divide its
leading coefficient M_K and could be removed from M_K, contradicting
minimality. Hence denominator costs at the specified primes cannot be
offset by numerator factors elsewhere. Every rational multiple of the
original determinant having integer coefficients is an integer multiple,
up to sign, of this primitive polynomial. Absolute values correctly handle
the sign of the monic leading normalization.

## 5. Independent computation and exact dependencies

[check_single_depletion_independent.py](../computations/check_single_depletion_independent.py)
checks the pole residues from their exact full harmonic expressions,
including the zero node and cases with s>p. It separately generates exact
uniform-distribution cumulants and Bernoulli moments, checks their equality
and local integrality, and constructs fully coupled matrices with singular
exceptional blocks (including the zero block) to test the coefficient-loss
argument. The scalar margins are checked as rational identities.

The author's checker was inspected and rerun: all 31 admissible reduced
matrix cases pass and two fixed exceptional primes are correctly skipped.
Its synthetic divisions correctly extract diagonal remainders, and its
Newton trace calculation operates only on an exceptional block of dimension
less than p. It does not silently divide by p in that step. The full
characteristic polynomial is assembled using the proved block factorization.

External mathematical inputs retained are the established original Hermite
pole/moment functional, positive Gram determinant identities, the prime
number theorem, and the previously independently certified canonical monic
Jensen bound. All new local cancellation and normalization arguments above
have been independently checked. No claim concerning the subsequent
Möbius extension is included in this review.
