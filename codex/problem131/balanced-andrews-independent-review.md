# Independent review: Andrews grid arithmetic and exponential normalization

Status: accepted. Source: [balanced-andrews-grid-arithmetic.md](balanced-andrews-grid-arithmetic.md).
This review includes the positive multisum, all-integer grid divisor, explicit
normalizer, and evaluation bounds. It does not claim a useful irrationality
inequality or historical novelty.

## The imported finite identity and its specialization

The root inspected the original page5 of Krattenthaler–Rivoal,
arXiv:math/0312148v2, equation(3.1), including a rendered page archived as
`sources/andrews-equation3-review.png`. This is Andrews' terminating identity,
not the subsequent nonterminating proposition with extra convergence hypotheses.
The polynomial note's substitution s=m, N=n, a=-2x-n,
b_i=c_i=1-x (i<=m), b_(m+1)=-x, c_(m+1)=1-x is correct.
The last denominator on the multiple-sum side becomes (1)_(k_m)=k_m!.
The increment parameters become -n-1, and the terminating parameter -n
supplies the final binomial. Their signs cancel because the increments sum
to k_m. Intermediate Pochhammer quotients are squared; the last has
(-x)_(k_m)(1-x)_(k_m), yielding precisely x/(x-k_m).
The prefactor times the zeroth hypergeometric summand is n!(x+1)_n^(2m).
Dividing by F(0) gives the stated normalized polynomial, with no missing n!.
Generic x avoids zero denominators. Finite rational continuation is legitimate;
at positive integer j the terms with k_m>=j vanish in the original Pochhammer
form, including k_m=j where cancellation of the displayed x/(x-k_m) must
be interpreted before specialization. The remaining terms are positive.

## Arithmetic lemma

After setting y=j-k, direct factorial cancellation gives equation(6).
At every prime-power level q the factorial valuation is
2 floor((n_0+y_0)/q)+[n_0<k_0]+[y_0=0]. It is nonnegative. If q divides
y+k, either y_0=0 or k_0=q-y_0; in the latter case one of the first two
terms supplies at least one. This cancels v_p(y+k) across all levels.
Thus the entire last block, rather than each of its fractions separately,
is integral. The remaining m-1 squared blocks cost d_n^(2m-2), because
v_p binom(n,k)<=floor(log_p n). This proves the all-index grid theorem.
The reduction to G_0=A/gcd(A,H) is ordinary integer divisibility and valid
whether or not this divisor is the exact fixed divisor.

## Full linear-form normalizer and its limitation

The monic trial has the correct pure-zeta coefficient (-1)^m by the prior
reviewed Lagrange pole identity. Substitution q(-j^2)=(A/H)u_j in its residues
shows that the proposed M clears each rational pole term separately. The
factor2j^2 cancels the two small pole corrections; the residual factor m
and the harmonic denominator are explicitly included. The factorial pair
has total argument2K and divides (2K)!. M itself is an integer since the
multinomial (4mn)!/(n!)^(4m) is integral. Bounding this multinomial by
(4m)^(4mn), then the remaining fixed number of factors, proves the stated
exponential upper bound for fixed m. The optional upper prime-band removal
uses the separately reviewed paired-residue theorem and the verified prime
support of all factors outside (2K)!.

The bound no longer has an n log n term, but its exponential constant has
not been shown to beat approximation decay. It is not a primitive denominator
formula. Nothing here proves a new irrationality result or an asymptotic
failure theorem for the primitive balanced approximants.

## Evaluation bounds and independent computation

The zero tuple gives the lower bound. The ratio Z_(k+1)/Z_k is exactly
(j-k-1)/(j+n-k); thus Z_k<=Z_0 and j/(j-k_m)<=j/(j-n) for j>n.
Summing the remaining increment binomials is Vandermonde and gives C.
This establishes both stated bounds with their precise ranges.

`check_andrews_arithmetic_independent.py` constructs F directly with fresh
Fraction raising operators, evaluates the multiple sum with dynamic
programming (not the author's tuple enumeration), and independently forms
every pole residue. All616 exact evaluation/bound comparisons and28 complete
normalizer checks pass for m=1,...,4, n=1,...,7. Each normalizer clears
individual rational terms and the exact zeta coefficient is (-1)^m.
These tests supplement the source specialization and all-index proof above.
