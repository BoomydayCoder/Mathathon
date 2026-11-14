# Where increasing the odd zeta index helps and hurts the arithmetic

Status: independently derived exact lemmas and sufficient normalization bounds. The input about Bernoulli denominators is the classical von Staudt–Clausen theorem, stated explicitly below. These estimates do not prove irrationality for all odd values. No novelty claim.

## Main finding: two opposing effects

Increasing s in the positive moment model has two opposing arithmetic effects:

* It cancels **more** large-prime Bernoulli denominators in the polynomial quotient.
* The direct denominator estimate for the harmonic pole values gets worse, by as much as two extra powers of lcm(1,...,K) for s→s+2.

This can be tested prime by prime. A useful search should optimize their interaction in the entire rational functional, rather than declaring every higher index uniformly harder.

## 1. Exact cancellation rule for polynomial moments

For odd s>=3 and e>=0, write m=2e+2. The functional from `odd-moment-generalization.md` has

    μ_s(t^e)=(-1)^e B_m * product_{a=1}^{s-2}(m+a)/(s-1)!
            =(-1)^e (B_m/m) binom(m+s-2,s-1).

These are factorial rearrangements, not estimates.

**Classical input (von Staudt–Clausen).** For positive even m, B_m + sum_{p-1|m} 1/p is an integer, where p runs over primes. Consequently v_p(B_m)=-1 when p-1 divides m, and v_p(B_m)>=0 otherwise.

**Exact rule.** If p>s-1 is prime, a denominator p occurs in μ_s(t^e) if and only if

    (p-1)|m  and no integer in [m+1,m+s-2] is divisible by p.

If it occurs, its exponent is exactly one. Proof: (s-1)! is a p-unit. In the first case the Bernoulli factor has valuation -1; the short product supplies a factor p exactly when the indicated interval contains a multiple of p. In the other case every factor is p-integral. QED.

When m=q(p-1), the cancellation condition is

    q mod p belongs to {1,...,s-2}.

For q<p this is particularly transparent: the first s-2 possible multiples m=(p-1),2(p-1),...,(s-2)(p-1) have their Bernoulli denominator canceled. The next multiple usually does not.

**Uniform corollary.** For every e<=E, the moment μ_s(t^e) is p-integral whenever

    p > max(s-1, 1+(2E+2)/(s-1)).

Indeed, if p-1|m, then q=m/(p-1)<s-1, hence 1<=q<=s-2<p. The preceding rule cancels p. All other cases were integral already.

This is an actual cutoff, not an asymptotic heuristic. In particular, moment denominators cannot contain an arbitrarily large Bernoulli prime once this short product is retained.

For the proposed proportions r=6, N~3K/40, h~37K/40, the largest polynomial quotient degree is E=rN+2h-2-K~1.3K. The large-prime cutoff is then approximately 1+2.6K/(s-1). At s=5 this is about 0.65K, and at s=7 about 0.433K. Thus, for large K, every polynomial quotient moment is integral at primes p>K for both s=5 and s=7. This does **not** imply integrality of the whole rational functional: simple-pole residue denominators remain.

The exact transition s→s+2 is

    μ_(s+2)(t^e)/μ_s(t^e)=(m+s-1)(m+s)/(s(s+1)).

It makes the tradeoff visible: two additional numerator factors cancel selected large primes; only primes dividing s(s+1) are added by this ratio's denominator. All moments are nonzero, so the quotient is defined.

## 2. Pole values have a different arithmetic cost

For j>=1,

    μ_X(1/(t+j²))=j^(s-1)X - j^(s-1)H_j^(s)-1/(s-1)+1/(2j).

Let a_p(j)=floor(log_p j). The harmonic sum has valuation at least -s a_p(j), hence the constant term has valuation at least

    min{ (s-1)v_p(j)-s a_p(j), -v_p(s-1), -v_p(2j) }.

The coefficient of X is integral. This lower bound is sufficient, not always sharp: cancellation between harmonic summands or between the three displayed terms may improve it.

Putting L_K=lcm(1,...,K), it follows immediately that

    (s-1)L_K^s μ_X(1/(t+j²)) belongs to Z[X], 1<=j<=K.

The factor 2 in 2j is absorbed because s-1 is even. This global bound worsens by L_K² when the index rises by two (apart from the fixed factors s-1 and s+1). Under the prime number theorem log L_K=K+o(K), this costs roughly exp(2K) per matrix entry and exp(2hK) in the entrywise determinant bound.

That last statement describes this **bound**, not the actual minimal denominator. It is precisely where an improvement might recover the higher-index analytic margin.

## 3. Entire rational functional: an explicit normalization

Set D_K(t)=product_{j=1}^K(t+j²), and take an integer polynomial A(t) of degree at most K+E. Divide by the monic integer polynomial D_K:

    A(t)/D_K(t)=Q(t)+sum_{j=1}^K c_j/(t+j²),
    Q in Z[t], degree Q<=E.

The residues are c_j=A(-j²)/D_K'(-j²), and the exact derivative formula is

    D_K'(-j²)=(-1)^(j-1) (K-j)!(K+j)!/(2j²).

To derive it, separate the products (l-j)(l+j), omit l=j, and simplify the factorials. The product (K-j)!(K+j)! divides (2K)!, since their quotient is an integer binomial coefficient. Therefore (2K)! c_j is an integer for every j. This bound deliberately ignores any helpful zero factors in A(-j²).

Define the exact moment-denominator integer

    T_(s,E)=lcm{denominator(μ_s(t^e)): 0<=e<=E},

with T_(s,-1)=1 when there is no polynomial quotient. It is computable using rational Bernoulli numbers, and the cancellation rule in Section 1 constrains its prime support.

**Normalization theorem.** With

    M=(2K)!(s-1)L_K^s T_(s,E),

we have M μ_X(A/D_K) in Z[X]. Proof: T clears each integer polynomial quotient coefficient's moment; (2K)! clears all residues; (s-1)L_K^s clears every pole value. Their product clears the sum. QED.

For a Hankel matrix with A=D_N^r t^(i+j), 0<=i,j<h, choose

    E=max(-1, rN+2h-2-K).

Then M^h det G(X) belongs to Z[X]. This gives a fully explicit sequence of integral polynomials and nonvanishing follows separately from the positive Gram representation. The only missing ingredient for irrationality is a sufficiently strong upper bound on their values.

A more economical sufficient factor is the lcm of the polynomial-part and pole-part factors, rather than their product. The product above keeps its proof transparent and overestimates prime overlap.

## 4. A coarse growth estimate, with its assumptions exposed

A uniform bound independent of s is

    T_(s,E) divides L_(2E+3)^2, E>=0.

Proof: use μ_s(t^e)=(-1)^e(B_m/m)binom(m+s-2,s-1). The Bernoulli denominator is squarefree and supported on primes <=m+1; it divides L_(m+1). Also m divides L_(m+1). Thus the denominator divides L_(m+1)^2, hence L_(2E+3)^2. This ignores the stronger large-prime cancellations proved above.

Using Stirling's formula and the prime number theorem, if E=cK+O(1), h=βK+O(1), with s,c,β fixed and c>=0, the explicit product choice satisfies

    log(M^h) <= 2β K² log K
                + β(2 log 2-2+s+4c)K² + o(K²).

For this formula use the coarse T bound, log L_x=x+o(x), and log(2K)!=2K log K+(2log2-2)K+o(K). The O(1) factors s-1 and degree rounding contribute only o(K²).

The K² log K term means an unnormalized determinant must be exceptionally small before this crude factor is useful. Determinant constructions normally insert factorial/basis normalizations precisely to cancel such structural growth. This estimate is a baseline bound, not a claim that optimal normalization has that cost.

**Effective sufficient condition.** If a candidate sequence with degree h has

    log det G(ζ(s)) + h log M <= -c_0 K²

for a fixed c_0>0 and all sufficiently large K, then the polynomial irrationality criterion proves ζ(s) irrational. More generally the left-hand side divided by h must tend to -infinity. All quantities in this condition refer to the same normalized polynomial; switching between differently normalized real and arithmetic estimates is invalid.

## 5. Falsifiable next mechanism

The promising question is not merely whether the integral stays positive: that is already proved for every odd s. It is whether residue cancellations in the *entire functional* save the extra ~2hK cost introduced by the harmonic bound when s rises by two.

A focused next lemma should partition primes by the residue classes of j and by the integer q=(2e+2)/(p-1). For polynomial terms the cancellation rule is exact. For pole terms the valuation bound in Section 2 supplies a testable target. One then asks whether a basis with zeros at selected pole classes forces cancellation between quotient and pole contributions before taking the determinant. Merely normalizing residues individually forfeits that possibility.

No conclusion about ζ(7) or all odd values follows from these estimates. They do isolate an explicit arithmetic contest and show that the polynomial quotient becomes better behaved at large primes as s increases, contrary to the naive expectation that every denominator cost must get worse.

## Exact computational check

`python codex/computations/check_odd_arithmetic.py` checks the exact Bernoulli cancellation rule for s=3,5,7,9 and several primes across q both below and above p; it also verifies the full normalization factor for all 33 previously saved Hankel determinants. All checks passed. A concrete switch is p=11,m=40: the moment has valuation -1 for s=3 and s=5, but valuation 0 for s=7 because m+4=44 supplies the canceled factor. These finite checks supplement the proofs; they do not replace the all-parameter arguments.
