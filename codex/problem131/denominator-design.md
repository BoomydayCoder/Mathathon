# Denominators as a design tool, and a smaller elementary lcm estimate

Status: theorems and proofs below independently developed in this session; their historical novelty has not been investigated. They explain and modestly streamline classical irrationality constructions. None proves a new higher odd-zeta irrationality result.

## 1. A period-six elementary lcm bound suffices for both ζ(2) and ζ(3)

Let D_n=lcm(1,...,n). Then

    log D_n <= c n+O(log²(n+1)),
    c=(log2)/2+3(log3)/4≈1.170532806781055.

This is weaker than the period-thirty estimate in `zudilin-explained.md`, but simpler: four factorial terms and a six-entry table replace five factorial terms and a thirty-entry table. The rate is still strictly adequate for both classical Beukers bounds.

**Proof.** For real x>=1 set ψ(x)=log D_floor(x), and

    L(x)=log(floor x)!-log(floor(x/2))!
         -2log(floor(x/3))!+log(floor(x/6))!,
    h(y)=floor y-floor(y/2)-2floor(y/3)+floor(y/6).

Counting the multiples of p,p²,... in each factorial yields

    L(x)=sum_{p prime,a>=1} h(x/p^a)log p.

Because 6-3-2·2+1=0, h(y+6)=h(y). On the six intervals [j,j+1), 0<=j<=5, the values are

    0,1,1,0,0,1.

Thus h>=0 everywhere, and h=1 on [1,3). Every prime power in (x/3,x] contributes to L(x), giving

    ψ(x)-ψ(x/3)<=L(x).

The elementary integral comparison for log t gives

    log(m!)=m log m-m+O(log(m+1)).

Substitute the four fixed dilations. Since 1-1/2-2/3+1/6=0, both the x log x and -x terms cancel. The surviving coefficient is

    L(x)=a x+O(log(x+1)),
    a=(log2)/2+2(log3)/3-(log6)/6
     =(log2)/3+(log3)/2.

Iterate ψ(x)<=L(x)+ψ(x/3) until the argument is less than one. The principal terms sum to (3/2)ax, while O(log x) error terms each contribute at most O(log x). Hence

    ψ(x)<=(3/2)a x+O(log²(x+1))=c x+O(log²(x+1)).

This proves the theorem using only factorial prime valuations, elementary integration, and a finite table. The O(log² n) error is o(n), which suffices for an exponential-rate contradiction. This theorem does not assert D_n<exp(cn) for every n.

### Why this weaker estimate is enough for ζ(3)

Let α=(sqrt2-1)^4. The needed strict inequality is 3c+log α<0. Raising its exponential to the fourth power, this is equivalent to

    2^6·3^9 < (1+sqrt2)^16.

The right side is 665857+470832sqrt2, while the left side is 1259712. Thus it suffices that 593855<470832sqrt2. Both sides are positive, and

    2·470832²-593855²=90701783423>0.

So the sign is certified by integer arithmetic. Numerically,

    3c+log α≈-0.0138959277350064.

Therefore any positive linear form with D_n³ denominator control and size O(poly(n)α^n) still produces the needed contradiction. The smaller negative margin affects quantitative convergence speed, not logical validity.

### Why it also suffices for ζ(2)

The classical ζ(2) integral has exponential factor

    β=((sqrt5-1)/2)^5=(5sqrt5-11)/2.

Its denominator exponent is two, so the requirement is 2c+log β<0. Squaring the exponential inequality reduces this to

    108β²<1,

or 6641<2970sqrt5. This follows from

    5·2970²-6641²=1619>0.

Thus the same six-entry arithmetic lemma supports both proofs. The ζ(2) integral and its denominator theorem are separate inputs; this paragraph verifies precisely the denominator-growth comparison, not a new integral identity.

### How the smaller factorial ratio was discovered

The search problem was: choose a short balanced combination of floors, nonnegative everywhere, equal to one on a long initial interval. Balanced means the coefficients of x cancel, so its factorial analogue has size exp(O(x)) rather than exp(O(x log x)). The usual central-binomial combination floor x-2floor(x/2) has period two and yields c=2log2≈1.386, too large. The combination floor x-floor(x/2)-floor(x/3)-floor(x/6) has period six but yields about 1.214, still too large. Reweighting the 1/3 term and changing the 1/6 sign gives the table above and c≈1.17053, just below the ζ(3) threshold ≈1.17516.

This is a controlled design search, not evidence of historical novelty or optimality among all elementary lcm arguments.

## 2. A pole-lattice theorem for building arithmetic approximations

Let n>=0 and c_0,...,c_n be integers. Set

    P(t)=sum_{k=0}^n c_k/(t+k),
    D=D_n.

For an integer m>=1, write its partial fractions as

    P(t)^m=sum_{k=0}^n sum_{j=1}^m a_{j,k}/(t+k)^j.

**Theorem.** For every k,j,

    D^(m-j) a_{j,k} ∈ Z.

For m>=2, sum_k a_{1,k}=0.

**Proof.** Near t=-k, put z=t+k. Then

    P(t)=c_k/z+sum_{a>=0} h_a z^a,
    h_a=(-1)^a sum_{l!=k} c_l/(l-k)^(a+1).

All nonzero differences |l-k| are at most n, so D^(a+1)h_a is an integer. To form a z^-j term in the mth power, choose p singular factors and m-p analytic factors. If their analytic exponents are a_1,...,a_(m-p), then

    sum a_i=p-j.

Their total denominator weight is

    sum(a_i+1)=p-j+m-p=m-j.

The multinomial coefficients and powers of c_k are integers, proving the assertion after summing all possibilities. The all-singular case j=m has weight zero. Finally P=O(t^-1), so P^m=O(t^-m); for m>=2 its coefficient of t^-1 vanishes, giving sum a_{1,k}=0. QED.

This theorem reveals what controls arithmetic: **differences between pole positions**, not the apparent size of roots in a logarithmic-derivative formula. A different finite integer pole set has the same theorem with D replaced by the lcm of its nonzero pairwise differences.

## 3. Derivatives explain both zeta isolation and denominator exponents

For an integer r>=1 define

    T_{m,r}=(-1)^r/r! sum_{t=1}^∞ (d/dt)^r[P(t)^m].

The differentiated series converges absolutely. Termwise partial fractions give

    T_{m,r}=sum_{j=1}^m binom(j+r-1,r)
                   sum_k a_{j,k}[ζ(j+r)-H_k^(j+r)].

Each rational constant term has denominator dividing D^(m+r): its coefficient uses at most D^(m-j), and the harmonic sum uses D^(j+r). Consequently D^(m+r) clears every coefficient and constant in this linear form. Some coefficients may already have much better denominators, as the pole-lattice theorem specifies.

For m>=2 the coefficient of ζ(r+1) is zero. In the particularly important case m=2,

    T_{2,r}=(r+1)(sum_k c_k²)ζ(r+2)-V,
    D^(r+2)V ∈ Z.

Thus **isolating a single higher zeta value algebraically is easy**: the square with r=3 gives a linear form in 1 and ζ(5). The difficult obligations are nonvanishing and sufficiently rapid decay after multiplication by D^5. An explanation that calls algebraic isolation alone the obstacle would miss this family.

For r=1 the theorem recovers Zudilin's arithmetic linear form in ζ(3), with D³ denominator control.

## 4. How to choose integer residues by interpolation

Let Q(t) be a polynomial of degree at most n, integer-valued at the integers -n,...,0. Define

    P(t)=n!Q(t)/[t(t+1)...(t+n)].

At t=-k its residue is

    c_k=(-1)^k binom(n,k) Q(-k),

which is an integer. Therefore every such Q is an admissible input to the pole-lattice theorem.

A particularly transparent sufficient choice is any integer-valued polynomial of degree <=n, expressed in the binomial basis binom(t,j) with integer coefficients. Another is

    Q(t)=binom(t-a-1,n),  a an integer,

which makes the numerator product_{j=1}^n(t-a-j). Choosing a=0 produces Zudilin's consecutive roots 1,...,n. All translated choices retain the basic denominator arithmetic; they need not retain good approximation or a short recurrence.

There is also a direct residue-design algorithm:

1. Fix the pole lattice 0,-1,...,-n to control differences and denominators.
2. Write P=sum c_k/(t+k), initially with rational unknowns c_k.
3. Require selected positive integer zeros r_i by imposing sum c_k/(r_i+k)=0.
4. Solve this finite rational linear system and scale to an integer residue vector.
5. Check how that scaling grows with n; arbitrary large scaling may destroy useful decay.
6. Only then investigate real-variable size, nonvanishing, and possible short recurrences.

If n distinct roots 1,...,n are prescribed, the numerator has degree at most n and is forced to be a scalar multiple of product(t-j). Thus the familiar binomial residues are the interpolation solution of natural zero conditions, rather than unexplained coefficients. Choosing their primitive integral normalization controls the arithmetic cost.

## 5. A precise warning about the naive ζ(5) derivative construction

Take the original P_n and R_n=P_n². For ζ(3), the summand -R_n'(t) vanishes at each integer t=1,...,n because R_n has double zeros there. For ζ(5), the analogous summand is -R_n'''(t)/6. Double zeros do not force a third derivative to vanish.

Indeed, for n>=1,

    P_n'(1)=(-1)^(n-1)/[n(n+1)],
    P_n''(1)=2P_n'(1)[-H_(n-1)-H_(n+1)],

so

    -R_n'''(1)/6
      =2[H_(n-1)+H_(n+1)]/[n²(n+1)²].

This is of order (log n)/n^4, not exponentially small. It is a concrete defect in the original initial-term removal mechanism. It does not alone prove that the entire signed series has that order: cancellation with other terms would have to be analyzed before claiming such an asymptotic. Numerical evidence for the full series is consistent with polynomial decay, but a rigorous global estimate is not asserted here.

More generally, raising P to the mth power and differentiating r times preserves all prescribed zeros only when m>r. For m=2,r=1 the arithmetic isolates ζ(3) and the zero mechanism works simultaneously. With m=2,r=3 one isolates ζ(5), but loses those zero terms. Increasing m repairs the zero multiplicity but introduces more zeta coefficients. Additional symmetries or transformations are then needed to eliminate them.

This is a useful design tension, not a theorem ruling out all elementary higher-zeta proofs.

## Reproducible checks and next steps

Run `python codex/computations/check_denominator_design.py` for the period-six table, exact exponential inequalities, finite pole-lattice examples, and the explicit first-root derivative. The general results are established by the proofs above; sampled examples merely check implementation conventions.

Next targets: quantify the full higher-derivative series without relying on cancellation guesses; classify root patterns with small primitive integer residues; and seek transformations that simultaneously improve zero multiplicity, zeta cancellation, and decay. These are concrete mathematical questions with separately testable arithmetic and analytic requirements.
