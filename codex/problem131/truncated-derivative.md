# Truncating the higher derivative: exponential decay returns, new denominators appear

Status: independent construction and proofs developed in this session. Historical novelty unresolved. This gives a precise repair of the analytic failure in `forced-recurrence-obstruction.md`, together with a rigorous arithmetic obstruction to retaining the original denominator bound. It does not prove irrationality of ζ(5).

## The construction

Use Zudilin's P_n(t)=product_{j=1}^n(t-j)/product_{j=0}^n(t+j), R_n=P_n². Define

    L_n=-(1/24) sum_{t=1}^∞ R_n'''(t),
    C_n=-(1/24) sum_{j=1}^n R_n'''(j),
    E_n=L_n-C_n=-(1/24) sum_{t=n+1}^∞ R_n'''(t).

The complete L_n cannot decay exponentially, as proved in the forced-recurrence note. The problem is that a third derivative no longer vanishes at a double zero. Removing precisely those first n terms is therefore a natural repair.

Write the usual partial fractions R_n=sum_k[a_k/(t+k)²+b_k/(t+k)], and A_n=sum a_k. Then

    L_n=A_n ζ(5)-B_n,
    B_n=sum_k[a_k H_k^(5)+(b_k/4)H_k^(4)],
    E_n=A_n ζ(5)-(B_n+C_n).

The coefficient A_n is still integral. The finite correction C_n is rational and explicit.

## 1. Exact correction formula

For 1<=j<=n, direct differentiation at the simple zero of P_n gives

    P_n'(j)=(-1)^(n-j) (j-1)!²(n-j)!/(n+j)!,
    P_n''(j)=2P_n'(j)[2H_(j-1)-H_(n-j)-H_(n+j)].

Since R_n'''(j)=6P_n'(j)P_n''(j),

    C_n=(1/2) sum_{j=1}^n P_n'(j)²
                    [H_(n-j)+H_(n+j)-2H_(j-1)].              (1)

This exact formula makes the newly exposed factorial denominators visible. It also shows why the correction is not simply a polynomial in n and H_n.

## 2. General arithmetic bound

Let D_n=lcm(1,...,n), D_0=1. Then

    4D_n²D_(2n)³ (B_n+C_n) ∈ Z.                            (2)

**Proof.** The original partial-fraction theorem gives 4D_n^5B_n integral. It remains to bound C_n.

For any prime p, Legendre's factorial formula expresses the exponent of p in the denominator of P_n'(j) as a sum, over q=p^a, of

    floor((n+j)/q)-2floor((j-1)/q)-floor((n-j)/q).

For q<=n this summand is at most two. To verify this, write j=bq+r with 0<=r<q and n=aq+s. If r=0 the summand is two; if r>0 it is floor((s+r)/q)-floor((s-r)/q), which is 0,1,or2. For n<q<=2n the summand is at most one, because all factorial arguments except n+j are below q. For q>2n it is zero. Consequently

    D_n D_(2n) P_n'(j) ∈ Z.

The harmonic bracket in (1) has denominator dividing D_(2n). Squaring P_n', multiplying by that bracket, and accounting for 1/2 gives

    2D_n²D_(2n)³ C_n ∈ Z.

Because D_n divides D_(2n), the same expression with factor4 also clears 4D_n^5B_n, proving (2). QED.

This is an upper bound, not a claim that every prime exponent is attained. It is nevertheless substantially worse than the original D_n^5 bound: new primes between n and 2n may enter.

## 3. An infinite family forces new prime denominators

**Theorem.** Let p be any odd prime other than 17, and let n=(p+1)/2. Then

    v_p(C_n)=-3,
    v_p(B_n+C_n)=-3.                                      (3)

In particular, no fixed constant multiple of D_n^5 can clear the corrected rational coefficient for every n.

**Proof.** Here p=2n-1>n. In (1), the only indices j for which the factorial (n+j)! contains p are j=n-1,n. For all earlier j, P_n'(j) and the harmonic bracket are p-integral. At those last two indices, P_n'(j) has valuation -1, and exactly one term 1/p occurs in H_(n+j). The other harmonic numbers have no p in their denominator. Thus the harmonic bracket has valuation -1, with p times that bracket congruent to 1 modulo p.

Put u=pP_n'(n-1) and v=pP_n'(n), both p-adic units. Their exact ratio is

    v/u=-(n-1)²/(2n) ≡ -1/4 (mod p),

because n≡1/2 modulo p. Hence

    p³ C_n ≡ (u²+v²)/2 ≡ 17u²/32 (mod p).

For odd p≠17 this is nonzero. All other terms are p-integral, so v_p(C_n)=-3. Since 4D_n^5B_n is integral and p>n is odd, B_n is p-integral. It cannot cancel a p^-3 denominator, proving the second assertion. There are infinitely many such primes, so multiplying D_n^5 by one fixed constant cannot repair the defect. QED.

The exceptional prime 17 really cancels at the displayed leading order; no claim about its exact valuation is needed. This theorem proves a denominator obstruction for the construction, not a lower bound on the total exponential size of all denominators.

## 4. Exponential decay of the truncated sum

**Theorem.** For every fixed θ with 1/16<θ<1,

    E_n=O(θ^n).

Equivalently, limsup |E_n|^(1/n)<=1/16. This upper-bound argument alone does not give a sign. Section 5 strengthens it to an exact asymptotic and eventual positivity; E_1 is negative, so positivity for every n would be false.

**Step 1: derive the forced recurrence.** Let

    (T u)_n=(n+1)³u_(n+1)-(2n+1)(17n²+17n+5)u_n+n³u_(n-1),
    S_n(t)=s_n(t)R_n(t),
    s_n(t)=4(2n+1)[-2t²+t+(2n+1)²].

Zudilin's exact certificate is T R_n(t)=S_n(t+1)-S_n(t). Differentiate three times and sum from t=n+1 to infinity. The neighboring tails E_(n+1) and E_(n-1) start at n+2 and n respectively, so keeping the two missing/extra terms gives

    (T E)_n=g_n,
    24g_n=S_n'''(n+1)+(n+1)³R_(n+1)'''(n+1)
                            -n³R_(n-1)'''(n).             (4)

The upper boundary derivative vanishes. The signs in (4) include the negative normalization -1/24 in E_n.

**Step 2: the forcing is exponentially small.** The elementary central-binomial lower bound

    binom(2n,n)>=4^n/(2n+1)

follows because the largest of 2n+1 binomial coefficients is at least their average. Therefore

    P_n(n+1)=1/[(2n+1)binom(2n,n)]<=4^-n,
    |P_(n+1)'(n+1)|=1/[(2n+2)(2n+1)binom(2n,n)]<=4^-n/(2n+2),
    P_(n-1)(n)<=4^(-(n-1)).

At the nonzero evaluation points, each logarithmic derivative of P through order three is bounded by a fixed constant times n+1: it is a sum of at most 2n+1 reciprocal powers whose denominators have absolute value at least one. The standard derivative formulas for a product then bound the first three derivatives of R by its value times O((n+1)^3). At the root of P_(n+1), use R'''=6P'P'' and the analogous logarithmic derivative after dividing out the root. The polynomial s_n and its first three derivatives at n+1 are O((n+1)^3). Substitution into (4) gives the deliberately generous uniform bound

    |g_n|<=K(n+1)^6 16^-n                              (5)

for an absolute constant K. Only the exponential base matters here; sharpening the polynomial is unnecessary.

**Step 3: the tail grows at most polynomially.** For integer t>=n+1,

    |P_n(t)|<=1/t,

since every (t-j)/(t+j) lies in [0,1]. The same logarithmic-derivative estimates give

    |R_n'''(t)|<=K'(n+1)^3/t².

Summing t>=n+1 shows |E_n|=O((n+1)^2). This excludes a hidden exponentially growing homogeneous component.

**Step 4: an elementary discrete maximum principle.** Fix θ in (1/16,1). For h_n=θ^n,

    (T h)_n/[n³h_n] → θ+1/θ-34<0.

The last inequality holds because the smaller positive root of z+1/z=34 is α=17-12sqrt2<1/16, and θ lies between that root and 1. Also T(2^n)/(n³2^n)→2+1/2-34<0. Thus for all sufficiently large n, both sequences have strictly negative T-values. By (5), a sufficiently large constant M makes

    T(Mθ^n)<=-|g_n|

for all n beyond some fixed N, and M can also be chosen to dominate |E_N|/θ^N.

For ε>0 set h_n=Mθ^n+ε2^n. Then h_N>=|E_N|, T h_n<=-|g_n|, and E_n-h_n tends to -infinity because E_n is polynomially bounded. If E-h had a positive maximum at an interior index n>N, then

    T(E-h)_n <= [(n+1)³+n³-(2n+1)(17n²+17n+5)](E_n-h_n)<0,

because the neighboring values do not exceed the positive maximum and the bracket is negative. But T(E-h)=g-T h>=0, a contradiction. Hence E<=h. Applying the same argument to -E gives |E_n|<=Mθ^n+ε2^n. Let ε decrease to zero at each fixed n. We obtain |E_n|<=Mθ^n for all n>=N, proving the theorem.

This proof establishes exponential decay without making numerical assumptions about cancellation in the infinite sum.

## 5. Exact forcing and sharp asymptotic: the tail is eventually positive

The upper bound above can be sharpened to a complete asymptotic:

    E_n ~ (π/2) (log n)/n³ ·16^-n.                         (6)

In particular E_n>0 eventually, and lim |E_n|^(1/n)=1/16. Thus truncation repairs both exponential size and eventual nonvanishing. The arithmetic obstruction persists.

**Exact forcing simplification.** Set

    r_n=R_n(n+1)=1/[(2n+1)²binom(2n,n)²],
    h_n=2H_n-H_(2n+1).

Then (4) simplifies to

    g_n=r_n/(16n²) [351n²+504n+176
                   -2n(287n²+335n+100)h_n].              (7)

Here is a reproducible derivation rather than an unexplained symbolic output. At t=n+1 the first three logarithmic derivatives of P_n are

    x=2H_n-H_(2n+1),
    y=H_(2n+1)^(2)-2H_n^(2),
    z=4H_n^(3)-2H_(2n+1)^(3).

Thus R_n'''/R_n=J(x,y,z), where J=8x³+12xy+2z. At t=n for P_(n-1), these three quantities become x+d1,y+d2,z+d3, with

    d1=-3/(2n)+1/(2n+1),
    d2=7/(4n²)-1/(2n+1)²,
    d3=-15/(4n³)+2/(2n+1)³.

Moreover n³R_(n-1)(n)/r_n=4n(2n+1)²=s_n(n+1). The root term in (4), divided by r_n, is 3(n+1)x-3/2. Finally s_n'(n+1)=-4(2n+1)(4n+3) and s_n''=-16(2n+1). Therefore 24g_n/r_n equals

    s_n(n+1)[J(x,y,z)-J(x+d1,y+d2,z+d3)]
      +3s_n'(n+1)(4x²+2y)+6s_n''x+3(n+1)x-3/2.

Expanding cancels x³, x², xy, y, and z, leaving exactly (7). The generic symbolic identity is included in the check script.

**Forcing asymptotic.** The standard Stirling estimate binom(2n,n)~4^n/sqrt(πn), and H_n=log n+O(1), give

    r_n~π/(4n)·16^-n,
    h_n~log n,
    g_n~-(287π/32)log n·16^-n.                            (8)

This sharper asymptotic uses Stirling's formula; the earlier exponential bound did not require it.

**Rigorous recurrence squeeze.** Put β=1/16 and f_n=β^n log n/n³ for n>=2. Direct ratios f_(n±1)/f_n give

    (T f)_n~(β+β^-1-34)log n·β^n
            =-(287/16)log n·β^n.

Hence g_n/(T f)_n→π/2. Fix δ>0, and choose η strictly between α=17-12sqrt2 and β, for example η=1/20. For all sufficiently large n, T(η^n)<0 and

    T[(π/2+δ)f]_n < g_n < T[(π/2-δ)f]_n.

Choose M sufficiently large to arrange the boundary comparisons at the first such index N. The upper barrier

    (π/2+δ)f_n+Mη^n+ε2^n

has T-value at most g_n and dominates E at N and eventually at infinity. The maximum principle from the preceding section bounds E above by this barrier. The lower barrier

    (π/2-δ)f_n-Mη^n-ε2^n

has T-value at least g_n and lies below E at the two boundaries; the minimum version of the same argument bounds E below it. Let ε decrease to zero, divide by f_n, and let n→infinity. Because η<β, the Mη^n/f_n terms vanish. Then let δ decrease to zero. This proves (6).

The constants 287 and π/2 are forced by the endpoint cancellation and the recurrence's value at the decay base β, not fitted to numerical data.

## 6. What the repair achieves and what it does not

Truncation repairs one genuine analytic defect: the complete derivative sequence cannot decay exponentially, whereas the truncated one provably does. However, the arithmetic deteriorates from D_n^5 to a bound involving D_(2n), and theorem (3) proves that this deterioration cannot simply be ignored.

Even granting the prime number theorem, the upper bound in (2) has logarithmic growth about 8n, whereas the established decay bound approaches exp(-n log16), with log16≈2.773. These bounds do not yield an irrationality contradiction. Their mismatch does not prove that sharper cancellation in denominators is impossible. Nor does theorem (3) alone establish an exponential lower bound on the full denominator.

A useful next target would therefore be a new arithmetic identity for B_n+C_n or a different correction that preserves small prime support. Merely deleting the initial terms is analytically effective but arithmetically expensive. Section 5 now supplies eventual nonvanishing, but the arithmetic gap remains.

## Reproducibility

Run `python codex/computations/check_truncated_derivative.py`. It checks exact rational denominator bounds, the prime-band family, and the moving-boundary recurrence at small n. These checks supplement the all-n proofs. The exponential theorem is proved by inequalities and a maximum principle, not inferred from observed numerical values.
