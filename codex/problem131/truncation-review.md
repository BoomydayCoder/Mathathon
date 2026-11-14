# Independent audit of truncation, and its exact leading asymptotic

Reviewed: `truncated-derivative.md`. Reviewer: independent integral/symmetric-square research agent.

**Audit result:** the moving-boundary recurrence, mixed denominator upper bound, prime-denominator obstruction and maximum-principle proof of exponential decay are valid. I reran `check_truncated_derivative.py`; all exact checks passed. The stronger asymptotic in Part II below is newly derived in this session and requests an additional independent review before being labelled independently checked. Historical novelty has not been assessed.

## I. Audit of the construction

### 1. Arithmetic upper bound

In the prime-power floor sum for P_n'(j), when j=bq+r and n=aq+s, the r=0 case is exactly 2; the r>0 case is floor((s+r)/q)-floor((s-r)/q), hence at most 2. For n<q≤2n, only the (n+j)! factorial can contribute, and it contributes at most 1. This proves D_nD_(2n)P_n'(j)∈Z. The harmonic bracket costs at most D_(2n), and the explicit factor 1/2 gives the stated bound for C_n. The common multiplier 4D_n²D_(2n)³ also clears B_n because D_n divides D_(2n). This is a valid upper bound, not a minimal-denominator formula.

### 2. New prime denominators

For odd prime p=2n-1, p>n and n+j≥p only when j=n-1 or n. At those indices P_n'(j) has valuation -1: the denominator factorial contains exactly one factor p and the numerator contains none. In the harmonic bracket only H_(n+j) contains a 1/p term. Therefore its valuation is exactly -1 and multiplying it by p gives residue 1.

The ratio of the two scaled derivatives is indeed

    pP_n'(n) / [pP_n'(n-1)] = -(n-1)²/(2n) ≡ -1/4 mod p.

Their squared contributions yield 17u²/32 modulo p, so v_p(C_n)=-3 unless p=17. Every other correction term is p-integral; B_n is also p-integral from its old denominator bound. Thus no cancellation with B_n can remove that denominator. The conclusion that no fixed multiple of D_n^5 clears all corrected coefficients is valid, since infinitely many primes of the form p=2n-1 occur when n is chosen from the primes, and a fixed integer has finitely many prime divisors.

This is not an exponential lower bound on the complete denominator. The original note correctly distinguishes the two assertions.

### 3. Moving-boundary recurrence

For E_n=-(1/24)Σ_{t≥n+1}R_n'''(t), summing the certificate at the central lower limit n+1 gives

    Σ_{t≥n+1} R_(n+1)'''(t)=-24E_(n+1)+R_(n+1)'''(n+1),
    Σ_{t≥n+1} R_(n-1)'''(t)=-24E_(n-1)-R_(n-1)'''(n).

These show that both endpoint signs in the proposed forcing term are correct:

    24g_n=S_n'''(n+1)+(n+1)³R_(n+1)'''(n+1)-n³R_(n-1)'''(n).

The upper boundary vanishes. The elementary factorial endpoint bounds and rough derivative estimates justify g_n=O(n^6 16^-n); their exact polynomial power is unimportant for the initial exponential theorem.

### 4. Maximum principle

For θ∈(1/16,1), θ+θ^(-1)-34<0, so T(θ^n)<0 eventually; T(2^n)<0 eventually as well. The exponentially decaying forcing can be dominated by -T(Mθ^n). The barrier h_n=Mθ^n+ε2^n dominates the polynomially bounded E_n at infinity and can dominate its value at the left boundary.

At a positive interior maximum of v=E-h, the positive off-diagonal coefficients of T imply

    (Tv)_n≤[(n+1)³+n³-(2n+1)(17n²+17n+5)]v_n<0.

But Tv=g-Th≥0, contradiction. The same argument applies to -E. Sending ε↓0 for each fixed n is legitimate because M is chosen independently of ε. This proves the stated exponential upper bound. No unproved claim about the sign of E_n is needed in that argument.

## II. A sharper asymptotic and eventual positivity

**Theorem (derived here).** With the same E_n,

    E_n ~ (π/2) (log n)/n³ ·16^-n.                    (1)

Consequently E_n>0 for all sufficiently large n and

    lim |E_n|^(1/n)=1/16.

This supplies the nonvanishing that was left open in the initial truncation note. It does not resolve its denominator-growth problem.

### 1. Simplify the exact forcing before estimating it

Let

    p_n=P_n(n+1)=n!²/(2n+1)!,
    h_n=2H_n-H_(2n+1),
    u_n=H_(2n+1)^(2)-2H_n^(2),
    v_n=4H_n^(3)-2H_(2n+1)^(3).

Thus h_n,u_n,v_n are the first three logarithmic derivatives of P_n at n+1. At that point,

    R_n=p_n²,
    R_n'=2h_n p_n²,
    R_n''=(4h_n²+2u_n)p_n²,
    R_n'''=(8h_n³+12h_n u_n+2v_n)p_n².

The multiplier s_n has values

    s_n(n+1)=4n(2n+1)²,
    s_n'(n+1)=-4(2n+1)(4n+3),
    s_n''(n+1)=-16(2n+1).

At the zero of P_(n+1),

    P_(n+1)'(n+1)=p_n/(2n+2),
    R_(n+1)'''(n+1)
       =12[p_n/(2n+2)]²[h_n-1/(2n+2)].

At the neighboring nonzero endpoint,

    P_(n-1)(n)=p_n·2(2n+1)/n,
    h_-=h_n-3/(2n)+1/(2n+1),
    u_-=u_n+7/(4n²)-1/(2n+1)²,
    v_-=v_n-15/(4n³)+2/(2n+1)³,
    R_(n-1)'''(n)
       =[p_n·2(2n+1)/n]²[8h_-³+12h_-u_-+2v_-].

Substitute these expressions into the reviewed recurrence forcing. Cubic powers of h, all u,v terms and the apparent higher powers in n cancel. The exact result is

    g_n=p_n² {
       -[(287n+335)/8+25/(2n)] h_n
       +351/16+63/(2n)+11/n² }.                       (2)

The accompanying script proves this rational identity symbolically with n,h,u,v indeterminate, and additionally checks the derivative definition at small integer n. Thus this cancellation is not inferred from numeric samples.

### 2. Extract the forcing asymptotic

The standard central-binomial asymptotic (from Stirling's formula or Wallis's product) is

    binom(2n,n) ~ 4^n/sqrt(πn).

Since p_n=1/[(2n+1)binom(2n,n)],

    p_n² ~ π/(4n) ·16^-n.

Also h_n=log n+γ-log2+o(1), so (2) yields

    g_n ~ -(287π/32) log n ·16^-n.                   (3)

Stirling's leading constant is an explicit prerequisite for the exact π/2 constant in (1). The earlier exponential bound remains elementary without it.

### 3. A candidate particular solution

Define for n≥2

    b_n=(π/2)(log n)/n³ ·16^-n.

Its neighboring ratios satisfy b_(n+1)/b_n→1/16 and b_(n-1)/b_n→16. Therefore

    (Tb)_n/(n³b_n)→1/16+16-34=-287/16.

Together with the definition of b_n, this gives

    (Tb)_n ~ -(287π/32) log n ·16^-n ~ g_n.          (4)

In particular Tb<0 eventually, and the residual r_n=g_n-(Tb)_n satisfies r_n=o(-Tb_n).

Matching leading terms alone would not establish (1), because a recurrence also has homogeneous solutions. The following comparison controls those possible additions.

### 4. Prove the relative error tends to zero

Fix θ_0 with α<θ_0<1/16, where α=17-12sqrt2. Then T(θ_0^n)<0 for all sufficiently large n. Put w_n=E_n-b_n, so Tw=r.

For any ε>0, choose N sufficiently large that

    |r_n|≤ε[-(Tb)_n] for n>N,
    (Tb)_n<0,   T(θ_0^n)<0.

Choose M≥0 with εb_N+Mθ_0^N≥|w_N|. The positive barrier

    h_n=εb_n+Mθ_0^n

satisfies Th≤-|r| for n>N. The already-proved exponential bound on E implies w_n→0, while h_n→0. If w-h had a positive value, its positive supremum would be attained at a finite index: a sequence tending to zero cannot approach a strictly positive supremum only at infinity. The left boundary is nonpositive, so such a maximum would be interior, contradicting the same maximum principle as above. Applying it also to -w gives

    |E_n-b_n|≤εb_n+Mθ_0^n   (n≥N).

Since θ_0<1/16, θ_0^n/b_n→0. Divide by b_n and take limsup to get at most ε. As ε is arbitrary, E_n/b_n→1. This proves (1), including eventual positivity.

## III. Relation to the full higher-log kernel

The companion `log-kernel-extension.md` gives the full sequence L_n asymptotic log n/n^4 (its J_n is 4L_n). The finite correction C_n contains that entire polynomial contribution; subtraction leaves the exponentially small E_n described above. Thus the construction exhibits a clean analytic/arithmetic tradeoff:

- keeping the initial double-zero terms preserves the old small-prime denominator support but leaves polynomial decay;
- deleting them restores exponential decay and eventual positivity, but exposes new primes up to 2n in the denominator.

This is a verified obstruction for this specific repair, not a proof that a different arithmetic reorganization cannot succeed.

## Reproducibility

Run `python codex/computations/check_truncation_asymptotic.py`. It checks the exact forcing simplification and evaluates E_n from exact rational coefficients at precision 4n+100 decimal digits. Sample ratios of E_n to the predicted leading term are approximately 0.7166, 0.9365, 0.9596, 0.9680, 0.9735 for n=5,10,20,40,80. Logarithmic corrections make convergence slow. The comparison argument, not this finite table, proves the limit.
