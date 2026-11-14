# The higher-log integral: a precise polynomial-decay obstruction

Status: independently derived all-n identities and asymptotic theorem below. This complements `forced-recurrence-obstruction.md`; it does not prove a general impossibility theorem for ζ(5). No historical novelty is claimed. All normalizations are explicit: the J_n here equals four times L_(n,2) in that companion note.

## 1. A natural higher-log kernel really does isolate ζ(5)

Keep the shifted Legendre polynomial

    L_n(x)=(1/n!)(d/dx)^n[x^n(1-x)^n]

and Zudilin's rational function

    R_n(t)=[(t-1)...(t-n)/(t(t+1)...(t+n))]².

The Mellin identity in `beukers-explained.md` gives

    J_n=∫_0^1∫_0^1 L_n(x)L_n(y)
                    (-log(xy))³/[6(1-xy)] dxdy
       =-(1/6) Σ_{t=1}∞ R_n'''(t).                     (1)

Every exchange is justified by boundedness of L_n and the integrable majorant (-log(xy))³/(1-xy). Its integral is 24ζ(5); after division by 6 it is 4ζ(5).

Writing R_n=Σ_k[a_k/(t+k)²+b_k/(t+k)], with a_k=binom(n,k)²binom(n+k,k)² and Σb_k=0, gives

    J_n=4A_n ζ(5)-V_n,
    A_n=Σ_k a_k,
    V_n=Σ_k[4a_k H_k^(5)+b_k H_k^(4)].                 (2)

Because a_k∈Z and D_n b_k∈Z, D_n^5 V_n∈Z. Thus a single ζ(5) value has been isolated with controlled denominators. That is a genuine positive feature of the construction, but it does not ensure sufficiently small linear forms.

## 2. The decisive asymptotic

**Theorem.** As n→∞,

    J_n = 2(H_(n-1)+H_(n+1))/[n²(n+1)²]
             + O(log n/n^8),                           (3)

and hence

    J_n ~ 4 log n/n^4.                                 (4)

In particular J_n is positive for all sufficiently large n, but it is not exponentially small. The obstruction is the decay rate, not merely lack of a positivity proof.

This theorem identifies the actual source of the error: the single former zero t=1 contributes the leading term. Terms at t=2,...,n and the entire tail t>n contribute only to the displayed much smaller remainder.

## 3. Derivatives at the former zeros

For 1≤j≤n write the unsquared rational function as

    (t-1)...(t-n)/(t(t+1)...(t+n))=(t-j)U_j(t).

Factorial products give

    U_j(j)=(-1)^(n-j) [(j-1)!]² (n-j)!/(n+j)!,
    U_j'(j)/U_j(j)=2H_(j-1)-H_(n-j)-H_(n+j).

Since R_n(t)=(t-j)²U_j(t)²,

    -(1/6)R_n'''(j)
      =2 U_j(j)² [H_(n-j)+H_(n+j)-2H_(j-1)].           (5)

At j=1 this becomes exactly the first term in (3). It is positive. For other j the bracket need not be positive, so a termwise positivity claim would be false: the sign changes near the right endpoint for large n.

To bound their total, define

    Q_j=binom(n,j)binom(n+j,j).

Then |U_j(j)|=1/(j²Q_j). The ratio

    Q_(j+1)/Q_j=(n-j)(n+j+1)/(j+1)²

crosses 1 at most once, from above to below. Therefore the minimum of Q_j over 2≤j≤n lies at j=2 or j=n. For n≥5,

    Q_n=binom(2n,n) ≥ Q_2=n(n-1)(n+1)(n+2)/4.

The n=5 case is 252≥210. Induction proves the inequality for larger n: the ratio of consecutive central binomial coefficients is 2(2n+1)/(n+1), while the ratio of consecutive Q_2 values is (n+3)/(n-1); the former exceeds the latter for n≥3.

It follows that Q_j≥c n^4 for 2≤j≤n, with an absolute c>0. The harmonic bracket in (5) is O(log n), uniformly in j. Since Σ_{j≥2}j^-4 converges,

    Σ_{j=2}^n |R_n'''(j)|=O(log n/n^8).                (6)

The n=1,...,4 cases do not affect the asymptotic statement.

## 4. A rigorous bound for the infinite derivative tail

Set N=n(n+1). For t>n, R_n(t)>0, and

    R_n(t)=t^-2 ∏_{j=1}^n[(1-j/t)/(1+j/t)]²
           ≤t^-2 exp(-2N/t),                          (7)

using log((1-u)/(1+u))≤-2u for 0≤u<1.

Split integer t>n into n<t<2n and t≥2n. On the first range, (7) gives R_n(t)≤t^-2e^(-n). Logarithmic differentiation of the original product shows

    |(log R_n)'|=O(log n),
    |(log R_n)''|+|(log R_n)'''|=O(1)

uniformly on those integer points: all nonzero denominator distances are positive integers, and their reciprocal-square/cube sums are bounded. Since

    R_n'''=R_n[(log R_n)'³
                 +3(log R_n)'(log R_n)''+(log R_n)'''],

the total over this finite range is exponentially small times a polynomial in log n, in particular O(n^-8).

On t≥2n use the cancellation in the first logarithmic derivative:

    (log R_n)'=-2/t+4Σ_{j=1}^n j/(t²-j²).

Since t²-j²≥3t²/4, differentiation gives the bounds

    |(log R_n)'| ≤ C(t^-1+N t^-2),
    |(log R_n)''|≤ C(t^-2+N t^-3),
    |(log R_n)'''|≤C(t^-3+N t^-4),

with an absolute constant C. Combining these with (7),

    |R_n'''(t)|≤C' exp(-2N/t)
          [t^-5+N t^-6+N²t^-7+N³t^-8].               (8)

For each fixed integer a≥2,

    Σ_{t≥1} t^-a exp(-2N/t)=O_a(N^(1-a)).              (9)

Indeed the function on the left is unimodal as a function of real t>0. Its sum is bounded by its integral plus twice its maximum. Substitution u=2N/t gives integral (a-2)!/(2N)^(a-1), and the maximum is O_a(N^-a). Applying (9) term by term in (8) yields

    Σ_{t>n}|R_n'''(t)|=O(N^-4)=O(n^-8).               (10)

Equations (5), (6) and (10) prove (3). Finally H_n=log n+O(1) proves (4). No unproved saddle-point claim or numerical fitting enters the argument.

## 5. Parameter deformation explains why the naive positive-integral argument breaks

For real h>-1 define

    K_n(h)=Σ_{t=1}∞ R_n(t+h)
          =∫∫ L_n(x)L_n(y)(xy)^h/(1-xy) dxdy.

The series representation shows K_n(h)≥0, since every term is a real square. It is analytic in h near 0, and its third derivative satisfies J_n=-K_n'''(0)/6. Positivity of a function does not constrain the sign or decay of its third derivative.

More concretely, the term t=1 is

    R_n(1+h)=h²/[n²(n+1)²]
        ·[1-2(H_(n-1)+H_(n+1))h+O(h²)].              (11)

At h=0 the double zero suppresses this term. Its first derivative is also zero, exactly what the ζ(3) proof needs. Its third derivative is not zero; it produces the leading term (3). Thus taking more derivatives does not simply preserve the small object used for ζ(3).

The proof of the positive Beukers triple integral also cannot be differentiated as if its clean h=0 form held unchanged for h≠0. In the Mellin deformation, x^h and y^h are present. Integration by parts differentiates these factors as well as the rational kernel, generating additional terms. Formula (11) is an exact, concrete witness of the missing contributions in any such naive manipulation.

## 6. A fixed-order extension

The same argument gives a broader diagnosis. For a fixed integer m≥2, use the normalization of `forced-recurrence-obstruction.md`:

    L_(n,m)=-(1/(2m)!) Σ_{t≥1} R_n^(2m-1)(t).

Then

    L_(n,m) ~ (4log n)^(2m-3) / [(2m)(2m-3)! n^4].    (12)

Here are the details needed beyond the m=2 case. Put q=2m-1. At t=1+h,

    R_n(1+h)=C_n h² exp(ell_1 h+ell_2 h²/2!+...+ell_q h^q/q!+...),
    C_n=1/[n²(n+1)²],
    ell_1=-2(H_(n-1)+H_(n+1))=-4log n+O(1).

For each fixed j≥2, ell_j=O_j(1), because it is a finite combination of reciprocal-j-power sums over nonzero integer distances. The coefficient of h^(q-2) in the exponential is consequently ell_1^(q-2)/(q-2)!+O_q((log n)^(q-4)) when q≥5; for q=3 it is exactly ell_1. Multiplying by -q!/(q+1)! gives the leading expression (12), with positive sign because q-2 is odd.

At zeros t=2,...,n the same product-derivative argument, together with the Q_j bound above, gives total O_q(n^-8(log n)^(q-2)). On t≥2n, repeated logarithmic differentiation gives |(log R_n)^(j)|≤C_j(t^-j+N t^(-j-1)). The product rule then bounds |R_n^(q)| by a fixed sum of terms N^a t^(-q-a-2)e^(-2N/t), 0≤a≤q. Formula (9) bounds their total by O_q(N^(-q-1))=O_q(n^(-4m)). The intermediate integer range is exponentially small as before. All remainder terms are negligible compared with the main term in (12).

Thus every fixed higher odd zeta value isolated by this unmodified higher-derivative construction has polynomial, rather than exponential, decay. The special success at m=1 comes from the double-zero/first-derivative match; it is not merely an accident of numerical optimization.

## 7. Implications for the project

One possible repair is to omit t=1,...,n before taking the derivative sum. That eliminates the polynomial contribution above, but changes the rational coefficient by explicit terms (5), whose factorial/binomial denominators need new analysis. This arithmetic cost is studied separately by the truncation agent; it must not be ignored.

Another repair is to increase the zero multiplicity of the rational function so that derivative order 2m-1 still vanishes at the initial integer zeros. That changes pole multiplicities as well, usually introducing several odd zeta values and new denominator costs. The theorem above explains exactly what must be repaired; it does not establish that either repair succeeds.

## Reproducible finite checks

Run `python codex/computations/check_log_kernel.py`. It evaluates the exact rational coefficient and then uses precision 3n+100 decimal digits for cancellation against ζ(2m+1). For n=5,10,20,40,80, the ratio J_n to its exact first-zero contribution decreases from approximately 1.000692 to 1.000000019. The normalized remainder (J_n−first-zero)n^8/log n stays of moderate size (about 1.69 to 3.45 on these samples), consistent with (3). Fixed-m leading asymptotics converge slowly because H_n−log n tends to a nonzero constant; finite data are supporting checks, not a substitute for the proof.
