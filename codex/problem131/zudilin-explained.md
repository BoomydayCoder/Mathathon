# Zudilin's elementary proof, explained and checked

**Scope and status.** This explains Wadim Zudilin, *An elementary proof of Apéry's theorem*, arXiv:math/0202159v1 (17 February 2002), retrieved directly from https://arxiv.org/pdf/math/0202159. The theorem is that ζ(3) is irrational. This is an exposition of an established proof, with a small independently derived improvement to its final tail estimate. No historical novelty is claimed. This file concerns the user's specified paper, not the still-unretrieved official competition statement. The companion integral exposition is a separate task.

All identities below have proofs or explicit verification instructions. Both telescoping certificates were independently checked symbolically in SymPy in this session; this is an exact rational identity check, not numerical sampling. A self-contained elementary denominator-growth argument is supplied at the end so the proof does not rely on the paper's unexplained assertion D_n<3^n.

## 1. Start from the contradiction we need

Suppose a real number ξ were p/q with integers p,q and q>0. If we could construct

    F_n = u_n ξ - v_n,
    u_n integral, D_n^3 v_n integral,
    0 < D_n^3 F_n → 0,

then q D_n^3 F_n would be a positive integer tending to zero. This is impossible. Here D_n=lcm(1,...,n), with D_0=1.

There are three logically independent requirements:

1. Arithmetic: coefficients have controlled denominators.
2. Analysis: the linear form is exponentially small.
3. Nonvanishing: it is not identically zero or zero along the chosen sequence.

A central insight of this proof is to use **two representations of the same number**. One makes arithmetic transparent; the other makes positivity and decay transparent. A recurrence proves that the representations agree. Trying to force every property out of a single formula creates unnecessary difficulty.

Set

    P_n(t) = product_{j=1}^n(t-j) / product_{j=0}^n(t+j),
    R_n(t) = P_n(t)^2,
    F_n = -sum_{t=1}^∞ R_n'(t).

A derivative of a double pole gives a triple pole, hence a ζ(3) sum. The zeroes at 1,...,n remove the beginning of the summation. The equally spaced poles make residues binomial coefficients. The extra denominator factor makes P_n=O(t^-1), and R_n=O(t^-2), which will cancel an unwanted ζ(2) term.

These are motivations for the design, not a claim that the construction is uniquely forced.

## 2. Lemma 1: arithmetic by squaring simple partial fractions

Define

    c_k = (-1)^(n-k) binom(n,k) binom(n+k,k),  0<=k<=n.

The partial fraction decomposition is

    P_n(t) = sum_{k=0}^n c_k/(t+k).

To prove it, multiply by t+k and substitute t=-k. The numerator is (-1)^n(n+k)!/k!, while the remaining denominator is (-1)^k k!(n-k)!. Their quotient is c_k. Subtracting these principal parts leaves a rational function without poles, hence a polynomial; it tends to zero, hence is zero.

Now square. The diagonal terms give c_k²/(t+k)². The off-diagonal terms are simplified using

    1/((t+k)(t+l)) = [1/(t+k)-1/(t+l)]/(l-k).

Thus

    R_n(t)=sum_{k=0}^n [a_k/(t+k)^2+b_k/(t+k)],
    a_k=c_k²,
    b_k=2 c_k sum_{0<=l<=n, l!=k} c_l/(l-k).

Since every nonzero |l-k| is at most n,

    a_k ∈ Z,       D_n b_k ∈ Z.

**Why this matters.** Direct logarithmic differentiation also gives

    b_k=2a_k(2H_k-H_{n+k}-H_{n-k}),

but that expression misleadingly exposes denominators up to 2n. Squaring the simple-pole expansion reveals that all denominators are differences of pole positions, hence at most n. This closes the arithmetic gap left open in the earlier auxiliary-function notebook.

Comparing coefficients of 1/t at infinity gives sum_k b_k=0, since R_n=O(t^-2). No residue theorem is necessary.

Differentiation and summation now yield

    F_n = sum_k [2a_k(ζ(3)-H_k^(3))+b_k(ζ(2)-H_k^(2))]
        = u_n ζ(3)-v_n,
    u_n = 2 sum_k a_k,
    v_n = 2 sum_k a_k H_k^(3)+sum_k b_k H_k^(2).

Here H_k^(j)=sum_{l=1}^k l^-j, and H_0^(j)=0. All series converge absolutely, and the k-sum is finite, so the rearrangement is legitimate. The ζ(2) coefficient vanishes because sum b_k=0.

Finally, D_n³ H_k^(3) and D_n² H_k^(2) are integers for k<=n. Consequently u_n and D_n³v_n are integers.

The first cases are

    F_0=2ζ(3),
    F_1=10ζ(3)-12.

For n=1, use a_0=1,a_1=4,b_0=-4,b_1=4. Notice the normalization: the earlier notebook used A_n=u_n/2 and B_n=v_n/2.

## 3. Lemma 2: why a recurrence is a reasonable thing to seek

The ratios of neighboring rational functions are simple:

    R_{n+1}(t)/R_n(t) = ((t-n-1)/(t+n+1))²,
    R_{n-1}(t)/R_n(t) = ((t+n)/(t-n))²,
    R_n(t+1)/R_n(t) = t⁴/((t-n)²(t+n+1)²).

This suggests looking for a linear combination of neighboring n-values that is a difference in t:

    a(n) R_{n+1}(t)+b(n) R_n(t)+c(n) R_{n-1}(t)
      = S_n(t+1)-S_n(t).

Such a difference disappears after summation, apart from boundary terms. This is creative telescoping. Its discovery can be algorithmic; its verification needs only polynomial algebra.

Zudilin's certificate is

    s_n(t)=4(2n+1)(-2t²+t+(2n+1)²),
    S_n(t)=s_n(t)R_n(t).

It gives the rational identity

    (n+1)³R_{n+1}(t) -(2n+1)(17n²+17n+5)R_n(t)+n³R_{n-1}(t)
       =S_n(t+1)-S_n(t).                         (T1)

Divide by R_n(t), insert the three ratios, and clear the denominator (t-n)²(t+n+1)². Both sides are now polynomials in n,t; their difference is the zero polynomial. This was checked exactly in SymPy. Apparent exceptions at zeroes of R_n are removed by the original rational identity; one must not evaluate a divided expression there without cancellation.

**Where do the strange coefficients come from?** They solve a compatibility problem: make the poles of neighboring n-values cancel into a t-difference with a low-degree multiplier. They are not selected merely to fit several terms of a numerical sequence. A companion derivation in [simplification.md](simplification.md) matches the double poles and the behavior at infinity to recover the quadratic multiplier and recurrence coefficients. That explains the discovery without requiring Zeilberger’s algorithm as a black box.

## 4. Lemma 3: boundary terms turn the certificate into a recurrence

Differentiate (T1), sum from t=1 to T, then let T→∞. Since R_n'=O(t^-3) and S_n'=O(t^-2), convergence and the upper boundary are controlled. The result is

    (n+1)³F_{n+1} -(2n+1)(17n²+17n+5)F_n+n³F_{n-1}=S_n'(1).

For n>=1, R_n has a double zero at t=1. Multiplication by the polynomial s_n preserves a zero of order at least two, so S_n'(1)=0. Thus

    (n+1)³F_{n+1} -(34n³+51n²+27n+5)F_n+n³F_{n-1}=0.      (A)

This is Apéry's recurrence.

The source's printed equation (7) has n³u_n as its last term. It must read n³u_{n-1}: the preceding certificate, the subsequent proof, and the initial values all establish that correction. Also, the initial conclusion is a recurrence for F_n. One should not infer separate recurrences for rational coefficients by comparing 1 and ζ(3) before irrationality has been proved. Zudilin's footnote gives a separate pole-coefficient argument for u_n, after which the recurrence for v_n follows; that extra assertion is unnecessary for the present irrationality proof.

## 5. The second representation: positivity is built in

Define Ball's rational function

    B_n(t)=(n!)²(2t+n)
       product_{j=1}^n(t-j) product_{j=n+1}^{2n}(t+j)
       / product_{j=0}^n(t+j)^4,
    G_n=sum_{t=1}^∞ B_n(t).

(The letter B here denotes a rational function, not the earlier notebook's rational coefficient B_n.)

For integers 1<=t<=n, B_n(t)=0. For t>n every factor is positive. Its degree gives B_n(t)=O(t^(-2n-3)), so the sum converges and G_n>0.

The reflected numerator factors and central factor 2t+n are characteristic of well-poised hypergeometric constructions. They furnish a function with symmetries that cancel unwanted zeta values. At this point, however, neither that terminology nor a partial-fraction expansion is needed: positivity follows directly from the displayed products.

The intended division of labor is now visible. F_n has good denominators; G_n has good signs. The remaining tasks are to bound G_n and prove F_n=G_n.

## 6. Lemma 4: why the small number (sqrt(2)-1)^4 appears

Write α=(sqrt(2)-1)^4=17-12sqrt(2), approximately 0.02943725.

For positive integers m,n, elementary bounds (1+1/j)^j<e<(1+1/j)^(j+1), multiplied over j=m,...,m+n-1, imply

    e^-n (m+n)^(m+n-1)/m^(m-1)
      < m(m+1)...(m+n-1)
      < e^-n (m+n)^(m+n)/m^m.                (P)

To see the cancellation explicitly, rewrite the single-j inequalities as

    (j+1)^j/j^(j-1) < e j < (j+1)^(j+1)/j^j,

then multiply: the intermediate powers telescope. This avoids invoking Stirling's formula.

Apply (P) to the three rising products in

    B_n(t)(t+n)^5/[(2t+n)(t+2n)]
       =(n!)² [(t-n)...(t-1)] [(t+n)...(t+2n-1)]
          /[t...(t+n-1)]^4.

Upper-bound the numerator products and lower-bound the denominator. For n>=1 and integer t>=n+1 the result is

    B_n(t) t⁴(t+n)/[(2t+n)(t+2n)(n+1)²]
      < (n+1)^(2n) t^(5t)(t+2n)^(t+2n)
           /[(t-n)^(t-n)(t+n)^(5(t+n))]
      = (1+1/n)^(2n) exp(n f(t/n)),

where

    f(τ)=5τ log τ+(τ+2)log(τ+2)
          -(τ-1)log(τ-1)-5(τ+1)log(τ+1),  τ>1.

The scaling τ=t/n is the crucial step: all exponential dependence is isolated into exp(n f(τ)). The largest value of f controls the whole series, just as a saddle point controls a Laplace integral.

Differentiation gives

    f'(τ)=log[τ⁵(τ+2)/((τ-1)(τ+1)⁵)].

The numerator minus denominator inside this ratio factors as

    -(τ+1/2)[2(τ+1/2)^4-5(τ+1/2)²-7/8].

For τ>1 there is exactly one zero, at

    τ_0=-1/2+sqrt(5/4+sqrt(2)).

The derivative is positive before this point and negative after it; hence it is the global maximum. At a stationary point one can cancel the large logarithmic terms by using f-τ f':

    exp(f(τ_0))=(τ_0+2)²(τ_0-1)/(τ_0+1)^5=α.

For an exact check of the last equality, τ_0²+τ_0=1+sqrt(2). Reduce

    (τ+2)²(τ-1)-(sqrt(2)-1)^4(τ+1)^5

modulo τ²+τ-(1+sqrt(2)); the remainder is zero. This algebraic reduction was checked independently.

Since (1+1/n)^(2n)<e², we have established

    0<B_n(t)<e²(n+1)² α^n (2t+n)(t+2n)/[t⁴(t+n)]   (t>n).   (B)

**A simpler final estimate.** Instead of the paper's estimates involving three zeta values, observe that for t>n,

    (2t+n)/(t+n)<=2,       t+2n<=3t.

The rational factor in (B) is therefore at most 6/t³. By monotonicity,

    sum_{t=n+1}^∞ t^-3 <= integral_n^∞ x^-3 dx=1/(2n²).

Thus, for n>=1,

    0<G_n<3e² ((n+1)/n)² α^n <=12e² α^n <96 α^n.       (B')

This removes both the polynomial prefactor and the final auxiliary-zeta calculations. It is a modest simplification of an established argument, not a claimed new irrationality proof. For n=0, G_0=2ζ(3), handled separately.

The source states the weaker bound 20(n+1)^4 α^n. Its t/n calculation is only valid for n>=1; n=0 must indeed be checked separately. Likewise the strict inequality D_0<3^0 is false, but the final contradiction only uses large n.

## 7. Lemmas 5 and 6: a second certificate, the same recurrence

Define

    Q_n(t)= -t^6-(8n-1)t^5+(4n²+27n+5)t^4
       +2n(67n²+71n+15)t³
       +(358n^4+339n³+76n²-7n-3)t²
       +(384n^5+396n^4+97n³-29n²-17n-2)t
       +n(153n^5+183n^4+50n³-30n²-22n-4),
    C_n(t)=B_n(t) Q_n(t)/[(2t+n)(t+2n-1)(t+2n)].

Zudilin's second certificate is

    (n+1)³ B_{n+1}(t)-(2n+1)(17n²+17n+5)B_n(t)+n³B_{n-1}(t)
       =C_n(t+1)-C_n(t).                         (T2)

A readable exact verification uses the following ratios:

    B_{n+1}/B_n=(n+1)²(2t+n+1)(t-n-1)(t+2n+1)(t+2n+2)
                  /[(2t+n)(t+n+1)^5],
    B_{n-1}/B_n=(2t+n-1)(t+n)^5
                  /[n²(2t+n)(t-n)(t+2n-1)(t+2n)],
    B_n(t+1)/B_n(t)=(2t+n+2)t^5(t+2n+1)
                  /[(2t+n)(t-n)(t+n+1)^5].

Divide (T2) by B_n, substitute these ratios and the multiplier defining C_n, and clear denominators. The resulting polynomial is zero. This was checked symbolically with n,t indeterminate, not by checking finitely many n.

Now sum (T2). For n>=1, C_n(t)=O(t^(-2n)), hence tends to zero at infinity. At t=1 its displayed denominator is nonzero and B_n(1)=0, so C_n(1)=0. Therefore G_n satisfies the same recurrence (A).

**Why this large polynomial is tolerable.** It is a certificate, analogous to a witness: its discovery may be difficult but its verification is finite algebra. For a conceptual presentation, its role belongs in the main argument and its expanded coefficients can go in an appendix. That shortens the exposition, but does not eliminate the mathematical verification.

## 8. Lemma 7: recurrence uniqueness identifies the representations

A second-order recurrence with nonzero leading coefficient determines every later term from two initial values. Here (n+1)³ never vanishes for n>=1.

For n=0, B_0(t)=2/t³, giving G_0=2ζ(3). For n=1,

    B_1(t)=-2/t⁴+2/(t+1)⁴+5/t³+5/(t+1)³-5/t²+5/(t+1)².

Summing, the fourth-power terms give -2, the third-power terms give 10ζ(3)-5, and the second-power terms give -5. Thus G_1=10ζ(3)-12.

These are precisely F_0,F_1. Induction through (A) proves F_n=G_n for every n>=0. No assumption about rational independence of zeta values has been used.

This is the recurrence's decisive conceptual role: it certifies equality between two different constructions, transferring arithmetic from one to the positivity and smallness of the other.

## 9. The final contradiction

Assume ζ(3)=p/q, with p,q positive integers. Then

    N_n=qD_n³F_n=D_n³u_n p-D_n³v_n q

is an integer. The equality F_n=G_n and (B') show

    0<N_n<96q D_n³ α^n  (n>=1).

Any denominator estimate log D_n<=c n+o(n), with

    3c+log α<0, equivalently c<(4/3)log(1+sqrt(2))≈1.17516478,

makes the right-hand side tend to zero. The paper uses D_n<3^n, with c=log3≈1.09861229. Its product 27α≈0.79480579 is less than 1.

For a self-contained substitute, the appendix below proves c≈1.10555043, still sufficient. Then eventually 0<N_n<1, a contradiction. Hence ζ(3) is irrational.

## 10. Elementary denominator growth, without the prime number theorem

This appendix supplies a sufficient estimate rather than proving the exact sharper inequality D_n<3^n quoted without proof in the paper.

For real x>=1, let ψ(x)=log lcm(1,...,floor x). Unique prime factorization gives

    ψ(x)=sum_{p prime, a>=1, p^a<=x} log p.

Let

    L(x)=log(floor x)!-log(floor(x/2))!-log(floor(x/3))!
         -log(floor(x/5))!+log(floor(x/30))!,
    w(y)=floor y-floor(y/2)-floor(y/3)-floor(y/5)+floor(y/30).

Legendre's counting formula, obtained by counting multiples of each p^a in a factorial, yields

    L(x)=sum_{p,a>=1} w(x/p^a) log p.

The function w is constant between consecutive integers and periodic with period 30. Its values at 0,...,29 are

    0,1,1,1,1,1,0,1,1,1,0,1,0,1,1,0,0,1,0,1,0,0,0,1,0,0,0,0,0,1.

Thus w is nonnegative everywhere and equals 1 on [1,6). All prime powers in (x/6,x] contribute log p to L(x), so

    ψ(x)-ψ(x/6)<=L(x).

The integral comparison for the increasing function log t gives

    log(m!)=m log m-m+O(log(m+1)).

Replacing m by floor(x/d), for any one of the fixed d=1,2,3,5,30, changes the principal expression by O(log(x+1)). Consequently

    L(x)=c_0 x+O(log(x+1)),
    c_0=(log2)/2+(log3)/3+(log5)/5-(log30)/30
       ≈0.9212920229340907.

The x log x and -x coefficients cancel because 1-1/2-1/3-1/5+1/30=0. Iterating ψ(x)<=L(x)+ψ(x/6) until x/6^j<1 gives

    ψ(x)<=(6/5)c_0 x+O((log(x+1))²).

Indeed the principal terms form a geometric series and there are O(log x) error terms, each O(log x). Therefore

    log D_n<=1.105550427521 n+O((log(n+1))²).

The exponential rate in the final integer bound is at most

    3(6/5)c_0+log α≈-0.208843065515445<0.

This establishes the needed decay by elementary factorial estimates and a finite table. For a wholly exact certification of the strict inequality, multiply the exponent by 25 and exponentiate. The desired inequality becomes 2^42 3^27 5^15 < (1+sqrt(2))^100. Since 1+sqrt(2)>12/5, it suffices to check the integer inequality 5^115 < 2^158 3^73, which holds by direct integer arithmetic. This comparison was independently evaluated using Python integers. Thus the numerical logarithm approximations are dispensable. This appendix is an elementary Chebyshev-style argument; no historical novelty is asserted.

## 11. Proof dependencies and what remains to explain

    simple-pole expansion → denominator bound ─────────────┐
                                                         │
    certificate T1 → recurrence for F ──┐                 │
                                       ├→ F=G ───────────┼→ positive small integer → contradiction
    certificate T2 → recurrence for G ──┘                 │
    two initial values ─────────────────┘                 │
    positive Ball series + saddle estimate ────────────────┤
    elementary lcm growth ─────────────────────────────────┘

**Which dependencies are indispensable?** Every version needs arithmetic control, nonvanishing, and decay faster than denominator growth. The recurrence is indispensable to Zudilin’s chosen comparison of F and G, but is not indispensable to irrationality itself. The companion [Beukers exposition](beukers-explained.md) proves a direct Mellin identity F_n=I_n, then transforms I_n into a positive triple integral. Its completed-square bound gives 0<I_n<=2ζ(3)α^n. Thus the combined route is simply: simple partial fractions → arithmetic; Mellin transform + integration by parts → positivity and decay; elementary lcm growth → contradiction. It bypasses both recurrence certificates, Ball’s series, and initial-value matching. The recurrence remains valuable for explaining Apéry’s coefficients and for computation.

The proof is complete after the explicit algebraic identities are verified. Its deepest remaining explanatory issue is discovery: why Ball's expression and the arithmetic expression should belong to the same recurrence. Zudilin points to a specialization of Bailey's very-well-poised hypergeometric identity as a structural explanation. That route may shorten the recurrence verification but introduces gamma functions and a contour integral. A serious simplification project should compare total conceptual and verification cost, not just the number of displayed lines.

A concrete elementary improvement already obtained here is (B'): the tail sum has an n-independent prefactor. The arithmetic cancellation is also much clearer when expressed through differences between pole positions. Neither change solves the much harder problem of adapting the construction to isolate ζ(5).

## References and verification scope

- W. Zudilin, *An elementary proof of Apéry's theorem*, arXiv:math/0202159v1, 2002, https://arxiv.org/pdf/math/0202159. The principal source for all seven lemmas.
- R. Apéry, *Irrationalité de ζ2 et ζ3*, Astérisque 61 (1979), 11–13, https://www.numdam.org/item/AST_1979__61__11_0/. Retrieved separately in this repository.
- Yu. V. Nesterenko, *A few remarks on ζ(3)*, Mathematical Notes 59 (1996), no. 6, 865–880 in the Russian original pagination: cited here through Zudilin; not independently retrieved for this file.
- F. Beukers, *A note on the irrationality of ζ(2) and ζ(3)*, Bulletin of the London Mathematical Society 11 (1979), 268–272: bibliographic identification from Zudilin; integral proof treated separately.

Exact symbolic checks performed in this session: T1, T2 after division by the central summand, and the algebraic maximum-value identity. The floor-weight table was independently enumerated. The earlier rational-number checks of u_n,v_n are useful consistency checks, but do not replace the general identities above. These exact checks are now persisted in `codex/computations/check_zudilin.py`; run `python codex/computations/check_zudilin.py` (requires SymPy). The script passed with SymPy 1.14.0.

Independent audit: see [independent-review.md](independent-review.md). The reviewing agent checked the lcm argument, initial values, short certificate, factorial bounds, and simplified Ball tail estimate and found no error in those reviewed portions. Its additional exact rational logarithm bounds are recorded in `codex/computations/check_independent_review.py`. This is independent mathematical review, not a claim of formal verification.
