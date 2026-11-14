# Designing irrationality proofs through three independent requirements

Status: general sufficient theorem proved below, together with a complete classical Legendre/logarithm family as a concrete example. The parameter search is finite experimental evidence, not a classification. No historical novelty claimed.

The useful discovery question is: **can one make weak identification, arithmetic, and recurrence structure separately easy, instead of forcing one formula to reveal everything?** The Apéry simplification suggests precisely this division.

## 1. Three requirements with different jobs

1. **Identify the constant with a weak analytic estimate.** Construct rational coefficients A_n,B_n for which F_n=A_n xi-B_n is bounded, or grows more slowly than A_n. Then B_n/A_n tends to xi. Exponential smallness of F_n is not needed at this stage.
2. **Protect arithmetic.** Prove that A_n is a positive integer and that a controlled integer d_n clears B_n. The size of d_n, not just the existence of rational coefficients, matters. Factorial denominator growth is usually too expensive for exponential approximation.
3. **Use a shared recurrence and a nonzero determinant.** If A and B satisfy the same second-order recurrence, their Casoratian determines all increments of B_n/A_n. Its nonzero sign supplies nonvanishing, and growth of A supplies smallness.

These requirements are independent. A beautifully decaying integral may have unusable denominators. An integral sequence with a recurrence may approximate an unidentified constant. A bounded linear form may be identically zero if the two solutions are proportional.

## 2. A quantitative theorem for power-weighted recurrences

Let d>=1 be an integer and let b(n) be a real polynomial of degree d with positive leading coefficient beta>2. Suppose rational sequences A_n,B_n satisfy

    (n+1)^d u_(n+1)=b(n)u_n-n^d u_(n-1), n>=1,          (R_d)

and assume A_n is a positive integer for every n>=0. Set

    lambda=(beta+sqrt(beta²-4))/2>1,
    w=A_0B_1-A_1B_0.

Assume w is nonzero. Then:

**Growth.** For every r with 1<r<lambda, A_(n+1)>=r A_n for all sufficiently large n. In particular A_n>=C_r r^n for a positive constant C_r.

**Approximation.** The ratios B_n/A_n converge to some real L, and

    L-B_n/A_n=w sum_{j=n}^infinity 1/[(j+1)^d A_j A_(j+1)].  (1)

Thus A_nL-B_n has the fixed nonzero sign of w. Once the ratio bound r holds,

    |A_nL-B_n|
      <= |w| r/[(r²-1)(n+1)^d A_n].                    (2)

**Identification.** If xi is a known real number and

    A_n xi-B_n=O(mu^n n^M),  with mu<lambda,

for some fixed M and mu>=0, then L=xi. A bounded representation is the particularly convenient case mu=1,M=0.

**Irrationality criterion.** Suppose positive integers d_n satisfy d_nB_n integral, and

    limsup d_n^(1/n)<=Delta<lambda.

Then L is irrational. Combined with the identification condition, this proves xi irrational.

The strict exponential threshold is

    lambda>Delta,

or equivalently beta>Delta+1/Delta when Delta>=1. This theorem is sufficient; equality of the exponential rates may require finer polynomial estimates and is not addressed here.

### Proof of growth without a black-box asymptotic theorem

Write

    r_n=A_n/A_(n-1),
    p_n=b(n)/(n+1)^d -> beta,
    q_n=(n/(n+1))^d ->1.

Then r_(n+1)=p_n-q_n/r_n. Since p_n-q_n tends to beta-1>1, once n is large enough, a ratio r_n>=1 implies r_(n+1)>1.

There must be such a ratio after any prescribed index: otherwise the positive integer sequence A_n would strictly decrease forever. Thus eventually all ratios are at least one.

Fix 1<r<lambda and put delta=beta-r-1/r>0. The function beta-x-1/x is decreasing on x>=1, so it is at least delta on [1,r]. As p_n->beta and q_n<=1, eventually

    r_(n+1)>=beta-1/r_n-delta/2.

While 1<=r_n<=r this increases r_n by at least delta/2. It reaches r in finitely many steps. Once r_n>=r, monotonicity of beta-1/x makes the same lower bound exceed r, so the bound persists. This proves growth.

For an explicit coefficient bound, write b(n) in powers of m=n+1 as

    b(n)=beta m^d+c_(d-1)m^(d-1)+...+c_0.

Then p_n>=beta-K/(n+1), with K=sum_j max(0,-c_j). The preceding argument works whenever n+1>=2K/delta. If a ratio >=1 is not already known there, it must occur within at most A_N subsequent indices, because strict integer decrease can last only finitely long. This cutoff is crude but effective.

### Proof of approximation and arithmetic criterion

Direct substitution in (R_d) gives

    W_n=A_nB_(n+1)-A_(n+1)B_n
        =(n/(n+1))^d W_(n-1)=w/(n+1)^d.

Divide by A_nA_(n+1), sum, and use eventual geometric growth. This proves convergence and (1). For j=n+k, bound (j+1)^d>=(n+1)^d and A_jA_(j+1)>=r^(2k+1)A_n². The geometric sum is

    sum_{k>=0}r^(-2k-1)=r/(r²-1),

giving (2). Identification follows by choosing r strictly between max(1,mu) and lambda, so F_n/A_n->0.

For irrationality choose r between max(1,Delta) and lambda. Then d_n|A_nL-B_n|->0. If L=p/q with integers p,q and q>0, the nonzero number

    q d_n(A_nL-B_n)=d_nA_np-q d_nB_n

is an integer of fixed sign. Its absolute value cannot tend to zero. Contradiction.

This proof shows why a nonzero Casoratian is an indispensable nonvanishing condition. If w=0, the two solutions are proportional and (1) is identically zero.

## 3. Concrete family: log((a+1)/a) for every positive integer a

This is a classical non-Apéry example, developed here to show the entire mechanism. It uses shifted Legendre polynomials, ordinary integrals, and a d=1 recurrence. It does not depend on a transcendence theorem.

For n>=0 define

    Q_n(x)=(1/n!)(d/dx)^n [x^n(1-x)^n]
          =sum_{k=0}^n (-1)^k binom(n,k)binom(n+k,k)x^k.

Fix an integer a>=1 and put

    xi_a=log((a+1)/a),
    F_n=integral_0^1 Q_n(x)/(a+x) dx,
    A_n=Q_n(-a)=sum_{k=0}^n binom(n,k)binom(n+k,k)a^k.

A_n is visibly a positive integer.

### Arithmetic from polynomial division

Since Q_n has integer coefficients and x+a is monic with integer coefficients,

    Q_n(x)=A_n+(x+a)V_n(x),   V_n in Z[x], deg V_n<=n-1.

Integration gives

    F_n=A_n xi_a-B_n,
    B_n=-integral_0^1 V_n(x)dx,
    D_n B_n in Z,   D_n=lcm(1,...,n).

For n=0 the quotient is zero and B_0=0. For n=1,

    A_1=2a+1,   F_1=(2a+1)xi_a-2,   B_1=2.

Thus w=A_0B_1-A_1B_0=2.

### Weak identification from integration by parts

Integrate by parts n times. Every boundary term vanishes because x^n(1-x)^n has zeros of order n at both endpoints. The two signs, from integration by parts and from differentiating 1/(a+x), cancel. Therefore

    F_n=integral_0^1 x^n(1-x)^n/(a+x)^(n+1) dx.

For 0<=x<=1 and a>=1, 0<=x(1-x)/(a+x)<=1. Consequently

    0<F_n<=integral_0^1 dx/(a+x)=xi_a.

This establishes boundedness and identifies the limit once A_n grows. A sharp integral maximum is optional; the recurrence will provide the needed quantitative approximation.

### Derive the recurrence rather than quote it

The coefficient formula gives the formal generating function

    sum_{n>=0}Q_n(x)z^n=[1-2(1-2x)z+z²]^(-1/2).        (3)

One direct verification expands the right side as

    (1-z)^(-1)[1+4xz/(1-z)²]^(-1/2)
      =sum_{k>=0}(-1)^k binom(2k,k)x^k z^k(1-z)^(-2k-1).

The coefficient of z^n is correct because

    binom(2k,k)binom(n+k,2k)=binom(n,k)binom(n+k,k).

Differentiate (3) with respect to z and compare coefficients to obtain

    (n+1)Q_(n+1)(x)=(2n+1)(1-2x)Q_n(x)-nQ_(n-1)(x).

At x=-a, this proves the recurrence for A:

    (n+1)A_(n+1)=(2a+1)(2n+1)A_n-nA_(n-1).            (4)

For F, divide the polynomial recurrence by a+x and integrate. Since

    (1-2x)/(a+x)=(2a+1)/(a+x)-2,

the only possible extra term is -2(2n+1)integral Q_n. This integral is zero for n>=1, directly from the Rodrigues formula and endpoint zeros. Thus F obeys (4), and B=A xi_a-F does too.

The leading coefficient beta=4a+2 gives

    lambda_a=2a+1+2sqrt(a(a+1))=(sqrt(a)+sqrt(a+1))².

At a=1 this is 3+2sqrt2, revealing the same quadratic field that appears in Apéry growth, now with the square rather than fourth power.

### A crude growth bound already suffices

Write t=2a+1>=3. Induction in (4) shows r_n=A_n/A_(n-1)>=t: r_1=t, and

    [(2n+1)t-n/t]/(n+1)>=t.

In particular r_n>=3, and then for n>=1,

    r_(n+1)>=[3(2n+1)-n/3]/(n+1)
             =(17n/3+3)/(n+1)>=4.

Therefore A_n>=3*4^(n-1) for n>=1. With w=2 and d=1, the tail formula gives

    0<A_n xi_a-B_n<=8/[15(n+1)A_n]  (n>=1).

The elementary period-6 factorial argument in `apery-asymptotics.md` proves

    log D_n<=c6 n+O(log²n),
    c6=log2/2+3log3/4<log4.

The strict inequality needs only exp(4c6)=2²3³=108<256=4^4. Hence D_n(A_nxi_a-B_n)->0 through positive quantities, proving **log((a+1)/a) is irrational for every integer a>=1**.

This proof supplies every gate explicitly: bounded rational-coefficient linear forms, integer A_n and lcm denominators, and a nonzero Casoratian.

## 4. Unexpected links that do actual work

- **Orthogonal polynomials and Diophantine approximation.** The integral of Q_n against constants vanishes, which removes the unwanted term when integrating its recurrence. The same endpoint zeros used for orthogonality also make the integral representation positive.
- **Polynomial division and denominator economy.** Division by x+a leaves an integer polynomial quotient. Integrating its monomials produces only denominators 1,...,n, explaining why an lcm suffices instead of n!.
- **Reciprocal characteristic roots and approximation.** The recurrence has asymptotic roots lambda_a and lambda_a^(-1). The Casoratian makes their complementary roles rigorous: the integer sequence grows and the approximation error shrinks.
- **The d=1 versus d=3 analogy.** In the logarithm example, integration of simple poles isolates log; in Zudilin, differentiation of double poles isolates zeta(3). Both constructions use endpoint cancellation and integer polynomial data. The power d controls the Casoratian factor (n+1)^(-d), while the denominator exponent has to be proved independently.

This is a structural analogy, not an assertion that increasing d automatically produces zeta(d).

## 5. Finite parameter searches: arithmetic is restrictive

To keep the search comparable, use the same negative backward term and normalized initial values U_0=1,U_1=b.

For d=2 search

    (n+1)² U_(n+1)=[a n(n+1)+b]U_n-n²U_(n-1),
    3<=a<=40, 1<=b<=20.

For d=3 search

    (n+1)³ U_(n+1)=(2n+1)[a n(n+1)+b]U_n-n³U_(n-1),
    2<=a<=40, 1<=b<=20.

Exact rational computation through U_13 found:

    d=2: no parameter pair retained integrality throughout the tested range;
    d=3: only (a,b)=(17,5), the Apéry pair, survived.

These statements describe only this finite box, this normalization, and this sign. They do not classify integral recurrences and do not rule out other families. In particular the familiar zeta(2) recurrence has a **positive** backward term, so it is outside this search and outside the fixed-sign theorem above. Its alternating Casoratian requires a corresponding modification of the positivity argument.

By contrast the d=1 Legendre construction proves an infinite integral family, not merely a finite pattern. This provides a useful baseline for future searches: a candidate should come with a structural source of integrality and a representation identifying its limit, rather than only many apparently integral terms.

## Verification and references

Run `python codex/computations/check_recurrence_design.py` for exact polynomial identities, logarithm-family recurrence/arithmetic checks, and the finite d=2,d=3 scan.

The Legendre/Rodrigues and generating-function identities are derived above. The organization is motivated by Zudilin, *An elementary proof of Apéry's theorem*, arXiv:math/0202159, and by the independently verified Casoratian theorem in `recurrence.md`. The logarithm family and its orthogonal-polynomial connection are classical; no historical originality is asserted. A focused primary-reference search for this particular presentation remains a bibliographic task, not a mathematical gap in the displayed proof.
