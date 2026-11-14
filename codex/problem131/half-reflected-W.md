# A positive half-integer companion to the cancellation sequence W

**Status.** Exact finite-field theorem, explicit positive rational formula, and reflected recurrence proved below. This gives a second family of finite-exception criteria for cubic denominators near the upper end p=2N of the prime band. It does not prove a density statement across that band or a global denominator-growth obstruction. No historical-novelty claim is made.

The starting objects and already proved identities are in `prime-band-arithmetic.md` and `W-recurrence.md`. Write

    P_n(t)=product_(r=1)^n(t-r)/product_(r=0)^n(t+r), R_n=P_n^2.

The rational constant W_n is the harmonic-sum expression in the prime-band note. Its denominator divides lcm(1,...,n)^2, and W_0=0,W_1=-13.

## 1. A direct positive formula

For an integer s, define the rational rising factorial

    (1/2)_s=product_(r=0)^(s-1)(1/2+r),                    s>=0,
    (1/2)_(-a)=1/product_(r=1)^a(1/2-r),                 a>=1.

The empty product is1. In particular,

    (1/2)_s=(2s)!/(4^s s!) for s>=0,
    (1/2)_(-a)=(-4)^a a!/(2a)!.

Set V_-1=V_0=0, and for j>=1 define

    V_j=sum_(t=0)^(2j-1)
          [((2j-1-t)!/t!)^2 ((1/2)_(t-j))^4].           (1)

Every summand is a positive rational number. Thus **V_j>0 for every j>=1**, without a sign inference from numerical data or modular residues. For j>=2 its denominator has no prime factor greater than2j-1, and V_1 is an integer.

The first two examples are especially transparent:

    V_1=16+1=17,
    V_2=1024/9+64+1/4+1/576=102545/576.                 (2)

The integer17 governing the exceptional prime in the p=2N-1 denominator family is the sum of the first two positive terms in this companion construction.

## 2. The finite-field sum extends beyond the earlier disjoint-root range

For every prime p>3 and every0<=n<p, put

    E_(n,p)=F_p minus {0,-1,...,-n}.

Then

    W_n=sum_(t in E_(n,p)) R_n(t)  (mod p).             (3)

The set in this assertion excludes **all specified denominator positions**. For large n some numerator and denominator factors cancel modulo p; it would then be inaccurate to describe E_(n,p) as the complement of the actual poles of the reduced rational function. The specified exclusion set is what the proof uses.

To verify(3), use the partial fractions

    R_n(t)=sum_(k=0)^n [a_k/(t+k)^2+b_k/(t+k)].

The a_k are integers; n! has all possible prime factors of the denominators of b_k. More precisely the explicit residue formula has only denominators l-k with0<=l,k<=n. Hence all coefficients are p-integral for n<p. The denominator positions remain distinct modulo p, whether or not some poles become removable.

For each k, the complete sums of z^-1 and z^-2 over F_p^* vanish when p>3. Subtracting the excluded positions -l, l!=k, gives

    sum_(t in E) (t+k)^-2=-H_k^(2)-H_(n-k)^(2),
    sum_(t in E) (t+k)^-1=-H_k+H_(n-k).

Multiplying by the partial-fraction coefficients proves precisely the defining harmonic expression for W_n. This proof needs n<p, not p>2n.

As ordinary integer representatives,

    E_(n,p)={1,...,p-n-1}.

Consequently, with h=(p-1)/2,

    W_n=0 (mod p) whenever h<=n<=p-1:                 (4)

every remaining evaluation point is one of the zeros1,...,n of P_n. For n=p-1 the set is empty, which causes no exception.

## 3. Wilson reflection produces the positive rational formula

Fix j>=1 and a prime p>max(3,2j+1), and put h=(p-1)/2, m=h-j. The only nonzero terms in(3) at n=m occur at

    t=m+r, 1<=r<=2j.

For these terms the factorial expression is

    P_m(m+r)=(m+r-1)!^2/[(r-1)!(2m+r)!].             (5)

All factorial arguments lie between0 and p-1. Write s=r-j-1. Since h=-1/2 modulo p,

    (h+s)!/h! = (1/2)_s (mod p),

including negative s, where this is a reciprocal product. Wilson's congruence supplies

    (h!)^2=(-1)^(h+1) (mod p),
    (2m+r)!=(p-1-(2j-r))!
             =(-1)^(r+1)/(2j-r)! (mod p).

Substituting into(5) and squaring removes both signs:

    R_m(m+r)=[(2j-r)!/(r-1)!]^2 ((1/2)_(r-j-1))^4 (mod p).

Summing the2j terms and reindexing r=t+1 proves the universal reflection congruence

    **W_((p-1)/2-j)=V_j (mod p).**                    (6)

The rational V_j is independent of p. All its denominators are prime to the primes allowed here. This is a finite factorial calculation, not a proposed analytic continuation of W_n to negative half-integer indices.

There is also a useful boundary case: if p=2j+1>3 is prime, the same argument remains valid with m=0. It gives V_j=0 modulo p because W_0=0. Thus this boundary prime, when it exists, necessarily divides the numerator of V_j. It is excluded from the main strict inequality in(6).

## 4. The reflected recurrence

The positive sequence(1) satisfies, for every j>=0,

    (2j+1)^3 V_(j+1)-4j(68j^2+3)V_j+(2j-1)^3 V_(j-1)
      =-(1148j^2+192j+17)/(2j-1)
          *binom(2j,j)^2/256^j.                     (7)

At j=0 this says V_1-V_-1=17, consistent with(2). For j>=1 the forcing on the right is negative.

Here is a short exact derivation using the proved W recurrence. For a prime p sufficiently large relative to j, substitute n=h-j=-j-1/2 modulo p into

    (n+1)^3 W_(n+1)
      -(2n+1)(17n(n+1)+5)W_n+n^3W_(n-1)
       =-(287n^2+239n+52)binom(2n,n)^2/[4(n+1)].

Use(6) for the three W values. The central binomial coefficient has the reflection

    binom(2(h-j),h-j)^2=binom(2j,j)^2/256^j (mod p).  (8)

Indeed binom(2h,h)^2=1 modulo p, and each backwards step contributes

    binom(2(h-s),h-s)/binom(2(h-s+1),h-s+1)
      =(2s-1)/(8s) (mod p).

Their product is binom(2j,j)/16^j before squaring. Substitution n=-j-1/2 and multiplication by-8 now gives(7) modulo every sufficiently large prime.

Both sides of(7) are fixed rational numbers for a fixed j. If their difference were nonzero, its numerator could be divisible by only finitely many primes. Euclid's infinitude of primes therefore proves the exact rational equality(7). This is a legitimate way to recover a rational identity from all sufficiently large finite fields; it assumes no distribution theorem about primes.

Conversely, the recurrence with V_-1=V_0=0 uniquely determines V, and the original proposed values follow:

    V_3=11537105377/3686400,
    V_4=155541579867829/2312110080.

## 5. An independent quantitative growth bound

Although positivity is immediate from(1), the recurrence also gives an elementary uniform estimate:

    V_(j+1)>10V_j for all j>=1.                      (9)

Write a_j=4j(68j^2+3), b_j=(2j+1)^3, c_j=(2j-1)^3, and let f_j be the positive magnitude of the forcing in(7). Since binom(2j,j)<=4^j,

    0<f_j<=1357j/16^j<=1357/16, j>=1.

Moreover

    a_j-c_j-10b_j=184j^3-108j^2-54j-9>=13j^3.

If V_j>=17 and0<=V_(j-1)<=V_j, the recurrence gives

    b_j(V_(j+1)-10V_j)
       >=13j^3 V_j-f_j>=221-1357/16>0.

The hypotheses hold at j=1 and propagate. Thus V_j>=17*10^(j-1), with strict inequality for j>=2. This is an independent sign proof and rules out an identically zero numerator in every fixed-offset family.

## 6. Exact denominator criteria near p=2N

Return to the third-derivative correction C_N of `prime-band-arithmetic.md`, and write B_N for its p-integral original rational coefficient when p>N. Fix j>=1 and take a prime p>max(3,2j+1). Put

    N=(p+2j-1)/2, equivalently p=2N-(2j-1).

Then N<p<=2N and p-N-1=(p-1)/2-j. The prime-band theorem and(6) give

    p^3 C_N = V_j/2 (mod p),
    v_p(B_N+C_N)=-3 iff p does not divide numerator(V_j).   (10)

If the numerator is divisible by p, the valuation is at least-2, rather than necessarily nonnegative. Since V_j>0, every fixed j has only finitely many exceptional primes. For example:

| j | p in terms of N | numerator factors of V_j | exceptions satisfying p>2j+1 |
|---|---|---|---|
|1|2N-1|17|17|
|2|2N-3|5*20509|20509|
|3|2N-5|7*1648157911|1648157911|
|4|2N-7|139*809*1163*1189333|139,809,1163,1189333|

All factors displayed in the table are prime, checked by exact trial division in the accompanying script. The factors5 and7 in rows2 and3 are precisely the boundary primes described after(6).

These are arithmetic statements about individual fixed-offset families. For one fixed N, finitely many fixed offsets account for only finitely many primes and at most O(log N) logarithmic denominator information. The positivity or growth of V_j does not prove modular nonvanishing when j varies with p and N. No weighted-density claim or global failure theorem for the truncated construction follows from(10).

## 7. Reproduction

Run `python codex/computations/check_half_reflected_W.py`. It checks the positive formula against the recurrence, the first explicit values and numerator factorizations, the extended finite-field sum for every n<p in a bounded prime range, half-band vanishing, the binomial reflection and V congruence, including boundary cases. It independently compares the W recurrence with the defining harmonic expression at small indices. The proofs above do not depend on extrapolating those finite checks.
