# A uniform theorem for truncated higher derivatives

Status: independently derived in this session; all-orders proof below with exact computational checks. A separate Zudilin/exposition agent audited Sections 1–7: normalization, common boundary factor, fixed-order asymptotics, recurrence comparison, denominator exponent and prime-adic congruence all passed. Historical novelty is unresolved. This extends the m=2 construction in `truncated-derivative.md`; it is not an irrationality proof for higher odd zeta values.

## 1. Definitions and result

Let

    P_n(t)=∏_{j=1}^n(t-j)/∏_{j=0}^n(t+j),   R_n=P_n²,
    q=2m-1,
    L_(n,m)=-(1/(2m)!) Σ_{t=1}∞ R_n^(q)(t),
    C_(n,m)=-(1/(2m)!) Σ_{j=1}^n R_n^(q)(j),
    E_(n,m)=L_(n,m)-C_(n,m)=-(1/(2m)!)Σ_{t=n+1}∞R_n^(q)(t).

Write R_n=Σ_k[a_k/(t+k)²+b_k/(t+k)], with a_k=binom(n,k)²binom(n+k,k)², Σb_k=0, a_k∈Z and D_n b_k∈Z. Then

    E_(n,m)=A_n ζ(2m+1)-(B_(n,m)+C_(n,m)),
    A_n=Σ_k a_k,
    B_(n,m)=Σ_k[a_k H_k^(2m+1)+(b_k/(2m))H_k^(2m)].

**Theorem.** For each fixed integer m≥2, as n→∞,

    E_(n,m) ~ π (2log n)^(2m-3)
                   /[(2m)(2m-3)! n³] ·16^-n.           (1)

In particular, each fixed-m tail is eventually positive and has exact exponential rate 1/16. Moreover,

    (2m)D_n²D_(2n)^(2m-1)(B_(n,m)+C_(n,m))∈Z.         (2)

For every odd prime p>2m with p≠17 and n=(p+1)/2,

    v_p(C_(n,m))=v_p(B_(n,m)+C_(n,m))=-(2m-1).        (3)

Thus truncation repairs exponential decay for all fixed odd zeta values, but introduces genuinely new primes into the rational coefficient. Increasing derivative order does not improve the decay base, while it increases the available denominator upper bound.

## 2. Put all moving-boundary terms into one rational factor

Let T be the Apéry recurrence operator,

    (Tu)_n=(n+1)³u_(n+1)-(2n+1)(17n²+17n+5)u_n+n³u_(n-1),

and s_n(t)=4(2n+1)[-2t²+t+(2n+1)²]. Differentiating the short certificate q times and keeping the moving lower boundary gives

    (T E_(_,m))_n=g_(n,m),
    (q+1)!g_(n,m)=
       [s_n(t)R_n(t)]^(q)|_(t=n+1)
       +(n+1)³R_(n+1)^(q)(n+1)-n³R_(n-1)^(q)(n).      (4)

The signs are the same as in the independently audited q=3 case; the differentiated order does not change the negative normalization of E.

Use a local variable x around t=n+1. The shift ratios give exact identities

    R_(n+1)(n+1+x)/R_n(n+1+x)=x²/(2n+2+x)²,
    R_(n-1)(n+x)/R_n(n+1+x)
       =[(2n+x)(2n+1+x)/(n+x)²]².

Therefore, writing [x^q] for Taylor coefficient extraction,

    g_(n,m)=1/(q+1) [x^q]{R_n(n+1+x) W_n(x)},          (5)

where the same rational function works for every derivative order:

    W_n(x)=s_n(n+1+x)
       +(n+1)³ x²/(2n+2+x)²
       -n³[(2n+x)(2n+1+x)/(n+x)²]².                  (6)

The simplification is structural: all moving-boundary contributions are encoded in W_n, and derivative order merely selects a Taylor coefficient.

## 3. The cancellation has exactly order two

Direct rational differentiation at x=0 gives

    [x^0]W_n=[x^1]W_n=0,
    [x²]W_n=-(287n²+335n+100)/(4n),
    [x³]W_n=(351n²+504n+176)/(4n²).                   (7)

For each fixed j≥2,

    [x^j]W_n=O_j(n^(3-j)).                            (8)

To verify the general bound, write x=ny. The third term of (6) is n³ times a rational function of y and 1/n analytic near y=0 with bounded derivatives of any fixed order. The second term is

    n³(1+1/n)³ y²/(2+2/n+y)²,

with the same property. The first term is a quadratic polynomial in x with coefficients of the corresponding orders n³,n²,n. Coefficient extraction in x therefore yields (8). The two vanishing coefficients are exact, not just asymptotic.

This explains the special m=1 case: q=1, so (5) is identically zero. Its truncated and full first-derivative sequences coincide because derivatives vanish at double zeros, and its recurrence is homogeneous. For every m≥2, q≥3, the nonzero quadratic coefficient (7) is detected.

## 4. Local logarithmic derivatives determine the forcing

Put

    r_n=R_n(n+1)=[n!²/(2n+1)!]².

For fixed k≥1, let ell_k be the kth derivative of log R_n at n+1. Then

    ell_1=2(2H_n-H_(2n+1))=2log n+O(1),
    ell_k=O_k(1) for each fixed k≥2.

The second bound follows from convergence of reciprocal-k-power sums. Thus

    R_n(n+1+x)=r_n exp(ell_1x+ell_2x²/2!+...).

For each fixed a≥0 its normalized coefficient has leading term

    [x^a](R_n(n+1+x)/r_n)
       ~ (2log n)^a/a!,                              (9)

with the interpretation 1 when a=0. This follows by expanding the exponential: the only term of logarithmic degree a is ell_1^a/a!, and all terms using ell_k with k≥2 have lower logarithmic degree. For a=1 the assertion is immediate.

In the convolution (5), the quadratic coefficient of W_n dominates. Indeed its size is asymptotic to -(287/4)n; the j≥3 terms are O(n^(3-j)(log n)^(q-j)), negligible compared with n(log n)^(q-2). Therefore

    g_(n,m) ~ -(287/4) n r_n
                   (2log n)^(q-2)/[(q+1)(q-2)!].      (10)

The standard central-binomial asymptotic gives

    r_n ~ π/(4n)·16^-n.

Substitution gives the exact dominant forcing

    g_(n,m) ~ -(287π/16)
                   (2log n)^(q-2)/[(q+1)(q-2)!]·16^-n. (11)

The π constant uses Stirling's formula or Wallis's equivalent leading estimate. The existence of an exponential bound requires only the elementary bound r_n≤16^-n and polynomial estimates on the fixed-order jets.

## 5. Control the recurrence solutions, rather than only matching coefficients

First, for fixed q, the tail E_(n,m) grows at most polynomially. For integer t>n,

    |P_n(t)|≤1/t,

and each logarithmic derivative through order q is O_q(n+1), since all denominator distances are nonzero integers. Repeated differentiation gives |R_n^(q)(t)|≤C_q(n+1)^q/t², hence

    E_(n,m)=O_q((n+1)^(q-1)).

Equation (5) and the local endpoint estimates give g_(n,m)=O_q(n^K16^-n) for some fixed K. The maximum-principle argument in `truncated-derivative.md` then applies unchanged: for every θ∈(1/16,1),

    E_(n,m)=O_(m,θ)(θ^n).                             (12)

For completeness, θ^n and 2^n have negative T-value eventually because their bases lie between the roots α=17-12sqrt2 and α^(-1) of z+z^(-1)=34. A sufficiently large M makes T(Mθ^n)≤-|g| and handles a fixed left boundary. The barrier Mθ^n+ε2^n dominates the polynomially bounded E at infinity; an interior positive maximum contradicts the negative diagonal dominance of T. Repeat for -E and send ε↓0. This proves (12), in particular E→0.

Now define the positive candidate b_n to be the right-hand side of (1), for n≥2. Since its neighboring ratios tend to 1/16 and 16,

    (Tb)_n/(n³b_n)→1/16+16-34=-287/16.

Equations (11) and the definition of b imply (Tb)_n~g_(n,m)<0. The residual r_n^*=g_(n,m)-(Tb)_n is o(-(Tb)_n).

Choose α<θ_0<1/16. For arbitrary ε>0 and all sufficiently large n,

    |r_n^*|≤ε[-(Tb)_n],   T(θ_0^n)<0.

The barrier εb_n+Mθ_0^n can be chosen to dominate |E-b| at the left boundary, and its T-value is at most -|r^*|. Both E-b and the barrier tend to zero by (12). Any positive value of their difference would imply an attained positive interior maximum, again impossible. Applying the argument to both signs gives

    |E_(n,m)-b_n|≤εb_n+Mθ_0^n.

Divide by b_n and let n→∞; θ_0^n/b_n→0. Since ε was arbitrary, E_(n,m)/b_n→1. This proves (1), not merely a formal asymptotic match.

## 6. Denominator control by Taylor coefficients

At an initial zero j, write

    R_n(j+x)=x² U_j(j)² V_j(x),
    U_j(j)=P_n'(j),
    V_j(x)=∏_{1≤l≤n,l≠j}(1+x/(j-l))²
                 /∏_{l=0}^n(1+x/(j+l))².

Every nonzero integer denominator in V_j has absolute value at most 2n. The coefficient [x^a]V_j therefore has denominator dividing D_(2n)^a: expand each denominator using (1+z)^(-2)=Σ_{k≥0}(-1)^k(k+1)z^k and collect total degree a.

The independently audited factorial valuation bound gives

    D_nD_(2n)U_j(j)∈Z.

Since

    C_(n,m)=-(1/(q+1))Σ_{j=1}^n U_j(j)²[x^(q-2)]V_j,

we conclude

    (q+1)D_n²D_(2n)^q C_(n,m)∈Z.

The old coefficient satisfies (q+1)D_n^(q+2)B_(n,m)∈Z, and D_n^(q+2) divides D_n²D_(2n)^q. This proves (2).

## 7. The same exceptional prime 17 persists at every fixed order

Take an odd prime p>q+1=2m, p≠17, and n=(p+1)/2. Only the endpoints j=n-1,n have a factor p among j,j+1,...,j+n. All other terms in the correction are p-integral.

At one of these two endpoints, factor

    V_j(x)=(1+x/p)^(-2) A_j(x),

where every Taylor coefficient of A_j is p-integral and A_j(0)=1. Put a=q-2, an odd integer. Then

    p^a[x^a]V_j(x) ≡ (-1)^a(a+1)=-(q-1) mod p.

All terms drawing positive degree from A_j are divisible by p after multiplication by p^a, so no omitted cross term affects this congruence. Since U_j(j) has valuation -1,

    p^q C_(n,m)
       ≡ (q-1)/(q+1) [u²+v²] mod p,
    u=pP_n'(n-1),  v=pP_n'(n).

As in the q=3 case, v/u≡-1/4 mod p. Thus

    p^q C_(n,m) ≡ [(q-1)/(q+1)]·17u²/16 mod p.        (13)

This is nonzero because p>q+1, p is odd and p≠17. It proves v_p(C)=-q. The old B is p-integral since p>n and p>2m, so v_p(B+C)=-q as well. This proves (3).

The persistent 17 comes from the squares of two neighboring residue ratios, 1+1/16; the derivative order only contributes (q-1)/(q+1). Thus an apparently isolated exception in the ζ(5) calculation reflects a stable local structure.

## 8. What this uniform theorem tells us

The full untruncated sequence has the asymptotic

    L_(n,m)~(4log n)^(2m-3)/[(2m)(2m-3)!n^4],

as derived independently in `log-kernel-extension.md`. Deleting the initial zeros removes this polynomial term and leaves the positive tail (1). The analytic repair therefore works uniformly.

However, even the prime number theorem applied to the denominator upper bound (2) gives logarithmic growth approximately

    2n+(2m-1)·2n=4mn,

whereas the decay exponent is always log16≈2.773. These available bounds do not produce an irrationality contradiction. The prime-family theorem establishes genuinely new denominator support, but does not prove the full denominator grows at that upper-bound rate. Sharper arithmetic cancellation or a modified construction remains possible and requires separate proof.

The lesson is concrete: increasing derivative order by itself has no analytic-rate benefit after truncation, while its arithmetic burden increases. A successful improvement would need a new arithmetic identity or a different family with a better balance.

## Verification

Run `python codex/computations/check_higher_truncated_tail.py`. It verifies the exact coefficients (7), the mixed denominator bound and prime-family valuations for finite samples, and numerical asymptotic ratios with exact rational coefficients. Finite tests supplement the general proofs above. The all-m recurrence comparison and arithmetic were independently audited as recorded above. Section 8 also cites the separately proved full-sequence asymptotic.
