# Shifting the odd zeta index: an exact operator and a genuine arithmetic cost

Status: independently derived identities and obstruction theorems. The operator and universal obstruction were also derived independently by the parent researcher. These results do not prove that Fauzan's determinant arithmetic cannot be extended; they rule out a particular cost-free transfer through a universal simple-pole-preserving differential operator.

## 1. Two integrations by parts give the exact shift

For odd s>=3 let

    w_s(y)=2y^s f^(s-1)(y)/(s-1)!,
    f(y)=1/(exp(2πy)-1),
    μ_s(R)=integral_0^infinity R(y²)w_s(y)dy.

Here μ_s denotes the actual positive integral functional, extended to rational functions with poles strictly off [0,infinity), and at most polynomial growth. On polynomials and simple poles at negative integer squares it is the single-zeta functional of `odd-moment-generalization.md`.

Put θ=t d/dt. Then

    μ_(s+2)(R)=μ_s(L_s R),
    L_s=[(2θ+s+1)(2θ+s+2)]/[s(s+1)]
       =[4t² d²/dt²+(4s+10)t d/dt+(s+1)(s+2)]/[s(s+1)].  (1)

Proof: write the left side as 2/(s+1)! times the integral of y^(s+2)R(y²)f^(s+1)(y). Integrate by parts twice. At zero, R(y²) is regular and both boundary products are O(y): f^s=O(y^(-s-1)), f^(s-1)=O(y^-s), while y^(s+2)R and its first derivative have orders s+2 and s+1. At infinity exponential decay beats every polynomial factor. Direct differentiation gives

    d²/dy²[y^(s+2)R(y²)]
       =y^s[(s+1)(s+2)R+(4s+10)tR'+4t²R''], t=y².

Comparing the normalizing factorials proves (1).

On a monomial t^e the operator multiplies by

    (2e+s+1)(2e+s+2)/[s(s+1)].

This recovers the exact polynomial-moment ratio already derived from Bernoulli numbers. No new large-prime denominators enter that polynomial ratio, apart from primes dividing the fixed integer s(s+1).

## 2. The same operator necessarily creates repeated poles

For a>0,

    L_s[1/(t+a)]
      =[(s-1)/(s+1)]/(t+a)
       +[2(2s-3)a/(s(s+1))]/(t+a)²
       +[8a²/(s(s+1))]/(t+a)³.                       (2)

The cubic-pole coefficient is nonzero. This is not merely an inconvenient presentation: the repeated-pole integral contains higher zeta values, as the triangular theorem in `compactification-obstruction.md` proves.

For example, s=5 and a=1 gives

    μ_7(1/(t+1))=(2/3)J_1+(7/15)J_2+(4/15)J_3,

where J_k=μ_5((t+1)^(-k)). The exact identities are

    J_1=ζ(5)-3/4,
    J_2=-2ζ(5)+(5/2)ζ(6)-1/4,
    J_3=ζ(5)-(35/8)ζ(6)+(15/4)ζ(7)-3/16.

Their indicated combination is ζ(7)-2/3, exactly the s=7 simple-pole formula. The ζ(5) and ζ(6) coefficients cancel, but evaluating the combination required ζ(7) in the extended functional. It cannot be obtained from the old affine-in-ζ(5) functional by pretending that repeated poles retain the old domain.

Thus (1) is a valid analytic transport identity, not a free reuse of all simple-pole prime-adic bounds.

## 3. No universal rational differential operator can avoid this pole growth

Let D be a finite-order differential operator with rational-function coefficients in Q(t). Suppose D maps the space of rational functions with arbitrary polynomial part and only simple poles among {-1,-4,-9,...} into itself.

**Pole-preservation lemma.** D must be multiplication by a polynomial in Q[t].

Proof: let p_q(t)d^q/dt^q be its nonzero highest-order term, q>=1. Exclude the finitely many j where a coefficient has a pole. Applied to 1/(t+j²), the coefficient of pole order q+1 at -j² is (-1)^q q!p_q(-j²). All lower-order derivatives contribute smaller pole orders. Preservation of simple poles therefore forces p_q(-j²)=0 for infinitely many j, impossible for a nonzero rational function. Descending in order leaves only multiplication by p_0(t). Applying this to1 shows p_0 has only allowed simple poles. If one of those poles were present, multiplication by the simple pole at the same point would make a double pole. Thus p_0 has no finite poles and is a polynomial. QED.

**No-shift theorem.** No such operator D can satisfy

    μ_s(DR)=μ_(s+2)(R)

for every polynomial R, and hence none can implement the index shift on the full simple-pole domain.

Proof: by the lemma D is multiplication by p(t)=sum_(j=0)^d c_j t^j. The positive moments are

    m_e=μ_s(t^e)=2(2e+s)! ζ(2e+2)/[(s-1)!(2π)^(2e+2)].

For fixed j, factorial ratios and ζ(2e+2)->1 give

    m_(e+j)/m_e ~ e^(2j)/π^(2j).

But the desired shift ratio is asymptotic to 4e²/[s(s+1)]. In the sum sum_j c_j m_(e+j)/m_e, the highest nonzero degree gives a unique dominant power e^(2d). Hence d must equal1. Comparing leading coefficients then forces

    c_1=4π²/[s(s+1)],

which is not rational. Contradiction. The zero polynomial cannot implement the nonzero moments. QED.

This theorem has precise limits. It requires one finite-order rational-coefficient differential operator preserving all prescribed simple poles. It does not rule out a K-dependent operator vanishing at finitely many relevant poles, a nonlocal transform, a new coefficient basis, or cancellations specific to a determinant.

## 4. A real denominator increase already occurs in the pole values

Write the constant coefficient of the simple-pole functional as

    h_s(j)= -j^(s-1)H_j^(s)-1/(s-1)+1/(2j).

Let p>s+1 be prime and choose j=p+1. Since j<2p, exactly one denominator in H_j^(s) is divisible by p, namely the term1/p^s. The factor j^(s-1) is a p-unit, and every other term, including both correction constants, is p-integral. Therefore

    v_p(h_s(p+1))=-s,
    v_p(h_(s+2)(p+1))=-(s+2).                         (3)

This proves an actual increase of two prime powers, not merely a worsening sufficient upper bound.

To formulate a uniform cost, let D_s^pole(K) be the least common denominator of h_s(1),...,h_s(K). For p in (K/2,K-1] and p>s+1,

    v_p(D_s^pole(K))=s,
    v_p(D_(s+2)^pole(K))=s+2.                         (4)

Indeed j<p is p-integral, j=p has valuation exactly -1 (the -1/p and +1/(2p) terms combine to -1/(2p)), and p<j<=K<2p has valuation exactly -s. The index p+1 exists because p<=K-1.

Consequently any positive integer M for which D_s^pole(K)M clears all the new pole constants must be divisible by

    product_(K/2<p<=K-1, p>s+1) p².

Under the prime number theorem this implies

    log M >= K+o(K), K->infinity with s fixed.         (5)

The statement concerns an **additional multiplier relative to the previous least common pole denominator**. It does not claim that the ratio of the two least common denominators has no cancellations at other primes.

If one transfers this scalar multiplier entrywise to an h-by-h matrix with h~βK, the resulting sufficient determinant normalization pays at least exp((β+o(1))K²). This shows why an O(K log K) cost cannot be assumed for the naive uniform-functional transfer.

However (5) is not a lower bound on the primitive denominator of Fauzan's determinant. Residue factors, special bases, and cancellations among entries may improve determinant valuations dramatically. Fauzan's local arithmetic arguments are designed to exploit precisely such structure. A new proof for s=7 would need to show that those coupled cancellations compensate for the genuine pole-value loss, rather than infer transfer from the easy polynomial-moment ratio.

## 5. Consequence for reuse of the ζ(5) normalizer

The analytic equilibrium budget can be shared by fixed odd s after lower-order weight changes, but the arithmetic transfer is not automatic:

* Polynomial moments shift by an explicit quadratic factor with only fixed rational denominators.
* Simple-pole values acquire demonstrably higher p-adic denominators for a whole band of large primes.
* The exact differential identity handles the shift by creating repeated poles, which lies outside the original simple-pole arithmetic theorem.
* No universal rational differential operator preserving the entire simple-pole domain can remove this issue.

These are concrete barriers to a particular reuse argument, not a claim that all odd zeta values resist determinant methods. The next defensible target is a K-dependent matrix-level identity or valuation comparison, with its actual prime costs tracked before any general irrationality conclusion.
