# Independent review of the proof for all positive integral powers of π

Reviewed: `pi-powers.md`, especially Section 3, against its displayed definitions. Reviewer: independent integral/symmetric-square research agent. This is a mathematical audit of the reconstructed argument, not verification of its historical attribution or a substitute for reading Niven's primary paper.

**Verdict:** no gap found in the displayed proof that π^r is irrational for every positive integer r. Its denominator exponent c^(mp) is sufficient; the nonzero-mod-p and complex-segment estimates are valid. The proof exposes a classical transcendence argument, not a novel elementary shortcut. The author's source-retrieval limitations should remain stated.

## A. The finite exponential relation

Under π^r=a/b, the polynomial h(X)=X^(2r)-(-1)^r(a/b)^2 does indeed have iπ as a root, since (iπ)^(2r)=(-1)^rπ^(2r). Its constant coefficient is nonzero, so no root is zero and its roots are distinct. Consequently the nonzero subset-sum list is nonempty: it contains all singleton sums, although repetitions elsewhere are allowed.

The product ∏(1+e^(α_j)) vanishes because of the single factor 1+e^(iπ). Expanding it over subsets yields k+Σe^(β_j)=0, with each subset counted once and repeated β counted with its multiplicity. In particular k is a positive integer (at least one from the empty subset). No claim that field automorphisms commute with exponentiation is needed or used; that would be invalid in a naive conjugation argument, and the proof correctly avoids it.

The full subset-sum polynomial is rational because its coefficients are symmetric polynomials with integer coefficients in the α_j. Every permutation of the α_j permutes the subsets. The fundamental theorem of symmetric polynomials then expresses its coefficients in the rational elementary symmetric functions of the α_j. Since exactly k of its roots are zero, division by X^k leaves a monic rational polynomial with nonzero constant coefficient. Clearing denominators by a positive integer c gives the displayed g∈Z[X], with degree m, leading coefficient c and g(0)≠0.

## B. Exact degree and factorial bookkeeping

Let H(X)=X^(p-1)g(X)^p and D=(m+1)p-1. Then deg H=D, so F_p=Σ_{ell=0}^D f_p^(ell) is exactly the finite antiderivative device required by F_p-F_p'=f_p. No top derivative is omitted.

At zero, H has a zero of order exactly p-1, since g(0)≠0. Thus the unique derivative below order p that can survive is order p-1, whose normalized value is c^(mp)g(0)^p. For ell≥p, ell!/(p-1)! includes a factor p. This proves that F_p(0) is an integer with the specified residue.

At each β_j, g^p has a zero of order at least p, even when β_j is a repeated root of g. Thus every derivative below order p vanishes. For ell≥p, writing derivatives using binomial coefficients gives exactly

    f_p^(ell)(β_j)/p
       =Σ_{v≥ell} H_v binom(v,ell) [ell!/(p(p-1)!)]
                    c^(mp) β_j^(v-ell).

The bracketed factor is integral. The maximal exponent v-ell is

    D-p=mp-1,

not mp+p-1; the derivative order at least p is essential. Since cβ_j is an algebraic integer,

    c^(mp)β_j^(v-ell)=c^(mp-v+ell)(cβ_j)^(v-ell)

is an algebraic integer; the remaining integer exponent is at least 1. Thus c^(mp) is sufficient, with one factor of c to spare. No missing growing denominator remains.

The proof that cβ_j is integral is also correct: c^(m-1)g(Y/c) is monic in Z[Y]. It has leading term Y^m, and the coefficient of Y^v is g_v c^(m-1-v), with nonnegative exponent for v≤m-1.

## C. From algebraic integers to ordinary integers

Each F_p(β_j)/p is an algebraic integer. Their sum is rational because F_p has rational coefficients and the power sums of the complete multiset of roots of the rational polynomial g are rational. This remains true with repeated roots. The closure of algebraic integers under addition and the fact Q∩overline(Z)=Z therefore give ΣF_p(β_j)∈pZ.

Consequently N_p=kF_p(0)+ΣF_p(β_j) is an ordinary integer and has nonzero residue k c^(mp)g(0)^p modulo any prime p not dividing kc g(0). Infinitely many primes satisfy this condition. Their infinitude is all that is needed; no estimate on prime gaps is used.

The underlying algebraic-integer facts are standard undergraduate algebra. For an olympiad audience, a short separate appendix should prove closure under sums/products (for example using a finite generated Z-module stable under multiplication), rational algebraic integers being integers (rational-root argument), and the elementary symmetric-polynomial theorem. These are prerequisites, not mathematical gaps in the stated undergraduate-level proof.

## D. Complex integration and smallness

The identity (e^(-z)F_p(z))'=-e^(-z)f_p(z) holds as an identity of entire functions. It can be integrated along z=tβ_j, 0≤t≤1, using the ordinary real-variable chain rule, so the argument does not need contour deformation or Cauchy's theorem.

The integrated identity has the correct sign:

    F_p(β_j)-e^(β_j)F_p(0)
       =-e^(β_j)∫_0^(β_j)e^(-z)f_p(z)dz.

Summing and substituting Σe^(β_j)=-k produces exactly N_p on the left. All β_j, c, m and g depend on the hypothesized algebraic relation, but remain fixed as p varies. The union of the finitely many segments is compact. Choose fixed M,G≥1 bounding |z|, |β_j| and |g(z)| as in the note. Then

    |f_p(z)|≤(c^m M G)^p/(p-1)!,

because |z|^(p-1)≤M^p. Segment lengths are at most M, and both exponential factors have modulus at most e^M. Therefore

    |N_p|≤m M e^(2M)(c^m M G)^p/(p-1)!→0.

Factorial growth dominates this fixed exponential: the ratio of consecutive C^p/(p-1)! terms is C/p→0. Together with N_p∈Z\{0}, this is a complete contradiction.

## E. Independent representative computation

For a small arithmetic stress test, use g(X)=4X²+1, roots ±i/2, m=2, c=4, k=2, and p=3,5,7. This corresponds to the root/subset-sum bookkeeping for the hypothetical relation π=1/2; it is only a test of the endpoint arithmetic, not of the exponential relation. Exact SymPy computation confirms F_p(0)∈Z, ΣF_p(β_j)∈pZ and N_p not divisible by p. The residues are respectively:

| p | F_p(0) mod p | ΣF_p(β_j) mod p | N_p mod p |
|---|---:|---:|---:|
| 3 | 1 | 0 | 2 |
| 5 | 1 | 0 | 2 |
| 7 | 2 | 0 | 4 |

The general algebraic argument above, rather than these finite checks, establishes the theorem.

## F. Interpretation and comparison

The proof does not obtain all powers from irrationality of π alone; that invalid inference is explicitly avoided. Instead it disproves every algebraic equation of the special form π^r=a/b by a construction powerful enough to handle a general algebraic equation for iπ. Hence it essentially establishes transcendence of π. Calling this an “explicit classical transcendence construction using undergraduate tools” is accurate; calling it a simpler direct irrationality proof for all powers would conceal its conceptual cost.

For ζ(3), the reviewer comparison is in `proof-comparison.md`. For the intended audience the Beukers proof has the least calculation once multiple integration is accepted, while the Mellin bridge gives the best explanation of its relation to Zudilin. The symmetric-square differential reduction is a valuable advanced connection but has not yet replaced the arithmetic and constant-identification steps of a full proof.

## G. Added Section 5: transcendence of exp(r) for nonzero rational r

The added extension is also valid. If exp(r) is algebraic and nonzero, its minimal polynomial has nonzero constant coefficient, and after clearing denominators supplies the required integer relation a_0+Σa_j exp(jr)=0. The retained jr are distinct nonzero rational endpoints. The same derivative calculation gives F_p(jr)∈p overline(Z); rationality now holds individually, so each value belongs to pZ. Primes avoiding a_0c g(0) make the weighted endpoint sum a nonzero integer. Summing the integral identities with integer weights a_j gives the same factorial decay, with a fixed additional factor Σ|a_j|. Negative r only reverses segment orientation and does not affect the absolute-value bound. This is a complete classical Hermite-type argument; the account accurately identifies modular nonvanishing as the extra device beyond the short positive-integral irrationality proof.
