# Fuss–Catalan valuations: forced carries and a tree quotient

Status: independently derived theorem and complete elementary proof below; finite computations independently check the formulas. This is **background research**, not a claimed solution of official problem #137: its exact statement has not been retrieved. The supplied URL returns a generic application shell, its API returns 403, and a separate web-tool open/search did not recover it. No Claude notes were read. Historical novelty is unresolved; the basic valuation identity is a direct consequence of classical factorial valuations/Kummer's theorem.

## 1. Exact theorem (checked proof)

For integers m>=1, n>=0 and prime p put

C_m(n) = ((m+1)n)! / (n! (mn)! (mn+1)),
a = v_p(mn+1).

Add n and mn in base p. A carry at position j means the carry leaving the digit of weight p^j; j starts at 0. Then

**Theorem.** The positions 0,...,a-1 all carry, and

v_p(C_m(n)) = number of carrying positions j>=a.

In particular C_m(n) is a p-adic unit precisely when the only carries are the a forced initial carries. If mn+1 is a power of p, C_m(n) is a p-adic unit.

### Proof

Let c_j be the carry entering position j, with c_0=0, and let x_j,y_j,z_j be digits of n,mn,(m+1)n. Thus x_j+y_j+c_j=z_j+p c_{j+1}, with each c_j in {0,1}. Summing over all digits gives

s_p(n)+s_p(mn)-s_p((m+1)n) = (p-1) sum_{j>=1} c_j.

To prove the factorial formula rather than assume it, count the multiples of p^k among 1,...,N:

v_p(N!) = sum_{k>=1} floor(N/p^k) = (N-s_p(N))/(p-1).

Subtracting three such formulas shows that v_p(binomial((m+1)n,n)) is exactly the number of carries. Division by mn+1 subtracts a.

If a>0 then mn is congruent to -1 modulo p^a: its first a digits are p-1. Also n is a unit modulo p (otherwise mn could not be -1 modulo p). The units addition therefore carries. Every subsequent position below a has a digit p-1 from mn and an incoming carry, so carries again. These are exactly a guaranteed carries, establishing the formula. For a=0 the assertion about initial carries is empty. For n=0, C_m(0)=1 and all counts are zero.

If mn+1=p^a with a>=1, mn=p^a-1 and 0<n<p^a. At position a, both summands have digit zero, the incoming carry is 1<p, and no further carry occurs. The a=0 case means n=0. This proves the corollary. QED.

An equivalent floor formula, useful for independent checking, is

v_p(C_m(n)) = sum_{k>a} [floor((m+1)n/p^k)-floor(n/p^k)-floor(mn/p^k)].

Each summand is zero or one. The omitted first a summands are one because mn modulo p^k equals p^k-1 and n modulo p^k is nonzero.

### Examples

|m|n|p|C_m(n)|a|all carry positions|valuation|
|--|--|--|--|--|--|--|
|1|7|2|429|3|0,1,2|0|
|2|4|3|55|2|0,1|0|
|2|7|3|7752|1|0,1|1|
|3|5|2|969|4|0,1,2,3|0|
|2|10|3|1430715|1|0,2|1|

The last row warns that surviving carries need not form one uninterrupted chain.

## 2. A genuinely combinatorial integrality proof, and its limitation

Set L=(m+1)n+1. Consider words with n letters of value m and mn+1 letters of value -1, hence total sum -1. Rotation acts freely: if a word were t repetitions of a shorter word, its total sum would be divisible by t, forcing t=1.

**Cycle lemma in this setting.** Each rotation orbit contains exactly one word whose first L-1 partial sums are nonnegative.

Proof: for partial sums S_0=0,...,S_L=-1, choose the first minimum among S_0,...,S_{L-1}, at index j. The word starting after j has nonnegative partial sums before its last step. Before wrapping this follows from minimality; after wrapping the partial sum is -1-S_j+S_k with k<j, and first-minimum plus integrality gives S_k>=S_j+1. Conversely any good starting index j must have S_k>=S_j for k>j (k<L) and S_k>=S_j+1 for k<j, making it that same first minimum. QED.

A good word is the preorder encoding of a plane full (m+1)-ary tree with n internal vertices: an internal vertex consumes one pending slot and creates m+1, net change m; a leaf consumes a slot, net change -1. Start with one slot. The partial-sum condition says a slot always remains until the last vertex. Reconstruction by filling the next pending preorder slot proves a bijection.

Thus the tree count is

binomial(L,n)/L = binomial((m+1)n,n)/(mn+1) = C_m(n).

This proves integrality without factorial valuations. It also explains division as forgetting the cyclic starting point, but it **does not yet give a combinatorial proof of the precise p-adic valuation**: knowing every word orbit has L elements does not by itself exhibit p-group orbits among trees corresponding to the surviving carries. A claim that the cycle lemma alone solves the valuation problem would leave this gap.

## 3. Concrete next constructions

1. Seek an explicit free action of a group of order p^{number of surviving carries} on an appropriate set of tree encodings, together with a quotient of size prime to p. The carrying positions depend only on m,n,p, so this asks for an action uniformly on all trees of a fixed size. An arbitrary partition exists from arithmetic divisibility, but would be circular; the action must be defined from tree structure.
2. Build a finite digit transducer for p-unit indices at fixed m,p. Keep multiplication carry u (0<=u<m), addition carry c in {0,1}, and a flag t recording whether all digits of mn seen so far equal p-1. On input digit d of n, compute e=(md+u) mod p, u'=floor((md+u)/p), c'=floor((d+e+c)/p), t'=t AND (e=p-1). The digit contributes c'-1_{t'} to v_p(C_m(n)). For valid inputs this contribution is nonnegative by the forced-carry theorem. Flush with zero digits until multiplication and addition carries vanish. This is a bounded-state, constructive route to regular-language classifications of p-unit indices. Its novelty and usefulness relative to official #137 remain unassessed.
3. Compare the two tree formulas: the carry count for n+mn loses v_p(mn+1), while that for n+(mn+1) loses v_p((m+1)n+1). A direct bijection between the corresponding surviving-carry descriptions could expose the role of pointing an extra leaf. At present this is a project direction, not a constructed bijection.

## Sources and provenance

The proofs above were developed independently before consulting mathematical research drafts. George N. Raney, *Functional composition patterns and power series reversion*, Transactions AMS 94 (1960), 441–451, DOI https://doi.org/10.1090/S0002-9947-1960-0114765-9 is a primary-source lead for the cycle-lemma framework. Crossref metadata was retrieved and verified; full text was not accessible through the web open attempted in this session, so no claim about a particular theorem number or exact formulation in that paper is made. Kummer's original paper and a primary source for factorial valuations should be retrieved before a final literature account. The digit-sum/factorial and cycle-lemma arguments are supplied in full above and do not depend on an unread citation.
