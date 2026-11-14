# A cycle-splicing bijection counts the typed skeletons

**Status.** Complete elementary bijection, independently checked by the valuation agent and by exhaustive small slot-matching computations. It replaces the multivariable Lagrange-inversion step in `prime-power-orbit-review.md`. Together with that note's group action and minimal-orbit classification, the argument is now combinatorial throughout its enumeration steps. This is a developed proof, not a historical-novelty claim; comparison with Konvalinka's primary paper remains necessary.

## 1. Alignment with the exact official problem

The retrieved statement `codex/sources/problem137-official.md` asks specifically for an exact valuation proof when the arity is q=p^e, preferably by tree symmetries and without factorial valuation formulas. It explicitly identifies exactness beyond a group-action divisibility bound as the main obstacle.

Our stable expanded-tree action is in precisely that scope. It proves every orbit has size divisible by p^D, classifies the orbits of size exactly p^D, and shows their number is a p-unit. The only previously algebraic enumeration step was the typed skeleton count. The following bijection removes that dependency. Factorials still occur as counts of bijections between finite sets; no formula for v_p(m!) or valuation of a binomial coefficient is invoked.

## 2. The counting lemma

Fix integers p>=2 and e>=1. Types are elements of Z/eZ. Let k_i and l_i be nonnegative integers satisfying

    k_i+l_i=delta_(i,0)+p k_(i-1),                    (1)

and suppose k_0>0. Consider plane rooted trees whose internal vertices of type i number k_i, whose leaves of type i number l_i, and in which every internal vertex has p ordered children, all of its next type i+1.

For a labelled-leaf version, prescribe individual label classes j, each confined to one type, with multiplicities d_j summing to l_i in that type. Then the number N of trees with those repeated leaf labels, and otherwise unlabelled internal vertices, is

    N=product_i (p k_i)!/[product_i k_i! product_j d_j!].    (2)

For type-only leaves take a single label class of multiplicity l_i in each type, obtaining the denominator product_i k_i! l_i!.

This lemma does not require p to be prime, and it allows some k_i=0. Primality and the digit restrictions enter the separate orbit argument. The zero-internal-vertex single leaf is handled separately and also agrees with the factorial expression.

## 3. Label the vertices and count arbitrary slot matchings

Give every internal vertex a distinct name within its type, and distinguish all copies of each leaf label. Every internal vertex has p distinct child slots, numbered1,...,p. A slot on a type-(i-1) internal vertex accepts a type-i child.

Fix one type-0 internal vertex r, and order the type-0 internal names so that r is largest. All vertices except r must be assigned a parent slot. By(1), the number of those vertices of type i is exactly p k_(i-1), the number of available slots for that type. Therefore the number of arbitrary assignments, each using every slot exactly once, is

    product_i (p k_(i-1))!=product_i (p k_i)!.         (3)

An assignment need not be a tree. Orient each assigned edge from the child to its parent. Every vertex except r has one outgoing parent edge; r has none. Hence there is one rooted-tree component containing r, and every other connected component has exactly one directed cycle with rooted branches attached.

Every directed cycle visits type0: moving along a parent edge subtracts one from the type, so a closed cycle has length a positive multiple of e and passes through every residue. Cycle vertices are internal, since leaves have no child slots. Thus each cycle has a uniquely largest type-0 internal name.

## 4. Splice the cyclic components into a single tree

For each cycle choose its largest type-0 vertex c. Order those vertices decreasingly,

    c_1>c_2>...>c_s.

Let sigma_i be the parent slot currently occupied by c_i. Initially the rooted component has root r. Perform the following operations in that order:

1. Remove the parent edge from c_i, leaving slot sigma_i empty. Its cyclic component becomes a rooted tree with root c_i.
2. Put the current root into slot sigma_i. This joins the two components into one rooted tree.
3. Declare c_i the new current root.

Equivalently, replace the edges from c_i to sigma_i by

    r -> sigma_1, c_1 -> sigma_2, ..., c_(s-1) -> sigma_s,

leaving c_s without a parent. With no cyclic components, do nothing.

Every moved child has type0, and each vacated slot accepts type0. Slot order, vertex types, and every vertex's number of children are preserved. At each step two distinct components merge, so the final result is a connected plane rooted tree with some type-0 internal vertex as root. Its root need not be r; allowing its name to vary is essential to the count.

## 5. The inverse: records along one path

Start with any connected labelled plane tree having the prescribed vertex types and child slots. Its root has type0. Consider the unique path from that root down to the distinguished vertex r. Read only the type-0 vertices on this path, and mark their successive record maxima. Write them in path order as

    b_1<b_2<...<b_s<r.

The root is the first record unless it is r, and r is the last record because it is the largest type-0 internal name in the entire tree. Between one record and the next, all type-0 names are smaller than the first record.

For each consecutive pair consisting of a record b_i and the next record (b_(i+1), or r at the end), do the following simultaneously:

1. Remove the edge attaching the next record to its parent slot sigma.
2. Attach b_i itself to that same slot sigma.

The path segment from b_i down to the parent of the next record closes into a directed cycle. Its largest type-0 vertex is b_i. The vertex r becomes parentless, and retains its own rooted-tree component. All branches off the distinguished path remain unchanged.

The reconstructed cycle maxima increase along the original root-to-r path, so their decreasing order is exactly the order used by the forward splicing. Each segment has no larger type-0 name than its initial record. Therefore the forward operation reconstructs the given tree exactly. Conversely, on a tree produced by splicing, the root-to-r path runs through the opened cycles in increasing order of their maxima; these and r are precisely its record maxima. The two operations are inverses.

Self-loops when e=1 cause no exception: opening a loop produces a one-vertex cycle segment, and the inverse closes it at the same slot. If some k_i=0, no directed cycle can pass through every type; all assignments are already trees, and the bijection is the identity. The balance equations then force k_0=1, in agreement with the root count.

Thus arbitrary assignments with the **fixed missing-parent name r** are in bijection with all labelled rooted trees whose type-0 root name is **allowed to vary**. This is a typed, slot-preserving form of the classical functional-digraph cycle-splicing construction; no claim that the general bijection idea is new is intended.

## 6. Forget the artificial names

For each plane tree with the intended repeated leaf labels, the number of assignments of distinct internal names is product_i k_i!. The number of assignments of distinct copy names within each leaf-label class is product_j d_j!. A plane tree has individually specified child positions, so there is no further automorphism quotient. Its root may receive any type-0 internal name, already included in k_0!.

Consequently the count in(3), divided by those naming factors, is exactly(2). This proves the typed skeleton formula by an explicit finite bijection. The factor1/k_0 that appeared during Lagrange inversion is now explained by the change from a fixed root name in arbitrary assignments to a variable root name among genuine trees.

## 7. Substitution into the prime-power valuation proof

For q=p^e and L=(q-1)n+1=sum_j d_j p^j, the minimal-orbit skeleton has leaf label j exactly d_j times. That leaf is allowed at depth i precisely when i+j=0 modulo e. Thus take

    l_i=sum_(j=-i mod e)d_j,
    (q-1)k_i=sum_(r=0)^(e-1) p^r l_(i-r)-p^i.

The earlier elementary flow argument proves these k_i are nonnegative integers and sum to D=(sum_j d_j-1)/(p-1). If D>0, k_0>0, so the bijection applies. If D=0, the unique compatible single leaf is immediate.

The stable p-group action has already proved that every minimum orbit has p^D elements. Their number is therefore

    M=N/p^D
      =[product_i product_(1<=a<=p k_i, p does not divide a) a]
           /product_j d_j!.

This exact product identity removes the multiples p,2p,...,p k_i from each slot-matching factorial. All remaining factors and all d_j! are p-units. Pairing inverse nonzero residues modulo p yields

    M = (-1)^D/product_j d_j! !=0 (mod p).

All larger orbits disappear after dividing the total count by p^D and reducing modulo p. Hence

    v_p(T_(p^e)(n))=(s_p((p^e-1)n+1)-1)/(p-1),

with the stronger leading-unit congruence stated in the review. Divisibility came from orbit sizes, and exactness came from counting the minimum orbits. The proof no longer imports multivariable Lagrange inversion or factorial valuation formulas.

## 8. Verification and scope

`check_cycle_splicing.py` exhaustively enumerates small slot assignments, checks that splicing produces a connected tree, verifies the inverse exactly, verifies that images are distinct, and separately enumerates every possible root name to check surjectivity. Cases include e=1 self-loops, multiple cyclic components, e>1, p=3, and a missing internal type. It does not use formula(2) to classify an assignment as a tree.

The independent depth-profile computation in `check_prime_power_orbit_review.py` supplies a different check of the final skeleton count and leading-unit congruence. These finite checks supplement the universal bijection, not replace it.

The exact official problem now matches this prime-power theorem. A final report should still separate the mathematical proof from historical assessment: compare the primary Konvalinka and Deutsch–Sagan constructions before describing the result as new, and seek further independent expert review of the complete combined argument.
