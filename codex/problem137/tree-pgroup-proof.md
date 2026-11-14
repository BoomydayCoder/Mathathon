# Tree p-groups: a combinatorial unit proof and exact valuations at prime arity

Status: complete independent proof audit/development of the parent agent's proposed finite tree action. The p-unit theorem works for q=p^e. The exact valuation theorem below is proved for q=p. A counterexample shows why minimum orbit sizes for cyclic q-child rotations do not immediately extend the exact proof to e>1. Historical novelty is unresolved; these are developed proofs, not novelty claims. The construction is in the orbit-method tradition associated with combinatorial congruence work on trees; a precise Deutsch–Sagan literature comparison remains to be retrieved and checked.

## 1. Make the finite p-group action precise

Let q=p^e, n>=0. A plane full q-ary tree with n internal vertices has (q-1)n+1 leaves. Encode an internal vertex by the word over {0,...,q-1} specifying its path from the root. The internal-vertex set is prefix-closed. Any root-to-leaf path has at most n internal vertices, so every such tree embeds in the complete q-ary tree of height n. This is a finite common ambient tree, independent of the shape being acted on.

At each ambient vertex of depth less than n, allow the cyclic permutation of its q children, carrying along the entire descendant subtrees. These are automorphisms of the finite ambient tree. They generate the iterated wreath product

G_n = C_q wr C_q wr ... wr C_q (n levels),

with the unambiguous recursive definition G_0=1, G_h=G_{h-1}^q semidirect C_q. Its order is q^{1+q+...+q^{n-1}}, a power of p. It permutes the set of embedded full q-ary trees with exactly n internal vertices: automorphisms preserve depth, prefix relations, vertex arities, and the number of internal vertices.

**Fixed-point lemma.** A tree fixed by every element of G_n is perfect: for some h, every vertex of depth below h is internal and every vertex at depth h is a leaf.

Proof: G_n is transitive on the ambient words at each fixed depth d, by induction using root rotation and independent transitive actions inside the q branches. An invariant internal-vertex set therefore contains either every vertex at a given depth or none. Prefix closure then says that the occupied internal depths are exactly 0,...,h-1. Conversely a perfect tree is visibly fixed by all ambient automorphisms. QED.

Thus there is exactly one fixed tree of size n=(q^h-1)/(q-1), and none at any other size. This includes n=0,h=0; the ambient group G_0 is trivial.

## 2. Complete combinatorial p-unit classification for q=p^e

Every orbit under a finite p-group has p-power cardinality. Modulo p, all non-singleton orbits disappear and singleton orbits are precisely the globally fixed trees. The cycle-lemma count from valuation.md identifies the number of trees as

T_q(n)=binom(qn,n)/((q-1)n+1).

Therefore

T_q(n) congruent 1 (mod p) if n=(q^h-1)/(q-1),
T_q(n) congruent 0 (mod p) otherwise.

This strengthens the earlier arithmetic p-unit classification by determining the actual nonzero residue, and proves it by a direct finite group action. It uses neither factorial valuations nor Kummer's theorem. The arithmetic proof in prime-power-arity.md is an independent consistency check.

## 3. Orbit recursion

A rooted tree orbit is specified by the cyclic necklace of the q child-orbit types, since the q descendant groups act independently and root rotation then acts on the tuple of child types. If child orbit sizes are p^{d_1},...,p^{d_q}, and the tuple has least cyclic period ell=p^j, the root orbit has size

p^{d_1+...+d_q+j}.

Indeed independent descendant actions give the product of child orbit sizes, and ell distinct tuples of child-orbit types occur under root rotations. This supplies both a proof and an exact enumeration algorithm; a tuple period depends on child *orbit types*, not merely child sizes or subtree vertex counts.

## 4. Exact valuation theorem for prime arity q=p

Let L=(p-1)n+1, s=s_p(L), and D=(s-1)/(p-1), an integer since digit sums preserve residue modulo p-1.

**Theorem.** Every orbit of full p-ary trees with L leaves has size divisible by p^D, and the number M of orbits having size exactly p^D satisfies

M congruent (-1)^D / product_j d_j! (mod p),

where L=sum_j d_j p^j, 0<=d_j<p. In particular M is nonzero modulo p, and

v_p(T_p(n))=(s_p((p-1)n+1)-1)/(p-1).

This is an exact valuation proof by tree orbits. The small factorials d_j! are units; division in the residue formula is in F_p.

### 4a. Lower bound and equality cases

Induct on tree size. A leaf has L=1, s=1, D=0 and orbit exponent zero. Let a non-leaf root have child leaf counts L_i and child exponents delta_i.

If all p child orbit types agree, write their common values as L_0,delta_0. The root orbit exponent is p delta_0 and L=pL_0, so s_p(L)=s_p(L_0). The induction bound gives

p delta_0 >= p(s_p(L_0)-1)/(p-1) >= (s_p(L)-1)/(p-1).

Equality holds exactly when s_p(L_0)=1 and delta_0=0. An orbit of exponent zero is a fixed tree, hence perfect by the fixed-point lemma. Thus equality in this repeated-child case produces a perfect tree and nothing else.

If the child orbit types do not all agree, primeness of p makes their cyclic tuple have period p. Then

delta=1+sum_i delta_i
>=1+sum_i (s_p(L_i)-1)/(p-1)
=(sum_i s_p(L_i)-1)/(p-1)
>=(s_p(sum_i L_i)-1)/(p-1)=D.

The final inequality is the elementary fact that addition cannot increase the sum of base-p digits: each carry lowers it by p-1. Equality holds iff every child is itself an equality case and the addition of its p leaf counts is carry-free. This proof uses digit arithmetic for ordinary addition, not a valuation formula.

### 4b. Count the equality orbits

The equality recursion has a useful unique description. Contract every maximal perfect subtree to a leaf, recording its height j as the leaf label. A perfect subtree of height j has p^j leaves. The carry-free condition says label j occurs exactly d_j times, the digit of L, so fewer than p times. The contracted skeleton is a full p-ary tree with s labeled leaves and D internal vertices. Conversely, take any such skeleton with label j repeated d_j times and replace each labeled leaf by the corresponding perfect tree. This gives an equality tree.

There is no accidental merger of skeleton nodes into a perfect subtree: that would require p identical perfect child blocks, and hence at least p occurrences of one label, forbidden by d_j<p. More generally no skeleton vertex has p equal child orbit types, because that would repeat every leaf label in those child subtrees p times, again forbidden. Consequently every skeleton orbit has size exactly p^D. The construction and contraction commute with the group actions, yielding a bijection of the relevant orbits.

The number of plane full p-ary skeletons with these leaf-label multiplicities is

[binom(pD,D)/s] * [s!/product_j d_j!]
= (pD)!/(D! product_j d_j!).

The first bracket is the cycle-lemma full-tree count, since s=(p-1)D+1; the second independently assigns the multiset of labels to the s ordered leaf positions. Divide by the common orbit size p^D to obtain

M=(pD)!/(p^D D! product_j d_j!).

The integers divisible by p among 1,...,pD have product p^D D!, so the quotient (pD)!/(p^D D!) is exactly the product of all the remaining integers, not an approximation involving valuations. In consecutive blocks of p each product is (p-1)! congruent -1 modulo p (Wilson's elementary inverse-pairing proof). Therefore M congruent (-1)^D/product_j d_j!, nonzero. All other orbits have sizes divisible by p^{D+1}, so dividing the total count by p^D and reducing modulo p gives M, proving the exact valuation. QED.

For p=2 all d_j are 0 or 1, and M=(2D-1)!! is odd (with (-1)!!=1). The result is the familiar exact Catalan formula v_2(C_n)=s_2(n+1)-1, now with explicit minimum orbits represented by rooted nonplane binary skeletons on the distinct labels corresponding to the 1-bits of n+1.

## 5. Why the cyclic q=p^e action alone is insufficient for exact valuations

At q=4,p=2,n=28, the tree leaf count is L=3*28+1=85, with binary digit sum 4. The independently proved arithmetic formula gives v_2(T_4(28))=3.

For cyclic four-child rotations, however, the minimum orbit exponent is 2. To see this without a huge enumeration, n-1=27 child internal vertices is odd. A tuple of period 1 or 2 would repeat child sizes four or two times, giving an even total; thus the root period is 4. Every root orbit consequently has exponent at least 2.

Equality requires all four child orbits to have exponent zero, hence all child trees perfect. Their allowed internal sizes below 28 are 0,1,5,21. The only four-term multiset with sum 27 is {0,1,5,21}: one needs 21, then 5, then 1, then 0. Its four entries are distinct. There are 4!/4=6 cyclic necklaces, each giving an orbit of size 4. Hence the minimum-orbit count is 6, divisible by 2, and its exponent 2 is strictly smaller than the true valuation 3.

Thus the tempting claim 'the valuation is the smallest orbit exponent' fails even in a simple prime-power-arity case. One must prove that the number of minimum orbits is a p-unit, or use a richer action.

A promising richer group replaces local C_{p^e} by a Sylow p-subgroup of S_{p^e}, obtained by expanding each q-ary vertex into a perfect p-ary tree of height e. In the q=4,n=28 example the local Sylow group has order 8, producing three minimum orbits of size 8 instead of six of size 4. The expanded trees form a stable subset of the full p-ary tree action: their leaves lie at depths divisible by e. A general proof that this restricted subset always has a p-unit number of minimum orbits remains open in these notes; it is not supplied by the prime-arity skeleton count without an additional depth-congruence argument.

## Literature and scope

This is a proof construction in the established group-action approach to combinatorial divisibility. A precise comparison with Deutsch–Sagan or related work must be made before presenting it as original. No literature search failure can establish novelty. The exact official #137 wording should still be checked before deciding whether prime arity plus prime-power-arity unit classification meets its intended scope.

## Later source verification and stronger extension

The primary Deutsch–Sagan paper was subsequently retrieved from arXiv math/0407326 and Section 2 read. Its Theorem 2.1/Lemma 2.3 give the binary wreath-action proof and minimum-orbit double factorial reconstructed above; this binary result is therefore explicitly attributed to Deutsch–Sagan (JNT117(2006),191–215, DOI10.1016/j.jnt.2005.06.005), not presented as new. See ../sources/deutsch-sagan-2006.pdf.

The previously open full-Sylow extension for q=p^e is now proved/audited in prime-power-sylow-audit.md using depth-constrained skeleton enumeration and the stated Lagrange–Good formula. That later proof resolves the extension gap identified in Section 5 here; the cyclic-group counterexample remains valid and useful.
