# A proved weighted extension: ordinary tree depth

Status: complete consequence of the proved unweighted Sylow/minimum-orbit theorem. This is a different weight model from Postnikov–Sagan's left-edge depth and Konvalinka/Gao–Gu's non-right-edge depth. Historical novelty unassessed.

Let p be prime, q=p^e, n>=0, L=(q-1)n+1=sum_j d_j p^j, D=(sum_j d_j-1)/(p-1). For a sequence of integers b_0,b_1,... assign a full q-ary tree T the weight

w_b(T)=product_{internal vertices v} b_{depth(v)},
W_q(n;b)=sum_T w_b(T),

where root depth is zero and depth counts all edges. Suppose there is a fixed integer u with p not dividing u and b_h=u modulo p for every h>=0. Then

**Theorem.** v_p(W_q(n;b))=D, and

p^{-D} W_q(n;b) congruent u^n (-1)^D / product_j d_j! (mod p).

In particular, for binary trees (p=q=2), *any* sequence of odd ordinary-depth weights preserves the exact Catalan 2-adic valuation. A nonconstant example at every prime is b_h=1+p h, or b_h=u+p h^2.

**Proof.** The finite Sylow action on expanded p-ary trees preserves depth; a macro q-tree internal vertex lies at micro depth e*depth(v), so macro depth is preserved too. Consequently w_b is exactly constant on each orbit O; call its value w_O. Since T has n internal vertices, w_O=u^n modulo p. Every orbit has cardinality at least p^D in p-power order, and there are M minimum orbits of cardinality p^D, with

M congruent (-1)^D/product_j d_j! modulo p,

as proved by typed-skeleton enumeration and cycle splicing. Therefore

W_q(n;b)/p^D
= sum_{minimum O} w_O + sum_{higher O} (|O|/p^D)w_O
= u^n M modulo p.

The last residue is nonzero, proving both exactness and the displayed leading unit. Negative weights cause no difficulty: all orbit sums are integers and the nonzero residue prevents complete cancellation. QED.

More generally, any integer weight constant on the Sylow orbits satisfies the divisibility bound p^D. Its exact valuation equals D precisely when the sum of its minimum-orbit weights is nonzero modulo p. A common nonzero residue on all those orbits is a simple sufficient condition. This broader criterion should not be mistaken for a statement about arbitrary weights that vary within an orbit.

## Why this does not automatically extend the cited path-weight model

Postnikov–Sagan weight a binary-tree vertex v by b(l_v), where l_v counts left edges from the root; their path model is Dyck ascent height. Exchanging left and right children changes l_v, so weights are generally not constant on their symmetry orbits. Their Theorem 2.1 controls the *average* weight by hypotheses on higher finite differences: b(0) odd and 2^{r+1} dividing Delta^r b(x) for every r>=1. These hypotheses are stronger than oddness because they must compensate for within-orbit weight variation.

For example take b_0=1,b_1=3. At n=2 our ordinary-depth weighted binary count is 2*b_0*b_1=6, of valuation one. The left-edge weighted count is b_0^2+b_0*b_1=4, of valuation two. Thus 'all weights odd' cannot be transferred from the present model to theirs.

Konvalinka's q-tree weight and Gao–Gu's subsequent refinement count non-right edges (all but one distinguished child direction); these too vary under child permutations and require orbit averaging. Our depth-weighted theorem is a rigorous separate extension, not a claim to subsume their stronger path-weight results.

## Reproducible computation

check_depth_weights.py computes the entire weighted count independently from the tree recursion

F_h(z)=1+z b_h F_{h+1}(z)^q,

truncated at the maximum size requested. It verifies the valuation and leading unit for several prime-power arities and nonconstant positive/negative depth weight sequences. The proof above, rather than the finite tests, establishes the general theorem.
