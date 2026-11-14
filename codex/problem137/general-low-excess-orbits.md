# The low-excess orbit recipe: all normalized digits below the first repeated-block mechanism

Status: complete general proof draft based on the proved Sylow action/cycle-splicing count. It unifies minimum, one-uncarry and excess-two counts. Independent review of this generalization is requested; the earlier exact valuation and odd-prime second-digit theorem already have independent reviews. No historical novelty claim.

Let p be prime, q=p^e, n>=0,

L=(q-1)n+1=sum_{i=0}^{h-1}d_i p^i, 0<=d_i<p,
D=(sum_i d_i-1)/(p-1), T_q(n)=binom(qn,n)/L.

Act by the full micro p-ary wreath p-group on the depth-constrained p-ary expansions of q-trees. Let M_k count its orbits of cardinality p^{D+k}.

**Theorem.** For every 0<=k<p-1, M_k is given by the finite carry-profile and visible-star inclusion–exclusion recipe below. Consequently

T_q(n)/p^D congruent sum_{k=0}^{p-2} p^k M_k modulo p^{p-1}.

Thus one combinatorial recipe determines all normalized digits through that modulus. It specializes to the leading digit for p=2, the reviewed two-digit formula for p=3, and at least the first four normalized digits for p>=5. The bound is structural: a repeated nonperfect child block can occur at excess exactly p-1, and an explicit family below shows that the recipe cannot simply be extended past that boundary.

## 1. Enumerate the raw profiles

For fixed k, enumerate finite nonnegative carry vectors c=(c_0,...,c_{h-2}) with sum_i c_i=k, taking c_{-1}=c_{h-1}=0. Set

a_i(c)=d_i+p c_i-c_{i-1}.                              (1)

Discard any vector for which a_i<0. Its raw skeleton has leaf label i repeated a_i times and

sum_i a_i p^i=L,
K=(sum_i a_i-1)/(p-1)=D+k

internal vertices. Conversely every finite nonnegative leaf profile expanding to L has a unique normalization carry vector satisfying (1). If the largest positive carry were c_j, then a_{j+1}=d_{j+1}-c_j>=0, so d_{j+1}>0. Therefore no carries beyond index h-2 need be considered.

For h>=2 there are at most binom(h+k-2,k) carry vectors before discarding negatives. For h=1 only k=0 is possible. At fixed p this is a finite polynomial-size search in the number of digits, with degree at most p-2.

## 2. The raw plane count, valid beyond digit multiplicities

For any retained profile a define

l_t(a)=sum_{i congruent -t mod e}a_i,
k_t(a)+l_t(a)=delta_{t0}+p k_{t-1}(a),

with t modulo e. The already proved cycle-splicing bijection gives

N(a)=product_t(p k_t(a))!/[product_t k_t(a)! product_i a_i!].    (2)

This counts plane labeled p-ary skeletons whose leaf depth plus label is divisible by e. It is valid for arbitrary a_i, including a_i>=p; those multiplicities do not invalidate the slot-matching bijection. The root-leaf case is included separately as N=1 when K=0. The prior flow proof supplies nonnegative integral k_t for all retained profiles.

## 3. Why absence of visible stars is sufficient in this range

For a micro subtree write its orbit exponent as delta and its excess as E=delta-(s_p(leaf count)-1)/(p-1). The proved root recursions are

E=sum E_child+carry count, when the p child orbit types are not all equal;
E=p E_child+(p-1)D_child, when they are all equal.

Hence a repeated-nonperfect vertex forces E>=p-1 at itself and at every ancestor. In a raw K=D+k skeleton, perfect leaf-label expansions contribute no local orbit factor, so delta<=K and E<=k. If k<p-1, a repeated-nonperfect vertex is impossible.

A remaining repeated-child vertex can only be part of a perfect merger. Such a merger exists iff the raw skeleton has a visible star: a vertex with p leaf children all carrying the same label. The nontrivial direction follows by inspecting the deepest skeleton vertex inside any alleged perfect merger: its children are skeleton leaves, and their perfect heights must agree. Thus no visible star means no merger at all; every skeleton internal vertex contributes p and delta=K=D+k.

Conversely, a genuine orbit of excess k<p-1 has no repeated-nonperfect vertex. Contract its maximal perfect subtrees. The resulting skeleton is star-free, has K=delta=D+k, and belongs to exactly one profile (1). Therefore counting star-free raw skeletons gives precisely all such orbits.

## 4. Full visible-star inclusion–exclusion

Fix a retained raw profile a. Select any subset of its visible stars, and let b_j count selected stars whose p children have label j. Because their children are leaves, distinct visible stars are disjoint; they cannot be nested. Necessarily

0<=b_j<=floor(a_j/p)<=c_j,                              (3)

where the last inequality follows from (1) and d_j<p. In particular the total number of selected stars is at most k. Set

a'_i=a_i-p b_i+b_{i-1}, b_{-1}=0.                     (4)

Simultaneously collapsing the selected stars produces a raw a'-profile skeleton with b_j distinguished leaves of label j+1, for each j. Conversely choose those leaves independently and expand them to p identically labeled leaf children. Selected leaves are distinct positions and no leaf is an ancestor of another, including when their labels are adjacent. Hence all created marked stars are simultaneously visible and disjoint. Unmarked stars may disappear or appear under collapse; the target is intentionally a raw skeleton, so this causes no loss of bijectivity.

The number of choices of marked leaves is

product_j binom(a'_{j+1}, b_j).

There is no missing b_j! factor: an inclusion–exclusion subset of stars is unordered within each height, and the binomial coefficient selects an unordered set of leaf positions. The bounds (3) also ensure a'_i>=0 and a'_{j+1}>=b_j. They automatically give c-b>=0, so the collapsed profile remains in the same finite carry family at lower total carry.

Therefore the exact number of star-free raw skeletons is

**F(a)=sum_{0<=b_j<=floor(a_j/p)} (-1)^{sum b_j} N(a') product_j binom(a'_{j+1},b_j),**    (5)

where a' is (4). Only indices j<h-1 need occur, since a_{h-1}<p. The number of summands is at most product_j(c_j+1)<=2^k. Equation (5) is ordinary inclusion–exclusion on the actual finite set of visible stars, justified by the explicit collapse/expansion bijection.

Each star-free skeleton orbit has exactly p^{D+k} elements, by Section 3. Distinct retained carry profiles are disjoint canonical contractions. Thus

**M_k=sum_{c:sum c=k, a(c)>=0} F(a(c))/p^{D+k}.**        (6)

Each quotient is a nonnegative integer. This proves the theorem and the normalized congruence by omitting higher orbit sizes, all divisible by p^{D+p-1}.

## 5. Agreement with the earlier formulas

For k=0 there are no stars and F(d)=N(d), recovering the exact minimum-orbit count.

For k=1 the vector has c_j=1. The only nonzero b is b_j=1, giving F=N(a)-(a_{j+1}+1)N(d)=N(a)-d_{j+1}N(d), the reviewed one-uncarry formula.

For k=2, the zero and single-star terms are those in excess-two-orbits.md. For c_j=2, selecting both stars gives binom(d_{j+1},2)N(d). For c_j=c_l=1 it gives d_{j+1}d_{l+1}N(d). If l=j+1 and the middle digit is zero, the upper star cannot be simultaneously visible; the corresponding binomial/pair count is zero or its bound (3) excludes it. Thus the general formula recovers the nested-merger-safe excess-two formula exactly.

## 6. The first excluded mechanism is real, for every prime

For q=p, take the nonperfect p-ary tree with two internal vertices: its root has one perfect child of height one and p-1 leaf children. Its child-placement orbit has size p, leaf count 2p-1, and D_child=1. Make a new root with p identical copies of that child orbit type. It has

n=2p+1, L=p(2p-1), D=1, delta=p,
E=delta-D=p-1.

After maximal-perfect contraction its leaf profile is a_0=p(p-1), a_1=p; it has K=p+1 internal vertices and total normalization carry p, rather than p-1. The uniform nonperfect root removes one local p factor, so delta=K-1. It is consequently missed by the carry-sum-(p-1) star-free recipe. At p=3 this is exactly the independently enumerated n=7 extra orbit of size27.

The boundary also occurs for every e>=2. Let h=e-2. Form a nonperfect child root with p-1 perfect branches of height h and one perfect branch of height h+e, then put p identical copies under a new root. All final leaf depths are e or2e, so this is a valid expansion of a q=p^e tree. The child has delta=1,D_child=1, and the whole tree again has delta=p,D=1,E=p-1. Its q-tree size is n=p^{e-1}+1. Thus the depth restriction does not remove the first excluded mechanism.

Past that threshold one must track repeated-nonperfect skeleton vertices and their automorphisms in addition to uncarry profiles and perfect-star collapse. The present theorem deliberately stops before that new phenomenon.

## Scope and verification

This is an explicit finite combinatorial recipe for all low-excess orbit counts; it does not claim a compact closed harmonic formula for every digit or historical novelty. Its factorials enumerate typed labelings via a proved bijection, while its p-divisibility conclusions arise from actual orbit sizes. Arithmetic normalized counts are used only as independent checks.

check_general_low_excess.py compares (5)–(6) with the reviewed k=1 and k=2 formulas, direct normalized counts modulo p^{p-1}, and independently generated Sylow orbit histograms. Boundary examples are checked separately and not represented as belonging to the theorem.
