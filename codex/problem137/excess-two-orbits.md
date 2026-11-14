# Excess-two orbits: a combinatorial normalized-unit formula modulo p^3

Status: complete proof draft building on the independently reviewed next-digit-orbits.md and the proved cycle-splicing enumeration. It gives an exact combinatorial count of the next-next orbit size, hence a finite digit-profile formula modulo p^3. It does not claim a simplified harmonic closed form or historical novelty. Independent audit of this new refinement remains appropriate.

Let p>=5 be prime, q=p^e, n>=0, L=(q-1)n+1=sum_i d_i p^i, and D=(sum_i d_i-1)/(p-1). The tree count is T_q(n). The group is the full micro p-ary wreath group on expanded q-trees. Let M_r be the number of its orbits of cardinality p^{D+r}. The complete theorem and the reviewed second-digit theorem already determine M_0 and M_1 exactly. Here we determine M_2 by explicit objects.

**Result.** Define the finite profile count N(a) and star-free count F(a) below. Then

M_2=sum_{c>=0, sum_i c_i=2, a(c)>=0} F(a(c))/p^{D+2},

and therefore

T_q(n)/p^D congruent M_0+p M_1+p^2 M_2 modulo p^3.

Every quotient in the sum is a nonnegative integer counting actual orbits. This is an orbit refinement, not a factorial-valuation derivation recast in combinatorial language.

## 1. Why only two uncarries occur

For a micro p-ary subtree put E=delta-(s_p(leaf count)-1)/(p-1), where p^delta is its orbit size. The reviewed excess recursions are

E=sum E_child + carry count, if child orbit types are not all equal;
E=p E_child+(p-1)D_child, if all p child orbit types agree.

A repeated-nonperfect vertex has excess at least p-1, and excess cannot decrease at ancestors. Since p>=5, no tree of excess two can contain such a vertex. Contract maximal perfect subtrees into leaves labeled by their heights. Every remaining internal vertex then has local orbit factor p, so the skeleton has exactly K=D+2 internal vertices.

Let a_i count its leaves labeled i. Then L=sum a_i p^i and K=(sum a_i-1)/(p-1). Normalizing a into the base-p digits d has total carry count two. Equivalently there is a finite nonnegative carry vector c with sum c_i=2 such that

a_i=d_i+p c_i-c_{i-1}, with c_{-1}=0.                 (1)

All indices beyond the highest digit are taken as zero. In a finite valid vector the largest positive c_i lies below a positive digit d_{i+1}, so only indices below the highest digit need be considered. The possible vectors are one entry c_j=2, or two distinct entries c_j=c_k=1, subject to a_i>=0. These include adjacent uncarries, which must not be omitted.

## 2. Count raw depth-constrained skeletons

For any nonnegative profile a with sum a_i p^i=L, define

l_t(a)=sum_{i congruent -t mod e}a_i,
k_t(a)+l_t(a)=delta_{t0}+p k_{t-1}(a),

where t is read modulo e. The cycle-splicing bijection proves the number of plane labeled skeletons with leaf labels counted by a and depth+label=0 modulo e is

N(a)=product_t (p k_t(a))!/[product_t k_t(a)! product_i a_i!].    (2)

This is a combinatorial enumeration valid for arbitrary multiplicities, not only a_i<p. The previously proved nonnegative integral solution for the k_t applies to these profiles as well. For a=d put N_0=N(d)=p^D M_0.

A raw skeleton of a profile (1) has K=D+2 internal vertices, so its expansion has orbit exponent at most K and excess at most two. Hence it cannot contain a repeated-nonperfect vertex, by Section 1. Its only possible zero local orbit factors arise from perfect mergers. A merger exists iff there is a visible leaf star: a vertex whose p children are all leaves carrying the same label. Inspecting the deepest internal skeleton vertex of a claimed perfect merger proves this equivalence.

Consequently the desired excess-two trees correspond precisely to raw skeletons having no visible identical-leaf star. For such a skeleton every internal vertex contributes p, giving exact orbit size p^{D+2}.

## 3. Inclusion–exclusion on visible stars, including nested merger cases

A visible star labeled j consumes p copies of leaf label j. From (1), floor(a_j/p)<=c_j. Therefore a raw skeleton has at most sum c_j=2 visible stars. Distinct visible stars are disjoint, since their children are leaves; a perfect merger created *after* a collapse is not a second initially visible star.

For every j with a_j>=p define

b^{(j)}=a-p e_j+e_{j+1}.

Collapse a distinguished visible j-star. The result is a raw b^{(j)}-profile skeleton with one distinguished leaf labeled j+1. Conversely expand any such distinguished leaf. This is a bijection even when collapsing creates a new visible star higher in the tree: the target is deliberately a raw skeleton, not a canonical contraction. Thus the number of pairs (raw skeleton, distinguished visible j-star) is

(a_{j+1}+1) N(b^{(j)}).                              (3)

For a pair of visible stars, collapsing both yields the digit profile d because it removes both carries. If c_j=2, the pair corresponds to two distinct chosen leaves of label j+1 in a d-profile skeleton; its count is binom(d_{j+1},2)N_0. If c_j=c_k=1 with j<k, the count is d_{j+1}d_{k+1}N_0. This includes k=j+1: if the shared middle digit is zero, the pair count is zero, while a nested merger after one collapse may still occur and is already handled by (3).

Set

Z(c)=binom(d_{j+1},2) if c_j=2;
Z(c)=d_{j+1}d_{k+1} if c_j=c_k=1, j<k.

Inclusion–exclusion therefore gives the exact star-free count

**F(a)=N(a)-sum_{j:a_j>=p}(a_{j+1}+1)N(b^{(j)})+Z(c)N_0.**    (4)

There is no triple term. The quotient F(a)/p^{D+2} is the exact number of excess-two orbits of that canonical profile. Different carry vectors give different canonical profiles and cannot overlap. Summing proves the result.

## 4. A concrete fully checked example

Take p=q=5,n=4. Then L=17=(32)_5, d=(2,3), D=1. The only two-carry profile is c_0=2, giving a=(12,1). The raw skeleton has three internal vertices and thirteen leaves. Formula (2) gives N(a)=455. Its marked-star term is 2 N(7,2)=2*180=360; its pair term is binom(3,2)N(2,3)=3*10=30. Hence

F(a)=455-360+30=125,
M_2=125/5^3=1.

Direct enumeration of these 455 labeled plane skeletons finds 125 with no visible stars, 300 with one and 30 with two. The total tree-orbit decomposition is

T_5(4)=285=2*5+6*25+1*125.

The normalized unit is therefore 285/5=57 modulo125, explained by the three orbit sizes rather than by a factorial valuation calculation.

## 5. Nested-merger and perfect-size boundary check

For p=q=5 and L=25 (n=6,D=0), the two adjacent carries give a=(5,4,0). There are N(a)=630 raw skeletons. Only five contain a visible star. Collapsing that star creates another perfect merger, but it remains one initially visible star; the marked count is 5*N(0,5,0)=5. The pair term is zero, since d_1=0. Thus F(a)=625 and M_2=625/25=25. Direct raw enumeration confirms 625 star-free and five one-star skeletons.

For the same L but q=25 (n=1), depth residues restrict this raw profile to the five forbidden skeletons, giving F(a)=0 and M_2=0. The exact orbit distributions thus depend on the arity/depth constraint even when the leaf count is the same. Both total counts are nevertheless one modulo125 after normalization; no general arity-independence theorem is asserted here.

The n=0 case has no valid two-carry vector and M_2=0. Zero internal types are handled by the already proved cycle-splicing enumeration. The restriction p>=5 is material. At p=q=3,n=7, L=15=(120)_3,D=1. Take a root with three identical child orbit types, each a minimal nonperfect tree with two internal vertices (three possible orientations). The child orbit has size3, so the root orbit has size3^3=27 and excess two. Its contracted skeleton has four internal vertices and carry sum three, outside the two-uncarry class. The two valid two-uncarry profiles give 1+13=14 orbits of size27, whereas direct enumeration gives15; the missing orbit is precisely this repeated-nonperfect construction. Thus the p=3 obstruction is an explicit extra orbit, not a hypothetical caveat. No p=3 version is claimed.

## 6. What this result does and does not simplify

This gives a finite digit algorithm: enumerate O(h^2) two-carry vectors for h digits, evaluate the already combinatorial typed counts, and apply (4). Together with the minimum and one-uncarry counts it computes the normalized count modulo p^3 and gives a precise structural interpretation of every retained orbit class.

A short closed harmonic expression is not derived in this note. Classical higher p-adic factorial-unit congruences (including Wolstenholme-type block-product refinements for p>=5) can provide an independent arithmetic simplification/check; their existence must not be confused with a historical novelty claim for the congruence. The mathematical content developed here is the actual excess-two orbit classification and its exact visible-star inclusion–exclusion.
