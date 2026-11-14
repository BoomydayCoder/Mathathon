# One more p-adic digit from the next tree orbits (odd primes)

Status: complete proposed refinement of the proved prime-power-arity Sylow theorem, with an explicit classification of the next orbit size. Independent audit passed; see `next-digit-independent-review.md`. Finite arithmetic and recursive orbit checks accompany it. Historical novelty is unresolved; the arithmetic congruence may be obtainable from classical factorial-unit congruences, and no claim of a new arithmetic formula is made.

Let p be an odd prime, q=p^e, n>=0, and write

L=(q-1)n+1=sum_j d_j p^j, 0<=d_j<p,
D=(sum_j d_j-1)/(p-1),
H_r=sum_{a=1}^r a^{-1} in F_p, H_0=0,
S=sum_{j>=0} d_{j+1} H_{d_j} in F_p,

with unlisted digits zero. Let T_q(n)=binom(qn,n)/L.

**Theorem.** The normalized count has the following second digit:

T_q(n)/p^D congruent ((p-1)!)^D / product_j d_j! * (1-p S) modulo p^2.

All factorials in the denominator are units modulo p^2. In p S, any lift of S from F_p gives the same residue. The formula is independent of e once L and its digits are specified.

The proof below refines the exact orbit enumeration rather than appealing to factorial valuations. It also gives an exact integer formula for the number of orbits with cardinality p^{D+1}.

## 1. Define the orbit excess and exclude repeated nonperfect blocks

Use the full micro p-ary wreath group on expanded q-trees, as in the complete proof. For a p-ary subtree with leaf count L_v and orbit exponent delta_v, put

D_v=(s_p(L_v)-1)/(p-1), E_v=delta_v-D_v>=0.

At a root whose p child orbit types are not all identical,

E_v=sum_i E_i + [sum_i s_p(L_i)-s_p(sum_i L_i)]/(p-1).

The second term is the total carry count, hence nonnegative. If the p child orbit types are identical, with common L_0,delta_0, then

E_v=p E_0+(p-1)D_0.

A repeated-child vertex which is not part of a perfect subtree therefore has excess at least p-1: if D_0>=1 this is immediate; if D_0=0 then either delta_0=0, giving a perfect child and hence a perfect parent, or E_0>=1, giving excess at least p. The excess cannot decrease as one moves to ancestors, by the two formulas above.

Because p is odd, p-1>=2. Thus an entire tree with excess one has no repeated-nonperfect vertex. Contract its maximal perfect subtrees to leaves labeled by their heights. At every remaining skeleton internal vertex, the p child orbit types are not all equal, so its local cyclic period is p. Consequently the skeleton orbit exponent equals its number K of internal vertices.

## 2. The only possible leaf profile is one uncarry

Let a_j be the number of leaves labeled j in the contracted skeleton. Their expansion gives L=sum_j a_j p^j and the full p-ary skeleton identity gives

K=(sum_j a_j-1)/(p-1).

For excess one, K=D+1, so sum_j a_j-s_p(L)=p-1. Normalizing the a_j into base-p digits performs exactly one elementary carry: if c_j are carries in this normalization, then sum_j a_j-s_p(L)=(p-1)sum_j c_j, and the c_j are nonnegative integers. Hence precisely one c_j is one and all others are zero. Equivalently there is one j with d_{j+1}>=1 and

a_j=d_j+p, a_{j+1}=d_{j+1}-1, a_i=d_i otherwise.       (U_j)

These are the only possible excess-one profiles. They are disjoint profiles for different j.

Conversely a skeleton with profile (U_j) can fail to be a canonical contraction only by having a vertex with p identical leaf children labeled j. Indeed any repeated child orbit type would repeat each label occurring in it p times. Every label except j occurs fewer than p times, so the repeated children would have to use label j only. If such a child were non-leaf, it would itself have at least p leaves, requiring at least p^2 copies of label j globally; but a_j<=2p-1<p^2. Therefore the only repeated type is a single labeled leaf. At most one such star can exist since a_j<2p.

If no such star exists, all internal vertices are nonuniform and no larger perfect merger can occur (inspect the deepest skeleton vertex of any claimed perfect merger). Its expansion is a tree of exact exponent K=D+1.

## 3. Count and subtract the forbidden stars

For any valid leaf profile a, let N(a) be the number of plane labeled skeletons satisfying the depth condition depth+label=0 modulo e. The complete cycle-splicing proof gives

N(a)=product_{i=0}^{e-1}(p k_i(a))! / [product_i k_i(a)! product_j a_j!],

where l_i(a)=sum_{j congruent -i mod e}a_j and

k_i(a)+l_i(a)=delta_{i0}+p k_{i-1}(a).

This count is valid even when some label multiplicity reaches p; the digit restriction was needed for orbit freeness, not for the labeled skeleton enumeration. For the digit profile d write N_0=N(d) and M=N_0/p^D, the exact minimum-orbit count.

A forbidden star in profile (U_j) contracts to a leaf labeled j+1, yielding a digit-profile skeleton with one distinguished leaf of that label. Conversely expand any distinguished label-(j+1) leaf into p leaves labeled j. This preserves depth+label modulo e and is a bijection, since the forbidden star is unique. The number of forbidden plane skeletons is exactly d_{j+1}N_0.

Every remaining skeleton has orbit size p^{D+1}. Therefore the exact number of next-size orbits is

**M_1=sum_{j:d_{j+1}>0} [N(U_j)-d_{j+1}N_0]/p^{D+1}.**       (1)

This is a nonnegative integer count, not just a congruence. In particular the normalized total satisfies

T_q(n)/p^D = M+p M_1 modulo p^2.                         (2)

## 4. Evaluate the next-orbit residue by unit products

For profile (U_j), the internal count increases by one in type t=-j-1 modulo e and stays unchanged elsewhere: this is exactly inserting the forbidden-star vertex before contraction. Let k=k_t(d), d=d_j and c=d_{j+1}. The exact ratio of factorial counts is

N(U_j)/N_0
= c * [(p(k+1))!/(p k)!] / [(k+1)(d+p)!/d!]
= c * [p product_{r=1}^{p-1}(pk+r)] / product_{r=1}^p(d+r).       (3)

This ratio is a p-adic unit, even when k+1 itself is divisible by p: the cancellation of p(k+1)/(k+1)=p is an exact identity before reduction. The denominator interval d+1,...,d+p contains the number p exactly once. Cancel it. Its remaining product is

(p-1)! product_{r=1}^d(1+p/r) = (p-1)!(1+p H_d) modulo p^2.

The numerator's unit product is

product_{r=1}^{p-1}(pk+r)=(p-1)!(1+p k H_{p-1})=(p-1)! modulo p^2,

because H_{p-1}=0 modulo p for odd p (pair r with p-r). Thus

N(U_j)/N_0 = c(1-p H_d) modulo p^2.

Substitute into (1), remembering that N_0=p^D M:

M_1 = -M sum_j d_{j+1}H_{d_j} = -M S modulo p.            (4)

To determine M itself modulo p^2, use its exact count formula from the complete proof:

M= product_i [(p k_i)!/(p^{k_i}k_i!)] / product_j d_j!.

Each bracket is exactly the product of p-units up to p k_i. Split it into k_i consecutive blocks of p; each block's unit product is (p-1)! modulo p^2 by the same H_{p-1}=0 identity. Since sum_i k_i=D,

M=((p-1)!)^D/product_j d_j! modulo p^2.                   (5)

Equations (2), (4), (5) prove the theorem. No factorial-valuation formula was used. QED.

## 5. Boundary cases and examples

* n=0 gives L=1,D=0,S=0,M=1,M_1=0, so the formula gives one.
* If L=q^h (a perfect-tree size), D=0 and all digits except one are zero. Every possible uncarry has d_j=0, so H_{d_j}=0 and S=0. In fact its skeleton is precisely the forbidden star, giving no next-size orbit. Thus T_q((q^h-1)/(q-1))=1 modulo p^2 for odd p.
* For p=q=3,n=2, L=5=(12)_3, D=1, S=H_2=0 modulo3. The theorem gives T_3(2)/3=1 modulo9, as T_3(2)=3.
* For p=q=3,n=3, L=7=(21)_3, D=1, S=2H_1=2. The theorem gives 2!/(1!2!)*(1-3*2)=4 modulo9, as T_3(3)=12 and 12/3=4.
* Oddness is essential in this derivation. At p=q=2,n=2, blindly using the displayed formula gives 3 modulo4, whereas C_2/2=1. Both the harmonic-block simplification and the exclusion of repeated-nonperfect excess-one blocks require separate treatment at p=2; that refinement is not claimed here.

No theorem about the next digits modulo p^3 is asserted. At that precision, several uncarry operations and repeated-nonperfect configurations already interact (especially p=3), so the present clean classification should not be extrapolated.
