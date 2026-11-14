# Independent audit of the excess-two orbit formula

Reviewed excess-two-orbits.md, its checker, the predecessor
next-digit-orbits.md, and the typed cycle-splicing counting statement.
**Verdict: accepted, with no proof gap found.** The restriction \(p\ge5\)
is essential to this particular two-uncarry classification.

This review does not claim historical novelty. The result is an exact
orbit classification and finite combinatorial formula for the normalized
count modulo \(p^3\), not yet a short harmonic closed form.

## 1. Excess and the raw-profile bound

At a nonuniform vertex the excess is the sum of child excesses plus the
nonnegative carry count. At a repeated-child vertex it is
\(pE_0+(p-1)D_0\). A repeated vertex whose expanded subtree is not perfect
therefore has excess at least \(p-1\): if \(D_0=0\), a nonperfect child
must instead have \(E_0\ge1\), giving excess at least \(p\).
Both recursions prevent excess from decreasing at ancestors.

The converse direction in Section 2 is valid and important. A raw
skeleton with \(K=D+2\) internal vertices expands with orbit exponent
\(\delta\le K\): its perfect inserted subtrees contribute zero, and
each raw internal vertex contributes at most one. Thus its excess is at
most two, before making any assumption about visible stars. For \(p\ge5\)
this excludes every repeated-nonperfect vertex in the expansion.

A noncanonical perfect merger necessarily contains a deepest internal
vertex of the raw skeleton. Its children are raw leaves, and perfection
forces their height labels to agree. Conversely such a visible star is
plainly a perfect merger. Thus star-free raw skeletons are exactly the
canonical excess-two skeletons, and each has orbit size \(p^{D+2}\).

## 2. Carry profiles and existence of the typed counts

For a valid profile the normalization recurrence is
\[
a_i+c_{i-1}=d_i+pc_i.
\]
It gives a unique nonnegative carry vector. Summing proves
\(\sum a_i-\sum d_i=(p-1)\sum c_i\), so the excess-two canonical profiles
are precisely the stated vectors with total carry two.

If \(j\) is the largest index with \(c_j>0\), nonnegativity of
\(a_{j+1}=d_{j+1}-c_j\) forces \(d_{j+1}>0\). This justifies the finite
enumeration range below the highest digit. Adjacent uncarries are valid
even when the intermediate digit is zero.

The typed counts do extend to arbitrary nonnegative multiplicities.
For the aggregate leaf counts \(l_t\), the solution
\[
k_t=\frac{\sum_{r=0}^{e-1}p^r l_{t-r}-p^t}{p^e-1}
\]
is integral by the fixed weighted leaf sum
\(L\equiv1\pmod{p^e-1}\). Its numerator is strictly greater than
\(-(p^e-1)\), since the weighted leaf sum is positive and
\(p^t\le p^e-1\); hence it is nonnegative. Summing the flow equations
gives the asserted total internal count.

When the total count is positive, \(k_0>0\): otherwise its weighted leaf
sum would be one, forcing a single leaf and no internal vertices. The
single-leaf case has the correct type because \(p^j\equiv1\pmod{p^e-1}\)
forces \(e\mid j\). Thus neither the cycle-splicing root condition nor
zero internal types cause an exception.

## 3. Inclusion–exclusion: initially visible stars, not later mergers

The inequality \(\lfloor a_j/p\rfloor\le c_j\) is immediate from
\(a_j=d_j+pc_j-c_{j-1}\). Consequently there are at most two initially
visible stars. Distinct initially visible stars cannot overlap or nest:
all their children are leaves.

Collapsing a distinguished \(j\)-star is bijective with a raw skeleton
of profile \(b=a-pe_j+e_{j+1}\) and a distinguished leaf labeled \(j+1\).
There are \(a_{j+1}+1\) choices of that leaf. The target is deliberately
raw, so the possibility that contraction creates a new star at an
ancestor causes no overcount or missing case.

For two initial stars, simultaneous contraction yields the digit profile.
If their labels agree, the two target leaves are an unordered pair,
giving the binomial factor. If their labels differ, the two target
leaves have different labels, giving the product factor.

This remains true for adjacent labels: the two chosen target leaves are
distinct locations in the digit-profile tree. Expanding each does not
make one expansion descend from the other. If the intermediate digit
vanishes there is no target leaf for that member of a pair, so the pair
term is zero; nested mergers after a single collapse are still counted
correctly by the first-order marked term.

There is no triple term. Different carry vectors give different
canonical profiles. Contraction commutes with the group action, so each
star-free profile is a union of whole orbits of size \(p^{D+2}\).
This proves nonnegativity and integrality of each quotient in the note,
rather than merely checking those properties arithmetically.

## 4. Explicit ternary counterexample

The \(p=3\) failure has a concrete structure. Let \(A\) be a ternary tree
with two internal vertices: its root has one perfect three-leaf child
and two leaf children. Its three plane placements form one orbit of
size three, so \(\delta(A)=1\).

Let \(T\) have a root whose three children all have orbit type \(A\).
It has seven internal vertices, fifteen leaves, and
\[
D(T)=(s_3(15)-1)/2=1,\qquad
\delta(T)=3\delta(A)=3,\qquad E(T)=2.
\]

Its canonical profile is \(a_0=6,a_1=3\): the skeleton has four internal
vertices, not \(D+2=3\). The normalized digit profile is
\((0,2,1)\), and the carry vector is \((2,1)\), totaling three.
This excess-two orbit is therefore genuinely absent from a two-uncarry
enumeration. It has 27 trees and contributes nine after division by
\(3^D\), exactly the discrepancy in the source checker.

## 5. Independent tests

The new check_excess_two_independent.py imports no project helper code.
It enumerates full ordered tree shapes, assigns repeated leaf labels
subject to depth restrictions, and counts actual initially visible
stars, including their labels and unordered pairs. Its raw enumeration
does not use the factorial count to generate the trees.

All nine tested cases pass. Some useful buckets, listed as
(zero stars, one star, two stars), are:

| \(p,e\) | profile \(a\) | buckets |
|---|---|---|
| \(5,1\) | \((12,1)\) | \((125,300,30)\) |
| \(5,1\) | \((5,4,0)\) | \((625,5,0)\) |
| \(5,2\) | \((5,4,0)\) | \((0,5,0)\) |
| \(5,1\) | \((8,5,0)\), adjacent | \((42500,2525,20)\) |
| \(5,1\) | \((8,0,5,0)\), separated | \((42500,2525,20)\) |
| \(5,1\) | \((11,0,2)\), same-label pair | \((1500,1200,30)\) |
| \(5,2\) | \((5,1,5,2)\), constrained depths | \((7500,0,60)\) |
| \(7,1\) | \((7,6,0)\) | \((12005,7,0)\) |
| \(7,2\) | \((7,6,0)\) | \((0,7,0)\) |

The same independent script enumerates all 7,752 ternary trees with
seven internal vertices, identifies local cyclic orbit types recursively,
and checks each orbit's actual cardinality against its recursively
computed exponent. It finds \(M_0=1,M_1=6,M_2=15\), with exactly one
excess-two orbit having the three-carry profile above.

The original checker also passed all 1,809 normalized congruence tests,
its two raw enumerations, and its complete recursive orbit histograms.
Those arithmetic congruence tests are diagnostics, not the proof.

## 6. Is this combinatorial?

Yes, given the previously established cycle-splicing bijection and
expanded-tree group action. Factorials occur only as exact counts of
finite slot assignments and repeated labels. Orbit cardinalities,
canonical contraction, and marked-star bijections prove the divisibility
and the normalized-unit formula. No factorial-valuation identity or
Kummer theorem is used in those deductions. Using binomial counts as
independent numerical checks does not change that proof structure.
