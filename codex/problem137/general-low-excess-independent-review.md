# Independent audit of the general low-excess recipe

Reviewed general-low-excess-orbits.md and check_general_low_excess.py.
**Verdict: accepted, with no correction required.** The proof gives exact
orbit counts for \(0\le k<p-1\) and hence the normalized count modulo
\(p^{p-1}\), including the \(p=2,k=0\) edge case.

The previously reviewed group action and cycle-splicing count are inputs.
The new marked-star inclusion–exclusion and sharp structural boundary
were checked independently below. Historical novelty remains unassessed.

## 1. The finite profile and selection boxes

The normalization identity
\[
a_i+c_{i-1}=d_i+pc_i
\]
uniquely recovers each nonnegative carry vector from its profile. Its
total carry is \(k\), because each carry decreases the digit sum by
\(p-1\). For the last positive \(c_j\), nonnegativity forces
\(a_{j+1}=d_{j+1}-c_j\ge0\), so \(d_{j+1}>0\).
Thus the stated range below the highest digit is exhaustive.

For \(h\ge2\), weak compositions into \(h-1\) slots give exactly the
stated upper bound \(\binom{h+k-2}{k}\). For \(h=1\), only the zero
carry vector exists. No negative-index binomial convention is needed.

The profile formula gives
\(\lfloor a_j/p\rfloor\le c_j\). For any selection vector in the
proposed box,
\[
a_i'=a_i-pb_i+b_{i-1}\ge0,\qquad
a_{j+1}'=a_{j+1}-pb_{j+1}+b_j\ge b_j.
\]
Moreover its carry vector is exactly \(c-b\), so the lower-profile
count has the same fixed weighted leaf sum and valid depth constraints.
The highest label cannot support a star since \(a_{h-1}<p\).

The number of box entries is bounded by
\[
\prod_j(\lfloor a_j/p\rfloor+1)
\le\prod_j(c_j+1)\le\prod_j2^{c_j}=2^k.
\]
The last inequality applies also to zero entries. Together with the
profile bound, this is the claimed fixed-prime polynomial-size search
in the number of digits.

## 2. The simultaneous collapse bijection

An initially visible star has only leaf children. Consequently two
distinct initially visible stars have disjoint vertex sets and neither
can contain the other. Simultaneously collapsing a selected set is
therefore unambiguous.

After collapse, the \(b_j\) selected stars of label \(j\) become an
unordered set of \(b_j\) leaf positions labeled \(j+1\). Conversely,
any such selected positions can be expanded simultaneously. Leaves of
different labels are different positions; no leaf is an ancestor of
another. This remains true for adjacent labels, even though a new leaf
created by one expansion might have the same label as another selected
leaf. The selections refer to original positions, not labels alone.

Hence the multiplicity is exactly
\(\prod_j\binom{a_{j+1}'}{b_j}\). There is no factorial for ordering
the selected stars and no additional multinomial coupling adjacent
heights.

Collapsing may create an unmarked star at an ancestor; expanding a
selected leaf may destroy such an ancestor star. This is harmless
because the target is a raw tree and the bijection remembers only the
chosen leaves. Its inverse returns precisely the original tree and
selected subset.

Depth residues are preserved directly. A star parent at depth \(d\)
whose children have label \(j\) satisfies \(d+1+j\equiv0\pmod e\).
The collapsed leaf at depth \(d\) has label \(j+1\), satisfying the same
condition. This also handles a selected root star.

For any given raw tree with \(v\) visible stars, summing signs over its
subsets gives \((1-1)^v\). Grouping these marked objects by \(b\)
therefore gives exactly the stated inclusion–exclusion, including
simultaneous selections of more than two stars.

## 3. Why the quotient counts whole orbits

For a raw profile of total carry \(k\), its expanded tree has
\(\delta\le K=D+k\). A repeated-nonperfect vertex would force excess
at least \(p-1\), so it is impossible for \(k<p-1\).
The deepest-vertex argument identifies every remaining perfect merger
with an initially visible identical-leaf star.

Thus a star-free raw tree in this range has exactly \(K\) nonuniform
internal vertices, orbit size \(p^{D+k}\), and is its own canonical
contraction. Conversely all genuine excess-\(k\) trees contract this
way. Canonical profiles are preserved by the group action and are
disjoint, proving each individual quotient is a nonnegative integer.
No factorial valuation is needed to establish that divisibility.

The raw profile count has the previously checked nonnegative integral
type solution for every weighted leaf profile. The root-leaf case,
zero internal types, and empty selection require no change.
In particular \(n=0\) gives \(M_0=1\) and no other retained profile.
For \(p=2\) the range is only \(k=0\), giving the leading normalized
digit modulo two, exactly as claimed.

## 4. Sharp structural boundary for every prime-power arity

For \(q=p\), the proposed child with two internal vertices has
\(\delta_{\rm child}=D_{\rm child}=1\). Repeating its orbit type at a
new root gives \(\delta=p,D=1,E=p-1\), while the canonical skeleton
has \(K=p+1\) and carry total \(p\). It is absent from the
carry-sum-\((p-1)\) recipe.

For \(e\ge2\), put \(r=e-2\). The child has \(p-1\) perfect branches of
height \(r\) and one of height \(r+e\). Its leaves number
\((p-1)p^r+p^{r+e}\), with digit sum \(p\), and its orbit exponent is
one because its child types are not all equal.

The repeated outer root has leaf depths exactly \(e,2e\), and hence is
an expanded \(p^e\)-ary tree. Its leaf count is
\[
(p-1)p^{e-1}+p^{2e-1}
=(p^e-1)(p^{e-1}+1)+1,
\]
giving the stated \(q\)-tree size. It again has
\(\delta=p,D=1,E=p-1\), canonical profile with
\(p(p-1)\) labels \(r\) and \(p\) labels \(r+e\), and total carry \(p\).
This checks all primes, including two.

A useful extra clarification is that raw carry-\((p-1)\) star-free
trees still cannot contain a repeated-nonperfect vertex. Such a vertex
would remove at least one local factor, giving
\(\delta\le K-1\) and \(E\le p-2\), contrary to \(E\ge p-1\).
Thus at the first boundary the old recipe counts a genuine subset:
it misses these new repeated-root orbits rather than introducing
spurious ones.

The examples establish a classification obstruction for every \(q\).
They alone need not imply that the omitted-orbit count is a unit modulo
\(p\) for every \(q\). The source checker separately verifies actual
next-digit mismatches for its stated \(q=p\) examples; no unproved
universal next-digit mismatch is needed for the theorem's sharp
structural scope.

## 5. Independent verification by a different counting method

The new check_low_excess_independent.py imports no project helper.
Its main check counts ordered trees by an ordered-child convolution,
tracking internal vertices, whether there is one exceptional leaf of
label one, and the number of initially visible stars. It does not use
the proposed inclusion–exclusion to generate those star counts.

For \(p=2,3,5,7,11\) and every \(0\le k<p-1\), it checks the raw profile
\(((p-1)(k+1),1)\), whose carry total is \(k\) and digit bound is one.
For each selection order \(b\), it compares the direct binomial moment
of the visible-star histogram against the collapse count, and then
checks full inclusion–exclusion and orbit divisibility.
All **23 profiles and 90 selection orders**, through order nine, pass.

The script also checks the explicit collapse and expansion maps on
18 selected subsets of trees with three stars, including adjacent
heights, distinct depth-constrained heights, and a chain where collapse
creates a new ancestor star. It verifies profiles, depth residues, and
inverse reconstruction directly.

Finally, it checks the repeated-block boundary symbolically at all
25 combinations \(p=2,3,5,7,11\), \(1\le e\le5\), including leaf
depths, internal counts, carry totals, and orbit exponents.

The original checker also passes all 641 normalized congruence tests,
its independently generated Sylow histograms, and its explicit
boundary mismatches. Those arithmetic checks support the proof; the
group-action and marked-object bijections remain the reason it is
combinatorial.
