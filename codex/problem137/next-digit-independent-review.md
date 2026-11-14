# Independent review of the odd-prime second-unit-digit theorem

**Verdict: the theorem and next-orbit classification pass independent mathematical review.** No gap was found in [next-digit-orbits.md](next-digit-orbits.md). The proof refines the established expanded-tree \(p\)-group action and the checked cycle-splicing count. It does not infer the answer from factorial valuations. The odd-prime hypothesis is essential.

Let \(p\) be odd, \(e\ge1\), \(q=p^e\), \(n\ge0\), and write

\[
L=(q-1)n+1=\sum_jd_jp^j,\qquad
D=\frac{\sum_jd_j-1}{p-1},\qquad
H_d=\sum_{a=1}^d a^{-1}\in\mathbb F_p,\quad H_0=0.
\]

The reviewed theorem is

\[
\boxed{\quad
\frac{T_q(n)}{p^D}\equiv
\frac{((p-1)!)^D}{\prod_jd_j!}
\left(1-p\sum_{j\ge0}d_{j+1}H_{d_j}\right)\pmod{p^2}.
\quad}
\]

The factorial prefactor must be evaluated modulo \(p^2\), not replaced by its reduction \((-1)^D/\prod d_j!\) modulo \(p\). Only the harmonic correction needs its class modulo \(p\). Historical novelty has not been established.

## 1. Excess cannot hide a repeated nonperfect block

For a micro \(p\)-ary subtree, let \(\delta\) be its orbit exponent,
\(D_v=(s_p(L_v)-1)/(p-1)\), and \(E_v=\delta-D_v\ge0\).
When the child orbit types are not all equal, direct subtraction gives

\[
E_v=\sum_iE_i+
\frac{\sum_i s_p(L_i)-s_p(\sum_iL_i)}{p-1}.
\]

When all child orbit types agree,

\[
E_v=pE_0+(p-1)D_0.
\]

Both formulas show that excess cannot decrease toward an ancestor. A uniform nonperfect vertex costs at least \(p-1\): if \(D_0\ge1\), this is immediate; otherwise nonperfectness forces \(E_0=\delta_0\ge1\), since a size-one orbit is perfect, giving cost at least \(p\).

Because \(p-1\ge2\), an excess-one tree has no such vertex anywhere. Contracting maximal perfect subtrees leaves a skeleton in which every internal vertex contributes exactly one factor \(p\). It has \(D+1\) internal vertices. Here “nonuniform” means the child types are not all equal; pairwise distinctness is unnecessary.

## 2. One extra internal vertex forces exactly one uncarry

If \(a_j\) counts contracted leaves of height \(j\), then

\[
\sum_ja_jp^j=L,\qquad
\sum_ja_j-\sum_jd_j=p-1.
\]

Let \(c_j\) be the carries when normalizing this profile into digits. These can in general be arbitrary nonnegative integers, but the digit equations imply

\[
\sum_ja_j-\sum_jd_j=(p-1)\sum_jc_j.
\]

Thus exactly one \(c_j\) is one. The profile is necessarily

\[
a_j=d_j+p,\qquad a_{j+1}=d_{j+1}-1,\qquad
a_i=d_i\ (i\ne j,j+1),
\]

where \(d_{j+1}>0\). Different \(j\) give disjoint profiles.

Only label \(j\) now occurs at least \(p\) times, and it occurs fewer than \(2p\) times. If a vertex had \(p\) equal child orbit types, those children could therefore use only label \(j\). A non-leaf child has at least \(p\) leaves, requiring \(p^2>2p-1\) copies globally, impossible. The only repeated block is a star of \(p\) leaves labelled \(j\).

There can be at most one such star, since two would be disjoint and require \(2p\) copies. A larger accidental perfect merger contains a bottommost such star. Consequently excluding precisely these stars leaves exactly the canonical excess-one skeletons; there is no additional symmetry correction.

## 3. The star subtraction preserves the type constraints

Let \(N(a)\) count plane labelled skeletons whose leaf of label \(j\) has depth \(-j\pmod e\). The cycle-splicing counting lemma allows arbitrary nonnegative multiplicities, including \(a_j\ge p\). The digit restriction was needed for the minimal orbit size, not for that enumeration formula.

Contracting the unique forbidden star gives a digit-profile skeleton with a distinguished label-\((j+1)\) leaf. Conversely, expand any such distinguished leaf. The old condition \(r+j+1\equiv0\pmod e\) becomes the identical condition \((r+1)+j\equiv0\pmod e\). Hence the exact forbidden count is \(d_{j+1}N(d)\), and

\[
M_1=\sum_{j:d_{j+1}>0}
\frac{N(U_j)-d_{j+1}N(d)}{p^{D+1}}
\]

is the exact number of next-size orbits. Each summand is a nonnegative integer.

The internal type vector increases by one at \(t=-j-1\pmod e\). This follows directly by star expansion or from the slot equations. When \(e=1\), parent and child types coincide and the net leaf change is \(p-1\); the same internal increment is valid.

Zero internal-type counts cause no problem: the counting bijection permits them. When \(D>0\), its required type-zero internal vertex exists by the previous flow argument. For \(D=0\), the digit skeleton is a single leaf and is treated separately: every possible uncarry produces only its forbidden star, so \(M_1=0\). No nonexistent internal root name is chosen.

## 4. The modular division and its precision are correct

Put \(k=k_t(d)\), \(d=d_j\), \(c=d_{j+1}\). The exact counting ratio is

\[
\frac{N(U_j)}{N(d)}
=c\,\frac{(p(k+1))!\,k!\,d!}
          {(pk)!\,(k+1)!\,(d+p)!}
=c\,\frac{p\prod_{r=1}^{p-1}(pk+r)}
             {\prod_{r=1}^{p}(d+r)}.
\]

Cancellation of \(k+1\) occurs over the rationals first. It need not be a unit. The denominator interval contains exactly one multiple of \(p\), namely \(p\), which cancels the displayed numerator factor. All remaining factors are units.

Their products satisfy

\[
\prod_{r=1}^{p-1}(pk+r)
\equiv(p-1)!(1+pkH_{p-1})
\equiv(p-1)!\pmod{p^2},
\]

\[
\prod_{\substack{1\le r\le p\\d+r\ne p}}(d+r)
=(p-1)!\prod_{a=1}^{d}\left(1+\frac pa\right)
\equiv(p-1)!(1+pH_d)\pmod{p^2}.
\]

Only \(H_{p-1}=0\pmod p\), from pairing \(a\) with \(p-a\), is used. Inversion therefore gives

\[
\frac{N(U_j)}{N(d)}\equiv c(1-pH_d)\pmod{p^2}.
\]

Since \(N(d)=p^DM\), the exact count reads
\(pM_{1,j}=M[N(U_j)/N(d)-c]\). Reducing modulo \(p^2\), then dividing by \(p\), yields
\(M_{1,j}\equiv-McH_d\pmod p\).
The loss of one digit of precision has been accounted for.

Likewise, split the exact minimum-orbit count into blocks of \(p\)-units:

\[
M=\frac{\prod_i\prod_{\substack{1\le a\le pk_i\\p\nmid a}}a}
         {\prod_jd_j!}
\equiv\frac{((p-1)!)^D}{\prod_jd_j!}\pmod{p^2}.
\]

All larger orbits disappear modulo \(p^2\) after dividing the total by \(p^D\), leaving \(M+pM_1\). These identities prove the theorem, including the negative sign in its harmonic correction.

## 5. Boundary cases and limits

- For \(p=3\), \(H_2=3/2\) is zero modulo \(3\), but equals \(6\) modulo \(9\). The extra factor \(3\) makes its block contribution vanish modulo \(9\). No stronger harmonic congruence is used.
- For \(n=0\), \(L=1,D=0,M=1\), and no uncarry exists.
- If \(D=0\) with \(n>0\), then \(L=q^h\). Every possible uncarry has \(d_j=0\), so its correction vanishes. Its only skeleton is the forbidden star, proving \(M_1=0\), not merely \(M_1\equiv0\pmod p\). The tree count is \(1\pmod{p^2}\).
- For \(p=q=3,n=3\), there is one orbit of size \(3\) and one of size \(9\). The formula gives \(12/3\equiv1-3\cdot2=4\pmod9\).
- At \(p=q=2,n=2\), the odd-prime formula would give \(3\pmod4\), while the actual normalized Catalan count is \(1\). Both the harmonic-block argument and the exclusion of repeated nonperfect excess-one blocks need a separate binary analysis.

No formula modulo \(p^3\) follows. Already at \(p=3\), a repeated nonperfect block can occur at excess two; it is excluded only from the layer considered here.

## 6. Independent executable checks

Run from the repository root:

    python codex/problem137/check_next_digit_independent.py

The script does not import the author's orbit-ID recursion or congruence functions. It enumerates actual plane \(q\)-ary trees, expands them into micro \(p\)-trees, and canonicalizes rotations. It verifies every full orbit cardinality, contracts perfect subtrees, and compares observed excess-one profiles and their exact plane counts with the uncarry-minus-star formula.

| Check | Passed coverage |
|---|---|
| Actual plane-tree enumeration | 16,145 trees, with \((p,e,n_{\max})=(3,1,7),(3,2,4),(5,1,5),(5,2,3),(7,1,4)\) |
| Structural classification | Full orbit cardinalities, one-uncarry profiles, and exact forbidden-star subtraction in every enumerated case |
| Additional congruences | 972 cases, \(p=3,5,7,11\), \(e=1,2,3\), \(0\le n\le80\) |
| Boundary cases | 48 perfect-tree cases and 65 nontrivial zero-internal-type cases; actual-tree enumeration also includes zero-type cases |
| Exact local ratios | 824 cases including \(p\mid(k+1)\), \(d=0\), \(d=p-1\), and \(p=3\) |
| Excluded phenomena | An actual repeated-nonperfect excess-two tree is encountered, and the \(p=2\) counterexample is checked |

These tests supplement the universal proof. They do not establish it by extrapolation. The arithmetic congruence and its orbit realization may have precedents; this review establishes validity and scope, not historical priority.
