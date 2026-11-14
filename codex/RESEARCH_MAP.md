# A 40-hour project: why Apéry's proof fits together

Build the project around one question: **what forces the auxiliary square,
the recurrence, and the reciprocal growth and decay?** A strong outcome is a
complete proof of \(\zeta(3)\)'s irrationality that the team can reconstruct,
plus one well-chosen extension and a combinatorial proof for prime-power
Fuss–Catalan numbers. The goal is an understandable chain of reasons, with
precise boundaries on generalization.

Allocate **32 focused hours to #131 and 8 to #137**. These are suggested project
effort blocks, not a prediction about autonomous runtime. Team members can work
in parallel, but should preserve the approximate 80/20 split. The exact
[#131 statement](sources/problem131-official.md) asks for recurrence motivation;
the [#137 statement](sources/problem137-official.md) asks for combinatorial
exactness of the valuation, beyond divisibility.

## The main proof story

**Shortest first reading:** the [elementary hybrid proof](problem131/short-elementary-hybrid-proof.md)
puts the rational arithmetic and positive integral together and finishes with
AM–GM and an elementary lcm bound. It gives a complete proof without using the
recurrence, so the team can first see exactly what must be achieved. Then use
the recurrence narrative below to address the official motivation question.

Use the [Apéry discovery guide](problem131/discovery-guide.md) as the spine.
Return to [Zudilin's seven lemmas](problem131/zudilin-explained.md) when checking
that the reconstructed story really proves the original assertions.

1. **Design the contradiction.** We need integral coefficients, a controlled
   rational denominator, and a positive nonzero error smaller than that
   denominator scale. Explain these three obligations before introducing a
   complicated expression.
2. **Force the square.** Double zeros at \(1,\ldots,n\), double poles at
   \(0,-1,\ldots,-n\), and decay \(O(t^{-2})\) exhaust the numerator degree:
   \[
   R_n(t)=
   \left[\frac{(t-1)\cdots(t-n)}{t(t+1)\cdots(t+n)}\right]^2.
   \]
   Partial fractions turn its integer residues into
   \(A_n=\sum_k\binom nk^2\binom{n+k}k^2\).
   This is uniqueness within the stated constraints, not a claim that every
   irrationality proof must use this square.
3. **Discover the recurrence on a quadratic grid.** The factors
   \(\phi_k(n)=\binom nk\binom{n+k}k\) are Newton polynomials in \(n(n+1)\).
   Transposing multiplication gives a nearest-neighbor operator. A finite
   balance identity closes its diagonal kernel:
   \[
   A_n+A_{n-1}=6\sum_k\phi_k(n)\phi_k(n-1).
   \]
   This is where \(17=6\cdot3-1\) and \(5=6\cdot1-1\) come from. The
   [two reviewed mixed telescopes](problem131/apery-mixed-numerator-independent-review.md)
   transfer the same balance to the full linear form: separate the square
   into \(P_nP_m\), prove its two-index recurrence and adjacent identity,
   then close the diagonal. Keep the nonzero source at index zero.
   The long diagonal certificate can stay in the appendix.
4. **Explain the reciprocal scales with one object.** For
   \(C=(1+x)(1+y)\), the positive Laurent polynomial
   \[
   \Lambda=(1+z)C\left(1+\frac{C}{zxy}\right)
   \]
   produces \(A_n\) as the constant term of \(\Lambda^n\) and the Beukers
   error as an integral of \(\Lambda^{-n-1}\). Its elementary minimum is
   \(17+12\sqrt2\). The same object explains growth and decay. The minimum
   alone is not a derivation of the recurrence; the Newton step supplies that.
   The [reciprocal-region bridge](problem131/reciprocal-region-beukers.md)
   makes this literal: a rational involution sends the positive octant to a
   bounded region, then a rational cube parametrization recovers Beukers'
   integral with its correct Jacobian. Draw the domains rather than reproduce
   the long geometric divergence certificate.
5. **Make nonvanishing transparent.** In the normalization
   \(A_0=1,A_1=5,B_0=0,B_1=6\), the Casoratian gives
   \[
   A_nB_{n+1}-A_{n+1}B_n=\frac6{(n+1)^3}.
   \]
   After separately identifying \(B_n/A_n\to\zeta(3)\), this yields
   \[
   A_n\zeta(3)-B_n
   =6A_n\sum_{j=n}^{\infty}
   \frac1{(j+1)^3A_jA_{j+1}}>0.
   \]
   A resistor-chain picture makes this a positive tail resistance. Finish
   with the cubic-lcm denominator estimate and the positive-integer
   contradiction.

The elementary route can be presented without a saddle-point theorem. The
integral/Laurent representation is an explanatory supplement and an independent
view of positivity. Do not prove every equivalent representation in the talk.

Attribute the Laurent model precisely: it is \(\phi(-x,-y,-z)\) in
Golyshev–Kerr–Sasaki, *Apéry extensions*, §5.4, which also gives the Beukers
relative period and its inhomogeneous differential equation. The
[Fano dictionary](problem131/fano-v12-apery-gamma-dictionary.md) records the
exact published match and a correction to a printed coordinate map.
The local elementary reconstructions do not make that geometric core new.

## Choose a small extension package

Choose **one** of the next two packages. The first emphasizes approximation
and positivity; the second emphasizes integer arithmetic. Both have complete
reviewed inputs, while neither proves irrationality of all odd zeta values.

**Option A: all-odd rational bounds from a positive minimum.**
The [rank-one construction and review](problem131/rank-one-zeta-independent-review.md)
turn a Hankel matrix with only one zeta-dependent entry into
\[
 E_{h,s}=(-1)^{(s-1)/2}(\zeta(s)-q_{h,s})
 =\min_{\substack{p\text{ monic}\\\deg p=h-1}}
 \int_0^\infty\frac{p(y^2)^2w_s(y)}
 {D_{2h+(s-3)/2}(y^2)}\,dy>0.
\]
Here \(D_K(t)=\prod_{j=1}^K(t+j^2)\), and \(w_s\) is the explicit positive
odd-zeta density in the linked proof.
The rational bounds are strictly monotone. Comparing the positive measures
shows that every fixed odd \(s\) shares the same lower and upper exponential
rates as \(s=3\), without assuming the rate has a limit. A
[reviewed repeated-root trial](problem131/rank-one-elementary-trial-independent-review.md)
proves \(\limsup E_{h,s}^{1/h}<25/12544<1/500\).

**Deliverable:** derive the constrained minimum and monotonicity on two pages,
then choose either the measure comparison or the one-variable trial potential
for a worked extension. Report denominator control as the missing arithmetic
step. A different, moving-pole Wilson family is already proved to equal
adjacent-Apéry quotients and inherits their arithmetic; identifying the
original \(D_{2h}\) determinant roots with that family remains conjectural.
Do not use the finite agreement to import an exact rate.

**Option B: cyclotomic moment integrality.**
The [cyclotomic depletion connection](problem131/cyclotomic-depletion-connection.md)
has a written proof and independent exact review/checks recorded in the linked
integer-depletion audit. Starting from the positive odd-zeta density \(w_s\),
the combination
\[
6w_s(y/6)-2w_s(y/2)-3w_s(y/3)+w_s(y)
\]
is still positive: inclusion–exclusion retains frequencies divisible by
neither \(2\) nor \(3\). Every polynomial moment is an integer. The unexpected
bridge is
\[
\frac{(q^6-1)(q-1)}{(q^2-1)(q^3-1)}
=q^2-q+1=\Phi_6(q).
\]
Its value at \(1\) is \(1\), so an integer-coefficient recurrence makes the
logarithmic derivatives of \(\Phi_6(e^z)\) integral. Those derivatives are
exactly the Bernoulli factors in the moments. For \(s=3\), the first three
moments are \(6,100,3822\).

**Deliverable:** a two-page proof, a one-page exact computation, and a paragraph
explaining the price: the real decay becomes slower, and integral polynomial
moments do not settle pole denominators or determinant normalization. The
broader Möbius/cyclotomic classification is a separate developing extension;
its review status should be checked before adding it to the submission.

**A common structural connection, if time permits.** The
[reviewed Euler-operator hierarchy](problem131/odd-hierarchy-independent-review.md)
expresses every \(\mu_s\) through the base functional \(\mu_3\):
\[
 \mu_s(g)=\frac1{(s-1)(s-2)/2}\,
 \mu_3\!\left[\binom{2\Theta+s}{s-3}g\right],
 \qquad \Theta=t\,\frac d{dt}.
\]
The cleared operator preserves integer polynomial coefficients, but need not
preserve pointwise positivity. This is a useful one-slide connection between
analysis and arithmetic. The subsequent hyperbolic multiple-orthogonality
construction now has a complete independently reviewed theorem for every
balanced index, with simple positive roots; its factorial constant-term
formula remains a finite-data conjecture. Neither normality nor the operator
identity supplies approximation denominators.

**Keep one honest barrier beside the success.**
The [reviewed monic-normalization obstruction](problem131/original-wide-monic-independent-review.md)
concerns the specified original, undepleted determinant with \(K=40n\),
\(N=3n\), numerator power six, and fixed odd \(s\ge13\). Its actual primitive
integer polynomial eventually has value larger than \(e^{K^2/20}\) in absolute
value at \(\zeta(s)\). Necessary coefficient denominators, a real lower bound,
and the PNT prove this; it is not merely an expensive sufficient lcm estimate.

**Deliverable:** one page separating “this actual normalization fails” from
“our available bound is insufficient.” State the family and \(s\ge13\)
restriction prominently. This does not decide \(s=7,9,11\), deformed
determinants, or irrationality of any higher zeta value. Its full technical
proof can remain an appendix; the logical distinction belongs in the main story.

**Analytic alternative if the team prefers contours.**
Use the [higher-period source comparison](problem131/rodrigues-primary-source-comparison.md)
and the [reviewed Eulerian period decomposition](problem131/eulerian-periods-independent-review.md).
The Rodrigues family, pure derivative forms, Barnes decomposition, and positive
central period are exact specializations of Nesterenko's 2003 results; the
polynomial is classical confluent Jacobi–Piñeiro. The independently checked
development makes the Laurent geometry and the contributing integration cycles
explicit. The central component is positive, but the full pure form can be
dominated by a complex conjugate pair.

Choose this **instead of either extension package**, not as a second
full analytic project. A good deliverable is the \(m=2\) diagram and a worked
parameter dictionary into the primary theorem. The finite \(m=2,3,4,5\) saddle
certificates have separate completed reviews, but reproducing all of them is
unnecessary for a 40-hour exposition. Small forms still need adequate arithmetic.

## The eight-hour combinatorial strand

Follow the [#137 discovery guide](problem137/discovery-guide.md), then use the
[reviewed low-excess recipe](problem137/general-low-excess-independent-review.md)
for one controlled extension.

For \(q=p^e\), expand each \(q\)-ary vertex into a perfect \(p\)-ary block.
Write \(T_q(n)=\binom{qn}{n}/((q-1)n+1)\). The full iterated \(p\)-group
records local rotations. If
\[
L=(q-1)n+1=\sum_i d_ip^i,\qquad
D=\frac{\sum_i d_i-1}{p-1},
\]
the smallest orbit size is \(p^D\). Contract maximal perfect subtrees, retain
their heights as leaf labels, and count the resulting depth-constrained
skeletons by the cycle-splicing bijection. Counting the smallest orbits modulo
\(p\) proves exactness, not just divisibility:
\[
v_p(T_q(n))=D,\qquad
p^{-D}T_q(n)\equiv\frac{(-1)^D}{\prod_i d_i!}\pmod p.
\]

For the extension, introduce **excess \(k\)** as the extra exponent in an
orbit of size \(p^{D+k}\). Below \(k=p-1\), repeated nonperfect child blocks
cannot occur. Carry profiles and inclusion–exclusion over visible identical-leaf
stars count these orbits exactly. At \(k=p-1\), a new mechanism really appears:
this is a structural boundary, not a failed computation.

**Deliverable:** two drawings of expansion/contraction, the forward and inverse
cycle-splicing maps on a small example, a proof of the leading-unit formula,
and one worked excess-one example for an odd prime. A small exact enumeration
should check the example. The general low-excess theorem can be stated with
its reviewed hypotheses; deriving every digit formula is optional.

Credit Konvalinka's arithmetic exactness and divisibility argument and
Deutsch–Sagan's binary orbit method. The current proof addresses the exactness
gap in the official prompt; historical novelty of its refinement is unresolved.

## A workable allocation

| Hours | Focus | Concrete result |
|---|---|---|
| 0–2 | Read the exact prompts; fix conventions and the irrationality criterion | One notation sheet and a proof dependency sketch |
| 2–9 | Auxiliary square, residues, denominator control | Team-written proof of the arithmetic lemmas |
| 9–16 | Newton grid, finite balance, and two mixed telescopes | Derive the recurrence and its source on a board without guessing |
| 16–22 | Laurent/integral bridge and Casoratian | One diagram connecting integers, positive error, and reciprocal scales |
| 22–27 | Choose rank-one variational bounds or cyclotomic integrality; contours are an alternative | One extension theorem and one precisely scoped limitation |
| 27–32 | Assemble and cross-check #131 | An 8–12-page main exposition and a rehearsed proof story |
| 32–40 | #137 action, bijection, exactness, one excess example | A 4–6-page combinatorial supplement with drawings and an exact check |

Have a second team member check at least the recurrence closure, the
normalization of the error, and the inverse bijection. These are more useful
than many extra numerical samples. A proof step should be marked “gap” if
the team cannot supply its justification, even when a supporting note reports
successful computational checks.

## The talk and the status labels

For a 25-minute talk, spend about 20 minutes on #131 and 5 on #137. Lead with
the positive-integer contradiction; explain the forced square; derive the
finite balance; show the reciprocal Laurent picture; close the proof. Use
one variational, cyclotomic, or complex-cycle slide to show what transfers and
what breaks. Finish with the tree action turning a base-\(p\) carry into an
orbit-size contribution.

Use four labels in the manuscript:

- **Classical input:** cite the primary theorem and the exact specialization.
  Apéry, Zudilin, Beukers, Golyshev–Kerr–Sasaki, Nesterenko, Jacobi–Piñeiro, and the earlier tree
  results belong here where applicable.
- **Independently developed and checked here:** a complete written derivation
  with its review scope. This label records verification, not historical priority.
- **Historical novelty unresolved:** the particular identity or formulation has
  not been attributed; an unsuccessful search is not evidence that it is new.
- **Proof draft or experiment:** an argument awaiting review, a conditional
  contour theorem whose hypotheses remain unchecked, or numerical evidence.
  None should silently become a theorem in the talk.

The [PNT-free zeta-five guide](problem131/pnt-free-zeta5-guide.md) is an advanced
optional appendix. It is a reviewed computer-assisted alternative derivation
retaining Fauzan's determinant and local arithmetic, with explicit finite
certificates. It is not an elementary proof, an independently invented
historically new construction, or a formalization of our changes. Reproducing
that entire route would crowd out the official #131 motivation goal. Keep it
as a research comparison unless the team deliberately reallocates the project.

The project succeeds if the team can explain and prove the selected results,
name the inputs honestly, and identify the next genuine lemma. Its success
does not depend on obtaining an all-odd irrationality theorem.
