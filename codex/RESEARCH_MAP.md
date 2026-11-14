# A focused 40-hour project on the motivation of Apéry's proof

The current scope is #131 and one selected extension. Further #137/Fuss–Catalan
work is stopped; its completed research remains in problem137/. Read the
[exact official statement](sources/problem131-official.md) and
[START_HERE.md](START_HERE.md) before choosing among the larger archive.

The central question is: **what forces the auxiliary function, the integer
coefficient, and the recurrence?** A good deliverable is a complete explanation
that the team can reproduce, plus one proved connection and one accurately
stated arithmetic gap. An all-odd irrationality breakthrough is not a realistic
requirement for completing this 40-hour exposition.

## The main proof story

Start with the positive-integer contradiction. A small rational approximation
is insufficient: both coefficients must become integers after a specified
multiplier, and the resulting nonzero linear form must tend to zero.

Then develop one kernel throughout:

\[
 K(xy)=-\frac{\log(xy)}{1-xy},\qquad
 I_n=\iint L_n(x)L_n(y)K(xy)\,dx\,dy.
\]

1. Expanding the monomial matrix makes its diagonal contain zeta3 and its
   off-diagonal entries rational. With the exact normalization
   I_n=2(A_n zeta(3)-B_n), this forces A_n=sum binomial(n,j)^2
   binomial(n+j,j)^2 and gives d_n³B_n integral, where d_n=lcm(1,...,n). See the
   [complete kernel proof](problem131/short-kernel-energy-zeta3-proof.md).
2. The kernel PDE yields (H+n(n+1))k_n=1. Explain the unit boundary flux;
   dropping that boundary term destroys the normalization. Polynomial division
   has a one-dimensional remainder, and its coefficient functional forces
   the SAME A_n. Positive energy and Delannoy growth give a small nonzero form.
3. For the official recurrence question, use the
   [two-index kernel derivation](problem131/kernel-commuting-grids-apery-recurrence.md).
   The two Legendre equations commute. Their compatibility gives the adjacent
   factor6; closing the diagonal explains17 and5. Derive the exceptional
   source at index0 explicitly. Do not replace this proof with coefficient fitting.
4. Recover Zudilin's rational square through the Mellin integral of L_n.
   Compare each ingredient with [the seven-lemma explanation](problem131/zudilin-explained.md).
   Use [Beukers' integral](problem131/beukers-explained.md) as a second view of
   positivity and reciprocal growth. A short proof can instead be read first
   in [the elementary hybrid note](problem131/short-elementary-hybrid-proof.md).

An elementary presentation can avoid recurrence asymptotics, the PNT, Bernoulli
numbers, complex contours and spectral theory in its main irrationality proof.
The recurrence still deserves its own derivation because it is the competition's
explicit motivation target. Do not prove every equivalent representation in one talk.

## Choose one extension

**Recommended accessible extension: projection and a random process.**
The [finite Christoffel argument](problem131/legendre-christoffel-cubic-tail.md)
constructs Q_k=(L_k+L_(k-1))/(1-x), with squared norm2/k. A Volterra integral
removes Q_k/k, so Pythagoras produces2/k³. The
[uniform-sampler theorem](problem131/uniform-gibbs-zeta-correlations.md) uses the
same modes: choose T uniform(x,1), then Y uniform(0,T). Two integrations give
eigenvalue1/k², and stationary logarithmic correlations give all odd zeta values.
The distinction to emphasize is that positivity and a representation do not
supply a denominator theorem. This package fits undergraduate analysis and
requires no general spectral theorem.

**Algebraic alternative: boundary repair, Schur bounds and Pell indices.**
The [design guide](problem131/schur-heun-apery-design-guide.md) begins with a
concrete two-endpoint defect. An anticommutator cancels it. Degree cutoffs
force a finite quotient; endpoint values determine a recurrence; rational
interpolation transfers it to the numerator. The resulting original upper
bound is (5B_h-B_(h-1))/(5A_h-A_(h-1)), now proved and twice audited.
At rank two the discriminant is8h²-12h+5, so a Pell equation selects rational
bounds. Choose this in place of the probability extension if the team prefers
linear algebra, orthogonal polynomials and exact computation. The general
rank-q statement is an optional appendix, not required for the presentation.

**Research alternative: one higher-zeta arithmetic construction.**
Use the [balanced Rodrigues theorem](problem131/balanced-generalized-wilson-rodrigues.md),
[grid arithmetic](problem131/balanced-andrews-grid-arithmetic.md), and one
positive remainder. State all contour/primary transformation inputs. Present
one proved normalizer and its exact analytic cost. The
[uniform single-pole limitation](problem131/balanced-diagonal-positive-rate.md)
is a useful case study: for fixed m>=2 and even n, moving the pole anywhere
in the allowed integer range a>n/2 cannot save the specified lcm multiplier,
while actual primitive cancellation remains open. Do not combine
this full package with both extensions above within40 hours.

## Time allocation and concrete output

| Hours | Work | Deliverable |
|---|---|---|
| 0–3 | Exact statement, Zudilin conventions, integer squeeze criterion | Notation sheet and proof dependency sketch |
| 3–10 | Monomial kernel arithmetic and Legendre orthogonality | Team-written proof of A_n and cubic-lcm clearing |
| 10–17 | Kernel PDE, endpoint flux, polynomial remainder and energy | Complete positivity/smallness proof with actual boundary estimates |
| 17–24 | Two commuting grids, adjacent identity and diagonal closure | Board derivation of the recurrence, including index0 |
| 24–28 | Mellin square and Beukers comparison; e^r/pi² design principle | One diagram and concise classical-source dictionary |
| 28–34 | One selected extension | One complete theorem and one clearly scoped limitation |
| 34–40 | Independent checks, attribution, manuscript and rehearsal | An8–12-page main exposition, appendices and25-minute talk |

Have a second team member check normalization, boundary signs and recurrence
closure. A numerical test is not a replacement for any of these proofs.

## Presentation and provenance

For25 minutes, spend about18 on the core proof/recurrence and7 on the chosen
connection. Start with the arithmetic target, explain the forced coefficient,
derive the recurrence, and close the contradiction. Leave the long symbolic
certificates and optional asymptotic constants in the appendix.

Use four explicit labels: classical input; independently developed and checked
here; historical novelty unresolved; proof draft/experiment. The first three
can coexist: a checked reconstruction can use classical objects while its exact
expository formulation has unresolved provenance. An unsuccessful search never
proves novelty. [CLAIM_LEDGER.md](CLAIM_LEDGER.md) and each linked audit specify
the exact scope of verification.

The [PNT-free zeta-five adaptation](problem131/pnt-free-zeta5-guide.md) is an
advanced comparison retaining Fauzan's construction, not the main40-hour project.
The [pi-power note](problem131/pi-powers.md) distinguishes the elementary pi²
argument from the advanced machinery needed in its uniform pi^r treatment.
The completed #137 notes are preserved as separate background and receive no
hours in this current plan.
