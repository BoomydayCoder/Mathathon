# Read this first: the checked research and the remaining gaps

Work is in the private **BoomydayCoder/Mathathon** repository, branch
**codex/independent-research**, entirely under codex/. Claude-owned files were
not edited. The live state and last save are in [HANDOFF.md](HANDOFF.md) and
[CHECKPOINT.md](CHECKPOINT.md); [REPRODUCE.md](REPRODUCE.md) gives a short set of exact checks. This guide prioritizes the current results over
older exploratory notes. “Checked” means a written proof and a separate
mathematical audit; it does not mean historical novelty or formal verification.

## The main result for the competition

The [exact official #131 statement](sources/problem131-official.md) asks for a
satisfying explanation of Apéry's proof and especially its recurrence. There
are now two complementary answers with complete checked arguments:

- [Zudilin explained](problem131/zudilin-explained.md) follows the primary paper's
  seven lemmas, their purposes and arithmetic. The [discovery guide](problem131/discovery-guide.md)
  shows how rational interpolation and a quadratic grid suggest its construction.
- The [one-kernel energy proof](problem131/short-kernel-energy-zeta3-proof.md)
  gives a complete elementary irrationality proof. Expanding one logarithmic
  kernel gives the lcm-cube arithmetic. Its differential equation and polynomial
  division force the Apéry binomial sum. Positive energy proves nonvanishing
  and smallness. The [same kernel derives the recurrence](problem131/kernel-commuting-grids-apery-recurrence.md):
  compatibility of two Legendre index equations forces the factor6, then
  17=6*3-1 and5=6*1-1. The index-zero source is retained.

For the integral supplement, use [Beukers explained](problem131/beukers-explained.md)
and the [exact kernel/resolvent bridge](problem131/log-kernel-legendre-resolvent-bridge.md).
For the shortest proof before learning the recurrence, the
[elementary hybrid proof](problem131/short-elementary-hybrid-proof.md) combines
rational arithmetic with a positive integral and an elementary lcm estimate.
These reconstruct classical objects; no priority claim is attached to the
shorter presentation.

## The strongest additional connections

**Boundary repair, Apéry and Pell equations.** The original positive Schur
construction is now proved to give

\[
 q_{h,3}=\frac{5B_h-B_{h-1}}{5A_h-A_{h-1}},\qquad h\ge1.
\]

Two independent audits accepted the all-index proof. Its formerly necessary
contour boundary step now has a finite-difference replacement. The
[design guide](problem131/schur-heun-apery-design-guide.md) explains how an
explicit boundary defect forces a symmetric difference operator and why its
degree cutoffs determine the recurrence. At rank two, the same mechanism
produces the field Q(sqrt(8h²-12h+5)); Pell indices give rational lower and
upper bounds. The [precise theorem](problem131/zeta3-rank-q-heun-reduction.md)
includes explicit Apéry formulas and inherited lcm-cube arithmetic.

The [sharp error theorem](problem131/rank-one-original-error-tail-and-rate.md)
identifies the optimized zeta3 error, with leading constant, and an exact
rational bracket width. The existing measure comparison then proves the same
exponential rate (17+12sqrt2)^(-2) for every **fixed** odd-zeta rank-one family.
It does not supply their missing denominator estimates.

**Orthogonal projection, probability and all odd zeta values.** The
[finite Christoffel argument](problem131/legendre-christoffel-cubic-tail.md)
explains the cubic tail by Pythagoras: the kth mode has coefficient1/k and
squared norm2/k. The [uniform-sampler theorem](problem131/uniform-gibbs-zeta-correlations.md)
then uses an elementary random process: choose T uniformly between x and1,
then choose the next point uniformly between0 andT. In stationarity, where
X_0 has density2(1-x), and with K(x)=-log(x)/(1-x),

\[
 \mathbb E[K(X_0)K(X_r)]=4\zeta(2r+3).
\]

The Euler deformation kernel is the resolvent of this process. The general
Gibbs/orthogonal-polynomial framework is classical and explicitly credited;
historical provenance of the particular formulation remains unresolved.
This is a positive representation, not an irrationality theorem.

**One proof-design principle for e^r, pi² and zeta3.** Read
[polynomial inverses and boundary pairings](problem131/polynomial-inverses-and-boundary-pairings.md).
For e^r the exponent is a nonzero rational number. The note explains why
classical finite derivative sums are inverses of nilpotent polynomial operators,
and why Apéry's degree-raising operator instead has a one-dimensional obstruction. The [pi-power note](problem131/pi-powers.md)
contains a direct elementary pi² proof and the explicit advanced route for
all nonzero rational powers pi^r. Arbitrary real powers are not asserted irrational.

## What was and was not achieved for higher odd values

The [balanced Rodrigues construction](problem131/balanced-generalized-wilson-rodrigues.md)
now has proved all-index polynomials, constant terms and endpoint values. Its
raising mechanism has an [explicit primary-source antecedent](problem131/balanced-rodrigues-primary-source-comparison.md);
the balanced orthogonality and arithmetic require their separate proofs.
Separate reviewed arguments prove binomial coefficient clearing, a positive
Andrews multisum, integer-grid divisibility, and associated rational numerator
normalizers. These are actual arithmetic theorems, not numerical patterns.
The [cubic specialization](problem131/cubic-balanced-adjacent-apery-identity.md)
is exactly an adjacent Apéry remainder, including its numerator. Its
[real kernel model](problem131/cubic-wilson-real-kernel-model.md) now constructs
the cubic Wilson polynomials without their contour/Rodrigues input and proves
an exact Dirichlet inner-product identity.

For each fixed odd value s=2m+1 with m>=2, [the positive single-pole analysis](problem131/balanced-diagonal-positive-rate.md)
proves an exact decay rate along even polynomial indices n and a uniform
limitation: the specified multiplier d_(max(n,a))^s prevents decay for
**every** allowed integer pole a>n/2. Here d_N=lcm(1,...,N); this statement does
not assume that the specified multiplier clears the associated numerator.
A [consecutive-pole construction](problem131/balanced-consecutive-pole-transport.md)
gives a second positive pure-zeta approximation and explicit arithmetic.
Neither family yet has an adequate primitive denominator theorem. Failure of
a supplied normalizer is not a proof that all smaller normalizers fail.

The [PNT-free zeta-five guide](problem131/pnt-free-zeta5-guide.md) is an advanced,
computer-assisted alternative retaining Fauzan's determinant and local
arithmetic construction. Its certificates and dependencies are documented.
It is not a new independent elementary proof. These reviews are internal to
this project; external acceptance of the preprint is not claimed. The pinned
existing Lean project was built and replayed separately, with an explicit PNT axiom; that verification
does not formalize our PNT-free changes. No proof that all odd zeta values are
irrational was obtained.

## A realistic next step

Use the [updated 40-hour map](RESEARCH_MAP.md) to turn a small selection into a
coherent manuscript. Keep the elementary kernel/recurrence story central,
then choose either the uniform-sampler connection or the Schur/Pell extension.
The main open research task is sharper **primitive arithmetic**, with actual
proofs of primewise cancellation. Novelty checking requires an exact primary
comparison, not a search that happened to find nothing. Further #137 work was
stopped at the user's request; its completed notes remain preserved.
