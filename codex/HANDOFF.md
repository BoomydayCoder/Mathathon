# Research handoff — 14 November 2026, 10:50 UTC

Continue autonomous research in `codex/` on `codex/independent-research` in the
exact private repository `BoomydayCoder/Mathathon`. User authorizes continued
work without permission pauses, remaining credits but no automatic resets.
Reporting target:14 November2026,08:00Los Angeles/16:00UTC. No usage-balance API
is available; do not invent a balance or guarantee runtime. No Claude-owned
files have been changed. Rough effort split remains80%#131/20%#137.

## Read first

1. [Discovery guide](problem131/discovery-guide.md): complete elementary zeta-three
   proof, reconstructed motivation, Newton-grid recurrence, positive reciprocal
   Laurent model, resistor chain, integral supplement.
2. [PNT-free zeta-five guide](problem131/pnt-free-zeta5-guide.md): strongest assembled
   simplification, with complete dependency map. Computer-assisted and retaining
   Fauzan's determinant/local arithmetic; not formally verified or claimed novel.
3. [Problem137 guide](problem137/discovery-guide.md), followed by
   [next-digit theorem](problem137/next-digit-orbits.md) and its independent review.
4. [Research results](problem131/research-results.md), [claim ledger](CLAIM_LEDGER.md),
   and the timestamped [checkpoints](CHECKPOINT.md). Individual proof/review files
   take precedence over older summaries.

## Strongest checked constructions and their limits

**The Apéry recurrence is designed, not merely verified.** Prescribed zeros,
poles and decay force the rational square. Its integral residues lie on a
quadratic Newton grid. Transposed multiplication and finite diagonal balance
explain17=6*3−1 and5=6*1−1; a rational certificate transfers the recurrence to
the full linear form. A Casoratian is constant current and the approximation
error is tail resistance. This is reconstructed mathematical motivation, not
historical testimony. Use F_n=2(A_n*zeta(3)−B_n), A_0=1,A_1=5,B_0=0,B_1=6.

**One Laurent polynomial supplies both sides.** The reviewed reciprocal model
has A_n=CT Lambda^n and F_n=integral_R^3 Lambda^(-n-1). One AM–GM minimum gives
17+12sqrt2 and its inverse. The same Hessian gives A_n F_n~1/(2sqrt2*n^3).
The real generating function is relative: L F=−12t, not zero. The remaining
full two-variable exact-form certificate is still unresolved; critical values
alone do not derive a recurrence. See apery-laurent-reciprocal.md and
laurent-legendre-relative-periods.md.

**PNT-free zeta-five normalization is complete and independently audited.**
A25-term superfactorial profile majorizes every prime valuation. The exact
standard-library checker verifies561 breakpoints and an infinite-tail period
identity. Small-prime lcm repair and lower-order factorial repair handle all
higher prime powers. Patched arithmetic exponent<1.3786, while the reviewed
rational equilibrium certificate gives analytic exponent<−1.3855. No prime
number theorem is used in this alternative.

The formerly retained40000th determinant power is now ONLY an arithmetic
intermediate. Gauss's lemma: if cF^q is integral, the primitive integer part P
of F satisfies |P(theta)|<=|cF(theta)^q|^(1/q). Thus final degree remains37n
when K=40n. Explicit prime-exponent FLOOR rounding gives the same result with
no added upper cost. See rounding-prime-exponents.md and its independent review.
The earlier PNT/plain-rank route and optional outer-Hermite improvement remain
valid separate routes, with their own normalizers; do not mix their constants.

**Why superfactorials appear.** Centered residue-floor counts have exactly a
Brownian-bridge covariance. Bernoulli B2 identities turn the covariance into
six elementary floor primitives, which are precisely superfactorial valuation
profiles. The floor-covariance derivation and short two-bound refined-rank lemma
are independently checked. The optional outer correction's full selected kernel
is an explicit Hermite ideal of dimensionK+3N−3p+1; this yields a proved
27/3200 arithmetic-budget saving, not an exact determinant valuation theorem.

**Higher odd constructions.** The monomial-kernel theorem isolates a single
zeta(s) from two integer polynomials. Multiplicity-lifted Rodrigues polynomials
repair unwanted early derivative contributions; for zeta5 a reviewed
exponentially small bound and actual new-prime denominator obstruction are proved.
The related positive higher-dimensional periods mix lower zeta values. Adaptive
adjacent-Legendre cancellation is now independently reviewed and eventually
positive, but available arithmetic remains too expensive; finite samples are
not a global failure theorem. Bernoulli period towers and Euler projectors give
positive Bose tails but lose rational-parameter arithmetic.

**Signed zeta-seven model: reviewed.** An ordinary signed atomic measure isolates
zeta7, with explicit lcm denominator control and a verified nonzero-Hankel input
from the primary paper. A proved lower bound shows that discarding all signs
and bounding by total variation cannot give the required small determinant
for this construction. It does not rule out exploiting signed cancellation.
See signed-single-zeta7-model.md and signed-zeta7-independent-review.md.

**Problem137: complete combinatorial exactness plus next digit.** Expanding
p^e-ary vertices into perfect p-ary blocks and counting minimum Sylow orbits
proves v_p(C)=(s_p(L)−1)/(p−1) and the leading unit. An explicit cycle-splicing
bijection replaces Lagrange inversion. The reviewed next-digit theorem for oddp
classifies one-uncarry skeletons and gives C/p^D modulo p^2; no p=2 extension
or modulo-p^3 theorem is claimed. Primary comparison credits Deutsch–Sagan and
Konvalinka; historical novelty remains unresolved.

## Actual formalization verification

The pinned existing Romik repository fc2fefa77076e22e7e80b4c101ded5a98d53e14d
was built locally:8995jobs completed; repository and independent ordinary-series
target checks succeeded. Printed axioms are propext,Classical.choice,Quot.sound,
and the explicit PNT axiom chebyshev_theta_asymptotic. A separate project kernel
replay of70source modules plus the root module completed exit0 at10:28:40UTC,
692.237seconds,maxRSS8294040KiB. This rechecks project declarations against
imports, not a fresh replay of all Mathlib or an independent kernel implementation.
It verifies that existing formalization with its stated axiom, NOT our PNT-free
variant. See lean-local-verification.md and verification/romik-lean/ logs.
No build/replay is currently running; toolchains/caches are ignored in codex/.cache.

## Active bounded work

- Connections independently audits the sharp oscillatory zeta5 Mellin–Barnes
  asymptotic in higher-rodrigues-laurent-geometry.md. Until review completion,
  its stronger asymptotic and all-index sign claims remain proof drafts.
- Zudilin-exposition develops the all-m Eulerian Barnes decomposition and bounded
  m3,4,5 saddle certification. Integral independently audits the all-m theorem.
  A saddle leaves the justified strip at m6: no arbitrary-m decay assertion.
- Odd-zeta independently audits the positive Euler-depleted moment construction.
  It removes odd-prime polynomial denominators, but pole divided differences
  explicitly regain expensive arithmetic. No asymptotic impossibility claim.
- Simplify audits root's integer-depletion-cumulants.md: any ell>=2, positive
  frequency deletion, discrete-uniform cumulants, proposed exact uniform
  denominator product_(p|ell)p^(1+v_p(s−1)). Parent exact checks pass;
  independent review is still pending at this timestamp.
- Manuscript updates proof/status indexes. Root maintains this handoff and
  checkpoints, saves independent derivations and commits/pushes codex-only work.

## Sources and reproduction

Exact official131/137 statements are saved with raw JSON and hashes under sources/.
The original Zudilin paper and primary references were used before AI drafts.
A single read-only Claude SUMMARY consultation at commitfc6da97 is logged in
sources/claude-read-log.md; no HANDOFF existed there. It supplied an outer-prime
lead, whose proof was independently developed and audited. No Claude files changed.

Standard-library scripts verify the scalar energy, endpoint and Barnes certificates;
other symbolic checks use SymPy/mpmath. Commands and dependencies are in README.md
and script headers. Finite computations supplement, not replace, proofs.
No all-odd irrationality theorem, historical novelty, or guaranteed runtime follows.
