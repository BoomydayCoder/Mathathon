# Research handoff

Continue independent research in `codex/` on `codex/independent-research`, with
approximately 80% effort on #131 and 20% on #137. The user authorizes continued
work without permission pauses, remaining credits but no automatic resets.
Reporting target: 14 November 2026, 08:00 Los Angeles / 16:00 UTC. No usage-balance
API is available; do not invent a remaining balance or guarantee runtime.
No Claude-owned files have been changed.

## Read first

- [Discovery guide](problem131/discovery-guide.md): complete elementary zeta-three
  proof, Newton-grid recurrence motivation, resistor chain, integral supplement.
- [Research results](problem131/research-results.md): constructions and limits;
  the newer proof/review links below take precedence where summaries lag.
- [Zeta-five guide](problem131/fauzan-equilibrium-guide.md),
  [assembly audit](problem131/fauzan-simplified-proof-audit.md), and
  [rational energy review](problem131/rational-energy-review.md).
- [Problem #137 guide](problem137/discovery-guide.md): complete combinatorial
  prime-power valuation proof, including counting bijection and inverse.
- [Claim ledger](CLAIM_LEDGER.md): its current-status addendum supersedes stale
  review/access entries. The parent maintains [CHECKPOINT.md](CHECKPOINT.md).

## Strongest checked results

**Recurrence motivation.** The square is forced by stated zero/pole/degree
constraints. Transposing multiplication on the quadratic Newton grid gives the
cubic operator; finite balance closes the adjacent diagonal, explaining
`17 = 6·3 − 1` and `5 = 6·1 − 1`. The rational certificate transfers the
recurrence to the full linear form. The Casoratian becomes constant current,
and approximation error becomes tail resistance. Identification with zeta three
and arithmetic normalization remain separate, proved steps.

**Shared coefficients and balance.**
[The reviewed endpoint-polynomial construction](problem131/shared-coefficients-exp-pi-zeta.md)
connects exponential and pi-squared proofs to Legendre coefficients and Apéry
residues. [Tangent](problem131/shared-coefficients-tangent.md) uses a Casoratian
where positivity fails. The elementary pi-squared argument alone does not
prove every nonzero rational power of pi irrational.
[The balance review](problem131/zeta-balance-independent-review.md) proves the
zeta-two recurrence with coefficients 11 and 3, the three low-degree cases
within the specified ansatz, and uniform mode/growth formulas for every fixed
positive pair of binomial exponents.

**Monomial kernel: independent review complete.**
[The new pairing theorem](problem131/monomial-kernel-design.md) isolates zeta(s)
from any two integer polynomials: diagonal monomials supply its coefficient;
off-diagonal terms are rational with controlled lcm denominators. Binomial and
Legendre partner choices recover the zeta-two and zeta-three coefficients.
[The dedicated review](problem131/monomial-kernel-review.md) checks every formula,
including finite-partial-sum cancellation when s=2. This supplies arithmetic,
not arbitrary-partner positivity or smallness.

**Zeta-five alternative and stronger rational certificate.**
The completed assembly retains Fauzan's determinant/local arithmetic engine,
replaces sixteen trial measures, and simplifies the arithmetic integration and
outer rank bound. It is a reviewed computer-assisted alternative derivation,
not a new independently invented construction, elementary proof, formal
verification, or external consensus.
The plain-rank option requires the **new** normalizer `m_tilde`, with outer
exponent `v_p(S_K) + 2 sum(weights) − rank`; do not retain the source's
stronger `min(rank, zero-count)` normalization.

The arithmetic budget is `2203747/1600000). The independently reviewed
[standard-library certificate](problem131/rational-energy-certificate.md)
proves energy `< -2771/2000`, leaving exact margin `13053/1600000`.
For `K = 40n`, the same degree-37n integer polynomials satisfy
`0 < Q_tilde_40n(zeta(5)) < exp(-13 n²)` eventually.
No explicit onset or irrationality measure is supplied. This removes Arb from
the scalar certificate, not the shared endpoint/potential inputs or the
advanced mathematics of the full proof.

**Elementary digit-character rule: review complete.**
For primes p >= 5, the digit polynomial is Delta^epsilon Q², with
`epsilon = (1 − (-6|p))/2` and `Delta = 1 − 34t + t²`.
[The independent review](problem131/apery-cm-independent-review.md) verifies the
ODE involution proof, including characteristic-p uniqueness and small primes.
Complex multiplication is optional interpretation with stated standard inputs.
The modular irrationality route and Wilson one-pole denominator obstruction
also have completed reviews.

**Problem #137: exactness by minimum orbits.**
Expand p^e-ary vertices to perfect p-ary blocks. Minimum orbit exponent is
`D = (s_p(L) − 1)/(p − 1)`, where `L = (p^e − 1)n + 1`.
[Cycle splicing](problem137/cyclic-type-cycle-splicing.md) counts depth-constrained
skeletons without Lagrange inversion or Matrix–Tree. The minimum-orbit count is
`(-1)^D / product(d_j!)` modulo p, proving exactness and the leading unit.
[The primary-source comparison](problem137/konvalinka-comparison.md) credits
Deutsch–Sagan's binary method and distinguishes Konvalinka's arithmetic
exactness and group-action divisibility. Novelty and weighted extensions remain
unresolved.

## Active bounded tasks

1. **Odd-zeta agent:** derive actual leading local-block nullvectors and test
   determinant valuation gains, rather than only improving sufficient entry bounds.
2. **Connections agent:** investigate half-reflected positive sums for W_n;
   distinguish signed identities, positivity, and nonvanishing modulo varying primes.
3. **Integral agent:** the monomial-kernel audit is complete. Constructive
   periodic-weight and mixed-partner extensions are being developed separately;
   their proof/review status is not inherited from the original theorem.

The common analytic budget for each fixed odd index is proved; matching
arithmetic normalizers are not. No all-odd irrationality theorem follows.
Specific-family obstructions do not exclude other families or cancellations.

## Sources, conventions, reproduction

Exact official [#131](sources/problem131-official.md) and
[#137](sources/problem137-official.md) statements are saved with raw JSON and
hashes; earlier HTML/access-blocked notes are superseded. Primary papers are
saved in `sources/`; [provenance](sources/odd-zeta-provenance.md) distinguishes
preprint claims and inspected formalization reports. Lean has not been built.

Use `F_n = 2(A_n zeta(3) − B_n)`, with A_0=1, A_1=5, B_0=0, B_1=6.
The early auxiliary-rational-function note uses half this normalization.
The recurrence's last term is indexed n−1. Reconstructed motivation is not a
historical account; no historical novelty is asserted.

Endpoint and rational-energy scripts use ordinary Python's standard library;
commands are in [README.md](README.md). Other checks use SymPy/mpmath.
Older Carlson/virial checks use python-flint/Arb and remain cross-checks.
Finite computations supplement, not replace, the general proofs.

## Newest continuation,09:45UTC

The common Laurent model is now independently reviewed and integrated into the
main discovery guide. A_n=CT Lambda^n, F_n=integral_R^3 Lambda^(-n-1), so one
AM–GM minimum gives both reciprocal exponential rates. Its same Hessian also
gives A_n F_n~1/(2sqrt(2)n^3). The real generating series has the boundary term
L F=-12t; a homogeneous all-period assertion would be incorrect.

The Bernoulli period tower and its projection extension are checked. They connect
higher odd constants, rational Taylor data and a positive Bose-moment tail, while
pinpointing the loss of rational-q arithmetic under logarithmic projection.
The multiplicity-lifted single-ζ5 forms have a reviewed exponential bound and
actual new-prime obstruction; no irrationality conclusion follows.

Active next tasks: stronger ζ5 outer-prime cancellation with an adversarial
review; exact higher-dimensional Laurent period mixing; multiple orthogonality
and root geometry; bounded #137 next-prime-digit refinement. A local Lean build
of the existing pinned Romik formalization is running under ignored codex/.cache;
verification logs/metadata are under codex/verification. No success claim yet.

Read-only Claude consultation is now documented in sources/claude-read-log.md.
Its published SUMMARY at fc6da97 supplied only the outer-prime investigation lead;
we are deriving and checking the mathematics independently. No Claude files changed.

09:56UTC update: outer-correction-adversarial-review.md accepts the ζ5 outer
lemma. `outer-gain-budget-update.md` defines the new normalizer and derives
budget2190247/1600000 and eventual exp(-26n²). Existing simpler route stays valid.
The pinned Romik Lean project and two final checks have now succeeded locally;
standard logic plus explicit PNT are the printed axioms. Separate kernel replay
is still active. See verification/romik-lean/ logs, not just source metadata.
The positive-period amplitude-mixing theorem has an independent review. The
#137 odd-prime modulo-p² refinement awaits manuscript-agent review.

## Read first after the tenth checkpoint

`problem131/pnt-free-zeta5-guide.md` is now the strongest assembled simplification:
25 explicit Barnes/superfactorial factors replace PNT averaging, with exact
standard-library verification and a complete independent normalizer audit.
The underlying proof remains computer-assisted and retains the source local
lemmas; it has not been formalized. The prior plain-rank and new outer-kernel
routes remain useful optional alternatives.

For the user’s central motivation request, `problem131/discovery-guide.md`
combines Newton/balance recurrence construction with the reviewed reciprocal
Laurent model. `floor-covariance-barnes-design.md` explains how Bernoulli
covariance leads to the new factorial construction. `bernoulli-period-tower.md`
and its projectors connect the higher-odd analytic constants and arithmetic.

Next bounded work: review adaptive adjacent-Legendre cancellation; finish optional
kernel replay with reliable capture; consider smaller factorial certificates using
the independently reviewed outer Hermite ideal; polish proof dependencies without
turning hypotheses or numerical experiments into claims. #137’s mod-p² refinement
is now reviewed. No all-odd irrationality theorem is claimed.
