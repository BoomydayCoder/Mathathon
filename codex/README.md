# Mathathon independent research

Work in progress on `codex/independent-research`. Read [CHECKPOINT.md](CHECKPOINT.md)
for the latest state and [HANDOFF.md](HANDOFF.md) for remaining work. Results derived
in this session are not claims of historical novelty.

## Three reading routes

**Understand and present Apéry's proof.** Start with the [discovery guide](problem131/discovery-guide.md).
It derives the rational square from its constraints, explains the recurrence
through a quadratic Newton grid and finite balance, and gives a complete elementary
irrationality proof. A single positive Laurent polynomial now connects its
integer coefficients to the reciprocal-power Beukers integral; one AM–GM
minimum explains their reciprocal scales. The resistor chain gives a further
exact interpretation of the error. For the original paper's sequence of steps, use
[Zudilin explained](problem131/zudilin-explained.md).

**Follow the zeta-five simplification without the prime number theorem.**
Start with the [PNT-free guide](problem131/pnt-free-zeta5-guide.md).
It connects residue-class covariance to Bernoulli polynomials and superfactorials,
constructs an explicit normalizer, and completes the rationality contradiction.
A 25-term rational certificate controls every prime individually; both the
small-prime and large-prime repairs are included. The arithmetic rate is below
1.379, while the reviewed equilibrium estimate supplies decay faster than 1.38
in the quadratic exponent. A reviewed primitive-part argument keeps the final
degree at \(37n\) when \(K=40n\); the certificate denominator \(40000\) is only
auxiliary bookkeeping. Explicit floor rounding gives the same degree without
any added upper arithmetic cost. This computer-assisted alternative retains Fauzan's
determinant and local arithmetic construction, including its refined rank bound.
For fuller analytic motivation, read the [equilibrium guide](problem131/fauzan-equilibrium-guide.md).
The earlier [plain-rank route](problem131/fauzan-equilibrium-guide.md#12-spend-a-little-analytic-margin-to-simplify-the-local-arithmetic-too)
remains available when simpler local arithmetic is preferred; it uses the PNT.
The [optional outer-kernel improvement](problem131/outer-gain-budget-update.md)
adds a local lemma for a larger margin.

**Explore results and promising boundaries.** The [research results](problem131/research-results.md)
collect precise theorem statements, proofs, obstructions, and links. The common
analytic budget works for every fixed odd index, but the arithmetic normalizations
needed for all-odd irrationality remain unproved. The negative results identify
specific constructions that fail; they are not global impossibility claims.

## Secondary project: exact prime-power Fuss–Catalan valuations

The [problem #137 discovery guide](problem137/discovery-guide.md) gives a standalone
combinatorial proof of the exact valuation and its leading unit modulo the prime.
Expand each vertex into a perfect prime-ary block, classify the smallest group
orbits, and count their depth-constrained skeletons by an explicit cycle-splicing
bijection with inverse. The proof uses no factorial-valuation formula or Lagrange
inversion. It addresses the exactness gap named in the [official statement](sources/problem137-official.md);
[the primary-source comparison](problem137/konvalinka-comparison.md) credits the
binary Deutsch–Sagan argument and distinguishes Konvalinka’s existing arithmetic
formula and tree-action divisibility proof. Historical novelty and weighted
extensions remain unresolved.

## Connections worth following

| Question | Read |
|---|---|
| Why precisely this square and this recurrence? | [Design classification](problem131/zudilin-design-classification.md); [Newton grid and finite diagonal balance](problem131/newton-symmetric-square-review.md) |
| Why does a growing integer solution create a small positive error? | [The same Laurent polynomial in positive and inverse powers](problem131/apery-laurent-reciprocal.md), with [independent review](problem131/apery-laurent-independent-review.md); [resistor chain](problem131/discrete-green-motivation.md) |
| Can advanced machinery make the structure clearer? | [Elliptic periods and modular branch cancellation](problem131/apery-modular-elliptic.md); [completed independent review](problem131/apery-modular-review.md) |
| How can the same moments control real size and prime denominators? | [Bernoulli local/global dictionary](problem131/bernoulli-local-global-dictionary.md); [local arithmetic reconstruction](problem131/fauzan-local-arithmetic-review.md) |
| Why do prime-allocation floors suggest a factorial product? | [Centered indicators, Brownian-bridge covariance, and superfactorials](problem131/floor-covariance-barnes-design.md); [PNT-free assembly](problem131/pnt-free-zeta5-guide.md) |
| What do positive matrices reveal—and where do they fail? | [Real-rooted zeta pencils](problem131/real-rooted-zeta-pencils.md); [Wilson-polynomial one-pole solution and exact denominator obstruction](problem131/one-pole-schur-structure.md) |
| What connects this to prime-adic Fuss–Catalan questions? | [Carry cancellation and digit factorization](problem131/apery-digit-arithmetic.md); [elementary character rule](problem131/apery-cm-character.md), now [independently reviewed](problem131/apery-cm-independent-review.md) |
| What transfers to other irrationality proofs? | [Shared exponential, pi-squared, and zeta coefficients](problem131/shared-coefficients-exp-pi-zeta.md); [zeta-two/zeta-three balance](problem131/zeta2-zeta3-balance-unification.md); [powers of pi](problem131/pi-powers.md) |

The new [monomial-kernel design](problem131/monomial-kernel-design.md) makes the
common arithmetic pairing explicit: diagonal monomials select the zeta constant,
while off-diagonal terms are rational. Its [dedicated review](problem131/monomial-kernel-review.md)
is complete; choosing partners with sufficient nonzero decay remains a separate
problem. Follow-up periodic-weight and mixed-partner constructions have separate
review status.

For the geometric interpretation, [the Legendre/relative-period note](problem131/laurent-legendre-relative-periods.md)
records the necessary boundary term: the real error generating function satisfies
an inhomogeneous equation. The Laurent critical values alone are not a derivation
of the full recurrence.

Detailed higher-derivative, truncation, and prime-band investigations are indexed
in [research results](problem131/research-results.md); their distinctions between
actual failure and insufficient estimates matter.

## What has been verified

The elementary zeta-three guide is a complete written proof. The advanced modular
route now includes a reviewed direct global covering argument, with standard
modular-form and covering-theory inputs stated explicitly.

For the zeta-five alternative, the [completed assembly audit](problem131/fauzan-simplified-proof-audit.md)
checks the shared positive moments, local arithmetic, real determinant transfer,
and final contradiction. The new [normalizer review](problem131/pnt-free-normalizer-independent-review.md)
checks the superfactorial replacement for prime averaging, including every higher
prime power and eventual coefficientwise integrality. Its exact sparse certificate
proves the arithmetic rate is below 6893/5000; the [rational energy review](problem131/rational-energy-review.md)
certifies the analytic magnitude exceeds 2771/2000. The
[primitive-part and rounding review](problem131/prime-exponent-rounding-independent-review.md)
removes the unnecessary determinant power. The
[new guide](problem131/pnt-free-zeta5-guide.md) assembles these dependencies.
The older plain-rank route and its exp(−13n²) bound remain valid with their own
normalizer and PNT input; the outer-kernel option gives exp(−26n²).
These are computer-assisted derivations retaining Fauzan's construction, not
claims of elementary proofs, historical novelty, formal verification, or external
acceptance.
No all-odd irrationality theorem is claimed. See [source provenance](sources/odd-zeta-provenance.md).

Separately, the pinned Romik Lean formalization **was built locally**, and its
ordinary-series statement and axiom checks passed. The completed kernel replay
checked all 71 project modules against their imported environments. The compiled
theorem retains the explicit prime number theorem axiom, along with Lean's
standard logical axioms. This was not a fresh replay of all Mathlib constants
and does not formalize our new equilibrium, PNT-free normalizer, or Python
certificates. The [verification report](problem131/lean-local-verification.md)
links the exact commit, successful build log, axiom output, and durable replay
result.

Symbolic scripts check identities; finite experiments alone do not prove all-index
statements. Reconstructed motivation explains how one might design a proof, not
what its historical authors thought. All work is under `codex/`; no Claude-owned
files have been changed. The exact [#131](sources/problem131-official.md) and
[#137](sources/problem137-official.md) statements were retrieved from the public
API at 08:46 UTC on 14 November 2026; raw JSON and SHA256 provenance are saved.
Earlier HTML snapshots are application shells, not statements. The prime-power
tree theorem matches #137's arity scope; the primary-source comparison with
Konvalinka is now recorded in [the comparison note](problem137/konvalinka-comparison.md).

## Reproduce

Most scripts use Python 3 with SymPy and mpmath; their headers state dependencies.
The [Newton closure check](computations/check_newton_diagonal_closure.py) and
[modular review check](computations/check_apery_modular_review.py) use exact identities.
The [independent Laurent check](computations/check_apery_laurent_independent_review.py)
verifies the reciprocal measure and covariance separately from the main derivation.
The sparse superfactorial certificate, its independent normalizer checks, and
the endpoint and rational-energy certificates use only the Python standard
library and local interval primitives. From the repository root:

```bash
python codex/computations/certify_equilibrium_endpoints.py
python codex/computations/check_rational_energy_primitives.py
python codex/computations/certify_equilibrium_energy_rational.py
python codex/computations/check_floor_covariance_barnes.py
python codex/computations/check_sparse_barnes_certificate.py
python codex/computations/check_pnt_free_normalizer_review.py
python codex/computations/check_prime_exponent_rounding_review.py
```

The sparse check verifies 561 exact breakpoints and propagates the inequality
to the infinite tail by an exact period identity. The independent normalizer
check tests the factorial identities and repair bounds. For the older PNT route,
run [check_table_free_arithmetic.py](computations/check_table_free_arithmetic.py).
The [rational energy review](problem131/rational-energy-review.md) explains every
remainder and outward-rounding bound; the saved JSON records exact dyadic endpoints.
Older [Carlson](computations/certify_equilibrium_energy.py) and
[virial](computations/certify_equilibrium_energy_elementary.py) checks use
`python-flint`/Arb and provide additional cross-checks. The rational certificate
shares the endpoint box and proved potential identities; it is not an independent
reconstruction of those inputs.

For the secondary project, [the cycle-splicing check](problem137/check_cycle_splicing.py)
exhaustively tests the bijection and its inverse; [the independent orbit check](problem137/check_prime_power_orbit_review.py)
counts skeleton depth profiles separately from the closed formula.
