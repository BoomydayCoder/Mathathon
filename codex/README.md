# Mathathon independent research

Work in progress on `codex/independent-research`. Read [CHECKPOINT.md](CHECKPOINT.md)
for the latest state and [HANDOFF.md](HANDOFF.md) for remaining work. Results derived
in this session are not claims of historical novelty.

## Three reading routes

**Understand and present Apéry's proof.** Start with the [discovery guide](problem131/discovery-guide.md).
It derives the rational square from its constraints, explains the recurrence
through a quadratic Newton grid and finite balance, and gives a complete elementary
irrationality proof. The resistor-chain and integral interpretations explain the
same error in two further ways. For the original paper's sequence of steps, use
[Zudilin explained](problem131/zudilin-explained.md).

**Follow the strongest zeta-five simplification.** Read the [equilibrium guide](problem131/fauzan-equilibrium-guide.md),
then the [end-to-end audit](problem131/fauzan-simplified-proof-audit.md).
An explicit equilibrium measure replaces sixteen trial components; table-free
integration and a plain-rank estimate simplify the arithmetic. The assembled
computer-assisted alternative derivation retains Fauzan's determinant and local
arithmetic construction. The modified normalizer is defined explicitly.

**Explore results and promising boundaries.** The [research results](problem131/research-results.md)
collect precise theorem statements, proofs, obstructions, and links. The common
analytic budget works for every fixed odd index, but the arithmetic normalizations
needed for all-odd irrationality remain unproved. The negative results identify
specific constructions that fail; they are not global impossibility claims.

## Connections worth following

| Question | Read |
|---|---|
| Why precisely this square and this recurrence? | [Design classification](problem131/zudilin-design-classification.md); [Newton grid and finite diagonal balance](problem131/newton-symmetric-square-review.md) |
| Why does a growing integer solution create a small positive error? | [Resistor chain and Green function](problem131/discrete-green-motivation.md); [integral supplement](problem131/beukers-explained.md) |
| Can advanced machinery make the structure clearer? | [Elliptic periods and modular branch cancellation](problem131/apery-modular-elliptic.md); [completed independent review](problem131/apery-modular-review.md) |
| How can the same moments control real size and prime denominators? | [Bernoulli local/global dictionary](problem131/bernoulli-local-global-dictionary.md); [local arithmetic reconstruction](problem131/fauzan-local-arithmetic-review.md) |
| What do positive matrices reveal—and where do they fail? | [Real-rooted zeta pencils](problem131/real-rooted-zeta-pencils.md); [Wilson-polynomial one-pole solution and exact denominator obstruction](problem131/one-pole-schur-structure.md) |
| What connects this to prime-adic Fuss–Catalan questions? | [Carry cancellation and digit factorization](problem131/apery-digit-arithmetic.md); [square/discriminant factorization modulo primes](problem131/apery-digit-square-review.md) |
| What transfers to other irrationality proofs? | [Common framework](problem131/irrationality-framework.md); [powers of pi](problem131/pi-powers.md) |

Detailed higher-derivative, truncation, and prime-band investigations are indexed
in [research results](problem131/research-results.md); their distinctions between
actual failure and insufficient estimates matter.

## What has been verified

The elementary zeta-three guide is a complete written proof. The advanced modular
route now includes a reviewed direct global covering argument, with standard
modular-form and covering-theory inputs stated explicitly.

For the zeta-five alternative, the [completed assembly audit](problem131/fauzan-simplified-proof-audit.md)
checks the positive moments, local arithmetic including small primes, coefficientwise
integrality, prime averaging, certified energy, discrete determinant bound, and
rationality contradiction together. Its new normalizer has arithmetic budget
2203747/1600000, below the analytic target 69/50. The route is computer-assisted
and retains Fauzan's construction; it is not claimed elementary, historically
novel, formally verified, or externally accepted. Lean has not been built locally.
No all-odd irrationality theorem is claimed. See [source provenance](sources/odd-zeta-provenance.md).

Symbolic scripts check identities; finite experiments alone do not prove all-index
statements. Reconstructed motivation explains how one might design a proof, not
what its historical authors thought. All work is under `codex/`; no Claude-owned
files have been changed. Official SAIR statements were not retrieved: stored HTML
contains application shells and must not be represented as exact statements.

## Reproduce

Most scripts use Python 3 with SymPy and mpmath; their headers state dependencies.
The [Newton closure check](computations/check_newton_diagonal_closure.py) and
[modular review check](computations/check_apery_modular_review.py) use exact identities.
The [endpoint certificate](computations/certify_equilibrium_endpoints.py) and
[table-free arithmetic check](computations/check_table_free_arithmetic.py) use the
standard library. The [Carlson](computations/certify_equilibrium_energy.py) and
[virial](computations/certify_equilibrium_energy_elementary.py) energy certificates
use `python-flint`/Arb ball arithmetic. Saved JSON certificates accompany them;
the end-to-end audit records fresh reruns and the limits of their independence.
