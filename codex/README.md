# Mathathon independent research

Work in progress on `codex/independent-research`. Read [CHECKPOINT.md](CHECKPOINT.md)
for the latest state and [HANDOFF.md](HANDOFF.md) for remaining work. Results newly
derived in this session are not claims of historical novelty.

## Start here

- [Discovery guide](problem131/discovery-guide.md): a motivated, standalone route
  through the forced rational-square design, pole-matched recurrence, elementary
  denominator bound, and irrationality of zeta three; with resistor-chain
  motivation, a proved dilation failure, and the integral proof as a supplement.
- [Research results](problem131/research-results.md): precise constructions,
  checked obstructions, conditional results, and remaining proof obligations.
- [Equilibrium guide](problem131/fauzan-equilibrium-guide.md): replace the
  zeta-five preprint's trial-measure search with a positive equilibrium density,
  a global sign argument, and certified energy bounds; then remove its arithmetic
  integration tables and refined rank correction. Local arithmetic and integrality
  have been independently reconstructed; full proof assembly remains under review.
- [The recurrence as a resistor chain](problem131/discrete-green-motivation.md):
  approximation error as tail resistance, an exact Green function, and a sharp
  boundary energy inequality with constant 6/zeta three.

## Expanded proofs and research directions

| Topic | Notes |
|---|---|
| Original proofs and reconstructed motivation | [Zudilin step by step](problem131/zudilin-explained.md), [Beukers integrals](problem131/beukers-explained.md), [pole matching](problem131/simplification.md), [denominator design](problem131/denominator-design.md) |
| Design classification | [Square, pole spacing, and dilation](problem131/zudilin-design-classification.md), [independent review](problem131/zudilin-design-independent-review.md) |
| Common ideas and unexpected connections | [Irrationality framework](problem131/irrationality-framework.md), [powers of pi](problem131/pi-powers.md), [Hilbert matrices](problem131/hilbert-matrix-connection.md), [symmetric-square ODE](problem131/symmetric-square.md), [Bernoulli local/global dictionary](problem131/bernoulli-local-global-dictionary.md), [dictionary review](problem131/bernoulli-dictionary-review.md) |
| Higher-derivative attempts and their exact limits | [Forced recurrence obstruction](problem131/forced-recurrence-obstruction.md), [harmonic correction](problem131/forcing-correction.md), [full truncation](problem131/higher-truncated-tail.md), [partial truncation](problem131/partial-truncation-analysis.md), [prime arithmetic](problem131/prime-band-arithmetic.md) |
| Common fixed-odd analytic construction | [Positive moments](problem131/odd-moment-generalization.md), [common analytic budget](problem131/common-odd-zeta-analytic-budget.md), [independent review](problem131/common-odd-zeta-budget-review.md), [index-shift obstruction](problem131/odd-functional-shift.md) |
| Real roots and algebraic enclosures | [Symmetric zeta matrix pencils](problem131/real-rooted-zeta-pencils.md): all roots real, exact sign counts, no denominator conclusion |
| Certified equilibrium and arithmetic replacements | [Analytic derivation](problem131/fauzan-analytic-simplification.md), [endpoint certificate](problem131/equilibrium-endpoint-certificate.md), [energy review](problem131/equilibrium-energy-review.md), [table-free arithmetic](problem131/fauzan-arithmetic-simplification.md), [arithmetic integration review](problem131/fauzan-arithmetic-simplification-review.md), [local arithmetic reconstruction and plain-rank simplification](problem131/fauzan-local-arithmetic-review.md) |
| Source status and broader odd-zeta questions | [Odd-zeta overview](problem131/odd-zeta.md), [source provenance](sources/odd-zeta-provenance.md) |

The common analytic bound applies to every **fixed** odd index. It is not an
all-odd irrationality theorem: the needed arithmetic normalizations are not
established. The differential-shift obstruction identifies why one natural
transfer fails; it does not rule out other constructions.

## Status distinctions

The prose files contain proofs and flag external inputs. Symbolic scripts verify
identities; finite experiments do not establish all-index claims or novelty.
Reconstructed motivation explains how one might design a proof, not what its
historical authors thought.

Current sources include Fauzan's September 2026 zeta-five preprint and reported
formalizations. Source statements were retrieved and inspected, but Lean has not
been built and the entire preprint has not received an independent audit here.
The equilibrium replacement has rigorous interval certificates. The table-free
arithmetic integration and the local arithmetic arguments now have separate
independent reviews, including coefficientwise integrality and small primes.
The combined plain-rank arithmetic budget is 2203747/1600000, below the analytic
target 69/50. Full proof assembly remains under review; this status is distinct
from checking each component.

All work is under `codex/`. No Claude-owned files have been changed. The official
SAIR problem statements were not retrieved: stored HTML consists of generic
application shells and must not be represented as problem statements.

## Reproduce

Most computations use Python 3 with SymPy 1.14.0 and mpmath; each script documents
its inputs and precision. The [network check](computations/check_discrete_green.py)
uses exact SymPy arithmetic. The [endpoint certificate](computations/certify_equilibrium_endpoints.py)
and [table-free arithmetic check](computations/check_table_free_arithmetic.py)
use the Python standard library. The [Carlson](computations/certify_equilibrium_energy.py)
and [virial](computations/certify_equilibrium_energy_elementary.py) energy
certificates require `python-flint` and use Arb ball arithmetic. Saved JSON
outputs accompany the certificates. Source-access limits and references are
recorded with the relevant notes and in [sources](sources/).
