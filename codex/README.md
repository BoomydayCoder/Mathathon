# Mathathon independent research

Work in progress, on codex/independent-research. Read CHECKPOINT.md for the latest
state and HANDOFF.md for remaining work. No historical novelty is claimed.

## Start here

1. problem131/zudilin-explained.md — original paper explained step by step,
   including motivations, exact certificates, and a self-contained lcm estimate.
2. problem131/beukers-explained.md — integral proof and the exact Mellin bridge.
3. problem131/simplification.md — discover the recurrence by matching poles.
4. problem131/irrationality-framework.md — e^r, pi, zeta2, Padé and Legendre links.
5. problem131/symmetric-square.md — derive the hidden second-order ODE.
6. problem131/denominator-design.md — denominator/derivative tradeoffs and a
   shorter elementary lcm proof.
7. problem131/forced-recurrence-obstruction.md — a proved reason the naive
   higher-derivative extension cannot decay exponentially for zeta5.
8. problem131/odd-moment-generalization.md and odd-hankel-experiments.md — a
   positive-moment construction for all odd values and exact small experiments.
9. problem131/odd-zeta.md — higher odd value connections and current-source status.

## Status distinctions

The prose files contain proofs and flag external theorem inputs. Symbolic scripts
verify identities; finite experiments do not establish all-n claims or novelty.
The authorial motivation is reconstructed mathematically, not claimed as a
historical account of what Zudilin or Beukers thought.

Current sources include Fauzan's September 2026 zeta5 preprint and reported
formalizations. We retrieved and inspected source statements but have not built
Lean or completed an independent audit of the entire preprint. Earlier chat
wording calling zeta5 open was outdated. See sources/odd-zeta-provenance.md.

All work is under codex/. No Claude-owned files have been changed. The official
SAIR problem statements were not retrieved: stored HTML is generic application
shells and must not be represented as problem statements.

## Reproduce

Python 3 with SymPy 1.14.0 and mpmath. Run the check_*.py files in computations/.
The exploration script documents numerical precision separately; its JSON stores
exact integer polynomial coefficients. References and source-access limits are
recorded with the relevant notes and in sources/.
