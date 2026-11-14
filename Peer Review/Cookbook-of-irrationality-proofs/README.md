# Cookbook of irrationality proofs (what exists so far)

This is **not a finished cookbook**. The general framework across $e$, $e^r$, $\pi^r$, ζ(2), ζ(3) was the
other agent's assignment (see `codex/problem131/irrationality-framework.md` on branch
`codex/independent-research`; not copied here). What Claude's side has is:

| document | content |
|---|---|
| `dependency-map.md` | the common scheme (arithmetic / non-vanishing and size / balance) and how each of four proofs instantiates it: Zudilin, the short elementary route, Beukers, Fauzan's ζ(5) |
| `catalan-status.md` | the Hankel recipe applied to other constants: a table of positive functionals with rational moments (ζ(2), log 2, π, ζ(3), Catalan's $G$, $L(2,\chi_{-3})$, β(4)), the exact cost/decay numbers, and why $G$ fails |
| `fauzan-status.md` | what the ζ(5) proof does, what is formalized, and the dictionary between its functional and the series language of the ζ(3) proofs |
| `literature-status.md` | what a limited arXiv/Zenodo search found on 14 Nov 2026 |

**What to check first**
1. `dependency-map.md`: is the three-demand scheme stated correctly for each proof?
2. `catalan-status.md` §1: the functionals (derivations are by a subagent; verified numerically to 32
   digits, not re-derived by the orchestrator) and the table in §2 (exact finite computations only).
3. Treat every statement about constants other than ζ(3) as numerical evidence, not proof.

The scripts behind `catalan-status.md` are on `research-claude` under
`research/claude/computations/catalan*_*.py`.
