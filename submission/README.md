# Mathathon submission: team "Still doing pure maths"

Problem: the motivation of Apéry's proof that ζ(3) is irrational and the intuition behind its
recurrence (reference: Zudilin, arXiv:math/0202159). Track: Algebra / Number Theory.

| Path | What it is |
|---|---|
| `report/paper.tex`, `report/paper.pdf` | The report (single self-contained `.tex` file and its PDF). |
| `slides/slides.pdf`, `slides/slides.tex` | The slides and their Beamer source (main talk, then backup slides after "Questions?"). |
| `supplementary/checks/` | SymPy/mpmath scripts that check the identities and numbers quoted in the report, with their saved outputs (`output_*.txt`). Run with `python3 check_*.py` (needs `sympy`, `mpmath`). |
| `supplementary/team-notes/` | The team's notes cited in the report: the motivated exposition of Beukers' proof (Codex), the guess-and-check maximum-principle note, the operator-picture note, and the line-by-line exposition of Zudilin's proof. All were produced with AI tools during the event; see the report's section on AI use. |
| `supplementary/ai-brief/INSTRUCTIONS.md` | The written brief given to the overnight AI agents. |
| `supplementary/SOURCES.md` | Annotated list of every source used, how each citation was checked, and the novelty search. |

The full working repository (all agent output and history) is linked, with its commit hash, in the
report's title block.
