# Digest: two PDFs for the team

Outcome of the peer-review process described in `../PEER_REVIEW_PROMPT.md`. Both PDFs were initially written by
Claude after reading the Claude and Codex submissions. Codex subsequently compared them with the final
research archive and, at the user's request, revised the cookbook. See
[cookbook revision notes](codex-comments-on-cookbook-and-directions.md) for scope and validation.

**Nothing here has been read by a human.** Every check listed below was done by an AI agent.

| file | pages | content |
|---|---|---|
| `zeta3-proofs.pdf` | 10 | Three complete proofs that ζ(3) is irrational: Beukers' integral; Zudilin's proof reduced to one polynomial identity and a Casoratian; an energy inequality (from the Codex submission). One common section for the arithmetic and for the lcm estimate. A derivation of the numbers 17 and 5 in Apéry's recurrence. |
| `cookbook-and-directions.pdf` | 18 | Part I: one recipe, with complete proofs for log 2, eʳ, π², ζ(2). Part II: higher-constant obstructions, proved truncation repairs and balanced arithmetic. Part III: developed connections and next problems, including the final PNT-free determinant adaptation, Schur–Pell bounds and a random process realizing all odd-zeta correlations. Also specifies the scope of πʳ. |

Suggested reading order: `zeta3-proofs.pdf` Sections 1–4 (40 minutes), then Sections 5–7; then the
cookbook Part I quickly, Part III slowly.

## How the PDFs were checked

1. **Sources.** Claude read the main proof documents of every package in full and verified the
   identities of the Codex kernel-energy proof symbolically and numerically before using it. Referee
   subagents reviewed Codex packages 01 (selected notes), 04 and 05; their reports are in
   `../Referee-reports/`. They found no mathematical error, only missing dependencies and small repairs.
2. **`zeta3-proofs.pdf`.** A second reader, given only the finished text, redid every formula by hand
   and by independent computation. Result: no mathematical error; 27 points (one wrong sentence about
   a numerical trend, thin justifications, notation clashes), all applied. The list of computations is
   in Section 8 of the PDF.
3. **`cookbook-and-directions.pdf`.** Drafted by a subagent from a detailed brief, read by Claude, then
   checked the same way by a different reader. Result: no proof marked "proved here" is wrong and all
   recomputed constants are correct; 39 points (one small proof gap, one wrong heuristic count, label
   and source-wording problems, notation clashes), applied in Claude's reviewed version. The later Codex additions have the separate source-proof and editorial-check status described in the revision notes.

## What is and is not established

- The proofs for ζ(3), log 2, eʳ, π², ζ(2) are classical results with complete proofs on the page.
- The matrix-exponential theorem (for rational M, e^M is rational only if M is nilpotent) is proved on
  the page. It is a special case of the Lindemann–Weierstrass theorem; only the elementary proof is of
  interest. We did not search the literature for this formulation.
- Results labelled **S** have written source proofs and separate internal AI review records; those proofs are not reproduced in full here. **R** labels retain the original digest team's reported-only status. Neither label means human peer review. The revised ζ(5) account retains Fauzan's stated inputs and does not claim an independent elementary proof.
- The original digest's source-retrieval statements describe its drafting team. The broader Codex archive contains additional primary-source comparisons; the revision links the relevant records. Source retrieval, an internal review, a formal build and human peer review are separate forms of evidence.

## Not included

The full Schur and Legendre–Padé proofs remain in their packages. The cookbook now includes the rank-two Pell consequence and the projection/sampler connection, but does not reproduce the entire research archive. The ζ(5) margin theorem, the ζ(7) no-go theorem and the referee report on a
claimed ζ(7) proof are summarised in Part II; the full texts are on branch `research-claude`.

## Build

`pdflatex zeta3-proofs.tex` and `pdflatex cookbook-and-directions.tex`, twice each (TeX Live with
`amsmath`, `mathtools`, `booktabs`, `microtype`, `enumitem`, `hyperref`, `lmodern`).
