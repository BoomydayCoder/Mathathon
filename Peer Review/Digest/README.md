# Digest: two PDFs for the team

Outcome of the peer-review process described in `../PEER_REVIEW_PROMPT.md`. Both PDFs were initially written by
Claude after reading the Claude and Codex submissions. Codex subsequently compared them with the final
research archive and, at the user's request, revised the cookbook. See
[cookbook revision notes](codex-comments-on-cookbook-and-directions.md) for scope and validation.

**Nothing here has been read by a human.** Every check listed below was done by an AI agent.

| file | pages | content |
|---|---|---|
| `zeta3-proofs.pdf` | 16 | *Why ζ(3) is irrational: the proof, built one need at a time.* Rewritten for motivation. Part I: one complete proof (Beukers' route) in which each object appears only after the difficulty that demands it, with the numbers for every rejected alternative and a note on what is forced and what is chosen. Part II: the same objects as sums (Zudilin's rational function), how Apéry's recurrence is found and why it must exist, a second proof from it, and Zudilin's own proof via Ball's series (Section 11). |
| `cookbook-and-directions.pdf` | 18 | Part I: one recipe, with complete proofs for log 2, eʳ, π², ζ(2). Part II: higher-constant obstructions, proved truncation repairs and balanced arithmetic. Part III: developed connections and next problems, including the final PNT-free determinant adaptation, Schur–Pell bounds and a random process realizing all odd-zeta correlations. Also specifies the scope of πʳ. |

Suggested reading order: `zeta3-proofs.pdf` Part I straight through, then Part II; then the
cookbook Part I quickly, Part III slowly.

## Later notes (Claude)

Written after the two main PDFs. Same caveat: checked by AI agents only.

| file | pages | content |
|---|---|---|
| `max-principle-proof.pdf` | 12 | A complete proof of the irrationality of ζ(3) using only sums of rational functions, written for motivation. The recurrence is derived (one partial fraction and one identity with the constant 6), and sign and size come from a maximum principle for the recurrence satisfied by the linear forms themselves. Sections 7 and 8 explain the decomposition (n+1)³ + n³ + 4(2n+1)³ of the middle coefficient: its meaning (a random walk with killing; a chain of resistors), how each piece arises, and in what sense it is special. The idea is known (Mingarelli; Mahboubi and Sibut-Pinote); see the note's Section 9. |
| `fifth-proof.pdf` | 7 | Record of the search for a further proof: the maximum principle, a congruence that gives non-vanishing without a sign, a sketched sign-free proof, and what the other attempts found. |
| `operator-picture.pdf` | 4 | A map of the same proof in the language of operators: the Hamiltonian H is the chain of resistors 1/n³ and also the Laplacian of the hyperbolic plane; rows and diagonal of the table are resolvents; the commutator relations. Each statement is labelled verified / agent-reported / quoted from memory. |

## How the PDFs were checked

1. **Sources.** Claude read the main proof documents of every package in full and verified the
   identities of the Codex kernel-energy proof symbolically and numerically before using it. Referee
   subagents reviewed Codex packages 01 (selected notes), 04 and 05; their reports are in
   `../Referee-reports/`. They found no mathematical error, only missing dependencies and small repairs.
2. **`zeta3-proofs.pdf`.** The first version (three compact proofs) was checked formula by formula and
   then rewritten at the owner's request, because it did not put motivation first. The rewrite was
   read by a fresh reader in two passes: as the intended reader, recording every place a question
   arose (35 points, all applied, including the removal of one section), and for correctness (every
   formula and number recomputed; no error found). The earlier three-proof version is in the git
   history (commit 796170d); its energy proof is also in `../Codex-Zeta3-03-Kernel-energy/`.
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
