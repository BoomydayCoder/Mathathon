# Goals: ChatGPT

Read `research/INSTRUCTIONS.md` first. It has the mission, the standards, and the rules for working alongside Claude. This file says what you lead on.

You write only in `research/chatgpt/`, on the branch `research-chatgpt`. Claude's work is in `research/claude/` on `research-claude`: read it every 15–30 minutes, never edit it.

## Your gate (do this first)

A complete, referee-checked exposition of **Beukers' integral proof** of the irrationality of ζ(2) and ζ(3) (Bull. London Math. Soc. 11 (1979) 268–272) in `exposition.md`. Every substantive step explained: purpose, motivation, prerequisites, and role in the final contradiction. Nothing else starts until this is done.

Claude's gate is Zudilin's proof. You still need to understand Zudilin's proof, but read Claude's exposition and check it before writing your own account of it.

## What you lead on

- [ ] **Why the integrals work.** Explain how Beukers' integrals and the Legendre polynomials could have been found, what orthogonality buys, and where the denominator bound and the size estimate come from.
- [ ] **The general framework (priority 2).** Compare proofs for \(e\), \(e^r\), \(\pi^r\), \(\zeta(2)\), \(\zeta(3)\). Organise the shared techniques, show how they fit together, and state exactly where the framework breaks.
- [ ] **Cross-area connections.** Test whether links to orthogonal polynomials, continued fractions, differential equations, algebraic geometry, modular forms or combinatorics can be used to discover or generalise a proof, not just observed.
- [ ] **Beyond ζ(5).** Extend integral and linear-form constructions towards \(\zeta(7)\), general odd \(n\), and other constants where irrationality is open. Record precisely what fails and why.

## What Claude leads on (do not duplicate; check and build on it)

- Zudilin's proof and the origin of Apéry's recurrence and its coefficients.
- Simplified and alternative proofs of Apéry's theorem (priority 3).
- The Fauzan ζ(5) preprint, its Lean formalizations, and attempts at a simpler ζ(5) proof.

## Shared

- [ ] The precise dictionary between Zudilin's series and Beukers' integrals. You bring the integral side; state it in your own exposition and compare with Claude's version.
- [ ] Each sync, check at least one substantial claim of Claude's and record the result in `PROGRESS.md` under "Notes for the other agent".

## Notes

(Add notes below this line as you work.)
