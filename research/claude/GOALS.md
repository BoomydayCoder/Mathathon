# Goals: Claude

Read `research/INSTRUCTIONS.md` first. It has the mission, the standards, and the rules for working alongside ChatGPT. This file says what you lead on.

You write only in `research/claude/`, on the branch `research-claude`. ChatGPT's work is in `research/chatgpt/` on `research-chatgpt`: read it every 15–30 minutes, never edit it.

## Your gate (do this first)

A complete, referee-checked exposition of **Zudilin's proof** (arXiv:math/0202159) in `exposition.md`. Every substantive step explained: purpose, motivation, prerequisites, and role in the final contradiction. Nothing else starts until this is done.

ChatGPT's gate is Beukers' integral proof. You still need to understand Beukers' proof, but read ChatGPT's exposition and check it before writing your own account of it.

## What you lead on

- [x] **Where Apéry's recurrence comes from.** Derive it from Zudilin's rational-function and series construction, and explain why the coefficients have their particular form. Show how someone could have found it.
- [x] **Better proofs of Apéry's theorem (priority 3).** Concrete, checked attempts at both an elementary demystified proof and a shorter proof using heavier machinery. Record the failures and what they reveal.
- [x] **The Fauzan ζ(5) proof.** Establish what the preprint and the two Lean formalizations actually prove and under which axioms. Explain the Hankel-determinant construction, and why it succeeds where an Apéry-style recurrence did not.
- [ ] **A simpler ζ(5) proof.** Starting from that understanding, try to shorten or restructure the argument, or to recast it in the series and recurrence language of the ζ(3) proofs.

## What ChatGPT leads on (do not duplicate; check and build on it)

- Beukers' integral proof and its links to Legendre polynomials and orthogonality.
- The general framework across \(e\), \(e^r\), \(\pi^r\), \(\zeta(2)\), \(\zeta(3)\), and the cross-area connections (priority 2).
- Extending integral and linear-form constructions to \(\zeta(7)\) and general odd \(n\), and to other constants.

## Shared

- [x] The precise dictionary between Zudilin's series and Beukers' integrals. You bring the series side; state it in your own exposition and compare with ChatGPT's version.
- [x] Each sync, check at least one substantial claim of ChatGPT's and record the result in `PROGRESS.md` under "Notes for the other agent".

## Notes

(Add notes below this line as you work.)

- Gate: `exposition.md` complete and referee-checked (§8 is the referee report). 2026-11-14.
- Recurrence: `recurrence-origin.md` (hand derivation) and exposition §7.
- Better proofs: `short-proof.md` (elementary, complete, no Ball series) and `hankel-zeta3.md` (a different proof via Hankel determinants; proof draft, two parts still sketched).
- Fauzan: `fauzan-status.md`. Sections 4-6 of the preprint were read by a subagent, not line by line by the orchestrator.
- Simpler zeta(5) proof: not achieved. Partial: the s=3 model (`hankel-zeta3.md` §7) and concrete simplification targets (`next-steps.md` item 2).
- The other agent works on branch `codex/independent-research`, folder `codex/` (not `research-chatgpt`).
