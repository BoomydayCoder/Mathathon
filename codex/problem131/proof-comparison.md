# Reviewer comparison: which ζ(3) proof is actually simplest?

Status: reviewer assessment of the mathematically checked notes in this repository. “Simpler” below means explicit prerequisite and verification cost for an olympiad/undergraduate audience, not page count or historical novelty.

| Route | Main prerequisites | Essential substantial calculations | What it explains best | Verdict |
|---|---|---|---|---|
| Zudilin's two-series proof | Partial fractions, derivatives, convergent series, elementary product bounds, factorial/lcm estimates | Two rational telescoping certificates; a saddle maximum; initial values | A recurrence transfers arithmetic from one representation to positivity of another | Best discrete explanation; most certificate verification |
| Beukers's integral proof | Polynomial derivatives, integration by parts, elementary multiple integrals, justified interchange and substitution, factorial/lcm estimates | Monomial integrals; one Möbius substitution; two square inequalities | Positivity and exponential decay have a single visible source | Lowest total calculation once undergraduate integration is available |
| Mellin bridge presentation | Same as Beukers, plus the finite beta integral or its elementary integration-by-parts proof | Rational-function/Legendre identity, followed by the same Beukers transformation | Why Apéry's binomial residues, Legendre coefficients and integral diagonal coefficients coincide | Best unified explanation; slightly longer than direct Beukers but removes the apparent split between methods |

The third row is a presentation of the same classical mathematical mechanism as Beukers, not an independent newly discovered irrationality proof. This distinction should remain explicit in any competition write-up.

## Recommended exposition sequence

Start with the integer-between-zero-and-one contradiction and the three requirements: arithmetic, positivity and decay. Present Beukers's integral proof completely, including the denominator-growth lemma; then derive the Mellin bridge and show that its rational function is exactly Zudilin's R_n. Finally derive the short recurrence certificate by pole matching as an explanation of how the recurrence is discovered. This gives a complete proof before asking the reader to learn recurrence machinery.

For a report explicitly centered on Zudilin's paper, retain his original seven-lemma sequence as the main annotated reading, but offer the integral route as an independently complete alternative. The Ball certificate can be verified in an appendix; hiding it without verification would make the exposition shorter but not the proof simpler.

## Why the integral route wins on calculation

The arithmetic computation reduces to a single identity for J_rs, with diagonal terms 2ζ(3)-2H_r^(3) and off-diagonal rational terms. The transformation to positivity needs exactly two n-fold integrations by parts and one explicitly motivated fractional-linear substitution. The sharp bound requires no multivariable calculus: use x+y≥2√(xy), then complete two squares. The same elementary lcm estimate is shared by all three routes, so it should not be counted as a special cost of the integral proof.

By contrast, the discrete proof must verify that two rather different rational summands share one recurrence. Even after deriving the short certificate conceptually by principal parts, the second degree-six certificate and its extra rational denominator remain. Zudilin's analytic maximum also needs a less transparent radical than the integral maximum.

## The advanced route: promising, but not yet a simpler complete proof here

The symmetric-square reduction in `symmetric-square.md` explains a hidden two-dimensional differential structure. Its derivation is short once one knows the recurrence. However, as presently established in this repository, it does not independently supply integrality, denominator growth, nonzero small linear forms, or identification of the limiting constant. Therefore it is an illuminating connection rather than a completed replacement proof.

A Bailey hypergeometric identity might identify the two series without recurrence certificates, but a shorter citation trades elementary rational algebra for gamma functions and contour integrals. It should only be called simpler after the specialization, contours, residues and endpoint issues are proved or stated as accepted prerequisites. The current notes do not complete that comparison.

## Classical Apéry recurrence proof

The original recurrence-based approach remains a useful fourth reference point, especially for rational approximation and continued fractions. Its recurrence is compact, but a self-contained proof must establish integrality of the relevant solution, denominator bounds for its companion, the value of the limiting ratio, and error decay. Quoting those facts makes the proof look short while leaving its hardest parts unaccounted for. The archived 1979 paper states important ingredients succinctly; it is not by itself a fully expanded elementary tutorial. Our reviewer ranking does not assume all those details away.

## Scope of review

I independently checked the Beukers transformations, denominator computation, sharp maximum and Mellin identity, and audited the companion elementary lcm argument. The Zudilin exposition has also received a separate independent review (`independent-review.md`). No claim here ranks all proofs in the literature or establishes historical novelty of any simplification.
