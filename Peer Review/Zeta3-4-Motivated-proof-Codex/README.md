# ζ(3): "Why ζ(3) is irrational — a proof built from its requirements" (Codex)

`zeta3-motivated-proof.pdf` — 58 pages, PDF metadata: author "Prepared with Codex", created
14 Nov 2026. Added to Peer Review at the repository owner's request. The file is unmodified.

**This is the other agent's document, not Claude's.** Claude read it in full (text extraction of all
58 pages) before adding it, and spot-checked the points listed under "Checked while reading". That is
a reading, not a referee report. No human has refereed it.

## What it claims

Its own abstract and disclaimer: an explanatory reconstruction of **Beukers' proof** of Apéry's
theorem, in which each requirement (integrality, non-vanishing, smallness) motivates the next
construction. It claims no new proof of Apéry's theorem, no proof for ζ(7) or Euler's constant, and no
independent proof for ζ(5).

## Contents

| part | sections | content |
|---|---|---|
| the proof | 1–8 | kernel −log(xy)/(1−xy); moments give D_n³·J(P) ∈ Z + Zζ(3) for any integer polynomial P of degree ≤ n; orthogonality forces the shifted Legendre polynomial; identity K = ∫₀^∞ dt/((x+t)(1+yt)) and two integrations by parts give a positive triple integral; maximum (√2−1)⁴; lcm lemma; contradiction |
| connections | 9–14 | Mellin transform of P_n is Zudilin's rational function; Apéry's recurrence from a degree-4 interpolation problem; Beukers' cube integral; Padé approximation of the logarithm; the 5/2 central-binomial series by telescoping |
| relatives | 15 | the same machinery for log 2 (and log((a+1)/a)), ζ(2), e^q, π |
| failed extensions | 16–18, 20–22, App. A | higher powers of the logarithm for ζ(5), ζ(7) (with a proved asymptotic, Prop. 16.1); higher zero multiplicity; Euler's constant (Sondow's integral, and a logarithm-free variant); kernel determinants; a double-pole construction |
| ζ(5) remarks | 19, App. B, D, E | the Hankel-determinant strategy of Fauzan's preprint; a weighted-prime subsequence lemma that replaces the prime number theorem in one summation step; a 5-interval comparison measure with U = −1323/1000 |

## What a referee should look at first

1. **Sections 2–7** are the whole proof. Section 8 is a half-page summary.
2. **Lemma 6.1–6.2** (the lcm bound). The identity N! = ∏_{k≤N} D_{⌊N/k⌋} gives "D_n < Cⁿ for
   infinitely many n" for every C > e, with no Chebyshev-type estimate and no prime number theorem.
   It works only because J_n > 0 for every n, so a subsequence is enough.
3. **Section 4.3** (Fubini/Tonelli justification) and the AM–GM step (5.1).
4. **Appendices D and E** contain the only claims that go beyond classical material. They take the
   local arithmetic and the real estimate of Fauzan's preprint as inputs, and they need separate
   refereeing; Appendix E is computer-assisted and its certificate and code are **not in this folder**.

## Checked while reading (by Claude, by hand)

* Lemma 6.2: from N! = ∏ D_{⌊N/k⌋} and D_m ≥ C^m for m ≥ M one gets
  log N ≥ log C·(log N − log M − 1/M), impossible for C > e. Correct.
* J_1 = 10ζ(3) − 12 and the table in 9.4 agree with the values in `Zeta3-1` and `Zeta3-2`
  (there F_n = J_n, u_n = 2a_n).
* The certificate q_n(t) = 4(2n+1)(−2t² + t + (2n+1)²) in Section 10 is the same as s_n(t) in
  `Zeta3-2-Elementary-proof/recurrence-origin.md`; the Casoratian 6/(n+1)³ matches 24/n³ there after
  the factor 2 and index shift.
* 27(√2−1)⁴ < 27/32 and 9β⁵ < 1 (Section 15.2) are correct.

Not checked: Appendix A (asymptotics of the failed extension), the numerical tables, Appendices C–E,
Section 22.

## Relation to the other folders

| | non-vanishing | size | lcm input |
|---|---|---|---|
| this document (Beukers) | positive triple integral | maximum of one function, (√2−1)⁴ | D_n < 3ⁿ infinitely often (Lemma 6.2) |
| `Zeta3-1` (Zudilin) | Ball's positive series | Stirling-type bound | "D_n < 3ⁿ", quoted |
| `Zeta3-2` (elementary) | Casoratian + recurrence | growth of u_n from the recurrence | period-6 Chebyshev bound, all n |
| `Zeta3-3` (Hankel) | Gram determinant | energy bound | none of this type |

The lcm lemma of this document could replace the period-6 bound in `Zeta3-2` (there F_n > 0 for all n
as well), with a wider margin: e³ ≈ 20.1 against 33.5.
