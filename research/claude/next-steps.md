# Next steps for a 40-hour team project (ranked)

Ranking is by (expected value) x (probability of a checkable result). Each item says what "done" means.

## 1. Finish and referee the Hankel-determinant proof of Apéry's theorem (`hankel-zeta3.md`)
Most promising: the core is already proved and exactly verified for K<=100.
- Write out Lemma L8 (small primes) and the analytic bound (d) in full; or replace (d) by the exact
  equilibrium computation (conjecture: log Delta_K(zeta(3)) = -(9/2 log 2) K^2 + o(K^2); the support of the
  equilibrium measure is [0,b] with 3b+1-sqrt(1+b^2)=2 lambda).
- Sharpen L8 to e_p = 3K floor(log_p K) (true in all data) so that the proved normaliser beats the decay
  at moderate K, giving an *effective* statement.
- Search the literature for an s=3 Hankel proof (Zudilin's and Krattenthaler-Rivoal-Zudilin's work on
  Hankel determinants of zeta values is the place to look) before any novelty claim.
- Done = a complete, human-refereed write-up, plus a Lean-checkable statement list.

## 2. Use the s=3 model to simplify Fauzan's zeta(5) proof
- Find a global (prime-independent) Newton-type basis adapted to the weight D_N^6/D_K. The two-sided
  Newton basis reproduces Fauzan's outer-range valuations for the largest primes only (loss 0.5 K^2 in
  the inner range). Target: one entrywise valuation bound + the weights lemma L5 replacing Prop. 4.1,
  Lemma 4.2 and Prop. 4.3.
- Replace the 684-cell energy certificate by the exact equilibrium measure. Our solver gives
  U_true = -1.3855 for Fauzan's field versus the proved -1.3670: the margin would double (1.3% -> 2.6%).
- Done = a zeta(5) proof whose arithmetic section is as short as §4 of `hankel-zeta3.md`.

## 3. Decide the fate of zeta(7) in the Hankel family
- Prove the cost formula s(s+4)/(2(s+2)) K^2 for W=1 (upper and lower bounds), and the s-independence
  of the decay. That would be a clean no-go theorem for the unmodified family.
- Look for functionals whose p-adic singular part per node is smaller than p^{-s} (e.g. combine several
  odd s with elimination as in FLSZ; or nodes on a sparser set than the squares; or weights with faster
  decay than e^{-2 pi y}).
- Done = either a no-go theorem with explicit hypotheses or a family with negative total at K<=80.

## 4. Formalize the short proof of Apéry's theorem (`short-proof.md`)
- One degree-4 polynomial identity, partial fractions, a Casoratian, and the period-6 lcm bound: a good
  Lean target (no PNT needed). Done = Lean statement `Irrational (riemannZeta 3)` with standard axioms.

## 5. Linear forms versus polynomials
- The 1x1 Hankel case gives linear forms with normalised exponent >= 0 (heuristic in NOTEBOOK). Make this
  a theorem: no choice of W gives an irrationality proof at h=1. It explains precisely why determinants
  (the K^2 scale) are needed.

## 6. Dictionary questions
- Is Ball's series the Mellin image of Beukers' triple integral (exposition §9)? Is there a
  Beukers-type multiple integral for the Hankel determinant Delta_K(zeta(3))? (Andréief gives a K-fold
  integral with a squared Vandermonde; compare with Sorokin/Vasilyev-type integrals.)
- The weight omega_K is a Wilson-type weight with a degenerate parameter (frontier subagent's remark,
  unchecked): do the orthogonal polynomials have a closed form, giving exact formulas for Delta_K?

## 7. By-products worth a short note
- Companion forms G_n = sum R_n(t) in 1, zeta(2) with Apéry's zeta(3) recurrence and right side -8(2n+1).
- Hand derivation of the recurrence; 17n^2+17n+5 = 4(2n+1)^2 + n(n+1) + 1.
