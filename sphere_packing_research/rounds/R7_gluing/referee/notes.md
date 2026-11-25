# Referee notes on R7_gluing (adversarial check, 07:50-08:01 UTC)

Tags: [L] checked against paper text, [P] proof checked line by line, [N] numerical (code in this directory), [H] heuristic, [GAP] gap found.
Code: `ref_lp.py` (LP written directly in the y-variables of eq. (6), brute force over supports; exact sympy resultants),
`ref_even.py` (exhaustive avoid-lcm search, n = 2 mod 4 allowed), `ref_count.py` (full recount of admissible families with two codes),
`ref_misc.py` (window table, filter sanity, cube-minus-top), `ref_final.py` (Parseval, counts). Outputs: `out_*.txt`.

## 1. Definitions vs. paper (arXiv:2610.12619 text file) [L]
* Sec 2.2 (lines 289-356): I = log prod_{{m,n} in N} |Res|, I_inc(S) = sum over pairs with {m,n} cap S nonempty. Apostol formula as quoted. OK.
  Apostol formula re-checked against exact integer resultants for all 190 divisor pairs of 60, 105, 45, 72, 75: 0 mismatches [N].
* Eq (6) (line 558), (8)-(11) (lines 773-800), Def 2.6 (lines 642-664), Prop 4.6 (lines 2038-2045): all quoted correctly by R7.
  S ranges over nonempty proper subsets; y_n only constrained by sum c_n y_n = 0 and |y_1| << log C.
* Correction 1: the paper's stated reason for 2 not in N is Sec 4.1 (all nontrivial components complex, formula (26)), not the
  criticality x_hi(2) = 1/2. R7's "this is why 2 is excluded" is an interpretation.
* Correction 2 [GAP, closed below]: Def 2.6 forbids only n = 2. R7's even-M pools excluded every n = 2 mod 4 (6, 10, 30, ...), which the paper allows.

## 2. Independent LP [N]
* LP in y-variables straight from eq. (6): margins agree with R7 for {1,3} (0.0719), {1,19,57} (0.2227), {1,3,5} (-0.4977),
  full cube of 105 (0.9147), and the four "best with 1 in N" sets for M = 105, 315, 525, 1155 (0.0182, 0.1132, 0.0029, 0.0180), plus 165, 195.
  With exact resultants instead of Apostol: identical.
* Full recount without any prefilter: the complete family of support-admissible N (lcm = M) for M = 45, 105, 315, 60
  has 7, 20, 74, 5 members with both codes, same sets, margins equal to 1e-14. None has C < phi(M); none avoids M.
* Minor: R7 says "all 30641 [sets] with C < 480 rejected" for M = 1155; the number of index sets with lcm = 1155 and C < 480 is 25148
  (30641 is the number tested before the third admissible set was found). Conclusion unaffected.

## 3. Proofs [P]
* gain = 1/mu(N), Parseval: correct (numerically re-verified for 5 sets incl. even/non-squarefree ones).
* (11) from (6),(8),(10): re-derived, correct. Cut / partition / vertex conditions: y (v) cancels exactly because the parts partition N. Correct.
* Per-prime window: every edge crossing {v_q(n) < e} | {v_q(n) = e} is a q-edge from n = m q^e down to m q^i, total weight per top element
  phi(m) q^{e-1} log q = phi(n) log q/(q-1). H_2(x)/(1-x) is increasing, so the admissible set is exactly x < x_hi(q). x_hi(q) < 2/q for q >= 3: argument correct.
  Table reproduced (x_hi(19): q*x_hi = 1.675, R7 has 1.676).
* Theorem A: correct, for ODD M only. For even M the same argument needs sum of x_hi with x_hi(2) = 1/2, and {2,3,5} (1.163), {2,3,7} (1.100),
  {2,3,11} (1.034), {2,3,13} (1.015) are NOT excluded by the counting argument.
* Theorem B: correct, and correctly described as far too weak for a no-go.
* The remark "every condition gets easier as N gets lighter / none bounds the gain from above" is a heuristic remark, not a statement with a proof.

## 4. Extension of the search (closing gap 2) [N]
Exhaustive over all N inside divisors(M) minus {2, M} with lcm(N) = M (so n = 2 mod 4, n > 2 allowed):
M in {12, 20, 24, 28, 30, 36, 40, 42, 60, 66, 78, 84, 90, 120, 126, 132, 150, 156, 168, 180, 210, 252, 270, 300, 330, 390, 462,
360, 420, 540, 600, 630, 660, 1050} and odd M in {105, 315, 525, 945, 1155}: in every case ZERO sets pass even the exact per-prime cut
conditions simultaneously (up to 4.16 million candidate sets per M, 22-element pools). Candidate counts for odd M (90, 1700, 1700, 28296, 31826)
coincide with R7's. Filter validated: with M put back into the pool it returns exactly the admissible families of Section 2.
So "no support-admissible N with C < phi(M)" also holds for these even M.

## 5. Cube minus top [N]
R7's slacks reproduced (y = 11, 47, 200, 1000, 10^4: -0.0617, -0.0355, -0.0175, -0.0044, +0.0003; closed form checked by brute force at M = 1155).
Correction: the per-prime cuts are all satisfied already from y between 2000 and 2500 (theta <= 0.147), not "about 10^4";
at y = 10^4 the smallest slack sits at the largest prime, not at q = 3. Vertex and singleton slacks not recomputed.

## 6. Does "feasible in y" give the hypotheses of Prop 4.6? [GAP]
* Not literally. (a) The "1 not required" searches (N = {M}, {105, 1155}, ...) are outside Def 2.6: no Z^a padding, only d = 0 mod C.
  (b) A fixed N never satisfies C >= log d. In the paper C >= log d is used only for (28) and r >> (log d)^4, both trivial for bounded C,
  so the conclusion very likely survives for fixed N (my reading of Sec 2.8 and 4.6; Sections 4.2-4.4 not re-read) [H]; otherwise the small sets
  matter only as coprime blocks in Lemma 3.3. (c) I/C << log C and |y_1| << log C are trivially true for fixed N.
* Consequence: the small-case numbers are evidence about blocks, not about growing non-product families.

## 7. What survives
* Survive: all definitions; reformulation; LP criterion; windows; Theorems A (odd M) and B; all exhaustive-search numbers re-tested; cube-minus-top numbers.
* Do not survive as stated / weakened: "Theorem A corollary" for even M; "each family fails at leading order" (family 1 is exactly borderline at
  leading order, fails by lower-order terms; families 2-4 unverified heuristics); "non-cyclic abelian groups do not help" (the argument given,
  faithful Q-irreducible => cyclic, does not apply to sums of several components, which is what GGS uses; not analysed).
* Still open: M/C <= K M/phi(M) for all admissible data. Evidence only covers moduli with at most 24 divisors; in the relevant regime
  (many primes, sum x_hi >> 1) nothing is proved beyond Theorem B (excess up to y^2).
