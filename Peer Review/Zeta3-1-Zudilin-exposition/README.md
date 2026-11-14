# ζ(3), proof 1: Zudilin's proof, expanded and refereed

**Document:** `exposition.md`.
**Claim:** a complete, checked exposition of W. Zudilin, *An elementary proof of Apéry's theorem*
(arXiv:math/0202159): every lemma with omitted calculations filled in, the purpose of each step, a
referee report on the paper (§8), where the recurrence comes from (§7), and the dictionary with
Beukers' integral proof (§9).

**What to check first**
1. §3.2 Steps 4–5: the maximum of $f$ and the exact value $e^{f(\tau_0)}=(\sqrt2-1)^4$ (proved by hand here;
   asserted without proof in the paper).
2. §5.1: the claim that the paper's "trivial bound $D_n<3^n$" is not trivial, and that the elementary
   period-6 bound closes the gap.
3. §8: the two misprints attributed to the paper.
4. §9: the Mellin-transform dictionary (our derivation; Beukers' 1979 paper itself was not retrieved).

**Depends on:** nothing outside the document except the cited prime estimate.
**Scripts:** `zudilin_lemma1.py`, `zudilin_lemma2_recurrence.py`, `zudilin_ball_series.py`, `dn_bound_check.py`.
Two FAIL lines in `output_zudilin_ball_series.txt` are errors in the brief given to the script writer
(a mis-transcribed factor and a constant stated without its scaling), not in the mathematics.
