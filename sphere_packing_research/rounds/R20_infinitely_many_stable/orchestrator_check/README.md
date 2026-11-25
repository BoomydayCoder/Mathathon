# R20 — orchestrator's spot check of the ideal-lattice model

What this is: a limited, independent check of the central claim of `../proof.md` (Theorem 1/2: the torus-fixed
sublattices of the SL_2(p) Weil lattice are explicit ideal lattices of subfields of Q(zeta_4p)).  It is **not** a
referee report on the proofs and **not** a replication of the t-circle scans.

Method (`chk.py`): only the lattice construction `build()` of `../ideal_model.py` is reused.  Everything else is
my own code, in exact arithmetic except for the enumeration tree (floating Gram-Schmidt with slack, every vector
re-evaluated exactly):

* exact tests on the integer Gram matrices `G` (lattice) and `E` (form twisted by eta): symmetry, evenness,
  `det G` (Bareiss), `eta = G^{-1}E` integral with `eta^2 = p` (exact rational solve);
* exact rational LLL, then enumeration of all vectors of norm 2 (roots), with exact `s = w.E.w/2` and
  `N(nu') = (p - s^2)/4` for each root.

Results (`out_p13_p37.txt`, `out_p61_p277_p661.txt`), rank-24 lattice `L_12(p)`, threshold `N(nu') >= 9`:

| p | even, det 1, eta^2 = p | roots (pairs) | (abs s, N(nu')) of the roots | agrees with |
|---|---|---|---|---|
| 13 | yes | 0 | — (so the lattice is the Leech lattice, by Niemeier's classification) | Rounds 5, 10 |
| 37 | yes | 72 | (3, 7): all 72 — **below the threshold 9** | Round 10 (the known violator, found there from the Weil representation) |
| 61 | yes | 360 | (1, 15): 288; (3, 13): 72 | Round 13 (E_8^3 type, abs s in {1,3}) |

(Further rows, if any, are in `out_p61_p277_p661.txt`.)

What was checked by hand (reading `../proof.md`): the algebra of Theorems 4 and 5 *given* the model
(ideal products, the identity `T_+- = 4 Tr_{K^+/F}(x_+)`, the AM-GM steps, the thresholds (A') and (A'')), and the
arithmetic of Proposition 5.2(ii).  Not checked: the proof of Theorem 1 itself (local analysis at p, the cited facts),
Corollary 4.2, Proposition 4.3, the exhaustive t-circle scans and their completeness.

Remark on wording: Section 5 of `../proof.md` shows that the Hermite-type route from (A) to (B)/(C) fails for
p > 73.  That is a statement about this method of proof; it is not a proof that (A) does not imply (B) or (C).
