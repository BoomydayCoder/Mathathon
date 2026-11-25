# R18. Integration referee report: R14 (Theorems D, E) + R15 (Theorem R)

(Saved by the orchestrator from the integration referee's final report, verbatim. The agent could not write the file
itself. Hard limit of the agent: 11 minutes.)

Read in this session: R14 `proof.md` (all), `referee_A.md`, `referee_B.md`; R15 `proof.md` (all), `referee.md`, `referee_B.md`; R5 `notes.md` (all); R6 `proof_revised.md` section 4 (Thm 4.1, Rem. 4.1a, Lemma 4.2) and section 6; R17 `notes.md` (check of Gargava Thm 16, all). No internet, no numerics. Not re-checked here: the analytic part of R14 (sections 2, 3, section 5 bookkeeping), the internal proof of Theorem R, Klartag's lemmas, the R5 computations.

## VERDICT: CONSISTENT_MODULO_MINOR_FIXES

Theorem D (Leech) does follow from (a) Gargava Thm 16 (+ Prop. 14) in ranks k, k-1, k-2 for (Delta^op, O^op), Delta_R = H x H; (b) Klartag's lemmas as adapted in R6; (c) the R5 facts including class number one; (d) standard tools. I found no gap at the interface. Every problem found is statement-level (wording, sides, hypothesis lists); none changes the mathematics.

## 1. Answers to the five questions

**Q1 (one probability space, one group family).** Yes. In every rank m take
G_m := { gamma in GL_m(Delta_R) : |det_R(y -> y gamma on Delta_R^{1 x m})| = 1 }, Gamma_m := G_m cap GL_m(O) = GL_m(O), X_m := Gamma_m\G_m with the right-invariant probability measure mu_m.
- Coherent (R15 0.3, family (B)): (G1) kernel of a character, norm 1 by definition; (G2) det_R of (x1 | x'') -> (x1 + x''n | x''h) equals det_R(R_h).
- For Leech (Delta_R = H x H) real determinants are positive, so G_m = {Nrd(g_1)Nrd(g_2) = 1}. This is exactly Gargava's total-norm-one group in every rank (R17 section 1: N(a) = det_R(a | D_R^m)^m, Gamma = GL_m(O)), transposed. So (L1)_m is literally [G, Thm 16 + Prop. 14] for every m >= 2 for the same family.
- Every use of the random gamma in R14 lives on (X_k, mu_k): Lemma 1.2 (one identity per subspace U, F(y) = phi(lambda M B_U y), j = dim_Delta U in {1..r}); Prop. 4.1 (B_1, B_2, B_3 are three Markov bounds on three first moments over the same space; X and X^nf are countable sums of nonnegative Borel functions of Gamma_k gamma; col(x) and Lambda^nf do not depend on the representative gamma).
- Section 5 uses a fixed lattice. Its last step multiplies by N = square root of a positive Hermitian matrix, N in GL_k(Delta_R), generally not in G_k. That is fine: the conclusion only asserts membership in {O^{1 x k} gamma' : gamma' in GL_k(Delta_R)} and uses the scale-invariant ratio Vol/covol; no mean value statement is applied to the new lattice.
- Measure preservation and covolume (R14 Lemma 1.2, last sentence): ||gamma|| = 1 on rows by definition of G_k, hence on E (r rows), so covol(Lambda) = lambda^n c_W^k = kappa_n. R14 cites [R6, Lemma 4.2] "from (L1)"; R15 has it as (G1)/Lemma 0.2. Both valid, no circularity.
- Statement match: R14 (H-Rogers)_j (proof.md lines 23-24) and R15 (R_{j,k}) are the same identity: same sum (Y in O^{j x k}, rk_Delta Y = j, F(Y gamma)), same normalisation covol(O^{j x k}) = 1, same F (Borel, values in [0, infinity], identity in [0, infinity]), same range 1 <= j < k, constant exactly 1. The only differences are R14's two slips: "G_k/Gamma_k" and "N" undefined.

**Q2 (sides).** Consistent.
- R14: W right Delta-space, Lambda_W = O^r columns (right O-module), L = O^{1 x k} gamma rows (left O-module). R15 uses the same (Delta, O), rows, scalars left, matrices right. R14 Lemma 1.1/1.2 need (h1) (right ideals: O^r cap U free); R15 Lemma 2.1 needs (h1^l) (left ideals: saturation of a row and its complement free). Both are needed in the merged theorem.
- Leech: R5's order acts on the LEFT, so Delta_template = Delta_R5^op and O_template = O_R5^op (same sets). O_R5 is a maximal quaternion order, so it contains O_F and is stable under x -> xbar = trd(x) - x; left and right class numbers agree, and h(O_R5) = 1 gives both (h1) and (h1^l) for O_template. Independent confirmation: R15 Lemma 4.1a (Lambda_24 unimodular for its standard rational G-invariant form, the form used in R14 section 6, c_W = 1). No new hypothesis for Theorem D.
- Gargava is applied to (Delta_template^op, O_template^op) = (Delta_R5, O_R5) literally: D = H tensor Q(sqrt13) inside End(W), O' = O(Lambda_24), D_R = H x H, column lattices g O'^k. Transposition (y gamma)^T = gamma^T ._op y^T maps rows O^{1 x k} gamma to columns gamma^T O'^k, Gamma_k\G_k to G^op/Gamma^op, and R_gamma on rows is the same real-linear map as L_{gamma^T} on columns, so norm conditions, Lebesgue normalisations and invariant measures correspond. Same family of lattices, correct algebra, correct order.

**Q3 (ranks).** No off-by-one. (R_{j,k}) needs (L1)_m for m = k, .., k-j+1, and k-j+1 >= 2 iff j <= k-1. For j = r = 3: ranks k, k-1, k-2, all >= 2 iff k >= 4. Chain at k = 4: (R_{3,4}) from (R_{2,3}) + (L1)_4 + [Gamma_3 lattice]; (R_{2,3}) from (R_{1,2}) + (L1)_3 + [Gamma_2 lattice]; (R_{1,2}) = (L1)_2. R15 Lemma 2.2's "k >= 3" holds at each step. Smallest k for the chain: k = r + 1 = 4. R14 Lemma 1.5(b) (k >= N_W^2 + 2 = 578 for Leech) and Remark 6.1 (k >= 10) do NOT constrain the chain: (MVF) is an identity in [0, infinity] for every k >= r+1; those two only give a k-uniform bound on Sigma_k. Theorem D itself holds for k >= k_0, non-effective (R14 referee A estimates n of order 10^14 or more).

**Q4 (does R14 need more than Theorem R?).** No.
- Tuples of rank < r are handled inside R14 by the bijection {X != 0} = disjoint union over U of {B_U Y : rk Y = dim U} (uses (h1)). Each piece is a full-rank j-frame sum with Borel F >= 0. No multiplicities, no lower-rank Rogers terms.
- Only first moments (three Markov inequalities). No second moment.
- The "fixed U" statement of Lemma 1.2 (used for Lambda^nf via Lemma 0.5) is Theorem R with j = dim U; nothing about a fixed sublattice.
- Integrability: all identities are in [0, infinity] (Tonelli over t, U, lattice points); test functions are bounded Borel; finiteness comes from R14 Lemma 1.4/1.5/Prop. 3.5.
- R14 Remarks 1.3 and 5.1 become superfluous (c_{j,k} = 1).

**Q5 (R14 referees' open issues not resolved by R15).**
- Resolved by R15: (H-Rogers) (A#1, B#1); side slip and meaning of N (B#2); left/right (h1) for Leech (B#3); list of ranks of (L1) (B#8).
- Not resolved by R15 but addressed in R17: the caveat on Gargava's proof for s >= 2 (A#4; R6 Rem. 4.1a(iv)). R17 repairs the false transitivity sentence in [G, Lemma 19] using the open orbit. After R17, (L1) for H x H rests on [G] as read plus Weil's reduction theory ([W58, Thm 2, Lemma 10.2]), cited in [G] without proof and read by nobody in the project. R17 has no referee report in its directory.
- Not resolved, a cited input: class number one via the Eichler mass formula (A#5, B#4). R14 referee B's computation (zeta_F(-1) = 1/6, mass 1/12 = 1/[O^x : O_F^x]) makes it exact given the mass formula. I re-added the character sum 576+576+7488+672+7488+672 = 17472 = 8 * 2184, so dim End_G(W) = 8 is exact.
- Presentation only: non-effective o(1) (A#3); Delta_R = K^s to be repeated in Theorem E (B#5); "k >= k_0" (B#6); Remark 5.1 (A#2); Lemma 3.4 (A#6).
- Nothing in either R14 report affects the truth of Theorem D beyond the cited inputs.

## 2. Interface problems (all minor), exact locations, fixes

1. **R14 proof.md line 23 ((H-Rogers)_j); inherited at line 66 (Lemma 1.2, "gamma random as in (H-Rogers)")**: "gamma Gamma_k on G_k/Gamma_k" is wrong, the summand is only left Gamma_k-invariant. Fix: Gamma_k\G_k with the right-invariant probability measure (R15 0.4). Not yet applied to R14's text.
2. **R14 line 23**: "N(gamma) = 1" with N undefined. Fix: |det_R(R_gamma on Delta_R^{1 x k})| = 1, Gamma_k = GL_k(O). For Leech this is Gargava's group (R17 sections 1, 4).
3. **R14 line 16 ((h1)) and the Theorem E/D statements (lines 26, 30)**: right ideals only. Fix: add (h1^l) to Theorem E; for Theorem D add one line (conjugation-stable maximal order, or R15 Lemma 4.1a). R15 Lemma 4.1a should state "B rational, Lambda_W B-unimodular" (R15 referee B, fix 4).
4. **R14 lines 30, 32, 157, 177 (dependencies)**: (L1) is needed in ranks k, .., k-r+1 (Leech: k, k-1, k-2), i.e. all m >= 2, including "Gamma_m is a lattice in G_m" ([G, Prop. 14]). R14 lists rank k only. Replace "(H-Rogers) [HYP]" by "Theorem R of R15".
5. **R15 proof.md status line, 0.4 and 4.5**: "(H-Rogers)_j of R14 is proved" holds only for the statement amended by items 1-3; the caveat of R6 Rem. 4.1a(iv) cited there should now point to R17.
6. **R15 section 4.2 / Lemma 4.2a against R14 line 23; only for Theorem E when some K_v = R**: Gargava's groups {N_m = 1} are then {det_R = 1} for odd m and {|det_R| = 1} for even m, which is not a coherent family as it stands. Fix: run both notes with the family G^pm above; (L1)_m for G^pm follows from Gargava's by Lemma 4.2a (index <= 2); R14 only needs (MVF) on some probability space of lattices of the family, so the j-frame transfer asked for by R15 referee B (fix 5) is unnecessary. Lemma 4.2a was checked by one referee only, and R17's note on sign slips in [G, Lemma 6/12] for real places and odd k is unresolved. Irrelevant for Leech and for K in {C, H}.
7. **R14 line 69 (Lemma 1.2, last sentence)**: replace "[R6, Lemma 4.2], from (L1)" by "definition of G_k (R15 Lemma 0.2)". Harmless either way.
8. **R14 line 32**: "(L1) checked in round 6, caveat applies" should be updated to R17.
9. **R14 Theorem E/D statements**: add "for k >= k_0(data), non-effective" and repeat "Delta_R = K^s" in Theorem E.

## 3. Consolidated statements

**Notation.** V(L) = sup Vol(B)/covol(L) over open origin-centred L-free balls; V_n = sup over lattices in R^n; kappa_m = Vol(B^m).

**Cited inputs.**
- **(L1)_m, m >= 2** [Gargava arXiv:2107.04844v3, Thm 16 and Prop. 14, for (Delta^op, O^op); statement read and proof repaired for s >= 2 in R17; depends on Weil's reduction theory, not read]: with G_m, Gamma_m = GL_m(O), X_m = Gamma_m\G_m as in Q1, Gamma_m is a lattice in G_m and for f in C_c(Delta_R^{1 x m}): the integral over X_m of sum_{0 != y in O^{1 x m}} f(y gamma) d mu_m equals the integral of f dx, covol(O^{1 x m}) = 1. (If some K_v = R: via R15 Lemma 4.2a.)
- **(L2)** Klartag arXiv:2504.05042 (Lemma 2.1, Lemma A.1, sections 2-4) as adapted in R6 `proof_revised.md` to a subspace Q containing Id (Prop. 1.1, Lemma 2.1, Cor. 2.2, Lemma 2.3, Prop. 3.1, Prop. 3.2, Lemma 3.3).
- **(L3)** Standard: Ito, Levy, Dambis-Dubins-Schwarz, reflection principle, GOE tail bound; Minkowski's second theorem; Stirling with remainder; double centraliser, Maschke; row rank = column rank and Gauss elimination over division rings; Weil's quotient integral formula, Borel fundamental domains, unimodularity of groups with lattices, push-forward of Haar measure under surjections of compact groups, Riesz representation.

**Theorem R (R15).** Delta a Q-division algebra with Delta_R = prod_v K_v; O an order with (h1^l): every nonzero left ideal principal; (G_m) a coherent family; 1 <= j < k; (L1)_m for m = k-j+1, .., k. Then for every Borel F: Delta_R^{j x k} -> [0, infinity], the integral over Gamma_k\G_k of sum_{Y in O^{j x k}, rk_Delta Y = j} F(Y gamma) d mu_k equals the integral of F dy (covol(O^{j x k}) = 1), in [0, infinity].

**Theorem E (template), consolidated.**
DATA (T): G a finite group; W an irreducible Q[G]-module, G in GL(W) faithful; Delta := End_G(W)^op (W a right Delta-space of dimension r); Delta_R isomorphic to K^s with one K in {R, C, H}, d = dim_R K; A = Q-span of G = End_Delta(W), dim_Q A = r^2 s d; Lambda_W in W a G-invariant lattice, W_R with a G-invariant Euclidean structure; O := {c in Delta : Lambda_W c in Lambda_W}; w = |G|, N_W = rsd, n = k N_W.
HYPOTHESES: (h1) every nonzero right ideal of O is principal; (h1^l) every nonzero left ideal of O is principal; (BS): with R_U := (dim W/dim U)^{dim_Q U} det(Lambda_W)^{f_U}/det(Lambda_W cap U), f_U = dim U/dim W, R_U <= 1 for every proper nonzero Delta-subspace U and R_U <= rho_0 < 1 for all but finitely many U.
INPUTS: (L1)_m for all m >= 2 for (Delta^op, O^op); (L2); (L3).
CONCLUSION: there is k_0 = k_0(data), non-effective, such that for k >= k_0 there is gamma' in GL_k(Delta_R) for which the G-invariant lattice Lambda_W tensor_O L, L = O^{1 x k} gamma' (realised as M O^{r x k} gamma' in Delta_R^{r x k} = R^n), has an empty open ball of relative volume (1 - o(1)) |G| n^2 / (4e dim_Q A), o(1) -> 0 as k -> infinity with the data fixed. Hence V_n >= (1 - o(1)) |G| n^2/(4e dim_Q A) for N_W | n. (For r = 1, (BS), (h1) and Theorem R are not needed; this is R6 Thm C'.)
Proof: R14 sections 0-5 with (H-Rogers)_j replaced by Theorem R for j = 2..r on Gamma_k\G_k (k >= r+1; Lemma 1.5 needs k >= N_W^2 + 2). (h1) and (h1^l) are equivalent if Lambda_W is unimodular for a rational G-invariant form (R15 Lemma 4.1a) or O is stable under an anti-automorphism of Delta. For K = R the theorem additionally uses R15 Lemma 4.2a and [G] with real places (see problem 6).

**Theorem D (Leech tensor family), consolidated.** Let Lambda_24 be the Leech lattice (standard unimodular form), g in Co_0 of order 13 and Delta = H tensor Q(sqrt13) in End(Lambda_24 tensor Q) as in R5, O = {x in Delta : x Lambda_24 in Lambda_24}, G = Aut_Delta(Lambda_24). ASSUME:
- (R5-a) [computer-assisted, R5] O is a maximal order; |G| = 2184 with R5's class/trace table (hence exactly: dim End_G(W) = 17472/2184 = 8, End_G(W) = Delta, W irreducible, dim_Q A = 72, r = 3, s = 2, d = 4, Delta_R = H x H); O^1 = SL_2(3);
- (R5-b) [computer-assisted, R5 Theorem 3] det(Lambda_24 cap U) >= 3^8 for every Delta-line U (hence for every Delta-plane, by unimodularity and x^* = xbar);
- (R5-c) [cited] h(O) = 1 (Eichler mass formula: mass = |zeta_F(-1)|/2 = 1/12 = 1/[O^x : O_F^x]); as O is stable under conjugation this gives (h1) and (h1^l) for O and O^op;
- (L1)_m for all m >= 2 for D = Delta (R5's algebra acting on the left, = Delta_template^op), order O, G_m = {(g_1, g_2) in GL_m(H)^2 : Nrd(g_1) Nrd(g_2) = 1}, Gamma_m = GL_m(O) [G Thm 16 + Prop. 14, repaired in R17; used in ranks k, k-1, k-2];
- (L2), (L3).
THEN there is k_0 (non-effective) such that for every k >= k_0, n = 24k, there is gamma' in GL_k(H)^2 for which the lattice Lambda_24 tensor_O L, L = O^{1 x k} gamma' (rows; a lattice in R^n invariant under SL_2(13) acting on the Lambda_24 factor), has an empty open ball of relative volume (1 - o(1)) (2184/(288e)) n^2 = (1 - o(1)) (91/(12e)) n^2. Hence V_n >= (91/(12e) - o(1)) n^2 for 24 | n, n -> infinity (91/(12e) = 2.7898).
Proof: Theorem E with (BS) from (R5-b): lines R_U = 3^8/det, in {1} or <= 6561/6562; planes R_U <= 3^8/2^16 < 0.1002; exceptional set finite (Lemma 1.5(a); the number 1638 is not used). The mean value chain is valid for k >= 4; Sigma_k is bounded for k >= 10 (Remark 6.1) or k >= 578 (Lemma 1.5); the theorem is purely asymptotic.

**Unverified in the whole chain (not interface issues).** Weil's reduction theorems behind [G, Thm 16]; the Eichler mass formula; the R5 computer runs; Klartag's paper. R17 is un-refereed.

Relevant files: rounds/R14_leech_tensor/proof.md, rounds/R15_rogers/proof.md, rounds/R17_gargava_check/notes.md, rounds/R5_leech13/notes.md, rounds/R6_quaternionic/proof_revised.md.
