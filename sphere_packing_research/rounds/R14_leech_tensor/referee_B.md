# 14회차 심사 보고 (referee2:algebra) — 하위 에이전트 보고 원문

(R14 증명 `proof.md` 에 대한 반박 담당 심사 에이전트의 구조화된 보고를 그대로 옮긴 것이다.)

## overall

Referee 2 (algebra, lattices, mean value formula); about 10 of 14 minutes used (09:36:38 to about 09:47 UTC). Object: /home/user/Mathathon/sphere_packing_research/rounds/R14_leech_tensor/proof.md, read in full, together with /home/user/Mathathon/sphere_packing_research/rounds/R5_leech13/notes.md and Section 4 of /home/user/Mathathon/sphere_packing_research/rounds/R6_quaternionic/proof_revised.md.

I found no mathematical error in the parts in my remit. Theorems D and E are correctly stated as CONDITIONAL theorems, and the hypothesis list is honest and complete as far as I can see. The constant 91/(12e) = 2.78975 is correctly derived from the stated inputs.

What I re-derived by hand and found correct:
- Set-up: left/right conventions, Lambda_W (x)_O L = M O^{r x k} gamma (column times row is O-balanced), G acting by unitary u_g and preserving every Lambda(L).
- Lemma 0.1 (Hermitian structure), Lemma 0.2 (invariant forms = Herm_k(K)^s, D_k = 2(2k^2-k) for Leech, trace form d r Re tr).
- Lemma 0.3: P(x(x)x) = R_K with K = x^*x/(dr); q = (1/(dr)) sum lambda^2/(sum lambda)^2 >= q_0 = 1/(r^2 s d) = 1/72.
- Lemma 0.4 (q >= (r/j) q_0 on E_U) and Lemma 0.5 (stabilisers).
- Lemma 1.1 (class number one implies freeness).
- Lemma 1.2: the parametrisation X = B_U Y of matrices with exact column space U is a disjoint decomposition; the Jacobian gives mu_U = kappa_n^{-f} (c_W^f/c_U)^k vol_{E_U}. This equals Vol_{fn}(B cap E_U)/covol(Lambda_U) in the covolume-1 normalisation, i.e. Vol/det^{k/2} for Leech.
- Lemma 1.4: T_U(V) <= 2 f^{-1/2} (pi n)^{-(1-f)/2} V^f R_U^{k/2}; the Stirling ratio was checked numerically (0.998 to 0.99999 for n = 240 to 24000).
- Lemma 1.5 (polynomial count of subspaces by Minkowski II; dyadic summation valid for k >= N_W^2 + 2) and the arithmetic of Remark 6.1 (3^24 X^6 lines; 2*3^36 (3/4)^{2k-6} for k >= 10; planes (9/16)^{2k}).
- Lemma 2.2 (orbit counting with non-free points at weight 1).
- The algebra of Lemma 3.4, including log t <= (t-1) - (t-1)^2/(2sj) on (0, sj] and the Gamma-function bounds.
- Cor. 3.6 exponent -5(1-f)/2; the bookkeeping of Prop. 4.1 and Section 5 (theta = 2 - 1/k).
- Section 6: (BS) for Leech follows from R5 Theorem 3. R_U = 3^8/det for lines; R_U <= 3^8/2^16 = 0.10011 for planes via det(L cap U) = det(L cap U^perp), where U^perp is a Delta-line because x^* = xbar. The number 1638 is not needed, only finiteness.

Not checked by me (outside my remit): the [V] steps taken from Klartag/R6 and the analytic uniformity in Prop. 3.5.

Main finding: the 'one missing lemma' (H-Rogers)_j for j >= 2 appears to be derivable, with constant exactly 1, from Gargava's Theorem 16 in ranks k-j+1,...,k (all >= 2 when j <= k-1) plus class number one. The route is the classical Siegel-Rogers primitive-frame unfolding; the sketch is in the first issue below. I derived it under time pressure and it is unrefereed. Until it is written out and checked, Theorems D and E must stay 'conditional on (H-Rogers)'. The literature agent found no division-algebra Rogers formula to cite, so a citation tag is not available. The statement of (H-Rogers) itself is correct and is also what Weil's Tamagawa-number argument would give.

Two inputs can be upgraded:
- dim End_G(W) = 8 follows EXACTLY from R5's integer trace table: sum_g tr(g)^2 = 17472 = 8 * 2184 (computed here). Hence End_G(W) = Delta, W is irreducible and dim_Q A = 72 without floating-point ranks.
- Class number one: I recomputed zeta_F(-1) = 1/6 for F = Q(sqrt13) (B_{2,chi} = 4), so the Eichler mass is 1/12. Since N(eps) = -1, nrd(O^x) lies in the squares of units, so [O^x : O_F^x] = |O^1|/2 = 12 given O^1 = SL_2(3) (R5 item 1). The mass is exhausted by the trivial class, so h(O) = 1, conditional only on the Eichler mass formula as a cited result.

## verdict

correct_modulo_minor_fixes

## issues

- **location**: Hypothesis (H-Rogers)_j, Remark 1.3, Remark 5.1, Section 7 item 1 (the main gap)
- **severity**: major
- **description**: Not an error, but the status is probably too pessimistic: (H-Rogers)_j for 2 <= j <= k-1 seems to follow with constant exactly 1 from (L1) in ranks m = k-j+1,...,k (each >= 2) and class number one. Referee's sketch (unrefereed, done in a few minutes):
(a) Primitive vectors. Every nonzero v in O^{1 x m} is v = c b with b primitive (the saturation Delta v cap O^{1 x m} is free of rank 1) and c in (O minus 0)/O^x. The measure F -> E sum_{b prim} F(b gamma) is G_m-invariant, dominated by Gargava's measure (so Radon), and carried by the single open orbit {y : y_v != 0 for all v}; hence it equals p_{1,m} Leb. Gargava then gives 1 = p_{1,m} sum_c [O : cO]^{-m}, i.e. p_{1,m} = 1/zeta_O(m), with zeta_O(s) = sum over right ideals of [O:I]^{-s}.
(b) Primitive j-frames. Unfold over the first row. Gamma_k is transitive on primitive vectors, with stabiliser Gamma_{k-1} semidirect O^{k-1} inside H_1 = G_{k-1} semidirect Delta_R^{k-1} (unimodular). Frames with first row e_1 are (e_1; (b, B')) with B' a primitive (j-1)-frame in O^{1 x (k-1)} and b in O^{j-1} arbitrary. Integrating over the torus Delta_R^{k-1}/O^{k-1} turns the sum over b into a Lebesgue integral, because a -> B'a is a surjective homomorphism of tori (B' has a right inverse over O). The induction hypothesis in rank k-1 then gives p_{j,k} = p_{1,k} p_{j-1,k-1} = prod_{i<j} 1/zeta_O(k-i).
(c) All rank-j matrices. Y = C B with C in (M_j(O) cap GL_j(Delta))/GL_j(O), corresponding to the finite-index right submodules C O^j of O^j, and int F(Cy) dy = [O^j : C O^j]^{-k} int F. The Hermite normal form (principal diagonal right ideals c_l O, off-diagonal entries modulo c_l O) gives sum_M [O^j:M]^{-k} = prod_{l<j} zeta_O(k-l). The product with (b) is 1.
This is the classical proof of Rogers' formula from Siegel's. It uses only Tonelli, unimodularity, and Gargava's theorem as a black box in each rank.
- **suggested_fix**: Write (a)-(c) out as a lemma '(H-Rogers)_j from (L1) in ranks k-j+1..k and (h1)' and have it refereed. Points needing care: the exact definition of N, G and Gamma in Gargava (see the next issue); transitivity of Gamma_k on primitive vectors; the measure normalisations in the two unfoldings (the same vol(Gamma_{H_1}\H_1) and the same measure on H_1\G appear in both, so no absolute volume is needed); left class number one for the row modules. Until this is refereed, keep (H-Rogers) as an explicit hypothesis. Do not cite it as literature: no division-algebra Rogers formula was found by the literature agent.

- **location**: (H-Rogers) statement and Lemma 1.2: the groups G_k = {N = 1}, Gamma_k, and the quotient
- **severity**: minor
- **description**: Two points.
(i) Side slip. The lattice O^{1 x k} gamma and the sum over Y gamma are invariant under gamma -> delta gamma with delta in Gamma_k, so the parameter space is Gamma_k\G_k with the right-invariant probability measure, not 'G_k/Gamma_k, gamma Gamma_k'. It is consistent after transposing to Gargava's column convention for (Delta^op, O^op).
(ii) The meaning of N is not pinned down. For Delta_R = H x H the group {N = 1} for the absolute norm is SL_k(H)^2 times a one-parameter central torus (t, 1/t), made compact by the unit eps; it is not SL_k(Delta_R). The proof only uses invariance and preservation of Lebesgue measure, so both conventions work. The identification with Gargava's theorem, and any derivation of (H-Rogers), must nevertheless use one convention throughout, including in the stabiliser G_{k-1}.
- **suggested_fix**: State N explicitly (absolute value of the real determinant of right multiplication on Delta_R^{1 x k}), write Gamma_k\G_k, and note that Gamma_k = GL_k(O) in the quaternionic case (the real determinant of an element of GL_k(O) is +1). This also gives transitivity of Gamma_k on primitive vectors.

- **location**: Standing data (T), hypothesis (h1); Section 6 '(Left versus right ... convention)'
- **severity**: minor
- **description**: (h1) is stated for right ideals of O only. That is exactly what Lemma 1.1/1.2 need (O^r cap U is a right module). The random lattices L = O^{1 x k} gamma are left modules, and any proof of (H-Rogers) (saturations of row submodules) needs left ideals principal as well. This is automatic: if all right ideals are principal then O is hereditary, and a left ideal J equals J^{**} = O c^{-1}, where {x : Jx in O} = cO. For a maximal quaternion order it follows from conjugation. It should be said. For Leech, R5's O acts on the left of Lambda_24, so the template's O is O_R5^op; class number one of a maximal quaternion order is left-right symmetric, so nothing breaks.
- **suggested_fix**: Add one line: (h1) implies the same for left ideals (duality argument or conjugation), and for Leech O_template = O_R5^op, which is isomorphic to O_R5 via conjugation.

- **location**: Section 6 (R5-a): 'numerical rank 8', 'numerical rank 72', 'class number one (from memory)'
- **severity**: minor
- **description**: These three inputs can be made exact or better founded, and the dependency list overstates the reliance on floating point.
(1) From R5's integer class/trace table (orders 1, 2, 3, 4, 6, 7, 12, 13, 14, 26 with traces 24, -24, 0, 0, 0, -4, 0, -2, 4, 2): sum_g tr(g)^2 = 17472 = 8 * 2184, so dim End_G(W) = 8 exactly. As End_G(W) contains the 8-dimensional division algebra Delta, it equals Delta; W is irreducible and dim_Q A = 72 by the double centraliser theorem.
(2) Class number: zeta_F(-1) = 1/6 for Q(sqrt13) (B_{2,chi_13} = 4, recomputed here), so the Eichler mass is 1/12. Because N(eps) = -1, all totally positive units are squares, so O^x = O_F^x O^1 and [O^x : O_F^x] = |O^1|/2 = 12 given O^1 = SL_2(3) (R5 item 1: every element of O^1 has order dividing 12 and commutes with g, so it lies in the order-24 group found there). The mass is exhausted by one class, so h = 1. The Eichler mass formula itself remains a cited result.
- **suggested_fix**: Replace 'numerical rank' by the exact character computation and record the mass-formula computation explicitly, with the mass formula tagged as cited literature.

- **location**: Theorem E statement / Standing data (T): 'Delta_R = K^s'
- **severity**: minor
- **description**: The template silently excludes all Delta whose real completion is not a product of copies of ONE division algebra K (indefinite quaternion algebras, mixed real/complex places, any place where Delta_R has a matrix factor). It is stated in (T) but not in the theorem. Lemma 0.1-0.4, D_k and q_0 = 1/(r^2 s d) all depend on it.
- **suggested_fix**: Repeat the restriction in the statement of Theorem E. Leech (H x H) is unaffected.

- **location**: Lemma 1.5(b) and Prop. 4.1 (summability of sum_U R_U^{k/2})
- **severity**: minor
- **description**: Correct, but the general lemma needs k >= N_W^2 + 2 (= 578 for Leech). Only Remark 6.1, which uses R5's finer statement det = N_U^4, gives the range k >= 10. I checked the arithmetic of Remark 6.1 and it is right. The range of k is irrelevant for the asymptotic theorem but should be visible. Also, the equality-case description of R5 Theorem 3 (exactly 1638 lines) depends on class number one, whereas the inequality det >= 3^8 does not; the proof correctly uses only the inequality plus finiteness.
- **suggested_fix**: State the range of k in Prop. 4.1 and Theorem E ('k >= k_0(data)'); keep the remark that 1638 is not used.

- **location**: Remark 1.3 (expected proof via Weil) and Remark 5.1 (unfolding up to a constant)
- **severity**: cosmetic
- **description**: Both remarks are consistent with what I derived. In Remark 5.1 the unfolding constant is in fact computable: c_{j,k} = p_{j,k} zeta_{O^j}(k) = 1 by the argument in the first issue. The sentence 'I could not derive j >= 2 from Gargava's theorem' should be replaced once that lemma is written. The 'robust form' with 1/c_{r,k} is then superfluous but harmless.
- **suggested_fix**: After the new lemma is refereed, delete or shorten Remarks 1.3 and 5.1.

- **location**: Dependency statement of Theorem D (end of file)
- **severity**: cosmetic
- **description**: If (H-Rogers)_{2,3} is obtained by the unfolding route, Theorem D needs Gargava's Theorem 16 for (Delta^op, O^op) in ranks k-2, k-1 and k (all >= 2 for k >= 4), not only in rank k. The R6 caveat (transitivity only on the conull open orbit for Delta_R = H x H) then applies in each of these ranks.
- **suggested_fix**: List (L1) as 'for all ranks m >= 2'.


## what_is_actually_established

ESTABLISHED (conditional). Theorem E (template) and Theorem D (Leech tensor family, V_n >= (91/(12e) - o(1)) n^2 for 24 | n, n = 24k, k -> infinity) hold under exactly these hypotheses:

(H1) (H-Rogers)_j for 2 <= j <= r (Leech: j = 2, 3; k >= r+1): the Siegel-Weil-Rogers formula with constant 1 (only '<= (1+o(1)) integral' is used) for Delta-independent j-tuples of the random O-lattice O^{1 x k} gamma. The statement is correct as formulated. It is NOT available as a citable result for orders in division algebras (literature agent). It is very probably derivable from (H2) in ranks k-j+1..k plus class number one by the classical primitive-frame unfolding (referee's sketch in the first issue: p_{1,m} = 1/zeta_O(m), p_{j,k} = product of p_{1,k-i}, Hermite normal form gives the reciprocal product, so the constant is exactly 1). That sketch is unrefereed, so (H1) stays a hypothesis for now.

(H2) Gargava arXiv:2107.04844 Theorem 16 for (Delta^op, O^op) with Delta_R = H x H, rank >= 2 (this is j = 1), with the R6 caveat that the cited proof was not re-verified for s >= 2.

(H3) The Klartag/R6 process machinery with Q = invariant forms ([V] steps: Prop. 2.1, 3.1, 3.2, Lemma 3.3) and the analytic uniformity of Prop. 3.5. Not in my remit; I checked only that the non-free replacements (Lemma 2.2, Prop. 2.3, the non-free bound in Prop. 4.1 via Lemma 0.5) are algebraically correct.

(H4) R5 computer-assisted facts: O maximal, O^1 = SL_2(3), |G| = 2184 with the stated class/trace table, and det(Lambda_24 cap U) >= 3^8 for every Delta-line (hence for every Delta-plane by unimodularity). Given these, the following are exact and no longer depend on floating point or memory: dim End_G(W) = 8 (character sum 17472/2184, computed here), dim_Q A = 72, and class number one (Eichler mass 1/12 = 1/[O^x : O_F^x], zeta_F(-1) = 1/6 recomputed here; the Eichler mass formula is a cited result).

(H5) Standard tools: Minkowski's second theorem, Stirling, double centraliser.

Verified by me in full:
- conventions and Lambda_24 (x)_O L = M O^{3 x k} gamma;
- invariant forms = Herm_k(H)^2, D_k = 2(2k^2-k), q_0 = 1/72, the formulas for P(x(x)x) and q, and q >= (r/j) q_0 on rank-j points;
- freeness from class number one;
- the decomposition of dependent tuples by Delta-subspaces U and the measure mu_U (T_U = Vol/det^{k/2}; borderline line: sqrt3 pi^{-1/3} V^{1/3} n^{-1/3});
- the counting and summability over U (lines and planes);
- (BS) for Leech with R_U <= 3^8/2^16 for planes and R_U in {1} or <= 6561/6562 for lines;
- the algebra of the angular lemma 3.4;
- the constant w q_0/(4e) = 2184/(288e) = 91/(12e) and theta = 2 - 1/k.

NOT established: an unconditional Theorem D. The numbers 1638 and 91 (borderline and fixed lines) are not used by the proof.
