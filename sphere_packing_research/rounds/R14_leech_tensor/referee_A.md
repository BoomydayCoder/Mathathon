# 14회차 심사 보고 (referee1:analysis) — 하위 에이전트 보고 원문

(R14 증명 `proof.md` 에 대한 반박 담당 심사 에이전트의 구조화된 보고를 그대로 옮긴 것이다.)

## overall

Referee 1 (analysis). I read /home/user/Mathathon/sphere_packing_research/rounds/R14_leech_tensor/proof.md in full, the round-6 template sections 0-4, and R5 Theorem 3, and re-derived the key estimates by hand. I found NO error in the analytic part. The document is honest: Theorems D and E are stated as CONDITIONAL on (H-Rogers)_j, j=2..r, and that is exactly their status. The main term (full-rank points, j=r=3) depends on that hypothesis, so nothing unconditional about the constant 91/(12e) is established.

Checked and confirmed:
- **Setup (Lemmas 0.1-0.5).** The commutant is right multiplication by Herm_k(K) per place; Tr R_N = d r Re tr N; P(x(x)x) = R_{x^*x/(dr)}; q = (1/(dr)) sum lambda^2/(sum lambda)^2 lies in [q_0, r s q_0] with q_0 = 1/(r^2 s d) = 1/72. On E_U, q = (q_0/f)(1+CV_j^2). P(x(x)x) is constant on G-orbits because P is the G-average. Stabiliser points have all columns in W^g.
- **Mean value formula (Lemmas 1.1, 1.2, 1.4, 1.5).** Freeness from principal right ideals is correct. The decomposition X = B_U Y with rank Y = j is correct (left row rank = right column rank). The Jacobian lambda^{fn} c_U^k and mu_U are correct. Stirling gives T_U(V) <= 2 f^{-1/2}(pi n)^{-(1-f)/2} V^f R_U^{k/2}, using fn = k dim_Q U. The Minkowski-II count and dyadic summation are correct.
- **Orbit counting and volume (Lemma 2.2, Prop 2.3).** N_t >= D_k - #orbits >= D_k - |contacts|/w - |nf contacts| is correct linear algebra.
- **Round-6 [V] steps.** I checked in the R6 text that none uses freeness or q >= q_0. Prop 1.1 is stated for an arbitrary subspace Q containing Id and an arbitrary lattice. Lemmas 2.1-2.3 hold for an arbitrary adapted projection. Prop 3.1 uses only pi_t = pi_t P and the upper bound [M]_t <= t q(x)|x|^4, valid for every lattice point of any rank or stabiliser. Prop 3.2 is a pointwise sum, so restriction to a subset is legitimate. Freeness enters R6 only in Lemma 1.3, correctly replaced by Lemma 2.2.
- **Angular lemma (Lemma 3.4).** I re-derived it: log t <= (t-1) - (t-1)^2/(2sj), the Bartlett radii chi^2_{d(k-i+1)}, the Gamma-ratio bounds, and the final (1-kappa)^{-2bsj}. An independent Monte Carlo (quaternionic 3 x 12 blocks, s=2, b=1/2, 40000 samples) gave E prod t^{-2b} = 1.959 against the exact Gamma formula 1.957, and E exp(b CV^2) = 1.106 against the bound 5.62.
- **Shell integral (Prop 3.5, Cor 3.6).** Polar coordinates in dimension m = fn, the substitution sigma = q t, and all three hypotheses of Lemma 3.3 hold uniformly (m sqrt S >= sqrt(8fz), m^3 S^2 <= 64 r^2 s^2 z^2/n, Y sqrt S = C_0 sqrt(nT)). Exponents: m^2 q T/8 = f z (1+CV^2), 8f/(q_0 m^2) = 8/(f q_0 n^2), and the n-exponent is -5(1-f)/2.
- **Good lattice (Prop 4.1).** Three Markov events with total probability <= 1 - e^{-rho} + O(n^{-1/(2r)}) < 1; E X^{nf} <= (w-1) sum_U X_U.
- **Optimisation (Section 5).** theta = 4D_k/(q_0 n^2) = 2 + (4-2d)/(dk) = 2 - 1/k; the third term is 2/p; E log det <= -2 log(w q_0 n^2/(4e)) + o(1); the error terms vanish. The constant 2184/(288e) = 91/(12e) = 2.78975 was checked numerically; w D_k/(2e n^2) = 2.7884 at k = 1000.
- **Leech hypotheses (Section 6).** (BS) follows from R5 Theorem 3 plus integrality of determinants; R_U^{k/2} = (9/N_U)^{2k} for lines. Class number one is consistent with my own evaluation of the mass: zeta_F(-1) = 1/6 for Q(sqrt13) by Siegel's formula, so the mass is 1/12 = 1/[O^x : O_F^x], because totally positive units are squares (N(eps) = -1) and O^x = O_F^x O^1.

Remaining caveats are inherited inputs and purely asymptotic character (slow, non-effective o(1)); see the issues.

## verdict

correct_modulo_minor_fixes

## issues

- **location**: Hypothesis (H-Rogers)_j, j = 2..r (Lemma 1.2, Remark 1.3, Remark 5.1); used in Prop 4.1 and Section 5
- **severity**: major
- **description**: Disclosed by the prover, but it is the crux. The mean value formula for Delta-independent j-tuples is not proved and not cited from a seen source. It carries the MAIN term (j = r = 3, full-rank points) as well as the planes (j = 2). Without it nothing about the constant 91/(12e) is established; only the lines (j = 1) are covered by Gargava's Theorem 16. My own assessment, from memory and unverified: it should follow from Weil's theory (Tamagawa number 1 for SL_m(D) and for the stabiliser SL_{k-j}(D) semidirect D^{(k-j) x j}, a single rational orbit for j < k, strong approximation for SL_k(D), k >= 2). One point is not mentioned in the note: Weil's formula lives on prod_v SL_k(K)/SL_k(O), whereas the hypothesis is stated on Gargava's larger group {N = 1}.
- **suggested_fix**: Keep the theorem labelled conditional. To close it, either find a citable statement (Weil, Adeles and Algebraic Groups, Ch. III; or an adelic Rogers formula extended to division algebras) or prove it by Siegel-Rogers induction on j. Add the reduction from {N = 1} to prod_v SL_k(K): the quotient is a compact torus (scalings (c, c^{-1}) modulo units of reduced norm eps^2) acting by Lebesgue-measure-preserving maps, so the formula averages over it.

- **location**: Remark 5.1 (robust form with constants c_{j,k})
- **severity**: minor
- **description**: This is only a sketch, and the stated condition is incomplete. Besides sup_k c_{j,k} < infinity for j < r, the argument needs c_{r,k} = o(log n): p = 1 - (c_{r,k}+2)e^{-rho} must be positive and tend to 1, and z must stay in [log n, 3 log n]. Otherwise rho has to be re-chosen. The unfolding, local finiteness and the sum over infinitely many Gamma-orbits are not proved.
- **suggested_fix**: State the extra condition on c_{r,k}, or drop the remark from the list of results; keep it tagged as a heuristic sketch.

- **location**: Prop 4.1 (event B_1) and Section 7.6: size of the o(1) for Leech
- **severity**: minor
- **description**: The theorem is purely asymptotic and the threshold is very large. The initial-ball Markov step needs 1638 * 2 sqrt3 (pi n)^{-1/3} (log n)^{-1/3} < 1/log n, i.e. n of order 10^14 or more by my estimate. With the cruder exponent 1/(2r) = 1/6 written in Prop 4.1 it would be about 10^33. The dominant relative error is 1/p - 1 = 3/log n, plus (log n)^{-1/4} from Lemma 3.3 and exp(about 48 log^2 n / k) from Lemma 3.4. Nothing is claimed or provable from this text for moderate n.
- **suggested_fix**: State explicitly that Theorem D is a k -> infinity statement with non-effective, slowly decaying o(1). Optionally use the exponent (1-f)/2 = 1/3 for lines in Prop 4.1, and count G-orbits rather than points in the Markov step to gain a factor w.

- **location**: Cited input (L1): Gargava Thm 16 for (Delta^op, O^op), Delta_R = H x H
- **severity**: minor
- **description**: Inherited caveat from round 6. The cited proof asserts transitivity on D_R^k minus 0, which holds only on the conull open orbit when s >= 2, and the proof was not re-verified in this session. Passing to the opposite algebra is harmless. Theorem D depends on this for the 1638 borderline lines and for the non-free points.
- **suggested_fix**: Keep it listed as a cited input with the caveat. It would be subsumed by a proof of (H-Rogers)_j for all j < k through Weil's formula.

- **location**: Section 6: inputs from round 5 (R5 Theorem 3; class number one; ranks 8 and 72; |G| = 2184)
- **severity**: minor
- **description**: Not re-verified here. They are computer-assisted, with floating-point ranks, and the class number comes from Eichler's mass formula quoted from memory. R5 Theorem 3 reads as sound: Hermite's constant gives N(a) >= 4, and N(a) = 4 is excluded by the exhaustive computation over minimal vectors; planes follow by unimodularity. Class number one agrees with my independent evaluation (mass 1/12), but that too rests on the mass formula from memory. Note that R5 says det >= 3^8 does not need class number one, whereas R14 does need it (hypothesis (h1), for freeness in Lemmas 1.1 and 1.2).
- **suggested_fix**: List these as explicit hypotheses of Theorem D (already done). If possible, verify class number one by a direct computation, e.g. show that O is right-norm-Euclidean or enumerate ideal classes, and recompute the ranks 8 and 72 in exact arithmetic.

- **location**: Lemma 3.4 / Section 7.7
- **severity**: cosmetic
- **description**: The lemma is correct; I re-derived it and checked the Gamma identity by an independent Monte Carlo. The bound is lossy (exp(about 8 b^2 s j/(dk)) against a true value of about exp(2.4 b/k)), which only affects the speed of convergence. The reference to the un-refereed R6 Lemma 7.2 for the quaternionic determinant is unnecessary: prod mu_i = prod rho_i^2 follows from the real representation of the triangular factor.
- **suggested_fix**: Make the lemma self-contained by adding one line on det_R(T) = prod rho_i^4 for lower triangular quaternionic T with real diagonal.

- **location**: Prop 3.2 and Section 7.3 ([V] steps)
- **severity**: cosmetic
- **description**: The prover asked for a re-check that no round-6 [V] step uses freeness. Confirmed against the R6 text: Prop 1.1, Lemmas 2.1-2.3, Cor 2.2, Prop 3.1 and Prop 3.2 use neither freeness nor a lower bound on q. Restricting Klartag's Prop 4.2 to a subset of lattice points is valid, since it is a pointwise sum plus the same exceptional event.
- **suggested_fix**: None needed; the status of these steps can be upgraded from 'should be re-checked' to 'checked by referee 1'.


## what_is_actually_established

A CONDITIONAL theorem, not an unconditional one.

**Established (my line-by-line check plus the refereed round-6 machinery).** IF the mean value formula (H-Rogers)_j holds for j = 2 and j = 3 (constant 1, or <= 1+o(1); k > 3) for the Haar-random O-lattice of rank k over Delta = H (x) Q(sqrt13), THEN for 24 | n, n = 24k -> infinity, there is an O-lattice L with Lambda_24 (x)_O L having an empty ball of relative volume (1-o(1)) (91/(12e)) n^2. Hence V_n >= (91/(12e) - o(1)) n^2 along 24 | n. Theorem E is the same statement for general data (T), (h1), (BS), with (H-Rogers)_j for 2 <= j <= r.

**Hypotheses the proof really needs** (all are stated in the note; I found no silent one):
1. (H-Rogers)_2 and (H-Rogers)_3. The j = 3 case carries the main term. Unproved here; believed to be Weil's theorem, source not seen.
2. Gargava's Theorem 16 for (Delta^op, O^op) with Delta_R = H x H, with the round-6 caveat on its proof.
3. Klartag's Lemma 2.1, Lemma A.1 and the proofs of his Sections 2-4 as adapted in round 6. I verified that those adaptations use no freeness and no lower bound on q.
4. The round-5 computer-assisted facts: End_G(W) = Delta, span of G of dimension 72, |G| = 2184, O maximal, det(Lambda_24 cap U) >= 3^8 for all Delta-lines and Delta-planes.
5. Class number one of O (Eichler mass formula, from memory; consistent with my independent evaluation, mass 1/12).
6. Standard tools: Ito, Levy, Dambis-Dubins-Schwarz, reflection principle, GOE tail bound, Minkowski's second theorem, Stirling.

**Proved in the note and confirmed by me, given the mean value formula:**
- the lower-rank measures mu_U and the bound T_U(V) <= 2 f^{-1/2}(pi n)^{-(1-f)/2} V^f R_U^{k/2};
- summability under (BS);
- non-free orbit counting N_t >= D_k - |contacts|/w - |non-free contacts|;
- exact q = (q_0/f)(1+CV_j^2) on E_U;
- the angular lemma E exp(b CV_j^2) <= (1-kappa)^{-2bsj};
- the shell integral on every E_U and X_U <= C R_U^{k/2} n^{-5(1-f)/2} (so borderline lines cost O(n^{-5/3}) and no restriction f <= 1/2 is needed for fixed data);
- the initial-ball bound e^{-rho} + O(n^{-1/(2r)});
- theta = 2 - 1/k and the constant w q_0/(4e) = 2184/(288e) = 91/(12e).

**Not established:**
- (H-Rogers)_j for j >= 2.
- The robust form of Remark 5.1 (sketch only, and it also needs c_{r,k} = o(log n)).
- Any effective statement: the argument as written becomes non-vacuous only for n of order 10^14 or more, with relative error about 3/log n.
- Anything when r or G grows with k.
