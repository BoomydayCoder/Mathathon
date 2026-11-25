# Referee report on R15 `proof.md` (Theorem R: (H-Rogers)_j from (L1))

**Verdict: CORRECT_MODULO_MINOR_FIXES.** I tried to refute the proof of the stated implication (sections 0-3) line by line and could not. No fatal error found. The theorem is conditional exactly as stated: on (L1)_m in ranks k, ..., k-j+1 (cited from Gargava, NOT verified here or there), on (h1^l) (LEFT ideals principal), and on classical Haar-measure facts. The minor fixes are added sentences and two amendments to the wording of R14, not changes of mathematics.

Not done: no numerical experiment (no Monte Carlo over SL_3(Z)\SL_3(R)); Gargava's paper was not opened (no internet); R14 was read only in its "Main results" block and Section 1 (Lemma 1.1, 1.2, Remark 1.3).

## 1. What is actually proved

Theorem R as stated: for a Q-division algebra Delta, an order O with every nonzero LEFT ideal principal, a family of closed subgroups G_m of GL_m(Delta_R) satisfying (G1) (normal in GL_m(Delta_R), |real det| = 1) and (G2) (stabiliser of e_1 is G_{m-1} semidirect Delta_R^{m-1}), Gamma_m = G_m cap GL_m(O), and 1 <= j < k: IF for each m in {k, ..., k-j+1} Gamma_m is a lattice in G_m and the rank-one mean value identity holds on Gamma_m\G_m for f in C_c, THEN for every Borel F >= 0 the mean over Gamma_k\G_k of the sum of F(Y gamma) over Y in O^{j x k} of Delta-rank j equals the Lebesgue integral of F (covol(O^{j x k}) = 1), as an identity in [0, infinity]. Constant exactly 1.

It does NOT prove (L1), and it does not prove (h1^l) for any particular order.

## 2. Checks, by the points of the task

**(1) Sides, "primitive".** Rows, scalars on the left, matrices on the right; the two actions commute, which is what makes Stab(a e_1) = Stab(e_1) for a in O\{0} (a is invertible in Delta_R) and rk(a e_1; Y') = rk(e_1; Y'). Verified. The lattice O^{1xk} gamma depends on Gamma_k gamma, so the space is Gamma_k\G_k; R14 writes "gamma Gamma_k in G_k/Gamma_k", on which the summand is not even well defined. This is a slip in R14, correctly diagnosed in 0.4; the transposition dictionary to Gargava's {g O'^k} for (Delta^op, O^op) is correct ((y gamma)^T = gamma^T ._op y^T, anti-isomorphism, same real-linear map). Lemma 2.1 (y = a e_1 delta, delta in GL_k(O)): I redid it. I = {c : cy integral} is a f.g. left O-submodule of Delta, tI = O c_0 by (h1^l), saturation M = O y', quotient Z-torsion-free hence free by the mirror image of R14 Lemma 1.1 (whose proof I read; it is symmetric), split, delta integral with integral inverse. Correct; it genuinely needs LEFT principality.

**(2) Transitivity.** Not used. The proof decomposes O^{1xk}\{0} into countably many Gamma_k-orbits and treats each one separately; delta in GL_k(O) need not lie in G_k, and this is handled by normality (G1): H_y = delta^{-1} H delta, Gamma_y = delta^{-1} Gamma_H delta because delta normalises both G_k and GL_k(O). Verified.

**(3) Unfolding and disintegration.** Lemma 2.6 is the standard unfolding plus Weil's quotient formula; I rechecked both equalities (fundamental domain bigsqcup delta F for Gamma'; inner integral over H' collapses to a fundamental domain D because exactly one epsilon in Gamma' moves hg into E; the inner integral is left-H'-invariant by right invariance of eta', H' unimodular). Lemma 2.2: u(n,h)u(n',h') = u(n + h n', h h'); dn dh is left Haar because ||h|| = 1 on COLUMNS (Lemma 0.2, checked on elementary and diagonal generators: both determinants are prod |c_i|^{dim K}) and right Haar because G_{k-1} is unimodular; Gamma_H = Gamma_{k-1} semidirect O^{k-1}; fundamental domain F' x P of mass 1. So the fibre Gamma_H\H is a torus bundle over Gamma_{k-1}\G_{k-1} with total mass 1 and the base carries exactly the invariant PROBABILITY measure of the rank k-1 family. For s >= 2: for family (A) (reduced norm one at every place) the stabiliser's Levi factor is prod_v SL_{k-1}(K_v); for family (B) ({||.|| = 1}, the group {N = 1} of R14/Gargava when no real place exists) it is {||h|| = 1} in rank k-1, because det_R of x = (x_1 | x'') -> (x_1 + x'' n | x'' h) is det_R(R_h). Both families are coherent; no mixing of the two occurs. Verified.

**(4) Constant, circularity.** Torus lemma 2.3 is correct (full row rank => right inverse over Delta => surjection of tori => push-forward of Haar is Haar). Fibre mean value = integral of Psi for EVERY orbit (Lemma 2.4/2.5; the factor a never enters). Then (3.1): the j-sum equals sum_i int_{H_i\G_k} f(y_i g) d nu_i with f the partial integral of F, and the j = 1 unfolding of f with the same (G_k, Gamma_k, H_i, eta_i), hence the same nu_i, gives the same expression, which (L1)_k evaluates as int f = int F. Not circular: (L1)_k is used only for the first row; (R_{j-1,k-1}) uses ranks k-1, ..., k-j+1. Constant exactly 1.

**(5) Tonelli / measurability.** All integrands are countable sums of nonnegative Borel functions; f and Psi_1 are Borel by Tonelli; identities in [0, infinity]. Lemma 0.5 (C_c to Borel via Riesz) is correct.

**(6) Match with R14.** Same group when no K_v = R ({N = 1} = {||.|| = 1} = G^+ = G^pm, coherent), same lattices O^{1xk} gamma, same normalisation covol(O^{j x k}) = 1, same range k > j; R14 Lemma 1.2 uses j <= r < k, so all ranks used are >= k-r+1 >= 2. Two amendments to R14 are needed: (a) read G_k/Gamma_k as Gamma_k\G_k; (b) R14's (h1) is for RIGHT ideals, Theorem R needs LEFT ideals. Lemma 4.1a (B-adjoint anti-involution; O(Lambda^#) = O(Lambda)^*) is correct, I rechecked it: for a B-unimodular Lambda_W with rational G-invariant B (Leech with its standard form, G inside Co_0) the involution preserves O and (h1) <=> (h1^l). For Theorem E of R14 in general, (h1^l) is an extra hypothesis.

**(7) Delta = Q, O = Z, j = 2, k >= 3, redone independently.** Orbits of SL_k(Z): y_a = a e_1, a >= 1. Stabiliser H = SL_{k-1}(R) semidirect R^{k-1} for every a. Fibre: sum over (b | Y''), Y'' != 0, of Psi(b + Y''.n | Y'' h); integrating n over [0,1)^{k-1} turns sum_b into int_R dt (Y'' != 0 integral); Siegel in rank k-1 >= 2 gives int Psi. So the fibre integral is f(a e_1 g), f(x) = int F(x, Z) dZ. H\SL_k(R) = R^k \ {0} with nu = c dx; Siegel's primitive-vector formula gives c = 1/zeta(k). Total: sum_a c a^{-k} int f = int f = iint F. This is the leading term of Rogers' formula for independent pairs, valid for j <= k-1. For j = k the induction would need rank 1, where the mean value identity is false, consistent with the classical divergence.

## 3. Fatal errors

None found.

## 4. Minor issues and fixes

1. Lemma 0.5: say explicitly that lambda is finite on compacta ((L1) applied to a C_c majorant of 1_K); otherwise Riesz uniqueness is not available. One sentence.
2. Lemma 2.6: state that Weil's formula is extended from C_c to Borel phi >= 0 by monotone classes (lcsc), and that the inner integral is a Borel function of H'g.
3. Lemma 2.1: "torsion-free" means Z-torsion-free (O-lattice); this is how R14 Lemma 1.1 uses it.
4. R14's statement of (H-Rogers) must be amended: Gamma_k\G_k instead of G_k/Gamma_k, and hypothesis (h1^l) added (automatic for Leech by Lemma 4.1a). Until then "(H-Rogers)_j of R14 is proved" is true only for the amended statement.
5. Lemma 4.2a (only when a real place exists; irrelevant for K = H): correct, but add why the stabiliser is unimodular (L semidirect V with ||h|| = 1 and L a closed normal subgroup of the unimodular group {||h|| = 1}), and say explicitly that it is applied to G^+_m inside G^pm_m (index <= 2).
6. The STATUS line should say "PROVED conditionally on (L1)" rather than "PROVED": (L1), including the lattice property and the caveat of R6 Rem. 4.1a(iv) on Gargava's proof for s >= 2 (which is the Leech case, Delta_R = H x H), is an unverified citation, and Gargava's definition of the norm-one group was not re-read. Theorem R holds for either reading, since both families are coherent.

## 5. Residual risk (outside the proof)

Everything now rests on (L1) in ranks k-r+1, ..., k for the group {N = 1} with Delta_R = H x H, quoted from Gargava and never re-derived in this project. Theorem R transfers that risk unchanged; it adds none.
