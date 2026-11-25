# R19 referee report (adversarial, 10 minute limit)

Verdict: **CORRECT_MODULO_MINOR_FIXES**. No fatal error and no gap found in Theorem A, Lemma A.6, Theorem B, Lemma 1.2'. Theorem E' / Corollary D' are correct *relative to* R14 sections 2-5 and the cited inputs; I only checked that the substitution is legitimate, not R14's analysis.

## What I checked (by hand, line by line; no machine numerics, time limit)

### Part A
- A.1: for c in a^{-1}L, delta = 1 mod a^2 M_m(O): c(delta-1) in aL. Cosets c+aL are Gamma''-stable. Gamma'' has finite index in Gamma' and in GL_m(O); Gamma' discrete, lattice. OK.
- A.3: nu_c is a Borel measure, right G_m-invariant (right invariance of mu''). Sum over all a^{2Dm} cosets = Siegel transform of a^{-1}L, which is GL_m(O)-invariant, so A.2 + (L1)_m (statement only, Borel sets) give sum_c nu_c = a^{Dm} Leb exactly. Hence nu_c <= a^{Dm} Leb: local finiteness comes from the black-box statement of (L1)_m alone; no new domination, Gargava's Lemma 18 is not re-used. Support on Omega_m: complement is Lebesgue-null. (U) on Omega_m = H_y\G_m (G_m lcsc, transitive for m >= 2, Lebesgue is an invariant Radon measure) gives nu_c = C_c Leb. OK.
- A.4: scaling identity eps^{Dm} int f_eps d nu_c = C_c int f; pointwise Riemann-sum limit a^{-Dm} int f for each gamma (covol(aL gamma) = a^{Dm} because ||gamma|| = 1); Fatou on the probability space gives C_c >= a^{-Dm} for every coset; exact total sum_c C_c = a^{Dm} with a^{2Dm} terms forces C_c = a^{-Dm}. Rigorous and non-circular (lower bound by Fatou, upper bound by exhaustion of the exact total).
- A.5: [N:aL] a^{-Dm} = covol(N)^{-1}. Final constant is covol(N)^{-1}. OK.
- Hand check (1,0)+2Z^2: (1/3) zeta(2)^{-1} * (3/4) zeta(2) = 1/4. OK.
- A.6: both inclusions up to finite index verified (a L_0 <= L <= a^{-1} L_0 follows from aL <= L_0 <= a^{-1}L).

### Part B
- B.1: I_1 = pi_1(M) is a right O-lattice, projective by (Her); section is x -> sx, s = (1,s'')^T (right Delta-linear); A unipotent, AM = e_1 I_1 + (0 x M_1); left multiplication by A commutes with G_k, preserves rank and Lebesgue measure and covolume. M L_0 = N_1 x M_1 L_0 (generators argument correct). N_1 = I_1 L_0 and M_1 L_0 are Gamma-stable (Gamma acts on the right, M on the left).
- B.2: orbit decomposition and unfolding correct; Phi_i is Gamma_i-invariant; Tonelli in [0, infinity].
- B.3: G_k normal in GL_k(Delta_R), so g_i H_i g_i^{-1} = Stab_{G_k}(e_1) = H for any g_i in GL_k(Delta); g_i Gamma_i g_i^{-1} = Stab_H(L_i). Left-module splitting L_i = (J e_1 + (0|L_i'')) u(t,1) correct (section x'' -> (x''t | x''), left Delta-linear); replacing g_i by u(t,1) g_i keeps e_1 g_i = y_i. Stab_H(L_i) = {u(n,h): h in Stab(L_i''), n in T_i}: verified in both directions using u^{-1} = u(-h^{-1}n, h^{-1}). T_i is a lattice, h T_i = T_i. H unimodular with Haar dn dh (left: n -> m + sigma n has Jacobian ||sigma||_col = 1; right: translation in n, right translation in h). Fundamental domain F'' x P_i, normalised Haar covol(T_i)^{-1} dn dmu''.
- B.4: Y'' T_i <= M_1 J (well defined), surjective since rk Y'' = j-1; push-forward of normalised Haar is normalised Haar. OK.
- B.5: Y' = Y_1 g_i^{-1} in M_1 L_i = M_1 J x M_1 L_i''; rank condition becomes rk Y'' = j-1; induction hypothesis B_{j-1,k-1}(M_1, L_i'') with group Stab(L_i'') (allowed, induction over all pairs). Jacobians ||g_i||^{-(j-1)} on both sides cancel. The fibre mean is covol(M_1 L_0)^{-1} int Psi for every orbit.
- B.6: same nu_i (defined from mu and eta_i) in the j-unfolding and the 1-unfolding; Theorem A for (N_1, Gamma). Product of covolumes. OK.
- Classical test O = Z, k = 3, j = 2, M = Z e_1 + 2Z e_2: orbits y = d e_1, g_i = diag(d,1,1), L_i = d^{-1}Z x Z^2, J = d^{-1}Z, T_i = d^{-1}Z^2, M_1 J = (2/d)Z, M_1 L_i'' = 2Z^2. Fibre mean = (d/2)(1/4) d^{-1} int Psi = (1/8) int Psi, independent of d; Theorem A for N_1 = Z^3 gives factor 1; total 2^{-3} int F, as required. Second test M = 2Z e_1 + Z e_2: N_1 = 2Z^3, orbits 2d e_1, fibre mean 2d * 1 * (2d)^{-1} = 1, Theorem A gives 1/8. Also M = Z(1,1)^T + Z(0,2)^T reduces by A to the first test. All consistent.

### Part C
- Read R14 lines 38-78 (0.1, Lemmas 0.1-0.5, 1.1, 1.2, 1.4, 1.5) and grepped the whole of R14 for free/(h1)/class number/principal/O^r/B_U: outside 0.1, Lemma 0.5 (X in O^{r x k}), Lemma 1.1, Lemma 1.2 and the notation M O^{r x k} gamma (line 126), "free" only means free G-orbits or lattice-free balls. Lemmas 0.1-0.4 concern W_R, E, Q only; Lemma 1.4 and 1.5 use c_U, c_W and Lambda_W as a Z-lattice only.
- Lemma 1.2': {X in Lambda_W^k : col X = U} = B_U {Y in M_U O^{1 x k} : rk Y = j}; Theorem B with L_0 = O^{1 x k}, Gamma = GL_k(O), j <= r < k; covariance covol(N)^{-1} int phi o T = covol(TN)^{-1} int phi gives (lambda^{fn} c_U^k)^{-1} vol_{E_U} = mu_U, identical to R14. Same probability space for all U. Lower-rank measures are unchanged. OK.
- Theorem E' hypothesis list compared with R18 consolidated Theorem E: (T) incl. Delta_R = K^s and the invariant Euclidean structure (via "as in R18"), O maximal, (BS) in R18's form, (L1)-(L3), k >= k_0 non-effective. Complete.
- Arithmetic: p(p^2-1)/(4e (p-1)^2/2) = p(p+1)/(2e(p-1)); p = 13 gives 91/(12e). OK.

## Minor issues (none affects the mathematics)
1. Section 0 says "O maximal" while Theorem B is claimed for hereditary orders: state Theorem A for arbitrary orders and Theorem B for hereditary orders explicitly in the text (the proofs support this).
2. B.1: say explicitly that covol((AM)L_0) = covol(M L_0) (A unipotent) when passing to the WLOG.
3. B.2/B.6: note that eta_i (hence nu_i) exists because Gamma_i is a lattice in H_i, which uses rank k-1 >= 2, i.e. j >= 2; the j = 1 comparison is only made in that situation.
4. Theorem E': the bracket "(L1)_m = Gargava for (Delta^op, O^op)" is literally true only without real places; for K = R one needs R15 Lemma 4.2a (index <= 2), as the author notes in D.4. Only ranks k-r+1 <= m <= k are used.
5. C.2: G lies in GL_r(Delta) stabilising Lambda_W (not in GL_r(O)); c_U must be read as covol of M(Lambda_W cap U) in M U_R (as in R14).
6. D.9 remark on non-free L (the identity (Lambda_W cap U) L_0 = Lambda_W L_0 cap U^k) is not proved and not used; keep it tagged as a remark.
7. Corollary D' relies on R14 section 6 deriving (BS) without freeness; I only checked this against R18's statement of (R5-b) (duality argument), not R14 section 6 itself.

## Not verified by me
- (L1)_m itself, (Her) (Reiner 21.4, 2.44; I agree from memory), (U), R15 Lemmas 0.2, 0.5, 2.2, 2.6 as cited, R14 sections 2-7 beyond the grep, the W_p data, (BS). No numerics were run.

## What is established
Relative to (L1)_m for (O^{1 x m}, GL_m(O), G_m): Theorem A (any order, m >= 2), Lemma A.6, Theorem B (hereditary O, (L1)_m for k-j+1 <= m <= k), hence Theorem R of R15 without (h1^l), and Lemma 1.2' (R14's mean value formula with the same lower-rank measures, O maximal, no class number hypothesis). Theorem E' and Corollary D' follow provided R14 sections 2-5 (and 6 for D') are correct, which they inherit from R14/R18. The Weil-lattice corollary is conditional on maximality and (BS) and gives no new unconditional instance.
