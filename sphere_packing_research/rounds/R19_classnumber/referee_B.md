# R19. Independent referee report (referee B) on `proof.md`

(Saved by the orchestrator from the referee agent's final report, verbatim; the agent was instructed not to write files.
Hard limit of the agent: 10 minutes. It read the version of `proof.md` available at about 10:21-10:25 UTC.)

VERDICT: CORRECT_MODULO_MINOR_FIXES (wording-level only). I found no mathematical error in Part B, in Lemma 1.2' / Part C, or in the use of (Her). I also re-derived Part A (A.1-A.6) and found it correct as written.

File refereed: rounds/R19_classnumber/proof.md. Not checked by me: the analytic correctness of R14 sections 2-5 (I read them only for uses of freeness), anything in R6/R17, and the Reiner theorem numbers.

## (1) Part B, step by step - every item holds

- **B.1 (line 59-61).** The section of pi_1 extends to x -> sx with pi_1(s) = 1, so M = sI_1 (+) (M cap ker pi_1). A fixes ker pi_1, so AM = e_1 I_1 (+) (0 x M_1). Y -> AY commutes with the right action, preserves rank, and has real determinant 1 on Delta_R^{j x k} (k copies of a unipotent column map), so F becomes F o A^{-1} and covol(ML_0) is unchanged.
- **Product structure.** ML_0 = e_1(I_1 L_0) (+) (0; M_1 L_0), because generators with x = 0 or m_1 = 0 are allowed. The blocks sit in complementary coordinates, so covol(ML_0) = covol(N_1) covol(M_1 L_0).
- **N_1 = I_1 L_0.** It is a full Z-lattice between c^2 O^{1xk} and c^{-2} O^{1xk}, and Gamma-stable since (xl)gamma = x(l gamma). It is only a module over the left order of I_1, which is fine because Theorem A needs no module structure. Gamma = Stab(L_0) is commensurable with GL_k(O) by A.6, so Theorem A applies to (N_1, Gamma).
- **B.2.** The orbit decomposition, the Gamma_i-invariance of Phi_i and the Weil unfolding are correct (G_k is unimodular as the kernel of a character of GL_k(Delta_R)). The j = 1 unfolding uses the same Gamma, orbits, eta_i and hence the same nu_i, so the comparison in B.6 is legitimate.
- **B.3 conjugation.** g_i H_i g_i^{-1} = Stab_{G_k}(e_1) by normality, and g_i Gamma_i g_i^{-1} = Stab_H(L_0 g_i^{-1}) exactly.
- **B.3 splitting.** The section is x'' -> (x''t | x'') with t in Delta^{(k-1)x1}, and (0|x'')u(t,1) = (x''t | x''), (a|0)u(t,1) = (a|0). So L_i u(t,1)^{-1} = Je_1 (+) (0|L_i''), and e_1 u(t,1) = e_1 keeps e_1 g_i = y_i. p(L_i) is unchanged by the replacement.
- **Stab_H(L_i).** Both inclusions hold: (a|x'')u(n,h) = (a + x''n | x''h) forces L_i''h in L_i'' and L_i''n in J; the inverse u(-h^{-1}n, h^{-1}) gives L_i''h^{-1}n = L_i''n in J. So it is exactly T_i semidirect Gamma_i''.
- **T_i.** c^2 O^{(k-1)x1} in T_i in c^{-2} O^{(k-1)x1}, so it is a full lattice, and hT_i = T_i. The fundamental domain F_i'' x P_i and the normalisation covol(T_i)^{-1} dn dmu_i'' are correct.
- **M_1 L_i.** M_1 J is a full lattice in Delta^{(j-1)x1}, and M_1 L_i = (M_1 J | 0) (+) (0 | M_1 L_i'') is direct (different columns), with product covolume.
- **B.4.** Y''T_i in M_1 J; Y'' has a right inverse over Delta, so n -> Y''n is R-linear onto and the torus map is a well-defined surjection carrying Haar probability to Haar probability. The displayed identity follows.
- **B.5.** The induction hypothesis is applied to (M_1, L_i''), with Gamma_i'' = Stab_{G_{k-1}}(L_i''), the probability measure, and 1 <= j-1 < k-1 (k-1 >= 2). The Jacobians ||g_i||^{-(j-1)} on both sides cancel, so the fibre value covol(M_1 L_0)^{-1} int Psi is independent of the orbit.
- **B.6.** Correct.

**Numerical tests of the bookkeeping.**
- O = Z, k = 3, j = 2, M = Ze_1 + 2Ze_2, L_0 = Z^3, orbit y_i = d e_1, g_i = diag(d,1,1): J = (1/d)Z, T_i = (1/d)Z^2, covol(M_1 J) = 2/d, covol(M_1 L_i'') = 4, ||g_i|| = d. Fibre constant = (d/2)(1/4)(1/d) = 1/8 for every d; total is 2^{-3} int F, matching Rogers.
- L_0 = Ze_1 + Ze_2 + (1/2)Ze_3, M = Z^2, orbit y_i = (1/2)e_3 with g_i the permutation-type matrix of ||g_i|| = 1/2: L_i = Z^3, fibre constant = 2 = covol(L_0)^{-1}. Globally, conjugating by g = diag(1,1,1/2) reduces to Rogers for Z^3 and gives 2^j int F = covol(L_0^j)^{-1} int F.
- Same L_0 with M = Ze_1 + 2Ze_2: both sides give 1/2 (covol(ML_0) = 2).

## (2) Part C - R14 sections 0-5 read line by line for freeness

Freeness, the form M O^{r x k} gamma, or class number one enter only in:
- 0.1 (coordinates);
- Lemma 1.1;
- the proof of Lemma 1.2 (O^r cap U = B_U O^j);
- pure notation "X in O^{r x k}" / "M O^{r x k} gamma" in Lemma 0.5, Prop. 4.1 and section 5.

Everything else is free of it:
- Lemma 0.1 uses only W_R and the real span of G.
- Lemmas 0.2-0.4 concern E, Q and E_U only.
- Lemma 0.5 needs only gX = X with g Delta-linear; it holds for X in Lambda_W^k.
- Lemma 1.4 is numeric in c_U, c_W.
- Lemma 1.5 uses Lambda_W as a Z-lattice plus Minkowski.
- Lemma 2.2 uses u_g Lambda = Lambda, A_t in Q, and trivial stabilisers outside Lambda^{nf}.
- Props 2.1, 2.3, 3.1-3.5 and Cor. 3.6 are for arbitrary lattices or for E_U.
- Prop. 4.1 uses (MVF), its fixed-U form, and Lemmas 0.5, 1.4, 1.5.
- Section 5, last step: sqrt(A_T) = R_N gives lambda M Lambda_W^k gamma N, still in the family and G-invariant.

**Lemma 1.2' (C.3) is correct.**
- {X in Lambda_W^k : col(X) = U} = {X in (Lambda_W cap U)^k : rk X = j}.
- M_U O^{1xk} is the set of matrices with all columns in M_U, and B_U M_U = Lambda_W cap U, so (Lambda_W cap U)^k = B_U(M_U O^{1xk}); rk(B_U Y) = rk Y.
- Gamma = GL_k(O) is the full stabiliser of O^{1xk}, as Theorem B requires; j <= r < k.
- The mean is (lambda^{fn} c_U^k)^{-1} int_{E_U} phi, which equals kappa_n^{-f} (c_W^f/c_U)^k vol_{E_U}: mu_U is unchanged, including U = W.

## (3) Hypothesis "O maximal"

- It is the right hypothesis for this proof. Only hereditarity is used: finitely generated torsion-free modules embed in free ones and are therefore projective.
- It is used on the correct sides: right modules in B.1 (I_1 = pi_1(M)), left modules in B.3 (L_i'' = p(L_i)).
- Maximality is left-right symmetric (O^op is maximal in Delta^op), so both are available. The pair (M_1, L_i'') in the induction is again (right O-lattice, left O-lattice) over the same O.
- The remark at line 102 is correct: the right order of Lambda_W O' is O', but c_W, c_U and hence (BS) change.

## Minor fixes (none affects validity)

1. **B.5, line 78-79:** say explicitly that the integrand is left Stab_H(L_i)-invariant before integrating over F_i'' x P_i (it is).
2. **B.3, line 68:** state that k-1 >= 2 because j >= 2; the forward reference from B.2 to B.3 is harmless.
3. **C.2 / C.4:** state once that every "O^{r x k}" and "G in GL_r(O)" in R14 (0.1, Lemma 0.5, Prop. 4.1, section 5) is to be read as Lambda_W^k and G in GL_r(Delta) with g Lambda_W = Lambda_W. Also that R14 line 23 (G_k/Gamma_k, N(gamma) = 1) is replaced by the R18 convention Gamma_k\G_k, ||gamma|| = 1.
4. **Theorem E' (line 100):** k_0 must include k >= max(r+1, N_W^2+2) (Lemma 1.5). If some K_v = R, the input (L1)_m for G_m = {||.|| = 1} comes from Gargava only via R15 Lemma 4.2a (R18 problem 6). It is listed as an input, so this is not a gap of R19.
5. **(Her), line 16:** Reiner numbering is from memory in the proof and unverified by me; the mathematical statement is standard.

## Strongest statement actually proved

- **Theorem A** as stated, from (L1)_m (Borel form), (U) and (Tr).
- **Theorem B** as stated, for every hereditary order O, every full right O-lattice M in Delta^{j x 1}, every full left O-lattice L_0 in Delta^{1 x k}, 1 <= j < k, Gamma = Stab_{G_k}(L_0), constant exactly covol(ML_0)^{-1}; from (L1)_m for m = k-j+1..k, (U), Weil's quotient formula and (Her).
- **Lemma 1.2'** (MVF with the same mu_U, including the fixed-U form) for O maximal, with no class number hypothesis.
- **Theorem E'** then holds exactly to the extent that R14 sections 2-5 are analytically correct and (L1)-(L3) hold; no step of R14 sections 0-5 other than those listed in (2) uses freeness.
