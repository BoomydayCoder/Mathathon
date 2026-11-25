# R17 -- Referee report on `notes.md` (Gargava, arXiv:2107.04844v3, Theorem 16, case $D_\mathbb R=\mathbb H\times\mathbb H$)

**Verdict: CORRECT_MODULO_MINOR_FIXES.** I tried to refute the notes and found no fatal error. The description of the source is accurate, the false step is correctly located, and the repair (R1)-(R4) is complete. (R5) is in fact a complete proof for $s=2$, not only a sketch. The minor fixes are in section 3. The result stays conditional on one reduction-theory input (section 4).

Referee work was done under a 10-minute limit: the source was read directly in the TeX (l. 296-321, 494-1010, 1476, 1948-1955), Appendix C and [W58] were not read.

## 1. Spot-check of the source (TeX of v3, read directly)

| Item | TeX lines | Notes' quotation |
|---|---|---|
| Thm 16, Rem. 17 | 718-728 | exact |
| Prop. 14, Rem. 15 | 679-681, 712-714 | exact |
| $G,\Gamma$ (norm one), Rem. 7 | 550-559 | exact; the condition is total norm one |
| Norm = determinant of left multiplication | 1476 | exact |
| Thm 4 with proof "See [W58, Theorem 2]", Rem. 5 | 307-321, 498-500 | exact |
| Lemma 6 with proof "See [W58, Lemma 10.2]" | 516-521 | exact |
| Lemma 12 | 635-663 | accurate, including the slips reported |
| Lemmas 18, 19 with proofs | 757-1010 | accurate |
| False sentence "G acts transitively on $D_\mathbb R^k\setminus\{0\}$" | 993, reused 1002 | exact |

PDF numbering (Prop. 14, Thm 16, Rem. 17, Lemmas 18, 19) confirmed with `pdftotext`.

The sentence at l. 993 is false for $s\ge2$: $G\subseteq\prod_vGL_k(\mathbb K_v)$ cannot move $(x,0)$ to a vector with two nonzero components. I read l. 613-964 looking for a second use of transitivity and found none: Lemma 12, Prop. 14, Lemma 18 and the proof of Thm 16 do not use it.

## 2. Check of the repair, item by item

1. **(R1) no mass on lower strata.** Correct. A nonzero $w\in D^k$ has a coordinate $w_i\in D^\times$; then $\sigma_v(\mathrm{Nrd}\,w_i)\neq0$ at every real place, so every component of $w$ is nonzero. $G$ preserves $\Omega$. Hence $\nu(V\setminus\Omega)=0$ exactly, not only up to a null set.
2. **(R2) transitivity for the group actually used.** Correct for $G$ (total norm one) and for $\prod_vSL_k(\mathbb K_v)$, and it needs $k\ge2$. I re-checked the determinant correction: $\mathrm{diag}(1,c,1,\dots)$ fixes $e_1$ and has real determinant $N_{\mathbb K_v/\mathbb R}(c)$, which takes every positive value. The module convention is consistent (columns, right $\mathbb K_v$-vector space, a nonzero vector extends to a basis).
3. **G-invariance and local finiteness of $\nu$.** Invariance: immediate from invariance of $dg$ and $\Gamma\mathcal O^k=\mathcal O^k$. $\Gamma=\mathrm{Stab}_G(\mathcal O^k)$ is right ($\det\pi(\gamma)=\pm1$ gives $\gamma^{-1}\in M_k(\mathcal O)$). Local finiteness is exactly source Lemma 18; see item 6.
4. **Uniqueness.** Correct as used. Lebesgue measure on $\Omega$ is $G$-invariant ($|\det_\mathbb Rg|=1$ under either normalisation of $\mathrm N$), so an invariant Radon measure on $G/G_w\cong\Omega$ exists and is unique up to a scalar (Folland 2.49); unimodularity of $G_w$ follows from that of $G$ and is only needed for the unfolding variant (R3). For $\mathbb H\times\mathbb H$, $G\cong SL_k(\mathbb H)^2\times\mathbb R_{>0}$ is visibly unimodular, so Prop. 51 is not needed.
5. **Constant, circularity.** No circularity. $C$ is fixed by $\varepsilon\to0$: $\int\varepsilon^n\Phi_{f_\varepsilon}=C\int f$ for all $\varepsilon$ (from $\nu=C\,\mathrm{Leb}$), the pointwise limit is $\int f$ (Riemann sums, covolume 1), and dominated convergence uses only the uniform bound of Lemma 18, whose proof does not use Lemma 19. Continuity of $\varepsilon\mapsto\varepsilon^n\Phi_{f_\varepsilon}(g\Gamma)$ on $(0,1]$ justifies the supremum over rational $\varepsilon$.
6. **Lemma 18 for $s\ge2$ (re-derived, given Lemma 12 and Cor. 10).**
   - $Y=\{a'an a'^{-1}\}$ is relatively compact in $G$: $a'$ is a real scalar matrix, $a_{ii}\in\omega_1$, and $\omega_1$ is relatively compact in $D^{(1)}_\mathbb R$ because $\mathrm N(\omega_0)$ is bounded away from $0$ and $\infty$ (this is where Rem. 5 is used).
   - Exponent of $y_j$ against $dy_j/y_j$: $d\big(\tfrac jk\#I-\#I_{\le j}+j(k-j)\big)$, $I\subseteq\{1,\dots,k\}$. It is $>0$: $j(k-j)\ge j\ge\#I_{\le j}$, with equality only if $j=k-1$ and $\#I_{\le j}=k-1$, and then $\#I\ge1$. Agrees with the notes.
   - Modular factor of Cor. 10 re-derived: for $g=\kappa a'an$, right Haar measure on $A_0N$ is $\prod_{i<j}(a'_{ii}/a'_{jj})^d\,da'\,da\,dn$. Agrees with the source.
   - The dominating function is well defined on $G/\Gamma$ as $H_f=\sup_\varepsilon\varepsilon^n\Phi_{f_\varepsilon}$, bounded on each piece of $\mathfrak S^*$.
7. **(R5).** Correct and complete for $s=2$: $G=Z_+\times G^{(1)}$ is a direct product, $p(\Gamma)$ is discrete (reduced norms of units are totally positive units of the centre), $\mathcal F_Z\times\mathcal F_1$ is a fundamental domain, so $\nu=\mathrm{avg}_z\,\nu_1(z^{-1}\cdot)$. If $\nu_1$ were not locally finite on $\Omega$ it would be infinite on every nonempty open subset (transitivity plus invariance), contradicting $\nu=\mathrm{Leb}$; hence $\nu_1=C_1\mathrm{Leb}$ and $C_1=1$. Dirichlet's theorem is not needed here.

## 3. Minor issues

1. **Dependency on [W58, Lemma 10.2] can be dropped for $\mathbb H\times\mathbb H$.** Lemma 6 is elementary: for $m\in M_k(\mathcal O)\cap GL_k(D)$ the row lattice $\mathcal O^{1\times k}m\subseteq\mathcal O^{1\times k}$ has index $|\det_\mathbb R m|$, which is bounded; there are finitely many sublattices of bounded index; and $\mathcal O^{1\times k}m=\mathcal O^{1\times k}m'$ gives $m m'^{-1}\in GL_k(\mathcal O)=\Gamma$ (all real determinants are positive here).
2. **State the remaining citation precisely.** It is [W58, Thm 2] *in the strengthened form of source Remark 5* ($\omega_0$ relatively compact in $D_\mathbb R^\times$). Remark 5 is asserted in the source without proof or reference. The notes say this in sections 2.2 and 5 but the summary line "rest on [W58, Thm 2]" hides it.
3. **Source slip not listed in the notes.** With $\mathrm N=\det(L_a\mid M_k(D_\mathbb R))$ one has $\mathrm N(\lambda)=\lambda^{dk^2}$, so the scaling in $\mathfrak S^*$ and in Lemma 18 should be $\mathrm N(b_i)^{1/dk^2}$, not $\mathrm N(b_i)^{1/dk}$. Only the constant $R_i$ changes.
4. **Side remark in notes section 2 is inexact for general $s$.** "$G$ is bigger than $\prod_vSL_k(\mathbb K_v)$ by the $(s-1)$-dimensional factor $Z_+$" ignores the phases of $\det_\mathbb C$ at complex places and signs at real places. Exact for $\mathbb H\times\mathbb H$.
5. **Label of (R5).** Tagged [P] in the heading and "sketch level" in the text. It can be tagged [P] for $s=2$ (item 7 above).
6. **Scope of the general statement.** $D_\mathbb R\cong\prod_v\mathbb K_v$ excludes factors $M_m(\mathbb K_v)$, $m>1$; the notes say so. The source's sentence is false there too.

## 4. Ingredients that remain citations

- **(C1) Reduction theory.** Source Lemma 12: $\mathfrak S^*\Gamma=G$ with $\mathfrak S^*$ a finite union of translates of $KA^{(1)}_{\omega_1}A^\mathbb R_{c_1}N_{c_2}$, $\omega_1$ relatively compact. It rests on [W58, Thm 2] plus Remark 5. Not read by the checker or by me. Two partial checks of mine: for $k=1$ the content of Remark 5 follows from Minkowski's convex body theorem (a nonzero $m\in\mathcal O$ with $|m_1|\le c\lambda$, $|m_2|\le c/\lambda$ has $\mathrm N(m)\ge1$, hence $|m_1|\ge\lambda/c$, $|m_2|\ge1/(c\lambda)$); and for general $k$ the statement has the shape of the Borel-Harish-Chandra reduction theorem $G_\mathbb R=\mathfrak S\,C\,\Gamma$ with compact $\omega\subseteq{}^0P_\mathbb R$, the anisotropic kernel being compact modulo its arithmetic subgroup because $D$ is a division algebra. I did not match constants or conventions.
- **(C2) Haar measure in $KA'AN$ coordinates** (source Cor. 10, Prop. 8, App. C). Standard; I re-derived the modular factor but did not read App. C.
- **(C3) Standard measure theory.** Uniqueness of invariant Radon measures on $G/H$; the orbit map of a transitive action of a second countable locally compact group is a homeomorphism; Riemann sums over lattices; quotient integral formula (only for (R3)).
- **Corroboration, not checked.** The $SL_k(\mathbb H)^2/SL_k(\mathcal O)$ form is what Weil's adelic Siegel formula gives ($\tau(SL_k(D))=1$, strong approximation for $k\ge2$). Stated from memory.

## 5. What a reader may rely on

Let $D$ be a totally definite quaternion algebra over a real quadratic field (so $D_\mathbb R=\mathbb H\times\mathbb H$), $\mathcal O$ any order, $k\ge2$. Take either

- $G=\{(g_1,g_2)\in GL_k(\mathbb H)^2:\mathrm{Nrd}(g_1)\mathrm{Nrd}(g_2)=1\}$, $\Gamma=GL_k(\mathcal O)$, or
- $G^{(1)}=SL_k(\mathbb H)^2$, $\Gamma^{(1)}=\{\gamma\in M_k(\mathcal O):\mathrm{Nrd}(\gamma)=1\}$,

with the invariant probability measure on the quotient and Lebesgue measure normalised by $\mathrm{covol}(\mathcal O^k)=1$. Then, **conditional on (C1)** and on the standard facts (C2), (C3),

$$\int\sum_{v\in g\mathcal O^k\setminus\{0\}}f(v)\,dg=\int_{\mathbb H^{2k}}f(x)\,dx$$

for every Borel $f\ge0$, in particular for indicators of balls, and for all $f\in C_c$. All lattices $g\mathcal O^k$ have covolume 1.

The printed proof of Theorem 16 must not be cited verbatim for $s\ge2$: l. 993 is false there. Cite it with the repair (replace $D_\mathbb R^k\setminus\{0\}$ by $\Omega=(\mathbb H^k\setminus0)^2$, using R1 and R2).
