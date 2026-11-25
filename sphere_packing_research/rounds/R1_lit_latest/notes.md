# R1_lit_latest — newest results around Klartag's bound V_n >= c n^2 (searched 2026-11-25)

Normalization: lattice L of covolume 1, V(L) = volume of largest open origin-centred lattice-free ball; density = V(L)/2^n; V_n = sup_L V(L).
Access labels: [FULL] = full text downloaded and read (pdftotext); [ABS] = arXiv abstract only; [INFER] = my own inference, NOT stated in any paper.

## 0. State of the art in one table (as of Nov 2026)

| Result | Bound on V_n | Dimensions | Source | Access |
|---|---|---|---|---|
| Klartag (Invent. Math. 244 (2026) 1251-1279) | >= c n^2 | all n | arXiv:2504.05042 | FULL (uploaded) |
| Rezvyakova | >= (1/(2e) - o(1)) n^2 ~ 0.1839 n^2 | n -> infinity | arXiv:2608.22590 (23 Aug 2026) | FULL |
| Abuya-Gargava-Zhao (AGZ) | >= c N^2 log log N | sparse sequence N = phi(m) * ceil((log phi(m))^4), m primorial | arXiv:2606.05105v2 (3-4 Jun 2026) | FULL |
| Gargava-Groebner-Simkin (GGS) | >= c d^2 log log d | EVERY d >= d_0 | arXiv:2610.12619 (8 Oct 2026) | FULL (intro, Sec. 2, parts of 3, Sec. 4, App. A) |
| GGS Thm A.1 (no stochastic process) | >= (1-o(1)) (d/e) (M/C) >> d log log d | every large d | same | FULL |

So the "log log n" paper the task refers to is AGZ (sparse sequence), and it has since been superseded by GGS (all large dimensions).
No result of order n^2 log n or better was found. No constant is given in AGZ / GGS (only "universal c > 0").

## 1. Abuya-Gargava-Zhao, "Stochastically evolving ellipsoids with symmetries", arXiv:2606.05105 [FULL, 21 pp.]

Authors: Elisha B. Abuya (pseudonym; Tel Aviv Univ.), Nihar Gargava (Paris-Saclay), Yufei Zhao (MIT). v1 3 Jun 2026, v2 4 Jun 2026. math.MG/NT/PR.
"Statement on AI use": the key new idea (let the module rank grow) and a proof were produced by GPT-5.5 Pro prompted by Zhao; the first two authors had independently tried rank 2 and got a bound strictly between Venkatesh and Klartag. Chat log is an ancillary file (not read).

**Theorem 1.1.** There is a universal c > 0 and an infinite sequence of dimensions N with Delta_N^lat >= c N^2 log log N 2^{-N}.

**Theorem 2.2 (arithmetic lattice ellipsoid).** There is a universal c > 0 such that for ALL sufficiently large integers m there are a mu_m-stable lattice Lambda in V_R (dim N = phi(m) r, r = ceil((log phi(m))^4)) of covolume vol(B_N) and an ellipsoid E with E cap Lambda = {0} and vol(E) >= c m N r vol(B_N).
(So V_N >= c m N r = c N^2 m/phi(m) for every N of this form; log log gain needs m/phi(m) ~ log log m, i.e. m = prod_{p<=x} p, by Mertens. log log N ~ log log m.)

Setup.
- K = Q(zeta_m), degree phi(m), s = phi(m)/2 conjugate pairs of embeddings; V = K^r, V_R = (C^r)^s, N = phi(m) r.
- Lattices: Lambda_g = g O_K^r, g in SL_r(K_infty) = prod_{j<=s} SL_r(C); Haar probability on SL_r(K_infty)/SL_r(O_K).
- Symmetry: mu_m acts by scalar multiplication, orthogonally and FREELY on V_R minus 0, so every nonzero orbit has exactly m points.
- Mean value (Thm 2.1): for r >= 2, E sum_{0 != v in Lambda_g} f(v) = covol^{-1} int f. (Venkatesh for rank 2; Gargava Math. Z. 2023 Thm 16; S. Kim, adelic Rogers formula, JLMS 2024.)
- Quadratic forms: Q_K = {A in Sym(V_R): A is K_infty-linear} = direct sum of s copies of Herm(C^r); dim_R Q_K = s r^2 = N r/2. Every A in Q_K is mu_m-invariant; Id in Q_K; closed under inversion.
- Process: exactly Klartag's contact-preserving process dA_t = Pi_{A_t}(dW_t) but with W a Brownian motion in Q_K and F_A = {B in Q_K: B[x]=0 for all contact points x}.

Key steps (numbers are the paper's).
- Lemma 2.3: ||Pi_K(x (x) x)||_HS^2 = (1/2) sum_j ||x_j||^4; beta(x)^2 := that / ||x||^4 = (1/2) sum_j p_j^2 in [1/(2s), 1/2] where p_j = ||x_j||^2/||x||^2 are the "block masses".
- Lemma 3.3 (the gain): dim F_{A_t} >= N r/2 - |contacts|/m, because a whole orbit imposes ONE linear condition on Q_K.
- Prop 3.5: for 0 < T <= 8 log(mNr)/(Nr): E log det A_T <= -N r T/4 + (1/2) int_0^T E(|contacts|/m) dt + O(1).
- Lemma 3.6: P(x is a contact at time t) <= 2 Phi((a_0 - |x|^{-2})/(beta(x) sqrt t)) — variance of the martingale A_t[x] is reduced by beta(x)^2 <= 1/2 since only the Q_K-projection of x (x) x is felt.
- Lemma 4.1 (radial): J_t(omega) <= exp((1 + C sqrt(N t)) N^2 beta(omega)^2 t / 8).
- Lemma 4.2 (new analytic step): for omega uniform on S^{N-1}, E exp(lambda s sum_j ||omega_j||^4) <= C e^lambda for 0 <= lambda <= c sqrt r. (Concentration of block masses; E[s sum ||omega_j||^4] = s(r+1)/(sr+1) = 1 + 1/r + o(1).)
- Prop 4.5: E_lattice K_t <= (C/m) exp(N r t/8). Stopping time T = (8/(Nr)) log(mNr/(800 C)). Then det A_T = O((mNr)^{-2}), vol(E) >= c m N r vol(B_N).

Why the rank must grow (Remark 4.3, verbatim content): for fixed r the exponential moment is at least exp((1+1/r+o(1)) lambda), so the contact bound is (1/m) exp((1+1/r) N r t/8), forcing an earlier stopping time and a volume gain of only O((mNr)^{r/(r+1)}) instead of mNr; e.g. r = 2 gives exponent 2/3. Growing r makes the coefficient 1+1/r -> 1. r = (log phi(m))^4 is chosen so that log(mNr) = o(sqrt r) (range of Lemma 4.2) and N^{3/2} r T^{3/2} = O(log^{3/2}/sqrt r) = o(1).

Statements about limits / open problems in AGZ: there is NO open-problem section. The only "barrier" statement is in the introduction: "On various occasions, Venkatesh has conjectured the upper bound Delta_N^lat <= N^2 (log N)^{O(1)} 2^{-N}. This is motivated by the fact that lattices in R^N have on the order of N^2 degrees of freedom." Remark 1.2 notes the lattices are the module lattices of ML-KEM/ML-DSA and that "short vectors seem to prefer certain directions more than others", more so for small rank.

## 2. Gargava-Groebner-Simkin, "A cyclotomic gluing construction for lattice sphere packings", arXiv:2610.12619 [FULL for the parts listed above, 39 pp.]

**Theorem 1.1.** There are an absolute c > 0 and d_0 such that for every integer d >= d_0: Delta_d^L >= c d^2 log log d 2^{-d}. (Abstract also says this improves the best general, not necessarily lattice, bound.)

**Proposition 4.6.** For every sufficiently large d for which admissible gluing data exist: Delta_d^L >> d^2 (M/C) 2^{-d}.
**Theorem 3.1.** Admissible gluing data with M/C >> log log d exist for every large d.
**Theorem A.1 (no stochastic process; Rogers/Minkowski successive-minima stretching, tau-equivariant).** Delta_d^L >= (1-o(1)) (d M/(e C)) 2^{-d}, hence >> d log log d 2^{-d} in EVERY large dimension (Venkatesh had this only on a sparse sequence).

Construction.
- Finite index set N of positive integers (1 in N, 2 not in N), F = prod_{n in N} Phi_n, order A = Z[X]/(F) inside O = sum Z[zeta_n] with index prod |Res(Phi_m, Phi_n)|. C = sum phi(n) = rank of A, M = lcm(N).
- Lattice module L = A^r + Z^a, d = C r + a, 0 <= a < C, r = floor(d/C). The a spare coordinates sit in the n=1 component (trivial action) — this is how every dimension is reached.
- Cyclic action tau = multiplication by X, order M. Orbit of v has size M_S = lcm(supp v); FULL-support vectors have orbit size M; proper-support vectors may have smaller orbits.
- Ensemble: G = prod_n SL_{r_n}(K_n (x) R), Haar probability on G/Gamma, plus a deterministic volume-preserving rescaling e^{y_n} of the components. Support-wise Siegel formula (Thm 2.3): E sum_{v in L_S^o} f(gv) = covol(L_S)^{-1} int_{V_S} f.
- "Gluing" (A is not the direct sum O) makes proper-support vectors arithmetically expensive: covol(L_S) = exp(r I_inc(S)) with I_inc from cyclotomic resultants (Apostol: |Res(Phi_m,Phi_n)| = l^{phi(min)} if max/min is a power of the prime l, else 1).
- Support exponent A(S) = (c_S/2) log(C/c_S) + (c_S/C) I - I_inc(S) - sum_{n in S} c_n y_n; admissibility requires A(S) <= -delta for all nonempty proper S, plus I/C << log C, |y_1| << log C, r >> C^3 log d, C >= log d, log M << log d. Then proper-support short vectors have expected number <= e^{-cr}.
- Parameters realised: M odd squarefree built from primes up to X = (log d)^{1/20} arranged in "divisor cubes"; C ~ phi(M), C = d^{o(1)}, r = d^{1-o(1)}, M/C ~ M/phi(M) ~ log X ~ log log d.
- Invariant forms: Q = R Id_{r+a} + sum over conjugate pairs of Herm(C^r); q = dim Q = 1 + (C-1) r^2/2 ~ d r.
- General angular lemma (Lemma 4.1), valid for ANY subspace Q of Sym(R^d) containing Id with dim q: E ||Pi_Q(omega (x) omega)||_HS^2 = (2q + d)/(d(d+2)), and E exp(lambda F^2) <= K e^lambda for lambda <= c sqrt((2q+d)/d), F^2 = d(d+2)/(2q+d) ||Pi(omega (x) omega)||^2.
- gamma = d(2q+d)/(8(d+2)) ~ q/4, T = gamma^{-1} log(M gamma/K_0); vol(E)/v_d >> (M gamma/K_0)^{q/(4 gamma)} with q/(4gamma) = 2q(d+2)/(d(2q+d)) = 1 - O(1/r); so vol(E) >> M q v_d and M q >> d^2 M/C.
- Contacts are counted with weight 1/M_supp(x) (Prop 4.3): dim F_A >= q - sum_{contacts} 1/M_supp(x).
- No open-problem section. AI statement: Simkin + ChatGPT produced a weaker form first.

## 3. Rezvyakova, "On the constant in Klartag's theorem...", arXiv:2608.22590 [FULL, 8 pp.]

**Theorem 1.** Delta_n >= (1/(2e) - o(1)) n^2 2^{-n} as n -> infinity.
Method: pure bookkeeping of Klartag's proof. Start ball radius 1 - rho/n with rho = log log n (so the initial ball is lattice-free with probability >= 1 - e^{-rho}); T = (16 log n + gamma)/n^2; sharp Laplace asymptotics int_0^T (n sqrt t/2) J(t) dt <= (1+o(1)) (8/n^2) exp(n^2 T/8); E log det A_T <= -4 log n - gamma/4 + 2 rho + 4(1+o(1)) e^{-rho} e^{gamma/8}; optimum gamma = 8 rho + 8 log(1/2) gives log det A_T <= -4 log n + 2 log(2e).
Not a numerical simulation; no claim of optimality of 1/(2e).

## 4. Aubrun, "Sphere Packings in Higher Dimension (after Boaz Klartag)", arXiv:2606.13313 [FULL for intro + Sec. 4.1-4.2; Bourbaki-seminar exposition, 17 pp.]

- Reformulation: C_Lambda = {A in Sym_n : <Ax,x> >= 1 for all nonzero x in Lambda} (this is the Ryshkov polyhedron); Klartag's process is a "Brownian exploration" of this closed convex set, a martingale that sticks to faces; "For t large enough, it is an extreme point of C_Lambda. The Brownian exploration process 'chooses' an extreme sphere packing." Minimisers of det on C_Lambda are extreme points (log det strictly concave).
- States "For large n, no explicit example is known of a lattice sphere packing in R^n with density greater than 2^{-n}."
- Barrier remark: "[Venkatesh, 2013] conjectures that, up to logarithmic corrections, the quantity 2^n delta_n^L is bounded above by a polynomial in n. We now know that such a polynomial must have degree >= 2."
- Calls AGZ "an AI-generated argument".

## 5. Klartag, "Poisson approximation of random lattices", arXiv:2605.10479 [FULL for intro, 17 pp.]

**Theorem 1.2.** If S in R^n is measurable, S cap (-S) empty, 0 < Vol(S) < c n, then d_TV(L cap S, Poisson cap S) <= C e^{-c' n} for a Haar-random unimodular lattice L.
Stated limitation: "The conclusion of Theorem 1.2 is false, in general, when Vol_n(S) >= n" (n lattice vectors in a half-ball have integer determinant; also E|L cap S|^n = infinity).
Relevance: the Poisson heuristic behind Klartag's n^2 is rigorous in total variation only up to volume ~ n; the n^2 theorem uses only the FIRST moment (Siegel) plus the martingale structure.

## 6. Other 2025-26 items found (complementary)

- Klartag's own paper: conjecture that Theorem 1.2 (lattice-free body of volume c n^2) holds for every origin-symmetric convex body K, known (Schmidt) only for volume <= c n. No paper proving it was found. [FULL, uploaded]
- Non-lattice packings: no improvement specific to non-lattice packings beyond CJMS (1/2 - o(1)) n log n was found; the best general bound is now the lattice bound c d^2 log log d (GGS). Related: Bukh-Gao-Liu-Pikhurko-Sun arXiv:2510.25685 (coverings; random-translate method fails below (1/2+o(1)) n ln n density) [ABS].
- Kissing numbers of lattices: Bennett-Hhan-Li-Stephens-Davidowitz arXiv:2610.14380: lattice kissing number >= 2^{Omega(sqrt n)}, improving the 60-year-old 2^{Theta(log^2 n)}; Laarhoven-Kominers arXiv:2609.13608: >= (1/(2 pi e^2)+o(1)) sqrt q e^{2 sqrt q} for prime powers q (Reed-Solomon / Bennett-Peikert type, efficient randomized or explicit constructions); Ducas-Nowakowski-Robertson arXiv:2611.22441 (Craig lattices). [ABS] Not related to the ellipsoid method.
- Lattice-specific upper bounds: Zubrilina arXiv:2610.35830 modifies the Cohn-Elkies LP using "a lattice cannot have too many short vectors"; improves lattice packing upper bounds for 11 <= n <= 16 and many n up to 91. [ABS] No asymptotic statement in the abstract.
- Exact small dimensions: Dutour Sikiric-van Woerden arXiv:2508.20719: all 2,237,251,040 perfect lattices in dimension 9 enumerated by Voronoi's algorithm; Lambda_9 is the densest lattice, gamma_9 = 2. [ABS]
- Vanoni-Guo-Torquato arXiv:2608.00920: conditional (on the Torquato-Stillinger realizability conjecture) polynomial improvements d^beta 2^{-d} for every fixed beta > 1 for (disordered, non-lattice) packings. [ABS]
- Gargava arXiv:2607.11695 (counterexamples in dims 8, 9 to Cartesian-product conjectures), Gargava-Serban-Viazovska-Viglino arXiv:2510.12893 (module lattices and shortest vectors; variance bounds). [ABS]
- NOT found: numerical simulation of Klartag's process; derandomized / explicit construction of a lattice with V >> 1; any proof that n^2 polylog is an upper barrier; any follow-up to GGS (searched arXiv API up to 2026-11-24).
- Could not access: Gil Kalai's blog post on Klartag's result (host not allowed in this session); the AGZ ancillary chat log.

## 7. My own inferences (NOT in the papers; unverified, to be checked by later rounds)

(I1) Meta-theorem shape. Klartag/AGZ/GGS all fit: take a Jordan subalgebra Q of Sym(R^d) containing Id (dimension q), a group of order M whose orbits impose one linear condition on Q each, and a lattice ensemble with a Siegel mean-value formula in which small-orbit vectors are negligible. Then vol(E)/v_d >~ (M gamma)^{q/(4 gamma)}, gamma = d(2q+d)/(8(d+2)); this is ~ M q when q >> d polylog. So the achievable volume is "(orbit size) x (number of free parameters of the invariant ellipsoid)".

(I2) Conjectural constant. Repeating Rezvyakova's optimisation with general (M, q): minimise -2u + (4/(M q)) e^{-rho} e^{u} + 2 rho over u = qT/4 gives log det A_T <= -2 log(M q/2) + 2, i.e. vol(E)/v_d >= (1-o(1)) M q/(2e). Check: M = 2, q = n(n+1)/2 gives n^2/(2e) (Rezvyakova). For AGZ (M = m, q = N r/2) this would give V_N >= (1-o(1)) (m/phi(m)) N^2/(4e), i.e. ~ (e^{gamma_Euler}/(4e)) N^2 log log N ~ 0.164 N^2 log log N. Needs the (1+o(1)) versions of AGZ Lemmas 4.1-4.2, which look available because r -> infinity.

(I3) Character formula for the gain of a general finite symmetry group G acting on V = R^d with generic orbit size |G|: |G| dim (Sym^2 V)^G = (1/2) sum_{g in G} (chi(g)^2 + chi(g^2)). The g = +-1 terms give Klartag's d(d+1); all gain comes from other elements with large |trace|. For mu_m on K^r: sum_g chi(g)^2 = r^2 sum_j c_m(j)^2 = r^2 m phi(m) = d^2 m/phi(m). Free modules over Z[G] (regular representation) give exactly d^2/2 x 1: no gain. For cyclic G and a union of faithful-ish cyclotomic components the gain is M/C <= M/phi(M) <= (e^{gamma_Euler}+o(1)) log log M, so cyclic symmetry cannot give more than log log. A quick calculation for G = (mu_p)^k acting componentwise on (K_p^r)^k with maximal mod-(1-zeta) gluing suggests the GGS admissibility condition fails (proper supports of one component have exponent ~ r[((p-1)/2) log k - (log p)/k] > 0), so the naive non-cyclic abelian generalisation does not work as is.

(I4) Deterministic counterpart: the restricted process ends at a "T-perfect" / G-perfect form in the sense of Berge-Martinet-Sigrist (Voronoi theory in a linear subspace T of Sym). For T = Hermitian forms over Z[zeta_m] this is the theory of perfect Hermitian forms over CM fields, where exact enumeration in small rank exists (not read; e.g. Parker-Yasaki arXiv:2609.31479, de la Torre Duran arXiv:2611.06300 seen only as titles). This is the natural place for the "small n, exact coordinates" experiments.

## 8. Files

- BASE/papers/abuya_gargava_zhao_2606.05105/2606.05105.{pdf,txt}
- BASE/papers/gargava_groebner_simkin_2610.12619/2610.12619.{pdf,txt}
- BASE/papers/rezvyakova_2608.22590/2608.22590.{pdf,txt}
- BASE/papers/aubrun_survey_2606.13313/2606.13313.{pdf,txt}
- BASE/papers/klartag_poisson_2605.10479/2605.10479.{pdf,txt}
- BASE/papers/bennett_kissing_2610.14380/2610.14380.{pdf,txt} (downloaded, only abstract read)
- BASE/papers/zubrilina_small_dims_2610.35830/2610.35830.{pdf,txt} (downloaded, only abstract read)
- BASE/papers/gargava_counterexamples_2607.11695/2607.11695.{pdf,txt} (downloaded, only abstract read)

(I5) Why "free action / V = D^r" seems forced (heuristic, unverified). If G acts on V = W^r with W absolutely irreducible of dimension d_W > 1 and End_G(W) = Q (e.g. S_k on its standard representation, where the naive gain |G|/d_W^2 = k!/(k-1)^2 looks huge), then away from primes dividing |G| every G-invariant lattice is a tensor product A (x) L' (A a G-lattice in W, L' a lattice in R^r). Tensor products contain pure tensors a (x) x of norm min(A) min(L'), and since min(A) < sqrt(d_W) x (covol A)^{1/d_W} by Minkowski, there are exponentially many (in r) pure tensors inside the ball of volume poly(d). So the first-moment method fails by an exponential factor. This is the same phenomenon as the "direct product has too many proper-support vectors" toy example in GGS Sec. 1.4. The cases without this obstruction are V = D^r with G a finite subgroup of D^x for a division algebra D (Venkatesh, Gargava, AGZ), where the gain is ~ |G|/dim_Q D, believed to be O(log log) (from memory of Venkatesh 2013 / Gargava 2023; not re-checked here). Any super-log-log gain from symmetry would need (a) a non-division-algebra symmetry together with a GGS-style "gluing at bad primes" argument that suppresses the pure-tensor / proper-support vectors, or (b) a genuinely different source of degrees of freedom.

Numerical sanity checks done for (I2), (I3): sum_j c_m(j)^2 = m phi(m) verified for m = 3, 5, 12, 30, 210; the optimisation in (I2) reproduces n^2/(2e) for (M, q) = (2, n^2/2); e^{gamma_Euler}/(4e) = 0.1638.
