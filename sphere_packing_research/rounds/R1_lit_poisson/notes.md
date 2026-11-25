# R1_lit_poisson — Poisson statistics and hole probabilities for Haar-random lattices

Written 2026-11-25 07:01–07:14 UTC (15-minute budget). Notation: L Haar-random in X_n (covolume 1), B_V = origin-centred ball of volume V, N_V = #(L ∩ B_V \ {0}), P_V = Prob(N_V = 0), V_n = sup{V : P_V > 0}.

Access labels: [FULL] = full text downloaded to BASE/papers/<dir>/ and the quoted statements were read there; [ABS] = abstract/intro only; [2ND] = statement taken from a secondary source that I read; [MEM] = from memory, not verified; [INFER] = my own inference, NOT in any paper.

Hosts that could NOT be reached from this session (so the originals were not read): math.stanford.edu (Venkatesh's note), ams.org (Schmidt 1958), projecteuclid.org (Schmidt 1963 + Correction). Rogers 1955/56 were not fetched either. Everything about those papers below is [2ND] or [MEM].

New folders created in BASE/papers/: kim_On_1611.02880, strombergsson_sodergren_1611.06332, sodergren_poisson_1001.3623, han_primitive_2407.00986, kim_adelic_rogers_2205.03138, vanoni_guo_torquato_2608.00920, torquato_stillinger_0508381 (schmidt_1958_admissible, schmidt_1963_MH, venkatesh_sp_note are empty: download failed). Already present and used: klartag_poisson_2605.10479, gsv_moments_2308.15275, aubrun_survey_2606.13313, uploaded_txt/klartag.txt.

---------------------------------------------------------------------------

## 1. Summary table: what is rigorously known about N_V and P_V

| Range of V | Statement | Source | Access |
|---|---|---|---|
| V fixed, n -> oo | N_V/2 -> Poisson(V/2) in distribution (any symmetric Borel set) | Rogers 1956 (PLMS 6) | 2ND (Kim Thm 1.1, GSV intro) |
| V fixed | lengths (volumes) of the shortest vectors -> Poisson process of intensity 1/2 on R_+, directions uniform | Sodergren 2011 (Math. Z. 269) | 2ND (Kim Thm 1.2, S-S Thm 1.1) |
| k-th moment, k <= 2 sqrt(n-3) | 2^k m_k(V/2) <= E N_V^k <= 2^k m_k(V/2) + E_{k,n} (V+1)^{k-1}, E_{k,n} <= 2*3^{ceil(k^2/4)} (sqrt3/2)^n + 21*5^{ceil(k^2/4)} 2^{-n} (m_k = k-th Poisson moment) | Rogers 1956, Thm 3 | 2ND (quoted verbatim in GSV p.2) |
| k <= n-1 | Rogers integral formula valid; sum over linearly independent k-tuples has mean exactly the Lebesgue integral (no error term). k = n: E N^n = +oo | Rogers 1955 / Schmidt 1957; divergence: Klartag eq. (6) | FULL (klartag_poisson, Sec. 1-2) |
| V = o(sqrt n) | Poisson for P(N_V/2 >= k) | Kim 2016 (Math. Z. 282) | 2ND (Kim 2020 intro) |
| V + k <= (n/4) log(4/3) (e.g. <= 0.05 n) | (1 - e^{-n}) Q(V/2,k) - e^{-cn} Q(V/2,k-1) < P(N_V/2 >= k) < (1 + e^{-n}) Q(V/2,k) + e^{-cn} Q(V/2,k-1) (Q = Poisson upper tail) — error is RELATIVE to neighbouring Poisson tails | Kim 2020 (IMRN), Thm 3.1 | FULL |
| V, k <= C n | P(N_V/2 >= k) = Q(V/2,k) + O(e^{-cn}) (additive) | Kim 2020, Prop 3.3 | FULL |
| V <= C n (C = 0.025 works) | P(L ∩ S' linearly dependent) < e^{-cn}; for S a ball this probability is >= V e^{-0.7 n} | Kim 2020, Prop 3.1 + Remark | FULL |
| V < C n, n >= 13 | P_V = e^{-V/2} + O(e^{-cn}) (additive error) | Schmidt 1958 (PAMS 9) | 2ND (Kim Thm 1.5; Klartag-Poisson p.3) |
| Vol(S) < c n, S ∩ (-S) = empty | d_TV(L ∩ S, Poisson process ∩ S) <= C e^{-c' n} (whole point process, not just counts) | Klartag, arXiv:2605.10479 (11 May 2026), Thm 1.2 | FULL |
| Vol(S) >= n | the TV statement is FALSE in general (n lattice vectors with integer determinant) | Klartag 2605.10479, p.3 | FULL |
| V = f(n) -> oo, f = O_eps(e^{eps n}) | (N_V - V)/sqrt(2V) -> N(0,1); functional CLT (Brownian motion); P(N_x <= 2N) - P(Poisson(x/2) <= N) -> 0 uniformly for min(x,N) <= f(n) | Strombergsson-Sodergren 2019 (arXiv:1611.06332), Thms 1.2, 1.3, 1.3', 1.5 | FULL |
| V = O(e^{cn}), c < c_k | k-th moment of (N_V - V)/sqrt(2V) -> Gaussian moment; diverges for c > c_k (k >= 3). c_k = ((k-1)/(k-2)) log(k-1) - log k; c_3 = 0.28768, c_4 = 0.26162, c_k ~ log k / k | S-S 2019, Thm 1.7 | FULL |
| primitive points, V fixed / subexponential | Poisson(…)/normal limits for primitive lattice point counts | J. Han, arXiv:2407.00986 | ABS |
| O_K-module lattices of rank t | k-th moments of (number of points)/omega_K -> Poisson(V/omega_K) for t > t_0(K,k) = O_K(k^3 log log k); uniform in K under a Weil-height (Bogomolov) lower bound; cyclotomic: V^2 + V omega_K <= E rho^2 <= V^2 + V omega_K (1 + C e^{-eps d (t - t_0)}), t >= 27 | Gargava-Serban-Viazovska arXiv:2308.15275, Thms 1, 2, Cor 3 | FULL (intro) |
| number fields | adelic Rogers formula, second moment formulas | S. Kim arXiv:2205.03138 (JLMS) | ABS |

IMPORTANT distinction. The CLT / "Poisson" results for V -> oo (S-S) are statements about the BULK of the distribution (fluctuations of size sqrt V around the mean V). They say nothing useful about the hole event {N_V = 0}, which sits sqrt(V/2) standard deviations into the tail. The only results that control the hole event to the Poisson scale e^{-V/2} are the Schmidt-sieve results, and they stop at V <= c n.

## 2. What exactly does Schmidt prove? (originals NOT accessed)

- Schmidt 1958, "The measure of the set of admissible lattices" [2ND via Kim Thm 1.5 and Klartag]: there are c, C > 0 with P(L ∩ S = empty) = e^{-V/2} + O(e^{-cn}) for symmetric Borel S of volume V < C n, n >= 13. Klartag (2504.05042, p.2): "Schmidt proved (2) [i.e. V_n >= cn] by considering random lattices and by analyzing large hole events; these are rare events that occur with a probability of only exp(-c~ n)."
  [MEM, unverified] The 1958 paper states the result in multiplicative form, measure = e^{-V/2}(1 - R) with |R| < 6 (3/4)^{n/2} e^{4V} + V^{n-1} n^{-n+1} e^{V+n} for V <= n - 1. The factor (3/4)^{n/2} is consistent with Klartag-Poisson eq. (11): E rho^{k-1}_k(L ∩ S) <= lambda^{k-1}/(k-1)! * [3^k (3/4)^{n/2} + 5^k 2^{-n}] ("Schmidt [12, Theorem 3]") and with Kim's hypothesis V + k <= (n/4) log(4/3).
- Schmidt 1963, "On the Minkowski-Hlawka theorem" (+ Correction) [2ND/MEM]: Klartag: "See Schmidt [28, 29] for a proof [of an admissible lattice for any convex body K] under the weaker assumption Vol_n(K) <= cn." [MEM] the theorem is for arbitrary Borel sets S, of the shape V(S) <= n*const - c (const = log sqrt2 for general S, twice that for symmetric S; I recall that the Correction changes the constant — NOT verified, do not quote the constant).
- So the answer to "what does Schmidt prove about P_V for V up to cn": P_V = e^{-V/2} up to an ADDITIVE error e^{-cn}; hence P_V > 0 as long as e^{-V/2} dominates, i.e. V < 2c n. It is a second-and-higher-moment (Bonferroni) argument, valid for every Borel set, and it is the only proof of V_n >> n that goes through an actual estimate of P_V.
- Barrier of the sieve (Kim 2020, Remark after Prop 3.1, FULL): P(L ∩ B_V' has a rank-1 dependency) > V/2^{n+1}, "we cannot explore with it a ball of volume greater than 1.4 n, as Q(V/2,k) will typically have size about e^{-0.7n} and thus will be clouded by the error term caused by linear dependence". (e^{-V/2} = 2^{-n} iff V = 2 n log 2 = 1.386 n.)
- Kim 2020 intro (FULL): "We suspect that the same would hold for V < 2n, cf. Section 3.1 of Venkatesh. ... the error in (1.1) must be strengthened to at least e^{-O(n log n)}, which seems impossible to finesse from the technique that we employ here"; "V ~ 2n is a transition point between some sort of randomness and structure. If V = 2.01 n, the vectors in L ∩ B(V) must satisfy at least e^{O(n)} linear relations on average"; "further shrinking the error term in (1.2) immediately leads to a better lower bound [on packing density]".

## 3. Bounds on P_V for n <= V <= n^2

Published, rigorous:
- Upper bounds. (i) Monotonicity + Schmidt: P_V <= P_{Cn} <= e^{-c n} for all V >= C n. (ii) Chebyshev with Rogers' second moment (Var N_V ~ 2V, n >= 3): P_V <= (2 + o(1))/V. (iii) [INFER from S-S Thm 1.7 + Markov on the k-th central moment] P_V <= (1 + o(1)) (k-1)!! (2/V)^{k/2} for fixed even k, V -> oo, V = O(e^{cn}), c < c_k. None of these is of Poisson strength e^{-V/2} once V >> n. I found NO published bound of the form P_V <= e^{-c V} for any V >> n.
- Lower bounds. P_V > 0 iff V < V_n (V <= V_n), and V_n >= c n^2 (Klartag), >= (1/(2e) - o(1)) n^2 (Rezvyakova 2608.22590, per R1_lit_latest). I found NO published quantitative lower bound on P_V for V >= C n.

Does Klartag's theorem/method give a quantitative lower bound on P_V at V = c n^2? Not directly [FULL reading of 2504.05042]: the proof (Prop 5.1, Lemma 5.2) shows E log det A_T <= C' - 4 log n at T = 16 n^{-2} log n, hence "with positive probability" over the PAIR (Haar lattice L, Brownian path) the final L-free ellipsoid is large. The lattice that has the round hole is A_T^{1/2} L, whose law is NOT Haar and whose density w.r.t. Haar is not controlled. So no Haar-measure estimate comes out of the paper.

[INFER — my own computation, needs checking] Neighbourhood bound. Let L_0 attain V_n (exists by Mahler compactness; it is extreme hence perfect, so its minimal vectors span R^n and near-identity automorphisms are trivial). If g in SL_n(R) has smallest singular value >= (V/V_n)^{1/n} then g L_0 has a hole of volume V. This set is bi-K-invariant; in the cone of positive quadratic forms it is {Q : det Q = 1, Q >= (1 - 2 eps) Q_0}, eps ~ delta/n, V = e^{-delta} V_n, a convex body of dimension N - 1 = n(n+1)/2 - 1 whose (cone-)Lebesgue measure is about (c delta / n^{3/2})^{N}. Dividing by Minkowski's volume of the reduced domain (2/(n+1)) prod_{k=2}^n pi^{-k/2} Gamma(k/2) zeta(k) = e^{(1/4 + o(1)) n^2 log n} gives
   P_V >= exp( -(1 + o(1)) n^2 log n )   for all V <= (1 - 1/log n) V_n.
Details of the estimate (so that it can be re-checked):
 - Body: after normalising Q_0 = I it is {I + S : tr S ~ 0, S >= -2 eps I} = (2 eps)^{N-1} * {traceless S >= -I}, and {S >= -I, tr S = 0} = n * (real density matrices) - I. Hilbert-Schmidt volume of the set of density matrices has volume radius ~ n^{-3/2} (Zyczkowski-Sommers formula in the complex case: log Vol = -(3/2) n^2 log n + O(n^2) in dimension n^2 - 1; real case analogous [MEM]). Hence volume radius n^{-1/2} for {traceless S >= -I} and (c delta n^{-3/2})^N for the body: exponent -(3/4 + o(1)) n^2 log n.
 - Normalisation: checked numerically (python) that log of Minkowski's volume = (n^2/4)(log n - 1/2 - log(2 pi e)) + ..., i.e. ratio to (n^2/4) log n is 0.15, 0.37, 0.52 at n = 50, 200, 1000 — the o(1) is of size 3.3/log n, so the statement is truly asymptotic.
 - Multiplicity of the map to X_n: not injective in general, but if Q and U^T Q U (U in GL_n(Z)) both lie in the body then U^T Q_0 U <= e^{2 delta} Q_0, so U sends n independent minimal vectors of Q_0 to vectors of norm <= e^{delta} r_0; there are at most (2 e^{delta} + 1)^n such vectors, so the multiplicity is <= e^{O(n^2)}, which is lower order.
Consequences (all INFER):
 (a) e^{-(1+o(1)) n^2 log n} <= P_V <= e^{-cn} is the whole rigorous knowledge for C n <= V <= c n^2.
 (b) Conditional barrier: IF the Poisson upper bound P_V <= e^{-(1-o(1)) V/2} held for all V up to 3 n^2 log n, THEN V_n <= (2 + o(1)) n^2 log n. Equivalently, in units of N = n(n+1)/2 degrees of freedom: Poisson mean V/2 <= (2 + o(1)) N log n. Klartag's V/2 ~ c N is one log n below this heuristic ceiling.
 (c) Contrapositive: if V_n is super-polynomial (e.g. exponential), then P_V must exceed the Poisson value e^{-V/2} by a factor e^{+V/2 - n^2 log n} for V near V_n, i.e. the hole event is hugely positively correlated ("structure").
 (d) Caveat on (b): the natural correlation for hole events of a GROUP is positive (absence of x and y makes absence of x + y more likely), so P_V >= e^{-V/2} is the plausible direction and P_V <= e^{-V/2} is exactly the direction with no supporting mechanism. I therefore regard n^2 log n as the ceiling of "Poisson-like" methods, not as a believable upper bound for V_n.

## 4. Polynomial vs exponential: published statements

- Venkatesh 2013 (IMRN) [NOT accessed; 2ND]: Aubrun's Bourbaki survey arXiv:2606.13313, p.2 (FULL): "[Venkatesh, 2013] conjectures that, up to logarithmic corrections, the quantity 2^n delta_n^L is bounded above by a polynomial in n. We now know that such a polynomial must have degree >= 2." A search snippet of Venkatesh's intro reads "However, it seems to me quite likely that, in the context of lattice packings, Minkowski's bound differs ..." (truncated). Kim 2020 refers to "the discussion in Section 3.1 of Venkatesh" for the suspicion that Poisson behaviour of L ∩ B(V) persists up to V < 2n. I could not read the exact wording or the exponent Venkatesh has in mind.
- Aubrun survey (FULL): "it is believed that the values of delta_n^L and delta_n differ for all sufficiently large n."
- Torquato-Stillinger 2006 (Exp. Math. 15; arXiv:math/0508381) [FULL, Conjecture 5.1]: "For sufficiently large d, a hard-core nonnegative tempered distribution g_2(r) that satisfies g_2(r) = 1 + O(|r|^{-d-eps}) is a pair correlation function of a translationally invariant disordered sphere packing in R^d at number density rho if and only if S(k) = 1 + rho h~(k) >= 0. The maximum achievable density is the terminal density phi_*." With their test function this gives the conjectural phi_max >= 2^{-(0.77865...+o(1)) d}, i.e. 2^d delta_d >= 2^{0.2213 d}. This is for DISORDERED (non-lattice) packings; they say "for sufficiently large d, lattice packings are likely not to be the densest".
- Vanoni-Guo-Torquato arXiv:2608.00920v3 (Sep 2026) [ABS + intro]: conditional on the T-S realizability conjecture, hyperuniform test functions give phi_max >~ d^beta 2^{-d} for every fixed beta > 1; they re-derive the rate 0.7786524795 from the Cohn-Elkies dual; a Gauss-Radau construction with floor((d-1)/4) shells gives rate 0.622556...; the unrestricted pair-correlation program attains the Cohn-Elkies rate 0.6044005.... They cite Klartag's c d^2 2^{-d} as motivation. All conditional, all non-lattice.
- Parisi-Zamponi (J. Stat. Mech. 2006; Rev. Mod. Phys. 2010) [MEM, not accessed]: replica theory predicts amorphous packings with 2^d phi ~ d log d (glass close packing / Kauzmann), threshold ~ 6.26 d. This matches the order of the Campos-Jenssen-Michelen-Sahasrabudhe theorem (1/2 - o(1)) d log d 2^{-d}. It concerns amorphous packings and says nothing about lattices.
- I found NO published heuristic or theorem asserting an upper barrier V_n <= n^2 polylog(n). The only items pointing that way are (i) Venkatesh's polynomial conjecture as reported by Aubrun, and (ii) my conditional statement 3(b).

## 5. Techniques that go beyond the first moment (for existence of holes)

1. Schmidt's sieve (1958/1963; Kim 2020; Klartag 2605.10479 Sec. 3 "Prop 3.1"): Bonferroni inequalities in which only rank-deficiency-0 and -1 tuples are needed, so that Rogers' formula has no (resp. small) error term for all k <= n - 1. Gives P_V = e^{-V/2} + O(e^{-cn}); ceiling V < 2 n log 2 because of rank-1 dependencies (probability ~ V 2^{-n}).
2. Rogers 1947 / Davenport-Rogers / Ball 1992: first moment + successive minima (deform the lattice, not measure-preserving). Ceiling ~ n.
3. Symmetry (Venkatesh 2013; Gargava; Gargava-Viazovska; GSV moments): first moment with orbits of size omega_K; ceiling n log log n. GSV show the moments are Poisson(V/omega_K) — i.e. the same "Poisson world" with a different pair-size, no extra structure.
4. Klartag 2025 and AGZ/GGS 2026: first moment (Siegel) applied along a stochastic deformation; the gain comes from N = n(n+1)/2 (resp. N_G) free parameters, not from higher moments.
5. Non-lattice: Krivelevich-Litsyn-Vardy, Jenssen-Joos-Perkins (hard sphere model), Campos-Jenssen-Michelen-Sahasrabudhe (independent sets in locally sparse graphs / nibble). [MEM] No lattice analogue of the nibble was found in the literature.
6. Covering analogue [MEM, not accessed]: Ordentlich-Regev-Weiss, "New bounds on the density of lattice coverings" (JAMS 2022): lattice covering density O(n^2) via finite-field Kakeya estimates for random Construction-A lattices — a structural (non-moment) input for a random-lattice extremal problem; note the same exponent 2.
7. No paper was found using Paley-Zygmund on a count of holes, exchangeable pairs, or Chen-Stein for random lattices in growing dimension. (Klartag's TV theorem is proved by Schmidt's sieve, not Stein's method.)

## 6. Open problems / barriers extracted

- Sharp asymptotics of P_V for C n <= V <= c n^2 are completely open: only e^{-(1+o(1)) n^2 log n} <= P_V <= e^{-cn}.
- Kim's programme: improve the additive error in Schmidt's theorem from e^{-cn} to e^{-C n log n}; this requires handling tuples of rank deficiency >= 2 (lattice structure), where Rogers' formula error terms are not small; n-th moments diverge.
- Klartag's conjecture (2504.05042, p.3): Theorem 1.2 (admissible lattice for volume c n^2) should hold for every origin-symmetric convex body, not only ellipsoids; known only for volume <= cn (Schmidt 1963).
- Is P_V >= e^{-V/2} (positive association of hole events) true in general? Is P_V <= e^{-cV} true for V up to n^2 log n? Either answer would be decisive for the heuristics in Section 3.

## 7. Ideas suggested by this survey (all INFER / untested)

- (I1) Prove a "tilted" lower bound P_V >= exp(-C V) for V <= c n^2 by running Klartag's process as a change of measure on X_n (Girsanov-type cost of the drift needed to keep contacts) — would quantify how typical the n^2-holes are and whether the entropy cost n^2 log n in 3(b) is attained.
- (I2) Entropy counting as a proof strategy for an UPPER barrier is hopeless without P_V <= e^{-cV}; but as a LOWER-bound strategy it suggests: any method that only uses Poisson-like independence cannot beat (2+o(1)) n^2 log n. A target "n^2 log n" theorem would be the natural endpoint of the Klartag-type approach (compare CJMS's n log n vs n for non-lattice packings).
- (I3) To get the missing log n one needs contacts to be absorbed at rate higher than one per degree of freedom, i.e. the ellipsoid should keep a hole while ~N log n Poisson points pass — impossible for an ellipsoid on one lattice (only N parameters), so the extra factor must come from choosing among e^{N log n} "independent" lattices/ellipsoids (a union/second-moment argument over deformations), or from symmetry (orbit size).
- (I4) Positive-association question: test numerically in small n whether P(L ∩ (S_1 ∪ S_2) = 0) >= P(L ∩ S_1 = 0) P(L ∩ S_2 = 0) for symmetric S_i.
