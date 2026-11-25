# R4_perceptron_theory -- the "random polytope / negative perceptron" model of the Ryshkov polyhedron

Start 07:23 UTC, end ~07:41 UTC (25 Nov 2026). Tags: [P] proved (complete argument here), [PM] proved as a statement about the
Gaussian/Poisson MODEL, [C] conditional, [H] heuristic, [N] numerical (code in this directory), [L] literature, [R] refuted/corrected.
Code: `lp_model.py` (LP in the Gaussian model), `lp_poisson_exact.py` (exact trace-LP for a Poisson process in R^n, row generation),
`check_formulas.py` (quadrature / Monte Carlo checks of every closed formula used below). Logs: `lp_run1.log`, `lp_run2.log`, `lp_run_exact.log`.

## 0. Verdict on the orchestrator's claims

| claim | verdict |
|---|---|
| C1 (trace-LP: v* = (2/e)N, rho^2 = 8, N active, V = 2N ~ n^2) | all constants re-derived and correct [PM for the AO computation; C for the w.h.p. statement]; LP numerics: V/N = 2.16, 2.07 (N = 210), 2.00 (N = 300) [N]; exact Poisson LP in R^n, n = 8..22, gives V/N = 1.86..2.58 (mean ~2.1) [N] |
| C2 (min-norm path: v = N rho^2 e^{-rho^2/8}/4, c = rho^2 N/8, V = rho^2 N/4 = 2c) | correct algebra given C1 [C] |
| C3 (sticky process freezes at V ~ N; rigorous analysis loses e) | correct; closed formula c(Q) = v(e^{Q/8} Phi(sqrt(Q)/2) - 1/2), V = v e^{Q/8}; the factor e is exactly the un-projected clock [H, formula N-checked]; agrees with R3_process 1.3/1.4 (independent re-derivation) |
| C4 (existence threshold V ~ c N log rho) | correct in structure and the constant can be named: V_exist(rho) = (2+o(1)) N log rho, i.e. (existence)/(LP) = log(rho/2)(1+o(1)), exactly MZZ's delta_s/delta_lin = log|kappa| with kappa = -rho/2 [L + PM pair formula + C second moment]. Two corrections [R]: (i) the plain second moment works only if v V / rho^2 -> 0 (tiny reference ball, rho^2 >= (16+o(1)) log N), not for v ~ 1; (ii) "decorrelation at distance O(1)" is true only in the weak sense needed for the entropy count; a fixed pair at distance delta is decorrelated only by the factor exp(-V(Phi(delta/4)-1/2)). |
| C4, "rho^2 ~ n gives n^2 log n in the model" | the Gaussian model is only valid for rho << n^{1/4} (quartic term of log det); within that range the model gives V up to ~ (1/4 - eps) n^2 log n, still order n^2 log n [H]. The un-linearised Poisson model has sup V = infinity (needles) [P], so the model MUST be cut off; the lattice cut-off is the geometry of SL_n(R)/SL_n(Z) (R3_process 3.1-3.3). |
| C5 (transfer) | first moment transfers only as far as P(L cap E = 0) <= e^{-cV}, which is NOT known (and false for V >> n^2 log n if V_n is larger); see 5. |

## 1. The model, re-derived

Lattice L of covolume 1, ball B_R of volume v = kappa_n R^n, ellipsoid A = R^{-2}(I + B_0), tr B_0 = 0, b = |B_0|_HS.

**1.1 Exact linear form of the constraints [P].** For x = R(1+eps)u, |u| = 1: x^T A x >= 1 iff u^T B_0 u >= -s with s := 1 - (1+eps)^{-2}
(exactly linear in B_0; s = 2 eps + O(eps^2)). So the Ryshkov polyhedron in these coordinates is {B_0 : u_x^T B_0 u_x >= -s_x for all x}.

**1.2 Intensity [P given the Poisson assumption].** #pairs with slack <= s = (v/2)(1+eps)^n = (v/2)(1-s)^{-n/2}. With sigma := n s this is
(v/2)(1 - sigma/n)^{-n/2} -> (v/2) e^{sigma/2}; intensity (v/4) e^{sigma/2} d sigma on the whole line (sigma < 0: the v/2 pairs inside B_R).
Relative error of the exponential form: exp(sigma^2/(4n) + ...), negligible iff sigma << sqrt(n) (R3_process 1.1 saw this factor).

**1.3 Variance [P].** For u uniform on S^{n-1}: E u_i^2 u_j^2 = 1/(n(n+2)) (i != j), E u_i^4 = 3/(n(n+2)), hence for symmetric B:
E (u^T B u)^2 = (2|B|^2 + (tr B)^2)/(n(n+2)); for traceless B_0: Var(u^T B_0 u) = 2 b^2/(n(n+2)). [N: MC 0.910522 vs 0.910502, n = 12.]

**1.4 Coordinates.** theta := sqrt(2) B_0 in the (N-1)-dimensional space of traceless forms with the HS norm, rho := |theta| = sqrt(2) b,
g_x := (n/sqrt 2)(u u^T - I/n). Then <theta, g_x> = n u^T B_0 u and the constraints are <theta, g_x> >= -sigma_x. Exact facts [P]:
E g g^T = (n/(n+2)) Id on traceless forms (isotropy: irreducibility of the O(n)-action), and |g|^2 = n(n-1)/2 = N - n deterministically.

**1.5 Volume [P].** V = v det(I+B_0)^{-1/2} >= v (concavity of log det, equality iff B_0 = 0); V = v exp(b^2/4 - tr B_0^3/6 + tr B_0^4/8 - ...)
= v exp(rho^2/8 + O(tr B_0^3)). For GOE-like B_0 the quartic term is ~ b^4/(4n) = rho^4/(16 n): the quadratic formula needs rho << n^{1/4}.
Exact first moment: E #{pairs in the ellipsoid} = V/2 (Siegel for lattices; for Poisson the number is Poisson(V/2)).

**1.6 MODEL ASSUMPTIONS (this is what "the model" means).**
 (M1) the pairs of L form a Poisson process of intensity 1/2 (radii and directions independent, directions uniform);
 (M2) the vectors g_x (which lie on a Veronese variety of dimension n-1 inside the sphere of radius sqrt(N-n) of R^{N-1}) are replaced by
      i.i.d. standard Gaussian vectors of R^{N-1}; (M3) V = v e^{rho^2/8}; (M4) intensity exactly (v/4)e^{sigma/2}.
 (M2) is a universality assumption: one-dimensional marginals n u^T B_0 u = sum b_i (n u_i^2 - 1) are asymptotically Gaussian iff
 |B_0|_op/|B_0|_HS -> 0; (M2)-(M4) are accurate only for rho << n^{1/4} and offsets sigma << sqrt(n). Model with only (M1) = "exact Poisson model".

## 2. C1/C2: the convex problem (Gordon / CGMT)

**Theorem used [L]** (Gordon 1985 Israel J. Math. 50; Gordon 1988 GAFA Sem. LNM 1317; Stojnic 2013 for perceptrons; CGMT: Thrampoulidis-Oymak-Hassibi
COLT 2015, Thrampoulidis-Abbasi-Hassibi IEEE-IT 2018). G in R^{M x d}, g in R^M, h in R^d standard Gaussian, S_w, S_u compact, psi continuous,
 Phi(G) = min_{w in S_w} max_{u in S_u} u^T G w + psi(w,u),   phi(g,h) = min_w max_u |w| g^T u + |u| h^T w + psi(w,u).
 (i) for every c: P(Phi < c) <= 2 P(phi <= c)  -- no convexity needed;
 (ii) if S_w, S_u are convex and psi is convex-concave: P(Phi > c) <= 2 P(phi >= c).
Here w = theta, u = lambda in {lambda >= 0, |lambda| <= 1} (convex), psi = -lambda^T sigma (bilinear), Phi = min_theta |(sigma + G theta)_-| = l2-distance
to feasibility (use -G ~ G). S_w = ball {|theta| <= rho}: (i)+(ii), two-sided. S_w = sphere {|theta| = rho}: only (i), i.e. only INFEASIBILITY
(upper bounds on V) is rigorous on the sphere; existence on the sphere is the non-convex problem of section 4.

**AO computation [PM].** For |theta| = r the minimum over the direction gives h^T theta = -r|h|, and the max over lambda gives
 phi = min_{r <= rho} ( |(sigma - r g)_-| - r|h| )_+ ,   sum_i (sigma_i + r Z_i)_-^2 has mean
 S(r) = int (v/4) e^{sigma/2} E(sigma + rZ)_-^2 d sigma = (v/4) * 16 * E e^{-rZ/2} = 4 v e^{r^2/8}        [N: quadrature agrees to 7 digits]
and variance (v/4) * 768 * e^{r^2/8} = 192 V (Campbell), so sd/mean = sqrt(12/V) -> 0 when V -> infinity; |h|^2 = d(1 + O(d^{-1/2})).
Hence phi > 0 w.h.p. iff 4 v e^{r^2/8} > r^2 N for all r <= rho, and phi = 0 w.h.p. iff 4 v e^{r^2/8} < r^2 N for some r <= rho.

**Model Theorem A [C: complete up to the routine uniform-in-r concentration; every displayed identity is proved].** In the model, for fixed eps > 0, N -> infinity:
 K_v cap {|theta| <= rho} is empty w.h.p. if v > (1+eps) max_{r <= rho} N r^2 e^{-r^2/8}/4 and non-empty w.h.p. if v < (1-eps) times the same.
 In particular K_v != empty iff v < v* = (2/e) N (1+o(1)), attained at r^2 = 8.
 Passing from "Phi <= c" to exact feasibility: apply (ii) to the system with all offsets lowered by eta (this is the same model with v e^{eta/2}
 instead of v); Phi <= eta for that system gives sigma_i + <g_i,theta> >= eta - eta = 0 for the original one.
 On the sphere: K_v cap {|theta| = rho} = empty w.h.p. if 4 v e^{rho^2/8} > (1+eps) rho^2 N, i.e. **V <= (1+eps) rho^2 N/4 for every point of K_v of norm rho** (one-sided, rigorous in the model).

**Truncation of the infinite constraint set.** Infeasibility: dropping constraints enlarges K_v, so it suffices to keep sigma_i <= S; the AO sum loses
only int_S^infty (...) which is a Gaussian tail in (S - r^2/2)/r; S = S(eps, rho) constant, M = Poisson((v/2)e^{S/2}) = O(N): proportional regime.
Feasibility: solve the truncated system (sigma_i <= S) and let theta^ be its minimum-norm point; the constraints with sigma_i > S are independent of
theta^ (Poisson on disjoint sets, independent g_i), so the number of them violated by theta^ is Poisson with mean <= (V/2) Phi(-(S/rho - rho/2)), which is
o(1) once S >= rho^2/2 + rho sqrt((2+eps) log N). Then M = N^{1+o(1)}; Gordon's comparison is non-asymptotic, so nothing changes.

**Consequences (C1, C2) [C].** Number of active constraints = #{lambda_i > 0} = #{sigma_i < r g_i}, mean (v/2) e^{r^2/8} = V/2; on the feasibility
boundary 4 v e^{r^2/8} = r^2 N this is c = r^2 N/8, and V = r^2 N/4 = 2c. At r^2 = 8: c = N (= number of unknowns (theta, v): a non-degenerate LP vertex,
a consistency check), v* = (2/e)N, V = e v* = 2N = n(n+1). Active offsets have density prop. to e^{sigma/2} e^{-sigma^2/(2 r^2)}: Gaussian with mean r^2/2 = 4
and sd r = 2.83, i.e. the active lattice points lie in the ball of volume ~ v* e^{2 +- 1.4}, a few times N: beyond the proven Poisson range (c n) by a factor n.

**Numerics [N].** LP "maximise tau = 2 log v s.t. <theta,g_i> - tau >= -sigma'_i" (HiGHS), truncation sigma <= 11 + (checked: 0 violations among the
~3.5x more constraints up to sigma <= 14). Prediction: v*/N = 0.7358, rho^2 = 8, V/N = 2.
| N | seeds | v*/N | rho^2 | V/N |
|---|---|---|---|---|
| 55 | 3 | 0.97, 0.86, 0.77 | 5.97, 6.36, 7.02 | 2.04, 1.91, 1.85 |
| 105 | 2 | 0.72, 0.76 | 7.04, 7.01 | 1.73, 1.83 |
| 210 | 2 | 0.76, 0.70 | 8.31, 8.70 | 2.16, 2.07 |
| 300 | 1 | 0.70 | 8.36 | 2.00 |
(rho^2 is poorly determined because x e^{-x/8} is flat at x = 8.) (N = 300: one seed, LP time 114 s.)
Exact Poisson model (assumption M1 only; true quadratic constraints, true volume det^{-1/2}; row generation makes the LP exact for the infinite process):
| n | N | V/N per seed | V/n^2 | v/N | rho^2 | eigenvalues of A/mean |
|---|---|---|---|---|---|---|
| 8 | 36 | 2.18, 2.16, 2.24, 2.58 | 1.21-1.45 | 0.55-0.97 | 5.3-11.1 | 0.22 .. 2.8 |
| 10 | 55 | 2.02, 1.86, 1.90, 2.21 | 1.02-1.22 | 0.47-0.81 | 6.3-12.6 | 0.22 .. 2.9 |
| 14 | 105 | 2.05, 1.86, 2.29 | 1.00-1.23 | 0.56-0.88 | 7.5-10.7 | 0.28 .. 2.7 |
| 18 | 171 | 2.02, 2.00 | 1.06 | 0.69-0.76 | 7.6-9.2 | 0.36 .. 2.2 |
| 22 | 253 | 2.23 | 1.17 | 0.80 | 7.95 | 0.39 .. 1.9 |
So already for n <= 22 the trace-LP optimum of a Poisson process has V ~ 2N ~ n^2 (the ellipsoids are far from round, eigenvalue ratio ~ 5-10, so this is
a statement about the exact Poisson model, not about the Gaussian linearisation). [In two runs (one Gaussian N = 210, one Poisson n = 18) fewer than N constraints were counted as active (0.98 N and 0.80 N; tolerances 1e-7 / 1e-9 on the slack): not investigated, probably a tolerance artefact.]
NOT tested here: real lattices (for a lattice the trace-LP over the Ryshkov polyhedron is Voronoi's algorithm restricted to one linear functional; see R2_voronoi_descent).

## 3. C3: fluid model of Klartag's sticky process, and where 2 and e are lost

Free Brownian motion in theta, projected off the active g's: with c contacts, E d rho^2 = (N - c) dt, and an untouched constraint
X_i = sigma_i + <theta_t, g_i> has d<X_i> = |pi_t g_i|^2 dt. Exact average [P, also R3_process 1.4]: for any subspace E of Sym_n,
E_u |P_E(u u^T)|^2 = (2 dim E + |P_E I|^2)/(n(n+2)), so for a direction independent of the contacts |pi_t g_i|^2 = (N - c)(1+o(1)).
Common clock Q := int (N - c) dt: rho^2 = Q, V = v e^{Q/8}, each X_i is a Brownian motion in the clock Q absorbed at 0, so (reflection principle, sharp)
 c(Q) = int_0^infty (v/4) e^{sigma/2} 2 Phi(-sigma/sqrt Q) d sigma = v ( e^{Q/8} Phi(sqrt(Q)/2) - 1/2 )   [N: quadrature agrees to 5 digits]
 => c = V Phi(sqrt(Q)/2) - v/2 ~ V: one contact pair per unit volume (a fixed ellipsoid of volume V contains only V/2 pairs: the path "sweeps twice what it keeps").
Freezing at c = N: **V_final = (N + v/2)/Phi(sqrt(Q_f)/2) ~ N = n^2/2** [H]. Data of R2_klartag_sim (n = 2..10): V_final/N = 1.21, 0.99, 0.93, 0.97, 1.08, 1.26,
1.24, 1.01, 1.03 and V(t)/k(t) -> 1: consistent, small n only.
Accounting in one line each (all in units of N):
 - LP vertex (section 2): V = 2N. Every one of the N active constraints is "used" (dual variable > 0).
 - sticky process, projected clock: V = N. Loss of 2 = reflection principle: constraints touched early stay active although the free endpoint would violate only half of them.
 - Klartag/Rezvyakova accounting: hits are bounded with the UN-projected clock (|pi_t(x x^T)|^2 <= |x|^4, the only bound available for a point whose
   direction is correlated with the history): E c(y) <= v e^{y}, y = N t/8, while log(V/v) = int (1 - c/N) dy >= y - (v/N)(e^y - 1); max at e^y = N/v:
   E log V >= log N - 1, i.e. **V >= N/e = n^2/(2e)** [this reproduces Rezvyakova's constant; literature-consistent, cf. R1_lit_latest (I2), R3_process 1.3].
   Loss of e = contacts accumulate at full speed while the volume grows at the reduced speed (1 - c/N).
 So n^2 (LP) -> n^2/2 (sticky, heuristic) -> n^2/(2e) (proved).

## 4. C4: the non-convex problem (points of K_v of large norm)

**4.1 Literature [L]** Montanari-Zhong-Zhou, arXiv:2110.15824 (text in BASE/papers/mzz_negative_perceptron_2110.15824/): unit sphere of R^d, n = delta d
constraints <theta, x_i> >= kappa, kappa < 0. Eq. (7), Thm 4.1 (lower bound, second moment-type functional Psi(q)) and Thm 4.2 (upper bound, Stojnic's lifted
Gordon bound): delta_s(kappa) = (1+o_kappa(1)) log|kappa| / Phi(kappa). Eq. (9): the LP algorithm succeeds up to delta_lin(kappa) = (1+o_kappa(1))/Phi(kappa).
Replica-symmetric (plain Gordon) upper bound: delta_RS = 1/E(kappa - G)_+^2 ~ kappa^2/(2 Phi(kappa)). They state a "multiplicative gap of order log|kappa|".
**Translation.** A fixed direction violates delta d Phi(kappa) constraints; in our problem a fixed theta of norm rho violates V/2; the local exponential
tilt e^{sigma/2} of our offsets equals the Gaussian tail slope at kappa = -rho/2. Hence
 LP: V/2 = N (V = 2N, = section 2);  existence: V/2 = N log(rho/2) (V = (2+o(1)) N log rho);  plain Gordon: V/2 = N rho^2/8 (= Theorem A, sphere).
So the non-constructive gain is the factor log(rho/2): ~ (1/2) log log N in Klartag's regime rho^2 ~ 8..16 log N, and (1/2) log n if rho = n^{a} were allowed.

**4.2 Exact pair formula for our offsets [PM].** For |theta| = |theta'| = rho, delta = |theta - theta'|, Y = -<theta,g>, Y' = -<theta',g>:
#constraints violated by theta or theta' is Poisson with mean (v/2) E e^{max(Y,Y')/2}; max = (Y+Y')/2 + |Y-Y'|/2 with independent parts of variances
rho^2 - delta^2/4 and delta^2 (for the difference), and E e^{|D|/4} = 2 e^{delta^2/32} Phi(delta/4). Therefore
 **P(theta in K_v) = e^{-V/2},  P(theta, theta' in K_v) = exp(-V Phi(delta/4))**,  V = v e^{rho^2/8}   [N: MC agrees within 2 s.e. at 4 parameter pairs].
Consequences: (a) P(theta' in K | theta in K) = exp(-V(Phi(delta/4) - 1/2)) ~ exp(-V delta/(4 sqrt(2 pi))): a FIXED neighbour is feasible only for delta <~ 1/V;
(b) a neighbour at distance delta violates only Poisson(V(Phi(delta/4)-1/2)) constraints, a fraction ~ delta/5 of V/2.

**4.3 First moment [H].** Z = normalised measure of K_v on the sphere of radius rho, E Z = e^{-V/2}. By (b), if a feasible point exists then a
constant fraction of its delta-neighbourhood violates <= m = V delta/5 constraints, and P(Poisson(V/2) <= m) = exp(-(V/2)(1 - O(delta log(1/delta)))).
With (C rho/delta)^N cells: no feasible point if V > (2+eps) N log rho (rho large, delta a small constant). Not a proof (the feasible point depends on the
data; one needs a uniform bound on #{i: 0 <= X_i(theta) <= x} over the sphere). The rigorous analogue for a common margin is MZZ Thm 4.2.

**4.4 Second moment [C: Laplace estimates sketched, identities proved].** E Z^2/(E Z)^2 = E exp(V Phi(-delta/4)), delta = rho sqrt(2(1-q)), q the overlap of two
uniform points (density prop. to (1-q^2)^{(N-4)/2}).
 - near-diagonal: need V Phi(-delta/4) < N log(rho/delta) for delta = O(1): holds iff V <= (2 - eps) N log rho (rho large);
 - bulk (q ~ 0): need V Phi(-rho/(2 sqrt 2)) = (2/(rho sqrt pi)) sqrt(v V) (1+o(1)) -> 0, i.e. **v V << rho^2**.
 If both hold, E Z^2/(EZ)^2 -> 1 and K_v meets the sphere w.h.p. The bulk condition fails for v ~ 1 (constraints with small sigma are violated by a positive
 fraction of all directions and make Z non-concentrated) [R: correction to C4]; it holds for v = o(rho^2/(N log rho)), which forces rho^2 >= (16+o(1)) log N.
 **Model statement:** for rho -> infinity with rho^2 >= 16 log N (and rho << n^{1/4} for the model to make sense), the largest V on the sphere of radius rho is
 (2+o(1)) N log rho; at rho^2 ~ 16 log N this is V = (1+o(1)) N log log N = (1/2 + o(1)) n^2 log log n -- the same order as AGZ/GGS, by a different mechanism
 (no symmetry; pure existence). For v ~ 1 a conditional second moment (conditioning on the points with sigma <= sigma_0) should give the same [H, not done].

**4.5 Limits of the model.** [P] In the exact Poisson model sup V = infinity: needles of length l and width r have ~ (l/r)^{n-1} independent directions, each empty
with probability e^{-V/2}. So "N log rho" must be cut off by hand; in the Gaussian model K_v is bounded (|theta| <= C sqrt N) and the formulae are valid for
rho << n^{1/4}, giving at most V ~ (1/4) n^2 log n (1+o(1)) [H]. For lattices all ellipsoids of volume V are equivalent under SL_n(R), and "exists an empty ellipsoid of
volume V for L" is the deterministic statement V <= V_n: the model can only be a model of a neighbourhood of a typical lattice in SL_n(R)/SL_n(Z).

## 5. What can be claimed, conjectured, and what is missing

**(a) Model theorems.** Theorem A (two-sided on balls, one-sided on spheres) [C]; pair formula 4.2 [PM]; S(r), c(Q) formulae [PM]; second-moment threshold 4.4 [C];
sup V = infinity in the exact Poisson model [P]. Numerical: trace-LP optimum of a Poisson process in R^n has V/N ~ 2 for n = 8..22 [N].

**(b) Conjectures for lattices.**
 (LP) For Haar-random L (and after Klartag's reduction, for every n) the minimiser of Tr A over the Ryshkov polyhedron of L has
      vol = (2 + o(1)) N = (1+o(1)) n^2, with trace-equivalent ball (2/e + o(1)) N. Would give V_n >= (1 - o(1)) n^2 (constant 1 vs proved 1/(2e)) [conjecture].
 (EX) V_n >= c n^2 log n [conjecture, from 4.4 extrapolated to rho = n^{a}]; together with the first-moment heuristic, V_n = n^2 (log n)^{1+o(1)} would be the
      prediction of the model, i.e. both the truth of Venkatesh's upper conjecture and its sharpness up to the power of log [heuristic].

**(c) Missing lemmas.**
 (LP-1) "Projected clock": for the lattice points x that matter (in the ball of volume ~ C N), |pi_t(x x^T)|^2 <= (1+o(1))(1 - c_t/N)|x|^4 on average over the
        Haar measure, jointly with the history of the process. For a FIXED subspace this is the identity in section 3; with the history it is a statement about
        correlations of ~N lattice points, i.e. needs Rogers-type k-point formulae with k ~ N or a martingale substitute. This would give n^2/2.
 (LP-2) To reach n^2 one must also release contacts (LP vertex vs sticky path): e.g. a reflected (non-sticky) process or a drift (R3_process 2(ii),(iii)); the
        first-moment input needed is the same as in (LP-1), plus control of how often a released point is re-hit.
 (LP-3) Primal-dual route: weak duality reads Tr A >= sum_x lambda_x whenever sum_x lambda_x x x^T = I, lambda >= 0 (eutaxy-type certificate). An approximate radial
        certificate lambda_x = f(|x|) is controlled by Rogers' second-moment formula (variance of a LINEAR statistic is Poisson-like for all volumes), but it only
        proves v <= (e/8) rho^2 N in model units (computed here: optimum cut-off sigma <= -2), i.e. an UPPER bound for the LP and only for near-round A. Duality
        does not produce feasible points; the primal side needs (LP-1).
 (EX-1) Pair estimate (C5): P(L avoids E u E') <= K P(L avoids E) P(L avoids E') exp((1+o(1)) vol(E cap E')/2) for ellipsoids of volume ~ N log rho
        inside a ball of volume V^2/v ~ N^3 (the range of the relevant offsets sigma ~ rho^2/2), together with the one-ellipsoid lower bound
        P(L avoids E) >= e^{-(1+o(1))V/2}. Known: total-variation Poisson only for volume <= c n (Klartag 2026). Note that by 4.4 the bulk pairs already need
        vol(E cap E') statistics to be exact to O(1), not to a factor.
 (EX-2) Upper side: P_Haar(V(L) >= V) <= e^{-cV} for V up to n^2 polylog n. By R3_process Lemma N, P(V(L) >= V_n/2) >= e^{-C n^2 log n}; so this single tail
        bound at V = V_n/2 would PROVE V_n <= C' n^2 log n. It is the precise form of the Poisson heuristic behind Venkatesh's conjecture.

**(d) Venkatesh's n^2 polylog.** In the model the largest empty ellipsoid has V = 2 x (entropy of the family of shapes) = 2 N log(rho_max/2); every bounded-eccentricity
cut-off rho_max = n^{a} gives Theta(n^2 log n). "Algorithmic" constructions (LP, Klartag's process, also with symmetry up to the factor m/phi(m) ~ log log n) sit at
Theta(N), a factor log(rho) below existence: this is the MZZ statistical-to-computational gap, and it suggests that process-based proofs cannot reach n^2 log n,
while a second-moment proof would need (EX-1). The log log n of AGZ/GGS and the log log N of 4.4 at rho^2 ~ 16 log N are numerically the same order but come from
different sources (symmetry reduces the number of constraints; here entropy of directions).

## 6. Not done / caveats
- Checked numerically after writing: the radial-certificate constant e/8 = 0.33979 (argmin S = -2) and the bulk asymptotic Phi(-rho/(2 sqrt 2)) ~ 2 e^{-rho^2/16}/(rho sqrt pi) (ratio 0.76, 0.86, 0.91 at rho = 4, 6, 8, converging slowly).
- Theorem A and 4.4 are written at the level "identities proved, uniform concentration sketched"; they are tagged [C], not [P].
- Universality (M2) for the Veronese vectors is untested beyond the exact-Poisson LP runs (which test the LP optimum only, n <= 22, 1-4 seeds per n).
- No lattice computation in this round; the LP conjecture (b) should be tested on real random lattices (n ~ 10-16) by minimising Tr(A) over the Ryshkov polyhedron.
- The constant in the first-moment bound 4.3 (2 N log rho) is heuristic; only its MZZ analogue is a theorem.
