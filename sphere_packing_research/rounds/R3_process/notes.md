# R3_process -- how far can Klartag-type evolving-ellipsoid arguments go without symmetry?

Agent notes, 2026-11-25 (budget ~18 min). Code: `meanfield.py` (same directory; all numbers below come from it).
Units (Klartag): lattice covolume kappa_n = Vol(B^n); relative volume of E_A is V = det(A)^(-1/2); D = n(n+1)/2 = dim Sym_n.
In these units a Poisson-like lattice has intensity 1 point (1/2 pair) per unit relative volume, and V_n = sup V.
Tags: [proved] [conditional] [conjecture] [numerical] [heuristic] [refuted] [literature].

## 1. Mean-field theory of Klartag's process

### 1.1 Contact count with the sharp reflection principle  [heuristic, formula checked numerically]
For a lattice point x = |x| w define its depth s = a_0 - 1/|x|^2 >= 0. x becomes a contact when the martingale
A_t w.w (start a_0, quadratic-variation rate |pi_t(w(x)w)|^2 <= 1) first reaches a_0 - s. Number of points with depth <= s is
(a_0 - s)^(-n/2) = V_0 exp(ns/2 + O(n s^2)), V_0 = a_0^(-n/2). With clock tau (= t if the variance is bounded by 1 as in Prop 4.1):
  K(tau) = int_0^inf V_0 (n/2) e^{ns/2} * 2 P(sqrt(tau) Z >= s) ds = 2 V_0 E[(e^{n sqrt(tau) Z/2} - 1) 1_{Z>0}]
         = V_0 ( 2 e^u Phi(sqrt(2u)) - 1 ),   u := n^2 tau / 8,
so the number of contact PAIRS is c(u) = V_0 ( e^u Phi(sqrt(2u)) - 1/2 ) ~ V_0 e^u.
Check: closed form vs quadrature with the exact density (n/2)(a_0-s)^(-n/2-1): n=1000: u=0.5: 1.774 vs 1.781; u=2: 13.44 vs 13.74;
u = 2 log n: 2.0e6 vs 4.8e6. The last discrepancy is the neglected factor exp(n s^2/4) at the saddle s* = 4u/n, i.e. exp(4u^2/n) = exp(16 log^2 n / n):
it tends to 1 but only for n >> 16 log^2 n  (finite-n loss; it disappears for the affine-invariant covariance of Remark 5.3, see 2(i)).

### 1.2 Volume equation  [heuristic]
Ito (Klartag Lemma 3.3 with lambda_i -> 1): d E log V = (1/4) N_t dt, N_t = dim F_t = D - c_t (contacts in general position).

### 1.3 Model K = the accounting of Klartag / Rezvyakova (variance bounded by 1, clock tau = t)  [literature-consistent]
 d log V / du = (2D/n^2)(1 - c(u)/D) ~ 1 - V_0 e^u / D  =>  log V(u) = log V_0 + u - V_0 e^u / D, maximal at e^u = D/V_0:
 V_max = D/e = n(n+1)/(2e), i.e. exactly Rezvyakova's constant 1/(2e) = 0.1839, independent of V_0.
 Numerics (meanfield.py (3)): max_t V/n^2 = 0.2024 (n=100), 0.1866 (1e3), 0.1839 (1e6); at Klartag's T_K = 16 log n/n^2 with V_0 = (1-1/n)^n: 0.1763 n^2.

### 1.4 Model R = self-consistent variance reduction  [heuristic; key identity proved]
Exact identity [proved; Monte Carlo check in meanfield.py (2), n=8, agreement to 3 digits]: for any subspace E of Sym_n and w uniform on S^{n-1},
  E |P_E(w(x)w)|^2 = (2 dim E + |P_E Id|^2) / (n(n+2)).
Proof: E[(w(x)w)(x)(w(x)w)] acts on Sym_n as (2 Id_Sym + Id (x) Id)/(n(n+2)) (4th moments of the sphere: E w_i w_j w_k w_l = (d_ij d_kl + d_ik d_jl + d_il d_jk)/(n(n+2)));
take the trace against P_E.
For E_t = span of c contact tensors in generic directions, |P_E Id|^2 ~ cn/(n-1+c) (-> n once c >> n: the trace of A_t freezes), so for a generic direction
  |pi_t(w(x)w)|^2 ~ N_t / D = 1 - c_t/D.
Hence every not-yet-touched point runs on the slower clock d tau = (1 - c/D) dt, while d log V = (1/4)(D - c) dt = (D/4) d tau. Therefore, in the clock u = D tau/4:
  V = V_0 e^u,  c = V_0 (e^u Phi(sqrt(2u)) - 1/2) ~ V      ("one contact pair per unit of relative volume": dV = dc),
and the process freezes (c = D) at
  **V_final ~ D + V_0/2 = n(n+1)/2 (1+o(1)), i.e. constant 1/2**  (not 1/(2e)).
In real time this is the logistic equation V' = (D/4) V (1 - V/D): V(t) = D / (1 + (D/V_0 - 1) e^{-Dt/4}).
At Klartag's T_K: V(T_K) = n^2/3 if V_0 = 1, and 0.212 n^2 for his V_0 = (1-1/n)^n ~ 1/e  (meanfield.py (3),(4)).
Numerics (model R): V_final/D = 1.67 (n=10), 1.27 (30), 1.10 (100), 1.014 (1e3), 1.000 (1e6): slow convergence, finite-n corrections are large.
Comparison with the simulation of agent R2_klartag_sim (their results_table.md, n = 2..10, real lattices, run to freezing):
 V_final/D = 1.21, 0.99, 0.93, 0.97, 1.08, 1.26, 1.24, 1.01, 1.03  [numerical, small n only] -- consistent with V_final ~ D, no trend; this is weak evidence
 because n <= 10 is far from the mean-field regime (and those frozen forms have more than D contact pairs: 60.6 vs 45 at n=9).
Why the factor 2 relative to a "typical vertex" (1.5): a Brownian path sweeps twice the volume it keeps (reflection principle).
What is missing for a proof of (1/2 - o(1)) n^2: a rigorous version of "|pi_t(x(x)x)|^2 <= (1+o(1)) (N_t/D)|x|^4 for the lattice points that matter".
For a FIXED subspace this is exactly GGS Lemma 4.1 / AGZ Lemma 4.2; here E_t depends on the lattice, and Klartag's per-point bound [M]_t <= t|x|^4 cannot see it.
In a Poisson model it would follow from the stopping-set property (unexplored points are independent of the past); for Haar lattices one needs a
two-point (Rogers k=2) input. [conjecture: Klartag's process frozen at tau_* has E log(V/n^2) -> log(1/2).]

### 1.5 Poisson model: law of the volume of a vertex  [conditional on the Poisson model; Mecke + scaling]
Vertices of the Ryshkov polyhedron of a Poisson process <-> D-tuples of pairs whose circumscribed ellipsoid is empty. By the Mecke formula and the
scaling x -> s x (V -> s^n V, product measure -> s^{nD}), the intensity of vertices of a given shape with volume in dV is prop. to V^{D-1} e^{-V/2} dV:
 V ~ Gamma(D, scale 2): mean 2D = n(n+1), s.d. 2 sqrt(D), and P(V >= 2 lam D) = P(Poisson(lam D) < D) ~ exp(-D(lam - 1 - log lam)).
So a typical vertex has V ~ n^2 (constant 1); reaching V = lam * n^2 costs probability exp(-D(lam-1-log lam)); with exp(kappa D log n) vertices the
first-moment extreme is V ~ 2D lam, lam - 1 - log lam = kappa log n, i.e. V ~ kappa n^2 log n (1+o(1)) (meanfield.py (5): kappa=1: 1.68 n^2 log n at n=100, 1.28 at n=1e6).

## 2. Modifications

### 2.0 A barrier for ALL sticky exploration processes in the Poisson model  [conditional: Poisson model + standard compensator time-change theorem, details not written]
Let the points be a Poisson process (1/2 pair per unit relative volume), and let (E_t) be ANY continuous ellipsoid-valued process, adapted to the
exploration filtration (it only knows the points it has touched, plus independent randomness), sticky (a touched pair stays on the boundary), started
from an empty E_0 of relative volume V_0. Let R_t = union_{s<=t} closure(E_s) (swept region). Every point of R_t \ E_0 is a contact, so
c_t = N(R_t \ E_0) and, R_t being an increasing stopping set with continuous compensator Lambda_t = relvol(R_t \ E_0)/2, c_t = N~(Lambda_t) for a unit-rate
Poisson process N~. Since contacts in general position satisfy c_t <= D, Lambda_T <= S_{D+1} ~ Gamma(D+1,1), and V_T <= relvol(R_T) = V_0 + 2 Lambda_T:
  **V_T - V_0 is stochastically dominated by 2 Gamma(D+1,1):  E V_T <= V_0 + n(n+1) + 2,  P(V_T >= V_0 + 2 lam D) <= P(Poi(lam D) <= D) <= exp(-D(lam-1-log lam)).**
(Numerical check of the last inequality: D=55, lam=2: 4.9e-9 <= 4.7e-8.) Consequences:
 - no sticky process -- whatever its covariance, drift, or adaptivity -- beats (1+o(1)) n^2 in the Poisson model, except with probability exp(-c D);
   n^2 log n would need probability exp(-D log n (1+o(1))), i.e. as many independent restarts as there are vertices;
 - Klartag's Brownian process realises half of the cap (1.4) because relvol(R_T) ~ 2 V_T.
Not proved for Haar lattices (no stopping-set independence). Lattices with symmetry evade it only by violating "general position" (an orbit of m
contacts costs one dimension), which is exactly the AGZ/GGS mechanism.

### (i) Time-dependent / anisotropic covariance, X -> pi_t(A_t^alpha X)  [heuristic: constant only, asymptotically not even that]
 - A deterministic or adapted scalar time change does not change the law of the path: nothing.
 - General covariance Q_t (PSD on Sym_n, range in F_t): gain rate (1/4) Tr(Q_t H_{A_t}), H_A(B) = A^{-1} B A^{-1}; hitting variance of direction w is <Q_t w(x)w, w(x)w>,
   whose spherical average is (2 Tr Q + <Q Id, Id>)/(n(n+2)) (same identity as 1.4). For traceless Q the ratio gain / average hitting variance is the same
   D/4 as for the isotropic process; variance in the direction Id is n/2 times more expensive. So an oblivious anisotropy gives no gain at leading order.
 - A_t^alpha: ||A_t - Id||_op = O(sqrt(log n / n)) -> 0, so the limit constant is unchanged; the affine-invariant choice dA = A^{1/2} pi(dW) A^{1/2}
   makes log(A w.w) an exact Brownian motion with drift and removes the exp(4u^2/n) finite-n loss of 1.1 (useful for simulations at n ~ 10-1000).
 - Non-oblivious covariance (steering away from visible lattice points) is capped by 2.0 anyway.

### (ii) Drift favouring volume growth  [heuristic: constant only, at most 1/2 -> 1]
 dA = pi_t(dW) - k pi_t(A^{-1}) dt. Accounting: d log V = (1 - c/D) du_noise + du_drift and c <= V_0 e^{u_noise + u_drift} * (2 Phi-factor, -> 1 for pure drift):
 a pure expansion sweeps exactly the volume it gains, so dc = dV/2 instead of dc = dV: at most a factor 2 (consistent with the cap 2D of 2.0).
 But |pi_t(Id)|^2 = n(n-1)/(n-1+c): the deterministic growth direction dies once c >> n, and the residual gradient flow on a face has normal velocity
 (n +- sqrt(2c))/(n+c) -- still mostly expanding while c << D, so efficiency ~1 until c ~ D/10.  Order n^2 unchanged.
 Proof obstacle: Prop 4.1 becomes P(sup(B_s + drift) >= depth), fine for bounded drift, but the gain term needs a lower bound on |pi_t(A^{-1})|^2.

### (iii) Releasing contacts (reflected BM on the Ryshkov polyhedron, Langevin with density det^{-beta})  [heuristic: this is NECESSARY for more than n^2, but not sufficient by any available argument]
 - Reflected BM: dA = dW + sum_x x(x)x dl^x; the local-time terms increase det, d log det = ... + sum_x (A^{-1}x.x) dl^x, and R_L widens towards large A
   (s A in R_L for s >= 1), so without drift the reflected process drifts to SMALL ellipsoids (entropy). The sign in Lemma 3.3 is lost.
 - Langevin with stationary density prop. to det(A)^{-beta} on R_L (integrable at infinity iff n beta > D): the stationary law is a Gibbs measure whose
   volume marginal is prop. to (density of states at V) * V^{2 beta - ...}. With the Gamma(D,2) density of states of 1.5 the Gibbs measure sits at
   lam = 1 + 2 beta/D, so beta ~ (D/2) log n would give n^2 log n -- IF the density of states is still positive there. Positivity of the density of
   states at V is the statement P_V > 0, i.e. the theorem to be proved. The dynamics gives an algorithm, not an existence proof.
 - Voronoi hill-climbing (release the contact with negative eutaxy coefficient rho_i; log V is convex along each edge so the volume increases along the
   whole edge): per step, Delta V ~ 2 Exp(1) |Tr B_i| / (n E(B_i w.w)_-) ~ c/sqrt(n) for a generic edge direction, and it stops at an extreme form.
   Generic local maxima have V = Theta(n^2) (conditioning on eutaxy is an exp(-cD) event, which moves a Gamma(D,2) variable by a constant factor only).

### (iv) Non-typical starting lattice  [heuristic]
 Conditioning on an empty ball of relative volume V_1 <= cn (probability e^{-V_1/2}, Klartag's Poisson theorem arXiv:2605.10479) only adds V_1/2.
 A multiplicative gain needs the point intensity in the whole shell up to the target volume to be reduced. Symmetry does this: for a finite group G
 acting freely with image algebra B in End(V), dim of invariant forms q = n^2/(2 dim_Q B) and orbits have size |G|, so V ~ |G| q = (n^2/2) |G| / dim_Q B.
 Cyclic G = mu_M: |G|/dim_Q B = M/C (GGS), = M/phi(M) ~ log log for a single cyclotomic field.

### (v) "Only the SPAN of the contact tensors matters"  [heuristic]
 N_t = D - rank{x(x)x}. Extra contacts are free iff x(x)x lies in the span of the others AND A x.x = 1 on the face; in a Poisson model this has
 probability 0; for a generic lattice face there is no reason for it (the parallelogram relation (x+y)^2 + (x-y)^2 = 2x^2 + 2y^2 is incompatible with all
 four being contacts). Symmetric lattices are exactly the systematic way of producing dependent contacts (one condition per orbit). So (v) = (iv).

## 3. The barrier question

### 3.1 Lemma N (neighbourhood lower bound on the hole probability)  [proved]
Notation: X_n = covolume-one lattices, mu Haar probability, P_V = mu{L : V(L) >= V}, V(L) = kappa_n lambda_1(L)^n.
**Lemma N.** For n >= 2 and 0 < eps < 1/3:  P_{(1-eps)^n V_n} >= n (eps/(2+2eps))^{n^2} kappa_{n^2} / (zeta(2) zeta(3) ... zeta(n)),
kappa_{n^2} = pi^{n^2/2}/Gamma(n^2/2+1). In particular (eps = 1/(2n))  **P_V >= exp(-2 n^2 log n - n^2) for all V <= V_n/2.**
Proof.
 (a) sup V(L) is attained at some L_0 (Mahler compactness, continuity of lambda_1). The minimal vectors of L_0 span R^n: otherwise they lie in a
     hyperplane y^perp; let g_s = e^{s} on y^perp and e^{-(n-1)s} on R y (det 1). Minimal vectors become longer by e^s; all other nonzero vectors have
     length >= (1+eta) lambda_1 for some eta > 0 (discreteness) and shrink by at most e^{-(n-1)s}; for small s > 0, lambda_1(g_s L_0) >= e^s lambda_1(L_0): contradiction.
 (b) U = {g in SL_n(R): ||g - 1||_op <= eps}. For g in U, sigma_min(g) >= 1 - eps, so lambda_1(g L_0) >= (1-eps) lambda_1(L_0), V(g L_0) >= (1-eps)^n V_n.
 (c) g -> g L_0 is injective on U: if g_1 L_0 = g_2 L_0, h = g_2^{-1} g_1 preserves L_0 and ||h - 1||_op <= 2 eps/(1-eps) < 1; for a minimal vector v,
     hv - v in L_0 and |hv - v| < |v| = lambda_1, so hv = v; minimal vectors span, so h = 1.
 (d) Let m(S) = Leb_{n^2}{t g : g in S, 0 <= t <= 1} (a Haar measure on SL_n(R)). Siegel: m(SL_n(R)/SL_n(Z)) = zeta(2)...zeta(n)/n [literature].
     By (c) and unimodularity, mu(U L_0) = m(U) n / (zeta(2)...zeta(n)).
 (e) Put delta = eps/(2+eps). If ||X - 1||_op <= delta then all singular values of X lie in [1-delta, 1+delta], det X > 0, t = det(X)^{1/n} in [1-delta, 1+delta],
     and g = X/t satisfies ||g - 1|| <= 2 delta/(1-delta) = eps. So {||X-1||_op <= delta} is contained in {t g: g in U, 1-delta <= t <= 1+delta}, whose
     Lebesgue measure is ((1+delta)^{n^2} - (1-delta)^{n^2}) m(U). The left side has measure >= delta^{n^2} kappa_{n^2} (||.||_op <= ||.||_HS).
     Hence m(U) >= (delta/(1+delta))^{n^2} kappa_{n^2} = (eps/(2+2eps))^{n^2} kappa_{n^2}.
 (f) eps = 1/(2n): (1-eps)^n >= 1/2 (Bernoulli), eps/(2+2eps) = 1/(4n+2) >= 1/(5n), Gamma(z+1) <= z^z, prod zeta(k) < e:
     log P >= -n^2 log(5n) - n^2 log n + (n^2/2) log(2 pi) - 1 >= -2 n^2 log n - 0.70 n^2 - 1 >= -2n^2 log n - n^2 (n >= 2).  QED
 Numerical evaluation of the exact bound (meanfield.py): n=8: -269.4 (simplified: -330.2); n=24: -3655.6 (-4237.1); n=100: -91828 (-102103).
 Remark [heuristic]: using SO(n)-invariance and the op-norm ball of Sym_n one expects the true local exponent to be about -(1+o(1)) n^2 log n.
**Corollary N' [proved].** If P_V <= exp(-cV) holds for V = V_n/2, then V_n <= (4 n^2 log n + 2 n^2)/c. Equivalently: if V_n >= n^{2+eps}, the hole
 probability is drastically non-Poissonian, P_V >= exp(-2n^2 log n - n^2) >> e^{-V/2} for 5 n^2 log n <= V <= V_n/2. A Poisson-type upper bound for the hole
 probability up to V_n/2 would PROVE Venkatesh's conjectured upper bound V_n <= n^2 polylog(n).

### 3.2 Proposition U (union bound over near-round ellipsoids)  [proved; the hypothesis needed to use it is NOT available]
**Prop U.** Let 0 < rho < 1 and L Haar-random of covolume kappa_n. For every V > 0,
  mu{ exists L-free A with ||A - Id||_op <= rho and det(A)^{-1/2} >= V } <= (1 + 4 rho n/(1-rho))^D * P_{V/sqrt(e)}.
Proof. Take an eta-net N of {||A-Id||_op <= rho} in (Sym_n, op), eta = (1-rho)/(2n), |N| <= (1+2rho/eta)^D. Given an L-free A choose A'' in N with
||A-A''|| <= eta and put A' = A'' + eta Id. Then A' >= A, so E_{A'} is contained in E_A and is L-free; A' <= A + 2 eta Id <= (1 + 1/n) A (as A >= (1-rho) Id),
so det A' <= e det A and relvol(E_{A'}) >= V/sqrt(e). For fixed A', P(E_{A'} L-free) = P_{relvol(E_{A'})} <= P_{V/sqrt e} by SL_n-invariance of Haar measure
and monotonicity. Union bound. QED
So statement 3(a) ("w.h.p. every L-free ellipsoid with ||A-Id|| <= rho has V <= C n^2 log n") follows IF P_V <= exp(-D log(1+4 rho n/(1-rho)) - omega(1)) at
V = C n^2 log n / sqrt(e), e.g. from a Poisson-type bound P_V <= exp(-cV).
What is available [literature]: Klartag's Poisson approximation (TV <= C e^{-c'n} for volume <= cn) gives only P_V <= C e^{-c'n} for V >= cn; Rogers' moments
(k <= n-1) give at best P_V <= (Cn/V)^{cn}. Both are exp(-O(n log n)), against the needed exp(-(n^2/2) log n). A moment problem with k <= n-1 moments can
never certify hole probabilities below ~ (k/V)^k, so Rogers' formulas cannot give 3(a). By Corollary N', P_V <= e^{-cV} is FALSE beyond V ~ n^2 log n
unless V_n = O(n^2 log n). Statement 3(a) is therefore [conjecture], equivalent in strength to a hole-probability upper bound that nobody has.
Note the first-moment threshold of Prop U depends on rho only through log(rho n): rho ~ 1/n gives n^2, rho ~ sqrt(log n/n) (Klartag's a-posteriori
range) already gives (n^2/2) log n, rho ~ 1 gives n^2 log n.

### 3.3 Is n^2 log n where generic methods stop?  [heuristic]
Yes, in the following sense. Lower bounds on V_n ARE statements P_V > 0. (1) Moment methods (Rogers, Schmidt, Klartag-Poisson) certify P_V ~ e^{-V/2} > 0
only up to V ~ n. (2) Sticky adaptive search (Klartag) reaches Theta(n^2) and, in the Poisson model, cannot exceed n(n+1) (2.0). (3) The count of
"shapes" (vertices of the Ryshkov polyhedron near Id, exp(Theta(n^2 log n)) of them; van Woerden's upper bound on perfect forms is of the same type)
times the Poisson hole probability e^{-V/2} crosses 1 at V ~ n^2 log n: this is simultaneously the first-moment threshold of Prop U and the scale of Lemma N.
Above n^2 log n any lattice must be non-Poisson for a structural reason.
Can P_V >= exp(-CV) be proved for V <= c n^2 by Klartag's method? Not as it stands: his output lattice S_T L is supported on strata of codimension c_T ~ D
(the contacts), and thickening the support back to positive Haar measure costs exp(-c_T log n) -- no better than Lemma N. What would be needed: a comparison
of the law of S_T L on each stratum with the natural (Haar-induced) measure of that stratum, i.e. an entropy bound for the sticky diffusion.

## 4. Most promising routes and the exact missing lemmas

A. (Generic lattices; order n^2 log n.) By 2.0 one must give up stickiness, so a first-moment (Siegel k=1) analysis of a single path is impossible in
   principle. Missing lemma [conjecture]: **"for V <= c n^2 log n and K the union of two ellipsoids with ||A-Id||_op <= 1/2 and relative volume V each,
   P_Haar(L cap K = {0}) <= exp(C n^2) * P(L cap E_A = {0}) * P(L cap E_A' = {0}) / (overlap correction as for Poisson)"**, i.e. a second-moment bound for
   the number of L-free ellipsoids. But the first moment of that count is nu(family) * P_V, which is positive iff P_V > 0: the second-moment method is
   circular here. Conclusion: none of (i)-(v) nor a union/second-moment argument can change the order n^2 for generic lattices with presently available
   inputs; the gain must come from structure.
B. (Structured lattices; the only mechanism with a proof today.) V ~ (n^2/2) |G|/dim_Q B. In the GGS gluing framework V >> d^2 M/C, C = sum_{n in N} phi(n), M = lcm(N).
   Exact missing lemma: **existence of admissible gluing data (support exponents A(S) <= -delta for all proper S) with M/C >= c log d** (currently log log d).
   Unverified back-of-envelope [heuristic, NOT checked against the paper's definitions]: averaging the admissibility condition over S and its complement gives
   the isoperimetric necessary condition log|Res(F_S, F_{S^c})| >= (C/2) H(c_S/C) + 2 delta; data with several maximal elements (e.g. Pg and Qg, lcm PQg, M not in N)
   are not obviously excluded and would give M/C ~ (g/phi(g)) PQ/(P+Q), which could exceed M/phi(M) by a factor up to ~ (g/phi(g)) log(g/phi(g)) -- i.e. possibly
   (log log d)^{2-o(1)} instead of log log d. This is a finite combinatorial optimisation worth handing to another agent.
C. (Constant.) Prove V_n >= (1/2 - o(1)) n^2 via the variance-reduction lemma of 1.4 (needs a Rogers two-point input); (1 - o(1)) n^2 would follow for a
   drift-dominated sticky process by 2(ii). Both are constant-factor improvements over Rezvyakova's 1/(2e).
