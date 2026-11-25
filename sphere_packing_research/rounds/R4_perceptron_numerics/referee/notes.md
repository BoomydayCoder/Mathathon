# Referee report on R4_perceptron_numerics (adversarial check)

Wall clock 07:42 - ~07:58 UTC.  Code: `ref_check.py` (independent re-implementation of both LPs, with exact dual
certificates), `check_consts.py` (quadrature / Monte Carlo of the constants, law of the Gordon auxiliary threshold),
`summarize.py` -> `referee_summary.txt`; raw data `refA_*.jsonl`, `refB_*.jsonl`, `K10_d*.jsonl` (the other agent's
`partK.py` re-run with other seeds and other values of the revelation parameter delta).  Every number below comes from these runs.
Tags: [proved], [model-proved] (proved inside the Gaussian/Poisson model), [sketch] (standard argument, not written in
full), [numerical], [heuristic], [conjecture], [refuted], [literature].

## 1. Constants re-derived independently
1. Offsets [proved, exact].  Put s := 2 log(|x|/R).  A pair +-x with |x| = R e^{s/2} lies on the sphere bounding the ball of
   volume v e^{ns/2}; with 1/2 pair per unit volume the number of pairs with n s <= sigma has mean (v/2) e^{sigma/2}, i.e.
   intensity (v/4) e^{sigma/2} d sigma.  The constraint x^T A x >= 1, A = R^{-2}(I+B_0), reads EXACTLY
   n log(1 + u^T B_0 u) >= -sigma (u = x/|x|); the model replaces n log(1+w) by n w (error n w^2/2 ~ |B_0|^2/n per constraint).
2. Variance [proved; correction].  For u uniform on S^{n-1} and symmetric traceless B: E (u^T B u)^2 = 2|B|_HS^2/(n(n+2))
   (from E u_i^4 = 3/(n(n+2)), E u_i^2 u_j^2 = 1/(n(n+2))), NOT 2|B|^2/n^2.  Monte Carlo (4e5 directions): n=10: 0.76733 vs
   0.76559 (the n^2-formula gives 0.91871); n=26: 0.94115 vs 0.94248 (1.01498).  So rho^2 = 2 |B_0|^2 n/(n+2), and the
   threshold rho^2 = 8 corresponds to |B_0|^2 = 4(n+2)/n = 4.80 (n=10), 4.50 (16), 4.31 (26), not 4.  This explains most of the
   "slow convergence of |B0|^2" in the other agent's table (5.10, 4.62, 4.39, 4.58, 4.29 for n = 10..26).
3. Balance equation [model-proved computation].  F(r) = int (1/4) e^{s/2} E(s+rZ)_-^2 ds = 4 e^{r^2/8}, count = (1/2) e^{r^2/8}
   (times v): quadrature gives 4.532594/4.532594 (r^2=1), 10.873127/10.873127 (8), 29.556224/29.556224 (16).  Maximiser of
   x e^{-x/8}/4: x = 8.00000, value 0.735759 = 2/e.  The arithmetic of (C1), (C2), (D1) is right.
4. Volume factor [proved exact / second order].  V/v = det(I+B_0)^{-1/2} = exp(|B_0|^2/4 - Tr B_0^3/6 + Tr B_0^4/8 - ...).
   EXACT first moment in the geometric Poisson model: E #pairs in E_A = (v/2) E_u (u^T(I+B_0)u)^{-n/2} = vol(E_A)/2 (polar
   coordinates), so P(E_A empty) = exp(-vol(E_A)/2) exactly; "V = v e^{rho^2/8}" is only its second-order form.  Size of the
   error: for a GOE-like B_0 with |B_0| = 2, n = 26: det(I+B_0)^{-1/2} = 3.2823 (MC of the polar formula 3.2519) versus
   e^1 = 2.7183, i.e. log(V/v) = 1.19 instead of 1 (the term Tr B_0^4/8 ~ |B_0|^4/(4n)).  At the LP optimum the spectrum of
   B_0 is skewed: Tr B_0^3 = 5.9 +- 1.7 (n=10), 2.3 +- 0.5 (n=16), 2.4 +- 0.5 (n=20).  Corrections are O(1/n) relative.
5. Trivial but useful [proved]: det A <= (Tr A/n)^n, hence V >= v always.  So a proof of "v_LP >= (2/e - eps) N" for lattices
   alone would already give V_n >= (1/e - o(1)) n^2, twice Rezvyakova's constant, without any control of the eccentricity.

## 2. Is Gordon's theorem used correctly?  (C1)/(C2)
Auxiliary problem.  Phi(G) = min_{|theta|<=rho} max_{lam>=0,|lam|<=1} lam^T(-sigma - G theta) = min |(-sigma-G theta)_+|_2.
With w = -theta the bilinear form is lam^T G w and Gordon's auxiliary optimisation is
  phi = min_{r<=rho} max_{beta<=1} beta ( |(r g - sigma)_+| - r |h| ) = min_{r<=rho} ( |(r g-sigma)_+| - r|h| )_+ ,
g in R^M, h in R^N standard Gaussian.  So AO-feasible iff exists r <= rho: X_r := sum_i (sigma_i - r g_i)_-^2 <= r^2 |h|^2;
E X_r = 4 v e^{r^2/8}, Var X_r = (v/4) e^{r^2/8} * int u^4 e^{-u/2} du = 192 v e^{r^2/8}, |h|^2 = N(1+O(N^{-1/2})).  [model-proved]
Points a referee must check, and my verdict:
 (a) convexity: {|theta|<=rho} and {lam>=0,|lam|<=1} are convex compact, objective bilinear + linear, sigma independent of G:
     both directions of the convex Gaussian min-max theorem apply (conditionally on sigma).  OK.
 (b) "Phi = 0" is not an event of the form {Phi > c}: use slack.  If Phi_v <= c then at volume v' = v e^{-c/2} every
     violation (which is <= c) is absorbed by the shift of the offsets, so K_{v'} cap B_rho is non-empty.  P(Phi_v > c) <=
     2 P(phi_v >= c) -> 0 when 4 v e^{r^2/8} < (1-eps) r^2 N for some r <= rho.  OK, costs a factor e^{-c/2} -> 1.
 (c) infinitely many constraints: the truncations K^{(S)} (offsets <= S) decrease to K; by compactness of B_rho,
     {K cap B_rho non-empty} = intersection over S of {K^{(S)} cap B_rho non-empty}; the AO values increase to the M = infinity
     value, which is finite and concentrated (variance above).  For infeasibility a single large S suffices.  OK.
 (d) unbounded theta (the LP has no norm constraint): for r >= R_0, X_r >= r^2 sum_{sigma_i<=S}(g_i - S/R_0)_+^2 ~ r^2 (v/2)
     e^{S/2} E(Z-S/R_0)_+^2 > 2 r^2 N for S = 2 log(16N/v), R_0 large; on [0,R_0] uniform LLN.  So P(K_v non-empty) =
     lim_rho P(K_v cap B_rho non-empty) <= 2 P(exists r >= 0: X_r <= r^2|h|^2) -> 0 when v >= (2/e+eps)N.  OK.
 Conclusion [sketch; I consider it a routine theorem in the abstract model, but no complete proof is written anywhere in
 this project]: for fixed eps and N -> infinity, K_v is empty w.h.p. if v >= (2/e+eps)N; if v <= (2/e-eps)N then K_v is
 non-empty and min{|theta|: theta in K_v} -> rho_-(v/N), the smaller root of x e^{-x/8}/4 = v/N.  NOT covered by Gordon:
 the number of active constraints on the min-norm path (c = rho^2 N/8), the law (D2) of the active offsets, and
 everything about the geometric model (the g_i = n(uu^T - I/n)/sqrt2 are rank-one, not Gaussian; universality is only
 supported by the numerics, and must fail for low-rank theta, where <theta,g> is chi-square-like).
 Finite-N law of the AO threshold (check_consts.py, cutoff c = 300): v_AO/N = 0.770 +- 0.011 (sd 0.211, N=100, 400
 samples), 0.752 +- 0.010 (sd 0.145, N=200, 200 samples), 0.774 +- 0.011 (sd 0.105, N=400, 100 samples).  So the AO itself
 is biased upwards by 2-5% at these N (max over r of a noisy function), and the per-sample fluctuation of v*/N is 15-20% at
 N = 100: the LP data test 2/e only at the 2-3% level.

## 3. Independent numerical re-runs (my own code, no shared functions)  [numerical]
Abstract LP (offsets generated from arrival times of a Poisson process; optimum CERTIFIED by an exact dual solution
lambda > 0 on the N+1 active rows, |dual gap| <= 1.1e-12; cutoff tested by checking the c=300 optimum against all offsets
up to c = 1500 (N=100) / 600 (N=200)):
| N | samples | v*/N (0.7358) | rho^2 (8) | v* e^{rho^2/8}/N (2) | active = N+1, certified | samples changed by offsets beyond c=300 | mean / sd active offsets (2 / 3.46) |
|---|---|---|---|---|---|---|---|
| 100 | 40 | 0.7725 +- 0.025 | 7.96 +- 0.25 | 2.054 +- 0.029 | 40/40 | 0 (min slack beyond cutoff 0.60) | 2.00 +- 0.08 / 3.41 +- 0.08 |
| 200 | 10 | 0.7695 +- 0.038 | 7.84 +- 0.32 | 2.032 +- 0.051 | 10/10 | 0 (min slack 0.52) | 2.05 +- 0.14 / 3.39 +- 0.11 |
The other agent: 0.7332 +- 0.023 (N=100), 0.7302 +- 0.024 (200).  Both sets are within 1.5 s.e. of 2/e; mine are above,
theirs below.  Verdict: confirmed at the 3% level; cutoff c = 300 is adequate (0 of 50 samples affected).
Geometric trace-LP (pair volumes from arrival times, V_max = 3000N (2000N for n=20); optimum certified by a dual
solution I = sum lambda_j x_j x_j^T with all lambda_j > 0 on exactly N active pairs):
| n | samples | PD, N active, certified | V_LP/N (2) | v_LP/N (0.7358) | \|B0\|^2 (4; 4(n+2)/n) | log(V/v) (1) | lam_max/lam_min | max E[#unsampled pairs inside] |
|---|---|---|---|---|---|---|---|---|
| 10 | 40 | 40/40 | 1.966 +- 0.042 | 0.619 +- 0.034 | 5.32 +- 0.48 (4.80) | 1.227 +- 0.074 | 11.0 +- 0.8 | 0.0095 |
| 16 | 10 | 10/10 | 2.003 +- 0.038 | 0.632 +- 0.035 | 4.70 +- 0.33 (4.50) | 1.165 +- 0.059 | 6.7 +- 0.5 | 0 |
| 20 | 4 | 4/4 | 2.036 +- 0.042 | 0.646 +- 0.046 | 4.87 +- 0.42 (4.40) | 1.155 +- 0.087 | 5.9 +- 0.4 | 0 |
Agreement with the other agent's table (1.976, 0.625, 5.10 at n=10).  Verdict: "V_LP ~ 2N in the Poisson model for n = 10..26"
is reproduced.  NOT established: the limits v_LP/N -> 2/e and log(V/v) -> 1; v_LP/N is 0.62-0.69 for all n <= 26 and
log(V/v) is 1.07-1.23; the value V_LP/N = 2 at n = 10 comes from two compensating finite-n errors of ~15% each.
Klartag process (their partK.py, new seeds): the output law must not depend on the revelation margin delta if the lazy
revelation is exact.  n = 10: delta = 0.1: V_K/N = 1.121 +- 0.016 (150 runs, 648 snapshots/run); delta = 2.0: 1.104 +- 0.021
(77 frozen runs, 12 snapshots/run; 1.113 +- 0.020 over all 83 runs if the 6 runs that hit my 20 s cap at rank 53-54 of 55
are included with their current volume; those have V/N = 1.01-1.53, mean 1.22, so dropping capped runs biases slightly
downwards; I stopped this batch at 83 of 150 runs for time); their delta = 0.4: 1.135 +- 0.023.  Consistent.
(The planned delta = 0.4 and n = 14 batches were not run: time.)  I also read
the code: the Loewner test at the end point of each straight sub-step is sufficient (A >= S is convex), new points are
filtered against all earlier snapshots and the initial ball: the scheme is exact up to the Euler step.  [numerical + code review]

## 4. Objections that survive (corrections to the report / to the orchestrator's picture)
(O1) [proved, lattice fact] For a lattice L = gZ^n the trace-LP is min <Q, G^{-1}> over the Ryshkov polyhedron
  {Q: Q[z] >= 1 for z in Z^n\0} (Q = g^T A g, G = g^T g).  Its optimum is a vertex, i.e. a PERFECT form Q*, and LP duality
  (the certificate above: G^{-1} = sum lambda_j z_j z_j^T, lambda >= 0, z_j minimal vectors) says exactly that the dual Gram
  matrix G^{-1} lies in the Voronoi domain of Q*.  Hence V_LP(L) = kappa_n det(Q*)^{-1/2} is the Hermite-type invariant of the
  perfect form whose Voronoi domain contains G^{-1}; the lattice version of (C1) is the statement "the perfect form whose
  Voronoi domain contains a Haar-random Gram matrix has V = (1+o(1)) n^2".  [Voronoi theory: literature]
(O2) [refuted for small n] The Poisson-model value V_LP ~ 2N "already at n = 10" cannot hold for lattices in low dimension:
  V_LP(L) <= V_n and the densest lattices have V_n/N = 1.21, 0.99, 0.99, 0.99, 1.14, 1.35, 1.80 for n = 2..8 (optimal,
  proved) and 1.66, 1.71, 1.88 for the best known lattices in n = 9, 10, 11 (2.60 at n = 12, 7.09 at n = 16, 107.9 at
  n = 24).  So for n <= 8 (and for n <= 11 unless denser lattices exist) NO lattice reaches 2N.  The Poisson numerics at
  n = 10..26 are therefore not predictions for lattices of the same dimension; Part C would have shown this.
(O3) [proved] For a fixed lattice the supremum of V over ALL empty ellipsoids equals V_n (every ellipsoid is a ball for
  another lattice), a finite constant independent of L, while for the Poisson process it is infinite.  Hence the Poisson
  model necessarily breaks down at some eccentricity, (C4)-type extrapolation to large rho has no support from this round, and
  the Klartag end points in the Poisson model (lam_max/lam_min ~ 500-2000) are in the regime where the model is least trusted.
(O4) (C4) upper bound.  The one-sided Gordon inequality on the sphere |theta| = rho is legitimate without convexity and
  gives V <= (1+o(1)) rho^2 N/4 [sketch, model].  The "union bound over a net with spacing O(1)" is NOT a proof: a net
  that captures every feasible theta needs spacing ~ 1/V (expected symmetric difference of the violated sets O(1)), or a
  margin argument (relax all offsets by a constant, random net point) that needs a uniform bound on the number of nearly
  active constraints; the naive version only gives V <= N log(C N rho^2).  Status: [heuristic / literature, unverified here].
(O5) Minor inconsistencies in the report: the reading guide says v_LP/N "0.62 -> 0.71" and |B0|^2 "5.1 -> 4.2", the table says
  0.6925 and 4.291; the pivot criterion should read "iff lam_j <= 0" (lam_j = 0 also increases, by strict concavity).
  "V_K -> N" is an extrapolation: the data 1.135, 1.119, 1.058, 0.998, 1.023 (n = 10..26) decrease and do not exclude a
  limit below 1.  The (C3) "reflection principle" derivation was not tested by anything in this round.

## 5. What I would state in a paper
(a) Theorems about the model: the two expectations (D1); the exact first moment P(E_A empty) = e^{-vol/2}; V >= v; the pivot
    criterion; (after writing out section 2) the threshold v* = (2/e)N and the min-norm radius rho_-(v/N) in the abstract
    Gaussian model; the one-sided bound V <= rho^2 N/4 on spheres; (O1) and (O3) for lattices.
(b) Numerical observations: N+1 active constraints and the law (D2); V_LP/N = 1.95-2.05 in the geometric Poisson model for
    n = 10..26; V_K/N = 1.0-1.14 with exactly N contacts for the Poisson Klartag process; min-norm path for v/N <= 0.5;
    no extreme form reached by the pivoting ascent.
(c) Conjectures for lattices only: V_LP(L) = (1+o(1)) n^2 for Haar-typical L as n -> infinity (false for n <= 8),
    v_LP(L) = (2/e+o(1)) N; Klartag's process ends at (1+o(1)) N; everything in (C4), (C5).

## 6. Proofs of the two lattice/Poisson facts used in (O1), (O3)
(O1) [proved, using Ryshkov/Voronoi: literature].  x = gz, x^T A x = Q[z] with Q = g^T A g, Tr A = Tr(Q g^{-1} g^{-T}) =
<Q, G^{-1}>, det A = det Q (covolume 1).  Every Q in the Ryshkov polyhedron R is positive definite and <Q,P> >=
lambda_min(P) Tr Q for P = G^{-1} > 0, so sublevel sets of the objective on R are compact; R is locally finite, so the
minimum is attained at a vertex Q*, i.e. a perfect form with minimum 1.  Locally the problem is a finite LP, and LP duality
gives: Q* optimal iff P is a non-negative combination of the active normals z z^T (Q*[z] = 1), i.e. P lies in the Voronoi
domain of Q*.  g -> g^{-T} preserves Haar measure, so G^{-1} is as random as G.
(O3) [proved].  (i) Lattice: E_A = {Ax.x<1} is free of non-zero points of L iff the unit ball is free of non-zero points of
A^{1/2}L (covolume det(A)^{1/2}), and every lattice is a rotation of some pL, p > 0 symmetric (polar decomposition); so
sup_A vol(E_A) = sup over all covolume-1 lattices of the largest empty ball = V_n.  (ii) Poisson process, n >= 2: fix V_0.
Let E_k be the ellipsoid with one semi-axis L_k along e_1 and the other semi-axes equal, of volume V_0, and R_k < L_k radii
with R_{k+1} >= L_k and R_k/L_k <= k^{-2}.  vol(E_k cap B(0,R_k)) <= vol(E_k cap {|x_1| <= R_k}) <= (3/2) V_0 R_k/L_k is
summable, so a.s. only finitely many E_k have a point in B(0,R_k) (Borel-Cantelli); the sets E_k \ B(0,R_k) lie in the
disjoint shells B(0,L_k) \ B(0,R_k), so the events {no point in E_k \ B(0,R_k)} are independent with probability >=
e^{-V_0/2}, and infinitely many occur a.s.  Hence a.s. some E_k is empty, for every V_0: sup V = infinity.
(For the slab bound: the section of E_k at height x_1 has (n-1)-volume <= that of the central section, and
V_0 = central section * L_k * 2 * int_0^1 (1-t^2)^{(n-1)/2} dt >= central section * 2 L_k * (2/3) for n-1 <= 2; in general the
constant 3/2 is replaced by 1/int_0^1 (1-t^2)^{(n-1)/2} dt, a constant depending on n only, which does not affect summability.)
