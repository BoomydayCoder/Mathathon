# R9: trace-LP on the Ryshkov polyhedron of GENUINE random lattices

Directory: `rounds/R9_trace_lp_lattices/`. Code: `tracelp.py` (LP by constraint generation, sticky walk with Q-isotropic
directions, greedy descent re-used from `R2_voronoi_descent/voronoi.py`), `lp_scan.py` (LP only, independent seeds, more n),
`validate.py`, `analyze.py`, `make_notes.py`. Raw data: `data/n{n}_{gm,gauss}.json` (full pipeline, with Q and every form),
`data/lpscan_n{n}_gm.json`, `data/validation.json`, `data/poisson_gamma_check.txt`. Full tables: `results_table.md`,
`summary.json`, `summary_lp_pooled.json`. Logs: `log1..5.txt`.

Normalisation: lattice Z^n with Gram matrix Q, det Q = 1; A is lattice-free iff min_{x != 0} A[x] >= 1;
V = kappa_n det(A)^(-1/2); N = n(n+1)/2; v = kappa_n (Tr(Q^-1 A)/n)^(-n/2).
Families: `gm` = Goldstein-Mayer lattice {x : x_n = a.x' mod p}, p = 100003 (close to Haar-random); `gauss` = Gram of an
i.i.d. N(0,1) basis (NOT Haar: more skewed lattices). Both LLL-reduced, det 1.

## 0. Checks (all numbers below come from the code in this directory)

* [numerical, validation] `validate.py`: for Q = A6, D6, E6, E7, E8, D4, A10 (extreme forms, so Q^-1 is interior to the
  Voronoi domain of Q) the trace-LP returns Q/min(Q) with the right V and number of minimal pairs: 7/7 OK.
* [numerical] Every reported form (LP vertex, walk end point, greedy end point) was re-verified by a final complete
  Fincke-Pohst enumeration: max |min A - 1| over all forms = 1.1e-11 (LP vertices: <= 2e-14); no V exceeds the best known
  V_n (max V/Vbest = 1.000000, attained at n = 6, 7, 8, 9, 10, 12 by E6, E7, E8, Lambda9, Lambda10, K12).
  0 LP failures, 0 walk/descent errors in all runs.
* The LP is solved exactly on a finite subset of constraints (8N shortest vectors of Q, then all violated vectors of the
  current optimum, 1-4 rounds); the final optimum is feasible for ALL integer vectors, hence optimal.

## 1. What the trace-LP is (proved, standard LP duality)

A* minimises Tr(Q^-1 A) over the Ryshkov polyhedron  iff  Q^-1 = sum_{x in Min(A*)} lambda_x x x^T with lambda_x >= 0,
i.e. iff Q^-1 lies in the Voronoi domain of the perfect form A*. So the trace-LP IS Voronoi's reduction of the dual form:
it returns the perfect form whose Voronoi cone contains Q^-1. The Ryshkov polyhedron is the same for every lattice; the
only randomness is the direction Q^-1 of the objective. For Haar-random lattices P(vertex = P) is proportional to the
invariant measure of the Voronoi domain of P, so fat (non-simplicial) domains with many minimal vectors are favoured.

## 2. Trace-LP vertex: pooled results (main runs + LP-only scan)

@@TABLE_LP@@

(mean+-standard error over lattices; "simple" = exactly N active pairs; full version incl. B0sq, eigenvalue ratio,
median/min/max, V_0, LP rounds: `results_table.md` Tables 1 and 5.)  Other measured quantities: mean ||B_0||_HS^2 = 2.9-5.4
at every n (3.6-4.3 for 9 <= n <= 26); mean eigenvalue ratio of Q^-1/2 A Q^-1/2 falls from 10-20 (n = 6-8) to ~5 (n = 20);
V_0/N (largest empty round ball) = 0.09-0.11 (n=6), 0.05 (n=8), 0.02-0.03 (n=10), 0.01 (n=14); the longest active vector
has Q[x] = 5-11 lambda_1^2 on average for n <= 10 (Table 1).

## 3. Comparison on the same lattices

@@TABLE_CMP@@

V_K = sticky random walk to a perfect form from the largest empty round ball, exact line search, direction uniform on the
Q-unit sphere of the free subspace, random sign (Klartag-like, unbiased) / "descent" sign (sign chosen so that det
decreases). V_ext = greedy edge descent (LP min Tr(A^-1 B) over the cone of active constraints) to a form with no
descending edge. n = 12: only 6 (gm) / 1-2 (gauss) walks finished in the time box; n = 14: 0 walks, and only 3 greedy
descents from the LP vertex (340, 382, 500 steps; the 500-step one hit the iteration cap and is NOT certified extreme;
V_ext/N = 1.53, 2.39, 1.93). At n <= 12 all greedy end points are certified (descent-LP value > 0).

## 4. Answers

1. [numerical] **n = 6, 7, 8: the LP lands on the root lattice.** E6 in 43/48, E7 in 23/24, E8 in 38/48 lattices
   (gm 17/24, gauss 21/24). These vertices are extremely degenerate (36, 63, 120 pairs for N = 21, 28, 36). The
   Poisson value V_LP/N = 2 is IMPOSSIBLE here because 2N exceeds the best known V_n for every n <= 11
   (Vbest/N = 1.14, 1.35, 1.80, 1.66, 1.71, 1.88 for n = 6..11). The sticky walk reaches E8 in only 17/48.
2. [numerical] **n = 9-13: not Poisson-like.** The LP vertex is never simple at n = 8-12 (0/172 lattices; 3/24 at
   n = 13); pairs/N = 1.2-1.7 on average with a heavy tail (n = 12: forms with 231-272 pairs, V/N = 1.69-1.91, in 4/20
   gm lattices; Lambda10 itself in 3/20 gauss lattices at n = 10). V_LP/N = 1.16-1.27, far below 2, and v_LP/N = 0.42-0.48
   instead of 2/e = 0.736.
3. [numerical] **n >= 14: the vertex becomes simple, as in the Poisson model, but V_LP/N is still far below 2 and
   grows slowly with n** (gm): 1.20, 1.28, 1.35, 1.45, 1.55 for n = 14, 16, 18, 20, 22 (about +0.04 per dimension);
   v_LP/N grows likewise (0.49 -> 0.59). The few lattices at n = 24, 26, 28 (5, 2, 3 samples) give 1.54, 1.55, 1.44, i.e.
   NO further growth; these rows are too small to tell levelling-off from noise. The spread of V_LP/N (std 0.04-0.07 for
   n = 16-28) is 2-3 times SMALLER than the Poisson spread 2/sqrt(N). The Gaussian family gives the same values at
   n = 16 and 20 (see table), so the result is not an artefact of the family. [conjecture / open] Whether V_LP/N -> 2 for
   genuine lattices as n -> infinity is NOT decided by these data: the values rise from 1.2 to about 1.5 between n = 14 and
   n = 22 and are still 25-30 % below the Poisson value at n = 22-28.
   Caveat for gm at large n: with p = 100003 the unscaled lattice has lambda_1 ~ 2.2 grid units for n >= 20
   (integer norms 4-6), so the family is visibly "granular" there and need not be close to Haar; the gauss rows and the
   p = 1000003 rows (label gm_p1000003, if present in the table) are the sensitivity checks.
4. [numerical] **LP vs sticky walk on the same lattice.** V_K/N (random sign) = 1.10, 1.26, 1.02-1.05, 1.06 for
   n = 6, 8, 10, 12, in line with round 2 (0.93-1.26) and with the Poisson-model value ~1. The LP vertex is better than or
   equal to the walk end point in most cases (n = 8: 30/48 strictly better, 14/48 equal; n = 10: 33/44 strictly better) but
   the mean gain is only 13-20 % at n = 10-12 (24-36 % at n = 8, the E8 effect; 0-2 % at n = 6), not the factor 2 of the
   Poisson model.
   Greedy descent adds 10-50 % at n = 10-12 (V_ext/N = 1.30-1.33 at n = 10, 1.48 (gm) / 1.80 (gauss) at n = 12 from the LP
   vertex; K12 reached in 2/10 gauss lattices, Lambda10 in 6/44 from the LP vertex and 5/44 from the walk end point).
   From n = 14 on the LP vertex is simple and the greedy descent needs hundreds of edge steps instead of 0-30.
5. [heuristic + numerical] **Why the Poisson model gives exactly 2 and lattices do not.** In the Poisson model the LP
   optimum is an empty ellipsoid determined by exactly N independent points in a scale-covariant way, so by the Mecke
   formula (Moller-Zuyev "Gamma-type" theorem) V_LP ~ Gamma(shape N, rate 1/2): mean 2N, std 2 sqrt(N). Check on the
   round-4 Poisson data (`data/poisson_gamma_check.txt`): std of V_LP/N = 0.277, 0.201, 0.121, 0.144 for n = 10, 14, 18,
   22 against 2/sqrt(N) = 0.270, 0.195, 0.153, 0.126. For a lattice the N contact vectors are integer combinations of n basis
   vectors (N >> n - 1, outside the range of Rogers' Poisson-like moment formulas): contact vectors x, y force x +- y to be
   lattice points with A[x+y] + A[x-y] = 4, i.e. |A(x,y)| <= 1/2 (a kissing configuration), and there is no independent
   scaling of N points. At small n this rigidity produces root-lattice/degenerate vertices; it relaxes only slowly with n.
6. [numerical] No new packing record: every form found has V <= best known (as it must for n <= 8; for n >= 9 all forms
   found are at or below the known lattices). The LP vertices at n = 16-22 have V = 1.28-1.55 N = 0.68-0.81 n^2, i.e. the
   Klartag scale c n^2 with a somewhat larger constant than the sticky walk (V_K ~ 1.0-1.1 N), far below Vbest/N = 7-35.

## 5. Caveats / gaps

* Sample sizes: 12-48 lattices per (n, family) for n <= 20, but only 9 (n = 22) and 5, 2, 3 for n = 24, 26, 28.
  Gaussian family only for n = 6, 8, 10, 12, 16, 20. The gm family with p = 100003 is granular for n >= 20 (see item 3);
  a run with p = 90000049 failed in the float LLL (Cholesky) and was discarded.
* Sticky-walk statistics at n = 12 rest on 6 + 1 walks, none at n = 14 (time box).
* The "Gamma law" in item 5 is a heuristic derivation checked numerically, not a written proof.
* The discrete sticky walk (move to the next contact along a random line) is not Klartag's continuous process; round 2
  compares the two.
* The column Vbest/N for n <= 14 uses the best known lattices of the task statement (E6, E7, E8, Lambda9, Lambda10, K11,
  K12, K13, Lambda14); for n >= 15 it uses laminated-lattice centre densities (Conway-Sloane table, quoted from memory:
  Lambda15..Lambda24, Lambda26, Lambda28) only as a sanity ceiling - for n = 26, 28 these are not claimed to be the best
  known lattices. All V found for n >= 13 are far below these values, so nothing depends on them.
* Sensitivity of the gm family to the modulus: p = 1000003 gives V_LP/N = 1.287+-0.012 (n = 16, 14 lattices) and
  1.417+-0.026 (n = 20, 7 lattices) against 1.278+-0.011 and 1.447+-0.012 for p = 100003, and 1.273+-0.012, 1.413+-0.023
  for the Gaussian family: all three agree within about 1.5 standard errors.
