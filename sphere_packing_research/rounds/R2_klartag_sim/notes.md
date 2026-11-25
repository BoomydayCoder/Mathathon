# R2: simulation of Klartag's stochastically evolving ellipsoid, n = 2..10

Files: `sim.py` (simulator), `analyze.py` (statistics), `data/*.jsonl` (one JSON record per run, incl. the full list of
contact events (t, V, #pairs, rank)), `results_table.md` / `summary.json` (tables), `klartag_sim.png`, `dt_sensitivity.txt`,
`data_dt/` (time-step test), `data_parts/`, `data_aborted_firstlaunch/` (discarded first launch, wrong time normalisation).

## Method (what was actually run)
- Lattice: Goldstein-Mayer family {x in Z^n : x_n = a.x' mod p}, p = 100003, a uniform, LLL-reduced (own float LLL),
  scaled to covolume kappa_n (`gm`); and L = scaled Z^n (`zn`).  lambda_1 found by own Fincke-Pohst enumeration.
- By the scaling symmetry (A -> cA <=> t -> c^2 t) the lattice is rescaled so that lambda_1 = 1 and A_0 = 1.001 Id
  (ball touches nothing).  All times are in this normalisation (a_0 ~ 1, Klartag's), V = kappa_n det(A)^(-1/2)/covol.
- dA = pi_F(dW): Euler steps dt = 1.25e-3/n; W standard Gaussian in the orthonormal coordinates of (Sym_n, Tr(AB));
  pi_F = projection onto the orthogonal complement of span{x x^T : x contact} (Gram-Schmidt basis, re-orthogonalised).
  Along each straight segment the first hitting time of every candidate lattice point is computed exactly; all points with
  A x.x <= 1+1e-9 at the event become contacts forever (ties = non-simple faces are handled); the rest of the Gaussian
  increment is re-projected and continued.  Candidates: all lattice points with A0 x.x <= 2, re-enumerated (LLL + Fincke-Pohst)
  before the relative Loewner change of A w.r.t. A0 can exceed 1/2 (this makes the candidate list provably complete).
  Contacts are re-projected to A x.x = 1 (min-norm correction) at every event/re-enumeration.
- Stop when rank span{x x^T} = n(n+1)/2 =: d.  The final form is checked by an independent LLL + enumeration.

## Validation
- n=2: 207/207 runs end at the hexagonal form V = 3.6276; n=3: 225/225 (random L) and 83/83 (Z^3) end at FCC V = 5.9238.
- n=4: only V = 9.8696 (D4, 12 pairs) and 8.8276 (A4, 10 pairs) occur = the two perfect forms.  n=5: exactly three values
  14.888 (D5, 20 pairs), 13.234 (A5^{+3}), 12.156 (A5) = the three perfect forms.  n=6: 23.869 (E6, 36 pairs), 20.671 (D6),
  17.653, 18.374, 17.858 and (Z^6 runs) 15.626 (A6): 6 of the 7 perfect forms.  n=8: 64.939 = E8 with 120 pairs.
- All 1878 recorded final forms are L-free in the independent check (min A x.x >= 1 - 5.2e-9).
- Time-step test (n=4, fraction ending in D4): dt: 0.468+-0.034 (220 runs); dt/4: 0.514+-0.084 (35 runs); 4dt: 0.463+-0.056 (80 runs).
- Problems: 1 run of 164 in n=8 aborted (Cholesky of a non-PD matrix, near a highly degenerate vertex) - excluded and counted.
  Internal consistency flags at re-enumeration (contact set rebuilt): n=4: 1 of 220 runs carry flags (flag counts [138]; their final V: [8.8276]; their final min A x.x >= 0.9999999948); n=8: 1 of 163 runs carry flags (flag counts [7]; their final V: [64.939]).  The final independent check passed for these runs.

## Observations (numerical; see results_table.md)
1. End point = a random perfect form, strongly biased to the dense, highly non-simple vertices.  Fraction of runs ending at
   the densest lattice: n=2: 1; n=3: 1; n=4: 0.47 (D4); n=5: 0.845 (D5); n=6: 0.67 (E6); n=7: 0.79 (E7); n=8: 0.33 (E8, 54/163);
   n=9: 0/100 (max 60.94 vs Lambda_9 74.64); n=10: 0/50 (max 75.38 vs Lambda_10 94.23).  Z^n start gives the same frequencies
   within errors (0.49, 0.85, 0.72, 0.34 for n=4,5,6,8).
2. "n(n+1)/2 contact pairs generically" fails in small n: the mean final number of contact pairs is 10.9, 19.2, 33.6, 56.8,
   71.4, 60.6, 76.1 for n=4..10 (d = 10, 15, 21, 28, 36, 45, 55); fraction with exactly d pairs: 0.53, 0.16, 0.04, 0.01, 0.01, 0, 0.
   The surplus arrives at the very last event (end of the last edge).
3. Volume per contact: E[V_final]/n^2 = 0.907, 0.658, 0.582, 0.579, 0.630, 0.718, 0.700, 0.559, 0.568 (n=2..10), i.e.
   E[V_final]/d = 1.21, 0.99, 0.93, 0.97, 1.08, 1.26, 1.24, 1.01, 1.03.  Along the path V(t)/k (k = rank of the contact set) is a
   nearly n-independent curve: about 1.75 at k=1, 1.2 at k=5, 1.13-1.19 at k=10, 1.0-1.13 at k=20, 0.96-1.02 for k>=36.
   So "swept volume ~ number of independent contact pairs" holds with constant ~1 (V ~ k), and the process delivers
   V_final ~ n(n+1)/2, i.e. c ~ 1/2 in V ~ c n^2, plus a bonus from the last edge when it ends in E6/E7/E8.
4. Growth: E[log V_final] = 1.289, 1.779, 2.230, 2.671, 3.118, 3.549, 3.755, 3.808, 4.033 (n=2..10).  Fits over n=2..10:
   a + b log n: a=-0.126, b=1.81, rms 0.110;  a + b n: a=0.815, b=0.350, rms 0.187.  Up to n=7 the growth looks exponential
   (it tracks the best lattices: mean/best = 1, 1, 0.944, 0.973, 0.950, 0.931), then it collapses to the n^2 law:
   mean/best = 0.690, 0.607, 0.603 for n = 8, 9, 10.
5. Time scale: median freezing time / T_K (T_K = 16 log n/n^2) = 0.25, 0.32, 0.45, 0.66, 1.09, 1.40, 1.33, 1.13, 1.14 for n=2..10;
   mean V(T_K)/mean V_final = 0.91, 0.82, 0.82, 0.95, 0.94 for n=6..10.  Most contacts arrive in a burst (n=10: rank 5 at
   t=0.05, rank 45 at t=0.14), the last few faces take most of the time (rank 54 at 0.28, rank 55 at 0.50).
6. 2 zeta(n) (Minkowski-Hlawka) is beaten by a factor growing like n^2/4: 4.3 (n=4), 11.1 (n=6), 22.3 (n=8), 28.4 (n=10).

## Caveats
Euler discretisation (bias not detected in the n=4 test, not tested for n>=6); n=8,9,10 have 163, 100, 50 runs only;
the Goldstein-Mayer lattices are only approximately Haar; identification of final forms is by (V, number of minimal pairs) only.
