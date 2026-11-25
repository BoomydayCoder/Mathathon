# R2: random perfect / extreme forms by sticky walk + greedy LP descent on the Ryshkov polyhedron (Z^n, n = 2..10)

All numbers below were produced by the code in this directory (see results_tables.md, structure_table.md, data/*.json).
Status of every statement: NUMERICAL OBSERVATION (floating point, tolerances 1e-8/1e-9), not proof.

## Method (voronoi.py, fp_enum.c)
* Forms A on Z^n normalized to min A[x] = 1; V = kappa_n det(A)^(-1/2) (= 2^n * density).
* Short vectors: own LLL (Gram version) + Fincke-Pohst enumeration in C (ctypes).
* Start: Wishart form G^T G, LLL-reduced, scaled to min 1.
* Sticky walk (discrete analogue of Klartag's process): B = Gaussian direction in F = {B : B[x]=0 on current
  minimal vectors} (orthonormal basis for the trace inner product); move A + tB to the first t where a new
  minimal vector appears (exact line search: t* = min (A[x]-1)/(-B[x]), iterated enumeration); repeat until rank
  {x x^T} = n(n+1)/2 (perfect form).  Two sign rules: "descent" (sign with Tr(A^-1 B) < 0, det decreases
  monotonically) and "random" (random sign, flipped only if B is PSD = infinite ray).
* Greedy descent: LP min Tr(A^-1 B) s.t. B[x] >= 0 (x in S), sum_x B[x] = 1 (HiGHS dual simplex).  The LP
  optimum is an extreme ray (checked: active contact set has rank N-1 in 100% of the steps; "nonedge" = 0), and
  by LP duality its value equals the max-min eutaxy coefficient; value > 0 <=> form is eutactic.  Move along the
  edge to the neighbouring perfect form; stop when the LP value is > 0.  Every final form was verified perfect
  and eutactic (independent eutaxy LP, min coefficient > 0), i.e. it is an extreme form (strict local optimum).
* After each move the form is "polished" (solve A[x]=1 on the integer contact set) and LLL-rebased.

## Validation (all passed)
* Enumeration: D4 12 pairs, D5 20, E8 120; V(D4)=9.8696, V(E8)=64.9394.
* n=2: 300/300 runs end at hexagonal (V=3.6276); n=3: 300/300 at FCC (5.9238);
  n=4: only D4 (9.8696, 12 pairs) and A4 (8.8276, 10 pairs) occur; n=5: exactly the 3 known perfect forms
  (V=14.888 D5, 13.234 A5^{+3}-type, 12.156 A5) ; n=6: 5 perfect classes seen (7 exist), 4 extreme (6 exist);
  n=7: 8 (descent) / 14 (random) distinct extreme forms seen (30 exist); E6, E7, E8, Lambda9 (det(4A)=512, 136 pairs),
  Lambda10 (det(4A)=768, 168 pairs), K10 (det 972, 138 pairs) all appear with the correct determinant and kissing number.
* No failed runs (0 exceptions) in 3403 samples (2030 descent rule, 1373 random rule).

## Main observations
1. Best known lattice is found as a local optimum with large frequency for n <= 8:
   descent rule: n=4 76%, n=5 96%, n=6 90%, n=7 89%, n=8 60% (already at the FIRST perfect form, no LP descent
   needed); n=9: Lambda9 1.5% at perfect stage, 23% after greedy descent; n=10: Lambda10 0% / 7%.
   random rule: n=8 29% -> 37%, n=9 0% -> 29%, n=10 4.5% -> 14%.
2. The sticky walk almost always ends directly on an extreme form for n <= 8 (mean number of LP descent steps
   0.00-0.12); for n = 9, 10 it needs 0.5-1.5 further edge steps on average.
3. Typical values.  Median V of first perfect form, random rule: n=8 37.8 (0.59 n^2), n=9 46.0 (0.57 n^2),
   n=10 53.0 (0.53 n^2); descent rule: n=9 49.0 (0.61 n^2), n=10 58.8 (0.59 n^2).  Median extreme V:
   n=9 51-52 (0.63 n^2), n=10 65-70 (0.65-0.70 n^2).  For n <= 7 the median is the best lattice itself.
   So in this range a "typical local optimum" has V ~ (0.5 .. 1.0) n^2, i.e. the same order as Klartag's c n^2,
   and it is below the best known lattice by a factor 1.4-1.8 for n = 8..10.
4. Growth fits of the median (n=3..10 or 5..10): power law exponent 1.8-2.2 with rms 0.11-0.28 (log2 units),
   exponential 2^(0.34..0.52 n) with rms 0.16-0.35.  The power law fits better in every case, but the range
   n <= 10 CANNOT discriminate polynomial from exponential growth (the data n<=8 are dominated by root lattices).
   Local slope of log2(median V_perfect), random rule: n=8->9: 0.28, n=9->10: 0.20 per dimension; n^2 would
   give 0.34, 0.30; Andreanov-Scardicchio's 2^(0.16 n) would give 0.16.  The slope is decreasing with n; the
   data are compatible with A-S (slope 0.16 asymptotically) and with a power law; samples at n=9,10 are small
   (48-66 and 22-28 runs) so the local slopes have errors of order 0.05-0.1.
5. Contacts: extreme forms found at n = 9, 10 have on average about 2 N pairs of minimal vectors (N = n(n+1)/2; minimum seen 1.16 N at n=9, 1.22 N at n=10); the
   best lattices have 3.0-3.3 N (E8: 120 = 3.33 N, Lambda9: 136 = 3.02 N, Lambda10: 168 = 3.05 N).
   V is strongly increasing in the number of contact pairs across the extreme forms found (see structure_table.md).
6. Trajectory (results_tables.md, last table): median V when the rank of the contact set reaches the fraction f
   of N grows concavely: e.g. n=10, descent: f=0.25: 22.5, 0.5: 36.9, 0.75: 46.6, 0.9: 51.4, 1.0: 58.8.  For
   n <= 8 the last 10% of the contacts produce a jump (landing on a root lattice: many contacts at once).
7. Rationality / structure of all extreme forms found: every final form is rational with a small denominator:
   q*A integral with q in {2,4,6,8,10,12,14} (q = minimum of the primitive integral form; 196 (mode,n,form) classes examined, 95 of them have q = 4).  q = 2 exactly for the
   root lattices (A_n, D_n, E_n); the best n=9,10 forms and most others have q = 4, i.e. they are integral
   lattices of minimum 4 whose minimal vectors have mutual inner products in {0, +-1, +-2} (and 4);
   det(4A): Lambda9 512, Lambda10 768, K10 972, others 576..1140.  In an LLL basis all minimal vectors have
   coordinates of absolute value <= 3 (<= 2 for the forms listed in structure_table.md with V near the top) and they generate Z^n (index 1) in all 196 cases examined (data/structure.json).
   => Suggestion for explicit constructions: search among integral lattices of minimum 4 (norm-4 "root-like"
   systems with inner products 0,+-1,+-2), e.g. laminated/Kappa-type chains, rather than generic real forms.

## Caveats
* The sampler is not uniform over perfect forms (it is the hitting distribution of a ray walk); "typical" is
  relative to this sampler.  Andreanov-Scardicchio use a random walk on the Voronoi graph; only slopes are comparable.
* n = 9, 10 have small samples; n >= 11 was not attempted (time budget ~18 min total).
* Distinct forms are separated by the invariant (V rounded to 1e-5, number of minimal pairs), which may merge
  non-isometric forms with equal det and kissing number (e.g. n=7: V=26.7273 occurs with 36 and with 42 pairs).
