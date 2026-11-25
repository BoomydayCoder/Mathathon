# R8: equivariant sticky random walk (discrete Klartag process) on G-invariant forms of O^k

Status of everything below: NUMERICAL (floating point, tolerances 1e-7..1e-9) unless tagged otherwise.
All numbers come from eq.py / summarize.py in this directory; raw runs: data/*.json (40 files, 1152 runs, 0 failed runs,
0 LP-descent errors); full tables: results_tables.md; best forms (integer Gram matrices q*A, Hermitian matrices): best_forms.json;
validation log: validate_out.txt; batch logs: logA..logE.txt.

## Method (eq.py)
* Lattice = fixed free module L = O^k in a Z-basis (coordinates Z^n). G = finite unit group of O acting diagonally
  (left multiplication), given by integer matrices. Q = {symmetric A : g^T A g = A for the generators g}, computed numerically
  as a null space (orthonormal basis for the trace inner product). Dimensions found = theory in every case:
  Eisenstein D=k^2, Gaussian k^2, Hurwitz k(2k-1), Z[zeta_m] k^2 phi(m)/2, icosians 2k(2k-1).
  Icosian order: 120 units generated from (1+i+j+k)/2 and (tau+sigma i+j)/2 embedded in H x H = R^8 by the two real places;
  Z-basis extracted from the units (all 120 units have integer coordinates; generator matrices integral: asserted).
* Start: Wishart matrix averaged over G (Reynolds operator), scaled to min 1.
* Sticky walk: B = Gaussian direction in the free subspace {B in Q : B[x]=0 on current minimal vectors}; move A+tB to the
  first t where a new minimal vector appears (exact line search with iterated Fincke-Pohst enumeration, C code fp_enum.c,
  LLL with warm start); contacts are kept; stop when the contact conditions have rank D (first Q-perfect form).
  Sign rules: "random" (flip only if the ray is infinite; closest analogue of Klartag's driftless motion) and "descent"
  (sign with Tr(A^-1 B) < 0).
* Greedy descent: LP min Tr(A^-1 B), B in Q, B[x] >= 0 on one representative per contact orbit, sum B[x] = 1; move along the
  optimal ray to the neighbouring Q-perfect form; stop when LP value >= 0 (Q-eutactic => Q-extreme = local maximum of the
  density among G-invariant forms). In all 1216 LP steps the optimal ray was an edge (active rank D-1).
* V = kappa_n / sqrt(det A) (covolume 1 in Z-basis coordinates, min A = 1).
* Every first-perfect and every final form was re-verified by an independent enumeration (cold LLL): no nonzero lattice
  vector with A[x] < 1-1e-7, and the number of vectors with A[x] <= 1+1e-7 equals the recorded contact number
  (column "verified" = True for all 40 data sets). max over all runs of V / V_best_known(n) = 1.000000 (never exceeded).

## Validation (validate_out.txt and the production data)
Eisenstein k=1 -> hexagonal 3.6276 (40/40); Gaussian k=1 -> Z^2 3.1416; Hurwitz k=1 -> D4 9.8696 (40/40);
Hurwitz k=2 -> E8 64.9394, 240 min. vectors (40/40); icosians k=1 -> E8 (20/20); Eisenstein k=3 -> E6 23.8686 (39/40);
Eisenstein k=4 -> E8 (32/40); Eisenstein k=6 and Gaussian k=6, Z[zeta_7] k=2, Z[zeta_12] k=3 -> K12 202.5643, 756 min. vectors,
det(4A)=729; Hurwitz k=4, icosians k=2, Z[zeta_8], Z[zeta_12] k=4 -> BW16 963.9143, 4320 min. vectors, det(4A)=256.
An early bug (null space computed from rounded rows -> walk stuck at t=0) was found on Hurwitz k=2 and fixed before any
production run.

## Answers to the key questions
(1) YES. Frequency with which the best known lattice is reached (first Q-perfect form / after LP descent), identified by
V (= determinant), kissing number and det(q*A):
* E8 (n=8): Hurwitz k=2 100%/100% (40 runs); icosians k=1 100%; Z[zeta_5] k=2 100%; Gaussian k=4 97%; Z[zeta_12] k=2 95%;
  Eisenstein k=4 80% (random rule), 93% (descent rule); Z[zeta_8] k=2 68%.  (Unsymmetric, R2: 29-37% random rule.)
* E6 (n=6): Eisenstein k=3 97%. Gaussian k=3: 0% (always V=20.67, 60 min. vectors = D6 by invariants).
* Lambda10 (n=10): Eisenstein k=5 7.5% / 10% (40 runs) - no better than the unsymmetric walk (4.5% / 14%); Gaussian k=5 0%.
* K12 (n=12): Z[zeta_7] k=2 (D=12) 33% / 33% (13 of 40); Z[zeta_12] k=3 5% / 5%; Eisenstein k=6 0% / 6% (1 of 18; descent
  rule 1 of 7); Gaussian k=6 0% / 6% (1 of 17); Hurwitz k=3: 0 of 80 (median = max = 170.91 = density of Lambda12, with
  624 or 648 minimal vectors, i.e. invariants of Lambda12^min / Lambda12^max); Z[zeta_5] k=3 0% (best 152.87, 560 vectors,
  det(4A)=1280); Z[zeta_8] k=3 0% (best 170.91).
* BW16 (n=16): icosians k=2 (D=12) 42% / 42% random rule (17 of 40), 67% / 67% descent rule (20 of 30); Hurwitz k=4 (D=28)
  12% / 12% random (5 of 40), 23% / 40% descent rule (30 runs); Z[zeta_12] k=4 0% / 27% (3 of 11); Z[zeta_8] k=4 0% / 16%
  (3 of 19); Eisenstein k=8 0 of 1.
* Lambda20 (n=20): Hurwitz k=5 (D=45): 0 of 10 at the first Q-perfect form, 2 of 10 after LP descent (random rule 1 of 3 after
  19 edge steps, descent rule 1 of 7); V=3382.5609, 17400 minimal vectors, det(4A)=64.
* Leech (n=24): icosians k=3 (D=30, w=120): random rule 2 of 11 (18%) already at the first Q-perfect form and 6 of 11 (55%)
  after LP descent; descent rule 4 of 13 (31%) / 5 of 13 (38%). V=32372.885, 196560 minimal vectors (98280 pairs, enumerated
  and verified), 4A integral with det(4A)=1.0000. In one run the last contact took the walk from V=2456 (rank 29) directly to
  the Leech lattice. Hurwitz k=6 (D=66): 0 of 2 (V_perf 907, 938; LP descent to 1381 and >= 1573, the second run was cut at
  200 edge steps and is NOT Q-extreme). Z[zeta_13] k=2 (D=24): 0 of 3 (V_ext 781..834).
  So the equivariant walk finds E8, K12, BW16, Lambda20 and Leech with frequencies of order 10-100% when the symmetry is
  large (w D/2 comparable to V_best), whereas the unsymmetric walk already misses Lambda9/Lambda10 most of the time.
  Small samples for n >= 20 (3-13 runs): frequencies there have large errors.

(2) Scaling. Median V at the first Q-perfect form (random rule) versus w D/2, 17 (family,k) points with n >= 10, w D/2 from 50
to 1800: ratio med V_perf / (w D/2) in [0.88, 1.59], geometric mean 1.06; log-log slope 0.94 (rms residual 0.13).
Examples: Hurwitz k=3,4,5: 0.95, 1.00, 0.99; (k=6, descent rule, 2 runs: 1.17); icosians k=2,3: 0.88, 1.02; Eisenstein k=5,6,8:
0.93, 0.96, 0.88; Z[zeta_12] k=3,4: 1.07, 0.99; Z[zeta_13] k=2: 1.16; Gaussian k=5,6: 1.15, 1.17. The outliers above 1.2 are
cases where the median is a "special" lattice (Z[zeta_7] k=2: 1.59, zeta_5 k=3: 1.23).  HEURISTIC reading: V_first-perfect
~ 1.0 * (w D/2) = (number of contact pairs forced by D independent orbits), the same law as V ~ N for Q = Sym_n (w=2, D=N).
Expressed through n: med V_perf / (n(n+1)/2) ranges from 1.0 (Z[zeta_8] k=4: 135.4/136; Gaussian k=5: 57.3/55) to 6.1 (icosians k=3: 1840/300) - the gain factor is the
predicted w D / n^2.  After LP descent the ratio V_ext/(w D/2) has geometric mean 1.5 and is dominated by special lattices
(up to 18 = Leech / 1800).

(3) Most efficient symmetry per dimension (median V_perf, random rule; = ordering by w D/2 in every dimension):
n=8: all families reach E8 with >= 68%; n=10: Eisenstein (69.4) > Gaussian (57.3); n=12: Hurwitz (170.9, but K12 never) >
Z[zeta_7] (133.4, K12 33%) > Z[zeta_12] (116.0) > Z[zeta_5] (111.0) > Eisenstein (104.1) > Gaussian (84.0) ~ Z[zeta_8] (83.5);
n=16: icosians (632) > Hurwitz (335) > Z[zeta_12] (190) > Eisenstein (169, 1 run) > Z[zeta_8] (135);
n=20: Hurwitz (534); n=24: icosians (1840) > Hurwitz (922, descent rule) > Z[zeta_13] (361).
For hitting the record lattice the best choice was: n=12 Z[zeta_7] (33%), n=16 icosians (42-67%), n=24 icosians (38-55%).

(4) Degeneracy. Fraction of first Q-perfect forms with exactly D contact orbits ("simple"): 0-20% in most families for
n <= 16 (mean number of orbits 1.2-2.5 D; e.g. Eisenstein k=6: 53.9 vs D=36, Hurwitz k=4: 70 vs 28, icosians k=2: 25 vs 12),
because the walk frequently lands on highly symmetric forms (E8 as Hurwitz form: 10 orbits vs D=6; BW16: 180 vs 28;
Leech as icosian form: 1638 orbits vs D=30). The generic (non-record) frozen forms in the largest cases are simple or nearly
simple: Hurwitz k=6: exactly 66 = D in 2 of 2 runs; Z[zeta_13] k=2: 24, 24, 25 vs D=24; Hurwitz k=5: 46-87 vs D=45.
So: typical frozen form ~ D orbits (V ~ w D/2), record lattices = extremely degenerate vertices.

## Explicit forms
* PROVED (elementary): every Q-perfect form is rational, because Q is cut out by integer linear equations and the form is the
  unique solution in Q of the integer linear system A[x]=1 (x in the contact set). Numerically q*A is integral with
  q = 2 (D4, E6, E8, D6), q = 4 (Lambda10, K12, Lambda12-type, BW16, Lambda20, Leech; det(4A) = 768, 729, 1024, 256, 64, 1).
* best_forms.json contains for the best form of every (family,k): q, the integer Gram matrix q*A in the Z-basis of O^k,
  kissing number, det(q*A), and for Eisenstein/Gaussian/Hurwitz the Hermitian matrix q*H (form = Re sum conj(x_i) x_j h_ij;
  reconstruction error <= 1e-15). Examples (components w.r.t. 1,i,j,k), copied from results_tables.md:
  E8 as Hurwitz 2x2: 2H = [[12, (1,-1,-4,-2)], [conj, 2]];
  BW16 as Hurwitz 4x4: 4H has diagonal (4,10,8,10) and off-diagonal entries with coordinates in {0,+-1,+-2,+-3,+-4};
  Lambda20 as Hurwitz 5x5: diagonal (6,6,4,8,12); K12 as Gaussian 6x6: 4H in M_6(Z[i]) with diagonal (4,16,12,10,8,4).
  These matrices are in the random basis produced by the walk (not reduced over O), so they are explicit but not canonical;
  the Leech form (icosians k=3) is stored as the 24x24 integer matrix 4A.

## Caveats / gaps
* Lattices are identified by invariants (V, number of minimal vectors, det of the integral form), not by an isometry test.
* Discrete sticky walk with exact line search and Reynolds-averaged Wishart start, not Klartag's continuous process.
* Only the free module O^k with scalar unit group; no non-free modules, no larger groups.
* n >= 20: 2-13 runs per case. Eisenstein k=8: 1 run. Hurwitz k=6: 2 runs, one truncated in the LP descent.
* Hurwitz k=3 never gave K12 (0 of 80). Not checked here whether K12 admits a Hurwitz-module structure at all.
* No theory in this round: the law V_perf ~ w D/2 is an empirical fit over w D/2 in [50, 1800].
