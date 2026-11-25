# R12_p61 - TESTER notes (tensor stability on group-theoretically structured Delta-subspaces)

Start 08:49:15 UTC, hard limit 30 min. Files: iso.py (main), extra.py (p=29 cross-check, second SL_2(5) class, sums),
extra2.py (joint components for subgroup chains), extra3.py (unimodular member inside structured sublattices),
logs log_*.txt, json iso_p*.json, search_p*.json, extra*_p*.json, unimod_sub_p*.json. weil.py/lat.c/lat.so = copies of builder B's code.

## Method (all R_U values below are EXACT integer computations unless said otherwise)
* Lattice: builder B's saved exact basis (Ared/B/E_p*.npy), re-verified here: Gram consistency exact, S,T act integrally (exact),
  S^2=-1,(ST)^3=-1,T^p=1 exact on Lambda, S,T isometries of B and E exact, float logdet(B) = (p-1) log D.
* Group elements: g in SL_2(F_p) -> integer matrix on Lambda via g = T^x (S T^-c S^-1) T^z; homomorphism property checked exactly
  on random products. Subgroups: split torus <d(g)>, non-split torus (order p+1 element), unipotent, SL_2(3)=<s,t|(st)^2=s^3=t^3=-1>,
  SL_2(5)=<s,t|(st)^2=s^3=t^5=-1> found by brute force in SL_2(F_p); defining relations re-verified exactly on the operators; group orders 24/120 by closure.
* Projectors: e_U = (chi(1)/|H|) sum_h Tr_{Q(chi)/Q} chi(h) h (class function by element order; cyclic groups: Ramanujan sums).
  Integer matrix Et = N e_U; checked exactly: Et^2 = N Et, B- and E-self-adjoint (so U is eta-stable and e_U is the ORTHOGONAL projector), trace = N dim U.
* KEY FORMULA (proved + checked): Lambda^# (w.r.t. B) = (eta/(r p)) Lambda  [exact check: X = r B^-1 E B^-1 is an integer matrix, B X B = r E exactly, |det X| = 1].
  Hence for every eta-stable U:  det_B(Lambda cap U) = D^{k/2} * g_U,  g_U := [Lambda : (Lambda cap U) + (Lambda cap U^perp)] (glue index), k = dim U, and
        R_U = (n/k)^k / g_U,   g_U = |rowspace(Et) mod N| = N^n / prod diag HNF_modN(Et)   (exact).
  Stable on U  <=>  g_U >= (n/k)^k.   (Same g for U and U^perp.)
* Independent certificate: for every component of rank <= 40 the saturated lattice Lambda cap U was ALSO computed as an integer kernel (LLL) and its exact
  Gram determinant equals g_U D^{k/2} (all 'exact=True' in logs) -> glue formula confirmed in 30+ cases.
* Cross-validation p=29 (log_extra_p29.txt): non-split C_30 eigen-order {6} (= fixed lattice of order-5 element, dim 8): R_U = 14.7579 = (49/25)^4;
  eigen-order {10} (= fixed lattice of order 3, dim 16): R_U = 11.7801. Both equal builder B's exact values -> the routine DOES detect instabilities.
* Exhaustive Gamma search: builder B's tsearch (LLL + Fincke-Pohst over t-grid 0.1, radius 2 sqrt(theta D) cosh(0.05)); complete => every w with Gamma(w) <= theta is found.

## Dimension remark
n = 2(p-1) = 120 for p=61, W_C|_{SL_2(5)} = 2 x (odd regular) => icosian component has dim 16 (Delta-rank 2), NOT 32; degree-4 comp. dim 32, degree-6 dim 72.
The predicted formula (625/256)^2 2^-e is the dim-16 formula; result: e = 8.

## Results p=13 (n=24) - validation: all R_U <= 1
split C_12: eigen-order {4} (=Fix C_3, dim 8) R=1 exactly (g=3^8); {12} (dim 16) R=0.1001.  SL_2(3): Hurwitz comp (dim 8) R=1 exactly; 2'+2'' (dim 16) R=0.1001.
Chain SL_2(3)>C_6: three joint components of dim 8, each R=1. Non-split C_14, unipotent C_13: W|_H Q-isotypic (free). Exhaustive Gamma search on all four
components: min Gamma = theta exactly, no violator.

## Results p=37 (n=72) - instability rediscovered
split C_36: {4} dim 8 R=1; {12} dim 16 R=0.1001; {36} dim 48 R=1.0e-3; {4,12}=Fix(C_3) dim 24 R=1 (g=3^24); other sums <<1.
SL_2(3): Hurwitz comp dim 24: R = 2^-8 (g=2^8 3^24); 2'+2'' dim 48: R=3.9e-6.  => NO isotypic component has R_U>1 at p=37 either!
But exhaustive Gamma search: Fix(C_3) (rank 24): min Gamma/theta = 7/9 (b=1026,e=5994; 1082 hits); Hurwitz component (rank 24): ALSO min Gamma/theta = 7/9
(b=2970,e=17982; 392 hits) -> violating Delta-LINES inside the components (R_line >= (9/7)^4 = 2.73). In the unimodular member: 72 pairs of norm 6, s=+-9 in Fix(C_3),
24 pairs of norm 6, s=+-9 in the Hurwitz component.  LESSON: isotypic-component R_U <= 1 is necessary but far from sufficient; the instability at p=37
lives in Delta-lines inside (H-stable) components, detected only by the vector search.

## Results p=61 (n=120, r=15, theta=14.754, D=13725)  - ALL structured tests passed
R_U exact (g = glue index):
* split C_60 atoms (covers all split cyclic subgroups, and by character theory all dicyclic subgroups of N(T_split) incl. Dic_3, Dic_5, Dic_15, Dic_30):
  {4}=Fix(C_15) dim 8 g=3^8 5^8 R=1 | {12} dim 16 g=3^8 5^16 R=0.1001 | {20} dim 32 g=3^32 5^8 R=3.2e-3 | {60} dim 64 R=1.0e-9
  sums: {4,12}=Fix(C_5) dim 24 g=5^24 R=1 | {4,20}=Fix(C_3) dim 40 g=3^40 R=1 | {4,60} 1.3e-20 | {12,20} 1.7e-17 | {12,60} 1.0e-5 | {20,60} 3.4e-8 |
  {4,12,20} 1.2e-8 | {4,12,60} 9.9e-10 | {4,20,60} 2.9e-9 | {12,20,60} 8.9e-7.
* SL_2(3): Hurwitz comp dim 40 g=2^16 3^40 R = 2^-16 = 1.5e-5 (e'=16) | 2'+2'' dim 80 R=1.5e-10.
* SL_2(5) (both conjugacy classes, identical): icosian comp dim 16 g=2^8 3^16 5^8 R = 5^8/2^24 = 0.02328 (e=8 >= 3: STABLE) |
  degree-4 dim 32 g=3^16 5^32 R=2.3e-12 | degree-6 dim 72 g=2^8 5^40 R=4.0e-15 | pairwise sums 5.4e-12, 7.1e-19, 6.8e-10.
* chains SL_2(5)>SL_2(3)>C_6, SL_2(5)>C_10, SL_2(5)>C_6 (joint components, extra2): max R = 1 (attained only at Fix(C_5) and Fix(C_3)), all others <= 1.5e-5.
* non-split C_62, unipotent C_61, C_4: W|_H Q-isotypic (computed) -> no proper component. Borel-type C_61:C_d, Dic_31, Q_8: Q-isotypic by character theory (NOT computed).
Exhaustive Gamma search (complete t-grid unless said): icosian comp (rank 16): no vector with Gamma <= theta (best seen 2.44 theta); degree-4 comp (rank 32): none
(best seen 6.78); Hurwitz comp (rank 40): none (best seen 1.44; 1.6e7 nodes); {4} rank 8: min = theta (border), no violator; {12} rank 16: none (1.333);
{4,12} rank 24: min = theta, no violator; {20} rank 32: none (1.44); {4,20}=Fix(C_3) rank 40: my run INCOMPLETE (66/74 t-points, node cap), no violator
(builder B's earlier run of the same lattice was complete: no violator).

## Open / caveats (honest)
* Tested = isotypic components (and sums/joint components for chains) for ALL conjugacy classes of subgroups of SL_2(61) (cyclic, dicyclic, Borel-type,
  SL_2(3), SL_2(5); no 2.S_4 since 61 = 5 mod 8). All R_U <= 1. This is a NECESSARY condition only.
* The vector test covers only Delta-lines that contain a vector with Gamma < theta; for non-principal O-ideal classes (class number of the maximal order
  is > 1 for p=37,61: mass = zeta_F(-1)/2 = 11/12 for p=61) a line can have R_U > 1 with min Gamma >= theta.  Not tested.
* H-stable Delta-subspaces inside components of multiplicity > 1 form continuous families (e.g. P^4(F) of SL_2(3)-stable Delta-lines in the Hurwitz comp);
  only their vectors with Gamma <= theta were excluded (exhaustively, rank <= 40). Degree-6 comp (72), 2'+2'' (80) etc. too large to enumerate.
* Generic Delta-subspaces of the full rank-120 lattice: untouched here.

## Task (4): even unimodular member (t = +-log eps), q = norm, s = (m~_+ - m~_-) sqrt(p)/2, Gamma = (q^2 - 4 s^2/p)/4
Exact check: Gram at t=+log eps integral (residual 2e-13), even diagonal, float logdet -4e-7 (p=61).
(a) EXHAUSTIVE inside structured sublattices (extra3.py, unimod_sub_p61.json), float FP enumeration at t=+log eps and t=-log eps:
  Fix(C_15) rank 8: no vector of norm <= 14 | icosian comp rank 16: min 12 (60 pairs, s=0) | {12} rank 16: none <= 14 |
  Fix(C_5) rank 24: min 10 (360 pairs: |s|=5 x288, |s|=15 x72; Gamma/theta = 1.667 / 1.444) | degree-4 comp rank 32: none <= 14 |
  Hurwitz comp rank 40: min 12 (1368 pairs, |s| in {0,6,12,18,24,30}, min Gamma/theta 1.44) | Fix(C_3) rank 40: min 12 (19800 pairs, same |s| set, 1.44).
  => smallest norm found in the 120-dim even unimodular lattice: 10 (so it is NOT extremal; min <= 10). No vector of norm <= 8 in any structured sublattice.
  p=37 control: norm 6, s=+-9 (72 pairs in Fix(C_3), 24 pairs in the Hurwitz comp), Gamma/theta = 7/9 -> rediscovered.
(b) HEURISTIC on the full lattice (extra4.py, unimod_full_orbit_p61.json): the 360 norm-10 pairs fall into 30 G-orbits of 22692 pairs each
  (24 orbits with |s|=5, 6 orbits with |s|=15; total 680760 pairs = 1361520 vectors of norm 10). Inner products of each orbit representative with all vectors
  of all 30 orbits: |ip| in {0,1,2,3,4,5} only (aggregate 4681860/4969200/858330/40800/720/840) -> no difference of two of these vectors has norm < 10
  (|ip|=5 gives further norm-10 vectors, |s| in {5,15}). Best Gamma/theta seen on the full lattice by this route: 1.444. No vector of norm <= 8 found.
  This is weak evidence only (pair differences of one family of minimal vectors); a genuine SVP search in dim 120 was NOT done.

## Conclusion p=61
Every group-theoretically structured necessary condition tested is satisfied: R_U <= 1 for all isotypic components / sums / chain-joint components of all
subgroup classes (equality exactly at Fix(C_15), Fix(C_5), Fix(C_3)); the decisive SL_2(5) icosian component has R = 5^8/2^24 = 0.0233 (2-adic glue exponent e = 8),
Hurwitz component R = 2^-16; exhaustive Gamma search finds no violator in any structured sublattice of rank <= 40. p=61 remains OPEN but un-refuted:
no instability found; not a proof of stability (see caveats above). Finished 09:08 UTC.

## Addendum 09:10 UTC - Fix(C_3) (rank 40, p=61) exhaustive search now COMPLETE in this round too
The first pass (log_p61.txt) stopped at 66/74 t-points on the time limit. Re-run in two pieces with per-point logging:
t-points 0..63 (log_extra5b_p61.txt, search_p61_fixC3_head.json): complete=True, 1.33e9 nodes, min Gamma = theta exactly (b=900,e=0), 0 violators, 13 border hits;
t-points 64..73 (log_extra5_p61.txt, search_p61_fixC3_tail.json): complete=True, 0 violators. => supersedes the 'INCOMPLETE' remark above; agrees with builder B.
