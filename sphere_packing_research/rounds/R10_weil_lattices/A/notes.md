# R10 Weil lattices -- Builder A notes (running)

Start 08:24:47 UTC, hard limit 09:06:47 UTC.

## Code
- `weil.py`  : generic pipeline (Weil rep odd part in all embeddings, orbit BFS, exact Z-basis of Lambda_0 = Z[G]v0,
               integer Gram matrices B, E, exact det, order test, fixed lattices, t-grid Fincke-Pohst, Gamma).
- `lat.c`    : floating LLL (modified Gram-Schmidt on actual vectors, exact int64 transformation) + BKZ with enumeration.
- `fp_enum.c`: Fincke-Pohst (copied from R2).
- `run1.py p`: full pipeline for one prime; `search_full.py`: BKZ search on the full lattice over the t-grid;
  `verify37.py p`: independent exact verification of a violator incl. explicit rank-8 sublattice O.w.

## How the Z-basis is obtained (exactness)
orbit vectors (floats, 1e-14) -> choose 4h independent orbit vectors M (pivoted QR) -> integer Gram G_M (rounded, residual 1e-12)
-> exact adjugate (fraction-free Gauss-Jordan, verified G_M*adj = det*I exactly) -> exact M-coordinates of ~4n orbit vectors
-> exact HNF modulo det(G_M) (python ints) -> chain of prime-index steps with exact modular solves + float LLL (unimodular, int64)
-> all orbit vectors have integer coordinates (float check, <1e-10), S and T map the basis to integer unimodular combinations,
and EXACT identity det(B) * index^2 = det(G_M).  Hence the saved basis is a Z-basis of Z[G]v0.

## IMPORTANT bug found and fixed during the round
The first version called the C routines through ctypes WITHOUT argtypes; numpy pointers were truncated to 32 bit and
several runs segfaulted (dmesg) and one run printed a bogus "min Gamma = 0" for p=61 C_3.  After declaring argtypes
(c_void_p) ALL runs were repeated; every number quoted below is from the re-runs (logs `log_p*.txt`, `log_C3_*.txt`,
`log_full13.txt`, `log_verify*.txt`).  All violator certificates are checked with exact Python integers anyway.

## Results (Lambda_0 = Z[G] v0, B = Tr_{K/Q}(u^* v), E = B(u, sqrt(p) v))
All five lattices: rank 2(p-1) (pivoted-QR gap 1e13), B and E integral (residual <1e-9), E B^-1 E = p B, O_F acts,
order O(Lambda_0) MAXIMAL (det of Tr nrd form on Lambda_0^H = (p/4)^4), det B = D^(p-1) EXACTLY with
  p=5: D=5 ; p=13: D=117=9*13 ; p=17: D=153=9*17 ; p=29: D=1421=49*29 ; p=37: D=2997=81*37 ; p=61: D=13725=225*61.
Gamma(w) = (p b^2 - e^2)/(4 p D), b = x B x^T, e = x E x^T (exact integers); violator <=> p b^2 - e^2 < (p-1)^2 D.

| p  | rank | min Gamma found / theta | status |
|----|------|-------------------------|--------|
| 5  | 8    | 1 (exhaustive, full lattice)                 | borderline (E_8), validation OK |
| 13 | 24   | 1 (exhaustive, FULL lattice, 19656 vectors enumerated, exact min N = THR = 16848) | exactly borderline, validation OK |
| 17 | 32   | 0.5625 = 9/16 on the g3-fixed line (v0)     | UNSTABLE, R_U = 9.98872 = (16/9)^4 (direct exact det of Lambda^H) |
| 29 | 56   | 5/7 (BKZ search, NOT claimed minimal)        | UNSTABLE, explicit rank-8 sublattice O.w: R_U >= 3.8416 = (7/5)^4 |
| 37 | 72   | 7/9 = Gamma 252/37 vs theta 324/37           | UNSTABLE, explicit rank-8 sublattice O.w: R_U >= 2.73261 = (9/7)^4 |
| 61 | 120  | 1 on C_15 line, C_5-fixed (rank 24), C_3-fixed (rank 40), all exhaustive | no violator found; full lattice NOT decided |

### p = 37 (the point of the round): UNSTABLE
- violator w (file `violator_p37.json`, integer coordinates x in the saved basis `BasC_p37_L0.npy`, Gram `B_p37_L0.npy`, `E_p37_L0.npy`):
  b = 1026, e = 5994, N = 37*1026^2 - 5994^2 = 3020976 < (p-1)^2 D = 3884112.  Gamma = 252/37 = 6.8108 < theta = 324/37 = 8.7568.
  m_plus = 1005.70, m_minus = 20.30 (cross-checked directly from the Minkowski embedding).
- independent certificate (`verify37.py`): the 8 vectors delta_k(w) (delta_k = Z-basis of O, realised through delta_u(g v0) = g u,
  G-equivariance checked to 3e-13) have INTEGER coordinates (1e-10) in Lambda_0, span a rank-8 lattice O.w containing w, with
  EXACT det_B(O.w) = 1270893532472362625841 = (N/16)^4, while r^8 D^4 = 3472858170158755180401.  R_U >= 2.7326 > 1.
- (7a) C_9 line (rank 8): exhaustive, min Gamma = theta (borderline, R_U = 1 exactly: det_B(Lambda^H) = r^8 D^4 exact identity).
- (7a) C_3-fixed lattice (rank 24): R_U = 1 exactly as a subspace (exact identity), BUT exhaustive t-grid Fincke-Pohst finds
  732 vectors below the bound, 143 distinct violators (+-pairs), min Gamma = 252/37 (ratio 7/9).  So the violating Delta-lines
  already live inside the C_3-fixed 24-dimensional sublattice, and 7/9 is the exact minimum of Gamma/theta on that sublattice.
- (7b) full lattice BKZ(20/30) over the t-grid: best found also ratio 7/9; no smaller Gamma found (not exhaustive).

### p = 61 (evidence only)
- lattice built and verified as above (index chain had a 32-digit prime; exact modular solves).  Order maximal (float det ratio 1.00000000).
- exhaustive on C_15 line, C_5-fixed (rank 24), C_3-fixed (rank 40, 2.9e8 nodes): min Gamma = theta exactly, 0 violators.
- full lattice (rank 120): BKZ-20/30 reaches |b1|^2 ~ 1000-5500 while the violation threshold is 900 = 1.04 * GH^2: generic
  reduction cannot decide; nothing found below theta (best candidates are the orbit of v0 with Gamma = theta).  NOT DECIDED.

## Discrepancies with the orchestrator's specification
1. With S, T as specified and kappa = +1: S^2 = -I and (ST)^3 = +I (not S^2).  kappa = -1 gives (ST)^3 = -I = S^2.  Both signs
   generate the same group (since -I = S^2), so nothing depends on it.  (Equivalent to T -> [[1,-1],[0,1]].)
2. "The form of determinant 1 with alpha = beta is E_8 / Leech" is NOT true for Lambda_0 with B = Tr(u^* v): det B = D^(p-1), D = p k^2,
   and alpha = beta (t = 0) is the member where the Gamma-minimal vectors (the orbit of v0) are shortest: normalised minimum
   2 sqrt(theta_p) (1.789 for p=5, 3.328 for p=13).  E_8 is the member at t = log(eps) (checked: det 1, min 2, 240 vectors).
   For p=13 the analogous check at t = log(eps) is in `log_unimod13.txt` (see below).  R_U is unaffected.
3. theta_p, the formula det_B(O w) = ((p/4) m+ m-)^4 and the calibration "split-torus fixed lattices are exactly borderline"
   were all CONFIRMED (exact integer identities for p=13, 29, 37; p=17 direct R_U = (16/9)^4).
4. A plain float approach for the Z-basis fails (index [Lambda_0 : Z-span of 4h orbit vectors] has huge prime factors, e.g. 2532462127963
   for p=37); exact HNF + prime-index chain with modified Gram-Schmidt LLL was needed.

## Late additions (09:00 UTC)
- p=13: determinant-1 member at t = log(eps) = 1.1948 has minimum 4 with 98280 pairs (= 196560 vectors): Leech (`log_unimod13.txt`).
  p=5: determinant-1 member at t = log(eps) = 0.4812 has minimum 2 with 120 pairs: E_8 (`log_unimod5.txt`).
- p=61: EXACT order test (det_B(Lambda^H) * 2^16 == p^4 * B(v0,v0)^8) True, exact calibration det_B(Lambda^H) == r^8 D^4 True (`log_verify61.txt`).
  `candidate_p61_NOT_a_violator.json` is just the best candidate of the BKZ search (an orbit vector, Gamma = theta); it is NOT a violator.
- p=37 deeper search (BKZ 30/40, 6 tours, 3 processes x 150 s, `log_searchdeep37_*.txt`): smallest Gamma found is still 252/37
  (ratio 7/9); smallest |b1|^2 over the t-grid = 286.0 (threshold 324).  Whether Gamma_min of the FULL rank-72 lattice is exactly
  252/37 is NOT proved (only proved on the C_3-fixed rank-24 sublattice, by exhaustive enumeration).
- p=29: the violator found (ratio 5/7, R_U >= 3.8416) is weaker than the one predicted from the non-split torus (R_U >= 11.8);
  the non-split-torus fixed spaces were not computed in this round.

## Exact vs floating point
EXACT (python integers): det B, D^(p-1) = det B, det(B)*index^2 = det(G_M), order test and torus calibration for p = 29, 37, 61,
  N = p b^2 - e^2 and the threshold comparison for every certificate, det_B(O.w) = (N/16)^4 for p = 29, 37, exact min N for
  p=13 full lattice, p=37 C_3, p=61 C_3 (`full13.py`: every enumerated vector evaluated with exact integers).
FLOATING POINT: the Weil matrices / orbit (1e-14), integrality of Gram matrices (rounded, residual <= 1e-8), integrality of
  coordinates (<= 1e-8), Fincke-Pohst enumeration (double Cholesky, slack 1e-9 in the bound, node flags 0 = complete),
  order test / calibration for p = 5, 13, 17 (ratio 1.00000000 printed), R_U(p=17) = 9.98872123.
NOT DONE: maximal-order enlargement code (never needed: O(Z[G]v0) was maximal in all six cases); non-principal lines
  (ideal-norm version); p=61 full lattice; p=37/29 global minimum of Gamma; counting violators modulo units/G.
