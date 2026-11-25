# R10 Weil lattices - Builder B notes

Start 08:24:48 UTC (hard limit 42 min). Independent implementation; directory A/ was not consulted.

## Method (exact-integer core)
* Model: odd part of the Weil rep on K^h, K=Q(zeta_p), S,T as in the spec; kappa = -1 (for p = 5,13,29,37,61) makes (ST)^3 = S^2 = -1.
* Explicit sublattice Lambda_B = Z[Borel] v0 = two copies of Z[zeta], Z-basis b_{j,c} = T^j w_c (j=1..p-1, c in {1,g}),
  w_c = 1_{c mu_m} - 1_{-c mu_m} (m=(p-1)/4, v0 = w_1). Gram_B(Lambda_B) = blockdiag(m(pI-J)) exactly, inverse (I+J)/(mp).
* S, T, eta = sqrt(p) act on these coordinates by N_S/(mp), M_T, M_eta. The integer matrices were obtained by rounding numerics
  (residual < 1e-12) and then ALL relations were verified in exact integer arithmetic (python ints): S^2=-1, (ST)^3=-1, T^p=1,
  S and T isometries of Gram_B, eta^2 = p, eta commutes with S,T, E = eta-twisted Gram symmetric.
* Lambda_0 = Z[G]v0 = smallest S,T-stable lattice containing Lambda_B. Since B is Z-valued on Lambda_0, (mp)Lambda_0 is inside Lambda_B,
  so the closure is an exact HNF computation modulo mp (weil.py: hnf_mod / Weil.build; 1-3 closure iterations).
  B, E (exact integer Gram matrices of B(u,v) and B(u,eta v)) are saved in an LLL-reduced basis: B_p*.npy, E_p*.npy;
  Ared_p*.npy = that basis in (mp)-scaled Lambda_B coordinates.
* Gamma(w) = (b^2 - e^2/p)/(4D), b = wBw, e = wEw. Violator test = exact integer inequality  p b^2 - e^2 < D (p-1)^2  (D = m^2 p integer).
* Fixed sublattices of split-torus subgroups: exact (HNF mod mp in orbit-difference coordinates). Words in S,T (non-split torus elements) and
  unit-twisted eigenlattices: LLL kernel, then exact verification of kernel membership and rank.
* Order: delta_u (u in Lambda^H) built explicitly as integer matrices on the saved basis (extras.py): delta_{v0}=1, commute with S,T exactly,
  preserve Lambda_0 exactly. Exact R_U of the Delta-line through any w: det(O w) (exact Bareiss; equals (N/16)^4, N = p b^2 - e^2, verified)
  and saturation index [Lambda cap Delta w : O w] via Smith diagonalisation: R_U = (theta/Gamma(w))^4 * index^2.
* Exhaustive searches: own LLL (exact int64 Gram, long double GSO) + Fincke-Pohst (lat.c) for every t on a grid of spacing 0.1 over
  [-log eps, log eps], bound 2 sqrt(M theta D) cosh(0.05)(1+1e-9) (M=1 for the plain violator test). Floating point (double) Cholesky; the
  Gamma test of every enumerated vector is exact.
* verify.py: independent end-to-end check: BFS orbit of v0 in the complex model (size |G|/|H| = 4p(p+1) => finite group), coordinates of all
  orbit vectors in the saved basis are integral and their Z-span is the whole saved lattice; numeric Gram from the embeddings equals B, E
  (1e-12); exact det(B) = (m^2 p)^(p-1); certificates recomputed from m_plus, m_minus directly.

## Results
All of p = 5, 13, 29, 37, 61: rank 2(p-1); det_B(Lambda_0) = (m^2 p)^(p-1) exactly, D = m^2 p; O = End_G(Lambda_0) is a MAXIMAL order
(det_B(Lambda_0^H) = p^4 m^16 exactly; no enlargement needed); the split-torus line Lambda^H is exactly borderline (R_U = 1);
theta_p D = 4 m^4 = m_+(v0) m_-(v0), i.e. v0 is a border vector. Units of reduced norm 1 in O: 24 (mod +-: 12, eight of order 3) for p=13,37,61.

* p=5: min Gamma = 0.8 = theta (exhaustive, rank 8). E_8 (even unimodular, min 2, 120 pairs) is the member at t = +-log(eps);
  at alpha=beta (t=0) the det-1 form is NOT integral, min 4/sqrt5 = 1.789 (60 pairs).  [spec discrepancy]
* p=13: min Gamma = theta_13 = 36/13 exactly; exhaustive on the full rank-24 lattice: no violator. Leech (even unimodular, min 4,
  98280 pairs = 196560 vectors) at t = +-log(eps); at t=0 the det-1 form has min 36/sqrt(117) = 3.328 (19656 pairs).
* p=29: unstable. Non-split order-5 element S T^5: fixed lattice = Delta-line with exact R_U = (49/25)^4 = 14.7579 (this line is NOT
  detected by single-vector Gamma: min Gamma on it = (100/49) theta, saturation index 16). Non-split order-3 element (ST)^2: rank 16,
  R_U = 11.780; it contains an explicit Gamma-violator b=294, e=1218, Gamma = (45/49) theta, exact line R_U = (49/45)^4 = 1.4058
  (violator_p29.json).
* p=37: **TENSOR-UNSTABLE (proved)**. The split C_3-fixed sublattice (rank 24, itself exactly borderline as a Delta-subspace) contains
  w with b = 1026, e = 5994: N = p b^2 - e^2 = 3020976 < D (p-1)^2 = 3884112, Gamma = (7/9) theta = 6.8108. Exact: O w is saturated,
  det_B(Lambda cap Delta w) = 188811^4 and R_U = (9/7)^4 = 6561/2401 = 2.7326. Certificate violator_p37_C3.json (re-verified by verify.py).
  Exhaustive on that rank-24 sublattice: min Gamma = (7/9) theta exactly; 72 violating vectors modulo +-eps^Z (all with the same N),
  12 border vectors. The unit-twisted eigenlattices ker(g_3 - delta_3) (rank 24) also contain violators with the same value 7/9.
  Split C_9 line: exactly borderline, no violator. Full-lattice LLL heuristics do NOT find these violators (best 1.69 theta): the
  heuristic has no power.
* p=61: NO violator found.
  - split C_15 line (rank 8), C_5-fixed (rank 24), C_3-fixed (rank 40): each exactly borderline (R_U = 1 exact determinants) and
    exhaustive enumeration over all 74 t-points completes with min Gamma = theta exactly (only border vectors).
  - unit-twisted eigenlattices ker(g_3 - delta_3) (rank 40, two different delta of order 3): exhaustive, no vector with Gamma <= theta.
  - exact line test in the C_5-fixed lattice: all 5964 vectors (mod +-eps^Z) with Gamma <= 4 theta generate lines with exact R_U <= 1.
  - full rank-120 lattice: LLL over the t-grid + pair sums, and + images under 30 random group elements: nothing below theta
    (best new vectors 2.35 theta). This is uninformative (same heuristic misses the known p=37 violators).

## Caveats
* "min Gamma >= theta" is NOT equivalent to line-stability: non-principal lines (Lambda cap U strictly bigger than O w) need the
  saturation index (p=29 order-5 line is such an example). Done exactly only where stated above.
* p=61 full lattice (rank 120) is far out of reach of exhaustive search; no BKZ was implemented. p=61 remains OPEN in general.
* Enumeration completeness relies on double-precision Fincke-Pohst (integer Gram data, 1e-9 slack), not interval arithmetic.
* p=17 (p = 1 mod 8) not implemented.
