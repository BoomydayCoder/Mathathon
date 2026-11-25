# R12_p61/verify - ADVERSARIAL VERIFIER notes (start 09:12 UTC, limit 12 min)
Scripts: v.py (p = 13, 37, 61), v2.py (follow-ups), v3.py (split C_60 atoms). Logs: log_v13/37/61.txt, log_v2_37/61.txt, log_v3_61.txt.
Own code: group words (M = T^x S T^c S T^z), presentation checks on operators, parallel closure, projectors, exact Bareiss det,
det mod 3 primes, rank mod q (primitivity). Reused from tester: setup() (builder B's saved B, E, S, T), hnf_mod, sat_basis (only to PRODUCE a basis).
Certificate per component (independent of the glue formula): C integer k x n, C(N - Et) = 0 exactly (in U), rank_q(C) = k for every prime q | det
(=> C is primitive => C spans Lambda cap U exactly, saturated), exact det(C B C^T), R_U as exact fraction.

## Verified computation (exact)
* p = 13, 37, 61: det(B) = D^(n/2) EXACTLY (Bareiss + 3 modular checks); eta = E B^-1 integral, eta^2 = p, commutes with S, T; (1+eta)/2 integral
  (Lambda is an O_F-module, so the t-period 2 log eps is legitimate); B X = r eta with X integral => Lambda^# = eta/(rp) Lambda. tr(S) = 0, tr(T) = -2.
* p = 61: Fix(C_3) dim 40 R = 1; Fix(C_5) dim 24 R = 1; Fix(C_15) dim 8 R = 1; SL_2(3) Hurwitz dim 40 R = 2^-16 (two SL_2-classes + the one inside SL_2(5));
  SL_2(5) icosian dim 16 R = 5^8/2^24 = 390625/16777216 (both classes); degree-4 dim 32 R = 3^16/2^64 (both classes);
  ico+deg4 dim 48 R = 5^8/2^56 = 5.42e-12 (=> deg-6 glue 2^8 5^40 confirmed via complement); eigen-order {12} dim 16 R = 3^8/2^16.
  All bases primitive; all agree with the tester's glue index.
* p = 37: Fix(C_3) dim 24 R = 1; Fix(C_9) dim 8 R = 1; Hurwitz dim 24 R = 2^-8. Saved violators: b=1026, e=5994, Gamma/theta = 7/9 exactly, lies in Fix(C_3);
  b=2970, e=17982, Gamma/theta = 7/9 exactly, lies in the Hurwitz component of the tester's SL_2(3) (gens #0), not in other conjugates (as expected).
* p = 13: Fix(C_3) = dim 8, R = 1; Hurwitz dim 8, R = 1.
## Logic
* U is a Delta-subspace: Et is a Z-combination of words in S, T, so it commutes with End_G(W); eta-commutation also checked exactly.
* Saturation: the tester's R_U are for the SATURATED lattice (glue formula is a theorem given the exact modularity; re-derived). No claim rests on a finite-index sublattice.
* icosian dim = 16 (W_C|H = 2 x odd-regular, 2*(4+4)); the task's "32" was wrong; tester's correction is right.
## NOT re-verified here
* Float Fincke-Pohst Gamma searches (logs inspected only: complete flags, truncation flag semantics in lat.c ok, grid covers a full period).
* Atoms {20}, {60}, remaining sums, chain-joint components, non-split/unipotent isotypy (only the tester's computation; character theory agrees).
* Identity of the saved lattice with Z[G] v0 (inherited from builder B).
* ADDENDUM 09:21: eigen-order {20} dim 32 certified: R = 5^24/2^64 = 0.00323 (primitive basis, exact det). {12,20} and {60} certificates (v3.py) still running at the deadline: NOT verified here. Via g(U) = g(U^perp): {4,12,60}, {4,20,60}, {12,60}, {20,60}, {12,20,60}, deg-6, 2'+2'', deg4+deg6, ico+deg6 follow from certified complements and reproduce the tester's values.
