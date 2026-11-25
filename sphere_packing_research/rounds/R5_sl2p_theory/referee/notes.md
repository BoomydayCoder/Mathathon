# Referee report on R5_sl2p_theory (07:42 - 07:58 UTC)

Code (all exact integer / rational arithmetic): `exact_enum.py` (circulant PSD enumeration, exact Gaussian elimination over Q),
`enum17_exact.py` -> `enum17_exact_out.txt`, `leech_check.py`, `group_check.py` (use the arrays of ../../R5_leech13 only as data whose
defining properties are re-verified: even unimodular basis with det 2^36, 196560 norm-4 vectors => Leech).

## Verdict: no claim refuted. Three wording/logic gaps, one upgrade, one strengthening.

Confirmed by re-derivation (by hand): Lemma A' (glue group L/(L_1+L_z) is killed by q and by 1-g, so F_q^c, c <= min(a,b); convexity of q^x),
Cor. A'' ((a,b) = (4(N-1),4N); character values checked by <eta,eta> = 1, <eta,1> = 0, p = 5 and the Leech traces -2, -4, 0),
Lemma B (trace-1 element x = -sum zeta^j; (W^H)^# = pi_H(W) = (1/h) N_H(W) = (1/h) W^H), Lemma 1 (det(wJ) = N(nu(w) nrd J)^4; planes for r = 3),
Lemma 2, Lemma 3 (a),(b), Lemma 4, mass values, the formula R = |H^0_Tate|/(1+sum rho)^d, the Siegel-type volume d omega_d^2 c^{d/2} log eps,
all quoted numbers (9.99, 11.8, 14.8, 12.5, 17.3, 13.0, 24.6; survivors list; 3.35 points per borderline line at k = 4; 1.380/5.52, 2.790/33.5).

Confirmed by own exact computation:
* p = 13 enumeration, exact: 1506 sequences with sum -2; PSD: s = 0: 48, +-2: 33, +-4: 18, |s| >= 6: 0 (of 90+90, 30+30, 6+6). The floating
  point caveat in the notes is removed.
* p = 17 (mu = 4), exact, |s| >= 6 only: 14 PSD sequences each for s = +-6, none for |s| >= 8 (as in enum17_out.txt).
* Leech: s-distribution over all 196560 minimal vectors = {0: 78624, +-2: 45864, +-4: 13104}, no s = +-6; E s^2 = 4, E s^4 = 41.6; 126 orbit types.
  Hence N(nu(w)) >= 9 for all w != 0 (N in {1,3} has a norm-2 vector in its unit orbit: nu = 1, (5-sqrt13)/2, 4-sqrt13 give |w|^2 = 2;
  N = 4 <-> minimal vector with s = +-6).
* O(L) = {x in H (x) Q(sqrt13) : xL in L} contains O_0 = Hurwitz (x) Z[(1+eta)/2] with index exactly 4 (4 of the 256 classes of (1/2)O_0/O_0
  preserve L), N(discrd O_0) = 4, so O(L) is a MAXIMAL order. => the hypothesis of the "conditional" p = 13 claim holds.
* <g, h> (h = stored Delta-linear automorphism) has order 2184, traces {24:1, -24:1, -2:168, 2:168, -4:468, 4:468, 0:910}, sum tr^2/|G| = 8
  (commutant = Delta exactly), element orders as in SL_2(13); every element commutes with i, j, eta.
* An order-3 element c: N_c(L) is primitive of rank 8 (Smith invariants 1), det = 6561 = 3^8, R_U = 1 exactly (Lemma B equality case).

Upgrade: "Leech is exactly borderline over Delta_13" is a computer-assisted THEOREM (lines: N(n_U) >= 9 without class number, by Hermite:
det >= 2^8 excludes N(n_U) in {1,3}; N(n_U) = 4 would give a line with det 2^8 containing a minimal vector w' with [M_U : Ow']^2 = N(nu(w'))^4/2^8
not an integer square for N(nu(w')) in {9,12,13}; planes by duality; equality on the s = +-4 lines and the sqrt3 E_8 lines).

Gaps / corrections:
1. Class number one for p = 13: "[O^x : O_F^x] >= 12, so mass 1/12 forces h = 1" is not a valid inference as written (a larger unit group
   would give a smaller term and need more classes). Fix: O^1 is a finite subgroup of H^x containing SL_2(3), hence 2T, 2O or 2I; 2O needs
   sqrt2, 2I needs sqrt5, neither in Q(sqrt13); so |O^1| = 24, the term is exactly 1/12 and h = 1. (Not needed for the borderline theorem.)
2. Lemma 3(c) "for no p >= 13 does the minimum suffice" only shows that the two general inequalities N >= mu^2/(4b^2), N >= mu^2/4 do not reach r^2.
   As a logical statement it is false for p = 13 (unimodular + min 4 => Leech => line condition holds).
3. Context: Delta = (-1,-1)_F only for p = 5 mod 8 (2 inert). For p = 1 mod 8 (17, 41, 73, ...) (-1,-1)_F is ramified at the two primes
   above 2 and is NOT End_G(W). The notes use "ramified exactly at the two real places", so no proof is affected.
4. "Tensor-unstable" (R_U > 1) is an obstruction to first-moment / Klartag arguments on that family (T_U exponentially large), not a theorem
   that every member has small V: a deterministic contradiction from Minkowski inside Lambda_U needs R_U > 4^{dim U}; otherwise one would need
   the O-lattice L itself to beat Minkowski-Hlawka exponentially.

Strengthening of Lemma A' (referee, proved): for ANY lattice L with an isometry g of prime order q (no unimodularity), U = fix(g):
L^# cap U = (pi_U L)^# and [pi_U L : L cap U] = q^c, c <= min(a,b) (same glue argument), so covol(L cap U) covol(L^# cap U) = q^c and
  R_U(L) R_U(L^#) = ((m/a)^a q^{-c})^2 >= ((m/a)^a q^{-min(a,b)})^2  (> 1 if a != b).
Since R_U is the same for L and LJ (same tensor family), Cor. A'' holds for every G-lattice whose dual is of the form LJ (in particular modular
lattices), and for an arbitrary G-lattice at least one of L, L^# is unstable. This reduces the scope caveat of Cor. A'' to G-lattices whose dual
is in a different family.
