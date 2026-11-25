# R5_leech13 - explicit integer verification of the p = 13 (Leech / SL_2(13)) computation

Wall clock: 07:23 - 07:42 UTC (25 Nov 2026). Code: stepA.py .. stepE.py (+ leechmem.py); outputs: stepA_out.txt .. stepE_out.txt.
Total CPU time of all scripts: under 1 minute.

VERDICT: the orchestrator's p = 13 computation is CONFIRMED in every point that was checked, and the
"borderline" statement is upgraded to a computer-assisted theorem (Theorem 3 below).

## Conventions
Leech lattice L in coordinates scaled by sqrt8 (integer vectors; minimal vectors have scaled norm 32);
all norms / determinants quoted below are in the unimodular normalisation (min 4, det 1).
Exactness: every group element has matrix entries in (1/8)Z, |entries| <= 1, so float64 products are exact; all
membership tests, Gram matrices, determinants (Bareiss), elementary divisors (integer row/column reduction) are done in
exact integer arithmetic. Floating point is used only (a) to propose candidates that are then verified exactly
(matrix inverses / solves, rounded to (1/8)Z and checked by exact identities), (b) numerical ranks (dimension counts
4, 8, 72), (c) the circulant eigenvalues in the PSD enumeration, (d) the "24 minimal vectors per line" count (tolerance 1e-8).

## Step 1-2 (stepA.py) [proved (computer)]
* Extended Golay code = extended QR code mod 23 (weight distribution 1, 759, 2576, 759, 1 checked).
* 196560 minimal vectors: 1104 of shape (+-4,+-4,0^22), 97152 of shape (+-2^8,0^16), 98304 of shape (-+3,+-1^23); all
  distinct, all satisfy the Leech congruences.
* Generators: x->x+1, x->2x, x->-1/x on PL(23), sign change on an octad, Conway's xi_T = (I - J/2 on each tetrad of
  the sextet of T={0,1,2,3}) * (sign change on T). Each is orthogonal and maps all 196560 minimal vectors into L.
* Random word -> M^12 =: g of order exactly 13 (8g integral, (8g)^13 = 8^13 I, g g^T = I, exact), tr g = -2,
  sum_{j=0}^{12} g^j = 0 (no eigenvalue 1).

## Step 3 (stepA.py, stepC.py) [proved (computer), exhaustive over all 196560 minimal vectors]
a_j = w.g^j w, s = 2(a_1+a_3+a_4)+2:

| s | -4 | -2 | 0 | 2 | 4 | +-6 |
|---|---|---|---|---|---|---|
| number of minimal vectors | 13104 | 45864 | 78624 | 45864 | 13104 | 0 |

= exactly the predicted counts. E s^2 = 4, E s^4 = 41.6 exactly. 126 orbit types (a_1..a_6) occur (full table in
stepA_out.txt; every type count is a multiple of 312 = |C(g)|); the PSD enumeration (numerical eigenvalues) allows 150
types with s in {0,+-2,+-4} (18,33,48,33,18) and none with s = +-6 (best minimal eigenvalue among the 180 candidates: -0.0325).

## Step 4 structure (stepB.py, stepC.py, stepD.py, stepE.py)
1. [proved (computer)] L is free of rank 2 over Z[g] = Z[zeta_13] with a basis b_1,b_2 of minimal vectors (det of the
   24x24 matrix (g^j b_i) = 2^36 exactly). Number of Z[g]-linear isometries = |C_{Co_0}(g)| = 312; the elements with
   h^24 = 1 form a group Q of order 24 with order statistics {1:1, 2:1, 3:8, 4:6, 6:8} (= SL_2(3)), Q-span of dimension 4
   (Hurwitz units 1,i,j,omega). C(g) acts freely on minimal vectors: 630 orbits of size 312.
2. [proved (computer)] eta = sum (j/13) g^j: eta^2 = 13, eta symmetric; theta = (1+eta)/2 and O_0 := Z<1,i,j,omega> (x) Z[theta]
   preserve L. Exactly 4 of the 256 classes (1/2) sum c_k e_k mod O_0 preserve L, so
   O(L) := {x in Delta : xL in L} has index 4 over O_0, i.e. O(L) is a MAXIMAL order of Delta = H (x) Q(sqrt13)
   (every order containing O_0 lies in the trd-dual of O_0, which is inside (1/2)O_0; discrd(O_0) = 2 O_F of norm 4).
3. [proved] For w in L put n = |w|^2, t = w.eta w (then t = n mod 4), nu(w) = (t + n sqrt13)/(4 eps) in O_F totally positive,
   N(nu(w)) = (13 n^2 - t^2)/16 = 13 m_+ m_-/4. Since x* = xbar on Delta, |x w|^2 = Tr_{F/Q}(2 beta nu(w) nrd(x)) and
   det(O w) = N(nu(w))^4 for the maximal order O (checked: det(O_0 w) = 16 N(nu)^4 on 4900 sampled vectors of norms 4,6,8,10
   and on representatives of all minimal vectors).
4. [proved (computer), exhaustive] For ALL minimal vectors w: M_w := L cap (Delta w) = O w and
   det M_w = N(nu(w))^4 = (13 - s^2/4)^4: 6561 = 3^8 (s=+-4; 26208 vectors), 20736 = 12^4 (s=+-2; 91728), 28561 = 13^4 (s=0; 78624).
5. [proved (computer), exhaustive via pairs w_1+w_2 of minimal vectors over the 630 orbit representatives; totals
   16773120 and 398034000 reproduced exactly] Distribution of N(nu) over ALL vectors of norm 6 and norm 8:
   norm 6: N=9: 26208 (t=+-18), N=17: 550368, N=23: 2515968, N=27: 5660928, N=29: 8019648;
   norm 8: N=16: 393120 (= 48 x 8190: further lines isometric to 2E_8, det 2^16, nu ~ 4), N=27: 5660928, N=36: 23862384, N=43: 57552768, N=48: 100284912,
   N=51: 135547776, N=52: 74732112.
6. Sampled lines (stepB_out.txt): 1500 norm-6, 1500 norm-8, 600 norm-10 vectors: det M_w in
   {6561, 20736, 28561, 65536, 83521, 279841, 531441, ...} = N(nu(w_0))^4 for the generator w_0 of the line;
   minimum found 6561 = 3^8, no line with det < 3^8.

### Theorem 3 (Leech is exactly borderline over Delta_13) [proved, computer-assisted]
For every Delta-line U (8-dimensional), det(L cap U) >= 3^8 = 6561, i.e. R_U = 3^8/det <= 1; for every Delta-plane U
(16-dim) R_U = (3/2)^16/det(L cap U) <= 3^8/2^16 = 0.1001. Equality R_U = 1 holds for exactly 1638 = 3 x 546 lines:
546 lines through the minimal vectors with s = +4 (24 on each), 546 through those with s = -4, and 546 lines isometric to
sqrt3 E_8 carrying the 26208 norm-6 vectors with t = +-18 (48 on each; nu = 3 or 3 eps^-2). (The three classes correspond to
the three ideals of norm 9 of O_F: P^2, P'^2, (3).)
Proof. M_U = I w for a left ideal I of the maximal order O (item 2), so det M_U = N(a)^4 with the integral O_F-ideal
a = nu(w) nrd(I), and nu(v) in a for all v in M_U. Ideal norms in O_F = Z[(1+sqrt13)/2] are 1,3,4,9,12,13,... (2,5,7 inert).
Hermite (gamma_8 = 2, min >= 4): det M_U >= 2^8, so N(a) >= 4, and N(a) = 4 (det = 256) would force min M_U = 4, i.e. U would be
a line through a minimal vector of L with determinant 256 - excluded by item 4 (all such lines have det 6561, 20736 or 28561).
Hence N(a) >= 9. (I is locally principal because O is maximal; index [Ow : Iw] = N(nrd I)^2; a is integral since a_p = nu(x_p w) O_{F,p}.)
Planes: U^perp is a line and det(L cap U) = det(L cap U^perp) (L unimodular). Equality cases: N(a) = 9 means a in
{P^2, P'^2, (3)} with totally positive generator unique up to eps^{2Z} (h_F = 1, N(eps) = -1), and every borderline line
contains a vector with nu exactly (7-sqrt13)/2, (7+sqrt13)/2 or 3 [uses class number one of O, see below, to have a
generator w_0 with nu(w_0) generating a; then scale by eps^k]; these are the minimal vectors with s = +-4 and the norm-6
vectors with t = 18; counts from items 3-5 and the unit count 24 = |O^1|. QED.
(The count 1638 uses h(O) = 1 [literature: Eichler mass formula, mass = |zeta_F(-1)|/2 = 1/12 = 1/[O^x : O_F^x] since
O^x contains O_F^x SL_2(3)]; the inequality det >= 3^8 does not.)
Independent proof of N(nu(w)) >= 9 for all w != 0 (orchestrator's unit argument, made precise): eps = (3+eta)/2 is in End(L),
|eps^k w|^2 = eps^{2k} m_+ + eps^{-2k} m_- >= 4 for all k; choosing k with eps^{4k} m_+/m_- in [eps^-2, eps^2) gives
4 <= sqrt(13 m_+ m_-), so N(nu) >= 4 with equality iff some eps^k w is a minimal vector with s = +-6; none exists (Step 3).

## Delta-linear automorphism group (stepB.py, stepC.py, stepD.py) [proved (computer)]
Aut_Delta(L) := {isometries of L commuting with i, j, eta}. Exhaustive search (images of a Delta-basis of three s=4 minimal
vectors, matched by the Delta-valued Hermitian products; orbit-stabiliser; in every failed search no triple with matching
Hermitian data exists, so failures are exact) run twice with different bases: stabiliser of a minimal vector trivial,
orbit 2184 (the 13104 vectors with s = 4 split into 6 orbits), so |Aut_Delta(L)| = 2184.
Element orders {1:1, 2:1, 3:182, 4:182, 6:182, 7:468, 12:364, 13:168, 14:468, 26:168} = those of SL_2(13); a group of
order 2184 with these statistics is SL_2(13) (unique involution, no normal Sylow subgroup or normal C_4/Q_8, so G/{+-1} is
simple of order 1092 = PSL_2(13), extension non-split). Traces: order 13: -2, order 26: 2, order 7: -4, order 14: 4,
orders 3,4,6,12: 0. dim_Q End_G(W) = 8 (so End_G(W) = Delta), dim_Q span(G) = 72 = 9*8 (A = M_3(Delta^op)),
sum_g (tr g/24)^2 = 2184/72 = 30.33, gain g = |G|/(2 dim A) = 15.17. An order-3 element has trace 0 and fixed
lattice of rank 8 and determinant exactly 3^8. The 91 subgroups of order 3 give 91 distinct such fixed lines (elements of order
2, 7, 13 have no eigenvalue 1, so point stabilisers are 3-groups), all of type (3) (stabilisers of minimal vectors are trivial):
91 of the 546 lines of type (3) are fixed lines of order-3 elements, the other 455 are not.

## What was NOT done / caveats
* Class number one of O is taken from the mass formula (memory), not recomputed; it is used only for the exact count 1638
  and for "M_U = O w_0" in general (verified directly on all lines through minimal vectors and on ~3600 sampled lines).
* Numerical ranks (4, 8, 72) and PSD eigenvalues are floating point (well separated).
* Nothing here addresses whether an equivariant Klartag argument actually works in the borderline case, nor p > 13.

## Addendum (stepF.py, stepF_out.txt) [proved (computer)]
No element h != 1 of G = Aut_Delta(L) = SL_2(13) fixes a minimal vector: G acts freely on the 196560 minimal vectors
(90 regular orbits). For an order-3 element c the trace map 1+c+c^2 sends minimal vectors to 0 (720), to norm-6 vectors
with t = +-18 (N(nu) = 9; 5184 each) or t = +-6 (N(nu) = 27), and to longer vectors: the fixed line is of type (3),
i.e. sqrt3 E_8 (240 norm-6 vectors = 48 with N(nu)=9 + 192 with N(nu)=27).
