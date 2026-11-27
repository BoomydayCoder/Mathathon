# R21 — Is "(p+1)/2 prime" really needed?  Exploration of weaker hypotheses for the symmetric tensor-family strategy

Round 21 (after the timed record; done by the orchestrator on request, 2026-11-27).  Nothing here was refereed.
Labels: **PROVED** (argument written here, by me, unchecked by anyone else), **CITED**, **COMPUTED** (exact computer
algebra), **NUMERICAL**, **HEURISTIC**.

## 0. Verdict

* Within the strategy "finite group G acting on W, random lattices Lambda_W (x)_O L, first-moment/Klartag argument",
  the arithmetic hypothesis **cannot be weakened inside the SL_2 family**: for SL_2(q) the only faithful rational
  characters passing the necessary conditions with gain > 2.4 are the Weil characters with q = p prime, p = 1 mod 4,
  (p+1)/2 prime (Section 3, PROVED + COMPUTED for q < 2000).  If (p+1)/2 is composite the gain drops from about p to
  at most about 2, whatever (mixed) representation is used (Section 4).  So almost-primes (Chen-type results) cannot
  replace primes.
* **No other family was found.**  A scan of the whole GAP character table library (2313 tables: ATLAS groups, covers,
  maximal subgroups, ...) finds 181 pairs (G, chi) passing the prime-order test; after the cyclic test only 25 have gain
  >= 1.9, and the only one above the icosian value 7.5 is SL_2(13) (gain 15.17) (Section 5, COMPUTED).  This is evidence,
  not a classification.
* The hypothesis is "strong" only in the sense of being unproved: such primes are plentiful (34 below 2000:
  5, 13, 37, 61, 73, 157, 193, ...; the Hardy-Littlewood heuristic predicts on the order of x/(log x)^2 up to x).  The real bottleneck is
  stability (BS), not the supply of primes.
* What can be weakened is the *use* of the hypothesis (Section 6): a constant-factor improvement needs one stable prime,
  not infinitely many; an order improvement needs infinitely many stable data from *any* admissible source, and none
  other than SL_2(p) is known.

## 1. Two general facts

Data: finite G, rational representation V (not necessarily irreducible), n = dim V, tau = its character,
rho = tau/n.  A = Q-span of G in End(V).

**Lemma C (gain formula; PROVED).**  For irreducible W: |G|/(2 dim_Q A) = (1/2) sum_{g in G} rho(g)^2.
(sum tau^2 = |G| dim_Q End_G(W) = |G| d, n = r d, dim_Q A = r^2 d.)  For general V the number of orbits that the
G-invariant quadratic forms can absorb is dim Sym^2(V)^G ~ (1/2|G|) sum tau^2, so the same expression should be the gain
(HEURISTIC: the template theorem of Round 14 is proved for isotypic V only).

**Prime-order rigidity (Round 5, PROVED there).**  If some G-lattice in V^k is to satisfy R_U <= 1 for U = Fix(x)
and its complement, x of prime order q, then x is *free* (no eigenvalue 1; tau(x) = -n/(q-1)) or *regular*
(V restricted to <x> is a multiple of the regular representation; tau(x) = 0).

So the gain is carried by elements with non-zero trace, and prime-order elements with fixed points contribute nothing.

## 2. Structure lemmas (all PROVED here, elementary)

**Lemma A.**  Let E = C_l x C_l act on V with every subgroup of order l free or regular.  Then at most one subgroup
is free, and if one is free the other l are regular.
*Proof.*  m_0 = multiplicity of the trivial character, m_i = total multiplicity of the characters with kernel C_i
(i = 0..l).  dim Fix(C_i) = m_0 + m_i is 0 or n/l, and n = m_0 + sum m_i.  If C_j is free, m_0 = m_j = 0, so the
remaining l numbers m_i <= n/l sum to n: all equal n/l.  QED

**Lemma B (free elements are Sylow-central).**  If c has prime order l and acts freely, then c is central in every
l-subgroup containing it; each Sylow l-subgroup contains at most one free subgroup of order l.
*Proof.*  Let P be an l-group containing c with c not central.  Choose x in N_P(C_P(c)) outside C_P(c).  Then c^x
lies in C_P(c), c^x is not in <c> (an l-element normalising <c> centralises it), so <c, c^x> = C_l x C_l has two free
subgroups, contradicting Lemma A.  Two central free subgroups commute, so they coincide.  QED
(COMPUTED check: in all 181 admissible pairs of Section 5 the free classes have class length prime to l.)

**Lemma D.**  If every prime-order power of g (order m) is free, all eigenvalues of g are primitive m-th roots and
rho(g) = mu(m)/phi(m).  In particular free elements of non-squarefree order have trace 0.

**Order-4 rigidity (PROVED for unimodular lattices).**  If g has order 4 and g^2 is a regular involution, then
tau(g) = 0.  (a = dim Fix g, b_2 = dim ker(g+1), a + b_2 = n/2.  The glue of Fix(g) is at most
|det(1-g | complement)| = 2^{b_2 + n/4}; R_U <= 1 needs it to be >= (n/a)^a, and a >= n/4.  These hold only for
a = n/4.  If a = 0 the subspace ker(g+1) has glue <= 2^{n/4} < 2^{n/2}.)

**Claim E (HEURISTIC in general; proved above for m = 4 and checked by hand for m = 6, 15, l^2).**  Under (BS), an
element with a regular prime-order power has trace 0; hence
gain = (1/2) sum over "totally free" g of 1/phi(ord g)^2  (squarefree orders only).
The scan flag `strong` tests exactly this pattern.

Consequence (PROVED from Lemmas B, C, D, given Claim E): a large gain needs *many conjugate free cyclic subgroups*:
the contribution of a prime l is (1/2) [G : N_G(C_l)]/(l-1), with C_l central in a Sylow l-subgroup.

## 3. The family SL_2(q), q odd (PROVED from the generic character table [CITED]; COMPUTED for all q < 2000)

Rational characters: DS(t) = sum of discrete series X_theta, theta of order t | q+1; HD = eta_1 + eta_2 (half discrete
series = odd Weil characters when q = 1 mod 4); PS(t), HP (principal series), St.  Prime-order elements: -1, unipotent,
split torus, non-split torus.  Result (`sl2_family.py`, `out_sl2_family.txt`):

* PS, HP, St are never admissible (unipotent or split-torus elements have a surplus of fixed vectors).
* q = p^f with f > 1: nothing (unipotent elements).
* DS(t), t > 2: admissible only if q + 1 is a power of 2 (Mersenne primes); gain about 2.
* HD: admissible iff q = p = 1 mod 4 and (p+1)/2 is prime; gain p(p+1)/(p-1).
  Reason: (psi + psi')|T' = twice (all odd characters of T' = C_{p+1} except the one of order 2).  An element of odd
  prime order l | M := (p+1)/2 fixes 2(M/l - 1) dimensions; regular would be 2(M/l - 1/l); free iff M = l.
  "Regular minus one character" is free only on a group of prime order.

Of the 35 pairs with gain >= 2.5 for q < 2000, 34 are of the last type (p = 5, 13, 37, 61, 73, 157, 193, 277, 313, 397, ...);
the remaining one is SL_2(3) (gain 3, fixed-point-free).

## 4. Composite M = (p+1)/2

**4.1 Same representation W_p (Round 5).**  Elements of order l | M have R_U >= ((M-1)/(M-l))^{4(M/l-1)}, a bounded
constant (<= e^4) but > 1; the first-moment term is proportional to R_U^{k/2}, exponential in the rank k, and the
method needs k >> r = (p-1)/4.  A bounded violation is therefore fatal; it cannot be absorbed.

**4.2 Mixed representations V = sum_theta w_theta X_theta (odd discrete series; PROVED/COMPUTED).**
Write -tau_V on T' minus centre as sum_theta c_theta theta, S = sum c_theta.  Then the gain carried by the p(p-1)/2
conjugates of T' is exactly

    G_T' = (2pM/(p-1)) ( sum c_theta^2 / S^2  -  1/M )   ~   p / K_eff  - 2,     K_eff = S^2 / sum c^2

(W_p: K_eff = 1, G_T' = p).  Freeness at l | M is impossible for composite M, and regularity at l is equivalent to:
the push-forward of c to the characters of C_l is uniform.  Rationality forces c to be constant on Galois classes.
* M = l1 l2 (l1 < l2): the admissible cone has two extreme rays and  max G_T' = 2p(l1-1)/((p-1)(l2-1)) <= 2
  (`composite_M.py`: exact, checked for p = 29, 41, 109, 173, 281, 409, 1009).
* More prime factors leave more room under the prime-order conditions alone (p = 2729, M = 3*5*7*13: 17.5 instead of
  2731), but:
* adding the two necessary conditions for the cyclic group C_M itself (fixed space a = 0 or a >= n/M; glue bound)
  leaves only the uniform weight, i.e. G_T' = 0, in the two-prime case (`composite_M_cyclic.py`, NUMERICAL grid + the
  monotonicity argument in the comments).  Total gain is then <= 1 + (p+1)/(p-1), about 2.
* Adding principal series or the Steinberg character puts every character of T' into V|T' and kills G_T' completely.

So for composite (p+1)/2 no representation of SL_2(p) gives more than a bounded gain.  Unconditional results of
Chen type ("p prime, (p+1)/2 a product of at most two primes") are useless here.

## 5. Other groups

**5.1 Three mechanisms (HEURISTIC organisation, each statement below is elementary).**
(I) *Regular-minus-one tori.*  A cyclic TI subgroup T on which a non-exceptional character is constant -1:
chi|T = (regular) - (one character).  Free iff |T/Z| is prime, and then chi(1) = |T/Z| - 1.  Examples: SL_2(p) with
(p+1)/2 prime; 2^{1+2a}_-.C_q with q = 2^a + 1 a Fermat prime (gain q); Weil characters of Sp_{2n}(p), SU_n(p)
(there other elements fail).  The primality is forced by the mechanism, and the group brings its own prime p: two
simultaneous prime conditions.
(II) *fpf-like tori* (chi|T consists of faithful characters, chi(1) = [N(T):T]): Jordan-bounded; SL_2(5) is typical;
gain O(log log).
(III) *Clifford spreading*: G = E.H with W|E irreducible (E extraspecial): gain = |H| x gain(E).  Then
End_G(W) is contained in End_E(W), which is Q, Q(i), the rational quaternions or Q(zeta_l): a small algebra.

**5.2 The discriminant obstruction (Round 3 no-go (1),(2); recalled).**  For r >= 2 every principal line must satisfy
det(Lambda cap Delta w) >= r^d.  For Delta = rational quaternions (maximal order O_H) this reads |w|^2 >= sqrt(2) r for
all w, impossible for large r and already false for E_8 (r = 2): a root gives R_U >= 4.  [PROVED given the standard
Hurwitz structure of E_8.]  So mechanism (III) with E a 2-group fails at lines, whatever H is.

**5.3 Scan of the GAP character table library (COMPUTED; `scan.g`, `scan3.g`, `scan_out.txt`, `scan3_lib.txt`).**
2320 non-duplicate tables, 7 skipped (more than 1200 classes: 2.O8+(7), O12-(3), U6(4), 2_1.O12+(3) and three large
2-local groups).  For every faithful irreducible chi: tau = sum of Galois conjugates, gain = |G|/(2 c chi(1)^2).
Tests: prime-order classes free or regular; `strong` = pattern of Claim E on all classes; free classes Sylow-central.

* 181 pairs pass the prime-order test with gain >= 1; all have Sylow-central free classes; 117 pass `strong`.
* Passing `strong`, gain >= 1.9 (25 pairs): 2.L2(13) and 2.L2(13).2 (15.17); 2.A5 and five relatives (7.5);
  NDG(2.HS,Q8) (3.75); 3x2.L3(2), (SL2(3)x7:3).2 (3.5); 2.(13:6xA4).2 (3.25); SL2(3), 2.S4 types (3);
  2.L3(2) = SL2(7) (2.33); 2.L2(31) (2.07).
* Failing `strong` (so excluded by order-4 rigidity) and also excluded by 5.2:
  2^{1+4}_-.A5 and 4 o 2^{1+4}.A5 (gain 60, on E_8), 2.2^4.2^4.S3 (24, real), 2.U4(2)N2 (12), 2(L2(7)x4).2 (9.3),
  various 2-groups (<= 8).
* Not admissible at all (prime-order test): 2.L2(17), 2.L2(25), 2.L2(29), 2.J2, 2.G2(4), 6.Suz, 2.Co1, 2.A7, 6.A7,
  6.L3(4), 6_1.U4(3), Sz(8), Sz(32), U3(q), S4(q), ... (every other table in the library).
* 2.L2(37), 2.L2(61) are not in the library; they are covered by Section 3.

Limits: irreducible faithful characters only (non-isotypic V only in Section 4); subgroups of the listed groups are
covered only as far as they have their own table; no classification theorem.

## 6. What can honestly be weakened

1. **Constant-factor improvement needs one prime, not infinitely many.**  p = 13 is done (2.79 n^2, 24 | n).  One more
   stable prime (p = 61 would give about 11.6 n^2) is a finite problem.
2. **Order improvement** needs gain -> infinity, hence infinitely many admissible stable data; by Sections 3-5 the only
   known source is SL_2(p), p and (p+1)/2 prime.  Products of data with disjoint free primes multiply gains
   (2 g_1 g_2) but still need the factors.
3. **Hybrid "Leech x cyclotomic"** (would need only p = 13): G = SL_2(13) x C_m (m odd, prime to 7*13) acting on Leech (x) Z[zeta_m] has gain
   15.17 m/phi(m).  But stability is not inherited: for m = 3 (C_3 inside the commutant order, Delta' = Q(sqrt 13, zeta_3),
   r' = 6) the line through a minimal vector with s = +-4 has det(O_K w) = 729 < 1296, so R_U >= 16/9 > 1.
   [PROVED given the Round-5 facts N(nu) = 9 for those vectors and O containing O_F[zeta_3].]  For m prime to 6 nothing
   is known; and a theorem uniform in m would be needed.  At best this changes the constant of n^2 log log n.
4. The congruence p = 5 mod 8 used in Rounds 13, 16, 20 is a limitation of the computations/model, not of the
   strategy: admissibility holds for all p = 1 mod 4 with (p+1)/2 prime (73, 193, 313, ... are untested).

## 7. Not done

No proof of Claim E in general; no proof that mechanisms (I)-(III) exhaust all possibilities; no scan of non-isotypic
representations beyond SL_2(p); the 7 skipped tables; (BS) itself was not tested for any new group because none
survived the necessary conditions.
