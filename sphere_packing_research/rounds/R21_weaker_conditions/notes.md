# R21 — Is "(p+1)/2 prime" really needed?  Exploration of weaker hypotheses for the symmetric tensor-family strategy

Round 21 (after the timed record; done by the orchestrator on request, 2026-11-27).  Nothing here was refereed.
**Revised on 2026-11-27 after a critic pass** (the critic compared the report text with this file and the raw outputs):
the scope "unimodular G-lattices" is now stated explicitly; the former heuristic "Claim E" is replaced by two proved
lemmas (prime-power regularity, fibrewise regularity), which also replace the LP/grid argument for composite (p+1)/2;
the scan was re-run with proved conditions only (`scan4.g`: same 117 survivors); several overstatements were corrected.
Labels: **PROVED** (argument written here, by me, unchecked by anyone else), **CITED**, **COMPUTED** (exact computer
algebra), **NUMERICAL**, **HEURISTIC**.

## 0. Verdict

Scope.  All necessary conditions below are for a **unimodular** G-lattice Lambda in the representation V that satisfies
R_U <= 1 for the subspaces U = ker f(g) (g in G, f a polynomial) and their complements.  For W_p this loses nothing
(all G-lattices Lambda I of the family have the same R_U, Round 20 Theorem 1(c)).  For a non-unimodular lattice in
some other representation only the entropy inequality of Round 3 is available for odd primes; that case was not scanned.

* **Inside the SL_2 family the arithmetic hypothesis cannot be weakened** (Sections 3-4, PROVED by me, not refereed;
  exact computation for q < 2000).  For SL_2(q) the only faithful rational characters passing the prime-order test with
  gain > 3 are the Weil characters with q = p prime, p = 1 mod 4, (p+1)/2 prime.  If M = (p+1)/2 is composite and not
  divisible by the squares of two distinct primes (so in particular if M is a product of at most two primes), then for
  *every* faithful rational representation the character vanishes on the non-split torus outside +-1, and for sums of
  odd discrete series the gain is at most 1 + (p+1)/(p-1), about 2.  So almost-primes (Chen-type results) cannot replace primes.
* **No other family was found** (Section 5, COMPUTED): in the GAP character table library (2313 tables) 181 pairs
  (G, chi) pass the prime-order test, 117 pass all proved tests, and the only one with gain above the icosian value 7.5
  is SL_2(13) (15.17).  This is evidence, not a classification.
* The hypothesis is "strong" only in the sense of being unproved: such primes are plentiful (34 below 2000;
  the Hardy-Littlewood heuristic predicts on the order of x/(log x)^2 up to x).  The real bottleneck is (BS).
* What can be weakened is the *use* of the hypothesis (Section 6).

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

## 2. Structure lemmas (PROVED here, by me; not refereed)

Standing hypothesis (U): Lambda unimodular G-lattice in V, n = dim V, and R_U <= 1 for U = ker f(g) and complements.
For a subspace U write Lambda_U = Lambda cap U and glue(U) = [Lambda : Lambda_U + Lambda_{U^perp}] = det(Lambda_U).

**L1 (prime-order rigidity, Round 5).**  x of prime order l is free, or regular with full glue: glue(Fix x) = l^{n/l},
pi_U(Lambda) = (1/l) Lambda_U, and Lambda_{U^perp} = (1-zeta) pi_{U^perp}(Lambda) as Z[zeta_l]-modules (x acting as zeta).

**L2.**  Let E = C_l x C_l act on V with every subgroup of order l free or regular.  Then at most one subgroup
is free, and if one is free the other l are regular.
*Proof.*  m_0 = multiplicity of the trivial character, m_i = total multiplicity of the characters with kernel C_i
(i = 0..l).  dim Fix(C_i) = m_0 + m_i is 0 or n/l, and n = m_0 + sum m_i.  If C_j is free, m_0 = m_j = 0, so the
remaining l numbers m_i <= n/l sum to n: all equal n/l.  QED

**L3 (free elements are Sylow-central).**  If c has prime order l and acts freely, then c is central in every
l-subgroup containing it; each Sylow l-subgroup contains at most one free subgroup of order l.
*Proof.*  Let P be an l-group containing c with c not central.  Choose x in N_P(C_P(c)) outside C_P(c).  Then c^x
lies in C_P(c), c^x is not in <c> (an l-element normalising <c> centralises it), so <c, c^x> = C_l x C_l has two free
subgroups, contradicting L2.  Two central free subgroups commute, so they coincide.  QED
(COMPUTED check: in all 181 admissible pairs of Section 5 the free classes have class length prime to l.)

**L4 (prime-power regularity).**  Let y have order l^i and let its power x of order l be regular.  Then V restricted
to <y> is a multiple of the regular representation; in particular tau(y) = 0.  (i = 2, l = 2: an element of order 4
whose square is a regular involution has trace 0.)
*Proof.*  Induction on i; i = 1 is the hypothesis.  Let z = y^l.  By induction V|<z> is regular.  A character lambda
of <y> of order l^t, t >= 2, has l Galois-conjugate extensions of lambda|<z> (lambda^{1 + l^{t-1} k}), so each has
multiplicity n/l^i.  The remaining characters (order 1 and l) have multiplicities a = dim Fix(y) and m_1 with
a + (l-1) m_1 = n/l^{i-1}.
(i) Suppose a > 0, U = Fix(y).  glue(U) is killed by 1-y on the complement side, so
glue(U) <= |det(1-y | U^perp)| = l^{m_1 + (i-1) n/l^i}, and R_U <= 1 needs glue(U) >= (n/a)^a.  The function
D(a) = a log(n/a) - (m_1 + (i-1) n/l^i) log l is concave on [n/l^i, n/l^{i-1}], D(n/l^i) = 0, D'(n/l^i) > 0 and
D(n/l^{i-1}) = (i-1)(l-1)(n/l^i) log l > 0; so D > 0 for a > n/l^i, a violation.  Also glue(U) <= l^{ia} gives a >= n/l^i.
Hence a = n/l^i and then m_1 = n/l^i: regular.
(ii) Suppose a = 0.  On U' = Fix(z) (dimension n/l^{i-1}) y has only eigenvalues of order l, so Phi_l(y) kills glue(U'):
glue(U') <= |det(Phi_l(y) | U'^perp)| = l^{(i-1)(l-1) n/l^i}  (prod over primitive l^t-th roots of Phi_l = l^{l-1}).
But z acts regularly with group <z> of order l^{i-1}, and R_{U'} <= 1 needs glue(U') >= (l^{i-1})^{n/l^{i-1}}
= l^{(i-1) l n/l^i}.  Contradiction.  QED

**L5 (fibrewise regularity).**  Let x be regular of prime order l and let h commute with x, the order of h prime to l.
Then all eigenspaces V_{zeta^j} of x (j = 0..l-1) are isomorphic as <h>-modules: V|<x> x <h> = reg x V_1.
In particular tau(x^a h^b) = 0 whenever x^a != 1.
*Proof.*  By L1 the glue group A = pi_U(Lambda)/Lambda_U is (1/l)Lambda_U / Lambda_U = Lambda_U / l Lambda_U and also
A = M'/(1-zeta)M' with M' = pi_{U^perp}(Lambda), a Z[zeta_l]-lattice of rank n/l in U^perp.  Both identifications are
h-equivariant.  Since |h| is prime to l, an F_l[h]-module is determined by its Brauer character, which for the
reduction of a lattice is the ordinary character: so the character of h on U = V_1 equals the character of h on the
zeta-eigenspace V_zeta of x.  A Galois automorphism fixing the |h|-th roots of unity and sending zeta to zeta^j gives
V_{zeta^j} = V_zeta (V is rational).  QED

**Corollary (which elements can carry gain).**  tau(g) = 0 whenever (a) some prime l divides ord(g) exactly once and
the power of g of order l is regular, or (b) ord(g) is a prime power and the power of order l is regular.
If every prime-order power of g is free, all eigenvalues are primitive and rho(g) = mu(m)/phi(m) (0 unless m is
squarefree).  Not covered: elements for which every regular prime divides the order at least twice and the order is
not a prime power.  (Expected to be 0 as well, by combining L4 and L5 inductively on the rescaled unimodular lattice
(Fix x, form/l), which inherits (BS); not written out.  In the library scan no class needed it.)

Consequence: the gain is carried by "totally free" elements, and by L3 a large gain needs *many conjugate free
cyclic subgroups*: the contribution of a prime l is (1/2) [G : N_G(C_l)]/(l-1), C_l central in a Sylow l-subgroup.

## 3. The family SL_2(q), q odd (PROVED from the generic character table [CITED]; COMPUTED for all q < 2000)

Rational characters: DS(t) = sum of discrete series X_theta, theta of order t | q+1; HD = eta_1 + eta_2 (half discrete
series = odd Weil characters when q = 1 mod 4); PS(t), HP (principal series), St.  Prime-order elements: -1, unipotent,
split torus, non-split torus.  Result (`sl2_family.py`, `out_sl2_family.txt`):

* PS, HP, St are never admissible (unipotent or split-torus elements have a surplus of fixed vectors).
* q = p^f with f > 1: nothing (unipotent elements).  (`sl2_family.py` does not test HD/HP for prime powers; by hand:
  for f odd tau(u) = -1 gives fixed dimension (q-p)/p, neither 0 nor (q-1)/p; for q a square eta_1, eta_2 are rational
  and one of the two unipotent classes has fixed dimension different from 0 and (q-1)/(2p).)
* DS(t), t > 2: admissible only if q + 1 is a power of 2 (q = 3, 7, 31, 127, ...); gain 3, 2.33, 2.07, ... -> 2.
* HD: admissible iff q = p = 1 mod 4 and (p+1)/2 is prime; gain p(p+1)/(p-1).
  Reason: (psi + psi')|T' = twice (all odd characters of T' = C_{p+1} except the one of order 2).  An element of odd
  prime order l | M := (p+1)/2 fixes 2(M/l - 1) dimensions; regular would be 2(M/l - 1/l); free iff M = l.
  "Regular minus one character" is free only on a group of prime order.

Of the 35 pairs with gain >= 2.5 for q < 2000, 34 are of the last type (p = 5, 13, 37, 61, 73, 157, 193, 277, 313, 397, ...);
the remaining one is SL_2(3) (gain 3, fixed-point-free).

## 4. Composite M = (p+1)/2

**4.1 Same representation W_p (Round 5).**  Elements of order l | M have R_U >= ((M-1)/(M-l))^{4(M/l-1)} > 1.
This lower bound is itself bounded (it is <= e^4), but the first-moment term is proportional to R_U^{k/2}, exponential
in the rank k, and the method needs k >> r = (p-1)/4.  A violation by a bounded factor is therefore still fatal.

**4.2 Any faithful rational representation V of SL_2(p), p = 1 mod 4 (PROVED from L1, L4, L5; not refereed).**
Let M be composite and b_l in T' of prime order l | M.  Every irreducible constituent of V restricted to T' contains a
character trivial on b_l: X_theta|T' = twice (all characters of the parity of theta) minus theta, theta^-1, and
there are M/l >= 3 characters of that parity trivial on b_l; principal series, Steinberg and the half characters
likewise.  So b_l is not free, hence regular (L1).  Now let b in T' be non-central, ord(b) = 2^e d with d | M, d > 1.
If some prime l divides d exactly once, tau_V(b) = 0 by L5.  Otherwise, if M is not divisible by the squares of two
distinct primes, d = l^i is a prime power and b = +-y: tau_V(y) = 0 by L4, and tau_V(-y) = 0 because -1 is either
the scalar -1 or a regular involution (then L5 applies with x = -1).  Hence

    tau_V = 0 on T' minus {+-1}        (M composite, not divisible by the squares of two distinct primes),

i.e. the non-split tori, which carry the whole gain ~p of W_p when M is prime, contribute nothing, whatever
representation is used.  For V a sum of odd discrete series characters (X_theta and psi + psi'), tau_V also vanishes
on the split torus outside +-1 (character table), so only the centre and the unipotent classes are left:

    gain(V)  <=  1 + (p+1)/(p-1)  ~  2.

(For representations with principal-series constituents the unipotent elements have fixed vectors, hence must be
regular; I did not bound the split-torus contribution in that case beyond what L4/L5 give.)
In particular for M = l1*l2 or M = l^2 ("(p+1)/2 a product of at most two primes") the gain drops from about p to
about 2.  **Chen-type results cannot replace the primality of (p+1)/2.**

**4.3 A weaker, independent computation (kept for the record).**  For V = sum_theta w_theta X_theta (odd discrete
series), write -tau_V on T' minus centre as sum c_theta theta, S = sum c.  The gain carried by T' is exactly
G_T' = (2pM/(p-1)) (sum c^2 / S^2 - 1/M)  (W_p: G_T' = p).  Using *only* the prime-order conditions (push-forward of c
to the characters of C_l uniform), the maximum for M = l1*l2 is 2p(l1-1)/((p-1)(l2-1)) (`composite_M.py`, exact vertex
enumeration, p = 29, 41, 109, 173, 281, 409, 1009); with more prime factors the prime-order conditions alone leave
more room (p = 389: 2.76; p = 2729, M = 3*5*7*13: 17.5).  Section 4.2 shows that the full conditions give 0 in all
these cases.  (`composite_M_cyclic.py` is an earlier numerical check of the two-prime case on a grid; it finds no
feasible point except the uniform one, which lies on the grid only for l1 = 5.  It is superseded by 4.2.)

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
(That 2^{1+2a}_-.C_q fails for this reason when a >= 2 is my inference from 5.2; it was not computed separately.)

**5.2 The discriminant obstruction (Round 3 no-go (1),(2); recalled).**  For r >= 2 every principal line must satisfy
det(Lambda cap Delta w) >= r^d.  For Delta = rational quaternions (maximal order O_H) this reads |w|^2 >= sqrt(2) r for
all w, impossible for large r and already false for E_8 (r = 2): a root gives R_U >= 4.  [PROVED given the standard
Hurwitz structure of E_8.]  So mechanism (III) with E a 2-group fails at lines, whatever H is.

**5.3 Scan of the GAP character table library (COMPUTED; `scan.g`, `scan3.g`, `scan4.g` and their outputs).**
2320 non-duplicate tables, 7 skipped (more than 1200 classes: 2.O8+(7), O12-(3), U6(4), 2_1.O12+(3) and three large
2-local groups).  For every faithful irreducible chi with gain = |G|/(2 c chi(1)^2) >= 1 (pairs with gain < 1 are not
examined): tau = sum of Galois conjugates.  Tests in the final run `scan4.g`, all PROVED necessary conditions under (U):
(R1) prime-order classes free or regular (L1); (R2) tau(g) = 0 if a prime l divides ord(g) exactly once and the power
of order l is regular (L5); (R3) tau(g) = 0 if ord(g) is a prime power with regular power of order l (L4).

* 181 pairs pass (R1); all have Sylow-central free classes (check of L3); 117 pass (R1)-(R3).  The same 117 pass the
  older heuristic flag `strong` of `scan3.g`, and no surviving pair has a class outside the reach of (R2)/(R3).
* Passing all tests, gain >= 1.9 (25 pairs): 2.L2(13) and 2.L2(13).2 (15.17); 2.A5 and five relatives (7.5);
  NDG(2.HS,Q8) (3.75); 3x2.L3(2), (SL2(3)x7:3).2 (3.5); 2.(13:6xA4).2 (3.25); SL2(3), 2.S4 types (3);
  2.L3(2) = SL2(7) (2.33); 2.L2(31) (2.07).
* Passing (R1) but failing (R2)/(R3) (most at classes of order 4 with regular square and non-zero trace), and also
  excluded by 5.2: 2^{1+4}_-.A5 and 4 o 2^{1+4}.A5 (gain 60, on E_8), 2.2^4.2^4.S3 (24, real), 2.U4(2)N2 (12),
  2(L2(7)x4).2 (9.3), various 2-groups (<= 8).
* Failing (R1) (among pairs with gain >= 1): 2.L2(17), 2.L2(25), 2.L2(29), 2.J2, 2.G2(4), 6.Suz, 2.Co1, 2.A7, 6.A7,
  6.L3(4), 6_1.U4(3), Sz(8), Sz(32), U3(q), S4(q), ... (every other table in the library).
* 2.L2(37), 2.L2(61) are not in the library; they are covered by Section 3.

Limits: irreducible faithful characters only; unimodular lattices only (scope in Section 0); subgroups of the listed
groups are covered only as far as they have their own table; no classification theorem.

## 6. What can honestly be weakened

1. **Constant-factor improvement needs one prime, not infinitely many.**  p = 13 is done (2.79 n^2, 24 | n).  One more
   stable prime (p = 61 would give about 11.6 n^2) is a finite problem.
2. **Order improvement** needs gain -> infinity, hence infinitely many admissible stable data; by Sections 3-5 the only
   known source is SL_2(p), p and (p+1)/2 prime.  Products of data with disjoint free primes multiply gains
   (2 g_1 g_2) but still need the factors.
3. **Hybrid "Leech x cyclotomic"** (would need only p = 13): G = SL_2(13) x C_m (m odd, prime to 7*13) acting on Leech (x) Z[zeta_m] has gain
   15.17 m/phi(m).  But stability is not inherited: for m = 3 (C_3 inside the commutant order, Delta' = Q(sqrt 13, zeta_3),
   r' = 6) the line through a minimal vector with s = +-4 has det(O_K w) = 729 < 1296, so R_U >= 16/9 > 1 (for s = +-2 and s = 0 the same sublattice determinants are 1296 and 1521, which by themselves show no violation).
   [PROVED given the Round-5 facts N(nu) = 9 for those vectors and O containing O_F[zeta_3].]  For m prime to 6 nothing
   is known; and a theorem uniform in m would be needed.  At best this changes the constant of n^2 log log n.
4. The congruence p = 5 mod 8 used in Rounds 13, 16, 20 is a limitation of the computations/model, not of the
   strategy: admissibility holds for all p = 1 mod 4 with (p+1)/2 prime (73, 193, 313, ... are untested).

## 7. Not done

No proof that mechanisms (I)-(III) exhaust all possibilities; the Corollary of Section 2 is not proved for elements in
which every regular prime divides the order at least twice (not needed in the scan); non-unimodular lattices; no scan of
non-isotypic representations beyond SL_2(p); the 7 skipped tables; (BS) itself was not tested for any new group because
none survived the necessary conditions.
