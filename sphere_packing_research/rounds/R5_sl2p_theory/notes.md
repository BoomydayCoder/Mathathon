# R5_sl2p_theory - theory + literature for the SL_2(p) tensor family (p = 1 mod 4)

Start 07:23 UTC, end ~07:44. Code/outputs in this directory: `arith.py` -> `arith_out.txt`, `checks.py` -> `checks_out.txt`,
`enum17.py` -> `enum17_out.txt`.
Tags: [P] proved here (complete argument below), [PC] proved by computer (exact integer arithmetic unless stated), [C] conditional,
[H] heuristic, [N] numerical, [L] literature, [L-mem] literature quoted FROM MEMORY (not re-read), [R] refuted, [O] open.

## 0. Verdict (short)

1. **New provable obstruction [P]** (Lemma A', Cor. A''): for an even/odd UNIMODULAR SL_2(p)-lattice W_p of rank 2(p-1), an element g of odd
   prime order q | (p+1)/2 with q < (p+1)/2 has a nonzero fixed space U = W^g with R_U >= ((Nq-1)/(Nq-q))^{4(N-1)} > 1, N = (p+1)/(2q).
   So **W_p is tensor-UNSTABLE whenever (p+1)/2 is composite**: p = 17 (R_U >= (4/3)^8 = 9.99), 29 (>= 11.8), 41 (>= 12.5), 53, 89, 97,
   101, 109, 113, ... The lower-rank term T_U is then >= R_U^{k/2}/poly(n): exponentially large.
   Survivors: **p with (p+1)/2 prime: 5, 13, 37, 61, 73, 157, 193, 277, 313, 397, 421, ...** (infinitely many only conjecturally: q and 2q-1 prime).
2. For the survivors all subgroup-fixed subspaces are harmless: unipotent, order-4 and non-split-torus elements act freely; odd-order
   subgroups of the split torus give R_U = 1 EXACTLY (Lemma B [P]): the family can be at best borderline, never strictly stable (p >= 13).
3. Lines: R_U <= 1 for all Delta-lines  <=>  N_{F/Q}(nu(w)) >= ((p-1)/4)^2 for all w != 0 (maximal order, class number one; Lemma 1 [P]).
   This is NOT implied by the minimum alone for any p >= 13 (Lemma 3 [P]); it is a statement about the Hermite-Humbert (norm) minimum.
   p = 13 (Leech): I re-derived the orchestrator's reduction and re-ran the enumeration: the only dangerous class is N = 4 <-> a minimal
   vector with s = +-6, which the circulant-PSD enumeration excludes (best min eigenvalue -0.032 < 0) [PC, floating point with margin].
   So **Leech is exactly borderline over Delta_13, conditional only on End_{Z[SL_2(13)]}(Leech) being a maximal order** (class number is 1 [P
   given Eichler's mass formula]). p = 37, 61, 73, ...: genuinely OPEN; heuristically a "generic" lattice with this structure has
   (log eps/pi)(pi e/sqrt p)^{p-1} dangerous vectors mod units: 1.6e5 (p=37), 246 (p=61), 2.3 (p=73), 2e-26 (p=157): no obstruction
   expected for large p, some structure needed for p = 37, 61.
4. So: provable obstruction for a positive-density set of p (all p with (p+1)/2 composite), conditional success for p = 13, trivial
   success p = 5 (r = 1), and **genuinely open for the thin (conjecturally infinite) set (p+1)/2 prime**. Nothing found here rules out
   V >= c n^2 p ~ n^{2.5} along such p; the missing pieces are listed in section 5.

## 1. Literature (what was and was not accessible)

Accessible: only web-search snippets. numdam.org is outside the session's allowed network domains (fetch failed); I did not obtain the
full text of ANY of the requested papers. Everything below marked [L-mem] is from memory and must be checked.
- Nebe, "Some cyclo-quaternionic lattices", J. Algebra 199 (1998) 472-498 [bibliographic data confirmed by search]. NOT read. [L-mem]: abstract
  "The faithful lattices of rank 2(p-1) of the groups SL_2(p) are described. For small primes p these and related lattices are investigated
  by computer. In particular a new extremal even unimodular lattice of rank 48 is constructed" (P_48n, group (SL_2(13) o SL_2(5)).2^2, a
  "related" lattice, not of the form W_p). Minima of W_p for p = 17, 29, 37: NOT verified. [L-mem] the catalogue lattices L_{56,2}(M),
  L_{56,2}(M~) (extremal, min 6, rank 56) come from this paper (p = 29).
- p = 41: **correction of the task's attribution** [L, search snippet of Watkins, "Another 80-dimensional extremal lattice", J. Theor. Nombres
  Bordeaux 24 (2012)]: "the unimodular lattice associated to the rank 20 quaternionic matrix group SL_2(F_41) (x) S~_3 in GL_80(Z) is a fourth
  example of an 80-dimensional extremal lattice ... As already noted by Nebe, this lattice has an additional 2-extension in its automorphism
  group." So extremality (min 8) of the SL_2(41)-lattice is Watkins 2012, not Bachoc-Nebe. Bachoc-Nebe (Crelle 494 (1998)) [L, snippet]:
  three extremal lattices of min 8, one 3-modular of rank 40, two unimodular of rank 80, built from the rank-20 Z[(1+sqrt-7)/2]-lattice with
  group 2.M_22.2. Third 80-dim extremal lattice: Stehle-Watkins (SL_2(79), quadratic residue code) [L, snippet].
- p = 37: no extremal 72-dimensional lattice was known before Nebe 2010 (Gamma_72, Aut = (PSL_2(7) x SL_2(25)):2 [L-mem], order not divisible
  by 37); hence the SL_2(37)-lattices of rank 72 studied in 1998 have minimum <= 6 [inference, H].
- p = 13: SL_2(13) < 2.G_2(4) < Co_0 [L-mem]; Leech is the unique even unimodular rank-24 lattice of min 4, so W_13 = Leech if W_13 is even
  unimodular without roots. Tits, "Four presentations of Leech's lattice": NOT read.
- Gross (JAMS 1990), Tiep (globally irreducible representations), Dummigan-Tiep, Scharlau-Tiep, Coulangeon (Acta Arith. 2000): NOT read.
  [L-mem] Gross' "globally irreducible" requires End_G(V) (x) R to be a division algebra (Q, imaginary quadratic, definite quaternion over Q);
  here End (x) R = H x H, so W_p is outside Gross' class - this is why Nebe treats them separately. [L-mem] Coulangeon: minima of tensor
  products of Hermitian lattices over imaginary quadratic fields / quaternion orders are min*min under "small minimum" hypotheses only; no
  general statement for non-commutative orders is known to me.
- Class numbers [P given Eichler's mass formula, computed in `checks.py`]: mass of the totally definite quaternion algebra over Q(sqrt p)
  unramified at finite places = zeta_F(-1)/2 = 1/60 (p=5), 1/12 (13), 1/6 (17), 1/4 (29), 5/12 (37), 2/3 (41), 7/12 (53), 11/12 (61), 11/6 (73).
  p = 13: the maximal order containing Hurwitz (x) O_F has O^1 >= SL_2(3), [O^x : O_F^x] >= 12, so mass 1/12 forces class number 1 [P].
  p = 17: h = 1 [L-mem, Kirschmer-Voight list]. A single class contributes 1/[O^x:O_F^x] = 2/|O^1| (reciprocal of an integer), so
  h >= 2 for p = 37 (5/12), 41 (2/3), 53 (7/12), 61 (11/12), 73 (11/6) [P]; p = 29 (mass 1/4): undecided here. For h >= 2 Lemma 1 must be
  used in its norm-ideal form (ii)-(iii).
- General lower bounds for min(W_p) (square-root bound as for QR codes): none known to me for p = 1 mod 4 [O]. The analogue for
  p = 3 mod 4 (lattices from QR codes, SL_2(p) of rank p+1) has the square-root bound of the code but those families are excluded (R3).

## 2. Theory for general p

### 2.0 Setup and character facts
G = SL_2(p), p = 1 mod 4, eta_1, eta_2 the faithful Weil characters of degree (p-1)/2. Standard table [L]: eta(1) = (p-1)/2, eta(-1) = -(p-1)/2,
eta(c) = (-1 +- sqrt p)/2 (order p), eta(a^l) = 0 (split torus, a^l != +-1), eta(b^m) = (-1)^{m+1} (non-split torus <b> of order p+1, b^m != +-1).
W (x) C = 2(eta_1 + eta_2), dim_Q W = 2(p-1) = 8r, r = (p-1)/4, Delta = quaternion algebra over F = Q(sqrt p) ramified exactly at the two real places.
g of order p: W|_<g> = Q(zeta_p)^2 (no fixed points), eta := sum_j (j/p) g^j acts as sqrt p on the eta_1-part V_+ and -sqrt p on V_- (Gauss sum),
is self-adjoint, central in A; O_F = Z[(1+eta)/2] acts on every G-lattice since (1+sqrt p)/2 in Z[zeta_p] [P].
For w in W_Z: m_+-(w) = |P_+- w|^2, w.eta w = sqrt p (m_+ - m_-) = 2s(w), s integer, m_+- = |w|^2/2 +- s/sqrt p.
Hermitian norm nu(w) in O_F totally positive with |w|^2 = 2 Tr_{F/Q}(beta nu(w)), beta = eps/sqrt p (eps fundamental unit, N(eps) = -1, beta >> 0,
(beta) = inverse different). Then m_+ m_- = 4 N(beta) N(nu) = (4/p) N(nu), i.e. **N(nu(w)) = (p/4) m_+ m_- = p|w|^4/16 - s^2/4** [P], and
s = |w|^2/2 mod 2. Dictionary: kappa := eps nu = (s + (|w|^2/2) sqrt p)/2.

### 2.1 Lemma 1 (line criterion) [P]
Assume O := End_{Z[G]}(W_Z) is a maximal order of Delta and W_Z is unimodular for x.y = Tr_{F/Q}(beta trd h(x,y)).
(i) (O, 2Tr_F(beta nrd)) is even unimodular of rank 8 (= E_8): its O_F-discriminant is the reduced discriminant squared = (1), and Tr_F(beta .) is
unimodular on O_F; det = d_F^4 N(beta)^4 = 1.
(ii) For a Delta-line U = w Delta: W_Z cap U = wJ, J = {d : wd in W_Z} an O-ideal, and det(W_Z cap U) = N(n_U)^4 with n_U := nu(w) nrd(J) an integral
ideal of O_F ("norm ideal of the line"); for every v in U cap W_Z, N(nu(v)) >= N(n_U), with equality for v = generator if J is principal.
(iii) R_U = r^8/det(W_Z cap U) (dim W/dim U = r, covol 1). Hence R_U <= 1 <=> N(n_U) >= r^2.
Consequences: [maximal order] (R_U <= 1 for all lines) => N(nu(w)) >= r^2 = ((p-1)/4)^2 for all w != 0 (because wO is contained in W_Z cap U, so
det(W_Z cap U) <= N(nu(w))^4); [maximal order + all one-sided ideals principal] the converse also holds. Equivalent form:
**(p/4) m_+(w) m_-(w) >= ((p-1)/4)^2 for all w.** R_U = 1 exactly on lines generated by w with N(nu(w)) = r^2.
For p = 13 (r = 3): planes U have U^perp a line and det(W cap U) = det(W cap U^perp) >= 3^8 > (3/2)^16, so the line condition implies the
condition for ALL proper Delta-subspaces [P].

### 2.2 Lemma 2 (necessary size of the minimum) [P]
AM-GM: m_+ m_- <= |w|^4/4, so the line condition for a minimal vector needs p mu^2/16 >= (p-1)^2/16, i.e. **mu >= (p-1)/sqrt p** (3.33, 3.88,
5.20, 5.92, 6.25 for p = 13, 17, 29, 37, 41: so mu >= 4, 4, 6, 6, 8), and minimal vectors must be balanced: s^2 <= p mu^2/4 - (p-1)^2/4, i.e.
|s| <= 4 (p=13, mu=4), 2 (p=17, mu=4), 7 (p=29, mu=6), 3 (p=37, mu=6), 16 (p=37, mu=8), 16 (p=41, mu=8). An extremal minimum (~p/6) is far
more than the necessary sqrt p for large p; the difficulty is with UNBALANCED vectors of larger norm.

### 2.3 Lemma 3 (unit / gamma_8 bound; it never suffices) [P]
(a) eps acts on W_Z; nu(w eps^j) = eps^{2j} nu(w): (m_+, m_-) -> (eps^{2j} m_+, eps^{-2j} m_-), product invariant. Choosing j with
m_+/m_- in [eps^{-2}, eps^2] gives |w eps^j|^2 <= sqrt(m_+ m_-)(eps + 1/eps), hence m_+ m_- >= (mu/(eps+1/eps))^2 and, with eps = (a + b sqrt p)/2,
eps + 1/eps = b sqrt p:  **N(nu(w)) >= mu^2/(4 b^2)**.
(b) Better and independent of b (maximal order): wO = (O, 2Tr(beta nu nrd)) has det N(nu)^4 and min >= mu, and gamma_8 = 2 gives
**N(nu(w)) >= mu^2/4** for all w != 0.
(c) So the minimum alone guarantees the line condition only if mu >= 2r = (p-1)/2 = n/4. For a unimodular lattice of rank n = 2(p-1),
mu <= 2[n/24]+2 (even case) < n/4 for all p >= 13 (and Minkowski: mu <= 0.234 n asymptotically). **For no p >= 13 does the minimum suffice**;
p = 5: mu = 2 = (p-1)/2, fine (r = 1, nothing to prove).
(d) Finer arithmetic criterion [PC, `arith.py`, exact]: list of classes (|w|^2, s) with N < r^2 whose whole unit orbit has norm >= mu
("dangerous classes compatible with mu"): p=13, mu=4: only N = 4, (4, +-6). p=17, mu=4: 12 classes, N in {1,4,8,9,13,...} e.g. (4,+-6) N=8,
(4,+-4) N=13, nu = 3 (N=9, norm 12 -> via O: norm 6). p=29, mu=6: 16 classes, smallest N = 9 (6,+-15). p=37: 54 (mu=6) / 40 (mu=8).
p=41, mu=8: 98. p=61, mu=12: 132. p=73, mu=14: 472. (Counts are w.r.t. the unit orbit only; multiplying by O removes some.)

### 2.4 p = 13 [PC + C]
Norms of O_F: 1, 3, 4, 9, 12, 13 (2, 5, 7 inert). N = 1: nu = eps^{2j}, some w eps^i has norm 2: excluded. N = 3: unit orbit contains norm 2:
excluded. N = 4: nu = 2 eps^{2j} <-> minimal vector with s = +-6 (kappa = 3 + sqrt13); NOT excluded by arithmetic (Leech does contain
sqrt2 E_8). Enumeration (`enum17.py`): a in {-2..2}^6, sum = -2, circulant Gram PSD: admissible s in {0,+-2,+-4} (18, 33, 48, 33, 18 sequences);
for s = +-6 the best minimal eigenvalue is -0.0325 < 0. Hence no minimal vector has s = +-6 and **N(nu(w)) >= 9 = r^2 for all w in Leech** [PC:
floating-point eigenvalues, margin 0.03; uses only |u.v| <= 2 for minimal vectors and g fixed-point-free of order 13].
Conclusion [C: End = maximal order]: Leech is exactly borderline (R_U = 1 on lines with N = 9: s = +-4 minimal vectors, nu = 3; R_U < 1 otherwise).
The orchestrator's design computation (counts 78624/45864/13104) was not re-verified here; the variance E s^2 = 4 is (Lemma 4).

### 2.5 Lemma 4 (design identity) [P]
If the minimal vectors form a spherical 4-design then avg (w.Aw)^2 = mu^2 (2 tr A^2 + (tr A)^2)/(n(n+2)); A = eta: tr = 0, tr eta^2 = pn, n+2 = 2p:
**avg s^2 = mu^2/4** (and avg s = avg s^3 = 0 for a 6-design). Extremal even unimodular lattices: 11-designs (n = 0 mod 24), 7-designs (n = 8 mod 24) [L].
- p = 17 (n=32, mu=4): line condition on minimal vectors needs |s| <= 2; variance 4 then forces |s| = 2 for EVERY minimal vector (N = 16 = r^2
  for all 146880 of them). Not contradictory with the 7-design moments (s = +-2 with probability 1/2 has the right 1st-3rd moments), and
  the circulant enumeration allows s in {0,+-2,+-4,+-6} (`enum17_out.txt`) - so failure on minimal vectors is NOT forced by design theory; but
  p = 17 fails anyway by Lemma A' (order-3 element, fixed lattice of rank 8 and det | 3^8, i.e. N(n_U) <= 9 < 16).
- p = 29 (n=56, mu=6): need |s| <= 7 (s odd), variance 9: not forced. p = 41 (n=80, mu=8): need |s| <= 16, variance 16: not forced.
  Both fail by Lemma A' instead.

### 2.6 Lemma A' (prime-order elements in unimodular lattices) [P]
L unimodular, g in Aut(L) of prime order q, L (x) Q = Q^a + Q(zeta_q)^b, m = a + b(q-1), a, b >= 1, L_1 = L^g, L_z = L cap ker Phi_q(g).
Then L_1^#/L_1 = pi_1(L)/L_1 = L/(L_1 + L_z) = pi_z(L)/L_z is killed by q and by (1 - g), so it is F_q^c with c <= min(a, b);
det L_1 = det L_z = q^c. Hence R_{W^g} = (m/a)^a q^{-c} >= 1, with equality iff a = b = c (L free/projective over Z[C_q]); if a != b then R > 1.
Proof of the inequality: a <= b: m/a >= q and c <= a. a > b: R >= ((1 + x(q-1))/q^x)^a with x = b/a in (0,1), and 1 + x(q-1) > q^x by strict
concavity of log. QED.  (Sharper than R3 Prop. 5, which did not use unimodularity.)
**Cor. A'' [P given the character table].** g in the non-split torus of odd prime order q | (p+1)/2: eta(g^j) = -1, so W|_<g> has
a = 4(N-1), b = 4N, N = (p+1)/(2q) (odd). If N > 1: a != b, a >= 8, and R_{W^g} >= (m/(aq))^a = ((Nq-1)/(Nq-q))^{4(N-1)} > 1.
Values [N]: p=17: 9.99; p=29: 11.8 (q=3), 14.8 (q=5); p=41: 12.5, 17.3; p=53: 13.0; p=89: 13.5, 21.1; p=97: 24.6. So T_U >= ~ R_U^{k/2}:
exponential in k. Scope: stated for unimodular W_p. If (as I believe but did not verify in the literature) every G-lattice in W is W_p J for an
O-ideal J, the families W_p J (x) L coincide with W_p (x) JL and the obstruction is unconditional for SL_2(p) [C].

### 2.7 Lemma B (split torus: exactly borderline) [P]
Let H < T (split torus) have odd order h | p-1. T normalises U_p = <g>, t g t^{-1} = g^{t^2}, so W_Z is a Z[zeta_p]-lattice on which H acts
semilinearly through a faithful action of H on Q(zeta_p) (tame). x := -sum_{j in R} zeta^j (R = representatives of F_p^x/H^2-image) has
sum_{h in H} h(x) = 1; for w in W_Z^H, w = N_H(x w): W_Z^H = N_H(W_Z). W|_H is a multiple of the regular representation (eta = 0 off +-1), so
dim W^H = 2(p-1)/h, det(W_Z^H) = [pi_H(W_Z) : W_Z^H] = h^{dim W^H} and **R_{W^H} = h^{dim}/det = 1 exactly**. (For p = 13, h = 3: fixed lattice of
rank 8 and determinant 3^8, consistent with sqrt3 E_8.) General formula [P]: for unimodular W_Z and any subgroup H,
R_{W^H} = |H^0_Tate(H, W_Z)| / (1 + sum_{h != 1} rho(h))^{dim W^H}, rho = chi_W/chi_W(1).
Classification of fixed spaces for (p+1)/2 = q prime [P given Dickson's list]: a subgroup with nonzero fixed space has odd order, contains
no unipotent, hence is cyclic in a torus; non-split: order q, free. So only Lemma B subspaces occur, all with R = 1.
Remark [C: needs W_Z locally free over Z[x]/(x^{(p-1)/2}+1)]: for the Q-isotypic components W_d of the split torus (eigenvalues of order d),
R_{W_d} = (rd(Q(zeta_d))/phi(d))^{4 phi(d)} <= 1, equality iff d is a power of 2 (which is the Lemma B subspace). No obstruction there.

### 2.8 Heuristic counts [H/N]
(a) Siegel-type count of dangerous vectors modulo units (region m_+ m_- < (p-1)^2/(4p), fundamental domain for eps^2):
Vol = d omega_d^2 c^{d/2} log eps, d = p-1, c = (p-1)^2/(4p), ~ (log eps/pi)(pi e/sqrt p)^{p-1}. Values: 1.2e4 (p=13), 7.5e4 (17), 2.1e5 (29), 1.6e5 (37),
1.3e5 (41), 2.5e3 (53), 246 (61), 2.3 (73), 3e-4 (89), 2e-26 (157). Threshold sqrt p = pi e, i.e. p ~ 73: same as R3 Prop. 4 (t rd_Delta = 2 sqrt p vs 2 pi e).
(b) Minimal vectors, Beta((p-1)/2,(p-1)/2) law for m_+/|w|^2 times kissing number: expected violators 1.2e3 (Leech; truth 0), 2.0e4 (p=17, mu=4),
8.4e4 (p=29), 6.7e3 (p=37 if extremal), 2.1e3 (p=41, K = 1250172000); p=37 with mu=6: 18.6% of directions violate. The Leech case shows that
integrality + positivity kill the Beta tail in highly structured lattices, so (b) is only an upper-bound heuristic.
(c) No modular-form obstruction [H]: theta_{W_p} is a Hilbert modular form of weight (p-1)/2 for SL_2(O_F); the line condition is the vanishing of
the coefficients a(nu), N(nu) < r^2: about (p^{1.5}/8) log eps classes, versus dim M_{(p-1)/2} ~ zeta_F(-1)(p/2)^2/2 ~ p^{3.5}: far fewer conditions
than dimensions (unlike classical extremality).

### 2.9 Asymptotics, Borel subgroup
Needed growth: mu >= (p-1)/sqrt p only (necessary), and N(nu) >= p^2/16 for all vectors (the real condition). Under the Borel subgroup B = U:T the
lattice is a rank-2 Z[zeta_p]-lattice with semilinear T-action; a B-orbit of w spans at most the Z[zeta_p]-line of w (rank p-1), which is NOT a
Delta-subspace (dim_F W over Q(zeta_p) is 2 while Delta-lines have Q-dimension 8), so Borel orbits do not produce small Delta-subspaces; the
Q(zeta_p)-lines only matter through Lemma B. A square-root bound for mu is not known to me and is not what is needed.

## 3. Conditional theorems

**Theorem C1 (first moment) [P modulo L: Siegel-Weil mean value formula for s-tuples (s < k) over a maximal order O of Delta, on the Steinitz
components that occur].** Let (G, W_Z) have A = M_r(Delta), W_Z an O-lattice, and R_U <= 1 for every proper Delta-subspace U, with sum_U c_U^{-(r+1)}
< infinity (true: height zeta, R3 Rem. 1.1). Then for n = k dim W, k -> infinity, there are lattices Lambda = W_Z (x)_O L with V(Lambda) >= |G|(1-o(1)).
Proof sketch: E#(Lambda \ 0 cap B) = V + sum_U T_U; borderline U (finitely many for each s): T_U = f^{-1/2} pi^{(f-1)/2} V^f n^{-(1-f)/2} -> 0;
non-borderline: sum_U T_U <= max_U' R_U'^{(k-r-1)/2} * (convergent sum at k = r+1) -> 0. A point with nontrivial stabiliser H lies in
Lambda_{W^H}, so P(some nonzero point) <= E#(free points)/|G| + E#(points in the union of the Lambda_U) <= V/|G| + o(1). QED.
The o(1) decays only like n^{-(1-f)/2} times the number of borderline subspaces: useless numerically (for Leech at k = 4 the borderline lines
alone contribute ~3 points each). p = 13 would give V_{24k} >= 2184(1-o(1)) - weaker than n^2 for n > 110; C1 is a sanity statement only.

**Theorem C2 (Klartag type) [C/H].** Same hypotheses, k >= C log n: V_n >= (1-o(1)) |G| D/(2e), D = dim of G-invariant forms = 2k(2k-1) for
Delta_R = H^2, n = 2(p-1)k:  **V >= (p(p+1)/(p-1)) (1 - 1/(2k)) n^2/(2e)**.
- p = 5 (icosians, E_8, |G| = 120, r = 1, n = 8k): V >= (7.5 n^2 - 30 n)/(2e) = **1.380 n^2 - 5.5 n** (7.5 x Rezvyakova's 0.184). No lower-rank terms;
  only the division-algebra mean value formula (Gargava arXiv:2107.04844, Thm 16, cited by AGZ) is needed: this case should be a routine
  adaptation of AGZ [C].
- p = 13 (Leech, |G| = 2184, r = 3, n = 24k, k >= 4): V >= (15.17 n^2 - 182 n)/(2e) = **2.790 n^2 - 33.5 n** (15.2 x) [C on: Leech borderline (2.4),
  the mean value formula with lower-rank terms, and the items below].
- p = 37, 61, 73 (if borderline): 39 x, 63 x, 75 x [O].
What must be re-proved (AGZ = arXiv:2606.05105, GGS = arXiv:2610.12619; statement numbers from the pdftotext files):
 1. AGZ Thm 2.1 (mean value over O_K^r) -> formula of R3 Thm 1 over the quaternion order with the terms sum_U c_U^{-k} int_{U_R^k}; template: GGS
    Thm 2.3 (mean value for vectors of prescribed support) - our Lambda_U play the role of GGS's proper-support vectors, and "R_U <= 1" the role
    of GGS's admissibility (Def. 2.6, Cor. 2.8). New inputs: Weil's formula for GL_k of a division algebra with non-principal ideal classes
    (p >= 29), uniform convergence of the height zeta at k >= r+1.
 2. AGZ Lemma 2.3 (projected rank-one norm): q(x) = |P(x (x) x)|^2/|x|^4 in terms of the two quaternionic k x k Gram matrices of x; GGS Lemma 4.1
    (E|P(w (x) w)|^2 = (2D+n)/(n(n+2))) is general and applies verbatim.
 3. AGZ Lemma 3.1, Prop. 3.2 (process in the invariant subspace), Lemma 3.4 (operator norm: GSE blocks instead of GUE), Prop. 3.5, Lemma 3.6:
    routine. AGZ Lemma 3.3 / GGS Prop. 4.3 (dimension loss with weights 1/orbit size): GGS's weighted version is needed because points of
    Lambda_{W^H} have orbits of size |G|/|H|.
 4. AGZ Lemma 4.2, 4.4 (angular block-mass moment) -> exponential moment of q for Wishart-type quaternionic Gram matrices; "sparse" vectors
    are now low-rank Grams (pure tensors have q about r times the typical value); requires k >= C log n as in R3 section 2d [H].
 5. AGZ Prop. 4.5, Lemma 3.7 + GGS Lemma 4.5 (contacts from proper-support vectors) -> contacts from each Lambda_U: heuristically
    T_U(V_0) e^{f z} ~ (wD)^f n^{-(1-f)/2} = O(n^{(5f-1)/2}) per borderline U, negligible against D ~ n^2/(p-1)^2 [H].
 6. AGZ Prop. 5.1 (good lattice): Markov over the family; unchanged.

## 4. Proved / not proved
Proved here: Lemma 1, 2, 3, 4, A', Cor. A'' (given the character table), Lemma B, classification of fixed spaces, C1 modulo the cited mean value
formula, the p = 13 reduction (PC). Conditional: Leech borderline (maximal order), C2. Heuristic: 2.8, constants in C2. Refuted: "p = 17, 29, 41
are candidates" (R3's weaker window test let 29, 41 pass; unimodularity kills them).

## 5. Gaps
- End_{Z[G]}(W_p) maximal? uniqueness of G-lattices up to O-ideals? minima of W_37, W_61, W_73 - all need Nebe's paper (not accessible here).
- For p = 37: is N(nu(w)) >= 81 for all w? With mu = 6 every minimal vector would need |s| <= 3. Needs the explicit lattice.
- Higher-dimensional Delta-subspaces that are not fixed spaces (2 <= s <= r-2) for p >= 37: no criterion beyond Lemma 1 applied to U and U^perp.
- Mean value formula over O with class number > 1; uniformity in p of the height-zeta tail (needed for k ~ p/4, the n^{2.5} regime).
- The floating-point enumeration for s = +-6 should be redone in exact arithmetic (margin 0.03 is large, so this is cosmetic).
