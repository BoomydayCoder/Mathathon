# R3_non_free — symmetric lattice families W_Z (x) L where the group does not act freely

Start 06:48 UTC, budget ~18 min. Code: `check_nonfree.py`, output: `check_output.txt` (all numbers quoted below come from it).
Tags: [P] proved here, [L] literature (cited, not re-proved), [C] conditional, [H] heuristic, [N] numerical, [R] refuted, [O] open/conjecture.

## 0. Setup

K number field (CM or totally real), d=[K:Q], O_K integers, K_R = K (x) R with Euclidean norm |x|^2 = sum_sigma |sigma x|^2
(all d embeddings), so covol(O_K) = |d_K|^{1/2}; rd_K := |d_K|^{1/d}.
W = K^r, W_Z = M O_K^r (M in GL_r(K_R); assume W_Z free, see Remark 1.2), Euclidean structure = Humbert form h=(h_sigma), h_sigma = M_sigma^* M_sigma.
L = O_K^k gamma, gamma in G_1 := {gamma in GL_k(K_R): |N det gamma| = 1}; "random L" = Haar probability on G_1/GL_k(O_K)
(or SL_k(K_R)/SL_k(O_K)).
Lambda = W_Z (x)_{O_K} L = M O_K^{r x k} gamma  in K_R^{r x k}, real dimension n = r k d, covol(Lambda) = c_W^k,
c_W := covol(W_Z) = |d_K|^{r/2} prod_sigma det(h_sigma)^{1/2}.
For a K-subspace U of W of dimension s: c_U := covol(W_Z cap U) (Euclidean, in U_R), f := s/r,
Lambda_U := (W_Z cap U) (x) L = sublattice of Lambda of dimension f n and covolume c_U^k   [P: the map y -> M B_U y gamma, B_U an O_K-basis of
W_Z cap U, is column-wise the injection K_R^s -> K_R^r with volume factor c_U/|d_K|^{s/2}, and right mult. by gamma has real determinant 1].
tau_U^2 := c_U^{2/(sd)} / c_W^{2/(rd)}.   V := Vol_n(B(R))/c_W^k  (so V is the "V(L)" normalisation of the task).

## 1. Mean value formula with lower-rank terms

**Theorem 1 [P modulo L].** Let k > r and f >= 0 Borel on K_R^{r x k}. Then
  E sum_{X in O_K^{r x k}, X != 0} f(M X gamma) = sum_{s=1}^{r} sum_{U in Gr_s(K^r)} c_U^{-k} int_{(MU_R)^k} f dvol_{skd}.
In particular for f = 1_{B(R)}:
  **E #((Lambda \ 0) cap B(R)) = V + sum_{U proper} T_U,   T_U = Vol_{fn}(R) / c_U^k.**
Proof. X != 0 has K-rank s in [1,r]; its column space U = col(X) is an s-dim K-subspace, and X in O_K^{r x k} with col(X) in U
iff X = B_U Y, Y in O_K^{s x k} (B_U basis of W_Z cap U, which is a saturated submodule), rank Y = s. So the sum splits as
sum_s sum_U sum_{Y in O_K^{s x k}, rank s} F_U(Y gamma), F_U(y) := f(M B_U y). Cited input [L]: Siegel's formula for linearly independent
s-tuples, s <= k-1: E sum_{Y rank s} F(Y gamma) = int_{K_R^{s x k}} F(y) dy (dy normalised with covol(O_K^{s x k}) = 1).
(K=Q: Siegel 1945 / Rogers 1955; general K: Weil 1965 "Sur la formule de Siegel..."; explicit module-lattice versions incl. higher moments:
Gargava-Serban-Viazovska arXiv:2303.15275 / Hughes 2023; s=1: Venkatesh 2013 Lemma.) Change of variables y -> M B_U y has Jacobian
(c_U/|d_K|^{s/2})^k w.r.t. Lebesgue measure with covol(O_K)=|d_K|^{1/2} per coordinate, giving c_U^{-k} int. QED
Remark 1.1. sum_U c_U^{-k} over Gr_s(K^r) is a height zeta function; it converges iff k > r (Schmidt/Thunder: #{U: c_U <= X} ~ c X^r). So k >= r+1,
hence n >= r(r+1)d always.
Remark 1.2. If W_Z cap U is not free (class number > 1) the U-term needs the Steinitz-class version of the formula; not checked [gap].

**Proposition 2 [P], checked [N] (ratios 0.998-0.9997 for n=200..800).** With V fixed as above,
  T_U = V^f (pi n)^{f/2} (pi f n)^{-1/2} ( (1/f) / tau_U^2 )^{f n/2} (1 + O(1/(fn))).
(Stirling on both Gamma factors.) Hence for V = poly(n):
- f tau_U^2 >= 1 + C log n / n for all proper U (and the height-zeta tail controlled) => lower-rank terms o(1);
- f tau_U^2 = 1 exactly ("borderline") => T_U = sqrt(1/f) pi^{(f-1)/2} V^f n^{-(1-f)/2}: polynomial;
- f tau_U^2 <= 1 - eps for one U => T_U >= exp(eps f n /2): exponential.
"Super-stability": c_U^{2/dim U} dim U >= c_W^{2/dim W} dim W for all proper K-subspaces U. The orchestrator's criterion is CONFIRMED.
(Compare Grayson-Stuhler semistability c_U^{1/dim U} >= c_W^{1/dim W}: super-stability is stronger by the factor dim W/dim U.)

**Proposition 3 (no-go) [P / C].** (a) [P, deterministic] Lambda_U is a lattice of dim fn, so V(Lambda_U) <= V_{fn} <= 2^{kappa fn}
(kappa = 0.401 KL, asymptotically). With R = lambda_1(Lambda_U) >= lambda_1(Lambda), Prop. 2 gives for EVERY L (random or not)
   V(W_Z (x) L) <= poly(n) * ( 2^{2 kappa} f tau_U^2 )^{n/2}.
So if f tau_U^2 < 2^{-0.802} = 0.5735 for one U, the whole family is exponentially worse than Minkowski-Hlawka.
(b) [C: needs a second-moment/Rogers bound P(V(L') >= v) <= C/v for the rank-k family attached to U] if f tau_U^2 < 1 then
V(W_Z (x) L) <= (f tau_U^2)^{n/2} poly(n) with probability -> 1. For K = Q, s = 1 this is Rogers' classical variance bound, so (b) is [P] there.
(c) In all cases with f tau_U^2 < 1 - eps the first moment E# is exponentially large, so Markov/first-moment and (a fortiori) any
Klartag-type argument started from a Haar-random member of the family fail.

## 2. Lines: Hermite-Humbert invariant and the root discriminant

For U = wK (w primitive): c_U = |d_K|^{1/2} (N h(w))^{1/2}, N h(w) := prod_sigma h_sigma(w,w); so the discriminant cancels:
  **tau_{wK}^2 = N h(w)^{1/d} / (N det h)^{1/(rd)}**   (Hermite-Humbert quotient, Icaza 1997), and super-stability on lines <=> this is >= r for all w.
AM-GM: N h(w)^{1/d} <= |w|^2/d, hence
  **Prop 4 [P].** line super-stability => gamma(W_Z) := lambda_1(W_Z)^2/c_W^{2/(rd)} >= rd / rd_K, hence rd_K >= rd/gamma_{rd}
  (gamma_m Hermite constant). As rd -> infinity: rd_K >= 2 pi e / 1.744 = 9.79 (KL), and >= 2 pi e = 17.08 if W_Z is only as dense as Minkowski-Hlawka.
  The orchestrator's root-discriminant requirement is CONFIRMED (c_0 in [1, 1.744]).
Non-commutative Delta (same proof, T(x) = Tr_{Delta/Q}(x x^*)/t, t^2 = dim_Z Delta, D = dim_Q Delta, rd_Delta := covol_T(O_Delta)^{2/D}):
  gamma(W_Z) >= dim_Q(W) / (t rd_Delta)   [P, AM-GM on the eigenvalues of h(w,w) in Delta_R].
**Corollary 4.1 (no-go for Delta = Q) [P].** If W is absolutely irreducible over Q with Delta = Q and r = dim W >= 2, then
f tau_line^2 = gamma(W_Z)/r <= gamma_r / r <= 2/pi < 1 (Minkowski gamma_r <= 4 omega_r^{-2/r} <= 2r/pi, using Gamma(1+x) <= x^x, x>=1).
So for every G-invariant W_Z in Q^r the family W_Z (x)_Z L has exponentially large first moment and, by Prop 3(b) (Rogers, K=Q), exponentially
small V with probability -> 1. This kills: Z^2/D_8, A_2/D_12, D_4/2_+^{1+4}, E_8/2_+^{1+6}, all Weyl groups, S_n on A_{n-1}, Barnes-Wall with the
real Clifford group, every monomial group over Q.

## 3. Fixed spaces of prime-order elements: entropy inequality

g in G of prime order q, W|_<g> = Q^a + Q(zeta_q)^b (Q-dims a and b(q-1), m = a + b(q-1)), U = W^g, U' = ker Phi_q(g) = U^perp (both Delta-subspaces).
[P] I := [W_Z : (W_Z cap U) + (W_Z cap U')] <= q^{min(a,b)}: the quotient embeds into pi_U(W_Z)/(W_Z cap U) (killed by q since q pi_U = N_g maps W_Z to
W_Z cap U; rank <= a) and into pi_U'(W_Z)/(W_Z cap U') (killed by 1-g, a quotient of a Z[zeta_q]-module of rank b by (1-zeta): rank <= b).
c_U c_U' = I c_W. Multiplying super-stability for U and U':  c_U c_U' >= c_W (m/a)^{a/2} (m/(m-a))^{(m-a)/2} = c_W exp((m/2) h(a/m)). Hence
  **Prop 5 [P].** super-stability => m h(a/m) <= 2 min(a,b) log q   (h = entropy in nats).  CONFIRMED.
  q = 2: h(x) > 2x log 2 on (0,1/2) by strict concavity, so a = b, I = 2^a (W_Z free over Z[<g>]) and c_U = c_U' = 2^{a/2} c_W^{1/2}, i.e. f tau^2 = 1
  EXACTLY for both U and U' (borderline, never strictly super-stable).
  [N] admissible windows for a/b: q=3: [0.709,1.283]; q=5: [0.460,1.497]; q=7: [0.343,1.587]; q=11: [0.228,1.672]. Equivalently the normalised
  character rho(g) = (a-b)/m must lie in about [-0.11,0.09] (q=3,5), [-0.10,0.08] (q=7): non-free prime-order elements must have |rho(g)| <~ 0.1,
  SMALLER than the free value |rho| = 1/(q-1). (Consequence: elements with fixed points contribute little to g = sum rho^2.)
Quantified failure [P]: if the inequality fails by Delta := (m/2)h(a/m) - min(a,b) log q > 0 then (f tau_U^2)^{dim U/2} (f' tau_U'^2)^{dim U'/2} <= e^{-Delta},
so T_U or T_U' is >= e^{Delta k /2}-ish (exponential in k).

## 4. Examples (explicit, `check_nonfree.py`)

| W_Z, G, Delta, r | U | tau_U^2 | needed 1/f | f tau^2 | verdict |
|---|---|---|---|---|---|
| Z^2, D_8, Q, 2 | line (1,0) | 1 | 2 | 0.500 | FAIL (Lambda = L+L), deterministic no-go |
| A_2, D_12, Q, 2 | min line | 1.1547 | 2 | 0.577 | FAIL |
| D_4, 2_+^{1+4}, Q, 4 | min line | 1.4142 | 4 | 0.354 | FAIL, deterministic |
| D_4, 2_+^{1+4}, Q, 4 | fixed plane of diag(1,1,-1,-1) | 1.4142 | 2 | 0.707 | FAIL; index 2 < 2^a = 4 |
| E_8, 2_+^{1+6}, Q, 8 | min line | 2 | 8 | 0.250 | FAIL, deterministic |
| O_K^p, Heisenberg p^{1+2}, Q(zeta_p), p | K e_1 | 1 | p | 1/p | FAIL, deterministic (p=3,5,7,...) |
| D_4 as Z[sqrt-2]-lattice, rank 2 | K-line through min vector | 2 | 2 | 1.000 | BORDERLINE (|d_K|=8 is exactly rd/gamma_4 squared) |
| O_E over O_K, E=Q(sqrt-3,sqrt5), K=Q(sqrt-15) | U = K | 2 | 2 | 1.000 | BORDERLINE |
| Craig A_{p-1}^{((p+1)/4)}, PSL_2(p), Q(sqrt-p), (p-1)/2 | min line | (p+1)/4 | (p-1)/2 | 0.667, 0.600, 0.556 (p=7,11,19) -> 1/2 | FAIL |
Heisenberg, ANY invariant O_K-lattice of rank p (Prop 4, Blichfeldt bound gamma_m <= (2/pi)Gamma(2+m/2)^{2/m}, rigorous):
p=3: tau^2 <= 1.59 < 3; p=5: <= 3.06 < 5; p=7: <= 5.40 < 7  => FAIL for p = 3,5,7 whatever W_Z [P].
p >= 11: Blichfeldt no longer excludes (12.65 > 11); at Minkowski-Hlawka density the line condition becomes satisfiable only for p >= 23 [N].
D_4 / E_8 with Q_8 resp. as modules over Hurwitz / icosian orders are r = 1 (free) - that is Vance's family, no proper subspaces.

## 5. Structural results on borderline families

**Prop 6 (unramified base change) [P for lines and subfields].** E/K Galois of degree r, unramified at all finite primes, W = E, W_Z = O_E,
h(x,y) = Tr_{E/K}(x ybar), G = mu_E x| Gal(E/K) (K-span = crossed product = M_r(K)). Then N det h = d_E/d_K^r = 1 and for every w in O_E \ 0
  tau_{wK}^2 = N_K(Tr_{E/K}(w wbar))^{1/d} >= r N_{E/Q}(w)^{2/(rd)} >= r   (AM-GM), equality for w in mu_E O_K^x;
for an intermediate field F, s = [F:K]: tau_F^2 = (r/s) rd_F / rd_E = r/s. So all lines are super-stable and the subfield subspaces (= fixed spaces
of subgroups of Gal) are exactly borderline. Numerically confirmed for K=Q(sqrt-15) < E=Q(sqrt-3,sqrt5): covol 15, c_K = 2 sqrt15, tau^2 = 2.
Gain: g = |G|/dim_Q A = r|mu_E|/(r^2 d) = |mu_E|/[E:Q]: IDENTICAL to the free family over E. So the extra Galois symmetry buys nothing for a
Klartag-type bound. For the first-moment method it would give V ~ |G| = r |mu_E| with n = k[E:Q], k >= r+1, i.e. V/n -> (r/(r+1)) m/phi(m),
up to twice Venkatesh's constant [H: requires bounding sum over all non-borderline U (height zeta with k = r+1, suppression only N(w)^{-k/r}
per line) - NOT done].
(i)-(iii) of the task: (i) large rd_K removes the line obstruction: by Prop 4 it is necessary and (random W_Z, Hlawka-Thunder) essentially sufficient
that rd_K >~ 2 pi e; (ii) a random O_K-lattice W_Z is line-super-stable when rd_K > 2 pi e (1+eps) [H] but has G = mu_K only (scalars), g = no gain;
symmetry forces short "structured" vectors (Prop 5: every non-central involution forces exact borderline, every Delta = Q factor kills);
(iii) W_Z (x) L symmetric: lines condition is exactly "split tensors w (x) l are not shorter than the generic radius".

## 6. Group-theoretic part

**Prop 7 [P].** W Q-irreducible with complex constituent chi, Schur index m_Q, e = [Q(chi):Q], chi_W = m_Q sum_{sigma} chi^sigma, rho = chi_W/chi_W(1).
sum_g rho(g)^2 = |G| <chi_W,chi_W> / chi_W(1)^2 = |G| m_Q^2 e / (m_Q e chi(1))^2 = |G|/(e chi(1)^2) = |G|/dim_Q A,
since A (x) C = prod_sigma M_{chi(1)}(C). Checked [N]: C_m (m/phi(m)), D_8 (2), D_12 (3), SL_2(5) on icosians (15).
**Correction [R].** "G acts freely on W \ 0 => A = Delta division, r = 1" is FALSE. Counterexample: G = SL_2(5) x C_7 on chi_2 (x) lambda:
fixed-point free, Q(chi) = K = Q(sqrt5, zeta_7) totally complex with even local degree at 2, so (-1,-1)_K splits, A = M_2(K), r = 2, dim_Q W = 24.
(Consistent with Amitsur: SL_2(5) x M is not a subgroup of a division ring for |M| > 1.) Free action does make G a Frobenius complement.
**Prop 8 [P modulo L (Zassenhaus structure theorem)].** For a fixed-point-free chi and any cyclic C <= G of order c: chi|_C is a sum of chi(1) faithful
characters; the stabiliser of chi in Gal permutes this multiset and acts freely on faithful characters, so its image in (Z/c)^x has order <= chi(1);
hence e >= phi(c)/chi(1) and
  g = |G|/(e chi(1)^2) <= ([G:C]/chi(1)) * c/phi(c).
Frobenius complements have a cyclic subgroup of index <= 24 chi(1) [L: Zassenhaus/Passman: normal Z-subgroup C_a x| C_b of index <= 24 resp.
SL_2(5) x Z-group of index <= 2; the Z-group has a cyclic subgroup of index = degree of its fpf characters]. So g <= 24 c/phi(c) = O(log log |G|), and
cyclic groups of primorial order give g = m/phi(m) ~ e^gamma log log m. **sup g = Theta(log log |G|) over fpf groups.**
Constant: G = SL_2(5) x C_N (gcd(N,30)=1) has g = 15 N/phi(N) ~ 4 e^gamma log log N, four times the cyclic constant - but there r = 2 (see correction),
so super-stability of an invariant rank-2 O_K-lattice, K = Q(sqrt5, zeta_N), must be checked [O]. Division-algebra case SL_2(3) x C_N: g = 6N/phi(N) but
N restricted to 2 of odd order mod N.

## 7. Candidate families beyond fpf: status

- PSL_2(p), p = 3 mod 4, Weil degree (p-1)/2 over Q(sqrt-p) (g ~ 2p): **NO-GO [P given Prop 5]**: the degree is odd, involutions of PSL_2(p) are
  non-scalar, so a != b for q = 2; deficiency Delta ~ log 2, T_U >= 2^{k/2}/poly. (Also the Craig lattice fails on lines, f tau^2 -> 1/2.)
- SL_2(p), p = 3 mod 4, Weil degree (p+1)/2 over Q(sqrt-p) (Leech p=23, P_48q p=47): unipotent u has (a,b) = (2,1), q = p, m = p+1:
  m h(2/m) = 2 log(p+1) - 2 log 2 + 2 + o(1) > 2 log p. **NO-GO [P given Prop 5]**, deficiency Delta -> 1 - log 2 = 0.31, lower-rank term >= e^{c k}.
- SL_2(p), p = 1 mod 4, faithful Weil degree (p-1)/2, Delta = definite quaternion algebra over Q(sqrt p) unramified at finite places, r = (p-1)/4,
  dim_Q W = 2(p-1) (E_8 for p=5 with r=1; Leech for p=13 reportedly): unipotents and order-4 elements act freely, torus elements have a = b or
  a = b - 1 (inside the windows), and the line condition only asks gamma(W_Z) >= (p-1)/sqrt p because t rd_Delta = 2 sqrt p is large.
  Detail [P given the standard Weil-character restriction, recalled not re-derived]: split torus, odd prime q | p-1: regular, a = b; non-split torus,
  odd prime q | p+1, N := (p+1)/(2q) (odd): (a,b) proportional to (N-1, N): free if N = 1, inside the window iff (N-1)/N >= lower end, which FAILS
  for q = 3, N = 3, i.e. p = 17 (0.667 < 0.709) and holds for p = 13 (free), p >= 29.
  Sanity check of Prop 5 on a known lattice: Leech with SL_2(23) (previous bullet): fixed lattice of an element of order 23 is 2-dimensional of
  determinant 23, f tau^2 = sqrt(23)/12 = 0.40 < 0.5735, consistent with the predicted failure.
  **Passes every necessary test derived here (p != 17) [O]**; g ~ 2p ~ 2 sqrt(2n) (n >= (p-1)(p+3)/2). Not verified: other Delta-subspaces (fixed spaces of
  non-prime-order subgroups, spans of several minimal vectors), existence of invariant lattices with the needed minimum for all p, the height-zeta
  tail at k = r+1, and of course any equivariant Klartag argument. Treat as the one concrete viable-looking candidate; a first-moment bound would
  give V ~ |G| ~ (2n)^{3/2} (< Klartag's n^2), an equivariant Klartag heuristic n^2 g ~ n^{2.5}. Too good not to be suspicious.
- Heisenberg p^{1+2} (g ~ 2: no gain anyway) and Clifford p^{1+2}.SL_2(p) over Q(zeta_p) (g ~ 2p^3 ~ 2n): excluded for p <= 7 by Prop 4; for larger p
  the standard O_K^p fails (f tau^2 = 1/p) and whether some invariant lattice is super-stable is [O]; SL_2(p) acts through elements with fixed points
  (permutation-like torus) which must satisfy Prop 5 with the right index.

## 8. Verdict

No universal no-go "r >= 2 => fails" exists: Prop 6 and D_4/Z[sqrt-2] are borderline families with r >= 2. Proven no-go statements: (1) Delta = Q
(Cor 4.1); (2) small root discriminant (Prop 4); (3) any prime-order element violating the entropy window, in particular any involution with
a != b or non-free Z[C_2]-structure (Prop 5); (4) deterministic exponential loss when f tau_U^2 < 0.5735 (Prop 3a). Proven structural fact: every
element with fixed points must have |rho(g)| <~ 0.1 and involutions with fixed points are at best exactly borderline, so the gain
g = sum rho^2 cannot come from elements with fixed points being "large"; whether sum over many small rho can still exceed log log is open
(SL_2(p), p = 1 mod 4 is the test case).
