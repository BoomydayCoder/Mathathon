# R1 literature notes: symmetric lattice families (Venkatesh / Vance / Gargava et al.)

Agent start 06:48 UTC, notes written ~06:58 UTC (2026-11-25). Normalisation as in the task: V(L) = volume of
largest lattice-point-free centred ball, covolume 1; density = V/2^n; c_d (Gargava) = V_n.

## 0. Access status (honest)

| paper | access |
|---|---|
| Vance, arXiv:1105.3779 | full text read (papers/vance_1105.3779/paper.txt) |
| Gargava, arXiv:2107.04844 | intro + Sec. 2, 2.3, 3 read in full; Sec. 1 (proof of mean value thm) only skimmed |
| Gargava-Serban, arXiv:2111.03684 | abstract + introduction read; body not read |
| Gargava-Viazovska, arXiv:2411.14973 | abstract + introduction read (Thm 1, Cor 2, Thm 7); proofs not read |
| Gargava-Serban-Viazovska, arXiv:2308.15275 | abstract + introduction read (Rogers thm, Thm 1, Thm 2, Cor 3); proofs not read |
| Moustrou, arXiv:1603.00743 | abstract + introduction read |
| Gaborit-Zemor, arXiv:math/0605098 | abstract + introduction + Sec. 2 start read |
| Venkatesh, IMRN 2013 | NOT ACCESSED: math.stanford.edu / math.ias.edu are outside the session's allowed network; not on arXiv. Everything about it below is from memory + how the other papers quote it. |
| NEW: Abuya-Gargava-Zhao, arXiv:2606.05105 (June 2026) | abstract + intro + Sec. 2.1-2.2, Remark 4.3 read |
| NEW: Gargava-Groebner-Simkin, arXiv:2610.12619 (Oct 2026) | abstract + intro + Appendix A start read |
| NEW: Bourbaki-style survey "Sphere packings in higher dimension (after Klartag)", arXiv:2606.13313 | only grepped |

## 1. STATE OF THE ART (important correction to the task's context)

* Abuya-Gargava-Zhao, "Stochastically evolving ellipsoids with symmetries", arXiv:2606.05105.
  Thm 1.1: universal c>0 and an infinite sequence of N with Delta^lat_N >= c N^2 loglog N 2^-N.
  Setup: K=Q(zeta_m), m primorial, module rank r = ceil((log phi(m))^4), N = phi(m) r, lattices g O_K^r,
  g in SL_r(K_R). Klartag's process restricted to the K_R-linear symmetric forms
  Q_K = (+)_{j=1..phi(m)/2} Herm(C^r), dim_R Q_K = phi(m) r^2/2 = N r/2.
  Mechanism (their words): "each new contact gives m lattice points for the price of one"; volume gain
  >= c m phi(m) r^2 = c (m/phi(m)) N^2.  So: (#invariant degrees of freedom) x (orbit size).
  Mean value formula used (Thm 2.1): for r>=2, int over SL_r(K_R)/SL_r(O_K) of sum_{0 != v in gO_K^r} f(v)
  = covol(O_K^r)^-1 int f.  Free action: "h-1 is invertible ... every nonzero mu_m-orbit has exactly m points".
  Remark 4.3 explains why r must grow: for fixed r the angular exponential moment is exp((1+1/r)lambda), giving
  an integral ~ (mNr)^{1/r}, unbounded (this is why the rank-2 attempt only interpolated between Venkatesh and Klartag).
  The intro records Venkatesh's conjectured UPPER bound Delta^lat_N <= N^2 (log N)^{O(1)} 2^-N ("lattices have ~N^2 degrees of freedom").
* Gargava-Groebner-Simkin, "A cyclotomic gluing construction for lattice sphere packings", arXiv:2610.12619.
  Thm 1.1: Delta^L_d >= c d^2 loglog d 2^-d for EVERY d >= d_0. Appendix Thm A.1: without stochastic methods,
  Delta^L_d >= (1-o(1)) (d M/(e C)) 2^-d (Rogers/Minkowski successive-minima stretching), M/C >> loglog d.
  Construction: tau acts with Q[tau] = prod_{n in N} K_n (product of several cyclotomic fields, NOT a field):
  a glued finite-index sublattice of a product of cyclotomic module lattices. THE ACTION IS NOT FREE IN THE
  UNIFORM SENSE: a vector of support S has orbit size M_S = lcm(n in S) < M. "Proper-support vectors may have smaller
  orbits, but the gluing estimates ensure that their total contribution is negligible"; density ~ 2^-d d^2 M/C with
  M = order of the cyclic symmetry, C = sum phi(n) = rank of the cyclotomic order. This is the first paper in this
  line that handles non-uniform orbit sizes (answer to (c): yes, but only for commutative semisimple algebras, and
  the gain is still M/C = O(loglog)).

## 2. Exact theorems and mean-value formulas (question (a))

**Vance 2011 (full text).** Hurwitz lattices = left modules over H=Z[i,j,(1+i+j+k)/2] in H^m (class number 1 => free).
- Thm 2.2 (Minkowski analogue): Lambda contains a quaternionic basis v_i with |v_i| = min_i(Lambda); if det=1 there is
  a Hurwitz lattice Lambda~ of det 1 with |Lambda~| = (prod min_i)^{1/m}.
- Thm 2.3 (Hlawka analogue): f bounded Riemann integrable compact support, m>=2, eps>0 => exists Hurwitz Lambda, det 1,
  sum_{u in Lambda\0} f(u) < int f + eps. Cor 2.4: same with primitive vectors and factor 1/zeta(4m).
  Proof is Davenport-Rogers style (explicit family), not Haar measure.
- Thm 3.1/3.2: Delta_{4m,L} >= 3 m zeta(4m) / (2^{4m-3} e (1-e^-m)), i.e. in our normalisation
  V_{4m} >= 24 m zeta(4m)/(e(1-e^-m)) ~ (6/e) n ~ 2.2 n with n = 4m. Group used: H^x, |H^x|=24, acting freely by scalar multiplication
  ("every Hurwitz lattice has at least 24 minimal vectors"). Uses Rogers' radial test function rho + successive minima.
- Thm 4.1: for an H^x-invariant convex body, MH improved by factor 12. Gaussian/Eisenstein analogues give only
  m zeta(2m)/2^{2m-2}, 3m zeta(2m)/2^{2m-1}: "the Hurwitz integer unit group is sufficiently larger".

**Gargava 2021 (Math. Z. 2023).**
- Thm 16 (Siegel for division algebras): D a Q-division algebra, O an order, G=SL_k(D_R), Gamma=SL_k(O), k>=2,
  f in C_c(D_R^k): int_{G/Gamma} sum_{v in gO^k\0} f(v) dg = int_{D_R^k} f dx.
- Prop 20/Thm 21: G_0 a finite subgroup of O^x => c_d >= #G_0 for d = 2 dim_Q D. Key sentence in the proof:
  "the G_0-orbit of any non-zero element of O^2 is of size #G_0 because O^2 and G_0 are made of elements of the
  division algebra D" -> sum is in {#G_0, 2#G_0, ...}, mean < #G_0 => empty. G_0 acts by right mult. v -> v g^-1, the
  ball is for a G_0-averaged form. "Tightening": may replace D by Q<G_0>.
- Thm 26 (Amitsur 1955, as quoted): finite subgroups of Q-division algebras are: D* (48, dim 16), I* (120, dim 20),
  G_{m,r} (order m ord_m r, dim phi(m) ord_m r), T* x G_{m,r} (order 24 m ord_m r, dim 4 phi(m) ord_m r; m odd, ord_m 2 odd).
- Thm 30: c_{8 phi(m)} >= 24 m for m = product of primes p<=x with ord_p 2 odd (density 7/24 by Hasse) =>
  c_d > C d (loglog d)^{7/24}. Prop 31: c_{4 phi(m)} >= 2m (m primorial, r=m-1, non-commutative dicyclic-type group)
  => c_d >= (1/2) d loglog d eventually.
- EXPLICIT BARRIER STATEMENT (end of Sec. 3, verbatim): "no sequence constructed using Theorem 21 can give us an
  asymptotic growth strictly better than O(d log log d). Indeed, looking at Theorem 26, we observe that
  |G_0|/dim_Q D can at most be 3m/phi(m) for some sequence of integers m."

**Gargava-Serban 2021.** Thm 1.1: A central simple DIVISION algebra over K, [A:K]=n^2, m=[K:Q], O an O_K-order,
G_0 < O^x finite, t>=2: lattice in dim n^2 m t with Delta >= |G_0| zeta(mn^2t) t / (2^{mn^2t} e (1-e^-t)) (Rogers-Vance
trick, rank t), effective via lifts of codes mod p (Thm 3.7 = finite Siegel mean value). Intro says quaternion
algebras over cyclotomic fields do NOT combine Vance x Venkatesh "due to a parity condition": only 3 (loglog n)^{7/24} n 2^-n.

**Gargava-Viazovska 2024.** Thm 1: K cyclotomic, B ball of volume V in K_R (trace form), Lambda Haar-random in the
Arakelov class group (unit-covolume IDEAL lattices, rank ONE): E #(B cap Lambda) = 1 + V + eps, |eps| <= c sqrt(V) Delta_K^-eta
(eta<1/12 via Petrow-Young subconvexity; Hecke integration formula). Cor 2: c_{phi(n)} >= n - c_1 e^{-c_2 phi(n) log n};
hence c_d >= d loglog d infinitely often (factor 2 over Venkatesh because rank 1 instead of rank 2: d=phi(n) vs 2phi(n)).

**Gargava-Serban-Viazovska 2023 (moments).** O_K-lattices of rank t, SL_t(K_R)/SL_t(O_K). Thm 1/2: n-th moment
omega_K^n m_n(V/omega_K) <= E[rho^n] <= omega_K^n m_n(V/omega_K) + E_{n,t,K}(V+1)^{n-1}, m_n = Poisson moment,
omega_K = # roots of unity, valid for t > t_0 = O(n^3 loglog n), error exponentially small in d(t-t_0) under a
Bogomolov (height gap) hypothesis (holds for cyclotomic towers). Cor 3: t>=27: V^2+V omega_K <= E rho^2 <= V^2+V omega_K(1+C e^{-eps d(t-t_0)}).
I.e. the number of omega_K-ORBITS in the ball is asymptotically Poisson(V/omega_K). Prop 58: without the height gap
the moments are strictly larger than Poisson.

**Moustrou 2016.** Effective Venkatesh: Thm 1: infinitely many n, Delta >= 0.89 n loglog n / 2^n constructible in
exp(1.5 n log n) operations; Construction A over Z[zeta_m]^2 -> (Z[zeta_m]/P)^2, all F_q-lines; mean value = finite
counting (Lemma 4). Explicitly: the cyclic group action "is in our case easier to deal with, because it is a free
action", whereas Gaborit-Zemor (double circulant codes, automorphism group of size n, n=2p, c ~ 0.06) had a
non-free permutation action and had to restrict to prime order.

**Venkatesh 2013 (FROM MEMORY, not re-read).** Thm 1: V_n >= 65963 n for n large; V_n >= (1/2) n loglog n for infinitely
many n (n = 2 phi(m), m primorial; lattices = rank-2 O_K-modules, K=Q(mu_m), Siegel-type mean value on
SL_2(K_R)/SL_2(O_K), balls for a mu_m-invariant form, orbits of size m => empty if V<m). I believe the general
statement there is for a Q-group/representation with a mean value formula and that the note mentions the Hurwitz
case of Vance; I cannot verify the exact wording.

## 3. Why the gain stops at loglog n (question (b))

Literature statement: Gargava (above) - with division algebras |G_0|/dim_Q D <= 3 m/phi(m) = O(loglog). GGS'26 phrase the
figure of merit as M/C (order of cyclic symmetry / rank of the cyclotomic order).

My own derivation (INFERENCE, not in the papers; representation-theoretic bookkeeping):
- First-moment method: V ~ |orbit| and N = r * dim_Q(simple module), r >= 2.
- Klartag-type method: V ~ c * |orbit| * dim(Sym^2 V)^G. If V_C = sum m_i rho_i then
  dim(Sym^2 V)^G ~ (1/2) sum m_i^2, N = sum m_i d_i. For V = S^r, S the simple Q[G]-module of a faithful
  irreducible chi with Schur index s: multiplicities m_i = s r over [Q(chi):Q] conjugates, hence
      V / N^2  ~  (c/2) * Phi(G,chi),   Phi(G,chi) := |G| / ( [Q(chi):Q] * chi(1)^2 ),
  INDEPENDENT of the Schur index and of r. For mu_m: Phi = m/phi(m). For Amitsur's G_{m,r}: Phi = M/phi(M).
- Fixed-point-free groups are a strictly larger class than Amitsur's list. Example I checked by hand:
  G = Z/19 x| Z/9 (order 171, r=7) acts freely on Q^36 \ 0; Amitsur's criterion (Gargava Thm 24) FAILS
  (p=19: gcd(3,(19-1)/3)=3), the enveloping algebra is M_3(F), [F:Q]=12, Schur index 1 (index divides 3 and
  is not 3). Numerically verified (numpy, 06:58 UTC): the 3-dim induced representation of the group of order 171
  has no non-identity element with eigenvalue 1. In general
  G = Z/p x| Z/q^2 with p = 1 mod q^2: |G| = p q^2 acting freely in dimension (p-1)(q-1): |G|/dim ~ q (up to n^{1/3})!
  BUT this does not help: the simple module is F^q (matrix size k=q), the lattice family is W (x)_{O_F} L with
  rank L = r, and vectors of full tensor rank exist only if r >= k; Rogers-type convergence needs r > k. Then
  N = r(p-1)(q-1) > q(p-1)(q-1) and |G|/N < pq/phi(pq). Equivalently Phi(G,chi) = M/phi(M) again (type I
  metacyclic f.p.f. groups: |G| = eM, chi(1)=e, [Q(chi):Q]=phi(M)/e). So for FREE actions the barrier
  Phi = O(loglog) seems robust (other Zassenhaus types differ by bounded index) - conjectural general statement.
- For NON-free actions Phi can be astronomically large: S_{k+1} on A_k (Phi=(k+1)!/k^2), Clifford group
  p^{1+2a}.Sp_{2a}(p) on Barnes-Wall type lattices (Phi ~ k^{2 log_p k}), monomial groups mu_m wr S_k.

## 4. Non-free actions / tensor families W (x) L (question (c))

Found in the literature: (i) Gaborit-Zemor (non-free cyclic permutation action, double circulant codes, only
constant gain); (ii) GGS'26 gluing (non-uniform orbit sizes for a product of fields, handled by forcing the
proper-support sublattices to have large covolume). I found NO paper using W (x)_{O_K} L with W a fixed highly
symmetric Hermitian lattice, nor M_r(K)-module families, for packing bounds (not searched exhaustively:
Coulangeon 2000, Gross 1990, Tiep, Nebe, Bachoc, Elkies/Shioda were NOT consulted in this 15-minute slot).

Obstruction (my derivation, INFERENCE; matches the classical fact min(W (x) L) <= min(W) min(L), Kitaoka/Coulangeon):
- Mean value for {W (x) gL}: GL_r(K) is not transitive on W (x) K^r \ 0; orbits = tensor rank / column space U in W_K.
  The first moment is a Rogers-type sum (GSV Thm 4 with n = rank W = k, t = r):
     E #(W(x)L cap B) = sum_{j<=k} sum_{W' < W primitive, rank j} vol_{jrd}(B_R) |d_K|^{jr/2} / covol(W')^r  (r > k).
  The j = k term is V; orbit of a full-rank tensor has size |G| (g fixes sum w_i (x) l_i iff g fixes all w_i).
- The rank-j terms are negligible (on the exponential scale) iff W is "super-stable":
     delta(W') >= (k/j) delta(W) for all O_K-sublattices W' of rank j < k,  delta := N(det h)^{1/(j d)}.
  For j=1: min_w N_{K/Q}(h_W(w))^{1/d} >= k * delta(W).
- Over K = Q this is IMPOSSIBLE for k >= 2: it needs gamma(W) >= k but gamma_k <= 2/pi + k/pi < k (Blichfeldt).
  Rigorous version: gamma(W (x)_Z L) <= gamma_k gamma_r <= 0.178 kr/(2 pi e) (Kabatiansky-Levenshtein) =>
  density(W (x) L) <= 0.178^{N/2} 2^-N poly: exponentially below Minkowski-Hlawka. So "fixed symmetric W tensor random
  L" over Z is dead, whatever the group.
- Over K with root discriminant rd_K: super-stability only needs gamma(W as Z-lattice of dim kd) >= kd/rd_K, which is
  weaker than MH-quality once rd_K > 2 pi e ~ 17.1 (random rank-k W satisfies j=1 if rd_K > (2 pi e)^{k/(k-1)} pi^{-1/(k-1)};
  k=2: rd_K > ~93). So over cyclotomic K of large conductor the tensor family is NOT obviously obstructed.
- Necessary conditions on a non-free element g of prime order p with fixed space of K-dimension j = theta k in a
  super-stable W:  e^{k H(theta)} <= [W : W^g (+) (W^g)^perp]^{2/d} <= min(p^{2j}, p^{2(k-j)/(p-1)}).
  Consequences: reflections are impossible for k >= 2; every non-central involution must have trace 0 (theta = 1/2,
  borderline equality); elements with small fixed spaces need large order. (Clifford groups satisfy the involution
  condition; Barnes-Wall itself fails j=1 since gamma(BW_n) ~ sqrt(n/2).)
- OPEN (candidate research direction): find Hermitian O_K-lattices W (K cyclotomic, large rd_K) with
  |Aut_{O_K}(W)| >> k^2 * m and super-stability; then V >= |Aut| (first moment) with N = k r d, r ~ k.

## 5. Higher moments (question (d))

GSV 2023: for rank t large (t > O(n^3 loglog n)) the number of omega_K-orbits in a ball of volume V is asymptotically
Poisson(V/omega_K) (all moments), uniformly in cyclotomic towers. Consequence (my inference): P(no nonzero point)
~ exp(-V/omega_K) > 0 for every V, so second-moment/Poisson information alone cannot push beyond V ~ omega_K * log(1/P):
with Haar measure, events of probability e^{-V/omega_K} exist but a second-moment (Paley-Zygmund) argument gives nothing
new versus first moment + Markov; one needs to exploit exponentially rare events (Klartag: probability ~ e^{-c n^2}?) or
change the measure. For Z-lattices: Rogers 1955-56 (Poisson(V/2) for pairs), Soedergren-Stroembergsson (CLT quoted in
Gargava-Viazovska footnote: (#-(1+V))/sqrt(2V) -> N(0,1) for log V = o(d)).
