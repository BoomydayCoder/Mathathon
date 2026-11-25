# R3: G-symmetric version of Klartag's stochastic-ellipsoid argument

Tags: [P]=proved here, [C]=conditional, [H]=heuristic, [N]=numerical, [L]=literature, [R]=refuted.
Scripts: `check.py` (identities, D, q), `check_ld.py` (large-deviation threshold). Both run in < 5 s.

## 0. Setting
G finite group acting orthogonally on R^n with character chi, L a G-invariant lattice.
S^G = G-invariant symmetric matrices, inner product Tr(AB), D = dim S^G, P = orthogonal projection Sym_n -> S^G,
P(B) = |G|^{-1} sum_g g B g^T.  Orbit size of generic lattice vectors: w (if G acts freely on R^n minus 0 and -1 in G, w=|G|).
Process: dA_t = pi_t(dW_t), W standard BM in S^G, pi_t = orthogonal projection (inside S^G) onto
F_t = {B in S^G : Bx.x = 0 for all contact points x}. E_{A_t} is G-invariant, so contact sets are unions of G-orbits.
The process freezes when F_t = 0, i.e. at a "G-perfect" form (Berge-Martinet-Sigrist terminology) [L].

## 1. Dimension count and the identity  (task item 1)

**(1a) [P] D = (1/|G|) sum_g (chi(g)^2 + chi(g^2))/2** (character of Sym^2). Checked numerically.

**(1b) [P] Structure.** W Q-irreducible, A = Q-span of G in End(W) = M_r(Delta), F = centre, [F:Q]=f, dim_F Delta = d^2,
dim_Q W = r d^2 f, dim_Q A = r^2 d^2 f. Commutant of G on W^k is M_k(Delta^op); over R it is
prod_v M_{N_v}(D_v), with (D_v,N_v) = (R,kd) at real split places, (H,kd/2) at real ramified places, (C,kd) at complex places.
Invariant symmetric forms = self-adjoint elements, hence
  D = sum_{v real split} kd(kd+1)/2 + sum_{v real ramified} kd(kd-1)/2 + sum_{v complex} (kd)^2
    = k^2 dim_Q(Delta)/2 + O(k d * #real places).
So the orchestrator's D ~ k^2 dim_Q(Delta)/2 is correct to leading order; exact value above.
Examples: cyclotomic (d=1, phi(m)/2 complex places): D = k^2 phi(m)/2 exactly [N-checked m=5,7,8,9,12,15];
Hurwitz (H at the real place, d=2): D = k(2k-1); icosians (two real ramified places): D = 2k(2k-1).

**(1c) [P] Identity.** chi_W = d * sum_{i=1}^f chi_i (Galois orbit of an absolutely irreducible chi_i of degree rd, Schur index d).
sum_g chi_W(g)^2 = |G| <chi_W,chi_W> = |G| d^2 f, chi_W(1)^2 = r^2 d^4 f^2, so
  sum_g rho(g)^2 = |G|/(r^2 d^2 f) = |G| / dim_Q(A).   (rho = chi_W/chi_W(1))
If the representation is fixed-point-free, W restricted to <g> is a multiple of Q(zeta_o), o = ord(g), so rho(g) = mu(o)/phi(o) and
  |G|/dim_Q A = sum_g mu(o(g))^2/phi(o(g))^2   (only elements of squarefree order contribute).
[N] verified: cyclic m (= m/phi(m), all m tested up to 210), Q8 (2), 2T=SL2(3) (6), 2O (6), 2I=SL2(5) (15).

**(1d) [P, correction of H2] Gain over Klartag is g = w D / n^2 -> |G|/(2 dim_Q A), not |G|/dim_Q A**
(with -1 in G, n = k r d^2 f, D ~ k^2 d^2 f/2). The factor 1/2 is forced by the trivial case G={+-1}: 2/1 would give "gain 2".
Cyclotomic, m even: g = m/(2 phi(m)) = prod_{p | m, p odd} p/(p-1); m odd: use mu_{2m}, same value m/phi(m).

## 2. Cyclotomic case in detail (task item 2)
K = Q(zeta_m), m >= 3, s = phi(m)/2 complex places, R^n = (C^k)^s, n = 2ks = k phi(m), G = mu(K), w = |G| = lcm(2,m).
G acts freely on nonzero vectors of K^k (K is a field). S^G = (H_sigma)_sigma, H_sigma in Herm_k(C); |B|^2 = sum_sigma 2 Tr(H_sigma^2); D = s k^2 = nk/2.

### (a) quadratic variation  [P]
M_t = A_t x.x - 1, dM = <dW, pi_t P(x (x) x)>, so [M]_t <= t |P(x(x)x)|^2.
General: |P(x(x)x)|^2 = <P(x(x)x), x(x)x> = |G|^{-1} sum_g <gx,x>^2.  q(x) := |P(x(x)x)|^2/|x|^4.
For x uniform on the sphere: E <gx,x>^2 = (chi(g)^2 + chi(g^2) + n)/(n(n+2)), hence
  **E q = (2D + n)/(n(n+2))**   (=1 for trivial G; check).
Cyclotomic: P(x(x)x) = (x_sigma x_sigma^*/2)_sigma, so
  **q(x) = (1/2) sum_sigma |x_sigma|^4 / |x|^4 = (1/2) sum_sigma u_sigma^2,  u_sigma = |x_sigma|^2/|x|^2,**
  u ~ Dirichlet(k,...,k) on the s-simplex.  Range: 1/phi(m) = k/n <= q <= 1/2.  E q = (k+1)/(n+2) ~ (1+1/k)/phi(m).
[N] direct group average vs formula agree to 1e-16; MC mean agrees with (k+1)/(n+2). The orchestrator's guess is correct.
Interpretation [P]: the excess 1/k = n/(2D) is exactly the Jensen gap of the "breathing" (trace) mode: log det A_t has martingale part
with quadratic variation ~ t|P(Id)|^2 = tn, so E Vol_t = e^{(2D+n)t/8} while the typical volume is e^{Dt/4}. (Removing the trace
direction from W gives E q_0 = 2(D-1)/(n(n+2)); not needed below.)

### (b) operator norm  [P up to standard Gaussian concentration]
Maurey's trick (Lemma 3.1) works verbatim inside S^G: A_t - a_0 Id = (W^1_t + W^2_t)/2 with W^i BMs in S^G.
W_t = (H_sigma(t)), H_sigma(t) = sqrt(t/2) * GUE_k (density exp(-Tr X^2/2)), so E||H_sigma(t)||_op <= sqrt(2kt), and
||.||_op is (1/sqrt2)-Lipschitz for the S^G norm, hence
  P( ||A_t - a_0|| >= sqrt(2kt) + r ) <= 2 s exp(-r^2/t).
With time horizon T = 4z/D = 8z/(nk) (z ~ log(wD) ~ 2 log n): ||A_T - a_0|| <~ sqrt(16 z/n) + sqrt(8 z u/(nk)) -> 0 with prob 1-e^{-u}, u = C log n.
Klartag's cruder event (radius C_0 sqrt(tn), prob 1 - e^{-n}, needed against the kissing bound 2(2^n-1) in Prop 4.2) has
C_0 sqrt(Tn) = C_0 sqrt(8z/k): it is < 1/2 iff k >= C z. For smaller k one must instead stop the process when lambda_min(A_t) <= a_0/2
(routine; P(stop before T) <= 2s exp(-c/T) -> 0).
Error term of Prop 3.4: int_0^T D sqrt((k+log n) t) dt ~ z^{3/2} sqrt((k+log n)/(nk)) -> 0.  OK for all k.

### (c) Ito / dimension count  [P]
delta_t = sum over an orthonormal basis f_b of F_t of Tr(A^{-1} f_b A^{-1} f_b) >= ||A_t||^{-2} dim F_t, so
  E log det A_T <= n log a_0 - (1/2) int_0^T E[ ||A_t||^{-2} N_t ] dt,   N_t = dim F_t.
For G-invariant B, B(gx).(gx) = Bx.x, so each G-orbit of contact points imposes one linear condition:
  **N_t >= D - |dE_t cap L| / w.**

### (d) contact count and final bound
Reflection principle (Prop 4.1 verbatim) [P]: P(x in dE_T) <= 2 P(Z >= (a_0 - |x|^{-2})/sqrt(T q(x^))).
Mean value formula for rank-k O_K-module lattices (Siegel/Weil; Venkatesh 2013) [L]: E sum_{x != 0} f(x) = (V_0/kappa_n) int f, where the
unit ball has volume V_0 * covolume. Free orbits: #(L cap ball minus 0) is a multiple of w, so P(initial ball not free) <= V_0/w: may take V_0 = w/2.
Radial integration (Klartag Lemma 4.3 with t -> t q(theta); exact Gaussian identity s_0 int_0^inf P(Z>=y)e^{s_0 y}dy = e^{s_0^2/2}Phi(s_0) - 1/2):
  E|dE_t cap L| <~ 2 V_0 * Psi(z),  Psi(z) := E_u exp( n^2 t q(u)/8 ) = E_u exp( z s sum_sigma u_sigma^2 ),  z := Dt/4 = nkt/8.
Since s sum u^2 >= 1, Psi(z) >= e^z; E[s sum u^2] = s(k+1)/(ks+1) ~ 1 + 1/k, bulk gives Psi ~ e^{z(1+1/k)}.
Bookkeeping [H]: E log det A_T <= -2 log V_0 - 2 z_T + (4 V_0/(wD)) e^{z_T} (if Psi ~ e^z). Optimum e^{z_T} = wD/(2V_0):
  E log det A_T <= -2 log(wD/2) + 2, i.e.  **V >= (1-o(1)) w D/(2e)**, independent of V_0. (Klartag: w=2, D=n^2/2: V ~ n^2/(2e) [H constant].)
With Psi ~ e^{z(1+1/k)}: V >= c V_0^{1/(k+1)} (wD)^{k/(k+1)} : loss factor (D)^{1/(k+1)} = "contacts ~ V^{1+1/k}"; bounded iff k >= c log n.

**Large deviations of q (new point, not in H2).** Sparse vectors (a fraction alpha of |x|^2 at ONE place) have probability ~ (1-alpha)^{k(s-1)}
and weight e^{z s alpha^2}: they dominate Psi unless z alpha^2 + k log(1-alpha) < 0 for all alpha in (0,1), i.e. unless
  **k > 0.407 z** (max of alpha^2/(-log(1-alpha)) is 0.407 at alpha ~ 0.715).
[N, rigorous lower bound] Psi(z)/e^z >= J := E_{u~Beta(k,k(s-1))} exp(zs(u^2+(1-u)^2/(s-1)) - z). `check_ld.py`, z = log(nk):
 s=100: k=2 (z/k=3.34): log J = 112;  k=4 (z/k=2.02): log J = 0.02.   s=400: k=2: log J = 822; k=4: 0.006.  s=20: k=2 (z/k=2.54): 6.0.
So for k=2 (Venkatesh's rank) the method with isotropic noise genuinely breaks: the reflection bound is essentially tight for such x
(early on pi_t ~ Id), the s blocks H_sigma fluctuate independently, and a vector living on one place sees variance s times larger. [H that it is not an artifact]
For k >= C log n everything is harmless. (My first heuristic used exponent (k-1)(s-1); the correct marginal exponent is k(s-1)-1. Corrected by the numerics.)

**Outcome [C]:** for n = k phi(m), k >= C log n:  V_n >= c * w D = c n^2 w/(2 phi(m)) >= c n^2 m/(2 phi(m)).
Optimal choice for given n: maximize prod_{p|m, p odd} p/(p-1) over m with phi(m) | n and n/phi(m) >= C log n. m = product of primes <= y gives
m/phi(m) ~ e^gamma log log n for every k in [C log n, n^{o(1)}] (indeed any k <= n^{1-eps}): V >= c n^2 log log n on these dimension sequences.
Heuristic constant (k >> log n, m even): V ~ (n^2/(2e)) * m/(2 phi(m)) ~ (e^{gamma-1}/4) n^2 log log n ~ 0.16 n^2 log log n.
NOT for all n: n odd has no m; n = 2^a gives only Fermat primes (gain <= ~2); n = 2p gives gain <= 3/2 (m = 6). Arbitrary n needs a new idea [open].

## 3. Target theorem and lemma list (task item 3)
**Theorem (target) [C].** There are universal c, C > 0: for m >= 3, k >= C log(k phi(m)), n = k phi(m), there is a lattice in R^n (a rank-k
Z[zeta_m]-module with a mu_m-invariant Euclidean structure) of covolume 1 with no nonzero point in a ball of volume c n^2 m/phi(m).
Hence V_n >= c' n^2 log log n for infinitely many n.
Lemmas:
 L1 process existence in S^G (Klartag Lemma 2.1, 2.2, Prop 2.3, Cor 2.4 with Sym_n replaced by S^G) - routine, NOT rechecked line by line here.
 L2 Maurey + operator norm (2b) - routine.   L3 Ito inequality (2c) - routine [P above].   L4 N_t >= D - #contacts/w - [P].
 L5 reflection bound with q(x) - [P].   L6 Psi(z) <= C e^z for z <= c k (Dirichlet exponential moment) - elementary but new; only sketched
    (bulk: Gamma representation; tails: exponent z alpha^2 + k log(1-alpha) < 0); this is the one step without a full written proof besides L1.
 L7 mean value formula for O_K-modules of rank k >= 2 - [L] (Venkatesh).   L8 confinement (kissing bound argument works verbatim once k >= Cz).
No step seems to need a new idea when k >= C log n. Presumably this is what the reported "Klartag x log log n" paper does (not checked by me).

## 4. Division algebras / other groups (task item 4, bookkeeping only)
[P] gain over Klartag g = wD/n^2 -> |G|/(2 dim_Q A) = (1/2) sum_g mu(o(g))^2/phi(o(g))^2 for a fixed-point-free G (then all orbits have size |G|).
 - mu_m: m/(2 phi(m)) (m even).  - 2T in Hurwitz (n=4k): 24(2k^2-k)/(16k^2) -> 3.  - 2O (n=8k): 3.  - 2I in icosians (n=8k): 120(4k^2-2k)/(64k^2) -> 7.5.
 - 2T x mu_m, gcd(m,6)=1, ord_m(2) odd (Hurwitz (x) Q(zeta_m) is then a division algebra): 3 m/phi(m), but such m only use primes with ord_p(2) odd
   (density 7/24), so growth ~ (log log n)^{7/24}: asymptotically worse than cyclic.
 - 2I x mu_m: character sum gives 7.5 m/phi(m) (4x the cyclic value at equal prime support) BUT icosian algebra (x) Q(zeta_m) = M_2(.) splits (r=2).
 [H, obstruction] For r > 1 and modules over a maximal order, L = (O^r, C) (x)_O (M, H) is a tensor product lattice; pure tensors c (x) M give
 ~ (r/min C)^{n/(2r)} lattice points in the unit-volume ball (exponentially many since min C < r by Hermite). So r > 1 needs non-maximal Z[G]-lattices
 and a mean value formula with rank-deficient terms: new idea needed. Within division algebras (Amitsur's list) the gain is O(log log |G|) [conjecture],
 only elements of squarefree order count; best constants: 2I (7.5), cyclic for growth.

## 5. Addenda (after numerics)
**(5a) [N] Bulk value of Psi.** `check_psi.py` (MC, 1e5 samples, z = log(nk)): Psi(z)/e^z = 2.508 (s=50,k=10; e^{z/k}=2.512), 1.690 (50,20; 1.699),
1.344 (50,40; 1.349), 2.888 (200,10; 2.885), 1.395 (200,40; 1.397), 1.190 (100,80; 1.192). So Psi(z) = e^{z(1+1/k)(1+o(1))} in the regime k >~ z:
the contact count is V^{1+1/k} as guessed, a bounded factor e^{z/k} once k >= c log n. (MC cannot see the rare sparse vectors; those are in check_ld.py.)

**(5b) [P-sketch] Elementary version of L6 for k >= C log^3 n.** sum u_sigma^2 <= max_sigma u_sigma, so Psi(z) <= E exp(z s max u_sigma)
<= e^{z(1+delta)} + sum_sigma E[e^{z s u_sigma}; s u_sigma > 1+delta]. u_sigma ~ Beta(k,k(s-1)), density of v = s u_sigma is <= C sqrt(k) e^{-k(v-1-log v)}.
With delta = 1/z the exponent zv - k(v-1-log v) is decreasing on v >= 1+delta once k >= z(z+1), and the tail sum is
<= C s sqrt(k) e^{z+1} e^{-k/(3z^2)} <= e^{z} once k >= C z^2 log(s k). Hence Psi(z) <= C e^z for k >= C z^2 log n ~ C log^3 n.
Since the conclusion c n^2 m/phi(m) does not depend on k (any k <= n^{1-eps} still gives log log m ~ log log n), the target theorem can be
stated with k >= C log^3 n and then L6 is fully elementary; the sharp range k >= C log n needs the finer large-deviation estimate (sketched only).

**(5c) [L]** Primes with ord_p(2) odd have density 7/24 (Hasse 1966: primes dividing some 2^j+1 have density 17/24). Used in section 4.

**(5d) Status summary.** Proved here: 1a-1d, 2a (q formula, E q), 2c, reflection bound, bookkeeping identities, necessity of k > 0.407 z for
Psi(z) <= C e^z (via the lower bound J; the constant 0.407 from the exact Beta exponent, confirmed numerically).
Conditional: the target theorem (needs L1 re-verified in S^G, L6 written in full, L7 quoted in the exact normalisation, L8).
Heuristic: constants (V ~ wD/(2e)); obstruction for r > 1; "sparse vectors are a real effect, not an artifact of the reflection bound".
Open: all n; whether anisotropic noise (non-identity covariance on S^G) or a non-Euclidean reference metric removes the k >= c log n requirement
(irrelevant for the n^2 log log n statement, relevant for small-rank modules such as Venkatesh's k = 2).
