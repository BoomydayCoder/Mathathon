# Referee notes on R3_G_klartag (adversarial verification)

Scripts (all < 1 min): `ref_dims.py` (independent D, q, E q), `ref_trunc.py` (1-coordinate sparse family, with/without cutoff),
`ref_full2.py` (rigorous upper bound for the truncated Dirichlet moment, full Dirichlet).
Tags: [P] proved here, [C] conditional, [N] numerical, [H] heuristic, [R] refuted.

## A. Claims that survive

A1 [P, re-derived]. sum_g rho(g)^2 = |G|/dim_Q A. chi_W = d * (f distinct abs. irreducible chi_i of degree rd), <chi_W,chi_W> = d^2 f,
chi_W(1) = r d^2 f. Fixed-point-free: all eigenvalues of g are primitive o(g)-th roots (else a nontrivial power has eigenvalue 1),
rationality gives chi(g) = dim W * mu(o)/phi(o). Hand check: 2T: 1+1+8/4+8/4 = 6; 2I: 1+1+20/4+20/4+24/16+24/16 = 15; Q8: 2. OK.

A2 [P + N]. D formula. Independent computation (nullspace of B -> gBg^T - B on Sym_n, explicit matrices, `ref_dims.py`):
mu_5,k=2: 8; mu_7,k=2: 12; mu_12,k=3: 18; mu_9,k=2: 12; mu_15,k=1: 4; mu_16,k=2: 16 (= k^2 phi(m)/2 in all cases);
Hurwitz k=1,2,3: 1, 6, 15 = k(2k-1); icosians (rational rep = rep + Galois conjugate) k=1,2: 2, 12 = 2k(2k-1). All match.

A3 [P + N]. q(x) = |P(x(x)x)|^2/|x|^4, E q = (2D+n)/(n(n+2)); cyclotomic q = (1/2) sum u_sigma^2 (asserted to 1e-10 on random x);
MC of E q matches in all cases above. EXTRA OBSERVATION [N, and P by the block formula]: for Hurwitz modules q is CONSTANT = 1/4
on the sphere (std 0 in MC, all k): a single archimedean place => no direction dependence at all. Direction dependence of q
(and hence any large-deviation lemma) only arises with several places (cyclotomic s>1, icosians s=2).

A4 [P]. gain = wD/n^2 -> |G|/(2 dim_Q A); Hurwitz 3-3/(2k), icosians 7.5-3.75/k, mu_m (m even) m/(2phi(m)). Bookkeeping only.
Consistency: G-Klartag/Vance = (24k^2/e)/(24k/e) = k = n/4 and Klartag/Rogers = (n^2/2e)/(2n/e) = n/4: same ratio. Trivial G gives Klartag.
At t=0 the Markov/orbit count gives V_0 < w = Venkatesh. No inconsistency found.

A5 [P]. Arithmetic of dimension sequences (odd n; n=2^a => Fermat primes, product = 2(1-2^-32); n=2p with k>=3 => m in {3,4,6}).

## B. Claim REFUTED as stated: "sparse vectors break the method unless k > 0.407 z; at k=2 the method genuinely breaks"

The mathematical statement about Psi(z) = E exp(z s sum u^2) (unrestricted) is true (their J lower bound reproduces: e.g. s=400,k=2: log J = 822).
But Psi is NOT the quantity the proof needs. In Klartag's Prop 4.2 / Lemma 4.3 the lattice sum runs only over the shell
D_t = {x : (a_0 - C_0 sqrt(tn)) |x|^2 < 1}, i.e. over points that can be reached at all while lambda_min(A_t) stays above a_0 - rho_max.
The unrestricted Gaussian integral e^{n^2 t q/8} has its saddle at rho* = n t q/2 = 4 z q/k; for a sparse direction (q ~ alpha^2/2) this is
rho* ~ 2 z alpha^2/k = O(1), far outside the shell (rho_max = sqrt(2kT) + sqrt(T u_0) ~ sqrt(z log n/n) -> 0, T = 8z/(nk)).
Correct per-direction bound:  I(theta) <= min( 2 e^{n^2 t q(theta)/8}, e^M ),  M = -(n/2) log(1 - rho_max/a_0) ~ (n/2) rho_max.

Stopping-time version (replaces Klartag's kissing-number step, so NO failure probability e^{-n} and no condition k >= Cz is needed) [P-sketch]:
tau = inf{t : ||A_t - a_0|| >= rho_max}. P(tau < T) <= 8 s e^{-u_0} (Maurey + Gaussian concentration + Levy maximal inequality).
Ito inequality for the stopped process: E log det A_{T^tau} <= n log a_0 - (1/2) E int_0^{T^tau} ||A_t||^{-2} N_t dt,
N_t >= D - #contacts_t / w, and  E #contacts_{t^tau} <= sum_{x in shell} P(x touched by time t) <= sum_{x in shell} 2 P(Z >= rho(x)/sqrt(t q(x)))
(x touched at a time <= tau forces rho(x) <= rho_max). int_0^{T^tau} N_t dt >= int_0^T (D - #contacts_{t^tau}/w) dt - D T 1{tau<T}.

Truncated moment lemma L6' [P-sketch + N]: Q := E_u min(e^{z s sum u^2}, e^M) <= e^{z(1+1/k)(1+o(1))} for EVERY k >= 1, provided z*rho_max -> 0.
Proof: u = g/S, g_sigma iid Gamma(k), S = sum g ~ Gamma(ks). With eps, g_0:
  Q <= e^M P(S < ks(1-eps)) + e^M s P(g > g_0) + ( E exp(c' min(g,g_0)^2) )^s,   c' = z/(s k^2 (1-eps)^2),
since on the complement z s sum u^2 <= c' sum min(g_sigma,g_0)^2. Using e^y <= 1 + y e^y, y <= c' g_0 g:
  E exp(c' min(g,g_0)^2) <= 1 + c' E[g^2 e^{c' g_0 g}] = 1 + c' k(k+1)(1 - c' g_0)^{-(k+2)},
and s c' k(k+1) = z(1+1/k)/(1-eps)^2. Take g_0 = M + log s + k log M + C (then e^M s P(g>g_0) <= 1), c' g_0 ~ z rho_max/k -> 0,
eps ~ sqrt(4M/n) ~ sqrt(2 rho_max) (then e^M P(S<ks(1-eps)) <= 1); loss exp(O(z eps + z^2 rho_max)) -> 1 (slowly: needs rho_max << z^-2).
Numerics (`ref_full2.py`, the right-hand side evaluated with mpmath): log(bound) - z(1+1/k) =
  k=2: 17.5 (s=1e6), 9.8 (1e7), 5.9 (1e8), 3.6 (1e9);  k=4: 10.8 (1e6), 3.8 (1e8);  k=16: 5.7 (1e6), 2.1 (1e8).
Bounded and decreasing, as predicted (slow convergence is due to the crude eps; a joint large-deviation treatment of S would sharpen it).
1-coordinate family with cutoff (`ref_trunc.py`): k=3: log J drops from 3840 to 0.002 (s=2000); k=2: 45407 -> 30.6 (s=1e4), still shrinking.
So the sparse-vector term is an ARTIFACT of summing the reflection bound beyond the confinement radius. It is asymptotically harmless for every k
(it does bite at moderate n for k=2,3: the cutoff needs rho_max < k/z, i.e. n k^3 >~ z^3 log n with poor constants).
Why the union bound overcounts: all far sparse points at place sigma_1 are touched only on the single event {||H_{sigma_1}(T)|| >= rho}.

## C. Improvement the notes mention but discard: trace-free noise removes the V^{1+1/k} bulk loss [C]

Run the Brownian motion in S^G_0 = {B in S^G : Tr B = 0}. Then [M]_t <= t q_0(x)|x|^4 with q_0 = q - 1/n EXACTLY (Id/sqrt n is a unit vector of S^G),
N_t >= (D-1) - #contacts/w. Exponent: n^2 t q_0/8 = z s sum u^2 - z/k (z = nkt/8), so with L6':
  E #contacts <= 2 V_0 e^{z(1+o(1))}   for every k >= 2,
instead of 2 V_0 e^{z(1+1/k)}. (The notes computed E q_0 = 2(D-1)/(n(n+2)) and wrote "not needed below"; it is exactly what kills the loss.)
Operator norm: ||W^0|| <= ||W|| + |Tr W|/n, extra term ~ sqrt(t/n), negligible. Hurwitz: q_0 n^2/(2 D_0) = 2k/(2k+1) < 1 (subcritical, like Klartag's n/(n+2)).
Consequence [C]: the restriction k >= C log n in the target theorem is not needed; conjecturally the method gives V >= c n^2 w/phi(m) for all
n = k phi(m), k >= 2 (including Venkatesh's n = 2 phi(m)), n large. Not written out line by line; same status as the notes' theorem.
Without trace removal and with the cutoff, k=2 gives V ~ V_0^{1/3} (wD)^{2/3} ~ m n^{2/3}: polynomially better than Venkatesh, NOT a breakdown.

## D. Remaining gaps in the notes' "conditional theorem" (unchanged by my checks)
- L1 (process in S^G, Klartag Sec. 2) not rechecked by them or by me. Looks routine: the L-free forms in S^G form a locally finite polyhedral set.
- L7 mean value formula: fine for free modules O_K^k, k >= 2, Haar on SL_k(K_R)/SL_k(O_K), ALL nonzero vectors (not primitive): lattice vectors lie
  in the open orbit {x_sigma != 0 for all sigma}. Not re-derived; normalisation taken from Venkatesh. k = 1 excluded (no such formula; norm form bounded below).
- E||GUE_k|| <= 2 sqrt k is only true up to a lower-order correction O(k^{-1/6}); irrelevant.
- Heuristic constants n^2/(2e), 0.16 n^2 log log n: consistent bookkeeping (I re-derived -2 z_T + (4V_0/(wD)) e^{z_T}, optimum wD/(2e)); not proved.
- r > 1 obstruction: the count (r/C(c))^{k/2} for the sublattice c (x) M re-derived for F = Q; fine as a heuristic.
- (5b) elementary L6: steps check out (sum u^2 <= max u; monotonicity needs k >= z(z+1); density bound needs s >> z^2). Superseded by L6'.

## E. Cleanest fully-checkable statement
Hurwitz modules: s = 1 place, q == 1/4 constant, so Klartag's proof goes through verbatim with (n^2/8, n(n+1)/2, 2) replaced by
(n^2 q/8 = n^2/32, D = k(2k-1), w = 24) and Vance's mean value formula: no large-deviation lemma, no stopping time beyond Klartag's own.
Heuristic outcome V ~ 24 k(2k-1)/(2e) ~ 3 n^2/(2e) for 4 | n (3x the heuristic Klartag constant). Icosians (two places, 8 | n): 7.5x, needs L6' with s = 2
(u ~ Beta, trivial).
