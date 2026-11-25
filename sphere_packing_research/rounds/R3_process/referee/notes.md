# Referee notes on R3_process (adversarial verification), 2026-11-25, ~14 min
Scripts: checks.py, primecut.py (same directory). Tags as in the project.

## A. Claims tagged proved
1. Lemma N [confirmed]. Re-derived every step. (c): ||h-1|| <= 2eps/(1-eps) < 1 iff eps<1/3, h fixes each minimal vector, minimal vectors span (stretching
   argument is correct: need 1+eta >= e^{ns}). (e): 2delta/(1-delta)=eps for delta=eps/(2+eps); delta/(1+delta)=eps/(2+2eps). (f) arithmetic ok (0.690 n^2+1<=n^2 for n>=2).
   Siegel constant: for n=2 I computed the cone volume of the reduced-basis fundamental domain by hand:
   pi * int_{-1/2}^{1/2} (1/4)(1-s^2)^{-1/2} ds = pi^2/12 = zeta(2)/2  (checks.py: 0.82247; brute-force MC 0.809 + truncated tail 0.0145). General n: literature.
   Exact bound reproduced: n=8: -269.4, n=24: -3655.6, n=100: -91828.2.
   Remark on exponent 1 instead of 2: plausible [heuristic]: polar decomposition, vol(O(n)) ~ exp(-(n^2/4)log n), Sym-part (c eps/sqrt n)^D with eps~1/n gives
   -(3/4) n^2 log n, total -(1+o(1)) n^2 log n; non-injectivity multiplicity <= (5^n)^n = e^{O(n^2)} (h is determined by images of n minimal vectors among vectors of length <= 2 lambda_1).
2. Corollary N' [confirmed]; the ">> e^{-V/2} for V >= 5 n^2 log n" clause needs n >= 8 (1 < 0.5 log n); "equivalently" is loose.
3. Prop U [confirmed]. Net size, A' >= A, A' <= (1+1/n)A, det A' <= e det A, SL_n-invariance all fine. Poisson first-moment threshold of the statement as written is
   sqrt(e) n^2 log(4 rho n)(1+o(1)) (constant sqrt e is a net artefact: eta = delta(1-rho)/n gives e^{delta/2}).
4. Spherical identity [confirmed]: hand derivation + MC with E = span of c random contact tensors, n=7 (c=1,5,12,20,27): MC vs formula agree to 3-4 digits.
   Caveats: |P_E Id|^2 ~ cn/(n-1+c) is poor near c ~ D (n=7,c=27: true 6.93, approx 5.73; limit n is right);
   exact: E|pi(w(x)w)|^2 = (2N + n - |P_E Id|^2)/(n(n+2)), slightly BELOW N/D (0.2632 vs 0.2857 at n=7, c=20).

## B. Conditional / heuristic claims
5. Mean-field numbers reproduce: closed form c(u) vs quadrature (n=1e5,u=2: 6.7224 vs 6.7210), 0.1763, 1/(2e), 1/3, 0.212. Rezvyakova arXiv:2608.22590 abstract: c = 1/(2e)-o(1) [literature, checked].
   Model R (constant 1/2) is NOT supported by the R2 data: at n=9,10 V_final/pairs = 0.747, 0.746 (model: 1), pairs/D = 1.35, 1.38 (model: 1), model predicts V/D=1.67 at n=10 vs 1.03.
6. Sticky barrier in the Poisson model: the chain c_t = N(R_t \ E_0), c_t <= D, Lambda_T < S_{D+1}, V_T <= relvol(R_T) is logically fine GIVEN the stopping-set time-change
   (not written) and general position (a.s. for Poisson). Tail inequality verified (D=55, lam=2: 4.95e-9 <= 4.68e-8; D=210).
   GAP for lattices: general position fails. R2 table: frac(pairs = D) = 0.000 at n=9 and n=10 (min 46 > 45, min 58 > 55; max 120), mean pairs at T_K = 54.2 > D = 45 with only
   46% of runs frozen. Reason: the process reaches a face boundary with probability 1, and faces of the Ryshkov polyhedron of a lattice are typically NOT simple at these n.
   Correct lattice statement [proved, trivial]: V_T <= relvol(R_T) and #lattice pairs in R_T \ E_0 = c_T; c_T is bounded only by the half kissing number of the final form.
   Also: in the Poisson model sup over all ellipsoids is +infinity (needles), so 2.0 is about the process, not about existence.
7. Vertex law Gamma(D,2): Mecke + scaling exponent re-derived (dx = s^{nD-1} ds dsigma, V = s^n => V^{D-1} dV). Root: n=100, kappa=1: lam=7.65, V=1.68 n^2 log n. ok (per shape only).
8. Modifications: (i) ratio n(n+2)/8 ~ D/4 verified. (ii) normal velocity mean (n-1)/(n-1+c), s.d. ~ sqrt(2c)/(n+c): verified from the identity.
   Added micro-lemma [proved]: once the contact vectors span R^n (c >= n) there is NO monotone expansion direction: B in F_t, B <= 0, B x_i.x_i = 0 => x_i in ker B => B = 0.
   (v) "dependent contacts have probability zero generically" is contradicted for real lattices at n <= 10 by the data cited in the same notes (see 6). Open whether pairs/D -> 1.
9. Symmetric formula q = n^2/(2 dim_Q B): checked for complex type (k^2 phi/2), quaternionic type (2k^2-k = n^2/8 - n/4, n=4k), real type (k(k+1)/2 with dim B = d^2).
   Consistency: Vance/Rogers = (24/4)/(2/1) = 3 ok; rank-2 cyclotomic: q = n, V ~ nm = (n^2/2) m/phi(m) ok (Markov alone gives V ~ m: Venkatesh); GGS intro states density of
   order 2^{-d} d^2 M/C ok. For G={+-1} the formula gives n^2, i.e. it is the "typical vertex" normalisation, a factor 2 above model R.

## C. GGS multi-top speculation (R3 section 4.B)
GGS eq. (11): A(S) = (c_S/2) log(C/c_S) - (1/2) W(S,S^c) - sum_{n in S} c_n v_n, sum_N c_n v_n = 0, w(m,n) = log|Res(Phi_m,Phi_n)|; Def 2.6 needs 1 in N, 2 not in N, c_n = phi(n).
(a) R3's isoperimetric condition is CORRECT: A(S)+A(S^c) <= -2 delta  <=>  W(S,S^c) >= (C/2) H(c_S/C) + 2 delta.
(b) Prime-cut lemma [proved given (11) and |Res(Phi_m,Phi_{m p^k})| = p^{phi(m)}, else 1 (spot-checked with sympy)]: for a prime p | M let S_p = {n in N : p | n}
    (nonempty, proper since 1 in N), x_p = c_{S_p}/C. Each n in S_p has at most one neighbour outside S_p, namely n/p^{v_p(n)}, with weight
    phi(n) log p /(p^{k-1}(p-1)) <= c_n log p/(p-1). Hence W <= x_p C log p/(p-1) and  H(x_p) < 2 x_p log p/(p-1).
    Since H(x)/x is decreasing this means x_p > x*(p): x*(3)=0.609, x*(5)=0.728, x*(7)=0.791, x*(11)=0.857, x*(101)=0.982, 1-x*(p) ~ 1.9/p.
    Corollary: any two primes p,q | M divide a common index n in N (x_p + x_q > 1); any set T of primes with sum_{p in T}(1-x*(p)) < 1 divides a common index.
    => R3's example "tops Pg and Qg, no index divisible by PQ, M/C ~ (g/phi(g)) PQ/(P+Q)" is REFUTED (inadmissible for all odd primes P,Q).
(c) Not excluded [heuristic]: a balanced fan of K >= p0/2 tops M/(p_i q_i) (p_i,q_i ~ p0 distinct): prime cuts ok; the cut "indices exclusive to top 1" needs
    M (log p0)/p0^3 >~ phi(M) (log p0)/(2 p0^2), i.e. p0 <~ 2 M/phi(M); then C ~ phi(M)/p0 and M/C ~ (M/phi(M))^2 ~ (log log d)^2. Full admissibility (all S, potentials v,
    I/C << log C, r >> C^3 log d) NOT checked. This is the concrete finite problem worth an agent.
