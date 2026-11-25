# R11 / T2 - devil's advocate on the "dichotomy conjecture" (SL_2(p) tensor family, V ~ p n^2 ~ n^{2.5})

Start 08:26 UTC, end ~08:50. Scripts: `t2_checks.py` -> `t2_out.txt`, `mc_q.py` -> `mc_q_out.txt`.
Tags: [P] proved here (short complete argument), [PS] proof sketch, [H] heuristic, [N] numerical, [L-mem] literature from memory, [O] open.
Notation: r = (p-1)/4, n = 8rk, w = |G| = p(p^2-1), D = D_k = 4k^2-2k, q_0 = 1/dim_Q A = 1/(8r^2), z = Dt/4.

## 0. Verdict in one paragraph

No FATAL structural obstruction found. But the claim AS STATED (k ~ r+1, n ~ p^2/2, V >~ p n^2) is NOT what the equivariant
Klartag argument (in the form of Klartag/Rezvyakova/AGZ/GGS/our R6) gives: in a tensor family the ratio n/(2D) equals r/k, and
this is exactly the relative size of the "breathing" (trace) mode of the noise. For k ~ r the argument yields
V ~ sqrt(V_0 w D) ~ 0.4 n^2, i.e. exponent 2.0, no gain at all (A3/A4). The gain p survives only when k/r -> infinity
(loss-free for k >= C r log n, i.e. n ~ C p^2 log p, V ~ n^{2.5}/sqrt(log n)), or if a trace-free variant of the process is proved.
Further: (i) new families of exactly-borderline or possibly unstable Delta-subspaces (isotypic components of SL_2(3), 2O, 2I
subgroups) - a new explicit necessary condition for p = 61, 421, 541, 661; (ii) "borderline for all Delta-subspaces" does not
by itself give E# = V(1+o(1)); a polynomial bound on the number of effectively-borderline subspaces is an extra hypothesis.

## A1. Main-term constant / zeta function of the O-module: NO OBSTRUCTION [P modulo Weil/Rogers]

The mean value formula (R3 Thm 1) is an exact identity: nonzero X in W_Z (x)_O L = r-tuples in L; the full-rank part is the sum
over ALL linearly independent r-tuples of L (not primitive ones), whose Siegel-Rogers-Weil mean is the integral with constant
EXACTLY 1 (needs r <= k-1). Decomposing the same term by the submodule M generated gives sum_M [W_Z:M]^{-k} x (density of
k-tuples generating M), and the two factors cancel. The Solomon zeta function itself [P, Morita at every finite place since
Delta is unramified there: O_v = M_2(O_{F,v})]:
  zeta_{O^r}(k) = sum_M [O^r:M]^{-k} = prod_{i=0}^{2r-1} zeta_F(2k-i),  convergent iff k >= r+1 (same threshold as the height zeta).
[N] k=r+1: 1.144 (p=13), 1.153 (37), 1.168 (61), 1.147 (157); k=r+2: 1.014; k=2r: 1.000. Uniform bound prod_{j>=3} zeta(j)^2 = 1.946.
So it is bounded uniformly in p and it is NOT a multiplier of the main term (c = 1). The L-side sublattices W_Z (x) lO and
W_Z (x) L' are already contained in the identity (rank-one tensors = points of Lambda_U for lines U); no extra terms.
Height-zeta tail for lines at k = r+1 under Schmidt asymptotics #{U: N(n_U) <= N} ~ C_r N^{2r}:
sum_U (r^2/N_U)^{2k} ~ C_r r^{4r} r/(k-r) = (generic number of lines below the borderline) x r/(k-r); by R5 2.8(a) this is
~ (pi e/sqrt p)^{p-1}: negligible for p >= 89 [H]. Caveat: R_U^{k/2} = (r^2/N)^{2k} does NOT suppress near-borderline lines
when k ~ r: all lines with N(n_U) <= r^2 + r^2/(2k) (~ r/2 values of the norm) count like borderline ones.

## A2. Number of borderline subspaces vs D_k: NO OBSTRUCTION from the forced ones; hypothesis needed for the others

Key invariance [P]: for every member of the family and every invariant ellipsoid E of relative volume V,
Vol(E cap U_R (x) L_R)/covol(Lambda_U) = T_U(V) (because det of the form on Lambda_U = det^f). So T_U(V) is the relevant count
during the whole process, and the final lattice needs V(Lambda_U) >= T_U(V_final) for each borderline U (deterministic).
Contact orbits from the G-class of a borderline U: ~ T_U(wD/2) |H_U|/|N_G(U)| (H_U pointwise, N_G(U) setwise stabiliser).
T_U(V) = f^{-1/2} pi^{(f-1)/2} V^f n^{-(1-f)/2} ~ n^{(5f-1)/2} p^f at V = wD/2.
- C_h-fixed spaces of split tori (Lemma B, f = 1/h, h | r odd; h = 3 always exists for p > 5): constraints ~ 1.1-1.7 per h [N].
- NEW: psi_0-isotypic component W_0(S) of each subgroup S = SL_2(3) (p^3/24 subspaces, f = 1/3): R = 2^{-e} (see A6);
  if e = 0 it is exactly borderline; constraints ~ T/24 = 0.4 (p=13) ... 31 (p=1009), versus D ~ p^2/4. Negligible.
- Borderline Delta-lines: T_line -> 1/sqrt(8 pi k); they fill D only if their number is >~ |G| D sqrt(8 pi k):
  1.2e4 (p=13,k=4; Leech has 1638), 3e7 (p=37), 1e9 (p=61), 4e11 (p=157), ~ p^{5.5}. No a priori bound known [O].
- A borderline class with f >= 1/2 and small normalizer WOULD matter (f=1/2: ~p^2/3 orbits ~ D; f=2/3: p^3 >> D, the process
  would stall at V ~ p^{3.5} < n^2). None is forced by the subgroup structure (checked: cyclic, dicyclic, Borel, SL_2(3), 2O, 2I).
Initial ball: P(point of Lambda_U in B_{V_0}) <= T_U(V_0) |H_U|/|N_G(U)|. SL_2(3)-type (if e=0) limits V_0 <~ 131 n, which is
below |G|/2 for p >~ 150. Borderline lines with trivial stabiliser limit (#G-orbits of such lines) << sqrt(k).
At V = |G|/4 the lower-rank terms are NOT o(V): sum over the p^3/24 SL_2(3)-subspaces of T_U / V = 0.35 (p=13) ... 1.6 (p=1009).

## A3. Projected variance: QUANTITATIVE RESTRICTION k >> r (as the argument stands); low-rank points harmless

[P] x = (x_v) in M_{k x r}(H)^2, S_v = x_v^* x_v with eigenvalues lambda_{v,i} (2r in total):
  q(x)/q_0 = 2r sum lambda^2/(sum lambda)^2 = 1 + CV^2(lambda).
[P] E q/q_0 = (2D+n)/(n(n+2)) / q_0 = 1 + (r-1/2)/k + O(1/(rk)); [N, mc_q.py] 1.587 (r=3,k=4), 1.845 (9,10), 1.907 (15,16),
1.961 (39,40); std ~ 0.65/r, minimum over 400 samples 1.92 for (39,40): q is CONCENTRATED at q_0(1+r/k), not at q_0.
[P] The Id-direction contributes exactly 1/n = q_0 r/k to q(x) for every x: the excess is the breathing mode, n/(2D) = r/k.
Consequence in the published bookkeeping (contacts ~ 2 V_0 E exp(z q/q_0) = V_0 e^{z(1+r/k)}):
  V_final = V_0^{r/(k+r)} (wD/2e)^{k/(k+r)}.
  k = r:   V = sqrt(V_0 wD/2e) ~ 0.43 n^2 (V_0 = w/4): exponent 2.0, NO GAIN.
  k = c r: exponent 2 + (c-1)/(2(c+1)) if V_0 ~ |G|; 2 + (c-2)/(2(c+1)) if V_0 ~ n. (c=3: 2.25 / 2.125; c=10: 2.41 / 2.36.)
  loss-free iff z r/k = O(1), i.e. k >= C r log n: n = 8rk ~ (C/2) p^2 log p, V ~ p n^2 ~ n^{2.5}/sqrt(C log n).
Same condition from Klartag's confinement step (needs C_0 sqrt(Tn) < 1/2, Tn = 8 z r/k) and from R6 Lemma 7.2 (kappa < 1).
[H] Removable?: noise restricted to trace-free invariant forms keeps the drift (D-1 directions) and has q_traceless/q_0 =
1 + CV^2 - r/k ~ 1 - 1/(2k) for typical x; then only fluctuations (std ~ 0.65/r, harmless when r >> z) and large deviations
(cost (1-alpha)^{~4kr} vs gain e^{2 r z alpha^2}: fine iff k >~ 0.2 z) remain. This variant is not in Klartag, Rezvyakova, AGZ, GGS
or R6 and must be proved (time-dependent gradient direction A_t^{-1}; alternatively truncate the breathing excursion).
Low-rank points [PS]: on Lambda_U, q/q_0 >= r/s = 1/f, so the Gaussian window is 1/f wider, but the density is d(V^f):
int f e^{fu - u^2/(4 z kappa)} du ~ e^{f^2 z kappa} = e^{f z (1+CV_U^2)}: contacts from Lambda_U = T_U(V_0 e^{z(1+s/k)}), i.e. the
same "volume reach" as the main term: no explosion; rank-1 points on borderline lines: e^{z/r(1+1/k)} = 1+o(1).

## A4. Uniformity in p: gaps, none fatal

(1) n/D = 2r/k = Theta(1) at k ~ r: breathing mode (A3), confinement event, angular lemma: all need k >> r log n or a new variant.
(2) Class number h(O) >= 2 for p >= 37 (R5): W_Z not free; mean value formula per Steinitz class [O, technical].
(3) Start: V_0 limited by borderline classes (A2), not only by |G|/2; irrelevant if (1) is solved (V_final independent of V_0),
    costly otherwise (factor V_0^{r/(k+r)}).
(4) Operator norm: T = 4z/D, GSE_k blocks scaled by (4r)^{-1/2}: ||A_T - a_0|| ~ sqrt(z/(rk)) -> 0. Fine.
(5) z ~ log(wD/2V_0) ~ 5 log p; fluctuation term z x std(q/q_0) ~ 3 log p / r: fine for p >= 61, not small for p = 13 (r = 3)
    unless k >> z (so Leech needs k >= C log n anyway, as in R5).

## A5. Upper bounds / no-go: NO OBSTRUCTION FOUND (not exhaustive)

Kissing number ~ wD ~ n^{2.5}: no contradiction. Deterministic requirement V(Lambda_U) >= T_U(V_final): for f = 1/3 this is
~ 0.34 p in dimension n/3 (far below the n^2-type capacity of that family); for lines ~ 1/sqrt(8 pi k). Minkowski/KL on Lambda_U
gives only exponential bounds (R3 Prop 3a) since f tau^2 = 1. Theta series / Hilbert modular constraints: not examined (time).

## A6. Sanity check and "disease in disguise"

Delta = Q (Barnes-Wall, Weyl groups): gain formula is astronomically large but lines have R_U = (r/gamma)^{..} > 1 always (R3 Cor 4.1);
for Delta_p the line condition is m_+ m_- >= (p-1)^2/(4p), satisfiable because t rd_Delta = 2 sqrt p. The first-moment identity
has ONLY Delta-subspace terms [P]: sublattices attached to subfields/suborders are not separate terms. So no hidden Delta = Q
disease in the first moment. But two things are hidden in "borderline-stable":
(a) NEW necessary conditions [P given: chi_W = 0 on split-torus elements; Lemma B => W_Z (x) Z_l projective over Z_l[H] for odd l
    when the Sylow-l subgroups of H lie in split tori]. For H < G with chi_W|_H supported on +-1 and W_Z locally free over
    Z[H]e_- the Q-isotypic component with simple algebra A_i (centre degree e_i, character degree chi_i(1), a_i = dim_Q A_i,
    Gamma_i = image of Z[H], m = 16r/|H|) has
        R_{W_i} = ( |disc_trd(Gamma_i)| / (e_i chi_i(1))^{a_i} )^m .
    Checks: cyclic C_h: (rd_d/phi(d))^{..} (= R5 2.7 remark); dicyclic N(T): same subspaces.
    At l = 2 W_Z is not projective; let 2^e be the actual 2-part of det(W_Z cap W_i):
    - H = SL_2(3) (all p), psi_0 (Hurwitz) component, f = 1/3: R = 2^{-e}, 0 <= e <= 4r/3: never unstable, borderline iff e = 0.
      (p = 13: these are 91 Delta-lines of Leech with det 3^8 2^e.)
    - H = 2I (p = 1 mod 5: p = 61, 421, 541, 661, ...), icosian component, f = 2/15, m = 2r/15: det = 5^{4m} 3^{8m} 2^e,
      R = (625/256)^m 2^{-e}: STABLE IFF e >= 1.288 m = 0.172 r. p = 61: need e >= 3 (range 0..16); if e = 0, R = 5.96 > 1 and
      W_61 is tensor-unstable. Other 2I components (degree 4, 6) are stable through their 5-part alone. [O: e not computed.]
(b) Borderline-stability does NOT imply E# = V(1+o(1)) at the relevant V without counting the borderline classes (A2).

## Overall estimate and corrected claim

Best case under the family's own heuristics: V <= (p(p+1)/(p-1)) n^2/(2e) with p <= sqrt(2n): 0.26 n^{2.5} (needs a trace-free
process at k = r+1); with the argument as it exists: V ~ n^{2.5}/sqrt(log n) (k ~ r log n); at k = c r: n^{2+(c-1)/(2(c+1))};
at k = r+1: n^{2.0}.

CORRECTED CONDITIONAL CLAIM. Suppose for infinitely many p = 1 mod 4 with (p+1)/2 prime:
 (H1) W_p unimodular with maximal endomorphism order and R_U <= 1 for all proper Delta-subspaces (including the explicit
      test e >= 0.172 r for the icosian component when p = 1 mod 5);
 (H2) the effectively-borderline subspaces (R_U^{k/2} >= n^{-A}) satisfy sum over G-classes T_U(wD) |H_U|/|N_G(U)| = o(D_k)
      and sum over G-classes T_U(V_0)|H_U|/|N_G(U)| <= 1/4 for the starting volume V_0;
 (H3) the equivariant Klartag argument holds for tensor families with r >= 2 (mean value formula with lower-rank terms and
      Steinitz classes, weighted orbit counting) in the regime k >= C r log n.
Then for n = 8rk with k = ceil(C r log p): V_n >= (1-o(1)) p n^2/(2e) >= c n^{5/2}/sqrt(log n) along this sparse sequence
(still contradicting V_n <= n^2 polylog n). The version with k = r+1, n ~ p^2/2 and V ~ 0.26 n^{5/2} requires in addition
 (H4) a trace-free (volume-preserving) version of the process valid when n/D = 2r/k is bounded below;
without (H4) the same heuristics give only V ~ V_0^{r/(k+r)} (wD)^{k/(k+r)} ~ n^2 at k = r+1.
