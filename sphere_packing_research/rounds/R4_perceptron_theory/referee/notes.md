# Referee report on R4_perceptron_theory  (25 Nov 2026, 07:38-07:55 UTC)
Tags: [P] proved here, [PM] proved as a statement about a model, [N] numerical (my own code, this directory), [H] heuristic, [R] refuted/corrected.
Code: ref_gauss_lp.py, ref_gauss_active.py (Gaussian model LP), ref_poisson_bruteforce.py (exact Poisson LP, brute-force pool + a-posteriori
exactness certificate, no row generation), ref_gamma_test.py (test of the exact chi-square law). Logs: *.log.

## 1. Constants re-derived independently (all agree with the reviewed notes)
- slack: x=R(1+eps)u, A=R^-2(I+B0): x^TAx>=1 <=> u^TB0u >= -(1-(1+eps)^-2) =: -s. Exact. #pairs with slack<=s = (v/2)(1-s)^(-n/2) -> (v/2)e^{sigma/2}; intensity (v/4)e^{sigma/2}. [P given Poisson]
- Var(u^TB0u)=2|B0|^2/(n(n+2)) (4th moments of the sphere); g=(n/sqrt2)(uu^T-I/n): <theta,g>=n u^TB0u, Var = rho^2 n/(n+2), |g|^2=N-n. [P]
- S(r)=int (v/4)e^{s/2}E(s+rZ)_-^2 ds = (v/4)*Gamma(3)/(1/2)^3*E e^{-rZ/2} = 4v e^{r^2/8}; Campbell variance (v/4)*4!/(1/2)^5*e^{r^2/8}=192V. [P]
- AO: phi = min_{r<=rho} ( |(sigma - r g)_-| - r|h| )_+ ; balance 4v e^{r^2/8} = r^2 N ; max_x x e^{-x/8}/4 at x=8: v*=(2/e)N ; V=e v*=2N ; active = (v/2)e^{r^2/8}=N. [PM]
- rho^2 -> 8 at the LP optimum needs BOTH directions: ball feasibility (lower bound on |theta*|) and the one-sided sphere bound (upper bound). Fine.
- pair formula exp(-V Phi(delta/4)), bulk term 2 sqrt(vV)/(rho sqrt(pi)), fluid c(Q)=v(e^{Q/8}Phi(sqrt Q/2)-1/2), Rezvyakova optimisation: algebra correct.
- MZZ (checked in the downloaded text, eq.(7), Thm 4.1/4.2 lines 609/632, lines 404-410): delta_s=(1+o(1))log|kappa|/Phi(kappa); gap "of order log|kappa|" to THEIR LP algorithm;
  they explicitly leave open whether other algorithms do better (abstract; lines 399-410, Section 6).

## 2. NEW exact theorem that replaces the numerics (and half of C1)  [P / PM]
**Theorem R (Gamma law).** (a) Exact Poisson model (pairs +-x, intensity 1/2 per unit volume, any n>=2): the minimiser A* of Tr A over
{A : x^TAx>=1 for all points} exists, is a.s. unique and positive definite, and V = det(A*)^{-1/2} satisfies **V ~ chi^2_{2N}** (V/2 ~ Gamma(N,1)),
N=n(n+1)/2. In particular E V = 2N = n(n+1) EXACTLY for every n, sd(V) = 2 sqrt N.
(b) Gaussian model (g_i iid of ANY law, offsets PPP (1/4)e^{s/2}ds): at the LP optimum (theta*,v*), the quantity
V := 2 Lambda{(s,g): s + <theta*,g> - 2 log v* < 0} (= v* e^{|theta*|^2/8} for Gaussian g) is exactly chi^2_{2N}.
*Proof of (a).* Existence: A=cI is feasible; feasible A are PSD (points in all directions at infinity); a singular minimiser would give a cylinder of
infinite volume, which contains points a.s.; near a PD form only finitely many constraints matter, so the minimiser is a vertex with N active pairs
(a.s. exactly N, a.s. unique). KKT (necessary and, by weak duality Tr A = sum lambda_x x^TAx >= sum lambda_x, sufficient): an N-subset S of the process
is the active set iff (i) the unique form A_S with x^TA_Sx=1 on S is PD, (ii) I = sum_{x in S} lambda_x xx^T with lambda>=0, (iii) the open ellipsoid E_S
contains no point. Exactly one S qualifies, so by the multivariate Mecke formula, for every Borel set W,
 P(V in W) = (1/N!) (1/2)^N int 1[(i),(ii)] 1[vol E_S in W] exp(-vol(E_S)/2) dx_1..dx_N .
(i),(ii) are invariant under x -> t x, and vol(E) scales by t^n. Polar coordinates x = t y on (R^n)^N with vol E(y)=1: dx = t^{nN-1}dt dmu(y),
w = t^n: P(V in dw) = C w^{N-1} e^{-w/2} dw. QED. (b) is the same with the translation s -> s+c in place of the dilation (Lambda scales by e^{c/2}).
This is the classical Gamma-type theorem for Poisson stopping sets (Miles; Moller-Zuyev; Zuyev 1999): the LP-optimal ellipsoid is a
scale-covariant stopping set determined by exactly N points.
*Checks [N].* n=4 (N=10): 400 samples, mean V=20.06+-0.30 (pred 20), var 35.9 (pred 40), KS p=0.65. n=6 (N=21): 149 samples, mean 40.9+-0.68 (42),
var 68.4 (84), KS p=0.31. n=8: 12 samples mean 69.1 sd 13.5 (72, 12). n=10: 6 samples mean 118.2 (110). n=12: 5 certified samples mean 162.4 (156).
Gaussian model: sd(V/N) observed 0.272 (N=55, 16 seeds), 0.176 (N=105, 8 seeds); predicted 2/sqrt N = 0.270, 0.195.
*Consequences.* (1) "V_LP = 2N" needs neither Gordon nor N -> infinity nor Gaussianity: it is exactly the Poisson void law on a stopping set with N
contact pairs. Gordon is needed only for the split V = v e^{rho^2/8} (v*=(2/e)N, rho^2=8). (2) The scatter in all LP tables (reviewed and mine) is the
chi^2 law, not solver noise; "V/N = 1.86..2.58" is exactly what chi^2_{2N}/N looks like. (3) The lattice transfer is NOT a universality question about
Veronese vectors but a Poisson question on a set of volume 2N >> n with N >> n marked points: Rogers' k-point formula exists only for k <= n-1.

## 3. The exact Poisson LP exceeds the true lattice optimum in low dimensions  [N + literature values of V_n]
V_n = 2^n * (max lattice density): V_4=9.870 (D4), V_6=23.87 (E6), V_8=64.94 (E8); 2N = 20, 42, 72.
My samples: n=4: 96.3% of 400 Poisson LP optima have V > V_4; n=6: 149/149 have V > V_6 (and 12/12 of the first batch, range 31.2-61.0);
n=8: 7/12 have V > V_8 (range 52.7-92.6). Also the reviewed table's n=8 values (V = 78-93) are all above V_8.
So for n<=8 NO lattice can have a trace-LP optimum of volume 2N; "V/n^2 = 1.0-1.45 for n=8..22" is a property of Poisson points, not evidence about lattices.
(For n>=12 the best known lattices exceed 2N, so the asymptotic conjecture is not contradicted.)

## 4. Corrections to the reviewed notes
- [R] Sect.2 "active offsets ... Gaussian with mean r^2/2=4 and sd r=2.83": wrong. Density of active offsets is prop. to e^{s/2}Phi(-s/r), i.e.
  sigma_act = r^2/2 + r Z - Exp(mean 2): mean 2, sd sqrt(12)=3.46. My runs: mean 1.86, sd 3.55 (N=55, cutoff 15, 12 seeds; 660 offsets); N=78: 1.50, 3.06 (6 seeds).
  So the active points lie in balls of volume v* e^{1 +- 1.7}, i.e. about 0.36N..11N (not "a few times N").
- [R] "truncation sigma<=11 ... zero violations beyond in every run": with cutoff 11 the expected number of constraints beyond the cutoff violated by the
  truncated optimum is (V/2)[Phi(-(S/r-r/2)) - e^{S/2-r^2/8}Phi(-S/r)] ~ 0.0023 N (0.24 at N=105, 0.7 at N=300). I checked up to cutoff+8 (not +3):
  violations in 2 of 8 runs at N=105 and 1 of 16 at N=55. So ~20-50% of the reported N>=105 optima are not optimal for the infinite system (bias on
  v* small, but "exact" is not justified). With cutoff 15: N=55: v/N=0.785+-0.078, rho^2=7.61+-0.45, V/N=1.971+-0.133; N=78: 0.917+-0.032, 6.89+-0.38, 2.17+-0.10.
- [R] "Gaussian model valid for rho << n^{1/4}": not sufficient. Needles B0 = a(I/n - uu^T) have rho^2 = 2a^2(1-1/n) < 2.2 and V/v -> infinity as
  a -> n/(n-1). Needles of aspect ratio Lambda have ~Lambda^{n-1} nearly disjoint directions, so empty ones exist up to V ~ 2(n-1) log Lambda:
  they beat the one-sided Gordon bound rho^2 N/4 (~ n^2/4 at rho^2 ~ 2) once Lambda > e^{n/8} and the LP value 2N once Lambda > e^{n/2}
  [N: n=100, lam_min(I+B0)=1e-12, aspect 1.0e6: 2(n-1)log(aspect)=2736 > rho^2 N/4=2550, rho^2=2.02]. So in the exact Poisson model the sphere bound is
  FALSE at rho^2 ~ 2 (exponentially thin needles only). The domain of validity of the Gaussian model must therefore include an operator-norm
  condition (e.g. lam_min(I+B0) >= e^{-o(n)}, or |B0|_op <= C|B0|_HS/sqrt n), not only rho << n^{1/4}. With such a condition I have no counterexample;
  without it the non-convex statements (sphere bound, 2N log rho) are statements about the Gaussian model only.
- [R] Sect.5(d) "algorithmic constructions sit at Theta(N) ... process-based proofs cannot reach n^2 log n": unsupported. MZZ prove a gap only for one
  LP algorithm and leave other algorithms open. In our problem iterating the LP from the new centre (= Voronoi-type ascent to an extreme form) is a
  polynomial-step procedure not covered by the model (the model is not invariant under change of base point).
- C3: the formula V = v e^{Q/8} needs rho^2 = Q_f = 8 log(N/v) << sqrt n, i.e. n >~ 10^4; for n<=10 (R2 data) only the robust relation
  "contacts ~ 2 x (expected points inside) = V" is being tested. By Theorem R's mechanism the SWEPT region of the sticky process (not its final
  ellipsoid) is the stopping set, consistent with "sweeps twice what it keeps".
- Theorem A: I found no misuse of Gordon (ball: convex-concave, two-sided; sphere: one-sided; random number of constraints is independent of G;
  truncation argument for feasibility is valid because the far constraints are independent of the truncated min-norm point; bound (V/2)Phi(-(S/rho-rho/2))
  re-derived). It is still a proof sketch (uniform-in-r concentration, Poisson number of rows) -> "conditional" is the right tag.

## 5. What I would state in a paper
(a) Theorems about models: Theorem R (chi^2_{2N}, exact, both models); the AO identities and, after writing out the concentration, Model Theorem A
    (v*=(2/e)N, rho^2=8, Gaussian model only); pair formula; needle statement.
(b) Numerical observations: Gaussian LP v*/N ~ 0.73-0.80, rho^2 ~ 7.4-8.1 at N = 55-171 (slow convergence, rho^2 biased low at small N);
    R2: sticky process on Z^n-type lattices V_final/N ~ 1 for n<=10.
(c) Conjectures for lattices: (LP) only as an asymptotic statement and with the warning of Sect.3; it must be tested on Haar-like random lattices
    (n = 10..16) before being stated. (EX) n^2 log n: I would NOT state it as a conjecture - it is the output of a cut-off-dependent model whose
    un-cut version is false (sup V = infinity); state it as a heuristic at most.

## 6. Addendum (numbers produced after the first draft)
- Theorem R(b) tested in the Gaussian model (cutoff 16): N=12, 250 samples: mean V=24.24+-0.47 (pred 24), var 55.0 (48), KS p=0.36, mean v*/N=0.892, mean rho^2=7.26;
  N=30, 80 samples: mean V=61.1+-1.2 (60), var 122.8 (120), KS p=0.62, mean v*/N=0.785, mean rho^2=8.08.
  => E V = 2N holds exactly at every N, while v*/N converges slowly from above to 2/e=0.7358 (0.89, 0.79, 0.76-0.80 at N=12, 30, 55-171).
- Expected number of constraints beyond the cutoff violated by the LP optimum (r^2=8), per unit N: cutoff 11: 0.00214; 14: 0.00005; 15: 0.00001.
- Gaussian LP, cutoff 11 (as in the reviewed runs), my seeds: N=55 (16): v/N=0.764+-0.041, rho^2=7.50+-0.51, V/N=1.929+-0.068; N=105 (8): 0.800+-0.050, 7.40+-0.45,
  1.994+-0.062; N=171 (3): 0.800+-0.006, 7.49+-0.74, 2.06+-0.19.
- Exact Poisson LP (brute force, certified): V/N = 2.06+-0.12 (n=6, 12), 1.92+-0.11 (n=8, 12), 2.15+-0.09 (n=10, 6), 2.08+-0.14 (n=12, 5 of 6 certified).
- Lattice maxima: V_n/(2N) = 0.605, 0.494, 0.493, 0.496, 0.568, 0.675, 0.902 for n=2..8 (A2, D3, D4, D5, E6, E7, E8).
- Not done: trace-LP on Haar-like random lattices (the decisive test of the lattice conjecture); R3 Lemma N and Rezvyakova's proof were not re-read.
