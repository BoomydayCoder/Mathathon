# Referee report: Codex-Zeta3-05-Schur-Apery
Reviewer: Claude (AI agent), not a human referee. Read in full: the final-status note, the design
guide, `rank-one-real-boundary-proof.md` ("RB"), `rank-one-original-apery-heun-rii-proof.md` ("HR"),
`rank-one-zeta-subspace-approximants.md` ("SUB"), `rank-one-loewner-quadratic-lattice.md` ("LOE"),
`rank-one-original-error-tail-and-rate.md` ("TAIL"). Skimmed: `zeta3-rank-q-heun-reduction.md`.
The package's own review notes were treated as claims.

## 1. Statements reviewed
Setting (the definition is **not in the package**, see §6): f(y) = 1/(e^{2πy}−1), w_3(y) = y³f''(y),
dμ_3(y²) = w_3(y)dy, D_M(t) = ∏_{j=1}^M (t+j²). The pole values are
∫dμ_3/(t+j²) = j²ζ(3) − j²H_j^{(3)} − ½ + 1/(2j). Replacing ζ(3) by X gives an affine functional
μ_{3,X} on rational functions with simple poles at −j².
- (SUB §2) The h×h Gram matrix G_h(X) = (μ_{3,X}(t^{i+j}/D_{2h})) depends on X only in its last
  entry, with coefficient −1. Its unique root q_h is rational, and
  q_h − ζ(3) = min over monic p of degree h−1 of ∫ p(y²)² w_3(y)/D_{2h}(y²) dy > 0.
- (RB 1) For even f, g with deg_{x²}f + deg_{x²}g ≤ M−1:
  ℓ_{M,X}(f𝒟_M g − g𝒟_M f) = −(M/(2(M!)²))[f(0)g(1) − f(1)g(0)].
- (HR T) q_h = (5B_h − B_{h−1})/(5A_h − A_{h−1}) for h ≥ 1, with A, B the two solutions of Apéry's
  recurrence, (A_0,A_1) = (1,5), (B_0,B_1) = (0,6).
- (TAIL) q_h − B_h/A_h = 6/(h³A_hC_h), C_h = 5A_h − A_{h−1}; B_h/A_h < ζ(3) < q_h;
  (q_h − ζ(3))/(ζ(3) − B_h/A_h) → 3+2√2; q_h − ζ(3) ~ 4π³(3−2√2)(17+12√2)^{−2h}.
  Rank two: the two roots of the rank-two pencil lie in Q(√(8h²−12h+5)).

## 2. Verdict
**(T), RB (1), TAIL (2)–(6): correct.** No error found. (T) is complete once two facts cited from
outside the package are supplied (§6 a, b); I supply (b) and checked (a) numerically. TAIL (7) and
rank q depend on files I did not referee. **Not a new irrationality proof**; the package agrees (§5).

## 3. What I re-derived or recomputed
By hand:
- RB: w_j = 2(−1)^{j−1}j²/((M−j)!(M+j)!); cancellation of the poles of 𝒟_M at 0 and ±½ on even
  polynomials; detailed balance α_jA_M(j) = α_{j+1}A_M(−j−1); ρ_j − ρ_{j+1} = (2j+1)/(2j³(j+1)³);
  the summation by parts with A_M(M) = A_M(−1) = 0; δ_j; K(−x−1) = −K(x); deg K ≤ 2M−3; the
  vanishing (2M−1)-th difference. All correct.
- HR §2: symmetry of H = X𝒟 + 𝒟X − 𝒟 − κX (the defects of (Xf,g) and (f,Xg) add to the defect of
  (f,g)); r_k = 2(k−n)(k−n+1); the eigen-argument; R_h(1) and R_h(4) from the values at x = 0, 1.
- HR §3: the old/new node bookkeeping (F_h(−t) = (t+(2h−1)²)(t+(2h)²); test degrees k+1 ≤ h−2 and
  k+2 ≤ h−1); W ∝ R_{h−1}.
- HR §4: T has 2h zeros and degree ≤ h+1; the leading coefficients give (12).
- HR §5: α_h = a_h s_h/s_{h+1}, γ_h = c_h s_{h−1}/s_{h+1}.
- TAIL (4) from A_{h−1}B_h − A_hB_{h−1} = 6/h³; TAIL (6): (λ−5)/(5−λ^{−1}) = 3+2√2.
- The Newton action 𝒟_Mφ_k = k(k+2−M)φ_k + k(k+1)²(k−M)φ_{k−1}, which HR §2 calls "already proved"
  without a proof in the package. Proof: φ_k(x+1)/φ_k(x) = x(x+k+1)/((x−k)(x+1)), so
  𝒟_Mφ_k = (kφ_{k−1}/(2x))[g(x) − g(−x)], g(x) = (x+1)²(x−M)(x+k), and
  g(x) − g(−x) = 2x[(k+2−M)x² + (k − M − 2kM)].
- The interpolation input of HR §1 (from LOE §2), by a shorter argument: Σ_j w_j g(−j²) = [t^{2h−1}]g
  for deg g ≤ 2h−1. At X = q_h, v_j = f_j(q_h)p_h(−j²) satisfies Σ_j w_j v_j π(−j²) = 0 for all
  deg π ≤ h−1. The annihilator of these π has dimension h and contains the value vectors of all
  polynomials of degree ≤ h−1, so v_j = N_h(−j²) with deg N_h ≤ h−1.
By machine:
- The pole values by quadrature at j = 1, 2, 5 (30 digits); μ_3 total mass ¼; J_3 = ½.
- q_h from the exact rational Gram matrix, h = 1,…,10: equal to (5B_h−B_{h−1})/(5A_h−A_{h−1}) each
  time (5/4, 577/480, 38687/32184, …). For the same h: R_h(1) = h²(h+1)/2,
  R_h(4) = h²(h+1)²(3h³+3h−2)/(16(2h−1)), 1/Q_h(0) = h²C_h/(12(2h)!).
- H_hQ_h = η_hQ_h, η_h = 1−2h−2h³(h+1), exactly for h ≤ 6; (H f)(0) = −2Mf(1) + (1−M)f(0).
- HR (10) by solving the three endpoint equations with h symbolic; HR (13) for 5u_h − u_{h−1}.
- RB (1) for M = 5, all degree pairs with r+s ≤ 4, X symbolic, random coefficients.
- Rank two: the discriminant has square class 13 at h = 2 and 41 at h = 3 (= 8h²−12h+5); at h = 11
  the two roots are exactly the rational numbers displayed in `zeta3-rank-q-heun-reduction.md`,
  with errors −7.1·10^{−34} and +4.3·10^{−29}.
- TAIL (7): at h = 10, (q_h − ζ(3))/(4π³(3−2√2)λ^{−2h}) = 1.114. Consistent, not a proof.
Not checked: rank-q proofs; "all-rank endpoint basis"; "remainder isometry"; the rate for other odd
s; all statements about cited papers (not retrieved by me).

## 4. Shortest correct proof of (T)
1. Pole values; G_h(X) = G_h(0) − Xee^T; variational formula (SUB §2).
2. Interpolation: N_h/Q_h takes the values j²q_h + β_j at j² (argument in §3 above).
3. Boundary identity RB (1): partial fractions, detailed balance, one vanishing finite difference.
4. H is symmetric with two consecutive degree cut-offs, so Q_h is an eigenfunction; evaluation at
   x = 0 and x = 1 gives R_h(1), R_h(4).
5. Three-term relation (HR 9) between consecutive measures; coefficients from z = 0, 1, 4.
6. The numerator obeys the same relation (root count); compare leading coefficients; rescale by
   h²/(12(2h)!); match two initial values. About four pages; §7 suggests a shorter route to try.

## 5. Is it a different proof of irrationality?
No. The linear form is D_h − C_hζ(3) = f_{h−1} − 5f_h with f_n = A_nζ(3) − B_n Apéry's form. Its
integrality (C_h ∈ Z, d_h³D_h ∈ Z), its size and Apéry's recurrence are imported; TAIL §4 says so.
The error is 3+2√2 times Apéry's. What differs is the *characterization*: an Apéry-type approximant
is the value of a positive least-squares problem, so rationality, sign and strict monotonicity are
immediate. Arithmetic is exactly what this description does not give. No step is simpler than in
the standard proofs.

## 6. Dependencies outside the package, and numerics-only claims
a. Definition of μ_3 and the pole values F_j(X) (used in RB, HR, LOE, SUB): not in the package (the
   density is in `Codex-Zeta3-01-…/independent-review.md`). Fix: add the derivation to SUB §1.
b. HR §2 "already proved Newton action": no proof in the package. Fix: the three lines in §3 above.
c. TAIL §3: constants C_A, C_f come from `apery-asymptotics.md` (package 01). Not refereed here.
d. TAIL §3, end (same rate for every odd s): needs `rank-one-universal-error-rate.md`, absent. Unverified.
e. Rank q ≥ 2: symmetry of H uses the all-degree identity of `rank-one-wilson-boundary-pearson.md`
   (a contour proof), not in the package. Unverified; my h = 2, 3, 11 computations agree with it.
f. SUB §6 table (q_{12,5}: 217-digit denominator, error 4·10^{−37}): numerics only, and so labelled.

## 7. Attribution and a remark the package does not make
Classical: rational Gauss quadrature and multipoint Padé (sources cited by the package, not
retrieved by me); Apéry's recurrence and arithmetic. Rearrangement: the boundary-defect repair
{r(X),𝒟} − (r(a)+r(b))𝒟 and the two-cutoff argument. We did not find (T) in the sources we read.

Remark (mine; proved from `apery-mixed-numerator-balance.md` (12), and checked to 30 digits for
h = 1, 2, 3). With P_n(t) = ∏_{j=1}^n(t−j)/∏_{j=0}^n(t+j) and F_{n,m} = −½Σ_{t≥1}((−1)^{n+m}P_nP_m)',
the balance f_h + f_{h−1} = 6F_{h,h−1} gives f_{h−1} − 5f_h = 6(F_{h,h−1} − F_{h,h}), that is

    (5B_h − B_{h−1}) − (5A_h − A_{h−1})ζ(3) = 6 Σ_{t≥1} d/dt [ t P_h(t)P_{h−1}(t)/(t+h) ].

So the Schur form is one rational-function series, and 5 = 6 − 1. A direct proof that the
least-squares problem of §1 produces this rational function would replace steps 3–6.

## 8. Motivation
The design guide states why the defect forces H, why M = 2h, why three endpoint values suffice. To
add: why 1/D_{2h}. With 2h poles the X-coefficient of the moments vanishes below degree 2h−2.

## 9. The idea most worth a mathematician's time
A positive extremal problem whose value is an Apéry-type error:

    min_{p monic, deg p = h−1} ∫_0^∞ p(y²)² w_3(y)/∏_{j=1}^{2h}(y²+j²) dy
        = (5B_h − B_{h−1})/(5A_h − A_{h−1}) − ζ(3) > 0,   w_3(y) = y³ (1/(e^{2πy}−1))''.

The same construction with w_s(y) = 2y^s f^{(s−1)}(y)/(s−1)! gives rational, strictly monotone,
one-sided approximants to every ζ(2m+1) with error ≤ J_s((h−1)!)⁴/((2h+m−1)!)². For s ≥ 5 the
denominators are observed to grow far too fast. Problems: (i) prove (T) directly from the series in
§7; (ii) explain the denominator growth for s = 5, where only arithmetic fails.
