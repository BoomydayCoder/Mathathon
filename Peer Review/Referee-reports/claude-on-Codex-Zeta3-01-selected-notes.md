# Referee report: three notes from Codex-Zeta3-01-Zudilin-and-recurrence

Reviewer: Claude (AI agent), not a human referee. Files read in full, and only these:
`denominator-design.md` ("DD"), `recurrence-design.md` ("RD"), `apery-mixed-numerator-balance.md`
("MB"). For comparison I read the certificate lines of `zudilin-explained.md` in the same package.
Zudilin's paper was **not retrieved**; statements about it are about the package's account of it.

Notation: P_n(t) = ∏_{j=1}^n(t−j)/∏_{j=0}^n(t+j), N = n(n+1), M = m(m+1), d_n = lcm(1,…,n).

## 1. Statements reviewed
DD §1. log d_n ≤ cn + O(log²(n+1)), c = ½log 2 + ¾log 3 = 1.170533; this suffices for ζ(3)
(3c + 4log(√2−1) = −0.013896 < 0) and for ζ(2) (108((√5−1)/2)^{10} < 1).
DD §2–3. If P = Σ_{k=0}^n c_k/(t+k) with c_k ∈ Z and P^m = ΣΣ a_{j,k}/(t+k)^j, then
d_n^{m−j}a_{j,k} ∈ Z, and (−1)^r/r! Σ_{t≥1}(P^m)^{(r)} is a linear form in ζ(r+1),…,ζ(r+m) whose
coefficients are cleared by d_n^{m+r}; for m = 2 it is (r+1)(Σc_k²)ζ(r+2) − V.
DD §5. For L_n = −(1/24)Σ_t (P_n²)'''(t) (a form A_nζ(5) − V_n/4), the Apéry operator T gives
(TL)_n = −(2n+1)[3 + 8n(n+1)(H_{n−1}+H_{n+1})]/(n²(n+1)²), so L_n does not decay exponentially.

RD §2. If (n+1)^d u_{n+1} = b(n)u_n − n^d u_{n−1} with b of degree d and leading coefficient β > 2,
A_n positive integers, w = A_0B_1 − A_1B_0 ≠ 0, then B_n/A_n → L,
L − B_n/A_n = w Σ_{j≥n} 1/((j+1)^d A_jA_{j+1}), and L is irrational if d_nB_n ∈ Z with
limsup d_n^{1/n} < λ = (β+√(β²−4))/2. RD §3: log((a+1)/a) is irrational for integers a ≥ 1.
RD §5: a finite search.

MB. With S_{n,m} = (−1)^{n+m}P_nP_m and F_{n,m} = −½Σ_{t≥1}S'_{n,m}(t):
(8)  (n+1)³S_{n+1,m} − (2n+1)(N+2M+1)S_{n,m} + n³S_{n−1,m} = −2(2n+1)[Q(t+1) − Q(t)],
     Q = (t−m−1)(t+m)S_{n,m};
(11) P_n² + P_{n−1}² + 6P_nP_{n−1} = T_n(t+1) − T_n(t),  T_n = (4/n²)(t+n)²P_n²;
(12) F_{n,n} + F_{n−1,n−1} = 6F_{n,n−1};  (13) A_n + A_{n−1} = 6Π_n(n−1);
(14) Apéry's recurrence for f_n = F_{n,n}, for A_n and for B_n, with 17 = 6·3 − 1, 5 = 6·1 − 1;
(1)  C_n(n−1) = (B_n + B_{n−1})/6.

## 2. Verdict
**All three notes: correct.** No error found. Two dependencies on files outside the package: MB (6)
(repaired in §4 below) and the identity in DD §5 (verified numerically only).

## 3. What I re-derived or recomputed
By hand:
- DD §1: the table of h on [0,6) is 0,1,1,0,0,1; a = (log 2)/3 + (log 3)/2; the factor 3/2; the two
  integer certificates 2·470832² − 593855² = 90701783423 and 5·2970² − 6641² = 1619. The two
  rejected combinations give 2log 2 = 1.386 and 1.214 as stated.
- DD §2: the weight count Σ(a_i+1) = m−j. DD §3: the coefficient C(j+r−1,r). DD §4: the residue.
  DD §5: P_n'(1), P_n''(1), −R_n'''(1)/6.
- RD §2: Casoratian W_n = w/(n+1)^d; the ratio argument (r_n ≥ 1 eventually, then r_n ≥ r); bound (2).
- RD §3: the n integrations by parts; B_1 = 2, w = 2; r_n ≥ 2a+1 and then r_n ≥ 4; 108 < 256.
- MB (8): all three ratios, and both sides after multiplying by (t−n)(t+n+1)/S; boundary values
  Q'(1) = 0 for n ≥ 1 and = 1 for n = 0; F_{0,m} = ζ(3) − H_m^{(3)}.
- MB (11): left side = 4(2t²−n²)r(t) with r = P_n²/(t−n)²; r(t+1)/r(t) = t⁴/((t−n+1)²(t+n+1)²).
- MB (12), including the sign of S_{n,n−1}; MB (13) from residues (−1)^{n+k}C(n,k)C(n+k,k).
- MB (14): 6·[mixed recurrence at m = n] with (n+1)³ + n³ = (2n+1)(N+1) gives the coefficient
  (2n+1)[(N+1) − 6(3N+1)] = −(2n+1)(17N+5).
- MB (4): it is identity (8) of `Codex-Zeta3-04/resolvent-legendre-pade-proof.md` with one index 0.

By machine:
- MB (8): the polynomial identity with n, m, t symbolic; as rational functions for 0 ≤ n,m ≤ 3.
  MB (11) for n ≤ 5. MB (6) for n ≤ 5 with y symbolic.
- F_{n,n} = A_nζ(3) − B_n, (12), (13) and F_{n,n−1} = Π_n(n−1)ζ(3) − (B_n+B_{n−1})/6 for n ≤ 4, to
  36 digits or exactly.
- DD §5: the formula for (TL)_n at n = 1,…,4, to 35 digits. L_n = 0.18, 0.039, 0.012, 0.0052.
- RD §5: the scan reproduced exactly. d = 2: no pair in 3 ≤ a ≤ 40, 1 ≤ b ≤ 20. d = 3: only (17,5).

Not checked: MB §2 beyond the algebra (the link to Prévost's associated numerator; not needed for
the recurrence); the proof of the DD §5 identity.

## 4. Gaps and repairs
1. MB §2, equation (6), is quoted from `newton-bispectral-motivation.md`, which is not in the package.
   MB §5 uses (6) for the recurrence of A_n. **Repair, two lines:** take the sum of the double-pole
   coefficients in (8). For S_{n,m} this sum is Σ_k C_{n,k}C_{m,k} = Π_n(m). The right side is
   Q(t+1) − Q(t) with Q bounded at infinity and with at most double poles, so its sum is 0. Hence
   (6) holds at every integer y = m ≥ 0, and both sides are polynomials in y.
2. DD §5: the identity for (TL)_n is quoted from `forced-recurrence-obstruction.md`, not in the
   package. Status here: verified numerically only (n ≤ 4, 35 digits). The deduction from it is
   correct: |TL_n| ≥ c(log n)/n, while L_n = O(ρ^n) would give TL_n = O(n³ρ^n).
3. RD §2: state "A_n > 0 for all n ≥ 0" as a hypothesis in the theorem line (it is used to divide).
4. MB "Theorem" uses C_n(n) = B_n "by Prévost's formula" before §5 proves it independently. Reorder:
   define B_n := C_n(n); §5 then shows it is the Apéry solution with B_0 = 0, B_1 = 6.

## 5. Different from, or simpler than, the standard proofs?
- DD §1 is not a proof of irrationality. It is a self-contained replacement for d_n < 3^n. It is
  weaker (e^{1.1705} = 3.2237 per step) but its proof is ten lines. The margin for ζ(3) is thin:
  e^{3c} = 33.50 against λ = 33.97. It is the simplest lcm estimate I know that is sufficient.
- RD §2 is the recurrence route with the analytic part reduced to "A_nξ − B_n is bounded". No
  asymptotic of the linear form is needed: the Casoratian gives |A_nL − B_n| ≤ C/(n^dA_n), and growth
  of A_n follows from integrality. This is simpler than the saddle-point estimate described in
  `zudilin-explained.md` §§5–7. The price is that the recurrence must be proved for both sequences.
- MB §§3–5 is a different *derivation of Apéry's recurrence*, not a different irrationality proof:
  the linear form f_n = −½Σ(P_n²)' is the rational-function one. According to
  `zudilin-explained.md`, Zudilin's certificate for P_n² is s_n(t)P_n(t)² with
  s_n = 4(2n+1)(−2t²+t+(2n+1)²). MB replaces it by two certificates, each checked by expanding one
  quadratic or quartic: Q = (t−m−1)(t+m)S_{n,m} and T_n = (4/n²)(t+n)²P_n². The length is about the
  same. The gain is explanatory: the off-diagonal family F_{n,m} has an eigenvalue-type recurrence
  in n with spectral parameter M, and one adjacent relation closes it on the diagonal.

## 6. Attribution
Classical: Chebyshev-type factorial-ratio bounds for d_n; the Casoratian argument; Legendre
approximations to logarithms; partial-fraction arithmetic (Nesterenko, Zudilin; not retrieved).
Rearrangement: the period-six choice; the three-requirements theorem RD §2. We did not find the pair
(8), (11) in the sources we read. The polynomials Π_n(y) = Σ_k C(n,k)C(n+k,k)C(y,k)C(y+k,k) are
classical (the package attributes them to Prévost and Wilson), so (6) is likely known.

## 7. Motivation
Stated in all three notes. DD explains how the period-six combination was found. RD separates
identification, arithmetic and recurrence. MB explains 17 and 5. One addition for MB: the reason
to introduce a second index m is that P_nP_m has the simple shift ratio
S(t+1)/S(t) = t⁴/((t²+t−N)(t²+t−M)), symmetric in the two quadratic variables N and M, so the
recurrence in n at fixed m has a degree-two certificate. On the diagonal m = n this structure is
hidden.

## 8. The idea most worth a mathematician's time
Apéry's recurrence from two one-line telescopes. Put S_{n,m} = (−1)^{n+m}P_nP_m. Then

    (n+1)³S_{n+1,m} − (2n+1)(N + 2M + 1)S_{n,m} + n³S_{n−1,m} = −2(2n+1)Δ_t[(t−m−1)(t+m)S_{n,m}],
    P_n² + P_{n−1}² + 6P_nP_{n−1} = Δ_t[(4/n²)(t+n)²P_n²],          Δ_t g = g(t+1) − g(t).

Apply −½Σ_{t≥1} d/dt (linear forms F_{n,m}) or "sum of double-pole coefficients" (integers Π_n(m)).
Both obey u_{n+1,m}(n+1)³ − (2n+1)(N+2M+1)u_{n,m} + n³u_{n−1,m} = 0 for n ≥ 1 and
u_{n,n} + u_{n−1,n−1} = 6u_{n,n−1}. Setting m = n: (2n+1)[(N+1) − 6(3N+1)] = −(2n+1)(17N+5).
Questions: is there a two-index family with such a pair of telescopes for the ζ(2) recurrence
(11n²+11n+3), and for the other integral recurrences of this shape? Does a three-index family
exist for ζ(5)-type forms, where DD §5 shows the one-index construction fails?
