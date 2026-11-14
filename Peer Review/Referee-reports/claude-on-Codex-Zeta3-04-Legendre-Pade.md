# Referee report: Codex-Zeta3-04-Legendre-Pade
Reviewer: Claude (AI agent), not a human referee. Read in full: `resolvent-legendre-pade-proof.md`
("RLP") and `half-shift-elementary-resolvent-proof.md` ("HS"). Not read: the two cubic notes,
`designing-the-apery-measure.md`. The package's own review notes were treated as claims.

## 1. Statement reviewed
With dν(t) = π² sech²(π√t) tanh(π√t) dt on t > 0, C_{n,k} = C(n,k)C(n+k,k),
Q_k(t) = (k!)^{-2} ∏_{r<k}(t+(r+½)²), Φ_n = Σ_k (−1)^k C_{n,k} Q_k, a = n+½:
- (RLP 2) ∫ dν/(t+a²) = 2ζ(3, a+½) for a > 0; (RLP 4) ν(t^j) = (−1)^j (2j+1) B_{2j}(½).
- (RLP 7) ν(Q_j Q_k) = 1/(j+k+1); hence ν(Φ_n p) = 0 for deg p < n and ν(Φ_n²) = 1/(2n+1).
- (RLP 11) Φ_n(−a²) = A_n = Σ_k C_{n,k}².
- (RLP 13–14) with N_n = ν((A_n − Φ_n)/(t+a²)), b_n = A_n H_n^{(3)} + N_n/2:
  0 < A_n ζ(3) − b_n = (1/(2A_n)) ∫ Φ_n²/(t+a²) dν ≤ 2/((2n+1)³ A_n).
- (RLP 15–18) N_n = Σ_{k=1}^n C_{n,k}² Σ_{m=1}^k (−1)^{m−1}/(m³ C_{n,m}); d_n³ b_n ∈ Z.
- (RLP 19) A_n ≥ λ^n/((n+1)(2n+1)²), λ = 17+12√2.
- (RLP §6) log d_n ≤ cn + O(log² n), c = ½log 2 + ¾log 3, 3c = 3.51160 < log λ = 3.52549.
- Conclusion: ζ(3) is irrational.

## 2. Verdict
**Correct.** No error and no gap that affects the conclusion; three small repairs (§7). The proof
is complete inside the two files, given Legendre's formula for v_p(n!) and Rodrigues' formula.

## 3. What I re-derived or recomputed
By hand (every line):
- HS §1: T_N(z)+1 = 2^{N−1}∏_{k=1}^N (z − cos((2k−1)π/N)); the substitution z = cosh 2u; the
  logarithmic derivative; the majorant sinh(2x)/(4(k−½)²); hence tanh(πy) = (2/π)Σ y/(y²+(j+½)²).
- HS §2: density = −y f''(y) with f = tanh(πy); ∂_a k_a = −h_a'; ∫k_a k_b = π/(2(a+b));
  ∫h_a'h_b' = π/(a+b)³; the bound ∫|h_a'h_b'| ≤ 1/(ab²).
- HS §3: M_e = (2e+1)π^{−2e}∫u^{2e} sech²u du, the exponential series, and the reduction to
  (2e+1)(−1)^e B_{2e}(½).
- RLP (5), (6); (8) by the root count (2j+2k+1 roots counted with multiplicity, leading coefficient
  2(j+k+1)); E_{jk}''(0) = 2(−1)^{j+k}(j!k!)²; hence (7).
- RLP (13), (14): the algebra from orthogonality, and ∫Φ²/(t+a²) ≤ a^{−2}·1/(2n+1).
- RLP (15): the divided-difference formula, the tail product, and the final simplification, which the
  file does not display: the j-th term is (j!)²(n−j−1)!/((j+1)(n+j+1)!) = 1/(m³C_{n,m}), m = j+1.
- RLP (16)–(17): the valuation count, including the p = 2 case m = 2^{ν_2}.
- RLP §5: Σ S_n z^n = (1−6z+z²)^{−1/2} from C_{n,k} = C(2k,k)C(n+k,2k); the coefficient bound;
  Cauchy–Schwarz.
- RLP §6: the period-six table 0,1,1,0,0,1; the linear coefficient (log 2)/3 + (log 3)/2; the factor
  3/2; 2⁶·3⁹ = 1259712 < (169/5)⁴ = 1305169.15….
By machine:
- (2) by quadrature at a = 0.3, 1, 2.5, 7: agreement to 30 digits. (4) for j ≤ 4: 30 digits.
- (7) exactly for 0 ≤ j,k ≤ 6 from the moments (4); (8) symbolically for j,k ≤ 4.
- Φ_n(−a²) = A_n, ν(Φ_n²) = 1/(2n+1), (15) against the direct definition of N_n, and
  b_n = Apéry's numerator (0, 6, 351/4, 62531/36, 11424695/288, …), all exactly for n ≤ 6.
- (14) numerically for n ≤ 6: the ratio of the true error to the bound is about 0.70.
- (16) for all 37 820 triples 1 ≤ m ≤ k ≤ n ≤ 60: integer and even in every case.
Not checked: RLP "optional Fourier verification" (not used); RLP §7 and (20); every statement about
the content of Prévost's paper (**not retrieved**; unverified).

## 4. Shortest correct proof extractable
1. tanh(πy) = (2/π)Σ_j y/(y²+(j+½)²)  (HS §1).
2. One integration by parts and ∫h_a'h_b' = π/(a+b)³ give ∫dν/(t+a²) = 2ζ(3,a+½)  (HS §2).
3. Moments (4), then ν(g) = L_x[(x+½)g(−(x+½)²)] with L(f(x+1)−f(x)) = f''(0)  (HS §3, RLP 5–6).
4. E_{jk}(x) = x²∏_{r≤j}(x²−r²)∏_{r≤k}(x²−r²) satisfies
   E(x+1) − E(x) = (j+k+1)(j!k!)²(2x+1)φ_jφ_k, so ν(Q_jQ_k) = 1/(j+k+1)  (RLP 7–8).
5. Transport Legendre: orthogonality, norm 1/(2n+1), A_n = Φ_n(−(n+½)²)  (RLP 9–11).
6. Error identity and bound (13)–(14).
7. Formula (15) and the valuation lemma (16).
8. A_n ≥ λ^n/((n+1)(2n+1)²); the lcm bound; contradiction.

## 5. Is it a different proof?
The approximants are the same: A_n and b_n are Apéry's numbers (checked exactly, n ≤ 6), so the
linear forms A_nζ(3) − b_n equal those of Beukers and of the rational-function proof. The package
itself attributes the route to Prévost (1996). The differences are in how the three properties are
obtained:

| Property | Beukers (integral) | Rational-function proof | This package |
|---|---|---|---|
| Linear form | triple integral with P_n(x)P_n(y) | −½Σ_t (P_n(t)²)' | (1/(2A_n))∫Φ_n²/(t+a²)dν |
| Non-vanishing | positive integrand | recurrence/positive series | positive integrand |
| Smallness | max of integrand = (√2−1)^{4n} | recurrence or asymptotics | **norm 1/(2n+1) and A_n large** |
| d_n³ | integrate monomials | partial fractions | formula (15) + valuation lemma (16) |

The step that differs is smallness: no bound of an integrand and no asymptotic of the linear form.
The error is at most 2/((2n+1)³A_n), so one only needs a *lower* bound for the integer A_n. This is
the Casoratian relation A_n(A_nζ(3)−b_n) ≈ const/n³, obtained without a recurrence. Legendre
polynomials enter through the isometry u^k ↦ Q_k, not through integration by parts.

## 6. Is any step simpler than in the standard proofs?
- Simpler: the explanation of A_n = Σ C(n,k)²C(n+k,k)² (the Legendre coefficient (−1)^kC_{n,k} times
  Q_k(−(n+½)²) = (−1)^kC_{n,k}); the smallness step; the growth bound from (1−6z+z²)^{−1/2}.
- Not simpler: the analytic input (steps 1–3 take about two pages); the denominator step, which is
  Apéry's formula with a valuation argument, harder than Beukers' or the partial-fraction argument.
  The whole proof is longer than Beukers'. It is more transparent on "why these numbers".

## 7. Repairs and dependencies
- RLP §2, (6): no proof is written ("The Bernoulli generating function also proves"). Add: from
  B(e^{zx}) = z/(e^z−1), B(f(x+1)−f(x)) = f'(0) and B(f(−x−1)) = B(f(x)); apply to f'.
- RLP §4, (15): "simplifies precisely to (15)" should display the one-line simplification in §3 above.
- RLP §4: ν_p denotes ⌊log_p n⌋ while ν is the measure. Rename.
- RLP §7 depends on scanned pages of Prévost's paper that are not in the package. Unverified here.
- Nothing in the irrationality proof rests on numerics only.

## 8. Attribution
Classical: Legendre polynomials; Apéry's formula (15) for b_n; the Padé/Stieltjes error identity; the
central Delannoy generating function; Chebyshev-type lcm bounds. Attributed by the package to
Prévost: the whole Wilson–Padé route (not verified by me). Rearrangement: the real proof of the
resolvent identity; the elementary proof of (7) by the difference identity (8); replacing the growth
asymptotic and the prime number theorem by (19) and the period-six bound. We did not find (7)–(8) in
this form in the sources we read, which is a statement about a small set of sources only.

## 9. Motivation
Stated well in RLP "Design route". To add: the measure is forced. One wants a positive ν whose
Stieltjes transform at a = n+½ is 2(ζ(3) − H_n^{(3)}); 2ζ(3,a+½) = ∫u²e^{−au}/(2 sinh(u/2)) du
determines ν, and ν is positive because −y(tanh πy)'' > 0.

## 10. The idea most worth a mathematician's time
The Hilbert-matrix Gram identity. For dν = π² sech²(π√t) tanh(π√t) dt,

    ν(Q_j Q_k) = 1/(j+k+1) = ∫_0^1 u^{j+k} du,   Q_k(t) = ∏_{r<k}(t+(r+½)²)/(k!)².

So u^k ↦ Q_k is an isometry from polynomials on [0,1] to polynomials in L²(ν). The image of the
shifted Legendre polynomial is the orthogonal polynomial Φ_n, with norm² 1/(2n+1), and
Φ_n(−(n+½)²) = Σ C(n,k)²C(n+k,k)². The proof is the single difference identity (8). Consequence:

    0 < A_nζ(3) − b_n = (1/(2A_n)) ∫ Φ_n(t)²/(t+(n+½)²) dν(t) ≤ 2/((2n+1)³A_n).

Question: for which (measure, basis) pairs is the Gram matrix a Hilbert or Cauchy matrix? Each such
pair transports a classical orthogonal family with its norm. (Cubic notes: not checked.)
