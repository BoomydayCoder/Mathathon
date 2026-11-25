# R7_gluing — how large can the gain M/C be in the GGS cyclotomic gluing framework?

Tags: [P] proved here, [L] literature (quoted from arXiv:2610.12619), [N] numerical (code in this directory), [H] heuristic, [O] open/conjecture, [R] refuted.
Code: `glue_search.py` (exact LP admissibility test, exhaustive search), `avoid_top.py` (index sets not containing lcm). Outputs: `out_noone.txt`, `out_one.txt`, `out_avoid.txt`.

## 1. Exact definitions from the paper [L]

* Sec 2.2: for S ⊂ N, c_S = Σ_{n∈S} c_n, c_n = φ(n); I = log Π_{{m,n}⊂N} |Res(Φ_m,Φ_n)|;
  I_inc(S) = Σ_{{m,n}⊂N, {m,n}∩S≠∅} log|Res(Φ_m,Φ_n)|. Apostol: |Res(Φ_m,Φ_n)| = ℓ^{φ(min)} if max/min is a positive power of a prime ℓ, else 1.
  covol(A) = e^I, covol(A∩K_S) = e^{I_inc(S)}; L = A^r ⊕ Z^a, d = Cr + a, 0 ≤ a < C, r = ⌊d/C⌋ (needs 1 ∈ N).
* Eq (6): A(S) := (c_S/2) log(C/c_S) + (c_S/C) I − I_inc(S) − Σ_{n∈S} c_n y_n, with Σ_{n∈N} c_n y_n = 0; A(N) = 0.
* Sec 3.1, eq (8)-(11): w(m,n) = log|Res|, ω_n = weighted degree, I_inc(S) = ½Σ_{n∈S} ω_n + ½W(S,S^c); with v_n = y_n − I/C + ω_n/(2c_n) (Σ c_n v_n = 0):
  **A(S) = (c_S/2) log(C/c_S) − ½ W(S,S^c) − Σ_{n∈S} c_n v_n**.   (11)
* Def 2.6 (admissible): 1 ∈ N, 2 ∉ N, C = Σ c_n, M = lcm(N); Σ c_n y_n = 0 and A(S) ≤ −δ for every nonempty proper S ⊊ N (δ > 0 absolute), together with
  I/C ≪ log C, |y_1| ≪ log C, r ≫ C^3 log d, C ≥ log d, log M ≪ log d.
* Prop 4.6: for admissible data, Δ_d ≫ d^2 (M/C) 2^{-d} (stochastic estimate applied with weights 1/M_supp(x); proper supports contribute e^{-cr} by Cor 2.8, so only full-support vectors with full orbit M matter). Thm 3.1: M/C ≫ log log d.
* GGS blocks D(q;N0) = {n | N0} ∪ {qN0}: C_loc = φ(qN0)(1 + P/(q−1)), P = N0/φ(N0) ≤ 5/4, primes of N0 ≫ q (log p ≥ 8(log q+1)). Hence in their construction C ≍ φ(M) and C > φ(M): gain < M/φ(M).

## 2. Reformulations

**2.1 [P] gain = 1/μ(N), any M.** μ_M(n) := φ(n)/M is a probability measure on the divisor lattice of M (Σ_{n|M} φ(n) = M), a product over primes (exponent e of p^a‖M has weight φ(p^e)/p^a). C/M = μ_M(N). So M/C = 1/μ(N); single field N={M} gives M/φ(M).

**2.2 [P] character/Parseval form.** tr(τ^k | Q(ζ_n)) = Ramanujan sum c_n(k), and (1/M)Σ_{k mod M} c_n(k)c_{n'}(k) = φ(n)δ_{nn'} for n,n' | M. Hence Σ_k tr(τ^k|A⊗Q)^2 = MC, i.e. M/C = Σ_{g∈⟨τ⟩} ρ(g)^2, ρ = tr/C. Checked numerically for N = {1,3,15}, {1,5,7,35,105}, {3,25,45} [N].

**2.3 [P] admissibility of the support part is an LP / cooperative-game core.** With u_n = c_n v_n, B(S) = (c_S/2)log(C/c_S) − ½W(S,S^c): need u(N)=0, u(S) ≥ B(S)+δ for all proper S. By LP duality (Bondareva–Shapley) this is feasible iff Σ_S λ_S (B(S)+δ) ≤ 0 for every balanced weight system (Σ_{S∋n} λ_S = 1 ∀n). Special cases (necessary conditions, the y/v terms cancel exactly because Σ_n c_n v_n = 0 and every n is covered once):
 * **Cut condition** (λ on {S,S^c}): W(S,S^c) ≥ (C/2) H_2(c_S/C) + 2δ, H_2 binary entropy (nats). Orchestrator's hint CONFIRMED (weights 1:1).
 * **Partition condition** (any partition 𝒫 of N): W(𝒫) (total weight between different parts) ≥ (C/2) H(ν_𝒫) + |𝒫|δ, ν = c/C. Singleton partition: **I ≥ (C/2) H(ν) + |N|δ**.
 * **Vertex condition** S={n}: ω_n ≥ (c_n/2)log(C/c_n) + ((C−c_n)/2) log(C/(C−c_n)) + 2δ. Squarefree, T ↔ n_T: Σ_{p∉T,T+p∈N} log p + Σ_{p∈T,T−p∈N} log p/(p−1) ≥ H_2(x_T)/(2x_T) ≥ ½log(1/x_T) + ½(1−x_T), x_T = c_T/C.
 NOTE the direction: every such condition is EASIER when C is small relative to the vertex weights, i.e. each one is a lower bound on the gain in terms of the shape of N; none bounds the gain from above by itself.

**2.4 [P] per-prime window, any M (odd q).** For a prime q | M let N_low(q) = {n∈N: v_q(n) < v_q(M)}, x_q = c(N_low)/C. Both sides nonempty when lcm(N)=M and x_q>0. Each n = m q^e in the top part has edges to m q^i (i<e) of total weight ≤ φ(m) q^{e−1} log q = φ(n) log q/(q−1). Hence W ≤ (1−x_q) C log q/(q−1) and the cut condition gives
  H_2(x_q) < 2(1−x_q) log q/(q−1)  ⟹  x_q < x_hi(q) < 2/q.
 (Same bound for the cut {q∤n}; in the squarefree case also W ≤ x_q C log q, giving x_q > x_lo(q) ≈ e/q^2.) This is R3's Prop 5 (m h(a/m) ≤ 2 min(a,b) log q) for the element τ^{M/q}.
 Proof of x_hi(q) < 2/q: at x=2/q, H_2 ≥ (2/q)(log q − log 2 + 1 − 2/q) > (2/q) log q (q−2)/(q−1) since 1 − log2 − 2/q + log q/(q−1) > 0 for q ≥ 3. Table [N]:

 | q | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 59 |
 |---|---|---|---|---|---|---|---|---|---|
 | x_hi(q) | 0.3909 | 0.2723 | 0.2092 | 0.1433 | 0.1239 | 0.0975 | 0.0882 | 0.0740 | 0.0305 |
 | q·x_hi | 1.173 | 1.362 | 1.465 | 1.576 | 1.610 | 1.658 | 1.676 | 1.703 | 1.797 |
 | q²·x_lo | 2.354 | 2.577 | 2.645 | 2.688 | 2.697 | 2.706 | 2.708 | 2.711 | 2.717 |
 (x_hi(2) = 1/2 exactly, attained with equality: this is why {1,2} is critical and 2 ∉ N.)

**2.5 [P] sharper "half-matching".** Squarefree: W_q = log q · c{T∈N: q∉T, T+q∈N}. So the matched fraction ρ_q of the mass missing q satisfies ρ_q ≥ H_2(x_q)/(2 x_q log q) (→ ≥ ½ for large q, x_q ≍ 1/q): at least about half of the mass of N missing q must have its q-up-neighbour in N. Summing over q: I/C ≥ ½ Σ_q H_2(x_q) ≥ ½ H(ν) (subadditivity), so the per-prime cuts together imply the singleton-partition condition.

## 3. Theorems on the gain

**Theorem A [P].** Let N satisfy the cut conditions (in particular: admissible), M = lcm(N) odd. If M ∈ N then C ≥ φ(M), so M/C ≤ M/φ(M). If M ∉ N then every n∈N lies in some N_low(q), so 1 ≤ Σ_{q|M} x_q < Σ_{q|M} x_hi(q) < 2 Σ_{q|M} 1/q.
 Corollary: **M/C > M/φ(M) is impossible unless Σ_{q|M} x_hi(q) > 1** (in particular unless Σ_{q|M} 1/q > 1/2). With the table: M needs ≥ 4 distinct prime factors, and with exactly 4 they must be {3,5,7,11} (Σ = 1.0157; {3,5,7,13} gives 0.9964).

**Theorem B [P] (weak universal bound, squarefree M).** With θ = φ(M)/M, μ(T) = θ Π_{p∉T} 1/(p−1) ≤ μ(N) for every T∈N; min ≤ mean under ν gives log(θ/μ(N)) ≤ Σ_p x_p log(p−1). Hence
  M/C ≤ (M/φ(M)) · Π_{p|M} (p−1)^{x_hi(p)} < (M/φ(M)) · exp(2 Σ_{p|M} (log p)/p).
 For M built from all primes ≤ y this is only (M/φ(M)) y^{2+o(1)}: NOT a log log no-go (the excess could a priori be a power of the largest prime).

**No proof of M/C ≤ K·M/φ(M) was found** [O]. Reason: all admissibility inequalities are of the type "entropy ≤ boundary weight" and get weaker as N gets lighter; the only upper bounds on the gain come from covering (lcm(N)=M) combined with the per-prime windows.

## 4. Exhaustive / random search [N]

Exact LP (maximise t s.t. u(S) − B(S) ≥ t, all 2^{|N|}−2 supports; admissible iff t* > 0), prefiltered by vertex, singleton-partition and all cut conditions. All index sets N ⊂ divisors(M) with lcm(N) = M:
 * M ∈ {15,21,35,45,75,105,135,165,195,225,315,525,1155} (up to 16 divisors, 64594 index sets for 1155): **every N with C < φ(M) fails**; best admissible gain = M/φ(M), attained only by N = {M} (no proper support). Next best: a·{1,p}-type sets, e.g. {105,1155} gain 2.1875 = 105/φ(105).
 * With 1 ∈ N required (as in Def 2.6): best gains are far below M/φ(M): M=105: 1.0825 (N={1,3,5,7,21,35,105}, t*=0.018); M=525: 1.1339 (N={1,5,7,25,35,105,175,525}, t*=0.003); M=1155: 1.1726 (N={1,5,7,11,35,55,77,231,385,1155}, t*=0.018) versus M/φ(M) = 2.19, 2.19, 2.41.
 * Index sets NOT containing M (lcm = M): M ∈ {105,225,315,525,945,1155,1365,2145} and, allowing 4 | n, M ∈ {60,120,420} (pools of ≤ 15 proper divisors, n ≡ 2 mod 4 excluded): **zero admissible sets**; in fact zero pass all cut conditions (only 3 sets for M=525 pass vertex + singleton-partition).
 * Random search, M = 15015 = 3·5·7·11·13: 2.6·10^6 random N avoiding M (|N| ≤ 13): 12 pass vertex+singleton conditions, 0 pass the cut conditions.
 * Sanity: {1,3}: t* = 0.0719 = ¼(log 3 − 2 log 1.5) as in GGS 3.2.1; {1,19,57} (GGS 3.2.2, b=3, p=19): t* = 0.223; full cube of 105: t* = 0.915; {1,3,5}: t* = −0.50.

## 5. Asymptotic candidate families that would beat M/φ(M) — all fail at leading order [H]

Write θ = φ(M)/M, relative mass m = μ(N)/θ (need m < 1), s = Σ 1/(p−1), Λ = Σ log p/(p−1) over the relevant primes; Mertens: Σ_{p≤y} log p/(p−1) = log y − γ + o(1), so primes in [X, X^a] have s ≈ log a, Λ ≈ (a−1) log X.
 1. Cube minus top, N = D(M)∖{M}: m = 1/θ − 1 < 1 needs θ > ½; vertex condition at M/p_max needs Λ ≥ ½ log p_max + ½ log(1/θ − 1) + ½, i.e. a ≥ 2, i.e. θ ≤ ½. Exactly borderline, lost by the constant ½. [R as an improvement]. (Numerics, `cube_minus_top.py` [N]: for M = product of odd primes ≤ y the minimal slacks of (level-1 vertex condition, per-prime cut, singleton-partition) are y=11: (−0.53, −0.062, −0.032); y=47: (−0.16, −0.036, +0.20); y=200: (+0.31, −0.018, +0.60); y=1000: (+0.89, −0.0044, +1.14); y=10^4: (+1.86, +0.0003, +2.05). The binding one is the q=3 cut, x_3 = 1/(3(1−θ)) ≤ x_hi(3) ⟺ θ ≤ 0.147. So for y ≳ 10^4 D(M)∖{M} passes all these necessary conditions; whether it is admissible is undecided, hence "lcm(N) ∈ N" is probably NOT forced in general [O]; but its gain 1/(1−θ) ≈ 1.14 is far below 1/θ ≈ 8.2.)
 2. Levels 1 and 2 with a graph G on the primes (N = {M/p} ∪ {M/pq : pq ∈ E}): vertex conditions allow s ∈ (log 2, 0.732) with m = s + s²/2 < 1, but the singleton-partition condition forces Σ_p (Λ_p − log p)/(p−1) ≳ 0, i.e. s ≥ 1 − o(1) and nearly complete G, so m ≈ 1.5 [R].
 3. Same plus pendant very large primes q (elements M/pq, M/q ∉ N), which repair the entropy condition cheaply (m ≈ 0.81–0.89): killed by the per-prime cut for p ∈ P_L (half-matching 2.5): needs s_L ≥ 1 + s_R [R].
 4. Idealised "J primes of size P, λ = J/P fixed" (orchestrator): windows j ≤ 2λ vs j > eλ exist only for j=1,2,3 with λ in (0.5,0.69), (1,1.15), (1.5,1.55); but real primes cannot realise this model (Σ_{p∈[P,cP]} 1/p → 0), and with Mertens weights the j=1 window closes (item 1) and j ≥ 2 top vertices M/(p1p2) with two large primes violate the vertex condition [H].
 Common mechanism: for a large prime p at least half of the ν-mass missing p must be matched upward along p, each matched pair multiplies mass by p−1, and Σ log p/(p−1) over primes below y is only log y − γ.

## 6. Verdict and what would be needed

 * [N]+[P partial] No admissible index set with M/C > M/φ(M) exists for any M tested, and Theorem A proves this whenever Σ_{q|M} x_hi(q) ≤ 1. [O] Conjecture: M/C ≤ K·M/φ(M) (probably K = 1) for all admissible data, so the cyclotomic gluing method is capped at d² log log d.
 * Since log M ≪ log d forces M/φ(M) ≤ (e^γ+o(1)) log log M ≪ log log d, a gain beyond log log d inside this framework would need N avoiding its lcm with μ(N) ≪ θ; Theorem B only excludes excess beyond exp(2Σ log p/p).
 * Modifications [H]: (i) non-cyclic abelian groups do not help (a faithful Q-irreducible representation of an abelian group forces it to be cyclic, and the Parseval identity M/C = Σ_g ρ(g)² shows the gain is the L²-mass of the character, bounded by fixed-space dimensions exactly through the windows of 2.4 = R3 Prop 5); (ii) one needs groups whose Q-irreducible constituents have dimension ≪ |G|/log log, i.e. non-commutative orders (fixed-point-free groups in division algebras give only bounded extra factors, cf. R3 Sec 6: SL_2(3)×C_N has g = 6N/φ(N)), or non-free actions whose lower-rank (R3 "super-stability") terms are controlled, which R3 Prop 5 restricts by the same entropy window.
