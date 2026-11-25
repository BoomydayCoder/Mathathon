# Referee notes on R3_non_free (07:16-07:36 UTC). Code: ref_check.py (ref_output.txt), nontensor_check.py (nontensor_output.txt)

Tags: [P] proved here, [N] numerical, [H] heuristic, [R] refuted, [G] gap.

1. Thm 1 (mean value with lower-rank terms): re-derived, Jacobian (c_U/|d_K|^{s/2})^k and covol(Lambda_U)=c_U^k correct. Conditional only on
   Siegel-Weil for K-independent s-tuples (s<=k-1), constant 1 for free O_K^k by strong approximation for SL_k. Non-free W_Z cap U: same formula
   expected from the adelic form (not written). [confirmed, conditional on literature]
2. Prop 2: re-derived: T_U = V^f Gamma(n/2+1)^f/Gamma(fn/2+1) tau^{-fn}; Stirling gives the stated form. [P][N ratios .997-.9997]
   Borderline: T_U ~ n^{a f-(1-f)/2} for V=n^a: f=1/2: n^{1/4} (V=n), n^{3/4} (V=n^2). So borderline terms are polynomially LARGE at Klartag scale
   unless f <= 1/5 (V=n^2) resp. f <= 1/3 (V=n).
3. Prop 3(a): exact relation V = T_U^{1/f} (f tau^2)^{n/2} (pi f n)^{1/(2f)} (pi n)^{-1/2}. The factor called poly(n) is (pi f n)^{1/(2f)}:
   NOT polynomial if f -> 0 ([N]: 10^{22.6} for n=1800, f=.05; 10^{37} for n=12000, f=.04). Correct statement:
   V(W_Z (x) L) <= [ (pi f n) (2^{2kappa+o(1)} f tau_U^2)^{fn} ]^{1/(2f)} (pi n)^{-1/2}; exponentially small iff (pi f n)(2^{2 kappa} f tau^2)^{fn} -> 0.
4. Prop 4, Prop 5, Heisenberg (Blichfeldt), Prop 6 (lines, subfields), Prop 7, counterexample SL_2(5)xC_7, Prop 8: re-derived, no error found.
   Entropy windows reproduced: q=3 [0.7084,1.2835], q=5 [0.4597,1.4970], q=7 [0.3424,1.5874], q=11 [0.2272,1.6723]; max |rho| = 0.121 (q=5).
   SL_2(p), p=1 mod 4, p<2000: only (p,q,N)=(17,3,3) violates. p=3 mod 4 deficiencies -> log 2 and 1-log 2 confirmed.
5. [R] SCOPE ERROR in Cor 4.1 / verdict (1) "Delta=Q, r>=2 always fails": true only for TENSOR-FORM lattices W_Z (x)_Z L (modules over the maximal
   order M_r(Z)). G-invariant lattices in W^k are modules over the image of Z[G], in general a non-maximal order, and need not be tensor products.
   Counterexamples [P, N: nontensor_check.py]:
   (a) G=D_8, W=Q^2: Lam_phi = {(x,y) in L1^2 : y = phi(x) mod 2 L1}, phi a free involution of L1/2L1 (k even). covol = 2^k c1^2.
       Lines: (1,0),(0,1): slice 2L1, tau^2 = 2; (1,+-1): slice index 2^{k/2}, tau^2 = 2; (a,b): tau^2 = 2(a^2+b^2) (a-b odd) or a^2+b^2 (a-b even).
       So f tau^2 >= 1 on ALL lines, =1 on 4 lines (borderline), versus 0.5 for Z^2 (x) L and the claimed bound gamma_2/2 = 0.577.
   (b) G=D_12, W=Q^2 (A_2): Lam_C = {X in A_2 (x) L1 : X mod theta in C}, theta = sqrt(-3), C a k/2-dim subspace of L1/3L1.
       covol = 3^{k/2} 3^{k/2} c1^2. Line through a+b w: tau^2 = 2N(a+bw) if theta does not divide, 2N/3 if it does. So tau^2 >= 2 on all lines,
       equality exactly on the 6 root lines (the fixed lines of the 6 reflections) -> borderline, not exponentially failing.
   Mean value: for L1 Haar-random with the level structure, rank-2 pairs give exactly V (Siegel-Rogers on SL_k(R)/Gamma_level, k>=3), rank-1 pairs
   give the line terms. So E# = V + O(V^{1/2} n^{-1/4}): no exponential term.
   These families are the complex-conjugation-symmetric subfamilies of the free mu_4 / mu_6 families (Lam is a Z[i]- resp. Z[w]-module), with
   max glue at the ramified prime; their gain is |G|/dim_Q A = 2 resp. 3 = gain of the free family. So the STRUCTURAL conclusion (no extra gain)
   survives, but the no-go must be restated. What does hold for every G-invariant lattice (any Z[G]-lattice Lam in W^k (x) R, any Delta):
   Prop 5 applied to Lam itself (it only uses the Z[<g>]-structure): an involution with a != b (e.g. any reflection: Weyl groups, S_n), or a
   prime-order element outside the window, forces an exponentially large term. D_4/2^{1+4}, E_8/2^{1+6}: involutions have a=b; non-tensor
   lattices are NOT excluded by the notes.
6. Prop 6 imprecision: for K=Q(sqrt-15), E=Q(sqrt-3,sqrt5): Q-span of G=mu_6 x| C_2 has dimension 4 (not 8), W=E is Q-reducible, sum rho^2 = 3,
   while |mu_E|/[E:Q] = 1.5 [N]. The 'gain' 1.5 is |G|/dim_Q(K-span) = |G| * (#parameters of K-Hermitian forms)/(n^2/2); fine as an effective gain
   of the family, but it is not the identity of Prop 7.
7. Consistency of the heuristic 'first moment gives V ~ |G|': reproduces Venkatesh (r=1,k=2: V<m, n=2phi(m)). Does NOT reproduce Vance
   (|G|=24 gives V<24 only); Rogers/Vance averaging over the k Delta-successive minima gives V ~ |G| k/e (2n/e Rogers, 24(n/4)/e = 6n/e Vance).
   So the benchmark should be |G| max(1,k/e). For the SL_2(p), p=1 mod 4 candidate (k~p/4, |G|~p^3, n~p^2/2) this would be ~ n^2/e, not (2n)^{3/2} [H].
8. Candidate SL_2(p), p=1 mod 4: 'passes all tests' is correct but the tests are weak. Line test must hold for ALL lattice vectors with the norm
   form N h(w)=h_1(w)h_2(w) (two real places of Q(sqrt p)), not only for minimal vectors. p=13 (Leech, r=3, c_W=1): need
   covol(w O_Delta) = (h_1 h_2)^2 p^2/16 >= 81, i.e. h_1 h_2 >= 36/13 with h_1+h_2=|w|^2. If h(w) in d^{-1}: h = N + j/sqrt13 (|w|^2=2N):
   norm 4 needs |j|<=4 (|j|<=7 possible), norm 6 |j|<=9 (10 possible), norm 8 |j|<=13 (14 possible), norm 10 |j|<=17 (18 possible).
   Volume of the norm body mod units: 12 omega_12^2 log(eta) Y^6 = 25.6 Y^6 -> ~1.1e4 expected violating classes at Y=36/13 for a random
   unimodular lattice [H]; the Eisenstein part of the Hilbert theta series is positive at all these nu. So p=13 is expected to FAIL on lines;
   generally the Gaussian count is (2 pi e/(2 sqrt p))^{p-1}: >1 for p<73. Not computed (needs explicit SL_2(13) action) [open].
