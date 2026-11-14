# Progress log (Claude)

## 2026-11-14 06:52 UTC
- Working on: gate (Zudilin exposition), first full draft done and referee-checked once.
- Since last entry: `exposition.md` (complete proof with expanded calculations, motivation, referee report §8, origin of the recurrence §7); `computations/zudilin_lemma1.py`, `zudilin_lemma2_recurrence.py`, `zudilin_ball_series.py`, `dn_bound_check.py` with saved outputs; `NOTEBOOK.md`. Fauzan preprint and both Lean repos retrieved and surveyed by a subagent (not yet written up).
- Status of claims: Zudilin's Lemmas 1–7 checked by hand and by script (symbolic for Lemmas 2 and 5). Proved by hand: e^{f(tau_0)}=(sqrt2-1)^4; tau_0 is the maximum. Numerical only: D_n<3^n (n<=10^5); uniqueness of the telescoping certificate (n<=8). Not proved here: a Chebyshev-type bound limsup D_n^{1/n}<3.2388 (the one external input).
- Notes for the other agent: research-chatgpt branch not yet checked at this entry. Two misprints in Zudilin's paper worth knowing: eq. (7) should read F_{n+1},F_n,F_{n-1}; in the last display of Lemma 4 the exponents should be 2/t^3+5n/t^4+2n^2/t^5. Also "trivial bound D_n<3^n" is not trivial: 4^n is not enough (64(sqrt2-1)^4=1.88>1); one needs Chebyshev/Hanson or PNT.

## 2026-11-14 07:10 UTC
- Working on: priority 3 (better proofs) and the Fauzan zeta(5) proof; a frontier subagent is testing a Hankel-determinant proof for zeta(3).
- Since last entry: `short-proof.md` (Apéry from one degree-4 identity + Casoratian, no Ball series; now fully elementary using the period-6 lcm bound); `recurrence-origin.md` (hand derivation of the recurrence and reading of its coefficients); `fauzan-status.md` (claims, Lean status with statements/axioms checked directly, series-language dictionary, why determinants); `dependency-map.md`; exposition §9 (Beukers–Zudilin dictionary via Mellin transform); scripts `short_proof_checks.py`, `fauzan_functional.py`, `dn_bound_check.py`.
- Status of claims: short proof: every step proved by hand, exact identities script-checked (n<=60). Fauzan functional series form: exact check on 23 test functions for s=3,5,7 (not a general proof beyond Fauzan's (3.1) for s=5). Inhomogeneous recurrence L[G]_n=-8(2n+1) for G_n=sum R_n(t): derived, not yet numerically checked. Hankel s-dependence: heuristic + small-K numerics only.
- Notes for the other agent: (1) Your branch is `codex/independent-research`, not `research-chatgpt`; we read it there. (2) Checked and found correct: period-6 lcm bound (denominator-design.md §1; re-derived by hand: table 0,1,1,0,0,1, a=log2/3+log3/2, c=1.17053); forced-recurrence formula for m=2 (numerically n<=7); odd-s moment and pole formulas (numerically and exactly, s=3,5,7); Beukers triple integral = -sum R_n'(t) for n=1,2. No errors found. (3) Result you can use: Zudilin's Lemmas 4-7 (Ball's series, degree-6 certificate) are unnecessary: from his Lemma 2, L[u]=0 and L[F]=0, so u_n F_{n-1} - u_{n-1} F_n = 24/n^3, and |F_n|<=2 gives F_n = u_n sum_{m>n} 24/(m^3 u_m u_{m-1}) > 0; see research/claude/short-proof.md on research-claude. With your lcm bound this needs u_n >= c*rho^n with rho > 2^{3/2} 3^{9/4} = 33.50, obtained by bootstrapping the recurrence (ratio -> 34 - 1/30). (4) FLSZ arXiv:2611.16156 §3 has an elementary Pochhammer-basis denominator bound for Hankel determinants; relevant to your odd-arithmetic-barrier.md.

## 2026-11-14 08:00 UTC
- Working on: Hankel-determinant method for odd zeta values (frontier subagent still running); checking the other agent's new results.
- Since last entry: `computations/companion_zeta2_forms.py` (inhomogeneous recurrence L[G]_n = -8(2n+1) for G_n = sum R_n(t), exact n<=30); in-progress `computations/hankel_*.py` and outputs from the frontier subagent (not yet refereed by the orchestrator).
- Status of claims: nothing new proved since last entry. Numerical only (unrefereed, from in-progress outputs): for the family G=[mu_X(t^{i+j}/D_K)], h=K, the measured (cost + log det)/K^2 at K=80 is -0.90 (s=3), +0.07 (s=5), +1.07 (s=7); with Fauzan's parameters -0.99 (s=3), -0.13 (s=5), +0.77 (s=7, K=40). Heuristic check by the orchestrator: 1x1 (linear-form) version has normalised exponent >= 0, so no proof at h=1.
- Notes for the other agent: Checked from your research-results.md (branch codex/independent-research), with independent exact/high-precision code: (3) E_{n,0,m} asymptotic, m=2,3, n<=1280: ratio = 1 + c/log n, c≈0.56 (m=2), 1.87 (m=3): consistent. (5) E_{n,n,m} ~ ... 16^{-n}: ratio = 1 - c/log n, consistent. (6) integrality exact for n<=40, m=2,3: correct; moreover the TRUE reduced denominators of beta_{n,n,m} grow like about e^{4m n} (7.8 n for m=2 at n=80), so the 16^{-n} decay cannot win: this is not an artefact of a weak bound. (7) v_p = -3 for all odd p<=61 except v_17 = -2 (m=2): correct. (8) ratio to F_n at the cut n/sqrt2 tends to -zeta(2), 2 zeta(4) with n^{-1/2} error: consistent. No errors found. Not checked: §4 (9)-(11), §5 (12)-(14), §6.

## 2026-11-14 09:40 UTC
- Working on: Hankel-determinant proof of Apéry's theorem; summary deliverables.
- Since last entry: `hankel-zeta3.md` (complete proof draft: statement, positivity, two-sided Newton basis, Lemmas L1-L8, analytic bound via Hadamard + inequality (H)); `computations/hankel_*.py` with outputs (exact checks K<=100; referee checks); `next-steps.md`; `SUMMARY.md`; GOALS ticked.
- Status of claims: Hankel proof for zeta(3): complete draft, three AI passes, exact verification K<=100, one interval-arithmetic step; not human-refereed, to be treated as a draft. s-dependence (zeta(5) marginal, zeta(7) fails by >0.7 K^2): numerical. Simpler zeta(5) proof: not achieved.
- Notes for the other agent: (1) Your open item "Gram/moment-matrix comparison for the positive odd-zeta moment matrices" (research-results.md §8): for s=3 we have a full proof draft that det[mu_X(t^{i+j}/D_K)]_{i,j<K}, normalised by an explicit Phi_K with log Phi_K = 6 int_{K/2}^K theta + O(K^{3/2}), is an integer polynomial that is positive and < exp(-(2.735 - 2.25 kappa)K^2) at zeta(3); see research/claude/hankel-zeta3.md on research-claude. Key device: rows prod_{j<=i}(t+j^2), columns prod_{j>K-i'}(t+j^2); entries are divided differences of m_j or mu of a product; v_p(entry) >= -3 - 2 N_p(i,i') for p > sqrt(2K)+3 by a reflection q -> -1-q argument. Your period-6 lcm bound is exactly what makes it PNT-free; we cite your file. (2) For s>=5 the same basis gives cost about s(s+4)/(2(s+2)) K^2 against an s-independent decay (9/2) log 2 K^2: +0.10 (s=5), +1.16 (s=7); consistent with your K<=6 table. (3) No errors found in the claims of yours that we checked (listed in SUMMARY.md).

## 2026-11-14 08:45 UTC
- Working on: second frontier attempt (prime-independent basis for Fauzan's zeta(5) determinant), in progress.
- Since last entry: wording fixes in `hankel-zeta3.md`. Correction: the two previous entries are stamped about an hour ahead of the container clock (they were written at about 07:50 and 08:15 UTC).
- Status of claims: unchanged.
- Notes for the other agent: none new.

## 2026-11-14 09:05 UTC
- Working on: wrap-up of the second frontier run.
- Since last entry: `zeta5-simplification.md`; `computations/hankel5_exact.py`, `hankel5_model.py` and outputs; addenda in `SUMMARY.md`, `next-steps.md`.
- Status of claims: negative-leaning, numerical: a single prime-independent basis does not reproduce Fauzan's arithmetic within the available margin; a small family might. Proposition H (entrywise bound for Hermite staircase bases) is the subagent's proof, not re-derived by the orchestrator. Possible improvement of Fauzan's (4.14) by about 0.02 K^2: exact evidence at 7 primes, unproved.
- Notes for the other agent: if you look at Fauzan's §4: his outer-range exponent (4.14) seems lossy for 2.47 < K/p < 3 (e.g. K=120, p=41: (4.14) gives -410, true v_p = -332); see research/claude/zeta5-simplification.md item 5.

## 2026-11-14 10:00 UTC
- Working on: low-burn mode (one modest item per hourly wake).
- Since last entry: `computations/equilibrium_s3.py` + output: exact solution of the equilibrium problem behind the s=3 Hankel determinant (support [0,3/4], explicit density, l = 2+6 log 2, E* = (9/2) log 2); remark added to `hankel-zeta3.md` §5.
- Status of claims: equilibrium problem: proved (subagent; 39 numerical cross-checks; orchestrator re-derived only the support equation). Link "log Delta_K <= -E* K^2 + o(K^2)" not written, so Theorem 1(d) is unchanged.
- Notes for the other agent: (1) Checked from your commits 315fa08..13a9247: resistor-chain identities for the Apéry recurrence (research-results.md §11): Casoratian 6/(n+1)^3, sum = zeta(3)/6 to 50 digits, quadratic-form identity on random f, Green kernel exactly for N=6: correct. Fuss-Catalan prime-power valuation statement (problem137 discovery-guide.md): brute force over 1320 cases, no failures; the bijective proof itself was not audited. (2) Your one-interval equilibrium simplification for Fauzan's field and our exact solution for the s=3, W=1 field use the same method (charge density of the field, closed-form sweep); for W=1 the answer is elementary: sigma(x) = (1/pi)[3 artanh(s/b) - artanh(s/c)], s = sqrt(b^2-x^2), b=3/4, c=5/4, E* = (9/2) log 2. (3) We have not checked your end-to-end zeta(5) derivation (fauzan-simplified-proof-audit.md). Our exact experiments suggest Fauzan's (4.14) is lossy for 2.47 < K/p < 3 (research/claude/zeta5-simplification.md item 5), which may help your normaliser. (4) Our work is on branch research-claude under research/claude/ (start with SUMMARY.md); there is no HANDOFF.md.

## 2026-11-14 11:00 UTC
- Working on: low-burn mode.
- Since last entry: `hankel-energy-bound.md` (subagent's write-up of log Delta_K(zeta(3)) <= -(9/2) log 2 K^2 + 5K log K + O(K)); pointer added in `hankel-zeta3.md` §5.
- Status of claims: `hankel-energy-bound.md` is an unrefereed draft (one step sketched; negativity lemma cited from Saff-Totik). Its bound is numerically consistent with exact values at K=20,60,100. Theorem 1 of `hankel-zeta3.md` unchanged.
- Notes for the other agent: Checked from commit 7ebb5d6 with independent code reading only your JSON data: the superfactorial normaliser arithmetic recomputes correctly (G_2 = 111/160, beta = 9377/80000, C = 148, H_* >= T on 68,800 grid points of [0,43], A = 1.37852, A + patch = 1.37858, K^2 log K cancellation confirmed numerically). IMPORTANT caveat: A = 1.3786 exceeds Fauzan's A_200 = 1.3496, so with Fauzan's proved real constant U = -1.3670 your PNT-free route has margin +0.0116 and FAILS; it closes only with your own real bound -1.38 (margin -0.0014), i.e. it stands or falls with your equilibrium-energy certificate and the "finite determinant transfer" in fauzan-simplified-proof-audit.md §7, which we did not audit. Also unchecked by us: that your T equals Fauzan's local cost (R_0, d), and his local lemmas. Our own energy solver gave -1.3855 for Fauzan's field (numerical), which is compatible with -1.38 but leaves only ~0.007 of room. Reciprocal Apéry model (CT Lambda^n = A_n, n<=5; integral = Beukers; min Lambda = (1+sqrt2)^4): correct.

## 2026-11-14 11:50 UTC
- Working on: low-burn mode; refereeing `hankel-energy-bound.md`.
- Since last entry: orchestrator refereed `hankel-energy-bound.md` line by line (Steps 1-13: scaling exponent, Riemann bound, identity I(mu)+l = E*, density bound, smoothing error integral, energy expansion, decay reserve, final bookkeeping: all correct). `hankel-zeta3.md` now states Theorem 1(d') with its dependencies; SUMMARY addendum.
- Status of claims: (d') log Delta_K(zeta(3)) <= -(9/2) log 2 K^2 + 5K log K + O(K): proved given (I3) Euler-Lagrange inequality off the support (subagent's proof; not re-derived by the orchestrator) and (I4) negativity lemma (cited). The Hankel proof of zeta(3) irrationality now has two independent analytic routes.
- Notes for the other agent: none new (a check of your commits a70ae89, 17a2f78 is running).

## 2026-11-14 12:05 UTC
- Working on: low-burn mode.
- Since last entry: check of the other agent's commits a70ae89 and 17a2f78.
- Status of claims: unchanged on our side.
- Notes for the other agent: Checked independently: (1) integer-depletion-cumulants.md §3 (least uniform denominator prod_{p|l} p^{1+v_p(s-1)}): exact lcm over e<160 for odd s=3..31, l=2..12 (165 cases) all equal the stated value; sharpness argument re-derived: correct. (2) coprime-depletion-primitive-obstruction.md eqs (2),(3),(6): X-coefficient (2^{s-1}-1)(3^{s-1}-1) r^{s-1}, leading coefficient with sign, monic reduction, the integral 4/3: correct; NOT checked: the local valuation (4) and the Jensen lower bound, so the obstruction theorem as a whole is unverified by us. No errors found. We note that you claim no new irrationality theorem; your PNT-free zeta(5) derivation remains conditional on your own real bound below -1.38 (see our 11:00 entry).

## 2026-11-14 12:55 UTC
- Working on: low-burn mode; closing the remaining unverified input of the sharp Hankel bound.
- Since last entry: orchestrator's own numerical and hand check of the Euler-Lagrange relations (NOTEBOOK, "Independently proved observations").
- Status of claims: `hankel-energy-bound.md` input (I3) now checked by the orchestrator (derivative formula numerically at 5 points to 8 digits, its negativity by hand); (I4) remains a cited classical lemma. Resume cadence reduced to every 3 hours.
- Notes for the other agent: your commits bafa595, db8494b not yet checked by us.

## 2026-11-14 14:35 UTC
- Working on: full-speed mode again (owner's instruction): three frontier subagents (zeta(5) margin theorem; Hankel method for Catalan/other constants; new constructions for zeta(7)), literature searches, check of the other agent's newest commits.
- Since last entry: `literature-status.md`; SUMMARY addendum.
- Status of claims: novelty downgrade: determinant-method proofs implying zeta(3) irrational already exist post-Fauzan (Tang; Lai arXiv:2610.06838), and Catalan's constant was claimed by OpenAI (24 Sep 2026, Lean). Our Hankel zeta(3) note is at best a simpler construction.
- Notes for the other agent: literature you may not have: Lai arXiv:2610.06838 (1, zeta(3), zeta(4) independent; lists Tang's 1, zeta(2), zeta(3), Kilani's p-adic results, OpenAI's Catalan proof); Niedbala Giraudin arXiv:2611.01708 (exponent of zeta(5) <= 171); Wachs arXiv:2609.22339 (gap in Sun's Catalan claim: the prime 2). No zeta(7) result found on arXiv.

## 2026-11-14 14:55 UTC
- Working on: full-speed mode; three frontier subagents still running (zeta(5) margin, Catalan/other constants, zeta(7) constructions).
- Since last entry: `apery-elementary.md` (one self-contained elementary proof of Apéry's theorem: linear forms, one identity verified by hand coefficient by coefficient, Casoratian, growth of u_n, period-6 lcm bound, conclusion).
- Status of claims: `apery-elementary.md`: every step hand-checked by the orchestrator; AI-checked only.
- Notes for the other agent: Refereed your short-elementary-hybrid-proof.md (commit bafa595): partial-fraction arithmetic exact for n<=7; I_n = Beukers triple integral numerically for n=1,2 and the Jacobian by hand; the AM-GM bound Lambda >= (3+2 sqrt2)^2 by hand and by sampling; the integer inequality for the lcm constant: all correct, proof complete as written. It is a different route from ours (yours: Beukers integral + AM-GM, no recurrence; ours in research/claude/apery-elementary.md: telescoping identity + Casoratian, no integrals); they share Lemma 1 and your period-6 lcm bound. Your rank-one approximants (rank-one-zeta-subspace-approximants.md): linearity in X, one-sidedness, monotonicity (s=3: h<=7; s=5,7: h<=5), factorial bound, and the Apéry match q_{h,3} = (5B_h - B_{h-1})/(5A_h - A_{h-1}) for h<=7 recomputed independently: correct; as you say, no irrationality consequence (denominator x error ~ 6e28 for s=5, h=5). Not checked: your second zeta(3) proof (resolvent/Padé), the log 2 proof, and the notes of db8494b, a3fa456, 1434d60 beyond their labels.

## 2026-11-14 15:25 UTC
- Working on: frontier lines; three subagents still running (Catalan/other constants, new zeta(7) constructions, referee of Anand's zeta(7) preprint).
- Since last entry: `zeta5-margin.md` (Theorem W + new constants + zeta(7) model no-go); `computations/zeta5margin_*.py` and outputs; `literature-status.md` update; `apery-elementary.md` now has motivation for every step.
- Status of claims: Theorem W: proved (subagent), re-derived by the orchestrator, exact verification K=80..200, matches a proved Lean statement. New constants A_200 = 1.329961 etc.: exact rational computation (computer-assisted). zeta(7) no-go: numerical/model.
- Notes for the other agent: (1) For your PNT-free zeta(5) normaliser: Fauzan's outer bound (4.14) can be replaced on K/3 < p < (69/170)K by the class-frame (water-filling) bound, which is valid for all p >= 5 with p^2 > 104n+6; this lowers the cost constant by exactly 3203/163200 = 0.019626 (see research/claude/zeta5-margin.md). With it your A = 1.3786 would become about 1.3590 < 1.3670, i.e. your route would close even with Fauzan's proved real constant — worth checking against your T(x). (2) A Zenodo preprint (P. Anand, 22920911) claims zeta(7) irrational with the same architecture and D_N^8; our exact computations (true totals +0.85 K^2 at K=80,120 for s=7) contradict it; a referee check is running. (3) Zudilin arXiv:1507.05697 §6 already records a Hankel-determinant criterion applying to Beukers' zeta(3) forms.

## 2026-11-14 15:45 UTC
- Working on: Catalan/other constants (subagent running); referee of Anand's zeta(7) preprint (subagent running).
- Since last entry: `zeta7-status.md`; `computations/zeta7_*.py` and outputs.
- Status of claims: all zeta(7) statements are negative and numerical (exact finite computations; lemmas L1-L6 by the subagent, not re-derived by the orchestrator).
- Notes for the other agent: zeta(7): 420 exact determinants in 8 families (numerators, mixtures zeta(7)+-zeta(5), other node sets, higher-order poles, one-sided derivative series, your signed FLSZ-type model with weights (45,-9,1), minors) all fail; best scale-free ratio lc-cost/decay = 1.218. Your signed single-zeta(7) model: exact Hankel determinants up to n=14 give cost/decay about 1.9-2.0. Conjecture C1 in research/claude/zeta7-status.md (sum over roots of max(0,-v_p) >= s * #nodes >= p) would turn your obstruction for s>=13 into one for s>=7 in the W=1 family.

## 2026-11-14 16:00 UTC
- Working on: Catalan/other constants (subagent still running).
- Since last entry: `anand-referee.md`; `computations/anand_zeta7_*.py` and outputs.
- Status of claims: the preprint's main estimate is impossible: elementary contradiction with integrality of the leading coefficient (logic checked by the orchestrator; constant 3.65 not recomputed), elementary lower bound for all n, exact computation n<=4. Error located: I_out printed -15.2819, exact value from the preprint's own definitions +453803/288000.
- Notes for the other agent: the Zenodo zeta(7) preprint (22920911) should not be relied on: see research/claude/anand-referee.md. Its local bound (59) is correct at all 66 applicable primes for K<=160; the error is the value of I_out (Table 5).

## 2026-11-14 16:20 UTC
- Working on: wrap-up of the afternoon's frontier runs; all subagents have reported.
- Since last entry: `catalan-status.md`; `computations/catalan_*.py` and outputs; SUMMARY state table; next-steps re-ranked.
- Status of claims: Catalan: functionals verified numerically (156 checks); failure of the method is exact finite-K evidence (K<=160), not an impossibility proof.
- Notes for the other agent: positive rational-moment functionals whose pole values involve only G exist (research/claude/catalan-status.md §1), but their Hankel determinants give totals +0.38..+0.52 per K^2 (K=40..160): decay 2 log 3 against cost -> 8/3 plus a K^2-order 2-adic term. The same family gives negative totals for zeta(2), log 2, pi. This may be useful for your log 2 / Catalan notes.

## 2026-11-14 15:20 UTC (container clock)
- Working on: second afternoon wave: Catalan mixed-node constructions; turning the zeta(7) no-go into theorems; independent recomputation of the zeta(5) margin constants (all running).
- Since last entry: nothing new recorded. Clock note: the entries above stamped 14:35-16:20 were written between about 14:20 and 15:10 by the container clock; my stamps ran ahead. Housekeeping: a stray scratch file `nmat.py` written by a subagent at the repository root (outside research/claude/, never committed) was moved to the scratch directory.
- Status of claims: unchanged.
- Notes for the other agent: your commits 1e747e1 and 4191392 are not yet checked by us.

## 2026-11-14 15:35 UTC
- Working on: Catalan mixed nodes and the independent zeta(5) constants check (both running).
- Since last entry: `zeta7-status.md` §5 (Theorems A, B, N); `computations/nogo_*.py` and outputs; remark in `hankel-zeta3.md`.
- Status of claims: Theorem B (lower bound, exact decay rate (9/2) log 2): re-derived by the orchestrator. Theorem A: subagent's proof, verified on exact Newton polygons, not re-derived. Theorem N (no-go for the W=1 family, odd s>=15): follows from A, B and PNT. s=7..13: open.
- Notes for the other agent: your "original-outer-monic-obstruction" (s>=29, per your ledger) can be compared with ours: for W=1 we get a no-go for all odd s>=15 from (i) at least n_p - m roots of p-adic valuation <= -s for primes K/(M+1) < p <= K, M = max{m: 2m^2-m+1 < s}, and (ii) the exact decay rate (9/2) log 2 with a Jensen lower bound using the explicit equilibrium density sigma(x) = (1/pi)[3 artanh(s/b) - artanh(s/c)], b=3/4, c=5/4. See research/claude/zeta7-status.md §5. An exceptional-prime phenomenon to be aware of: s=7, K=120, p=43 | num H_2^(7) = 129 breaks the root count (cost inequality still holds).

## 2026-11-14 15:50 UTC
- Working on: waiting for the independent recomputation of the zeta(5) margin constants.
- Since last entry: `catalan-status.md` §6 (mixed nodes: negative, with cancellation and obstruction lemmas); `computations/catalan2_*.py` and outputs.
- Status of claims: Catalan mixed-node constructions fail (exact determinants up to h=64; totals +2.2 to +8 per h^2); the two lemmas are the subagent's, not re-derived.
- Notes for the other agent: for Catalan with positive weights: cancelling pi^2 in a Hankel matrix of size h needs at least 2h pi^2-carrying nodes, and a non-negative pi^2-free test function must decay exactly like 1/t (research/claude/catalan-status.md §6). Relevant to your log2-endpoint-catalan-arithmetic note.

## 2026-11-14 16:00 UTC
- Working on: all subagents of the second wave have reported.
- Since last entry: `zeta5-margin.md` corrected (two table rows added, U not rounded) and §6 added; `computations/zeta5margin_independent_check.py` and output.
- Status of claims: Theorem W and the new constants confirmed by an independent recomputation for n=1,2,3 (different method, no shared code). Exact margin at M=200: -3488452892524542562969/94195881588024216000000 = -0.037034.
- Notes for the other agent: none new.

## 2026-11-14 16:20 UTC
- Working on: third wave (re-optimising the zeta(5) construction with Theorem W and exact energies; extending the no-go to s=7); both running.
- Since last entry: `SUMMARY.md` rewritten as one consolidated document; check of the other agent's commits 1e747e1..50f6ed1.
- Status of claims: unchanged on our side.
- Notes for the other agent: Checked independently (exact arithmetic): (1) q_{h,3} = (5B_h - B_{h-1})/(5A_h - A_{h-1}) as the unique root of the 2h x 2h interpolation determinant for h=1..6; the increment formula and the bracket B_h/A_h < zeta(3) < q_h (h=2..11); error rate (17+12 sqrt2)^{-2h}: correct. Caveat: your measure mu_3/D_{2h} was not rebuilt from scratch. (2) balanced Rodrigues integrality F_{m,n}/n! in Z[x^2] and binom((m+1)n+m, n) p_{m,n} in Z[t], endpoint values, for m<=3, n<=5: correct (the binomial is not always sharp). No errors found. Your zeta(5) status is unchanged (cost 1.3786 against Fauzan's 1.3670): our Theorem W (research/claude/zeta5-margin.md) lowers Fauzan's cost constant by 3203/163200 and may close your route with his proved real constant. None of your new notes overlaps or contradicts our results on the zeta(5) margin, the s>=15 no-go, the exact decay rate, or the refutation of the Zenodo zeta(7) preprint.
