# Research notebook (Claude)

## Established results (from sources, checked by us)
- Zudilin, arXiv:math/0202159: retrieved, read in full, checked line by line. Proof correct. See `exposition.md` (§8 is the referee report: two misprints, one unproved-but-true assertion "D_n<3^n is trivial", two unproved steps in Lemma 4 which we prove by hand).
- Script verification (all in `computations/`, outputs saved as `output_*.txt`):
  - `zudilin_lemma1.py`: 93/93 PASS (partial fractions, A2k, sum A1k=0, u_n=2a_n, D_n^3 v_n integral, F_n numerics n<=6).
  - `zudilin_lemma2_recurrence.py`: 48/48 PASS (Lemma 2 symbolic; S_n'(1)=0; L[u]=L[v]=0 n<=10; uniqueness of the certificate; footnote telescoping).
  - `zudilin_ball_series.py`: 93/95; the two FAILs are errors in the orchestrator's brief, not in the mathematics (a mis-transcribed factor (tau+1/2), and a constant claim stated without the (n+1)^2 scaling). Lemma 5 symbolic PASS; F~_n = F_n exact for n<=8 via partial fractions.
  - `dn_bound_check.py`: D_n<3^n for n<=10^5 (numerical evidence only).
- Fauzan preprint (doi:10.5281/zenodo.22826419) retrieved (32 pp). Both Lean repos cloned. Subagent survey received; NOT yet checked by the orchestrator against the files. Headline (unverified by me): Hankel determinant of size 37n built from a functional mu_X; positivity at X=zeta(5) from a positive weight; arithmetic via p-adic valuations + PNT; decay exp(-139 n^2/5); margin about 1.3%. danromik repo: `Irrational zeta5` with PNT as the one axiom; mo271 repo: PNT imported as a theorem from PrimeNumberTheoremAnd, final theorem about `riemannZeta 5`, weaker constants.

## Independently proved observations (ours)
- Exact hand evaluation e^{f(tau_0)} = (sqrt2-1)^4 in Lemma 4 (exposition §3.2 Step 5).
- A1k = 2 A2k (2H_k - H_{n+k} - H_{n-k}) (exposition §1.2).
- The cubes (n+1)^3, n^3 in Apéry's recurrence come from evaluating the cleared telescoping identity at t=-n-1 and t=n (exposition §7.2).
- The certificate/recurrence pair of Lemma 2 is unique up to scaling within the ansatz S=s(t)R_n(t), deg s<=6 (computed for n=1..8; not proved for general n).

- `short-proof.md`: Apéry's theorem from Lemma 1 + the single degree-4 identity + Casoratian; no Ball series. All steps checked by hand; `short_proof_checks.py` 15/15 PASS (Casoratian exact n<=60, bounds).
- `recurrence-origin.md`: complete hand solution of the telescoping linear system; 17n^2+17n+5 = 4(2n+1)^2 + n(n+1) + 1.
- By-product (derivation only): G_n = sum_{t>=1} R_n(t) (forms in 1, zeta(2)) satisfies L[G]_n = -8(2n+1). Numerical check pending.
- Fauzan's functional for every odd s: mu_X(R) = tau_X(x^s R(-x^2)), and mu_{zeta(s)}(R) = -(1/(s-1)!) sum'_{x>=1} g^{(s-1)}(x) + boundary term, g = x^s R(-x^2). `fauzan_functional.py` 34/34 PASS (exact for s=3,5,7 on 23 test functions). The all-s moment/pole formulas were found first by the other agent (codex/problem131/odd-moment-generalization.md); we confirmed them independently.

## Proof drafts
- `hankel-zeta3.md`: Hankel-determinant proof of zeta(3) irrationality. Complete draft. Passes: frontier subagent (construction, L1-L8, sketch of (d)); orchestrator (L1-L7 re-derived by hand); referee subagent (L8 and (d) re-proved from scratch, uniform argument for (d)); orchestrator (checked the referee's proofs line by line). Computer-assisted step: inequality (H) on [1,60]. External: Hermite's formula; a Chebyshev-type bound with kappa < 1.2155 (period-6 bound suffices). Exact checks: K<=100. NOT human-refereed.

## Conjectures
- log Delta_K(zeta(3)) = -(9/2 log 2 + o(1)) K^2. The equilibrium problem itself is now SOLVED exactly (E* = 9/2 log 2, explicit density; `computations/equilibrium_s3.py`); what remains conjectural is only the standard link between the determinant and the energy (upper bound not written).
- True arithmetic cost of Delta_K (s=3, W=1) is (9/4+o(1))K^2; for s>=5 it is s(s+4)/(2(s+2)) K^2.
- lcm(1..K)^{3K} Delta_K in Z[X].
- omega_K is (1/4pi) x Wilson weight with parameters (1,1,1,-K) (identity checked by the orchestrator by hand); are there closed forms for the orthogonal polynomials or for Delta_K?

## Unresolved gaps
- CLOSED: elementary lcm bound. The other agent's period-6 argument (codex/problem131/denominator-design.md §1) gives log D_n <= 1.17053 n + O(log^2 n); we checked the proof by hand. Sufficient for zeta(3) (threshold 1.17516).
- Bailey's identity (16) (direct proof of Lemma 7) not verified.
- Beukers' paper (doi:10.1112/blms/11.3.268) not yet attempted.

## Independently proved observations, continued
- Equilibrium support for the s=3, W=1 field (orchestrator, by hand; standard theory of the singular integral equation assumed, positivity of the density and the inequality off the support NOT checked): symmetrising, the even measure sigma of mass 2*lambda on [-b,b] satisfies PV int dsigma(y)/(x-y) = 3pi/2 - arctan x for 0<x<b; the bounded solution has mass (2/pi) int_0^b y(3pi/2 - arctan y)/sqrt(b^2-y^2) dy = 3b + 1 - sqrt(1+b^2), using int_0^b sqrt(b^2-y^2)/(1+y^2) dy = (pi/2)(sqrt(1+b^2)-1). Hence 3b+1-sqrt(1+b^2) = 2 lambda, b = 3/4 at lambda = 1. This confirms the subagent's closed form for the support. The energy value (9/2) log 2 remains numerical (6 digits).
- Euler-Lagrange relations for the s=3 equilibrium measure, checked by the orchestrator (12:50 UTC): numerically (25 digits) mass = 2, l = 2+6 log 2, equality on the support, int V dsigma = 6 log 2 - 4, and the off-support derivative formula (2U-V)'(x) = 2[3 arcsin(b/x) - arctan(c/sqrt(x^2-b^2)) - pi] at x = 0.8, 1, 1.5, 3, 10; and by hand: that expression is 0 at x=b and its derivative has the sign of x^2(c-3b) - 3b < 0, so it is negative for x>b. Hence input (I3) of `hankel-energy-bound.md` is now checked by the orchestrator except that the derivative formula itself is verified numerically, not symbolically.
- zeta(5) simplification attempt: see `zeta5-simplification.md` (negative-leaning; numerical).

## Exploratory numerics (not proofs)
- Simplest Hankel family G=[mu_X(t^{i+j}/D_K)] (no D_N factor), K<=6: log|P(zeta(s))| for the primitive integer polynomial P at (K,h)=(6,5): -33.6 (s=3), -15.8 (s=5), +9.8 (s=7). Consistent with: arithmetic cost grows with s, real decay does not. Far from the asymptotic regime.

## Assumptions
- From 12:50 UTC the hourly resume task was slowed to every 3 hours to limit usage (the owner's 75% weekly cap cannot be measured from inside the session).
- The other agent's branch is `origin/codex/independent-research` (folder `codex/`), not `research-chatgpt` as INSTRUCTIONS.md says; we read that branch instead.
- The owner asked to stop at 75% of weekly usage; this session cannot read the usage meter, so it cannot enforce that. Mitigation: Opus for routine work, the xhigh Fable subagent only for frontier work.
- The owner's mid-run instruction: use a Fable subagent at extra-high reasoning ONLY for absolute frontier work; Opus subagents for coding/reading; orchestrator at medium otherwise.
- Trivial git commits/pushes and PROGRESS entries are done by the orchestrator directly (one shell command) because spawning a subagent for them costs more usage than it saves; Opus subagents are used for reading ChatGPT's branch once it has substantial content.
- An hourly scheduled task "Resume Mathathon overnight research" pokes this session so work resumes after a usage-limit pause (set up at the owner's request).
- Subagents cannot write report files outside the repo (harness refused); their reports come back as text.
