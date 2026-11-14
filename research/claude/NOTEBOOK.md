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
(none yet; a Hankel-determinant proof of zeta(3) irrationality is being attempted by the frontier subagent)

## Conjectures
(none yet)

## Unresolved gaps
- CLOSED: elementary lcm bound. The other agent's period-6 argument (codex/problem131/denominator-design.md §1) gives log D_n <= 1.17053 n + O(log^2 n); we checked the proof by hand. Sufficient for zeta(3) (threshold 1.17516).
- Bailey's identity (16) (direct proof of Lemma 7) not verified.
- Beukers' paper (doi:10.1112/blms/11.3.268) not yet attempted.

## Exploratory numerics (not proofs)
- Simplest Hankel family G=[mu_X(t^{i+j}/D_K)] (no D_N factor), K<=6: log|P(zeta(s))| for the primitive integer polynomial P at (K,h)=(6,5): -33.6 (s=3), -15.8 (s=5), +9.8 (s=7). Consistent with: arithmetic cost grows with s, real decay does not. Far from the asymptotic regime.

## Assumptions
- The other agent's branch is `origin/codex/independent-research` (folder `codex/`), not `research-chatgpt` as INSTRUCTIONS.md says; we read that branch instead.
- The owner asked to stop at 75% of weekly usage; this session cannot read the usage meter, so it cannot enforce that. Mitigation: Opus for routine work, the xhigh Fable subagent only for frontier work.
- The owner's mid-run instruction: use a Fable subagent at extra-high reasoning ONLY for absolute frontier work; Opus subagents for coding/reading; orchestrator at medium otherwise.
- Trivial git commits/pushes and PROGRESS entries are done by the orchestrator directly (one shell command) because spawning a subagent for them costs more usage than it saves; Opus subagents are used for reading ChatGPT's branch once it has substantial content.
- An hourly scheduled task "Resume Mathathon overnight research" pokes this session so work resumes after a usage-limit pause (set up at the owner's request).
- Subagents cannot write report files outside the repo (harness refused); their reports come back as text.
