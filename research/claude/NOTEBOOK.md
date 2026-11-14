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

## Proof drafts
(none yet)

## Conjectures
(none yet)

## Unresolved gaps
- Self-contained elementary proof of limsup D_n^{1/n} < (sqrt2+1)^{4/3}=3.2388 (Chebyshev/Hanson) not written out; attributions quoted from memory.
- Bailey's identity (16) (direct proof of Lemma 7) not verified.
- Beukers' paper (doi:10.1112/blms/11.3.268) not yet attempted.

## Assumptions
- The owner's mid-run instruction: use a Fable subagent at extra-high reasoning ONLY for absolute frontier work; Opus subagents for coding/reading; orchestrator at medium otherwise.
- Trivial git commits/pushes and PROGRESS entries are done by the orchestrator directly (one shell command) because spawning a subagent for them costs more usage than it saves; Opus subagents are used for reading ChatGPT's branch once it has substantial content.
- An hourly scheduled task "Resume Mathathon overnight research" pokes this session so work resumes after a usage-limit pause (set up at the owner's request).
- Subagents cannot write report files outside the repo (harness refused); their reports come back as text.
