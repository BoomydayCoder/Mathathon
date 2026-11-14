# The Fauzan ζ(5) proof: what is claimed, what is formalized, how it works

Sources (all retrieved in this session, stored outside the repo):

* A. Fauzan, *ζ(5) is irrational*, preprint, 17 Sep 2026, doi:10.5281/zenodo.22826419 (v1, 32 pp).
* `github.com/danromik/zeta5-irrationality` (HEAD `fc2fefa`) and `github.com/mo271/zeta5` (HEAD `7fe7367`).
* Fischler–Lai–Sprang–Zudilin, *Many more odd zeta values are irrational*, arXiv:2611.16156 (12 Nov 2026),
  found via the other agent's source list; it uses Fauzan's Hankel idea.

Provenance of statements below: **[me]** = read or checked by the orchestrator directly in the
files; **[sub]** = reported by a reading subagent from the files and not re-checked line by line;
**[ours]** = our own analysis.

## 1. What the preprint claims **[me: §1–§3.1 read; sub: the rest]**

* **Theorem 1.1.** $\zeta(5)\notin\mathbb Q$. **Corollary 1.2.** $|\zeta(5)-a/b|>b^{-260}$ for $b\ge b_0$ ($b_0$ not specified).
* **Theorem 2.1 (engine).** For $K=40n$, $N=3n$, $h=37n$ there are $Q_n\in\mathbb Z[X]$ with
  $\deg Q_n=37n$ and $0<Q_n(\zeta(5))<\exp(-139n^2/5)$ for large $n$.
* Deduction: if $\zeta(5)=a/b$, then $b^{37n}Q_n(a/b)$ is a positive integer less than
  $\exp(37n\log b-139n^2/5)\to0$. Note: **polynomials of growing degree, not linear forms**; the
  quadratic decay beats the linear cost $37n\log b$, and no height bound is needed.

### The construction

A $\mathbb Q$-linear functional $\mu_X$ on rational functions of $t$ with simple poles at $t=-j^2$:
$$\mu_X(t^e)=\frac{(-1)^eB_{2e+2}(2e+3)(2e+4)(2e+5)}{24},\qquad
\mu_X\Bigl(\frac1{t+j^2}\Bigr)=j^4\bigl(X-H_j^{(5)}\bigr)-\frac14+\frac1{2j}.$$
With $D_m(t)=\prod_{j\le m}(t+j^2)$:
$$G_K(X)=\Bigl[\mu_X\Bigl(\frac{D_N(t)^6\,t^{i+j}}{D_K(t)}\Bigr)\Bigr]_{0\le i,j<h},\qquad Q_{K,M}=m_{K,M}\,S_K\,\det G_K(X),$$
$S_K=\frac{(K!)^{2h}4^{h-1}}{(N!)^{12h}\prod_{i<h}((2i)!)^2}$, and $m_{K,M}$ a rational number defined prime by prime.

Three pillars:

1. **Positivity (Prop. 2.2).** $\mu_{\zeta(5)}(R)=\int_0^\infty R(y^2)w(y)\,dy$ with
   $w(y)=\frac{(2\pi)^4}{12}y^5\sum_{\ell\ge1}\ell^4e^{-2\pi\ell y}>0$. So $G_K(\zeta(5))$ is a Gram
   matrix and $\det>0$. This replaces "the linear form is non-zero".
2. **Real size (§6).** By Andréief, $\det G_K(\zeta(5))=\frac1{h!}\int\prod_{i<j}(t_i-t_j)^2\prod_i\frac{D_N^6}{D_K}(t_i)\,d\mu$,
   bounded by a weighted logarithmic-energy problem: $\log F_K(\zeta(5))\le UK^2+O(K\log K)$,
   $U=-2733991/2000000=-1.3670$. Computer-assisted (Appendix A: 684 cells; Appendix B: 143 pieces).
3. **Arithmetic (§3–5).** $p$-adic lower bounds for the valuation of *the whole determinant*
   (not of its entries), in three ranges of primes, summed with the prime number theorem:
   $\limsup K^{-2}\log m_{K,M}\le A_M$, $A_{200}\approx1.3496$.

**Margin:** $A_{200}+U\approx-0.0174$, about 1.3% of either term.

### Weak points noted in the danromik audit **[sub]** (audit is AI-written, not human-refereed)

* p. 10: "$b_a\le6\alpha x+3$ from (4.4)" is unproved and does not follow from (4.4); both Lean
  projects prove a correct replacement. **[me: confirmed that the README of danromik lists this as deviation 4.]**
* p. 12: a divided-difference step needs the $-\frac14+\frac1{2j}$ tail; not shown in the paper.
* Several numerical margins are extremely thin ($6\cdot10^{-8}$ in (6.4)).
* Appendices A, B are computer-assisted without published code; all thresholds are ineffective.

## 2. What is formalized **[me: statements, axiom and sorry greps checked directly]**

| | danromik/zeta5-irrationality | mo271/zeta5 |
|---|---|---|
| final statement | `theorem zeta5_irrational : Irrational zeta5` (`Zeta5/Interface.lean:437`), `zeta5 := ∑' v, 1/((v:ℝ)+1)^5` | `theorem irrational_five : ∃ x, Irrational x ∧ riemannZeta 5 = x` (`Apery/Main.lean:26`) |
| prime number theorem | **axiom** `chebyshev_theta_asymptotic : Tendsto (θ x / x) atTop (nhds 1)` (`Axioms.lean:46`), the only `axiom` in the repo | **theorem** `theta_asymptotic := chebyshev_asymptotic`, imported from the PrimeNumberTheoremAnd project (`Apery/PNT.lean:17`) |
| `sorry` | none in code (one occurrence inside a comment) | one, the intended stub `Challenge.lean:34` |
| route | follows the paper, with 13 listed deviations (README) | own simplifications: weaker rank bound, $M=400$, direct partial fractions, proves only decay $\approx5.6n^2$ |
| build | claims a clean build + `leanchecker` replay (Lean 4.34.0) **[sub]**; not rebuilt by us | no build log in repo **[sub]**; not rebuilt by us |

Caveats: (i) we did not build either project; (ii) for mo271 the PNT dependency contains `sorry`s in
files that the README says are not reached **[sub]**; we did not verify non-reachability; (iii) the
danromik statement is about $\sum_{v\ge1}v^{-5}$, the mo271 one about Mathlib's `riemannZeta 5`.

**Relation of the formal statements to the paper.** Both formalize Theorem 1.1 only (not
Corollary 1.2). Both *repair* the paper in places rather than transcribe it (danromik: 13 places,
of which nos. 4 and 5 are genuine gaps in the printed proof; mo271: different constants and a
weaker Lemma 4.2). So the correct summary is: *the theorem is formally verified modulo PNT
(danromik) or with PNT imported (mo271, subject to caveat ii); the paper as printed has at least
two gaps that the formalizations had to fill.*

**Independent human uptake.** FLSZ (arXiv:2611.16156) **[me: abstract and §1 read]** state that the
"principal idea was a new use of Hankel determinants", and build on it; their own proof uses an
elementary Pochhammer-basis denominator bound ("no sophisticated $p$-adic estimates") and a general
non-vanishing section. They do not claim to have refereed Fauzan's paper.

## 3. The construction in the series language of the ζ(3) proofs **[ours; (3.1) is Fauzan's]**

Fauzan's (3.1) says $\mu_X(R)=\tau_X\bigl(x^5R(-x^2)\bigr)$ where $\tau(P)=\mathsf L(P''')/24$,
$\mathsf L(x^k)=B_k$, and $\tau_X\bigl(\frac1{x-r}\bigr)=H^{(5)}_{d(r)}-X$. Since
$B_k=-k\zeta(1-k)$, $\mathsf L(P)$ is the zeta-regularised value of $-\sum_{x\ge1}P'(x)$. So, morally,

$$\mu_{\zeta(5)}(R)=-\frac1{4!}\sum_{x\ge1}\frac{d^4}{dx^4}\Bigl[x^5R(-x^2)\Bigr]\quad\text{(regularised)}.$$

Compare Zudilin/Nesterenko: $F_n=-\sum_{t\ge1}\frac{d}{dt}R_n(t)$ with $R_n$ having **double** poles.
The parallel and the difference:

| | Apéry–Nesterenko–Zudilin ($\zeta(3)$) | Fauzan ($\zeta(5)$) |
|---|---|---|
| summand | $\frac{d}{dt}$ of a function with double poles at integers | $\frac{d^4}{dx^4}$ of an **odd** function with **simple** poles at $\pm j$ |
| why only one zeta value | $\zeta(2)$ killed by $\sum\operatorname{Res}=0$ | simple poles differentiated 4 times give only order 5; nothing to kill |
| non-vanishing | Ball's positive series, or Casoratian | positive weight $w$ (the summand is a Laplace-type transform of a positive function) |
| object | one linear form per $n$ | $h\times h$ determinant of linear forms: a polynomial of degree $h$ in $\zeta(5)$ |
| smallness | $\lambda^n$, from zeros of $R_n$ | $e^{-cK^2}$, from the squared Vandermonde (energy), i.e. from orthogonality |
| arithmetic | $D_n^3$ entrywise | $p$-adic valuation of the determinant as a whole |

(Numerical verification of this dictionary, and of its analogue for every odd $s$, is in
`computations/fauzan_functional.py`; the all-odd-$s$ moment formulas were found independently by
the other agent, `codex/problem131/odd-moment-generalization.md` on branch
`origin/codex/independent-research`, and spot-checked by us numerically.)

## 4. Why an Apéry-style argument was not what worked **[ours; heuristic, flagged]**

* Linear forms in $1,\zeta(5)$ alone are easy to write down: any $\mu_{\zeta(5)}(R)$ is one. The
  classical obstruction is quantitative: denominators of order $D_K^5\approx e^{5K}$ against forms
  that cannot be made smaller than about $e^{-cK}$ with $c<5$ by the known constructions; and the
  well-poised series that *are* small enough bring $\zeta(3)$ along (Ball–Rivoal).
* The determinant changes the scale. A Gram determinant of $h\approx K$ functions is as small as
  $e^{-cK^2}$ because it is the product of the squared norms of the orthogonalised functions.
  Entrywise denominators would cost $e^{c'K^2}$ with $c'$ much larger than $c$; the point of §3–5 is
  that the determinant is $p$-adically far better than its entries (rank and vanishing arguments
  per prime). So the proof is won on the $K^2$ scale by a 1.3% margin.
* **Dependence on $s$ (heuristic, under test by a subagent).** On the $K^2$ scale the real side is
  independent of $s$ (the factor $y^s$ in the weight is invisible), while the arithmetic cost grows
  with $s$ through the denominators of $H_j^{(s)}$. This predicts: comfortable margin for $s=3$
  (so there should be a *Hankel-determinant proof of Apéry's theorem*, possibly with elementary
  arithmetic), razor-thin for $s=5$ (as observed), failure of the same family for $s=7$.
  Small exact experiments of the other agent ($K\le6$) are consistent with this ordering
  ($-\log|P(\zeta(s))|/\deg$: $7.30$, $3.17$, $-1.68$ for $s=3,5,7$ at $K=6$) but are far from the
  asymptotic regime.

## 5. Open items

* Read §4–§6 of the preprint line by line (only §1–§3.1 read directly so far).
* Outcome of the $s=3$ Hankel experiment and the $s$-dependence bookkeeping: see `NOTEBOOK.md`.
