# The Hankel method applied to Catalan's constant and other low-weight constants

**Status: a valid, very simple functional exists, but the method does not prove $G\notin\mathbb Q$.**
Source: frontier subagent (scripts `computations/catalan_*.py`, outputs `output_catalan_*.txt`). The
orchestrator derived the first functional independently in outline; nothing else was re-derived.
Exact finite computations; not a proof of impossibility. Context: irrationality of $G$ was claimed by an
OpenAI note of 24 Sep 2026 with a much more elaborate determinant (see `literature-status.md`); we did
not verify that note.

## 1. Positive functionals with rational moments and pole values in $\mathbb Q+\mathbb QG$

All have $\mu_{X_0}(R)=\int_0^\infty R(y^2)w(y)\,dy$ with $w>0$.

| | weight | moments $\mu(t^e)$ | pole values |
|---|---|---|---|
| (I) half-integer nodes $a_j=j+\frac12$, $X=4G$ | $\pi y^2\cosh(\pi y)/\sinh^2(\pi y)$ | $(2^{2e+2}-1)|B_{2e+2}|$ | $(-1)^ja_j\bigl(X-\sum_{v<j}\frac{(-1)^v}{(v+1/2)^2}\bigr)-\frac1{2a_j}$ |
| (II) integer nodes, $X=G$ | $\frac{\pi y}4\sinh(\pi y)/\cosh^2(\pi y)$ | $(2e+1)|E_{2e}|/2^{2e+3}$ | $(-1)^n\bigl(X-\sum_{k<2n}\chi_4(k)/k^2\bigr)$ |
| $L(2,\chi_{-3})$ | $-\frac2{3\sqrt3}y\frac{d}{dy}\frac{\sinh\pi y}{\sinh3\pi y}$ | rational | $X-\sum_{k<3n}\chi_{-3}(k)/k^2$ |

No $\zeta(2)$ or $\pi$ appears: the alternation removes it (with the kernel $1/(e^{2\pi y}-1)$ at nodes
$n-\frac34$ one gets $\pi^2+8G$ instead). 156 numerical checks to $\ge32$ digits
(`catalan_functional_check.py`). The Hankel determinant $\Delta_K=\det[\mu_X(t^{i+j}/D_K)]$ has degree $K$ and
is positive at the constant.

## 2. Numbers (per $K^2$; total $=\log|P_K(X_0)|$ for the primitive integer polynomial; a proof needs total $<0$)

| constant [functional] | $K$ | cost | $\log\Delta$ | total |
|---|---|---|---|---|
| **$G$** [II] | 40 / 80 / 120 / 160 | 2.58 / 2.66 / 2.69 / 2.72 | −2.198 | **+0.38 / +0.46 / +0.49 / +0.52** |
| **$G$** [I] | 40 / 80 / 120 | 2.46 / 2.59 / 2.64 | −2.07 … −2.15 | **+0.39 / +0.47 / +0.50** |
| $L(2,\chi_{-3})$ | 20 … 100 | 2.95 … 3.15 | −3.12 | −0.18 … +0.03 (changes sign near $K=80$) |
| $\zeta(3)$ [standard] | 100 | 2.13 | −3.04 | −0.91 |
| $\zeta(2)$ [standard] | 120 | 2.25 | −3.08 | −0.83 |
| $\zeta(2)$ [sinh kernel] | 100 | 1.43 | −2.15 | −0.73 |
| $\log2$ | 100 | 1.49 | −2.20 | −0.71 |
| $\pi$ | 100 | 1.80 | −2.18 | −0.38 |
| $\beta(4)$ | 60 | 4.44 | −2.01 | +2.43 |

So the simple Hankel family appears to reprove irrationality of $\zeta(2)$, $\log2$, $\pi$ (numerical
only; no proofs written), and fails for $G$, $\beta(4)$. $L(2,\chi_{-3})$ is borderline.

## 3. Why $G$ fails

* **Decay.** For kernels $\sim e^{-\gamma y}$ with unit node density the equilibrium problem has the same
  closed form as in `equilibrium_s3.py`; for $\gamma=\pi$ the energy is $2\log3=2.197$ (numerical, 30
  digits), confirmed by exact determinants ($-2.1976$ at $K=160$). Weaker than the $\frac92\log2=3.12$ of
  the ζ(3) kernel, because the kernel decays half as fast.
* **Cost.** The partial sums $\sum\chi_4(k)/k^2$ run up to $k\approx2K$, so primes up to $2K$ enter.
  Proved for functional (II): $v_p(\Delta_K)\ge-3(K-\frac{p-1}2)$ for $K<p<2K$; equality observed at every
  such prime for $K\le160$. With the observed pattern for $p\le K$ the cost is $\frac83K^2$ under PNT,
  against decay $2.197$. **"Small weight means wide margin" is false here.**
* **No reflection gain**: the lemma $m_a\equiv m_{p-a}$ of the ζ(3) proof fails (opposite signs at $n$, $p-n$).
* **The prime 2 costs order $K^2$** (about $0.12K^2$), as in the flaw found by Wachs in an earlier claim.
* Variants (same-parity nodes, Fauzan-type numerators, rectangular $h=\lambda K$, shifted windows; a
  94-configuration scan): best total $+0.24$ to $+0.33$; none negative.

## 4. Relation to the OpenAI proof (as read by the subagent; unverified)

Theirs: signed determinants of size $48N$ with entries in $\mathbb Q+\mathbb QG+\mathbb Q\zeta(2)$, the $\zeta(2)$ part
cancelled by high Taylor contact, non-vanishing by a separate argument, margin about $10^{-4}n^2$. Ours:
one positive weight, no $\zeta(2)$, non-vanishing for free, but short by about $0.5K^2$ out of $2.2K^2$.
The comparison suggests that using both parity classes of nodes (hence $\zeta(2)$ and its cancellation)
is necessary.

## 5. Next steps

1. Mixed nodes with one positive weight: half-integer poles give $G$, integer poles give $\pi^2$; impose
   vanishing of the $\pi^2$ part by divided-difference conditions. Untested.
2. $L(2,\chi_{-3})$ needs only about $0.1K^2$: compute the exact asymptotic cost and try numerators.
3. Write the easy Hankel proofs for $\zeta(2)$, $\log2$, $\pi$ along `hankel-zeta3.md`.

## 6. Mixed nodes and π²-cancellation with positive weights: negative, with a structural reason

Second frontier run (`computations/catalan2_*.py`, outputs `output_catalan2_*.txt`). Not re-derived by the
orchestrator; determinants exact, lemmas by the subagent.

* **Cancellation lemma.** If poles at some nodes carry $\pi^2$, the $\pi^2$-part of $\mu(t^e\Psi)$ vanishes for
  $e=0,\dots,E$ iff the residues there are $q(-a^2)/D_S'(-a^2)$ with $\deg q\le|S|-E-2$. A Hankel matrix of size
  $h$ needs $E=2h-2$, hence **at least $2h$ extra $\pi^2$-carrying nodes**, each bringing the usual
  denominators.
* **Obstruction lemma.** A non-negative $\pi^2$-free test function $\Psi$ must satisfy $t\Psi(t)\to c>0$: it
  decays exactly like $1/t$ and cannot have the $1/D_K$ shape that produces the decay. So cancellation is
  paid in nodes (cost) and buys nothing analytically.
* **Numbers** (per $h^2$; baseline single-class Hankel: $+0.45$ at $h=64$): two quarter-classes with
  alternating residues, positive: $+5.5\to+6.1$; half-integer plus $2h$ integer nodes: $+4.0\to+4.5$; pole
  rows plus polynomial rows (best $\pi^2$-cancelling family): $+2.3\to+2.6$, increasing. Polynomials in two
  unknowns $(G,\pi^2)$ are dead by logic.
* **What the OpenAI construction has that this class lacks** (as read by the subagent; unverified): a
  Cauchy kernel $1/(1-ts)$ in which cancellation of $\zeta(2)$ is per row and costs no nodes ("Taylor
  contact"), and a Frobenius self-similarity of the moment array that keeps the cost near $1.94n^2$ with a
  small 2-adic part; its determinants are signed and non-vanishing is proved separately.

**Conclusion for this line:** positive-weight Hankel determinants with simple poles, with or without
mixed nodes, do not reach $G$. Item 5.1 above is closed negatively.
