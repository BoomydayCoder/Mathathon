# Referee report on a claimed proof that ζ(7) is irrational

Preprint: P. Anand, "ζ(7) is irrational", Zenodo 22920911 (Sept 2026), retrieved 14 Nov 2026.
Referee work: frontier subagent of this session (scripts `computations/anand_zeta7_*.py`, outputs
`output_anand_zeta7_*.txt`). The orchestrator checked the logic of §2 (a) and (b) below by hand but did
**not** recompute the constants or re-read the preprint. AI-refereed only. We make no statement about
whether ζ(7) is irrational.

## Verdict

**The main estimate of the preprint cannot be true, and its error is located.** The preprint's local
($p$-adic) estimates were not found to be wrong; the failure is a numerical constant in the summation
over primes.

## 1. What the preprint claims

Same architecture as Fauzan's ζ(5) proof: $K=40n$, $N=3n$, $h=37n$;
$\mu_X(1/(t+j^2))=j^6(X-H_j^{(7)})+\frac1{2j}-\frac16$; $G_K=[\mu_X(D_N^8t^{i+j}/D_K)]_{i,j<h}$; $F_K=S_K\det G_K$;
$Q_{K,M}=m_{K,M}F_K$. Claims: $Q_{K,M}\in\mathbb Z[X]$ and $0<Q_{40n,M}(\zeta(7))<\exp(-19000n^2)=\exp(-11.875K^2)$,
from an arithmetic constant $A_{200}=-13.7529$ (so the normaliser $m\approx e^{-13.75K^2}$ is *tiny*) and
a real-side constant $U=+1.76$.

## 2. Why this is impossible

**(a) Internal contradiction (elementary).** If $Q_{K,M}=m\,F_K\in\mathbb Z[X]$, its leading coefficient
$m\cdot\operatorname{lead}(F_K)$ is a non-zero integer, so $m\ge1/|\operatorname{lead}F_K|$. By the preprint's own formulas
$\log|\operatorname{lead}F_K|\le3.65K^2+O(K\log K)$, hence $\log m\ge-3.65K^2$. This contradicts the claimed
$\log m\le(-13.75+o(1))K^2$. So the integrality statement and the size of the normaliser cannot both hold.

**(b) The bound is false for every $n\ge1$ (elementary lower bound).** Any integer polynomial
proportional to $\det G_K$ has absolute value at ζ(7) at least $\det G_K(\zeta(7))/|\operatorname{lead}|$; bounding the
Gram determinant from below by Legendre-type minimal norms on an interval gives
$\log|Q(\zeta(7))|\ge-7.171K^2-K(2+\log2K)>-11.875K^2$.

**(c) Exact computation, $n=1,\dots,4$.** With the preprint's exact matrix (the positive integral
representation was confirmed to 36 digits):

| $K$ | $\log F_K(\zeta7)/K^2$ | true normaliser $\log m/K^2$ | $\log P_K(\zeta7)$ | claimed |
|---|---|---|---|---|
| 40 | −0.655 | +1.452 | +1274 | −19000 |
| 80 | −0.707 | +1.579 | +5576 | −76000 |
| 120 | −0.725 | +1.609 | +12728 | −171000 |
| 160 | −0.734 | +1.620 | +22661 | −304000 |

$P_K$ is the primitive integer polynomial; every integer polynomial proportional to the determinant
is an integer multiple of it, so its value at ζ(7) is at least $e^{1274},\dots,e^{22661}$: far above 1.

## 3. Where the error is

* The outer integral $I_{\rm out}$ (eq. (122), Table 5) is printed as $-15.2819$. Recomputed exactly from
  the preprint's own definitions (58), (59), (62) it is $\frac{453803}{288000}=+1.5757$. Table 5 agrees with
  (59) only on one sub-interval; e.g. on $[1,\frac{37}{20}]$ the preprint itself sets the local bound to $0$,
  so the integrand should be $\frac{37}{20}-y\ge0$, but the table prints $-\frac{239}{20}-2y$ (that row alone
  contributes $-12.58$ instead of $+0.36$).
* With the corrected value the arithmetic constant is about $+1.72$ at $M=200$ instead of $-13.75$, and
  the total exponent is positive: $+3.48$ with the preprint's $U$, about $+0.96$ with the true real constant
  ($-0.76$ from our equilibrium solver). No contradiction with integrality arises.
* The normaliser evaluated directly from the preprint's definitions at an admissible pair
  $(K,M)=(320000,40)$ gives $\log m/K^2=+1.907$, against the claimed $\le-13.6$.
* Other slips: a constant "$A_0$" is digit-for-digit Fauzan's ζ(5) constant; (82) has $25$ where $16$ is
  correct; the theorem states $-19000n^2$ where its proof derives $-17000n^2$.
* Checked and found **correct**: the functional and its positive integral representation; the outer
  local bound (59) at all 66 applicable primes for $K\le160$ (equality at 50 of them); the crude small-prime
  bound; the table for $\int_3^{20}$.

## 4. Is it repairable?

No, as far as the evidence goes: the local bounds appear asymptotically sharp, and then the true limit
of $\log P_K(\zeta(7))/K^2$ is about $+0.92$ (exact values $0.797,0.871,0.884,0.885$ for $n=1..4$), in
agreement with our independent computations (`zeta7-status.md`, `zeta5-margin.md` §5). The family fails
for ζ(7).

## 5. Scope

Rigorous for all $n$: (a) and (b). Exact for $n\le4$: (c). Asymptotic "$+0.92$": numerical extrapolation.
Not verified: the preprint's §3, its inner-range proposition (not testable at feasible $K$), its
interval-arithmetic certificate, and its irrationality-measure corollary.
