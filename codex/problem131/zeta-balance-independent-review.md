# Independent review and a uniform growth theorem for binomial balance

Status: the identities and recurrence in `zeta2-zeta3-balance-unification.md`
are correct. No mathematical correction was needed. The degree classification
has exactly the stated limited scope. A stronger uniform statement about modes
and exponential growth is proved in Section 3 below. These are classical
sequences; neither this audit nor the explanatory organization establishes
historical novelty.

The original script passes. The independent
`codex/computations/check_balance_independent.py` checks the rational basis
actions with a symbolic, arbitrary basis index, both two-variable equations
including zero-index boundaries, the adjacent relations, the three growth
factors, and the exact degree formulas.

## 1. Basis transposition and the adjacent closure

The operator in the source is sought in the two-shift class
\(a(n)f(n+1)+b(n)f(n-1)\), with no central term. Its action at \(k=0\)
gives \(a+b=1\); the action at \(k=1\) gives
\[
a(n)(n+1)(n+2)+b(n)n(n-1)=3n(n+1)+2.
\]
Solving yields \(a=(n+1)^2/(2n+1)\) and \(b=-n^2/(2n+1)\), as claimed.
If an arbitrary central coefficient were allowed, two conditions would not
uniquely determine all three coefficients; the specified two-shift class is
what makes this discovery step precise.

For the all-index verification, divide by \(\phi_k(n)\) where allowed and use
\[
\frac{\phi_k(n+1)}{\phi_k(n)}=\frac{n+k+1}{n-k+1},\quad
\frac{\phi_k(n-1)}{\phi_k(n)}=\frac{n-k}{n+k},\quad
\frac{\phi_{k-1}(n)}{\phi_k(n)}=
\frac{k^2}{(n-k+1)(n+k)}.
\]
The required identity becomes the rational identity
\[
\frac{(n+1)^2(n+k+1)/(n-k+1)-n^2(n-k)/(n+k)}{2n+1}
=2k+1+\frac{2k^3}{(n-k+1)(n+k)}.
\]
Clearing denominators proves it identically; continuation of the polynomial
basis handles the removable exceptional values. The analogous calculation for
\(\psi_k\) is
\[
\frac{(m+1)^2k}{m+1-k}-mk
=k(k+1)+\frac{k^3}{m+1-k}.
\]
This verifies the second operator too.

In the kernel calculation, regard the sums as running over all \(k\ge0\),
with the integer binomial coefficients making them finite. This avoids
incorrectly retaining the cutoff \(\min(n,m)\) while shifting \(n\).
Index transposition then proves both two-variable equations. At a zero index,
the backward term is simply absent, so no undefined negative-index kernel is
needed.

For the adjacent closure, put \(x=k/N\) and use the common weights \(w_k\).
The four summand ratios are
\[
C_N:1+x,\qquad U:1-x^2,\qquad V:1-x,\qquad C_{N-1}:(1-x)^2.
\]
These are polynomial identities after cancellation, so include \(k=N\).
The balance identity says \(\sum w_k(1-x-x^2)=0\). Thus
\[
2U+V-C_N=2\sum w_k(1-x-x^2)=0,
\qquad 2V-U=C_{N-1}
\]
with the latter equality termwise. Solving gives exactly the two fractions
with denominator 5. Substitution in the first two-variable equation yields
\[
(n+1)^2C_{n+1}=(11n^2+11n+3)C_n+n^2C_{n-1}.
\]
The \(n=0\) case is consistent with \(C_0=1,C_1=3\). There is no implicit
identification of an off-diagonal kernel value with a diagonal one.

## 2. The low-degree classification is exact

For fixed \(N\ge1\), set
\(P_{r,s}(k)=(N-k)^r(N+k)^s-k^{r+s}\), with positive integers \(r,s\).
Its degree and leading coefficient are precisely

| Case | Degree | Leading coefficient |
|---|---:|---:|
| \(r\) odd | \(r+s\) | \(-2\) |
| \(r\) even, \(r\ne s\) | \(r+s-1\) | \((s-r)N\) |
| \(r=s\) even | \(2r-2\) | \(-rN^2\) |

The last row follows from \((N^2-k^2)^r-k^{2r}\). Since \(N>0\), the
listed coefficients do not vanish. Therefore degree at most two occurs exactly
for \((1,1),(2,1),(2,2)\). This classifies the degree of this particular single
balance polynomial. It neither classifies all low-order recurrences nor proves
that every higher-degree balance resists some different closure mechanism.

## 3. A uniform theorem: the same balance equation determines the growth rate

**Theorem.** Fix positive integers \(r,s\). Define
\[
w_{N,k}=\binom Nk^r\binom{N+k-1}k^s,\quad
S_N=\sum_{k=0}^Nw_{N,k},\quad
D_N=\sum_{k=0}^N\binom Nk^r\binom{N+k}k^s.
\]
There is a unique \(x_*\in(0,1)\) satisfying
\[
(1-x_*)^r(1+x_*)^s=x_*^{r+s}.                    \tag{A}
\]
The modes of \(w_{N,k}\) lie at \(Nx_*+O(1)\), and
\[
\log S_N=N\log\frac{(1+x_*)^s}{(1-x_*)^r}+O(\log N),\qquad
\log D_N=N\log\frac{(1+x_*)^s}{(1-x_*)^r}+O(\log N).    \tag{B}
\]
The constants in the error terms may depend on the fixed pair \((r,s)\).

**Proof of the mode statement.** The consecutive ratio is
\[
R_N(k)=\frac{(N-k)^r(N+k)^s}{(k+1)^{r+s}}.
\]
For \(0\le k<N\), its continuous logarithmic derivative is
\[
-\frac r{N-k}+\frac s{N+k}-\frac{r+s}{k+1}
\le-\frac r{N-k}-\frac r{k+1}<0,
\]
because \(N+k\ge k+1\). Thus strict monotonicity holds for **all** positive
\(r,s\); the condition \(s\le r\) used in the original sufficient argument
is unnecessary. A ratio equal to one can give two adjacent integer modes, but
cannot create a distant second peak.

After \(k=Nx\), the equation \(R_N(k)=1\) becomes
\[
(1-x)^r(1+x)^s=(x+1/N)^{r+s}.
\]
The limiting logarithmic ratio has negative derivative, tends to \(+\infty\)
at zero and \(-\infty\) at one, and hence has the unique root in (A).
Its derivative is nonzero there. Uniform convergence near that root, or the
mean value theorem applied to the logarithmic equation, gives
\(x_N=x_*+O(1/N)\). Rounding the crossing to an integer proves the mode
statement. The finitely many small \(N\), including the tie at \(N=1\), do
not affect the asymptotic claim.

**Proof of the exponential rate.** Uniformly for \(0\le k\le N\), the
factorial estimate \(\log j!=j\log j-j+O(\log(j+1))\), with
\(0\log0=0\), gives
\[
\log w_{N,k}=NF(k/N)+O(\log N),
\]
where
\[
F(x)=-r(1-x)\log(1-x)+s(1+x)\log(1+x)
       -(r+s)x\log x.
\]
Its derivative is the logarithm of the ratio in (A), and
\[
F''(x)=-\frac r{1-x}+\frac s{1+x}-\frac{r+s}{x}
=-\frac r{x(1-x)}-\frac s{x(1+x)}<0.
\]
Thus its unique maximum is at \(x_*\). The largest summand bounds \(S_N\)
from below, and \((N+1)\) times that summand bounds it from above. Choosing
an integer nearest \(Nx_*\) in the lower bound proves
\(\log S_N=NF(x_*)+O(\log N)\).

Using (A) in \(F(x_*)\) cancels the terms multiplied by \(x_*\), giving
\[
F(x_*)=-r\log(1-x_*)+s\log(1+x_*).
\]
Finally, the diagonal summand differs from the edge summand by
\(((N+k)/N)^s\), between 1 and \(2^s\). Hence
\(S_N\le D_N\le2^sS_N\), proving the same rate for \(D_N\).
This proves the theorem without a Gaussian approximation or a recurrence
asymptotics theorem.

## 4. The three classical scales follow together

Applying the theorem to the low-degree cases gives:

| \((r,s)\) | Positive balance root \(x_*\) | Exponential factor | Diagonal sequence |
|---|---|---|---|
| \((1,1)\) | \(1/\sqrt2\) | \((1+\sqrt2)^2\) | \(P_N(3)\) |
| \((2,1)\) | \((\sqrt5-1)/2\) | \(((1+\sqrt5)/2)^5\) | \(C_N\) |
| \((2,2)\) | \(1/\sqrt2\) | \((1+\sqrt2)^4\) | \(A_N\) |

For the middle case, \(1-x_*=x_*^2\) and \(1+x_*=1/x_*\), so the factor is
\(x_*^{-5}\). For the other two cases, the ratio
\((1+x_*)/(1-x_*)=(1+\sqrt2)^2\). Thus the same exact balance relation has
two roles: its low degree permits the adjacent-value closure, and its scaled
root locates the dominant summands and determines the exponential growth.

This theorem concerns the integer denominator sequences. It does not identify
companion limits with logarithms or zeta values, prove their denominator bounds,
or construct the necessary small nonzero linear forms. Those remain separate
steps in the respective irrationality arguments.
