# Independent audit of the Catalan endpoint arithmetic for log 2

**Verdict: accepted as a complete irrationality proof**, conditional only
on the explicitly cited and already independently reviewed positive
deformation/isometry theorem. I found no missing integer factor, boundary
case, or growth assumption in
[log2-endpoint-catalan-arithmetic.md](log2-endpoint-catalan-arithmetic.md).
The logarithmic irrationality result and Wilson family are classical;
this review makes no historical novelty claim.

The independent checker is
[check_log2_endpoint_independent.py](../computations/check_log2_endpoint_independent.py).
It reconstructs the endpoint moments from the Newton-basis integrals,
rather than reusing the author's Bernoulli-moment implementation.

## 1. Endpoint normalization and the denominator evaluation

The cited deformation review proves that
\[
 d\nu(y^2)=16y/\sinh(2\pi y)\,dy
 \]
is a probability measure, including the endpoint with no additional atom.
Its Hilbert Gram identity has
\(\delta_k=(1/2)_k(3/2)_k\), and the factorial expression
\[
 \delta_k=\frac{(2k)!(2k+1)!}{16^k(k!)^2}
 \]
is correct. Thus the transported shifted-Legendre polynomial has norm
\(1/(2n+1)\).

At \(s_n=-(n+1/2)^2\),
\[
 D_k(s_n)=(-1)^k\frac{(n+k)!}{(n-k)!}.
 \]
Multiplying by its coefficient in \(\Phi_n\) cancels the sign and gives
precisely \(16^k\binom{n+k}{2k}^2/(2k+1)\).
This proves the displayed formula for \(A_n\), including \(k=0\).

The binomial identity
\[
 \frac{2n+1}{2k+1}\binom{n+k}{2k}
 =\binom{n+k}{2k}+2\binom{n+k}{2k+1}
 \]
is exact, with the last binomial zero at \(k=n\). Multiplying it
by the second binomial factor proves \(Jw_{n,k}\in\mathbb Z\),
and hence \(JA_n\in\mathbb Z\), without cancellation between terms.

## 2. Newton quotient, sign and Catalan–beta factor

The divided-product identity is correct. Its tail at the special point is
\[
 \prod_{r=m}^{k-1}(r-n)(r+n+1)
 =(-1)^{k-m}\frac{(n-m)!(n+k)!}{(n-k)!(n+m)!}.
 \]
The outer minus sign from \(A_n-\Phi_n\) combines with the coefficient
\((-1)^k\) to produce \((-1)^{m-1}\). Substituting
\(\nu(D_{m-1})=\delta_{m-1}/m\) gives exactly equation (6).

Using the factorial expression for \(\delta_{m-1}\) extracts a factor
\(16^{-(m-1)}\). After division of the entire numerator by \(16\),
the remaining power is \(16^{k-m}\), an integer because \(m\le k\).
There is no remaining power-of-two denominator hidden in equation (7).

The beta factor is
\[
 \int_0^1x^{n-m}(1-x)^{2m-1}\,dx
   =\frac{(n-m)!(2m-1)!}{(n+m)!}.
 \]
Multiplication by
\(\operatorname{Cat}_{m-1}=(2m-2)!/[m((m-1)!)^2]\)
is exactly \(H_{n,m}\). Both exponents in the integral are
nonnegative, including the case \(m=n\).
Catalan integrality follows from the displayed binomial difference;
at \(m=1\) it is simply \(\operatorname{Cat}_0=1\), with the
negative lower-index binomial interpreted as zero.

The integrand is an integer polynomial of degree \(n+m-1\).
Termwise integration introduces only denominators at most
\(n+m\le2n<J\). Hence \(\ell_J H_{n,m}\) is an integer.
Combining this with \(Jw_{n,k}\in\mathbb Z\) and \(16^{k-m}\in\mathbb Z\)
proves
\[
 J\ell_J\,\mathcal N_n/16\in\mathbb Z
 \]
term by term. The proof does not assume cancellation among the signed
summands, nor extrapolate a finite divisibility pattern.

## 3. Harmonic conversion and the excluded index

The partial fractions and finite harmonic identity are correct:
\[
 \frac{T_n+4}{8}=H_{2n+1}-H_n-\frac1{2(2n+1)}.
 \]
For \(n\ge1\), \(J\ge3\) is odd and \(\ell_J\) is divisible by
both \(2\) and \(J\); coprimality then gives \(2J\mid\ell_J\).
Therefore the harmonic factor is cleared by this same single lcm.
Together with the numerator result this proves
\(J\ell_JB_n\in\mathbb Z\).

At \(n=0\), the definitions give \(A_0=1,\mathcal N_0=0,B_0=1/2\).
Since \(J=\ell_J=1\), the second normalizer really fails.
The restriction \(n\ge1\) is necessary, not merely convenient.

The first case is exactly
\[
 A_1=19/3,\quad \mathcal N_1=8/3,\quad
 B_1=79/18,\quad B_1/A_1=79/114.
 \]
It agrees with the direct polynomial quotient and with the harmonic
conversion independently.

The endpoint target \(8\log2-4\) also follows directly from the finite
identity: \(H_{2n+1}-H_n\to\log2\) by elementary integral comparison,
and \(1/(2J)\to0\). This is consistent with the already reviewed
deformed resolvent and requires no special-function evaluation.

## 4. Positivity and the exact final scale

The positive Padé identity gives
\[
 A_nZ_{1/2}-b_n
 =\frac1{2A_n}\int\frac{\Phi_n(t)^2}{t+(n+1/2)^2}\,d\nu(t)>0.
 \]
The norm bound gives its upper bound \(2/(J^3A_n)\).
Since
\[
 A_nZ_{1/2}-b_n=8(A_n\log2-B_n),
 \]
the scaled logarithmic error is bounded by
\(\ell_J/(4J^2A_n)\), exactly as written. The factors two from
the resolvent and eight from the target have both been accounted for.

The comparison with the ordinary Apéry denominator is termwise:
\[
 \frac{(k!)^2}{\delta_k}
 =\prod_{r=1}^k(1-1/(4r^2))^{-1}\ge1.
 \]
Thus the already proved positive Delannoy lower bound applies without
any new recurrence or asymptotic assertion.

## 5. Elementary lcm induction and irrationality

Let \(r=\lceil M/2\rceil\). A prime power \(p^a\in(r,M]\)
has floor contribution
\[
 \lfloor M/p^a\rfloor-\lfloor r/p^a\rfloor
       -\lfloor\lfloor M/2\rfloor/p^a\rfloor=1.
 \]
All other floor contributions to a binomial valuation are nonnegative.
There is at most one such power for each prime, because the ratio of
successive powers is at least two. It follows that
\(\ell_M/\ell_r\) divides the central binomial coefficient, proving
the claimed lcm divisibility.

For \(M\ge2\), Pascal's identity bounds each binomial coefficient
by the sum \(2^{M-1}\) of the previous row. Applying the induction
hypothesis at \(r<M\) gives
\[
 \ell_M\le4^{r-1}2^{M-1}\le4^{M-1}.
 \]
The last inequality is equivalent to \(2r\le M+1\), valid for
both parities. The base case \(M=1\) is exact.
No prime distribution estimate is used.

Since \(J=2n+1\), the lcm bound is \(\ell_J\le16^n\).
Combining it with the denominator lower bound gives precisely
\[
 0<J\ell_J(A_n\log2-B_n)
 \le\frac{n+1}{4}\left(\frac{16}{17+12\sqrt2}\right)^n\longrightarrow0.
 \]
If \(\log2=p/q\), \(q>0\), multiplying by \(q\) gives a positive
integer: \(J\ell_JA_n\) and \(J\ell_JB_n\) are integers.
Its convergence to zero is impossible. The irrationality assembly
is therefore complete.

## 6. Reproducible checks and scope

The independent Fraction-only checker passes:

- 169 endpoint Gram identities reconstructed from the Newton integrals;
- 13 direct polynomial quotients, norms, and rational log-error
  certificates, including \(n=0\);
- 37,820 termwise Catalan–beta integrality checks through \(n=60\);
- both normalizers, the harmonic identity, the first approximation and
  the exact \(n=0\) failure;
- every lcm-induction ingredient through \(M=1000\).

The logarithm certificates use the rational bounds
\[
 2\sum_{j=0}^{K-1}\frac1{(2j+1)3^{2j+1}}
 <\log2<
 2\sum_{j=0}^{K-1}\frac1{(2j+1)3^{2j+1}}
 +\frac9{4(2K+1)3^{2K+1}},
 \]
not numerical logarithm evaluation. The author's checker also passed.

The Catalan connection is an exact factorial identity and an
integer-polynomial integration mechanism. It does not prove a general
arithmetic theorem for intermediate deformation parameters, nor does
irrationality at this endpoint imply irrationality of any parameter
limit.
