# Independent review: integer depletion and discrete-uniform cumulants

**Status: independently proved and checked.** The uniform polynomial denominator theorem in `integer-depletion-cumulants.md` is correct for every odd integer \(s\ge3\) and integer \(\ell\ge2\). This review checks the positive density, generating series, prime support, constructive sharpness, real rescaling and Hurwitz constants. No new irrationality theorem or historical novelty claim follows.

## 1. Positive density and its moments

Write \(f(y)=(e^{2\pi y}-1)^{-1}\). Since \(s-1\) is even,
\[
 w_s(y)=\frac{2y^s}{(s-1)!}f^{(s-1)}(y)
 =\frac{2(2\pi)^{s-1}}{(s-1)!}y^s
   \sum_{n\ge1}n^{s-1}e^{-2\pi ny}>0.
\]
In \(\ell w_s(y/\ell)\), the terms indexed by multiples of \(\ell\) are exactly those of \(w_s(y)\). Removing them proves the stated positive series for the depleted density. The classical Riemann-sum limit
\(x^s\sum n^{s-1}e^{-nx}\to\Gamma(s)\)
gives \(w_s(0+)=1/\pi\), hence depleted limit \((\ell-1)/\pi\). The remaining term \(n=1\) gives exponential rate \(2\pi/\ell\).

Changing variables proves
\[
 \nu_{s,\ell}(R)=\ell^2\mu_s(R(\ell^2t))-\mu_s(R(t)).
\]
Positive-term integration gives the undepleted moment directly:
\[
 \mu_s(t^e)=\frac{2(2e+s)!}{(s-1)!(2\pi)^{2e+2}}\zeta(2e+2)
 =(-1)^e B_k\frac{(k+s-2)!}{(s-1)!k!},\quad k=2e+2.
\]
Thus both normalizations and the sign in the stated depleted moment formula are correct.

## 2. Cumulants remove all primes outside ell

For a discrete uniform variable on \(\{0,\ldots,\ell-1\}\),
\[
 M(z)=\frac{e^{\ell z}-1}{\ell(e^z-1)},\qquad
 \frac{M'(z)}{M(z)}
 =\frac{\ell-1}{2}+
   \sum_{k\ge2}\frac{(\ell^k-1)B_k}{k!}z^{k-1}.
\]
The odd terms beyond the linear cumulant vanish. Consequently, for even \(k\ge2\),
\(\kappa_k=(\ell^k-1)B_k/k\).

There is no hidden logarithm-denominator problem. If \(m_j=M^{(j)}(0)\), then
\[
 \kappa_n=m_n-\sum_{j=1}^{n-1}\binom{n-1}{j-1}\kappa_jm_{n-j}.
\]
Since every \(m_j\in\mathbb Z[1/\ell]\), this recurrence proves the same for every cumulant using integer coefficients only. The identity
\[
 \frac{k}{s-1}\binom{k+s-2}{s-2}=\binom{k+s-2}{s-1}
\]
then proves
\[
 \nu_{s,\ell}(t^e)=(-1)^e\kappa_k\binom{k+s-2}{s-1}
 \in\mathbb Z[1/\ell].
\]
In particular, a prime dividing \(s-1\) but not \(\ell\) really does cancel.

## 3. Exact least uniform denominator

For \(p\mid\ell\), \(\ell^k-1\) is a unit and von Staudt–Clausen gives \(v_p(B_k)\ge-1\). The binomial is integral. This proves the upper bound
\[
 C_{s,\ell}=\prod_{p\mid\ell}p^{1+v_p(s-1)}.
\]
To prove it is least, put \(a=v_p(s-1)\). Lucas's theorem is being applied in its precise no-carry form to \(\binom{k+s-2}{s-2}\).

For \(p=2\), oddness of \(s\) ensures \(a\ge1\). The \(a\)-th digit of \(s-2\) is zero, so \(k=2^a\) causes no carry and is even. For odd \(p\), write \(s-1=p^au\). The digit at position \(a\) of \(s-2\) is \(u-1\pmod p\), between 0 and \(p-2\). This remains true when \(a=0\). Choose \(b>a\) with \(p^b>s-2\) and take
\[
 k=p^a+(p-2)p^b.
\]
There are no carries; the digit sum is \(p-1\), so \(p-1\mid k\) and \(k\) is even. Von Staudt–Clausen gives exactly \(v_p(B_k)=-1\). Therefore the moment has valuation exactly \(-1-a\).

Different primes may attain their worst denominators at different moments. This is sufficient for a *uniform* denominator: every required prime power is individually forced. The proof does not require a single moment attaining all of them.

## 4. Rescaling: an actual denominator-growth corollary

The normalization must be stated. With
\[
 \sigma_{s,\ell}(x)=\ell^{-1}w_{s,\ell}(\ell x)
 =w_s(x)-\ell^{-1}w_s(\ell x)>0,
\]
its \(e\)-th moment equals \(\ell^{-k}\nu_{s,\ell}(t^e)\), \(k=2e+2\). This is not the unnormalized pushforward measure; the explicit prefactor accounts for the extra \(\ell^{-2}\). The parent's revised statement specifies it correctly.

The resulting denominator growth is exact along explicit subsequences. For odd \(p\mid\ell\), allow \(b\to\infty\) in the preceding construction. For \(p=2\), instead take \(k=2^a+2^b\), with \(b>a\) and \(2^b>s-2\). These remain even no-carry witnesses. Hence infinitely many normalized rescaled moments satisfy
\[
 \boxed{v_p\!\left(\int_0^\infty x^{2e}\sigma_{s,\ell}(x)\,dx\right)
 =-k\,v_p(\ell)-1-v_p(s-1).}
\]
Thus the restored decay rate comes with genuinely growing prime-power denominators in these polynomial moments, not merely a growing sufficient normalizer. This is a statement about this particular rescaling and moment sequence, not a general determinant or irrationality obstruction.

## 5. Pole identity and endpoint audit

For every real \(a>0\), the underlying formula is
\[
 \mu_s\!\left(\frac1{t+a^2}\right)
 =a^{s-1}\zeta(s,a)-\frac1{s-1}-\frac1{2a}.
\]
The sign of the last term matters. Here is a direct reduction to Hermite's formula. Put \(g(y)=y^s/(y^2+a^2)\). Integration by parts \(s-1\) times gives
\[
 \mu_s\!\left(\frac1{t+a^2}\right)
 =\frac2{(s-1)!}\int_0^\infty g^{(s-1)}(y)f(y)\,dy.
\]
At zero, every boundary product is \(O(y)\):
\(g^{(j)}=O(y^{s-j})\) and
\(f^{(s-2-j)}=O(y^{-s+1+j})\).
At infinity exponential decay dominates all polynomial factors. Polynomial division of \(g\), followed by differentiation, proves
\[
 \frac{g^{(s-1)}(y)}{(s-1)!}
 =a^{s-1}\operatorname{Im}(a-iy)^{-s}.
\]
Hermite's integral formula now yields the displayed identity, including both constants.

The dilation identity cancels the constant \(-1/(s-1)\), leaving
\[
 \nu_{s,\ell}\!\left(\frac1{t+j^2}\right)
 =(j/\ell)^{s-1}\zeta(s,j/\ell)-j^{s-1}\zeta(s,j)
 -\frac{\ell-1}{2j}.
\]
This holds for real \(j>0\). Statements about half-integer or integer Hurwitz arguments additionally assume positive integer \(j\). The parent was alerted to make that scope explicit. General rational Hurwitz arguments need not simplify to this single zeta value; no independence claim is needed here.

## 6. Independent reproducible checks

`codex/computations/check_integer_depletion_independent_review.py` uses only the Python standard library. It constructs cumulants from exact uniform raw moments, independently of Bernoulli numbers, and checks:

- 17,980 polynomial moments against both formulas and the claimed uniform denominator;
- 1,120 constructive Lucas witnesses, including \(a=0\);
- 438 actual rational sharp valuations, with repeated prime factors in \(\ell\);
- the arbitrary-large no-carry subsequences for the rescaling corollary;
- 15 exact rational-function derivative identities underlying the Hermite pole formula.

All checks pass. The unbounded claims above rest on the proofs, not the finite tests.
