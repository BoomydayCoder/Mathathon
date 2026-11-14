# Positive frequency depletion and cyclotomic logarithmic derivatives

**Status: proved; independent exact checks pass.** The coprime double-depletion construction was proposed by the parent. The cyclotomic interpretation below was developed independently in the review. It is an exact structural connection, not a new irrationality theorem or a historical novelty claim.

## The construction

Fix odd \(s\ge3\), and let
\[
 w_s(y)=\frac{2(2\pi)^{s-1}}{(s-1)!}y^s
 \sum_{n\ge1}n^{s-1}e^{-2\pi ny},\qquad
 \mu_s(R)=\int_0^\infty R(y^2)w_s(y)\,dy.
\]
For coprime integers \(a,b\ge2\), define
\[
 W_{s;a,b}(y)=ab\,w_s(y/(ab))-a\,w_s(y/a)-b\,w_s(y/b)+w_s(y).
\]
Then \(W_{s;a,b}>0\), and every polynomial moment
\(\int y^{2e}W_{s;a,b}(y)\,dy\) is an integer.

The positivity is elementary inclusion-exclusion: with all exponentials written at scale \(ab\), the remaining frequency coefficient is
\[
 1-\mathbf1_{a\mid n}-\mathbf1_{b\mid n}+\mathbf1_{ab\mid n}
 =\mathbf1_{a\nmid n,\ b\nmid n}.
\]
The term \(n=1\) remains, so positivity is strict and the exponential decay rate is \(2\pi/(ab)\).

## Why a cyclotomic polynomial proves integrality

Put
\[
 R_{a,b}(q)=\frac{(q^{ab}-1)(q-1)}{(q^a-1)(q^b-1)}.
\]
Since \(\gcd(a,b)=1\), factorization of \(q^n-1\) gives
\[
 R_{a,b}(q)=\prod_{\substack{d\mid ab\\d\nmid a,\ d\nmid b}}\Phi_d(q)
 \in\mathbb Z[q],\qquad R_{a,b}(1)=1.
\]
In particular its degree is \((a-1)(b-1)\). Define
\[
 c_k=\left.\frac{d^k}{dz^k}\log R_{a,b}(e^z)\right|_{z=0}.
\]
These numbers are integers. Indeed, every derivative \(r_k\) of \(R_{a,b}(e^z)\) at zero is an integer, \(r_0=1\), and differentiation of \(R'=R(\log R)'\) yields
\[
 c_n=r_n-\sum_{j=1}^{n-1}\binom{n-1}{j-1}c_jr_{n-j}.
\]
This proves integrality inductively without Bernoulli denominator theorems. The polynomial may have negative coefficients; no probability interpretation is required here.

For the uniform-distribution generating function
\(M_a(z)=(e^{az}-1)/(a(e^z-1))\), one has
\[
 R_{a,b}(e^z)=M_a(bz)/M_a(z).
\]
Expanding its logarithmic derivative with the Bernoulli generating function gives, for even \(k\ge2\),
\[
 c_k=\frac{(a^k-1)(b^k-1)B_k}{k}.
\]
On the other hand, direct integration gives
\[
 \mu_s(t^e)=\frac{(-1)^eB_k}{s-1}\binom{k+s-2}{s-2},\quad k=2e+2.
\]
Dilation therefore proves the exact identity
\[
 \boxed{\int_0^\infty y^{2e}W_{s;a,b}(y)\,dy
 =(-1)^e c_k\binom{k+s-2}{s-1}\in\mathbb Z_{>0}.}
\]
The positive integral also explains the alternating signs of these even logarithmic derivatives.

## The smallest example is Phi_6

For \(a=2,b=3\),
\[
 R_{2,3}(q)=q^2-q+1=\Phi_6(q),\qquad
 \log R_{2,3}(e^z)=z+\log(2\cosh z-1).
\]
The first even derivatives are \(c_2=2,c_4=-10,c_6=182\). Thus the first three moments are
\[
 2s,\qquad 10\binom{s+2}{s-1},\qquad
 182\binom{s+4}{s-1}.
\]
For \(s=3\), these are \(6,100,3822\). Their integrality follows from one small integer polynomial with value 1 at 1.

## What arithmetic was actually improved

The original undepleted functional has no fixed uniform polynomial denominator. For any prime \(p>s\), choose \(k=p(p-1)\). Then \(p-1\mid k\), so \(v_p(B_k)=-1\), while \(p\nmid(s-1)\) and Lucas's theorem gives
\[
 \binom{k+s-2}{s-2}\equiv1\pmod p.
\]
Hence the original moment has valuation exactly \(-1\) at arbitrarily large primes.

Single depletion by an integer \(\ell\) already changes this to the fixed optimal denominator
\(\prod_{p\mid\ell}p^{1+v_p(s-1)}\).
The coprime double depletion removes even that fixed denominator through an exact positive operation. Multiplying the single-depleted functional by its fixed denominator would also have produced integral moments, so the double construction does **not** automatically improve the leading quadratic arithmetic budget in a growing determinant.

Restoring the original decay rate gives the positive density
\[
 \Sigma_{s;a,b}(y)=\frac1{ab}W_{s;a,b}(ab,y)
 =w_s(y)-\frac1a w_s(ay)-\frac1b w_s(by)+\frac1{ab}w_s(ab,y).
\]
Its \(e\)-th moment is the preceding integer divided by \((ab)^k\). Thus the prime support remains restricted, but denominator exponents may grow. Pole evaluations, determinant normalization and the real bound must still be audited separately. This note makes no claim that these costs are harmless.

## Checks and references

The standard-library script `codex/computations/check_integer_depletion_independent_review.py` independently computes the integer quotient polynomial, its exponential derivatives and logarithmic derivatives, and checks 5,100 double-depleted moments for coprime pairs up to 12. All checks pass. The proofs above establish the general statements.

See `integer-depletion-cumulants.md` for the parent construction and exact single-depletion denominator theorem, and `integer-depletion-independent-review.md` for the independent all-prime, rescaling and Hurwitz-constant audit. The external ingredients here are the standard cyclotomic factorization, Bernoulli generating series, Lucas's theorem and von Staudt–Clausen; the last two are used only to explain the undepleted denominator contrast, not to prove double-depletion integrality.
