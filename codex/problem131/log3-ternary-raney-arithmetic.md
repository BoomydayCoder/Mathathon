# The one-third parameter: ternary forests in the logarithmic numerator

Status: proved finite-sum factorization and arithmetic normalizers; exact finite experiment through index 15. Independent audit pending. No irrationality proof for log 3 is claimed, and no historical novelty is asserted. The positive Padé model is the established [Wilson deformation](logarithm-zeta3-wilson-deformation.md).

Fix \(\tau=1/3\), \(\sigma=1/9\). The multiplication formula, or simply the quotient of consecutive terms, gives
\[
 \delta_k=(2/3)_k(4/3)_k=\frac{(3k+1)!}{27^k k!}.
\]
Let
\[
 W_{nk}=\frac{C_{nk}^2(k!)^3}{(3k+1)!}
 =\frac{(n+k)!^2}{(n-k)!^2 k!(3k+1)!},\qquad
 A_n=\sum_{k=0}^n27^kW_{nk},\quad
 L_n=\mathop{\mathrm{lcm}}_{0\le k\le n}\operatorname{den}(W_{nk}).
\]
Write \(\mathcal N_n\) for the deformation's polynomial-difference-quotient numerator. Then
\[
 \boxed{\frac{\mathcal N_n}{27}
 =\sum_{k=1}^n\sum_{m=1}^k(-1)^{m-1}27^{k-m}W_{nk}
 R_{3,2}(m-1)\,\mathrm B(n-m+1,2m).} \tag{1}
\]
Here
\[
 R_{3,2}(m-1)=\frac{(3m-2)!}{m(m-1)!(2m-1)!}
 =\frac2{3m-1}\binom{3m-1}{m-1}\in\mathbb Z.
\]
These are the numbers of ordered pairs of full plane ternary trees with a total of \(m-1\) internal vertices: \(1,2,7,30,\ldots\). The usual forest cycle lemma proves the formula by encoding internal vertices by steps +2 and leaves by −1; among the \(3r+2\) rooted cyclic positions, two produce the valid two-root forest order. Equivalently, the generating function is \(T(x)^2\), where \(T=1+xT^3\), and the elementary forest coefficient formula gives \([x^r]T^2=2\binom{3r+2}{r}/(3r+2)\).

To prove (1), use the finite divided-product identity from the [Catalan endpoint proof](log2-endpoint-catalan-arithmetic.md), which is independent of the parameter. Before substituting the parameter it gives
\[
 \mathcal N_n=\sum_{k=1}^n C_{nk}^2\frac{(k!)^2}{\delta_k}
 \sum_{m=1}^k(-1)^{m-1}\frac{\delta_{m-1}(n-m)!}{m(n+m)!}.
\]
Substitute the displayed \(\delta\) formula and factor
\[
 \frac{(3m-2)!(n-m)!}{m(m-1)!(n+m)!}
 =R_{3,2}(m-1)\mathrm B(n-m+1,2m).
\]
This is the precise link with Fuss–Catalan/Raney combinatorics: the associated numerator contains a two-ternary-tree count, rather than merely an arithmetically integral factorial ratio.

## One lcm clears the inner integral

Let \(d_j=\mathrm{lcm}(1,\ldots,j)\). The beta integral
\[
 \mathrm B(n-m+1,2m)=\int_0^1x^{n-m}(1-x)^{2m-1}dx
\]
integrates an integer polynomial of degree \(n+m-1\le2n-1\). Hence \(d_{2n}\mathrm B\) is an integer, and (1) proves
\[
 L_n A_n\in\mathbb Z,\qquad
 L_nd_{2n}\mathcal N_n/27\in\mathbb Z. \tag{2}
\]
The endpoint target is \(Z_{1/3}=\frac{27}{2}(\log3-1)\). Put
\[
 T_n=\sum_{r=1}^n\frac9{r(9r^2-1)},\quad
 \mathcal B_n=A_n(1+2T_n/27)+\mathcal N_n/27.
\]
The lower approximant to log 3 is \(q_n=\mathcal B_n/A_n=1+2b_n/(27A_n)\). Partial fractions give the exact tail identity
\[
 \frac2{27}T_n=-\frac23H_n+
 \sum_{r=1}^n\left(\frac1{3r-1}+\frac1{3r+1}\right).
\]
For \(n\ge1\), \(d_{3n+1}\) clears its denominators, including 3. Therefore
\[
 \boxed{L_nd_{3n+1}A_n,\quad L_nd_{3n+1}\mathcal B_n\in\mathbb Z.} \tag{3}
\]
This is an exact normalizer using the minimal common outer denominator. It does not assert that (3) is the primitive common normalizer after summing.

## A simple universal factorial normalizer

The following explicit integer clears **every** outer summand:
\[
 \boxed{M_n=\frac{(3n+1)!}{(n!)^3},\qquad M_nW_{nk}\in\mathbb Z
 \quad(0\le k\le n).} \tag{4}
\]
In particular \(L_n\mid M_n\), so \(M_n\) can replace \(L_n\) in (2)–(3). To prove (4), count multiples of each prime power \(d=p^a\) in the factorials. Put \(x=n\bmod d,y=k\bmod d\). The contribution to the factorial valuation is
\[
 \left\lfloor\frac{3x+1}d\right\rfloor-
 \left\lfloor\frac{3y+1}d\right\rfloor+
 2\left\lfloor\frac{x+y}d\right\rfloor-
 2\left\lfloor\frac{x-y}d\right\rfloor.
\]
If \(x\ge y\), the first difference is nonnegative and the last floor is zero. If \(x<y\), the last contribution is +2, while the first difference is at least −2 because both floors lie between 0 and 2. Thus every prime-power contribution is nonnegative. This proof establishes a useful universal normalizer but need not be sharp enough for irrationality. In particular \(M_n\) has exponential scale \(27^n\) up to polynomial factors.

## Exact finite costs

The script independently checks (1) against direct polynomial moments and divided quotients for indices 1 through 6. At all fifteen indices it checks (2)–(4) and certifies the actual error using the positive rational series
\(\log3=2\sum_{j\ge0}(1/2)^{2j+1}/(2j+1)\), truncated at 220 terms with its geometric tail bound.

| n | outer LCD L_n | primitive denominator digits | error | primitive error |
|---:|---:|---:|---:|---:|
| 1 | 6 | 3 | 0.000127440183261 | 0.0168221041905 |
| 2 | 70 | 6 | 1.20901157523e-7 | 0.0985271893119 |
| 3 | 42 | 8 | 1.08974016177e-10 | 0.00128645918398 |
| 4 | 18018 | 13 | 9.65308998834e-14 | 0.537190567298 |
| 5 | 1716 | 14 | 8.4830542561e-17 | 0.00259052732888 |
| 6 | 184756 | 21 | 7.42282875059e-20 | 26.9667593477 |
| 7 | 125970 | 18 | 6.47832982975e-23 | 4.4069372837e-5 |
| 8 | 965770 | 27 | 5.64454180428e-26 | 9.45377716899 |
| 9 | 3120180 | 27 | 4.91242362245e-29 | 0.0122342281874 |
| 10 | 8415125460 | 37 | 4.27175784885e-32 | 46924.6513041 |
| 11 | 165002460 | 37 | 3.71239322501e-35 | 141.986796362 |
| 12 | 13082337900 | 39 | 3.2247868486e-38 | 29.8433170023 |
| 13 | 413126460 | 41 | 2.80022047069e-41 | 0.997980002953 |
| 14 | 242780649660 | 49 | 2.43085931029e-44 | 82103.4337761 |
| 15 | 19352080770 | 50 | 2.10973407533e-47 | 523.343138926 |

The first approximant is 145/132. Primitive error is below 1 at indices 1–5, 7, 9 and 13, and above 1 at the other seven tested indices. Large arithmetic fluctuations preclude a monotonicity conclusion. The top summand reciprocal, n!(3n+1)!/(2n)!², has exponential scale (27/16)^n but is generally not an integer and does not clear all outer summands; it is therefore not an established common normalizer. Neither this finite table nor the coarse universal normalizer proves the asymptotic arithmetic estimate needed for a log 3 irrationality argument.

Reproduce with `python codex/computations/explore_ternary_log3_arithmetic.py`; exact approximants, normalizers and rational error certificates are in `codex/computations/ternary_log3_arithmetic.json`.
