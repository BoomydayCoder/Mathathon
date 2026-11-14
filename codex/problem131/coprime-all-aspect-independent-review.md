# Independent audit of the all-aspect-ratio depletion obstruction

**Verdict.** `coprime-depletion-all-aspect-ratios.md` and its exact interval certificate pass independent review. The certificate proves the asserted analytic inequality for every real \(0\le\alpha<1\). Combined with the established local-prime and Jensen arguments, it proves the stated primitive-value obstruction for each fixed rational aspect ratio and each fixed odd \(s\ge3\). No mathematical correction was needed.

The result retains the integer-square pole grid, numerator power six, and dimension \(h=K-N\). It does not establish uniform asymptotics for a sequence \(\alpha_K\to1\), and does not address different numerator powers or polynomial combinations. The arithmetic assembly uses the prime number theorem.

## 1. The varying-ratio arithmetic cost

The previously audited local theorem gives bad-prime cost
\(s[K-\max(N,\lfloor p/6\rfloor)]\) for \(2K<p\le6K\), outside fixed exceptions. The leading coefficient remains a unit in this band for any \(N<K\), because its variable factors have prime divisors below \(K+N<2K\). Thus the monic coefficient-denominator lcm is an integer with these exact valuations; no compensation from other primes is possible.

For fixed \(\alpha=N/K\), the prime number theorem gives the rate \(sB(\alpha)\). Direct integration independently gives
\[
 B(\alpha)=\begin{cases}
 4/3,&0\le\alpha\le1/3,\\
 (1-\alpha)(3\alpha+1),&1/3\le\alpha<1.
 \end{cases}
\]
The threshold occurs at \(u=6\alpha\) inside the prime-density integral. The floor error changes the prime sum by only \(O(K)\). Keeping the constant \(4/3\) beyond \(\alpha=1/3\) would be wrong, and the note correctly does not do so.

After division by \(\lambda=1-\alpha\), the minimum odd-zeta cost is the continuous increasing function \(4/(1-\alpha)\) before \(1/3\), and \(9\alpha+3\) afterward. Its two values match at \(1/3\). Thus evaluation at a cell's left endpoint is a valid lower bound even when the cell crosses that threshold.

## 2. Jensen's formula and the endpoint ratio

The varying trial interval \(b=1/2+2\alpha\) stays positive and bounded for \(0\le\alpha\le1\). For every fixed \(\alpha<1\), the earlier Jensen proof applies without alteration. In particular the exact powers of \(K\) cancel between determinant and leading coefficient, leaving
\[
 R(\alpha)=\lambda[6J_\alpha(b)-J_1(b)-\pi b]
 +\lambda^2[2\log b+2\log2-3]-5I(\alpha).
\]
The expressions for \(J\) and \(I\) in the note agree with independently integrated rectangle logarithms and the leading-coefficient constant in the previous review.

At \(\alpha=0\), \(N=0\), the leading cross-product is empty and \(I(0)=0\). The remaining scaled factor \(\prod_{r=1}^K(r/K)^{s-1}\) has logarithm \((s-1)[\log K!-K\log K]=O_s(K)\), so this endpoint does not introduce a missing quadratic term. All moment and Gram arguments remain valid.

The tail replacement for \(\alpha\ge1/2\) is in the correct direction. Put
\[
 f_\alpha(x)=(x+\alpha)\log(x+\alpha)
 -(x-\alpha)\log(x-\alpha)-2\alpha.
\]
This function is continuous at \(x=\alpha\) with the usual \(0\log0=0\) convention. Its derivative in \(x\) is nonnegative on \((\alpha,1]\), so
\[
 I(\alpha)/(1-\alpha)\le f_\alpha(1).
\]
Moreover \(d f_\alpha(1)/d\alpha=\log(1-\alpha^2)\le0\). On a tail cell \([A,B]\), therefore,
\[
 -5I(\alpha)/(1-\alpha)\ge-5f_A(1).
\]
This bound contains no division by \(1-B\). It is valid when \(B=1\), for every actual parameter \(\alpha<1\) in that cell. The endpoint 1 is used only to close the analytic interval cover, not to assert a determinant of positive dimension there.

## 3. Audit of the original rational interval engine

Every acceptance comparison in `check_coprime_all_aspect_ratios.py` uses Fractions. The displayed decimal minimum is produced only afterward. I checked the following enclosure details:

- Addition, negation, multiplication by all endpoint products, and division by a strictly positive interval are inclusion-preserving. Correlated appearances of \(\alpha\), \(b\), and \(\lambda\) can widen bounds but cannot invalidate them.
- After logarithm range reduction to \([1,2]\), the omitted atanh series is positive and bounded by the stated geometric tail. Negative multiples of the interval for \(\log2\) reverse endpoints correctly.
- The arctangent next-term enclosure is valid also for negative reduced arguments. In the addition reduction, \((x-1/2)/(1+x/2)\in[-1/2,1/3]\) for \(0\le x\le1\), and the principal-angle identity is valid throughout. Larger positive arguments use the reciprocal identity with the enclosed π value.
- Machin's formula and these alternating enclosures give a rigorous π interval.
- For \(H(v)=v^2\log v/2-v^2/4\), the derivative \(v\log v\) proves that the only interior minimum is at 1. Taking endpoint bounds and including \(H(1)=-1/4\) when appropriate encloses the entire range. The value at zero is the correct continuous value 0.
- The rewritten \(J_t\) uses \(\arctan(b/t)=\pi/2-\arctan(t/b)\) and extends continuously to \(t=0\). It avoids a hidden division by zero at the leftmost cell.
- Non-tail cells lie in \([0,1/2]\), so their division by \(\lambda\) has denominator at least \(1/2\). Tail cells use the separate monotone bound above.

The adaptive process begins with two genuine intervals and only replaces an interval by its two exact rational halves. The final explicit adjacency test additionally proves there are no gaps or overlaps of interiors. The depth and trial limits are termination safeguards, not assumptions about the sign. All 76 cells were accepted by strict inequalities.

## 4. Independent verification of every saved cell

I saved the generator's exact partition as `codex/computations/coprime_all_aspect_intervals.json`. The separate verifier `check_all_aspect_intervals_independent.py` uses its endpoints and tail flags only. It ignores all reported lower bounds and recomputes them with a separate interval implementation, 60-term logarithm bounds, and 50-term arctangent bounds using different angle reductions. It uses no floating-point arithmetic.

The verifier checks:

1. every cell has positive length and the cells form an exact unbroken cover of \([0,1]\);
2. every tail flag is appropriate to its cell, including the final endpoint;
3. all interval function bounds and the relevant cost lower bound are recomputed;
4. each resulting lower bound is strictly larger than \(1/5\).

All 76 cells passed. In fact every independently computed ratio lower bound exceeds the exact rational number \(2019/10000>1/5\). This extra reserve is not required for the theorem, but verifies that the strict threshold is not a floating-point artifact.

Therefore the finite computation rigorously proves
\[
 3B(\alpha)+R(\alpha)>\frac{1-\alpha}{5}
 \qquad(0\le\alpha<1).
\]
It is a verified interval inequality, not a point-sampling argument or a claim that the chosen trial interval is optimal.

## 5. Quantifiers in the resulting no-go theorem

Fix a rational \(\alpha<1\) and a fixed odd \(s\ge3\), and take \(K\) through multiples of the denominator of \(\alpha\). The local-prime and real limits then give
\[
 \liminf\frac{\log P_K(\zeta(s))}{K^2}
 \ge sB(\alpha)+R(\alpha)
 \ge3B(\alpha)+R(\alpha)>(1-\alpha)/5.
\]
The strict gap implies the claimed eventual bound on the primitive value. Every nonzero rational scalar clearing of the original determinant is an integer multiple of its primitive polynomial, so it cannot make the values tend to zero.

The analytic inequality is uniform over the real interval, but the errors in the determinant and prime asymptotics have only been invoked for fixed parameters. Hence it does not by itself justify a varying \(\alpha_K\to1\) statement. The note explicitly preserves this distinction, and its quantifiers are correct.

## 6. Reproduction

Run:

```
python codex/computations/check_coprime_all_aspect_ratios.py
python codex/computations/check_all_aspect_intervals_independent.py
```

Both passed. The generator reports 76 accepted cells, 150 evaluations, and maximum depth seven. The independent verifier recomputes the bounds on the persisted partition and reports its strictly rational safety margin.
