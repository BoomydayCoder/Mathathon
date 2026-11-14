# Independent audit of the primitive coprime-depletion obstruction

**Verdict.** The theorem in `coprime-depletion-primitive-obstruction.md` is valid for the stated family: for each fixed odd \(s\ge3\), its primitive integer determinant, signed positively at \(\zeta(s)\), eventually exceeds \(\exp(K^2/3)\). This is an obstruction to **nonzero rational scalar normalization of that family**, not merely a mismatch between sufficient bounds. The monic normalization is what upgrades the bad-prime information to an actual lower bound on primitive values.

I independently derived the Jensen limit and its constant, checked the leading-coefficient and prime-number-theorem assembly, and verified the exact local-prime argument used as input. The separate local audit by the connections agent covers that prime theorem in further detail. This result explicitly uses the prime number theorem. It says nothing against other pole grids, numerator powers, dimension ratios, or combinations that change the polynomial family.

## 1. The monic step closes the arithmetic loophole

Let \(E=D_K/D_N\). After polynomial division, the surviving pole at \(-r^2\), \(N<r\le K\), has residue
\[
 \rho_r=\frac{D_N(-r^2)^5}{E'(-r^2)}.
\]
Its coefficient of \(X\) is \(C_sr^{s-1}\), where \(C_s=(2^{s-1}-1)(3^{s-1}-1)\). The top-degree matrix is a square Vandermonde congruence of that diagonal matrix. The squared Vandermonde cancels the product of the derivatives \(E'\), including the sign, giving
\[
 \ell_K=(-1)^{h(h-1)/2}C_s^h
       \prod_{r=N+1}^K r^{s-1}D_N(-r^2)^5.
\]
It is a nonzero integer and proves that the degree is exactly \(h\). Its absolute value is
\[
 C_s^h\prod_{r=N+1}^K r^{s-1}
             \prod_{j=1}^N(r^2-j^2)^5.
\]

For the monic rational polynomial \(F_K=\Delta_K/\ell_K\), every rational scalar that makes it an integer polynomial is itself an integer, by inspecting its leading coefficient. Let \(L_K\) be the least positive such integer, equivalently the lcm of its coefficient denominators. Then \(L_KF_K\) is primitive: a common prime factor would divide its leading coefficient \(L_K\) and allow a smaller clearing integer. Thus the primitive value of the original determinant is exactly
\[
 P_K(\zeta(s))=L_K|F_K(\zeta(s))|.
\]
The absolute value is necessary, since the monic normalization may reverse the sign. The original Gram determinant is strictly positive.

This argument prevents numerator factors at unrelated primes from offsetting the new denominator cost. Such an offset would have been possible when reasoning directly about an arbitrary rational scalar multiplying the unnormalized determinant.

## 2. Local input and its prime-density consequence

For \(2K<p\le6K\), all node differences \((r-j)(r+j)\) are units modulo \(p\). The pole constants have valuation \(-s\) exactly when \(r>p/6\), outside the listed fixed prime exceptions. In the upper subrange \(p/3<r<p/2\), the leading coefficient is
\[
 3^{s-1}-6^{s-1}(1+2^{-s})
 =-\frac{3^{s-1}(2^s-1)}2,
\]
so the exclusion of divisors of \(2^s-1\) is precisely what is required. The other nonzero subrange has coefficient \(-6^{s-1}\). All polynomial-quotient moments are integral at these primes, because their denominators are supported at 2 and 3.

Conjugating by the unit Vandermonde leaves an integral constant matrix plus an affine diagonal. Multiplying each bad row by \(p^s\) gives, modulo \(p\), a row with only its diagonal constant nonzero. Expanding along these rows leaves a good block whose highest \(X\)-coefficient is a product of units. Thus the Gauss valuation and the specified coefficient valuation are exactly \(-sb\), not merely at least that value.

Apart from the fixed factor \(C_s\), every prime divisor of \(\ell_K\) is at most \(K+N<2K\). The leading coefficient is consequently a unit at every prime used in the band. Since \(N/K=3/40<1/3\), the monic clearing integer satisfies
\[
 v_p(L_K)=s\bigl(K-\lfloor p/6\rfloor\bigr)
\]
throughout the band, after the fixed exceptional primes have fallen below it.

The prime number theorem gives
\[
 \sum_{2K<p\le6K}\left(K-\frac p6\right)\log p
 =K^2\int_2^6(1-u/6)du+o(K^2)
 =\frac43K^2+o(K^2).
\]
The floor error is \(O(K)\), since \(\sum_{p\le6K}\log p=O(K)\). Other primes have nonnegative valuations in the positive integer \(L_K\), hence
\[
 \liminf\frac{\log L_K}{K^2}\ge\frac{4s}{3}.
\]
There is no cancellation across primes left in this statement.

## 3. Independent Jensen derivation and complete power accounting

Write \(\alpha=3/40\), \(\lambda=1-\alpha\), and \(b=1/2\). In fact the depleted density itself contains the positive first-frequency term, so
\[
 \sigma_s(y)\ge c_s^* y^s e^{-2\pi y},
 \qquad c_s^*=2(2\pi)^{s-1}/(s-1)!.
\]
The comparison with the original density in the main proof is also valid, but this direct lower bound suffices.

Andreief's formula, restriction to \(0<y_i<Kb\), and \(y_i=Kx_i\) give a positive multiple of an expectation over independent uniform \(x_i\in(0,b)\). Jensen applies to its logarithm: \(\log x_i\), \(\log|x_i-x_j|\), and \(\log(x_i+x_j)\) are integrable under this distribution. Coincident points and the endpoint zero have measure zero.

The exact powers of \(K\) in this lower integral are
\[
 \underbrace{2h(h-1)}_{\text{Vandermonde}}+
 \underbrace{(10N-2h)h}_{D_N^5/E}+
 \underbrace{h}_{dy}+\underbrace{sh}_{y^s}
 =10Nh+(s-1)h.
\]
These equal the powers extracted from the exact leading coefficient. Thus there is no residual \(K^2\log K\) term. The factor \(1/h!\), fixed density constants, interval length, the remaining \(x_i^s\), and the scaled factors \((r/K)^{s-1}\) cost only \(O_s(K\log K)\), which is negligible here.

For independent uniform \(X,Y\) on \((0,b)\),
\[
 \mathbb E\log|X^2-Y^2|=2\log(2b)-3=-3.
\]
The pair contribution is \(h(h-1)\) times this expectation, not half that value: the Vandermonde is squared. This gives \(-3\lambda^2\) in the limit.

Define the rectangle average
\[
 J(a)=\frac1b\int_0^b\int_0^a\log(x^2+u^2)du\,dx
 =a\log(a^2+b^2)-3a+\frac{a^2}{b}\arctan(b/a)+b\arctan(a/b).
\]
The limiting external-field contribution is
\(\lambda[6J(\alpha)-J(1)-\pi b]\), since \(2\pi\mathbb EX=\pi b\). This independently reproduces the main proof's \(V\).

## 4. The logarithmic Riemann limits are legitimate

For the numerator and denominator node sums, first average over the uniform variable. The resulting function of the node coordinate is
\[
 g(u)=\frac1b\int_0^b\log(x^2+u^2)dx
 =\log(b^2+u^2)-2+\frac{2u}{b}\arctan(b/u),
\]
with continuous extension \(g(0)=2\log b-2\). Its Riemann sums therefore converge ordinarily. The apparent corner singularity of the two-variable logarithm causes no issue after this averaging.

For the leading coefficient, a convenient uniform verification avoids invoking a Riemann-sum theorem for an unbounded function. Use the exact identity
\[
 |D_N(-r^2)|=\frac{(r+N)!}{r(r-N-1)!}
\]
and the uniform integer Stirling bound
\(\log(q!)=q\log q-q+O(\log(q+2))\), with \(0\log0=0\). For \(v=r/K\), uniformly over \(N<r\le K\),
\[
 \log|D_N(-r^2)|-2N\log K
 =K[(v+\alpha)\log(v+\alpha)
       -(v-\alpha)\log(v-\alpha)-2\alpha]+O(\log K).
\]
The bracket is continuous at \(v=\alpha\), and the summed error is only \(O(K\log K)\). It equals \(\int_0^\alpha\log(v^2-u^2)du\). Thus the leading-coefficient constant is exactly
\[
 L=5\int_\alpha^1\int_0^\alpha\log(v^2-u^2)du\,dv.
\]

An independently convenient primitive is
\(F(z)=z^2\log z/2-3z^2/4\), for which \(F''(z)=\log z\). Splitting the logarithm into \(\log(v-u)+\log(v+u)\) gives
\[
 L=5[F(1+\alpha)-F(1-\alpha)-F(2\alpha)].
\]
This equals the main note's expression using its different primitive \(H\), because \(F=H-z^2/2\). The finite correction \(-2\alpha(1-\alpha)\) in that expression is therefore correct.

Combining the pieces independently gives
\[
 \liminf\frac{\log|F_K(\zeta(s))|}{K^2}
 \ge -3\lambda^2+\lambda[6J(\alpha)-J(1)-\pi b]-L
 =:J_*.
\]
The constant is independent of each fixed \(s\), since all \(s\)-dependent prefactors contribute below the quadratic scale.

## 5. Exact constants and the actual primitive-value conclusion

My independent rational interval calculation uses the direct primitive \(F\) above, a symmetric atanh-series error bound for logarithms without dyadic reduction, and 60-term alternating arctangent bounds. It yields
\[
 -3.577448593858<J_*<-3.577448593857,
 \qquad J_*>-18/5.
\]
The independent component values agree with the parent's certificate: the pair contribution is exactly \(-2.566875\), the external contribution is approximately \(-1.56878207356247\), and the subtraction of the leading-coefficient constant contributes approximately \(+0.558208479705224\). Only rational inequalities are used for acceptance.

Consequently
\[
 \liminf\frac{\log P_K(\zeta(s))}{K^2}
 \ge\frac{4s}{3}+J_*>rac{4s}{3}-\frac{18}{5}
 \ge\frac25>\frac13.
\]
This proves the asserted eventual \(\exp(K^2/3)\) lower bound. The threshold may depend on the fixed \(s\); no explicit uniform threshold in \(s\) is claimed.

A nonzero rational scalar multiple of a primitive integer polynomial is integral only if that scalar is an integer, by Bezout applied to its coefficients. Therefore every nonzero scalar clearing of the original determinant has value of absolute size at least \(P_K(\zeta(s))\). This proves the intended scalar-normalization obstruction. It does not prohibit operations that produce a different polynomial family.

## 6. Reproduction

Run `python codex/computations/check_coprime_obstruction_independent.py`. It checks the strict Jensen margin with Fraction intervals, 234 exact leading-coefficient/sign/scaling cases, and the exact prime-density integral \(4/3\). All checks passed. The parent's `certify_coprime_depletion_obstruction.py` was rerun and gives the same enclosing interval. These checks support the all-parameter proofs above; they do not replace the prime number theorem or the separately proved local determinant theorem.
