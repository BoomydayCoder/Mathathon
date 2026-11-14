# Independent audit of the signed single-zeta(7) model

**Verdict.** The construction and stated obstruction in `signed-single-zeta7-model.md` pass independent review. The linear functional isolates \(1,\zeta(7)\); its displayed denominator bound is valid; the general signed-measure nonvanishing proposition applies with the stated finite window; and the total-variation lower bound has the correct constants and quantifiers. No mathematical correction was needed. This does not prove irrationality of \(\zeta(7)\), and the obstruction applies only to the specified absolute-value estimate and normalization.

I checked the statement and proof of Proposition 5.4, the growth condition, and the residue argument in Sections 5.1–5.2 of the primary source stored as `codex/sources/fischler-lai-sprang-zudilin-2026.txt`. Its large-parameter main theorem is not an input here. The nonvanishing proposition is an explicit external theorem dependency; the arithmetic and variation estimates below are checked directly.

## 1. Exact isolation of a single zeta value

The weights satisfy
\[
 45-9\cdot2^3+3^3=0,\quad45-9\cdot2^5+3^5=0,
 \quad45-9\cdot2^7+3^7=1080.
\]
Their divided-difference formula is also correct: \(w_bb^3=1080/\prod_{c\ne b}(b^2-c^2)\).

For even \(n\), the numerator fractions are stable under \(\theta\mapsto1-\theta\), and direct reversal of the Pochhammer factors gives \(R_n(-t-n)=-R_n(t)\). The multiplier
\(P(t+n/2)+P(-t-n/2)\) is invariant under this same reflection. The resulting rational function is therefore odd about \(-n/2\). Its coefficients obey
\[
 c_{\ell,n-k}=(-1)^{\ell+1}c_{\ell,k},
\]
so the total coefficients at even pole orders vanish. Its degree is at most \(-2\) for the allowed polynomial degrees, hence \(\sum_kc_{1,k}=0\) as well.

For \(\ell\ge2\), grouping the dilated lattice gives exactly
\[
 \sum_{j=1}^b\sum_{\nu\ge0}(\nu+j/b+k)^{-\ell}
 =b^\ell[\zeta(\ell)-H_{bk}^{(\ell)}].
\]
The simple-pole case must be interpreted after combining the finitely many \(k\)-terms or by a common finite cutoff: its divergent harmonic part cancels because \(\sum_kc_{1,k}=0\). No \(\zeta(1)\) or logarithmic dilation constant survives. Applying the weights then removes the order-three and order-five terms. This proves the claimed affine expression in \(\zeta(7)\), with no tacit assumption about its irrationality.

## 2. The measure and the nonvanishing theorem

Because \(n/2\) is an integer, the shifted positive lattices have exactly the periodic weights 37 on integers, \(-9\) on half-integers, and 1 on the two nonintegral third-lattices. Reflection makes the measure even, although it remains signed. The atom at a negative point has weight \(-\widehat w(x)\widehat R_n(x)\); the oddness of \(\widehat R_n\) gives precisely the reflected positive atom.

At infinity its absolute weights are \(O(|x|^{-4n-7})\). Thus multiplication by \((1+|x|)^{D-1}\), for \(D=4n+6\), leaves a summable \(O(|x|^{-2})\) tail. All moments used are absolutely defined.

The displayed polynomials \(A_n,B_n\) both have degree 12, and their quotient is exactly \(\widehat R_n(z+1)/\widehat R_n(z)\). The extra endpoint factors are necessary for the extension by zero through the gap. At a transition into the right tail, \(A_n\) vanishes: the integer endpoint is handled by \(z-n/2\), and each fractional endpoint by one factor \(z+\theta-n/2\). At a transition from the left tail into the gap, the analogous factors of \(B_n\) vanish. Off the supporting lattices, both weights are zero. This checks the recurrence for every real argument, not just for tail points.

If a rational solution \(r\) of the associated inhomogeneous equation existed, the primary lemma would imply that \((f-r)/\widehat R_n\) is period one. For sufficiently large integer \(T\), neither \(T+1/3\) nor \(-T+1/3\) is a pole of the rational function \(r/\widehat R_n\); also \(\widehat R_n\) is finite and nonzero there. The corresponding residues of \(f/\widehat R_n\) are +1 and −1. Their arguments differ by the integer \(2T\), contradicting periodicity. This is a valid nonrationality proof independent of the arithmetic nature of \(\zeta(7)\).

All hypotheses of Proposition 5.4 now hold: \(d=12\le D\) for even \(n\ge2\), and
\[
 2(M+d-1)-2<D\iff M\le2n-8.
\]
The conclusion is therefore a nonzero Hankel determinant at some \(M\le m\le M+11\). For \(n=2,4\) this choice of \(M\) is empty, not a claim of a nonzero determinant in an inadmissible window. For growing even \(n\), choose \(M\) near \(cn\) for \(0<c<2\), or \(M=2n-8\) for the endpoint ratio 2. The stated conclusions follow. They do not assert nonvanishing at every size.

The \(n=0\) example is correctly labeled separately: only finitely many moments are available, and its vanishing second and third determinants do not contradict the proposition's hypotheses.

## 3. The denominator bound counts every source of fractions

For each rational shift \(a/d\), at a prime not dividing \(d\), the product of \(n\) consecutive terms in the associated arithmetic progression contains at least \(v_p(n!)\) factors of \(p\), counted with multiplicity. At the exceptional primes, \(v_2(n!)\le n\) and \(v_3(n!)\le n/2\) prove the displayed factors \(4^n\) and \(3^{3n/2}\) sufficient. Evenness of \(n\) makes the latter exponent integral. This verifies integer-valuedness at every integer, including negative integers.

The Taylor-coefficient lemma is valid at every integer: expand an integer-valued polynomial of degree at most \(n\) in its integral binomial basis, shift by Vandermonde, and use the product formula for \(\binom uj\). An order-\(r\) Taylor coefficient has denominator dividing \(D_n^r\). Products preserve the bound at their *total* Taylor order, even when the product degree exceeds \(n\).

At a pole \(-k\), regularizing \(R_n\) produces an integer constant factor \(\pm\binom nk^7\), the three integer-valued numerator factors, and inverse powers whose nonzero integer denominators have magnitude at most \(n\). The shifted multiplier of \(P\) has the same Taylor property, since \(n/2-k\) is integral. Thus
\[
 D_n^{7-\ell}c_{\ell,k}\in\mathbb Z.
\]
The rational tails add only harmonic denominators up to \(bk\le3n\) and order \(\ell\). Consequently \(D_{3n}^7\) clears every rational tail. The coefficient of \(\zeta(7)\) is already integral, since it equals \(1080\sum_kc_{7,k}\).

Each \(U_i=(t)_i/b_i\) is a product of shifted binomial polynomials of degrees at most \(n\). Its leading coefficient is \(1/b_i\); hence changing the moment basis contributes exactly \(1/\prod_i b_i^2\) to the determinant. This proves the displayed determinant normalization, with no importation of the primary paper's arithmetic bounds for its different denominator set.

## 4. The total-variation obstruction is a lower bound, not a failed upper estimate

For the positive measure \(|\mu_n|\), the Cauchy–Binet expansion has nonnegative terms. The single unordered atom configuration
\(x_j=n/2+j/2\), \(1\le j\le m\), appears with coefficient one: its \(m!\) orderings cancel the \(1/m!\) in Andreief's integral. Its absolute lattice weights are at least one. Its Vandermonde square is
\[
 2^{-m(m-1)}\prod_{i=1}^{m-1}(i!)^2.
\]
Since \(i!/b_i\) is an integer, cancellation of these factorials against the normalizer gives exactly
\[
 J^{\rm abs}_{n,m}\ge D_n^{7m}2^{-m(m-1)}\prod_{j=1}^mR_n(j/2).
\]
There is no missing factorial from either the determinant expansion or the basis change.

Logarithmic differentiation gives \(R_n'(t)<0\) for \(t>0\), because each of the three shifted numerator sums is smaller than one copy of the sevenfold denominator sum. Uniform Stirling asymptotics for \(t=yn\), with \(y\) in a compact positive interval, give the explicit continuous rate
\[
 \frac1n\log R_n(yn)\longrightarrow
 \log108-4[(y+1)\log(y+1)-y\log y].
\]
At \(y=1\), this is \(\log(27/64)\). Thus monotonicity bounds every selected factor below by \(R_n(n)\) when \(m\le2n\); if \(m/n\to2\) and \(m\le2n+3\), use \(R_n(m/2)\), with the same limiting rate. This treats the endpoint quantifier correctly.

For the lcm bound, every prime-power term in
\[
 v_p\binom{2r}r=\sum_{a\ge1}(\lfloor2r/p^a\rfloor-2\lfloor r/p^a\rfloor)
\]
is either zero or one. There are at most \(v_p(D_{2r})\) nonzero possible terms, so \(\binom{2r}r\mid D_{2r}\). Taking \(r=\lfloor n/2\rfloor\) and the largest-term bound in the binomial sum proves
\(\log D_n\ge n\log2-O(\log n)\), without the prime number theorem.

For every sequence of admissible sizes along even \(n\) with \(m/n\to c\in(0,2]\), these estimates yield
\[
 \liminf\frac{\log J^{\rm abs}_{n,m}}{n^2}
 \ge7c\log2-c^2\log2+c\log(27/64)
 =c\log54-c^2\log2
 \ge c\log(54/4)>0.
\]
All discarded terms are \(o(n^2)\). The positivity is uniform in the sense stated for each fixed positive limiting ratio; no assertion is made for \(m/n\to0\).

The actual proven arithmetic uses \(D_{3n}\), making the associated absolute majorant still larger. Thus replacing the signed determinant by its total-variation determinant cannot give a small normalized value for this construction at these size ratios. This is a real lower bound on the proposed majorant, not merely an insufficient upper bound. It does not give a lower bound on the signed determinant, and does not exclude a substantially smaller primitive arithmetic normalization.

## 5. Independent exact checks

`python codex/computations/check_signed_zeta7_independent.py` uses Fraction arithmetic and truncated local products, independently of the original SymPy series implementation. It checks 32 basis-entry pairs at \(n=2,4,6,8\), including products of degree \(4n+4\); each individual local Taylor denominator; cancellation of all unwanted zeta coefficients; the three displayed moments; and the exact Vandermonde, lcm, and basis-factor identities. All checks passed. The original `check_signed_zeta7_model.py` was rerun as well.

The general conclusions rest on the arguments above and the explicitly identified primary nonvanishing theorem, not on those finite computations.
