# Combined independent audit: all fixed coprime pairs and fixed aspect ratios

**Verdict.** `general-coprime-depletion-band.md` is correct. Combining it with the independently verified all-aspect-ratio inequality proves the following actual primitive-value obstruction:

> Fix odd \(s\ge3\), coprime integers \(a,b\ge2\), and a rational \(0\le\alpha<1\). Set \(q=ab\), \(N=\alpha K\), and \(h=K-N\), taking \(K\) through multiples of the denominator of \(\alpha\). For the normalized depleted density and the integer-square-node determinant with numerator \(D_N^6\), its primitive integer polynomial, signed positively at \(\zeta(s)\), eventually satisfies
> \[
> P_K(\zeta(s))>\exp((1-\alpha)K^2/5).
> \]

Thus nonzero rational scalar normalization cannot produce values tending to zero for any of these fixed choices. This uses the prime number theorem. It does not cover varying depletion parameters, a moving ratio approaching one, different numerator powers, different grids, or polynomial combinations changing the family.

## 1. The step-value bound has the necessary sign and uniformity

Write \(m=s-1\). If \(0<x<1/2\), the undepleted term \(H_{\lfloor x\rfloor}^{(s)}\) is zero. When \(qx<1\), every term is zero. When \(qx\ge1\), the positive leading harmonic sum contains its first term, giving
\[
 \frac{R_s(x)}{q^m}\ge1-(a^{-m}+b^{-m})\zeta(s)>79/144.
\]
The last strict inequality is justified for all allowed parameters: \(m\ge2\), coprimality rules out \(a=b=2\), hence \(a^{-m}+b^{-m}\le13/36\); and the elementary decreasing-integral estimate gives \(\zeta(s)\le\zeta(3)<5/4\). Thus the sign is genuinely uniform, not sampled or extrapolated from small pairs.

For \(k=\lfloor qx\rfloor\), the identities \(\lfloor ax\rfloor=\lfloor k/b\rfloor\) and \(\lfloor bx\rfloor=\lfloor k/a\rfloor\) hold exactly. The range is \(1\le k\le\lfloor(q-1)/2\rfloor\). Therefore all nonzero step values belong to the finite positive rational list in the note.

The exceptional set correctly includes prime divisors of **both numerators and denominators** of that list, as well as \(2qC\). Consequently every surviving step value is a p-adic unit. Positivity ensures that no list value is zero, so this exclusion really is finite. It is a sufficient exception set, not a claim of minimality or absence of special prime cancellations.

## 2. The exact local determinant theorem survives the generalization

Separating multiples of \(p\) in each harmonic sum is an exact identity, even if multiples of \(p^2\) occur. The nonmultiple remainder is p-integral. The final rational pole term is p-integral because \(p\nmid2q\) and \(r<p/2\). Thus a pole is bad exactly when \(qr>p\), with constant valuation exactly \(-s\); equality \(qr=p\) is excluded by \(p\nmid q\). The number of surviving bad poles is precisely
\(K-\max(N,\lfloor p/q\rfloor)\).

All affine \(X\)-coefficients are units, since \(p\nmid C\). The polynomial moments are p-integral because they belong to \(\mathbb Z[1/q]\). The unit Vandermonde and residue factors remain valid under \(p>2K\). After the bad rows are multiplied by \(p^s\), expansion modulo \(p\) leaves a nonzero highest-degree coefficient from the complementary good block. This proves both the Gauss valuation and the stated individual coefficient valuation. The case in which every row is bad uses the empty complementary determinant, equal to one, and causes no exception.

The exact leading coefficient is the stated nonzero integer. Its variable prime factors lie below \(K+N<2K\), and its fixed factor \(C\) is excluded. Therefore it is a unit at every prime in this theorem. Monic normalization turns the required clearing scalar into a positive integer, so other primes cannot offset the band contribution. This is the essential distinction between an actual primitive obstruction and an expensive sufficient denominator estimate.

## 3. Prime density and comparison with the (2,3) cost

The piecewise integration is correct:
\[
 B_q(\alpha)=\begin{cases}
 (q-2)^2/(2q),&\alpha\le2/q,\\
 (1-\alpha)[q(1+\alpha)/2-2],&\alpha>2/q.
 \end{cases}
\]
Both formulas agree at the breakpoint. For fixed parameters the floor error contributes \(O(K)\), the finite exception set eventually leaves the band, and the prime number theorem yields the denominator rate \(sB_q(\alpha)\). The statement with \(N/K\to\alpha\) is also valid for this arithmetic bound: the max function is Lipschitz in \(\alpha\), giving only \(o(K^2)\) change.

Coprimality and \(a,b\ge2\) imply \(q\ge6\). On \([2,6]\), \(u/q\le u/6\), so the integrand defining \(B_q\) is at least the integrand for \(B_6\). Its remaining integral on \([6,q]\) is nonnegative. Hence \(B_q(\alpha)\ge B_6(\alpha)\) for every real \(0\le\alpha<1\). This is a comparison of the established asymptotic band costs; it does not assert an ordering of full coefficient-denominator lcms.

## 4. The monic real estimate transfers without a new equilibrium calculation

The retained positive frequencies are those divisible by neither \(a\) nor \(b\). For composite dilation factors this should not be confused with requiring frequencies to be relatively prime to \(q\); that stronger assertion is unnecessary. The first frequency always remains, giving strict positivity and the original exponential tail.

The normalized density ratio to \(w_s\) tends to \((1-1/a)(1-1/b)>0\) at zero and to one at infinity. It is positive and continuous in between, so any two fixed coprime choices have ratios bounded above and below by positive constants independent of \(K,N,h\).

The corresponding positive Gram matrices therefore bound one another in matrix order; conjugation by a positive square root gives determinant ratios between \(c^h\) and \(C_0^h\). Their leading coefficients differ exactly by \((C_{s;a,b}/C_{s;2,3})^h\). Taking absolute monic values consequently changes the logarithm by only \(O_{s,a,b}(h)=O(K)\), for every permitted \(N<K\). The previously established Jensen lower rate \(R(\alpha)\) transfers unchanged on the quadratic scale.

Combining these facts with the independently certified inequality from `coprime-all-aspect-independent-review.md` gives
\[
 \liminf\frac{\log P_K(\zeta(s))}{K^2}
 \ge sB_q(\alpha)+R(\alpha)
 \ge3B_6(\alpha)+R(\alpha)>(1-\alpha)/5.
\]
This proves the combined theorem above for every fixed rational aspect ratio, extending the fixed \(3/40\) conclusion explicitly written in the general-pair note.

## 5. Checks and scope

I reran `check_general_coprime_depletion_band.py`: 492 exact positive step values, 84 exact band integrals, and 41 full determinant-prime cases passed, including the specified coefficient and monic denominator valuations. The separate independent all-aspect verifier recomputed all 76 interval cells without floating-point arithmetic. The all-parameter proofs above, rather than those finite tests, establish the general conclusion.

There is no claim that primitive polynomials from different pairs coincide, that their coefficients or denominator lcms are termwise ordered, or that fixed exceptional sets remain harmless when \(a,b,s\) grow with \(K\). The constants in the real comparison and the asymptotic onset may depend on all fixed parameters. These restrictions are necessary and are preserved in the combined theorem.
