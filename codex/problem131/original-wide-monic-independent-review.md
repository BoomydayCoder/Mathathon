# Independent audit of the wider-band obstruction for the original determinant

**Verdict.** The new Section 7 of `original-outer-monic-obstruction.md`, together with its monic-basis construction in Section 1, is valid. For every fixed odd \(s\ge13\), the actual primitive integer polynomial of the specified original, undepleted determinant family eventually has value larger than \(\exp(K^2/20)\) in absolute value. The proof establishes necessary coefficient-denominator costs; it is not an inference from an expensive sufficient normalizer.

The parameters remain \(K=40n\), \(N=3n\), \(h=K-N\), and numerator power six. The conclusion uses the prime number theorem and the independently certified Jensen lower bound. It does not decide the remaining cases \(s=7,9,11\), or exclude other determinant designs.

## 1. The rational basis change is legitimate after making the polynomial monic

The polynomials \(f_j=D/(t+j^2)\), \(N<j\le K\), form a rational basis: evaluation at the nodes gives a nonsingular diagonal matrix. Their change of basis need not be p-adically unimodular. This creates no problem, since its squared determinant multiplies both the determinant polynomial and its leading coefficient, and cancels in their ratio.

In this basis, off-diagonal rational integrands are integer polynomials. The diagonal integrand has integer polynomial quotient and residue
\(R_j=W(-j^2)D'(-j^2)\). Thus the affine matrix has diagonal \(X\)-coefficient \(C_j=j^{s-1}R_j\), and the monic ratio is exactly
\[
 F(X)=\det(XI+C^{-1}B).
\]
The common polynomial-quotient degree bound \(E=K+4N-2\) is correct for diagonal and off-diagonal entries.

## 2. Polynomial moments and every row valuation in K/4<p≤K

Here is a direct verification of the Bernoulli cancellation needed for all polynomial entries. The moment is
\[
 \mu_s(t^e)=\frac{(-1)^e B_{2e+2}}{s-1}\binom{2e+s}{s-2}.
\]
For \(p>s-1\), a possible Bernoulli denominator occurs only when \(2e+2=\ell(p-1)\). If
\(p>1+(2E+2)/(s-1)\), then \(1\le\ell\le s-2\). In the binomial numerator \(\prod_{i=1}^{s-2}(2e+2+i)\), the factor with \(i=\ell\) is divisible by \(p\); its factorial denominator is a unit. It cancels the possible single Bernoulli p-factor. All moments through degree \(E\) are therefore p-integral.

For \(s\ge13\), the cutoff is at most \(13K/60+5/6<K/4\) for \(K\ge40\). The other local hypotheses hold throughout the band once the fixed exceptional primes are excluded. Put \(q=\lfloor K/p\rfloor\le3\). An ordinary square class has at most \(2q+2\) nodes and at most one removed node. Since \(p^2>2K\), an ordinary same-class node difference has valuation exactly one; both factors cannot be divisible by \(p\) for an ordinary node. Consequently
\[
 v_p(C_j)=\ell_a+4\delta_a-1\le2q+5\le11<s.
\]
At a zero-class node \(j=ap\), the numerator weight is a unit, while each other zero node supplies exactly two p-factors. Thus
\[
 v_p(C_{ap})=(s-1)+2(q-1)=s+2q-3.
\]
The factors \(a^2-b^2\), with \(a,b\le3\), are units for the primes under consideration.

For ordinary nodes, the scaled pole diagonal reduces to \(-H_{\lfloor j/p\rfloor}^{(s)}\). At a zero node it reduces to \(-H_a^{(s)}+1/(2a^s)\); the \(-1/(s-1)\) pole constant vanishes after scaling. Therefore only the polynomial entries in the at most three zero rows can remain nonintegral in \(A=p^sC^{-1}B\).

Multiplying those rows by \(p^d\), \(d=(2q-3)_+\), repairs them. The respective total costs are \(E_0=0,2,9\) for \(q=1,2,3\). The repaired determinant is
\[
 \mathcal Q(Z)=p^{E_0+sh}F(Z/p^s)\in\mathbb Z_p[Z],
\]
with leading coefficient exactly \(p^{E_0}\). Its Gauss content \(v\) therefore satisfies \(0\le v\le E_0\), regardless of cancellation in the zero-row block.

## 3. The polynomial-ideal lemma handles the cancellation correctly

For each nonzero residue \(\epsilon\), every entry in its ordinary-row group belongs to \((p,Z+\epsilon)\). Determinant multilinearity consequently gives
\(\mathcal Q\in(p,Z+\epsilon)^{n_\epsilon}\). No complementary-block invertibility or particular residual minor is needed.

In the coordinate \(Y=Z+\epsilon\), membership in \((p,Y)^n\) implies that the coefficient of \(Y^j\), for \(j<n\), is divisible by \(p^{n-j}\). Integral translation preserves Gauss valuation: both the substitution and its inverse have integral coefficients. After division by the content \(p^v\), reduction modulo \(p\) therefore retains a zero of multiplicity at least \((n-v)_+\).

This proof is valid even when the polynomial degree or the multiplicities exceed \(p\). It does not differentiate modulo \(p\), divide by factorials, or use Newton identities with unjustified division by a multiple of \(p\).

For fixed \(s\), the finitely many rational numbers \(H_1^{(s)},H_2^{(s)},H_3^{(s)}\) are positive and distinct. Excluding primes in their numerators, denominators, and differences makes their reductions nonzero and distinct. These fixed primes eventually lie below \(K/4\). The nonzero groups contain
\(b_{\rm ord}=K-p+1-q\) ordinary nodes in total, in at most \(q\) groups. Hence the primitive nonzero reduction of \(\mathcal Q\) has at least
\[
 R=\max(0,b_{\rm ord}-qE_0)
\]
nonzero roots counted with multiplicity. This is a statement about the entire repaired polynomial, not just one principal block.

## 4. Retained roots force an actual coefficient denominator

Write \(F(X)=\sum_{r=0}^ha_rX^{h-r}\), \(a_0=1\). The corresponding coefficient of \(\mathcal Q\) is \(p^{E_0+sr}a_r\). A nonzero polynomial of degree at most \(h\) with at least \(R\) nonzero roots has zero-root multiplicity at most \(h-R\). Thus some surviving coefficient in the primitive reduction has index \(r\ge R\). For that coefficient,
\[
 v_p(a_r)=v-E_0-sr\le-sR.
\]
The least integer clearing the monic polynomial must therefore contain at least \(sR\) p-factors. This argument establishes a necessary denominator even if every obvious constant or minor vanishes modulo \(p\).

Substituting \(E_0=q(2q-3)_+\) gives the displayed bound
\[
 v_p(M_K)\ge s\max\{0,K-p+1-q-q^2(2q-3)_+\}.
\]
The loss from \(K-p\) is bounded independently of \(K,p\) (the stated bound 30 is safely conservative). Thus its prime-weighted contribution is only \(O_s(K)\).

## 5. The primitive-growth conclusion and its quantifiers

The prime number theorem now gives
\[
 \liminf\frac{\log M_K}{K^2}
 \ge s\int_{1/4}^1(1-u)du=9s/32.
\]
The previously reviewed Jensen calculation applies to the original positive density as well as to its fixed normalized depletion. Removing the fixed depletion factor from the leading coefficient changes only \(O_s(K)\), and the original first-frequency density lower bound is exactly the one used in that calculation. It supplies the lower rate greater than \(-1789/500\).

For \(s=13\),
\[
 9s/32-1789/500=313/4000>1/20.
\]
The margin increases with \(s\). Therefore for every fixed odd \(s\ge13\), sufficiently large \(K=40n\) satisfy
\(|M_KF(\zeta(s))|>\exp(K^2/20)\). The onset may depend on \(s\) and its finite harmonic exception set. The leading coefficient of \(F\) is one, so \(M_KF\) is primitive; every nonzero rational scalar clearing of the original determinant is an integer multiple of this primitive representative. Other primes cannot compensate for the denominator cost.

## 6. Cancellation-sensitive exact checks

I inspected the full checker `check_original_wide_monic.py`, including its precision handling in Newton identities. It reserves \(v_p(h!)\) extra digits and subtracts \(v_p(k)\) when dividing at step \(k\); each earlier coefficient has at least the currently available precision. The subsequent rescaling reserves enough further precision to determine every coefficient of \(\mathcal Q\) modulo \(p^{E_0+1}\). Thus the use of Newton identities is legitimate even when \(h>p\).

I independently reran its two \(K=40,s=13\) cases. At \(p=13\), the repaired determinant has actual Gauss content 6 despite repair cost 9, and a detected coefficient denominator exponent 367. At \(p=17\), content 1 occurs with repair cost 2, and the detected exponent is 313, exceeding the theorem's lower bound 234. These checks exercise genuine cancellation rather than a diagonal special case.

The separate `check_original_ideal_repair_independent.py` constructs exact rational matrices with deliberately dependent repaired rows. Four cases, of degrees 6, 11, and 27 and primes 3, 5, and 7, verify the content/multiplicity lemma and the coefficient-denominator conversion even with degree greater than characteristic. Their contents are respectively 0, 2, 4, and 5. It also checks the exact integral \(9/32\) and the \(s=13\) margin. All passed. These artificial matrices test the general lemma; they are explicitly distinguished from the actual-family checks above.
