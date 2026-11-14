# A higher-multiplicity polynomial family repairs the zeta(5) boundary terms

**Status.** The polynomial, Mellin, arithmetic, cancellation and exponential-bound theorems below are proved and independently reviewed in `multiplicity-lifted-independent-review.md`. For \(\zeta(5)\) the new forms are exponentially small and infinitely often nonzero. The available denominator upper bound is too large to finish an irrationality argument; no assertion that the *actual* denominator growth rules it out is made. A precise new-prime obstruction is proved instead. This is a concrete modification of the polynomial family, not a fixed harmonic correction to the original Apéry forms. Historical novelty has not been investigated.

## 1. Impose the desired Mellin zeros before choosing the polynomial

Fix positive integers \(n,m\). Define
\[
 J_{n,m}(x)=\sum_{k=0}^{mn}a_kx^k,
 \qquad a_k=(-1)^k\binom{mn}{k}\binom{n+k}{k}^{m}.
\tag{1}
\]
All coefficients are integers and \(J_{n,m}(0)=1\). These are iterated Rodrigues polynomials:
\[
 J_{n,m}=\mathcal R_n^m[(1-x)^{mn}],
 \qquad\mathcal R_n f=\frac1{n!}\frac{d^n}{dx^n}(x^nf).
\tag{2}
\]
Indeed \(\mathcal R_n(x^k)=\binom{n+k}{k}x^k\).

Their Mellin transforms have prescribed multiple zeros:
\[
 \boxed{F_{n,m}(t):=\int_0^1x^{t-1}J_{n,m}(x)\,dx
 =(-1)^{mn}\frac{(mn)!}{(n!)^m}
 \frac{\prod_{j=1}^n(t-j)^m}{\prod_{k=0}^{mn}(t+k)}.}
\tag{3}
\]
The integral initially requires \(\Re t>0\); the identity extends rationally. One direct proof compares residues at \(t=-k\): the residue of the right side is exactly (1), and both sides tend to zero at infinity. This also fixes the sometimes easy-to-miss sign \((-1)^{mn}\).

Thus, for \(1\le j\le n\) and \(0\le r<m\),
\[
 \int_0^1x^{j-1}(\log x)^rJ_{n,m}(x)\,dx=0.
\tag{4}
\]
This is simultaneous orthogonality for the weights \(1,\log x,\ldots,(\log x)^{m-1}\). It supplies a concrete multiple-orthogonality construction rather than merely changing the exponent of a logarithmic kernel.

### Minimal degree within this multiplicity design

If a nonzero polynomial \(P\) of degree at most \(d\) has a Mellin transform with zeros of order at least \(m\) at \(1,\ldots,n\), then \(d\ge mn\). Indeed
\[
 \int_0^1x^{t-1}P(x)\,dx=\sum_{k=0}^d\frac{p_k}{t+k}
\]
has numerator of degree at most \(d\) over the common denominator \(\prod_{k=0}^d(t+k)\); its numerator has at least \(mn\) prescribed roots counted with multiplicity.

At degree \(mn\), normalizing \(P(0)=1\) makes the numerator, its scalar, and hence \(P\) unique. It is exactly \(J_{n,m}\). For an asymmetric pairing whose product of Mellin transforms has zeros of total order at least \(2m\) at each of these \(n\) points, the same count gives \(\deg P+\deg Q\ge2mn\).

These are minimal-degree assertions for **imposed zero multiplicities**. They do not rule out a different design making a high derivative vanish through cancellation while retaining lower-order zeros.

## 2. Exact single-zeta forms and elimination of the root contributions

Use the previously proved monomial kernel
\[
 \mathcal I_s(P,Q)=\frac1{(s-2)!}\int_0^1\!\int_0^1
 \frac{P(x)Q(y)(-\log(xy))^{s-2}}{1-xy}\,dx\,dy.
\]
Put
\[
 \mathscr L_{n,m}=\mathcal I_{2m+1}(J_{n,m},J_{n,m}).
\]
Then
\[
 \boxed{\mathscr L_{n,m}=A_{n,m}\zeta(2m+1)-B_{n,m},
 \quad A_{n,m}=2m\sum_{k=0}^{mn}a_k^2>0,}
\tag{5}
\]
where \(B_{n,m}\) is rational and
\[
 D_{mn}^{,2m+1}B_{n,m}\in\mathbf Z.
\tag{6}
\]
The explicit rational part, with \(s=2m+1\), is
\[
 B_{n,m}=(s-1)\sum_k a_k^2H_k^{(s)}
 -2\sum_{k<\ell}a_ka_\ell
 \frac{H_\ell^{(s-1)}-H_k^{(s-1)}}{\ell-k}.
\tag{7}
\]
The denominator claim follows term by term, as in `monomial-kernel-design.md`.

Let \(R_{n,m}=F_{n,m}^2\). Expansion of the geometric kernel and Mellin differentiation give
\[
 \mathscr L_{n,m}=-\frac1{(2m-1)!}
 \sum_{t=1}^\infty R_{n,m}^{(2m-1)}(t).
\tag{8}
\]
Since \(R_{n,m}\) has a zero of order \(2m\) at every \(t=1,\ldots,n\),
\[
 \boxed{R_{n,m}^{(2m-1)}(j)=0\quad(1\le j\le n).}
\tag{9}
\]
The sum in (8) is therefore exactly the tail \(t>n\), with **no rational root correction subtracted afterward**. This removes the particular source of polynomial boundary contributions in the unmodified higher-derivative family. The recurrence itself changes; no claim is made that these forms obey the old homogeneous Apéry recurrence.

For \(m=2\), the first polynomials are
\[
 J_{1,2}=1-8x+9x^2,
 \qquad J_{2,2}=1-36x+216x^2-400x^3+225x^4.
\]
The first exact form is
\[
 \mathscr L_{1,2}=584\zeta(5)-\frac{9689}{16}.
\tag{10}
\]

## 3. A rigorous exponential bound for the zeta(5) case

In this section \(m=2\), write \(F=F_{n,2}\), \(R=F^2\), and \(\mathscr L_n=\mathscr L_{n,2}\). Set
\[
 a_n=\frac1{n^3\binom{3n}{n}}.
\]
This is the coefficient of \((t-n)^2\) in \(F(t)\), as follows immediately from (3).

**Theorem.**
\[
 \boxed{\mathscr L_n=O\left(n^{-7/2}\left(\frac{16}{729}\right)^n\right).}
\tag{11}
\]
The proof supplies an explicit bound, not a numerical fit.

### Cotangent summation and the endpoint cancellation

Let \(K(z)=\pi\cot(\pi z)\). Close an upward vertical contour \(\Re z=n\) to the right. The orientation is clockwise. The rational function has no poles there, and the residue of \(KR'''\) at \(n\) is zero, since \(R'''(n)=0\). In fact the apparent singularity is removable. Using large half-integer right edges and horizontal edges tending to infinity justifies the residue sum: for fixed \(n\), \(R'''(z)=O(|z|^{-5})\).

After three integrations by parts on the vertical line,
\[
 \mathscr L_n=-\frac1{12\pi i}
 \int_{n-i\infty}^{n+i\infty}K'''(z)R(z)\,dz.
\tag{12}
\]
All apparent singularities on this line are removable because \(R\) has a fourth-order zero at \(n\), exactly matching the fourth-order pole of \(K'''\). There are no hidden endpoint terms.

For \(z=n+iy\),
\[
 |K'''(z)|=2\pi^4\left(2\operatorname{csch}^2(\pi y)
 +3\operatorname{csch}^4(\pi y)\right).
\tag{13}
\]
The rational product (3) and Euler's standard hyperbolic-sine product give
\[
 |F(n+iy)|\le
 \frac{a_n|y\sinh(\pi y)|}{\pi}
 \left(1+\frac{y^2}{9n^2}\right)^{-(2n+1)/2}.
\tag{14}
\]
Indeed the numerator ratio is
\(y^2\prod_{j=1}^{n-1}(1+y^2/j^2)\), bounded by
\(y\sinh(\pi y)/\pi\); every denominator distance uses an integer between \(n\) and \(3n\), so its normalized modulus is at least \((1+y^2/(9n^2))^{1/2}\).

Squaring (14) and multiplying by (13) cancels the growing hyperbolic sine. Since \(y^2\operatorname{csch}^2(\pi y)\le1/\pi^2\),
\[
 |\mathscr L_n|\le\frac{\pi a_n^2}{6}
 \int_{-\infty}^{\infty}
 \left(2y^2+\frac3{\pi^2}\right)
 \left(1+\frac{y^2}{9n^2}\right)^{-2n-1}dy.
\tag{15}
\]
This is already a finite explicit bound for every \(n\ge1\). Evaluating the two beta integrals gives
\[
 |\mathscr L_n|\le\frac{\pi^{3/2}a_n^2}{6}
 \left[
 27n^3\frac{\Gamma(2n-1/2)}{\Gamma(2n+1)}
 +\frac{9n}{\pi^2}\frac{\Gamma(2n+1/2)}{\Gamma(2n+1)}
 \right].
\tag{16}
\]
The bracket is \(O(n^{3/2})\). Stirling's formula gives
\[
 a_n\sim\frac{2\sqrt\pi}{\sqrt3}
 n^{-5/2}\left(\frac4{27}\right)^n,
\]
which proves (11).

This estimate establishes the analytic repair. It does not assert that \(16/729\) is the sharp rate. It also does not imply a sign: the forms actually change sign in the numerical examples below.

## 4. A new-prime obstruction for the rational part

The new family retains an explicit lcm bound, but it does not retain the old \(D_n\) prime support. This can be proved for the *actual rational number* \(B_{n,m}\), not just for its sufficient normalizer.

Let \(s=2m+1\), and suppose \(p\) is an odd prime satisfying
\[
 p>n,\qquad p\le mn<2p.
\]
Write \(r=mn-p\), and use the coefficients \(a_k\) from (1). Lucas's theorem gives
\[
 a_{p+u}\equiv-a_u\pmod p\quad(0\le u\le r),
\]
and every coefficient with \(r<k<p\) vanishes modulo \(p\). These conclusions remain true if the second binomial coefficient has a carry: both corresponding low and high coefficients then vanish.

Multiply (7) by \(p^s\). The diagonal terms contribute
\((s-1)\sum_{u=0}^r a_{p+u}^2\). For off-diagonal terms, only \(\ell-k=p\) can survive modulo \(p\): other denominator differences are units, leaving an extra factor of \(p\). Each surviving harmonic difference contains exactly the term \(1/p^{s-1}\). Therefore
\[
 \boxed{p^sB_{n,m}\equiv(s+1)\sum_{u=0}^r a_u^2\pmod p.}
\tag{17}
\]
All terms have valuation at least \(-s\), so a nonzero right side proves exact valuation \(-s\).

For \(m=2\), choose a prime \(p=2n-1\ge5\). Then \(r=1\), \(n\equiv1/2\pmod p\), and
\[
 a_0\equiv1,\qquad a_1\equiv-9/4\pmod p.
\]
Thus
\[
 \boxed{p^5B_{n,2}\equiv\frac{291}{8}\pmod p.}
\tag{18}
\]
For every such prime other than 97,
\[
 \boxed{v_p(B_{n,2})=-5.}
\tag{19}
\]
Since \(p>n\), **no power of \(D_n\)** clears this rational part. There are infinitely many such \(n\), because every odd prime \(p\ge5\), \(p\ne97\), gives \(n=(p+1)/2\).

More generally, at a prime \(p=mn-1>n\), (17) becomes
\[
 p^{2m+1}B_{n,m}\equiv
 \frac{2(m+1)\bigl(m^{2m}+(m+1)^{2m}\bigr)}{m^{2m}}
 \pmod p.
\tag{20}
\]
This identifies a finite exceptional set for each fixed \(m\); no density assertion about primes in this progression is required for the conditional congruence.

One isolated prime factor per index does **not** prove an exponential lower bound on the whole denominator. Accordingly, (19) is not misrepresented as an impossibility theorem for this family.

## 5. Infinitely many forms are nonzero

There is a useful nonvanishing argument despite the sign changes. If \(\zeta(5)\) is irrational, every \(A_{n,2}\zeta(5)-B_{n,2}\) is nonzero because \(A_{n,2}>0\). If it is rational, say \(a/b\), choose infinitely many primes \(p=2n-1\ge5\) with \(p\nmid b\) and \(p\ne97\). Then \(A_{n,2}\zeta(5)\) is \(p\)-integral, while \(B_{n,2}\) has valuation \(-5\), so equality is impossible.

Thus, unconditionally, the sequence in (11) is **infinitely often nonzero**. Under a hypothetical rationality assumption, the argument directly supplies infinitely many nonzero terms. What is missing for an irrationality contradiction is a sufficiently cheap denominator normalization, not an analytic decay estimate or an eventual-sign assertion.

## 6. Evidence and the remaining arithmetic question

The exact script computes \(A\) and \(B\) first and then evaluates the real form at two high precisions. Representative values are:

| \(n\) | \(\mathscr L_{n,2}\), approximately |
|---:|---:|
| 1 | \(3.3090037280\cdot10^{-3}\) |
| 2 | \(3.7687198207\cdot10^{-5}\) |
| 3 | \(-9.1006735188\cdot10^{-8}\) |
| 5 | \(6.4836064368\cdot10^{-12}\) |
| 10 | \(-2.1481583716\cdot10^{-22}\) |
| 20 | \(-6.2878008509\cdot10^{-43}\) |

The rigorous available clearing factor is \(D_{2n}^5\), whose logarithm is \(10n+o(n)\). Its growth exceeds the decay supplied by (11), whose logarithmic rate is \(\log(729/16)\approx3.819\). That comparison shows the **proved upper bounds do not complete the argument**. It does not establish the actual primitive-denominator growth or exclude a sharper analytic/arithmetic treatment.

Next useful tasks are to classify the sum-of-squares residue in (17) over entire prime bands, and to understand the oscillatory sharp asymptotic of (12). For the general \(m\)-fold family, (1)–(9) and (17)–(20) are proved, but the particular exponential estimate in §3 is only established for \(m=2\). Its hyperbolic cancellation is especially favorable there; blindly repeating the same bound for larger \(m\) would not give an integrable majorant.

Reproduce the exact Mellin orthogonality, arithmetic residues and high-precision values with `python codex/computations/check_multiplicity_lifted_zeta_family.py`.
