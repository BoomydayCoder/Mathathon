# Independent review of adaptive total-boundary cancellation

**Status: no proof gap found** in `adaptive-boundary-cancellation.md`. The exact tuning, eventual existence and boundedness, contour comparison, sharp asymptotic, and stated sufficient arithmetic normalization check out. The numerical arithmetic samples remain bounded evidence; they do not establish an all-index denominator obstruction.

The sharp conclusion depends on the already proved and independently reviewed unweighted tail asymptotic in `truncated-derivative.md`. This review checks the new comparison and its normalization, rather than treating that earlier theorem as numerical evidence.

## 1. Exact cancellation and the adjacent polynomial partner

At a positive integer root \(j\) of
\[
 P_n(t)=\frac{\prod_{r=1}^n(t-r)}{\prod_{k=0}^n(t+k)},
\]
put \(u_j=P_n'(j)\) and \(\delta_j=H_{n-j}+H_{n+j}-2H_{j-1}\). The identities
\[
 P_n''(j)=-2u_j\delta_j,\quad
 R_n''(j)=2u_j^2,\quad R_n'''(j)=-12u_j^2\delta_j
\]
are correct. With \(b=n+1\) and \(w=(t-a)/(t+b)\), the product derivative at the root is
\[
 (R_nw)'''(j)=R_n'''(j)w(j)+3R_n''(j)w'(j).
\]
The omitted product-rule terms vanish because \(R_n(j)=R_n'(j)=0\). Since \(w'= (a+b)/(t+b)^2\), summation yields exactly
\[
 \sum_{j=1}^n(R_nw)'''(j)=-12D+6(a+b)(2E+G).
\]
Thus \(a_n=2D/W-b\), \(W=2E+G\), is the unique tuning whenever \(W\ne0\). The cancellation concerns the **sum** of the root terms, not each individual root term; the proof never silently strengthens it to the latter.

The identity
\[
 P_nw=(1-d_n)P_n+d_nP_{n+1},\quad d_n=(a_n+b)/(2b),
\]
is algebraically exact. Because the Mellin transforms of the adjacent Legendre polynomials carry opposite signs, the partner must be
\(Q_n=(1-d_n)L_n-d_nL_{n+1}\), as stated. Multiplying the two Mellin transforms gives \(R_nw\) with no extra sign. The monomial kernel then produces a pure \(\zeta(5)\) form with rational coefficients.

At \(n=1\), the exact values \(a_1=7/10\), \(d_1=27/40\), and
\(20Q_1=-7+68x-81x^2\) check. The displayed coefficient \(-572\) of \(\zeta(5)\) is also correct.

## 2. Uniform control of every root weight

The claim that the first root dominates needs a uniform bound across \(2\le j\le n\), including the opposite endpoint \(j=n\). Such a bound is available.

Let
\[
 g_j=\binom nj\binom{n+j}j=\frac{(n+j)!}{j!^2(n-j)!}.
\]
Then
\[
 \frac{g_{j+1}}{g_j}=
 \frac{(n+j+1)(n-j)}{(j+1)^2}
 =\frac{n(n+1)-j(j+1)}{(j+1)^2}
\]
is strictly decreasing with \(j\). Hence the sequence first increases and then decreases, and its minimum on \([2,n]\) is at an endpoint.

For \(n\ge5\),
\[
 g_n=\binom{2n}n\ge g_2
 =\frac{n(n-1)(n+1)(n+2)}4\ge\frac{n^4}4.
\]
The first inequality holds at \(n=5\), and persists by induction: the successive central-binomial ratio is \(2(2n+1)/(n+1)\), whereas the successive ratio of \(g_2\) is \((n+3)/(n-1)\), which is smaller for \(n\ge5\).

Because \(|u_j|=1/(j^2g_j)\),
\[
 \sum_{j=2}^n u_j^2
 \le\frac{16}{n^8}\sum_{j=2}^\infty j^{-4}
 \le\frac{5}{3n^8}.
\tag{1}
\]
Here \(\sum_{j=2}^\infty j^{-4}\le1/16+\int_2^\infty x^{-4}dx=5/48\). Also \(|\delta_j|\le4H_{2n}\), so
\[
 \sum_{j=2}^nu_j^2|\delta_j|
 \le\frac{20H_{2n}}{3n^8}=O(n^{-8}\log n).
\tag{2}
\]
This explicitly justifies the uniform root estimates used in the draft. Merely bounding every term and then multiplying by \(n\) would lose a power; the convergent \(j^{-4}\) sum is the reason that loss does not occur.

At \(j=1\), \(u_1^2=1/[n^2(n+1)^2]\) and \(\delta_1=2\log n+O(1)\). Since \(j+b\ge n+3\) for \(j\ge2\), (1)–(2) give
\[
 D=u_1^2\delta_1+O(n^{-8}\log n),
\]
\[
 W=u_1^2\left(\frac{2\delta_1}{n+2}+\frac1{(n+2)^2}\right)
 +O(n^{-9}\log n).
\]
Their main terms are respectively of orders \(n^{-4}\log n\) and \(n^{-5}\log n\), so both relative errors are \(O(n^{-4})\). Thus their ratio has absolute error \(O(n^{-3})\), even after the subtraction of \(n+1\):
\[
 a_n=1-\frac{n+2}{2\delta_1(n+2)+1}+O(n^{-3})
 =1-\frac1{4\log n}+O((\log n)^{-2}).
\]
This proves eventual \(D>0,W>0\), \(0<a_n<1\), and \(d_n\to1/2\). No uniform all-small-index assertion is required.

## 3. Contour orientation, integration by parts, and the multiplier

For \(c=n+1/2\), the added pole at \(-n-1\) lies strictly outside the right contour. The upward vertical line is closed clockwise to the right, so the sum of residues carries a minus sign. With \(K(z)=\pi\cot\pi z\),
\[
 -\frac16\sum_{t>n}(R_nw)'''(t)
 =-\frac1{12\pi i}\int_{c-i\infty}^{c+i\infty}K'''(z)R_n(z)w(z)\,dz.
\]
The integration-by-parts boundary terms vanish by rational decay before applying reflection. Half-integer right edges and horizontal edges tending to infinity justify the contour closure. No assumptions about a sign of the summands are needed.

The exact derivative is
\[
 K'''=4\pi^4\csc^2\pi z-6\pi^4\csc^4\pi z.
\]
With \(R_n=\sin^2(\pi z)U_n^2/\pi^2\) and \(\sin^2(\pi(c+iy))=\cosh^2(\pi y)\), this gives
\[
 T_n(w)=-\frac\pi{12}\int_{\mathbb R}
 (4-6\operatorname{sech}^2\pi y)U_n(c+iy)^2w(c+iy)\,dy.
\tag{3}
\]
Both the minus sign and the factor \(\pi/12\) are correct.

## 4. Uniform beta bound and the sharp comparison constant

The factorization
\[
 U_n(z)=B(z,n+1-z)B(z,n+1)
\]
is exact. On \(\Re z=c=n+1/2\), the second beta integral satisfies
\(|B(c+iy,n+1)|\le B(c,n+1)\). For the first factor,
\[
 |\Gamma(c+iy)|\le\Gamma(c),\qquad
 |\Gamma(1/2-iy)|^2=\pi\operatorname{sech}(\pi y).
\]
Consequently
\[
 |U_n(c+iy)|^2\le U_n(c)^2\operatorname{sech}(\pi y).
\tag{4}
\]
The bound is uniform in \(n\) and \(y\); it is not an asymptotic valid only on bounded vertical segments.

For the constant, duplication gives
\[
 U_n(c)=\frac{\pi}{2^{2n+1/2}}
 \frac{\Gamma(n+1/2)^2}{\Gamma(n+3/4)\Gamma(n+5/4)}
 \sim\frac\pi{\sqrt2\,n}4^{-n}.
\]
Thus \(U_n(c)^2\sim\pi^2 16^{-n}/(2n^2)\), exactly as stated.

If \(|a_n|\le A\), then \(w\) is uniformly bounded on the line, and
\[
 w(c+iy)-\frac12
 =\frac{-1/2-2a_n+iy}{2(c+b+iy)}.
\]
Since \(c+b\ge2n\), its modulus is bounded by a constant depending only on \(A\), times \((1+|y|)/n\). The first moment of \(\operatorname{sech}(\pi y)\) is finite. Combining this with (3)–(4) proves the stated uniform estimate
\[
 T_n(w)=\frac12T_n(1)+O(n^{-3}16^{-n}).
\tag{5}
\]
The first-root analysis supplies the needed boundedness of the actual tuning \(a_n\), so there is no circular dependence.

The earlier reviewed tail theorem is
\[
 E_n=-\frac1{24}\sum_{t>n}R_n'''(t)
 \sim\frac\pi2\frac{\log n}{n^3}16^{-n}.
\]
Since \(T_n(1)=4E_n\), (5) yields
\[
 \boxed{\mathcal L_n=2E_n+O(n^{-3}16^{-n})
 \sim\pi\frac{\log n}{n^3}16^{-n}.}
\]
The error lacks the logarithm, so it is lower order. This establishes eventual positivity and the exact exponential rate, not just an upper bound.

## 5. The arithmetic statements have the advertised limited scope

The claim \(D_nD_{2n}u_j\in\mathbf Z\) can be proved directly at each prime power. With \(q=p^a\), the factorial-floor contribution to \(v_p(u_j)\) is
\[
 2\left\lfloor\frac{j-1}q\right\rfloor
 +\left\lfloor\frac{n-j}q\right\rfloor
 -\left\lfloor\frac{n+j}q\right\rfloor.
\]
For \(q\le n\), this is at least \(-2\): write \(j-1=qa+r\), \(n-j=qb+s\), so the loss is \(\lfloor(2r+s+2)/q\rfloor\le2\). For \(n<q\le2n\), the loss is at most one. These are precisely the valuations supplied by \(D_nD_{2n}\).

Squaring, then accounting for \(\delta_j\)'s harmonic denominator and \(j+b\le2n+1\), proves that
\[
 K_n=D_n^2D_{2n+1}^4
\]
clears both \(D\) and \(W\). Therefore \((n+1)|K_nW|\) is a sufficient integer multiplier for \(Q_n\), whenever \(W\ne0\). This bound includes the cost of the adaptive parameter.

The reported primitive linear-form normalization is also correct. If the partner has first been scaled so that the zeta coefficient is \(A\in\mathbf Z\), and the rational part is reduced as \(p/q\), multiplying by \(q\) gives coefficients \(qA,p\). Since \(\gcd(p,q)=1\), their gcd is \(\gcd(A,p)\), exactly the additional cancellation used by the script.

The new prime factors in the first few tuning denominators are actual finite examples. They neither prove infinitely many new primes nor give an asymptotic lower bound. The draft correctly separates this evidence from its proved analytic success.

## 6. Verification

The existing exact checker was rerun successfully. The independent checker `check_adaptive_boundary_independent_review.py` additionally verifies the uniform endpoint comparison, root-weight arithmetic, and the rational identity for \(w-1/2\). These finite checks supplement the proofs above.

No change to the construction is required. The open issue is the primitive arithmetic cost of the tuning, not cancellation, existence, real decay, or eventual nonvanishing.
