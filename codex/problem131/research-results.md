# Research results: constructions, obstructions, and unexpected links

This is a compact record of independently derived results from the ongoing investigation. “Proved” means that a complete argument is written in the linked note; review status is stated where relevant. It does **not** mean historically novel. These results do not constitute a new irrationality proof for a higher odd zeta value.

The main exposition is [the discovery guide](discovery-guide.md). The results here test its mechanism by changing the derivative, the summation range, the arithmetic normalization, and the representation of the auxiliary function.

## 1. One notation for the derivative constructions

Let

\[
P_n(t)=\frac{\prod_{j=1}^n(t-j)}{\prod_{j=0}^n(t+j)},
\qquad R_n(t)=P_n(t)^2,
\qquad D_n=\operatorname{lcm}(1,\ldots,n),\quad D_0=1.
\]

For an integer \(m\ge1\) and \(0\le k\le n\), define

\[
E_{n,k,m}=-\frac1{(2m)!}\sum_{t=k+1}^{\infty}R_n^{(2m-1)}(t).
\tag{1}
\]

Partial fractions show

\[
E_{n,k,m}=A_n\zeta(2m+1)-\beta_{n,k,m},
\qquad
A_n=\sum_{j=0}^n\binom nj^2\binom{n+j}{j}^2,
\tag{2}
\]

where \(\beta_{n,k,m}\) is rational. Thus a single desired odd zeta value is already isolated. The real question is whether a nonzero value in (2) remains small after clearing its rational denominator.

Write

\[
F_n=-\sum_{t=1}^{\infty}R_n'(t),
\qquad \rho=(\sqrt2-1)^4.
\]

Then \(E_{n,k,1}=F_n/2\) for every \(k\le n\), because first derivatives vanish at the initial double zeros. This cutoff independence is special to the first derivative.

## 2. A proved failure of the unchanged higher-derivative construction

**Theorem.** For every fixed integer \(m\ge2\),

\[
E_{n,0,m}\sim
\frac{(4\log n)^{2m-3}}
 {(2m)(2m-3)!\,n^4}.
\tag{3}
\]

In particular the full higher-derivative sum is eventually positive, but has polynomial decay rather than exponential decay.

The source of the leading term is exactly the former zero at \(t=1\). Taking one derivative preserves that zero; taking three or more detects its Taylor coefficients. The other initial zeros contribute a smaller error, and a uniform derivative-tail estimate controls the infinite remainder. This is an asymptotic theorem for the **whole signed sum**, not an inference from one nonzero summand. Its proof and independent review are in [the logarithmic-kernel analysis](log-kernel-extension.md) and [the review](odd-zeta-independent-review.md).

There is also a shorter recurrence obstruction. For

\[
(Tu)_n=(n+1)^3u_{n+1}
 -(2n+1)(17n^2+17n+5)u_n+n^3u_{n-1},
\]

Zudilin’s certificate gives

\[
(TE_{\cdot,0,2})_n
=-\frac{(2n+1)\{3+8n(n+1)(H_{n-1}+H_{n+1})\}}
 {n^2(n+1)^2}.
\tag{4}
\]

The right side is not exponentially small. Applying \(T\) to an exponentially small sequence would give an exponentially small sequence times a cubic factor, so (4) already rules out that behavior. See [the forcing theorem](forced-recurrence-obstruction.md).

**A stronger correction obstruction.** No expression

\[
Q_n=\sum_{j=0}^{d}a_j(n)H_n^j,
\qquad a_j\in\mathbb Q(n),\quad d<\infty,
\]

can make \(E_{n,0,2}-Q_n\) exponentially small. The proof uses extreme-pole propagation for a rational-function difference equation, then the independence of powers of \(H_n\) over rational functions of \(n\). This includes subtracting the contributions of any **fixed** number of initial zeros. It excludes neither growing-length sums nor other kinds of correction. The precise proof is in [forcing-correction.md](forcing-correction.md).

The design lesson is therefore specific: algebraic isolation of \(\zeta(5)\) succeeds, but the boundary cancellation responsible for the successful \(\zeta(3)\) recurrence does not survive.

## 3. Full truncation repairs analysis and exposes arithmetic

**Theorem.** For every fixed integer \(m\ge2\),

\[
E_{n,n,m}\sim
\frac{\pi(2\log n)^{2m-3}}
 {(2m)(2m-3)!\,n^3}\,16^{-n}.
\tag{5}
\]

Thus deleting all initial zero contributions restores exponential decay and eventual positivity. The rate \(1/16\) is exact, and does not improve as the fixed derivative order increases.

This is proved through a forced recurrence with a moving lower boundary. All its boundary terms factor through one rational function whose constant and linear Taylor coefficients vanish. Its first surviving quadratic coefficient explains the logarithmic exponent \(2m-3\). A discrete maximum-principle comparison turns the formal leading balance into a genuine asymptotic. The argument was independently reproduced and reviewed; see [the uniform theorem](higher-truncated-tail.md) and [the truncation review](truncation-review.md).

The arithmetic also has a proved statement:

\[
(2m)D_n^2D_{2n}^{,2m-1}\,\beta_{n,n,m}\in\mathbb Z.
\tag{6}
\]

Moreover, for every odd prime \(p>2m\), \(p\ne17\), with \(n=(p+1)/2\),

\[
v_p(\beta_{n,n,m})=-(2m-1).
\tag{7}
\]

These primes lie above \(n\), so the original \(D_n\)-based denominator support genuinely changes. The exceptional \(17\) comes from the squared residue ratio \(1+1/16\), and persists at every fixed derivative order.

Equation (6) is an **upper bound on a sufficient denominator multiplier**. Equation (7) is an **exact valuation on a specified prime family**. Neither establishes the full asymptotic reduced-denominator growth. In particular, the failure of the upper-bound decay budget is not an impossibility theorem for every sharper normalization.

## 4. A saddle cut reconnects every odd derivative to the original \(\zeta(3)\) form

**Theorem.** Fix \(m\ge1\). If integers \(k_n\) satisfy
\(k_n=n/\sqrt2+O(1)\), then

\[
\frac{E_{n,k_n,m}}{F_n}
 =\frac{(-1)^{m+1}(2\pi)^{2m-2}}{(2m)!}
   +O_m(n^{-1/2}).
\tag{8}
\]

Thus every fixed-order saddle-cut family has exact exponential rate \(\rho\), and its eventual sign is \((-1)^{m+1}\). The first constants are

\[
\frac12,\qquad-\zeta(2),\qquad2\zeta(4),
\qquad-\frac32\zeta(6).
\]

This is stronger than observing matching numerical decay: the leading error is an explicit even-zeta multiple of the original \(\zeta(3)\) linear form. The all-fixed-\(m\) identity and its normalization were independently checked; see [partial-truncation-analysis.md](partial-truncation-analysis.md), Section 8, and [the review](odd-zeta-independent-review.md), Section 10. No uniformity for \(m\) growing with \(n\) is asserted.

For \(m=2\), the exact identity behind the theorem is particularly transparent. Put \(c=k+1/2\) and

\[
U_n(z)=\frac{\Gamma(z)^2\Gamma(n+1-z)}{\Gamma(n+1+z)}.
\]

Then

\[
E_{n,k,2}=-\zeta(2)F_n+
 \frac\pi8\int_{\mathbb R}\operatorname{sech}^2(\pi y)
                  U_n(c+iy)^2\,dy,
\tag{9}
\]

and

\[
\bigl|E_{n,k,2}+\zeta(2)F_n\bigr|\le\frac14U_n(c)^2.
\tag{10}
\]

The complex square in (9) is essential; it is not an absolute square. The bound follows from the beta factorization
\(U_n(z)=B(z,z)B(n+1-z,2z)\), valid in the required strip. Gamma reflection and cotangent contour summation explain the even powers of \(\pi\) in (8).

The saddle is also visible in exact factorial arithmetic:

\[
U_n(j)^2=
\frac1{j^4\binom nj^2\binom{n+j}{j}^2}.
\tag{11}
\]

The small correction envelope is reciprocal to the large binomial summands of \(A_n\). This is why the same \(n/\sqrt2\) saddle appears in both places.

A central cut is not automatically the cut of smallest absolute error. For the third derivative, exact harmonic monotonicity gives one discrete minimum and at most two sign crossings. These notions must remain distinct; [the shape analysis](partial-truncation-structure.md) treats them separately.

## 5. Movable cuts admit an exact prime-by-prime test

For the third derivative, define the rational correction

\[
C_{n,k}=-\frac1{24}\sum_{j=1}^{k}R_n'''(j),
\qquad \beta_{n,k,2}=\beta_{n,0,2}+C_{n,k}.
\]

A safe general bound is

\[
4D_n^2D_{n+k}^{,3}\beta_{n,k,2}\in\mathbb Z.
\tag{12}
\]

More precisely, if \(p>3\) is prime and \(n<p\le n+k\), put \(a=p-n\). Then

\[
p^3C_{n,k}\equiv
 \frac12\sum_{j=a}^{k}P_{a-1}(j)^2\pmod p.
\tag{13}
\]

All denominators in the scaled congruence are prime to \(p\). The corrected coefficient has valuation exactly \(-3\) if and only if the finite-field sum is nonzero. At a newly crossed prime boundary \(p=n+k\), the sum contains one nonzero square, so cancellation is impossible at that step. Later terms can cancel it.

For the full cut, the finite-field sum becomes a rational constant \(W_{p-n-1}\), independent of the large degree \(n\) except through that complementary index. The first values are

\[
W_0=0,\qquad W_1=-13,\qquad
W_2=-\frac{905}{4},\qquad W_3=-\frac{167965}{36}.
\]

These constants also satisfy the real-lattice identity

\[
2A_m\zeta(2)+W_m
 =\sum_{t\in\mathbb Z\setminus\{-m,\ldots,0\}}R_m(t)>0.
\tag{14}
\]

This links denominator valuations, finite-field sums, and a positive \(\zeta(2)\) linear form built from the same rational function. Proofs are in [prime-band-arithmetic.md](prime-band-arithmetic.md) and [partial-truncation-arithmetic.md](partial-truncation-arithmetic.md).

> **Unresolved arithmetic theorem.** To obstruct the full-cut small-integer argument asymptotically using only the cubic primes, one would need a surviving-prime logarithmic mass exceeding \((\log16/3+\varepsilon)n\). The exact criterion does not yet prove this. Fixed-distance prime families and finite high survival counts do not establish the required weighted-density bound.

## 6. Two exact links to other parts of mathematics

**Inverse Hilbert matrices.** Let \(H_N=(1/(i+j+1))_{0\le i,j<N}\), and let \(c_j\) be the coefficient vector of the shifted Legendre polynomial, padded to length \(N\). Then

\[
H_N^{-1}=\sum_{j=0}^{N-1}(2j+1)c_jc_j^{\mathsf T},
\qquad
\operatorname{tr}(H_N^{-1})
 =\sum_{j=0}^{N-1}(2j+1)A_j.
\tag{15}
\]

This follows directly by changing from monomials to the orthogonal Legendre basis. It proves inverse integrality and explains why the same exponential base \((1+\sqrt2)^4\) controls Apéry coefficient growth and Hilbert-matrix ill-conditioning. See [the exact matrix theorem](hilbert-matrix-connection.md).

**Rational orthogonality.** With

\[
M_n(t)=\int_0^1x^{t-1}\mathcal P_n(x)\,dx,
\]

the shifted Mellin transform is exactly

\[
M_n(z+\tfrac12)=
\frac{(-1)^n}{z+n+\tfrac12}
 \prod_{j=1}^{n}\frac{z-(j-\tfrac12)}{z+(j-\tfrac12)}.
\tag{16}
\]

The factors have modulus one on the imaginary axis. After multiplication by \(\sqrt{2n+1}\), they form an orthonormal rational system for the boundary inner product

\[
\langle f,g\rangle=
\frac1{2\pi}\int_{\mathbb R}f(iy)\overline{g(iy)}\,dy.
\]

The substitution \(x=e^{-u}\), with the half-density \(e^{-u/2}\), and Fourier Plancherel identify this rational orthogonality exactly with Legendre orthogonality on \((0,1)\). Thus the auxiliary product is simultaneously an interpolation object, a Mellin transform, and a reflected-pole orthogonal function. The precise identities and an explicit rational Gram-matrix comparison are in [rational-orthogonality.md](rational-orthogonality.md).

These connections explain the construction and suggest alternative estimates. They do not automatically supply integer coefficients or useful denominator bounds.

## 7. Distinguish three kinds of negative evidence

The project has examples of genuinely different logical strengths:

- **Proved asymptotic failure:** the unchanged higher-derivative forms have the polynomial asymptotic (3). No hidden cancellation can turn that particular family exponentially small.
- **Proved finite failure after exact normalization:** for the examined two-dimensional positive weighted family, exact primitive reduction and positive rational partial sums prove \(|U_n\zeta(5)-V_n|>1\) for \(n=1,\ldots,25\). This does not prove failure for every later index. See [weighted-family-arithmetic.md](weighted-family-arithmetic.md).
- **Insufficient bounds:** an available denominator upper bound may exceed the proved decay budget. That says the estimates do not close the argument; it does not prove the actual reduced denominator is that large.

A broader positive-weight space has already improved some finite examples, so the second and third categories must not be promoted to impossibility claims. Conversely, a single primitive error below one proves no irrationality theorem without an infinite nonvanishing sequence tending to zero.

## 8. Concrete next lemmas

The next useful work has explicit mathematical targets:

1. Derive a recurrence or modular description for \(W_m\) strong enough to control the changing divisibility diagonal \(p=n+m+1\). This would upgrade the exact local criterion toward an arithmetic rate theorem.
2. Prove reduced-denominator bounds at the saddle cuts, rather than relying only on the sufficient multiplier (12). The analytic rate is already known sharply enough to specify the needed arithmetic advantage.
3. Quantify the rational Gram comparison for the positive odd-zeta moment matrices. The proposed estimate separates an exact rational determinant from a one-variable weighted maximum; its cancellation on the relevant large-parameter scale remains unproved.
4. Find a growing-dimensional positive auxiliary space whose primitive arithmetic error tends to zero. Fixed-dimensional coefficient tuning can be arithmetically equivalent to rescaling, and fixed positive cones have their own lower-gap obstructions.

The original [Zudilin paper](https://arxiv.org/abs/math/0202159) remains the primary starting point. Current claims concerning Fauzan’s 2026 \(\zeta(5)\) work are documented separately in [source provenance](../sources/odd-zeta-provenance.md); this project has not independently certified that complete preprint or locally built its Lean formalizations. No claim here depends on treating an unreviewed research proposal as an established theorem.
