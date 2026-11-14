# Research results: constructions, obstructions, and unexpected links

This is a compact record of independently derived results from the ongoing investigation. “Proved” means that a complete argument is written in the linked note; review status is stated where relevant. It does **not** mean historically novel. The completed zeta-five assembly is a reviewed computer-assisted alternative derivation retaining Fauzan’s construction, not a newly invented irrationality construction or evidence of external consensus. No all-odd irrationality proof has resulted.

The main exposition is [the discovery guide](discovery-guide.md). The results here test its mechanism by changing the derivative, the summation range, the arithmetic normalization, and the representation of the auxiliary function.

**Current status, 14 November 2026.** The [PNT-free zeta-five assembly](pnt-free-zeta5-guide.md) now keeps degree \(h\): independently reviewed floor rounding or the primitive-part lemma removes the unnecessary determinant power. The adaptive boundary repair, signed zeta-seven model and its total-variation obstruction, and Euler-depleted functional have completed independent reviews; Sections 21–24 record their exact limits. Higher-Rodrigues geometry and Eulerian contour developments have their separate status in Section 25. The pinned existing Lean formalization also completed its local build and 71-module project replay, with an explicit PNT axiom; it does not formalize these new constructions.

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

The original [Zudilin paper](https://arxiv.org/abs/math/0202159) remains the primary starting point. Current claims concerning Fauzan’s 2026 \(\zeta(5)\) work are documented separately in [source provenance](../sources/odd-zeta-provenance.md). The [assembled alternative derivation](fauzan-simplified-proof-audit.md) passed end-to-end mathematical review. Separately, the [pinned existing Lean formalization](lean-local-verification.md) passed its local build, ordinary-series statement/axiom checks, and replay of all 71 project modules. Its external mathematical axiom is the usual PNT; replay was against imported environments, not a fresh replay of all Mathlib constants. This does not Lean-verify the new equilibrium construction or PNT-free certificate. No claim here depends on treating an unreviewed research proposal as an established theorem.

## 9. An analytic equilibrium construction replaces Fauzan’s sixteen-component measure

This result improves a specific analytic component of Fauzan’s preprint. Its analytic proof alone does not validate the arithmetic; those components have their own reviews.

For the source parameters \(\alpha=3/40\), \(\lambda=37/40\), the external field has derivative

\[
V'(x)=\int_0^\infty\frac{q(u)}{x+u^2}\,du,
\qquad
q(u)=\begin{cases}
-3,&0<u<\alpha,\\
3,&\alpha<u<1,\\
2,&u>1.
\end{cases}
\]

This single negative-to-positive sign change makes it possible to construct the actual equilibrium measure on one interval \([a,b]\). Sweeping each charge onto that interval gives the density

\[
d\rho(x)=-\frac{\sqrt{(x-a)(b-x)}}{2\pi}
 \left[\int_0^\infty
 \frac{q(u)\,du}{(x+u^2)\sqrt{(u^2+a)(u^2+b)}}\right]dx.
\tag{17}
\]

Two endpoint equations fix the mass \(\lambda\) and make the density vanish at both edges. Subtracting the integrand’s value at the charge-sign-change point proves that the bracket is negative. The same sign proves that \(2U_\rho-V\) increases before the support, is constant on it, and decreases afterwards. This establishes positivity and the **global** potential inequality analytically; it does not rely on a grid of potential values.

An exact integer-interval certificate locates the unique endpoint pair in

\[
0.00007119508<a<0.00007119517,
\qquad
0.75394957<b<0.75394967.
\]

Two validated energy evaluations then prove

\[
\mathcal J(\rho):=\lambda\sup_{x\ge0}(2U_\rho(x)-V(x))
 -I(\rho)+C_*< -\frac{69}{50}=-1.38.
\tag{18}
\]

One evaluates the density through Carlson integrals; the other uses an equilibrium virial identity and elementary off-support potentials. They use different integral representations but share the certified endpoint box and Arb arithmetic library. The respective enclosures are approximately \(-1.3855\pm0.0000964\) and \(-1.3855\pm0.0000295\). The mathematical identities and certificate reductions were independently checked within the session; no historical novelty is asserted.

**Conditional quantitative improvement.** The preprint’s arithmetic constant is

\[
A_{200}=\frac{127125602969131786927559}
 {94195881588024216000000},
\qquad
1600\left(\frac{69}{50}-A_{200}\right)>48.
\]

Retaining the preprint’s separate moment, integrality, and arithmetic-normalization arguments, this stronger analytic comparison implies for its existing polynomial sequence

\[
0<Q_{40n,200}(\zeta(5))<e^{-48n^2}
\]

for all sufficiently large \(n\), improving the source’s stated exponent \(139/5=27.8\). This concerns an eventual rate, not a certified starting index or unchanged finite-error constants. This quantitative comparison retains the source’s original refined normalizer. The completed assembled alternative below uses a different plain-rank normalizer and the weaker but sufficient eventual exponent 4.

The [equilibrium guide](fauzan-equilibrium-guide.md) explains the construction and exact scope. The [endpoint proof](equilibrium-endpoint-certificate.md), [Carlson energy certificate](../computations/certify_equilibrium_energy.py), and [virial energy certificate](../computations/certify_equilibrium_energy_elementary.py) make the local analytic claim reproducible. The main gain is that the sign of the field’s charge distribution replaces sixteen trial components and many separate potential comparisons.

## 10. A common fixed-odd analytic budget, with an explicit arithmetic gate

The equilibrium improvement extends analytically to each **fixed** odd \(s\ge3\). For the positive moment determinant \(\Delta_{K,s}\) with the source proportions \(K=40n\), \(N=3n\), \(h=37n\), and its fixed factorial normalizer \(S_K\), the [reviewed theorem](common-odd-zeta-analytic-budget.md) proves

\[
\deg\Delta_{K,s}=h,\qquad \Delta_{K,s}(\zeta(s))>0,
\qquad
\limsup K^{-2}\log[S_K\Delta_{K,s}(\zeta(s))]\le-\frac{69}{50}.
\]

The change of odd index affects polynomial weight factors and hence only \(O_s(K\log K)\) terms, leaving the leading external field unchanged. Uniform control at the origin and infinity is included in [the independent review](common-odd-zeta-budget-review.md).

This is not an all-odd irrationality theorem. The missing condition is a positive rational normalization \(\nu_{K,s}\) with

\[
\nu_{K,s}S_K\Delta_{K,s}(X)\in\mathbb Z[X],
\qquad \limsup K^{-2}\log\nu_{K,s}<69/50.
\]

The exact analytic index shift uses
\(\mathcal L_s=(2\theta+s+1)(2\theta+s+2)/[s(s+1)]\), \(\theta=t\,d/dt\). It turns a simple pole into a combination including a nonzero cubic pole. Moreover, no universal finite-order rational differential operator both preserves the entire prescribed simple-pole domain and implements the shift. Actual pole constants lose two additional prime powers in an explicit prime family. These are [proved obstructions to a particular free transfer](odd-functional-shift.md), not impossibility results for determinant-specific arithmetic cancellation.

## 11. The recurrence is a resistor chain, and \(\zeta(3)\) is its sharp boundary constant

The Casoratian has a concrete energy interpretation. Put

\[
c_n=(n+1)^3A_nA_{n+1},\qquad r_n=c_n^{-1}.
\]

Writing a recurrence solution as \(u_n=A_nv_n\) gives exactly Kirchhoff’s law

\[
c_n(v_{n+1}-v_n)=c_{n-1}(v_n-v_{n-1}).
\]

For the companion solution \(B_n\), the current is \(6\). The independently established limit \(B_n/A_n\to\zeta(3)\) therefore identifies

\[
\frac{\zeta(3)}6=\sum_{n\ge0}r_n,
\qquad
A_n\zeta(3)-B_n=6A_n\sum_{j=n}^\infty r_j.
\]

The approximation error is the resistance of the omitted tail. Positive resistance explains its sign; growing conductance explains its smallness.

More strongly, for every finitely supported real sequence \(f\), with \(b_n=(2n+1)(17n^2+17n+5)\),

\[
\sum_{n\ge0}b_nf_n^2-2\sum_{n\ge0}(n+1)^3f_nf_{n+1}
=\sum_{n\ge0}c_n
 \left(\frac{f_n}{A_n}-\frac{f_{n+1}}{A_{n+1}}\right)^2
\ge\frac6{\zeta(3)}f_0^2.
\]

The constant is sharp, approached by explicit rational finite sequences whose voltage drops are proportional to the edge resistances. Positive definiteness also proves positivity of every Schur-complement denominator in the continued fraction, resolving the earlier concern about its nested denominators.

The [independent audit](discrete-green-independent-review.md) also gives the entire finite Green kernel. If \(J_N\) is the tridiagonal matrix of this quadratic form on indices \(0,\ldots,N\), then

\[
(J_N^{-1})_{ij}=A_iA_j\sum_{k=\max(i,j)}^N r_k.
\]

Its entrywise limit is \(A_i(A_j\zeta(3)-B_j)/6\) for \(i\le j\). The difference between this limit and any finite leading block is positive semidefinite and rank one. Thus the growing and decaying Apéry solutions are exactly the two factors of the Green kernel. This identity makes no assumption about a self-adjoint extension of an infinite operator.

The [network derivation](discrete-green-motivation.md) is an exact identity, not merely an analogy. Its finite matrix and minimizer identities pass [exact checks](../computations/check_discrete_green.py). It still requires the separate rational-function argument to identify the total resistance with \(\zeta(3)\), and separate arithmetic to prove irrationality.


## 12. Arithmetic-envelope integration without the source's tables

For \(\alpha=3/40\), \(\lambda=37/40\), the [table-free arithmetic calculation](fauzan-arithmetic-simplification.md) proves, **conditional on the source's local valuation envelope and Proposition 5.2**, that its fixed-cutoff normalizer satisfies

\[
\limsup_{\substack{K\to\infty\\40\mid K}}
 K^{-2}\log m_{K,1000}
 \le\frac{439979}{320000}<\frac{11}{8}.
\]

The construction cancels the main quadratic terms before integration, bounds the residual floor-function fluctuations by a two-value variance estimate, and evaluates the fractional-part averages using Euler's constant. Convexity supplies a sufficient elementary bound for that constant; five positive \(\operatorname{atanh}\)-series terms supply a rational logarithm bound. The remaining outer integral uses positive-part primitives. No 143-piece inner or 11-piece outer integration table is needed.

Combining with the certified analytic budget gives the exact positive margin

\[
\frac{69}{50}-\frac{439979}{320000}
 =\frac{1621}{320000}>\frac1{200}.
\]

The [exact rational script](../computations/check_table_free_arithmetic.py) passes, but it checks the concluding algebra, not the source's local prime-adic lemmas. Those premises and integrality are separate from this calculation; they have now been independently reconstructed in [the local arithmetic review](fauzan-local-arithmetic-review.md), including the small-prime bound. The alternative assembly, with the modified normalizer below, has now [passed end-to-end review](fauzan-simplified-proof-audit.md). The larger cutoff also retains the source's very large threshold \(K\ge200{,}000{,}000\). This is a simpler sufficient asymptotic arithmetic estimate for the source's zeta-five construction, not a proven arithmetic transfer to all odd indices. See [the guide](fauzan-equilibrium-guide.md#11-a-table-free-arithmetic-budget-for-the-sources-zeta-five-construction) for the combined interpretation.


## 13. The rational square is forced within its minimal interpolation ansatz

For \(n\ge1\), a rational function with only poles \(0,-1,\ldots,-n\) of order at most two, decay \(O(t^{-2})\), and double zeros at \(1,\ldots,n\), is necessarily \(CR_n(t)\). The proof is a degree count: the numerator has degree at most \(2n\), already exhausted by its prescribed roots. A nonzero function with these zeros and decay needs at least \(n+1\) double poles. This is a classification under explicit constraints, not a classification of all irrationality constructions.

For any \(n+1\) distinct integer pole positions, the lcm \(L\) of pairwise differences is divisible by \(D_n\), by the pigeonhole principle modulo each prime power. Consecutive positions attain equality. The least universal multiplier for the simple-pole coefficients of a square with arbitrary integral residues is more precisely \(L/\gcd(L,2)\); consecutiveness minimizes that too. Special residue cancellation and the shifted harmonic sums remain additional costs.

The actual reflection is reciprocal, \(R_n(t)R_n(-t)=t^{-4}\), not additive reflection about \(-n/2\). It is equivalent to the finite half-plane inner-function structure of \(tP_n(t)\).

A deformation shows why these properties are not enough. For each fixed integer \(d\ge2\),

\[
-\sum_{t\ge1}\frac{d}{dt}\bigl[d^{-2}R_n(t/d)\bigr]
\sim\frac4{d^3}\left[\frac{\Gamma(1/d)}{\Gamma(1-1/d)}\right]^2
 n^{-4/d}\log n.
\]

Dilation preserves integral residues and zeta isolation but introduces unmatched lower-boundary terms in the telescoper. The asymptotic proves actual polynomial decay, not merely the insufficiency of an exponential estimate. The [classification proof](zudilin-design-classification.md) and [independent review](zudilin-design-independent-review.md) establish all these statements; the [discovery guide](discovery-guide.md) integrates their motivation into the elementary proof.

## 14. Plain rank plus table-free integration still fits the analytic budget

The [independent local arithmetic reconstruction](fauzan-local-arithmetic-review.md) checks inner and outer primes, the Bernoulli distribution mechanism, small-prime control, and coefficientwise integrality. It also permits the simpler outer bound

\[
v_p^{\mathrm G}\det(A+p^{-1}L)\ge2\sum_iw_i-\operatorname{rank}L
\]

for the actual local matrices, whose weights satisfy \(w_i\le0\). This replaces the refined rank/zero-weight correction. Its extra prime-averaged cost is exactly \(9/640\).

Using the fixed cutoff \(M=10000\), the [independently reviewed table-free integration](fauzan-arithmetic-simplification-review.md) gives the relaxed normalization budget

\[
A_{\mathrm{plain}}\le\frac{2203747}{1600000},
\qquad
\frac{69}{50}-A_{\mathrm{plain}}\ge\frac{4253}{1600000}>0.
\]

Dropping the refined correction changes the normalizer. Define \(\widetilde m_{K,M}=\prod_{p\le2h}p^{-\widetilde L_p}\), retaining the source’s small-prime and inner exponents but replacing each outer exponent by \(v_p(S_K)+2\sum_iw_{p,i}-r_p\); for \(p>K\), the determinant exponent is zero. The exact piecewise definition and local bases are in [the assembly audit](fauzan-simplified-proof-audit.md). Put

\[
\widetilde Q_K=\widetilde m_{K,10000}S_K\Delta_{K,5}.
\]

The **completed end-to-end review** establishes coefficientwise integrality, degree \(37n\), strict positivity at \(\zeta(5)\), and eventually

\[
0<\widetilde Q_{40n}(\zeta(5))<e^{-4n^2}.
\]

A hypothetical rational value \(a/b\) makes \(b^{37n}\widetilde Q_{40n}(a/b)\) a positive integer tending to zero. This is a complete reviewed computer-assisted alternative derivation of the preprint’s claimed theorem, retaining its determinant and local arithmetic engine. It is not a newly invented construction, an elementary proof, formal verification, or external acceptance. The enormous sufficient local threshold \(K\ge20{,}000{,}000{,}000\) and unspecified asymptotic onset make this a proof-design result rather than a practical approximation algorithm. No higher-odd normalization follows.

## 15. Positive moments force real roots and locate the target between them

Fix odd \(s\ge3\), integers \(K>N\ge0\), \(h=K-N\), and even \(r\ge2\). Let \(D_m(t)=\prod_{j=1}^m(t+j^2)\), and form

\[
G(X)_{ab}=\mu_{s,X}\!\left(\frac{D_N(t)^r t^{a+b}}{D_K(t)}\right),
\qquad 0\le a,b<h.
\]

Then \(\det G(X)\) has degree exactly \(h\), all its roots are real, and none equals \(\zeta(s)\). The number below \(\zeta(s)\) is the number of odd integers in \(\{N+1,\ldots,K\}\); the number above it is the number of even integers. Multiplicities are allowed.

The proof uses two exact facts. First \(G_0=G(\zeta(s))\) is a positive definite Gram matrix, so \(G_0^{-1/2}[X]G\,G_0^{-1/2}\) is real symmetric. Second the coefficient matrix is a Vandermonde congruence of the diagonal residue matrix

\[
c_j=\frac{j^{s-1}D_N(-j^2)^{r-1}}
 {(D_K/D_N)'(-j^2)},\qquad \operatorname{sign}c_j=(-1)^{j-1}.
\]

These residues are nonzero. Their signs determine the inertia, and a positive eigenvalue produces a root below the target. When both parities occur, \(G(X)\) is positive definite exactly between the nearest root on either side of \(\zeta(s)\). This produces algebraic enclosures whose rational polynomial coefficients can be computed without knowing the target value.

For example, \(N=0,K=2,s=3\) yields

\[
\frac{17-\sqrt5}{16}<\zeta(3)<\frac{17+\sqrt5}{16}.
\]

The [matrix-pencil proof](real-rooted-zeta-pencils.md) applies standard spectral theory to the constructed moments; it does not assert distinct roots, a shrinking-rate theorem, or adequate denominator control. Its proof and exact example received editorial checking; no separate completed review is recorded here. No historical novelty is claimed.

## 16. Bernoulli numbers connect the same rational data to two geometries

For fixed odd \(s\ge3\), set \(L(x^k)=B_k\), with \(B_1=-1/2\), and

\[
\tau_s(g)=\frac{L(g^{(s-2)})}{(s-1)!},
\qquad
\tau_s(e^{tx})=\frac{t^{s-1}}{(s-1)!(e^t-1)}.
\]

The generating function makes a finite difference into one boundary derivative:

\[
\tau_s(g(x+1)-g(x))=\frac{g^{(s-1)}(0)}{(s-1)!}.
\]

Adjoining simple integer poles and setting \(\tau_{s,X}(1/x)=-X\) forces their harmonic values by this boundary rule and reflection. Pulling back through \(t=-x^2\), with the factor \(x^s\), gives exactly the positive real moment functional \(\mu_{s,X}\). Its correction terms \(-1/(s-1)+1/(2j)\) are consequences of polynomial division and the two paired harmonic indices.

On polynomials the same Bernoulli functional is the p-adic limit

\[
L(g)=\lim_{m\to\infty}\frac1{p^m}\sum_{a=0}^{p^m-1}g(a).
\]

Residue-class subdivision gives the arithmetic distribution law. For primes \(p>s\), its rational extension has parameter change \(Y=p^sX+C_{p,s}\), with \(C_{p,s}\in p^s\mathbb Z_p\). This constant is p-adic, not asserted rational. Continuity is established on restricted analytic power series; no boundedness claim on all continuous functions is needed.

The [dictionary](bernoulli-local-global-dictionary.md) and [independent review](bernoulli-dictionary-review.md) prove these identities and their domains. This is an exact bridge between real positivity and local divisibility, not an analogy. It still leaves the difficult task of designing a determinant basis with sufficient local valuation gains for each higher odd index.


## 17. Newton interpolation and finite balance derive the numbers 17 and 5

Write \(\phi_k(n)=\binom nk\binom{n+k}{k}\). As a polynomial in \(z=n(n+1)\), it is the normalized Newton basis on the grid \(j(j+1)\). Multiplication by \(z\) raises its index. Transposing that rule leads to the nearest-neighbor operator

\[
(\mathcal Lf)(n)=\frac{(n+1)^3(f_{n+1}-f_n)+n^3(f_{n-1}-f_n)}{2(2n+1)},
\qquad
\mathcal L\phi_k=k(k+1)\phi_k+k^2\phi_{k-1}.
\]

Its coefficients are discovered by imposing this rule on \(\phi_1\) and \(\phi_2\); the contiguous-ratio calculation then proves it for all \(k\). The symmetric kernel \(K(n,m)=\sum_k\phi_k(n)\phi_k(m)\) satisfies \(\mathcal L_nK=m(m+1)K\), with \(A_n=K(n,n)\).

The nontrivial diagonal closure is

\[
A_n+A_{n-1}=6K(n,n-1).
\]

For \(w_k=[n\phi_k(n)/(n+k)]^2\), the finite balance relation
\(w_k(n^2-k^2)^2=(k+1)^4w_{k+1}\) gives
\(\sum k^2w_k=(n^2/2)\sum w_k\). Expanding the two adjacent Newton factors proves the closure. Substitution in the kernel equation gives the Apéry middle coefficient
\((2n+1)[17n(n+1)+5]\), with \(17=6\cdot3-1\) and \(5=6\cdot1-1\).

This proof does not assume the Apéry recurrence and does not illegally replace an off-diagonal kernel value by a diagonal one. It is explained in [the Newton note](newton-bispectral-motivation.md), proved and audited in [the closure review](newton-symmetric-square-review.md), and integrated into [the discovery guide](discovery-guide.md). Its birth–death interpretation also identifies the exact root-mean-square summation index \(n/\sqrt2\), but does not by itself prove concentration. Integrality and the zeta-valued limiting ratio still require their separate arguments.

## 18. The modular explanation is now checked through the global continuation step

Let \(\mathcal A(t)=\sum A_nt^n\), \(H(z)={}_2F_1(1/3,2/3;1;z)\),
\(T(x)=x(1-9x)/(1-x)\), and \(Z(x)=27x^2(1-x)/(1-3x)^3\). The exact analytic-germ identity is

\[
\mathcal A(T(x))=\frac{1-x}{(1-3x)^2}H(Z(x))^2.
\]

A contour integral realizes \(H(27w^3)\) as a period of the elliptic curve
\(Y^2=X^2(1-wX)^2-4w^2X\). Its square explains the third-order differential equation and cubic recurrence coefficients. Eta products then provide the integral generating series; a triple primitive provides cubic lcm denominator control.

The decisive analytic step is Fricke symmetry. A direct three-moment calculation identifies the constant as \(\zeta(3)\) and removes the nearest branch value \(\rho=(\sqrt2-1)^4\). The covering argument proves continuation to the larger disk \(|t|<\rho^{-1}\); it is not inferred from finite coefficient matching. [The completed modular review](apery-modular-review.md) checks the global branching, nonpolynomiality, arithmetic, and constant normalization in [the full construction](apery-modular-elliptic.md), retaining explicit standard modular-form and covering-theory inputs.

The same review proves a restricted weight-six failure: in the specified four-dimensional level-six Eisenstein span, the compatible fifth primitive has a nonremovable extra moment. Even the best scalar correction leaves a cubic local defect. This excludes that particular branch-cancellation attempt, not all modular approaches to \(\zeta(5)\).

## 19. Wilson orthogonality solves the one-pole model—and proves its arithmetic failure

For the zeta-three positive functional, take the size-\(m\) matrices
\(\mu_{3,X}(t^{i+j}/(t+1))\) and
\(\mu_{3,X}(t^{i+j+1}/(t+1))\). Their common positive-definite interval has exact endpoints

\[
L_m=H_m^{(3)}+\frac{m-1}{2m^3},\qquad
U_m=H_m^{(3)}+\frac1{2m(m+1)}.
\]

The Schur-complement minimizers are normalized Wilson polynomials. Their classical orthogonality norms make the endpoint increments telescope, proving these formulas for every size. In particular,

\[
\zeta(3)-L_m\sim U_m-\zeta(3)\sim\frac1{4m^4}.
\]

This alone would only show slow approximation. The stronger fact concerns the **actual reduced denominators** \(d_m^-\), \(d_m^+\). For every prime \(m/2<p<m\), \(p\ge3\),

\[
v_p(L_m)=v_p(U_m)=-3.
\]

The unique \(1/p^3\) term cannot be canceled by the endpoint corrections. The prime number theorem gives \(\log d_m^\pm\ge(3/2+o(1))m\). Therefore both cleared errors tend to \(+\infty\). No smaller integer multiplier, including one exploiting all rational cancellation, can make either endpoint sequence prove irrationality by a positive-integer-below-one argument. The [all-size Wilson proof and exact obstruction](one-pole-schur-structure.md) distinguish this verified failure from a merely pessimistic denominator upper bound. Other pole families are not ruled out.

## 20. Carry cancellation produces digit products and an algebraic series modulo p

For positive integers \(a,b\), let
\(B_n^{a,b}(T)=\sum_{k=0}^n\binom nk^a\binom{n+k}{k}^bT^k\).
For every prime \(p\), \(n\ge0\), and \(0\le r<p\),

\[
B_{pn+r}^{a,b}(T)=B_n^{a,b}(T^p)B_r^{a,b}(T)
\quad\text{in }\mathbb F_p[T].
\]

Lucas factorization separates the digits unless the low digits of \(n+k\) carry; in that case the second binomial vanishes modulo \(p\), on both sides of the identity. Consequently \(A_n\equiv\prod_iA_{n_i}\pmod p\), where \(n_i\) are the base-\(p\) digits. For example, a zeta-three Apéry number is a \(5\)-unit exactly when every base-five digit belongs to \(\{0,2,4\}\).

With \(P_p(t)=\sum_{r=0}^{p-1}A_rt^r\), the generating series satisfies

\[
\mathcal A(t)=P_p(t)\mathcal A(t^p),\qquad
P_p(t)\mathcal A(t)^{p-1}=1\quad\text{in }\mathbb F_p[[t]].
\]

Thus it becomes algebraic after reduction modulo \(p\). A sharper [independently proved identity](apery-digit-square-review.md), for \(p\ge5\), is

\[
P_p(t)=Q_p(t)^2\quad\text{or}\quad
P_p(t)=(1-34t+t^2)Q_p(t)^2.
\]

Here \(Q_p(0)=1\), its roots are simple, and it is coprime to \(1-34t+t^2\). The proof reduces an exact characteristic-zero differential first integral and analyzes root multiplicities; deriving only a zero derivative in characteristic \(p\) would be insufficient. The rule deciding which case occurs is proved separately in the [completed character review](apery-cm-independent-review.md): the extra quadratic factor occurs exactly when \(\left(\frac{-6}{p}\right)=-1\).

[The digit-factorization proof](apery-digit-arithmetic.md) connects this to the Fuss–Catalan project: carries are removed by a quotient there, while carry terms vanish under summation here. These are complete finite-digit mechanisms, not claims of historical novelty or stronger prime-power congruences.

## 21. Rational valuation weights need not increase the polynomial degree

The [reviewed superfactorial normalizer](pnt-free-normalizer-independent-review.md) and exact 25-term certificate replace prime-density averaging by a pointwise valuation majorant. For \(K=40n\), its repaired arithmetic rate is below \(6893/5000<1.379\), against a separately certified analytic magnitude exceeding \(1.38\). Both small-prime and large-prime repairs remain necessary.

Rational coefficients in the majorant do not force a power \(40000\) in the final degree. If \(F=gP\), with \(g>0\) rational and \(P\) primitive integral, and \(cF^q\in\mathbb Z[X]\) for \(c>0\), Gauss's lemma gives \(cg^q\in\mathbb Z_{>0}\). Hence, whenever \(F(\theta)>0\),

\[
0<P(\theta)\le[cF(\theta)^q]^{1/q}.
\]

Equivalently, once the repairs prove
\(v_p(EF)+b_p\ge0\) for rational \(b_p\), the integer \(v_p(EF)\) implies

\[
E\left(\prod_p p^{\lfloor b_p\rfloor}\right)F\in\mathbb Z[X],
\qquad
\log\prod_p p^{\lfloor b_p\rfloor}\le\sum_p b_p\log p.
\]

Thus floor rounding adds no upper arithmetic cost, including at negative valuations. The [independent review](prime-exponent-rounding-independent-review.md) proves both versions. The [complete guide](pnt-free-zeta5-guide.md) consequently ends with degree at most \(37n\), preserving the same quadratic decay budget. This is a simplification of normalization bookkeeping, not a new arithmetic theorem for all odd zeta values.

An [independently audited eleven-factor option](eleven-factor-normalizer.md) now uses the reviewed outer Hermite-kernel gain pointwise. Its repaired arithmetic cost is below \(6923/5000=1.3846\), so it requires the stronger analytic magnitude \(2771/2000=1.3855\). The exact remaining margin is \(9/10000\). The list is shorter but its rational coefficients are less pleasant and it adds a local lemma; it does not supersede the primary 25-term route. Its common denominator \(3545916\) is auxiliary, with final degree still \(h\).

## 22. An adaptive adjacent-Legendre partner cancels the entire unwanted boundary

Let \(P_n,R_n\) be as in Section 1, and set \(w_{n,a}(t)=(t-a)/(t+n+1)\). There is a uniquely determined rational \(a_n\) for all sufficiently large \(n\), bounded as \(n\to\infty\), such that

\[
\sum_{j=1}^n(R_nw_{n,a_n})'''(j)=0.
\]

The multiplier changes the polynomial partner only to a rational combination of two adjacent Legendre polynomials. The resulting full sum is a pure zeta-five linear form and equals its tail:

\[
\mathcal L_n=-\frac16\sum_{t\ge1}(R_nw_{n,a_n})'''(t)
=-\frac16\sum_{t>n}(R_nw_{n,a_n})'''(t)
\sim\pi\frac{\log n}{n^3}16^{-n}.
\]

This cancels the sum of boundary contributions, not every contribution separately. The [construction](adaptive-boundary-cancellation.md) and [completed independent review](adaptive-boundary-independent-review.md) prove the tuning, uniform contour comparison, exact constant, and eventual positivity. The sufficient arithmetic multiplier counts the extra denominators of the tuning parameter. Finite primitive-arithmetic samples show a serious cost but do not prove an all-index obstruction. This is an analytically successful repair with unresolved asymptotic arithmetic.

## 23. A signed single-zeta-seven model has nonzero determinants, but total variation is too costly

For even \(n\ge2\), the [signed model](signed-single-zeta7-model.md) uses

\[
R_n(t)=108^n(n!)^4
\frac{(t+1/2)_n(t+1/3)_n(t+2/3)_n}{(t)_{n+1}^7}
\]

and lattice-dilation weights \((45,-9,1)\). Reflection removes even zeta values, while
\(\sum_{b=1}^3w_bb^3=\sum w_bb^5=0\) and \(\sum w_bb^7=1080\) remove the two lower odd values. The resulting ordinary signed measure has moments in \(\mathbb Q+\mathbb Q\zeta(7)\).

The hypotheses of Fischler–Lai–Sprang–Zudilin's general nonvanishing Proposition 5.4 are checked directly. For integers \(1\le M\le2n-8\), some Hankel determinant with \(M\le m\le M+11\) is nonzero. Their large-\(s\) main theorem is not specialized to \(s=7\); the general proposition is the explicit external input. This does not imply nonvanishing at every size.

There is also a genuine obstruction to discarding the signs. For the stated normalized total-variation determinant \(J^{\rm abs}_{n,m}\), along admissible sizes with \(m/n\to c\in(0,2]\),

\[
\liminf\frac{\log J^{\rm abs}_{n,m}}{n^2}
\ge c\log54-c^2\log2>0.
\]

A single positive atom configuration proves this lower bound, using only an elementary lcm lower bound. The [independent review](signed-zeta7-independent-review.md) checks its normalization, endpoint \(c=2\), coefficient denominators, and nonvanishing window. It rules out this total-variation bound with the displayed multiplier; it does not bound the signed determinant from below or rule out a smaller primitive normalization. No irrationality theorem for \(\zeta(7)\) follows.

## 24. Euler depletion improves polynomial denominators but not automatically the local arithmetic

For every fixed odd \(s\ge3\), the positive functional with density

\[
\widetilde w_s(y)=2w_s(y/2)-w_s(y)>0
\]

satisfies

\[
\nu_s(t^e)=(2^{2e+2}-1)\mu_s(t^e),\qquad
2^{1+v_2(s-1)}\nu_s(t^e)\in\mathbb Z
\quad(e\ge0).
\]

The power of two is uniformly sharp. Euler/Genocchi identities remove every odd prime from these polynomial denominators. The pole moments instead isolate \(\eta(s)=(1-2^{1-s})\zeta(s)\):

\[
\nu_s((t+j^2)^{-1})
=(-1)^{j-1}j^{s-1}
\left(\eta(s)-\sum_{k=1}^j\frac{(-1)^{k-1}}{k^s}\right)+\frac1{2j}.
\]

The [completed Euler review](euler-depletion-independent-review.md) checks these identities, sharpness, exact divided-difference obstructions, and determinant scaling in [the construction](euler-depleted-positive-moments.md), with 1,686 exact checks. The polynomial Euler distribution law still has factor \(p^{1-s}\): losing the Bernoulli averaging factor costs one additional derivative. Extending this law to rational poles remains a separate task. Moreover, explicit paired poles have genuine coefficientwise losses, and the real weight decays at rate \(\pi\) rather than \(2\pi\). The finite determinant experiments do not establish an improved asymptotic normalizer.

The general integer-depletion construction now has its own [completed review](integer-depletion-independent-review.md). For odd \(s\ge3\) and integer \(\ell\ge2\), the positive density \(\ell w_s(y/\ell)-w_s(y)\) has the **least uniform polynomial denominator**

\[
C_{s,\ell}=\prod_{p\mid\ell}p^{1+v_p(s-1)}.
\]

The [cumulant derivation](integer-depletion-cumulants.md) identifies
\((\ell^k-1)B_k/k\) with the \(k\)-th cumulant of a discrete uniform variable on \(\{0,\ldots,\ell-1\}\). An integer-coefficient cumulant recurrence puts it in \(\mathbb Z[1/\ell]\); explicit no-carry witnesses establish sharpness at every prime dividing \(\ell\). Restoring the original real decay by the specified rescaling creates genuinely growing moment denominators along explicit subsequences. This is a result about those moments and that rescaling, not a lower bound for arbitrary determinant normalizers.

## 25. Higher Laurent geometry: distinguish the construction from the contour theorem

The [primary-source comparison](rodrigues-primary-source-comparison.md) establishes exact containment of the underlying rational-function family and pure derivative forms in Nesterenko's 2003 construction. His Barnes representation also contains the finite Fourier decomposition, and his multiple-integral theorem contains the positive central period. The Rodrigues polynomial is a classical confluent Jacobi–Piñeiro specialization. These are explicit parameter dictionaries, not merely resemblances. The Laurent coordinates, critical-value identities and sharp certified rates remain historically unattributed; that is not a novelty claim.

For \(m\ge1\), write \(C=\prod_{i=1}^{2m}(1+x_i)\), \(P=\prod_i x_i\). The higher-Rodrigues construction introduces

\[
\Lambda_m=(1+z)^mC\left(1+\frac C{zP}\right)^m,\qquad
\operatorname{CT}\Lambda_m^n
=\sum_{k=0}^{mn}\binom{mn}{k}^2\binom{n+k}{k}^{2m}.
\]

Its nonzero critical points reduce to
\((m-x)(1+x)^m=\pm x^{m+1}\), with value
\([m(1+x)/(m-x)]^{2m}\). The [geometry note](higher-rodrigues-laurent-geometry.md) and [completed independent \(m=2\) review](higher-rodrigues-asymptotic-independent-review.md) establish the sharp \(\zeta(5)\) contour asymptotic. If \(\lambda\) is the positive minimum of \(\Lambda_2\), the selected conjugate pair has modulus \(\lambda^{-1/2}\), rather than the reciprocal real minimum \(\lambda^{-1}\). The full form has sharp limsup rate \(\lambda^{-1/2}\) and infinitely many positive and negative terms. Saddle selection follows from the contour proof, not merely from listing algebraic critical values.

The [Eulerian decomposition](eulerian-rodrigues-period-decomposition.md) now has a [completed independent core review](eulerian-periods-independent-review.md). For every fixed \(m\), it expresses the pure form as a finite sum of contour components and proves a conditional component saddle theorem and an unconditional central-component asymptotic. A positive Mellin convolution additionally proves the central component positive at every index. The other components are not thereby controlled, and no uniform estimate as \(m\) grows is asserted.

The separate [finite-saddle audit](finite-rodrigues-saddles-independent-review.md) has now passed for \(m=2,3,4,5\). Exact rational certificates establish all required strip roots, logarithm branches, a unique dominant nonreal conjugate pair, and modulus below one. Combined with the component theorem, they prove

\[
\mathscr L_{n,m}=
\frac{2n^{-m-1/2}}{(2m-1)!}
\Re\{\mathcal C_m(t_\star)\eta_\star^n\}
+O(n^{-m-3/2}|\eta_\star|^n),
\quad
\limsup|\mathscr L_{n,m}|^{1/n}=|\eta_\star|<1,
\]

with infinitely many positive and negative terms. This includes the pure forms for \(\zeta(7),\zeta(9),\zeta(11)\), as well as \(\zeta(5)\); it does not assert that every individual term is nonzero.

**Review boundary:** these finite certificates do not cover \(m=6\) or all multiplicities. A numerical candidate outside the permitted contour strip at \(m=6\) is not a theorem that the full form grows. The proved decay still leaves insufficient arithmetic estimates; a sufficient denominator bound is not an actual primitive-denominator lower bound. No all-odd irrationality theorem follows.
