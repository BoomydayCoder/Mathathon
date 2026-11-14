# Optimizing a cut under a denominator budget

Status: independent audits and all-alpha statements proved below. The obstruction concerns one specified sufficient multiplier, not the least possible denominator and not the irrationality of zeta(5).

## 1. Independent audits

**discovery-guide.md.** Checked the residue normalization, factor of two in B_n, short certificate's pole conditions, separate A recurrence, bounded F estimate, Casoratian sign, growth bound, period-six factorial lemma, and exact comparison 16875>16384. No mathematical error found. Recurrences involving n-1 are used for n>=1; stating that range explicitly near the certificate would improve presentation. The integral supplement's Mellin sign and positive integral bound also check.

**optimal-cut-higher-derivatives.md.** The contour orientation and integration-by-parts sign agree with the independently audited m=2 case. With y=csc²(pi z),

    D²(y^j)=pi²[-4j²y^j+2j(2j+1)y^(j+1)],

which gives the coefficient recurrence and C_(m,1)=(-1)^m4^(m-1). Thus the F_n coefficient is exactly

    (-1)^(m+1)(2pi)^(2m-2)/(2m)!.

The remaining terms have integrable sech weights and are O_m(U_n(c)²), giving relative O_m(n^-1/2) at the saddle. The m=1 normalization is F_n/2, m=2 gives -zeta(2)F_n, and m=3 gives (pi^4/45)F_n. Performing contour summation on R_n before introducing U_n's right-side gamma poles is essential and correctly done. No error found.

**hilbert-matrix-connection.md.** Orthogonality gives C H C^T=D^-1, hence H^-1=C^T D C. The inverse-entry integrality and trace identity follow. Rodrigues integration gives norm 1/(2j+1) with the stated signs. The trace/operator-norm inequalities and 1<=||H||<=N correctly yield the exponential condition-number rate. No sharp spectral theorem is silently assumed. No error found.

The accompanying script independently checks the cotangent derivative polynomials and Hilbert identities. These audits establish the displayed mathematics, not historical novelty.

## 2. Fixed-cut analytic rate

For 0<alpha<1 put k_n=floor(alpha n). The contour theorem proves

    E_(n,k)=A_n zeta(5)-(B_n+C_(n,k)),
    lim |E_(n,k_n)|^(1/n)=exp[-phi(alpha)],
    phi(alpha)=2[(1+alpha)log(1+alpha)
                  -(1-alpha)log(1-alpha)-2alpha log alpha].

For alpha!=1/sqrt2 the correction integral dominates F_n and is asymptotic to

    U_n(k_n+1/2)² * b/[4sinh b],
    b=log[(1-alpha²)/alpha²].

It is eventually positive. At alpha=1/sqrt2 with bounded rounding error, E/F_n tends to -zeta(2), so the form is eventually negative with the same exponential rate. Nonvanishing is thus established before taking logarithms.

The maximal analytic decay occurs at alpha=1/sqrt2, where phi=4log(1+sqrt2). Arithmetic optimization has a different objective.

## 3. Optimize the stated sufficient multiplier

The safe arithmetic multiplier is

    M_(n,k)=4D_n²D_(n+k)³,
    M_(n,k)(B_n+C_(n,k)) in Z.

The proved elementary factorial bound log D_N<=c6 N+O(log²N), with c6=log2/2+3log3/4, gives the sufficient upper-bound budget

    log(M|E|)/n <= c6(5+3alpha)-phi(alpha)+o(1).

For any c>0 define G_c(alpha)=c(5+3alpha)-phi(alpha). Since

    phi'=2log[(1-alpha²)/alpha²],
    phi''=-4/[alpha(1-alpha²)]<0,

G_c is strictly convex and has unique minimum at

    alpha_c=[1+exp(3c/2)]^(-1/2).

The identity phi-alpha phi'=2log[(1+alpha)/(1-alpha)] gives

    min G_c=5c-2log[(1+alpha_c)/(1-alpha_c)].

For c=c6, alpha_c≈0.383819 and the minimum is ≈4.23454>0. Thus the elementary upper-bound estimates never prove a successful squeeze.

If one additionally uses the prime number theorem, c=1 describes the actual exponential size of this specified multiplier. Its optimum is alpha≈0.427113, with residual exponent≈3.17456. This benchmark is informative but unnecessary for the elementary failure theorem below.

The budget-optimal cut is therefore to the left of the best analytic saddle. Making the real error smaller can make normalization harder.

## 4. Elementary theorem: this multiplier actually fails at every fixed cut

An oversized upper estimate alone shows only that the estimate is insufficient. Here we can prove more about the specified multiplier.

**Lemma.** D_N>=2^N/(N+1) for every integer N>=1.

**Proof.** Legendre's formula expresses v_p binom(N,K) as a sum over p^a<=N of terms floor(N/p^a)-floor(K/p^a)-floor((N-K)/p^a), each zero or one. Thus every binomial coefficient in row N divides D_N. One coefficient is at least the average 2^N/(N+1). QED.

Consequently

    M_(n,k)>=4*2^(5n+3k)/[(n+1)²(n+k+1)³],

and the exact analytic rate gives

    liminf log(M_(n,k_n)|E_(n,k_n)|)/n
       >=(5+3alpha)log2-phi(alpha).

**Theorem.** This lower bound is greater than log2 for every alpha in (0,1). Hence the named multiplier times the error grows exponentially.

**Proof.** Minimize G_c at c=log2. Its minimizer is

    r=[1+2^(3/2)]^(-1/2)<3/5,

because 1+2sqrt2>3>25/9. Hence (1+r)/(1-r)<4, and

    min G_(log2)=5log2-2log[(1+r)/(1-r)]>log2.

This is an exact inequality, not a numerical-grid conclusion. QED.

At k=0 the untruncated form decreases only polynomially. At k=n the full-tail decay base is 1/16 while the multiplier's elementary lower rate is 8log2. These endpoint choices also fail with the named multiplier.

**Scope.** This rules out the direct integer squeeze using M=4D_n²D_(n+k)³ at every fixed proportional cut. It does not show that the reduced denominator of B_n+C_(n,k) is so large. A smaller denominator or a common-divisor cancellation involving A_n could change the problem. No reduced-denominator asymptotic or general impossibility theorem is claimed.

## 5. A simple concrete cut: halfway

For k=floor(n/2), phi(1/2)=3log3=log27. Thus the form is eventually positive with decay base 27^-n. Under the c=1 benchmark its gap is 13/2-3log3≈3.20416, close to the optimum≈3.17456.

Leading constants show a rounding effect. For odd n, k=(n-1)/2 and c=k+1/2=n/2 exactly. Stirling gives

    U_n(c)²~(16pi²/3)n^-2 27^-n.

The contour factor has b=log3 and sinh b=4/3. Therefore

    E_(n,(n-1)/2)~pi²log3 * n^-2 27^-n  (n odd).

For even n and k=n/2, c=n/2+1/2; this fixed half-step shifts U_n(c)² by the limiting factor exp(-log3)=1/3. Hence

    E_(n,n/2)~(pi²log3/3)n^-2 27^-n  (n even).

Quoting a single leading constant for both parities would be incorrect. This rational cutoff is analytically clean, but its denominator cost still defeats the specified multiplier.

## 6. Next design requirement

Moving the cut alone cannot rescue this standard lcm multiplier. Further progress requires improved arithmetic, altered coefficients, or a different approximation object. Numerical optimization should record the reduced denominator and error together. The successful main zeta(3) proof and the all-m analytic limit do not remove this separate arithmetic requirement.

Run check_cut_budget_optimization.py for exact Hilbert/cotangent checks, numerical sanity checks of the proved convex inequality, and half-cut asymptotic diagnostics.
