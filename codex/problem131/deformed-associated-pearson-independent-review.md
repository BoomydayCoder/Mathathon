# Independent audit: Pearson reciprocity and the associated deformed recurrence

Status: **accepted as a full proof**, after independent hand verification and exact free-parameter checks. This review concerns `deformed-associated-pearson-symmetry.md` together with `generic-sigma-diagonal-recurrence.md`. It uses the previously reviewed positive measure, polynomial moments and sampled resolvent on $0\le\sigma\le1/4$. It does not assert a new historical recurrence or any new irrationality theorem.

## Functional and second-kind calculation

The primitive definition of $\mathcal L_\sigma$ is independent of the primitive's constant. Its monomial values are
\[
\ell_d=\frac2{d+1}\sum_{h=1}^{\lfloor(d+1)/2\rfloor}
 {d+1\choose2h}B_{d+1-2h}\sigma^{h-1}.
\]
Here $B_1=-1/2$ is essential. In particular $(\ell_0,\ell_1,\ell_2,\ell_3)=(0,1,-1,(1+\sigma)/2)$. Bernoulli translation and reflection give the stated difference and antisymmetry identities. The Pearson rule follows from $v(0)=v(\pm\tau)=0$; at $\tau=0$, $(x^3h)''(0)=0$ proves the limiting case directly.

The claimed eigen-equation is correct. The invariant polynomial operator $D=L/[2(2x+1)]$ satisfies
\[
D\phi_k=k(k+1)\phi_k+(k^2-\sigma)\phi_{k-1}.
\]
One can verify this without appealing to spectral theory: reflection gives invariance, zeros at $0,\ldots,k-2$ and their reflected partners give divisibility by $\phi_{k-1}$, and the leading coefficient and value at $k-1$ determine the two coefficients. The displayed finite-sum coefficient ratio then proves the eigen-equation for every $n$; $k=0$ is treated separately.

The formal Cauchy series is legitimate coefficientwise: every coefficient is a polynomial in the functional variable. The Pearson flux constant is **2**, since the remainder after polynomial division has functional value $\mathcal L_\sigma(2x+y+1)=2$. The adjoint kernel identity has constant **minus 1**. Combining them gives exactly
\[
L\Theta_R-\Theta_{LR}=2\{R(x-1)-R(x+1)\}+\mathcal L_\sigma R.
\]
As a sensitive non-invariant test, $R=x$ gives $-3=-4+1$; omitting either constant breaks the identity. For invariant $\Pi_n$, $\mathcal L_\sigma\Pi_n=0$, and applying $\Theta$ to $(2x+1)\Pi_n$ introduces no extra term. Thus the inhomogeneous equation for $\Theta_n$ is correct.

The identification $N_{nm}=\Theta_n(m)$ is exact, not an asymptotic argument. Under the quadratic substitution its integrand is the polynomial quotient
\[
(x+\tfrac12)\frac{\Pi_n(x)-\Pi_n(m)}{(x-m)(x+m+1)}.
\]
Partial fractions give half the sum of the two linear quotients. Antisymmetry and $\Pi_n(-x-1)=\Pi_n(x)$ make their functional values equal. This also works at $n=0$ and $m=0$.

## Sources, symmetry and compatibility

For $m\ge1$, the tail satisfies $v(m+1)\Delta_+f_m=-1$ and $v(m)\Delta_-f_m=1$. The product-rule source is $K_{n,m-1}-K_{n,m+1}$ and is canceled by half the associated-polynomial source. At $m=0$ the product contributes $-\Pi_n(1)$ and the second-kind contribution is $-[\Pi_n(-1)-\Pi_n(1)]$, leaving **minus 1**. Consequently the $m$-grid source is $-\delta_{m0}$, with no negative-index remainder required.

The transported Legendre recurrence, integrated against the resolvent, gives the $n$-grid source $-\pi_n\nu_\sigma(\Phi_n)=-\delta_{n0}$. Orthogonality is applicable to these exact polynomial moments. In the $m=0$ column this grid is a constant-flux equation, yielding $F_{n0}=Z_\sigma-T_n=F_{0n}$. Transposition now has the same $n$-recurrence, source, and initial row. Since $v(n+1)>0$ throughout the positive interval, the index-zero equation determines the first step and forward recursion determines all subsequent values. Thus symmetry does not require an unmentioned boundary condition at infinity.

Commuting the two grid operators is legitimate including the sources: $L1=0$ and $\lambda_j\delta_{j0}=0$ remove all source terms. The commutator with $\lambda_n=n(n+1)$ is
\[
[L,\lambda]F_n=2(n+1)v(n+1)F_{n+1}-2nv(n)F_{n-1}.
\]
Substituting $m=n-1$ gives the two displayed off-diagonal eliminations and
\[
v(n)(F_{nn}+F_{n-1,n-1})=2n(3n^2-\sigma)F_{n,n-1}.
\]
At $n=1$, the apparent exceptional source is multiplied by $n-1$ and vanishes; $v(0)=0$ removes the negative-index term. The adjacent identity therefore starts at **1**, as required.

## Diagonal and initial condition

The resulting coefficients $\alpha_n,\beta_n,\gamma_n$ in the draft are correct, with $\gamma_0=0$ defined explicitly. The remainder equation is
\[
\alpha_nF_{n+1,n+1}-\beta_nF_{nn}+\gamma_nF_{n-1,n-1}=-\delta_{n0}.
\]
The denominator satisfies its homogeneous equation independently (finite telescoping, or the same compatibility argument without sources). Subtracting this equation from $F_{nn}=A_nZ_\sigma-b_n$ proves the numerator equation without assuming irrationality of $Z_\sigma$. In particular
\[
b_0=0,\qquad \alpha_0b_1=1,\qquad
b_1=\frac{2(3-\sigma)}{(1-\sigma)^2}.
\]
The numerator is homogeneous only for $n\ge1$. At $\sigma=0$, $b_1=6$; this normalization should not be confused with notes using a remainder twice as large.

The generic finite telescope in the companion note is also correct: its endpoint at $k=0$ vanishes and its terminating endpoint vanishes. The transported operator $\tfrac12(uJ+Ju)-\sigma u+\tfrac14$, $J=-D[u(1-u)D]$, has the claimed Newton-basis action. Rational continuation of the finite identities is valid wherever their denominators are finite; positive-measure arguments are asserted only on the reviewed interval.

## Reproducible checks and scope

Run:

```
python codex/computations/check_deformed_pearson_independent.py
python codex/computations/check_deformed_pearson_symmetry.py
python codex/computations/check_generic_sigma_diagonal_recurrence.py
```

All passed. The independent checker keeps $\sigma$ symbolic. It checks arbitrary-polynomial Pearson/reflection/commutator identities through degree 14, the full rational adjoint kernel, directly constructed associated numerators, symmetry, both grid directions (by symmetry after a direct sourced grid check), adjacent balances and the diagonal source/initial normalization. The author's rational checker additionally tests six parameters, degrees through 18, and grids through index 6. These finite checks support the all-index proof audited above; they do not replace it.

The denominator recurrence is the exact classical Askey–Wilson 1984 specialization documented in `koecher-askey-wilson-primary-dictionary.md`. This review accepts the independently derived numerator reconstruction and does not treat the denominator's historical identification as evidence for it. No remaining mathematical gap was found in these two recurrence proofs. Arithmetic denominator growth, irrationality outside the original Apéry specialization, and priority of this particular proof remain separate questions.
