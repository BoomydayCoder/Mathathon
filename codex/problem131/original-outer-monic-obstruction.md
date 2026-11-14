# Actual denominator costs for the original odd-zeta determinant

**Status: proved local coefficient theorem and an actual scalar-normalization obstruction for fixed odd s>=29.** This concerns the original, undepleted Hermite functional, numerator power six, and K=40n,N=3n,h=37n. It is not inferred from a sufficient normalizer. The s=7,9,...,27 cases are not ruled out by the bound proved here. No general irrationality obstruction or novelty claim is asserted.

## 1. A rational basis change makes the monic arithmetic visible

Put
\[
 D_j(t)=\prod_{r=1}^j(t+r^2),\quad D=D_K/D_N,\quad W=D_N^5,
\]
\[
 \Delta_{s,K,N}(X)=\det_{0\le i,j<h}\mu_{s,X}(Wt^{i+j}/D),\quad
 F_{s,K,N}(X)=\Delta_{s,K,N}(X)/[X^h]\Delta_{s,K,N}(X).
\]
The leading coefficient is the nonzero integer
\[
 \ell_K=(-1)^{h(h-1)/2}\prod_{j=N+1}^K j^{s-1}D_N(-j^2)^5.
\]
Use the integer polynomials
\[
 f_j(t)=D(t)/(t+j^2),\qquad N<j\le K,
\]
as a basis of the h-dimensional polynomial space. Their rational change-of-basis determinant is nonzero; it need **not** be a p-adic unit, because it cancels between a determinant and its leading coefficient.

In this basis the X-coefficient matrix is diagonal:
\[
 C_j=j^{s-1}R_j,\qquad R_j=W(-j^2)D'(-j^2).
\]
For j!=k, the rational integrand W f_j f_k/D is the integer polynomial
\(WD/((t+j^2)(t+k^2))\). For j=k, monic division gives
\[
 Wf_j^2/D=Q_{jj}(t)+R_j/(t+j^2),\quad Q_{jj}\in\mathbb Z[t].
\]
Write Q_jk for the off-diagonal polynomial too, and set
\[
 P_{jk}=\mu_s(Q_{jk}),\qquad
 c_j=-j^{s-1}H_j^{(s)}-\frac1{s-1}+\frac1{2j}.
\]
Then the constant matrix is B=P+diag(R_j c_j), and exactly
\[
 \boxed{F_{s,K,N}(X)=\det(XI+C^{-1}B).}                 \tag{1}
\]
Every Q_jk has degree at most E=K+4N-2.

## 2. The outer-band matrix has an explicit reduction

Assume K=40n,N=3n, fixed odd s>=7, and a prime p with
\[
 K/2<p\le K,\qquad p>s-1,\qquad p^2>2K.
\]
The last two hold eventually throughout the band for fixed s. All P_jk are p-integral. Indeed the previously proved exact Bernoulli cancellation rule gives integrality through degree E whenever
\[
 p>1+\frac{2E+2}{s-1}.
\]
For s>=7 and K>=40 this follows from p>K/2, since the right side is at most 13K/30+2/3.

Let a in {1,...,(p-1)/2} index an ordinary squared residue class. There are at most four original nodes j<=K congruent to +/-a. Let their number be ell_a. The condition 2N<p guarantees at most one removed node in this class; denote its presence by delta_a. Since p^2>2K, each ordinary same-class difference has valuation exactly one. Consequently, for a surviving ordinary node,
\[
 v_p(C_j)=5\delta_a+(\ell_a-\delta_a-1)
 =\ell_a+4\delta_a-1\le7.                            \tag{2}
\]
At the unique zero-class node j=p, all node differences and W(-p^2) are units, so
\[
 v_p(C_p)=s-1.                                      \tag{3}
\]
The harmonic constants satisfy
\[
 \overline{p^s c_j/j^{s-1}}=
 \begin{cases}0,&j<p,\\-1/2,&j=p,\\-1,&j>p.\end{cases} \tag{4}
\]
For j>p only the harmonic summand p^(-s) contributes; at j=p, the extra 1/(2j) contributes +1/2 after the scaling. The constant -1/(s-1) vanishes after scaling in all cases.

If s>=9, (2)-(3) give v_p(C_j)<s everywhere. Thus
\[
 A=p^s C^{-1}B\in M_h(\mathbb Z_p)
\]
reduces to the diagonal matrix in (4): the entire polynomial-part matrix vanishes modulo p after row scaling. Therefore
\[
 \boxed{\overline{p^{sh}F_{s,K,N}(Z/p^s)}
 =Z^{p-N-1}(Z-1/2)(Z-1)^{K-p}.}                     \tag{5}
\]
This is a coefficientwise identity over F_p, not a statement merely about one value of the determinant.

## 3. An exact coefficient valuation forces a denominator

Put b=K-p+1. The coefficient of Z^(h-b) on the right of (5) is (-1)^b/2, a unit. Hence
\[
 \boxed{v_p([X^{h-b}]F_{s,K,N})=-sb,\qquad s\ge9.}    \tag{6}
\]
Other coefficients could force still larger denominators; (6) is not being misreported as the exact Gauss valuation of the entire polynomial.

Let M_K be the least common multiple of all coefficient denominators of the monic polynomial F. It is a positive integer, and M_K F is primitive integral. Thus (6) proves the actual lower bound
\[
 \boxed{v_p(M_K)\ge s(K-p+1),\quad K/2<p\le K.}        \tag{7}
\]
There can be no cancellation of this cost by numerator content at other primes after monic normalization.

The prime number theorem, used only here to summarize the exact finite theorem, gives
\[
 \boxed{\liminf_{K\to\infty}\frac{\log M_K}{K^2}
 \ge s\int_{1/2}^{1}(1-u)\,du=\frac{s}{8}.}           \tag{8}
\]
The fixed exceptional primes leave the band; the +1 and endpoint rounding contribute o(K^2).

## 4. Actual failure of every scalar normalization for s>=29

The independently certified uniform-trial Jensen calculation in
`coprime-depletion-primitive-obstruction.md` and its independent review gives the monic real lower exponent greater than -1789/500=-3.578. The same lower bound applies directly to the original density here: the original-density first-exponential lower bound is precisely the one used there, and the depleted leading-coefficient factor C_s^h changes only O_s(K) in the logarithm. For every fixed odd s, the y^s and j^(s-1) factors cancel their K-log powers and affect only lower-order terms. Equivalently one may repeat that identical uniform-trial calculation with the original density and C_s=1.

Consequently
\[
 \liminf_{K\to\infty}\frac1{K^2}
 \log|M_K F_{s,K,N}(\zeta(s))|
 \ge\frac{s}{8}-\frac{1789}{500}.                    \tag{9}
\]
For s>=29 the right side is at least 47/1000>0. In particular, for all sufficiently large K=40n,
\[
 \boxed{|M_KF_{s,K,N}(\zeta(s))|>\exp(K^2/40),\quad
 s\ge29\text{ fixed and odd}.}
\]
Every rational scalar making the original determinant polynomial integral produces an integer multiple, up to sign, of this primitive polynomial. Thus **no scalar normalization of this specified determinant family can yield small nonzero integral-polynomial values** for these odd indices. This conclusion uses an actual coefficient denominator and a real lower bound; it is stronger than failure of an available sufficient arithmetic budget.

The result does not rule out s=7,9,...,27, different proportions, altered numerators, or other irrationality constructions.

## 5. The precise residual block at s=7

For s=7, the preceding clean diagonal proof still holds whenever
\[
 p>(K+N)/2.
\]
Indeed a removed class can have four nodes only if its removed index a satisfies 2p-a<=K, hence a>=2p-K. No such a<=N exists in this narrower band.

In the remaining short interval K/2<p<=(K+N)/2, put
\[
 d=K+N-2p+1>0.
\]
Exactly d removed classes have four nodes, indexed by 2p-K<=a<=N. Each has three surviving nodes p-a,p+a,2p-a. Let J be this set of 3d nodes. These and only these rows have v_p(C_j)=7. The matrix A=p^7C^(-1)B is still integral; all rows outside J reduce diagonally as in (4). On J the surviving residual block is explicitly
\[
 \boxed{R_{jk}=\delta_{jk}\epsilon_j+
 \left(\frac{C_j}{p^7}\right)^{-1}\overline{\mu_7(Q_{jk})},
 \quad j,k\in J,}                                  \tag{10}
\]
where epsilon_j=0 for j<p and -1 for j>p. There is no zero-class node in J.

Rows outside J allow an exact block-triangular expansion, regardless of the couplings from J to its complement. Thus the characteristic reduction is its outside diagonal product times det(ZI+R). Of the b=K-p+1 outside-and-inside nonzero diagonal entries in (4), precisely 2d occur inside J. Therefore it has at least b-2d nonzero roots counted with multiplicity, giving
\[
 \boxed{v_p(M_K)\ge7(b-2d).}                        \tag{11}
\]
More precisely, if r is the greatest index of a nonzero coefficient of det(ZI+R), written in descending powers, then one coefficient of F has valuation exactly -7(b-2d+r).

This produces the unconditional s=7 band summary
\[
 \liminf\frac{\log M_K}{K^2}
 \ge7\left(\frac18-\frac{\alpha^2}{2}\right)
 =\frac{2737}{3200},\qquad \alpha=3/40.
\]
It is far too small, by itself, to prove primitive growth. Determining whether the concrete residual determinant det R vanishes modulo p is a well-defined remaining arithmetic question. It is not answered by the real positivity of the original functional.

## 6. Independent exact checks and limits

`codex/computations/check_original_outer_monic.py` independently constructs exact Bernoulli moments, the integer basis leading factors, and the finite-field residual matrix. At K=40,N=3, it verifies (6) for s=7,9,11 and p=23,29,31,37. At K=80,N=6,s=7, the residual blocks for p=41 and p=43 have sizes 15 and 3 and are nonsingular modulo their primes. They force denominator exponents at least 315 and 273, respectively; the clean band case p=47 forces exponent at least 238. These residual nonsingularity observations are finite checks, not an all-prime theorem.

The proof of (8)-(9) uses the prime number theorem and the separately certified Jensen lower bound explicitly. No general s>=7 failure theorem is claimed, and no inference is drawn from a sufficient normalizer alone.
