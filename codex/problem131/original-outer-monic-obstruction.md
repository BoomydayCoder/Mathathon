# Actual denominator costs for the original odd-zeta determinant

**Status: proved and independently reviewed local coefficient theorem and actual scalar-normalization obstruction for fixed odd s>=13.** See `original-wide-monic-independent-review.md` and its separate ideal-repair checker. This concerns the original, undepleted Hermite functional, numerator power six, and K=40n,N=3n,h=37n. It is not inferred from a sufficient normalizer. The s=7,9,11 cases are not ruled out by the bound proved here. Sections 1–6 first establish the simpler s>=29 result; Section 7 strengthens it to s>=13. No general irrationality obstruction or novelty claim is asserted.

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

The proof of (8)-(9) uses the prime number theorem and the separately certified Jensen lower bound explicitly. No general s>=7 failure theorem is claimed, and no inference is drawn from a sufficient normalizer alone. Section 7 below strengthens the threshold to 13.

## 7. A bounded zero-class repair improves the threshold to s>=13

The parent proposed extending the argument to p>K/4; an independent derivation gives the same result. The only exceptional rows are the at most three zero-class nodes. Their finite p-adic cost can be removed without losing a linear number of forced coefficient denominators.

### 7.1 All local hypotheses in the wider band

Fix odd s>=13, K=40n, N=3n, and a prime
\[
 K/4<p\le K.
\]
For all sufficiently large K (depending on fixed s), p>s-1, p^2>2K, 2N<p, and
\[
 p>1+(2E+2)/(s-1),\qquad E=K+4N-2.
\]
Indeed the last right side is at most 13K/60+5/6<K/4 for K>=40. Thus the entire polynomial-part matrix P is p-integral.

Put q=floor(K/p), so 1<=q<=3. An ordinary squared-residue class has at most 2q+2<=8 nodes. At most one was removed, so (2) now gives
\[
 v_p(C_j)\le 2q+5\le11<s\quad(p\nmid j).             \tag{12}
\]
At a zero-class node j=ap, 1<=a<=q, the other q-1 zero nodes each contribute exactly two to v_p(D'(-j^2)). All other factors are units. For sufficiently large p, the factors a^2-b^2 are units since a,b<=3; therefore
\[
 v_p(C_{ap})=s-1+2(q-1)=s+2q-3.                     \tag{13}
\]
In A=p^sC^(-1)B the ordinary rows are integral, with reduction diagonal. Their diagonal residues are
\[
 \epsilon_j=-H_{\lfloor j/p\rfloor}^{(s)}\pmod p,
 \qquad p\nmid j.                                 \tag{14}
\]
The off-diagonal entries vanish modulo p by (12). For the zero rows, the pole contribution is integral too:
\[
 p^s c_{ap}/(ap)^{s-1}
 \equiv-H_a^{(s)}+\frac1{2a^s}\pmod p.
\]
Only their polynomial-part entries can have negative valuation, bounded below by -(2q-3). Multiply each zero row by p^d, where
\[
 d=(2q-3)_+,\qquad E_0=qd\le9.
\]
The resulting polynomial matrix with variable Z is integral. Its determinant is
\[
 \mathcal Q(Z)=p^{E_0}\det(ZI+A)
 =p^{E_0+sh}F(Z/p^s)\in\mathbb Z_p[Z].              \tag{15}
\]
Its leading coefficient is p^(E_0), so its Gauss valuation v satisfies
\[
 0\le v\le E_0.                                    \tag{16}
\]
No Schur-complement invertibility or residual-minor nonvanishing has been assumed.

### 7.2 A polynomial-ideal lemma controls every cancellation

For a lift epsilon in Z_p of a nonzero residue, let n_epsilon be the number of ordinary rows having residue epsilon in (14). Every entry in each of these rows of the repaired matrix belongs to
\[
 I_\epsilon=(p,Z+\epsilon)\subset\mathbb Z_p[Z]:
\]
the diagonal entry is Z+epsilon modulo p, and every off-diagonal entry is zero modulo p. Multilinearity of the determinant implies
\[
 \mathcal Q\in I_\epsilon^{n_\epsilon}.              \tag{17}
\]

**Elementary lemma.** If Q in Z_p[Z] belongs to (p,Z+epsilon)^n and v=v_p^G(Q), then the reduction of p^(-v)Q is divisible by
\((Z+\overline\epsilon)^{(n-v)_+}\).

**Proof.** Set Y=Z+epsilon. Membership in (p,Y)^n means that the coefficient of Y^j has valuation at least n-j for j<n. Translation by an integral constant preserves Gauss valuation: substitution in either direction has integral binomial coefficients. Thus after division by p^v and reduction, every Y^j coefficient with j<n-v vanishes. QED.

Applying the lemma and (16), the nonzero polynomial
\(\overline{p^{-v}\mathcal Q}\) is divisible by
\[
 \prod_{\epsilon\ne0}(Z+\epsilon)^{(n_\epsilon-E_0)_+}.       \tag{18}
\]
For fixed s, the rational numbers H_1^(s),H_2^(s),H_3^(s) are positive and distinct. Apart from finitely many primes dividing their numerators, denominators or pairwise differences, their reductions are nonzero and distinct. All these fixed exceptions are eventually outside p>K/4. Consequently there are at most q distinct nonzero groups, containing in total
\[
 b_{\rm ord}=K-p+1-q
\]
ordinary nodes. Formula (18) therefore forces at least
\[
 R=\max(0,b_{\rm ord}-qE_0)                         \tag{19}
\]
nonzero roots, counted with multiplicity, in the primitive reduced polynomial.

This argument handles all bottom-row dependence on Z, all couplings, and arbitrary cancellation inside their determinant. It needs no assumption that a particular finite-field minor is nonzero.

### 7.3 Convert the retained roots into an actual coefficient denominator

Write F(X)=sum_(r=0)^h a_r X^(h-r), a_0=1. The coefficient of Z^(h-r) in (15) is
\[
 p^{E_0+sr}a_r.
\]
A polynomial of degree at most h having R nonzero roots counted with multiplicity has some nonzero coefficient with r>=R: its zero-root multiplicity is at most h-R. Hence some such r survives in \(\overline{p^{-v}\mathcal Q}\), and
\[
 v_p(a_r)=v-E_0-sr\le-sR.
\]
Thus the **actual** monic coefficient denominator satisfies
\[
 \boxed{v_p(M_K)\ge
 s\max(0,K-p+1-q-q^2(2q-3)_+),\quad K/4<p\le K.}   \tag{20}
\]
The loss is bounded by 30, independently of K and p. Even primes close to K, where the displayed lower bound can be zero, affect only the same lower-order error.

The prime number theorem and partial summation now give
\[
 \boxed{\liminf\frac{\log M_K}{K^2}
 \ge s\int_{1/4}^1(1-u)\,du=\frac{9s}{32}.}         \tag{21}
\]
The constant loss per prime contributes O_s(K), not a quadratic term. This use of the prime number theorem is explicit; the finite local theorem does not require it.

Combining with the same independently certified Jensen bound as Section 4 yields
\[
 \liminf\frac1{K^2}\log|M_KF(\zeta(s))|
 \ge\frac{9s}{32}-\frac{1789}{500}.
\]
For s=13 the right side is
\[
 \frac{117}{32}-\frac{1789}{500}
 =\frac{313}{4000}=0.07825>0.
\]
Therefore for **every fixed odd s>=13**, all sufficiently large K=40n satisfy
\[
 \boxed{|M_KF(\zeta(s))|>\exp(K^2/20).}              \tag{22}
\]
As before, every rational scalar producing an integral polynomial gives an integer multiple of this primitive representative. This proves an actual obstruction for the original undepleted determinant family at these indices, not merely for a proposed normalizing formula. The remaining fixed odd cases s=7,9,11 require additional ideas; the present argument does not decide them.

### 7.4 Full characteristic-polynomial checks

The separate standard-library checker `codex/computations/check_original_wide_monic.py` constructs the full rational moment matrix modulo a sufficiently high power of p. It computes its complete characteristic polynomial using Newton identities, explicitly increasing precision to account for divisions by characteristic indices divisible by p. It then reverses the auxiliary uniform scaling, extracts the exact Gauss content through the guaranteed bound E_0, and verifies the retained-root multiplicities of the primitive reduction. Thus it checks the nontrivial cancellations in the zero block, rather than merely the easy diagonal reductions.

All three cases pass:

| K,N,s,p | Repair E_0 | Actual Gauss content of Q | One forced coefficient denominator exponent | General bound (20) |
|---|---:|---:|---:|---:|
| 40,3,13,13 | 9 | 6 | 367 | 0 |
| 40,3,13,17 | 2 | 1 | 313 | 234 |
| 80,6,13,23 | 9 | 6 | 757 | 364 |

The exponents in the fourth column are exact valuations of particular detected coefficients, not a claim to have found the least denominator of the entire polynomial. The finite examples support the all-parameter proof above; the general theorem does not infer valuations or root multiplicities from the samples.
