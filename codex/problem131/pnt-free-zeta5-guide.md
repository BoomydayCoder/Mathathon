# A prime-number-theorem-free route through the zeta-five determinant

This guide assembles a complete computer-assisted alternative derivation of the irrationality theorem claimed in [Aabir Fauzan's preprint, version 1](https://doi.org/10.5281/zenodo.22826419). It retains his rational functional, determinant, and local arithmetic construction. The changes are an explicit equilibrium measure for the real estimate and an explicit superfactorial product for the arithmetic estimate. The latter removes the prime number theorem (PNT) and every asymptotic prime-counting input.

The mathematical components have independent written reviews, linked below. This is a coherent assembly of those components, not a new independent invention of the determinant, a formally verified proof, or evidence of external acceptance of the preprint. Historical novelty has not been established. “PNT-free” does not mean the entire proof is elementary or short: the local interpolation lemmas, positive moment identity, potential argument, and finite exact certificates remain substantial.

The guiding question is concrete: **can we make a positive, extremely small determinant into an integer polynomial without estimating the distribution of primes?** The answer supplied by the construction below is yes. Residue-class overlaps produce Bernoulli polynomials; these become exact superfactorial valuations; a finite majorant then controls all primes individually.

## 1. Design the contradiction before choosing the construction

It suffices to produce polynomials \(P_K\in\mathbb Z[X]\), for arbitrarily large \(K\), such that

\[
\deg P_K=O(K),\qquad
0<P_K(\zeta(5))\le e^{-cK^2}
\]

for a fixed \(c>0\). If \(\zeta(5)=a/b\), then
\(b^{\deg P_K}P_K(a/b)\) is a positive integer. Its logarithm is at most
\(O(K)\log b-cK^2\), which tends to minus infinity. This contradicts the least positive integer \(1\).

There are three separate obligations: coefficientwise integrality, positivity, and sufficient decay. Integer-valuedness on integers would not suffice; a small value might be zero; and a denominator bound without a size estimate would not suffice. The determinant construction gives each obligation a distinct proof.

Fix throughout

\[
40\mid K,\qquad
N=\alpha K,\quad h=\lambda K,\qquad
\alpha=\frac3{40},\quad\lambda=\frac{37}{40}.
\]

The final polynomial will be the primitive integer part of a degree-\(h\) rational determinant. A common denominator \(Q=40000\) appears in the arithmetic certificate, but a short application of Gauss's lemma removes this power from the final polynomial. Its degree remains at most \(h\).

## 2. The determinant has rational coefficients and a positive real meaning

Let

\[
D_m(t)=\prod_{j=1}^{m}(t+j^2),\qquad
H_j^{(5)}=\sum_{r=1}^{j}\frac1{r^5}.
\]

Define a rational functional by polynomial division and simple partial fractions:

\[
\mu_X(t^e)=
\frac{(-1)^eB_{2e+2}(2e+3)(2e+4)(2e+5)}{24},
\qquad e\ge0,
\]

\[
\mu_X\!\left(\frac1{t+j^2}\right)
=j^4\bigl(X-H_j^{(5)}\bigr)-\frac14+\frac1{2j},
\qquad j\ge1.
\]

The Bernoulli convention is \(B_1=-1/2\). Set

\[
G_K(X)_{ij}=
\mu_X\!\left(\frac{D_N(t)^6t^{i+j}}{D_K(t)}\right),
\quad 0\le i,j<h,\qquad
\Delta_K(X)=\det G_K(X).
\tag{1}
\]

Every entry is affine in \(X\), so \(\Delta_K\in\mathbb Q[X]\) and
\(\deg\Delta_K\le h\). The scalar used on both sides of the proof is

\[
S_K=
\frac{(K!)^{2h}4^{h-1}}
{(N!)^{12h}\prod_{i=1}^{h-1}((2i)!)^2},
\qquad F_K=S_K\Delta_K.
\tag{2}
\]

Why these moments? At \(X=\zeta(5)\), they are integration against the positive weight

\[
w(y)=\frac{y^5}{12}\frac{d^4}{dy^4}\frac1{e^{2\pi y}-1}
=\frac{(2\pi)^4y^5}{12}\sum_{\ell\ge1}\ell^4e^{-2\pi\ell y}>0.
\tag{3}
\]

Termwise integration and Euler's even-zeta formula give the polynomial moments. For a simple pole, four integrations by parts and Hermite's Hurwitz-zeta integral give

\[
\int_0^\infty\frac{w(y)}{y^2+j^2}\,dy
=j^4\zeta(5,j)-\frac1{2j}-\frac14
=j^4(\zeta(5)-H_j^{(5)})+\frac1{2j}-\frac14.
\tag{4}
\]

The boundary terms vanish; the second equality uses
\(\zeta(5,j)=\zeta(5)-H_j^{(5)}+j^{-5}\). Thus the sign of the last term is positive, as in the rational definition.

For every nonzero polynomial \(r\) of degree below \(h\),

\[
\int_0^\infty
\frac{D_N(y^2)^6r(y^2)^2}{D_K(y^2)}w(y)\,dy>0.
\]

Consequently \(G_K(\zeta(5))\) is positive definite and
\(\Delta_K(\zeta(5))>0\). This proves nonvanishing before any asymptotic estimate. The [full assembly audit, Section 3](fauzan-simplified-proof-audit.md#3-the-functional-gram-matrix-and-degree-agree-exactly), checks the moment calculation and also proves exact degree \(h\); the degree upper bound is enough here.

The scalar \(S_K\) is purposeful. An integer-valued Newton basis makes the arithmetic manageable, and its triangular change-of-basis determinant produces precisely (2). Retaining this same scalar in the real asymptotics is essential.

## 3. The real estimate: an equilibrium measure instead of many trial measures

The reviewed analytic input is

\[
\boxed{\quad
\limsup_{\substack{K\to\infty\\40\mid K}}
\frac{\log F_K(\zeta(5))}{K^2}
\le-\frac{69}{50}.
\quad}
\tag{5}
\]

Here is the construction and why it applies. Andréief's determinant identity writes (1) as an \(h\)-fold integral with a squared Vandermonde. After \(y=K\sqrt t\), the large exponential part is governed by

\[
V(t)=2\pi\sqrt t+
\int_0^1\log(t+u^2)\,du
-6\int_0^\alpha\log(t+u^2)\,du.
\tag{6}
\]

For a measure \(\rho\) of mass \(\lambda\), write
\(U_\rho(t)=\int\log|t-x|\,d\rho(x)\) and
\(I(\rho)=\int U_\rho\,d\rho\). We seek

\[
2U_\rho(t)-V(t)\le m_0\quad(t\ge0),
\]

with equality on its support. Then the zero-mass energy inequality
\(I(\sigma-\rho)\le0\), for measures of equal mass, bounds every competing configuration by
\(\lambda m_0-I(\rho)\).

The field derivative has a particularly simple signed charge density:

\[
V'(t)=\int_0^\infty\frac{\chi(u)}{t+u^2}\,du,\qquad
\chi(u)=
\begin{cases}
-3,&0<u<\alpha,\\
3,&\alpha<u<1,\\
2,&u>1.
\end{cases}
\tag{7}
\]

Choose \(0<a<b\), let \(\omega\) be unit arcsine measure on \([a,b]\), and put

\[
C_u=\sqrt{(u^2+a)(u^2+b)},\qquad
d\omega_u(x)=\frac{C_u}{x+u^2}\,d\omega(x),
\]

\[
\rho=\lambda\omega+
\frac12\int_0^\infty\chi(u)(\omega_u-\omega)\,du.
\tag{8}
\]

Each \(\omega_u\) has mass one. The grouped difference in (8) is necessary: the signed charge itself has infinite mass, whereas this difference is integrable.

Choose the endpoints by

\[
\int_0^\infty\frac{\chi(u)}{C_u}\,du=0,\qquad
\lambda+\frac12\int_0^\infty\chi(u)
\left(\sqrt{\frac{u^2+a}{u^2+b}}-1\right)du=0.
\tag{9}
\]

These conditions remove both endpoint singularities of the provisional density. Define, for \(x>0\),

\[
K_\rho(x)=\int_0^\infty
\frac{\chi(u)}{(x+u^2)C_u}\,du.
\]

Subtracting \(1/(x+\alpha^2)\) inside this integral, using the first equation in (9), shows \(K_\rho(x)<0\): on both sides of \(\alpha\), the charge and the difference have opposite signs. Therefore

\[
d\rho(x)=
-\frac{\sqrt{(x-a)(b-x)}}{2\pi}K_\rho(x)\,dx>0.
\tag{10}
\]

The same sign, through the arcsine Stieltjes transform, shows that
\(2U_\rho-V\) increases up to \([a,b]\), is constant on it, and decreases afterwards. Its value at \(0\) follows by continuity; the potentially divergent expression \(K_\rho(0)\) is not used.

The [endpoint certificate](equilibrium-endpoint-certificate.md) proves existence and uniqueness in the rational rectangle

\[
a\in[0.00007119508,0.00007119517],\qquad
b\in[0.75394957,0.75394967].
\tag{11}
\]

The [rational energy certificate](rational-energy-certificate.md), with [independent review](rational-energy-review.md), proves

\[
\lambda m_0-I(\rho)+C_*<
-\frac{2771}{2000}<-\frac{69}{50},
\tag{12}
\]

where

\[
C_*=-2\lambda+12\alpha\lambda(1-\log\alpha)
+3\lambda^2-2\lambda^2\log(2\lambda).
\]

These are outward interval bounds with explicit quadrature remainders, not values from a numerical root finder.

There is one important bridge between this continuous statement and the finite determinant. Point masses have infinite self-energy. The [reviewed transfer argument](fauzan-simplified-proof-audit.md#7-from-a-continuous-measure-to-the-finite-determinant) replaces them by small circles, producing an \(O(K\log K)\) error. It also retains an integrable tail factor on the unbounded domain; a supremum alone would not control an integral over infinite volume. Uniform Riemann-sum estimates and Stirling for the same \(S_K\) then cancel the \(K^2\log K\) terms and prove (5).

The [equilibrium guide](fauzan-equilibrium-guide.md) gives the full analytic explanation. No prime-distribution theorem occurs in this part.

## 4. State precisely the local arithmetic that remains necessary

For a rational polynomial \(A\), let \(v_p^{\mathrm G}(A)\) be the minimum \(p\)-adic valuation of its coefficients. It is a coefficientwise notion and satisfies
\(v_p^{\mathrm G}(AB)=v_p^{\mathrm G}(A)+v_p^{\mathrm G}(B)\).

We retain the following local inputs, reconstructed from the primary source in the [local arithmetic review](fauzan-local-arithmetic-review.md) and checked for uniformity in the [PNT-free normalizer review](pnt-free-normalizer-independent-review.md).

**Small-prime lemma.** For every prime \(p\le5K\), put
\(\ell_p=\lfloor\log_p(5K)\rfloor\). Then

\[
v_p^{\mathrm G}(F_K)\ge
-6h\ell_p-hv_p(24)\ge-9K\ell_p.
\tag{13}
\]

The mechanism is an integer-valued binomial basis, a root-ball estimate, and
\(\mathcal L\binom{x}{j}=(-1)^j/(j+1)\) for the Bernoulli functional
\(\mathcal L(x^j)=B_j\). It includes the primes \(2,3\) and all higher prime powers.

**Large-prime profile lemma.** For each fixed \(M\ge40\), there is a constant \(E_M\) such that, for sufficiently large admissible \(K\) and every prime \(K/M<p\le5K\),

\[
v_p^{\mathrm G}(\Delta_K)\ge-p\,T(K/p)-E_M.
\tag{14}
\]

For \(p>K\), the stronger assertion \(v_p^{\mathrm G}(\Delta_K)\ge0\) holds. Primes \(p>5K\) also contribute no denominator.

Here is the exact profile in (14). Put \(J(v)=\lfloor v\rfloor v-\lfloor v\rfloor(\lfloor v\rfloor+1)/2\), and

\[
\begin{aligned}
\Gamma(x)={}&-\frac{197}{40}x
+\frac12J(23x/10)-\frac12J(2x)+\frac12J(3x/10)\\
&-3J(37x/40)+3J(43x/40)-9J(3x/20).
\end{aligned}
\tag{15}
\]

Write \(u_+=\max(u,0)\). For \(1/3<y<1\), set

\[
R_0(y)=
\begin{cases}
8-9y-8\alpha-5\min(\alpha,1-2y)-5(1+\alpha-3y)_+,
&y<1/2,\\
7(1-y)-6\min(\alpha,1-y)-6(1+\alpha-2y)_+
+(1+4\alpha-2y)_+,&y\ge1/2,
\end{cases}
\]

\[
d(y)=
\begin{cases}
(1+4\alpha-3y)_+-(1+\alpha-3y)_+,&y<1/2,\\
0,&y\ge1/2.
\end{cases}
\]

Then

\[
T(x)=
\begin{cases}
0,&0\le x\le1,\\
x[R_0(1/x)-d(1/x)],&1<x<3,\\
-\Gamma(x),&x>3.
\end{cases}
\tag{16}
\]

At \(x=3\), take the larger one-sided value.

The inner range uses a balanced Chinese-remainder interpolation basis. Assigning rows to residue classes forces zeros at the other classes; the determinant adds the resulting row valuations. Counting those rows gives \(p\Gamma(K/p)+O_M(1)\). The error is uniform, including points where the preferred allocation changes.

The outer range has fewer colliding poles. The matrix becomes \(A+p^{-1}L\), with \(L\) integral and of actual rank at most \(r\). Its assigned weights are nonpositive half-integers, and \(v_p^{\mathrm G}(A_{ij})\ge w_i+w_j\). Two elementary bounds give the refined estimate:

\[
v_p^{\mathrm G}\det(A+p^{-1}L)
\ge2\sum_iw_i-\min(r,z),
\tag{17}
\]

Here \(z\) counts zero weights. First, expand using complementary minors: terms with \(k>r\) vanish, and deleting nonpositive weights cannot hurt the remaining estimate. This gives \(2\sum_iw_i-r\).

Second, replace each zero weight by \(-1/2\), leaving all other weights unchanged. The new weights \(w_i'\) satisfy \(w_i'\le-1/2\), so both \(A\) and \(p^{-1}L\) satisfy the entry bounds \(w_i'+w_j'\): for the correction, use \(v_p(L_{ij})-1\ge-1\ge w_i'+w_j'\). Every determinant term then has valuation at least
\(2\sum_iw_i'=2\sum_iw_i-z\).
Taking the better of these two bounds proves (17), with no optimization over sorted minors. The [independent two-bound proof](refined-rank-two-bounds.md) explains the saving \(d\) in (16).

This route **retains (17)**. Dropping its zero-weight saving gives the separate plain-rank proof; the present 25-term certificate must not be applied to that different, larger target without a new verification.

## 5. Discover the right arithmetic object: covariance, Bernoulli, superfactorial

Why should any short factorial product match the complicated allocation?
The [floor-covariance derivation](floor-covariance-barnes-design.md), with [independent symbolic review](floor-covariance-independent-review.md), starts with

\[
\ell_u(z)=\lfloor u-z\rfloor+\lfloor u+z\rfloor+1,
\qquad A_u(z)=\ell_u(z)-2u,\quad 0<z<1/2.
\]

If \(r=\{u\}\), \(\rho_u=\min(r,1-r)\), and
\(\varepsilon_u=1\) for \(r\le1/2\), \(-1\) otherwise, then almost everywhere

\[
A_u(z)=\varepsilon_u
\left(\mathbf1_{z<\rho_u}-2\rho_u\right).
\tag{18}
\]

Thus the fluctuations of a floor count are centered interval indicators. Multiplying them gives

\[
\begin{aligned}
\int_0^{1/2}A_u(z)A_v(z)\,dz
&=\varepsilon_u\varepsilon_v
\bigl(\min(\rho_u,\rho_v)-2\rho_u\rho_v\bigr)\\
&=\frac{B_2(\{u-v\})-B_2(\{u+v\})}{2}.
\end{aligned}
\tag{19}
\]

The middle expression is the covariance kernel of a Brownian bridge on an interval of length \(1/2\). No probability theorem is needed: the equality follows from the overlap length of two intervals. The probabilistic name explains the shape, while the calculation supplies the proof.

Here \(B_2(t)=t^2-t+1/6\). This ordinary quadratic polynomial packages the overlap cases.

Expanding the row-allocation expression into these means and covariances, and using

\[
B_2(\{v\})=v^2-v+\frac16-2J(v),
\tag{20}
\]

gives (15). In particular \(197/40=5-\alpha\) is derived, not fitted.

Now define the superfactorial

\[
\mathcal B(n)=\prod_{j=0}^{n-1}j!.
\]

Count the contribution of each prime power to each factorial:

\[
\boxed{\quad
v_p(\mathcal B(n))
=\sum_{r\ge1}\sum_{k\ge1}(n-kp^r)_+
=\sum_{r\ge1}p^rJ(n/p^r).
\quad}
\tag{21}
\]

The product ends at \((n-1)!\); ending at \(n!\) would change the formula. Equation (21) is the discovery step: the floor profile is already the valuation of a concrete arithmetic object.

We should therefore search for a pointwise majorant
\(\sum a_jJ(c_jx)\ge T(x)\), rather than sum \(T(K/p)\) over primes and estimate their density. The Brownian-bridge interpretation motivates this search; it is not an additional proof assumption.

## 6. A finite rational certificate proves the majorant for every real \(x\ge0\)

Use the 25 nonzero pairs \((c_j,a_j)\) in the exact
[coefficient file](../computations/sparse_barnes_certificate_100.json), and put
\(H_*(x)=\sum_j a_jJ(c_jx)\).
The file also displays one zero coefficient, which has no effect. It is part of the explicit specification of the construction.

The [standard-library checker](../computations/check_sparse_barnes_certificate.py) proves

\[
H_*(x)\ge T(x)\quad(x\ge0),\qquad
\frac12\sum_j a_jc_j^2=-G_2,\qquad
G_2=10\alpha\lambda=\frac{111}{160},
\tag{22}
\]

\[
Q=\mathop{\rm lcm}_j\operatorname{den}(a_j)=40000,\qquad
W_1=\sum_j|a_j|c_j=\frac{2053577}{40000}.
\tag{23}
\]

Every \(c_j\) is a positive multiple of \(1/40\), and every nonzero frequency is at most \(9/8\). Thus \(c_jK\) is an integer whenever \(40\mid K\).

The check is not a mesh experiment. Both profiles are affine between their rational breakpoints. The script verifies the inequality at all 561 necessary points on \([0,43]\), including transitions. For the infinite tail, exact quadratic cancellation gives

\[
[H_*+\Gamma](x+40)-[H_*+\Gamma](x)
=40\beta,\qquad
\beta=\frac{9377}{80000}>0.
\tag{24}
\]

Every \(x\ge3\) reduces to \([3,43]\) by subtracting a nonnegative multiple of \(40\); (24) then propagates the inequality. This is a finite proof of an all-\(x\) assertion.

A floating-point optimization suggested the coefficients. It is not a logical dependency: the rational list, exact breakpoint comparisons, and exact periodic identity suffice to verify the result.

## 7. Construct an integral intermediate, including both repairs

Set

\[
\mathcal T_K=\prod_j\mathcal B(c_jK)^{Qa_j}\in\mathbb Q_{>0},
\qquad
M=10^7,\quad C=148=\lceil2W_1+45\rceil,
\]

\[
H_K=\lfloor K/M\rfloor,\quad
D_m=\operatorname{lcm}(1,\ldots,m),\quad
r_K=\lceil\sqrt K\rceil.
\]

Define the repair factor and powered intermediate

\[
E_K=D_{H_K}^{\,3CK}((5K)!)^{r_K},
\]

\[
\boxed{\quad
\mathcal W_K(X)=E_K^Q
\mathcal T_K\,\Delta_K(X)^Q.
\quad}
\tag{25}
\]

The \(D_m\) here is an integer lcm, not the polynomial \(D_m(t)\) in Section 2.
All superfactorial exponents \(Qa_j\) are integers; negative exponents produce rational numbers. We do not take fractional powers in a \(p\)-adic field. This intermediate makes the valuation verification convenient; Section 9 extracts an integer polynomial of the original degree.

**Large primes.** Fix \(M\) before allowing \(K\) to grow. If \(p>H_K\), then \(p>K/M\), and eventually \(p^2>5K\). Equation (21) therefore becomes exactly

\[
\frac1Qv_p(\mathcal T_K)=pH_*(K/p).
\]

By (14) and (22),
\(v_p^{\mathrm G}(\mathcal T_K\Delta_K^Q)\ge-QE_M\).
For \(p\le5K\), its valuation in \((5K)!\) is at least one. Eventually \(r_K\ge E_M\), so the factorial repair clears this bounded deficit at every large prime simultaneously. For \(p>5K\), all factors are already integral.

**Small primes, including higher powers.** For \(0\le m\le5K\), (21) gives

\[
\left|v_p(\mathcal B(m))-\frac{m^2}{2(p-1)}\right|
\le2m\ell_p.
\tag{26}
\]

Indeed, write \(J(v)=v^2/2-v/2+\{v\}(1-\{v\})/2\), sum over \(p^r\le m\), and bound the remaining geometric tail of the quadratic term. Elementary factorial digit-sum bounds applied to (2) also give

\[
\left|v_p(S_K)+\frac{G_2K^2}{p-1}\right|
\le36K\ell_p.
\tag{27}
\]

The quadratic terms in (26)–(27) cancel exactly by (22). Combining these estimates with (13) yields

\[
v_p^{\mathrm G}(\mathcal T_K\Delta_K^Q)
\ge-Q(2W_1+45)K\ell_p
\ge-QCK\ell_p.
\tag{28}
\]

For \(K\ge200M^2\), one has \(5K\le H_K^2\). If \(p\le H_K\), put
\(u=\lfloor\log_pH_K\rfloor\ge1\). Then
\(\ell_p\le2u+1\le3u\), and \(v_p(D_{H_K})=u\).
The lcm factor in (25) clears (28).

Every coefficient now has nonnegative valuation at every prime, proving

\[
\mathcal W_K\in\mathbb Z[X]
\quad\text{for all sufficiently large admissible }K.
\tag{29}
\]

Gauss valuation is multiplicative, so raising the determinant to \(Q\) is legitimate coefficientwise. Positivity follows from Section 2 and the positive scalar, and \(\deg\mathcal W_K\le Qh\). This degree belongs to the arithmetic intermediate, not to the final primitive polynomial.

The threshold \(K\ge200M^2=2\cdot10^{16}\) is sufficient for the displayed elementary size comparisons, but is not an explicit starting index for the whole theorem: the uniform constant \(E_M\) and subsequent asymptotics can require a larger \(K\). The unbounded repair \(r_K\) avoids needing its numerical value.

## 8. Price the actual product by Stirling, without a prime sum

Summing Stirling's formula gives

\[
\log\mathcal B(m)=\frac{m^2}{2}\log m-\frac{3m^2}{4}
+O(m\log m).
\]

For the same scalar (2),

\[
\log S_K=-G_2K^2\log K+C_*K^2+O(K\log K).
\]

The exact balance in (22) cancels the \(K^2\log K\) term, leaving

\[
\frac1Q\log\frac{\mathcal T_K}{S_K^Q}
=A K^2+O(K\log K),
\qquad
A=\frac12\sum_j a_jc_j^2\log c_j+\frac32G_2-C_*.
\tag{30}
\]

The factorial repair costs only
\(Qr_K\log(5K)!=O(QK^{3/2}\log K)=o(QK^2)\).

For the lcm repair, the elementary central-binomial argument proves

\[
\log D_m\le m\log4+O((\log m)^2).
\tag{31}
\]

To see this, each prime power in \((\lfloor m/2\rfloor,m]\) contributes its prime to \(\binom m{\lfloor m/2\rfloor}\), except possibly the single value \((m+1)/2\). Thus
\(\log D_m\le\log D_{\lfloor m/2\rfloor}+m\log2+\log m\);
iterate. No asymptotic statement about primes is being used.

Writing \(\mathcal W_K\) as the scalar in (25) divided by \(S_K^Q\), times \(F_K^Q\), its arithmetic cost per \(QK^2\) is therefore at most

\[
A+\frac{3C\log4}{M}.
\]

The sparse checker bounds every logarithm by a rational atanh series with an explicit positive remainder, after dyadic range reduction. It certifies the stronger bound

\[
A+\frac{3C\log4}{M}
<\frac{6893}{5000}=1.3786
<\frac{1379}{1000}=1.379.
\tag{32}
\]

The diagnostic value is about \(1.3785819071\); decimals do not decide any inequality. Combining the convenient bound \(1.379\) with (5) gives

\[
\limsup\frac{\log\mathcal W_K(\zeta(5))}{QK^2}
<\frac{1379}{1000}-\frac{69}{50}
=-\frac1{1000}.
\tag{33}
\]

In particular, for all sufficiently large admissible \(K\),

\[
0<\mathcal W_K(\zeta(5))<
\exp\!\left(-\frac{QK^2}{2000}\right).
\tag{34}
\]

We deliberately use a weaker decay constant to leave visible room for every lower-order term.

## 9. Remove the power by Gauss's lemma, then contradict rationality

There is a useful general principle. Write a nonzero rational polynomial as

\[
F=gP,\qquad g\in\mathbb Q_{>0},\quad
P\in\mathbb Z[X]\ \text{primitive}.
\]

Here primitive means that the greatest common divisor of the integer coefficients is \(1\); the sign remains in \(P\). Suppose \(c>0\) is rational and \(cF^Q\in\mathbb Z[X]\). Gauss's lemma makes \(P^Q\) primitive, so its rational scalar \(cg^Q\) must be a positive integer. In particular \(cg^Q\ge1\), and at any real \(\theta\) with \(F(\theta)>0\),

\[
0<P(\theta)=\frac{F(\theta)}g
\le\bigl(cF(\theta)^Q\bigr)^{1/Q}.
\tag{35}
\]

Apply this to \(F=\Delta_K\) and \(c=E_K^Q\mathcal T_K\). Its primitive part \(P_K\) has degree at most \(h\), and (34) gives

\[
\boxed{\quad P_K\in\mathbb Z[X],\qquad
\deg P_K\le h,\qquad
0<P_K(\zeta(5))<e^{-K^2/2000}.\quad}
\tag{36}
\]

Thus clearing the rational certificate's denominator need not multiply the final degree. The [independent rounding and primitive-part review](prime-exponent-rounding-independent-review.md) verifies this step, including the positivity of the content scalar.

Suppose now that \(\zeta(5)=a/b\), with integers \(a\) and \(b>0\). By (36),

\[
b^{h}P_K(a/b)\in\mathbb Z_{>0}.
\]

It is smaller than
\(\exp(\lambda K\log b-K^2/2000)\), which tends to zero. This proves the contradiction. No coefficient-height estimate is needed, \(b\) is fixed, and positivity excludes a zero polynomial value.

There is also an explicit normalization without taking a determinant power. Put

\[
b_p=\sum_j a_jv_p(\mathcal B(c_jK)),\qquad
m_K^-=\prod_p p^{\lfloor b_p\rfloor}.
\]

Only finitely many exponents are nonzero. The fully repaired proof in Section 7 says
\(v_p(E_K)+v_p^{\mathrm G}(\Delta_K)+b_p\ge0\).
The first two terms have an integer sum \(z_p\). Consequently
\(z_p\ge\lceil-b_p\rceil=-\lfloor b_p\rfloor\), even when \(b_p<0\). Hence

\[
\widehat P_K=E_Km_K^-\Delta_K\in\mathbb Z[X],
\qquad \deg\widehat P_K\le h.
\]

Rounding down adds no upper cost:

\[
\log m_K^-
\le\sum_p b_p\log p
=\sum_j a_j\log\mathcal B(c_jK).
\tag{37}
\]

This direct formula gives the same quadratic estimate. It requires the repairs to have been included first; rounding an uncorrected profile is not sufficient. The primitive-part proof (35) is the shortest route to the contradiction, while (37) gives a finite explicit normalizer.

The complete chain is:

~~~mermaid
flowchart TD
  A["Rational functional and Fauzan determinant"] --> B["Positive moment identity"]
  A --> C["Reviewed local interpolation and refined rank lemmas"]
  B --> D["Equilibrium measure and global potential inequality"]
  D --> E["Exact endpoint and energy certificates"]
  E --> F["Finite determinant transfer: real rate at most -1.38"]
  C --> G["Required prime valuation profile T"]
  G --> H["25-term rational majorant and periodic propagation"]
  H --> I["Superfactorial product plus both prime repairs"]
  I --> J["Integral intermediate and arithmetic rate below 1.379"]
  J --> L["Gauss primitive-part lemma: degree h"]
  F --> K["Positive integer under rationality tends to zero"]
  L --> K
~~~

The [PNT-free normalizer audit](pnt-free-normalizer-independent-review.md) checks (21) and (25)–(32), including all prime powers, rational exponents, the constant \(148\), and the uniform large-prime repair. The [rounding review](prime-exponent-rounding-independent-review.md) checks (35)–(37). The [earlier end-to-end audit](fauzan-simplified-proof-audit.md) supplies the shared moment, local-functional, and real-transfer proofs; its prime-averaging stage is replaced here by Sections 5–8. Its plain-rank normalizer is not silently reused.

## 10. Reproduce, compare routes, and keep the limits explicit

From the repository root, the main certificates use only Python's standard library and the local interval primitives:

    python codex/computations/check_floor_covariance_barnes.py
    python codex/computations/check_sparse_barnes_certificate.py
    python codex/computations/check_pnt_free_normalizer_review.py
    python codex/computations/check_prime_exponent_rounding_review.py
    python codex/computations/certify_equilibrium_endpoints.py
    python codex/computations/check_rational_energy_primitives.py
    python codex/computations/certify_equilibrium_energy_rational.py

The proofs explain why these finite checks imply the asserted inequalities; experiments at many parameter values alone would not. Standard Hermite-integral, Bernoulli, integer-valued-polynomial, and determinant identities remain mathematical inputs. The covariance identities are elementary; the entire proof has not thereby become elementary.

There are three distinct arithmetic routes:

| Route | What it simplifies | What it retains |
|---|---|---|
| This PNT-free construction | Removes prime-density averaging; retains degree \(h\) by primitive parts or floor rounding | Refined outer rank estimate, 25-term rational certificate, both prime repairs; \(40000\) is only an auxiliary denominator |
| [Plain-rank route](fauzan-equilibrium-guide.md#12-spend-a-little-analytic-margin-to-simplify-the-local-arithmetic-too) | Removes the zero-weight optimization and long integration tables | PNT, its own normalizer \(\widetilde m\), degree \(h\) |
| [Optional outer-kernel improvement](outer-gain-budget-update.md) | Gives a stronger arithmetic margin using a finite-field kernel | Additional reviewed local lemma, PNT, another explicitly changed normalizer |

The [Bernoulli period tower](bernoulli-period-tower.md), with [independent review](bernoulli-period-independent-review.md), offers a related structural lesson: the same Bernoulli algebra appears in real periods and discrete arithmetic. It is useful motivation, not a hidden dependency of this proof.

The common real determinant budget extends to every fixed odd index, as explained in the [equilibrium guide](fauzan-equilibrium-guide.md). The local functional exponents do not automatically extend with the same cost. This guide proves no all-odd irrationality theorem, no practical starting index, and no new irrationality measure. The remaining research objective is to design and prove the corresponding arithmetic normalization for higher odd values—not to infer it from the shared analytic estimate.
