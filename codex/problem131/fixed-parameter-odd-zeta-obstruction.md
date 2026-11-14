# A fixed numerator power and fixed aspect ratio cannot handle all odd zeta values

**Status: independently reviewed structural obstruction; see fixed-parameter-obstruction-independent-review.md.** This generalizes the exact monic-coefficient method of `original-outer-monic-obstruction.md`. The parent proposed the generalization; the proof below checks the numerator-power-one and zero-aspect-ratio endpoints and supplies an explicit rational threshold. It uses the prime number theorem. No historical novelty or general impossibility claim is made.

## 1. Precise theorem

Fix an integer r>=1 and a rational number 0<=alpha<1. Let K tend to infinity through multiples of the denominator of alpha, and put
\[
 N=\alpha K,\qquad h=K-N,\qquad D_j(t)=\prod_{u=1}^j(t+u^2).
\]
For each fixed odd s>=3 form the original, undepleted determinant
\[
 \Delta_{s,K}(X)=\det_{0\le i,j<h}
 \mu_{s,X}\left(\frac{D_N(t)^r t^{i+j}}{D_K(t)}\right).
\]
Here the Hermite functional has positive density
\[
 w_s(y)=\frac{2(2\pi)^{s-1}}{(s-1)!}y^s
 \sum_{n\ge1}n^{s-1}e^{-2\pi ny},
\]
polynomial moments
\[
 \mu_s(t^e)=\frac{(-1)^e B_{2e+2}}{s-1}
 \binom{2e+s}{s-2},
\]
and pole values obtained by replacing zeta(s) with X in
\[
 \mu_s((t+j^2)^{-1})
 =j^{s-1}\zeta(s)-j^{s-1}H_j^{(s)}-\frac1{s-1}+\frac1{2j}.
\]
Let P_(s,K) be a primitive integral representative of Delta, with sign chosen positive at zeta(s).

**Theorem. There is a finite explicit threshold S(r,alpha) such that for every fixed odd s>=S(r,alpha), there is c_(s,r,alpha)>0 with**
\[
 \boxed{P_{s,K}(\zeta(s))>\exp(c_{s,r,\alpha}K^2)
 \quad\text{for all sufficiently large }K.}
\]
Consequently no nonzero rational scalar normalization of this fixed family can supply small integral-polynomial values for all odd indices. Rational basis changes within its full h-dimensional polynomial space are also scalar changes of the determinant, so they are covered. Numerator powers, aspect ratios or constructions that change with s are **not** ruled out.

The quantifiers are important: s is fixed before K tends to infinity. No estimate uniform for s growing with K is asserted.

## 2. An explicit threshold involving only rational arithmetic

Choose any rational
\[
 \max(1/2,\alpha)<\beta<1.
\]
For a completely prescribed choice one may take
\(\beta=(1+\max(1/2,\alpha))/2\).
Define
\[
 \lambda=1-\alpha,\quad \kappa=1+(r-2)\alpha>0,
\]
\[
 C=\lambda(4r\alpha+9/4)+3\lambda^2,
\]
\[
 T=\max\left\{3,2r-1,1+\frac{2\kappa}{\beta},
                    \frac{2C}{(1-\beta)^2}\right\},\qquad
 S=2\left\lfloor\frac{T+1}{2}\right\rfloor+1.         \tag{1}
\]
This S is the smallest odd integer strictly greater than T. All quantities in (1) are rational before the final floor, so it is an explicit finite threshold without numerical transcendental comparisons. It is deliberately coarse; optimizing beta or using the exact real constant below improves it substantially.

For any fixed odd s>=S, one may take any positive c less than
\[
 \frac12\left[\frac{s(1-\beta)^2}{2}-C\right].       \tag{2}
\]
The bracket is strictly positive by construction.

## 3. The leading coefficient and monic normalization

After cancelling D_N from the denominator, put
\[
 D=D_K/D_N,\qquad W=D_N^{r-1}.
\]
The X-coefficient matrix factors through the h surviving nodes. Its determinant is
\[
 \ell_K=(-1)^{h(h-1)/2}
 \prod_{j=N+1}^K j^{s-1}D_N(-j^2)^{r-1}
 \in\mathbb Z\setminus\{0\}.                        \tag{3}
\]
This proves degree exactly h for every r>=1, including r=1. The positivity of the Gram matrix at zeta(s) holds for every r>=1 too, because D_N(y^2)>0 on the integration domain; evenness of r is not required.

Set F_K=Delta/ell_K, a monic rational polynomial, and let M_K be the least positive integer clearing all its coefficient denominators. Then M_KF_K is primitive, and
\[
 P_{s,K}(\zeta(s))=M_K|F_K(\zeta(s))|.              \tag{4}
\]
Every integral rational scalar multiple of Delta is a nonzero integer multiple of this primitive polynomial. Monic normalization prevents numerator content at other primes from offsetting the denominator lower bounds that follow.

## 4. A prime band with an exact coefficient valuation

Take primes beta K<p<=K. Eventually p>N, K<2p, p>s-1 and p^2>2K. Use the rational polynomial basis
\[
 f_j=D/(t+j^2),\qquad N<j\le K.
\]
As in the fixed-parameter example, its X-matrix is diagonal with entries
\[
 C_j=j^{s-1}W(-j^2)D'(-j^2).
\]
Its constant matrix B is a polynomial-part matrix plus the diagonal pole contributions. Every polynomial quotient has degree at most
\[
 E=(r-1)N+h-2=K+(r-2)N-2=\kappa K-2.               \tag{5}
\]
If E<0 the polynomial part is zero. Otherwise the Bernoulli cancellation lemma makes it p-integral whenever
\[
 p>1+(2E+2)/(s-1).
\]
The condition (s-1)beta>2kappa in (1) guarantees this eventually throughout the band.

For an ordinary nonzero squared-residue class, let ell<=4 be the number of original nodes, and delta<=2 the number removed among 1,...,N. At a surviving node,
\[
 v_p(C_j)=(r-1)\delta+(\ell-\delta-1)
 =\ell+(r-2)\delta-1
 \le B_r:=\max(3,2r-1).                             \tag{6}
\]
The bound also covers r=1: in that case the expression is ell-delta-1<=3. At the unique zero-class node j=p,
\[
 v_p(C_p)=s-1.                                     \tag{7}
\]
For s>B_r all rows of p^s C^(-1) times the polynomial-part matrix vanish modulo p. The harmonic constants give exactly the diagonal residues
\[
 0\ (j<p),\qquad -1/2\ (j=p),\qquad -1\ (j>p).
\]
The rational change of basis cancels from the monic determinant; it need not be unimodular. Hence
\[
 \overline{p^{sh}F_K(Z/p^s)}
 =Z^{p-N-1}(Z-1/2)(Z-1)^{K-p}.                     \tag{8}
\]
In particular, with b_p=K-p+1,
\[
 \boxed{v_p([X^{h-b_p}]F_K)=-s b_p,\qquad
        v_p(M_K)\ge s(K-p+1).}                     \tag{9}
\]
These are actual coefficient valuations. Other coefficients can require still larger denominators; no equality for the full Gauss valuation is needed.

Using the prime number theorem and partial summation only to aggregate (9),
\[
 \boxed{\liminf\frac{\log M_K}{K^2}
 \ge s\int_\beta^1(1-u)\,du
 =\frac{s(1-\beta)^2}{2}.}                         \tag{10}
\]
Fixed exceptional primes and the +1 term have lower-order total cost. The exact finite-prime theorem (9) is independent of the prime number theorem.

## 5. A finite real lower constant, independent of fixed s

Apply Andreief's identity to the positive Gram determinant and restrict every integration variable to y=Kx with 0<x<b, choosing b=1/2. The first exponential term gives
\(w_s(y)\ge c_s y^s e^{-2\pi y}\).
Jensen's inequality under independent uniform x variables gives a real lower bound. All logarithms, including the pair-distance singularities, are integrable.

Define
\[
 J_t=\frac1b\int_0^b\int_0^t\log(x^2+u^2)\,du\,dx
 =t\log(b^2+t^2)-3t+\frac{t^2}{b}\arctan(b/t)
                         +b\arctan(t/b),
\]
with J_0=0, and
\[
 I_\alpha=\int_\alpha^1\int_0^\alpha\log(x^2-u^2)\,du\,dx.
\]
The latter is finite and nonpositive. If desired it has the explicit form
\[
 I_\alpha=H(1+\alpha)-H(2\alpha)-H(1-\alpha)
                         -2\alpha(1-\alpha),
 \quad H(v)=\tfrac12v^2\log v-\tfrac14v^2,\ H(0)=0.
\]
The pair expectation at b=1/2 is
\(\mathbb E\log|x^2-y^2|=-3\).

The K-log powers must be compared before taking a monic ratio. From the rational weight, Vandermonde, density power and differentials they total
\[
 (2rN-2K)h+2h(h-1)+(s+1)h
 =2(r-1)Nh+(s-1)h.
\]
These are exactly the powers in (3). The remaining leading-coefficient constant is (r-1)I_alpha. The fixed-s density constants and remaining s-dependent terms contribute only O_s(K log K), hence vanish after division by K^2.

The resulting bound is
\[
 \boxed{\liminf\frac1{K^2}\log|F_K(\zeta(s))|
 \ge L_{r,\alpha}:=
 \lambda[rJ_\alpha-J_1-\pi/2]-3\lambda^2
                  -(r-1)I_\alpha.}                \tag{11}
\]
No equilibrium-measure existence or minimization is required.

For completeness, the limit of the averaged logarithmic Riemann sums is valid even at u=0: their averaged integrand is continuous there and has derivative bounded by pi/b. The leading-coefficient limit follows also from
\[
 |D_N(-j^2)|=(j+N)!/[j(j-N-1)!]
\]
and uniform Stirling error O(log K) per factorial. This includes alpha=0, where D_N=1 and I_alpha=J_alpha=0. At r=1 the entire (r-1) term disappears; there is no singular limiting argument.

A crude rational bound is enough for existence. Since log(x^2+u^2)>=2log x,
\[
 J_\alpha\ge-2\alpha(\log2+1)>-4\alpha
\]
when alpha>0; the weak inequality suffices at alpha=0. Also
\(J_1\le\log(5/4)<1/4\), \(\pi/2<2\), and \(I_\alpha\le0\). Therefore
\[
 \boxed{L_{r,\alpha}>-\lambda(4r\alpha+9/4)-3\lambda^2=-C.} \tag{12}
\]
The strictness remains valid at alpha=0 because the other two inequalities are strict.

## 6. Finish the obstruction and state its limits

Combining (4), (10) and (11),
\[
 \liminf\frac1{K^2}\log P_{s,K}(\zeta(s))
 \ge\frac{s(1-\beta)^2}{2}+L_{r,\alpha}
 >\frac{s(1-\beta)^2}{2}-C>0.
\]
The final inequality is precisely the last condition built into (1), and proves the theorem.

Thus a fixed power and fixed aspect ratio cannot supply one unchanged full-pole determinant recipe for all odd zeta values. This is a statement about **actual primitive values**, stronger than failure of a sufficient normalization budget. It says nothing against allowing r or alpha to depend on s, selecting a different subspace dimension, adding genuine polynomial combinations, changing the functional, or using a different proof mechanism.

## 7. Exact checks

`codex/computations/check_fixed_parameter_obstruction.py` verifies the local valuation and harmonic reductions at nine parameter choices. They include r=1, alpha=0, alpha=9/10, and six powers r=1,...,6 with alpha=3/5, where two removed nodes share a squared residue class. In the latter cases the ordinary leading-factor valuation attains 2r-1 for r>=2, checking the feature absent from the earlier small-alpha proof. Exact Bernoulli moments verify the required polynomial-part integrality in each case.

The same script constructs S and its positive rational margin for five parameter pairs. For example, the prescribed simple beta gives S(6,3/40)=203; this is a coarse existence threshold, not a replacement for the sharper independently reviewed s>=13 theorem for that particular pair. All checks pass. The unbounded theorem follows from the proof, not the finite cases.

## 8. Quantitative design consequence: numerator power must grow at least linearly

Fix a rational A with 0<=A<1, and consider any of these full-pole families with alpha<=A. The parameters r and alpha may now be chosen separately for each fixed odd s. Define
\[
 \beta_A=\frac{1+\max(1/2,A)}2,\qquad
 c_A=\frac{(1-\max(1/2,A))^2}{8}>0.
\]
Then every parameter choice satisfying
\[
 \boxed{c_A s>r+21/4}                               \tag{13}
\]
has exponentially growing primitive values as K tends to infinity.

Indeed 4alpha(1-alpha)<=1 implies
\[
 C=4r\alpha\lambda+(9/4)\lambda+3\lambda^2
 \le r+21/4.
\]
The prime-band cost for beta_A is c_A s. Moreover c_A<=1/32 and beta_A>=3/4, so (13) implies s>32r+168. This automatically gives s>max(3,2r-1) and
\((s-1)\beta_A>2\kappa\), since kappa<=r for r>=1. Every hypothesis of the preceding proof is therefore satisfied.

Consequently, a parameter choice that could yield small primitive values must obey the necessary condition
\[
 \boxed{r\ge c_A s-21/4.}                           \tag{14}
\]
For all aspect ratios alpha<=1/2 this specializes to
\[
 r\ge s/32-21/4.
\]
The constants are deliberately weak, but the **linear growth requirement is an actual consequence of coefficient denominators**, not a sufficient-budget heuristic. Thus keeping the aspect ratio bounded away from one forces numerator powers to grow at least linearly along any potentially successful sequence of arbitrarily large odd indices. Alternatively the aspect ratios may approach one; the present uniform statement does not exclude that possibility.

As throughout, each s and its selected parameters are fixed before K tends to infinity. This does not address simultaneous growth of s with K, altered subspace dimensions, or different determinant/function constructions.
