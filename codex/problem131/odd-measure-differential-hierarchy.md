# One positive measure generates every odd-zeta measure by Euler operators

Status: independently derived from the defining densities, with exact algebra
checks and independent positive-density quadrature. The integration-by-parts
endpoint and integer-normalizer arguments are proved below. Independent review
is still required. This is a structural connection and an arithmetic transfer
rule, not an irrationality theorem or a historical novelty claim.

The hierarchy was proposed by the root agent. The Stieltjes/pole-order and
sharp integer-normalizer consequences below were developed independently in
this investigation.

## 1. The two-integration-by-parts hierarchy

For odd integers s≥3, put
\[
 f(y)=\frac1{e^{2\pi y}-1},\quad
 w_s(y)=\frac{2y^s}{(s-1)!}f^{(s-1)}(y),\quad
 \mu_s(g)=\int_0^\infty g(y^2)w_s(y)\,dy.
\]
These positive weights and their Hermite/Hurwitz interpretation are established
in `odd-moment-generalization.md`. An admissible rational g has no pole on
[0,∞); polynomial growth at infinity is allowed. Such g and its derivatives
are bounded near zero, and all integrals below converge absolutely.

**Theorem 1.** With Θ=t∂ₜ,
\[
 \boxed{\mu_{s+2}(g)=\mu_s(T_sg),\qquad
 T_s=\frac{(2\Theta+s+1)(2\Theta+s+2)}{s(s+1)}.}       \tag{1}
\]
Equivalently,
\[
 T_sg=\frac{4t^2g''+(4s+10)tg'+(s+1)(s+2)g}{s(s+1)}.
\]

Proof: integrate the defining integral for μₛ₊₂ by parts twice, moving the
two outer derivatives of f onto h(y)=yˢ⁺²g(y²). Direct differentiation gives
\[
 h''(y)=y^s\{4t^2g''+(4s+10)tg'+(s+1)(s+2)g\}_{t=y^2}.
\]
The factorial ratio is (s−1)!/(s+1)!=1/[s(s+1)]. At zero,
\(f^{(k)}(y)=O(y^{-k-1})\), h=O(yˢ⁺²), and h′=O(yˢ⁺¹). Both boundary
products h f⁽ˢ⁾ and h′ f⁽ˢ⁻¹⁾ are O(y). The post-integration integrand is
bounded near zero since h″f⁽ˢ⁻¹⁾=O(1). At infinity exponential decay beats
all these rational/polynomial factors. These estimates justify ordinary IBP
on [ε,R], followed by the improper limits.

All the operators are polynomials in the same Θ and commute. For s=2m+3,
\[
 \boxed{\mu_s(g)=\mu_3(P_m(\Theta)g),\qquad
 P_m(\Theta)=\frac2{(2m+2)!}
                  \prod_{j=4}^{2m+3}(2\Theta+j).}       \tag{2}
\]
For m=0 the product is empty and P₀=1. In particular these are not unrelated
positive moment models: they are adjoint differential transforms of one
another.

## 2. A sharp integral normalization of the differential transform

The factorial expression hides a much smaller exact uniform integer constant.
Set
\[
 N_s=\frac{(s-1)(s-2)}2=(m+1)(2m+1).
\]
Then
\[
 \boxed{P_m(\Theta)=\frac1{N_s}
                    \binom{2\Theta+s}{s-3}.}           \tag{3}
\]
The binomial here is its polynomial definition. On tᵉ it acts diagonally:
\[
 P_m(\Theta)t^e=\frac{\binom{2e+s}{s-3}}{N_s}t^e.
\]

**Theorem 2.** Nₛ is the least positive integer N such that
\(NP_m(\Theta)\mathbb Z[t]\subseteq\mathbb Z[t]\).
Moreover, for any integer polynomial D,
\[
 N_sP_m(\Theta)\mathbb Z[t,D^{-1}]
                  \subseteq\mathbb Z[t,D^{-1}].         \tag{4}
\]
Thus the raising operator has a fixed exact integral normalization independent
of the number of poles or polynomial degree.

Proof of sufficiency and the localization statement: write k=s−3=2m. The
integer-valued polynomial in e has a Newton expansion
\[
 \binom{2e+s}k=\sum_{r=0}^k b_{s,r}\binom er,
\quad
 b_{s,r}=[u^k](1+u)^s(2u+u^2)^r\in\mathbb Z_{\ge0}.      \tag{5}
\]
Indeed expand (1+u)²ᵉ=(1+2u+u²)ᵉ and extract the coefficient of uᵏ.
Since \(\binom\Theta r=t^r\partial_t^r/r!\),
\[
 \boxed{N_sP_m(\Theta)
        =\sum_{r=0}^{s-3}b_{s,r}\frac{t^r}{r!}\partial_t^r.} \tag{6}
\]
Divided derivatives preserve the integer localization: their values are the
coefficients of hʳ in g(t+h), and D(t+h)⁻¹ is a power series in h with
coefficients in \(\mathbb Z[t,D(t)^{-1}]\). This proves (4) without paying a
factorial denominator for each derivative.

For minimality, the integers \(\binom{2e+s}k\), e≥0, have gcd1. For each odd
prime p, choose L with pᴸ>k and choose a sufficiently large odd n≥s with
n≡k mod pᴸ. Lucas's congruence gives binom(n,k)≡1 mod p; such n has the form
2e+s. For p=2 choose n=2ᴸ−1≥s, k<2ᴸ; all its binary digits are1, so again
binom(n,k) is odd. Therefore for every p dividing Nₛ there is an e for which
Pₘ(e) has exactly the full denominator p^vₚ(Nₛ). N must contain every such
prime power. This proves sharpness. No prime asymptotics are involved.

For example N₃=1, N₅=6, N₇=15, N₉=28. The first transform does **not**
preserve integer polynomials without its normalizer: P₁(0)=5/3 and
P₁(1)=7/2, forcing both3 and2.

## 3. The Stieltjes transform is raised by derivatives

Define the shifted Stieltjes transform
\[
 F_s(z)=\mu_s\!\left(\frac1{t+z}\right),
 \qquad z\in\mathbb C\setminus(-\infty,0].
\]
The shift z is positive at the integer-square samples; the rational test
function has its pole at t=−z. The finite measure and locally uniform kernel
bounds justify differentiation on this slit domain. Homogeneity gives
\[
 (\Theta_t+\Theta_z)(t+z)^{-1}=-(t+z)^{-1}.
\]
Consequently (1) gives the clean Stieltjes recurrence
\[
 \boxed{F_{s+2}(z)=
 \frac{(2\Theta_z-s+1)(2\Theta_z-s)}{s(s+1)}F_s(z).}     \tag{7}
\]
In ordinary derivatives this is
\[
 F_{s+2}=
 \frac{4z^2F_s''+(6-4s)zF_s'+s(s-1)F_s}{s(s+1)}.
\]
In particular,
\[
 \boxed{F_5(z)=\tfrac12F_3(z)-\tfrac z2F_3'(z)
                         +\tfrac{z^2}3F_3''(z).}        \tag{8}
\]
Equivalently, without differentiating an integral sign,
\[
 \boxed{F_5(z)=\mu_3\!\left(
 \frac{1/2}{t+z}+\frac{z/2}{(t+z)^2}
                     +\frac{2z^2/3}{(t+z)^3}\right).}  \tag{9}
\]
This is an exact constructive bridge from ζ5's positive measure to ζ3's.
The bridge changes the rational kernel from simple poles to higher poles.

## 4. Exact pole-order theorem and coefficient control

**Theorem 3.** For z≠0 and s=2m+3,
\[
 \boxed{P_m(\Theta_t)\frac1{t+z}
  =\frac1{N_s}\sum_{r=0}^{2m}
                    c_{m,r}\frac{z^r}{(t+z)^{r+1}},
 \qquad c_{m,r}\in\mathbb Z,}                          \tag{10}
\]
with
\[
 c_{m,0}=2m+1,\qquad c_{m,2m}=4^m.
\]
Thus the pole order is **exactly s−2**, with highest-pole coefficient
\[
 \frac{4^m z^{2m}}{N_s}\ne0.                           \tag{11}
\]

Proof: apply (6) to the kernel. Its rth divided derivative is
(−1)ʳ/(t+z)ʳ⁺¹. Replace tʳ by ((t+z)−z)ʳ and collect powers. This proves
integrality of all c. The top coefficient is bₛ,₂ₘ=2²ᵐ=4ᵐ; the constant
coefficient is the value at Θ=−1,
\(\binom{s-2}{s-3}=s-2=2m+1\).

The first coefficient polynomials \(\sum c_{m,r}h^r\) are
\[
 3+3h+4h^2\quad(s=5),
\]
\[
 5+5h+5h^2+4h^3+16h^4\quad(s=7).
\]
Starting at s=9 some coefficients in this particular pole expansion are
negative; positivity of the original measure should not be inferred from
coefficient signs in every transformed representation.

If g has rational partial fractions
\[
 g=A(t)+\sum_{j=1}^K\frac{a_j}{t+j^2},
\]
then NₛPₘg has the same pole locations and orders at most s−2. At a pole
with aⱼ≠0 its order is exactly s−2 and its leading coefficient is
\(4^m j^{4m}a_j\). If M clears the coefficients of A and all aⱼ, then
M Nₛ clears **all** coefficients in the resulting higher-pole partial fractions.
There is no additional K-dependent coefficient-denominator cost from the
Euler transform itself. A common-denominator expression may use D_K^(s−2),
but its structured numerator cancels the apparent extra pole-separation
costs in the partial fractions.

This is a useful integrality bridge, but not a normalization theorem for
zeta-valued determinants. The values of the higher-pole moments themselves
have a different arithmetic content, as the next section makes explicit.

## 5. Why the derivatives contain ζ5 rather than just more ζ3 samples

The classical Hermite formula gives, for every real a>0,
\[
 \boxed{F_s(a^2)=a^{s-1}\zeta(s,a)-\frac1{s-1}-\frac1{2a}.} \tag{12}
\]
The hierarchy independently agrees with this formula: with Dₐ=a∂ₐ,
\[
 F_{s+2}(a^2)=
 \frac{(D_a-s+1)(D_a-s)}{s(s+1)}F_s(a^2).
\]
The two factors remove the power aˢ⁻¹ before differentiating the Hurwitz zeta
twice, using ∂ₐζ(s,a)=−sζ(s+1,a). They also send −1/(s−1) to −1/(s+1)
and preserve −1/(2a). Thus the elementary endpoint corrections are part of
the same differential hierarchy.

For j a positive integer, (12) gives the familiar rational samples
\[
 F_s(j^2)=j^{s-1}(\zeta(s)-H_j^{(s)})
                           -\frac1{s-1}+\frac1{2j}.
\]
These are rational linear forms in the indicated single ζ(s). However,
knowing the discrete values F₃(j²) is not the same as knowing the derivatives
in (8). Put \(Z_k=\zeta(k,a)\). Direct differentiation of (12) at s=3 gives
\[
 \mu_3((t+a^2)^{-2})=-Z_3+\frac{3a}2Z_4-\frac1{4a^3},
\]
\[
 \mu_3((t+a^2)^{-3})=-\frac9{8a}Z_4+\frac32Z_5
                                             -\frac3{16a^5}. \tag{13}
\]
Substitution in (9) cancels the Z₃ and Z₄ terms exactly, leaving
\(a^4Z_5-1/4-1/(2a)\). At a=j these Hurwitz values are ordinary zeta
values minus rational harmonic sums. The third-order pole is precisely where
the ζ5 term appears in this transfer. No linear-independence assertion about
ζ3, ζ4, or ζ5 is needed.

More generally let \(J_k(a)=\mu_3((t+a^2)^{-k})\). Differentiating k−1 times
shows that it is a rational combination of ζ(3,a),…,ζ(k+2,a) and rational
functions of a. Its highest Hurwitz-zeta coefficient is exactly
\[
 \boxed{[\zeta(k+2,a)]J_k(a)
                  =\frac{k(k+1)}{2^k}a^{3-k}\ne0.}      \tag{14}
\]
This denotes the coefficient in the displayed derivative expansion, not a
claim that these special functions or their integer values are independent.
The highest contribution comes from differentiating ζ(3,a) at every step;
all other product-rule terms have lower zeta index. Combining (11) and (14)
for k=2m+1 yields aˢ⁻¹ as the coefficient of ζ(s,a), exactly as in (12).

There is a precise **kernel** obstruction to transferring this through values
alone: a finite combination of polynomials and simple kernels at negative
squares cannot equal the kernel in (10), since it has a nonzero pole of order
s−2. Equivalently there is no such finite value-only identity that is valid
for every underlying measure. This does not rule out unexpected numerical
relations for the particular constants ζ3 and ζ5; that would require a wholly
different arithmetic argument. At noninteger rational a, (12) contains Hurwitz
values and need not reduce to the same single ordinary zeta value.

## 6. What transfers constructively, and what does not

The hierarchy gives exact algorithms to transport a rational test function,
its polynomial moments, and its coefficient normalizers to the base measure.
It also explains why the common analytic decay scale of the odd positive
measures can coexist with different arithmetic problems: the target index is
encoded in higher pole jets and their cancellation, not solely in the location
of the simple-pole samples.

The differential transform is not a pointwise positivity-preserving change
of variables. For instance
\[
 P_1(\Theta)(t-1)^2=6t^2-7t+5/3,
\]
which is −3/8 at t=7/12, although (t−1)²≥0 everywhere. Its μ₃ integral is
nevertheless μ₅((t−1)²)>0 by Theorem1. Thus a positive μₛ Gram matrix pulls
back to a globally positive differential bilinear form; it does not become
a μ₃ Gram matrix by simply applying Pₘ to each basis vector. The operator
is neither multiplicative nor a positivity-preserving kernel replacement.
For s=5 its Leibniz rule explicitly introduces p″q, p′q′, and pq″ terms.

Accordingly (3)–(11) provide a constructive coefficient-integrality bridge,
not an automatic transfer of the already optimized ζ3 integer approximants.
One still has to control the higher-pole moment values, the cancellation of
unwanted lower zeta components, and the actual normalization of the proposed
linear forms or determinants. The hierarchy identifies the work needed and
does not resolve it merely by differentiation.

## Reproduction

```sh
python codex/computations/check_odd_measure_differential_hierarchy.py
python codex/computations/check_odd_measure_differential_hierarchy.py --numeric
```

The checker verifies the general symbolic IBP operator, its Stieltjes adjoint,
the integer Newton decomposition and normalizer samples for s=3,…,23, the
pole expansion endpoints, the general Hermite raising identity, (13), the
ζ5 cancellation, and the positivity counterexample. Optional quadrature uses
the original positive densities to test four resolvent samples and three
rational-function transfers at40-digit precision. Finite normalizer samples
corroborate but do not replace the all-prime proof of sharpness above.

Primary background: the defining moment construction and Hermite integral
formula are recorded with source context in `odd-moment-generalization.md`
and `integer-depletion-independent-review.md`. This note derives the
hierarchy and its consequences directly; historical priority has not been
investigated.
