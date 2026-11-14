# Independent review of the common fixed-odd analytic budget

Status: reviewed the theorem in `common-odd-zeta-analytic-budget.md` and strengthened its proof with explicit uniform estimates. No error was found in the fixed-s analytic conclusion. The arithmetic gate remains unproved for general odd s. This review does not convert a family of small rational polynomials into integer polynomials without a normalization theorem.

## 1. Scope and prerequisites

The theorem fixes one odd integer s>=3 before K grows through positive multiples of40. It uses the proved positive-moment functional from `odd-moment-generalization.md`, the independently certified equilibrium comparison with α=3/40 andλ=37/40, Andréief's identity, and elementary factorial asymptotics.

The claimed conclusion is

    limsup K^-2 log[S_K Δ_(K,s)(ζ(s))] <= -69/50.

All constants in lower-order estimates may depend on the fixed s. There is no claim uniform when s grows with K.

## 2. Positivity and exact polynomial degree

At ζ(s), the matrix is the Gram matrix of1,t,...,t^(h-1) against the positive measure obtained by multiplying w_s(y)dy by D_N(y²)^6/D_K(y²). Its density is positive for every y>0. Every polynomial moment exists because w_s has exponential decay. Therefore the matrix is positive definite and its determinant is strictly positive.

Each entry is affine in X. After cancellation of D_N against D_K, the surviving poles are exactly -j² for N<j<=K. Their number K-N equals h. The X-coefficient matrix is

    sum_(j=N+1)^K w_j v_j v_j^T,
    v_j=(1,-j²,...,(-j²)^(h-1))^T,
    w_j=j^(s-1)D_N(-j²)^6/D_K'(-j²).

Every w_j is nonzero, and the h column vectors form a Vandermonde matrix on distinct nodes. Its determinant is nonzero, so the determinant polynomial has degree exactly h. Positivity of the individual w_j is unnecessary for this algebraic leading-coefficient argument; positivity of the actual Gram measure proves the real determinant's sign separately.

## 3. A global weight bound with an explicit constant

The geometric expansion gives

    w_s(y)=2(2π)^(s-1)y^s/(s-1)!
                  sum_(l>=1) l^(s-1)e^(-2πly).

For l>=1,

    l^(s-1) <=l(l+1)...(l+s-2)
             =(s-1)!binom(l+s-2,s-1).

Summing the binomial series and using1-e^-x>=x/(1+x) yields

    0<w_s(y)<=(1/π)(1+2πy)^s e^(-2πy)
                <=2(2π)^(s-1)(1+y)^s e^(-2πy).

This controls all y>0, including the origin, without relying only on a large-y asymptotic. The extension to y=0 follows from w_s(0)=1/π. This explicit proof was added to the main note.

The extra polynomial factor K^s after y=Ksqrt(t) occurs once per particle, so its logarithmic contribution is sh logK=O_s(KlogK), not K² times a new field.

## 4. Uniform Riemann-sum error, including t=0

For1<=m<=K set

    e_m(t)=sum_(j=1)^m log(t+(j/K)²)
            -K integral_0^(m/K)log(t+u²)du.

The function of u is increasing, so e_m(t)>=0. For t>0, differentiation gives a right-Riemann-sum error for the decreasing function1/(t+u²), so e_m'(t)<=0. The finite t=0 limit is exactly

    e_m(0)=2log(m!)-2m logm+2m.

Elementary factorial bounds imply

    0<=e_m(t)<=2logK+2  for every t>=0.

The singular value log0 occurs only at the integration endpoint and is integrable; no grid evaluation uses it. This closes the potential uniformity gap in a proof that merely appeals to smooth Riemann sums away from0.

For the rational products, the correction is6e_N-e_K<=12logK+12. This contributes only O(KlogK) after summing h particles.

## 5. Explicit K-power accounting

After y_i=Ksqrt(t_i):

| factor | power of K |
|---|---|
| squared Vandermonde | 2h(h-1) |
| D_N^6/D_K over h particles | 2h(6N-K) |
| Jacobians | h |
| polynomial weight upper bound | sh |
| numerator Riemann-sum bound | 12h |

The total is

    2h(h+6N-K)+(s+11)h.

All factors such as C_s^h,2^-h,e^(12h), and1/h! have logarithms O_s(KlogK) or smaller. Thus the coefficient of K²logK is

    2λ(λ+6α-1),

independent of s.

The source normalizer satisfies

    log S_K=(2λ-12αλ-2λ²)K²logK+C*K²+O(KlogK),

where

    C*=-2λ+12αλ(1-logα)+3λ²-2λ²log(2λ).

For the last term of S_K, the estimate

    sum_(i=1)^(h-1)log((2i)!)
       =h²log(2h)-(3/2)h²+O(hlogh)

follows from Stirling's formula and integral comparison. It gives the displayed3λ² and -2λ²log(2λ) contributions. The K²logK coefficients cancel exactly.

## 6. The infinite domain can be handled without discarding the polynomial

There is no need to absorb the entire polynomial weight into an unspecified exponential constant. Retain it and apply the equilibrium bound to V(t)-sqrt(t)/K. The fixed support and the reviewed tail estimate give

    sup_(t>=0)[2Uρ(t)-V(t)+sqrt(t)/K] <=M+sqrt2/K

for K>=2. Bounded-density circle regularization produces only O(KlogK) error, as detailed in `equilibrium-energy-review.md`.

The remaining one-variable factor has exact finite integral

    Z_s=integral_0^infinity t^-1/2(1+sqrt(t))^s e^-sqrt(t)dt
       =2sum_(j=0)^s binom(s,j)j!.

It is integrable at0 because of t^-1/2, and at infinity because of the exponential. The h-fold integral therefore costs only Z_s^h, whose logarithm is O_s(K). For s=5, Z_5=652 exactly, matching the source. This more explicit argument was inserted into the main note.

Combining the configuration estimate, these factors, and S_K leaves

    log[S_KΔ_(K,s)(ζ(s))]
       <=[λM-I(ρ)+C*]K²+O_s(KlogK).

The certified comparison coefficient is less than -69/50. This establishes the stated limsup bound for every fixed odd s.

## 7. The exact arithmetic gate

The analytic result alone gives positive values of a rational polynomial. To deduce irrationality one still needs positive rational normalizers m_(K,s) such that

    m_(K,s)S_KΔ_(K,s)(X) ∈ Z[X],
    limsup K^-2 log m_(K,s)<69/50.

If these are proved, the degree h=37K/40 and positivity yield the familiar integer contradiction under ζ(s)=a/b: multiply the polynomial value by b^h. Its logarithmic denominator cost is O(K), whereas the net decay is a negative constant times K².

The gate is a sufficient condition, not a theorem asserting that suitable normalizers exist for all odd s. The s=5 arithmetic claim belongs to the source's separate valuation/integrality argument and is not re-proved by this analytic theorem. For other s it requires new work. In particular Bernoulli denominators and local pole costs cannot be dropped merely because the leading real field is unchanged.

## Review outcome

Accepted the fixed-s analytic theorem with its stated prerequisites. Expanded the proof rather than changing its conclusion: exact degree, uniform t=0 handling, K powers, explicit global weight constant, and a direct unbounded-domain integral are now recorded. No new irrationality theorem follows until the arithmetic gate is met.
