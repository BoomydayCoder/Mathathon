# Integral polynomial moments do not rescue the depleted determinant

Status: independently reviewed twice; see coprime-primitive-obstruction-independent-review.md
and coprime-depletion-primitive-independent-review.md. Exact scalar certificate saved.
The local prime theorem is proved separately in
[coprime-depletion-pole-arithmetic.md](coprime-depletion-pole-arithmetic.md).
This is an actual obstruction for scalar normalization of one specified family,
not merely failure of a sufficient bound. It uses the prime number theorem;
that input is unrelated to the PNT-free positive result in the main zeta5 guide.

## 1. The tempting construction and the precise negative result

Fix an odd integer s>=3 and let

    sigma_s(y)=w_s(y)-w_s(2y)/2-w_s(3y)/3+w_s(6y)/6 > 0.

This is the normalized coprime frequency depletion. Its exponential tail has
the same rate2*pi as w_s. Its polynomial moments have denominators supported
only at2 and3. Its integer-square pole moments are rational affine functions
of zeta(s), so it still defines a single-zeta rational functional sigma_X.

Set K=40n, N=3n, h=K-N, D_j(t)=product_(r=1)^j(t+r^2), and

    Delta_K(X)=det[ sigma_X(D_N(t)^6 t^(i+j)/D_K(t)) ]_(0<=i,j<h).

The determinant is positive at X=zeta(s). Let P_K be its primitive integer
part, with the sign chosen to be positive there. Thus Delta_K=g_K P_K for
a positive rational g_K and gcd of the coefficients of P_K is one.

**Theorem. For each fixed odd s>=3, eventually**

    P_K(zeta(s)) > exp(K^2/3).                          (1)

In particular, no nonzero rational scalar clearing the coefficients of this specific
Delta_K can produce values tending to zero. Every such integral scalar
multiple is a nonzero integer multiple of P_K. This theorem does not exclude
other pole grids, numerator powers, aspect ratios, or additional polynomial
combinations that alter the construction.

## 2. Monic normalization prevents an arithmetic bookkeeping mistake

Put C_s=(2^(s-1)-1)(3^(s-1)-1). The exact leading coefficient is

    ell_K=(-1)^(h(h-1)/2) C_s^h
           * product_(r=N+1)^K r^(s-1) D_N(-r^2)^5.    (2)

It is a nonzero integer. To see this, polynomial division leaves only poles
r=N+1,...,K with residues D_N(-r^2)^5 / E'(-r^2), where E=D_K/D_N.
The X coefficient at each pole is C_s r^(s-1). Squared Vandermonde factors
cancel product E'(-r^2), with the sign in (2). This also proves degree h.

Let F_K=Delta_K/ell_K be monic. Let L_K be the least positive integer clearing
all coefficients of F_K. Then L_K F_K is primitive integral: a common prime
dividing its coefficients would divide its leading coefficient L_K and allow
that prime to be removed from the clearing integer. Hence

    P_K(zeta(s)) = L_K |F_K(zeta(s))|.                 (3)

This distinction matters. Individual bad-prime valuations of an unnormalized
rational polynomial do not by themselves give a lower bound on its primitive
value: common numerator factors elsewhere may offset them. A monic rational
polynomial has leading coefficient one, so its clearing factor is an integer.
No negative prime exponents can offset its denominator cost.

## 3. Actual denominators in a new prime band

The separate exact local theorem gives, for all sufficiently large primes
2K<p<=6K outside the fixed finite divisors of C_s(2^s-1),

    v_p^G(Delta_K)=-s [K-max(N,floor(p/6))].           (4)

Every factor of ell_K apart from the fixed C_s has prime divisors at most
K+N<2K; hence ell_K is a unit at these primes. Since N/K=3/40<1/3,

    v_p(L_K)=s[K-floor(p/6)] .                        (5)

The prime number theorem and partial summation therefore give

    liminf_(K->infinity) log L_K/K^2
      >= s integral_2^6 (1-u/6) du = 4s/3.           (6)

The O(1) floor error contributes O(K), since sum_(p<=6K)log p=O(K).
Fixed exceptional primes eventually lie below2K. Other primes can only add
to log L_K, because L_K is a positive integer.

## 4. An elementary lower bound for the real determinant

A lower bound is needed here: the familiar upper smallness estimates cannot
prove a construction fails. We use only a uniform trial distribution and
Jensen's inequality, without an equilibrium measure.

Write alpha=3/40, lambda=37/40. There exist constants0<c_s<C_s' such that

    c_s w_s(y) <= sigma_s(y) <= C_s' w_s(y), y>0.

Indeed the ratio is positive and continuous, tends to1/3 at zero and to1 at
infinity. In particular sigma_s(y)>=c_s c_s^* y^s exp(-2*pi*y), since the
first exponential term lower-bounds w_s. Constants depend on fixed s only.

Andreief's Gram determinant identity gives

    Delta_K(zeta(s)) = 1/h! integral_(0,infinity)^h
        product_i [D_N(y_i^2)^6/D_K(y_i^2) sigma_s(y_i)]
        product_(i<j)(y_i^2-y_j^2)^2 dy_1...dy_h.      (7)

Restrict the integral to0<y_i<Kb, with b=1/2, and write y_i=K x_i.
Take x_i independent uniform on(0,b). The logarithm of every factor is
integrable; coincident x_i form a null set, and log-distance singularities
are integrable. Jensen gives log E(exp Z)>=E Z for the logarithm Z of this
positive integrand. The factorial, interval-length, fixed-density and y^s
terms contribute at most O_s(K log K) after separating the powers of K.

Define

    I_t(x)=integral_0^t log(x^2+u^2)du
          =t log(x^2+t^2)-2t+2x atan(t/x),
    J_t=1/b integral_0^b I_t(x)dx
       =t log(b^2+t^2)-3t+(t^2/b)atan(b/t)+b atan(t/b).

The last formula follows by elementary integration; its limiting value is
used at t=0. Uniform integration in x makes the logarithmic Riemann sums
continuous even at u=0. Thus their limits are J_alpha and J_1.
The pair interaction is also explicit:

    E log|x-y|=log b-3/2,
    E log(x+y)=log b+2log2-3/2.

Consequently the scale-free determinant constant is bounded below by

    V=lambda*(6 J_alpha-J_1-pi*b)
         +lambda^2*(2log b+2log2-3).                  (8)

At b=1/2 the second bracket is simply-3. More explicitly the powers of K
from the rational weight, Vandermonde, dy and y^s sum to
10Nh+(s-1)h; this is exactly the power of K in the leading coefficient (2).
Hence no K^2 log K term remains in the monic ratio.

The constant from the leading coefficient is

    L=5 integral_alpha^1 integral_0^alpha log(x^2-u^2)du dx
      =5[H(1+alpha)-H(2alpha)-H(1-alpha)
                                      -2alpha(1-alpha)],
    H(v)=v^2 log(v)/2-v^2/4, H(0)=0.                 (9)

This limit can also be obtained directly by Stirling from
abs(D_N(-r^2))=(r+N)!/[r*(r-N-1)!]. The isolated r^(s-1) and C_s factors
contribute only O_s(K log K), with their logarithmic scale already separated.
Combining (7)--(9) proves

    liminf log|F_K(zeta(s))|/K^2 >= V-L.             (10)

## 5. The constants leave a comfortable strict gap

The standard-library rational certificate bounds every logarithm by its
atanh series and every arctangent by its alternating series, using Machin's
formula for pi. It gives approximately

    V=-4.135657073562474...,
    L=-0.558208479705225...,
    V-L=-3.577448593857249... > -18/5.                (11)

Only the final strict rational inequality is used. There is no floating-point
input in its acceptance checks. Therefore (3),(6),(10) imply

    liminf log P_K(zeta(s))/K^2 >= 4s/3+V-L
                                    >4s/3-18/5 >=2/5.

This proves the eventual bound (1), with1/3<2/5.

The conclusion is about the actual primitive integer polynomial, not just an
expensive sufficient normalizer. Removing polynomial Bernoulli denominators
has moved an unavoidable cost into primes between2K and6K. For this fixed
positive determinant, that cost already exceeds an elementary lower estimate
for the real decay, before accounting for primes2 and3.
