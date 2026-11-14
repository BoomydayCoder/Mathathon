# Frequency depletion, discrete-uniform cumulants, and exact denominator support

Status: independently reviewed; see integer-depletion-independent-review.md.
This generalizes the separately reviewed ell=2 Euler model. No new irrationality theorem or
historical novelty is claimed.

## 1. A positive construction for every integer ell>=2

For odd s>=3 use the positive density

    w_s(y)=c_s y^s sum_(n>=1) n^(s-1) exp(-2*pi*n*y),
    c_s=2(2*pi)^(s-1)/(s-1)!.

Define

    w_(s,ell)(y)=ell*w_s(y/ell)-w_s(y)
       =c_s ell^(1-s) y^s sum_(n>=1, ell does not divide n)
                    n^(s-1) exp(-2*pi*n*y/ell).

Cancellation of the terms with n divisible by ell proves the equality and
strict positivity. Its limit at zero is (ell-1)/pi, and its exponential tail
has rate2*pi/ell. For nu_(s,ell)(R)=integral R(y^2)w_(s,ell)(y)dy,

    nu_(s,ell)(R)=ell^2 mu_s(R(ell^2 t))-mu_s(R(t)).       (1)

For e>=0 put k=2e+2. The Bernoulli moment formula gives

    nu_(s,ell)(t^e)=(-1)^e (ell^k-1)B_k/(s-1)
                             * binom(k+s-2,s-2).        (2)

So a positive analytic operation multiplies Bernoulli numbers by exactly the
factor that cancels their denominators away from the primes dividing ell.

## 2. The unexpected elementary object is a uniform distribution

Let U be uniformly distributed on {0,1,...,ell-1}. Its exponential generating
function is

    M(z)=E(exp(zU))=(exp(ell*z)-1)/(ell*(exp(z)-1)).

Define cumulants simply as derivatives at zero of log M(z). No probability
theorem is needed. The Bernoulli generating series implies, for even k>=2,

    kappa_k=(ell^k-1)B_k/k.                            (3)

Since all derivatives of M at zero lie in Z[1/ell], M(0)=1, and differentiation
of log M uses only integer sums and products of these derivatives, every
kappa_k belongs to Z[1/ell]. One can see this without formal logarithm
denominators by recursively differentiating M'=M*(log M)'. Thus (2) becomes

    nu_(s,ell)(t^e)=(-1)^e kappa_k
                              * binom(k+s-2,s-1),      (4)

and belongs to Z[1/ell]. This removes even primes dividing s-1 but not ell;
that cancellation is less transparent in (2).

Equations (2) and (4) serve different purposes. Formula (4) fixes the prime
support. Formula (2), combined with von Staudt--Clausen, bounds the exponents
uniformly in e.

## 3. Exact optimal uniform denominator

**Theorem.** The least positive integer C such that

    C nu_(s,ell)(t^e) belongs to Z for every e>=0

is

    C_(s,ell)=product_(p prime, p divides ell) p^(1+v_p(s-1)).       (5)

In particular, it depends on the prime divisors of ell, not their multiplicity.
For ell=2 this recovers2^(1+v_2(s-1)); for ell=3 the denominator is
3^(1+v_3(s-1)); for ell=6 it is their product.

**Upper bound.** By (4), no prime outside ell occurs. For p dividing ell,
ell^k-1 is a p-adic unit. Von Staudt--Clausen says v_p(B_k)>=-1 for every
even k>=2. The binomial coefficient in (2) is integral. Thus its p-adic
valuation is at least -1-v_p(s-1), proving the bound.

**Sharpness at each prime.** Put a=v_p(s-1). We need an even k>=2 with
p-1 dividing k and binom(k+s-2,s-2) a p-adic unit. Lucas's theorem says
this last condition is equivalent to addition of k and s-2 in base p having
no carries.

For p=2 take k=2^a. Because s is odd, a>=1. The a-th digit of s-2 is zero
and k has no other nonzero digit, so there is no carry.

For odd p, write s-1=p^a u with p not dividing u. The a-th digit of s-2
is u-1 modulo p, hence lies between0 and p-2. Choose b>a so large that
p^b>s-2 and set

    k=p^a+(p-2)p^b.

Its digit at a is1, and its only other nonzero digit is p-2 at b. Therefore
addition to s-2 has no carry. Its digit sum is p-1, so k is divisible by p-1;
as p is odd this also ensures k is even. Von Staudt--Clausen now gives
v_p(B_k)=-1, while all other factors in (2) except s-1 are p-adic units.
Hence v_p(nu(t^e))=-1-a for e=(k-2)/2.

Each prime's maximal exponent is attained, possibly at a different e. Any
uniform denominator must contain them all, proving minimality of (5). QED.

## 4. Why this is not a free improvement in a zeta proof

The positive density now decays at rate2*pi/ell rather than2*pi. The normalized rescaled density

    sigma_(s,ell)(x)=w_s(x)-ell^(-1)w_s(ell*x)>0

has the original exponential rate2*pi and its e-th moment equals
ell^(-2e-2)nu_(s,ell)(t^e). Thus this normalization reintroduces denominators
growing with e.

There is also a precise pole cost. For j>0,

    nu_(s,ell)(1/(t+j^2))=mu_s(1/(t+(j/ell)^2))
                                      -mu_s(1/(t+j^2)).

Using Hermite's formula for real a>0 gives

    (j/ell)^(s-1) zeta(s,j/ell)-j^(s-1) zeta(s,j)
                                           -(ell-1)/(2j).       (6)

For ell=2 and positive integer j, the integer or half-integer Hurwitz value reduces to a rational multiple of
zeta(s), explaining the especially simple alternating-zeta model. For a
general ell, Hurwitz values at j/ell need not reduce to the same single zeta
value. Restricting j to multiples of ell restores integer Hurwitz arguments
and a single zeta value, but changes the pole spacing and its arithmetic.
We do not assert a general independence theorem for those Hurwitz values.

Even at ell=2 the previously established divided-difference examples show
that odd-prime denominators can return through the pole terms. Formula (5)
concerns all polynomial moments exactly; it is not a denominator theorem for
the whole rational-function space used by the determinant.

### The rescaling cost is attained, not just an upper bound

For each p dividing ell, there are arbitrarily large even k for which the
valuation in the sharpness proof is exactly -1-v_p(s-1). For odd p, simply let
b tend to infinity in the displayed construction. For p=2, use
k=2^a+2^b with b>a and2^b>s-2. The same no-carry argument works. Therefore
the normalized rescaled moment has exact valuation

    -k*v_p(ell)-1-v_p(s-1)

on an infinite subsequence. In particular no fixed denominator clears all its
moments. This is a statement about this rescaling, not every possible basis.

## 5. Two coprime depletions give positive INTEGER moments

Here is a stronger structural consequence. Let a,b>=2 be coprime integers.
Compose the two dilation-difference operations and set

    W_(s;a,b)(y)=ab*w_s(y/(ab))-a*w_s(y/a)-b*w_s(y/b)+w_s(y).

Because a and b are coprime, inclusion-exclusion in the exponential expansion
gives exactly

    W_(s;a,b)(y)=c_s(ab)^(1-s)y^s
      * sum_(n>=1, a does not divide n, b does not divide n)
                  n^(s-1) exp(-2*pi*n*y/(ab)) > 0.       (7)

Its limit at zero is(a-1)(b-1)/pi. Its e-th moment is

    M_e=(a^k-1)(b^k-1) mu_s(t^e),  k=2e+2.             (8)

**Corollary. Every M_e is an integer, for every odd s>=3.**
Indeed (a^k-1)mu_s(t^e) belongs to Z[1/a] by the theorem already proved;
multiplying by the integer b^k-1 preserves that ring. Reversing a,b puts
M_e also in Z[1/b]. Coprimality gives Z[1/a] intersect Z[1/b]=Z. QED.

Thus a strictly positive continuous measure on the positive half-line can
simultaneously encode higher-zeta pole values and have every polynomial moment
integral. Its ordinary polynomial Hankel determinants are positive integers:
Gram positivity supplies strict positivity and (8) supplies integrality.
Consequently those particular determinants cannot tend to zero. Any useful
small-determinant construction still needs rational pole factors or a rescaling,
whose arithmetic must be analyzed separately.

For a=2,b=3, formula (8) gives
M_0=(4-1)(9-1) B_2 binom(s,s-2)/(s-1)=2s.
The exact value2s checks the normalization directly. The first few higher
moments are independently computable integers; no asymptotic inference is
made from finite data.

For positive integer r, poles at j=6r involve only integer Hurwitz arguments.
Their zeta(s) coefficient is

    r^(s-1)(2^(s-1)-1)(3^(s-1)-1),

an integer. This follows by applying the original pole formula at r,2r,3r,6r.
The rational harmonic part and collisions between the pole locations still
require new denominator estimates; polynomial integrality does not supply them.
For other pole locations, additional rational-argument Hurwitz values can occur.

## 6. What this contributes to proof design

The same construction has three exact descriptions: frequency removal in a
positive integral, multiplication of Bernoulli moments by ell^k-1, and
cumulants of a finite uniform distribution. The first guarantees positivity;
the second supplies a uniform denominator exponent; the third proves the
prime support without complicated factorial cancellation. This separates
three jobs that can look accidental in one formula.

The useful research question is whether a pole basis can exploit this sharp
polynomial arithmetic without paying more through its collisions and scale.
The Euler case has explicit unfavorable examples, not a general impossibility
theorem. Larger ell introduce additional Hurwitz values, another constraint
that must be handled rather than omitted.
