# Consecutive-pole transport and a second pure odd-zeta approximation

Status: accepted in [the independent review](balanced-consecutive-pole-independent-review.md). This extends the accepted
balanced Rodrigues and positive moving-pole identities. It supplies explicit
positive rational approximations and an exponential-size integer normalizer;
no irrationality inequality or historical novelty is claimed.

## 1. A closed family of kernels under the quadratic difference operator

Write

    Dcal g(x)=[g(x+1/2)-g(x-1/2)]/(2x),
    g_(a,K)(x)=1/product_(j=0)^(K-1)[(a+j)^2-x^2], K>=1.

Then, as an identity of rational functions,

    Dcal g_(a,K)=K g_(a-1/2,K+1),
    Dcal^n g_(a,K)=(K)_n g_(a-n/2,K+n).                 (1)

Indeed, after putting the first difference over the denominator of
h=g_(a-1/2,K+1), the two numerators are respectively

    (a+K-1/2-x)(a-1/2+x),
    (a-1/2-x)(a+K-1/2+x).

Their difference is 2Kx. Iteration proves the second identity. Equivalently,
the family consists of gamma quotients Gamma(a+x)Gamma(a-x) divided by
Gamma(a+K+x)Gamma(a+K-x), and Dcal shifts its two endpoints outward.

Let p=p_(m,n) be the monic balanced polynomial, let n=2r be even, and put

    D=mn+m+1,
    S_(m,n)(t)=t product_(j=1)^(r-1)(t+j^2)^(2m+2)
                        *(t+r^2)^(2m+1).

The same contour adjoint used in `balanced-positive-moving-pole.md` now gives,
for every real a>r and positive integer K,

    mu_(2m+1)[p(t)/product_(j=0)^(K-1)(t+(a+j)^2)]
      = (K)_n/(D)_n * mu_(2m+1)[S_(m,n)(t)/
                      product_(j=0)^(K+n-1)(t+(a-r+j)^2)] >0.        (2)

There are no additional contour issues: after s half-shifts the nearest
kernel pole is at a-s/2. For each required strip shift s<=n-1, the assumption
a>n/2 keeps it strictly outside the strip. The gamma-weight strip and tail
estimates are exactly those already proved for the single-pole kernel.
Since n is even, the periodic multipliers return to their original values,
and the quotient of the shifted gamma weights is S_(m,n). Formula (1) then
replaces the factor n! from the single-pole case by (K)_n.

## 2. Matching the block length removes the Rodrigues prefactor

Choose the integer block length K=D=mn+m+1 and an integer a>r, and set
b=a+K-1. For a<=j<=b define the rational partial-fraction coefficient

    c_j=p(-j^2)/product_(ell=a, ell!=j)^b(ell^2-j^2).

Since deg p=mn=K-m-1, expansion at infinity gives

    sum_(j=a)^b c_j j^(2m)=(-1)^m.                      (3)

In detail, p(t)/product_(j=a)^b(t+j^2)=t^(-m-1)+O(t^(-m-2)); comparison
with sum c_j/(t+j^2) proves (3), while all earlier coefficients vanish.
Using the exact resolvent evaluation

    mu_(2m+1)[1/(t+j^2)]
       =j^(2m)(zeta(2m+1)-H_j^(2m+1))-1/(2m)+1/(2j),

put

    q_(m,n,a)=(-1)^m sum_(j=a)^b c_j
                  [j^(2m)H_j^(2m+1)+1/(2m)-1/(2j)].

Then q is rational and (2) becomes

    (-1)^m[zeta(2m+1)-q_(m,n,a)]
      =mu_(2m+1)[S_(m,n)(t)/product_(j=a-r)^(b+r)(t+j^2)]>0.          (4)

Thus the approximations are below zeta for even m and above it for odd m.
The prefactor has disappeared because (K)_n=(D)_n. At the smallest permitted
integer pole a=r+1, the positive remainder simplifies further to

    (-1)^m[zeta(2m+1)-q_(m,n,r+1)]
      =mu_(2m+1)[S_(m,n)(t)/product_(j=1)^((m+1)n+m+1)(t+j^2)].       (5)

This differs from the previous many-pole construction using p^2: the numerator
before transport is p itself, and positivity comes from its Rodrigues
adjoint. The special block length isolates one odd zeta value by degree alone.

## 3. An explicit integer normalizer

Let

    C=binom((m+1)n+m,n), d_h=lcm(1,...,h), H=d_n^(2m-2).

The proved grid theorem says

    H R_(m,n)(j) in Z,
    p(-j^2)=(-1)^(mn)(n!)^(2m) R_(m,n)(j)/C.

The following quantity is a positive integer and clears q in (4):

    M=m C * (2b)!/[(n!)^(2m)((a-1)!)^2] * H * d_b^(2m+1),
    M q_(m,n,a) in Z.                                  (6)

To see first that M is integral, the factorial arguments in its denominator
sum to 2mn+2a-2, whereas 2b=2mn+2a+2m. Thus the factorial ratio in (6) is
an integer multinomial coefficient times (2m+2)!.

For the coefficient clearing, extend the partial fraction to the full grid
1,...,b by inserting the numerator

    p(t) product_(ell=1)^(a-1)(t+ell^2).

The reciprocal full-grid derivative equals

    1/product_(ell=1,ell!=j)^b(ell^2-j^2)
       =(-1)^(j-1) 2j^2 binom(2b,b-j)/(2b)!.

Also, at integer j>=a,

    product_(ell=1)^(a-1)(j^2-ell^2)
       =((a-1)!)^2 binom(j-1,a-1)binom(j+a-1,a-1).

Consequently M c_j is, up to sign, 2m j^2 d_b^(2m+1) times an integer:
all remaining factors are binomial coefficients and H R_(m,n)(j). This clears
both constants 1/(2m),1/(2j) and the harmonic denominator in the formula
for q, proving (6).

For a=O(n), the factorial ratio is exponential in n up to a polynomial
factor, because the factorial arguments balance at order n. The binomial
factor C and both lcms also have at most exponential growth. This is only
an exponential normalizer theorem: neither its size nor (4) implies that the
normalized positive remainder tends to zero. No primitive-denominator bound
or successful irrationality estimate is claimed here.

The derivations of (3) and (6) are algebraic and remain valid for odd n as
well. Positivity and the explicit S-transport in (2), (4), and (5) are asserted
only for even n.

The exact minimal-pole decay rate and the scope of the fifth-zeta normalizer
limitation are proved separately in `balanced-consecutive-pole-rate.md`.
