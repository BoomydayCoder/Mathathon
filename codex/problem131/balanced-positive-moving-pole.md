# An even balanced index gives a positive moving-pole remainder

Status: the all-index positive transport and associated arithmetic have passed
the [independent audit](balanced-moving-pole-independent-review.md). This is a different approximation from the
squared-polynomial many-pole trial. It does not prove irrationality or an
asymptotic primitive-error estimate.

Fix m>=1, take an even index n=2r>=2, and let p=p_(m,n) be the monic balanced
polynomial. Put N=2m+1 and D=mn+m+1. Define the integer polynomial

    S_(m,n)(t)=t product_(k=1)^(r-1)(t+k^2)^(2m+2)
                      (t+r^2)^(2m+1).                  (1)

For any real a>r, the exact positive identity is

    mu_(2m+1)[p(t)/(t+a^2)]
      = [n!/(D)_n] mu_(2m+1)[
           S_(m,n)(t)/product_(k=-r)^r(t+(a+k)^2)] >0.    (2)

## 1. Proof by moving the divided differences

Use the gamma weights w_s(x), the polynomial F and the operators Dcal,B_s
from `balanced-generalized-wilson-rodrigues.md`. Since n is even, mn is even,
so p(-x^2)=F(x)/(D)_n without an additional sign. Apply the contour adjoint n
times with the test function

    g(x)=h_j(x)/(a^2-x^2),
    h_j(x)=sin(pi x)^(2(m-j)), j=1,...,m.

The previously proved periodic-multiplier rule gives

    Dcal^n g(x)=h_j(x+n/2) Dcal^n[1/(a^2-x^2)]
              =h_j(x) Dcal^n[1/(a^2-x^2)].

The second equality is where the even index is essential: n/2=r is an
integer, and h_j has period one. The ordinary divided difference of the
resolvent z ->1/(a^2-z) is the reciprocal product of its pole distances.
Therefore

    Dcal^n[1/(a^2-x^2)]
       =n!/product_(k=-r)^r[a^2-(x+k)^2].                (3)

At x=iy, pairing the two linear factors and reversing k in one of them
converts this denominator exactly into

    product_(k=-r)^r[y^2+(a+k)^2].

The gamma recurrence for the integer shift n/2=r gives

    w_n(iy)/w_0(iy)
      =(iy)_r(-iy)_r[(1+iy)_r(1-iy)_r]^(2m+1)
      =S_(m,n)(y^2).                                    (4)

Combining (3)-(4) with the adjoint identity proves (2) first for each
hyperbolic component omega_j, and then for their positive linear combination
mu_(2m+1).

There is one new analytic point compared with polynomial orthogonality.
After s divided differences, the poles of the rational test factor have
real parts whose absolute values are at least a-s/2. At every one of the n
contour-deformation steps s<=n-1, this is strictly larger than 1/2 because
a>n/2. Hence the test functions introduce no poles into the deformation
strip. At infinity they have only polynomial growth, so the previously
established exponential tail bounds remain valid. Endpoint integrability
also holds because every a+k is strictly positive. The integrand on the
right of (2) is strictly positive for t>0, proving the strict inequality.
No claim at the boundary a=r is made here.

## 2. Explicit rational lower approximants to one odd zeta value

For an integer a>r, use the established resolvent formula

    mu_(2m+1)[1/(t+a^2)]
       =a^(2m)(zeta(2m+1)-H_a^(2m+1))-1/(2m)+1/(2a).

Polynomial division gives the rational quantity

    Q_(m,n,a)=H_a^(2m+1)+1/(2m a^(2m))-1/(2a^(2m+1))
       - mu_(2m+1)[(p(t)-p(-a^2))/(t+a^2)]/[a^(2m)p(-a^2)].            (5)

All moments in the final term are rational polynomial moments. The balanced
polynomial has even degree mn and positive roots, so p(-a^2)>0. Equation (2)
therefore gives the exact positive error

    0<zeta(2m+1)-Q_(m,n,a)
      = mu_(2m+1)[S_(m,n)(t)/product_(k=-r)^r(t+(a+k)^2)]
           /[a^(2m)(n!)^(2m) R_(m,n)(a)],               (6)

where R_(m,n)(a)=p(-a^2)/p(0)>0. The factorial cancellation in (6) uses the
proved constant term p(0)=(n!)^(2m+1)/(D)_n for even n.

The positive Andrews sum supplies the general bound

    0<zeta(2m+1)-Q_(m,n,a)
      <= mu_(2m+1)[S_(m,n)/product_(k=-r)^r(t+(a+k)^2)]
              /[a^(2m)((a+1)_n)^(2m)].                 (7)

This is an explicit positive error representation for a single odd zeta
value. It contains no unproved root or sign assumption, and no lower odd
zeta values occur in (5).

## 3. Scope of the arithmetic question

The integer-grid theorem controls the zeta coefficient R_(m,n)(a) after
multiplication by d_n^(2m-2). The associated rational numerator in (5) still
requires its own arithmetic analysis. The many-pole normalizer from
`balanced-andrews-grid-arithmetic.md` does not automatically apply: the
polynomial quotient in (5) has degree mn-1, and its Bernoulli moments must be
retained. The positive remainder alone does not show that a primitive integer
linear form is small.

For orientation only, exact rational computations with m=2 and a=n give
log10 of the primitive positive form B(zeta(5)-Q), where Q=A/B is reduced,
approximately

    n=4: 1.09844; n=6: 3.05766; n=8: 3.85043;
    n=10: 5.52512; n=12: 8.05054.

These finitely many positive values do not establish an all-index obstruction.
At n=2,a=4 the primitive value is below one, so even the finite behavior is
not a uniform failure statement. Numerical diagnostics here are not used to
prove (2), (5), (6), or any irrationality claim.

## 4. Integer-valued polynomials clear the associated numerator

The grid theorem supplies enough structure to prove an all-index associated
numerator bound as well. This arithmetic statement holds for every n>=1,
including odd n; the positive error statement above remains restricted to
even n.

Put d=mn, H=d_n^(2m-2), and

    f(x)=H p_(m,n)(-x^2)/p_(m,n)(0).

The proved grid theorem says f belongs to Int(Z), the ring of rational
polynomials taking integer values at every integer. It has degree 2d, is
even, and f(0)=H. For a positive integer a define

    P_a(x)=x^(2m)[f(x)-f(a)]/(x-a),
    L=2d+2m-1, J=max(a,L), M=(2m)! d_J^(2m+1).           (8)

Then there is an explicitly computable integer B such that

    M H mu_(2m+1)[(p(t)/p(0))/(t+a^2)]
       = M a^(2m) f(a) zeta(2m+1) - B.                 (9)

In particular both affine coefficients in (9) are integers. For even n and
a>n/2, the left side is strictly positive by (2) and p(0)>0.

### The single-quotient Bernoulli identity

Let BB(x^k)=B_k, with B_1=-1/2. The known rational moment formula gives, for
any polynomial g,

    mu_(2m+1)(g)=BB D^(2m-1)[x^(2m+1)g(-x^2)]/(2m)!.

The Bernoulli generating function also gives, for every polynomial Q,

    BB(Q(-x))=BB(Q(x))+Q'(0).                            (10)

If g(t)=[p(t)/p(0)-p(-a^2)/p(0)]/(t+a^2), then evenness of f implies

    H x^(2m+1)g(-x^2)=-[P_a(x)-P_a(-x)]/2.

Apply D^(2m-1), then (10). Because this derivative order is odd and
P_a has a factor x^(2m), the result is

    H mu_(2m+1)(g)
      =-BB P_a^(2m-1)/(2m)!-[f(a)-H]/(2a).               (11)

The final term follows from
P_a^(2m)(0)=(2m)![f(a)-H]/a. Combining (11) with the rational pole formula
cancels the f(a)/(2a) boundary term and proves the exact identity

    H mu_(2m+1)[(p(t)/p(0))/(t+a^2)]
       =a^(2m)f(a)(zeta(2m+1)-H_a^(2m+1))
         -f(a)/(2m)+H/(2a)-BB P_a^(2m-1)/(2m)!.         (12)

Thus (9) uses

    B=M {a^(2m)f(a)H_a^(2m+1)+f(a)/(2m)-H/(2a)
                       +BB P_a^(2m-1)/(2m)!}.           (13)

The reflection step is useful arithmetically: only one linear polynomial
quotient remains, rather than two successive divisions by x-a and x+a.

### Three elementary integer-valued-polynomial facts

First, if F belongs to Int(Z) and deg F<=h, then its Newton expansion about
an integer a has integer coefficients:

    F(a+x)=sum_(k=0)^h c_k binom(x,k), c_k in Z.

Since binom(x,k)/x=binom(x-1,k-1)/k for k>=1,

    d_h [F(x)-F(a)]/(x-a) belongs to Int(Z).              (14)

Multiplication by x^(2m) preserves integer values, so d_(2d)P_a belongs to
Int(Z) and has degree at most L.

Second, on polynomials of degree at most L the shift identity exp(D)=1+Delta
implies the finite nilpotent-operator identity

    D=log(1+Delta)=sum_(k=1)^L (-1)^(k+1)Delta^k/k.

Every Delta^k preserves Int(Z), so d_L D preserves Int(Z). Iteration shows
that d_L^r D^r preserves Int(Z) for each r>=0.

Third, the Bernoulli generating function gives

    BB[(1+z)^x]=log(1+z)/z,
    BB binom(x,k)=(-1)^k/(k+1).                          (15)

Hence d_(h+1) BB F is an integer whenever F belongs to Int(Z) and deg F<=h.

Apply these facts successively to P_a and r=2m-1. Its differentiated degree
is at most L-r=2d. Consequently

    d_(2d) d_L^(2m-1) d_(2d+1) BB P_a^(2m-1) is integral.              (16)

Each lcm index in (16) is at most J, and the total exponent is 2m+1.
Therefore d_J^(2m+1) clears that Bernoulli value. The same lcm power clears
the harmonic sum in (13) and the denominator a; the fixed factorial (2m)!
clears 2m and the displayed factorial. This proves B is integral and finishes
(9). No bound on individual Bernoulli-number denominators was needed.

For m=2, the explicit multiplier is 24 d_(max(a,4n+3))^5 and the initial grid
normalization is H=d_n^2. This is a proved exponential bound when a=O(n),
not a claim of arithmetic optimality. Its size currently overwhelms the
available error estimates. The positive integer relation (9) therefore does
not supply an irrationality contradiction.
