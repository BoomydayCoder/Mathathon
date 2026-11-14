# A positive moment functional for every odd zeta value

Status: independent extension of the construction read in Fauzan (2026), Section
2.4. Exact analytic theorem below, conditional only on the classical Hermite
integral formula stated explicitly. It does NOT prove irrationality. Historical
novelty not assessed. Arithmetic normalization and decay are unresolved.

Let s>=3 be odd, f(y)=1/(exp(2*pi*y)-1), and

    w_s(y)=2*y^s*f^(s-1)(y)/(s-1)!.

**Theorem.** w_s is positive on y>0, tends to 1/pi at zero, and decays
exponentially times a polynomial at infinity. Define a rational functional on
polynomials and simple poles at negative integer squares by

    mu_X(t^e)=(-1)^e B_(2e+2)*(2e+s)!/((s-1)!*(2e+2)!), e>=0,
    mu_X(1/(t+j^2))=j^(s-1)*(X-H_j^(s))-1/(s-1)+1/(2j), j>=1.

Then for every rational function R in this domain,

    mu_(zeta(s))(R)=integral_0^infty R(y^2)*w_s(y) dy.

Thus all these odd zeta values have a positive single-variable moment model;
their isolation in such a functional is not by itself the main obstruction.

**Proof.** The geometric series for f, differentiated s-1 times on compact
subsets of (0,infinity), gives positive summands because s-1 is even:

    w_s(y)=2*(2*pi)^(s-1)/(s-1)! * y^s * sum_(l>=1) l^(s-1) exp(-2*pi*l*y).

At zero use f(y)=1/(2*pi*y)-1/2+O(y); its differentiated Laurent expansion
gives w_s(y)->1/pi. At infinity the l=1 term bounds the tail up to a constant.
Integrating monomials termwise by Tonelli gives

    integral y^(2e)*w_s(y)dy
      =2*(2e+s)!/((s-1)!*(2*pi)^(2e+2))*zeta(2e+2).

Euler's even-zeta formula gives exactly the displayed Bernoulli moment.

For a>0, polynomial division of y^s/(y^2+a^2), followed by s-1 derivatives,
kills the polynomial part and gives

    D^(s-1)[y^s/(y^2+a^2)]
        =(s-1)!*a^(s-1)*Im((a-i*y)^(-s)).

To verify the sign, put s=2m+1: the remainder is (-a^2)^m*y/(y^2+a^2),
and D^(2m) Im((a-i*y)^(-1))=(-1)^m*(2m)!*Im((a-i*y)^(-s)).
The two signs cancel. Integrate by parts s-1 times. Each boundary term at
zero is O(y): in the term g^(k) f^(s-2-k), g=y^s/(y^2+a^2), the powers
are y^(s-k) and y^(-(s-1-k)). At infinity exponential decay dominates.

Hermite's formula, valid here for real a>0 and integer s>1, is

    zeta(s,a)=1/(2*a^s)+a^(1-s)/(s-1)
              +2*integral_0^infty f(y)*Im((a-i*y)^(-s))dy.

Consequently the pole integral is

    a^(s-1)*zeta(s,a)-1/(2*a)-1/(s-1).

For integer a=j, use zeta(s,j)=zeta(s)-H_j^(s)+j^(-s).
This yields the stated functional, and linearity proves the result. QED.

**Corollary (positive Hankel candidates).** For integers K>N>=0 and even r>=0,
put D_m(t)=product_(j=1)^m(t+j^2). For any h>=1 define

    G(X)_(i,j)=mu_X(D_N(t)^r*t^(i+j)/D_K(t)), 0<=i,j<h.

All entries are affine in X. The rational functions have only simple poles after
cancellation (including r=0). At X=zeta(s), G is a positive definite Gram matrix:
its quadratic form at a nonzero polynomial q is the integral of
D_N(y^2)^r*q(y^2)^2*w_s(y)/D_K(y^2), strictly positive. Thus det G(zeta(s))>0.
The determinant is a rational polynomial in X of degree at most h.
No claim of degree exactly h is needed for this corollary.

**The missing inequality.** Clearing all coefficients produces some integer
polynomial P. Irrationality follows only if one constructs a sequence with
-log|P_n(zeta(s))|/max(1,deg P_n)->infinity. Positivity alone and a small
unnormalized determinant do not imply this. Prime-by-prime normalization costs
can overwhelm real decay. Fauzan's s=5 analysis addresses precisely that issue;
this generalization does not transfer those estimates automatically.

**Why this connection matters.** The same positive-moment idea links Beukers
integrals, Gram matrices, orthogonal polynomials, and determinant smallness.
Higher degree polynomials can replace linear rational approximations: under a
rationality assumption a/b, b^degree*P(a/b) is an integer. If degree is O(n)
and decay is exp(-c*n^2), every fixed denominator is eventually defeated.

Primary inputs: A. Fauzan, 'zeta(5) is irrational', v1, 17 September 2026,
DOI https://doi.org/10.5281/zenodo.22826419 (Section 2.4 read locally);
Hermite's Hurwitz-zeta integral, DLMF 25.11.29;
Euler even-zeta formula, DLMF 25.6.2. Check source-status audit for what was
retrieved and what formal claims have not been independently compiled.
