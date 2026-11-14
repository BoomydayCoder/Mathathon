# A binomial coefficient clears every balanced-polynomial denominator

Status: all-index proof accepted in the [independent arithmetic audit](balanced-integrality-independent-review.md), following the accepted
Rodrigues construction in `balanced-generalized-wilson-rodrigues.md`. No
historical novelty or irrationality conclusion is asserted.

Fix m,n>=1, N=2m+1, D=mn+m+1, and let F_(m,n) be the Rodrigues polynomial
from that note. The new arithmetic statement is

    F_(m,n)(x)/n! belongs to Z[x^2],                       (1)

and consequently

    binom((m+1)n+m,n) p_(m,n)(t) belongs to Z[t].           (2)

Thus the denominator of the monic polynomial is bounded by a single explicit
binomial coefficient. The half-integer shifts in the raising operators do
not produce an additional power of two in this normalization.

## 1. Integer polynomial divided differences

For distinct nodes z_0,...,z_n, denote the ordinary divided difference by

    P[z_0,...,z_n]=sum_(k=0)^n P(z_k)/product_(j!=k)(z_k-z_j).

If P is a polynomial, this rational expression extends to a polynomial in
all the nodes. More specifically,

    z^ell[z_0,...,z_n]=0                         if ell<n,
    z^ell[z_0,...,z_n]=h_(ell-n)(z_0,...,z_n)    if ell>=n,                (3)

where h_q is the complete homogeneous symmetric polynomial, defined by

    product_(k=0)^n (1-z_k u)^(-1)=sum_(q>=0)h_q u^q.

Its coefficients are integers. Formula (3) follows, for example, by partial
fractions in u for distinct nodes, and then by polynomial continuation.
Therefore P in Z[z] and all nodes in Z[x] imply P[z_0,...,z_n] in Z[x],
including coincident-node specializations. This is the arithmetic mechanism.

## 2. Two explicit divided-difference representations

For even n=2r define

    z_k=(x+r-k)^2,  k=0,...,n,
    P_(m,n)(z)=z^m product_(a=0)^(r-1)
                         (a^2-z)((a+1)^2-z)^N.            (4)

For odd n=2r+1 define instead

    z_k=(x+n/2-k)^2-1/4
       =x^2+(n-2k)x+((n-2k)^2-1)/4,
    P_(m,n)(z)=(-z)^N product_(a=0)^(r-1)
                  (a(a+1)-z)((a+1)(a+2)-z)^N.            (5)

Every node and both versions of P have integer coefficients. In the odd case
n-2k is odd, so its squared value minus one is divisible by four. Empty
products are one. In both cases deg P=(m+1)n+m.

The exact identity is

    F_(m,n)(x)/n!
       =(-1)^n x^(-2m) P_(m,n)[z_0,...,z_n].             (6)

The right side is a polynomial despite its displayed x^(-2m); this follows
from the derivation below and the already established polynomial Rodrigues
construction. Thus the numerator divided difference vanishes to order at
least 2m at x=0. Combining (3)-(6), division by the monomial x^(2m) retains
integer coefficients. Evenness of F then proves (1).

For example, when n=1, P(z)=(-z)^N and the two nodes are x^2+x,x^2-x.
Equation (6) reduces to

    F_(m,1)(x)=[(x+1)^N-(x-1)^N]/2,

which has integer coefficients and constant term one. For m=n=2,

    F_(2,2)(x)/2!=28x^8+287x^6+654x^4+311x^2+16.

## 3. Derivation from the Rodrigues formula

Use the notation w_s, E and Dcal of the Rodrigues note. The elementary
operator identity

    E(2x g)=2x Dcal g

implies E^n v=2x Dcal^n[v/(2x)]. Repeatedly applying the Pearson identity
E(w_(s+1)f)=-w_s B_s f gives

    F_(m,n)(x)=(-1)^n E^n w_n(x)/w_0(x).                  (7)

For arbitrary function values at the points

    u_k=x+n/2-k,  k=0,...,n,

iteration of Dcal gives n! times the divided difference on the quadratic
nodes u_k^2. This statement is an algebraic identity wherever the nodes are
distinct and the function values exist; no global evenness of the function
is needed. To prove it, compare the n-node lists after the +1/2 and -1/2
shifts. They overlap in all interior nodes, and their two extreme squared
nodes differ by 2(n+1)x. The divided-difference recursion therefore cancels
the denominator 2x in Dcal and introduces the factor n+1. The n=1 statement
is the definition, so induction proves the assertion. Equivalently,

    Dcal^n g(x)=sum_(k=0)^n (-1)^k binom(n,k)
                 (2x+n-2k)g(u_k)/(2x-k)_(n+1).          (8)

A common shift of all squared nodes changes neither the divided difference
nor its weights; this allows the subtraction of 1/4 in (5).

It remains to evaluate w_n(u_k)/(2u_k). The reflection formula gives, as a
meromorphic identity,

    w_0(x)=4 pi^N x^N cos(pi x)/sin(pi x)^N.              (9)

If n=2r, then u_k-x=r-k is an integer. Since N is odd, the trigonometric
ratio in (9) has period one, and gamma recurrence gives

    w_n(u_k)/(2u_k)
      =w_0(x)/(2x^N) (u_k^2)^m
          product_(a=0)^(r-1)(a^2-u_k^2)((a+1)^2-u_k^2)^N
      =w_0(x)/(2x^N) P_(m,n)(z_k).                       (10)

If n=2r+1, first calculate

    w_1(u)=-4 pi^N u(1/4-u^2)^N sin(pi u)/cos(pi u)^N.

Now u_k-x is an integer plus 1/2. The half-period substitution in this
formula, followed by r gamma shifts, yields

    w_n(u_k)/(2u_k)
      =w_0(x)/(2x^N) (1/4-u_k^2)^N
        product_(a=0)^(r-1)
          ((a+1/2)^2-u_k^2)((a+3/2)^2-u_k^2)^N
      =w_0(x)/(2x^N) P_(m,n)(z_k).                       (11)

Here z_k=u_k^2-1/4 turns the half-integer squared differences into exactly
the integer factors in (5). All identities may first be taken at generic x;
meromorphic continuation handles special points afterward.

Substitute (10) or (11) into the finite sum (8), then use (7) and
E^n=2x Dcal^n(1/(2x)). The common factor w_0(x)/(2x^N) cancels and leaves
precisely (6). As both sides are rational functions and the left is an even
polynomial, the identity extends through x=0 and every coincident-node value.

## 4. A uniform exponential coefficient-denominator bound

Let T_(m,n) be the least positive integer for which T_(m,n)p_(m,n) has
integer coefficients. The monic normalization

    p_(m,n)(-x^2)=(-1)^(mn)F_(m,n)(x)/(D)_n

and (1) prove the divisibility, not merely an inequality,

    T_(m,n) divides C_(m,n):=binom((m+1)n+m,n).            (12)

For fixed m, set K_m=(m+1)^(m+1)/m^m. The binomial theorem with probabilities
1/(m+1) and m/(m+1) gives

    binom((m+1)n,n) <= K_m^n.

Also

    C_(m,n)/binom((m+1)n,n)
       =product_(a=1)^m [((m+1)n+a)/(mn+a)]
       <=((m+1)/m)^m.

Therefore

    T_(m,n) <= ((m+1)/m)^m K_m^n.                         (13)

In particular, for the two-component family relevant to the fifth odd index,

    T_(2,n) divides binom(3n+2,n),
    T_(2,n) <= (9/4)(27/4)^n.                            (14)

This controls the polynomial coefficient denominators for every n. It does
not claim a primitive denominator bound for any particular zeta linear form:
those forms introduce additional rational moments and pole denominators, which
must be analyzed separately. The result supplies one previously missing
arithmetic ingredient, not an irrationality proof.

## Subsequent arithmetic advance

The [Andrews grid theorem](balanced-andrews-grid-arithmetic.md), now
[independently audited](balanced-andrews-independent-review.md), proves
an all-integer negative-square divisor and an explicit exponential normalizer
for the pure-zeta rank-one trial forms. Its exponential constant is not shown
to be effective for irrationality. This resolves the formerly missing
factorial cancellation, while leaving the decisive rate comparison open.
