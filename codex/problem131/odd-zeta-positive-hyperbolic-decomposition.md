# A positive hyperbolic decomposition of the odd-zeta measures

Status: elementary proof draft, independently reviewable; no historical novelty
claim. This refines the common-density construction and explains one structural
difference between the zeta(3) measure and its higher-odd counterparts. It does
not prove that the latter have no explicit orthogonal polynomials.

Put f(u)=(exp(2u)-1)^(-1), c(u)=csch(u)^2 and s=2m+1. The positive measure used
in our rank-one approximation construction has density

    w_(2m+1)(y)=2*pi^(2m)*y^(2m+1)/(2m)! * f^(2m)(pi*y).

The following exact expansion holds for every integer m>=1:

    f^(2m)(u)=coth(u) * sum_(j=1)^m a_(m,j) c(u)^j,
    a_(m,j)=(2j)!/2 * 4^(m-j) h_(m-j)(1^2,2^2,...,j^2).

Here h_d is the complete homogeneous symmetric polynomial, with h_0=1.
In particular all coefficients are positive integers. The formula follows from

    (d/du)^2[c^j*coth(u)]
      =4j^2 c^j*coth(u)+2(j+1)(2j+1)c^(j+1)*coth(u).

The initial derivative is f''=c*coth(u), so

    a_(m+1,j)=4j^2 a_(m,j)+2j(2j-1)a_(m,j-1).

After division by (2j)!/2, this is precisely the usual recurrence for complete
homogeneous polynomials on the nested alphabet 4,16,...,4j^2. This proves the
closed formula by induction, without needing any special-function identity.
Equivalently, for each fixed j,

    sum_(m>=j) a_(m,j) z^(m-j)
       = ((2j)!/2) / product_(k=1)^j (1-4k^2 z).

This is a formal power-series identity, so no convergence condition is needed.
The coefficients are a scaled central-factorial triangle; that name supplies
context, not a novelty claim.

Since w_3=pi^2*y^3*coth(pi*y)*csch(pi*y)^2, division gives a positive finite
mixture:

    w_(2m+1)/w_3
       = [2*pi^(2m-2)/(2m)!]
         * sum_(j=1)^m a_(m,j) y^(2m-2)csch(pi*y)^(2j-2).

The first three cases are

    w_3/w_3 = 1,
    w_5/w_3 = (pi^2*y^2/3)*(1+3*csch(pi*y)^2),
    w_7/w_3 = pi^4*y^4*(2/45+(2/3)*csch(pi*y)^2+csch(pi*y)^4).

All summands are positive on y>0. The last summand has limit one at y=0;
every earlier summand tends to zero. At infinity only the first summand has
the largest polynomial growth. These two ends give a direct explanation of
why w_s/w_3 is comparable to (1+y^2)^((s-3)/2), the analytic fact underlying
rank-one-universal-error-rate.md.

## Why the Wilson case is structurally special

The standard gamma identities give the exact real-axis identity

    w_3(y)=1/(4*pi) * |Gamma(1+iy)^3 Gamma(iy)/Gamma(2iy)|^2.

Also |Gamma(1+iy)|^2=pi*y/sinh(pi*y). Thus the j-th summand in the mixture
is a positive constant times

    y^(2(m-j))
      * |Gamma(1+iy)^(2j+1) Gamma(iy)/Gamma(2iy)|^2.

More precisely its constant is a_(m,j)*pi^(2m-2j)/[2*pi*(2m)!].
For m=1 there is exactly the four-gamma Wilson weight. For m>1 the same
measure is a positive sum of a polynomial modification of that weight and
weights with six, eight, ..., 2m+2 gamma factors. This gives a concrete
special-function explanation for why direct substitution into the ordinary
four-parameter Wilson formula stops being sufficient.

This observation alone does not establish an obstruction to an explicit formula,
a denominator bound, or irrationality. A positive mixture can sometimes have
unexpected exact structure. The next constructive question is whether the
mixture admits a useful system of multiple orthogonal polynomials or a coupled
Padé problem whose arithmetic is simpler than that of the scalar Gram matrix.
