# An unexpected uniform limit: optimal cuts reconnect every odd zeta form to zeta(3)

Status: independent derivation from the contour identity developed in
partial-truncation-analysis.md. General theorem below has been independently audited; see
odd-zeta-independent-review.md Section10 and partial-truncation-analysis.md. Historical novelty unassessed. This is an analytic theorem,
not an irrationality result for higher odd values.

Use R_n=P_n² as before, F_n=-sum_(t>=1)R_n'(t), and for integer 0<=k<=n put

    E_(n,k,m)= -1/(2m)! sum_(t=k+1)^infinity R_n^(2m-1)(t).

For m>=2 these isolate zeta(2m+1) with a rational coefficient, but the cut changes
the denominator. For m=1, all removed first derivatives vanish, so E_(n,k,1)=F_n/2.

**Theorem.** Fix m>=1. If k_n=n/sqrt2+O(1) is an integer sequence, then

    E_(n,k_n,m)/F_n -> (-1)^(m+1)*(2*pi)^(2m-2)/(2m)!.

The convergence has error O_m(n^-1/2). In particular all these optimally cut
forms have the Apéry exponential rate rho=(sqrt2-1)^4, and their eventual sign
alternates with m. For m=2 the ratio is -zeta(2); for m=3 it is pi^4/45.

## Exact contour decomposition

Set c=k+1/2 and

    U_n(z)=Gamma(z)^2*Gamma(n+1-z)/Gamma(n+1+z).

Reflection gives R_n(z)=sin(pi*z)^2 U_n(z)^2/pi². Though U_n has poles to the
right, the original R_n has none there; perform contour summation on R_n first.
Closing the upward vertical line Re z=c clockwise to the right gives

    E_(n,k,m)=1/[(2m)!*2*pi*i] integral_(c-i∞)^(c+i∞)
                                        pi*cot(pi*z)*R_n^(2m-1)(z) dz.

This follows from residues at the integers greater than c. Rational decay of
R_n and its derivatives justifies the limiting contour, with sides chosen away
from integer poles of cot. Integrating by parts 2m-1 times introduces a minus
sign; boundary terms vanish by rational decay and boundedness of cot and its
derivatives on the vertical line. This yields

    E= -1/[(2m)!*2*pi*i] integral pi*[D^(2m-1)cot(pi*z)]*R_n(z) dz.

Define integer coefficients C_(m,j) by

    D^(2m-1)cot(pi*z)=pi^(2m-1)*sum_(j=1)^m C_(m,j)csc(pi*z)^(2j).

They exist because D²(y^j)=pi²[-4j²y^j+2j(2j+1)y^(j+1)] for
y=csc(pi*z)^2. Thus C_(1,1)=-1 and

    C_(m+1,j)=-4j² C_(m,j)+2(j-1)(2j-1) C_(m,j-1),

with missing indices zero. In particular C_(m,1)=(-1)^m4^(m-1).

On z=c+iy, csc(pi*z)^2=sech(pi*y)^2. Substituting the reflection identity gives

    E_(n,k,m)= -pi^(2m-2)/(2m)! * (1/(2*pi))
       *sum_(j=1)^m C_(m,j)
           integral_R sech(pi*y)^(2j-2)*U_n(c+iy)^2 dy.       (1)

The j=1 integral divided by 2pi is F_n. This is the same contour argument with
one derivative: removing the first k integer points makes no change because
R_n' vanishes at all its first n roots. The remaining terms have integrable
sech weights. Complex squares in (1) must not be replaced by absolute squares.

## Bound the correction, and locate the best cut

The beta factorization

    U_n(z)=B(z,z)*B(n+1-z,2z)

and beta-integral triangle inequalities give |U_n(c+iy)|<=U_n(c) for
0<c<n+1. Therefore the sum of terms j>=2 in (1) is O_m(U_n(c)^2).

Uniform Stirling estimates at c=n/sqrt2+O(1) give

    U_n(c)^2=O(n^-2*rho^n).

For completeness the exponential is obtained by writing x=c/n:

    log U_n(c)^2
      =-n*phi(x)+O(log n),
    phi(x)=2[(1+x)log(1+x)-(1-x)log(1-x)-2xlogx].

The unique maximum of phi is at1/sqrt2, with value -logrho. The factorial
square-root factors give n^-2 in U_n(c)^2. The sharper computation is in
apery-asymptotics.md and partial-truncation-analysis.md.

The already proved Casoratian and binomial asymptotics imply

    F_n~C_F*rho^n*n^-3/2,
    C_F=2^(3/4)*pi^(3/2)/(1+sqrt2)^2>0.

Divide (1) by F_n: the j=1 term has coefficient
-C_(m,1)*pi^(2m-2)/(2m)!, and every other term is O_m(n^-1/2).
This proves the theorem.

## Motivation and limitation

A cut near n/sqrt2 balances the factorial products in the removed-root
correction. The same saddle controls Apéry's integer coefficients. The contour
identity explains why this agreement is structural: its leading term is exactly
a multiple of the original zeta(3) form, with an explicitly smaller correction.
The cotangent derivatives force the even powers of pi in the ratio.

This improves the analytic rate from16^-n for the full cut to rho^n near the
balanced cut. It does not preserve the original denominator bound. The
partial-truncation arithmetic note proves new-prime contributions and gives
exact finite tests; no successful irrationality squeeze is claimed here.
