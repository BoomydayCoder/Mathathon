# Independent review: the cubic balanced form is an adjacent Apéry form

Verdict: accepted after correcting one explanatory sign typo in the Newton
conversion. Root reviewed cubic-balanced-adjacent-apery-identity.md. Its
ordinary-Wilson orthogonality input remains the separately accepted Rodrigues
theorem; this identification does not remove that input's contour dependency.

The corrected conversion is q_(0,j)(x)=(j!)² Q_j(-x²), without an extra
(-1)^j. The coefficient(-n)_j already supplies the alternating sign. With
this correction the raising identity and its finite coefficient induction
give p_n/p_n(0)=U S_n exactly, including n=0. The m=1 orthogonality and mass
then determine the functional identity on the complete polynomial basis S_n.
Evaluation at a=n+1 rearranges the four binomial factors term by term into
the already established mixed Apéry coefficient. This alone would not prove
the numerator, so the later transform is necessary.

T_0=H+xD has the displayed monomial action and UT_0=tU. Its ordinary formal
adjoint is H-xD-1. The ladder x(L_a+L_(a-1))'=a(L_a-L_(a-1)) gives the
stated derivative identity for the kernel averages; combining their two H
source equations yields the adjoint equation with source1. Thus h_a has
source(1-x)/2. The finite Christoffel coefficient gives its zero value at1,
and the two cubic tail identities give integral h_a=F(a), including the
1/(2a) correction. Polynomial division and the exact Green bracket therefore
match the original rational numerator, without an irrationality assumption.
At0 the x²h' and xh factors tend to zero; at1 the remaining xrh term vanishes
because h_a(1)=0. All displayed boundary conditions are sufficient.

In the final mixed-integral calculation, the first integration has no y=1
term since P(1)=0. The second has exactly the boundary term
4(-1)^(a-1)/a³. The ladder D_a'=-aP/(1-x) gives -(J-F), with this sign.
Multiplication by a²/8 cancels the rational correction -integral C_a/4
because integral C_a=2(-1)^(a-1)/a. Symmetry gives F-G=4I_(a,a-1).
This verifies the full normalization a²I_(a,a-1)/2 and hence the claimed
(a²/6)(epsilon_a+epsilon_(a-1)). The n=0 value is zeta(3)-1.

The form is therefore positive for every n, not only the even indices allowed
by the general higher-zeta transport. Its Apéry exponential and arithmetic
follow from an exact identity. The sign correction is incorporated in the
proof and explicitly tested in the author's checker; it does not alter the
final theorem or the previously computed polynomial family.
