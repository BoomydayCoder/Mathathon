# Why the polynomial choice matters: even the zeta(3) kernel is indefinite

Status: independently derived elementary theorem. Historical novelty unassessed.
This strengthens a conceptual warning in the higher-log exploration: positivity
in Beukers's proof comes from the specific Legendre construction and its
integration-by-parts identity, not from generic positivity of the logarithmic
kernel as a quadratic form.

For integer q>=0 define the quadratic form on real polynomials

    Q_q(P)=integral_0^1 integral_0^1
             P(x)P(y)*[-log(xy)]^q/[q!*(1-xy)] dx dy.

The integrals converge absolutely. A pointwise positive kernel need not give a
positive quadratic form when P changes sign.

**Theorem.** For every integer q>=1, Q_q is indefinite already on the
2-dimensional polynomial space span{1,x}. In particular this holds for q=1,
the zeta(3) kernel, as well as q=3, the zeta(5) kernel.

**Proof.** Put s=q+2 and z=zeta(s). The matrix in basis {1,x} is

    [[(s-1)*z, 1], [1, (s-1)*(z-1)]].

The diagonal entries follow by differentiating the Mellin product 1/t² q times
and summing. The off-diagonal entry is 1: the product 1/[t(t+1)] is
1/t-1/(t+1), so its qth derivative divided by (-1)^q*q! telescopes after
summing t>=1. Thus its determinant is

    d_s=(s-1)²*zeta(s)*(zeta(s)-1)-1.

For s=3, a decreasing-integral tail gives

    zeta(3)<sum_(j=1)^5 j^-3+integral_5^infinity x^-3 dx
           <193/160.

Hence 4*zeta(3)*(zeta(3)-1)<4*193*33/160²=25476/25600<1.
For s=4,

    zeta(4)-1 <1/16+1/81+integral_3^infinity x^-4 dx
              =1/16+2/81<1/10.

Therefore 9*zeta(4)*(zeta(4)-1)<99/100<1. Finally, for s>=4,

    zeta(s+1)-1 <=(zeta(s)-1)/2,
    zeta(s+1)<zeta(s).

The positive quantity (s-1)²*zeta(s)*(zeta(s)-1) consequently decreases as s
increases from 4: its ratio to the previous term is less than
(s/(s-1))²/2<=8/9. Thus d_s<0 for every s>=3.
The symmetric 2-by-2 matrix has negative determinant and positive diagonal,
so one positive and one negative eigenvalue. QED.

An explicit rational example for q=1 is P(x)=1-(5/2)x:

    Q_1(P)=(29*zeta(3)-35)/2<0,

using the same bound193/160<35/29. In contrast the shifted Legendre polynomial
L_1(x)=1-2x gives Q_1(L_1)=10*zeta(3)-12>0, whose sign follows from the
positive integral transformation or the recurrence argument.

**Design consequence.** It is unsafe to replace Legendre polynomials by
arbitrary integer polynomials and retain positivity merely because the kernel
is pointwise positive. Orthogonality, endpoint vanishing, and the exact
transformation are substantive ingredients. The failure is already present at
zeta(3), so it is not exclusively a higher-zeta pathology.

For q=0 the situation changes: 1/(1-xy)=sum_(k>=0)x^k y^k is a positive
semidefinite kernel, since its quadratic form is sum_(k>=0)(integral x^kP(x)dx)².
Inserting a total logarithm destroys that Gram representation. This is another
precise reason why a seemingly small change in an auxiliary kernel can change
the proof mechanism.
