# A short elementary proof with the arithmetic and positive integral visible together

Status: assembled from classical, independently checked ingredients in this
workspace. This is a presentation of Apéry/Beukers/Zudilin methods, not a claim
of a historically new proof. It uses no prime number theorem, recurrence
asymptotics, or complex contour integration. The longer discovery guide explains
the recurrence; this shorter route makes its use unnecessary for irrationality.
The [independent whole-proof review](short-hybrid-independent-review.md) checks
the arithmetic, both improper-integral transformations, and the strict final
comparison.

The design is to find one number with two descriptions. Rational poles expose
its denominators; a positive integral exposes its sign and smallness.

## 1. Orthogonality designs the rational square

Let

    L_n(x)=1/n! * (d/dx)^n [x^n(1-x)^n]
          =sum_(k=0)^n (-1)^k binom(n,k)binom(n+k,k)x^k.

These are shifted Legendre polynomials, with integer coefficients. Integrating
by parts n times, initially for t>n, gives

    integral_0^1 x^(t-1)L_n(x)dx=(-1)^n P_n(t),
    P_n(t)=product_(j=1)^n(t-j)/product_(j=0)^n(t+j).       (1)

Both sides are rational functions of t, so their equality then holds for all
real t>0. This formula is also the orthogonality statement: the first n integer
Mellin arguments vanish. Pairing two such transforms produces

    R_n(t)=P_n(t)^2,

whose first n zeros are double. Thus the rational square is the Mellin image
of two orthogonal polynomials; the zero/pole constraints also force it directly.

Define the absolutely convergent number

    I_n=-sum_(t=1)^infinity R_n'(t).                     (2)

## 2. Its partial fractions give exactly the needed arithmetic

The residues of P_n at its consecutive integer poles are

    c_k=(-1)^(n-k)binom(n,k)binom(n+k,k) in Z.

Therefore

    R_n(t)=sum_(k=0)^n [a_k/(t+k)^2+b_k/(t+k)],
    a_k=c_k^2,
    b_k/2=c_k sum_(j!=k)c_j/(j-k).                      (3)

Put D_n=lcm(1,...,n), with D_0=1. Then a_k is integral and D_n b_k/2 is
integral. Since R_n(t)=O(t^(-2)), comparison of the coefficient of t^(-1)
at infinity gives sum b_k=0. Differentiating and summing (3) consequently yields

    I_n=2(A_n zeta(3)-B_n),
    A_n=sum a_k in Z,
    B_n=sum [a_k H_k^(3)+(b_k/2)H_k^(2)],
    D_n^3 B_n in Z.                                    (4)

Here H_k^(r)=sum_(j=1)^k j^(-r). The possible zeta(2) term is eliminated by
sum b_k=0. Pole differences j-k explain the denominator D_n directly.

## 3. The same number is a positive triple integral

Let K(v)=-log(v)/(1-v), with its continuous value at v=1. The elementary
Stieltjes identity

    K(v)=integral_0^infinity dt/[(1+t)(v+t)], v>0,        (5)

follows by partial fractions. Expanding K(xy) geometrically, or integrating
its monomials first, and using (1) gives

    I_n=integral_(0,1)^2 L_n(x)L_n(y)K(xy) dxdy.          (6)

The logarithm differentiates the Mellin product: it is precisely the derivative
in (2). Bounded polynomial factors and the integrable log kernel justify the
interchanges; alternatively finite partial sums and absolute bounds suffice.

Now (5) turns (6) into a positive expression by two integrations by parts.
First integrate the Rodrigues derivative in x. Next put t=y u, cancelling
all the y powers before integrating the Rodrigues derivative in y. The result is

    I_n=integral x^n(1-x)^n y^n(1-y)^n u^n
             /[(1+yu)^(n+1)(x+u)^(n+1)] dxdydu,

with x,y in(0,1) and u>0. Finally put u=(1-z)/(yz). This proves

    I_n=integral_(0,1)^3
      [x(1-x)y(1-y)z(1-z)]^n/(1-z+xyz)^(n+1) dxdydz.   (7)

All Rodrigues boundary terms vanish because x^n(1-x)^n and its first n-1
derivatives vanish at both endpoints. For fixed other variables each step
is ordinary integration by parts. Here is the absolute-convergence justification
between the two integrations, where one Legendre factor is still signed.
Put \(M_n=\max_{0\le x\le1}|L_n(x)|\). The first integration gives the
nonnegative quantity
\[
J_n(y,t)=\int_0^1\frac{x^n(1-x)^ny^n}{(xy+t)^{n+1}}\,dx
=\int_0^1\frac{L_n(x)}{xy+t}\,dx
\le M_n\int_0^1\frac{dx}{xy+t}.
\]
Thus its integral against \(|L_n(y)|/(1+t)\) is bounded by
\(M_n^2\iint K(xy)\,dx\,dy=2M_n^2\zeta(3)<\infty\).
The substitution \(t=yu\) preserves this absolute convergence, so the
reordering needed for the second integration by parts is justified.
The resulting positive inner integral is bounded by the corresponding
absolute integral before that step; Tonelli then proves finiteness of the
positive triple integral. The final substitution is legitimate, and (7)
proves \(I_n>0\). For \(n=0\), no integrations by parts are needed and the
same convergence bounds apply.

## 4. One AM–GM inequality gives the exponential decay

Make the odds substitution x=X/(1+X), y=Y/(1+Y), z=Z/(1+Z). Define

    Lambda=(1+X)(1+Y)(1+Z)
               *[(1+X)(1+Y)+XYZ]/(XYZ).

The Jacobian in (7) gives exactly

    I_n=integral_(0,infinity)^3 Lambda^(-n-1) dX/X dY/Y dZ/Z.       (8)

Put p=(1+X)(1+Y), q=p^2/(XY), r=sqrt(XY), A=sqrt(p). Then

    Lambda=pZ+(p+q)+q/Z >= (sqrt(p)+sqrt(q))^2
          =(A+A^2/r)^2.

Since X+Y>=2sqrt(XY), A>=1+r. Hence

    Lambda >= (3+2r+1/r)^2 >= (3+2sqrt2)^2
            =lambda:=17+12sqrt2.                    (9)

Equality is attained at X=Y=1/sqrt2 and Z=1+sqrt2. Since I_0=2zeta(3),
which follows at once from (6), equations(8)--(9) imply

    0<I_n<=2zeta(3) lambda^(-n).                     (10)

The same Laurent polynomial has A_n as the constant term of Lambda^n;
that additional identity explains why one often sees reciprocal growth and
decay in the recurrence. It is not needed for the present contradiction.

## 5. An elementary lcm estimate finishes the proof

We need a little better than D_n<=4^n. Put psi(x)=log D_floor(x), and define

    T(x)=log(floor(x)!)-log(floor(x/2)!)
          -2log(floor(x/3)!)+log(floor(x/6)!).

Counting multiples of prime powers in factorials gives

    T(x)=sum_(p prime,a>=1) h(x/p^a)log p,
    h(y)=floor y-floor(y/2)-2floor(y/3)+floor(y/6).

The function h has period6 and takes the values0,1,1,0,0,1 on the six unit
intervals. Thus h>=0 and h=1 on[1,3), so

    psi(x)-psi(x/3)<=T(x).

The elementary factorial estimate log(m!)=m log m-m+O(log(m+2)) gives

    T(x)=[(log2)/3+(log3)/2]x+O(log(x+2)).

Summing over x,x/3,x/9,... proves

    log D_n <= c n+O(log^2(n+2)),
    c=(log2)/2+3(log3)/4.                            (11)

The strict inequality exp(3c)<lambda is elementary. Indeed sqrt2>7/5 gives
lambda>169/5, and

    (exp(3c))^4=64*3^9=1259712
       <169^4/5^4=815730721/625.

Combining(10)--(11), D_n^3 I_n tends to zero while remaining positive.
If zeta(3)=a/b with a,b positive integers, (4) says

    b D_n^3 I_n=2(bD_n^3 A_n zeta(3)-bD_n^3B_n)

is a positive integer. A positive integer cannot tend to zero. Contradiction.

## What each step contributes

Orthogonality produces the zeros. Mellin transforms produce the rational product.
Taking two transforms produces the square, and the logarithmic kernel produces
the derivative. Consecutive poles control denominators, while decay at infinity
cancels the lower zeta value. A Stieltjes integral and two carefully ordered
integrations by parts give positivity. The odds substitution exposes the same
Laurent kernel as the integer sequence, and AM–GM supplies smallness. The lcm
estimate proves that the arithmetic cost fits inside that decay.

The recurrence remains important for explanation, generalization, asymptotics
and efficient computation. The [discovery guide](discovery-guide.md) derives it
from the quadratic Newton grid; the geometric divergence construction gives an
independent advanced route. The exact q-deformation in
[deformed-laurent-period-reduction.md](deformed-laurent-period-reduction.md)
explains why the lower-zeta cancellation is special to q=1.
