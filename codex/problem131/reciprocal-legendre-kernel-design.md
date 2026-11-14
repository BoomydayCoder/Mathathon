# Why this logarithmic kernel: reciprocal Legendre geometry and a uniqueness lemma

Status: root-reviewed elementary supplement to the [kernel-energy proof](short-kernel-energy-zeta3-proof.md). Exact symbolic checks are supplied in the companion checker. No historical novelty is asserted. This is motivation within a specified design, not uniqueness among all irrationality proofs.

The operator and logarithmic kernel need not be guessed independently. Legendre's equation outside the interval, carried back by the reciprocal coordinate, gives exactly the operator H used in the proof. Requiring a product kernel that intertwines the two equations, stays regular at 1, and has unit endpoint flux then determines K(v)=-log(v)/(1-v) uniquely.

## 1. Reciprocal coordinates explain the extra factor x

The ordinary shifted Legendre operator on (0,1) is
\[
 J=-D_x[x(1-x)D_x].
\]
On the exterior interval (1,infinity), use its positive differential expression
\[
 \mathscr A=-D_z[z(z-1)D_z].
\]
The reciprocal change of variable preserves the L^2 norm when accompanied by one factor of z:
\[
 (Vf)(z)=\frac{f(1/z)}z,\qquad
 \int_1^\infty|Vf(z)|^2dz=\int_0^1|f(x)|^2dx.
\]
Direct differentiation proves
\[
 \boxed{VHV^{-1}=\mathscr A,\qquad
 H=-D_x[x^2(1-x)D_x]+x.}                                  \tag{1}
\]
Thus the cubic coefficient x^2(1-x) is the reciprocal image of the ordinary quadratic Legendre coefficient. The added potential x comes from differentiating the norm-preserving factor. This statement concerns explicit differential expressions and the displayed change of variables, not an assumed self-adjoint extension.

## 2. The product-kernel requirements have a unique answer

Seek a twice continuously differentiable function k on (0,1) such that
\[
 (H_x+J_y)k(xy)=0 \quad(0<x,y<1).                         \tag{2}
\]
With v=xy, direct calculation reduces the left side to
\[
 x[(v^2-v)k''(v)+(3v-1)k'(v)+k(v)]
 =x\frac d{dv}\left(v\frac d{dv}[(v-1)k(v)]\right).
\]
Therefore every solution has the form
\[
 k(v)=\frac{c\log v+d}{v-1}.
\]
Regularity at v=1 forces d=0. The normalization
\[
 \lim_{v\downarrow0}v k'(v)=-1
\]
then forces c=1. Consequently
\[
 \boxed{k(v)=K(v)=-\frac{\log v}{1-v}.}                   \tag{3}
\]
This is the exact uniqueness statement. The endpoint normalization means that integrating the Legendre equation against K supplies the source +1 in (H+n(n+1))k_n=1. It fixes the normalization used by the positive energy proof, instead of requiring a later adjustment.

## 3. Under reciprocity it is the divided difference of log

Apply V in the x variable to the product kernel:
\[
 \boxed{\frac1zK(y/z)=\frac{\log z-\log y}{z-y}
          =\int_0^\infty\frac{du}{(z+u)(y+u)}.}           \tag{4}
\]
The elementary integral follows by partial fractions. The kernel is therefore a logarithmic divided difference on the two complementary Legendre intervals. Equation (2) becomes
\[
 (\mathscr A_z+J_y)\frac{\log z-\log y}{z-y}=0.
\]

There is also a short general explanation. Let p(s) be a quadratic polynomial, let D=-partial_s[p(s)partial_s], and let f be a smooth function. For its divided difference F(z,y)=(f(z)-f(y))/(z-y), writing h=p f' gives
\[
 (D_z-D_y)F
 =-\frac{h'(z)+h'(y)}{z-y}
       +2\frac{h(z)-h(y)}{(z-y)^2}.                       \tag{5}
\]
The terms containing f(z)-f(y) cancel because p is quadratic. If h has degree at most two, the right side vanishes by the exact trapezoidal identity. For p(s)=s(s-1), f(s)=log s gives h(s)=s-1. Since D_y=-J_y on (0,1), (5) proves the claimed kernel equation again.

This explains both useful descriptions of the kernel: differentiating the rational Mellin square gives log(xy), while reciprocal Legendre symmetry gives the divided difference of log. Equation (4) connects these descriptions directly and also supplies the Stieltjes representation used to obtain the Beukers cube.

## 4. Mellin transformation explains the consecutive rational poles

For a polynomial f, let Mf(s)=integral_0^1 x^(s-1)f(x)dx. Ordinary integration by parts, initially for Re(s)>1, gives
\[
 \boxed{M(Jf)(s)=s(s-1)Mf(s)-(s-1)^2Mf(s-1),}           \tag{6}
\]
\[
 \boxed{M(Hf)(s)=s^2Mf(s+1)-s(s-1)Mf(s).}             \tag{7}
\]
These extend as rational identities. If f=L_n, the eigenvalue n(n+1) in (6) forces
\[
 \frac{ML_n(s)}{ML_n(s-1)}
       =\frac{(s-1)^2}{(s-n-1)(s+n)}.
\]
The rational solution with ML_n(s)~(-1)^n/s at infinity is exactly
\[
 \boxed{ML_n(s)=(-1)^n
   \frac{\prod_{j=1}^n(s-j)}{\prod_{j=0}^n(s+j)}.}        \tag{8}
\]
To justify uniqueness, the quotient of two nonzero rational solutions is a rational function with period one; such a function has no pole (otherwise all integer translates would be poles), hence is a periodic polynomial and is constant. The stated behavior at infinity fixes that constant. Alternatively, polynomial degree bounds, the Legendre orthogonality zeros, and the consecutive Mellin poles already force the same product.

Thus the rational interpolation product, the Newton multiplication rule, and the logarithmic integral are three coordinate descriptions of the same Legendre structure. The actual irrationality contradiction still needs the explicit cubic-lcm arithmetic and growth comparison proved in the standalone note.
