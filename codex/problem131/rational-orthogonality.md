# Prescribed poles, Blaschke products, and a rational Gram comparison

Status: exact elementary rational-function theorems reconstructed independently. Classical complex analysis (residues and maximum modulus) is used. These are candidate analytic tools for simplifying zeta determinant estimates; no successful new irrationality estimate is asserted.

## 1. An exact extremal principle explains the Apéry product

For positive numbers a_1,...,a_n put

    B_n(z)=product_(j=1)^n (z-a_j)/(z+a_j).

It has prescribed zeros in the right half-plane, prescribed reflected poles in the left half-plane, and

    |B_n(iy)|=1, y real;    |B_n(z)|<=1, Re z>=0.

The first equality follows factor by factor; the second follows from |z-a|<=|z+a| for Re z>=0 and a>0.

**Rational extremal theorem.** Let F be a rational function with no poles in the closed right half-plane, bounded at infinity, with zeros at each a_j counted with multiplicity. Suppose |F(iy)|<=1 for real y. Then

    |F(z)|<=|B_n(z)|, Re z>=0.

Proof: F/B_n has removable singularities at the a_j and no other poles in that half-plane. Its boundary modulus is at most1. At infinity B_n→1; bounded rational F has a finite limit of modulus at most1, since the same limit is approached along the imaginary axis. Maximum modulus on expanding right half-disks therefore gives |F/B_n|<=1. Equality is attained by F=B_n. QED.

Thus B_n is not merely a useful product guessed from factorials: it is the largest possible interior modulus under these rational boundary and zero constraints. Smallness at interior real points is forced by the prescribed zeros, without paying any boundary-modulus cost.

For a_j=j,

    R_n(t)=[B_n(t)/t]^2

is exactly Zudilin's auxiliary rational function. For t>n, the elementary inequality log((1-u)/(1+u))<=-2u gives

    |B_n(t)|<=exp(-n(n+1)/t),
    R_n(t)<=t^-2 exp(-2n(n+1)/t).

This recovers the suppression mechanism used in the derivative-tail estimates. The boundary pole at0 means B_n(z)/z is not itself an ordinary square-integrable Hardy function. One must not suppress this qualification. The regularized B_n(z)/(z+ε), ε>0, has squared boundary norm 1/(2ε) under the normalization below; the limit ε↓0 diverges.

## 2. An explicit rational orthonormal system

For boundary functions use

    <f,g>_H=(1/(2π))integral_(-infinity)^infinity f(iy)overline(g(iy))dy.

For any positive a_1,a_2,... define

    e_j(z)=sqrt(2a_j)/(z+a_j)
                 product_(k<j)(z-a_k)/(z+a_k).

Then <e_i,e_j>_H=δ_ij. Repeated a_j are allowed.

**Proof of norms.** Every product factor has modulus1 on the boundary, so

    ||e_j||²=(1/(2π))integral 2a_j/(y²+a_j²)dy=1.

**Proof of orthogonality.** All coefficients are real, so conjugation at iy replaces z by -z. If j>i, the shared product factors cancel in e_j(z)e_i(-z). The remaining factor (z-a_i)/(z+a_i) cancels the denominator a_i-z, up to a minus sign. What remains is a rational function with poles only in the left half-plane and decay O(z^-2). Its integral along the imaginary axis is zero by closing a contour in the right half-plane. QED.

When the a_j are distinct, each e_j lies in the span of 1/(z+a_1),...,1/(z+a_j), by simple partial fractions. Thus the product is also the explicit outcome of Gram–Schmidt for a prescribed-pole rational space. This is an exact link between the Apéry product shape and rational orthogonality; it does not claim that Apéry's arithmetic proof follows from Hardy-space orthogonality alone.

## 3. Exact rational Gram matrices

The elementary Cauchy identity is

    <1/(z+a),1/(z+b)>_H=1/(a+b), a,b>0.

One proof uses partial fractions on the imaginary axis; another uses a single residue. Differentiating with respect to a and b, with dominated differentiation, gives

    <1/(z+a)^p,1/(z+b)^q>_H
       =binom(p+q-2,p-1)/(a+b)^(p+q-1), p,q>=1.       (1)

Every proper rational function with rational coefficients and poles at negative rational points therefore has an exactly computable rational Gram matrix under this boundary norm. Repeated poles are not an analytic problem.

For the simple basis φ_j=1/(z+a_j), distinct positive a_j, the Cauchy determinant gives

    det(<φ_i,φ_j>_H)
       = product_(i<j)(a_i-a_j)^2
          /[product_i(2a_i) product_(i<j)(a_i+a_j)^2]. (2)

This can also be derived by the preceding explicit orthogonalization. Equations (1) and (2) supply exact norms and determinants, rather than numerical guesses for a weighted capacity.

## 4. A usable comparison theorem for the zeta Hankel matrix

Fix odd s>=3 and the positive weight w_s from `odd-moment-generalization.md`. Write

    P_m(z)=product_(j=1)^m(z+j),
    D_m(y²)=|P_m(iy)|².

Consider the single-zeta matrix

    G[i,j]=integral_0^infinity y^(2i+2j)
                 D_N(y²)^r/D_K(y²) w_s(y)dy,
    0<=i,j<h,

where K>N>=0 and r>=0 is an integer. Define real rational functions

    f_i(z)=(-1)^i z^(2i) P_N(z)^r/P_K(z).

At z=iy their pairwise products satisfy

    f_i(iy)overline(f_j(iy))
      =y^(2i+2j) D_N(y²)^r/D_K(y²).

The equality includes the signs: (-1)^i(iy)^(2i)=y^(2i).

The f_i may be improper, which prevents direct use of the unweighted boundary norm. Choose any positive rational a and any integer

    L>=max(0, rN+2h-1-K).

Then g_i(z)=f_i(z)/(z+a)^L is proper with O(1/z) decay or better. Its poles remain at negative rational points. Let

    H[i,j]=<g_i,g_j>_H,
    C=sup_(y>=0) w_s(y)(y²+a²)^L.

Exponential decay of w_s makes C finite. All H entries are rational and can be computed by (1). The functions g_i are linearly independent, so H is positive definite.

**Comparison theorem.** In positive-semidefinite matrix order,

    0<G<=π C H,

and consequently

    0<det G<=(π C)^h det H.                            (3)

Proof: for every complex vector v, the quadratic form of G is the half-line integral of w_s times |sum_i v_i f_i(iy)|². Because the f_i(iy) share a common rational factor times real monomials, the corresponding absolute square is even in y. Substituting f_i=(z+a)^L g_i yields

    v*Gv=(1/2)integral_R w_s(|y|)(y²+a²)^L
                               |sum_i v_i g_i(iy)|²dy
          <=π C v*Hv.

Conjugating the inequality by H^(-1/2) bounds every eigenvalue by πC and proves the determinant estimate. QED.

**Important distinction:** the new repeated pole at -a is introduced only in an analytic comparison. It is not inserted into the original rational moment functional whose coefficients involve ζ(s). Thus it does not introduce ζ(s+1),ζ(s+2),... into the arithmetic. This avoids the fixed-compactification obstruction proved in `compactification-obstruction.md`.

## 5. A completely explicit, though crude, scalar bound

Set

    c_s=2(2π)^(s-1)/(s-1)!,
    A_s=sum_(ell>=1) ell^(s-1) exp(-πell).

A_s is a convergent positive constant and can also be evaluated as a rational function of exp(-π) by repeated differentiation of a geometric series.

For 0<y<=1, compare each term ell^(s-1)exp(-u ell) with the integral over [ell,ell+1] after multiplying by exp(u), where u=2πy. This gives

    sum ell^(s-1)e^(-u ell) <= e^u (s-1)!/u^s,
    w_s(y)<=e^(2π)/π.

For y>=1, the inequality 2ell*y>=y+ell gives

    w_s(y)<=c_s A_s y^s exp(-πy).

Using y²+a²<=(1+a²)y² for y>=1 and maximizing y^(s+2L)e^(-πy) over y>0, we obtain

    C<=(1+a²)^L max{
           e^(2π)/π,
           c_s A_s[(s+2L)/(πe)]^(s+2L)
       }.                                              (4)

This makes (3) an elementary fully explicit upper bound. It is intentionally conservative; direct calculus on the true weight may substantially reduce C. The exact H determinant retains the rational orthogonality structure lost by the earlier total-mass compact bound.

## 6. How this becomes a falsifiable irrationality test

Let d clear all rational affine coefficients of the original single-zeta matrix G(X). Then d^h det G(X) is an integer polynomial, and positivity makes its value at ζ(s) nonzero. Combining arithmetic with (3), a sufficient condition is

    h log d+h log(π C)+log det H <=-c K²

along a family with h=O(K), for some c>0. A weaker condition with the left side divided by h tending to minus infinity also suffices.

Every item here is independently meaningful:

* d belongs to the original simple-pole arithmetic, where the Bernoulli cancellation and harmonic estimates apply.
* C is a one-variable weighted maximum depending on the analytic damping choice a,L.
* H is an exactly computable rational Gram determinant with prescribed poles.

The inequality has **not** been established for a family proving ζ(5), ζ(7), or every odd value. In the Fauzan proportions, L is of order K, so C can have large K log K terms. It is necessary to compute their cancellation against log det H and the true arithmetic scale; merely observing that H is small is insufficient.

This comparison is a concrete candidate simplification: replace a many-variable energy estimate by exact rational orthogonality plus an optimized scalar weight bound. Its success or failure can be determined by explicit estimates. It also explains a structural coincidence: the same reflected-pole products underlying Zudilin's auxiliary function form a rational orthogonal system for the denominator factors of the Hankel construction.

## 7. Exact Legendre–Malmquist equivalence

A sharper connection, identified by the parent researcher and independently checked here, uses half-integer poles. Let

    L_n(x)=(1/n!)(d/dx)^n[x^n(1-x)^n],
    M_n(t)=integral_0^1 x^(t-1)L_n(x)dx.

The Mellin identity proved elsewhere in the project is

    M_n(t)=(-1)^n product_(j=1)^n(t-j)/product_(j=0)^n(t+j).

Set t=z+1/2. Then exactly

    M_n(z+1/2)=(-1)^n/(z+n+1/2)
                  product_(j=1)^n (z-(j-1/2))/(z+(j-1/2)).

This is (-1)^n/sqrt(2n+1) times the (n+1)st rational orthonormal function of Section2 for the pole sequence a_j=j-1/2. Thus Legendre polynomials and the half-integer prescribed-pole system are not merely analogous: the shifted Mellin transform takes one to the other exactly.

To see why the half shift appears, substitute x=e^(-u):

    M_n(z+1/2)=integral_0^infinity e^(-zu)
                          [e^(-u/2)L_n(e^(-u))]du.

The bracketed functions have real L² inner product

    integral_0^infinity e^(-u)L_n(e^(-u))L_m(e^(-u))du
       =integral_0^1 L_n(x)L_m(x)dx.

Boundary values of their Laplace transforms are Fourier transforms. Plancherel's identity equates this inner product with the Hardy boundary inner product of the shifted Mellin transforms. The rational orthogonality therefore yields

    integral_0^1 L_n(x)L_m(x)dx=δ_nm/(2n+1),

including the correct norm and the harmless sign (-1)^n. Conversely, classical Legendre orthogonality yields that rational orthogonality under the same transform.

This explains in one exact diagram the recurring ingredients:

    polynomial orthogonality on (0,1)
      -- x=e^(-u), half-density e^(-u/2) --> L²(0,infinity)
      -- Laplace/Fourier boundary transform --> rational Hardy orthogonality.

The integer zeros and poles in the Apéry rational function are the shifted image of this half-integer rational orthogonal system. Squaring the Mellin transform explains the two Legendre factors, while differentiating in its parameter inserts the logarithm. This establishes a concrete common origin for the auxiliary product, Legendre orthogonality, and the rational-function boundary estimates. Arithmetic integrality and the final smallness bound remain additional requirements.
