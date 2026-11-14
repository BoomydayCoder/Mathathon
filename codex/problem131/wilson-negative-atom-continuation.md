# Explicit negative atoms beyond the half-shift Wilson interval

Status: root-proposed formula verified and proved below using Bernoulli
polynomials; exact rational checks pass. Separate independent review is pending.
The input is the already proved positive measure and Bernoulli moments in
[logarithm-zeta3-wilson-deformation.md](logarithm-zeta3-wilson-deformation.md).
This is a self-contained derivation of the continuation, not a historical
novelty claim about classical Wilson measures.

## 1. The continued positive Hamburger measure

For |e|<=1/2 let nu_e=nu_|e| be the probability measure on t=y²>=0 with

\[
 d\nu_e(y^2)=
 \left(\frac{\sin\pi e}{\pi e}\right)^2
 \frac{2\pi^2 y\tanh\pi y}{\cosh^2\pi y-\sin^2\pi e}\,dy,
\]

with the established limits at e=0 and |e|=1/2. Its moments are

\[
 M_k(e)=\frac{(-1)^k}{(k+1)e^2}
 [B_{2k+2}(1/2+e)-B_{2k+2}(1/2)].                     \tag{1}
\]

At e=0 use the polynomial limit. For every k, M_k is a polynomial in e²,
so (1) defines the continued moment functional for arbitrary real e.

Fix tau>1/2 and set

\[
 K=\lfloor\tau+1/2\rfloor,\quad \epsilon=\tau-K\in[-1/2,1/2),
 \qquad v_j=\tau-j-1/2\quad(0\le j<K).
\]

Then the continued functional has the explicit positive probability measure
on the real line

\[
 \boxed{\widetilde\nu_\tau=
 \frac{\epsilon^2}{\tau^2}\nu_{|\epsilon|}
   +\sum_{j=0}^{K-1}\frac{2v_j}{\tau^2}\delta_{-v_j^2}.}       \tag{2}
\]

A zero coefficient is omitted. In particular, at half-integer tau the last
v_j is zero and contributes NO atom at zero. If epsilon=0 the continuous
part is absent. All nonzero atomic weights are strictly positive, and their
nodes are distinct and strictly negative. Thus this is a positive Hamburger
measure; it is not a positive measure supported entirely on [0,infinity).

For 1/2<tau<1 this simplifies to exactly one negative atom:

\[
 \widetilde\nu_\tau=
 \frac{(1-\tau)^2}{\tau^2}\nu_{1-\tau}
 +\frac{2\tau-1}{\tau^2}\delta_{-(\tau-1/2)^2}.        \tag{3}
\]

The continuous density in (2) equals the original displayed hyperbolic
formula with e replaced by tau: sin²(pi tau)=sin²(pi epsilon), and its
normalizing coefficient changes by epsilon²/tau². Simply continuing that
hyperbolic density alone therefore misses the negative atoms and the missing
probability mass.

## 2. Bernoulli shift proves every moment and the mass

The two required identities are

\[
 B_d(x+1)-B_d(x)=d x^{d-1},\qquad
 B_{2r}(1/2+\epsilon)=B_{2r}(1/2-\epsilon).
\]

Telescope the first over K consecutive shifts from 1/2+epsilon to 1/2+tau.
For d=2k+2,

\[
 B_d(1/2+\tau)-B_d(1/2+\epsilon)
 =d\sum_{j=0}^{K-1}(\tau-j-1/2)^{d-1}.
\]

Insert this into (1), use its evenness, and obtain

\[
 M_k(\tau)=\frac{\epsilon^2}{\tau^2}M_k(|\epsilon|)
 +\sum_{j=0}^{K-1}\frac{2v_j}{\tau^2}(-v_j^2)^k.       \tag{4}
\]

When epsilon=0 the first summand means zero, with no division by zero.
This proves (2) as a representing measure directly, rather than by analytic
continuation of an integral outside its original domain.
For k=0 the normalization also follows from the finite identity

\[
 \epsilon^2+2\sum_jv_j
 =(\tau-K)^2+K(2\tau-K)=\tau^2.
\]

## 3. The full Hilbert Gram identity survives except for integer degeneration

Let

\[
 D_k(t)=\prod_{r=0}^{k-1}(t+(r+1/2)^2),\qquad
 d_k(\tau)=\prod_{r=1}^k(r^2-\tau^2).
\]

The already proved Bernoulli identity is polynomial in tau²:

\[
 \widetilde\nu_\tau(D_jD_k)=\frac{d_j(\tau)d_k(\tau)}{j+k+1}. \tag{5}
\]

Alternatively, it follows on the original open parameter interval and hence
everywhere by equality of polynomials. For tau not an integer all d_k are
nonzero. Thus Q_k=D_k/d_k defines an isometry on the entire polynomial space,

\[
 \widetilde\nu_\tau(Q_jQ_k)=\frac1{j+k+1}.             \tag{6}
\]

Some d_k can now be negative; that does not invalidate the Gram identity or
positivity. The continued measure has a nonzero positive continuous part at
all noninteger parameters, so its polynomial-square form is positive definite.

At the positive integer tau=N, formula (2) instead becomes

\[
 \boxed{\widetilde\nu_N=
   \sum_{r=0}^{N-1}\frac{2r+1}{N^2}\delta_{-(r+1/2)^2}.}      \tag{7}
\]

It has exactly N positive-weight atoms. Its moment matrix has rank N;
more precisely the moment matrix of degrees 0,...,d has rank min(d+1,N), by
the Vandermonde matrix at these distinct nodes. Equation (6) remains valid
for 0<=j,k<N. There is no full infinite isometry: d_N=0, and the multiplication
coefficient a_(N-1)=N²-tau² vanishes. The annihilating polynomial D_N vanishes
at every atom, so (5) with a factor D_j, j>=N, vanishes on both sides.
This describes a finite quotient, not a division by a zero Newton denominator.

## 4. A real differential explanation of the first atom

The Newton multiplication rule transports to

\[
 T_\tau=-D[x^2(1-x)D]+(1-\tau^2)x-1/4
\]

on polynomials in L²[0,1]. For tau>1/2 the function f(x)=x^(tau-1) is in L²,
and direct differentiation gives the formal identity

\[
 T_\tau f=-(\tau-1/2)^2f.                             \tag{8}
\]

Indeed the x^tau coefficient cancels, while the x^(tau-1) coefficient is
-tau(tau-1)-1/4. The normalized squared overlap with the constant function is

\[
 \frac{(\int_0^1f\,dx)^2}{\int_0^1f^2dx}
 =\frac{(1/\tau)^2}{1/(2\tau-1)}
 =\frac{2\tau-1}{\tau^2},                            \tag{9}
\]

exactly the weight of the most negative atom in (2). This is a consistency
and motivation calculation. The existence and weight of the measure's atom
were already proved algebraically; no claim that f belongs to any selected
closed-operator domain is needed to establish them.

### Every atom has an explicit Jacobi counterpart

There is a complete finite dictionary, not only the first atom. For every
integer j>=0 with j<tau-1/2, put

\[
 s=\tau-j-1,\quad \alpha=2s+1>0,\quad
 v_j(x)=P_j^{(\alpha,0)}(1-2x),\quad f_j(x)=x^s v_j(x).
\]

The Jacobi differential equation is

\[
 x(1-x)v_j''+[2s+2-(2s+3)x]v_j'
                      +j(2\tau-j)v_j=0.
\]

Substitution therefore gives

\[
 T_\tau f_j=-(s+1/2)^2f_j
          =-(\tau-j-1/2)^2f_j.                        \tag{9a}
\]

Both required integrals have especially simple exact values:

\[
 \boxed{\int_0^1 f_j dx=1/\tau,\qquad
        \int_0^1 f_j^2dx=1/(2\tau-2j-1).}            \tag{9b}
\]

Here is a proof without a spectral norm formula. The shifted Jacobi Rodrigues
formula is

\[
 v_j(x)=\frac{x^{-\alpha}}{j!}
       D^j[x^{j+\alpha}(1-x)^j].
\]

Integrate it against x^s. After j integrations by parts, differentiating
x^(s-alpha)=x^(-s-1), the result is

\[
 \int x^s v_j dx=\frac{(s+1)_j}{j!}\,\mathrm B(s+1,j+1)
               =\frac1{s+j+1}=\frac1\tau.
\]

For the squared norm, write v_j=v_j(0)+xw, where deg w<j (w=0 for j=0).
Rodrigues orthogonality gives integral x^alpha v_j w=0. Hence

\[
 \int x^{\alpha-1}v_j^2dx
 =v_j(0)\int x^{\alpha-1}v_jdx
 =\frac{(\alpha+1)_j}{j!}\,\mathrm B(\alpha,j+1)
 =\frac1\alpha.
\]

The second equality again follows by integrating Rodrigues by parts, now
against x^(-1). Every zero-endpoint boundary term vanishes: its smallest
power is x^(s+1) in the first calculation and x^alpha in the second. Both
exponents are positive. At one, the factor (1-x)^j supplies the required
vanishing orders.

Equations (9a–b) give the squared normalized overlap
(2tau-2j-1)/tau², precisely the corresponding atom weight in (2).
The functions for distinct j are mutually orthogonal: Green's integration
identity has zero boundary terms, since at zero their exponent is at least
s_i+s_j+1>0, and their formal eigenvalues are distinct. At integer tau=N
these are N ordinary polynomials of degree at most N-1, giving a direct
finite-dimensional realization of (7). At noninteger tau, these calculations
remain explicit integrals and differential identities; no closed-domain
assertion is imported implicitly.

## 5. Why the boundary is sharp even at the polynomial-functional level

One representing measure with negative atoms does not, by itself, exclude a
different representing measure on [0,infinity). The following direct test
avoids that logical gap, and avoids any moment-determinacy or spectral theorem.

For 1/2<tau<1 put b=1-tau in (0,1/2) and take finite polynomials

\[
 f_L(x)=\sum_{k=0}^L\frac{(b)_k}{k!}(1-x)^k.
\]

On 0<x<1 these converge to f=x^(-b)=x^(tau-1), with
0<f_L<=f and |f_L'|<=|f'|, because the binomial-series coefficients are positive.
Both f² and x²(1-x)f'² are integrable at zero precisely because tau>1/2.
Dominated convergence therefore gives

\[
 \langle T_\tau f_L,f_L\rangle\longrightarrow
 \int x^2(1-x)f'^2+\int[(1-\tau^2)x-1/4]f^2
 =-\frac{(\tau-1/2)^2}{2\tau-1}<0.                   \tag{10}
\]

For the equality one can substitute elementary powers directly, or integrate
(8): its boundary term is (tau-1)x^(2tau-1)(1-x), which vanishes at both ends.
For tau>=1, the constant polynomial is already a negative witness since
<T_tau 1,1>=1/4-tau²/2<0. Consequently for EVERY tau>1/2 some polynomial p
has continued moment functional \widetilde nu_tau(t p(t)²)<0, using the
isometry (6). At integer N>=1 the constant witness suffices and is also
visible directly from M_1(N)=1/4-N²/2.

Thus the continued moments have no representing positive measure supported
on [0,infinity) beyond the original interval. They do retain the explicit
positive Hamburger representation (2). The obstruction concerns support,
not positivity of the ordinary polynomial-square Gram matrix.

Checks: [check_wilson_negative_atoms.py](../computations/check_wilson_negative_atoms.py)
uses exact rational Bernoulli values to test all moments and Gram identities,
integer ranks, finite negative witnesses, and the overlap constants. Finite
checks supplement the all-parameter proof, rather than replacing it.
