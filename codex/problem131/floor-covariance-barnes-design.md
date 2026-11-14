# Why the prime-allocation floors turn into superfactorials

Status: elementary identities proved below and independently computable with
rational arithmetic. This derivation supplies the missing motivation for the
Barnes-profile formula used in `barnes-normalizer-without-pnt.md`. It is a
reconstruction of the mathematics, not a historical-priority claim.

## 1. A floor count is a centered interval indicator

For u>=0 and 0<z<1/2, let

    ell_u(z)=floor(u-z)+floor(u+z)+1,
    A_u(z)=ell_u(z)-2u.

Write r={u}, rho=min(r,1-r), and epsilon=1 if r<=1/2, epsilon=-1 otherwise.
Apart from endpoints, which do not affect integrals,

    A_u(z)=epsilon [1_(z<rho)-2rho].                  (1)

This follows by considering whether r lies below or above1/2. In particular,
integral_0^(1/2) ell_u(z) dz=u and integral A_u=0.

Multiplying the two indicators gives the exact covariance

    integral A_u A_v dz
      =epsilon_u epsilon_v [min(rho_u,rho_v)-2rho_u rho_v]
      =[B2({u-v})-B2({u+v})]/2,                     (2)

where B2(t)=t^2-t+1/6. The first equality is immediate from interval lengths;
the second is the same quadratic identity in the four possible endpoint
orders. Setting v=u gives

    integral A_u^2 dz=[B2(0)-B2({2u})]/2.            (3)

The expression min(a,b)-2ab is the covariance kernel of a Brownian bridge on
an interval of length1/2. No probability theorem is required here: the equality
comes from multiplying centered indicators. This explains why a Bernoulli
polynomial naturally records overlap of the prime-residue classes.

## 2. Derive the allocation formula without a piecewise integration table

Let 0<alpha<1, lambda=1-alpha, H=1+2alpha. Define

    T=floor(2Hx), xi=Hx-T/2,
    q=floor(2x), eta=x-q/2,
    b(z)=3 ell_(alpha x)(z).

The source's inner row-allocation expression is

    Gamma_alpha(x)=integral_0^(1/2) (T-b)(T+b-ell_x-5) dz
                      +xi(2T-q-5)+(xi-eta)_+.

For alpha=3/40 and the source's permitted x range this is the actual valuation
envelope. The following identity is algebraic for all the displayed parameters;
we do not infer valid row allocations outside their proved range.

Expanding the product, using integral ell_u=u, and collecting the quadratic
and linear terms gives

    Gamma_alpha(x)=10alpha lambda x^2-5lambda x
                  +3 Cov(x,alpha x)-9 Var(alpha x)+E,
    E=-2xi^2+2xi eta+(xi-eta)_+.

The positive part in E has only two cases. They give

    E=[-B2({2Hx})+B2({2x})-B2({4alpha x})+B2(0)]/4.

Insert (2)--(3). Finally use the exact relation

    B2({v})=v^2-v+1/6-2J(v),
    J(v)=floor(v)v-floor(v)(floor(v)+1)/2.

Every quadratic term cancels, leaving the compact general identity

    Gamma_alpha(x)=-(5-alpha)x
      +(1/2)J(2(1+2alpha)x)-(1/2)J(2x)+(1/2)J(4alpha x)
      -3J((1-alpha)x)+3J((1+alpha)x)-9J(2alpha x).      (4)

At alpha=3/40 this is exactly the six-superfactorial profile in the new
normalizer. The coefficient197/40 is 5-alpha, not a fitted decimal.

## 3. Why J suggests a factorial product

Let B(N)=product_(j=0)^(N-1) j!, for integer N>=0. Count each prime power
separately:

    v_p(B(N))=sum_(a>=1) sum_(k>=1)(N-kp^a)_+
             =sum_(a>=1) p^a J(N/p^a).               (5)

Thus J is not merely a convenient way to rewrite a floor function. It is the
exact prime-power valuation profile of a superfactorial. The source's local
allocation already has the right shape to be compared directly with a finite
product of such objects.

A pointwise majorant sum_j a_j J(c_j x)>=required_cost(x), plus an exact
quadratic cancellation, can therefore control every sufficiently large prime
without estimating how many primes lie in an interval. Stirling estimates the
size of the factorial product directly. The separately proved small-prime
repair handles the omitted higher powers; that repair is essential.

The finite majorant is still a nontrivial design step, supplied by the exact
certificate in `check_sparse_barnes_certificate.py`. Equations (1)--(5) explain
how to look for it: residue overlaps give a covariance, the covariance gives
Bernoulli polynomials, and those polynomials give exact factorial valuations.
This replaces prime-density averaging with an explicit arithmetic object.

The standalone rational checker `check_floor_covariance_barnes.py` integrates
the original floor expression directly, using its z breakpoints, and compares
it with (4) for varied alpha and x. It also checks (2), (3), and (5); these finite
checks supplement the all-parameter algebra above.
