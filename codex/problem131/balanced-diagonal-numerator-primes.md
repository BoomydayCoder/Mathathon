# Large primes cancel from the diagonal associated numerator

Independent review: accepted in [balanced-diagonal-independent-review.md](balanced-diagonal-independent-review.md).

Status: all-index proof drafted for independent review. The only additional
classical arithmetic input is the primewise consequence of the von
Staudt–Clausen theorem: for even k>=2, B_k has at most a simple p in its
denominator, and such a p can occur only if p-1 divides k.
No irrationality or optimal-denominator assertion is made. The prime cancellation
is an application of the existing general polynomial-moment threshold in
`odd-inner-prime-shift.md`, section 6; the derivative proof is included here
to make its use in the new associated normalizer explicit.

Fix m,n>=1, set s=2m+1, a=n+1, d=mn and

    R(x)=p_(m,n)(-x^2)/p_(m,n)(0).

The Andrews diagonal sum proves R(a) is an integer. Define

    alpha=a^(2m)R(a) in Z,
    P(x)=x^(2m)[R(x)-R(a)]/(x-a),
    b=alpha H_a^(s)+R(a)/(2m)-1/(2a)
                       +BB P^(2m-1)/(2m)!.              (1)

Then the exact normalized moving-pole form is

    mu_s[(p(t)/p(0))/(t+a^2)]=alpha zeta(s)-b.            (2)

For even n it is strictly positive. The arithmetic results below hold for
all n.

## 1. Prime support from the derivative order

The coefficient theorem F/n! in Z[x^2] and the constant F(0)=(n!)^(2m+1)
show that R has coefficients in Z[1/n!]. Consequently, for every prime p>n,
R and P have p-integral coefficients; division by the monic polynomial x-a
introduces no new denominators. Also deg P<=L=2mn+2m-1.

Put r=2m-1. A monomial c x^(k+r) of P contributes

    c (k+1)(k+2)...(k+r) B_k                             (3)

to BB P^(r), where k<=L-r=2mn. Consider a prime p>n+1. The k=0 term is
integral at p; the k=1 term has no denominator at this odd prime. For even
k>=2, the only possible Bernoulli p-denominator requires k=ell(p-1).
Because p-1>=n+1,

    1<=ell<=floor(2mn/(p-1))<=2m-1=r.

The factor k+ell=ell p occurs in the derivative product in (3), so it cancels
that possible simple p-denominator. Every monomial contribution is p-integral.
Thus

    BB P^(2m-1) is p-integral for every p>n+1.             (4)

At a prime p=n+1, all coefficients of P are still p-integral, so the same
classical Bernoulli theorem gives the weaker but sufficient estimate

    v_p(BB P^(2m-1))>=-1.                                (5)

No claim that this one remaining possible denominator always survives is
made. The argument is a derivative-order cancellation, not a finite-prime
experiment.

The other terms of (1) have no variable primes above a. When p=a is prime,
the harmonic term alpha H_a^(s) has valuation at least -1, since alpha
contains a^(2m) and H_a^(2m+1) has at worst p^(-2m-1). The correction 1/(2a)
also costs at most one p beyond the fixed factor 2. Hence

    (2m)! b is p-integral for p>n+1,
    v_p((2m)! b)>=-1 when p=n+1 is prime.                (6)

The fixed factorial is retained to cover the normalization in (1).

## 2. A smaller explicit diagonal normalizer

Let

    L=2mn+2m-1,
    Lambda_(L;n)=product_(p<=n, p prime) p^floor(log_p L),
    U=(2m)!(n+1) d_n^(2m-2) Lambda_(L;n)^(2m+1).         (7)

Then U is an integer and

    U b belongs to Z.                                    (8)

At primes p<=n, this follows from the previously proved Int(Z) clearing
bound: (2m)!d_n^(2m-2)d_L^(2m+1) clears b, and Lambda_(L;n) retains exactly
the p-adic part of d_L at those primes. At p>n, use (6); the single factor
n+1 clears its sole possible variable prime above n, and (2m)! covers the
fixed normalization. This proves (8) prime by prime.

For even n, therefore,

    U alpha zeta(2m+1)-U b>0

has integer coefficients and an exact positive-integral representation.
This is a clearing statement, not an assertion that its positive value tends
to zero.

## 3. The larger lcm index now costs only a subexponential correction

Write Lambda_(L;n)=d_n E_(m,n). If a prime occurs in E_(m,n), then some
power p^e lies in (n,L] while p<=n, so e>=2 and p<=sqrt(L). Each such prime
contributes at most log L to log E. The elementary bound on their number by
sqrt(L) therefore gives

    log E_(m,n)<=sqrt(L) log L.

Consequently, for fixed m,

    log U=(4m-1) log d_n+O_m(sqrt(n)log(n+1)).             (9)

This improves the leading exponential cost compared with retaining the
whole lcm d_(2mn+2m-1). Any elementary bound log d_n<=c n+o(n) immediately
gives log U<=(4m-1)c n+o(n). In particular the previously proved period-six
lcm estimate can be used without invoking the prime number theorem.

For m=2, the bound is

    U=24(n+1)d_n^2 Lambda_(4n+3;n)^5,
    log U=7 log d_n+O(sqrt(n)log(n+1)).                   (10)

Its growth is still too large for the currently available error estimates.

## 4. A sharper finite pattern is kept separate

Exact computations for m=1,...,4 and n=1,...,15 satisfy

    (2m)!(n+1)d_n^(2m+1) b in Z.                         (11)

This would remove the extra factor d_n^(2m-2), as well as the small-prime
correction in (7). It is **not proved here**. The large-prime cancellation
(4)-(6) does not control the remaining exponents at primes p<=n, and must
not be presented as proving (11). A genuine small-prime cancellation lemma
or a transformed formula for the associated numerator is still needed.

The positive Andrews sum proves the diagonal zeta coefficient is integral;
it does not automatically prove (11) for the rational coefficient. The
classical Krattenthaler–Rivoal denominator theorem is not being applied
without an identification of its linear forms with (1).
