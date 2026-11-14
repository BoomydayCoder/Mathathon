# Compactifying the odd-zeta moment model: an exact obstruction and a repair

Status: independently derived theorems, using the Hermite moment identity proved in `odd-moment-generalization.md` and classical transcendence of π. No new irrationality theorem is claimed. This investigates whether the elementary compact-capacity criterion could replace the harder determinant estimates.

## 1. The attractive change of variables

Fix an odd integer s>=3 and let μ_s be the positive measure on t>=0 defined by

    integral R(t)dμ_s(t)=integral_0^infinity R(y²)w_s(y)dy.

Its mass is s/12. Its polynomial moments and its simple poles at negative integer squares are rational affine expressions in ζ(s). It is tempting to push μ_s forward through

    x=t/(1+t),    0<=x<1,

and apply `determinant-capacity-recipe.md` to a fixed measure on [0,1]. But x^d=t^d/(1+t)^d introduces repeated poles at -1. Their zeta content is different from that of simple poles.

## 2. Exact repeated-pole formula

For a>0 define

    J_k(a)=integral (t+a²)^(-k)dμ_s(t), k>=1.

The simple-pole identity is

    J_1(a)=a^(s-1)ζ(s,a)-1/(2a)-1/(s-1).

Differentiation under the integral is justified locally in a>0 by domination: the denominator stays bounded away from zero and the weight has finite mass and exponential decay. Hence

    J_k(a)=(-1)^(k-1)/(k-1)! * [(1/(2a))d/da]^(k-1) J_1(a).

Use ∂_a ζ(u,a)=-u ζ(u+1,a). The result is a rational-function combination of

    1, ζ(s,a), ζ(s+1,a), ..., ζ(s+k-1,a).

The coefficient of the last value is exactly

    (s)_(k-1) a^(s-k) / [2^(k-1)(k-1)!],

where (s)_r=s(s+1)...(s+r-1). To obtain that term, every differentiation must hit the zeta factor; differentiating any other factor loses at least one rise in its index. This proves that the displayed top coefficient cannot cancel.

For example,

    J_2(a)=-(s-1)a^(s-3)ζ(s,a)/2
            +s a^(s-2)ζ(s+1,a)/2 -1/(4a³).

At a=1 all coefficients are rational and ζ(u,1)=ζ(u).

## 3. A triangular span theorem gives an unconditional obstruction

Let ν_s be the fixed pushforward measure on [0,1], and M_d=integral x^d dν_s. Expanding (1-(1+t)^(-1))^d gives

    M_0=s/12,
    M_d=s/12+sum_(k=1)^d (-1)^k binom(d,k)J_k(1), d>=1.

Therefore the coefficient of ζ(s+d-1) in M_d is

    (-1)^d (s)_(d-1)/[2^(d-1)(d-1)!],

which is nonzero. Lower k cannot contribute that value. By triangular elimination,

    span_Q{1,M_1,...,M_d}
      =span_Q{1,ζ(s),ζ(s+1),...,ζ(s+d-1)}.             (1)

This is an equality of actual vector spaces and requires no unproved independence conjecture: both inclusions follow by solving a triangular system with nonzero rational diagonal.

**Obstruction theorem.** The first four nonconstant moments of this fixed compactified measure cannot all be rational affine expressions in any one real constant α.

Proof: by (1) their span together with 1 contains 1,ζ(s+1),ζ(s+3). Since s is odd, Euler's formula writes the latter two as nonzero rational multiples of π^(s+1),π^(s+3). These three numbers are linearly independent over Q: a dependence would give a nonzero rational polynomial vanishing at the transcendental number π. The span consequently has dimension at least 3, whereas Q+Qα has dimension at most 2. QED.

Thus the single-constant compact-capacity criterion does not apply to the fixed pushforward even at the 3-by-3 Gram stage, which uses moments through degree 4. This is a proof of an obstruction to that precise shortcut, not a general prohibition on compact measures in zeta arguments.

For s=5 the first examples are

    M_1=7/6-ζ(5),
    M_2=5/3-4ζ(5)+(5/2)ζ(6).

Already assuming ζ(5) rational would leave M_2 transcendental. At this stage one cannot legitimately clear its coefficients using the hypothetical denominator of ζ(5) alone. The full four-moment theorem is stronger: it rules out replacing ζ(5) by some other single constant to save the fixed construction.

## 4. A repair that preserves distinct simple poles

Return to a simple-pole Gram family

    G[i,j]=integral D_N(t)^r t^(i+j)/D_K(t) dμ_s(t),
    0<=i,j<h,

where D_m(t)=product_(j=1)^m(t+j²), K>N>=0 and r is a nonnegative integer. Its entries remain rational affine in ζ(s), and the matrix is positive definite.

Define an h-dependent positive measure ρ_h on [0,1] by pushing forward

    (1+t)^(2h-2) D_N(t)^r/D_K(t) dμ_s(t).

It has finite mass because the weight grows at most polynomially in t while μ_s has exponential decay in sqrt(t). For 0<=d<=2h-2 its x-moment is

    integral x^d dρ_h
      =integral t^d(1+t)^(2h-2-d) D_N(t)^r/D_K(t)dμ_s(t).

The numerator is an integer polynomial; the denominator still has only distinct simple poles. Thus precisely the moments needed for an h-by-h Gram determinant remain rational affine in the single value ζ(s). The statement stops at degree 2h-2; higher moments again introduce negative powers of 1+t.

Now choose the polynomial basis

    b_i(x)=x^i(1-x)^(h-1-i), 0<=i<h.

At x=t/(1+t), b_i(x)=t^i/(1+t)^(h-1). Therefore its Gram matrix under ρ_h is exactly G. Its coefficient matrix relative to 1,x,...,x^(h-1) is triangular with diagonal 1, so the determinant equals the ordinary monomial Gram determinant of ρ_h.

This is a concrete way to compactify while retaining the simple-pole arithmetic: vary the measure with the determinant size, and use the Bernstein-type polynomial basis above. It avoids falsely assigning a single-constant representation to all moments of one fixed measure.

## 5. Why the simplest mass-and-support bound still fails

Let M_h=ρ_h([0,1]). The Chebyshev argument now gives the valid inequality

    0<det G <= M_h^h 4^(h-1) 4^(-h(h-1)).             (2)

If d_h clears the rational affine coefficients of its first 2h-2 compact moments, a sufficient irrationality condition is

    limsup [log d_h+log M_h]/h < log 4.

Arithmetic bounds from `odd-arithmetic-barrier.md` supply explicit d_h candidates because the compact moment numerators are integer polynomials. But the mass M_h is not uniformly bounded, and in the natural families it is enormous.

Here is a rigorous diagnosis. Suppose E=rN+2h-2-K>=0. On y in [K,K+1] with K>=1 and t=y²,

    D_N(t)^r >= t^(rN),
    (1+t)^(2h-2) >= t^(2h-2),
    D_K(t) <= 2^K t^K.

The first positive exponential term in w_s gives

    w_s(y)>=c_s y^s exp(-2πy),
    c_s=2(2π)^(s-1)/(s-1)! >0.

Integrating only over that unit interval yields

    M_h >= c_s 2^(-K) K^(2E+s) exp(-2π(K+1)).         (3)

For h~βK and E~cK with β,c>0, (3) implies

    log M_h/h >= (2c/β)log K-O_s(1) → infinity.

Hence the sufficient condition based only on total mass and the full support [0,1] cannot hold for these raw compactified weights, since d_h>=1. This includes N=0,h=K (E=K-2) and the Fauzan proportions r=6,N~3K/40,h~37K/40 (E~1.3K).

This is **failure of the bound (2) as a useful estimate**, not a lower bound on det G. A large mass can coexist with a very small Gram determinant because the measure is concentrated and its polynomial directions nearly dependent. The orthogonal-polynomial norms, or an equivalent logarithmic-energy estimate, capture that dependence; a uniform norm bound on the whole interval discards it.

## 6. What to try next

The two precise barriers suggest two separate repairs:

* Keep the h-dependent measure so that the required moments retain distinct simple poles and a single target zeta value.
* Replace the crude full-interval/mass estimate by weighted monic competitors adapted to the concentration near x=1. Alternatively try an h-dependent map x=t/(a_h²+t), tracking both the basis determinant and every rational denominator factor introduced by a_h.

The compactification itself is exact and elementary. It exposes rather than removes the difficult approximation estimate. A successful simplification of Fauzan would need explicit competitors whose weighted norms beat the proved arithmetic cost, not just an interval of capacity 1/4.
