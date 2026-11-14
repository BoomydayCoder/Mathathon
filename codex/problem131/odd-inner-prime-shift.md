# The fixed-odd-index cost in the inner prime allocation

Status: exact algebraic shift and a proof draft of the corresponding local
extension for fixed odd s≥5. This builds on the source's local basis method and
does not provide a completed general-s outer prime analysis or an irrationality
proof for ζ(7). The distinction between the envelope and the actual determinant
is essential throughout.

## 1. Why the index enters a local prime calculation

Let L(x^k)=B_k and define, on polynomials,

\[
 \tau_s(P)=\frac{L(P^{(s-2)})}{(s-1)!}.
\]

On simple integer poles set τ_{s,X}(1/(x-r))=H_{d(r)}^{(s)}−X,
where d(r)=r for r≥0 and d(r)=−r−1 otherwise. Polynomial division and partial
fractions extend this functional. The generalized moment functional satisfies

\[
 \mu_{s,X}(R)=\tau_{s,X}(x^sR(-x^2)).
\]

For fixed odd s and a prime p>s, the distribution formula has factor
p^{−(s−1)}:

\[
 \tau_{s,X}(g)=p^{-(s-1)}\sum_{a=0}^{p-1}
      \tau^{\rm ext}_{s,Y}(g(a+px)),\qquad Y=p^sX+C_{p,s},
 \quad C_{p,s}\in p^s\mathbb Z_p.
\]

Here the polynomial identity follows from Bernoulli multiplication after s−2
differentiations. To check the constant, the far-pole series is

\[
 \tau_s^{\rm an}\left(\frac1{x-z}\right)
 =-\frac1{s-1}\sum_{k\ge0}\binom{k+s-2}{s-2}B_kz^{-k-s+1}
 \quad(v_p(z)<0).
\]

In the sum over z=−a/p, the k=0 term is divisible by p^s because
∑_{a=1}^{p−1}a^{−(s−1)}≡0 mod p when p>s. The k=1 term already has p^s;
k≥2 terms have at least p^{k+s−2} by von Staudt–Clausen. The identity at 1/x
determines Y. Reflection and the translation identity
τ_s(g(x+1)−g(x))=g^{(s−1)}(0)/(s−1)! extend it to every integer pole,
just as in the source's s=5 argument.

The analytic completion still loses at most one p-power. Indeed
τ_s(x^d)=d(d−1)⋯(d−s+3)B_{d−s+2}/(s−1)!, and p>s makes the factorial
denominator a unit. For s≥5 and d≤p+1, its Bernoulli index is at most p−2,
so it is p-integral. Thus the source's leading-degree restriction d≤p+1
continues to suffice. Higher expansion terms are absorbed by their p-powers.

## 2. Ordinary and zero classes behave differently

Keep the source's determinant geometry r=6, α=3/40, λ=37/40 and its inner
range K/M<p≤K/3, with M fixed. At an ordinary class a≠0, the factor x^s is a
unit. Extracting near factors and applying the distribution formula gives

\[
 \nu_i(a)+\nu_j(a)+6\ell_N(a)-\ell_K(a)-(s-1). \tag{1}
\]

At the zero class x=pz, x^s contributes p^s, which cancels the index-dependent
distribution cost, leaving the same bound as for s=5:

\[
 2\nu_i(0)+2\nu_j(0)+12m_N-2m_K+1. \tag{2}
\]

The zero-class leading numerator degree increases from 5+2ν_i+2ν_j+12m_N to
s+2ν_i+2ν_j+12m_N. Since s and M are fixed while p→∞ with K, this remains
below p+1 eventually. Ordinary leading degrees are unchanged: the constant
term of (a+pz)^s is a unit and its other terms carry p-powers.

Use the same local dimensions L_a and CRT basis as the source, and define

\[
 w^{(s)}_{a,i}=i+3\ell_N(a)-\frac{\ell_K(a)+s-1}{2},\quad a>0,
\]
\[
 w^{(s)}_{0,i}=\min\left(2i+6m_N-m_K+\frac12,
       \min_{c>0}\left[Z_c-\frac{\ell_K(c)+s-1}{2}\right]\right).
\]

The ordinary-source comparisons are identical to the source's after subtracting
(s−5)/2 from both ordinary half-weights. The zero-source comparisons remain
valid because all assigned weights weakly decrease for s≥5. Hence the same
determinant argument gives the lower envelope γ_{p,s}^{in}=2∑w^{(s)}.
If L₀ is the reserved zero-block dimension, the exact comparison is

\[
 -(s-5)h\le\gamma_{p,s}^{in}-\gamma_{p,5}^{in}
 \le-(s-5)(h-L_0). \tag{3}
\]

This local proof draft uses the same CRT, far-pole completion and determinant
valuation arguments as the source. It does not replace their independent
audit; it identifies explicitly every place where the fixed index changes.

## 3. Exact limiting formula and why reallocation alone cannot help

Let T=⌊2Hx⌋, ξ=Hx−T/2, q=⌊2x⌋ and n₊=(2x−q)/2, with H=1+2α.
Replacing the source's constant 5 by s gives

\[
 \Gamma_s(x)=\int_0^{1/2}(T-b)(T+b-\ell-s)\,dz
       +\xi(2T-q-s)+(\xi-n_+)_+,
 \quad b=3\ell(\alpha x,z).
\]

Since ∫b=3αx and T/2+ξ=Hx,

\[
 \boxed{\Gamma_s(x)=\Gamma_5(x)-(s-5)\lambda x.} \tag{4}
\]

The factorial normalization S_K does not depend on s, so
R_s(x)=R_5(x)+(s−5)λx. This also follows from (3), since L₀=O_M(1) and
h/p→λx. Thus the inner prime sum changes by exactly

\[
 (s-5)\lambda\int_3^M x^{-2}\,dx
 =(s-5)\lambda\left(\frac13-\frac1M\right). \tag{5}
\]

More generally, in any allocation governed by the same ordinary row inequality
(1), each ordinary row has the same half-weight loss (s−5)/2. Reallocating a
fixed total of h−O_M(1) ordinary rows among classes therefore leaves this leading
loss unchanged. This is an obstruction to **that particular improvement**.
It does not rule out sharper local functionals, extra residue cancellation,
another determinant, or allocations outside the stated assumptions.

## 4. The small-prime term confirms a cutoff-independent cost

For fixed s, the source's small-prime proof adapts to

\[
 v_p^G(F_{K,s})\ge -(s+1)h\lfloor\log_p(5K)\rfloor
                       -h\,v_p((s-1)!)
\]

for sufficiently large K. The reason is that the polynomial functional loses
at most (s−1) logarithmic factors, while the quotient bound loses at most two;
a residue loses one logarithmic factor and H_j^{(s)} loses at most s. The
fixed factorial only affects finitely many small primes and contributes o(K²).
The pullback degree grows by the fixed amount s−5, so the bound by 5K persists
eventually. This is a sufficient bound, not a minimal denominator formula.

Consequently the small-prime contribution becomes (s+1)λ/M instead of 6λ/M.
Together with (5), the net extra small-plus-inner cost is

\[
 \boxed{(s-5)\lambda/3.}
\]

At s=7 it equals 37/60≈0.616667. The independently certified analytic exponent
is about −1.38551 for every fixed odd s in this common geometry, whereas the
successful s=5 arithmetic budget has only a small spare margin. Transferring
the same small/inner envelope therefore spends far more than that margin.

This calculation alone does **not** prove failure of ζ(7) irrationality or of
the actual determinant sequence. The outer-prime arithmetic also changes,
and the displayed envelopes need not be sharp. In particular higher s cancels
more Bernoulli denominators; see
[`odd-arithmetic-barrier.md`](odd-arithmetic-barrier.md). A complete outer
analysis or a different source of local cancellation is still needed before
making a conclusion about the total normalized exponent.

## 5. A stronger obstruction for the direct outer-weight extension

The source's outer residue argument also extends transparently for fixed s.
After removing the nonintegral polynomial moments, its ordinary row weights are

\[
 w^{(s)}_{a,i}=\begin{cases}
 \min(0,i+3\delta-(\ell+s-1)/2),&i<\ell-2,\\
 0,&i\ge\ell-2.
 \end{cases}
\]

Here ℓ is the class size and δ indicates a removed pole. The only change in
the bad residue estimate is that a harmonic value may lose s powers of p
instead of 5. The two good poles a and p−a remain congruent modulo p for odd s:
H_{p−a}^{(s)}≡H_{a−1}^{(s)}, their X coefficients agree since s−1 is even,
and the 1/(2j) terms cancel the remaining difference. These statements require
p>s+1, automatic eventually for fixed s. The one- and two-pole zero-class
costs remain 1 and 4: j^{s−1}H_j^{(s)} has valuation at least −1 for j=p,2p,
and the node difference contributes the same two powers as before.

All ordinary residue costs are therefore nondecreasing with s. For example,
for removed classes with ℓ=2,3,4,5,6, the s=5 costs are 0,1,2,4,6, whereas
the s=7 costs are 0,3,6,9,12. For an unremoved class the cost is
(s+2)(ℓ−2). These are costs in a sufficient valuation estimate, not exact
valuations of the determinant.

Higher s improves the polynomial correction. Write m=(s−1)/2. The Bernoulli
denominator cancellation makes all polynomial moments integral below
e₀=m(p−1)−1. Thus the correction rank is bounded by

\[
 r_{p,s}=(K+4N-mp+m)_+.
\]

To see this, the quotient degree is i+j+6N−K; the first potentially nonzero
correction lies at i+j≥K−6N+e₀. The initial zero rows and columns give rank
at most 2h−1−(K−6N+e₀)=K+4N−mp+m. For s=5 this recovers the source's bound.

Even granting **complete removal** of the s=5 plain-rank correction, its
maximum possible outer saving is only

\[
 \int_{1/3}^{\infty}(1+4\alpha-2y)_+\,dy
 =\frac{(1+4\alpha-2/3)^2}{4}=\frac{361}{3600}.
\]

This observation gives a precise limitation of the direct extension. Let
A_s(M) denote the prime-sum arithmetic bound obtained with the ordinary
weights above, the inner allocation already analyzed, and any outer rank
improvement whose penalty is nonnegative. Then

\[
 A_s(M)\ge I_{\rm out,5}^{\rm plain}
  +\frac{6\lambda}{M}+\int_3^M R_5(x)x^{-3}\,dx
  +\frac{(s-5)\lambda}{3}-\frac{361}{3600}.
\]

Here I_{out,5}^{plain}=129101/96000. The displayed inequality is a **lower
bound on an upper-bound expression**, not a lower bound on actual normalizer
growth. Its purpose is to test whether this proposed estimate could suffice.

There is a short entirely rational lower bound. Since Q≥−1/2,
F≤6λ and Q≤13/8,

\[
 \frac{6\lambda}{M}+\int_3^M R_5x^{-3}
 \ge J-\frac1{36}-\frac{13}{16M^2}.
\]

The elementary upper bound γ<H₃−log3 gives
J>λ(4/3+10α+12α logα). Moreover atanh z<z/(1−z²) gives
logα>−2(4/15)−6(3/8)=−167/60. Thus J>−9361/24000.
For all fixed odd s≥7 and M≥40 we conclude

\[
 \boxed{A_s(M)>\frac{1662179}{1152000}
                 =1.442863715\ldots>\frac75.}
\]

Consequently this direct family of arithmetic envelopes cannot meet the shared
certified analytic budget 1.38, even if every outer rank penalty is waived.
Rearranging the same inner row allocation cannot change this conclusion.
This is a concrete obstruction to a proof strategy, not to ζ(7) irrationality
or to sharper arithmetic for these determinants. It suggests looking for new
local residue cancellation, or changing the determinant parameters, before
investing in more precise integration of this same envelope.

The rational comparisons and outer class cost tables are checked by
`python codex/computations/check_odd_inner_shift.py`.
