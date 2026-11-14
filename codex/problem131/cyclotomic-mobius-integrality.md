# Cyclotomic frequency depletion: an exact integrality classification and tangent moments

**Status:** complete proof below, with exact and high-precision checks in
`../computations/check_cyclotomic_mobius_integrality.py`; independent review
passed in cyclotomic-mobius-independent-review.md. The starting Möbius proposal came from the parent agent and was
checked rather than assumed. No irrationality or historical novelty claim.

Fix an odd s≥3, and use the already established positive density

\[
 w_s(y)=c_s y^s\sum_{n\ge1}n^{s-1}e^{-2\pi ny},\quad
 c_s=\frac{2(2\pi)^{s-1}}{(s-1)!},\quad
 \mu_s(R)=\int_0^\infty R(y^2)w_s(y)\,dy.
\]

Its polynomial moments and real pole values are

\[
 \mu_s(t^e)=(-1)^e\frac{B_k}{k}\binom{k+s-2}{s-1},\quad k=2e+2,
 \qquad
 \mu_s((t+a^2)^{-1})=a^{s-1}\zeta(s,a)-\frac1{2a}-\frac1{s-1}.
\]

## 1. Remove precisely the frequencies sharing a prime factor with m

For an integer m≥2 define

\[
 W_{s,m}(y)=\sum_{d\mid m}\mu(m/d)d\,w_s(y/d),\qquad
 \nu_{s,m}(R)=\int_0^\infty R(y^2)W_{s,m}(y)\,dy,
\]

where μ without a subscript is the Möbius function. Putting every exponential
on the common scale m gives

\[
 \boxed{W_{s,m}(y)=c_s m^{1-s}y^s
 \sum_{\substack{n\ge1\\(n,m)=1}}n^{s-1}e^{-2\pi ny/m}>0.} \tag{1}
\]

The coefficient of a frequency n is Σ_{q|(n,m)}μ(q), exactly the indicator
of coprimality. The density tends to φ(m)/π at zero and has exponential rate
2π/m at infinity. All nonnegative polynomial moments converge.

Let Jₖ(m)=Σ_{d|m}μ(m/d)dᵏ=mᵏ∏_{p|m}(1−p^(−k)) be the Jordan totient.
Dilation gives

\[
 \boxed{\nu_{s,m}(t^e)=(-1)^e\frac{J_k(m)B_k}{k}
                 \binom{k+s-2}{s-1},\quad k=2e+2.} \tag{2}
\]

## 2. Cyclotomic logarithmic derivatives classify integrality

Write Φₘ for the cyclotomic polynomial and define

\[
 \kappa_{m,k}=\left.\frac{d^k}{dz^k}\log\Phi_m(e^z)\right|_{z=0}.
\]

Möbius factorization Φₘ(q)=∏_{d|m}(qᵈ−1)^{μ(m/d)}, with cancellation of
the log z terms, gives

\[
 \kappa_{m,1}=\phi(m)/2,\quad
 \kappa_{m,k}=J_k(m)B_k/k\ (k\ge2\text{ even}),\quad
 \kappa_{m,k}=0\ (k>1\text{ odd}). \tag{3}
\]

**Theorem.** For each fixed odd s≥3, every polynomial moment of νₛ,ₘ is an
integer **if and only if m is composite**. Every such moment is strictly
positive. In fact, for composite m, every logarithmic derivative κₘ,ₖ is an
integer, including the first derivative.

Here is an elementary proof of the positive direction that does not use
Bernoulli denominator theorems.

If m is not a prime power, Φₘ(1)=1. All exponential derivatives
rₖ=(d/dz)ᵏΦₘ(eᶻ)|₀ are integers and r₀=1. Differentiating M′=M(log M)′ gives

\[
 \kappa_k=r_k-\sum_{j=1}^{k-1}\binom{k-1}{j-1}\kappa_jr_{k-j},
\]

so κₖ∈ℤ by induction.

If m=pᵃ with a≥2, use Φ_{pᵃ}(q)=Φ_p(q^{p^{a-1}}). The normalized exponential
series M_p(z)=Φ_p(eᶻ)/p has M_p(0)=1 and every derivative in p^(−1)ℤ.
The same recursion proves

\[
 p^k\kappa_{p,k}\in\mathbb Z\quad(k\ge1):
\]

each product in the recursion has denominator dividing p^(j+1), with j+1≤k.
Since κ_{pᵃ,k}=p^{(a−1)k}κ_{p,k}, all these derivatives are integral. This
also handles p=2 without a separate parity argument. Formula (2) is now an
integer binomial coefficient times (−1)ᵉκₘ,ₖ. Positivity follows from (1).

The fact Φₘ(1)=1 outside prime powers, and Φ_{pᵃ}(1)=p, follows directly by
letting q→1 in the Möbius product; the surviving product is
∏_{d|m}d^{μ(m/d)}. Thus no unmentioned arithmetic theorem enters this direction.

For prime m=p the density is exactly the previously studied single depletion
p w_s(y/p)−w_s(y). Its sharp uniform denominator is

\[
 \boxed{p^{1+v_p(s-1)}>1,} \tag{4}
\]

so the converse holds. For completeness, the sharpness argument is explicit.
Put h=v_p(s−1). For p=2 choose k=2ʰ. For odd p choose b>h with pᵇ>s−2 and
k=pʰ+(p−2)pᵇ. In both cases k is even, p−1 divides k, and adding k to s−2
in base p creates no carry. Von Staudt–Clausen gives v_p(Bₖ)=−1 and Lucas
makes binom(k+s−2,s−2) a p-adic unit. Rewriting (2) with denominator s−1
shows that the moment has valuation exactly −1−h. The upper bound and
absence of other primes are proved in
[integer-depletion-cumulants.md](integer-depletion-cumulants.md).

## 3. The smallest composite modulus gives tangent-number moments

For m=4,

\[
 \Phi_4(q)=q^2+1,\qquad
 \log\Phi_4(e^z)=\log2+z+\log\cosh z.
\]

Let T₂ᵣ₋₁ be the positive tangent numbers defined by

\[
 \tan z=\sum_{r\ge1}T_{2r-1}\frac{z^{2r-1}}{(2r-1)!}.
\]

Differentiating log cosh and using tanh z=−i tan(iz) gives
κ₄,₂ᵣ=(−1)^(r−1)T₂ᵣ₋₁. Consequently

\[
 \boxed{\nu_{s,4}(t^e)=T_{2e+1}\binom{2e+s}{s-1}\in\mathbb Z_{>0}.} \tag{5}
\]

The differential equation (tan z)′=1+tan²z proves the tangent numbers are
integers by an elementary coefficient recursion. Their first values are
1,2,16,272,… . For s=3 the moments begin 3,20,336,9792,… . This is a
particularly small positive continuous moment model with integer moments.

## 4. Prime-power integrality is a dilation, not a new normalized density

Let r=rad(m) and h=m/r. Inclusion-exclusion gives the exact identities

\[
 W_{s,m}(y)=hW_{s,r}(y/h),\qquad
 \nu_{s,m}(R)=h^2\nu_{s,r}(R(h^2t)). \tag{6}
\]

After restoring the original exponential decay,

\[
 \boxed{\sigma_{s,m}(y)=m^{-1}W_{s,m}(my)
 =\sum_{q\mid m}\frac{\mu(q)}q w_s(qy),} \tag{7}
\]

which depends only on the prime divisors of m. Its e-th moment is
m^(−k)νₛ,ₘ(tᵉ)=∏_{p|m}(1−p^(−k))μₛ(tᵉ).

In particular Wₛ,₄(y)=2Wₛ,₂(y/2), while
σₛ,₄(y)=w_s(y)−w_s(2y)/2 is **exactly the old normalized Euler density**.
Raw tangent integrality is useful and elementary, but after this normalization
the moments are the integers in (5) divided by 4^(2e+2). There is no new real
weight after normalization. The existing coprime 2-and-3 determinant
obstruction does not apply automatically to this one-prime family; the
separate Euler pole and collision calculations do apply after the indicated
change of scale.

## 5. Ramanujan sums give the exact pole dictionary

For every real a>0, scaling the original Hermite pole formula gives

\[
 \nu_{s,m}((t+a^2)^{-1})=
 a^{s-1}\sum_{d\mid m}\mu(m/d)d^{1-s}\zeta(s,a/d)
                         -\frac{\phi(m)}{2a}. \tag{8}
\]

The constant −1/(s−1) cancels because Σ_{d|m}μ(m/d)=0. For an integer
j≥1, define the integer Ramanujan sum

\[
 c_m(v)=\sum_{d\mid(m,v)}d\mu(m/d)
       =\sum_{\substack{a\bmod m\\(a,m)=1}}e^{2\pi iav/m}.
\]

Expanding each Hurwitz tail by residue classes gives

\[
 \boxed{\nu_{s,m}((t+j^2)^{-1})=
 j^{s-1}\sum_{\ell\ge j}\frac{c_m(\ell-j)}{\ell^s}
                         -\frac{\phi(m)}{2j}.} \tag{9}
\]

Thus cyclotomic logarithmic derivatives govern polynomial moments, while
Ramanujan sums govern poles. Both come from the same primitive roots of unity.
Formula (9) can equally be written as a finite sum of Liₛ at primitive mth
roots of unity, plus a rational tail. It need not reduce to a single zeta value.
No independence assumption about these constants is being made.

On the grid j=mr all Hurwitz arguments are integral. With
Dₘ,ₛ=∏_{p|m}(1−p^(s−1)), one obtains

\[
 \boxed{\nu_{s,m}((t+(mr)^2)^{-1})=
 r^{s-1}\left[D_{m,s}\zeta(s)
   -\sum_{q\mid m}\mu(q)q^{s-1}H_{qr}^{(s)}\right]
 +\frac{\phi(m)}{2mr}.} \tag{10}
\]

This retains a single ζ(s), but denominator and pole-spacing costs must still
be evaluated; integral polynomial moments do not clear these rational tails.

## 6. For modulus four, pole parity chooses between ζ(s) and π^s

Let η(s)=(1−2^(1−s))ζ(s), Aᵣ=Σ_{l≤r}(−1)^(l−1)/lˢ, χ₄ be the primitive
character modulo4, β(s)=Σ_{l≥1}χ₄(l)/lˢ, and Bⱼ=Σ_{l≤j}χ₄(l)/lˢ.
The dilation (6) gives, at even nodes j=2r,

\[
 \boxed{\nu_{s,4}((t+(2r)^2)^{-1})=
 (-1)^{r-1}r^{s-1}(\eta(s)-A_r)+\frac1{2r}.} \tag{11}
\]

This is exactly the old Euler pole value at node r. At odd nodes,
c₄(v) is 0 for odd v and 2(−1)^(v/2) for even v. Equation (9) therefore gives

\[
 \boxed{\nu_{s,4}((t+j^2)^{-1})=
 2\chi_4(j)j^{s-1}(\beta(s)-B_j)+\frac1j,
 \quad j\text{ odd}.} \tag{12}
\]

For odd s=2r+1, β(s)=(−1)ʳE₂ᵣπˢ/[4^(r+1)(2r)!], with E₂ᵣ the usual
Euler numbers. So the same positive tangent-moment measure encodes ζ(s) at
even nodes and a rational multiple of πˢ at odd nodes. This is a concrete
connection between two irrationality targets, not a proof that their mixture
has the required rational arithmetic.

Restricting to even nodes preserves the pure ζ(s) target but returns, by an
exact dilation, to the Euler model and its known pole-collision costs. Allowing
odd nodes introduces β(s). Either design still needs a separate determinant
normalization theorem and a real decay estimate after normalization.

## 7. Primary-source context

The Jordan/Bernoulli identity (3) is classical cyclotomic-derivative structure.
A directly applicable primary reference is Andrés Herrera-Poyatos and Pieter
Moree, [*Coefficients and higher order derivatives of cyclotomic polynomials:
old and new*, arXiv:1805.05207v2](https://arxiv.org/abs/1805.05207), Theorem1.3
(and the fuller Theorem5.8). Its left side is the Stirling-number expansion of
(q d/dq)ᵏ log Φₘ(q) at q=1, hence exactly our exponential logarithmic derivative.
Theorem3.3 gives the ordinary derivative version, related to Lehmer's 1966
work; §3.1 treats the Ramanujan-sum logarithmic expansion at zero.

That paper was retrieved and inspected after the proof above. This note
applies those classical quantities to the positive zeta-moment density and
its pole arithmetic. It does not claim historical originality for the
classification or this application. The elementary proof of composite
integrality and the exact scale/parity distinctions are self-contained.

Archived source: `../sources/herrera-moree-cyclotomic-180505207.pdf`, SHA256
`73ea070ea7cec72270cea512794e9e816c08c8213b9f69a3c30047895792704c`;
text extraction at the same stem with `.txt`, retrieved 2026-11-14 UTC.


## 8. The tangent-moment sequence is already primitive

For every fixed odd s≥3, a finite initial segment already has greatest common
 divisor one:

\[
 \boxed{\gcd_{0\le e\le(s-1)/2}\nu_{s,4}(t^e)=1.} \tag{13}
\]

Here is a short proof that does not require factoring large tangent numbers.
For every odd prime p and h≥1,

\[
 T_{p^h}\equiv(-1)^{h(p-1)/2}\pmod p. \tag{14}
\]

To see this, work over the quadratic extension of Fₚ containing i²=−1 and
let δ=(1+y²)d/dy. In the Möbius coordinate v=(y−i)/(y+i), one has δv=2iv,
so δ=2iv d/dv. The Euler derivation E=v d/dv satisfies Eᵖ=E on rational
functions: its pth power is a derivation in characteristic p and agrees with
E on the generator v. Consequently

\[
 \delta^p=(2i)^{p-1}\delta=(-1)^{(p-1)/2}\delta,
 \qquad
 \delta^{p^h}=(-1)^{h(p-1)/2}\delta.
\]

The tangent differential equation gives Tᵣ=(δʳy)|_{y=0} for odd r, proving
(14). All evaluations here are of polynomials in y; the extension field and
Möbius coordinate only establish the derivation identity.

Now the zeroth moment is s, so a prime dividing every moment would be an odd
prime p|s. Write s=pʰu with p∤u and choose e=(pʰ−1)/2. Formula (5) gives

\[
 \nu_{s,4}(t^e)=T_{p^h}\binom{s+p^h-1}{p^h}.
\]

Lucas's congruence gives the binomial factor congruent to u modulo p: the
base-p digit at position h in s+pʰ−1 is u modulo p, and the lower h digits
are all p−1. Equation (14) makes the tangent factor a unit modulo p as well.
Thus p does not divide this moment. Since e≤(s−1)/2, every prime divisor of
the zeroth moment is excluded within the stated finite segment, proving (13).

This rules out reducing the integral moment sequence by a fixed integer
factor greater than one. It is not a lower bound for determinant denominators
or an obstruction to index-dependent cancellations in a proposed proof.
