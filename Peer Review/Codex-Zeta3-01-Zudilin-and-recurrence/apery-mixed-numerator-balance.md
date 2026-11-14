# Two short rational telescopes explain Apéry’s recurrence and mixed numerator

Status: complete elementary proof, independently reviewed in
[apery-mixed-numerator-independent-review.md](apery-mixed-numerator-independent-review.md). This
closes the mixed-numerator conjecture recorded in
[rank-one-primary-source-context.md](../Codex-Zeta3-05-Schur-Apery/rank-one-primary-source-context.md), §7.
It does **not** by itself identify a truncated Hankel Schur root with the
resulting difference quotient. No historical novelty claim is made.

Let

\[
 \phi_k(x)=\binom xk\binom{x+k}k
 =\frac1{(k!)^2}\prod_{j=0}^{k-1}[x(x+1)-j(j+1)],
 \qquad \Pi_n(x)=\sum_{k=0}^n\phi_k(n)\phi_k(x).
\]

The subscript n here is Prévost's subscript2n. Set A_n=Π_n(n).
Let ℬ be the Bernoulli functional ℬ(xʲ)=Bⱼ, with B₁=−1/2, and set
ℒ(f)=ℬ(f′). Define

\[
 \Theta_n(y)=\mathcal L_x\left(
             \frac{\Pi_n(x)-\Pi_n(y)}{x-y}\right),\qquad
 C_n(m)=H_m^{(3)}\Pi_n(m)+\frac12\Theta_n(m),\quad m\ge0.
\]

Prévost's formula identifies C_n(n)=B_n, the usual Apéry rational numerator
(B₀=0,B₁=6). The proof below independently identifies the entire mixed family
with rational parts of rational-derivative sums.

**Theorem.** For every n≥1,

\[
 \boxed{C_n(n-1)=\frac{B_n+B_{n-1}}6.} \tag{1}
\]

Consequently the finite-difference quotient is exactly

\[
 \boxed{\frac{C_n(n)-C_n(n-1)}{\Pi_n(n)-\Pi_n(n-1)}
       =\frac{5B_n-B_{n-1}}{5A_n-A_{n-1}}.} \tag{2}
\]

The denominator does not vanish: Π_n(n)>Π_n(n−1), since each Newton-basis
summand is nondecreasing there and the k=1 term strictly increases.

**A shorter route to the recurrence itself.** Readers interested primarily in
why Apéry's cubic recurrence has coefficients17 and5 can go directly to
§§3–5. They need only the mixed rational products PₙPₘ and two low-degree
polynomial identities. The Bernoulli-functional work in §§1–2 is needed to
identify Prévost's associated numerator, not to derive the recurrence for the
original irrationality linear forms. No complex integral or large diagonal
creative-telescoping certificate is required for that derivation.

## 1. Two elementary Bernoulli identities replace an orthogonality import

The generating function ℬ(eᶻˣ)=z/(eᶻ−1) proves, for every polynomial f,

\[
 \mathcal B(f(x+1)-f(x))=f'(0),\qquad
 \mathcal B(f(-x-1))=\mathcal B(f(x)).
\]

Differentiating inside the functional gives

\[
 \mathcal L(f(x+1)-f(x))=f''(0),\qquad
 \mathcal L(f(-x-1))=-\mathcal L(f(x)). \tag{3}
\]

Since φₖ(−x−1)=φₖ(x), it follows that ℒ(φₖ)=0. A less obvious but equally
useful identity is

\[
 \boxed{\mathcal L(x\phi_k(x))=\frac{(-1)^k}{k+1}.} \tag{4}
\]

To prove it, introduce the even polynomial

\[
 f_k(x)=\frac{x^2\prod_{j=1}^k(x^2-j^2)}{((k+1)!)^2}.
\]

Its difference is exactly

\[
 f_k(x+1)-f_k(x)=\frac{2x+1}{k+1}\phi_k(x).
\]

Both sides have the same roots and leading coefficient, or this can be
checked directly from their factored products. Now
fₖ″(0)=2(−1)ᵏ/(k+1)², so (3) proves(4).

Therefore ℒ(Πₙ)=0 for all n, and for n≥1,

\[
 \mathcal L(x\Pi_n)=\sum_{k=0}^n
             \phi_k(n)\frac{(-1)^k}{k+1}=0. \tag{5}
\]

For the final equality, the sum is the integral on[0,1] of
Lₙ(u)=Σₖ(−1)ᵏφₖ(n)uᵏ=(1/n!)Dⁿ[uⁿ(1−u)ⁿ]. That integral vanishes by
one integration of the nth derivative and the endpoint zeros. At n=0,
ℒ(xΠ₀)=1 instead. No Wilson orthogonality theorem is needed for(5).

## 2. The associated numerator obeys the same difference equation

The elementary Newton-basis calculation in
[newton-bispectral-motivation.md](https://github.com/BoomydayCoder/Mathathon/blob/03cc75d694220d0fba888c637e573a29799e6d5b/codex/problem131/newton-bispectral-motivation.md) gives

\[
 (n+1)^3\Pi_{n+1}(y)=
 (2n+1)[n(n+1)+2y(y+1)+1]\Pi_n(y)-n^3\Pi_{n-1}(y). \tag{6}
\]

It is a polynomial identity in y. Take the divided difference in x,y and
apply ℒₓ. For multiplication by x(x+1), the extra term is
ℒ[(x+y+1)Πₙ(x)]. By(5) this vanishes for n≥1. Consequently Θₙ satisfies
exactly(6) for n≥1, with initial values Θ₀=0, Θ₁=2. Indeed Π₀=1 and
Π₁(x)=2x(x+1)+1, whose divided difference is2(x+y+1).

For each fixed nonnegative integer m, Cₙ(m) also satisfies(6) in n≥1,
with

\[
 C_0(m)=H_m^{(3)},\qquad
 C_1(m)=[2m(m+1)+1]H_m^{(3)}+1. \tag{7}
\]

## 3. A rational telescoping certificate identifies the mixed linear forms

Set

\[
 P_n(t)=\frac{\prod_{j=1}^n(t-j)}{\prod_{j=0}^n(t+j)},\qquad
 S_{n,m}(t)=(-1)^{n+m}P_n(t)P_m(t),\qquad
 F_{n,m}=-\frac12\sum_{t=1}^{\infty}S_{n,m}'(t).
\]

These sums converge absolutely, since S=O(t⁻²), S′=O(t⁻³). They are
symmetric in n,m. Define

\[
 Q_{n,m}(t)=(t-m-1)(t+m)S_{n,m}(t).
\]

Then the following rational identity holds for n≥0 (omit the n³ term when
n=0):

\[
\begin{aligned}
 &(n+1)^3S_{n+1,m}(t)
 -(2n+1)[n(n+1)+2m(m+1)+1]S_{n,m}(t)
 +n^3S_{n-1,m}(t)\\
 &\hspace{2em}=-2(2n+1)[Q_{n,m}(t+1)-Q_{n,m}(t)]. \tag{8}
\end{aligned}
\]

Here is a short verification usable without creative-telescoping software.
The ratios in n are −(t−n−1)/(t+n+1) and −(t+n)/(t−n), while

\[
 \frac{S_{n,m}(t+1)}{S_{n,m}(t)}=
 \frac{t^4}{(t-n)(t+n+1)(t-m)(t+m+1)}.
\]

Put N=n(n+1), M=m(m+1), and multiply (8) by
(t−n)(t+n+1)/Sₙ,ₘ(t). Its right side is −2(2n+1) times

\[
 t^4-(t^2+t-N)(t^2-t-M)
 =(N+M+1)t^2+(M-N)t-NM.
\]

Expanding the three quadratic terms on the left gives exactly the same
polynomial. Cancellation proves the rational identity even at apparent zeros
of the displayed ratios.

Differentiate(8) and sum at positive integers. The derivative of Q tends to
zero at infinity. For n≥1, Q′(1)=0: if m≥1, PₙPₘ has a double zero at1;
if m=0, its simple zero is supplemented by the factor t−1 in Q.
Thus Fₙ,ₘ satisfies the homogeneous recurrence(6) for n≥1.

At n=0, Q′₀,ₘ(1)=1 for every m≥0. For m=0 this follows from
Q₀,₀=(t−1)/t. For m≥1,
Pₘ′(1)=(−1)^(m−1)/[m(m+1)] and the other factors give the same value1.
Accordingly the boundary equation is

\[
 F_{1,m}-[2m(m+1)+1]F_{0,m}=-1. \tag{9}
\]

Since F₀,₀=ζ(3), symmetry and the recurrence at spectral parameter0 imply

\[
 (m+1)^3(F_{0,m+1}-F_{0,m})
 =m^3(F_{0,m}-F_{0,m-1})=-1,
\]

starting with F₀,₁−F₀,₀=−1. Hence

\[
 F_{0,m}=\zeta(3)-H_m^{(3)}.
\]

Comparing these initial values, (9), and(7), uniqueness of the second-order
recurrence proves the exact all-n,m identification

\[
 \boxed{F_{n,m}=\Pi_n(m)\zeta(3)-C_n(m).} \tag{10}
\]

In particular Cₙ(m)=Cₘ(n); this symmetry follows because both F and Π are
symmetric. The proof has never compared coefficients of an assumed-irrational
constant. It uses a rational telescoping identity and equal initial values.

## 4. A second, especially short telescope proves the adjacent balance

For n≥1, put

\[
 T_n(t)=\frac4{n^2}(t+n)^2P_n(t)^2.
\]

The decisive identity is

\[
 \boxed{P_n(t)^2+P_{n-1}(t)^2+6P_n(t)P_{n-1}(t)
                     =T_n(t+1)-T_n(t).} \tag{11}
\]

Indeed set r(t)=Pₙ(t)²/(t−n)². The left side is4(2t²−n²)r(t), while
Tₙ(t)=4(t²−n²)²r(t)/n² and
r(t+1)/r(t)=t⁴/[(t−n+1)²(t+n+1)²]. Thus the difference is
4[t⁴−(t²−n²)²]r(t)/n², exactly the left side.

Both Tₙ′(1)=0 and lim_(t→∞)Tₙ′(t)=0. Differentiating and summing(11)
therefore gives

\[
 \boxed{F_{n,n}+F_{n-1,n-1}=6F_{n,n-1}.} \tag{12}
\]

The same rational identity proves the denominator balance without any
irrationality assumption. The residues of Pₙ at−k are
(−1)^(n+k)φₖ(n). Hence the sum of double-pole coefficients on the left
of(11) is Aₙ+Aₙ₋₁−6Πₙ(n−1). Translation preserves the sum of double-pole
coefficients of a rational function, so the right side has sum zero. Thus

\[
 A_n+A_{n-1}=6\Pi_n(n-1). \tag{13}
\]

Finally insert(10) into(12) and use(13). The ζ(3) term cancels by a
separately proved rational identity, leaving exactly(1). Equation(2) follows.

## 5. A direct elementary derivation of the Apéry recurrence

Let fₙ=Fₙ,ₙ=−(1/2)Σₜ₌₁∞(Pₙ(t)²)′. Apply the mixed recurrence proved by(8)
at m=n, and then substitute the adjacent balance(12):

\[
 (n+1)^3\frac{f_{n+1}+f_n}{6}
 -(2n+1)[3n(n+1)+1]f_n
 +n^3\frac{f_n+f_{n-1}}6=0.
\]

Since (n+1)³+n³=(2n+1)[n(n+1)+1], this is exactly

\[
 \boxed{(n+1)^3f_{n+1}
 -(2n+1)(17n^2+17n+5)f_n+n^3f_{n-1}=0,
 \qquad n\ge1.} \tag{14}
\]

Thus17=6·3−1 and5=6·1−1: the number6 comes from the adjacent rational
finite difference(11), and the3 comes from the two-variable eigenvalue
coefficient at m=n. These constants arise from the interaction of two simple
identities, rather than being guessed in a large diagonal recurrence.

The boundary term must not be dropped. At n=0 the mixed equation has the
nonzero right side−1 in(9). It gives F₀,₁=ζ(3)−1 and
F₁,₁=5F₀,₁−1=5ζ(3)−6. Consequently f₀=ζ(3), f₁=5ζ(3)−6, and(14) is
asserted for n≥1, exactly as required. The positive Beukers integral is2fₙ;
its normalization therefore satisfies the same recurrence.

The denominator sequence Aₙ also satisfies(14), by applying the elementary
polynomial recurrence(6) and its adjacent closure(13). Subtracting the two
recurrences from fₙ=Aₙζ(3)−Bₙ then proves the recurrence for Bₙ, without
assuming ζ(3) irrational or separating coefficients of an unknown constant.
This recovers B₀=0,B₁=6 as well.

## 6. Interpretation and remaining scope

The same number6 appears in denominator and numerator because both arise
from the same rational finite difference(11). Summing its double-pole
coefficients gives the denominator identity; summing its derivative at the
positive integers gives the linear-form identity. The Bernoulli functional
is the bridge to Prévost's associated numerator, and its endpoint correction
is the explicit−1 in(9).

This proves the source-guided mixed numerator conjecture. It remains a
separate task to show that the rank-one Hankel root q_(n,3) equals the
finite-difference quotient(2). An exact proof of that further identification
must account for the varying measure and cannot substitute a formal
negative-parameter Wilson polynomial for its true orthogonal polynomial.

Run `python codex/computations/check_apery_mixed_numerator_balance.py`.
The symbolic rational certificates hold with indeterminate n,m,t; finite
checks separately test the Bernoulli identities, rational parts of the mixed
sums, associated-polynomial identification, and the claimed diagonal balance.
