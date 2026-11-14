# Euler depletion: a positive η functional with dyadic polynomial moments

Status: the positive model, pole identities, sharp uniform polynomial
denominator, and scaling identities below are proved. The determinant
comparisons are finite exact computations with numerical evaluations, not
an asymptotic theorem. This is an independent development of the parent's
proposed Euler depletion, using the previously proved Hermite moment model.
No historical novelty or new irrationality theorem is asserted.

## 1. Removing the even frequencies preserves positivity

For odd \(s\ge3\), recall

\[
w_s(y)=c_s y^s\sum_{n\ge1}n^{s-1}e^{-2\pi ny},
\qquad c_s=\frac{2(2\pi)^{s-1}}{(s-1)!}.
\]

Define

\[
\boxed{\widetilde w_s(y)=2w_s(y/2)-w_s(y)
=c_s2^{1-s}y^s\sum_{\substack{n\ge1\\n\ {\rm odd}}}
n^{s-1}e^{-\pi ny}>0.}                                  \tag{1}
\]

Indeed, the even-frequency terms in the first expression cancel the whole
second expression, term by term. At zero this density tends to \(1/\pi\);
at infinity it is

\[
c_s2^{1-s}y^s e^{-\pi y}(1+O_s(e^{-2\pi y})).
\]

Set \(\nu_s(R)=\int_0^\infty R(y^2)\widetilde w_s(y)\,dy\).
For polynomial and simple negative-square poles, substitution \(y=2u\)
gives the exact operator identity

\[
\nu_s(R)=4\mu_s(R(4t))-\mu_s(R(t)).
\]

Thus

\[
\boxed{\nu_s(t^e)=(2^{2e+2}-1)\mu_s(t^e).}                \tag{2}
\]

The price visible immediately in (1) is the change of real exponential
decay from \(2\pi\) to \(\pi\).

## 2. The polynomial denominator is a fixed power of two

Put \(k=2e+2\), and let

\[
G_k=2(1-2^k)B_k
\]

be the even-index Genocchi number. Formula (2) becomes

\[
\nu_s(t^e)=
\frac{(-1)^{e+1}G_k}{2(s-1)}
 \binom{k+s-2}{s-2}.                                    \tag{3}
\]

Von Staudt–Clausen implies that \(G_k\) is an odd integer: every odd
prime in the Bernoulli denominator divides \(2^k-1\), while the factor two
cancels the unique factor two in that denominator. Consequently the
denominator in (3) divides \(2(s-1)\), already a substantial sharpening
of the immediate bound \(2(s-1)!\).

In fact every odd prime dividing \(s-1\) cancels too. Define the Euler
polynomial values at zero by

\[
\frac2{e^z+1}=\sum_{n\ge0}E_n(0)\frac{z^n}{n!}.
\]

Their coefficient recurrence is

\[
E_0(0)=1,\qquad
E_n(0)=-\frac12\sum_{j=0}^{n-1}\binom njE_j(0)\quad(n\ge1),
\]

so all \(E_n(0)\) are dyadic rationals. The generating function gives
\(G_k=kE_{k-1}(0)\), and the other binomial form of (3) is

\[
\nu_s(t^e)=
\frac{(-1)^{e+1}}2 E_{k-1}(0)
 \binom{k+s-2}{s-1}\in\mathbb Z[1/2].
\]

Combining the two denominator statements proves the sharp uniform result

\[
\boxed{2^{\,1+v_2(s-1)}\nu_s(t^e)\in\mathbb Z
\quad\text{for every }e\ge0.}                            \tag{4}
\]

For sharpness, let \(a=v_2(s-1)\ge1\), choose \(k=2^a\), and
\(e=2^{a-1}-1\). The binomial coefficient in (3) is odd: addition of
\(k\) and \(s-2\) in base two has no carries, since the \(a\)-th bit
of \(s-2\) is zero and \(k\) has no other nonzero bit. Equivalently,
Lucas's parity criterion gives
\(\binom{k+s-2}{s-2}\equiv1\pmod2\).
Because \(G_k\) is odd, equality is attained in the denominator bound.
For \(s=5\) the optimal uniform denominator is \(8\); for \(s=7\) it is \(4\).

This removes the odd-prime *polynomial* correction. It says nothing yet
about denominators of pole moments or divided differences.

### Bernoulli averaging becomes one extra Euler derivative

Let \(\mathcal E(v^n)=E_n(0)\), and for polynomials define

\[
\tau_E(P)=-\frac{\mathcal E(P^{(s-1)})}{2(s-1)!}.
\]

The preceding binomial identity proves exactly

\[
\nu_s(R(t))=\tau_E\!\left(v^sR(-v^2)\right)
\quad\text{for polynomial }R.
\]

For every odd positive integer \(p\), the finite geometric sum gives

\[
\frac2{e^z+1}
=\sum_{a=0}^{p-1}(-1)^ae^{az}\frac2{e^{pz}+1}.
\]

Equating coefficients proves the Euler distribution relation

\[
\mathcal E(P(v))=\sum_{a=0}^{p-1}(-1)^a
                    \mathcal E(P(a+pv)).
\]

After the \(s-1\) derivatives, this becomes

\[
\boxed{\tau_E(P)=p^{1-s}\sum_{a=0}^{p-1}(-1)^a
                                    \tau_E(P(a+pv)).}
\]

The scale \(p^{1-s}\) is the same as for the original Bernoulli
derivative functional. Bernoulli distribution has an averaging factor
\(1/p\) and \(s-2\) derivatives; Euler distribution has no averaging
factor but needs \(s-1\) derivatives. Thus the disappearance of the
polynomial Bernoulli denominator does not automatically improve the local
distribution exponent. These identities are proved here on polynomials;
an extension to rational poles must separately respect their defining
values and is not being assumed.

## 3. The poles isolate the alternating zeta value

Write

\[
\eta(s)=(1-2^{1-s})\zeta(s),\qquad
A_j^{(s)}=\sum_{k=1}^{j}\frac{(-1)^{k-1}}{k^s}.
\]

The pole formula is

\[
\boxed{\nu_s((t+j^2)^{-1})
=(-1)^{j-1}j^{s-1}\bigl(\eta(s)-A_j^{(s)}\bigr)
 +\frac1{2j}.}                                          \tag{5}
\]

To verify all constants, the original Hermite formula for every real \(a>0\)
is

\[
\mu_s((t+a^2)^{-1})
=a^{s-1}\zeta(s,a)-\frac1{2a}-\frac1{s-1}.
\]

Scaling gives the difference of this expression at \(a=j/2\) and \(a=j\):

\[
j^{s-1}\bigl(2^{1-s}\zeta(s,j/2)-\zeta(s,j)\bigr)
-\frac1{2j}.
\]

Separating even and odd terms in these Hurwitz tails gives

\[
2^{1-s}\zeta(s,j/2)-\zeta(s,j)
=j^{-s}+(-1)^{j-1}(\eta(s)-A_j^{(s)}),
\]

which proves (5). In particular, the original term \(-1/(s-1)\) cancels;
it must not be retained in the new pole formula.

Define the affine rational functional \(\nu_{s,X}\) by (2) on polynomials
and (5) with \(\eta(s)\) replaced by \(X\) on poles. Its positive integral
representation holds at \(X=\eta(s)\). If an integer polynomial criterion
could be met, irrationality of this target would be equivalent to
irrationality of \(\zeta(s)\).

## 4. Two precise failures of the old pole argument

### All integer nodes: the coefficient signs break a congruence

For an odd prime \(p\), the colliding nodes \(j\) and \(p-j\), both below
\(p\), have opposite parity. Their \(X\)-coefficients in (5) are opposite
modulo \(p\):

\[
(-1)^{p-j-1}(p-j)^{s-1}
\equiv-(-1)^{j-1}j^{s-1}\pmod p.
\]

Their difference is a nonzero unit modulo \(p\). The node difference
\((p-j)^2-j^2=p(p-2j)\) has valuation one. Thus the divided difference's
\(X\)-coefficient has valuation \(-1\). Both pole values individually
are integral, but they are not congruent coefficientwise.

This is exactly where the original two-good-pole argument fails. It also
permits denominator primes greater than the largest integer node index
\(K\), when two available nodes add to such a prime.

### Odd nodes: the second colliding pole is no longer good

Restricting to odd nodes makes all \(X\)-coefficients positive, but \(j\)
and \(p-j\) cannot both occur. The next reflected odd node is \(2p-j\),
whose alternating harmonic sum contains a term with denominator \(p^s\).

For a completely explicit example use \(j=1\) and \(j'=2p-1\).
The constant term at \(j=1\) is \(-1/2\). Among \(1,\dots,2p-1\) the
only multiple of \(p\) is \(p\), whose sign in the alternating harmonic
sum is positive. Hence the constant term of the pole value at \(j'\)
has exact \(p\)-adic valuation \(-s\). Their squared-node difference is
\(4p(p-1)\), of valuation one. The divided difference therefore has
constant coefficient of valuation exactly

\[
\boxed{-s-1.}
\]

Thus the sign repair changes the local geometry from two good poles
below \(p\) to only one. Removing the polynomial correction does not
license reuse of the old row-weight estimates.

This is not restricted to a pole cancelled by the numerator: for any
surviving odd \(j<p\), the pair \(j,2p-j\), when present, has the same
valuation \(-s-1\). The second alternating harmonic sum again contains
exactly one multiple of \(p\), and the squared-node difference is
\(4p(p-j)\). Thus one may choose \(j>2N-1\) in the odd-node determinant.

At the prime two there is another explicit warning. For odd \(j\ge3\),
the unique largest power of two at most \(j\) is the unique term of
smallest valuation in \(A_j^{(s)}\). Therefore the constant term in (5)
has exact valuation

\[
v_2\!\left(-j^{s-1}A_j^{(s)}+\frac1{2j}\right)
=-s\lfloor\log_2j\rfloor.
\]

Moreover every difference of two distinct odd squares is divisible by
eight. The fixed denominator in (4) cannot be substituted for these pole
and interpolation costs.

## 5. Odd-node scaling restores the real field but exposes dyadic costs

Let

\[
E_m(t)=\prod_{j=1}^{m}(t+(2j-1)^2),\qquad
\widehat D_m(t)=\prod_{j=1}^{m}(t+(j-\tfrac12)^2).
\]

Then \(E_m(4t)=4^m\widehat D_m(t)\). Define the positive density and
functional

\[
v_s(y)=w_s(y)-\tfrac12w_s(2y)>0,\qquad
\sigma_s(R)=\int_0^\infty R(y^2)v_s(y)\,dy.
\]

Positivity follows from \(\widetilde w_s(2y)=2v_s(y)\).
At infinity \(v_s(y)\sim w_s(y)\), so it has the original \(2\pi\)
exponent. But

\[
\sigma_s(R)=\mu_s(R(t))-\tfrac14\mu_s(R(t/4)),
\qquad
\boxed{\sigma_s(t^e)=4^{-e-1}\nu_s(t^e).}                 \tag{6}
\]

Thus restoring the real decay scale converts the fixed dyadic denominator
into one that grows with the polynomial degree. More precisely, (3) gives

\[
v_2(\sigma_s(t^e))
=v_2\binom{2e+s}{s-2}-1-v_2(s-1)-2e-2
=-2e+O_s(\log(e+2)).
\]

This is an exact cost for individual moments, not by itself a lower bound
on the normalization of a determinant, where cancellations may occur.

For the source-like choice \(h=K-N\) and numerator power six, let
\(\Delta_{\nu,\mathrm{odd}}\) be the determinant with entries
\(\nu_s(E_N(t)^6t^{i+j}/E_K(t))\), and let
\(\Delta_{\sigma,\mathrm{half}}\) use \(\widehat D\) and \(\sigma_s\).
Substitution \(y=2u\) in each entry gives

\[
\nu_s(E_N^6t^{i+j}/E_K)
=4^{6N-K+i+j+1}\,
 \sigma_s(\widehat D_N^6t^{i+j}/\widehat D_K),
\]

and hence

\[
\boxed{\Delta_{\nu,\mathrm{odd}}
=4^{5Nh}\Delta_{\sigma,\mathrm{half}}.}                  \tag{7}
\]

This scalar is explicit and must be included in any normalization.
Primitive integer polynomial normalization is unchanged by a nonzero
rational scalar, so it does not itself create a better approximation.

At the leading real-energy level the odd-node family does recover the
original external field: with \(y=2K\sqrt t\), the nodes
\((2j-1)/(2K)=(j-1/2)/K\) have the same limiting uniform distribution
as the original \(j/K\), and the density contributes \(2\pi K\sqrt t\).
The midpoint log products differ from the limiting integrals by at most
\(2+2\log K\), uniformly for \(t\ge0\). Here is a direct bound, including
the logarithmic endpoint. Write \(f(u)=\log(t+u^2)\). On all cells after
the first, monotonicity bounds the absolute midpoint error by the endpoint
variation; the sum telescopes to
\(f(m/K)-f(1/K)\le2\log m\le2\log K\).
For the first cell put \(T=K^2t\). Its error is
\(\log(T+1/4)-\int_0^1\log(T+v^2)\,dv\). Since
\[
v^2\le\frac{T+v^2}{T+1/4}\le4\qquad(0<v\le1),
\]
the absolute value of this integral difference is at most two.
This identifies the common limiting field, but a complete normalized
determinant bound still needs its scalar normalization and the new local
arithmetic.

## 6. Exact small determinants: the arithmetic saving is not automatic

The script explore_euler_depleted_hankel.py evaluates three families:
original integer nodes with \(\mu_s\), integer nodes with \(\nu_s\), and
odd nodes with \(\nu_s\). It uses \(1\le K\le6\), \(N=0\), and also
\(N=1\) when \(K>1\),
and \(s=5,7\), for 66 exact determinant polynomials. All three use the
same target variable \(Z=\zeta(s)\), substituting
\(X=(1-2^{1-s})Z\) in the depleted models.

Each polynomial is converted to a primitive integer polynomial by clearing
its exact coefficient denominator and dividing its integer content.
All rational multipliers and integer coefficients are saved in
explore_euler_depleted_hankel.json. Evaluations at 250 decimal digits are
diagnostics, not certified error enclosures.

For \(K=6,N=0\), the results are:

| \(s\) | family | \(v_2\) of primitive multiplier | \(\log P(\zeta(s))\) |
|---:|---|---:|---:|
| 5 | original | 41 | −19.031930 |
| 5 | depleted, all nodes | 41 | 19.142352 |
| 5 | depleted, odd nodes | 79 | 68.885261 |
| 7 | original | 49 | 10.098049 |
| 7 | depleted, all nodes | 49 | 37.645231 |
| 7 | depleted, odd nodes | 101 | 119.221511 |

For \(s=5,K=6,N=0\), the exact raw coefficient denominators are

\[
\begin{array}{ll}
\text{original:}&2^{41}3^{19}5^8,\\
\text{depleted, all nodes:}&2^{41}3^{18}5^8\,7^3\,11,\\
\text{depleted, odd nodes:}&2^{79}3^{26}5^{19}7^{13}11.
\end{array}
\]

These explicitly show that cancellation of polynomial Bernoulli
denominators need not reduce the determinant's denominators. The pole
terms and collided-node divisions matter. None of these small experiments
is a theorem of eventual failure; \(N=0,1\) also does not test an
asymptotic sequence with \(N/K=3/40\).

The same script separately verifies 729 polynomial moments, including
sharpness of the dyadic denominator, and 360 pole identities by expanding
the two Hurwitz tails independently.

## 7. What remains worth trying

The rigorous new object is a positive single-\(\eta(s)\) functional with
no odd-prime polynomial denominators, uniformly in degree. To exploit it,
one needs a local interpolation construction adapted to alternating
harmonic poles, or a different pole set and normalization that handles
the dyadic collision cost. Neither the original arithmetic budget nor
its divided-difference congruences transfer unchanged. This is a concrete
design problem with explicit formulas, rather than an irrationality claim.

Primary analytic input: Hermite's Hurwitz-zeta formula, DLMF 25.11.29,
as proved/applied in odd-moment-generalization.md. The new density identity,
Euler-polynomial denominator argument, and scaling relations are derived
here from that model and the standard von Staudt–Clausen theorem.
