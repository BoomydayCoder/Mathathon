# Independent review: Bernoulli antidifferences for the balanced system

**Status: accepted for the parameterization, moment reductions, and
normalizations. The proposed factorial constant-term formula remains a
conjecture.** This reviews
[hyperbolic-balanced-bernoulli-construction.md](hyperbolic-balanced-bernoulli-construction.md).
The separate [all-balanced normality audit](hyperbolic-balanced-normality-independent-review.md)
is complete, but does not prove explicit coefficients or denominator bounds.
No historical novelty or irrationality claim is made.

## 1. The moment identity, including signs and shifts

Fix integers \(m,n\ge1\), set \(d=mn\), \(N=2m+1\), and let \(p\) have
degree at most \(d\). Define \(R(x)=x^Np(-x^2)\).
For \(1\le j\le m\) and any integer \(k\ge0\), the claimed identity is

\[
 \mu_{2j+1}(t^{m-j+k}p(t))
 =\frac{(-1)^{m-j+k}}{(2j)!}
 \mathcal B\!\left(D^{2j-1}[x^{2k}R(x)]\right).
\]

It follows by checking \(p=t^e\). Put \(\ell=m-j+k+e\). The sign on the
right is \((-1)^{m-j+k+e}=(-1)^\ell\); its derivative coefficient is

\[
 \frac{(2m+2e+2k+1)!}{(2j)!(2\ell+2)!},
\]

and its Bernoulli index is \(2\ell+2\). These are precisely the defining
moment \(\mu_{2j+1}(t^\ell)\). Thus the identity is valid before imposing
any vanishing condition. All exponents in these polynomial moments are
nonnegative.

The sign convention \(B_1=-1/2\) gives
\[
 \mathcal B(H(x)-H(x-1))=H'(-1),
\]
not \(H'(0)\). This shift is handled correctly in the construction.

## 2. Antidifference and the exact first-layer parameterization

The difference map \(Q(x)\mapsto Q(x)-Q(x-1)\) lowers each nonconstant
monomial degree by one with nonzero leading coefficient. It is surjective
onto polynomials, with constants as its kernel. Consequently \(R\) has a
unique polynomial antidifference satisfying \(Q(0)=0\).

Since \(R(0)=0\), also \(Q(-1)=0\). Oddness of \(R\) shows that
\(\widetilde Q(x)=Q(-x-1)\) has the same difference and the same normalization
as \(Q\). Hence \(\widetilde Q=Q\). With \(u=x+1/2\), this invariance is
evenness in \(u\), so the invariant ring is exactly
\(\mathbb Q[x(x+1)]\) when \(p\) has rational coefficients.

Commuting differentiation with finite difference gives
\[
 \mathcal B(D^{2j-1}R)=Q^{(2j)}(-1)=Q^{(2j)}(0).
\]
Thus the first \(m\) moment conditions, namely \(k=0\), are precisely the
even derivative vanishings of orders \(2,4,\ldots,2m\).
Reflection also gives \(Q(x-1)=Q(-x)\). Since
\(R=Q(x)-Q(-x)\) is divisible by \(x^N\), its odd derivatives of orders
\(1,3,\ldots,2m-1\) vanish at zero. Together with \(Q(0)=0\), all
derivatives through order \(N-1\) vanish.

Reflection now supplies the same multiplicity at \(-1\). The two factors
are coprime, so
\[
 Q(x)=[x(x+1)]^N V(x(x+1)).
\]
Conversely, any such \(Q\) has the required reflection and endpoint zeros;
its odd part divided by \(x^N\) is even and therefore defines a polynomial
\(p(-x^2)\). It satisfies every \(k=0\) condition.
The degree bound is
\[
 \deg Q\le2d+2m+2,\qquad
 \deg V\le d-m=m(n-1).
\]
Both directions of the parameterization are therefore proved, including
polynomials of smaller degree. No normality theorem was used here.

## 3. Monic normalization and constant term

Put \(D=d+m+1\). If \(p\) has degree exactly \(d\) and is monic, then \(R\)
has leading coefficient \((-1)^d\) and degree \(2D-1\). An antidifference
of degree \(2D\) multiplies its leading coefficient by \(2D\). Thus
\[
 [z^{m(n-1)}]V(z)=\frac{(-1)^d}{2D}.
\]
There is no factorial factor missing in this normalization.

Since \(z=x(x+1)=x+O(x^2)\), the coefficient of \(x^N\) in \(Q\) is
\(V(0)\). The integer \(N\) is odd, so its coefficient in
\(Q(x)-Q(-x)\) is twice that value. Therefore
\[
 p(0)=2V(0).
\]
For \(n=1\), \(V\) is a constant, and these formulas recover the previously
reviewed explicit degree-\(m\) family. The case \(m=n=1\) gives
\(p(t)=t-1/3\), correctly orthogonal to constants for \(\mu_3\).
The stated \(m=n=2\) example and its scale \(1/14\) are also consistent.

The conjecture for the monic normalization \(\widehat V\) is algebraically
equivalent to the conjecture for \(p(0)\): writing
\((m+1)n+m=D+n-1\), the factor \(D\) cancels from the binomial denominator.
This equivalence is proved; the values conjectured on either side are not.

## 4. The next moment layer retains a real constraint

With the parameterized \(Q\), the identity
\[
 x^2R=\Delta_-[x^2Q]-(2x-1)Q(x-1)
\]
is exact. Apply \(\mathcal B D^{2j-1}\). The difference term becomes
\((x^2Q)^{(2j)}(-1)=0\), since \(2j\le2m<N\).
For \(H(x)=(2x+1)Q(x)\), the other term is \(H(x-1)\), and
\[
 \mathcal B(D^{2j-1}H(x-1))
 =\mathcal B(D^{2j-1}H(x))-H^{(2j)}(-1).
\]
The last term is again zero by the same multiplicity. This proves the
negative sign in the note's equation (7), without dropping a boundary term.

Arbitrary choices of \(V\) do not generally satisfy this next layer. In
particular the proposed ordinary-difference polynomial
\(t^2-7t/4+1/4\) has zeroth \(\mu_3\) moment zero but first moment
\(1/40\). Thus the parameterization does not accidentally prove a Rodrigues
formula or the higher orthogonality conditions.

After monic normalization, \(V\) has \(m(n-1)\) free coefficients, matching
the remaining \(k=1,\ldots,n-1\) equations. The separately reviewed
all-balanced normality theorem ensures that this reduced affine system has
exactly one solution. It does not evaluate that solution's constant term.

## 5. Checks, source scope, and unresolved claim

The original checker passed all seventeen listed systems: \(1\le m,n\le4\)
and \((m,n)=(5,2)\). It checked rational equations, Sturm counts, exact
antidifferences, normalizations, and the finite constant-term pattern.
Only its final JSON write was suppressed during this review; no mathematical
test was altered.

The [independent checker](../computations/check_balanced_bernoulli_independent.py)
uses nine arbitrary \(V\)'s, rather than only fully solved polynomials.
It checks reflection, recovery, monic normalization, the constant term, all
first-layer conditions, and the next-layer identity. It also checks the
moment formula on arbitrary \(p\), where neither side generally vanishes.
There are 54 such checks, and the constructed examples retain nonzero
higher constraints. This guards against a vacuous comparison of two zero
moments. All checks passed.

The archived Beckermann–Coussement–Van Assche primary source, §3.2 and
Theorem 3.3, varies one gamma parameter while keeping three fixed and
imposes its stated restrictions on parameter differences. Its usual
positive Wilson weights have the same exponential tail. The present
hyperbolic components have different tail exponents. Therefore the note
correctly declines to treat its system as an immediate substitution into
that theorem. This is not a claim excluding a subtler transform or an
existing result elsewhere.

The formula
\[
 p_{m,n}(0)\stackrel{?}{=}
 \frac{(-1)^{mn}(n!)^{2m}}{\binom{(m+1)n+m}{n}}
\]
remains finite exact evidence. General normality, positivity of roots, and
the proved identity \(p(0)=2V(0)\) do not establish it. A determinant
evaluation or a genuine raising relation remains necessary.
