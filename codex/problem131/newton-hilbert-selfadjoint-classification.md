# Why a quadratic Newton grid and cubic recurrence arise from the Hilbert Gram matrix

Status: the elementary classification and its spectral dictionary passed the
[independent review](newton-hilbert-classification-independent-review.md);
exact checks pass. The spectral family is classical:
a primary-source parameter match with Ismail–Koelink is given in section7.
No historical novelty is asserted. The elementary classification under the
specified Newton/Gram assumptions is proved directly, not imported from that
spectral theorem.

The design mechanism is short: prescribe the simplest integral Gram matrix;
require multiplication by the spectral variable to be symmetric; this forces
a quadratic Newton lattice. Legendre polynomials then make that multiplication
tridiagonal, and the forced quadratic coefficient becomes a cubic coefficient.
The usual Apéry binomial sums appear by evaluating this change of basis on its
own Newton lattice.

## What the design forces, and what it leaves free

Start with the Legendre integral on \([0,1]\). Its monomials have Gram
matrix \(1/(j+k+1)\), so preserving that matrix would let us transport
Legendre orthogonality to a new polynomial family. To keep the new
polynomials explicit, ask for Newton products: multiplication by \(t\)
should add just one factor,
\[
 tQ_j=a_jQ_{j+1}-\lambda_jQ_j.
\]
These are the design assumptions, rather than conclusions about every
irrationality proof.

There is a crucial consistency test. If the transported inner product is to
come from a single real moment functional, multiplication by \(t\) must be
symmetric:
\[
 \langle tQ_j,Q_k\rangle
 =L(tQ_jQ_k)=\langle Q_j,tQ_k\rangle.
\]
Putting the Hilbert entries into this equality turns an integral requirement
into a relation between \(a_j\) and \(\lambda_j\). Comparing with just the
constant and linear basis elements, \(k=0,1\), determines every later
coefficient. Section 1 proves that the answer has to be quadratic.
The converse matters too: these forced coefficients really do produce a
positive moment functional, not merely a consistent-looking recurrence.

The remaining choices have different roles:

| Input or additional requirement | Exact consequence |
|---|---|
| Hilbert Gram matrix and symmetric Newton multiplication | \(a_j=a_0+cj(j+2)\), \(\lambda_j=\lambda_0+cj(j+1)\) |
| Evaluation by the Legendre binomial product | \(a_0=c\), hence \(a_j=c(j+1)^2\) |
| Choose unit positive scale and the zero spectral edge | \(c=1,\lambda_0=1/4\), giving the half-square grid |
| The separately proved adjacent balance | Diagonal closure with the coefficients \(17\) and \(5\) |

Why do cubes then appear? In the Legendre basis, multiplication by \(x\)
has upper coefficient \(-(n+1)/(2(2n+1))\).
The symmetric differential operator multiplies it by the forced quadratic
\(a_n=(n+1)^2\). Their product is exactly
\[
 -\frac{(n+1)^3}{2(2n+1)}.
\]
Evaluation on the Newton grid gives the two-index Apéry recurrence; the
[reviewed mixed telescopes](apery-mixed-numerator-independent-review.md)
supply the additional adjacent balance and transfer it to the full linear
form. This separates the origin of the cubic outer terms from the extra
identity needed for \(17=6\cdot3-1\) and \(5=6\cdot1-1\).

Here “symmetric” means the polynomial inner-product identity just displayed.
Essential selfadjointness and the precise spectral support are stronger
analytic assertions, attributed to Ismail–Koelink in §7. The quarter shift
is not forced by symmetry or by the integer binomial evaluations. This
design narrative reconstructs a mathematical route to the formulas; it
makes no claim about their historical discovery.

## 1. Classification for the uniform Hilbert Gram matrix

Let Q₀=1 be a real polynomial Newton basis satisfying
\[
 tQ_j=a_jQ_{j+1}-\lambda_jQ_j,\qquad a_j\ne0\quad(j\ge0).
\]
Assume its inner products are
\[
 \langle Q_j,Q_k\rangle=\frac1{j+k+1}.
\]
**Theorem1.** Multiplication by t is symmetric if and only if
\[
 \boxed{a_j=a_0+cj(j+2),\qquad
 \lambda_j=\lambda_0+cj(j+1),\qquad c=(a_1-a_0)/3.}      \tag{1}
\]
The three initial constants are arbitrary real numbers subject to the stated
nonvanishing of every aⱼ. In particular c=0 is a valid, separately interpretable
case; there is no implicit assumption c>0.

Proof: comparing ⟨tQⱼ,Qₖ⟩ and ⟨Qⱼ,tQₖ⟩ gives
\[
 (a_j-a_k)(j+k+1)=(\lambda_j-\lambda_k)(j+k+2).          \tag{2}
\]
The k=0 equation gives
λⱼ−λ₀=(aⱼ−a₀)(j+1)/(j+2). Insert it into the k=1 equation and put
c=(a₁−a₀)/3. After multiplication by j+2, the coefficient of aⱼ−a₀ is
(j+2)²−(j+1)(j+3)=1. The result is aⱼ−a₀=cj(j+2), and then λ follows.
This works for every j≥0 without division by j or j−1. Conversely substituting
(1) into (2) proves every pairwise symmetry condition.

### Converse positivity is exact, but support requires more

For any real sequences (1) with all aⱼ nonzero, construct Q recursively and
map UQⱼ=xʲ. This is a linear bijection between real polynomial spaces. Define
⟨p,q⟩=∫₀¹(Up)(Uq)dx. It is strictly positive definite and has the prescribed
Hilbert Gram matrix. Symmetry from (2) gives
⟨tʳ,tˢ⟩=⟨1,tʳ⁺ˢ⟩. Hence with L(p)=⟨p,1⟩ one obtains
\[
 \boxed{L(pq)=\langle p,q\rangle,\qquad L(p^2)>0
                                  \quad(p\ne0).}       \tag{3}
\]
Thus it is a positive polynomial moment functional; the standard Hamburger
moment theorem supplies a positive real-line representing measure. Symmetry
alone does not place that measure on [0,∞), and this proof makes no generic
support or uniqueness assertion.

If c=0, then Qⱼ=((t+λ₀)/a₀)ʲ and the representing measure is explicitly the
pushforward of uniform x∈[0,1] under t=a₀x−λ₀. Its support is the interval
with endpoints −λ₀ and a₀−λ₀; nonnegative support is equivalent to
λ₀≤min(0,a₀). If some aⱼ vanishes, the full Newton-basis assumption fails;
possible finite-dimensional degenerations are not silently included in this
theorem.

## 2. The forced operator is a Sturm–Liouville operator

Under U, multiplication by t becomes
\[
 \boxed{T=-c\frac d{dx}\left[x^2(1-x)\frac d{dx}\right]
                          +a_0x-\lambda_0.}            \tag{4}
\]
Indeed its action on xʲ is
\(a_jx^{j+1}-\lambda_jx^j\), exactly (1). Integration by parts is legitimate
on polynomials since x²(1−x) vanishes at both endpoints. This independently
explains the symmetry and turns the Newton multiplication rule into a real
differential operator on a familiar interval.

For a₀=c>0 and λ₀=c/4,
\[
 \boxed{\langle Tf,f\rangle=
 c\int_0^1x^2(1-x)\left(f'+\frac f{2x}\right)^2dx
                  +\frac c4\int_0^1xf^2dx\ge0.}       \tag{5}
\]
The cross-term integration has boundary x(1−x)f², which vanishes even when
f(0)≠0. The square is integrable at0 because its x⁻² factor is canceled by
x². Expansion gives precisely c∫x²(1−x)f′²+c∫(x−1/4)f².

This selects a natural zero-edge normalization. It is **not** a consequence
of (2), which leaves λ₀ free. Smaller λ₀ adds a positive constant to T.
The elementary finite-polynomial argument in
[newton-hilbert-critical-quarter.md](newton-hilbert-critical-quarter.md),
accepted in [its independent review](newton-critical-quarter-independent-review.md),
proves that c/4 is the maximal allowed shift for positivity of multiplication
on polynomial squares in this particular a₀=c>0 family. It therefore excludes
nonnegative-support measures above that shift without using a spectral theorem.
This is a maximal shift, rather than a minimal shift.

## 3. Integer Legendre-binomial evaluations force the remaining scale choice

The Newton products are
\[
 Q_k(t)=\prod_{r=0}^{k-1}\frac{t+\lambda_r}{a_r}.
\]
Suppose one additionally wants the exact evaluation rule
\[
 \boxed{Q_k(-\lambda_n)=(-1)^k\binom nk\binom{n+k}k.}   \tag{6}
\]
Its k=1,n=1 instance gives −2c/a₀=−2, so a₀=c≠0. Conversely, if a₀=c,
then aᵣ=c(r+1)², and
\[
 \lambda_r-\lambda_n=c(r-n)(r+n+1),
\]
which proves (6) by taking the product. For k>n both sides are zero. The
shift λ₀ cancels completely; the binomial rule does not fix it.

With the convenient positive scale c=a₀=1 and zero-edge shift λ₀=1/4,
\[
 \boxed{Q_k(t)=\frac{\prod_{r=0}^{k-1}[t+(r+1/2)^2]}{(k!)^2}.} \tag{7}
\]
Thus the half-integer-square grid is not guessed inside this design problem:
the quadratic form is forced by Hilbert symmetry, a₀=c is forced by the
Legendre binomial evaluation, and the remaining quarter shift is a selected
spectral normalization.

## 4. Legendre tridiagonalization gives cubic coefficients directly

Let Lₙ(x)=Pₙ(1−2x), so
\[
 L_n=\sum_{j=0}^n(-1)^j\binom nj\binom{n+j}j x^j,
\qquad\int_0^1L_nL_mdx=\delta_{nm}/(2n+1).
\]
The standard Legendre operator J=−D[x(1−x)D] satisfies JLₙ=n(n+1)Lₙ.
An exact anticommutator identity is
\[
 \boxed{T=\frac c2(xJ+Jx)+(a_0-c)x+c/2-\lambda_0.}      \tag{8}
\]
Since multiplication by x is tridiagonal in the Legendre basis, so is T.
Using the elementary Legendre recurrence gives
\[
 TL_n=-\frac{(n+1)a_n}{2(2n+1)}L_{n+1}
 +\left(\frac{cn(n+1)}2+\frac{a_0}2-\lambda_0\right)L_n
 -\frac{na_{n-1}}{2(2n+1)}L_{n-1}.                      \tag{9}
\]
For n=0 omit the last term. The cubic coefficients in the chosen case are
therefore products of a linear Legendre coefficient with the forced quadratic
aₙ=(n+1)²; no recurrence fitting is involved.

Put Rₙ(t)=U⁻¹Lₙ. For (7),
\[
 \boxed{tR_n=-\frac{(n+1)^3}{2(2n+1)}R_{n+1}
       +\frac{2n(n+1)+1}4R_n
       -\frac{n^3}{2(2n+1)}R_{n-1}.}                   \tag{10}
\]
These polynomials have norm1/(2n+1) for L and are finite terminating Wilson
polynomials, as identified independently in section7.

### The Apéry sum is a diagonal overlap on the Newton grid

Evaluating Rₙ at −λₘ and using (6) gives
\[
 K(n,m)=R_n(-(m+1/2)^2)
 =\sum_{j=0}^{\min(n,m)}
 \binom nj\binom{n+j}j\binom mj\binom{m+j}j.           \tag{11}
\]
In particular K(n,n)=Aₙ, exactly the Apéry integer sum. Equation (10) yields
the all-n mixed recurrence
\[
 (n+1)^3K(n+1,m)
 -(2n+1)[n(n+1)+2m(m+1)+1]K(n,m)
 +n^3K(n-1,m)=0.                                      \tag{12}
\]
Combined with the separately proved adjacent balance
K(n,n−1)=(Aₙ+Aₙ₋₁)/6, this gives the usual Apéry recurrence immediately.
Indeed at m=n, the middle coefficient is (2n+1)[3n(n+1)+1], and eliminating
the two adjacent overlaps leaves (2n+1)[17n(n+1)+5].
The classification/tridiagonalization explains the quadratic grid, the cubic
outer terms, and the mixed recurrence; the adjacent balance remains an
explicit additional identity, not a consequence silently assumed from Gram
symmetry alone.

## 5. The same classification holds for every beta Gram matrix

Let A>0, B>A and κ>0. Prescribe
\[
 \langle Q_j,Q_k\rangle
 =\kappa\,\mathrm B(j+k+A,B-A)
 =\kappa\int_0^1x^{j+k+A-1}(1-x)^{B-A-1}dx.
\]
**Theorem2.** Symmetric Newton multiplication is equivalent to
\[
 \boxed{a_j=a_0+cj(j+B),\quad
 \lambda_j=\lambda_0+cj(j+A),\quad c=(a_1-a_0)/(B+1),}   \tag{13}
\]
again assuming every aⱼ≠0.

The ratio of successive beta moments is (j+k+A)/(j+k+B), so symmetry says
\[
 (a_j-a_k)(j+k+A)=(\lambda_j-\lambda_k)(j+k+B).
\]
For d=aⱼ−a₀, the k=0 and k=1 equations reduce to
\[
 \frac{B-A}{j+B}[d-cj(j+B)]=0.
\]
The factor B−A>0 is exactly the nondegeneracy needed. This proves (13);
substitution proves the converse. Positivity follows by the same weighted
isometry argument as in (3).

With w=xᴬ⁻¹(1−x)ᴮ⁻ᴬ⁻¹, the operator is
\[
 T=-\frac c{w}D[x^2(1-x)wD]+a_0x-\lambda_0.             \tag{14}
\]
For the Jacobi operator J=−w⁻¹D[x(1−x)wD],
\[
 \boxed{T=\frac c2(xJ+Jx)+(a_0-cB/2)x+cA/2-\lambda_0.}  \tag{15}
\]
Thus every beta version is tridiagonal in its Jacobi basis. For reference,
its weighted square completion is
\[
 \begin{aligned}
 \langle Tf,f\rangle_w={}&c\int x^2(1-x)
                         (f'+Af/(2x))^2w\,dx\\
 &+\int\left[(a_0-cA(2B-A)/4)x+cA^2/4-\lambda_0\right]f^2w\,dx.
 \end{aligned}                                        \tag{16}
\]
The boundary x(1−x)wf² vanishes because A>0 and B−A>0. This is a useful
sufficient positivity test when c>0; no exhaustive support classification is
being inferred from it.

## 6. The beta(2,2) companion is the positive Wilson-one model

For A=2,B=4, a₀=4,c=1,λ₀=1 one has
\[
 a_j=(j+2)^2,\quad\lambda_j=(j+1)^2,\quad
 Q_k(t)=D_k(t)/((k+1)!)^2.
\]
With κ=1/2, its Gram matrix is
\[
 \boxed{\langle Q_j,Q_k\rangle
              =\tfrac12\mathrm B(j+k+2,2).}            \tag{17}
\]
The representing measure is precisely t dμ₃(t), the positive Wilson measure
with all four parameters1, with the previously established normalization.
Its mass is1/12, agreeing with (17). Thus the Legendre/all-half model and the
Jacobi(1,1)/all-one model sit in the same elementary Gram/Newton design.
The identification uses the exact operator and Wilson parameters in §7;
the conventional all-one Wilson measure is twice \(t\,d\mu_3\), so its
mass fixes the factor \(1/2\). It does not rely on the finite Gram checks alone.
This explains the derivative companion used in the moving-pole upper bounds.
It does not identify the separate original D₂ₕ optimized norm with that
moving-pole norm; the unresolved bridge remains unresolved.

## 7. Exact classical Wilson identification and primary-source comparison

For c≠0 there is a direct finite-series identification, not just a resemblance
of recurrence coefficients. Factor
\(a_j=c(j+\eta)(j+B-\eta)\), so η(B−η)=a₀/c. Set
\(Y^2=(t+\lambda_0-cA^2/4)/c\). Then
\[
 Q_j(t)=\frac{(A/2+iY)_j(A/2-iY)_j}{(\eta)_j(B-\eta)_j}.
\]
The normalized Jacobi polynomial
\({}_2F_1(-n,n+B-1;A;x)\), transported through U⁻¹, is
\[
 {}_4F_3\left(\begin{matrix}-n,n+B-1,A/2+iY,A/2-iY\\
                         A,\eta,B-\eta\end{matrix};1\right).
\]
These are Wilson polynomials with parameter multiset
\[
 (A/2,A/2,\eta-A/2,B-\eta-A/2),
\]
up to their conventional scalar normalization. This algebraic identity is
valid whenever the Newton denominators are nonzero; positivity of a Wilson
integral for arbitrary parameters is a separate matter.

The primary paper retrieved and read is M.E.H. Ismail and E. Koelink,
*Spectral properties of operators using tridiagonalisation*,
[arXiv:1108.5716](https://arxiv.org/abs/1108.5716), §3, equation(3.1) and
Theorem3.1. Local PDF/text copies are in `codex/sources/ismail-koelink-11085716.*`.
Here is the exact parameter dictionary, including the operator scale:

    paper coordinate ξ=1−2x;
    α_paper=A−1, β_paper=B−A−1, δ_paper=η−A;
    T_ours=−(c/2)T_paper−λ₀.

Changing variables in their displayed differential operator proves this
identity term by term; the checker verifies it exactly. In particular:

- the uniform chosen case has α=β=δ=0 and Wilson parameters(1/2)^4;
- the beta(2,2) companion has α=β=1,δ=0 and Wilson parameters1^4.

Their theorem supplies the essentially selfadjoint spectral realization in
these positive-parameter cases. The all-half case has no discrete extra masses.
For c=1,λ₀=1/4 its probability spectral measure in y=√t is

    2π² y sech²(πy)tanh(πy)dy,   y>0.

The gamma identities identify this exactly with the ordinary Wilson all-half
weight; its mass is1. More generally a₀=c>0 gives the pushforward

    t=cY²+c/4−λ₀.

Its support has lower edge c/4−λ₀, confirming the selected zero-edge shift.
The maximal-shift statement in section2 now has an independently reviewed
elementary polynomial proof; the classical theorem remains the source for the
spectral realization and density, not a required input to maximality. Gram
symmetry alone still does not select that shift. No claim that this is a new
Sturm–Liouville spectral family is appropriate. The elementary two-anchor classification and
its design interpretation are the independent derivation recorded here.

## Reproduction and remaining work

`check_newton_hilbert_classification.py` checks the general classification,
weighted operator/anticommutator/square identities, the exact primary-source
map, the general uniform Legendre tridiagonal action through degree8, the
half-grid Hilbert moments and mixed binomial evaluations, and the beta(2,2)
companion Gram. All checks pass. The all-j classification and positivity proof
are analytic/algebraic arguments above; finite checks do not replace them.

The next mathematical tasks concern arithmetic and optimized norms beyond
these classical isometries. The original D₂ₕ Schur-root/Apéry identification
has not been proved by this classification. The useful result here is a clear,
reproducible answer to why quadratic Newton factors and cubic mixed recurrence
coefficients are natural once the uniform beta integral is chosen.
