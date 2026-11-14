# Bernoulli polynomials repair powers of a Stieltjes function

Status: complete analytic/algebraic proof from the root agent's proposal,
with exact and numerical checks, independently reviewed in
[bernoulli-wilson-independent-review.md](bernoulli-wilson-independent-review.md). No
historical novelty, irrationality theorem, or useful denominator bound
is asserted. This refines the jump mechanism in
[balanced hyperbolic normality](hyperbolic-balanced-normality.md).

## 1. The normalized analytic function and its jump

For \(z>0\), put
\[
 F(z)=\pi^{-2}\operatorname{arsinh}(1/\sqrt z)^2.
\]
Its single-valued continuation to \(\mathbb C\setminus[-1,0]\) is specified by
\[
 F(z)=\pi^{-2}\int_0^1\frac{\rho(v)}{z+v}\,dv,\qquad
 \rho(v)=\operatorname{arcosh}(1/\sqrt v).
 \tag{1}
\]
An elementary integration-by-parts proof appears in the linked normality
note: differentiating the integral gives the derivative of
\(\operatorname{arsinh}(1/\sqrt z)^2\), and both vanish at infinity.
Thus the normalization and branches in (1) are unambiguous.

Define \(G(v)=\rho(v)^2/\pi^2=F(v/(1-v))\).
It decreases continuously from infinity to zero on \(0<v<1\).
The boundary values are
\[
 F_\pm(-v)=(\sqrt{G(v)}\mp i/2)^2.
 \tag{2}
\]
In particular the upper boundary value has negative imaginary part.

For a real polynomial \(P\), its Wilson divided difference is the polynomial
\[
 (\mathcal DP)(t)=
 \frac{P((\sqrt t+i/2)^2)-P((\sqrt t-i/2)^2)}{2i\sqrt t}.
 \tag{3}
\]
Although square roots appear here, the numerator is odd in \(\sqrt t\).
The result is in \(\mathbb R[t]\) and
\(\mathcal D(t^k)=kt^{k-1}+\text{lower terms}\).
This is the classical Wilson operator, not a claim that the polynomials
constructed below are themselves Wilson orthogonal polynomials.

## 2. Bernoulli polynomials invert the jump

For \(k\ge1\), define
\[
 R_k(t)=(-1)^k\left[
 B_{2k}(1/2-i\sqrt t)-B_{2k}(1/2)\right].
 \tag{4}
\]
Reflection \(B_{2k}(1-x)=B_{2k}(x)\) shows this is a real rational
polynomial in \(t\), monic of degree \(k\), with constant term zero.
In fact
\[
 R_k(t)=\sum_{j=1}^k c_{k,j}t^j,\qquad
 c_{k,j}=\binom{2k}{2j}|B_{2k-2j}(1/2)|>0.
 \tag{5}
\]
Here \(B_0(1/2)=1\) and
\(\operatorname{sign}B_{2r}(1/2)=(-1)^r\) for \(r\ge1\).
The latter sign fact also follows from the positive representation below,
so no even-zeta evaluation is needed to establish it here.
The first examples are
\[
 R_1=t,\qquad R_2=t^2+\tfrac12t,\qquad
 R_3=t^3+\tfrac54t^2+\tfrac7{16}t.
\]

The essential identity is
\[
 \boxed{\mathcal D R_k(t)=kt^{k-1}.}                \tag{6}
\]
To prove it set \(t=x^2\). The two Bernoulli arguments in the numerator
become \(1-ix\) and \(-ix\). The classical identity
\(B_{2k}(w+1)-B_{2k}(w)=2kw^{2k-1}\) gives
\[
 R_k((x+i/2)^2)-R_k((x-i/2)^2)=2ikx^{2k-1}.
\]
Dividing by \(2ix\) proves (6), including the polynomial continuation at zero.
Triangularity of \(\mathcal D\) also proves uniqueness of \(R_k\)
among polynomials with zero constant term satisfying (6).

One compact generating function is
\[
 \boxed{\sum_{k\ge1}R_k(t)\frac{u^{2k}}{(2k)!}
   =\frac{u}{2\sin(u/2)}\bigl(\cosh(u\sqrt t)-1\bigr).}
 \tag{7}
\]
This follows directly from
\(\sum B_n(a)s^n/n!=se^{as}/(e^s-1)\), taking the even part
and substituting \(s=iu\). It is a formal identity, and converges
for \(|u|<2\pi\) at each fixed \(t\).

## 3. Exact positive Stieltjes representation

**Theorem.** For every \(k\ge1\) and \(z\notin[-1,0]\),
\[
 \boxed{R_k(F(z))=\frac{k}{\pi}
     \int_0^1\frac{G(v)^{k-1/2}}{z+v}\,dv.}        \tag{8}
\]

More generally, every real polynomial \(P\) with \(P(0)=0\) satisfies
\[
 \boxed{P(F(z))=\frac1\pi\int_0^1
       \frac{\sqrt{G(v)}\,(\mathcal DP)(G(v))}{z+v}\,dv.}
 \tag{9}
\]

Proof: (2)–(3) give the upper-minus-lower jump
\[
 P(F_+(-v))-P(F_-(-v))
   =-2i\sqrt{G(v)}(\mathcal DP)(G(v)).
\]
The function vanishes as \(O(1/z)\) at infinity. At the cut endpoint
zero it is \(O(1+|\log|z||^{2\deg P})\); at \(-1\) it is bounded
with a square-root local correction. These follow either from the
continued logarithm
\(\operatorname{arsinh}(z^{-1/2})
 =\log((1+\sqrt{1+z})/\sqrt z)\), continued from positive \(z\),
or from (1) and the local boundary formulas.
Hence the small endpoint contour arcs are bounded by
\(O(\epsilon|\log\epsilon|^{2\deg P})\) and \(O(\epsilon)\).
They vanish. The density in (9) is integrable: at zero it has at most
logarithmic-power growth; at one it is \(O(\sqrt{1-v})\).
Cauchy's formula around the cut therefore gives (9), with the sign shown.
Substitution of (6) gives (8).

The proof of (6) and (8) did not use Bernoulli signs.
At infinity the left side of (8) has coefficient
\((-1)^{k+1}\binom{2k}{2}B_{2k-2}(1/2)/\pi^2\)
in front of \(1/z\). The right side has strictly positive total mass.
Taking \(k=r+1\) therefore proves
\(\operatorname{sign}B_{2r}(1/2)=(-1)^r\) directly, and hence all
the coefficient signs in (5).

The unmodified square is already a counterexample to positivity:
\(\mathcal D(t^2)=2t-1/2\). Thus \(F^2\) has a negative density where
\(G<1/4\), whereas \(F^2+F/2\) has the positive density \(2G^{3/2}/\pi\).

This diagonalizes the jump used in the balanced-normality proof.
Replacing the basis \(1,F,\ldots,F^{n-1}\) by
\(1,R_1(F),\ldots,R_{n-1}(F)\) is a fixed rational triangular change
of basis. It preserves the degree bounds on every coefficient polynomial.
The jump of \(\sum_k P_k(z)R_k(F(z))\) is simply
\(-2i\sqrt G\sum_{k\ge1}kP_k(-v)G^{k-1}\).
The degree drop is now visible without a separate triangular jump
polynomial calculation. The contour and zero-count arguments are
still needed; this change of basis does not by itself prove normality.

## 4. The entire polynomial positivity cone

For \(P(0)=0\), formula (9) gives the exact equivalence
\[
 \boxed{P(F)\text{ is the Stieltjes transform of a positive measure on }[0,1]
 \iff(\mathcal DP)(t)\ge0\quad(t\ge0).}             \tag{10}
\]
The forward implication uses uniqueness of the Cauchy transform, or
Stieltjes inversion. Formula (9) already supplies its unique finite signed
measure, so its continuous density must be nonnegative. Because \(G\)
maps \((0,1)\) onto \((0,\infty)\), this is precisely the condition in
(10); continuity includes zero. The reverse implication is immediate.

This is constructive. For \(Q(t)=\sum_{j=0}^d q_jt^j\), set
\[
 \mathcal D^{-1}_0Q(t)
      =\sum_{j=0}^d\frac{q_j}{j+1}R_{j+1}(t).
 \tag{11}
\]
Then \(P(0)=0,\ \mathcal DP=Q\), with no undetermined integration constant.
Choose arbitrary rational polynomials \(A,B\) and take
\[
 Q(t)=A(t)^2+tB(t)^2,\qquad P=\mathcal D^{-1}_0Q.
 \tag{12}
\]
This gives a rational polynomial \(P\) for which \(P(F)\) is positive
Stieltjes. No sign restrictions on the individual coefficients of \(P\)
or its Bernoulli-basis coefficients are needed.
For example,
\[
 P=\tfrac13R_3-R_2+R_1
   =\tfrac13t^3-\tfrac7{12}t^2+\tfrac{31}{48}t,
 \qquad \mathcal DP=(t-1)^2.
\]

Over **real coefficients**, (12) describes the full cone.
Here is an elementary factorization proof. A nonnegative polynomial on
\([0,\infty)\) has even multiplicity at every positive root; such factors
are squares. A negative-root factor \(t+a\), \(a>0\), equals
\((\sqrt a)^2+t\). The zero-root factor is \(t=0^2+t\,1^2\).
For a conjugate pair \(u\pm iv\), \(v\ne0\), write \(r=\sqrt{u^2+v^2}\):
\[
 (t-u)^2+v^2=(t-r)^2+2(r-u)t.
\]
Every factor is of the required form, and this form is closed under
multiplication:
\[
 (A^2+tB^2)(C^2+tD^2)=(AC-tBD)^2+t(AD+BC)^2.
\]
A positive leading constant is absorbed into the squares. This proves
the real-coefficient assertion. It does **not** assert that every
rational nonnegative \(Q\) admits this two-square representation over
\(\mathbb Q\); the rational construction above is a sufficient recipe.

## 5. Moments, rationality, and a second positive moment functional

Let \(M_{k,m}\) be the \(m\)-th moment of the measure in (8).
The substitution \(v=\operatorname{sech}^2r\), followed by integration
by parts, gives
\[
 M_{k,m}=\frac{k(2k-1)}{(m+1)\pi^{2k}}
       \int_0^\infty r^{2k-2}\operatorname{sech}^{2m+2}r\,dr.
 \tag{13}
\]
The boundary terms vanish by exponential decay at infinity and
the power \(r^{2k-1}\) at zero. In particular comparison at infinity
in (8) gives
\[
 M_{k,0}=\frac{c_{k,1}}{\pi^2}
 =\frac{\binom{2k}{2}|B_{2k-2}(1/2)|}{\pi^2}.
 \tag{14}
\]

The elementary base expansion is
\[
 F(z)=\pi^{-2}\sum_{m\ge0}(-1)^m a_m z^{-m-1},
 \quad
 a_m=\frac{4^m}{(m+1)(2m+1)\binom{2m}{m}}.
 \tag{15}
\]
Indeed \(a_m=\int_0^1v^m\rho(v)\,dv
 =(m+1)^{-1}\int_0^\infty\operatorname{sech}^{2m+2}r\,dr\);
putting \(s=\tanh r\) evaluates the latter beta integral.
Consequently all \(M_{k,m}\) are explicit polynomials in \(\pi^{-2}\)
with rational coefficients. A finite formula is
\[
 M_{k,m}=\sum_{j=1}^{\min(k,m+1)}
 (-1)^{j-1}c_{k,j}\pi^{-2j}
 \sum_{\substack{r_1+\cdots+r_j=m+1\\r_i\ge1}}
                  a_{r_1-1}\cdots a_{r_j-1}.       \tag{16}
\]
They are not, in general, rational numbers. Even after dividing by
total mass, the first moment for \(k=2\) is \(1/3-2/\pi^2\).

There is nevertheless a distinct **rational** positive moment functional
in the polynomial-design variable \(t\). Taking the coefficient of \(1/z\)
in (9), then putting \(t=(r/\pi)^2\), yields
\[
 \boxed{P'(0)=\int_0^\infty Q(t)\,d\nu(t),\quad
 d\nu(t)=\pi^2\operatorname{sech}^2(\pi\sqrt t)
                     \tanh(\pi\sqrt t)\,dt,\quad Q=\mathcal DP.}
 \tag{17}
\]
This is a probability measure. From (11),
\[
 \int_0^\infty t^j\,d\nu(t)
       =(2j+1)|B_{2j}(1/2)|\in\mathbb Q.
 \tag{18}
\]
Thus Bernoulli half-values have a direct positive moment interpretation.
It implies positive definite Hankel matrices and rational ordinary
orthogonal polynomials for \(\nu\). No useful common denominator bound
or irrationality construction is being inferred.

For each nonzero \(Q\ge0\), the measure in (9) is positive on every
subinterval apart from finitely many zero-density points. Its moment
Hankel matrices are strictly positive definite:
\(\int |q(v)|^2\,d\mu(v)>0\) for every nonzero polynomial \(q\).
Therefore the ordinary \([N-1/N]\) Padé denominators at infinity exist uniquely and
their zeros lie in \((-1,0)\), by the usual orthogonal-polynomial
sign-change proof on \((0,1)\). This is an analytic consequence,
not rational-coefficient Padé arithmetic for the measures in (8).

## 6. The quadratic repair gives exactly Euler's ζ(2) tails

Write \(H_m^{(2)}=\sum_{\ell=1}^m\ell^{-2}\), with \(H_0^{(2)}=0\).
There is a particularly simple exact formula:
\[
 \boxed{M_{2,m}
       =\frac{a_m}{2\pi^4}\bigl(\pi^2-6H_m^{(2)}\bigr).}
 \tag{19}
\]
To prove it set
\[
 U_m=\int_0^\infty\operatorname{sech}^{2m+2}r\,dr,\qquad
 V_m=\int_0^\infty r^2\operatorname{sech}^{2m+2}r\,dr.
\]
For \(q=2m+2\), the identity
\[
 (\operatorname{sech}^q r)''
 =q^2\operatorname{sech}^q r-q(q+1)\operatorname{sech}^{q+2}r
\]
and two integrations by parts give
\[
 U_{m+1}=\frac q{q+1}U_m,\qquad
 V_{m+1}=\frac{q^2V_m-2U_m}{q(q+1)}.
\]
Thus
\[
 \frac{V_{m+1}}{U_{m+1}}
 =\frac{V_m}{U_m}-\frac1{2(m+1)^2}.                \tag{20}
\]
We have \(U_0=1\). Equations (13)–(14) at \(k=2,m=0\) give
\(V_0=\pi^2/12\), without assuming Euler's evaluation of \(\zeta(2)\).
Hence \(V_m/U_m=\pi^2/12-H_m^{(2)}/2\). Using
\(U_m=(m+1)a_m\) in (13) proves (19).

Comparison with \(R_2(F)=F^2+F/2\) and (15) gives the exact rational
convolution identity
\[
 \boxed{\sum_{i=0}^{m-1}a_i a_{m-1-i}
                    =3a_m H_m^{(2)}\quad(m\ge1).} \tag{21}
\]
This identity is proved for every \(m\), not inferred from finite checks.

Positivity of \(V_m\) gives \(6H_m^{(2)}<\pi^2\).
Moreover \(V_m/U_m\to0\): after normalizing
\(\operatorname{sech}^{2m+2}r\,dr\) to mass one, its mass and
second moment concentrate at \(r=0\).
For detail, for every fixed \(\delta>0\),
the integral outside \([0,\delta]\), with factor \(r^2\), is at most
\(\operatorname{sech}^{2m}\delta\int_\delta^\infty r^2\operatorname{sech}^2r\,dr\),
whereas \(U_m\ge c/\sqrt{m+1}\) by integrating on
\([0,1/\sqrt{m+1}]\). The contribution inside is at most \(\delta^2\).
First let \(m\to\infty\), then \(\delta\downarrow0\).
Equation (20) consequently recovers the classical identity
\(\zeta(2)=\pi^2/6\). This derivation uses the explicitly checked \(R_2\),
so it does not presuppose an even-zeta formula to establish the needed
general Bernoulli signs.

## 7. The measure is classical Wilson, and Apéry is its diagonal

There is an exact identification beyond the name of the divided difference.
Put \(t=y^2\) in (17). Gamma reflection gives
\[
 |\Gamma(1/2+iy)|^2=\frac{\pi}{\cosh\pi y},\qquad
 |\Gamma(2iy)|^2=\frac{\pi}{2y\sinh2\pi y}.
\]
Consequently
\[
 \boxed{d\nu(y^2)=
  2\pi^2y\,\operatorname{sech}^2(\pi y)\tanh(\pi y)\,dy
  =\frac1{2\pi}\left|
       \frac{\Gamma(1/2+iy)^4}{\Gamma(2iy)}
                \right|^2dy.}                    \tag{22}
\]
This is the classical Wilson probability measure with all four parameters
equal to \(1/2\). The weight and its orthogonality are given in
Cheng–Chiang, §6.1, with original attribution to J. A. Wilson (1980).
They are also recorded with explicit measure conventions in
[DLMF 18.25.2–18.25.5](https://dlmf.nist.gov/18.25.ii).

The standard Wilson polynomial in this specialization is the terminating sum
\[
 W_n(t)=(n!)^3
 \sum_{k=0}^n
 \frac{(-n)_k(n+1)_k}{(k!)^4}
       \prod_{j=0}^{k-1}\bigl(t+(j+1/2)^2\bigr).
 \tag{23}
\]
Its monic normalization is
\(q_n(t)=(-1)^n n!W_n(t)/(2n)!\).
Classical Wilson orthogonality therefore supplies explicit rational
orthogonal polynomials for the Bernoulli moments (18). For example
\(q_1=t-1/4\), \(q_2=t^2-3t/2+11/48\).
This identification is established classical mathematics, not a new
family of orthogonal polynomials.

For integers \(m,n\ge0\), evaluate at a **negative** half-integer square:
\[
 \prod_{j=0}^{k-1}\bigl((j+1/2)^2-(m+1/2)^2\bigr)
           =(-m)_k(m+1)_k .
\]
Thus, purely by the finite sum (23),
\[
 \boxed{\frac{W_n(-(m+1/2)^2)}{(n!)^3}
 =\sum_{k\ge0}\binom nk\binom{n+k}k
                  \binom mk\binom{m+k}k=K(n,m),}   \tag{24}
\]
where the sum terminates at \(\min(n,m)\).
This is exactly the kernel in
[our bispectral derivation](newton-bispectral-motivation.md). In particular
\[
 \boxed{A_n=\frac{W_n(-(n+1/2)^2)}{(n!)^3}.}        \tag{25}
\]
The orthogonality variable is positive, while this evaluation point is
negative and moves with \(n\). Hence (25) does not justify applying an
ordinary fixed-variable three-term recurrence directly on the diagonal.
The separate finite-grid closure argument in the linked note remains needed.

For completeness, the classical Wilson norm in this normalization is
\[
 \int_0^\infty W_n(t)^2\,d\nu(t)=\frac{(n!)^6}{2n+1},
 \qquad
 \int q_n(t)^2\,d\nu(t)
       =\frac{(n!)^8}{(2n+1)((2n)!)^2}.            \tag{26}
\]
It follows by substituting the four half-parameters in DLMF 18.25.5,
and dividing its unnormalized weight by \(2\pi\), exactly as in (22).
Therefore the Bernoulli Hankel determinants have the closed product
\[
 \det\bigl[(2i+2j+1)|B_{2i+2j}(1/2)|\bigr]_{i,j=0}^{N-1}
    =\prod_{j=0}^{N-1}\frac{(j!)^8}{(2j+1)((2j)!)^2}. \tag{27}
\]
This uses the standard Gram determinant/product-of-monic-norms identity.
The finite kernel identity (24) is proved directly here; the all-degree
Wilson orthogonality/norm theorem is explicitly imported from its sources.

An elementary alternative proof of these particular orthogonality and norm
identities is now available in the
[resolvent–Legendre proof](resolvent-legendre-pade-proof.md), §§1–2.
It proves \(\nu(Q_jQ_k)=1/(j+k+1)\) by a finite polynomial-difference
identity and transports the Legendre inner product. That later proof
does not require the general Wilson orthogonality theorem.

## 8. Source context and verification

The Bernoulli difference/reflection identities are classical and follow
from their generating function. The operator (3) is exactly Definition 2.1
of K. H. Cheng and Y.-M. Chiang, *Nevanlinna theory of the Wilson
divided-difference operator*, Ann. Acad. Sci. Fenn. Math. 42 (2017),
175–209, DOI 10.5186/aasfm.2017.4211.
[Primary journal PDF](https://www.acadsci.fi/mathematica/Vol42/vol42pp175-209.pdf).
The paper is cited for the operator and its established context, not as
a source for the particular Stieltjes representation proved here.
Its §6.1 also supplies the Wilson weight and its orthogonality.
Its PDF/text and [SHA-256 provenance](../sources/cheng-chiang-wilson-2017.json)
are saved. No unsuccessful search is evidence of novelty.

Original attribution: J. A. Wilson, *Some hypergeometric orthogonal
polynomials*, SIAM Journal on Mathematical Analysis 11 (1980), 690–701,
DOI [10.1137/0511064](https://epubs.siam.org/doi/10.1137/0511064).
The title, author, journal and pages are given by the retrieved
Cheng–Chiang paper; the DOI was checked against the SIAM publisher's
search result. The original paper was not separately retrieved here.
The norm normalization was checked against the displayed formulas in
NIST DLMF §18.25, rather than assumed from the word “Wilson.”
An [exact transcription of that retrieval](../sources/dlmf-wilson-norm-retrieval.md)
records the measure and norm conventions without claiming an archived HTML file.

[Checker](../computations/check_bernoulli_wilson_stieltjes.py):
12 exact Bernoulli/Wilson identities, eight generating coefficients,
54 moments compared with an independent sech-power recurrence,
the nontrivial positivity-cone example, and direct high-precision
Stieltjes quadratures. The expanded checks also include 100 exact
Euler-tail convolution identities, six rational-density moments,
117 Wilson-grid evaluations, exact orthogonality and norms through
degree eight, and six Hankel determinant products.
The analytic endpoint and uniqueness arguments
are proofs in this note, not conclusions from quadrature.
