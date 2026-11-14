# Independent audit of the geometric Apéry operator and real boundary

**Verdict:** after the correction to the discriminant estimate in Section 5 of
[geometric-divergence-proof.md](geometric-divergence-proof.md), the algebraic
certificates, torus integration, rational lift, improper real-chain limits,
and boundary normalization are valid. No unresolved mathematical gap was found.
In particular the boundary calculation proves

\[
\mathcal L_t H(t)=-12,\qquad \mathcal L_uF(u)=-12u.
\]

The previously claimed inequality \(a_*-b_*\ge17+12\sqrt2\) was false. At
\(x=y=R\), its left side is \((1+1/R)^2\to1\). This review accepts the
replacement \(a_*-b_*>1\), not the old assertion. All downstream estimates
below are checked using the repaired bound alone.

This is an exact computer-assisted geometric reduction with an analytic proof
of its boundary limits. It is not a short hand derivation of the large
certificates, a minimal-order theorem, or a historical novelty claim.

## 1. The torus branch is single valued where it is used

Write \(s=xy\), \(p=(1+x)(1+y)\), \(r=1+x+y\), and

\[
\Lambda=pz+p+p^2/s+p^2/(sz),\qquad
f=sz(1-t\Lambda)=az^2+bz+c.
\]

Thus \(dx\,dy\,dz/f\) is exactly
\(d\log x\,d\log y\,d\log z/(1-t\Lambda)\), with no missing \(s\) or \(z\)
factor. Fix nonzero radii for all three torus circles and a compact product
of annuli around them. For sufficiently small \(|t|\), the geometric series
converges uniformly there and \(D/s^2\) is uniformly close to \(1\).
Define the branch by

\[
\sqrt D=s\left(D/s^2\right)^{1/2},
\]

using the binomial branch near \(1\). It is single valued and analytic on the
two annuli and satisfies \(\sqrt D=s\) at \(t=0\). There is no claim that this
branch extends across \(x=0\) or \(y=0\).

Constant-term extraction in \(z\), or the convergent central-binomial series,
then gives

\[
\operatorname{CT}_z(1-t\Lambda)^{-1}=s/\sqrt D.
\]

The remaining form is \(dx\,dy/\sqrt D\), with normalized integral
\((2\pi i)^{-2}\) times the double contour integral. The same branch defines
\(D^{5/2}\) in the primitives \(U/D^{5/2}\) and \(V/D^{5/2}\).
Each is single valued along the corresponding closed circle. Integrating its
ordinary derivative gives zero, without requiring analyticity throughout the
disk bounded by that circle. Compactness also justifies parameter derivatives.

The torus argument only needs this two-variable algebraic primitive. The
later rational \(q_v\) formulas can have extra denominators at \(p=0\), which
may occur on unit circles. They need not be integrated on those circles.
Their real-chain use is on \(x,y,t>0\), where \(a=-tsp\ne0\).

## 2. What the coefficient calculation does establish

The recursive identity
\(\theta^kD^{-1/2}=N_kD^{-k-1/2}\) follows by one differentiation. Clearing
denominators in the proposed divergence gives exactly the coefficient equation
in the construction. There are 405 primitive coefficients and eight unknown
operator coefficients.

I read and reran the discovery program with file output alone suppressed.
It returned 1,180 equations, augmented matrix width 414, rank 383, and
30 free primitive coefficients. Its pivot checks establish that no operator
coefficient depends on a free variable. The result is

\[
\mathcal L=\theta^3
-t(2\theta+1)(17\theta^2+17\theta+5)+t^2(\theta+1)^3.
\]

The program does not insert these eight coefficients before row reduction.
Its conversion of coefficient entries to integers is legitimate:
\(D_t\) has even integral coefficients, so the recursion leaves every \(N_k\)
integral; all displayed coefficient-system entries are integers.

The saved \(U\) also satisfies the cleared polynomial identity directly.
This existence verification is independent of the discovery algorithm.
Uniqueness is only within the stated degree/symmetry ansatz; neither global
minimal order nor a conceptual explanation of that particular search space
follows from its rank.

## 3. The quadratic connection and the rational lift

For a generic parameter derivative write \(A=a_v,B=b_v,C=c_v\),
\(g=D_v/(2D)\), and

\[
R_v=(A/a-g)z+(B-bg)/(2a),\qquad q_v=R_v/f.
\]

The connection identity is equivalent to

\[
-f_v=-g f+R_v'f-R_vf_z.
\]

Comparing its \(z^2,z,1\) coefficients gives respectively \(-A,-B,-C\).
For the constant coefficient use \(D_v=2bB-4(Ac+aC)\).
This proves the identity for an arbitrary quadratic, before substituting the
Laurent coefficients.

With \(S_k=N_k/D^k\), differentiation gives
\(\theta S_k+tS_kh_t=S_{k+1}\). Therefore

\[
T_{k+1}=\theta T_k+tS_kq_t
\]

really proves \(\theta^k(f^{-1})=S_kf^{-1}+\partial_zT_k\).
Applying the \(x,y\) connection to the two primitive factors yields the
rational three-variable divergence, including the minus signs in
\(B=T_{\mathcal L}-(U/D^2)q_x-(V/D^2)q_y\).

At \(z=0\), the independently derived recursion gives

\[
T_1=\frac{2s}{D},\qquad
T_2=-\frac{3stD_t}{D^2},\qquad
T_3=-\frac{4s(tD_t+t^2D_{tt})}{D^2}
     +\frac{15st^2D_t^2}{2D^3}.
\]

These agree with the construction. The rational boundary identity
\(W_0=P/(pD^3)=\partial_x(E/(pD^2))+\partial_y(E^\top/(pD^2))\)
passes a separate polynomial-ring check. Here \(E^\top(x,y,t)=E(y,x,t)\).
No differentiation of the ambiguously named polynomial \(B_x\) is used.

## 4. The repaired positive-domain bound

Let

\[
a_*=p+p^2/s,\qquad b_*=2\sqrt{p^3/s}.
\]

Writing \(A=\sqrt p\), \(v=\sqrt s\), the inequality
\(x+y\ge2\sqrt{xy}\) gives \(A\ge1+v\). Hence

\[
\frac{A^2}{v}-A
=\frac{A(A-v)}v\ge1+\frac1v>1,
\]

and \(a_*-b_*=(A^2/v-A)^2>1\). It follows that
\(a_*^2-b_*^2\ge a_*\). Since
\(s^2(a_*^2-b_*^2)=p^2r^2\), for any fixed \(t>2\),

\[
D=t^2s^2(a_*^2-b_*^2)-2ts^2a_*+s^2
\ge t(t-2)p^2r^2>0.
\tag{1}
\]

The construction uses \(t>\lambda=17+12\sqrt2>2\), so (1) applies.
The positive minimum \(\Lambda\ge\lambda\) also gives
\(f=sz(1-t\Lambda)<0\) throughout the positive octant. The rational lift has
no pole on that domain.

All constants below may depend on this fixed \(t\), or uniformly on a compact
parameter interval contained in \((\lambda,\infty)\). Such dependence is harmless
for proving the differential identity.

## 5. Boundary \(E\): the corner limit is dominated uniformly

The exact support audit gives
\(\deg_xE\le8,\deg_yE\le7,\deg_{x,y}E\le13\), and \(x\mid E\).
By (1),

\[
\left|\frac E{pD^2}\right|\le C_t\frac{|E|}{p^5r^4}.
\]

At \(x=\varepsilon\le1\), the integral over all \(y>0\) is \(O(\varepsilon)\):
the numerator has that factor and at most \(y^7\), while the denominator
is bounded below by a fixed multiple of \((1+y)^9\).

At \(x=R\ge1\), a monomial \(x^iy^j\) of total degree at most 12 contributes
\(O(R^{-1})\) from \(0<y<1\). From \(1<y<R\), its bound is

\[
C_tR^{i-9}\int_1^R y^{j-5}\,dy.
\]

For \(j<4\), use \(i\le8\); for \(j=4\) there is one logarithm; for
\(j>4\), use \(i+j\le12\). Every case is \(O(\log R/R)\).

Only

\[
E_{13}=12t^3x^6y^5
\bigl[t(x+y)^2-2(x^2+y^2)\bigr]
\]

survives. After \(y=Rv\), for \(0\le v\le1\), (1) gives the uniform bound

\[
R\left|\frac{E_{13}(R,Rv)}{pD^2}\right|
\le C_t\frac{(Rv)^5}{(1+Rv)^5(1+v)^2}\le C_t.
\tag{2}
\]

This explicitly handles the degenerating corner \(v=0\). For every fixed
\(v>0\), the exact highest homogeneous parts give the limiting integrand

\[
\frac{12t[t(1+v)^2-2(1+v^2)]}
     {[t(1+v)^2-4v]^2}.
\]

Its pointwise value at \(v=0\) need not agree with the prelimit value:
that one point has measure zero, and (2) supplies dominated convergence.
The denominator stays positive since \(t>2\).
The antiderivative
\[
\frac{6t(v^2-1)}{t(1+v)^2-4v}
\]
has increment 6 from 0 to 1. The transposed edge contributes another 6.

## 6. \(W_0\) is absolutely integrable

The saved numerator has
\(\deg_xP,\deg_yP\le11\) and total degree at most17. Thus
\[
|W_0|\le C_t|P|/(p^7r^6).
\]
On \(x\ge y\ge1\), a monomial is bounded by
\(x^{i-13}y^{j-7}\).
If \(j<6\), its \(y\)-integral is bounded; if \(j=6\), it is logarithmic;
if \(j>6\), it is \(O(x^{j-6})\). The remaining \(x\)-integral converges
using respectively \(i\le11\) or \(i+j\le17\).
The other large quadrant is symmetric. Where one variable is bounded,
the other has exponent at most \(-2\); near the origin (1) prevents a pole.

Consequently the square exhaustion computes an actual absolutely convergent
integral, not merely a conditionally assigned corner value:
\[
\int_{(0,\infty)^2}W_0\,dx\,dy=12.
\]

## 7. Integrate in \(z\) first, then exhaust \(x,y\)

For fixed positive \(x,y,t\), every \(q_v=O(1/z)\) as \(z\to\infty\).
Parameter differentiation preserves this bound, so \(B(\infty)=0\);
at \(z=0\), \(f=c=-tp^2\ne0\). On each compact rectangle
\(\varepsilon\le x,y\le R\), these estimates and their needed derivatives
are uniform. The \(f^{-1}\) integrands and their \(x,y,t\) derivatives are
\(O(z^{-2})\) at infinity and bounded at zero. This justifies differentiating
their convergent \(z\)-integrals. The integrated \(z\)-flux is precisely
\(-W_0\).

The remaining flux involves the exact positive identity
\[
\int_0^\infty\frac{dz}{xyz\Lambda}
=\int_0^\infty\frac{dz}{p(z+1)(sz+p)}
=\frac{\log(p/s)}{pr}.
\]
For fixed \(t>\lambda\),
\(|1-t\Lambda|^{-1}\le C_t/\Lambda\). Hence the \(x\)-flux is bounded by
\[
C_t\frac{|U|\log(p/s)}{p^5r^5},
\qquad
\log(p/s)=\log(1+1/x)+\log(1+1/y).
\]

At \(x=\varepsilon\), \(x\mid U\) and \(\deg_yU\le8\) give
\(O(\varepsilon\log(1/\varepsilon))\) after integration over all \(y>0\).
At \(x=R\), the segment \(y<1\) contributes \(O(1/R)\).
For \(1<y<R\), use \(\log(p/s)\le1/R+1/y\). A monomial is bounded by
\[
C_tR^{i-10}
\left(\frac1R\int_1^R y^{j-5}dy+\int_1^R y^{j-6}dy\right).
\]
The conditions \(i\le9,j\le8,i+j\le14\) give \(O(\log R/R)\).
The \(y\)-flux is identical. These estimates are uniform in the opposite
cutoff, so the positive-octant limit does not exchange uncontrolled boundary
limits simultaneously.

Finally, all needed \(t\)-derivatives of
\((1-t\Lambda)^{-1}d\log x\,d\log y\,d\log z\)
are dominated, locally uniformly in \(t\), by
\(C_t\Lambda^{-1}d\log x\,d\log y\,d\log z\).
This integrable measure is the already established zeroth Beukers period;
its finiteness suffices here, not its numerical evaluation.
Differentiation under the full real integral is therefore valid. Combining
the vanishing \(x,y\)-fluxes with the negative \(z\)-flux proves
\(\mathcal L_tH=-12\).

## 8. Reciprocal normalization and generating functions

Directly from the kernels,
\[
H(t)=-t^{-1}F(1/t),\qquad
F(u)=\int\frac{d\log x\,d\log y\,d\log z}{\Lambda-u}.
\]
The operator identity
\[
\mathcal L_u[u^{-1}h(1/u)]
=-u(\mathcal L_th)(1/u)
\]
is correct. The additional minus sign in \(F(u)=-u^{-1}H(1/u)\) then gives
\(\mathcal L_uF(u)=-12u\), not \(+12u\) or a constant.

Since \(\Lambda\ge\lambda\) and the zeroth reciprocal period is finite,
the geometric expansion for \(F\) is absolutely dominated for \(|u|<\lambda\).
It defines the analytic germ \(\sum I_nu^n\). The identity initially obtained
on \(0<u<1/\lambda\) therefore holds as a germ at zero and throughout that
disk. Its first coefficient is \(I_1-5I_0=-12\); higher coefficients give
the homogeneous Apéry recurrence. This checks the normalization against the
exceptional boundary index without using those coefficients to derive 12.

## 9. Reproduction and remaining limitations

Fresh runs passed:

- The discovery program: 1,180 equations, rank383, 30-dimensional primitive
  freedom, and the eight displayed operator coefficients.
- The original exact checker, including both saved polynomial certificates.
- [check_geometric_independent.py](../computations/check_geometric_independent.py):
  a generic quadratic connection, the independently derived boundary recursion,
  reciprocal operator conjugation, polynomial-ring certificate checks, full
  support bounds, highest homogeneous pieces, and the edge integral.

The independent checker does not turn analytic estimates into consequences
of a computer assertion. Sections 1 and 4–8 above supply those arguments.
All checks concern the repaired construction. The initial false discriminant
claim is neither retained nor used.

The remaining problem is exposition: compressing the 231-term and 188-term
primitives, or replacing them by a short natural geometric calculation.
The present certificates close the exact-form and boundary gap, while the
Newton-balance derivation remains a more accessible explanation of the
recurrence for the team's main talk.
