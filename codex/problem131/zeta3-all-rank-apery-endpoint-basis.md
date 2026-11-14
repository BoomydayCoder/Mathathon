# Every original zeta(3) rank has a basis of rank-one polynomials

Status: accepted in [the independent audit](zeta3-all-rank-independent-review.md). The inputs are the
accepted original Schur/Apéry identity, strict monotonicity of its q_h, and
the fixed Heun quotient from `zeta3-rank-q-heun-reduction.md`. No historical
novelty claim is made. This note does not replace any already accepted proof.

Let 1<=q<=h, M=2h-q+1, n=h-1, a=h-q=M-h-1 and l=a+1=M-h. Let W_(h,M)
be the q-dimensional space of polynomials of degree<=n that are orthogonal
to every polynomial of degree<=a-1 for the original positive ell_M measure
in z=-t. All these orthogonality equations are rational and independent of
the affine zeta variable X: their maximum degree is n+a-1=M-3, below its
first nonzero coefficient degree M-2.

Use the accepted rank-one polynomials R_j(0)=1 and quantities

\[
 L_j=[z^{j-1}]R_j=\frac{j^2C_j}{12(2j)!},\quad
 q_j=D_j/C_j,\quad C_j=5A_j-A_{j-1},\quad D_j=5B_j-B_{j-1}.
\]

The original rank-one nullvector property is needed in its full form:

\[
 \ell_{2j,q_j}(R_j P)=0\qquad(\deg P\le j-1).      \tag{1}
\]

This includes the last Gram row because q_j is its determinant root.

## 1. A concrete q-element basis

For every integer j=l,...,h define

\[
 F_j^{(M)}(z)=
 \begin{cases}
 R_j(z),&2j\ge M,\\
 \displaystyle\prod_{r=2j+1}^M(r^2-z)R_j(z),&2j<M.
 \end{cases}                                                    \tag{2}
\]

Then these q polynomials form a basis of W_(h,M).

First check membership with all degree bounds. In the upper branch,
deg F_j=j-1<=h-1. For a test P of degree<=a-1, changing ell_M to ell_(2j)
multiplies P by product_(r=M+1)^(2j)(r²-z). Its degree is at most

\[
 a-1+2j-M=2j-h-2\le j-2,
\]

so rank-one orthogonality applies. In the lower branch the degree is
M-j-1<=h-1 because j>=M-h=l. Cancelling its extra denominator factors
leaves ell_(2j)(R_jP), and a-1<=j-2. Thus membership holds in both branches,
including the equality 2j=M.

For q>=2, precisely the two extreme polynomials have degree n:

\[
 F_h=R_h,\qquad
 F_l=\prod_{r=2l+1}^M(r^2-z)R_l,
\]

with leading coefficients L_h and epsilon L_l, where epsilon=(-1)^(q-1).
All interior F_j have degree at most n-1.

### Independence from strict monotonicity

Fix M and induct by increasing h, so q increases by two. When M=2m,
the initial q=1 space at h=m is spanned by R_m. When M=2m-1, start from
the empty q=0 space at h=m-1 (the whole degree<=h-1 space is then required
to be orthogonal to itself). The first nonempty step has q=2.

At the induction step the interior polynomials j=l+1,...,h-1 are exactly
the preceding basis of W_(h-1,M). That smaller space consists of degree<=h-2
polynomials orthogonal to all degrees<=a. Suppose
alpha F_h+beta F_l belongs to this smaller space. Its leading coefficient
must vanish:

\[
 \alpha L_h+\beta\epsilon L_l=0.                 \tag{3}
\]

The affine moment computation in the next section gives, with
s=(-1)^(M-1)=(-1)^q,

\[
 \ell_{M,X}(F_hz^a)=sL_h(X-q_h),\qquad
 \ell_{M,X}(F_lz^a)=s\epsilon L_l(X-q_l).
\]

Orthogonality to z^a at X=zeta(3), together with (3), therefore forces
s alpha L_h(q_l-q_h)=0. Since l<h at every nontrivial induction step,
strict monotonicity gives q_l>q_h. Thus alpha=beta=0. This proves that the
two new extremes are independent modulo the old space, completing the
basis proof. The index l is always at least one in our range q<=h.

## 2. Only the two extremes pair with the first missing test degree

For every j=l,...,h,

\[
 \ell_{M,q_j}(F_jz^a)=0.                          \tag{4}
\]

In the upper branch, conversion to ell_(2j) gives a test of degree

\[
 a+2j-M=2j-h-1\le j-1;
\]

in the lower branch cancellation leaves the test z^a with a<=j-1.
Thus the FULL nullvector property (1), rather than only ordinary lower-degree
orthogonality, proves (4).

The product F_jz^a has degree at most n+a=M-2. The affine coefficient
functional vanishes below degree M-2 and equals s=(-1)^(M-1) on the monic
polynomial of that degree. Consequently

\[
 \boxed{\begin{aligned}
 \ell_{M,X}(F_hz^a)&=sL_h(X-q_h),\\
 \ell_{M,X}(F_lz^a)&=s\epsilon L_l(X-q_l),\\
 \ell_{M,X}(F_jz^a)&=0\quad(l<j<h).
 \end{aligned}}                                                   \tag{5}
\]

The first two lines are for q>=2; when q=1 there is only one basis vector
and the already accepted rank-one formula applies. Equation (5) is an exact
affine identity for all X. No limiting or numerical argument enters it.

## 3. Every rank-q root is a fractional combination of two Apéry endpoints

Assume q>=2. Let eta be one of the distinct eigenvalues of the accepted
integer q-by-q Heun quotient, and let its eigenpolynomial be expanded in
(2):

\[
 F_\eta=\sum_{j=l}^h\gamma_j(\eta)F_j.
\]

The coefficients can be chosen in Q(eta), because the basis and the Heun
operator are rational. Each eigenpolynomial has degree EXACTLY n. Indeed,
the quotient raising coefficient is 2(k-n)(k-a+1), which is nonzero for
a<=k<n; an eigenpolynomial of any such smaller degree would be raised.
Hence

\[
 \gamma_hL_h+\gamma_l\epsilon L_l\ne0.            \tag{6}
\]

Let r_eta be the associated affine Gram root. Simultaneous diagonalization
with the Heun quotient makes F_eta a nullvector of the Gram form at X=r_eta
on W_(h,M). It is orthogonal to the lower space for every X as well. It is
therefore a nullvector on the entire degree<=n space at this X, including
against z^a. Substituting (5) yields

\[
 \boxed{r_\eta=
 \frac{\gamma_h L_hq_h+\gamma_l\epsilon L_lq_l}
      {\gamma_h L_h+\gamma_l\epsilon L_l}.}       \tag{7}
\]

In terms of the classical Apéry pairs, with s_j=j²/[12(2j)!], this is

\[
 \boxed{r_\eta=
 \frac{\gamma_hs_hD_h+\gamma_l\epsilon s_lD_l}
      {\gamma_hs_hC_h+\gamma_l\epsilon s_lC_l},
       \qquad l=h-q+1.}                          \tag{8}
\]

The denominator is nonzero by (6). Individual extreme coefficients are
NOT asserted to be nonzero; if one vanishes, the formula is still valid.
Repeated roots for different Heun eigenvalues are likewise not excluded
at q>=3.

This result gives a constructive route: build the finite rank-one basis
using the accepted Heun or RII recurrence, represent the fixed Heun operator
in it, and take the two extreme coordinates of an eigenvector. Every
original rank-q root then uses just the two Apéry endpoints h and h-q+1.
The rank-two closed formula in the earlier note is the case where the two
coordinates can be written immediately from evaluations at z=0 and z=1.
The coefficients gamma retain the remaining finite spectral and basis
calculation; (8) is not an assertion of a new two-term recurrence in h or
a uniform arithmetic size bound for those coefficients. The statement
supplies no new denominator bounds for arbitrary algebraic rank-q eigenfields.
