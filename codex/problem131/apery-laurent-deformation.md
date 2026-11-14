# A positive Laurent deformation: why the Apéry point loses two critical points

Status: direct algebraic classification and analytic consequences proved below;
exact symbolic checks in `check_apery_laurent_deformation.py`. This note extends
the reciprocal Laurent model, with no historical novelty claim. The special
number of critical values does not by itself prove the order of a recurrence.

For \(q>0\), let
\[
\Lambda_q(x,y,z)=
\frac{(1+x)(1+y)(1+z)((1+x)(1+y)+qxyz)}{xyz}.
\]
Write \(A=(1+x)(1+y)\), \(B=xy\), and define
\[
A_n(q)=\operatorname{CT}\Lambda_q^n
=\sum_{k=0}^n\binom nk^2\binom{n+k}k^2q^{n-k}.              \tag{1}
\]
The coefficient-extraction proof is the same binomial-theorem calculation as
for \(q=1\), with the unused term carrying the factor \(q^{n-k}\).

## 1. The reciprocal construction persists throughout the family

The odds substitution \(a=x/(1+x),b=y/(1+y),c=z/(1+z)\) gives
\[
h_q(a,b,c)=\frac{a(1-a)b(1-b)c(1-c)}{1-c+qabc}
=\Lambda_q^{-1},
\]
\[
\frac{da\,db\,dc}{1-c+qabc}
=\Lambda_q^{-1}\frac{dx\,dy\,dz}{xyz}.
\]
Thus the positive integral
\[
I_n(q)=\int_{(0,1)^3}\frac{h_q^n}{1-c+qabc}\,da\,db\,dc
=\int_{(0,\infty)^3}\Lambda_q^{-n-1}\frac{dx\,dy\,dz}{xyz}
\tag{2}
\]
is paired with the integer-coefficient polynomial \(A_n(q)\). For every fixed
\(q>0\), \(I_0(q)\) is finite. Indeed its cube denominator is at least
\(\min(1,q)(1-c+abc)\), whose reciprocal has finite integral.

Put \(d=\sqrt q\) and \(v=\sqrt{1+d}>1\). The same two AM–GM steps give
\[
\min\Lambda_q=(1+v)^4=:\lambda_q,\qquad
x=y=\frac1v,\quad z=\frac1{v-1}.                            \tag{3}
\]
For clarity, the inequality is
\[
\Lambda_q=qAz+\frac{A^2}{Bz}+qA+\frac{A^2}{B}
\ge A\left(d+\sqrt{A/B}\right)^2.
\]
With \(r=\sqrt B\) and \(\sqrt A\ge1+r\),
\[
\sqrt{\Lambda_q}\ge d(1+r)+(1+r)^2/r
=(1+d)r+(2+d)+1/r
\ge2+d+2\sqrt{1+d}=(1+v)^2.
\]
Equality conditions give exactly (3). In particular,
\[
A_n(q)\le\lambda_q^n,\qquad
0<I_n(q)\le I_0(q)\lambda_q^{-n}.
\]
Reciprocal exponential scales are therefore a property of the whole positive
family, not a uniqueness criterion for the Apéry point.

## 2. All nonzero torus critical points

Treat \(x,y,z\) as nonzero complex variables and exclude \(\Lambda_q=0\).
As in the undeformed case,
\[
x\partial_x\Lambda_q-y\partial_y\Lambda_q=(x-y)F_A,
\qquad BF_B=-\frac{A^2(1+z)}{Bz}\ne0.
\]
The critical equations force \(x=y=u\). The \(z\)-equation becomes
\[
z=\pm\frac{1+u}{du}.
\]
On the branch labelled by \(\epsilon=\pm1\), the remaining function is
\[
\Lambda_q\left(u,u,\epsilon\frac{1+u}{du}\right)
=\left((1+\epsilon d)u+(2+\epsilon d)+\frac1u\right)^2.
\tag{4}
\]
For nonzero values, its derivative vanishes precisely when
\[
u^2=\frac1{1+\epsilon d}.                                  \tag{5}
\]
If \(1+\epsilon d=0\), there is no finite solution of this equation. Apart
from that exception, the two values on the branch are
\[
\left(\sqrt{1+\epsilon d}+1\right)^4,
\qquad
\left(\sqrt{1+\epsilon d}-1\right)^4.                       \tag{6}
\]
Changing the square-root sign only exchanges these two values.

The product of each pair is \(q^2\), and its sum is
\(2(q+8)+16\epsilon\sqrt q\). Consequently the generic critical-value
polynomial is
\[
\boxed{\quad
\left[V^2-2(q+8)V+q^2\right]^2-256qV^2.
\quad}                                                     \tag{7}
\]
This formula records values, with multiplicity, when all four critical points
remain in the torus. Exceptional parameters must be checked before interpreting
its factors as actual critical values.

## 3. Exact exceptional-parameter classification for q>0

For \(q>0\), the two branch quadratics cannot share a root: subtracting them
would force \(V=0\), whereas their constant term is \(q^2\ne0\).
The discriminant within branch \(\epsilon\) is
\[
64(1+\epsilon d)(d+2\epsilon)^2.
\]
The plus branch always has two distinct values. The minus branch has a repeated
value exactly at \(d=1\) or \(d=2\). It follows that:

- If \(q>0\) and \(q\ne1,4\), there are four critical points with four
  distinct nonzero values.
- At \(q=4\), there are still four critical points, but only three distinct
  values. The minus branch consists of \(u=\pm i\), both with value \(-4\).
  The factorization is
  \((V^2-56V+16)(V+4)^2\).
- At \(q=1\), the minus branch has no finite critical point. Its two points
  escape to the boundary as \(u\to\infty,z\to-1\). The two remaining
  critical values are \(17\pm12\sqrt2\). Formula (7) specializes to
  \((V^2-34V+1)(V-1)^2\), but **the factor \((V-1)^2\) does not represent
  torus critical points at this parameter**.

Thus \(q=1\) is the unique positive parameter with exactly two nonzero torus
critical values. This is a concrete geometric simplification of the Apéry point.
It is not sufficient, on its own, to infer a third-order Picard–Fuchs equation,
a three-term recurrence, or an irrationality theorem.

The limiting parameter \(q=0\) is outside this classification. There
\(A_n(0)=\binom{2n}n^2\), but \(\partial_z\Lambda_0\) is nonzero whenever
\(\Lambda_0\ne0\), so there are no nonzero torus critical points. The
positive minimum becomes an infimum at \(z\to\infty\), and the real
integral (2) diverges. Hence that boundary sequence does not retain the same
finite auxiliary-integral construction.

## 4. A uniform Hessian and the product of the two asymptotic constants

The covariance/Hessian at the unique positive critical point has a compact
expression in \(v=\sqrt{1+\sqrt q}\):
\[
\Sigma_q=\begin{pmatrix}
\dfrac{2v}{(v+1)^2}&\dfrac{v-1}{(v+1)^2}&\dfrac{v-1}{v(v+1)}\\
\dfrac{v-1}{(v+1)^2}&\dfrac{2v}{(v+1)^2}&\dfrac{v-1}{v(v+1)}\\
\dfrac{v-1}{v(v+1)}&\dfrac{v-1}{v(v+1)}&\dfrac{2(v-1)}{v^2}
\end{pmatrix},
\]
\[
\det\Sigma_q=\frac{4(v-1)}{v(v+1)^3}>0.
\tag{8}
\]
All coefficients of the Laurent polynomial are positive on its support for
\(q>0\); zero and the signed coordinate steps remain present. The Fourier and
real Laplace proofs in `apery-laurent-reciprocal.md` therefore apply unchanged,
with \(q\) fixed:
\[
A_n(q)\sim\frac{\lambda_q^n}{(2\pi n)^{3/2}\sqrt{\det\Sigma_q}},
\qquad
I_n(q)\sim\frac{(2\pi)^{3/2}\lambda_q^{-n-1}}
{n^{3/2}\sqrt{\det\Sigma_q}}.
\]
In particular,
\[
\boxed{\quad
A_n(q)I_n(q)\sim
\frac{\sqrt{1+\sqrt q}}{4\sqrt q}\,n^{-3}.
\quad}                                                     \tag{9}
\]
This recovers \(1/(2\sqrt2)\) at \(q=1\). It also shows why the small
parameter limit is not uniform: the product constant diverges as \(q\downarrow0\),
just as the real integral loses integrability at the boundary.

For positive integer \(q\), integrality and Lucas congruences persist because
the exponent support remains in the same cube. Their persistence, like that of
the reciprocal rates, does not solve the arithmetic identification of the
companion constant. The parent's separate deformation-arithmetic investigation
identifies the new period values and explains why \(q=1\) is also special there.
