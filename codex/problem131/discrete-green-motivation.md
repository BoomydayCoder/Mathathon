# Apéry's recurrence as an electrical network and a sharp inequality

Status: complete elementary derivation below, newly developed in this session.
The ingredients are standard ground-state transformation and resistance in series;
no historical novelty is claimed. Identification with zeta(3) uses the separately
proved bounded rational-function representation in discovery-guide.md. The network
alone does not identify its resistance with zeta(3).

## 1. Why a growing solution controls a shrinking error

Write
\[
b_n=(2n+1)(17n^2+17n+5),\qquad
(n+1)^3u_{n+1}=b_nu_n-n^3u_{n-1}.
\]
Let \(A_0=1,A_1=5\), \(B_0=0,B_1=6\). The binomial formula proves
\(A_n\) is a positive integer; its recurrence gives \(A_{n+1}\ge5A_n\).

For a general solution put \(u_n=A_nv_n\). Subtract the equation for \(A\)
multiplied by \(v_n\). The zeroth-order term disappears, leaving
\[
c_n(v_{n+1}-v_n)=c_{n-1}(v_n-v_{n-1}),\qquad
c_n=(n+1)^3A_nA_{n+1}.
\tag{1}
\]
This is Kirchhoff's law on the chain \(0,1,2,\ldots\), with edge conductance
\(c_n\), or resistance \(r_n=1/c_n\). Current has the same value across every edge.
For \(u=B\), the first voltage drop is \(B_1/A_1=6/5\) and \(c_0=5\),
so the current is exactly 6. Consequently
\[
\frac{B_n}{A_n}=6\sum_{j=0}^{n-1}r_j,\qquad
L:=\lim_n\frac{B_n}{A_n}=6\sum_{j=0}^\infty r_j.
\tag{2}
\]
Convergence follows from \(A_n\ge5^n\).

The previously proved weak estimate \(A_n\zeta(3)-B_n=O(1)\) identifies
\(L=\zeta(3)\). Thus
\[
\boxed{\quad \zeta(3)=6\sum_{n\ge0}
\frac{1}{(n+1)^3A_nA_{n+1}}.\quad}
\tag{3}
\]
This familiar reduction-of-order identity now has an immediate motivation:
the rational approximations measure successively longer pieces of a chain.
Their error is the resistance still omitted:
\[
A_n\zeta(3)-B_n=6A_n\sum_{j=n}^\infty r_j>0.
\tag{4}
\]
Rapid growth of the conductances explains smallness; positive resistance explains
nonvanishing. Neither conclusion needs a delicate signed-summand estimate.

## 2. A ground-state identity with no spectral theory required

For a real sequence \(f\) with finite support in the nonnegative integers, set
\[
\mathcal Q(f)=\sum_{n\ge0}b_nf_n^2
-2\sum_{n\ge0}(n+1)^3f_nf_{n+1}.
\]
Then
\[
\boxed{\quad
\mathcal Q(f)=\sum_{n\ge0}c_n
\left(\frac{f_n}{A_n}-\frac{f_{n+1}}{A_{n+1}}\right)^2.
\quad}
\tag{5}
\]
To prove it, expand the squares. The cross terms are correct because
\(c_n/(A_nA_{n+1})=(n+1)^3\). For \(n\ge1\), the coefficient of \(f_n^2\)
is
\[
\frac{n^3A_{n-1}+(n+1)^3A_{n+1}}{A_n}=b_n.
\]
At \(n=0\) it is \(A_1/A_0=5=b_0\). This checks the boundary as well.

The operation \(f\mapsto f/A\) is called a ground-state transformation.
It turns the unfamiliar cubic recurrence into the ordinary energy of a resistor
chain. Formula (5), rather than the name, is all that is needed below.

## 3. Zeta(3) is a sharp constant in an elementary discrete inequality

For every finitely supported real sequence,
\[
\boxed{\quad
\mathcal Q(f)\ge\frac{6}{\zeta(3)}f_0^2,
\quad}
\tag{6}
\]
and the constant \(6/\zeta(3)\) cannot be increased.

Indeed put \(v_n=f_n/A_n\). Telescoping and weighted Cauchy--Schwarz give
\[
f_0^2=\left(\sum_{n\ge0}(v_n-v_{n+1})\right)^2
\le\left(\sum_{n\ge0}r_n\right)
\left(\sum_{n\ge0}c_n(v_n-v_{n+1})^2\right).
\]
Use (3) and (5). To check sharpness rather than just an upper estimate, for
\(N\ge0\) choose
\[
v_n=\frac{\sum_{j=n}^N r_j}{\sum_{j=0}^N r_j}\quad(0\le n\le N),
\qquad v_n=0\quad(n>N),\qquad f_n=A_nv_n.
\]
Then \(f_0=1\) and
\[
\mathcal Q(f)=\frac1{\sum_{j=0}^N r_j}
=\frac{6A_{N+1}}{B_{N+1}}\longrightarrow\frac6{\zeta(3)}.
\tag{7}
\]
Thus the best constant is approached by explicitly rational finite sequences.
This supplies a variational characterization:
\[
\frac{\zeta(3)}6=
\sup_{f\ne0,\,\mathrm{finite\ support}}\frac{f_0^2}{\mathcal Q(f)}.
\tag{8}
\]
No infinite-dimensional operator theorem is used.

## 4. Why the continued fraction has positive denominators

Let \(J_N\) be the tridiagonal matrix indexed by \(0,\ldots,N\), with diagonal
\(b_n\) and neighboring entries \(-(n+1)^3\). Identity (5) proves it is
positive definite: if its energy were zero, all voltage differences would vanish,
and the zero voltage past \(N\) forces every voltage to vanish.

Minimizing its energy subject to \(f_0=1\), or completing squares in a positive
definite matrix, gives
\[
(J_N^{-1})_{00}=\sum_{j=0}^N r_j
=\frac{B_{N+1}}{6A_{N+1}}.
\tag{9}
\]
Successive Schur complements express the same entry as
\[
\frac1{b_0-\displaystyle\frac{1^6}{b_1-
\displaystyle\frac{2^6}{b_2-\cdots-\displaystyle\frac{N^6}{b_N}}}}.
\tag{10}
\]
Every pivot is positive because each remaining Schur complement is positive
definite. This fills the earlier tail-positivity gap flagged in recurrence.md,
without guessing a separate numerical lower bound for each tail denominator.
Alternatively ordinary determinant expansion gives
\(\det J_N=((N+1)!)^3A_{N+1}>0\), with positive-definiteness additionally
following directly from (5).

## 5. What this does and does not simplify

The recurrence has three distinct roles which this viewpoint separates:

- The positive integral sequence \(A_n\) builds the conductances.
- The initial companion value \(B_1=6\) fixes the current.
- The weak original rational-function estimate identifies the total resistance.

After those steps, positivity, convergence, error bounds, continued fractions,
and the sharp inequality are consequences of elementary energy calculations.
Arithmetic denominator control remains a separate essential input to the
irrationality proof. This network does not explain by itself why Apéry's particular
\(b_n\) has integral solutions, nor why a higher odd zeta should arise from a
similar chain. Those are the construction problems, not consequences of a
formal electrical analogy.

A concrete search principle follows: a positive rational tridiagonal matrix family
with an explicitly integral ground-state solution may provide good approximations
to its limiting boundary Green function. One must independently identify that
constant and prove economical denominators. Positivity alone is not an
irrationality criterion.

Reproduce finite identities with computations/check_discrete_green.py.
