# The refined rank lemma from two elementary bounds

Status: independently checked proof of the parent's simplification. It
retains the original source estimate without sorted-minor optimization.
The optional residue-field-rank extension below is also proved. No symmetry
or invertibility hypothesis on the whole matrix is needed.

## Statement

Let \(v_p^G\) be the Gauss valuation on \(\mathbb Q_p(X)\): for a polynomial
it is the minimum valuation of its coefficients, and for a quotient it is
the difference of those valuations. Let \(A,L\) be \(h\)-by-\(h\) matrices
over this field. Suppose

\[
G=A+p^{-1}L,\qquad
w_i\in\tfrac12\mathbb Z_{\le0},\qquad
v_p^G(A_{ij})\ge w_i+w_j,\qquad v_p^G(L_{ij})\ge0.
\]

If \(\operatorname{rank}_{\mathbb Q_p(X)}L\le r\), and
\(z=\#\{i:w_i=0\}\), then

\[
\boxed{v_p^G(\det G)\ge2\sum_{i=1}^h w_i-\min(r,z).}       \tag{1}
\]

The source has \(L\) independent of \(X\), so its ordinary
\(\mathbb Q_p\)-rank provides such an \(r\).
In fact (1) remains valid with the potentially smaller rank of
\(\overline L\) over \(\mathbb F_p(X)\); the additional argument for this
extension is in Section 3.

## 1. Two estimates prove the lemma

Work temporarily in \(\mathbb Q_p(\sqrt p)(X)\), with the extension
valuation normalized by \(v_p(p)=1\). Set
\(D=\operatorname{diag}(p^{-w_i})\). Then \(DAD\) is integral in Gauss
valuation, while every entry of \(D(p^{-1}L)D\) has valuation at least
\(-1\), because \(w_i\le0\). Expand \(\det(DGD)\) by multilinearity in
its columns. A term using \(k\) columns from the second matrix has valuation
at least \(-k\); if \(k>r\), it vanishes by exact column dependence.
Therefore

\[
v_p^G(\det G)\ge2\sum_iw_i-r.                            \tag{2}
\]

Independently replace every zero weight by \(-1/2\), leaving negative
weights unchanged, and call the result \(w_i'\). The half-integer hypothesis
implies \(w_i'\le-1/2\) for every \(i\). Hence

\[
v_p^G(A_{ij})\ge w_i+w_j\ge w_i'+w_j',\qquad
v_p^G(p^{-1}L_{ij})\ge-1\ge w_i'+w_j'.
\]

Scaling rows and columns by \(p^{-w_i'}\) thus makes the whole matrix
integral, giving

\[
v_p^G(\det G)\ge2\sum_iw_i'
                   =2\sum_iw_i-z.                     \tag{3}
\]

The larger of (2) and (3) is (1). Extending the coefficient field has not
changed the valuation of the original determinant. If its entries are
polynomials in \(X\), the result is exactly a coefficientwise assertion.

## 2. Interpretation and application

The correction costs at most one power of \(p\) per independent column:
that is (2). Separately, every negative half-integer weight already supplies
enough row and column scaling to absorb the correction entrywise. Only the
zero-weight rows need an additional half unit on each side: that is (3).
Use whichever estimate is cheaper.

All source outer weights are nonpositive half-integers, and the correction
has the required actual rank bound. Thus the source's refined estimate
follows without optimizing complementary minors. Counting zero weights
remains necessary when a numerical formula for \(z\) is wanted. The older
plain-rank normalization remains valid, but it is no longer needed merely
to shorten the proof of this lemma.

The half-integer grid is used in (3), to ensure every negative weight is
at most \(-1/2\). The rank estimate (2) needs only nonpositive weights and
a valued extension containing the required scaling factors.

## 3. A stronger version using rank modulo \(p\)

Let \(\mathcal O\) be the valuation ring of the Gauss valuation on
\(\mathbb Q_p(X)\). It has uniformizer \(p\) and residue field
\(\mathbb F_p(X)\). Suppose

\[
\operatorname{rank}_{\mathbb F_p(X)}\overline L=\rho.
\]

For \(0<\rho<h\), permute rows and columns to place a nonzero
\(\rho\)-minor of \(\overline L\) in the upper left corner, and write

\[
L=\begin{pmatrix}U&V\\W&Z\end{pmatrix}.
\]

Here \(U\) is invertible over \(\mathcal O\): its determinant has Gauss
valuation zero. Define

\[
L_0=\begin{pmatrix}U&V\\W&WU^{-1}V\end{pmatrix}.
\]

This matrix is integral and has exact rank \(\rho\). The Schur complement
\(Z-WU^{-1}V\) reduces to zero, by the assumed residue-field rank.
Consequently

\[
L=L_0+pL_1,\qquad L_1\in M_h(\mathcal O).
\]

Undo the permutations to obtain the same decomposition in the original
coordinates. The cases \(\rho=0\) and \(\rho=h\) use \(L_0=0\) and \(L_0=L\),
respectively. Now

\[
G=(A+L_1)+p^{-1}L_0.
\]

Since all \(w_i+w_j\le0\), adding the integral matrix \(L_1\) preserves
every bound on \(A\). Apply Sections 1–2 with actual rank \(\rho\).
This proves

\[
\boxed{v_p^G(\det G)\ge2\sum_iw_i-\min(\rho,z).}
\]

For polynomial \(L(X)\), the pivot inverse may introduce rational functions
of \(X\); this is harmless in the Gauss valuation ring. A nonzero pivot
minor modulo \(p\) has Gauss valuation zero even if it is not a constant
polynomial. The original determinant is still polynomial and the final
valuation statement is still coefficientwise.

Thus residue-field rank does suffice, but not by pretending the original
large minors vanish exactly. One first constructs a lift of small actual
rank and absorbs its integral remainder into \(A\). This distinction is
essential when comparing the short multilinearity proof with modular
kernel arguments elsewhere in the notes.

## 4. In the source outer matrix, the residue-field rank is already exact

The stronger version above does not by itself give another improvement for
the source outer correction. Its original rank bound is attained modulo
\(p\), as follows from a general monic-division observation.

Let \(D,W\) be monic over a field, with degrees \(h,w\), respectively.
Let \(\nu\) be a linear functional on polynomials satisfying

\[
\nu(t^e)=0\ (0\le e<e_0),\qquad \nu(t^{e_0})=c\ne0.
\]

Form the \(h\)-by-\(h\) matrix

\[
M_{ij}=\nu\!\left(\operatorname{quotient}
                   \frac{Wt^{i+j}}D\right),\qquad 0\le i,j<h.
\]

Put \(r=\max(0,h+w-e_0-1)\), and assume \(r\le h\).
If \(r=0\), every quotient has degree below \(e_0\), so \(M=0\).
If \(r>0\), the first \(h-r\) rows and columns vanish, because the largest
possible quotient degree in any such row or column is below \(e_0\).
In the remaining block, write \(i=h-r+u\), \(j=h-r+v\).
For \(u+v<r-1\) its entry is zero; for \(u+v=r-1\), the quotient
is monic of degree \(e_0\), so the entry is \(c\). Thus this block is
anti-triangular, with determinant

\[
(-1)^{r(r-1)/2}c^r\ne0.
\]

Consequently \(\operatorname{rank}M=r\), independently of every lower
coefficient of \(W,D\) and every higher moment of \(\nu\).

For the source ζ(5) functional over \(\mathbb F_p\), \(p\ge7\),

\[
\nu(t^e)=
\begin{cases}
\displaystyle
\frac{(-1)^{e+1}(2e+3)(2e+4)(2e+5)}{24},
 &p-1\mid2e+2,\\[4pt]
0,&p-1\nmid2e+2.
\end{cases}
\]

The first three positive multiples of \(p-1\) are cancelled by the three
displayed factors. Hence the first nonzero moment is at

\[
e_0=2p-3,\qquad \nu(t^{e_0})=-\tfrac14\ne0.
\]

Take \(W=D_N^5\), \(w=5N\), \(D=D_K/D_N\), \(h=K-N\). This gives

\[
\operatorname{rank}_{\mathbb F_p}\overline L
 =\max(0,K+4N-2p+2)=r_p.
\]

The source condition \(5N\le2p-2\) ensures \(r_p\le h\).
Combined with the independently proved upper bound on actual rank, both
actual and residue-field ranks are exactly \(r_p\). Unimodular CRT basis
changes preserve these statements.

For completeness, the same calculation works for odd \(s\ge3\) and
\(p>s\), using the moment functional defined in
outer-correction-Hermite-ideal.md:

\[
\nu_s(t^e)=
\begin{cases}
\displaystyle\frac{(-1)^{e+1}\prod_{j=3}^{s}(2e+j)}{(s-1)!},
 &p-1\mid2e+2,\\[4pt]
0,&p-1\nmid2e+2.
\end{cases}
\]

The first \(s-2\) positive multiples are cancelled. The first surviving
index and its value are

\[
e_0=\frac{(s-1)(p-1)}2-1,\qquad
\nu_s(t^{e_0})=-\frac1{s-1}.
\]

Here \(e_0\) is odd and the product of the factors modulo \(p\) is
\(-(s-2)!\), which verifies the sign. Provided
\(5N\le(s-1)(p-1)/2\), the same monic-division lemma gives exact rank

\[
r_s=\max\!\left(0,K+4N-\frac{(s-1)(p-1)}2\right).
\]

The additional ζ(5) Hermite-ideal gain therefore comes from the position of
the correction kernel relative to the row weights, not from an unexploited
reduction of its global rank.
