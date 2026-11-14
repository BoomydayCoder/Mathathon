# An explicit Hermite ideal for the outer ζ(5) correction

Status: proved finite-field construction, with the same local-arithmetic
dependencies as `outer-correction-adversarial-review.md`. This replaces the
selected-row/target-column count in `outer-polynomial-correction-kernel.md`
by an explicit polynomial space. It gives the **same**, not an additional,
gain \(27/3200\). No historical novelty is asserted. The research lead and
its attribution remain as recorded in that earlier note.

## 1. The polynomial that explains 49/40

Retain \(K\in40\mathbb Z_{>0}\), \(N=3K/40\), \(h=K-N\), the source's outer
prime hypotheses, \(2p<K<3p\), and

\[
d=K+3N-3p+1>0,\qquad m=(p-1)/2.
\]

In particular \(p\ge7\), \(2N<p\), \(p^2>2K\), and \(5N\le2p-2\).
All polynomials in the next three sections are over \(\mathbb F_p\).
Put \(D=D_K/D_N\), \(W=D_N^5\), and define

\[
\boxed{F(t)=t^2\prod_{a=1}^{N}(t+a^2)^2
                  \prod_{a=N+1}^{m}(t+a^2)^6.}                 \tag{1}
\]

Its degree is

\[
\deg F=2+2N+6(m-N)=3p-1-4N,
\qquad h-\deg F=d.                                        \tag{2}
\]

Consequently

\[
\mathcal I=F\mathbb F_p[t]_{<d}
 \subset \mathbb F_p[t]_{<h}                              \tag{3}
\]

has dimension exactly \(d\), with explicit independent basis

\[
F,tF,\ldots,t^{d-1}F.
\]

The coefficient \(49/40=1+3(3/40)\) now comes directly from the degree
budget \(K+3N-3p+1\). No computation of a correction-matrix rank is needed.

## 2. These vectors lie in the correct CRT row space

Let \(\ell_a\) count \(1\le j\le K\) with \(j\equiv\pm a\pmod p\).
For \(a>0\), \(\ell_a\in\{4,5,6\}\); at zero, \(\ell_0=2\).
At a removed node \(1\le a\le N\), the multiplicity in \(D\) is
\(\ell_a-1\); elsewhere it is \(\ell_a\).

There is no removed class of size four in the positive-\(d\) range.
Indeed, writing \(v=K-2p\), the assumption \(d>0\) gives

\[
v-N> K/3-3N-2/3=13K/120-2/3>0.
\]

Every \(a\le N<v\) therefore has \(\ell_a\ge5\).

The reduced source CRT rows in class \(a\) have the form

\[
g_{a,i}=P_a(t)(t+a^2)^i,\qquad P_a=D/(t+a^2)^{\ell_a-\delta_a},
\]

where \(\delta_a=1\) exactly in removed classes. In the local quotient
at \(a\), \(P_a\) is a unit. Thus the span of its rows of local orders
at least two is exactly the ideal of elements divisible by \((t+a^2)^2\)
in that quotient. Multiplication by this unit does not change the ideal.

Every \(f\in\mathcal I\) has local order at least two in every removed
class. It is zero in every other CRT quotient: order at least six at an
unremoved ordinary node, and at least two at zero. Hence its CRT coordinates
are supported precisely in the permitted subspace \(C\): removed classes'
old zero-weight rows, with orders \(2,3\) for \(\ell=5\), and \(2,3,4\) for
\(\ell=6\). The word "supported" allows some of these coordinates to vanish.

This is an exact statement about the whole CRT decomposition, not just the
vanishing of \(f\) at its own class. CRT is an isomorphism on polynomials
of degree below \(h=\deg D\); therefore the \(d\) displayed polynomials remain
independent as vectors of coordinates in \(C\).

## 3. They annihilate every critical column

Since \(D\mid Wf\), set \(H_f=Wf/D\), a polynomial. Monic division shows
that reduction of the polynomial quotient of \(Wfg/D\) is exactly \(H_fg\).
The actual leading correction functional is

\[
\overline{p\mu(Q)}=
\frac1{24}\sum_{a\in\mathbb F_p}
\left.\frac{d^3}{dv^3}\bigl[v^5Q(-v^2)\bigr]\right|_{v=a}. \tag{4}
\]

At every ordinary node, local order at least four in \(t+a^2\) kills
its contribution. At zero the factor \(v^5\) kills it automatically.
We check \(Q=H_fg\) for **all** old zero-weight columns and all old weight
\(-1/2\) columns, the only critical sectors in the subsequent lift.

At a removed node of size five, \(H_f\) has order at least

\[
5+2-(5-1)=3.
\]

An own-class old zero column has order at least two, and the own-class
old \(-1/2\) column has order one. Either product has order at least four.
Other-class columns have full local order four, so also cause no problem.

At a removed node of size six, \(H_f\) has order at least two. Own-class
critical columns have order at least two; there is no own-class old
\(-1/2\) column. Other-class columns have local order five. Again all
products have order at least four.

At an unremoved node of size \(\ell\in\{4,5,6\}\), \(H_f\) has order
at least \(6-\ell\). Own-class old zero columns have orders
\(\ell-2,\ell-1\), giving product order at least four. There are no
own-class old \(-1/2\) columns. Every other critical column has local
order \(\ell\), giving product order at least six.

Thus

\[
L(f,g)=0\pmod p\quad
(f\in\mathcal I,\ g\text{ any old zero or old }-1/2\text{ row}). \tag{5}
\]

The statement includes couplings between different CRT classes and all
local unit factors. Formula (1) was chosen to force exactly these orders:
six at unremoved nodes compensates for a denominator of order \(\ell\)
and a critical column of order \(\ell-2\).

### The kernel is exactly this ideal

In fact no additional independent kernel vector is available **inside this
selected space \(C\)**. Suppose \(f\in C\), and set \(H_f=Wf/D\). The same
order checks show that its correction against a fixed unremoved zero column
can receive a contribution only from that column's own node. At all other
unremoved nodes that column has order at least four; at removed nodes its
full local order combines with \(H_f\)'s order to kill the derivative.

At the remaining node write \(u=t+a^2\). The combined contribution at
\(v=\pm a\) of the column of local order \(i\) has the form

\[
[u^3]\,\eta_a(u)H_f(u)u^i,\qquad \eta_a(0)\ne0,
\]

where the unit includes the actual CRT factor and the change of local
coordinate. Both signs contribute equally, not oppositely. The explicit
formula in Section 5 also proves the nonzero constant coefficient.

For \(\ell=4\), the order-three column forces \(H_f(0)=0\); the
order-two column then forces its linear coefficient to vanish. For
\(\ell=5\), the order-three column forces \(H_f(0)=0\), and the
order-four column is automatically zero. For \(\ell=6\), both zero
columns are automatically zero. Thus the full condition is exactly

\[
\operatorname{ord}_{-a^2}H_f\ge6-\ell_a
\quad\Longleftrightarrow\quad
\operatorname{ord}_{-a^2}f\ge6
\qquad(a>N).
\]

Together with \(C\)'s already prescribed orders two at removed nodes and
at zero, this says precisely \(F\mid f\). Therefore

\[
\boxed{\ker(L_{C,Z\cup T}^{\mathsf T})=\mathcal I,
\qquad\dim\ker=d.}
\]

In particular the target-column rank bound in the earlier proof is attained
for every admissible prime, not just the seven computed cases. This proves
optimality of that particular selected-row kernel construction. It does
**not** give an upper bound on the full determinant valuation or rule out
different rows, higher-order cancellations, or other normalizations.

## 4. Why these finite-field vectors give a full \(p^d\)

Here is the lifting argument, included to make clear that (5) is not merely
a self-block nullity claim. The source matrix is \(G=A+p^{-1}L\), where
the Gauss valuations of \(A_{ij}\) are at least \(w_i+w_j\), and \(L\) is
integral and independent of \(X\). Replace every old zero weight by
\(-1/2\), obtaining weights \(w_i'\). All new weights are at most \(-1/2\).
The scaled matrix \(B_{ij}=p^{-w_i'-w_j'}G_{ij}\) is integral over
\(\mathbb Q_p(\sqrt p)\), coefficientwise in \(X\).

For a selected row, against an old zero or old \(-1/2\) column,

\[
B_{ij}=pA_{ij}+L_{ij},
\]

and \(A_{ij}\) is integral: its old lower bound is either zero or
\(-1/2\), and its coefficients have integer \(p\)-adic valuations.
Against a negative half-integer column of weight \(-k-1/2\), \(k\ge1\),
the scaled \(A\) and \(L\) contributions are automatically divisible by \(p\).
Against an integer-weight column they are automatically divisible by
\(\sqrt p\).

All the vectors in (3) have coordinates only in \(C\), whose new weights
are the identical value \(-1/2\). Their independence modulo \(p\) allows
them to be completed to an invertible matrix over \(\mathbb F_p\) on \(C\).
Lift its entries to \(\mathbb Z_p\): its determinant is a unit, so the
result is a unimodular transformation **inside this one common-weight
space**. Apply it to rows and columns simultaneously. The chosen \(d\)
rows now have entries divisible by \(p\) against every half-integer-weight
column, and by \(\sqrt p\) against every integer-weight column.

Divide each of the \(d\) selected rows and columns by \(\sqrt p\).
Selected-selected entries were divisible by \(p\), so remain integral;
all other affected entries remain integral as well. The determinant has
therefore gained at least \(d\) in its \(p\)-adic valuation. There is no
assertion of a rank drop over \(\mathbb Q_p\), nor any need for the
complementary block to be invertible.

Finally the source's rank bound exceeds the number \(z\) of old zero rows
in this range. In the notation \(u=(K+N-3p+1)_+\), they are

\[
r=K+4N-2p+2,\qquad z=p+u,\qquad r-z=d+N+1-u>0
\]

because \(u\le N\). Hence the baseline from replacing all zero weights is
the source's original refined bound \(\gamma_p^{\rm out}=2\sum w_i-z\).
We have proved

\[
\boxed{v_p^G(\Delta_K)\ge\gamma_p^{\rm out}+K+3N-3p+1.}
\]

The usual prime-number-theorem conversion gives \(27/3200\), as in the
earlier note. All other local-prime ranges, the real estimate, and source
integrality constructions remain independent dependencies.

## 5. The finite-field jet form behind the construction

This section is a general algebraic observation, not a claim of a further
determinant improvement. For odd \(s\ge3\) and \(p>s\), the analogous
polynomial moment is
\(\mu_s(Q)=\mathcal L_B(D_v^{s-2}[v^sQ(-v^2)])/(s-1)!\),
where \(\mathcal L_B(v^j)=B_j\). Its correction is

\[
\overline{p\mu_s(Q)}=
\frac1{(s-1)!}\sum_{a\in\mathbb F_p}
D_v^{s-2}\big[v^sQ(-v^2)\big](a).                          \tag{6}
\]

At a nonzero pair \(\pm a\), the two contributions are equal. Set
\(u=t+a^2=-2az-z^2\), and let \(z=z(u)\in u\mathbb F_p[[u]]\) be the
formal inverse. If locally \(Q=u^\kappa\eta(u)f(u)g(u)\), with
\(\kappa\ge0\) and \(\eta(0)\ne0\), the pair's contribution is

\[
[u^{L-1}]\,\omega(u)f(u)g(u),\qquad L=s-1-\kappa,          \tag{7}
\]

where, if \(L>0\),

\[
\omega(u)=\frac2{s-1}(a+z(u))^s z'(u)
                    \left(\frac{u}{z(u)}\right)^{s-1}\eta(u),
\quad
\omega(0)=\frac2{s-1}a^s(-2a)^{s-2}\eta(0)\ne0.           \tag{8}
\]

More explicitly, let \(b(u)^2=a^2-u\) with \(b(0)=a\). Then

\[
z=b-a,\qquad
\omega(u)=-\frac{[b(u)(a+b(u))]^{s-1}}{s-1}\eta(u),
\qquad \omega(0)=-\frac{2^{s-1}a^{2s-2}}{s-1}\eta(0).
\]

This exhibits the coordinate-change factor as a fixed scalar times a square
unit, because \(s-1\) is even.

To prove this, the derivative at \(a\) is \((s-2)!\) times the coefficient
of \(z^{s-2}\), and formal residue substitution gives

\[
[z^{s-2}]A(z)B(u(z))
=[u^{s-2}]A(z(u))z'(u)(u/z(u))^{s-1}B(u).
\]

When \(\kappa\ge s-1\), the contribution is zero. Otherwise (7) is a
nondegenerate form on \(\mathbb F_p[u]/u^L\). Since two is invertible,
\(\omega/\omega(0)\) has a unique square root with constant coefficient
one; multiplying both arguments by that square root reduces the form to

\[
\omega(0)[u^{L-1}]fg.
\]

Its matrix is a nonzero scalar times the anti-identity. It is a direct
sum of hyperbolic planes if \(L\) is even; if \(L\) is odd, one additional
one-dimensional form remains. Its Witt index is \(\lfloor L/2\rfloor\),
and

\[
(u^jA)^\perp=u^{L-j}A,\qquad A=\mathbb F_p[u]/u^L.
\]

Thus powers of the local parameter give explicit isotropic ideals. This
is the Hermite-interpolation connection: prescribe enough zero derivatives
at all nodes, and the surviving global polynomial space is a principal
ideal whose dimension is a degree deficit.

**Important qualification.** The effective length need not be \(s-1\):
the amplitude can vanish. At a removed class of size \(\ell\) in the
present determinant, \(W/D\) has order \(\kappa=6-\ell\), so for ζ(5)
the lengths are three at \(\ell=5\), and four at \(\ell=6\).
Also \(W/D\) has poles at unremoved nodes. Formula (7) cannot be applied
there without first cancelling them. On the denominator-clearing space
\(\{f:D\mid Wf\}\), each argument has local order at least \(\ell\)
at such a node, leaving amplitude order \(\ell\). This distinction
prevents an erroneous extension to ζ(7).

## 6. What the same idea does, and does not, give for ζ(7)

Keep the same \(r=6,N=3K/40\) denominator family and the band above, but
replace the third derivative by a fifth derivative. Clearing the denominator
alone no longer kills couplings between selected rows: unremoved sizes
four and five leave effective jet lengths two and one. In particular,
vanishing of a removed class's own self-block does **not** imply that
different selected rows pair to zero globally.

A sufficient uniform repair is the explicit ideal generated by

\[
F_7(t)=t^2\prod_{a=1}^N(t+a^2)^3
                 \prod_{a=N+1}^{m}(t+a^2)^8.               \tag{9}
\]

At removed size five its product with an old \(-1/2\) column has order
\(5+3-4+2=6\); at removed size six its product with an old zero column
has order \(5+3-5+3=6\). At an unremoved node, the smallest own-class
zero-row order is \(\ell-2\), so the product order is
\(8-\ell+(\ell-2)=6\). Hence all critical fifth derivatives would vanish.
Its CRT support lies in old zero rows, just as required.

But

\[
\deg F_7=4p-2-5N,\qquad h-\deg F_7=K+4N-4p+2<0.          \tag{10}
\]

For the strict inequality, \(K\le3p-1\), \(4N=3K/10\), and \(K\ge40\)
give \(p\ge17\) and
\(K+4N-4p+2\le-\tfrac1{10}p+\tfrac7{10}<0\).
Thus this sufficient ideal contains **no** nonzero polynomial of degree
below \(h\). This is a limitation of this construction, not an impossibility
theorem for other kernels or for irrationality of ζ(7). Moreover the source's
ζ(7) rank bound can be smaller than the number of zero rows, so a prospective
kernel must also recover that baseline difference before furnishing a gain.

## 7. Reproducible checks

`python codex/computations/check_outer_Hermite_ideal.py` constructs (1),
the actual reduced source CRT polynomials, and the Bernoulli correction.
It verifies that every displayed vector lies in \(C\), that they are
independent, that the denominator cancels, and that all critical columns
(including old \(-1/2\) rows) are annihilated. Finite checks complement
the proof; they do not replace its all-prime or coefficientwise arguments.
The seven cases from the earlier matrix check pass, totaling 5,880 critical
correction entries. The same script also checks 2,970 local jet identities
and the nondegenerate anti-triangular Gram forms directly over finite fields.
