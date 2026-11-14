# Independent review: the negative-atom continuation

**Verdict:** accepted without correction. This independently audits
[wilson-negative-atom-continuation.md](wilson-negative-atom-continuation.md),
using the already reviewed measure on \(|\tau|\le1/2\) and its Bernoulli
moment formula. The continuation, Jacobi formulas, and polynomial support
obstruction are proved for their stated full parameter ranges. Finite exact
checks corroborate them. These are explicit derivations within classical
Wilson theory, not a claim of historical novelty or an irrationality result.

## 1. Folding and all boundary cases

The continued moments are the polynomials
\[
 M_k(\tau)=\frac{(-1)^k}{(k+1)\tau^2}
 \bigl[B_{2k+2}(1/2+\tau)-B_{2k+2}(1/2)\bigr].
\]
At zero the numerator has a double zero, so the polynomial limiting value
is well defined. For \(\tau>1/2\), set
\(K=\lfloor\tau+1/2\rfloor\), \(\epsilon=\tau-K\), and
\(v_j=\tau-j-1/2\). Telescoping Bernoulli polynomials from
\(1/2+\epsilon\) to \(1/2+\tau\) gives exactly
\[
 M_k(\tau)=\frac{\epsilon^2}{\tau^2}M_k(|\epsilon|)
       +\sum_{j=0}^{K-1}\frac{2v_j}{\tau^2}(-v_j^2)^k.
\]
The sign \((-1)^k\), factor two, and number of shifts all agree. The
evenness about \(1/2\) justifies replacing \(\epsilon\) by its absolute
value. Every \(v_j\) is nonnegative. The last is zero precisely at a
half-integer parameter, and then its weight is zero: no extra atom at
zero is present. At an integer parameter \(\epsilon=0\), the continuous
component has coefficient zero, irrespective of the nonzero limiting
base density. It must be omitted rather than evaluated by division by zero.

The mass calculation is
\[
 \epsilon^2+2\sum_{j=0}^{K-1}v_j
   =(\tau-K)^2+K(2\tau-K)=\tau^2.
\]
Thus the displayed measure is a positive probability measure on the real
line. Its nonzero atoms are distinct and strictly negative. Its continuous
density agrees with the periodic hyperbolic expression for noninteger
\(\tau\); that density alone has mass \(\epsilon^2/\tau^2\), not one.
This proves a representing measure by exact moments, with no invalid
analytic continuation under an integral sign.

## 2. Gram identities and integer degeneration

Both sides of
\[
 \widetilde\nu_\tau(D_jD_k)
   =\frac{d_j(\tau)d_k(\tau)}{j+k+1},\qquad
 D_k(t)=\prod_{a=0}^{k-1}(t+(a+1/2)^2),\quad
 d_k(\tau)=\prod_{a=1}^k(a^2-\tau^2),
\]
are polynomials in \(\tau^2\). Their established equality on the original
open interval therefore proves the identity at every real parameter.
For noninteger positive \(\tau\), all \(d_k\ne0\), and
\(Q_k=D_k/d_k\) gives a genuine full polynomial isometry with the
Hilbert Gram matrix \(1/(j+k+1)\). Negative signs of individual
\(d_k\) do not affect this conclusion. The positive continuous component
ensures positive definiteness in every finite degree.

At \(\tau=N\in\mathbb N\), the measure is instead
\[
 \sum_{r=0}^{N-1}\frac{2r+1}{N^2}\delta_{-(r+1/2)^2}.
\]
Its moment matrix through degree \(d\) has rank \(\min(d+1,N)\):
the weighted evaluation matrix has distinct nodes and positive weights,
so Vandermonde rank gives the result. \(D_N\) vanishes at every node,
and \(d_N=0\). Consequently the Hilbert identity remains meaningful
for \(j,k<N\), but dividing by \(d_N\) would be invalid. The manuscript
correctly uses a finite quotient at these parameters.

## 3. Every Jacobi function and overlap

The transported operator acts on a power by the particularly useful formula
\[
 T_\tau x^a=-(a+1/2)^2x^a+
             ((a+1)^2-\tau^2)x^{a+1}.
\]
For an integer \(j\ge0\) with \(j<\tau-1/2\), put
\(s=\tau-j-1\), \(\alpha=2s+1>0\), and
\(f_j=x^sP_j^{(\alpha,0)}(1-2x)\). The shifted Jacobi equation
then gives
\[
 T_\tau f_j=-(s+1/2)^2f_j.
\]
In particular the constant term and the highest-power cancellation have
the correct signs: \((s+j+1)^2-\tau^2=0\). The strict restriction on
\(j\) is exactly \(s>-1/2\), needed for square integrability. At a
half-integer boundary, the excluded \(s=-1/2\) function corresponds
precisely to the omitted zero-weight atom.

The shifted Rodrigues formula is
\[
 P_j^{(\alpha,0)}(1-2x)
   =\frac{x^{-\alpha}}{j!}
       \frac{d^j}{dx^j}[x^{j+\alpha}(1-x)^j].
\]
Integrating by parts against \(x^s\) gives
\((s+1)_j\mathrm B(s+1,j+1)/j!=1/\tau\). For the norm,
write \(v_j=v_j(0)+xw\), where \(\deg w<j\). Jacobi orthogonality
with weight \(x^\alpha\) removes the \(w\) term. A second Rodrigues
integration gives
\[
 \int_0^1x^{\alpha-1}v_j(x)^2\,dx
 =\frac{(\alpha+1)_j}{j!}\mathrm B(\alpha,j+1)
 =\frac1\alpha.
\]
Thus
\[
 \int_0^1f_j\,dx=\frac1\tau,\qquad
 \int_0^1f_j^2\,dx=\frac1{2\tau-2j-1},\qquad
 \frac{|\langle1,f_j\rangle|^2}{\|f_j\|^2}
       =\frac{2\tau-2j-1}{\tau^2}.
\]
All integration-by-parts boundary terms vanish. At zero their smallest
powers are \(x^{s+1}\) or \(x^\alpha\), both with positive exponent;
at one the differentiated Rodrigues factor retains the required zero.
For distinct indices, the Green boundary term is
\(O(x^{s_i+s_j+1})\) at zero and vanishes at one. Since the eigenvalues
are distinct, the functions are mutually orthogonal. At integer \(N\)
these are \(N\) nonzero ordinary polynomials in the invariant space of
degree at most \(N-1\), so they give its full orthogonal eigenbasis.

At noninteger parameters this review accepts the stated differential and
integral identities. Neither the proof nor this review silently assigns a
domain to an unbounded closed operator; the atom representation was already
established without that extra spectral assertion.

## 4. Why no alternative positive half-line measure is possible

For \(1/2<\tau<1\), set \(b=1-\tau\in(0,1/2)\) and
\[
 f_L(x)=\sum_{k=0}^L\frac{(b)_k}{k!}(1-x)^k.
\]
Positive coefficients give \(0<f_L\le x^{-b}\) and
\(|f_L'|\le b x^{-b-1}\), with pointwise convergence of both.
The dominating functions for the potential and kinetic integrals are
constant multiples of \(x^{-2b}\) and \((1-x)x^{-2b}\). They are
integrable because \(2b<1\). Hence the polynomial quadratic forms
converge to
\[
 -\frac{(\tau-1/2)^2}{2\tau-1}<0.
\]
The limit may also be checked directly by integrating powers; its Green
boundary factor \((\tau-1)x^{2\tau-1}(1-x)\) vanishes at both ends.
Some finite \(f_L\) must therefore have negative quadratic form.

To make the support contradiction explicit, the isometry sends
\(x^k\) to \(Q_k(t)\), and
\[
 tQ_k=((k+1)^2-\tau^2)Q_{k+1}-(k+1/2)^2Q_k.
\]
Thus it intertwines multiplication by \(t\) with \(T_\tau\).
The finite witness becomes a polynomial \(p\) satisfying
\(\widetilde\nu_\tau(tp^2)<0\). Any positive measure on
\([0,\infty)\) with the same moments would give a nonnegative value,
an immediate contradiction. For \(\tau\ge1\), including the degenerate
integer parameters, the simpler witness \(p=1\) already works:
\(M_1(\tau)=1/4-\tau^2/2<0\).

This argument correctly excludes *every* positive half-line representation,
not merely the particular representation found in the manuscript. Ordinary
Hamburger positivity survives; the obstruction is support.

## 5. A finite binomial form of the integer Hilbert identity

One useful explicit consequence requires no additional analytic theorem.
For \(0\le k,r<N\), direct evaluation gives
\[
 q_{N,k}(r):=Q_k(-(r+1/2)^2)
 =\begin{cases}
 \displaystyle\frac{N(r+k)!(N-k-1)!}{(r-k)!(N+k)!},&k\le r,\\
 0,&k>r.
 \end{cases}
\]
Equivalently,
\[
 q_{N,k}(r)=\frac{N}{2k+1}
       \frac{\binom{r+k}{2k}}{\binom{N+k}{2k+1}}.
\]
The finite atomic Gram identity therefore becomes the exact binomial sum
\[
 \boxed{\sum_{r=0}^{N-1}\frac{2r+1}{N^2}
       q_{N,j}(r)q_{N,k}(r)=\frac1{j+k+1},
       \qquad 0\le j,k<N.}
\]
Its evaluation matrix is triangular with nonzero diagonal, directly
exhibiting the finite isometry. For polynomials \(f,g\) of degree below
\(N\), this is a finite quadrature identity for their *transported*
coefficient transforms. It should not be mistaken for ordinary quadrature
that evaluates \(f\) and \(g\) at a single common set of scalar nodes:
in general \(q_{N,k}(r)\ne q_{N,1}(r)^k\).

## 6. Reproducible independent checks

[check_wilson_atoms_independent.py](../computations/check_wilson_atoms_independent.py)
uses standard-library exact fractions and does not import the author's code.
It reconstructs Jacobi coefficients from the triangular differential
operator, independently compares them with Rodrigues coefficients, and checks:

- 192 continued moment identities, including parameters on both sides of
  integers, half-integers, and close to \(1/2\);
- 35 Jacobi functions, their differential equations, norms, overlaps, and
  mutual orthogonality;
- 385 finite binomial Hilbert identities, integer triangular ranks, and
  half-integer omission of the zero node;
- two finite negative polynomial quadratic forms, calculated directly from
  their bilinear coefficient matrix.

All pass. The author's checker also passes: 156 moments, 432 Gram identities,
21 rank checks, three finite negative witnesses, and 32 Jacobi functions.
Finite checks supplement the all-parameter arguments above.
