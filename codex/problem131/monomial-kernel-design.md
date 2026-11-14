# Designing the zeta forms by their monomial kernel

Status: elementary general theorem proved below. This reorganizes the previously
proved rational-function arithmetic into a transparent polynomial pairing.
Historical novelty is not claimed.

## 1. The diagonal selects the constant; off-diagonal terms are rational

For an integer s>=2 and polynomials P,Q with integer coefficients and degrees<=n,
put
\[
\mathcal I_s(P,Q)=\frac1{(s-2)!}\int_0^1\!\int_0^1
\frac{P(x)Q(y)(-\log(xy))^{s-2}}{1-xy}\,dx\,dy.
\]
Let \(P=\sum p_kx^k\), \(Q=\sum q_kx^k\), and
\(D_n=\operatorname{lcm}(1,\ldots,n)\), with D0=1. Then
\[
\boxed{\quad
\mathcal I_s(P,Q)=(s-1)\left(\sum_{k=0}^np_kq_k\right)\zeta(s)-B,
\qquad D_n^sB\in\mathbb Z.
\quad}
\tag{1}
\]
Every coefficient and the denominator exponent can be seen in a single monomial
calculation. Write \(H_k^{(r)}=\sum_{j=1}^k j^{-r}\), H0=0. Then
\[
\mathcal I_s(x^k,y^\ell)=
\begin{cases}
(s-1)(\zeta(s)-H_k^{(s)}),&k=\ell,\\[2mm]
\displaystyle\frac{H_\ell^{(s-1)}-H_k^{(s-1)}}{\ell-k},&k\ne\ell.
\end{cases}
\tag{2}
\]
Indeed expand(1-xy)^(-1), and use
\[
\frac{(-\partial_t)^{s-2}}{(s-2)!}
\frac1{(t+k)(t+\ell)}
=\begin{cases}
(s-1)/(t+k)^s,&k=\ell,\\
\big[(t+k)^{-(s-1)}-(t+\ell)^{-(s-1)}\big]/(\ell-k),&k\ne\ell.
\end{cases}
\]
Summation over t>=1 proves(2). The integral is absolutely convergent: bounded
polynomials reduce it to \(\mathcal I_s(1,1)=(s-1)\zeta(s)<\infty\).
Geometric expansion and differentiation can therefore be justified first for
monomials by positive integrands and then by a finite linear combination.
The off-diagonal denominator uses at most D_n^(s-1) from the harmonic sum and
one more D_n from ell-k. The diagonal uses D_n^s. This proves(1).

Thus algebraically isolating one higher zeta value is easy in this setting. The
research problem is to choose P,Q so that the same form is nonzero and shrinks
quickly enough after this denominator is cleared.

## 2. Two deliberate partner choices give the classical denominator sequences

Let
\[
L_n(x)=\frac1{n!}\frac{d^n}{dx^n}[x^n(1-x)^n]
=\sum_{k=0}^n(-1)^k\phi_k(n)x^k,
\quad\phi_k(n)=\binom nk\binom{n+k}k.
\]
For s=2 choose \(P=L_n\), \(Q=(1-y)^n\). Formula(1) gives the zeta coefficient
\[
C_n=\sum_k\phi_k(n)\binom nk.
\]
Integration by parts n times in x, with all endpoint terms zero, gives
\[
\boxed{\quad
(-1)^n\mathcal I_2(L_n,(1-y)^n)
=\int_0^1\!\int_0^1
\frac{[x(1-x)y(1-y)]^n}{(1-xy)^{n+1}}\,dx\,dy>0.
\quad}
\tag{3}
\]
So the binomial partner has an analytic purpose: it supplies the remaining
endpoint factor(1-y)^n after differentiating the kernel produces y^n.
The mixed binomial kernel and its11,3 recurrence are consequently attached to a
specific positive approximation, not selected only because its numbers look nice.

For s=3 choose \(P=Q=L_n\). The coefficient is
\[
2A_n=2\sum_k\phi_k(n)^2.
\]
Its Mellin transform is
\[
\int_0^1x^{t-1}L_n(x)\,dx
=(-1)^n\frac{\prod_{j=1}^n(t-j)}{\prod_{j=0}^n(t+j)}.
\]
Hence \(\mathcal I_3(L_n,L_n)=-\sum_{t\ge1}R_n'(t)\), exactly Zudilin's
normalization. The Beukers triple-integral identity proves its positivity and
exponential bound; alternatively the recurrence and boundedness argument proves
nonvanishing and decay. Those identities are fully proved in the existing guides.

The square in A_n is therefore the diagonal coefficient of a polynomial kernel.
The Mellin transform turns that same pairing into a square of rational functions.
This is an exact bridge between orthogonal polynomials, partial fractions and
integer recurrence solutions.

## 3. Why keeping the square while changing the logarithm can fail

Using L_n twice in(1) isolates every zeta(s), but does not automatically produce
a useful irrationality form. For s=2 it is a sum of positive squares of Mellin
transforms; that alone does not provide a sufficiently fast decay. For s>=3 the
bilinear form is not even positive semidefinite on arbitrary polynomials.
For example its matrix on1,x is
\[
\begin{pmatrix}
(s-1)\zeta(s)&1\\
1&(s-1)(\zeta(s)-1)
\end{pmatrix}.
\]
For s>=3 this determinant is negative, as proved in kernel-positivity-obstruction.md.
Pointwise positivity of the displayed two-variable kernel does not make its
quadratic form positive when P changes sign.

The unmodified higher odd-derivative family is now known to decay only
polynomially, with exact boundary forcing; see forced-recurrence-obstruction.md
and log-kernel-extension.md. Thus formula(1) should be used as the arithmetic
part of a construction, not mistaken for a complete proof for every odd index.

## 4. One endpoint polynomial supplies several different proof mechanisms

The primitive polynomial \(x^n(1-x)^n\) has high-order zeros at both endpoints.
Different operations on it expose different properties:

- Its nth derivative is L_n, whose coefficients are the integer phi_k(n).
- Its Mellin transform creates prescribed consecutive zeros and poles.
- Its endpoint derivatives against exp give the reversed Bessel/Padé polynomials.
- The imaginary-axis specialization of that same exponential identity gives the
  positive sine integral proving irrationality of pi².
- Pairing its coefficients with ordinary binomial coefficients gives the zeta(2)
  denominator sequence; pairing with themselves gives the zeta(3) sequence.

These are linked by actual transformations, not an analogy based on the word
'auxiliary'. What changes between proofs is how one proves sign/nonvanishing and
how the real decay competes with the arithmetic normalization. For exp and pi²,
a factorial beats a fixed exponential. For the zeta forms, a much tighter
comparison between two exponential rates is required.

Read shared-coefficients-exp-pi-zeta.md for the elementary exp/pi² identities and
zeta2-zeta3-balance-unification.md for the coefficient-level recurrence recipe.
