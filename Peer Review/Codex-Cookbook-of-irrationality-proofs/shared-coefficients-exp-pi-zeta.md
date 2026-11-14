# The same integer coefficients behind exponential, pi-squared and zeta proofs

Status: complete elementary exponential and pi-squared arguments below, with exact
polynomial checks. The constructions are classical Hermite/Padé constructions;
no historical novelty is claimed. The point is a concrete coefficient-level link
to the binomial kernels already used for Apéry's recurrences.

## 1. Start with the coefficients, not a list of unrelated constants

Let
\[
\phi_k(n)=\binom nk\binom{n+k}k.
\]
These are the quadratic-grid Newton coefficients and Zudilin's signed simple-pole
residues. Consider the two integer polynomials
\[
P_n(z)=\sum_{k=0}^n\phi_k(n)k!z^{n-k},\qquad Q_n(z)=P_n(-z).
\tag{1}
\]
Equivalently,
\[
Q_n(z)=\sum_{j=0}^n\frac{(2n-j)!}{(n-j)!j!}(-z)^j.
\]
The binomial form makes integrality immediate. The first pairs are
\[
P_0=Q_0=1,\quad P_1=2+z,\quad Q_1=2-z,
\]
\[
P_2=12+6z+z^2,\quad Q_2=12-6z+z^2.
\]
They are reversed Bessel/diagonal exponential Padé polynomials in a convenient
integer normalization. Those names are optional: the defining sums and the next
identity prove everything needed.

## 2. Endpoint derivatives produce those coefficients exactly

Set
\[
f_{n,z}(x)=\frac{x^n(z-x)^n}{n!}.
\]
Its only nonzero endpoint derivatives at0 have orders n+k, and
\[
f_{n,z}^{(n+k)}(0)=(-1)^k\phi_k(n)k!z^{n-k}.
\tag{2}
\]
This follows by expanding(z-x)^n. Reflection x↦z-x gives
\(f^{(j)}(z)=(-1)^j f^{(j)}(0)\). Repeated integration by parts, terminating
because f has degree2n, therefore yields the exact identity
\[
\boxed{\quad
Q_n(z)e^z-P_n(z)
=\frac{(-1)^n}{n!}\int_0^z x^n(z-x)^n e^x\,dx.
\quad}
\tag{3}
\]
The integral can be along the straight segment in the complex plane. Alternatively
one can first prove it for real z and use equality of entire functions. No
asymptotic expansion or unnamed approximation theorem is needed.

The zeros at both endpoints have the same job as zeros in Zudilin's rational
function: they remove unwanted boundary terms. Here the surviving endpoint
derivatives are integral because their coefficients are exactly(2).

## 3. Irrationality of e^r for every nonzero rational r

First let r=a/b>0 with positive integers a,b. The numbers
\(b^nQ_n(r)\) and \(b^nP_n(r)\) are integers. Formula(3) has a nonzero sign,
and
\[
0<|Q_n(r)e^r-P_n(r)|
\le\frac{r e^r}{n!}\left(\frac{r^2}{4}\right)^n.
\tag{4}
\]
The bound uses0<=x(r-x)<=r²/4 on[0,r]. Multiplication by b^n still leaves
an expression tending to0, since any fixed exponential divided by n! tends to0.
If e^r=c/d were rational, then
\[
db^n(Q_n(r)e^r-P_n(r))
\]
would be a nonzero integer tending to0, a contradiction. For negative rational r,
rationality of e^r would imply rationality of its reciprocal e^(-r), which has
just been excluded.

This is an irrationality proof, not a claim that this elementary argument proves
the full Lindemann theorem.

## 4. The same identity proves irrationality of pi^2

Put z=i*pi in(3), using e^(i*pi)=-1. Along x=i*u,
\[
Q_n(i\pi)e^{i\pi}-P_n(i\pi)
=\frac{i}{n!}\int_0^\pi u^n(\pi-u)^n e^{iu}\,du.
\]
The cosine part vanishes by the reflection u↦pi-u. Thus
\[
\boxed{\quad
P_n(i\pi)+Q_n(i\pi)
=\frac1{n!}\int_0^\pi u^n(\pi-u)^n\sin u\,du>0.
\quad}
\tag{5}
\]
Only even powers survive in P_n(z)+P_n(-z), so its left side is an integer
polynomial in pi² of degree floor(n/2). If pi²=a/b with integers a,b>0,
\[
b^{\lfloor n/2\rfloor}[P_n(i\pi)+Q_n(i\pi)]
\]
is a positive integer. But(5) bounds it above by
\[
\frac\pi{n!}\,b^{\lfloor n/2\rfloor}
\left(\frac{\pi^2}{4}\right)^n\longrightarrow0,
\]
a contradiction. Irrationality of pi follows as well.

Complex notation makes the connection to the exponential transparent. The same
identity can be obtained using real integration by parts against sin u, so no
complex-analysis theorem is essential.

This proves pi², not every rational power pi^r by the same argument. For a
nonzero rational exponent r, the standard general statement follows from the
transcendence of pi: if pi^r were algebraic, an integer power relation would make
pi algebraic. The stronger theorem and its machinery are discussed separately
in pi-powers.md; irrationality of pi² alone would not justify that conclusion.

## 5. A second-order recurrence and the shared coefficient dictionary

Both P_n and Q_n satisfy, for n>=1,
\[
R_{n+1}(z)=2(2n+1)R_n(z)+z^2R_{n-1}(z).
\tag{6}
\]
It follows directly by comparing coefficients in(1). The two initial conditions
supply different solutions. Their Casoratian is
\[
Q_nP_{n+1}-Q_{n+1}P_n=2z(-z^2)^n.
\]
Again a nonzero determinant records the difference between two solutions, as in
the Apéry recurrence. Here the real integral already gives the needed sign and
factorial smallness, so a separate recurrence estimate would be redundant.

The coefficient dictionary is now concrete:

| Construction | Operation on the same phi_k(n) |
|---|---|
| log((a+1)/a) | Sum phi_k(n) a^k, a shifted Legendre value |
| e^r | Sum phi_k(n) k! (±r)^(n-k), from endpoint derivatives |
| pi² | Take the even part of the preceding polynomial at i*pi |
| zeta(2) | Pair phi_k(n) with binom(n,k) and take the diagonal |
| zeta(3) | Pair phi_k(n) with itself and take the diagonal |

This does not say every constant follows from one universal theorem. It identifies
an actual shared piece of integer algebra. The choice of partner or functional
changes the analytic representation, the recurrence and the denominator cost.

For exponential and pi² proofs, the numerator bound is a fixed exponential
while n! sits in the denominator of the integral. That large separation makes
the contradiction easy. For Apéry's zeta forms, both real decay and denominator
growth are exponential, so their exact rates and cancellation structure matter.
This explains why the same-looking endpoint zeros can support proofs of very
different difficulty.

Reproduce(2), the endpoint form of(3), recurrence(6), and the even-power property
with computations/check_shared_coefficient_family.py.

## 6. Tangent, logarithms and hyperbolic tangent

The same polynomials also give Lambert's tangent continued fraction and an
irrationality proof for tan r whenever r is nonzero rational. Writing
P_n(iθ)=C_n(θ)+iS_n(θ), their centered integral bounds
C_n tan(θ/2)−S_n by a fixed exponential divided by n!, while

    C_n S_(n+1)−C_(n+1)S_n=θ^(2n+1)

ensures that two consecutive errors cannot both vanish. This last point is
essential when the centered cosine integrand changes sign. The complete proof,
including convergence of the finite continued fractions, is in
[shared-coefficients-tangent.md](https://github.com/BoomydayCoder/Mathathon/blob/03cc75d694220d0fba888c637e573a29799e6d5b/codex/problem131/shared-coefficients-tangent.md).

Two immediate corollaries of the exponential result are irrationality of log q
for positive rational q≠1, and irrationality of tanh r for nonzero rational r.
They follow respectively from e^(log q)=q and
 e^(2r)=(1+tanh r)/(1−tanh r). These statements do not imply irrationality of
each separate sinh r or cosh r by the same reasoning.

Independent review of the original exponential/pi-squared proof, signs,
normalizations and pi-power scope is recorded in
[shared-coefficients-review.md](shared-coefficients-review.md). Its further
all-parameter recurrence calculation supplements the original finite checks.
