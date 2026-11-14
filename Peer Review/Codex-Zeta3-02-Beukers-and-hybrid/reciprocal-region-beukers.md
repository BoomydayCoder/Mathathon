# The reciprocal involution sends the positive octant to a Beukers region

Status: proved here by explicit rational changes of variables. The Laurent
polynomial and its Beukers period are published mathematics; this note makes
the real-domain interpretation of the corrected reciprocal map explicit.
The proposed region and reciprocity were supplied by the parent research
agent; the inverse, cube parametrization, orientation and convergence below
were independently derived by the simplify agent.

Put
\[
 C=(1+x)(1+y),\quad Q=C+xyz,\quad
 \Lambda=\frac{C(1+z)Q}{xyz},\quad
 T=\left(-\frac z{1+z},-\frac C Q,-\frac x{1+x}\right).
\]
The [Fano dictionary](https://github.com/BoomydayCoder/Mathathon/blob/03cc75d694220d0fba888c637e573a29799e6d5b/codex/problem131/fano-v12-apery-gamma-dictionary.md) proves
\(T^2=1\), \(\Lambda\circ T=\Lambda^{-1}\), and
\(T^*(dx\,dy\,dz/(xyz))=\Lambda^{-1}dx\,dy\,dz/(xyz)\).
All are rational identities on the open set on which the expressions exist.

## 1. The exact real image

Define
\[
 R=\{(a,b,c)\in(0,1)^3:(1-a)(1-b)<abc\}.
\]
The map \(F=-T\) is an orientation-preserving smooth bijection from the
positive octant to \(R\). Indeed,
\[
 (a,b,c)=\left(\frac z{1+z},\frac C Q,\frac x{1+x}\right),\qquad
 (x,y,z)=\left(\frac c{1-c},
 \frac{(1-a)(1-b)}{abc-(1-a)(1-b)},\frac a{1-a}\right).       \tag{1}
\]
The displayed inverse is positive precisely on \(R\), and direct
substitution verifies both compositions. Thus this is a statement about
entire open real domains, not merely a formal birational identity.

Set
\[
 \Phi(a,b,c)=\Lambda(-a,-b,-c)
 =\frac{(1-a)(1-b)(1-c)[abc-(1-a)(1-b)]}{abc}>0.             \tag{2}
\]
Then
\[
 \Phi\circ F=\Lambda^{-1},\qquad
 F^*\frac{da\,db\,dc}{abc}
   =\Lambda^{-1}\frac{dx\,dy\,dz}{xyz}.                   \tag{3}
\]
The multiplier in (3) is positive, so the ordinary Jacobian of \(F\)
is positive. In contrast, \(T\) has negative ordinary Jacobian into
the negative octant: the negative target product in its logarithmic
volume form supplies the other sign. Working with \(F\) avoids that
orientation ambiguity.

Consequently, for every integer \(n\ge0\),
\[
 I_n:=\int_{(0,\infty)^3}\Lambda^{-n-1}\frac{dx\,dy\,dz}{xyz}
 =\int_R\Phi^n\frac{da\,db\,dc}{abc}.                       \tag{4}
\]
The usual change-of-variables theorem first applies on compact subdomains;
nonnegative monotone exhaustion extends it to these open domains, even
before finiteness has been proved.

## 2. A rational cube parametrization

For \(0<u,v,w<1\), write \(D=1-w+uvw\). The map
\[
 (a,b,c)=\left(w,\frac{1-w}{D},u\right)                    \tag{5}
\]
is a bijection from the open cube to \(R\). Its inverse is
\[
 u=c,\qquad w=a,\qquad v=\frac{(1-a)(1-b)}{abc}\in(0,1).
\]
Its determinant is \(uw(1-w)/D^2>0\), while \(abc=uw(1-w)/D\).
Thus
\[
 \frac{da\,db\,dc}{abc}=\frac{du\,dv\,dw}{D},\qquad
 \Phi=\frac{u(1-u)v(1-v)w(1-w)}D.                           \tag{6}
\]
Equations (4)–(6) recover exactly the Beukers integral:
\[
 I_n=\int_{(0,1)^3}
 \frac{[u(1-u)v(1-v)w(1-w)]^n}{(1-w+uvw)^{n+1}}
 \,du\,dv\,dw.                                           \tag{7}
\]
For \(n=0\), integrating in \(w\) and expanding a nonnegative series gives
\[
 I_0=\int_0^1\!\int_0^1\frac{-\log(uv)}{1-uv}\,du\,dv
 =\sum_{k\ge0}\frac2{(k+1)^3}=2\zeta(3).                   \tag{8}
\]
Tonelli justifies the interchange. In particular, the logarithmic volume
of \(R\) is finite, despite its approach to coordinate boundaries.

## 3. Minimum–maximum duality

The unique positive minimum of \(\Lambda\) is
\(\lambda=17+12\sqrt2\), attained at
\((x,y,z)=(1/\sqrt2,1/\sqrt2,1+\sqrt2)\).
One elementary verification is to differentiate at this point and observe
that \(\Lambda(e^X,e^Y,e^Z)\) is a positive sum of exponentials whose
exponent vectors span three-dimensional space. Its Hessian is therefore
positive definite everywhere; a stationary point is its unique global
minimum. No assertion of compactness is needed for this argument.

The bijection \(F\) therefore proves
\[
 0<\Phi\le\rho:=17-12\sqrt2,
\]
with unique maximum at
\[
 (a,b,c)=(1/\sqrt2,1/\sqrt2,\sqrt2-1).
\]
In cube coordinates this is
\((u,v,w)=(\sqrt2-1,\sqrt2-1,1/\sqrt2)\). Thus
\[
 0<I_n\le2\zeta(3)\rho^n,\qquad n\ge0.                    \tag{9}
\]
This also proves convergence of all the integrals above. The exponential
smallness and the reciprocal positive saddle are literally the same
optimization problem after a rational change of variables.

## 4. Scope

The real chain \(R\) is bounded and open, with boundary; it is not the
closed complex torus defining \(\operatorname{CT}\Lambda^n\). The identity
\(T^*(\Lambda^n\Omega)=\Lambda^{-n-1}\Omega\) explains the exponent shift
pointwise but does not make the two period sequences equal or remove the
inhomogeneous boundary term. No new irrationality theorem or historical
novelty is asserted.

The mirror polynomial and reciprocal-map context are in Golyshev–Kerr–Sasaki,
*Apéry extensions*, arXiv:2009.14762v3, §5.4. The correction to that version's
printed map is localized and checked in the companion dictionary. The
integral (7) is Beukers' published integral. See the dictionary's archived
primary references and source manifest.

Reproducible verification: `python codex/computations/check_reciprocal_region.py`.
It checks full determinants independently, both inverse maps, the cube
pullback, critical point, strict-convexity exponent support, and the
printed-map rational counterexample. The convergence and domain arguments
are the proofs above, not consequences of finite numerical tests.
