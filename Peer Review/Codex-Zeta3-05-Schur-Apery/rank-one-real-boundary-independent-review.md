# Independent review: finite boundary cancellation

Verdict: accepted. Root independently reviewed
[rank-one-real-boundary-proof.md](rank-one-real-boundary-proof.md).
It replaces the contour dependency in the accepted original Schur–Apéry proof.
The claim is explicitly restricted to deg_(x²) f+deg_(x²) g<=M-1.

For p of degree<=M-1, p(i sqrt(t))/D_M(t) is proper. Its residue at -j²
is p(j)/D_M'(-j²), with the derivative weight exactly as displayed. This
checks both the sign and the admissible degree range of the finite model.
The adjacent alpha-weight ratio is
-(j+1)^4(M-j)/(j^4(M+j+1)); multiplying A(-j-1) gives A(j).
The harmonic update gives rho_j-rho_(j+1)=(2j+1)/(2j³(j+1)³).
The absent end edges have A(M)=A(-1)=0. Summation by parts therefore gives
the stated delta_j coefficient, with negative sign for odd j.

Because f and g are even, K(-x-1)=-K(x). Also delta_(-j-1)=-delta_j;
the two products pair with the same sign. The leading terms of the two
products defining K cancel, so deg K<=2(r+s)-1<=2M-3. Multiplication by
2(2M-1)! turns delta into alternating binomial coefficients on precisely
2M consecutive integer nodes. Their sum against K vanishes. Removing the
j=0 term produces -delta_0 K(0), which is the required -M/[2(M!)²] defect.
The M=1 case and zero/constant inputs cause no exception.

For the Heun symmetry argument, deg f,deg g<=n and M=2n+2. The raised
pairs (Xf,g) and (f,Xg) have total degree<=2n+1=M-1, so the restricted
identity applies exactly where needed. The two resulting boundary wedges
sum to the original wedge because X(0)+X(1)=1. Subtracting the D term
cancels this boundary defect. This proves the required symmetry on the finite
space; no assertion about a global operator domain is necessary.

The rest of the all-index bridge already passed two independent audits.
This algebraic replacement therefore yields the promised entirely real proof
from the previously established real resolvent identity. It does not itself
claim that the rational functional is positive for arbitrary X; positivity
is used only at X=zeta(3). The author's 364 exact finite cases and free-symbol
balance checks supplement the written proof, rather than extend its range.
