# Second independent audit of the finite real boundary identity

Status: **accepted**, independently of the [root review](rank-one-real-boundary-independent-review.md). This separately checks [rank-one-real-boundary-proof.md](rank-one-real-boundary-proof.md), including the restricted degree range needed in the already accepted [original Schur/adjacent-Apéry theorem](rank-one-original-apery-heun-rii-proof.md). The check does not use the unrestricted contour identity.

## Finite grid and boundary sign

For an even polynomial p of degree at most M-1 in x², p(i sqrt(t))/D_M(t) is proper. Its residue at -j² is exactly p(j)/D_M'(-j²), because p is even. Direct factorization gives

    1/D_M'(-j²)=2(-1)^(j-1)j²/[(M-j)!(M+j)!].

The difference operator preserves even polynomial degree: its apparent poles at x=0 cancel between the reflected terms, and at x=+-1/2 the relevant difference vanishes by evenness. Its coefficients are O(x²), while the reflected first differences cancel their top odd term. Thus the partial-fraction formula applies to f Dg-g Df whenever deg_(x²) f+deg_(x²) g<=M-1.

Writing alpha_j=j²/D_M'(-j²), rho_j=F_j(X)/j², and K(j)=f(j)g(j+1)-g(j)f(j+1), direct finite summation pairs neighboring terms. The coefficient of K(j) is

    alpha_j A_M(j) rho_j-alpha_(j+1) A_M(-j-1) rho_(j+1).

The exact detailed-balance identity makes this alpha_j A_M(j)(rho_j-rho_(j+1)). The supplied rational pole data give

    rho_j-rho_(j+1)=(2j+1)/[2j³(j+1)³].

Consequently the coefficient is delta_j=(-1)^j/[2(M-j-1)!(M+j)!], with the sign stated in the proof. The outer terms vanish because A_M(M)=A_M(-1)=0. In particular, no undefined sample at zero enters the finite-grid summation.

## Finite difference and permitted degrees

Evenness gives K(-j-1)=-K(j), while delta_(-j-1)=-delta_j. Its polynomial degree is at most 2(r+s)-1, since the top product terms cancel. This is at most 2M-3, below the order 2M-1 of the finite difference over the 2M nodes -M,...,M-1. The sum of delta_j K(j) over all nodes is therefore zero, and reflection makes it twice the sum over j=0,...,M-1. Removing its zero-node term proves the boundary coefficient -delta_0=-M/[2(M!)²]. The M=1 case allows only constant inputs and both sides vanish.

For the Heun step, M=2n+2 and f,g have degree at most n in x². Applying the identity to (Xf,g) and (f,Xg) uses degree sums at most 2n+1=M-1; applying it to (f,g) uses at most 2n. These are precisely within the proved range. The skew forms of XD and DX are respectively the boundary forms for (Xf,g) and (f,Xg). Their sum is the original boundary form because X(0)=0 and X(1)=1. The subtractive -D term cancels it. Multiplication by X is symmetric. This proves every symmetry statement actually used by the all-h bridge, including n=0/h=1.

Thus the all-h theorem no longer needs the complex-contour input. The unrestricted all-degree identity is correctly kept separate.

## Executed exact checks

The supplied checker completed with 364 finite-grid boundary identities covering every permitted pair through M=12, plus symbolic detailed-balance and source identities. The mathematical review above independently checks the reduction and degree bounds; the finite calculation is corroboration, not an all-index substitute.
