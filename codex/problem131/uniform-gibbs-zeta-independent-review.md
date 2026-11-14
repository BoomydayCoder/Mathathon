# Independent audit of the uniform-sampler zeta correlations

Status: **accepted** for [uniform-gibbs-zeta-correlations.md](uniform-gibbs-zeta-correlations.md). This is an exact correlation identity and interpretation of the already constructed kernel family; it is not an irrationality theorem or a historical novelty claim.

## Transition normalization and Hilbert-space factors

Integrating out the intermediate uniform T gives p(x,y)=-log(max(x,y))/(1-x). Multiplying by the stationary density 2(1-x) gives the symmetric joint density -2log(max(x,y)). Both total mass and stationarity follow either by direct integration or by the uniform triangle model. In the latter, T has density 2t and X,Y are conditionally independent uniforms on (0,t), so the claimed positive quadratic form is justified, including its normalization.

For each monomial, E[Y^j|T]=T^j/(j+1) and averaging T on (x,1) gives (1+...+x^j)/(j+1)^2. Thus P preserves polynomial degrees and scales the leading coefficient in degree k-1 by 1/k². Symmetry and the already proved orthogonality of Q_k make P Q_k=Q_k/k² an all-index consequence. Jensen and stationarity make P a contraction on L²((1-x)dx), not only on bounded or polynomial functions.

The accepted finite projection identity K=sum_(k<=n)Q_k/k+G_n has ||G_n||²=2 sum_(k>n)k^-3. This tends to zero, so it supplies the exact expansion needed for K without any assertion of a complete eigenbasis. Applying bounded P^r and taking the inner-product limit is valid. The contribution of one mode is (1/k)^2 k^(-2r)(2/k)=2k^(-2r-3). The normalized stationary measure contributes one additional factor 2, giving E[K(X0)K(Xr)]=4 zeta(2r+3). Its mean is 2, so subtracting 4 gives the stated covariance. For r=0 the same argument reduces to the checked mean-square integral. The unbounded logarithm lies in L², hence all correlations exist by Cauchy-Schwarz.

The finite-tail correlation formula follows by applying the same argument after removing the first n orthogonal modes. Positivity does not imply a favorable arithmetic normalization or transfer irrationality to another zeta value; the note correctly makes neither inference.

## Euler resolvent, and a direct simplification

The stated spectral proof of the resolvent is valid: the coefficient of Q_k is k/(k²-sigma), and the separately accepted Euler-kernel expansion has the same coefficient. The operator series converges in norm for 0<=sigma<1.

There is also a shorter proof using only the scalar ODE already established for the Euler integral. Put h(x)=(1-x)K_tau(x), with sigma=tau². Its ODE is exactly

\[
 (x h'(x))'=-\sigma K_\tau(x).
\]

The proved unit logarithmic flux implies lim_(x->0) x h'(x)=-1. Integrating once gives

\[
 xh'(x)=-1-\sigma\int_0^xK_\tau(y)dy.
\]

Since K_tau is finite at one, h(1)=0. A second integration therefore gives

\[
 h(x)=-\log x+\sigma\int_x^1\frac{dt}{t}\int_0^tK_\tau(y)dy.
\]

The two conditional-uniform integrations defining P are precisely

\[
 Pf(x)=\frac1{1-x}\int_x^1\frac{dt}{t}\int_0^tf(y)dy.
\]

Hence K_tau=K_0+sigma P K_tau pointwise on (0,1), and also in L²_w. The logarithmic endpoint bound justifies both integrations; for sigma>=0 Tonelli also applies directly. The contraction estimate now identifies K_tau=(I-sigma P)^(-1)K_0 without using a complete spectral system or an unbounded operator domain. This simplification preserves the exact normalization because the integration constant is the unit boundary flux.

Finally, direct differentiation of the same integral formula proves (J+1)Pf=f on polynomials, where J=-(1-x)^(-1)D[x(1-x)^2D]. It is an identity on explicitly given functions; the note does not overstate it as an unproved global-domain theorem.

The primary paper is used only for the classical Gibbs/orthogonal-polynomial framework and provenance discussion. No mathematical claim above depends on an unverified identification of this exact triangular sampler in that paper.
