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


## Sampler determinant supplement

Section 5 is also accepted. The finite determinant follows directly from the triangular monomial action, with diagonal entries 1/k² for k=1,...,N. Completeness of the orthogonal polynomial modes is justified by polynomial density for the finite weight measure. This identifies the whole bounded operator, not merely its restriction to a possibly proper invariant subspace. Its positive eigenvalues have finite sum, so the trace-class assertion and the absolutely convergent determinant product are valid. Euler's sine product gives the displayed entire expression, with no choice of square-root branch in its power series.

The identity for A_n is just substitution of the finite product into its already established Newton weights. The two uniform limits retain exactly the closed disk |sigma|<=1/4 from the separately accepted asymptotic theorem. In particular the determinant interpretation does not justify extending those error asymptotics to a zero or pole elsewhere, and the note makes no such extension. The new connection adds no denominator-growth or irrationality assertion.


The added dilogarithm illustration in Section 3 also has the stated normalization. Under the stationary triangle law, the latent density is 2t, so the unnormalized pairing is integral t Af(t)Ag(t)dt. Differentiating zeta(2)-Li_2(1-t) gives -log(t)/(1-t), and its value at t=0 is zero. Thus the integral of its square divided by t is exactly <K,PK>_w=2 zeta(5). Near zero the primitive is O(t(1+|log t|)), and near one it is bounded, so both endpoints converge. This is a consequence of the already reviewed correlation formula.


## Harmonic-truncation limitation and the second Legendre factor

Section 6 is accepted for each fixed integer r>=0, s=2r+3, and n>=1. The earlier exact tail correlation has the factor 2, so one half of <G_n,P^rG_n>_w is precisely the harmonic tail. Since x^-s is decreasing, summing its values at integers k>=n+1 dominates its integral from n+1 to infinity, giving the stated lower bound 1/[(s-1)(n+1)^(s-1)]. The companion upper integral bound from n to infinity also confirms its polynomial order.

The multiplier d_n^s really does clear the harmonic sum because every k<=n divides d_n. The accepted elementary bound d_n>=2^(n-1)/n, raised to the fixed power s and multiplied by the positive tail lower bound, gives exactly the displayed exponential-over-polynomial lower bound. It tends to positive infinity. This conclusion concerns this specified clearing multiplier; it neither lower-bounds the reduced denominator nor excludes further cancellation under a smaller one.

The final comparison preserves the distinct normalizations: I_(n,0)=2(zeta(3)-H_n^(3)) is the boundary row, while I_(n,n)=2(A_n zeta(3)-B_n) is the diagonal linear form. The accepted energy estimate and the exponential lower growth of A_n make the latter exponentially small. Thus the explanation of why the second Legendre factor matters is accurate and does not claim that positivity alone proves the improved rate or an irrationality theorem.
