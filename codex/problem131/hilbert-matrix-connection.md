# Apéry numbers measure the growth of inverse Hilbert matrices

Status: independently derived exact theorem below, based on classical Legendre
orthogonality. No historical novelty asserted. This is a concrete link to
numerical linear algebra, not a replacement for the arithmetic part of the
irrationality proof.

Let H_N be the N-by-N Hilbert matrix

    (H_N)_(i,j)=1/(i+j+1), 0<=i,j<N.

Let c_j be the coefficient vector of the shifted Legendre polynomial L_j,
padded with zeros to length N:

    (c_j)_k=(-1)^k binom(j,k)binom(j+k,k), 0<=k<=j.

The Apéry number is A_j=c_j^T c_j.

**Theorem.** For N>=1,

    H_N^(-1)=sum_(j=0)^(N-1) (2j+1)c_j c_j^T,
    trace(H_N^(-1))=sum_(j=0)^(N-1)(2j+1)A_j.

In particular the inverse Hilbert matrix has integer entries, and

    A_n=[trace(H_(n+1)^(-1))-trace(H_n^(-1))]/(2n+1)

for n>=1 (with trace(H_0^(-1)) interpreted as0 if n=0).

**Proof.** H_N is the Gram matrix of1,x,...,x^(N-1) in L²[0,1]. Rodrigues'
formula and integration by parts prove

    integral_0^1 L_j L_k=0 (j!=k),
    integral_0^1 L_j²=1/(2j+1).

For the norm, integrate the Rodrigues derivative j times onto L_j; its leading
coefficient is (-1)^j binom(2j,j), and the remaining beta integral is
(j!)²/(2j+1)!. These factors give the displayed norm.
Let C be the triangular matrix with rows c_j^T, and D=diag(1,3,...,2N-1).
Then C H_N C^T=D^(-1). Inverting and rearranging gives H_N^(-1)=C^T D C.
Taking traces yields the result. All c_j and diagonal entries of D are integers,
proving inverse integrality directly. QED.

## The same exponential constant controls numerical ill-conditioning

The independently established binomial asymptotic is

    A_n~C_A*(1+sqrt2)^(4n)*n^(-3/2), C_A>0.

Thus trace(H_N^(-1)) has exponential rate lambda=(1+sqrt2)^4. Let || ||_2
be the Euclidean operator norm. Since H_N^(-1) is positive definite,

    trace(H_N^(-1))/N <= ||H_N^(-1)||_2 <= trace(H_N^(-1)).

Also1<=||H_N||_2<=N: the lower bound uses its (0,0) entry, while the upper
bound follows from the Frobenius norm since every entry has absolute value<=1.
Hence the spectral condition number satisfies

    lim_(N->infinity) condition_2(H_N)^(1/N)=(1+sqrt2)^4.

No sharp Hilbert-matrix spectral theorem is needed for this exponential-rate
statement. Polynomial factors disappear under Nth roots.

## What the connection explains

The integer coefficients in Apéry's rational function are also the change of
basis that orthogonalizes a matrix of rational moments. Their Euclidean norm
measures how large the inverse moment problem becomes in the monomial basis.
The same coefficient growth that causes numerical instability supplies the
large integer solution needed for a small irrationality linear form through
the Casoratian.

This turns an apparent coincidence into an exact identity: the familiar
(1+sqrt2)^4 growth of Hilbert-matrix conditioning and Apéry's recurrence are
both encoded by the coefficient vectors of Legendre polynomials. It suggests
looking at other rational moment matrices, but their inverse integrality and
companion small forms must be proved separately.
