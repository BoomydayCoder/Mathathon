# Root independent audit: the original Schur–Apéry identity

Verdict: accepted. Reviewed proof:
[rank-one-original-apery-heun-rii-proof.md](rank-one-original-apery-heun-rii-proof.md).
This is the original positive measure mu_3/D_(2h), not the earlier moving-pole
substitute. The all-h quotient, constant term and endpoint consequences are
proved by the written argument. Historical priority remains unresolved.
A second extra-high agent independently accepted the proof in
[rank-one-heun-rii-second-independent-review.md](rank-one-heun-rii-second-independent-review.md).

## The interpolation input does not assume the desired identity

The original Gram matrix has G(X)=G(0)-Xee^T. If v is the monic minimizer
coefficient vector, its lower components satisfy the original orthogonality
and its final component is1. Thus G(q)v=0 at the original Schur root.
Passing to the two Cauchy polynomial bases gives the Loewner nullvector.
In the right basis its entries are Q_h(b^2)/C_R'(b^2): the signs from
p_h(-b^2)=(-1)^(h-1)Q_h(b^2) and from C'(-b^2) cancel. Their sum is1
by the monic leading coefficient. The Lagrange numerator then interpolates
f_j(q)=j^2q+beta_j at every node. Since Q_h has only negative roots,
none of the required positive-node evaluations vanish. Consequently
S_h=(N_h-qzQ_h)/Q_h(0) has degree<=h and leading coefficient -qL_h.
This rederives the necessary input directly from the original minimization;
no conjectural Apéry formula is embedded in it.

## Boundary-defect cancellation and the two cutoffs

The previously independently reviewed Wilson defect is a wedge of evaluations
at0 and1. Applying it to (Xf,g) and (f,Xg) gives the same wedge, because
X(0)=0 and X(1)=1. Hence {X,D}-D is genuinely symmetric; the additional
multiplication term is also symmetric. This checks the central Heun step.
The raising coefficient is lambda_k+lambda_(k+1)-kappa
=2(k-n)(k-n+1), with n=h-1. Both degree n and degree n-1 cutoffs vanish.
Thus the operator preserves both nested polynomial spaces. Orthogonality,
symmetry, and positive definiteness imply that the true minimizer is an
eigenpolynomial. The second cutoff ensures that a lower leading term cannot
alter its eigenvalue; there is no coefficient assumption at that step.
The Newton action gives eta_h=1-2h-2h^3(h+1), including eta_1=-5.

Direct expansion at zero yields Hf(0)=-2M f(1)+(1-M)f(0), so R_h(1)
has the stated value. In the Newton basis, the coefficient equation at phi_0
uses only e_0 and e_1; its lower coefficient is16(1-2h), nonzero for h>=1.
This gives the exact e_1/e_0 and hence R_h(4). The constant h=1 case has
e_1=0 and agrees with both endpoint formulas.

## The varying-measure recurrence and its boundary cases

Choose the linear interpolant at the TWO OLD nodes (2h-1)^2 and(2h)^2.
They are distinct and R_h is nonzero there. Dividing the polynomial difference
by their product gives W of degree<=h-2. The key measure conversion is

    ell_(2h-2)(W q)=ell_(2h)(difference*q).

The old factors are exactly those removed from the denominator. In the
R_(h+1) term, converting to ell_(2h+2) multiplies by the TWO NEW factors;
it does not divide by them. The degree of the test polynomial rises by2,
so k<=h-3 gives k+2<=h-1, within the next orthogonality range. The other
term needs only k+1<=h-2. Therefore W lies in the one-dimensional orthogonal
complement spanned by R_(h-1). For h=2 it is simply a constant and the
argument remains valid with no orthogonality tests. No undefined R_0 is used.

The endpoint system has determinant
3 h^6(h+1)^4(12h^2-15h+5)/8, strictly positive for h>=2. The displayed
three rational coefficients follow by solving it and are unique. Thus they
are not fitted to a finite list of original polynomials.

The interpolation-numerator recurrence has degree<=h+1. At the first2h-2
nodes all three interpolants apply; at the last two old nodes the third
term vanishes. There are2h distinct zeros, more than h+1 for every h>=2.
The h=2 case is four zeros against degree<=3. Hence the numerator recurrence
is identically zero, with the same coefficients and no missing forcing term.

## Scaling, initial data and conclusion

Leading coefficients of the polynomial and numerator give the same scalar
recurrence for L_h and q_h L_h. With s_h=h^2/[12(2h)!], the two factorial
ratios give exactly alpha_h and gamma_h in(13). Independently eliminating
u_(h+1) and u_(h-2) from Apéry's recurrence gives the same equation for
5u_h-u_(h-1). Its determinant -12h(12h^2-15h+5) is nonzero for h>=2.
Both original initial pairs h=1,2 agree with the separately known A,B pairs.
Uniqueness proves L=s(5A-A_previous) and qL=s(5B-B_previous). Dividing is
legitimate because L>0. This proves the original quotient theorem for all h.
The constant-term formula follows by inversion, and substituting it and the
endpoint ratio into the already proved Riccati normalization gives the claimed
72(h+1)/[h(5A_h-A_(h-1))^2] exactly.

## Independent reproduction

`check_original_schur_apery_independent.py` uses fresh standard-library Fraction
Gaussian elimination on original rational pole moments; it imports none of
the author's Gram generators or recurrence data. Apéry A,B are independently
computed by binomial and harmonic partial-fraction sums. All10 original
minimizers, both endpoint formulas, quotient identities and8 full polynomial
RII relations pass. A separate symbolic solve of the endpoint linear system
verifies its determinant and the free-index Apéry recurrence transformation.
These finite checks supplement the all-h proof audited above.
