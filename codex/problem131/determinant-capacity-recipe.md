# From a positive integral to a determinant: an elementary capacity criterion

Status: independently reconstructed general theorem, proved below. No novelty
claim. Motivated by Fauzan's positive Hankel construction, but the fixed compact
support hypotheses here differ substantially from that paper's n-dependent
rational weight on an unbounded domain. This is a design criterion, not a proof
for a new zeta value.

## The theorem

Let mu be a finite positive measure of mass M>0 with infinite support contained
in [a,b], b>a. Suppose for every k>=0 its moment has the form

    integral x^k dmu(x)=u_k*alpha+v_k,   u_k,v_k rational.

For each h choose a positive integer d_h clearing all u_k,v_k for
0<=k<=2h-2. If

    limsup_(h->infinity) log(d_h)/h < log(4/(b-a)),

then alpha is irrational.

The hypothesis in particular forces b-a<4, since d_h>=1. Rescaling a measure
may shrink its support but also changes moment denominators. Both effects must
be counted; changing units cannot manufacture irrationality.

## Proof by monic polynomial norms

Define G_h(X)=(u_(i+j)X+v_(i+j))_(0<=i,j<h), and P_h(X)=d_h^h det G_h(X).
Then P_h is an integer polynomial of degree at most h. At X=alpha the matrix
is the Gram matrix of 1,x,...,x^(h-1). Its quadratic form is integral q(x)^2 dmu.
A nonzero polynomial has finitely many roots, so infinite support makes this
strictly positive. Thus P_h(alpha)>0.

Gram-Schmidt with monic polynomials p_j changes basis by a unit triangular
matrix. Therefore

    det G_h(alpha)=product_(j=0)^(h-1) integral p_j(x)^2 dmu(x).

Each monic orthogonal p_j minimizes its squared norm among monic degree-j
polynomials: write any competitor as p_j+r, deg r<j, and use orthogonality.
For j>=1 use the rescaled Chebyshev competitor

    q_j(x)=2*((b-a)/4)^j*T_j((2x-a-b)/(b-a)).

Its leading coefficient is one because T_j has leading coefficient 2^(j-1),
and its absolute value is at most 2*((b-a)/4)^j on [a,b]. Since p_0=1,

    0<det G_h(alpha)<=M^h*4^(h-1)*((b-a)/4)^(h(h-1)).

Consequently

    log P_h(alpha)
       <=h log d_h+h log M+(h-1)log4+h(h-1)log((b-a)/4)
       <=-c*h²+O(h)

for some c>0 and all sufficiently large h. If alpha=A/B with B>=1 integer,
then B^h P_h(A/B) is a positive integer; its logarithm tends to minus infinity,
a contradiction. QED.

## Why this is a different approximation mechanism

A linear-form proof spends all its effort making one number small. A Gram
determinant multiplies h orthogonal residual norms. The jth monic norm can be
exponentially small in j, so the determinant is exponentially small in h².
Its degree in the target constant is only h. Under rationality, clearing the
unknown denominator costs B^h, leaving the quadratic decay decisive.

The factor (b-a)/4 is the logarithmic capacity of an interval, but the proof
above needs only Chebyshev polynomials and Gram-Schmidt. This gives a concrete
bridge from elementary approximation to the logarithmic-energy viewpoint in
more sophisticated determinant proofs.

## Exact example: why the arithmetic cost cannot be ignored

For Lebesgue measure on [0,1], moments are 1/(k+1), rational. Take u_k=0,
v_k=1/(k+1). Then d_h must be divisible by lcm(1,...,2h-1), with exponential
rate 2 by the prime number theorem, greater than log4. The sufficient criterion
fails, as it must: this construction contains no irrational constant at all.
This comparison uses PNT only to illustrate the asymptotic failure, not in the
theorem's proof. Exact determinants are positive rationals (Hilbert matrices),
showing that positivity and determinant smallness alone never suffice.

For zeta moment constructions the measure generally depends on h, so M, support,
and the rational functions used as the basis also vary. Those changes require
new estimates. The argument identifies what must be preserved; it does not
silently transfer the compact fixed-measure theorem to Fauzan's setting.

## A broader version, with the same proof

The support interval may depend on h if the moment data still define an h-by-h
positive Gram matrix and suitable monic competitors exist. It is enough to
exhibit bounds N_(h,j) on their squared norms such that

    h log d_h+sum_(j=0)^(h-1) log N_(h,j) <=-c*h²

for c>0. Tracking these norms directly is often more useful than quoting
capacity: it makes the analytic-versus-arithmetic budget explicit and suggests
how to optimize the auxiliary weight.

## The Hilbert example gives an exact arithmetic boundary

There is an even sharper lesson requiring no prime number theorem. For Lebesgue
measure on [0,1], the h-by-h Gram matrix is the Hilbert matrix
H_(i,j)=1/(i+j+1). Shifted Legendre polynomials have leading coefficient
(-1)^j binom(2j,j) and squared norm 1/(2j+1). The norm follows by j integrations
by parts in the Rodrigues formula, followed by the beta integral.

Thus the monic orthogonal norms give the exact identity

    det H=product_(j=0)^(h-1) 1/[(2j+1)*binom(2j,j)^2].

Its reciprocal is an integer. The exact least positive integer multiplier
turning this positive rational determinant into an integer therefore turns it
into 1, never into a positive integer below 1. The analytic smallness and the
true denominator cost match exactly. This is a clean example of how a strong
positive determinant estimate can be mathematically real yet arithmetically
useless for irrationality.
