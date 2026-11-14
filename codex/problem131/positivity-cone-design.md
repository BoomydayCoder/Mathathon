# The sharp positivity cone for quartic fifth-zeta weights

Status: proved cone and duality certificates; exact rational verification for n=1,...,5. No all-n sign theorem or irrationality result is claimed. This note explains why unrestricted coefficient freedom does not automatically give useful approximations.

## 1. Domain and notation

Use H_n(t) and v=t+n/2 from weighted-family-extension.md, and weights

    (2t+n) q(v²) H_n(t),    degree q<=2.

Here impose the stronger continuous positivity condition

    q(z)>=0 for z>=L_n=9n²/4.

This corresponds to real t>=n, and suffices for positive summands at integer t>=n+1. The earlier note used the weaker half-line z>=(3n/2+1)²; its numerical cone endpoints therefore differ. Both choices are valid. The formulas below work for any rational L; the tables use L=9n²/4 as specified.

Shift y=z-L and write

    q(z)=a y²+b y+c.

The summed rational function has coefficient functionals

    S(q)=Z(q) zeta(5)+ell(q) zeta(3)-R(q),
    ell(q)=ell_a a+ell_b b+ell_c c,

with Z,R also rational linear functionals. Imposing ell(q)=0 isolates zeta(5). All functionals are explicitly computed by the principal-part formulas already proved in the companion notes.

## 2. Exact 2x2 semidefinite representation

Let K be the cone of real quadratics nonnegative for y>=0. Then

    q in K iff q(y)=[y,1] M [y,1]^T+d y,
    M positive semidefinite, d>=0.                       (1)

To prove this, q>=0 on the half-line is equivalent to a,c>=0 and b>=-2sqrt(ac). If b>=0, choose M=diag(a,c), d=b. If b<0, choose M=[[a,b/2],[b/2,c]], d=0; its determinant is nonnegative. The reverse implication of (1) is immediate.

For a linear functional f(q)=f_a a+f_b b+f_c c, nonnegativity on K is equivalent to

    [[f_a,f_b],[f_b,f_c]] positive semidefinite,
    f_b>=0.                                             (2)

Indeed pairing with (1) gives the matrix trace pairing plus d f_b. Necessity follows by testing every positive semidefinite M and every d>=0. Thus K* is described by one 2x2 semidefinite constraint and one scalar inequality.

A useful dual certificate for a rational or algebraic candidate slope rho is: find eta such that

    f=rho Z-R-eta ell satisfies (2).                     (3)

Then for every admissible positive q with ell(q)=0,

    rho Z(q)-R(q)>=0.

For Z(q)>0 this implies R/Z<=rho; for Z(q)<0 it implies R/Z>=rho. This is an explicit sufficient certificate, valid without a duality theorem. Optimizing R on the normalized slices Z=1 or Z=-1, ell=0, q in K is a 2x2 semidefinite optimization problem. Under the usual strict-feasibility hypothesis, conic duality makes these certificates sharp; the concrete family below is solved directly, so it does not rely on an unproved strong-duality assertion.

## 3. Closed boundary equations and generalized eigenvalues

For a nondegenerate two-dimensional section ell(q)=0, boundary rays can arise in three ways:

1. a=0, with b,c>=0;
2. c=0, with a,b>=0;
3. q(y)=(y-r)² up to a positive factor, r>=0.

The first two give rational boundary rays by one linear equation. The third requires

    ell_c r²-2ell_b r+ell_a=0.                            (4)

Hence the only potentially irrational ray parameters are roots of one quadratic. Retain only real nonnegative r, and discard a ray with Z(q)=0 when assigning a finite slope. Degenerate sections are handled directly by the semidefinite representation (1).

There is a compact generalized-eigenvalue form for the square boundary. Assume ell_c!=0 and define

    C=[[0,-ell_a/ell_c],[1,2ell_b/ell_c]].

This is multiplication by r in the quadratic algebra given by (4), on the basis 1,r. Reduce the quadratic expressions R(1,-2r,r²) and Z(1,-2r,r²) modulo (4):

    R_0=r_a-r_c ell_a/ell_c,
    R_1=-2r_b+2r_c ell_b/ell_c,
    Z_0=z_a-z_c ell_a/ell_c,
    Z_1=-2z_b+2z_c ell_b/ell_c.

Every finite square-boundary slope rho satisfies

    det[(R_0 I+R_1 C)-rho(Z_0 I+Z_1 C)]=0.               (5)

This is a 2x2 generalized eigenvalue condition with rational entries. Its algebraic roots must still be checked against the admissibility conditions on r. The rational a=0 and c=0 boundaries must also be included; equation (5) is not a substitute for those cases.

## 4. The actual first five cases simplify to a rational sector

The exact coefficient computation for n=1,...,5 gives

    ell_a ell_c>0,    ell_a ell_b<0.

Set

    u=-ell_a/ell_b>0,    v=-ell_c/ell_b>0.

The cancellation equation becomes

    b=u a+v c.                                          (6)

Nonnegativity on the half-line necessarily gives a,c>=0. Conversely, these two inequalities and (6) give b>=0, making every coefficient of a y²+b y+c nonnegative. Therefore the positive canceled cone is *exactly*

    q=a q_1+c q_2,   a,c>=0,
    q_1=y²+u y,    q_2=1+v y.                            (7)

This is a complete cone description for each verified n, not merely a sufficient inner approximation. The curved square boundary is absent: ell(1,-2r,r²) has a fixed nonzero sign for r>=0 under the displayed sign pattern.

The coefficient computation also gives

    Z_1=Z(q_1)>0,    Z_2=Z(q_2)<0.

Both sums S(q_1),S(q_2) are strictly positive, because their summands are nonnegative and infinitely many are positive. Define

    r_-=R(q_1)/Z_1,    r_+=R(q_2)/Z_2.

Then positivity proves

    r_- < zeta(5) < r_+.                                (8)

Moreover the complete set of slopes R/Z from this cone is

    (-infinity,r_-] union [r_+,infinity),

with endpoints attained by q_1 and q_2. To verify sharpness, substitute (7): when Z=aZ_1+cZ_2>0, the ratio is at most r_-; when Z<0 it is at least r_+. Letting Z approach zero from either side gives the infinite ends of the intervals.

Thus the positive cone excludes a sharp open interval containing the target, even though unrestricted rational weights can realize every rational coefficient pair. This is the missing distinction between algebraic representability and a positive approximation construction.

## 5. The n=1 case in a few exact numbers

At L=9/4, the cancellation functional is

    ell=-880a+370b-156c,

so

    b=(88/37)a+(78/185)c.

Its positive canceled cone is generated by

    q_1=y²+(88/37)y,    q_2=1+(78/185)y.

Their coefficient ratios give the exact sharp gap

    659/640 < zeta(5) < 557/536,

and every positive canceled weight has slope at most 659/640 or at least 557/536. No numerical value of zeta(5) is used: strict inequalities follow from the corresponding positive series.

For the weaker domain L=25/4, the lower endpoint becomes 23063/22272, as in weighted-family-extension.md. The different lower endpoints reflect different imposed positivity domains, not a contradiction.

## 6. Exact rational intervals for n=2,...,5

| n | r_- | r_+ | numerical width |
|---|---|---|---|
| 1 | 659/640 | 557/536 | 9.49160e-3 |
| 2 | 1127751817/1087590144 | 62525737/60299008 | 9.06857e-7 |
| 3 | 1006821768888589/970966168047360 | 669057626547713/645230705022720 | 1.27196e-10 |
| 4 | 12229035509065451320139/11793527030602742538240 | 2030658482618565161729/1958341333372622069760 | 1.81686e-14 |
| 5 | 356315607505039857895754777/343626261075222339456000000 | 8281273639772626474979412613/7986355460827281327283200000 | 2.65053e-18 |

Each interval is certified by exact rational coefficients, the verified signs in section4, and strict positivity of the two infinite series. The decimal width is only a summary. Narrow intervals alone do not prove irrationality: the denominators are also growing rapidly, and the relevant integer linear forms require their own arithmetic estimates.

## 7. An explicit rational lower bound for every primitive form in each cone

Let t_0=n+1, y_0=(t_0+n/2)²-L, and w_0=(2t_0+n)H_n(t_0)>0. Since each endpoint sum dominates its first positive term,

    zeta(5)-r_- >= w_0 q_1(y_0)/Z_1,
    r_+-zeta(5) >= w_0 q_2(y_0)/(-Z_2).

The right sides are positive rational numbers. Let delta_n be their minimum. For every q in the cone, the decomposition (7) gives

    S(q)>=delta_n [aZ_1+c(-Z_2)]>=delta_n |Z(q)|.

After primitive integer normalization with U>0,

    |U zeta(5)-V|>=delta_n U>=delta_n.                   (9)

This proves a coefficient-uniform positive gap at each fixed n. It does not give a uniform-in-n obstruction: delta_n becomes small as n grows. The exact rational delta_n values are saved in the computational certificates.

## 7a. Discrete positivity introduces no extra freedom in these five cases

The actual sum only needs q(z_t)>=0 at z_t=(t+n/2)², integer t>=n+1. In general, a quadratic could be negative between consecutive sampled points while remaining nonnegative at every sample. Here the cancellation condition rules out that extra freedom.

Take L_*=(3n/2+1)², so y=0 is the first summation node. The same sign pattern ell_a ell_c>0, ell_a ell_b<0 is verified exactly for n=1,...,5. If q is nonnegative at all summation nodes, then a>=0 (from the unbounded sequence of nodes) and c=q(L_*)>=0. Cancellation forces b=u a+v c>=0. Thus q is automatically nonnegative on the entire half-line y>=0.

Therefore, for these five n, the discrete positivity cone equals the continuous cone for L=L_*. The sharp upper endpoint remains unchanged, and the sharp lower endpoints are

| n | sharp lower endpoint for actual summandwise positivity |
|---|---|
| 1 | 23063/22272 |
| 2 | 53042811667/51153825024 |
| 3 | 43073743711003271/41539773140951040 |
| 4 | 68531629954079712478393/66091036346697600000000 |
| 5 | 5277853143068723705136714343/5089894755820270113715200000 |

The q_1 extreme ray now vanishes at the first summation node. Its next term, t=n+2, supplies a strictly positive rational gap certificate; q_2 remains positive. These additional exact certificates are saved under `discrete_domain` in the JSON. This closes the potential loophole of allowing narrow negative intervals between the summation nodes, within the verified cases.

## 8. Reproduction and limitations

Run `python codex/computations/check_positivity_cone.py`. It imports the exact coefficient routines from weighted_family_extension.py and writes check_positivity_cone.json containing L, the cancellation row, u,v, both extreme polynomial weights, their rational coefficient pairs, exact interval endpoints, and first-term rational gap certificates.

The sign pattern and opposite Z signs have been checked exactly for n=1,...,5. A general proof for all n remains open in this investigation. The semidefinite and generalized-eigenvalue formulations are general; the particularly simple rational sector is presently a conditional theorem with these five exact instances.

The next task is arithmetic, not further fixed-n coefficient tuning: derive how the primitive endpoint denominators and the rational gap delta_n scale together. The cone analysis explains why a large coefficient search cannot bypass the arithmetic cost simply by exploiting the extra dimension.
