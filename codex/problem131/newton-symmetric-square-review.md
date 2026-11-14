# Independent review and a short closure of the Apéry diagonal

**Status.** The two-parameter symmetric-square theorem and the Newton/bispectral identities were independently checked. No error found. A further exact adjacent-diagonal identity below supplies a direct elementary derivation of the constants17 and5 from the bispectral equation. It does not assume the Apéry recurrence, and it makes no historical-novelty claim.

Reviewed notes: `symmetric-square-design-family.md` and `newton-bispectral-motivation.md`.

## 1. Symmetric-square family: accepted with the stated domains

For Delta=1-2at+t^2, the proposed coefficients are

    p=1/t+Delta'/(2Delta), q=(t-2b)/(4tDelta).

Direct expansion of the generating differential equation gives precisely

    y'''+3p y''+(p'+2p^2+4q)y'+(2q'+4pq)y=0.

Differentiating y=u^2 and substituting u''=-pu'-qu proves this identity for the square of any solution. At zero the second-order indicial factor is n^2, so u(0)=1 selects one formal series and forces u'(0)=b/2. Its square has the same constant and linear terms as the unique third-order normalized formal series; the n^3 recurrence then proves equality.

The integrating factor t sqrt(Delta) is correct. The singular exponents0,0 at zero,0,1/2 at either simple root of Delta, and1/2,1/2 at infinity follow from the displayed p,q. The exclusions a=+/-1 in the singularity table are necessary, though not in the rational identity. A branch of sqrt(Delta) must be specified for the self-adjoint analytic form; its normalized formal branch always exists.

Thus the symmetric-square property alone imposes no restriction on a,b and cannot select17,5 or imply integrality. The examples(1,1) and(2,1) check exactly.

## 2. Newton operator: accepted with boundary qualifications

The product identity

    phi_k(n)=product_(j=0)^(k-1)[n(n+1)-j(j+1)]/(k!)^2

gives its symmetry under n -> -n-1 and the raising identity

    (k+1)^2 phi_(k+1)=[n(n+1)-k(k+1)]phi_k.

Inserting the two contiguous n-ratios into L and cancelling gives

    L phi_k=k(k+1)phi_k+k^2 phi_(k-1).

Although the ratios have apparent poles at particular integer n, the resulting identity is polynomial after cancellation. At n=0, the backward edge has coefficient0, so no value of an arbitrary function at-1 is required. For k=0 the lowering term is0, not an actual use of phi_-1.

When proving the eigenfunction identity at fixed m, write K(n,m)=sum_(k=0)^m phi_k(n)phi_k(m). This is a single polynomial in n; it agrees with the min(n,m) truncation for nonnegative integers. Applying L to a truncated sum whose upper limit is mistakenly held fixed at its value at the current n would miss a boundary term. The polynomial interpretation avoids that error.

The weighted summation-by-parts formula also holds exactly. Its nonpositive quadratic form is a statement on finitely supported functions (or an appropriate completed domain), whereas these polynomial eigenfunctions are not square-summable with weight2n+1. Positive eigenvalues on those polynomials do not violate self-adjoint negativity on that Hilbert-space domain.

## 3. The missing finite-grid identity

For every integer n>=1,

    A_n+A_(n-1)=6 K(n,n-1),                         (10)

where A_n=K(n,n). This identity closes the first strip adjacent to the diagonal. It can be proved directly from the Newton factors, without any known recurrence for A_n.

Set f_k=phi_k(n), g_k=phi_k(n-1), for0<=k<=n. Then

    g_k=(n-k)/(n+k) f_k,
    f_(k+1)/f_k=(n-k)(n+k+1)/(k+1)^2.

The first formula remains true at k=n, where g_n=0. Define

    T_k=k^4 f_k^2/[n^2(n+k)^2], 0<=k<=n+1,

using f_(n+1)=0. Thus T_0=T_(n+1)=0. For0<=k<=n the raising ratio gives

    T_(k+1)=(n-k)^2 f_k^2/n^2.

Consequently

    4(T_k-T_(k+1))
      =4[2k^2-n^2]f_k^2/(n+k)^2
      =f_k^2+g_k^2-6 f_k g_k.                      (11)

All denominators in this argument are nonzero for n>=1. Summing(11) from k=0 to n gives(10), including the endpoint where the shorter sum terminates.

A boundary-safe binomial form is

    T_k=[binom(n-1,k-1)binom(n+k-1,k-1)]^2,

with T_0=0 and the usual vanishing binomial convention at k=n+1. This also makes clear that the telescoping quantities themselves are integers.

### Why six is natural here

For an unknown constant c,

    f_k^2+g_k^2-c f_k g_k
      =[(2-c)n^2+(2+c)k^2]f_k^2/(n+k)^2.

The most direct telescoping ansatz uses k^4 to cancel the fourth-power factorial ratio in f_(k+1)^2/f_k^2, and(n+k)^2 to cancel its shifted numerator. Its difference is proportional to2k^2-n^2, because

    k^4-(n^2-k^2)^2=n^2(2k^2-n^2).

Matching these two quadratic coefficients forces2+c=-2(2-c), hence c=6. Thus six emerges from a simple difference of squares, rather than from guessing the final Apéry recurrence.

### A finite weighted-average interpretation

The parent investigation supplied a useful interpretation of the same telescoper. For n>=1 put

    c_k=n phi_k(n)/(n+k)=binom(n,k)binom(n+k-1,k),
    w_k=c_k^2, 0<=k<=n, and w_(n+1)=0.

The contiguous ratio says

    w_k(n^2-k^2)^2=(k+1)^4 w_(k+1).

Summing and shifting the finite index gives

    sum_(k=0)^n w_k(n^2-k^2)^2=sum_(k=0)^n k^4 w_k.

The k^4 terms cancel after expansion. Therefore the exact weighted second moment is

    sum k^2 w_k=(n^2/2)sum w_k.                    (13)

Writing phi_k(n)=c_k(1+k/n) and phi_k(n-1)=c_k(1-k/n), equation(13) immediately gives

    A_n+A_(n-1)=3sum w_k,
    K(n,n-1)=(1/2)sum w_k,

and hence(10). In probability language, the displayed edge relation is detailed balance for a finite birth-death chain with upward rates(n^2-k^2)^2 and downward rates k^4, with rates at the two outside edges equal to zero. No Markov-chain theorem is needed: the finite index shift is the entire proof.

The same ratio also locates the maximal weights without calculus. The inequality w_(k+1)>=w_k is equivalent to

    2k^2+2k+1<=n^2.

Thus the mode lies at n/sqrt(2)+O(1), and the exact root-mean-square index from(13) is n/sqrt(2). This explains the appearance of that location in a saddle analysis. Neither the second moment alone nor the mode calculation is being asserted as a full concentration or asymptotic theorem.

## 4. Deriving17 and5, with no illicit diagonal substitution

The bispectral equation, evaluated at the genuine two-variable point(n,n), says

    (n+1)^3 K(n+1,n)+n^3 K(n-1,n)
       =(2n+1)[3n(n+1)+1]A_n.                     (12)

By(10) and symmetry, its two off-diagonal terms are

    K(n+1,n)=(A_(n+1)+A_n)/6,
    K(n-1,n)=(A_n+A_(n-1))/6.

Substituting these actual identities, rather than replacing the entries by A_(n+/-1), gives

    (n+1)^3 A_(n+1)
      ={6(2n+1)[3n(n+1)+1]-(n+1)^3-n^3}A_n
         -n^3 A_(n-1).

Finally

    (n+1)^3+n^3=(2n+1)[n(n+1)+1],

so the middle coefficient is exactly

    (2n+1)[17n(n+1)+5].

This proves the Apéry recurrence for n>=1. At n=0, K(1,0)=1 and(10) at n=1 give A_1+A_0=6; since A_0=1, A_1=5. Equivalently the defining sum gives A_1=1+4=5. There is no undefined A_-1 requirement.

This argument identifies a concrete route to the characteristic integers:

    quadratic Newton grid -> cubic difference operator
       -> adjacent-diagonal closure with coefficient6
       ->17=6*3-1 and5=6*1-1.

The two-variable identity alone did not supply the closure; the elementary telescoping relation(11) is the additional ingredient. Unlike the generic symmetric-square condition, these identities use the actual integer-residue/binomial structure of the Apéry sequence.

## 5. Exact verification and limits of the conclusion

The supplied symbolic scripts for both original notes passed independently. The new script `check_newton_diagonal_closure.py` proves the rational summand identity symbolically and checks boundary conventions, adjacent diagonals and the recurrence by exact integer arithmetic. These finite checks supplement the general argument(11)–(12).

This explains the homogeneous recurrence for A_n. Identifying the companion solution's limiting constant as zeta(3), establishing its denominator bound, and proving sufficiently small nonzero linear forms remain separate steps of the irrationality proof. The closure identity does not extend those claims to arbitrary a,b, and no such extension is asserted.

## 6. The fixed-m symmetric square is explicitly Legendre

There is also an exact link between the two reviewed notes. Write P_m for the usual Legendre polynomial normalized by P_m(1)=1. Then

    sum_(n>=0) K(n,m)t^n
       =(1-t)^-1 P_m((1+t)/(1-t))^2
       =[sum_(j=0)^m binom(m,j)^2 t^j]^2/(1-t)^(2m+1).       (14)

Here is a proof from the established differential family, not a numerical guess. Put z=(1+t)/(1-t), u=(1-t)^-1/2 v(z), and b=2m(m+1)+1. Substitution into the a=1 second-order operator gives exactly

    u''+p u'+q u
      =-[ (1-z^2)v''-2zv'+m(m+1)v ]/[t(1-t)^(5/2)].

Thus v=P_m makes u a solution; its formal value at zero is1. Uniqueness and the symmetric-square identity prove the first equality in(14). Rodrigues' formula and Leibniz's rule give

    P_m(z)=sum_j binom(m,j)^2 ((z-1)/2)^j ((z+1)/2)^(m-j),

which proves the second equality. The square root is the normalized formal branch, or the analytic branch near t=0.

Formula(14) gives another elementary explanation of the fixed-m family's integrality and polynomial growth: its generating function has a single pole of order2m+1 at1 and a polynomial numerator with positive integer coefficients. It also ties the quadratic Newton basis to ordinary Legendre orthogonality. It does not give the diagonal generating function by putting m=n, since m is fixed throughout the coefficient extraction.
