# Exact positive tails, a rational bracket, and the sharp optimized error rate

Status: complete corollaries of the independently accepted original
Schur/Apéry identity in `rank-one-original-apery-heun-rii-proof.md` and the
reviewed classical Apéry asymptotics in `apery-asymptotics.md`. Independent
audit of the normalizations is accepted in `rank-one-original-error-independent-review.md`. No historical novelty or new
arithmetic input is claimed.

Write

\[
 C_h=5A_h-A_{h-1},\quad D_h=5B_h-B_{h-1},\quad
 q_h=D_h/C_h,\quad E_h=q_h-\zeta(3)>0,
\]

where q_h is now the PROVED original mu_3/D_(2h) Schur root. The symbol D_h
in this note denotes the adjacent numerator, not the product polynomial used
in the measure. Also put d_h=12h²-15h+5.

## 1. Exact differences and a positive tail series

The RII proof gives L_h=h²C_h/[12(2h)!] and Y_h=q_hL_h satisfying the same
second-order recurrence V_(h+1)=a_hV_h-c_hV_(h-1), with

\[
 c_h=\frac{(h-1)d_{h+1}}{4h^2(h+1)(2h-1)(2h+1)d_h},\qquad h\ge2.
\]

Their Casoratian W_h=Y_hL_(h-1)-Y_(h-1)L_h consequently satisfies
W_(h+1)=c_hW_h. Since W_2=-23/96 and d_2=23, the factors telescope to

\[
 \boxed{W_h=-\frac{d_h}{2(h-1)(2h)!(2h-2)!}\quad(h\ge2).}   \tag{1}
\]

Division by L_hL_(h-1) gives the exact decreasing increments

\[
 \boxed{q_{h-1}-q_h=
 \frac{72d_h}{h^2(h-1)^3C_hC_{h-1}}>0\quad(h\ge2).}        \tag{2}
\]

The original variational construction already proved q_h->zeta(3), so
ordinary telescoping yields an exact positive error series:

\[
 \boxed{E_h=72\sum_{k=h+1}^{\infty}
   \frac{12k^2-15k+5}{k^2(k-1)^3C_kC_{k-1}}.}       \tag{3}
\]

This is the optimized original error, not the earlier nonoptimal trial norm.
The series follows from an exact finite telescoping identity plus the
established limit, so no unproved interchange is involved.

## 2. A simple rational two-sided bracket

The classical Apéry Casoratian is
A_(h-1)B_h-A_hB_(h-1)=6/h³. Therefore

\[
 \boxed{q_h-\frac{B_h}{A_h}=
              \frac{6}{h^3A_hC_h}.}               \tag{4}
\]

Together with the known positive lower error, this gives for every h>=1

\[
 \boxed{\frac{B_h}{A_h}<\zeta(3)<q_h,
  \qquad\text{bracket width}=\frac6{h^3A_hC_h}.}    \tag{5}
\]

At h=1 this is 6/5<zeta(3)<5/4 and the width is 1/20, fixing the index and
normalization. The original Schur construction thus supplies the upper end
of an entirely rational bracket whose width is known without evaluating zeta.

## 3. The sharp exponential rate and exact leading constant

Let

\[
 \lambda=17+12\sqrt2=(1+\sqrt2)^4,\qquad \rho=\lambda^{-1}.
\]

The reviewed asymptotics give

\[
 A_h\sim C_A\lambda^h h^{-3/2},\quad
 f_h:=A_h\zeta(3)-B_h\sim C_f\lambda^{-h}h^{-3/2},
\]
\[
 C_A=\frac{(1+\sqrt2)^2}{2^{9/4}\pi^{3/2}},\qquad
 C_f=\frac{\pi^{3/2}}{2^{1/4}(1+\sqrt2)^2}.
\]

In particular C_f/C_A=4pi³/lambda. The original upper error is exactly

\[
 E_h=\frac{f_{h-1}-5f_h}{5A_h-A_{h-1}}.
\]

The limiting coefficient is nonzero because lambda>5, so substitution gives

\[
 \boxed{\frac{E_h}{\zeta(3)-B_h/A_h}\longrightarrow
 \frac{\lambda-5}{5-\lambda^{-1}}
 =3+2\sqrt2=\sqrt\lambda.}                       \tag{6}
\]

Consequently the optimized norm has the full sharp asymptotic

\[
 \boxed{E_h\sim4\pi^3(3-2\sqrt2)\,(17-12\sqrt2)^{2h},}
                                                               \tag{7}
\]

and in particular

\[
 \boxed{\lim_{h\to\infty}E_h^{1/h}
                 =(17+12\sqrt2)^{-2}.}            \tag{8}
\]

There is no power of h in (7): the h^(-3/2) powers in the small linear form
and the large denominator cancel. The earlier h^(-3)16^(-h) statement was
an upper estimate from a particular trial polynomial, not this optimum.

The already reviewed comparison theorem in `rank-one-universal-error-rate.md`
now automatically identifies the same exponential rate for each fixed odd-zeta
rank-one family: its upper and lower comparison factors with E_h have only
polynomial growth/decay in h. Thus their h-th roots tend to the same value
in (8). This is a direct corollary of that existing comparison, not a new
higher-zeta construction; it supplies no leading constant or denominator
control for those other families.

## 4. What arithmetic transfers to the original family

The classical Apéry denominator theorem says A_h is integral and
lcm(1,...,h)^3 B_h is integral. Hence C_h is integral and the same lcm cube
clears D_h. Also

\[
 0<D_h-C_h\zeta(3)=f_{h-1}-5f_h<f_{h-1}.
\]

The existing elementary Apéry growth and lcm estimates make the cleared
positive quantity tend to zero. Under rationality of zeta(3), multiplication
by its fixed denominator would produce a positive integer tending to zero.
This transfers the established Apéry arithmetic to the ORIGINAL Schur
family. It is not an independent new proof of that arithmetic theorem.
