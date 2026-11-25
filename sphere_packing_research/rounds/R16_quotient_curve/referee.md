# R16 referee report: Proposition "quotient-curve dictionary" (report/sec_p61.tex, prop:quotient, lines 38-53)

VERDICT: CORRECT_MODULO_MINOR_FIXES  (no item refuted; three statements are used without proof, all provable)

Conventions (rounds/R11_mordell_weil/T1/notes.md, M2-M4): E_0/F_p supersingular, pi_0^2=-p, M_p=Hom_{Fbar_p}(J,E_0),
b(phi,psi)=deg(phi+psi)-deg phi-deg psi, pi_J = p-Frobenius of J=Jac C_p (pi_J^2=p), q_F=(u+v sqrt p)/2, u=deg phi,
v=b(phi pi_J,phi)/(2p), beta = O_F-bilinear form with Tr beta=b, Q_t=e^t b'_+ + e^-t b'_-, b'=Tr(lambda beta), lambda=eps/sqrt p.

## (a) item (0): CORRECT
* diag(a,a^-1): x->a^2 x, y->y/(a^-1)^((p+1)/2)=a^((p+1)/2) y; (a^((p+1)/2))^2=a^(p+1)=a^2, a^(2p)=a^2. OK.
* m odd => m | (p-1)/2; Y -> a^((p-1)/2+m) Y = Y. Y^2=x^(m-1)(x^p-x)=X^(N+1)-X. OK.
* C_m acts faithfully (a^2=1, a of odd order => a=1), k(X,Y) is invariant and [k(x,y):k(X,Y)]<=m, so k(X,Y) is the
  whole invariant field and deg pi = m exactly.
* X^(N+1)-X = X(X^N-1) is separable (N | p-1), N is even (in fact 4 | N since p=1 mod 4, m odd), so N+1 is odd, the
  affine model is smooth with one point at infinity, genus N/2=(p-1)/(2m); g_D is EVEN and >= 2 (g_D=1 never occurs).
* Riemann-Hurwitz: p-3 = m(N-2) + 2(m-1); exactly two branch points (0,0) and infinity, both totally ramified; no other
  point has a non-trivial stabiliser. OK.

## (b) item (1), bijectivity: CORRECT (at the level of lattices, index 1)
* G acts on M_p by psi -> psi o h_* (or its inverse; same invariants). With iota=iota_{P_0}, phi=psi o iota:
  phi(hP)=psi(h_*[P-P_0])+psi([hP_0-P_0]); hP_0=P_0 kills the constant, so "psi invariant" <=> "phi o h=phi"
  (converse because iota(C) generates J). The fixed-point argument is valid.
* phi factors uniquely through the categorical quotient C/C_m=D: phi=phi' o pi, phi'(pi P_0)=0, and psi=psi' o pi_*
  (pi_*: Jac C -> Jac D Albanese push-forward; this is what "pi^*" must mean). Injective because pi_* is surjective.
  The proof text omits the (trivial) injectivity and the definition of pi^*; add one line.
* Without a fixed point the statement would be false in general (index = order of the unramified abelian part of pi);
  here pi is totally ramified so pi^*: Jac D -> Jac C is injective and ker pi_* is connected. Consistent.

## (c) Fix(C_m): norm image versus full invariants -- TRUE BUT NOT PROVED IN THE TEXT (minor gap, line 45/51)
The proposition is proved for the full invariant lattice M_p^{C_m}. The scans build "Fix(C_m)" from norm images of
samples (cf. the p=541 index-2 incident in sec:r13verify). The two agree, but this needs an argument:
Norm = sum_h h_* = pi^* pi_* on Jac C, so Norm(M_p) = pi^*( Hom(Jac C,E_0) o pi^* ), and restriction
Hom(Jac C,E_0) -> Hom(Jac D,E_0) along the closed immersion pi^*: Jac D -> Jac C is surjective because Jac C is
superspecial (End(E_0) maximal => hereditary => every abelian subvariety with connected complement is a direct factor:
Jac C = Jac D x K).  Hence H^0-hat(C_m, M_p)=0 and Norm(M_p)=M_p^{C_m}.  Add this sentence.
(Rank check: trace of h on H^1 is 2-#Fix = 0, so M_p (x) Q is Q[C_m]-free and rank Fix = 2(p-1)/m = 4 g_D.)

## (d) twisted forms: CORRECT
pi is defined over F_p, so pi_* pi_{J_C} = pi_{J_D} pi_*; hence pi_{J_D}^2 = p (pi_* surjective), pi^* is O_F-linear,
deg(psi' o pi_*) = psi' pi_* pi^* psi'^vee = m deg psi', b_C(pi^* . , pi^* .) = m b_D, v = m v', q_F(pi^* phi') = m q'_F(phi'),
beta_C = m beta_D, and Q_t = m Q'_t for every t (same lambda, same t).  Only pi_J of the curve enters v, so the fact that
E_0 has Frob_{p^2} = -p while C_p has +p (T1, M2(ii)) is irrelevant as long as Hom is taken over Fbar_p on both sides
and the SAME E_0 is used on both sides.

## (e) integrality of q' (item (3), line 47 "q' 이 0 이 아닌 O_F 의 원소이면"; proof line 51): HYPOTHESIS HOLDS, NOT JUSTIFIED IN THE TEXT
From (d) alone one only gets q' in (1/m) O_F.  Integrality needs T1 (M3) applied to D, which works verbatim because
 (i) D and pi_D are defined over F_p and pi_D^2 = p;
 (ii) (1+pi_D)/2 is in End(Jac D): the Weierstrass points 0, mu_N, infinity are F_p-rational since N | p-1
     (or: (1+pi_C)/2 preserves the connected kernel of pi_*);
 (iii) phi pi_D = pi_0 phi^{sigma^-1}, so phi pi_D phi^vee is in pi_0 O_B cap Z = pZ, giving v' in Z; (ii) gives u'=v' mod 2;
 (iv) D is superspecial (pi^* embeds H^0(D,Omega) Cartier-equivariantly into H^0(C_p,Omega), where Cartier = 0) and
     g_D >= 2, so Jac D = E_0^{g_D}, M_D^# = P^-1 M_D, and b'_D = Tr(lambda beta_D) is even unimodular of rank 4 g_D.
Also u'^2 > p v'^2 for phi' != 0 (deg(phi'(a+b pi_D)) > 0), so q' is totally positive and N(q') is an integer >= 1.
So the g_D=2 claim N(q') >= 1 is a theorem; state (i)-(iv) explicitly.  "unimodular 원소" in item (3) also relies on (iv).

OVERCLAIM (line 47, "랭크 8 고정격자가 모든 p 에서 정확히 경계인 이유"): N(q') >= 1 only proves "never violating" on
Fix(C_{(p-1)/4}).  "Exactly borderline" needs a vector with q' a totally positive unit; this is numerical (and follows
from the degree-2 bielliptic map of T1 (M5) only when its target curve is the E_0 in use, or under a class-number-one
type hypothesis for the rank-one module M_D over the maximal order).  Reword.

## (f) arithmetic: CORRECT
* threshold: ((p-1)/4)^2 / m^2 = (g_D/2)^2; g_D/2 is the Delta-rank of M_D, so it is genuinely the same line-stability
  condition (T1 M4) for the unimodular family (M_D, Q'_t).
* With eps=(a+b sqrt p)/2, a^2-p b^2=-4: norm 2l = Tr(lambda 2q') = a v'+b u', s' = Tr(eps q') = (a u'+p b v')/2,
  and p l^2 - s'^2 = u'^2 - p v'^2 = 4 N(q').  So N(q') = (p l^2 - s'^2)/4, identical to N(nu)=(p norm^2-4 s^2)/16
  with norm=2l; scaling norm=m norm', s=m s' is consistent.  p=37, m=3: (6,9) -> (2,3), 63/9=7, 81/9=9. OK.
* roots: s' odd, s'^2<p; condition s'^2 <= p-36.  p=37: |s'|=3 gives N=7 <9; q'=(19-3 sqrt37)/2, i.e. a degree-19 map
  Y^2=X^13-X -> E_0 with v'=-3 (degree 57, v=-9 on C_37; (57^2-37*81)/4=63).  The EXISTENCE of these vectors is the
  numerical result of sec:weilcomp transported through thm:MW (tag PL), not a consequence of the proposition.

## Side remark (line 53): correct, and can be strengthened to "superspecial".

## Strongest statement actually proved
For p=1 mod 4, m odd | p-1, E_0/F_p supersingular: pi: C_p -> D_{p,m} is the quotient by C_m (degree m, totally ramified
exactly over 0 and infinity, g_D=(p-1)/(2m) even); psi' -> psi' o pi_* is an O_F-linear bijection
Hom(Jac D,E_0) -> Hom(Jac C_p,E_0)^{C_m} = Norm_{C_m}(Hom(Jac C_p,E_0)) with beta_C = m beta_D, hence q_F = m q'_F and
Q_t = m Q'_t for all t; D is superspecial, q'_F takes totally positive values in O_F, (M_D, b'_D) is even unimodular of
rank 4 g_D; N(q_F) >= ((p-1)/4)^2 on the fixed lattice <=> N(q'_F) >= (g_D/2)^2, which holds automatically for g_D=2.
The link to W_p (and all statements about which values N(q') actually occur, e.g. 7 at p=37, equality for g_D=2) remains
conditional on thm:MW / numerical.
(No computer check was run in the time limit; all verifications above are by hand.)
