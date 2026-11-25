# R11 / T1 — literature + derivation check of the "Mordell-Weil reformulation"
(start 08:26 UTC, limit 26 min; files: heur.py, heur_out.txt, downloads/ (untrusted, only text-extracted))

## 0. Access caveat
Blocked hosts in this session (DNS refused): www.ams.org, resolve.cambridge.org, www.numdam.org,
www.quantamagazine.org, www.math.rwth-aachen.de.  So Gross (JAMS 1990), Dummigan (JLMS 1997), Nebe (J. Alg. 1998),
Watkins (JTNB 2012) could NOT be read; only search snippets / arXiv.  Everything about their content is marked
[memory] or [snippet].  Reached: arxiv.org (2606.05105 html; 2109.09993 pdf).

## 1. Literature (what is actually established here)
* Venkatesh conjecture — SEEN (arXiv:2606.05105, Abuya-Gargava-Zhao, "Stochastically evolving ellipsoids with
  symmetries", Introduction, via fetch): "On various occasions, Venkatesh has conjectured the upper bound
  Delta_N^lat <= N^2 (log N)^{O(1)} 2^{-N}. This is motivated by the fact that lattices in R^N have on the order of
  N^2 degrees of freedom."  No written source by Venkatesh is cited for it (oral conjecture); cite it as
  "Venkatesh, as reported in [AGZ, arXiv:2606.05105, Sect. 1]".  Their Thm 1.1: Delta_N^lat >= c N^2 loglog N 2^{-N}
  along an infinite sequence of N.  (Klartag arXiv:2504.05042 abstract: ellipsoid of volume c n^2; no Venkatesh mention there.)
* Gross, "Group representations and lattices", J. AMS 3 (1990) 929-960 — existence/bibliography seen; content NOT seen.
  [memory] it treats Mordell-Weil lattices Hom(Jac C, E) for curves with large automorphism group (SL_2(q), SU_3(q)).
  I could not verify that "y^2=x^p-x" and "p=13 gives Leech" are stated there.
* Tiep, "Globally irreducible representations of finite groups and integral lattices", Geom. Dedicata 64 (1997) —
  [snippet of abstract] a globally irreducible rep related to the Weil rep of degree (p^n-1)/2 of Sp_2n(p), p = 1 mod 4,
  "new series of even unimodular lattices of rank 2(p^n-1)" (n=1 is the SL_2(p) family).  Consistent with Sect. 2 below.
* Dummigan, "Algebraic cycles and even unimodular lattices", J. London Math. Soc. (1997) — landing page only [snippet];
  content NOT seen.  [memory, low confidence] realizes these lattices via cycles/Mordell-Weil and gives lower bounds
  for the minimum of the type coming from intersection theory (i.e. of Castelnuovo-Severi strength, not sqrt(p)).
* Nebe, "Some cyclo-quaternionic lattices", J. Algebra 199 (1998) 472-498 — bibliographic record seen; abstract [memory]:
  "The faithful lattices of rank 2(p-1) of the groups SL_2(p) are described. For small primes p these and related
  lattices are investigated by computer. In particular a new extremal even unimodular lattice of rank 48 is constructed."
* Watkins, "Another 80-dimensional extremal lattice", JTNB 24 (2012) — [snippet] "the unimodular lattice associated to
  the rank 20 quaternionic matrix group SL_2(F_41) (x) S~_3 in GL_80(Z) is a fourth example of an 80-dimensional
  extremal lattice" (min 8).  Dimension 80 = 2(41-1): an SL_2(41)-invariant even unimodular lattice in W_41 has min 8.
  (Needed at t=0 for p=41: min >= 40/sqrt(41)=6.25, OK; but p=41 is unstable anyway, see 3(M5).)  Not verified that it
  is the maximal-order lattice.
* Minima for p = 17, 29, 37, 61: NOT found in any accessible source.  [memory, unverified] dim 32 (p=17) min 4,
  dim 56 (p=29) min 6 (extremal), dim 72 (p=37) min 6.  Treat as unknown until Nebe's paper is read.
* arXiv:2109.09993 (seen, text-extracted): Thm 5.3 quoting [BFM94]: even unimodular lattices of dim < 48 from the
  ideal-lattice construction on quaternion algebras have min 2 — a different (rank-1 over the algebra) construction; irrelevant.

## 2. Derivation check (all PROVED below unless marked)
Notation: C=C_p, J=Jac(C), g=(p-1)/2, r=(p-1)/4, p = 1 mod 4, E_0/F_p supersingular with Frobenius pi_0 (pi_0^2=-p),
O_B=End(E_0) maximal order in B_{p,infty}, P=pi_0 O_B, M=M_p=Hom_{Fbar_p}(J,E_0), b(phi,psi)=deg(phi+psi)-deg phi-deg psi.

(M1) CORRECT.  (ax^p+b)(cx+d)-(ax+b)(cx^p+d)=(ad-bc)(x^p-x), so y -> y/(cx+d)^{(p+1)/2} works; -1 acts as iota since
 (p+1)/2 is odd.  Order-p element: uniformizer t=x^{(p-1)/2}/y at infinity, v(tau t - t)=3, Lefschetz 3 = 2 - Tr => Tr=-1. OK.
 Extra check that H^1 = psi+psi' and not a discrete-series rep of degree p-1: a generator s of the non-split torus has no
 fixed point on C (it swaps the two points over each of its fixed x's because s^{(p+1)/2}=iota and (p+1)/2 is odd), so
 Tr(s)=2, impossible for -(theta+theta^{-1})(s) with theta^2 != 1.
(M2) CORRECT with two precisions.
 (i) Point counts: x in F_{p^2}\F_p gives z=x^p-x with z^{p-1}=-1, z^{(p^2-1)/2}=(-1)^{(p+1)/2}=-1: non-square.  So
   #C(F_p)=#C(F_{p^2})=p+1, Tr(Frob^2)=p(p-1)=p*2g, hence pi^2=p on J.  C is MINIMAL over F_{p^2} (Frob_{p^2}=+p); it is
   the quadratic twist over F_{p^2} of the quotient w^2=v^p+v of the Hermitian curve (maximal), so "quotient of the
   Hermitian curve" holds only after twisting / over F_{p^4}.  Superspecial: Cartier operator kills x^i dx/y, 0<=i<=g-1
   (needs j = i+(p+1)/2 mod p with 0<=j<=(p-1)/2: impossible).  So J ~ E^g over Fbar_p.
 (ii) E must have Frob_{p^2}=+p (not the base change of E_0/F_p, which has -p); work over Fbar_p to avoid this.
   O_F=Z[(1+pi)/2] acts: pi=1 on J[2] since all Weierstrass points are F_p-rational.  pi is Rosati-symmetric (pi^2=p).
   O_F (x) O_B has reduced discriminant p O_F and is NOT maximal; but phi -> pi_0^{-1} phi pi_J (= Galois conjugation
   sigma^{-1} on Hom) is an automorphism xi of M with xi^2=-1, xi = pi_0^{-1} (x) sqrt p in Delta, and
   O_F (x) O_B + xi(...) is maximal at sqrt p.  Hence End_G(M) is a MAXIMAL order and contains a unit xi of order 4
   not in F (so [O^x:O_F^x] is even).
(M3) CORRECT, and now proved: write b = Tr_{F/Q}(beta), beta the O_F-bilinear form.  Facts: (a) M^# = P^{-1}M (J=E^g,
 principal polarization, d(O_B)=p^2), det b = p^{2g}=p^{p-1};  (b) M pi_J = pi_0 M = P M (pi_0 phi = phi^sigma pi_J);
 (c) b(phi alpha, phi) is even for every Rosati-symmetric alpha (it is the reduced trace of the symmetric element
 phi alpha phi^vee of End(E_0), i.e. twice an integer).  From (a),(b): beta(M,M) in O_F and M is beta-unimodular over O_F;
 from (c) with alpha=1,(1+pi)/2: beta(phi,phi)=2 q_F(phi), q_F(phi)=(u+v sqrt p)/2 in O_F, u=deg phi, v=b(phi pi,phi)/(2p),
 u = v mod 2.  With lambda=eps/sqrt p (totally positive, generates the inverse different, N(lambda)=1/p since N(eps)=-1):
   b' := Tr_{F/Q}(lambda beta) is EVEN UNIMODULAR on M, G-invariant, and 2deg = Tr(beta) = sqrt p * Q_{t=-log eps}
 in the circle Q_t = e^t b'_+ + e^{-t} b'_-.  So (M_p,2deg) is the antipode (t = log eps mod 2 log eps) of the
 unimodular lattice, up to the factor sqrt p; it is p-modular (via pi_0).  p=5: (M_5,2deg) 5-modular, min 4, deg 2;
 p=13: min deg 6 (given round-5 line stability + existence, see M5).  [The identification "even unimodular rank 24
 with irreducible SL_2(13) => Leech" is round 5's computation.]
(M4) CORRECT up to normalisation.  For a principal line U=Delta w:  det(W_Z cap U)=N_{F/Q}(q_F(w))^4, so
   R_U <= 1  <=>  N_{F/Q}(q_F(w)) >= r^2 = ((p-1)/4)^2      [m_+ m_- = 4 N(q_F)/p].
 In degrees: deg(phi)^2 - <phi,phi Frob>^2/p >= ((p-1)/2)^2 with <phi,phi>=deg phi; the orchestrator's "/(4p)" is right
 only if <.,.> denotes b (b(phi,phi)=2deg).  AM-GM gives deg >= (p-1)/2 = g.  Castelnuovo-Severi gives only
 deg >= (g+1)/2, i.e. deg >= (p+3)/4; applied to all eps^k phi it yields the unconditional bound
 N(q_F) >= ((p+3)/(4 b sqrt p))^2 (eps=(a+b sqrt p)/2) ~ p/(16 b^2), versus the needed p^2/16.
 General rank j: det(W_Z cap U) = N_{F/Q}(Mdet H_U)^4 (Moore determinant of the hermitian Gram matrix over the maximal
 order), so stability <=> N_{F/Q}(Mdet H_U) >= (r/j)^{2j} for all free primitive rank-j sublattices.
 p=13 check: Leech norm-4 vectors must have (deg,v) in {(7,-1),(10,-2),(13,-3),(16,-4),(19,-5)} (N(q_F)=9,12,13,12,9);
 the types (4,0) and (22,-6) (N=4) are exactly the excluded violators; degree-6 maps (q_F=3) are norm-6 Leech vectors.
(M5) First part CORRECT: p = 5 mod 8, m=(p-1)/4 odd, (X,Y)=(x^m, y x^{(m-1)/2}) maps C_p -> D: Y^2=X^5-X, degree m, Galois
 with group = fourth powers in the split torus; D is bielliptic (X -> (1-X)/(1+X) is an edge-axis involution of the
 octahedron), E has CM by Z[sqrt -2], supersingular since (-2/p)=-1.  Total degree (p-1)/2.
 Last sentence needs CORRECTION.  Classification (H < Aut(C_p)=2.PGL_2(p), iota not in H; then H cap SL_2(p) has odd
 order; branch points of C_p/H -> P^1/Hbar = Hbar-orbits on P^1(F_p) plus orbits where the stabiliser swaps the fibre):
   - cyclic split C_k: genus >= 2 (=2 iff k=(p-1)/4, p=5 mod 8);   - cyclic non-split odd C_k: genus (p+1)/(2k)-1, never 1;
   - Borel p:k: genus 0;
   - dihedral D_k, k odd | p-1 (involutions x -> c/x, c non-square): (p-1)/(2k)+2 branch points; genus 1 iff k=(p-1)/4, p=5 mod 8;
   - C_k x <j> and D_k with k odd | (p+1)/2: (p+1)/(2k)+1 branch points; genus 1 iff k=(p+1)/6.
 So: if 3 | p+1 there ARE Galois elliptic quotients of degree (p+1)/3 < (p-1)/2 (p=17: degree 6 < 8; p=29: 10 < 14;
 p=41: 14 < 20) — these p have (p+1)/2 composite, consistent with round 5.  If (p+1)/2 is prime and p>5 the claim holds:
 only dihedral groups of order (p-1)/2, and only for p = 5 mod 8.  For p = 1 mod 8 with (p+1)/2 prime (p=73, 193, ...)
 there is NO Galois elliptic quotient at all.  Also j: x -> c/x gives C_p/j of genus r=(p-1)/4.
(M6) CORRECT (phi iota = -phi after translation; descend to R: P^1 -> E/+-1; converse by normalised fibre product).
 Count: 4d >= p+1 recovers Castelnuovo-Severi; parameter count 2d+2 versus 2d+(p+1)/2 conditions gives no information.

## 3. (M7) re-derived, and planes  (heur.py, heur_out.txt) — HEURISTIC
Violating vectors mod eps^Z in a covolume-1 lattice of R^h x R^h, h=p-1: volume of {m_+m_- < c, 0 <= log(rho_+/rho_-) < 2log eps}
 = h omega_h^2 c^{h/2} log eps, c=(p-1)^2/(4p).  Stirling: (log eps/pi)(pi e/sqrt p)^{p-1}.  CONFIRMED (exact / Stirling):
   p=13: 1.153e4 / 1.185e4;  p=37: 1.584e5 / 1.60e5;  p=61: 246.1 / 247.5;  p=73: 2.343 / 2.354;  p=157: 1.80e-26.
 Per |G|: 5.28 (13), 15.3 (17), 8.6 (29), 3.13 (37), 1.9 (41), 1.08e-3 (61), 6.0e-6 (73), 4.7e-33 (157).
 Correction: vectors of one line come in packets of size 2[O^x:O_F^x] >= 4 (xi above), so the expected number of
 G-orbits of violating LINES is N_1/(|G| [O^x:O_F^x]) <= half the orchestrator's figures.
 The model already fails conservatively at p=13 (predicts ~0.4-5 orbits; truth: none, only borderline).
Planes (rank 2), mean-value heuristic by double counting flags (line in plane) with the rank-1 formula applied to the
 quotient lattice: #planes with N(Mdet H) < Y  ~ A_2 Y^{2r},  A_2 = (r-1) A_1(r) A_1(r-1) / (2u A_1(2)),
 A_1(s) = 4s omega_{4s}^2 (4/p)^{2s} log eps, u = 2[O^x:O_F^x] (u=2 used; zeta-factors ~1 dropped).  Threshold Y=(r/2)^4:
   expected violating planes: p=13: 183 (0.17 orbits); p=29: 0.18; p=37: 1.7e-5 (7e-10 orbits); p=61: 3e-23; p=73: 6e-34; p=157: 7e-133.
 Planes (and higher rank, by the factor j^{-4jr}) are heuristically irrelevant for p >= 37; lines are the only issue.

## 4. Cleanest open problem
For p = 1 mod 4 with (p+1)/2 prime (first open case p=37, n=72; then 61, 73, 157, ...):
 (L)  N_{F/Q}(q_F(phi)) >= ((p-1)/4)^2 for all non-constant phi: C_p -> E, i.e.
      deg(phi)^2 - <phi, phi Frob>^2/p >= ((p-1)/2)^2;       necessary: no map C_p -> E of degree < (p-1)/2,
      i.e. (M6) has no solution with d < (p-1)/2;  necessary at t=0: min(W_p) >= 6 (p=37), 8 (61), 10 (73), 14 (157).
 (S)  N_{F/Q}(Mdet H_U) >= (r/j)^{2j} for all rank-j Delta-sublattices, 2 <= j <= r-1.
Known: p=5 (trivial, r=1), p=13 (round 5 computer proof).  Unconditional: only deg >= (p+3)/4.  Expected: (L) true for
p >= 61 unless a structured low-degree map exists; p=37 genuinely undecided (heuristic ~1 orbit).  Equality in (L) is
attained for p = 5 mod 8 (the degree-(p-1)/2 maps through D), so (L) is sharp there; for p = 1 mod 8 no borderline
vector is known.

## 5. Cross-checks done
* NUMERICAL: #C_p(F_{p^2}) = p+1 for p = 5, 13, 17, 29 by brute force (scratch script) — confirms pi^2 = p.
* Consistency with round 5 (R5_leech13/notes.md): their nu(w) has N(nu) = 13 m_+ m_-/4 = my N(q_F), threshold 9 = r^2. OK.
  Their norm-6 Leech vectors with N=9, t = w.eta w = +-18 (26208 of them) are exactly w_0 and eps^{-1} w_0 with
  q_F(w_0)=3, i.e. (deg, v) = (6,0) [t = a*u + 13*b*v = 18 for eps=(3+sqrt13)/2] and (33,-9).  Hence [derived, not
  independently recomputed]: C_13 has exactly 13104 = 6|G| maps to E of degree 6 with <phi,phi Frob>=0 (mod translation),
  lying in 546 Delta-lines (24 = |O^1| per line).  Only 91 of these lines (one per split torus) are pointwise fixed by
  an element of order 3, i.e. come from C_13 -> D -> E of (M5); the other 455 lines carry degree-6 maps NOT explained
  by Galois quotients.  So the borderline locus is larger than the "structured" part even at p=13.
* Mass check: zeta_F(-1) for p=37 is 5/6 (Siegel formula), mass 5/12, so h(O) > 1 possible for p=37: the "ideal version"
  of the line criterion (non-principal lines) is genuinely needed from p=37 on.
