# R10 referee notes (adversarial check of builders A and B)

Files: ref.py (exact Gram checks, independent S/T action on A's saved complex basis, commutant/order, O.w determinants),
enum_ref.py (independent t-grid Fincke-Pohst using the R2 LLL/enumeration code), chk_F.py (primitivity / exact R_U of
torus-fixed sublattices, p=17 certificate), logs log_ref_p*.txt, log_enum_*.txt, ref_p*.json.

## Verified by the referee (exact integer arithmetic unless stated)
* A and B Gram pairs (B,E) for p=5,13,29,37,61 (A also 17): symmetric, B positive definite (all leading minors >0, Bareiss),
  det B = D^(p-1) with D = 5,117,153,1421,2997,13725; M=B^-1E integral, M^2=p, hence E B^-1 E = pB; (I+M)/2 integral (O_F acts).
  A and B have identical det for every common prime. Isometry A ~ B NOT tested.
* A's saved basis BasC: with my own S_a,T_a (kappa=-1: S^2=-1,(ST)^3=-1,T^p=1) the images of the basis are integral
  (defect <2e-11) and the resulting integer matrices are exact isometries of B and E and satisfy the relations exactly.
  Gram from the embedding equals the saved B,E to 4e-11. So A's lattice is an SL_2(p)-stable lattice of rank 2(p-1).
* Commutant: 8 integer matrices X_j commuting exactly with S,T, closed under multiplication, containing sqrt(p) and (1+sqrt p)/2;
  discriminant det(Tr_{F/Q} trd(x_i xbar_j)) = p^4 for p=5,13,17,29,37,61 => End_G(Lambda_0) is a MAXIMAL order.
* det_B(O w) = (N/16)^4 = ((p/4) m+ m-)^4 exactly, N = p b^2 - e^2, for the violators (p=29,37) and 3 random vectors per prime.
* p=37 violator (A and B files): b=1026, e=5994, N=3020976 < 3884112=(p-1)^2 D, Gamma/theta = 7/9; my own O.w has
  det_B = 1270893532472362625841 < r^8 D^4 = 3472858170158755180401, R_U >= (9/7)^4 = 2.7326.  p=37 is UNSTABLE (for Lambda_0).
* p=29: A's vector N=795760 (5/7), B's vector N=1023120 (45/49), threshold 1114064: both violators.
* p=17: A's violator_p17.json coordinate vector is STALE (gives b=294,e=306, not a violator in the saved basis);
  but basis vector 0 has b=36,e=0,N=22032<39168 (9/16) and det_B(F_H)=3595305184641 gives R_U>=9.98872. p=17 unstable: confirmed.
* Torus-fixed sublattices (A's F files): primitive, fixed by d(a) (which equals the word T^a S T^(1/a) S T^a S up to sign),
  exact R_U = 1 for p=13 (C3), 29 (C7), 37 (C9, C3), 61 (C15, C5, C3).
* Independent enumeration (double precision FP, grid 0.1 over one period, bound 2 sqrt(theta D) cosh(0.05)(1+1e-9)):
  p=5 full: 60 pairs all border; p=13 full rank 24: 19656 pairs all border, 0 violators; p=37 C9 line: 12 border;
  p=37 C3 (rank 24): 72 distinct violating +-pairs (36 with e=+5994, 36 with e=-5994) + 12 border, min Gamma/theta = 7/9;
  p=61 C15 line, C5-fixed (24), C3-fixed (40; 2.9e8 nodes): 12 border pairs, 0 violators.

## Corrections / caveats
* A's "143 distinct violating +-pairs" in the p=37 C3 lattice is not reproduced: 72 (agrees with B).
* A's violator_p17.json is inconsistent with the saved basis (see above).
* "order maximal" and "split-torus line has R_U=1" are the same determinant identity (det_B(Lambda^H) = p^4 m^16 = r^8 D^4).
* p=61 is NOT decided; only three torus-fixed sublattices are clean. Full rank 120 unsearched in any meaningful sense.
* All statements concern the single G-lattice Lambda_0 = Z[G]v0 (any invariant form on it, incl. the even unimodular one).
