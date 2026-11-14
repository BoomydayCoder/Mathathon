/- Independent statement/axiom check for pinned Romik formalization.
   Run with the pinned checkout's `lake env lean` after a successful build. -/
import Zeta5

example : _root_.Irrational
    (_root_.tsum (fun v : _root_.Nat =>
      (1 : _root_.Real) / ((v : _root_.Real) + 1) ^ 5)) :=
  Zeta5.zeta5_irrational

#print Zeta5.zeta5
#print Zeta5.Axioms.chebyshev_theta_asymptotic
#print axioms Zeta5.zeta5_irrational
