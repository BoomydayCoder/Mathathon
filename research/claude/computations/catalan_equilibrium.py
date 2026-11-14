#!/usr/bin/env python3
"""
catalan_equilibrium.py -- the equilibrium problem governing log Delta_K(X0)/K^2 for kernels ~ e^{-gamma y}
with nodes of density 1 (generalises equilibrium_s3.py, which is gamma = 2 pi).

  minimise over probability measures rho on [0,inf):  E(rho) = int V drho - iint log|x^2-y^2| drho drho,
  V(x) = gamma x + int_0^1 log(x^2+u^2) du.          Put A = 1 + gamma/pi.
Ansatz (same derivation as in equilibrium_s3.py: V'(x) = gamma + 2 arctan(1/x) = (A pi) - 2 arctan x - ... ):
  support [0,b], b = 2A/(A^2-1), c = sqrt(1+b^2) = (A^2+1)/(A^2-1);
  symmetrised density (mass 2) sigma(x) = (1/pi)[A artanh(s/b) - artanh(s/c)], s = sqrt(b^2-x^2).
Checked here NUMERICALLY (30 digits): mass = 2; Euler-Lagrange equality V - 2U = ell on the support;
strict inequality V - 2U > ell off the support (grid); E* = (1/4) int V dsigma + ell/2.
Results:  gamma = 2 pi (A=3): E* = (9/2) log 2 = 3.1192 (zeta kernels);
          gamma =   pi (A=2): E* = 2 log 3    = 2.1972 (sinh / cosh kernels: Catalan, log 2, pi);
          other A: values printed (no closed form claimed; the guess (A^2/2) log((A+1)/(A-1)) is FALSE for
          A = 3/2, 5/2, 4 although it holds at A = 2, 3).
STATUS: CHECKED NUMERICALLY ONLY (no symbolic proof written for A != 3).
"""
import mpmath as mp

mp.mp.dps = 30


def solve(A):
    A = mp.mpf(A)
    b = 2 * A / (A * A - 1)
    c = mp.sqrt(1 + b * b)
    gam = (A - 1) * mp.pi

    def sig(x):
        x = abs(x)
        if x >= b:
            return mp.mpf(0)
        s = mp.sqrt(b * b - x * x)
        return (A * mp.log((b + s) / x) - mp.atanh(s / c)) / mp.pi

    def V(x):
        x = abs(x)
        return gam * x + mp.log(1 + x * x) - 2 + 2 * x * mp.atan(1 / x)

    def U(x):
        pts = sorted(set([-b, -b / 2, 0, b / 2, b] + ([x] if -b < x < b else [])))
        return mp.quad(lambda y: mp.log(abs(x - y)) * sig(y), pts)

    mass = 2 * mp.quad(sig, [0, b / 2, b])
    IV = 2 * mp.quad(lambda x: V(x) * sig(x), [0, b / 2, b])
    ell = V(b) - 2 * U(b)
    on = max(abs(V(x) - 2 * U(x) - ell) for x in (b / 7, b / 3, b / 2, mp.mpf('0.9') * b))
    off = min(V(x) - 2 * U(x) - ell for x in (b * mp.mpf('1.01'), b * mp.mpf('1.1'), b * 1.5, b * 3, b * 10, b * 100))
    dens_pos = all(sig(b * mp.mpf(k) / 200) > 0 for k in range(1, 200))
    return dict(A=A, b=b, mass=mass, ell=ell, on=on, off=off, E=IV / 4 + ell / 2, dens_pos=dens_pos)


if __name__ == '__main__':
    ok = True
    for A in (mp.mpf(3) / 2, 2, mp.mpf(5) / 2, 3, 4):
        r = solve(A)
        closed = {2: 2 * mp.log(3), 3: mp.mpf(9) / 2 * mp.log(2)}.get(int(r['A']) if r['A'] == int(r['A']) else -1)
        good = abs(r['mass'] - 2) < 1e-20 and r['on'] < 1e-18 and r['off'] > 0 and r['dens_pos'] and (closed is None or abs(r['E'] - closed) < 1e-18)
        ok &= bool(good)
        print(f"A={mp.nstr(r['A'], 4)} gamma/pi={mp.nstr(r['A'] - 1, 4)} b={mp.nstr(r['b'], 8)} mass-2={mp.nstr(r['mass'] - 2, 3)} "
              f"EL-equality defect={mp.nstr(r['on'], 3)} min(V-2U-ell) off support={mp.nstr(r['off'], 5)} density>0:{r['dens_pos']}")
        print(f"     E* = {mp.nstr(r['E'], 25)}   closed form: {mp.nstr(closed, 25) if closed is not None else 'none claimed'}   {'PASS' if good else 'FAIL'}")
    print("2 log 3 =", mp.nstr(2 * mp.log(3), 25), "  (9/2) log 2 =", mp.nstr(mp.mpf(9) / 2 * mp.log(2), 25))
    print("ALL PASS" if ok else "FAILURES")
