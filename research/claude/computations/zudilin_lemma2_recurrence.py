#!/usr/bin/env python3
"""Zudilin math/0202159: Lemma 2 (creative-telescoping identity (6)), Lemma 3 (recurrence),
discovery of the recurrence by undetermined coefficients, and the footnote (telescoping of A2k).
All checks exact (sympy rationals / fractions)."""
from fractions import Fraction as Fr
from math import comb, lcm
import sympy as sp

results = []


def check(name, ok, extra=""):
    results.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {extra}" if extra else ""))


n, t = sp.symbols('n t')


def s_poly(nn, tt):
    return 4 * (2 * nn + 1) * (-2 * tt ** 2 + tt + (2 * nn + 1) ** 2)


def R(m, tt):
    if m < 0:
        raise ValueError
    num = sp.prod([tt - j for j in range(1, m + 1)]) if m else sp.Integer(1)
    den = sp.prod([tt + j for j in range(0, m + 1)])
    return (num / den) ** 2


def H(k, j):
    return sum((Fr(1, l ** j) for l in range(1, k + 1)), Fr(0))


def A2(m, k):
    return comb(m + k, m) ** 2 * comb(m, k) ** 2 if 0 <= k <= m else 0


def A1(m, k):  # closed form verified in zudilin_lemma1.py
    return 2 * A2(m, k) * (2 * H(k, 1) - H(m + k, 1) - H(m - k, 1))


def uv(m):
    u = 2 * sum(A2(m, k) for k in range(m + 1))
    v = 2 * sum((A2(m, k) * H(k, 3) for k in range(m + 1)), Fr(0)) + sum((A1(m, k) * H(k, 2) for k in range(m + 1)), Fr(0))
    return Fr(u), v


print("=" * 72)
print("(a) Reduced identity of Lemma 2, symbolic in n and t")
print("=" * 72)
lhs = (n + 1) ** 3 * ((t - n - 1) / (t + n + 1)) ** 2 - (2 * n + 1) * (17 * n ** 2 + 17 * n + 5) + n ** 3 * ((t + n) / (t - n)) ** 2
rhs = s_poly(n, t + 1) * (t ** 2 / ((t - n) * (t + n + 1))) ** 2 - s_poly(n, t)
check("(n+1)^3((t-n-1)/(t+n+1))^2 - (2n+1)(17n^2+17n+5) + n^3((t+n)/(t-n))^2 == s_n(t+1)(t^2/((t-n)(t+n+1)))^2 - s_n(t)",
      sp.cancel(sp.together(lhs - rhs)) == 0)
# sanity: the ratios used in the reduction
for m in range(1, 6):
    ok = sp.cancel(R(m + 1, t) / R(m, t) - ((t - m - 1) / (t + m + 1)) ** 2) == 0 \
        and sp.cancel(R(m - 1, t) / R(m, t) - ((t + m) / (t - m)) ** 2) == 0 \
        and sp.cancel(R(m, t + 1) / R(m, t) - (t ** 2 / ((t - m) * (t + m + 1))) ** 2) == 0
    check(f"  ratio formulas R_(n+1)/R_n, R_(n-1)/R_n, R_n(t+1)/R_n(t) correct for n={m}", ok)
# full identity (6) directly for n = 1..6
for m in range(1, 7):
    S = lambda tt: s_poly(m, tt) * R(m, tt)
    d = (m + 1) ** 3 * R(m + 1, t) - (2 * m + 1) * (17 * m ** 2 + 17 * m + 5) * R(m, t) + m ** 3 * R(m - 1, t) - (S(t + 1) - S(t))
    check(f"  identity (6) holds exactly as rational functions, n={m}", sp.cancel(sp.together(d)) == 0)

print()
print("=" * 72)
print("(b) Double zeros at t=1 and S_n'(1) = 0, n = 1..6")
print("=" * 72)
for m in range(1, 7):
    Rm = R(m, t)
    Sm = s_poly(m, t) * Rm
    r0, r1 = sp.simplify(Rm.subs(t, 1)), sp.simplify(sp.diff(Rm, t).subs(t, 1))
    s0, s1 = sp.simplify(Sm.subs(t, 1)), sp.simplify(sp.diff(Sm, t).subs(t, 1))
    check(f"n={m}: R_n(1)=R_n'(1)=0, S_n(1)=S_n'(1)=0", r0 == 0 and r1 == 0 and s0 == 0 and s1 == 0)
    # order of zero of R at t=1 is exactly 2
    num = sp.numer(sp.together(Rm))
    mult = sp.roots(sp.Poly(num, t)).get(1, 0)
    check(f"n={m}: multiplicity of zero of R_n at t=1 is exactly 2", mult == 2, f"(mult={mult})")
# also n=0: S_0'(1) != 0 explains why the recurrence starts at n=1
S0p = sp.diff(s_poly(0, t) * R(0, t), t).subs(t, 1)
print(f"  info: S_0'(1) = {S0p} (nonzero; Lemma 2/3 need n>=1)")

print()
print("=" * 72)
print("(c) Apery recurrence L[X]_n = (n+1)^3 X_(n+1) - (2n+1)(17n^2+17n+5) X_n + n^3 X_(n-1) = 0")
print("=" * 72)
U, V = {}, {}
for m in range(0, 12):
    U[m], V[m] = uv(m)
for m in range(1, 11):
    Lu = (m + 1) ** 3 * U[m + 1] - (2 * m + 1) * (17 * m ** 2 + 17 * m + 5) * U[m] + m ** 3 * U[m - 1]
    Lv = (m + 1) ** 3 * V[m + 1] - (2 * m + 1) * (17 * m ** 2 + 17 * m + 5) * V[m] + m ** 3 * V[m - 1]
    check(f"n={m}: L[u]_n = 0 and L[v]_n = 0 (exact)", Lu == 0 and Lv == 0)
# the paper prints (7) with "n^3 u_n" in the last term; show that this literal form fails
Lbad = 2 * U[2] - 3 * 39 * U[1] + 1 * U[1]
print(f"  info: paper's eq. (7) as printed reads '... + n^3 u_n'; with u_n (not u_(n-1)) at n=1: {Lbad} != 0 -> typo, must be n^3 u_(n-1)")
check("literal (7) with n^3 u_n fails (confirms typo)", Lbad != 0)

print()
print("=" * 72)
print("(d) Discovery by undetermined coefficients, fixed numeric n = 1..8")
print("    a R_(n+1)(t) + b R_n(t) + c R_(n-1)(t) = s(t+1)R_n(t+1) - s(t)R_n(t),  s polynomial in t")
print("=" * 72)


def solve_space(m, sdeg, use_c=True, use_a=True):
    a, b, c = sp.symbols('a b c')
    ss = sp.symbols(f's0:{sdeg + 1}')
    s = sum(ss[i] * t ** i for i in range(sdeg + 1))
    unknowns = ([a] if use_a else []) + [b] + ([c] if use_c else []) + list(ss)
    terms = [(a if use_a else 0, R(m + 1, t)), (b, R(m, t)), (c if use_c else 0, R(m - 1, t)),
             (s.subs(t, t + 1), -R(m, t + 1)), (s, R(m, t))]
    dens = [sp.denom(sp.together(r)) for _, r in terms]
    L = sp.lcm(dens) if len(dens) > 1 else dens[0]
    L = sp.lcm_list(dens)
    expr = sum(coef * sp.cancel(r * L) for coef, r in terms)
    P = sp.Poly(sp.expand(expr), t)
    eqs = P.coeffs()
    M = sp.Matrix([[sp.diff(e, x) for x in unknowns] for e in eqs])
    ns = M.nullspace()
    return unknowns, ns


dims = {}
for m in range(1, 9):
    unknowns, ns = solve_space(m, 2)
    dims[(m, 'deg2')] = len(ns)
    ok = len(ns) == 1
    msg = f"dim = {len(ns)}"
    if ok:
        vec = ns[0]
        vec = vec / vec[0] * (m + 1) ** 3
        sol = dict(zip([str(x) for x in unknowns], vec))
        b_ok = sol['b'] == -(2 * m + 1) * (17 * m ** 2 + 17 * m + 5)
        c_ok = sol['c'] == m ** 3
        s_target = sp.Poly(s_poly(m, t), t).all_coeffs()[::-1]
        s_ok = [sol['s0'], sol['s1'], sol['s2']] == list(s_target)
        ok = b_ok and c_ok and s_ok
        msg += f"; normalized a=(n+1)^3={vec[0]}: b={sol['b']}, c={sol['c']}, s=({sol['s0']}) + ({sol['s1']})t + ({sol['s2']})t^2"
    check(f"n={m}: deg s <= 2 solution space is 1-dim and reproduces b, c, s_n", ok, msg)

print()
print("Solution-space dimensions for variants (n = 1..8):")
variants = [("deg s<=0", 0, True), ("deg s<=1", 1, True), ("deg s<=2, c=0 (first order)", 2, False),
            ("deg s<=2", 2, True), ("deg s<=3", 3, True), ("deg s<=4", 4, True), ("deg s<=6", 6, True)]
table = {}
for name, d, use_c in variants:
    row = []
    for m in range(1, 9):
        _, ns = solve_space(m, d, use_c=use_c)
        row.append(len(ns))
    table[name] = row
    print(f"  {name:32s}: {row}")
check("deg s<=1 admits only the zero solution for n=1..8", all(x == 0 for x in table["deg s<=1"]))
check("first-order (c=0, deg s<=2) admits only zero solution for n=1..8", all(x == 0 for x in table["deg s<=2, c=0 (first order)"]))
check("deg s<=3 and deg s<=4 still give 1-dim space (no new solutions)",
      all(x == 1 for x in table["deg s<=3"]) and all(x == 1 for x in table["deg s<=4"]))
# also c=0 with higher degree s
row = [len(solve_space(m, 6, use_c=False)[1]) for m in range(1, 9)]
print(f"  {'deg s<=6, c=0':32s}: {row}")

print()
print("=" * 72)
print("(e) Footnote: multiply (6) by (t+k)^2 and set t = -k")
print("=" * 72)
print("  (t+k)^2 R_m(t)|_{t=-k} = A2k^(m);  (t+k)^2 S_n(t+1)|_{t=-k} = s_n(1-k) A2,(k-1)^(n)  [R_n(t+1) has its double pole at t+1 = -(k-1), i.e. t = -k]")
print("  so the k-th coefficient identity is")
print("  (n+1)^3 A2k^(n+1) - (2n+1)(17n^2+17n+5) A2k^(n) + n^3 A2k^(n-1) = s_n(1-k) A2,(k-1)^(n) - s_n(-k) A2k^(n) = B_(k-1) - B_k,")
print("  with B_k := s_n(-k) A2k^(n)  (B_(-1) = B_(n+1) = 0).")
allok = True
for m in range(1, 7):
    for k in range(0, m + 2):
        L = (m + 1) ** 3 * A2(m + 1, k) - (2 * m + 1) * (17 * m ** 2 + 17 * m + 5) * A2(m, k) + m ** 3 * A2(m - 1, k)
        B = lambda kk: s_poly(m, -kk) * A2(m, kk)
        if L != B(k - 1) - B(k):
            allok = False
            print(f"   mismatch n={m}, k={k}: L={L}, B_(k-1)-B_k={B(k-1)-B(k)}")
    # also directly from rational function: coefficient extraction
    S = s_poly(m, t) * R(m, t)
    for k in range(0, m + 2):
        expr = ((m + 1) ** 3 * R(m + 1, t) - (2 * m + 1) * (17 * m ** 2 + 17 * m + 5) * R(m, t) + m ** 3 * R(m - 1, t)) * (t + k) ** 2
        val = sp.cancel(expr).subs(t, -k)
        rhsv = sp.cancel((S.subs(t, t + 1) - S) * (t + k) ** 2).subs(t, -k)
        if val != rhsv or val != (s_poly(m, 1 - k) * A2(m, k - 1) - s_poly(m, -k) * A2(m, k)):
            allok = False
            print(f"   rational-function mismatch n={m}, k={k}")
check("footnote coefficient identity, n=1..6, k=0..n+1 (both from formula and from (6) directly)", allok)
sumok = all(sum((m + 1) ** 3 * A2(m + 1, k) - (2 * m + 1) * (17 * m ** 2 + 17 * m + 5) * A2(m, k) + m ** 3 * A2(m - 1, k)
                for k in range(0, m + 2)) == 0 for m in range(1, 7))
check("summing over k telescopes to 0 => L[u]_n = 0, n=1..6", sumok)

print()
print(f"SUMMARY: {sum(results)}/{len(results)} checks PASS" + ("" if all(results) else "  -- SOME FAILED"))
