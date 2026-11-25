# Referee re-implementation (independent of glue_search.py).
# Uses eq. (6) of arXiv:2610.12619 DIRECTLY in the y-variables:
#   A(S) = (c_S/2) log(C/c_S) + (c_S/C) I - I_inc(S) - sum_{n in S} c_n y_n,   sum_n c_n y_n = 0,
#   I_inc(S) = sum over pairs {m,n} in N with {m,n} cap S nonempty of log|Res(Phi_m,Phi_n)|.
# margin(N) = max_y min_{S proper nonempty} ( -A(S) );  support-admissible iff margin > 0.
import math, sys, itertools
import numpy as np
from scipy.optimize import linprog
from sympy import totient, factorint, cyclotomic_poly, resultant, Symbol, Integer

X = Symbol('x')

def w_apostol(m, n):
    if m == n: return 0.0
    a, b = min(m, n), max(m, n)
    if b % a: return 0.0
    f = factorint(b // a)
    if len(f) != 1: return 0.0
    return float(totient(a)) * math.log(next(iter(f)))

def w_exact(m, n):
    """log|Res(Phi_m,Phi_n)| from the exact integer resultant (sympy)."""
    R = resultant(cyclotomic_poly(m, X), cyclotomic_poly(n, X), X)
    return math.log(abs(int(R)))

def margin_y(N, exact=False, return_y=False):
    """LP in the paper's own y-variables, brute-force over all supports (no vectorisation tricks)."""
    k = len(N)
    c = [int(totient(n)) for n in N]
    C = sum(c)
    wf = w_exact if exact else w_apostol
    W = [[0.0] * k for _ in range(k)]
    for i in range(k):
        for j in range(i + 1, k):
            W[i][j] = W[j][i] = wf(N[i], N[j])
    I = sum(W[i][j] for i in range(k) for j in range(i + 1, k))
    rows, rhs = [], []
    for mask in range(1, (1 << k) - 1):
        S = [i for i in range(k) if mask >> i & 1]
        Sc = [i for i in range(k) if not mask >> i & 1]
        cS = sum(c[i] for i in S)
        Iinc = I - sum(W[i][j] for a_, i in enumerate(Sc) for j in Sc[a_ + 1:])
        const = 0.5 * cS * math.log(C / cS) + cS / C * I - Iinc
        # A(S) = const - sum_{S} c_n y_n <= -t   <=>   -sum_S c_n y_n + t <= -const
        row = [0.0] * (k + 1)
        for i in S: row[i] = -float(c[i])
        row[k] = 1.0
        rows.append(row); rhs.append(-const)
    res = linprog(c=[0.0] * k + [-1.0], A_ub=np.array(rows), b_ub=np.array(rhs),
                  A_eq=np.array([[float(x) for x in c] + [0.0]]), b_eq=[0.0],
                  bounds=[(None, None)] * k + [(None, 50.0)], method='highs')
    if res.status != 0: return (float('nan'), None) if return_y else float('nan')
    if return_y: return res.x[-1], res.x[:k]
    return res.x[-1]

if __name__ == "__main__":
    # 1. Apostol formula vs exact resultants
    bad = 0; npairs = 0
    for M in (60, 105, 45, 72, 75):
        D = [d for d in range(1, M + 1) if M % d == 0]
        for a, b in itertools.combinations(D, 2):
            npairs += 1
            if abs(w_apostol(a, b) - w_exact(a, b)) > 1e-9:
                bad += 1; print("MISMATCH", a, b, w_apostol(a, b), w_exact(a, b))
    print(f"Apostol check: {npairs} pairs, mismatches = {bad}")
    # 2. LP margins in y-variables, compare with R7 values
    tests = [([1, 3], 0.0719), ([1, 19, 57], 0.223), ([1, 3, 5], -0.50),
             ([1, 3, 5, 7, 15, 21, 35, 105], 0.915),
             ([1, 3, 5, 7, 21, 35, 105], 0.018),
             ([1, 5, 7, 35, 45, 63, 105, 315], 0.113),
             ([1, 5, 7, 25, 35, 105, 175, 525], 0.003),
             ([1, 5, 7, 11, 35, 55, 77, 231, 385, 1155], 0.018),
             ([105, 1155], None), ([1, 5, 11, 33, 55, 165], 0.435), ([1, 3, 13, 39, 65, 195], 0.021)]
    for N, r7 in tests:
        t = margin_y(N)
        te = margin_y(N, exact=True) if len(N) <= 8 else float('nan')
        C = sum(int(totient(n)) for n in N); M = math.lcm(*N)
        print(f"N={N}: C={C} M={M} M/C={M/C:.4f} M/phi(M)={M/int(totient(M)):.4f}  margin(y-LP)={t:+.4f}  margin(exact resultants)={te:+.4f}  R7 reported={r7}")
