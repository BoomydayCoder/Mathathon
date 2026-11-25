# Check Lemma 0.3: P(x (x) x) onto Q = {R_N : N in Herm_k(H)} has |P|^2/|x|^4 = (1/(d r)) sum lambda^2/(sum lambda)^2,
# lambda = eigenvalues of x x^* (r x r quaternionic Hermitian). One place, K = H, r = 3, k = 2 (E = R^24), and dim Q = 2k^2-k.
import numpy as np
rng = np.random.default_rng(0)
r, k = 3, 2
def chi(q):  # q: array (..., 4) real quaternion coords (a0,a1,a2,a3) -> 2x2 complex
    a = q[..., 0] + 1j*q[..., 1]; b = q[..., 2] + 1j*q[..., 3]
    return np.stack([np.stack([a, b], -1), np.stack([-b.conj(), a.conj()], -1)], -2)
def to_c(xr, rows, cols):  # real vector -> (2rows x 2cols) complex matrix
    q = xr.reshape(rows, cols, 4); c = chi(q)             # rows, cols, 2, 2
    return c.transpose(0, 2, 1, 3).reshape(2*rows, 2*cols)
def to_r(C, rows, cols):
    c = C.reshape(rows, 2, cols, 2).transpose(0, 2, 1, 3)
    a = c[..., 0, 0]; b = c[..., 0, 1]
    return np.stack([a.real, a.imag, b.real, b.imag], -1).reshape(-1)
n = 4*r*k
def R(Nr):  # real matrix of x -> x N
    Nc = to_c(Nr, k, k); M = np.zeros((n, n))
    for i in range(n):
        e = np.zeros(n); e[i] = 1
        M[:, i] = to_r(to_c(e, r, k) @ Nc, r, k)
    return M
# span of Hermitian N: symmetrise random N
basis = []
for _ in range(40):
    Nr = rng.standard_normal(4*k*k); Nc = to_c(Nr, k, k); Nh = (Nc + Nc.conj().T)/2
    M = R(to_r(Nh, k, k)); basis.append(M.reshape(-1))
B = np.array(basis); U, S, Vt = np.linalg.svd(B, full_matrices=False)
rank = int((S > 1e-9).sum()); Qb = Vt[:rank]
Nh0 = (lambda c: (c + c.conj().T)/2)(to_c(rng.standard_normal(4*k*k), k, k)); M0 = R(to_r(Nh0, k, k))
print("dim Q =", rank, " expected 2k^2-k =", 2*k*k-k, "; R_N symmetric for Hermitian N:", np.allclose(M0, M0.T),
      "; Tr R_N = d r Re tr N:", np.isclose(np.trace(M0), 4*r*np.trace(Nh0).real/2))
for trial in range(4):
    x = rng.standard_normal(n)
    if trial == 3:  # rank-one point (a "line" point): column * row
        col = to_c(rng.standard_normal(4*r), r, 1); row = to_c(rng.standard_normal(4*k), 1, k); x = to_r(col @ row, r, k)
    xx = np.outer(x, x).reshape(-1); P = Qb.T @ (Qb @ xx)
    q = (P @ P)/(x @ x)**2
    xc = to_c(x, r, k); lam = np.linalg.eigvalsh(xc @ xc.conj().T)[::2]
    pred = (1/(4*r))*(lam**2).sum()/lam.sum()**2
    print(f"trial {trial}: q = {q:.10f}  formula = {pred:.10f}  q0(s=1) = {1/(4*r*r):.6f}  q/q0 = {q*4*r*r:.4f}")
