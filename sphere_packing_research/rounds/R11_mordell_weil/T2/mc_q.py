import numpy as np
rng = np.random.default_rng(1)
def quat_to_complex(a, b):  # x = a + b j, a,b complex (k x r) -> complex 2k x 2r
    return np.block([[a, b], [-b.conj(), a.conj()]])
def sample_ratio(r, k, N=400):
    out = []
    for _ in range(N):
        num = 0.0; den = 0.0
        for v in range(2):
            a = rng.normal(size=(k, r)) + 1j*rng.normal(size=(k, r))
            b = rng.normal(size=(k, r)) + 1j*rng.normal(size=(k, r))
            X = quat_to_complex(a, b)
            S = X.conj().T @ X
            num += 0.5*np.trace(S @ S).real
            den += 0.5*np.trace(S).real
        out.append(2*r*num/den**2)
    out = np.array(out)
    return out.mean(), out.std(), out.min()
for (r, k) in [(3, 4), (3, 12), (9, 10), (9, 27), (15, 16), (15, 150), (39, 40)]:
    m, s, mn = sample_ratio(r, k)
    print(f"r={r} k={k}: E q/q0 = {m:.4f} (pred 1+(r-1/2)/k = {1+(r-0.5)/k:.4f}), std {s:.4f}, min {mn:.4f}; trace part r/k = {r/k:.4f}; traceless mean-1 = {m-1-r/k:+.4f}")
