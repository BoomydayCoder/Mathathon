# Monte Carlo sanity check of Lemma 3.4 (angular lemma, j x k quaternionic blocks, s = 2 places)
# and of E CV^2 ~ (j - 1/2)/k * (sj... ) ; Leech: j = 3, d = 4, s = 2.
import numpy as np
rng = np.random.default_rng(1)
def cv2_samples(j, k, s, nsamp):
    mus = []
    for v in range(s):
        A = rng.standard_normal((nsamp, j, k)) + 1j*rng.standard_normal((nsamp, j, k))
        B = rng.standard_normal((nsamp, j, k)) + 1j*rng.standard_normal((nsamp, j, k))
        top = np.concatenate([A, B], axis=2)
        bot = np.concatenate([-B.conj(), A.conj()], axis=2)
        X = np.concatenate([top, bot], axis=1)            # 2j x 2k complex rep of y_v
        S = X @ X.conj().transpose(0, 2, 1)
        ev = np.linalg.eigvalsh(S)                         # each eigenvalue of y y^* twice
        mus.append(ev[:, ::2])
    mu = np.concatenate(mus, axis=1)                       # nsamp x (s j)
    t = (s*j) * mu / mu.sum(axis=1, keepdims=True)
    return ((t - 1.0)**2).mean(axis=1)
for (j, k) in [(3, 50), (3, 200), (3, 800), (1, 200), (2, 200)]:
    s, d = 2, 4
    cv = cv2_samples(j, k, s, 20000)
    n = 24*k
    line = f"j={j} k={k}: E CV^2 = {cv.mean():.5f}  (k*E = {k*cv.mean():.3f}), min CV^2 = {cv.min():.5f}"
    for b in [2.0, 2*np.log(n)/3*j/1.0/ (3.0/ j) if False else j/3*2*np.log(n)]:
        kappa = (d*(j-1) + 4*b + 2)/(d*k)
        bound = (1-kappa)**(-2*b*s*j) if kappa < 1 else float('inf')
        line += f" | b={b:.2f}: E exp(b CV^2) = {np.exp(b*cv).mean():.4f} <= bound {bound:.4f}"
    print(line)
