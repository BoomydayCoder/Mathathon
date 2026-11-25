"""(2)+(4): independent re-verification of every stored final form A_final in data/*.json.
usage: verify_forms.py tag file1 file2 ...   -> writes forms_<tag>.json"""
import numpy as np, json, sys, os, math, time
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vlib

HERE = os.path.dirname(os.path.abspath(__file__))
DISC = {'zeta5': 125, 'zeta7': 7 ** 5, 'zeta8': 256, 'zeta12': 144, 'zeta13': 13 ** 11}


def is_power(fr, e):
    def ip(m):
        if m < 0:
            return False
        r = math.isqrt(m)
        if r * r != m:
            return False
        if e == 2:
            return True
        r2 = math.isqrt(r)
        return r2 * r2 == r
    return ip(fr.numerator) and ip(fr.denominator)


def arith_check(name, k, detA):
    """necessary arithmetic condition for a rational O-Hermitian form on O^k"""
    if name == 'gauss':
        return is_power(detA, 2)
    if name == 'eis':
        return is_power(detA * Fraction(4, 3) ** k, 2)
    if name == 'hur':
        return is_power(detA * 4 ** k, 4)
    if name == 'ico':
        return is_power(detA * Fraction(16, 25) ** (2 * k), 4)
    return is_power(detA / Fraction(DISC[name]) ** k, 2)


def process(fn):
    r = json.load(open(fn))
    name, k, n, w, D = r['name'], r['k'], r['n'], r['w'], r['D']
    els, wm, n0, extra = vlib.base_group(name)
    gsI = [np.kron(np.eye(k, dtype=np.int64), g) for g in els]
    QB = vlib.Qbasis_reynolds(els, k)
    Vb = vlib.Vbest(n)
    kap = vlib.kappa(n)
    res = []
    for i, x in enumerate(r['recs']):
        t0 = time.time()
        A = np.array(x['A_final'])
        has_ext = 'V_ext' in x
        Vst = x['V_ext'] if has_ext else x['V_perfect']
        pst = x['pairs_ext'] if has_ext else x['pairs_perfect']
        q = vlib.find_q(A)
        d = dict(file=os.path.basename(fn), idx=i, n=n, q=q, V_stored=Vst, pairs_stored=pst, has_ext=has_ext)
        if q:
            Z = np.rint(q * A).astype(np.int64)
            d['symmetric'] = bool(np.array_equal(Z, Z.T))
            d['invariant'] = bool(all(np.array_equal(g.T @ Z @ g, Z) for g in gsI))
            detZ = vlib.bareiss_det(Z)
            d['detZ'] = str(detZ)
            V = kap * math.sqrt(q ** n / detZ)
            U = vlib.my_lll(Z)
            Zr = U.T @ Z @ U
            Y, nodes = vlib.enum(Zr.astype(float), q + 0.5)
            X = Y @ U.T
            nrm = np.einsum('ij,jk,ik->i', X, Z, X)
            d['min_exact'] = int(nrm.min())
            d['below'] = int((nrm < q).sum())
            S = X[nrm == q]
            d['arith_ok'] = bool(arith_check(name, k, Fraction(detZ, q ** n)))
            d['even'] = bool(all(int(Z[j, j]) % 2 == 0 for j in range(n)))
        else:
            gsf = [g.astype(float) for g in gsI]
            d['invariant'] = bool(max(np.abs(g.T @ A @ g - A).max() for g in gsf) < 5e-9)
            sgn, ld = np.linalg.slogdet(A)
            V = kap * math.exp(-0.5 * ld)
            U = vlib.my_lll(A)
            Uf = U.astype(float)
            Ar = Uf.T @ A @ Uf
            Y, nodes = vlib.enum((Ar + Ar.T) / 2, 1 + 1e-6)
            X = Y @ U.T
            Xf = X.astype(float)
            nrm = np.einsum('ij,jk,ik->i', Xf, A, Xf)
            d['min_float'] = float(nrm.min())
            d['below'] = int((nrm < 1 - 1e-6).sum())
            S = X
        d['V_mine'] = V
        d['V_relerr'] = abs(V / Vst - 1)
        d['pairs_mine'] = int(len(S))
        d['ratio_best'] = V / Vb
        Sf = S.astype(float)
        Phi = np.einsum('mi,dij,mj->md', Sf, QB, Sf, optimize=True)
        sv = np.linalg.svd(Phi, compute_uv=False)
        d['rank'] = int((sv > 1e-8 * sv[0]).sum())
        d['D'] = D
        # eutaxy-type sanity: nothing
        d['time'] = time.time() - t0
        # famous-lattice certificates
        if q and abs(V / Vb - 1) < 1e-9:
            cert = None
            if n == 8 and q == 2 and detZ == 1 and d['even']:
                cert = 'E8 (even unimodular, rigorous)'
            if n == 24 and q == 4 and detZ == 1 and d['even'] and d['min_exact'] == 4:
                cert = 'Leech (even unimodular, min 4, rigorous)'
            if n in (12, 16) and q == 4 and d['even']:
                p = 3 if n == 12 else 2
                inv = vlib.frac_inverse(Z)
                pm = all((p * inv[a][b]).denominator == 1 for a in range(n) for b in range(n))
                pe = all((p * inv[a][a]).numerator % 2 == 0 for a in range(n))
                if n == 12 and detZ == 729 and pm and pe and d['min_exact'] == 4:
                    cert = 'K12 (even, det 3^6, 3*dual even integral, min 4)'
                if n == 16 and detZ == 256 and pm and pe and d['min_exact'] == 4:
                    cert = 'BW16 (even, det 2^8, 2*dual even integral, min 4)'
            if n == 20 and q == 4 and detZ == 64 and d['min_exact'] == 4:
                cert = 'Lambda20 invariants (det 64, min 4)'
            if n in (2, 4, 6, 10):
                cert = 'record density n=%d' % n
            d['cert'] = cert
        res.append(d)
    return res


if __name__ == '__main__':
    tag = sys.argv[1]
    allres = []
    for fn in sys.argv[2:]:
        t0 = time.time()
        try:
            rs = process(fn)
        except Exception as e:
            print('ERROR', fn, repr(e), flush=True)
            continue
        allres += rs
        nb = sum(1 for d in rs if d['below'] > 0)
        ni = sum(1 for d in rs if not d['invariant'])
        npair = sum(1 for d in rs if d['pairs_mine'] != d['pairs_stored'])
        nr = sum(1 for d in rs if d['rank'] != d['D'])
        nq0 = sum(1 for d in rs if d['q'] == 0)
        na = sum(1 for d in rs if d['q'] and not d['arith_ok'])
        mx = max(d['ratio_best'] for d in rs)
        me = max(d['V_relerr'] for d in rs)
        certs = {}
        for d in rs:
            if d.get('cert'):
                certs[d['cert']] = certs.get(d['cert'], 0) + 1
        print(f"{os.path.basename(fn):28s} recs={len(rs):3d} shorter_vec={nb} noninvariant={ni} pairs_mismatch={npair} "
              f"rank!=D={nr} q_unresolved={nq0} arith_bad={na} maxVrelerr={me:.1e} max V/Vbest={mx:.9f} certs={certs} t={time.time()-t0:.1f}s", flush=True)
        json.dump(allres, open(os.path.join(HERE, f'forms_{tag}.json'), 'w'))
