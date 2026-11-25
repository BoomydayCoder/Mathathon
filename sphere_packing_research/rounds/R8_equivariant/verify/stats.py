"""(3) recompute summary statistics from raw data (independent of summarize.py) and compare with results_tables.md"""
import json, glob, os, math, sys, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vlib
HERE = os.path.dirname(os.path.abspath(__file__))
R8 = os.path.dirname(HERE)
rows = {}
tot_runs = tot_fail = tot_steps = tot_nonedge = 0
maxratio = 0
nfiles = 0
bad_verify = 0
desc_err = 0
for f in sorted(glob.glob(os.path.join(R8, 'data', '*.json'))):
    r = json.load(open(f))
    recs = r['recs']
    nfiles += 1
    tot_fail += r['fails']
    if not recs:
        continue
    n, w, D = r['n'], r['w'], r['D']
    Vb = vlib.Vbest(n)
    Vp = np.array([x['V_perfect'] for x in recs])
    Ve = np.array([x.get('V_ext', np.nan) for x in recs])
    op = np.array([x['pairs_perfect'] / (w / 2) for x in recs])
    tot_runs += len(recs)
    tot_steps += sum(x.get('descent_steps', 0) for x in recs)
    tot_nonedge += sum(x.get('nonedge', 0) for x in recs)
    desc_err += sum(1 for x in recs if 'descent_error' in x)
    for x in recs:
        vp = x['verify_perfect']
        if not (vp[0] == 0 and vp[1] == x['pairs_perfect'] and abs(vp[2] - 1) < 1e-8):
            bad_verify += 1
        if 'verify_ext' in x:
            ve = x['verify_ext']
            if not (ve[0] == 0 and ve[1] == x['pairs_ext'] and abs(ve[2] - 1) < 1e-8):
                bad_verify += 1
        if 'V_ext' in x and x['V_ext'] < x['V_perfect'] * (1 - 1e-9):
            print('WARNING V_ext < V_perfect', f, x['V_ext'], x['V_perfect'])
    maxratio = max(maxratio, np.nanmax(np.concatenate([Vp, Ve])) / Vb)
    key = (r['name'], r['k'], r['mode'])
    assert key not in rows or True
    d = dict(fam=r['name'], k=r['k'], n=n, w=w, D=D, mode=r['mode'], runs=len(recs), fails=r['fails'],
             Vp_min=Vp.min(), Vp_med=float(np.median(Vp)), Vp_mean=Vp.mean(), Vp_max=Vp.max(),
             Ve_med=float(np.nanmedian(Ve)), Ve_max=float(np.nanmax(Ve)), orb_mean=op.mean(), orb_min=op.min(), orb_max=op.max(),
             simple=float(np.mean(np.abs(op - D) < 1e-9)), ratio_med=float(np.median(Vp)) / (w * D / 2),
             hitp_n=int((np.abs(Vp / Vb - 1) < 1e-6).sum()), hite_n=int((np.abs(Ve / Vb - 1) < 1e-6).sum()),
             n_ext=int(np.isfinite(Ve).sum()), file=os.path.basename(f))
    if key in rows:   # several seeds: merge is not done by summarize.py; keep separately
        key = key + (r['seed'],)
    rows[key] = d

print('files', nfiles, 'runs', tot_runs, 'fails', tot_fail, 'descent errors', desc_err, 'LP steps', tot_steps, 'nonedge', tot_nonedge,
      'records failing stored verify', bad_verify, 'max V/Vbest %.9f' % maxratio)
print()
hdr = '| family | k | n | w | D | mode | runs | V_perf min/med/mean/max | V_ext med/max | orbits mean [min,max] | frac D | med Vp/(wD/2) | hits perf | hits ext (of n_ext) |'
print(hdr)
for key in sorted(rows):
    d = rows[key]
    print(f"| {d['fam']} | {d['k']} | {d['n']} | {d['w']} | {d['D']} | {d['mode']} | {d['runs']} | {d['Vp_min']:.2f} / {d['Vp_med']:.2f} / {d['Vp_mean']:.2f} / {d['Vp_max']:.2f} | "
          f"{d['Ve_med']:.2f} / {d['Ve_max']:.2f} | {d['orb_mean']:.1f} [{d['orb_min']:.0f},{d['orb_max']:.0f}] | {d['simple']:.2f} | {d['ratio_med']:.3f} | "
          f"{d['hitp_n']}/{d['runs']} | {d['hite_n']}/{d['n_ext']} | {d['file']}")

# compare with their table
their = {}
for line in open(os.path.join(R8, 'results_tables.md')):
    c = [s.strip() for s in line.strip().strip('|').split('|')]
    if len(c) == 20 and c[0] not in ('family', '---'):
        their[(c[0], int(c[1]), c[5], int(c[6]))] = c
mism = 0
matched = 0
for key, d in rows.items():
    tk = (d['fam'], d['k'], d['mode'], d['runs'])
    if tk not in their:
        print('NOT IN THEIR TABLE:', tk)
        mism += 1
        continue
    c = their[tk]
    mine = [f"{d['Vp_min']:.2f} / {d['Vp_med']:.2f} / {d['Vp_mean']:.2f} / {d['Vp_max']:.2f}", f"{d['Ve_med']:.2f} / {d['Ve_max']:.2f}",
            f"{d['orb_mean']:.1f} [{d['orb_min']:.0f},{d['orb_max']:.0f}]", f"{d['simple']:.2f}", f"{d['ratio_med']:.3f}",
            f"{d['hitp_n']/d['runs']:.2f} / {d['hite_n']/max(d['n_ext'],1):.2f}"]
    th = [c[8], c[9], c[10], c[11], c[12], c[15]]
    if mine != th:
        mism += 1
        print('MISMATCH', tk, mine, th)
    else:
        matched += 1
print('table rows matched:', matched, 'mismatched:', mism, 'rows in their table:', len(their))

# scaling fit
sel = [d for d in rows.values() if d['mode'] == 'random' and d['n'] >= 10]
x = np.log([d['w'] * d['D'] / 2 for d in sel])
y = np.log([d['Vp_med'] for d in sel])
p = np.polyfit(x, y, 1)
res = y - np.polyval(p, x)
rat = np.array([d['ratio_med'] for d in sel])
print('scaling fit: points', len(sel), 'slope %.3f intercept %.3f rms %.3f' % (p[0], p[1], res.std()),
      'ratio min %.3f geomean %.3f max %.3f' % (rat.min(), math.exp(np.log(rat).mean()), rat.max()),
      'wD/2 range', min(d['w'] * d['D'] / 2 for d in sel), max(d['w'] * d['D'] / 2 for d in sel))
# bootstrap-ish uncertainty of slope: leave-one-out
sl = [np.polyfit(np.delete(x, i), np.delete(y, i), 1)[0] for i in range(len(x))]
print('leave-one-out slope range %.3f .. %.3f' % (min(sl), max(sl)))
# fit using means instead of medians
ym = np.log([d['Vp_mean'] for d in sel])
pm = np.polyfit(x, ym, 1)
print('fit with MEAN V_perf: slope %.3f; ratio mean V_perf/(wD/2): min %.3f geomean %.3f max %.3f' % (
    pm[0], min(d['Vp_mean'] / (d['w'] * d['D'] / 2) for d in sel),
    math.exp(np.mean([math.log(d['Vp_mean'] / (d['w'] * d['D'] / 2)) for d in sel])), max(d['Vp_mean'] / (d['w'] * d['D'] / 2) for d in sel)))
for d in sorted(sel, key=lambda d: d['w'] * d['D']):
    print('   %-7s k=%d n=%2d wD/2=%5d runs=%2d med Vp=%8.2f ratio=%.3f  medVp/N=%.3f' % (d['fam'], d['k'], d['n'], d['w'] * d['D'] // 2, d['runs'], d['Vp_med'], d['ratio_med'], d['Vp_med'] / (d['n'] * (d['n'] + 1) / 2)))
ye = np.log([d['Ve_med'] for d in sel]) - x
print('V_ext/(wD/2): geomean %.3f max %.3f' % (math.exp(ye.mean()), math.exp(ye.max())))

# Leech trajectories (icosian k=3)
for f in sorted(glob.glob(os.path.join(R8, 'data', 'ico_k3_*.json'))):
    r = json.load(open(f))
    Vb = vlib.Vbest(24)
    for i, x in enumerate(r['recs']):
        hp = abs(x['V_perfect'] / Vb - 1) < 1e-6
        he = abs(x.get('V_ext', 0) / Vb - 1) < 1e-6
        if hp or he:
            print(os.path.basename(f), 'run', i, 'Leech at perfect' if hp else 'Leech after descent', 'pairs_perf', x['pairs_perfect'], 'pairs_ext', x.get('pairs_ext'),
                  'last traj', x['traj'][-2:], 'descent_steps', x.get('descent_steps'), 'V_perf %.2f' % x['V_perfect'])
