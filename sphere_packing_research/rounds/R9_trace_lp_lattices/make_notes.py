#!/usr/bin/env python3
"""Compose notes.md = prose + compact tables generated from summary.json / summary_lp_pooled.json."""
import json, os, math
here = os.path.dirname(os.path.abspath(__file__))
S = json.load(open(os.path.join(here, 'summary.json')))
Pp = json.load(open(os.path.join(here, 'summary_lp_pooled.json')))
rows = sorted(Pp.values(), key=lambda r: (r['n'], r['kind']))
T1 = ['| n | family | #lat | Vbest/N | V_LP/N mean+-se | std | Poisson std 2/sqrt(N) | v_LP/N | pairs/N | simple vertices |', '|---|---|---|---|---|---|---|---|---|---|']
for r in rows:
    N = r['n'] * (r['n'] + 1) // 2
    se = ('%.3f' % r['VLP_se']) if r['VLP_se'] is not None else '-'
    sd = ('%.3f' % r['VLP_std']) if r.get('VLP_std') is not None else '-'
    T1.append(f"| {r['n']} | {r['kind']} | {r['nlat']} | {r['Vbest_over_N']:.3f} | {r['VLP_mean']:.3f}+-{se} | {sd} | {2/math.sqrt(N):.3f} | {r['vLP_mean']:.3f} | {r['pairs_over_N_mean']:.3f} | {r['nsimple']}/{r['nlat']} |")
T2 = ['| n | family | V_LP/N | V_K/N random sign [#] | V_Kd/N descent sign [#] | V_ext/N from LP [#] | V_ext/N from walk [#] | best known reached: LP / walk / greedy-LP / greedy-walk |', '|---|---|---|---|---|---|---|---|']
for k, s in sorted(S.items(), key=lambda t: (t[1]['n'], t[1]['kind'])):
    T2.append(f"| {s['n']} | {s['kind']} | {s['VLP']} | {s['VK']} [{s['nK']}] | {s['VKd']} [{s['nKd']}] | {s['VextLP']} [{s['nextLP']}] | {s['VextK']} [{s['nextK']}] | {s['fLP']} / {s['fK']} / {s['fextLP']} / {s['fextK']} |")
prose = open(os.path.join(here, 'notes_prose.md')).read()
txt = prose.replace('@@TABLE_LP@@', '\n'.join(T1)).replace('@@TABLE_CMP@@', '\n'.join(T2))
open(os.path.join(here, 'notes.md'), 'w').write(txt)
print('notes.md written', len(txt))
