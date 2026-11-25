import json, sys, glob, math
import numpy as np
def load(fn):
    out=[]
    for l in open(fn):
        l=l.strip()
        if l.startswith('{'):
            try: out.append(json.loads(l))
            except Exception: pass
    return out
def ms(x):
    x=np.array(x,float)
    return "%.4g +- %.2g (sd %.3g, n=%d)"%(x.mean(), x.std(ddof=1)/math.sqrt(len(x)) if len(x)>1 else float('nan'), x.std(ddof=1) if len(x)>1 else float('nan'), len(x))
for fn in sorted(glob.glob('refA_*.jsonl')):
    R=load(fn)
    if not R: continue
    print('==',fn,'N=',R[0]['N'],'c1,c2=',R[0]['c1'],R[0]['c2'])
    print(' v*/N   ',ms([r['vstar_over_N'] for r in R]),' pred 0.7358')
    print(' rho^2  ',ms([r['rho2'] for r in R]),' pred 8')
    print(' V/N    ',ms([r['V_over_N'] for r in R]),' pred 2')
    print(' log(v*/N) ',ms([math.log(r['vstar_over_N']) for r in R]),' pred %.4f'%math.log(2/math.e))
    print(' act sig mean',ms([r['act_sig_mean'] for r in R]),' sd',ms([r['act_sig_sd'] for r in R]))
    print(' nact==N+1 in',sum(r['nact']==r['N']+1 for r in R),'of',len(R),'; dual-certified',sum(r['cert'] for r in R),'; max|dual gap|',max(abs(r['dual_gap']) for r in R))
    print(' samples whose c1-optimum is violated by an offset beyond c1:',sum(r['nviol_beyond_c1']>0 for r in R),'; min slack beyond c1 over samples',min(r['min_slack_beyond_c1'] for r in R))
    print(' max act sigma (shifted)',max(r['act_sig_max'] for r in R),' cutoff1 shifted min',min(r['cutoff1_shifted'] for r in R))
for fn in sorted(glob.glob('refB_*.jsonl')):
    R=load(fn)
    if not R: continue
    print('==',fn,'n=',R[0]['n'],'N=',R[0]['N'],'VmN=',R[0]['VmN'])
    print(' PD in',sum(r['pd'] for r in R),'of',len(R),'; nact==N in',sum(r['nact']==r['N'] for r in R),'; dual-certified',sum(r['cert'] for r in R))
    R=[r for r in R if r['pd']]
    print(' V_LP/N ',ms([r['V_over_N'] for r in R]),' pred 2')
    print(' v_LP/N ',ms([r['v_over_N'] for r in R]),' pred 0.7358')
    print(' |B0|^2 ',ms([r['B0sq'] for r in R]),' pred 4 (4(n+2)/n = %.3f)'%(4*(R[0]['n']+2)/R[0]['n']))
    print(' tr B0^3',ms([r['trB0cube'] for r in R]))
    print(' log(V/v)',ms([math.log(r['V_over_N']/r['v_over_N']) for r in R]),' pred 1')
    print(' lam ratio',ms([r['lam_ratio'] for r in R]))
    print(' missed_mean max',max(r['missed_mean'] for r in R),'; circumscribed ball inside sample in',sum(r['circ_over_N']<r['sampled_over_N'] for r in R),'of',len(R))
    print(' max active W/N',max(r['act_W_max_over_N'] for r in R),'; act logW mean',ms([r['act_logW_mean'] for r in R]),' sd',ms([r['act_logW_sd'] for r in R]))
for fn in sorted(glob.glob('K*_d*.jsonl')):
    R=load(fn)
    if not R: continue
    F=[r for r in R if r['status']=='frozen']
    print('==',fn,'n=',R[0]['n'],'delta=',R[0]['delta'],'runs',len(R),'frozen',len(F))
    if F:
        print(' V_K/N ',ms([r['V_over_N'] for r in F]),'; median %.4g'%np.median([r['V_over_N'] for r in F]))
        print(' npairs==N in',sum(r['npairs']==r['N'] for r in F),'; snapshots',ms([r['n_snapshots'] for r in F]),'; points',ms([r['n_points'] for r in F]))
        print(' |B0|^2',ms([r['B0sq'] for r in F]),'; median lam ratio %.4g'%np.median([r['lam_ratio'] for r in F]),'; min q',min(r['min_q'] for r in F))
