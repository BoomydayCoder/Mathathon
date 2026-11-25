import eq, numpy as np, time, sys
specs = [('eis',1,3),('gauss',1,3),('hur',1,3),('ico',1,5),('hur',2,10),('eis',3,10),('eis',4,10),('eis',6,3),('hur',4,2),('ico',2,2)]
for name,k,nr in specs:
    t0=time.time()
    r=eq.run_family(name,k,nr,40,'random',1,outdir='data_test')
    Vs=sorted(set(round(x['V_perfect'],4) for x in r['recs'])); Ve=sorted(set(round(x.get('V_ext',-1),4) for x in r['recs']))
    print('   Vperf set',Vs[-6:],'Vext set',Ve[-6:], 'verify', [x.get('verify_ext') for x in r['recs'][:2]], 'time/run %.2f'%((time.time()-t0)/max(1,len(r['recs']))), flush=True)
