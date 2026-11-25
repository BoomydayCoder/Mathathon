import sys, time
sys.path.insert(0,'.'); sys.path.insert(0,'../R10_weil_lattices/orchestrator')
import ideal_model as IM, proto as P, numpy as np
for p,N in [(397,12),(1237,12),(3733,12),(421,20),(397,36)]:
    t0=time.time(); Gi,Ei,Bz,F=IM.build(p,N,verbose=False); t1=time.time()
    mx=max(abs(int(x)) for x in Gi.flat)
    G2,U=P.lll_gram_fraction([[int(x) for x in row] for row in Gi]); t2=time.time()
    print(p,N,'build %.1fs'%(t1-t0),'max entry digits',len(str(mx)),'fraction LLL %.1fs'%(t2-t1),'diag',sorted(set(G2[i][i] for i in range(len(G2)))),flush=True)
