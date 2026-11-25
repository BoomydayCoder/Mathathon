import sys, collections
sys.path.insert(0,'.'); sys.path.insert(0,'../R10_weil_lattices/orchestrator')
import ideal_model as IM, proto as P, numpy as np
V = IM.V
for p,N in [(5,4),(13,12),(37,12),(61,12),(61,20)]:
    out = {}
    for var in ('L1','L0'):
        Gi,Ei,Bz,F = IM.build(p,N,verbose=False,variant=var)
        n = 2*N
        G1,U1 = P.lll_gram_fraction([[int(x) for x in row] for row in Gi]); U1 = np.array(U1,dtype=object).T
        Gr = U1.T@Gi@U1; Er = U1.T@Ei@U1
        det = IM.exact_det(Gr)
        Xs,vals,flag = V.short_vectors(np.array(Gr,dtype=float), 4.01, maxc=600000, maxnodes=2_000_000_000)
        cnt = collections.Counter()
        for x in np.array(Xs,dtype=object):
            nrm = int(x@Gr@x); s = int(x@Er@x)//2
            cnt[(nrm,abs(s))] += 1
        out[var] = (det, all(int(Gr[i,i])%2==0 for i in range(n)), tuple(sorted(cnt.items())))
    print(p, N, 'L1:', out['L1'][0], out['L1'][1], ' L0:', out['L0'][0], out['L0'][1], ' same (norm,|s|) table:', out['L1'][2]==out['L0'][2], flush=True)
    if out['L1'][2]!=out['L0'][2]: print('   L1', out['L1'][2][:8]); print('   L0', out['L0'][2][:8])
