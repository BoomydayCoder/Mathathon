# Generates best_known_lattices.csv and prints derived tables.  Own script; reads no downloaded file.
import math, csv
s2, s3 = math.sqrt(2), math.sqrt(3)
MEM = "memory: SPLAG Table 1.2/6.1 (Nebe-Sloane site unreachable from this session)"
CE  = "; cross-checked numerically with Cohn-Elkies math/0110009 Table 3 (best packing known 2001)"
KAL = "; 1/delta cross-checked with Kallus arXiv:1305.1310 Table I"
rows = {}
def add(n, d, name, kiss, src): rows[n] = (d, name, kiss, src)
# n = 1..24 : laminated / K_n lattices (exact closed forms)
lam = {1:(1/2,"Z=Lambda1",2),2:(1/(2*s3),"A2=Lambda2",6),3:(1/(4*s2),"A3=D3=Lambda3",12),4:(1/8,"D4=Lambda4",24),
5:(1/(8*s2),"D5=Lambda5",40),6:(1/(8*s3),"E6=Lambda6",72),7:(1/16,"E7=Lambda7",126),8:(1/16,"E8=Lambda8",240),
9:(1/(16*s2),"Lambda9",272),10:(1/(16*s3),"Lambda10",336),11:(1/(18*s3),"K11",432),12:(1/27,"K12 (Coxeter-Todd)",756),
13:(1/(18*s3),"K13",918),14:(1/(16*s3),"Lambda14",1422),15:(1/(16*s2),"Lambda15",2340),16:(1/16,"Lambda16=BW16",4320),
17:(1/16,"Lambda17",5346),18:(1/(8*s3),"Lambda18",7398),19:(1/(8*s2),"Lambda19",10668),20:(1/8,"Lambda20",17400),
21:(1/(4*s2),"Lambda21",27720),22:(1/(2*s3),"Lambda22",49896),23:(1/2,"Lambda23",93150),24:(1.0,"Lambda24 (Leech)",196560)}
nonlat = {10:"0.03906 (P10c)",11:"0.03516 (P11a)",13:"0.03516 (P13a)",18:"0.07508",20:"0.13154",22:"0.33254"}
for n,(d,name,k) in lam.items():
    src = MEM + CE
    if 9 <= n <= 20: src += KAL
    if n <= 8: src += "; PROVEN optimal lattice (classical)"
    if n == 9: src += "; PROVEN optimal lattice (Dutour Sikiric-van Woerden arXiv:2508.20719)"
    if n == 24: src += "; PROVEN optimal lattice (Cohn-Kumar)"
    if n in nonlat: src += "; NOTE a denser NON-lattice packing is known: delta=" + nonlat[n]
    add(n, d, name, k, src)
# n = 25..36
add(25, 1/s2, "Lambda25", 196656, MEM+CE)
add(26, 1/s3, "Lambda26 (also T26)", 196848, MEM+CE)
add(27, 1/s2, "B27 (Bachoc) [name from memory]", "", MEM+CE+"; lattice attribution from memory (laminated Lambda27 has delta=0.5, kissing 197142)")
add(28, 1.0, "B28 (Bachoc) [name from memory]", "", MEM+CE+"; lattice attribution from memory (Lambda28: 0.5, kissing 197736)")
add(29, 1/s2, "B29 [name from memory]", "", MEM+CE+"; lattice attribution from memory (Lambda29: 0.5, kissing 198506)")
add(30, 1.0, "Q30 (Quebbemann) [name from memory]", "", MEM+CE+"; lattice attribution from memory (Lambda30: 0.57735, kissing 200046)")
add(31, 3**15/2**23.5, "Q31 (Quebbemann)", "", MEM+CE+"; closed form 3^15/2^23.5 inferred (matches 1.2095)")
add(32, 3**16/2**24, "Q32 (Quebbemann)", 261120, MEM+CE+"; closed form 3^16/2^24 (BW32 and Lambda32 have delta=1)")
add(33, 3**16.5/2**25, "Q33 [name from memory]", "", MEM+CE+"; closed form 3^16.5/2^25 inferred (matches 2.2220)")
add(34, 3**16.5/2**25, "Q34 [name from memory]", "", MEM+CE+"; value 2.2220")
add(35, 2*s2, "B35 [name from memory]", "", MEM+CE)
add(36, 2**18/3**10, "Ks36 (Kschischang-Pasupathy)", 234456, MEM+CE+"; closed form 2^18/3^10; kissing number from memory")
# n = 37..47 : laminated lattices via lambda_{24+k} = 2^{-k} lambda_k  (delta = 2^{k/2} delta(Lambda_k))
for n in range(37, 48):
    k = n - 24
    lamdelta = {11: 1/32, 12: 1/32, 13: 1/32}.get(k, lam[k][0])   # Lambda11,12,13 have delta=1/32 (K11,K12,K13 are denser, not laminated)
    d = 2**(k/2) * lamdelta
    add(n, d, "Lambda%d (laminated)" % n, "",
        "LOWER BOUND ONLY, from memory/derived: laminated lattice rule delta(Lambda_{24+k})=2^{k/2} delta(Lambda_k) (SPLAG Ch.6 Table 6.1; anchors Lambda32=1, Lambda40=16, Lambda48=4096 checked); the Nebe-Sloane table may list a denser lattice in this dimension (could not be fetched)")
add(48, 1.5**24, "P48p, P48q, P48n (extremal even unimodular, min 6)", 52416000, "exact from definition: (sqrt6/2)^48 = (3/2)^24; names/kissing from memory (SPLAG)")
add(56, 1.5**28, "extremal even unimodular, min 6 (e.g. Nebe L56,2(M))", 15590400, "exact from (min 6, det 1); existence and kissing number from memory")
add(64, 3.0**16, "Ne64 (Nebe, extremal 3-modular, min 12)", "", "memory (Nebe-Sloane table): delta=3^16; BW64 has delta=2^16, kissing 9694080")
add(72, 2.0**36, "Gamma72 (Nebe, extremal even unimodular, min 8)", 6218175600, "exact from (min 8, det 1); existence Nebe 2012; kissing from memory")
add(80, 2**40.14, "MW80 (Mordell-Weil; Shioda/Elkies) [log2 delta=40.14 from memory]", "", "memory (SPLAG Table 1.3) log2 delta = 40.14, LOW PRECISION; rigorous alternative: extremal even unimodular L80, M80 (Bachoc-Nebe) have delta = 2^40 exactly, kissing 1250172000")
add(128, 2**97.40, "MW128 (Mordell-Weil; Elkies) [log2 delta=97.40 from memory]", 218044170240, "memory (SPLAG Table 1.3) log2 delta = 97.40, LOW PRECISION; BW128 has delta = 2^64 exactly, kissing 1260230400")
def ballvol(n): return math.pi**(n/2)/math.gamma(n/2+1)
with open("best_known_lattices.csv","w",newline="") as f:
    w = csv.writer(f)
    w.writerow(["n","center_density","log2_center_density","name","kissing_number","source","packing_density","V_eq_2^n_density","V_over_n","V_over_n2","V_over_n3"])
    for n in sorted(rows):
        d,name,k,src = rows[n]
        dens = d*ballvol(n); V = dens*2.0**n
        w.writerow([n, "%.10g"%d, "%.6f"%math.log2(d), name, k, src, "%.6e"%dens, "%.6g"%V, "%.5g"%(V/n), "%.5g"%(V/n**2), "%.5g"%(V/n**3)])
print("n  delta  V=2^n*density  V/n  V/n^2  V/n^3")
for n in sorted(rows):
    d = rows[n][0]; V = d*ballvol(n)*2.0**n
    print("%3d %12.6g %12.6g %9.4g %9.4g %9.4g" % (n, d, V, V/n, V/n**2, V/n**3))
# Andreanov-Scardicchio Table II: <e> = -(1/d) <ln phi> of random perfect lattices (random walk on Voronoi graph), best known
AS = {8:(0.180572,0.171465),9:(0.23352,0.21396),10:(0.26828,0.23857),11:(0.29347,0.25511),12:(0.31505,0.25041),13:(0.33106,0.27178),
14:(0.34277,0.27386),15:(0.35405,0.27218),16:(0.36507,0.26370),17:(0.37205,0.27833),18:(0.38322,0.28489),19:(0.39000,0.28903)}
print("\nAndreanov-Scardicchio typical perfect lattice: d, <e>, phi_typ=exp(-d<e>), V_typ=2^d phi_typ, V_typ/d^2, V_typ/d^3, V_best(from e_best)")
xs=[];ys=[]
for d,(e,eb) in AS.items():
    V = 2.0**d*math.exp(-d*e); Vb = 2.0**d*math.exp(-d*eb)
    print("%2d %.5f %.4e %8.2f %7.4f %8.5f %9.1f" % (d,e,math.exp(-d*e),V,V/d**2,V/d**3,Vb))
    xs.append(math.log(d)); ys.append(math.log(V))
def fit(xs,ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    b=sum((x-mx)*(y-my) for x,y in zip(xs,ys))/sum((x-mx)**2 for x in xs); return b, math.exp(my-b*mx)
print("power-law fit V_typ = C d^nu, d=8..19:  nu=%.3f C=%.4f" % fit(xs,ys))
print("power-law fit V_typ = C d^nu, d=12..19: nu=%.3f C=%.4f" % fit(xs[4:],ys[4:]))
print("power-law fit V_typ = C d^nu, d=14..19: nu=%.3f C=%.4f" % fit(xs[6:],ys[6:]))
print("\nKallus-Marcotte-Torquato fit for isostatic extreme (jammed) lattices: V = 2^d<phi> = 0.0567 d^3.01, d=15..24")
for d in (15,19,20,24,32,48,64,100): print(d, "%.1f"%(0.0567*d**3.01), " ratio to AS typical perfect:" , ("%.2f"%(0.0567*d**3.01/(2.0**d*math.exp(-d*AS[d][0]))) if d in AS else "-"))
