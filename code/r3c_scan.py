import itertools, json, numpy as np
from multiprocessing import Pool
from r3c_kin import *
MG=[0.0,0.1,1.0]; PHI=[0.1,0.2,0.4]; DEL=[0.10,0.27]; KON=[1e5,1e6,1e7]; KOBS=[0.3,0.62,1.3]; CONC=[1.0,5.0,25.0]; KDIN=[10.0,100.0,3600.0]
TS=list(range(22,47,2))
def job(p):
    mg,phi,d,kon,kobs,c,kdi=p
    J=calibrate(kobs,kon,c)
    l1=lam_one(42,J,kobs,kon,c); l2=lam_two(42,J,kobs,kon,c,phi,mg,d,kdi)
    # lowest temperature at which the two-phase system still grows at least as fast as one-phase does at 42 C
    Tlow=None
    for T in TS:
        if lam_two(T,J,kobs,kon,c,phi,mg,d,kdi)>=l1: Tlow=T; break
    return dict(p=p,J=J,l1=l1,l2=l2,speed=l2/l1,Tlow=Tlow)
if __name__=="__main__":
    combos=list(itertools.product(MG,PHI,DEL,KON,KOBS,CONC,KDIN))
    with Pool(5) as pool: res=pool.map(job,combos,chunksize=20)
    json.dump(res,open("r3c_scan.json","w"))
    sp=np.array([r["speed"] for r in res]); print(len(res),"combos")
    print("median",np.median(sp),"IQR",np.percentile(sp,[25,75]),"min",sp.min(),"max",sp.max())
    print(">1.5x:",(sp>1.5).mean(),"  <1x:",(sp<1).mean())
    for mg in MG:
        s=np.array([r["speed"] for r in res if r["p"][0]==mg]); print("Mg_in",mg,"mM  median",round(np.median(s),2),"range",round(s.min(),2),round(s.max(),1))
    for k in KDIN:
        s=np.array([r["speed"] for r in res if r["p"][6]==k and r["p"][0]==0.1]); print("kdin_max",k,"(Mg_in 0.1) median",round(np.median(s),2))
