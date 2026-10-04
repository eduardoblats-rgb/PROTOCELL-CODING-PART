import json, numpy as np
from multiprocessing import Pool
from population import *
M=200; SEEDS=[11,22,33]
SIZES=[10,20,40,80,160]; ERRS=[0.5,1,2,4,8,12,16,24]
def comp_task(a):
    N,e,seed=a; rng=np.random.default_rng(seed)
    ts,ff,_=run_compart(M,N,np.full(M,kd_for_ratio(4)),None,e/100,2.0,200,0.1,rng)
    return float(ff[(ts>=150)].mean())
def mix_task(a):
    N,e,seed=a; rng=np.random.default_rng(seed)
    ts,ff=run_mixed(M,N,kd_for_ratio(4),e/100,2.0,200,0.1,rng)
    return float(ff[(ts>=150)].mean())
def fig5():
    out={}
    kd4,kd13=kd_for_ratio(4),kd_for_ratio(1.3)
    for lab,kdx in (("4x",kd4),("1.3x",kd13)):
        curves=[]
        for s in SEEDS:
            rng=np.random.default_rng(s); kdt=np.array([KD1]*(M//2)+[kdx]*(M//2))
            ts,ff,ft=run_compart(M,20,kdt,None,0.0,1.0,40,0.5,rng,rec_every=1.0); curves.append(ft)
        out[lab]=dict(t=ts.tolist(),f=np.mean(curves,0).tolist())
    pc=[];pm=[]
    for s in SEEDS:
        rng=np.random.default_rng(s); ts,ff,_=run_compart(M,20,np.full(M,kd4),None,0.02,1.5,120,0.1,rng); pc.append(ff)
        rng=np.random.default_rng(s); ts,fm=run_mixed(M,20,kd4,0.02,1.5,120,0.1,rng); pm.append(fm)
    out["par"]=dict(t=ts.tolist(),comp=np.mean(pc,0).tolist(),mixed=np.mean(pm,0).tolist())
    return out
if __name__=="__main__":
    with Pool(5) as p:
        tasks=[(N,e,s) for N in SIZES for e in ERRS for s in SEEDS]
        r=p.map(comp_task,tasks)
        H=np.array(r).reshape(len(SIZES),len(ERRS),len(SEEDS)).mean(2)
        mt=[(40,e,s) for e in ERRS for s in SEEDS]
        rm=np.array(p.map(mix_task,mt)).reshape(len(ERRS),len(SEEDS)).mean(1)
    res=dict(fig5=fig5(),sizes=SIZES,errs=ERRS,heat=H.tolist(),mixed40=rm.tolist())
    json.dump(res,open("fig56.json","w"))
    print((100*H).round(0)); print((100*rm).round(1))
