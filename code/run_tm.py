import json, numpy as np
from multiprocessing import Pool
from melt_lib import native_pairs, fnat, tm
from seqs import RIBOZYMES
SALT=0.015; Ts=list(range(0,101,4)); DELTA=0.27
def job(a):
    name,seq=a
    s,native=native_pairs(seq,SALT)
    out={}
    for d in (0.0,DELTA):
        out[d]=[fnat(seq,native,T,SALT,d) for T in Ts]
    return name,len(seq),len(native),out
if __name__=="__main__":
    with Pool(5) as p: res=p.map(job,list(RIBOZYMES.items()))
    json.dump({n:{"len":l,"npairs":k,"F":{str(d):F for d,F in o.items()}} for n,l,k,o in res},open("tm_result.json","w"))
    for n,l,k,o in res:
        t0=tm(Ts,o[0.0]); t1=tm(Ts,o[DELTA])
        print(f"{n:42s} pairs {k:3d}  Tm outer {t0}  Tm inner {t1}  dTm {None if t0 is None or t1 is None else round(t0-t1,1)}  F25 out {o[0.0][6]:.2f} in {o[DELTA][6]:.2f}")
