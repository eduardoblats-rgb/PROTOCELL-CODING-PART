import json, numpy as np
from multiprocessing import Pool
from melt_lib import native_pairs, fnat, tm
from seqs import RIBOZYMES
Ts=list(range(0,101,4))
S_LOW=0.015                # inner, Mg-poor: 15 mM monovalent
S_HIGH=1.0                 # outer, Mg-rich: >=3 mM Mg2+ ~ 1 M NaCl (capped at 1.0 M)
S_MG=0.015+333*0.0001      # inner with 0.1 mM Mg2+ (monovalent equivalent)
def job(a):
    name,seq=a
    _,native=native_pairs(seq,S_LOW)
    F={}
    F["out_low"]=[fnat(seq,native,T,S_LOW,0) for T in Ts]
    F["out_high"]=[fnat(seq,native,T,S_HIGH,0) for T in Ts]
    for d in (0.10,0.27):
        F[f"in_mg_{d}"]=[fnat(seq,native,T,S_MG,d) for T in Ts]
    for d in (0.10,0.15,0.27,0.45):
        F[f"in_low_{d}"]=[fnat(seq,native,T,S_LOW,d) for T in Ts]
    return name,F
if __name__=="__main__":
    with Pool(5) as p: res=dict(p.map(job,list(RIBOZYMES.items())))
    json.dump(res,open("fig23.json","w"))
