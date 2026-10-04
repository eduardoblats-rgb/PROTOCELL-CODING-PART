"""Melting temperatures of (a) full-length ribozyme:complement duplexes and (b) the R3C product-enzyme complex,
in the outer (Mg-rich, 1.0 M equivalent) and inner (0.1 mM Mg = 48 mM, + delta per base pair) phases."""
import RNA, numpy as np, json
from scipy.optimize import brentq
from seqs import RIBOZYMES
import r3c_kin as K
COMP={"A":"U","U":"A","G":"C","C":"G"}
def revcomp(s): return "".join(COMP[c] for c in reversed(s))
def dG_full(seq,Tc,salt,delta):
    md=RNA.md(); md.temperature=float(Tc); md.salt=float(salt)
    comp=revcomp(seq); n=len(seq)
    fc=RNA.fold_compound(seq+"&"+comp,md)
    return fc.eval_structure("("*n+")"*n)+delta*n
def tm_from(dGfun,c):
    f=lambda T: dGfun(T)-0.0019872*(T+273.15)*np.log(c/2)
    try: return brentq(f,5,200)
    except ValueError: return None
res={}
C=1e-6
for name,seq in RIBOZYMES.items():
    out={}
    for lab,salt,d in (("outer (1.0 M eq.)",1.0,0.0),("outer (15 mM, no Mg)",0.015,0.0),("inner (48 mM, d 0.10)",K.salt_in(0.1),0.10),("inner (48 mM, d 0.27)",K.salt_in(0.1),0.27)):
        out[lab]=tm_from(lambda T: dG_full(seq,T,salt,d),C)
    res[name]=out
    print(name,len(seq),{k:(None if v is None else round(v,1)) for k,v in out.items()})
# R3C product-enzyme complex
kobs,kon,c=0.62,1e6,5.0
J=K.calibrate(kobs,kon,c)
cc=c*1e-6
for lab,salt,d in (("outer 1.0 M",1.0,0.0),("inner 48 mM d 0.10",K.salt_in(0.1),0.10),("inner 48 mM d 0.27",K.salt_in(0.1),0.27),("inner 15 mM d 0.27",0.015,0.27)):
    t=tm_from(lambda T: K.S(K.ALL,T,salt,d)+J,cc); print("R3C product complex Tm",lab,None if t is None else round(t,1))
for h in K.HELICES if hasattr(K,"HELICES") else []: pass
from r3c import HELICES
for h in HELICES:
    n=len(HELICES[h][0])
    for lab,salt,d in (("outer 1.0 M",1.0,0.0),("inner 48 mM d 0.27",K.salt_in(0.1),0.27)):
        t=tm_from(lambda T: K.hdG(h,T,salt)+d*n,cc)
        print(h,n,"bp",lab,None if t is None else round(t,1))
