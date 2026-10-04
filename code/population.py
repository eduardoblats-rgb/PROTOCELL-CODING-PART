"""Stochastic-corrector model of RNA protocells (tau-leaping + Moran replacement).
Per protocell: Ef free replicases, Pf free parasites, CE / CP replicase-bound copies
(functional / parasite). Copying -> bound complex -> release (rate kd) -> free molecules.
"""
import numpy as np
from scipy.optimize import brentq
A=3.0      # catalysis rate per free replicase (1/h)
KD1=0.5    # product-release rate, one phase (1/h)
DEG=0.02   # degradation of free molecules (1/h)
def lam(a,kd): return (-(a+kd)+np.sqrt((a+kd)**2+4*a*kd))/2
def kd_for_ratio(f,a=A,kd1=KD1):
    return brentq(lambda k: lam(a,k)/lam(a,kd1)-f,kd1,1e6)
def step(S,rng,dt,kd,err,w):
    Ef,Pf,CE,CP=S
    tot=Ef+w*Pf
    ppt=np.where(tot>0,w*Pf/np.maximum(tot,1e-12),0.0)
    nc=np.minimum(rng.poisson(A*Ef*dt),Ef)
    nP=rng.binomial(nc,ppt); nE=nc-nP
    m=rng.binomial(nE,err)
    Ef=Ef-nc; CE=CE+nE-m; CP=CP+nP+m
    rE=np.minimum(rng.poisson(kd*CE*dt),CE); rP=np.minimum(rng.poisson(kd*CP*dt),CP)
    CE=CE-rE; CP=CP-rP
    Ef=Ef+2*rE+rP; Pf=Pf+rP
    pd=1-np.exp(-DEG*dt)
    Ef=Ef-rng.binomial(Ef,pd); Pf=Pf-rng.binomial(Pf,pd)
    return [Ef,Pf,CE,CP]
def nmol(S): return S[0]+S[1]+2*(S[2]+S[3])
def frac_func(S):
    E=S[0]+S[2]+S[3]+S[2]; P=S[1]+S[3]
    return E.sum()/max(E.sum()+P.sum(),1)
def run_compart(M,Ndiv,kd_types,type_frac,err,w,T,dt,rng,rec_every=1.0):
    """kd_types: array of kd per cell (cell-type specific). Returns times, frac_func, frac_type2."""
    N0=Ndiv//2
    S=[np.full(M,N0,dtype=np.int64),np.zeros(M,np.int64),np.zeros(M,np.int64),np.zeros(M,np.int64)]
    kd=kd_types.copy().astype(float)
    ts=[];ff=[];ft=[];t=0.0;nxt=0.0
    nsteps=int(round(T/dt))
    for k in range(nsteps+1):
        if t>=nxt-1e-9:
            ts.append(t); ff.append(frac_func(S)); ft.append(np.mean(kd>KD1+1e-9)); nxt+=rec_every
        S=step(S,rng,dt,kd,err,w)
        N=nmol(S)
        for i in np.where(N>=Ndiv)[0]:
            d=[rng.binomial(S[c][i],0.5) for c in range(4)]
            j=rng.integers(M-1); j=j+1 if j>=i else j
            for c in range(4): S[c][i]-=d[c]; S[c][j]=d[c]
            kd[j]=kd[i]
        t+=dt
    return np.array(ts),np.array(ff),np.array(ft)
def run_mixed(M,Ndiv,kd,err,w,T,dt,rng,rec_every=1.0):
    """No compartments: one pool of M*Ndiv/2 molecules, binomially diluted 1/2 when it reaches M*Ndiv."""
    n0=M*(Ndiv//2)
    S=[np.array([n0],dtype=np.int64),np.zeros(1,np.int64),np.zeros(1,np.int64),np.zeros(1,np.int64)]
    ts=[];ff=[];t=0.0;nxt=0.0
    for k in range(int(round(T/dt))+1):
        if t>=nxt-1e-9: ts.append(t); ff.append(frac_func(S)); nxt+=rec_every
        S=step(S,rng,dt,np.array([kd]),err,w)
        if nmol(S)[0]>=M*Ndiv: S=[rng.binomial(S[c],0.5) for c in range(4)]
        t+=dt
    return np.array(ts),np.array(ff)
