"""Two-phase kinetics of R3C cross-replication (see r3c.py for helices).
Growth rate of a cross-replicating pair with catalysis a and product release kd:
  lam = (-(a+kd) + sqrt((a+kd)^2 + 4 a kd))/2
Outer phase: Mg-rich, catalysis + release. Inner phase: enzyme inactive, release only.
Time-weighted two-phase rates: a_eff=(1-phi)*a_out ; kd_eff=(1-phi)*kd_out + phi*kd_in.
Single free parameter J (kcal/mol) is calibrated so that one-phase growth = 0.98 /h at 42 C, 25 mM MgCl2 (Lincoln & Joyce 2009).
"""
import numpy as np
from functools import lru_cache
from scipy.optimize import brentq
from r3c import helix_dG, NBP, A_SET, B_SET, ALL
RT=lambda Tc: 0.0019872*(Tc+273.15)
KD_MAX=3600.0   # /h, ceiling on release rate (1 per second)
SALT_OUT=1.0    # 25 mM Mg2+ >> 3 mM ~ 1 M Na+ -> capped at 1.0 M
def salt_in(mg_mM): return 0.015+333.0*mg_mM*1e-3
@lru_cache(maxsize=None)
def hdG(h,Tc,salt): return helix_dG(h,float(Tc),float(salt))
def S(names,Tc,salt,delta=0.0): return sum(hdG(h,Tc,salt)+delta*NBP[h] for h in names)
def lam(a,kd): return (-(a+kd)+np.sqrt((a+kd)**2+4*a*kd))/2
THETA_REF=0.9
def _Jsub(names):
    # offset so that theta=THETA_REF for a 5 uM substrate at the reference condition
    rt=RT(42); Kd=5e-6*(1-THETA_REF)/THETA_REF
    return rt*np.log(Kd)-S(names,42,SALT_OUT)
_JA=None
def rates(Tc,salt,delta,J,kcat_min,kon,conc_uM,active=True,kdmax=KD_MAX):
    rt=RT(Tc); c=conc_uM*1e-6
    dGp=S(ALL,Tc,salt,delta)+J
    kd=min(kdmax, kon*3600*np.exp(dGp/rt))          # kd = kon*Kd ; Kd=exp(dG/RT) M
    dGa=S(A_SET,Tc,salt,delta)+_Jsub(A_SET); dGb=S(B_SET,Tc,salt,delta)+_Jsub(B_SET)
    th=lambda dG: c/(c+np.exp(dG/rt))
    # kcat_min = measured observed ligation rate (/min) at 42 C, 5 uM substrates, 25 mM Mg (Lincoln & Joyce 2009);
    # it already includes substrate occupancy, so scale by occupancy relative to that reference condition.
    k42=kcat_min*60.0; EA=15.0
    kc=k42*np.exp(-EA/0.0019872*(1/(Tc+273.15)-1/315.15))
    a=kc*th(dGa)*th(dGb)/THETA_REF**2 if active else 0.0
    return a,kd
def lam_one(Tc,J,kcat,kon,conc,delta=0.0):
    a,kd=rates(Tc,SALT_OUT,0.0,J,kcat,kon,conc); return lam(a,kd)
def lam_two(Tc,J,kcat,kon,conc,phi,mg_in,delta,kdin_max=KD_MAX):
    ao,ko=rates(Tc,SALT_OUT,0.0,J,kcat,kon,conc)
    _,ki=rates(Tc,salt_in(mg_in),delta,J,kcat,kon,conc,active=False,kdmax=kdin_max)
    return lam((1-phi)*ao,(1-phi)*ko+phi*ki)
def calibrate(kcat,kon,conc,target=0.98):
    return brentq(lambda J: lam_one(42,J,kcat,kon,conc)-target,-25,15)
