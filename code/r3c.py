"""R3C cross-replicating ligase (Lincoln & Joyce 2009, Fig. 1B): helix thermodynamics + two-phase kinetics.
Helices read from Fig. 1B (top strand = substrate, 5'->3'; bottom = enzyme E', written 5'->3').
"""
import RNA, numpy as np
R=0.0019872
HELICES={ # name: (substrate strand 5'->3', enzyme strand 5'->3')
 "A1":("GUUCAUGU","GCAUGAAU"),   # 8 bp, G.U at both ends
 "A2":("GGUU","GACC"),           # 4 bp
 "A3":("GAAU","GUUU"),           # 4 bp, G.U wobbles at the junction
 "B1":("GACC","GGUU"),           # 4 bp
 "B2":("GCAACUU","AAGUUGU"),     # 7 bp
}
NBP={k:len(v[0]) for k,v in HELICES.items()}
def helix_dG(name,Tc,salt):
    a,b=HELICES[name]; n=len(a)
    md=RNA.md(); md.temperature=float(Tc); md.salt=float(salt)
    fc=RNA.fold_compound(a+"&"+b,md)
    return fc.eval_structure("("*n+")"*n)
def sumdG(names,Tc,salt,delta=0.0):
    return sum(helix_dG(h,Tc,salt)+delta*NBP[h] for h in names)
A_SET=["A1","A2","A3"]; B_SET=["B1","B2"]; ALL=A_SET+B_SET
if __name__=="__main__":
    for h in HELICES: print(h,NBP[h],round(helix_dG(h,42,1.0),2),round(helix_dG(h,42,0.015),2))
    print("sum all 42C 1.0M",sumdG(ALL,42,1.0),"0.015M",sumdG(ALL,42,0.015))
