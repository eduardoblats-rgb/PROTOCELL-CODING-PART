import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=json.load(open("fig56.json")); BL,OR,GR="#2a78d6","#ee6a33","#1aae7d"
f5=R["fig5"]
fig,ax=plt.subplots(1,2,figsize=(14,4.8))
a=ax[0]
a.plot(f5["4x"]["t"],100*np.array(f5["4x"]["f"]),color=GR,lw=2.6); a.plot(f5["1.3x"]["t"],100*np.array(f5["1.3x"]["f"]),color=BL,lw=2.6)
a.text(9,95,"Two-phase 4× faster",color=GR,fontweight="bold"); a.text(14,62,"Two-phase only 1.3× faster",color=BL,fontweight="bold")
a.set_xlabel("Time (h)"); a.set_ylabel("Two-phase protocells (%)"); a.set_ylim(0,105); a.set_xlim(0,40)
a.set_title("Competition: populations start 50% of each type",loc="left",fontsize=11)
b=ax[1]; p=f5["par"]
b.plot(p["t"],100*np.array(p["comp"]),color=GR,lw=2.6); b.plot(p["t"],100*np.array(p["mixed"]),color=OR,lw=2.6)
b.text(55,70,"With compartments",color=GR,fontweight="bold"); b.text(50,12,"Without compartments (fully mixed)",color=OR,fontweight="bold")
b.set_xlabel("Time (h)"); b.set_ylabel("Functional RNA (% of total)"); b.set_ylim(0,105); b.set_xlim(0,120)
b.set_title("Parasites: 2% of copies are parasites, copied 1.5× better",loc="left",fontsize=11)
for x in ax: x.spines[['top','right']].set_visible(False)
fig.suptitle("Two-phase protocells displace one-phase protocells; compartments restrain parasites",x=0.01,ha="left",fontsize=13)
fig.text(0.01,0.005,"Stochastic-corrector model: 200 protocells, division at 20 molecules, product-inhibited replicases (catalysis 3 h⁻¹, release 0.5 h⁻¹ in one phase). Mean of 3 simulations.",fontsize=8,color="#555")
fig.tight_layout(rect=(0,.03,1,.94)); fig.savefig("fig5.png",dpi=200)
# Fig 6
H=100*np.array(R["heat"]); sizes=R["sizes"]; errs=R["errs"]
fig,ax=plt.subplots(1,2,figsize=(14,5),gridspec_kw={"width_ratios":[1.35,1]})
a=ax[0]; im=a.imshow(H,cmap="Blues",vmin=0,vmax=100,aspect="auto",origin="lower")
a.set_xticks(range(len(errs))); a.set_xticklabels([str(e) for e in errs]); a.set_yticks(range(len(sizes))); a.set_yticklabels(sizes)
for i in range(len(sizes)):
    for j in range(len(errs)): a.text(j,i,f"{H[i,j]:.0f}",ha="center",va="center",fontsize=9,color="white" if H[i,j]>55 else "#333")
a.set_xlabel("Copy error (% of copies that become parasites)"); a.set_ylabel("Protocell size (molecules at division)")
a.set_title("With compartments: smaller protocells tolerate more error",loc="left",fontsize=11)
cb=fig.colorbar(im,ax=a); cb.set_label("Functional RNA (%)")
b=ax[1]; i40=sizes.index(40)
b.plot(errs,H[i40],"-o",color=BL,lw=2.4); b.plot(errs,100*np.array(R["mixed40"]),"-o",color=OR,lw=2.4)
b.text(12,52,"With compartments",color=BL,fontweight="bold"); b.text(8,8,"Without compartments",color=OR,fontweight="bold")
b.set_xlabel("Copy error (%)"); b.set_ylabel("Functional RNA (%)"); b.set_ylim(0,100); b.spines[['top','right']].set_visible(False)
b.set_title("Protocell of 40 molecules",loc="left",fontsize=11)
fig.suptitle("Error threshold: there is a maximum tolerated error, and compartments raise it",x=0.01,ha="left",fontsize=13)
fig.text(0.01,0.005,"Stochastic-corrector model; parasite copy advantage ×2. Mean functional RNA over 150–200 h, mean of 3 simulations.",fontsize=8,color="#555")
fig.tight_layout(rect=(0,.03,1,.94)); fig.savefig("fig6.png",dpi=200)
