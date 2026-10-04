import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=json.load(open("fig23.json")); Ts=list(range(0,101,4))
names=list(R); short=[n.split(" (")[0]+"\n("+n.split("(")[1].split(",")[0].split(")")[0]+", "+n.split(", ")[-1].replace(")","") if False else n for n in names]
BL,OR,GR,YE="#2a78d6","#ee6a33","#1aae7d","#f0a30a"
# Fig 3
conds=[("Equal Mg in both phases · prebiotic polycation (δ 0.10)",BL,"out_low","in_low_0.1"),
       ("Equal Mg in both phases · arginine-like polycation (δ 0.27)",OR,"out_low","in_low_0.27"),
       ("Mg-rich outer phase · prebiotic polycation (δ 0.10)",GR,"out_high","in_mg_0.1"),
       ("Mg-rich outer phase · arginine-like polycation (δ 0.27)",YE,"out_high","in_mg_0.27")]
fig,ax=plt.subplots(figsize=(10,6.2))
vals={}
for k,n in enumerate(names):
    for j,(lab,c,o,i) in enumerate(conds):
        v=100*max(np.array(R[n][o])-np.array(R[n][i])); vals[(n,j)]=v
        y=k*5+(3-j)
        ax.barh(y,v,color=c,height=.8,label=lab if k==0 else None); ax.text(v+.8,y,f"{v:.0f}",va="center",fontsize=8,zorder=6,bbox=dict(fc="white",ec="none",pad=0.8))
ax.set_yticks([k*5+1.5 for k in range(len(names))]); ax.set_yticklabels([n.replace(" (","\n(") for n in names],fontsize=9)
ax.invert_yaxis(); ax.axvline(40,ls="--",color="#444"); ax.text(40.5,-1.3,"≈ useful-switch threshold",fontsize=9,color="#444")
ax.set_xlabel("Maximum difference in native fold between outer and inner phase (percentage points)")
ax.legend(fontsize=8.5,frameon=False,loc="upper center",bbox_to_anchor=(0.45,-0.12),ncol=2); ax.spines[['top','right']].set_visible(False)
ax.set_title("Folding contrast between phases",loc="left",fontsize=13)
fig.text(0.01,-0.02,"ViennaRNA 2.7.2. Mg²⁺ treated as monovalent-equivalent salt (3 mM Mg²⁺ ≈ 1 M NaCl). Mg-rich outer phase ≥ 3 mM ≈ 1.0 M (cap); Mg-poor inner phase 0.1 mM ≈ 48 mM; equal-Mg bars use 15 mM in both phases.",fontsize=8,color="#555")
fig.tight_layout(rect=(0,.04,1,1)); fig.savefig("fig3.png",dpi=200,bbox_inches="tight")
for n in names: print(n,[round(vals[(n,j)]) for j in range(4)])
