import json, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from melt_lib import tm
R = json.load(open("fig23.json")); Ts = list(range(0, 101, 4)); names = list(R)
BL, OR, WIN, INK = "#2a78d6", "#e8601c", "#f2c14e", "#1d2330"
plt.rcParams["font.family"] = "DejaVu Sans"
fig, axs = plt.subplots(2, 3, figsize=(13.5, 8.2), sharey=True)
for a, n in zip(axs.flat, names):
    o, i = R[n]["out_low"], R[n]["in_low_0.27"]
    t0, t1 = tm(Ts, o), tm(Ts, i)
    a.axvspan(t1, t0, color=WIN, alpha=.35, lw=0)                     # switch window
    a.fill_between(Ts, R[n]["in_low_0.15"], R[n]["in_low_0.45"], color=OR, alpha=.15, lw=0)
    a.plot(Ts, o, color=BL, lw=3); a.plot(Ts, i, color=OR, lw=3)
    a.axhline(.5, color="#999", lw=1, ls=":")
    for t, c in ((t0, BL), (t1, OR)):
        a.plot([t], [.5], "o", color=c, ms=9, mec="white", mew=1.5, zorder=5)
    a.text(3, 0.1, f"ΔTm = {t0-t1:.0f} °C", ha="left", va="center", fontsize=14, fontweight="bold", color=INK)
    nm, rest = n.split(" (", 1)
    a.set_title(f"{nm}\n{rest.rstrip(')')}", fontsize=12, color=INK, pad=8)
    a.set_ylim(-0.03, 1.05); a.set_xlim(0, 100)
    a.spines[["top", "right"]].set_visible(False); a.tick_params(labelsize=11)
for a in axs.flat: a.set_xlabel("Temperature (°C)", fontsize=12)
for a in axs[:, 0]: a.set_ylabel("Fraction folded (native)", fontsize=12)
# explanatory panel
ax = axs[1, 2]; ax.axis("off")
h = [Line2D([0], [0], color=BL, lw=3), Line2D([0], [0], color=OR, lw=3),
     Patch(fc=OR, alpha=.25), Patch(fc=WIN, alpha=.5),
     Line2D([0], [0], marker="o", color="w", mfc="#777", ms=9)]
l = ["Outer phase (reference)", "Inner phase (δ = 0.27 kcal/mol per bp)",
     "Inner phase, δ range 0.15–0.45", "Switch window: folded outside,\nunfolded inside",
     "Melting temperature, Tm (50 % folded);\nΔTm = Tm outer − Tm inner"]
ax.legend(h, l, loc="center left", fontsize=12, frameon=False, labelspacing=1.2, handlelength=1.8)
ax.set_title("How to read", fontsize=12, color=INK, loc="left", pad=8)
fig.suptitle("The inner phase lowers ribozyme melting temperature, opening a switch window",
             x=0.01, ha="left", fontsize=15, fontweight="bold", color=INK)
fig.text(0.01, 0.005, "ViennaRNA 2.7.2; 15 mM monovalent salt, no Mg²⁺ (conditions of Choi et al. 2022). "
         "δ: penalty per base pair in the inner phase, calibrated to Choi et al.", fontsize=9.5, color="#555")
fig.tight_layout(rect=(0, .03, 1, .94), h_pad=2); fig.savefig("fig2.png", dpi=200, facecolor="white")
