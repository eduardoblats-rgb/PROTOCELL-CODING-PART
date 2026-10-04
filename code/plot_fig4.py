import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
import r3c_kin as K
BL, INK = "#2a78d6", "#1d2330"; GS = ["#8fd3b5", "#2fae82", "#0b6e4d"]
kobs, kon, c, phi, mg, d = 0.62, 1e6, 5.0, 0.2, 0.1, 0.27
J = K.calibrate(kobs, kon, c); Ts = np.arange(30, 47, 0.5)
plt.rcParams["font.family"] = "DejaVu Sans"
fig, ax = plt.subplots(1, 2, figsize=(14.5, 6.4), gridspec_kw={"width_ratios": [1.2, 1], "wspace": .28})
# ---- A
a = ax[0]
a.semilogy(Ts, [K.lam_one(T, J, kobs, kon, c) for T in Ts], color=BL, lw=3.4)
labs = ["Two phases, slow exchange (10 /h)", "Two phases, medium exchange (100 /h)", "Two phases, fast exchange (3,600 /h)"]
for kd, col in zip((10, 100, 3600), GS):
    a.semilogy(Ts, [K.lam_two(T, J, kobs, kon, c, phi, mg, d, kd) for T in Ts], color=col, lw=3.2)
a.text(30.3, 36, labs[2], color=GS[2], fontsize=11.5, fontweight="bold", va="bottom")
a.text(30.3, 3.9, labs[1], color=GS[1], fontsize=11.5, fontweight="bold", va="bottom")
a.text(30.3, 0.72, labs[0], color="#3d9172", fontsize=11.5, fontweight="bold", va="top")
a.text(30.4, 3e-3, "One phase", color=BL, fontsize=13, fontweight="bold", ha="left", va="center", rotation=0)
a.axvline(42, color="#999", ls=":", lw=1.4)
for kd, col, txt in ((3600, GS[2], "28×"), (100, GS[1], "10×"), (10, GS[0], "2.4×")):
    y = K.lam_two(42, J, kobs, kon, c, phi, mg, d, kd)
    a.plot([42], [y], "o", color=col, ms=9, mec="white", mew=1.5, zorder=5)
    a.text(42.35, y*1.0, txt, color=INK, fontsize=12, fontweight="bold", va="center", bbox=dict(fc="white", ec="none", pad=1.5), zorder=6)
a.plot([42], [0.98], "o", color=BL, ms=11, mec="white", mew=1.5, zorder=5)
a.annotate("Measured: ~1 per hour\n(doubling in ~1 h, 42 °C)", xy=(42, 0.98), xytext=(36.2, 0.05), fontsize=11.5, color=INK,
           arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.3), ha="center")
a.set_ylim(1e-3, 150); a.set_xlim(30, 46.3)
a.yaxis.set_major_locator(FixedLocator([1e-3, 1e-2, 1e-1, 1, 10, 100]))
a.yaxis.set_major_formatter(FixedFormatter(["0.001", "0.01", "0.1", "1", "10", "100"])); a.yaxis.set_minor_locator(NullLocator())
a.set_xlabel("Temperature (°C)", fontsize=13); a.set_ylabel("Replication growth rate (per hour, log scale)", fontsize=13)
a.tick_params(labelsize=11.5); a.spines[["top", "right"]].set_visible(False)
a.set_title("A   Two phases keep replicating at lower temperature", loc="left", fontsize=14, fontweight="bold", color=INK, pad=12)
# ---- B
R = json.load(open("r3c_scan.json")); groups = [[r["speed"] for r in R if r["p"][6] == k] for k in (10, 100, 3600)]
b = ax[1]; bp = b.boxplot(groups, positions=[1, 2, 3], widths=.55, patch_artist=True, showfliers=False,
                          medianprops=dict(color="k", lw=2), whiskerprops=dict(lw=1.4), capprops=dict(lw=1.4), boxprops=dict(lw=1.4))
for p, col in zip(bp["boxes"], GS): p.set_facecolor(col); p.set_alpha(.9)
rng = np.random.default_rng(0)
for i, g in enumerate(groups, 1):
    b.scatter(i + rng.uniform(-.22, .22, len(g)), g, s=6, color="#333", alpha=.18, zorder=3)
    b.text(i + .36, np.median(g), f"{np.median(g):.1f}×" if np.median(g) < 10 else f"{np.median(g):.0f}×", fontsize=13, fontweight="bold", va="center", color=INK)
b.set_yscale("log"); b.set_ylim(.8, 90)
b.yaxis.set_major_locator(FixedLocator([1, 3, 10, 30, 60])); b.yaxis.set_major_formatter(FixedFormatter(["1×", "3×", "10×", "30×", "60×"])); b.yaxis.set_minor_locator(NullLocator())
b.axhline(1, color="#888", ls=":", lw=1.4); b.text(0.55, 1.07, "no gain", fontsize=10.5, color="#666")
b.set_xticks([1, 2, 3]); b.set_xticklabels(["Slow\n(10 /h)", "Medium\n(100 /h)", "Fast\n(3,600 /h)"], fontsize=12)
b.set_xlim(.5, 3.75); b.tick_params(axis="y", labelsize=12)
b.set_xlabel("Speed of exchange between the two phases", fontsize=13); b.set_ylabel("Speed-up over one phase (at 42 °C)", fontsize=13)
b.spines[["top", "right"]].set_visible(False)
b.set_title("B   The gain depends mostly on exchange speed", loc="left", fontsize=14, fontweight="bold", color=INK, pad=12)
fig.text(0.01, 0.945, "Two phases accelerate R3C cross-replication by relieving product inhibition", fontsize=16, fontweight="bold", color=INK)
fig.text(0.01, 0.012, "A: reference parameters (substrates 5 µM, 20 % of time in the inner phase, 0.1 mM inner Mg²⁺, δ = 0.27 kcal/mol per bp). Dot at 42 °C = calibration to Lincoln & Joyce (2009), 25 mM MgCl₂;\n"
         "multipliers = gain over one phase at 42 °C.  B: 486 parameter combinations per column; box = quartiles, line = median, dots = individual combinations.\n"
         "The model has no thermal inactivation of the enzyme, so it is not valid far above 42 °C.", fontsize=9, color="#555", va="bottom")
fig.subplots_adjust(left=.07, right=.98, top=.86, bottom=.22)
fig.savefig("fig4.png", dpi=200, facecolor="white")
