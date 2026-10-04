"""Original schematic of the R3C cross-replicating ligase (Fig. 1B). Drawn from scratch."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

TEAL, TEAL_L = "#1b7f72", "#dcecea"
SLATE, SLATE_L = "#4a5563", "#eceff3"
GOLD = "#b8801f"
INK = "#1d2330"
plt.rcParams["font.family"] = "DejaVu Sans"

fig, ax = plt.subplots(figsize=(10, 6.2))
ax.set_xlim(-6, 106); ax.set_ylim(-9.5, 11.5); ax.axis("off")

yE, hE = 0.0, 2.0          # enzyme bar
yS, hS = 6.0, 2.0          # substrate bars
def bar(x0, x1, y, h, fc, ec):
    ax.add_patch(FancyBboxPatch((x0, y), x1-x0, h, boxstyle="round,pad=0,rounding_size=0.9",
                                fc=fc, ec=ec, lw=2.2, zorder=2))

bar(0, 100, yE, hE, TEAL_L, TEAL)       # enzyme E'
bar(3, 43, yS, hS, SLATE_L, SLATE)      # substrate A
bar(57, 97, yS, hS, SLATE_L, SLATE)     # substrate B

# helices: (x_start, n_bp, substrate label side)
helices = [(5, 8), (17, 4), (62, 4), (74, 7)]
for x0, n in helices:
    for i in range(n):
        x = x0 + i + 0.5
        ax.plot([x, x], [yE+hE, yS], color=TEAL, lw=2.4, solid_capstyle="round", zorder=3)
    # bracket + label above substrate bar
    yb = yS + hS + 0.6
    ax.plot([x0+0.2, x0+0.2, x0+n-0.2, x0+n-0.2], [yb-0.35, yb, yb, yb-0.35], color=TEAL, lw=1.6)
    ax.text(x0+n/2, yb+0.5, f"{n} bp", ha="center", va="bottom", fontsize=13, fontweight="bold", color=TEAL)

# third helix of A sits at the junction; its two terminal pairs are G·U wobbles (gold)
x0 = 36
for i in range(4):
    x = x0 + i + 0.5
    ax.plot([x, x], [yE+hE, yS], color=(GOLD if i in (0, 3) else TEAL), lw=2.4,
            ls=("-" if i in (1, 2) else (0, (2, 1.2))), zorder=3, solid_capstyle="round")
yb = yS + hS + 0.6
ax.plot([x0+0.2, x0+0.2, x0+3.8, x0+3.8], [yb-0.35, yb, yb, yb-0.35], color=TEAL, lw=1.6)
ax.text(x0+2, yb+0.5, "4 bp", ha="center", va="bottom", fontsize=13, fontweight="bold", color=TEAL)
ax.text(x0-1.2, 3.9, "G·U at\nboth ends", ha="right", va="center", fontsize=10.5, color=GOLD, fontweight="bold",
        bbox=dict(fc="white", ec="none", pad=1.5), zorder=4)

# substrate / enzyme names
ax.text(23, yS+hS/2, "Substrate A", ha="center", va="center", fontsize=15, fontweight="bold", color=INK, zorder=4)
ax.text(80, yS+hS/2, "Substrate B", ha="center", va="center", fontsize=15, fontweight="bold", color=INK, zorder=4)
ax.text(50, yE+hE/2, "Enzyme E′  (cross-replicating ligase)", ha="center", va="center",
        fontsize=15, fontweight="bold", color=INK, zorder=4)

# 5'/3' ends
for x, y, t, ha in [(1.6, yS+hS/2, "5′", "right"), (98.4, yS+hS/2, "3′", "left")]:
    ax.text(x, y, t, ha=ha, va="center", fontsize=13, color=SLATE)

# ligation arrow and chemistry labels (clear area above the gap)
ax.add_patch(FancyArrowPatch((44.2, yS+hS+0.7), (55.8, yS+hS+0.7), connectionstyle="arc3,rad=-0.6",
             arrowstyle="-|>", mutation_scale=22, lw=2.6, color=GOLD, zorder=5))
ax.text(50, yS+hS+2.9, "ligation", ha="center", va="bottom", fontsize=15, fontweight="bold", color=GOLD)
ax.text(50, 5.0, "A 3′-OH", ha="center", va="center", fontsize=12, color=GOLD, fontweight="bold")
ax.text(50, 3.1, "B 5′-ppp", ha="center", va="center", fontsize=12, color=GOLD, fontweight="bold")

# catalytic core
ax.plot([50, 50], [yE, -2.2], color=TEAL, lw=2.4)
ax.add_patch(FancyBboxPatch((40, -6.2), 20, 4.0, boxstyle="round,pad=0,rounding_size=1.0",
                            fc="white", ec=TEAL, lw=2.4, zorder=2))
ax.text(50, -4.2, "catalytic core", ha="center", va="center", fontsize=11.5, color=TEAL, fontweight="bold")

# legend / take-home line
ax.text(50, -8.4, "Product A+B (= enzyme E) stays bound to E′ through 5 short helices (27 bp); "
        "releasing it is the rate-limiting step",
        ha="center", va="center", fontsize=11.5, color=INK)

fig.savefig("fig1B.png", dpi=220, bbox_inches="tight", facecolor="white")
