"""
Predicts the self-consistent local pH of each coacervate phase (inner/template
and outer/catalytic) from its titratable composition, coupling acid-base
equilibria to Donnan partitioning of protons against a bulk pool.

Physics
-------
A phase holds fixed (polyelectrolyte) titratable groups plus mobile pool ions.
Mobile monovalent cations partition into the phase by a Donnan factor lambda,
monovalent anions by 1/lambda, divalent cations by lambda**2. Protons follow
the monovalent cations, so the LOCAL pH is:

        pH_in = pH_bulk - log10(lambda)

The net FIXED charge of the phase depends on that local pH (carboxyls
deprotonate, amines protonate). Electroneutrality inside the phase,

        lambda*Cc + 2*lambda**2*Cd - Ca/lambda + Zfix(pH_in) = 0,

is solved for lambda (Brent). A net-anionic phase needs lambda>1 (cations,
incl. H+, enriched) -> pH drops; but as pH falls toward the carboxyl pKa the
carboxyls reprotonate and REDUCE the net charge -> the drop is buffer-limited
and self-consistent. That coupling is the whole point.
"""
import numpy as np
from scipy.optimize import brentq

# ---- bulk pool (freshwater geothermal; mol/L) -------------------------------
POOL_pH = 6.0
Cc = 0.020    # background monovalent cations (Na+K+NH4), mol/L
Ca = 0.020    # background monovalent anions (Cl+HCO3), mol/L
Cd = 0.0043   # background divalent cations (Mg+Ca+Fe), mol/L


def z_fixed(pH, comp):
    """Net fixed charge density (mol/L, signed) at a given local pH."""
    # carboxyls (Asp/Glu): -1 when deprotonated
    q_cooh = -comp["cooh"] / (1.0 + 10.0 ** (comp["pka_cooh"] - pH))
    # amines (putrescine/lysine): +1 when protonated
    q_amine = +comp["amine"] / (1.0 + 10.0 ** (pH - comp["pka_amine"]))
    # RNA phosphates: fully ionized (pKa ~1) -> -1 each above pH ~2
    q_phos = -comp["phos"]
    # guanidinium (if any): pKa ~12.5, essentially always +1
    q_guan = +comp.get("guan", 0.0) / (1.0 + 10.0 ** (pH - 12.5))
    return q_cooh + q_amine + q_phos + q_guan


def electroneutrality(lam, comp):
    pH_in = POOL_pH - np.log10(lam)
    mobile = lam * Cc + 2.0 * lam ** 2 * Cd - Ca / lam
    return mobile + z_fixed(pH_in, comp)


def solve_phase(comp):
    lam = brentq(electroneutrality, 1e-4, 1e4, args=(comp,), xtol=1e-12)
    pH_in = POOL_pH - np.log10(lam)
    return {
        "pH": pH_in,
        "lambda": lam,
        "mono_cation_M": lam * Cc,
        "divalent_mM": lam ** 2 * Cd * 1000.0,
        "z_fixed_mM": z_fixed(pH_in, comp) * 1000.0,
    }


# ---- compositions (mol/L of groups in the dense phase) ----------------------
# Real coacervates are nearly charge-balanced (polyanion ~ polycation); only a
# modest NET charge remains, and that net charge is the single Donnan knob that
# sets BOTH the local pH offset and the divalent enrichment.
OUTER_NEUTRAL = dict(phos=0.20, cooh=0.02, pka_cooh=4.5, amine=0.205,
                     pka_amine=10.0, guan=0.0)          # net ~ 0
INNER_CATIONIC = dict(phos=0.20, cooh=0.00, pka_cooh=4.5, amine=0.36,
                      pka_amine=10.0, guan=0.0)         # polycation excess, net +
INNER_ANIONIC = dict(phos=0.20, cooh=0.30, pka_cooh=4.5, amine=0.08,
                     pka_amine=10.0, guan=0.0)          # acidic/RNA excess, net -


def report(label, comp):
    r = solve_phase(comp)
    print(f"{label:26s} pH {r['pH']:5.2f}   "
          f"divalent {r['divalent_mM']:7.1f} mM   "
          f"net fixed {r['z_fixed_mM']:+6.0f} mM   lambda {r['lambda']:5.2f}")
    return r


if __name__ == "__main__":
    print(f"Bulk pool: pH {POOL_pH}, monovalent {Cc*1000:.0f} mM, "
          f"divalent {Cd*1000:.1f} mM  (outer-phase reference)\n")
    ro = report("Outer (near-neutral)", OUTER_NEUTRAL)
    print()
    print("Two ways to make the inner phase destabilizing:\n")
    rc = report("Inner A: polycation-excess", INNER_CATIONIC)
    ra = report("Inner B: acidic/anion-excess", INNER_ANIONIC)
    print()
    print(f"  A (alkaline, Mg-poor):  pH gap {rc['pH']-ro['pH']:+.2f} vs outer, "
          f"divalent {rc['divalent_mM']:.1f} vs {ro['divalent_mM']:.1f} mM")
    print(f"  B (acidic,  Mg-rich) :  pH gap {ra['pH']-ro['pH']:+.2f} vs outer, "
          f"divalent {ra['divalent_mM']:.1f} vs {ro['divalent_mM']:.1f} mM")

    # ---- the core coupling: net charge sets pH AND divalents together --------
    print("\nThe unavoidable coupling (scan of net fixed charge on the inner "
          "phase):")
    print("  net fixed[mM]   local pH   divalent[mM]   duplex effect")
    for amine in (0.40, 0.33, 0.26, 0.205, 0.14, 0.05):
        comp = dict(phos=0.20, cooh=0.02, pka_cooh=4.5, amine=amine,
                    pka_amine=10.0, guan=0.0)
        r = solve_phase(comp)
        if r["pH"] > ro["pH"] + 0.3:
            eff = "alkaline + Mg-poor -> destabilizes"
        elif r["pH"] < ro["pH"] - 0.3:
            eff = "acidic + Mg-rich -> mixed"
        else:
            eff = "~neutral, ~uniform Mg"
        print(f"   {r['z_fixed_mM']:+7.0f}       {r['pH']:5.2f}      "
              f"{r['divalent_mM']:7.1f}     {eff}")
