"""Make all report figures for the DOX/membrane project.
Run from ~/Desktop/dox_md:   python report/make_figures.py
Figures are written to report/figures/ (PNG, 300 dpi)."""
import os, sys, warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Desktop/dox_md")
A, B = os.path.join(ROOT, "system_A"), os.path.join(ROOT, "system_B")
OUT = os.path.join(ROOT, "report", "figures"); os.makedirs(OUT, exist_ok=True)

CA, CB = "#2a78d6", "#eb6834"
C_WATER, C_P, C_TAIL, C_CHOL, C_DOX = "#1baf7a", "#eda100", "#6250d6", "#e87ba4", "#222222"
INK, INK2, GRID = "#222222", "#6b6b6b", "#e6e6e3"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "lines.linewidth": 2,
                     "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight"})
made, missing = [], []
def need(*paths):
    miss = [p for p in paths if not os.path.exists(p)]
    missing.extend(miss); return not miss
def save(fig, name):
    fig.savefig(os.path.join(OUT, name)); plt.close(fig); made.append(name)
def csv(p): return np.genfromtxt(p, delimiter=",", names=True)
def smooth(y, n):
    if n <= 1: return y
    k = np.ones(n); return np.convolve(y, k, mode="same") / np.convolve(np.ones_like(y), k, mode="same")
def bound_states(t, dmin, cut=5.0, minlen=1.0):
    dt = np.median(np.diff(t)); s = dmin < cut
    def runs(s):
        out, st = [], 0
        for i in range(1, len(s) + 1):
            if i == len(s) or s[i] != s[st]: out.append([bool(s[st]), st, i]); st = i
        return out
    while True:
        r = runs(s); short = [k for k, (v, a, b) in enumerate(r) if (b - a) * dt < minlen]
        if not short or len(r) == 1: break
        k = min(short, key=lambda k: r[k][2] - r[k][1])
        nb = [r[j] for j in (k - 1, k + 1) if 0 <= j < len(r)]
        s[r[k][1]:r[k][2]] = max(nb, key=lambda x: x[2] - x[1])[0]
    return runs(s)

trA, trB = os.path.join(A, "dox_track_A_full.csv"), os.path.join(B, "dox_track_B.csv")

if need(trA, trB):
    fig, axs = plt.subplots(2, 1, figsize=(7.5, 5.2), sharex=True)
    for ax, f, c, lab in ((axs[0], trA, CA, "A: fluid DPPC (323 K)"), (axs[1], trB, CB, "B: DPPC + 50% cholesterol (303 K)")):
        d = csv(f); t, h, pp = d["t_ns"], np.abs(d["dz_A"]), d["PP_A"] / 2
        ax.axhspan(0, pp.mean(), color="#f1f0ec", zorder=0, lw=0)
        ax.axhline(pp.mean(), color=INK2, lw=1, ls="--")
        ax.text(149, pp.mean() - 1.5, f"phosphate plane ({pp.mean():.1f} Å)", color=INK2, fontsize=8, va="top", ha="right",
                bbox=dict(fc="#f1f0ec", ec="none", pad=1.5))
        ax.text(1, 2, "inside the membrane", color=INK2, fontsize=8)
        ax.plot(t, h, color=c, lw=0.4, alpha=0.35)
        ax.plot(t, smooth(h, 50), color=c, lw=1.8)
        ax.set_ylim(0, 50); ax.set_ylabel("|DOX⁺ height| from\nbilayer centre (Å)")
        ax.set_title(lab, loc="left", fontsize=10, color=INK)
    axs[1].set_xlabel("time (ns)"); axs[1].set_xlim(0, 150)
    fig.text(0.99, -0.03, "thin line: every 10 ps; thick line: 0.5 ns running mean", ha="right", fontsize=7, color=INK2)
    save(fig, "fig1_depth_vs_time.png")

if need(trA, trB):
    fig, ax = plt.subplots(figsize=(7.5, 1.9))
    for y, f, c, lab in ((1, trA, CA, "A"), (0, trB, CB, "B")):
        d = csv(f); t = d["t_ns"]; dt = np.median(np.diff(t))
        for st, a, b in bound_states(t, d["mindist_P_A"]):
            if st: ax.broken_barh([(t[a], (b - a) * dt)], (y - 0.32, 0.64), color=c, lw=0)
        ax.broken_barh([(t[0], t[-1] - t[0])], (y - 0.32, 0.64), facecolor="none", edgecolor=GRID, lw=0.8)
    ax.set_yticks([1, 0]); ax.set_yticklabels(["A: fluid DPPC", "B: DPPC/chol"]); ax.grid(False)
    ax.set_xlim(0, 150); ax.set_ylim(-0.6, 1.6); ax.set_xlabel("time (ns)")
    ax.set_title("Filled = DOX⁺ in contact with the membrane (within 5 Å of a phosphorus; blips < 1 ns ignored)",
                 loc="left", fontsize=9, color=INK)
    save(fig, "fig2_binding_timeline.png")

dA, dB = os.path.join(A, "dens_A_ext.csv"), os.path.join(B, "dens_B_late.csv")
if need(dA, dB):
    fig, axs = plt.subplots(1, 2, figsize=(8.5, 3.4), sharey=True)
    for ax, f, lab in ((axs[0], dA, "A: fluid DPPC, 100–150 ns"), (axs[1], dB, "B: DPPC/chol, 118–136 ns")):
        d = csv(f); z = d["z_A"]; m = z >= 0
        def fold(col): return (d[col] + d[col][::-1])[m]
        groups = [("water_O", "water", C_WATER), ("phosphate_P", "phosphate", C_P), ("lipid_tails", "lipid tails", C_TAIL)]
        if "cholesterol_OH_O3" in d.dtype.names: groups.append(("cholesterol_OH_O3", "chol. OH", C_CHOL))
        groups.append(("DOX", "DOX⁺", C_DOX))
        for j, (col, name, c) in enumerate(groups):
            y = fold(col); y = y / y.max()
            ax.plot(z[m], y, color=c, lw=2.4 if name == "DOX⁺" else 1.8, label=name)
            if name == "water": continue
            i = np.argmax(y); ax.annotate(name, (z[m][i], 1.0), xytext=(0, {"phosphate": 34, "DOX⁺": 19}.get(name, 5)),
                                          textcoords="offset points", ha="center", fontsize=8, color=INK,
                                          arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))
        ax.set_xlim(0, 50); ax.set_ylim(0, 1.3); ax.set_xlabel("distance from bilayer centre (Å)")
        ax.set_title(lab, loc="left", fontsize=10, color=INK)
    axs[0].set_ylabel("density (each normalised to its peak)")
    axs[1].legend(loc="center right", fontsize=8)
    save(fig, "fig3_density_profiles.png")

rA, rB1, rB2 = (os.path.join(A, "rdf_A_ext.csv"), os.path.join(B, "rdf_B.csv"), os.path.join(B, "rdf_B_late.csv"))
if need(rA, rB1, rB2):
    fig, axs = plt.subplots(1, 3, figsize=(9.5, 3.0), sharex=True)
    sets = [(rA, CA, "-", "A 100–150 ns"), (rB1, CB, "-", "B 30–49 ns"), (rB2, CB, ":", "B 118–136 ns")]
    for ax, col, title in zip(axs, ("phosphate_O", "ester_O", "water_O"), ("phosphate O", "ester O", "water O")):
        for f, c, ls, lab in sets:
            d = csv(f); ax.plot(d["r_A"], d[col], color=c, ls=ls, label=lab)
        ax.set_title(f"NH₃⁺ N – {title}", loc="left", fontsize=10, color=INK); ax.set_xlim(1.5, 10)
        ax.set_xlabel("r (Å)")
    axs[0].set_ylabel("g(r)"); axs[0].axvline(2.65, color=INK2, lw=1, ls="--")
    axs[0].text(3.4, axs[0].get_ylim()[1] * 0.9, "← paper: 2.65 Å", fontsize=8, color=INK2)
    axs[2].legend(fontsize=8)
    save(fig, "fig4_rdf.png")

sA, sB = os.path.join(A, "scd_A.csv"), os.path.join(B, "scd_B.csv")
if need(sA, sB):
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    for f, c, lab in ((sA, CA, "A: fluid DPPC"), (sB, CB, "B: DPPC + 50% cholesterol")):
        d = csv(f); n = d["carbon"]
        ax.plot(n, d["chain2"], color=c, marker="o", ms=4, label=f"{lab}, sn-2")
        ax.plot(n, d["chain3"], color=c, marker="s", ms=4, ls="--", mfc="white", label=f"{lab}, sn-1")
    ax.set_xlabel("carbon number along the tail"); ax.set_ylabel("order parameter  −S$_{CD}$")
    ax.set_ylim(0, 0.5); ax.set_xticks(range(2, 17, 2)); ax.legend(fontsize=8)
    save(fig, "fig5_order_parameter.png")

eA, eB = os.path.join(A, "energy_A.xvg"), os.path.join(B, "energy_B.xvg")
if need(eA, eB):
    a = np.loadtxt(eA, comments=("#", "@")); b = np.loadtxt(eB, comments=("#", "@"))
    k = 4.184; n = 100
    fig, axs = plt.subplots(2, 1, figsize=(7.5, 5.0), sharex=True)
    series = [(a[:, 0] / 1000, a[:, 1] / k, a[:, 2] / k, CA, "A: DOX⁺–DPPC"),
              (b[:, 0] / 1000, (b[:, 1] + b[:, 3]) / k, (b[:, 2] + b[:, 4]) / k, CB, "B: DOX⁺–(DPPC + cholesterol)")]
    for t, el, lj, c, lab in series:
        axs[0].plot(t, smooth(el, n), color=c, label=lab); axs[1].plot(t, smooth(lj, n), color=c, label=lab)
    axs[0].set_ylabel("electrostatic\n(kcal/mol)"); axs[1].set_ylabel("van der Waals\n(kcal/mol)")
    axs[1].axhline(-39, color=CA, lw=1, ls="--"); axs[1].text(2, -36, "paper, fluid DPPC: −39", fontsize=8, color=INK2)
    axs[0].legend(fontsize=8, loc="lower left"); axs[1].set_xlabel("time (ns)"); axs[1].set_xlim(0, 150)
    fig.text(0.99, -0.03, "1 ns running means; short-range (≤ 1.2 nm) part only", ha="right", fontsize=7, color=INK2)
    save(fig, "fig6_interaction_energy.png")

hb = {
    "A 60–100 ns":   (0.75 + 0.36, 0.15 + 0.10, 1.23 + 0.98),
    "A 100–150 ns":  (0.86 + 0.27, 0.43 + 0.33, 0.82 + 0.71),
    "B 30–49 ns":    (0.99 + 0.21, 0.01 + 0.00, 1.25 + 1.29),
    "B 118–136 ns":  (0.17 + 0.16, 0.00 + 0.01, 1.87 + 1.27)}
fig, ax = plt.subplots(figsize=(7, 3.2)); x = np.arange(3); w = 0.2
styles = [(CA, None), (CA, "////"), (CB, None), (CB, "////")]
for i, ((lab, v), (c, h)) in enumerate(zip(hb.items(), styles)):
    ax.bar(x + (i - 1.5) * w, v, w * 0.9, color=c if h is None else "white", edgecolor=c, hatch=h, lw=1.2, label=lab)
    for xb, vb in zip(x + (i - 1.5) * w, v):
        ax.text(xb, vb + 0.04, f"{vb:.2f}", ha="center", fontsize=7, color=INK)
ax.set_xticks(x); ax.set_xticklabels(["phosphate O", "ester O", "water"]); ax.grid(axis="x", visible=False)
ax.set_ylabel("H-bonds donated by DOX⁺\nper frame"); ax.legend(fontsize=8, ncol=2)
save(fig, "fig7_hbonds.png")

D = [("A in water (16.7–30.5 ns)", 50.21, CA), ("B in water (48.9–68.9 ns)", 34.69, CB), ("B in water (75–91 ns)", 95.62, CB),
     ("B bound (32–49 ns)", 4.89, CB), ("B bound (118–136 ns)", 5.98, CB),
     ("A inserted (60–90 ns)", 2.13, CA), ("A inserted (90–120 ns)", 1.74, CA), ("A inserted (120–150 ns)", 0.81, CA)]
fig, ax = plt.subplots(figsize=(6.5, 3.4)); y = np.arange(len(D))[::-1]
for yi, (lab, v, c) in zip(y, D):
    ax.errorbar(v, yi, xerr=[[0.4 * v], [0.4 * v]], fmt="o", ms=8, color=c, ecolor=c, elinewidth=1, capsize=3)
    ax.text(v * 1.55, yi, f"{v:.1f}", va="center", fontsize=8, color=INK)
ax.set_yticks(y); ax.set_yticklabels([d[0] for d in D], fontsize=8); ax.set_xscale("log"); ax.set_xlim(0.3, 300)
ax.set_xlabel("lateral diffusion coefficient D (10⁻⁷ cm²/s, log scale)")
fig.text(0.99, -0.10, "bars: ±40% (one trajectory piece per value)", ha="right", fontsize=7, color=INK2)
save(fig, "fig8_diffusion.png")

print("saved to", OUT)
for m in made: print("  ", m)
if missing:
    print("\nskipped figures because these files were not found:")
    for m in sorted(set(missing)): print("  ", m)
