import sys, warnings
import numpy as np
import MDAnalysis as mda
warnings.filterwarnings("ignore")

gro, xtc = sys.argv[1], sys.argv[2]
t0, t1 = float(sys.argv[3]), float(sys.argv[4])          # analysis window in ns
stride = int(sys.argv[5]) if len(sys.argv) > 5 else 1    # use every Nth saved frame
tag = sys.argv[6] if len(sys.argv) > 6 else "scd"

u = mda.Universe(gro, xtc.split(","))
lip = u.select_atoms("resname DPPC")
dt = u.trajectory.dt
f0 = int(round(t0 * 1000 / dt)); f1 = min(u.trajectory.n_frames, int(round(t1 * 1000 / dt)))
frames = list(range(f0, f1, stride))
print(f"{len(lip.residues)} lipids; using {len(frames)} frames from {t0} to {t1} ns")

def hnames(chain, n):
    h = [f"H{n}R", f"H{n}S"] if chain == 2 else [f"H{n}X", f"H{n}Y"]
    if n == 16:
        h.append("H16T" if chain == 2 else "H16Z")
    return h

carbons = range(2, 17)
sel = {}
for c in (2, 3):
    for n in carbons:
        C = lip.select_atoms(f"name C{c}{n}")
        Hs = [lip.select_atoms(f"name {h}") for h in hnames(c, n)]
        assert len(C) == len(lip.residues) and all(len(h) == len(C) for h in Hs), f"atom selection problem at chain {c} carbon {n}"
        sel[(c, n)] = (C, Hs)

S = {k: np.zeros(len(frames)) for k in sel}
for i, fr in enumerate(frames):
    ts = u.trajectory[fr]
    box = ts.dimensions[:3]
    for k, (C, Hs) in sel.items():
        cp = C.positions
        vals = []
        for H in Hs:
            v = H.positions - cp
            v -= box * np.round(v / box)                 # minimum image: bonds can straddle the box edge
            cos2 = v[:, 2] ** 2 / (v ** 2).sum(axis=1)
            vals.append(-(3 * cos2 - 1) / 2)             # S_CD as usually quoted (positive for ordered chains)
        S[k][i] = np.mean(vals)

def stat(a):
    blocks = np.array_split(a, 5) if len(a) >= 5 else [a]
    bm = np.array([b.mean() for b in blocks])
    return a.mean(), (bm.std(ddof=1) / np.sqrt(len(bm)) if len(bm) > 1 else 0.0)

rows = []
print(" carbon   chain-2 S_CD     chain-3 S_CD")
for n in carbons:
    m2, e2 = stat(S[(2, n)]); m3, e3 = stat(S[(3, n)])
    rows.append((n, m2, e2, m3, e3))
    print(f"  C{n:<3d}   {m2:6.3f} +/- {e2:.3f}   {m3:6.3f} +/- {e3:.3f}")
plateau = np.mean([(r[1] + r[3]) / 2 for r in rows if 3 <= r[0] <= 8])
print(f"plateau (mean of carbons 3-8, both chains): {plateau:.3f}")
np.savetxt(f"{tag}.csv", np.array(rows), delimiter=",", header="carbon,chain2,err2,chain3,err3", comments="")
try:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    r = np.array(rows)
    plt.figure(figsize=(5, 3.6))
    plt.errorbar(r[:, 0], r[:, 1], r[:, 2], marker="o", label="chain 2")
    plt.errorbar(r[:, 0], r[:, 3], r[:, 4], marker="s", label="chain 3")
    plt.xlabel("carbon number"); plt.ylabel("|S_CD|"); plt.ylim(0, 0.5); plt.legend(); plt.tight_layout()
    plt.savefig(f"{tag}.png", dpi=200)
except Exception as e:
    print("plot skipped:", e)
