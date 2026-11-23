import sys, warnings
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
warnings.filterwarnings("ignore")

gro, xtc = sys.argv[1], sys.argv[2]
t0, t1 = float(sys.argv[3]), float(sys.argv[4])
stride = int(sys.argv[5]) if len(sys.argv) > 5 else 1
tag = sys.argv[6] if len(sys.argv) > 6 else "rdf"

u = mda.Universe(gro, xtc.split(","))
N = u.select_atoms("resname LIG and name N")
refs = {"phosphate O":  u.select_atoms("resname DPPC and name O11 O12 O13 O14"),
        "ester O":      u.select_atoms("resname DPPC and name O21 O22 O31 O32"),
        "water O":      u.select_atoms("resname TIP3 and name OH2")}
if len(u.select_atoms("resname CHL1")):
    refs["cholesterol OH"] = u.select_atoms("resname CHL1 and name O3")
assert len(N) == 1, "did not find a single DOX nitrogen named N"
for k, g in refs.items(): assert len(g) > 0, f"empty selection: {k}"

dt = u.trajectory.dt
f0 = int(round(t0*1000/dt)); f1 = min(u.trajectory.n_frames, int(round(t1*1000/dt)))
frames = list(range(f0, f1, stride))
rmax, dr = 12.0, 0.1
edges = np.arange(0, rmax+dr/2, dr); r = 0.5*(edges[1:]+edges[:-1])
shell = 4/3*np.pi*(edges[1:]**3 - edges[:-1]**3)
cnt = {k: np.zeros(len(r)) for k in refs}; rho = {k: 0.0 for k in refs}
for fr in frames:
    ts = u.trajectory[fr]; V = np.prod(ts.dimensions[:3])
    for k, g in refs.items():
        d = distance_array(N.positions, g.positions, box=ts.dimensions)[0]
        cnt[k] += np.histogram(d, edges)[0]
        rho[k] += len(g)/V
n = len(frames)
res = [r]
print(f"{n} frames, {t0}-{t1} ns   (NH3+ nitrogen vs each oxygen group)")
for k in refs:
    g = cnt[k]/n/(shell*rho[k]/n); res.append(g)
    i = np.argmax(g[r > 1.5]) + np.argmax(r > 1.5)
    cut = 4.0 if k != "water O" else 3.6
    cn = (cnt[k]/n)[r < cut].sum()
    print(f"  {k:15s} first peak at {r[i]:.2f} A  (g = {g[i]:.1f});  mean # within {cut} A = {cn:.2f}")
np.savetxt(f"{tag}.csv", np.column_stack(res), delimiter=",", header="r_A,"+",".join(refs), comments="")
print(f"saved {tag}.csv")
