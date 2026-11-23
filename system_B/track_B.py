import sys, warnings
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
warnings.filterwarnings("ignore")

gro, xtc = sys.argv[1], sys.argv[2]
every = int(sys.argv[3]) if len(sys.argv) > 3 else 100   # print every N frames

u = mda.Universe(gro, xtc.split(","))
P = u.select_atoms("resname DPPC and name P")
lipids = u.select_atoms("resname DPPC CHL1")
dox = u.select_atoms("resname LIG")
dox_O = u.select_atoms("resname LIG and name O*")
chol_O = u.select_atoms("resname CHL1 and name O*")      # cholesterol has one oxygen: its 3-OH
print(f"{len(P)} P atoms, {len(chol_O)} cholesterol O, {len(dox)} drug atoms ({len(dox_O)} O), {u.trajectory.n_frames} frames")
assert len(P) > 0 and len(dox) > 0 and len(chol_O) > 0 and len(dox_O) > 0, "a selection came back empty"

def cmean(z, L):
    th = 2 * np.pi * z / L
    return (np.angle(np.exp(1j * th).mean()) % (2 * np.pi)) * L / (2 * np.pi)
def wrap(d, L):
    return (d + L / 2) % L - L / 2

rows = []
try:
    for ts in u.trajectory:
        L = ts.dimensions[2]
        zc = cmean(lipids.positions[:, 2], L)
        rel = wrap(P.positions[:, 2] - zc, L)
        thick = rel[rel > 0].mean() - rel[rel < 0].mean()
        zd = cmean(dox.positions[:, 2], L)
        dz = wrap(zd - zc, L)
        dP = distance_array(dox.positions, P.positions, box=ts.dimensions).min()
        dC = distance_array(dox_O.positions, chol_O.positions, box=ts.dimensions)
        rows.append((ts.time / 1000.0, dz, thick, dP, dC.min(), (dC.min(axis=0) < 3.5).sum()))
except Exception as e:
    print("stopped reading at frame", len(rows), "(", type(e).__name__, ") - probably the file is still being written")

r = np.array(rows)
np.savetxt("dox_track_B.csv", r, delimiter=",", header="t_ns,dz_A,PP_A,mindist_P_A,mindist_CholO_A,n_CholO_3p5", comments="")
print(" t(ns)   DOX height (A)   P-P (A)   closest P (A)   closest Chol-O (A)   Chol-OH within 3.5 A")
for x in r[::every]:
    print(f"{x[0]:6.1f}   {x[1]:+8.1f}   {x[2]:8.1f}   {x[3]:8.1f}   {x[4]:12.1f}   {int(x[5]):12d}")
t, dz, pp, dP, dC, nC = r.T
h = np.abs(dz)
print(f"\nclosest phosphorus < 5 A in {100*(dP < 5).mean():.0f}% of frames")
print(f"a DOX oxygen within 3.5 A of a cholesterol OH in {100*(nC > 0).mean():.0f}% of frames")
print(f"below the phosphate plane in {100*(h < pp / 2).mean():.0f}% of frames; smallest |height| {h.min():.1f} A at {t[h.argmin()]:.1f} ns")
print("by 10 ns window (contact with P | with Chol-OH | mean |height| A):")
for a in np.arange(0, t[-1], 10):
    m = (t >= a) & (t < a + 10)
    if m.sum() > 5:
        print(f"  {a:4.0f}-{a+10:<4.0f}  {100*(dP[m] < 5).mean():4.0f}%  {100*(nC[m] > 0).mean():4.0f}%  {h[m].mean():6.1f}")
