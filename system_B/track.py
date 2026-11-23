import sys, warnings
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
warnings.filterwarnings("ignore")

gro, xtc = sys.argv[1], sys.argv[2]
every = int(sys.argv[3]) if len(sys.argv) > 3 else 100   # print every N frames

u = mda.Universe(gro, xtc)
P   = u.select_atoms("resname DPPC and name P")
lipids = u.select_atoms("resname DPPC CHL1")
dox = u.select_atoms("resname LIG")
print(f"{len(P)} P atoms, {len(dox)} drug atoms, {u.trajectory.n_frames} frames")

def cmean(z, L):                 # mean of positions on a periodic axis
    th = 2 * np.pi * z / L
    return (np.angle(np.exp(1j * th).mean()) % (2 * np.pi)) * L / (2 * np.pi)

def wrap(d, L):                  # wrap a displacement to [-L/2, L/2)
    return (d + L / 2) % L - L / 2

rows = []
try:
    for ts in u.trajectory:
        L = ts.dimensions[2]
        zc = cmean(lipids.positions[:, 2], L)             # robust: centre of all lipid atoms
        rel = wrap(P.positions[:, 2] - zc, L)
        thick = rel[rel > 0].mean() - rel[rel < 0].mean()   # P-P thickness
        zd = cmean(dox.positions[:, 2], L)
        dz = wrap(zd - zc, L)
        dmin = distance_array(dox.positions, P.positions, box=ts.dimensions).min()
        rows.append((ts.time / 1000.0, dz, thick, dmin))
except Exception as e:
    print("stopped reading at frame", len(rows), "(", type(e).__name__, ") - probably the file is still being written")

rows = np.array(rows)
np.savetxt("dox_track.csv", rows, delimiter=",", header="t_ns,dz_A,PP_A,mindist_P_A", comments="")
print(" t(ns)   DOX height from bilayer centre (A)   P-P thickness (A)   closest DOX-P distance (A)")
for r in rows[::every]:
    print(f"{r[0]:6.1f}   {r[1]:+8.1f}   {r[2]:8.1f}   {r[3]:8.1f}")
