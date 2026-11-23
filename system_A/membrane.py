import sys, warnings
import numpy as np
import MDAnalysis as mda
warnings.filterwarnings("ignore")
gro, xtc = sys.argv[1], sys.argv[2]
t0, t1 = float(sys.argv[3]), float(sys.argv[4])
stride = int(sys.argv[5]) if len(sys.argv) > 5 else 20
u = mda.Universe(gro, xtc.split(","))
P = u.select_atoms("resname DPPC and name P")
lip = u.select_atoms("resname DPPC CHL1")
nd, nc = len(u.select_atoms("resname DPPC").residues), len(u.select_atoms("resname CHL1").residues)
def cmean(z, L):
    th = 2*np.pi*z/L
    return (np.angle(np.exp(1j*th).mean()) % (2*np.pi))*L/(2*np.pi)
def wrap(d, L): return (d + L/2) % L - L/2
dt = u.trajectory.dt
f0 = int(round(t0*1000/dt)); f1 = min(u.trajectory.n_frames, int(round(t1*1000/dt)))
A, T = [], []
for ts in u.trajectory[f0:f1:stride]:
    Lx, Ly, Lz = ts.dimensions[:3]
    z = wrap(P.positions[:, 2] - cmean(lip.positions[:, 2], Lz), Lz)
    T.append(z[z > 0].mean() - z[z < 0].mean())
    A.append(Lx * Ly)
A, T = np.array(A), np.array(T)
def blk(x):
    b = np.array_split(x, 5); m = np.array([bb.mean() for bb in b])
    return x.mean(), m.std(ddof=1) / np.sqrt(5)
per_leaf = (nd + nc) / 2
(a, ea), (t, et) = blk(A / per_leaf), blk(T)
print(f"{len(A)} frames, {t0}-{t1} ns | {nd} DPPC + {nc} cholesterol")
print(f"  area per molecule (all lipids incl. cholesterol): {a:6.2f} +/- {ea:.2f} A^2")
if nc:
    (ad, ead) = blk(A / (nd / 2))
    print(f"  box area per DPPC (ignores cholesterol):         {ad:6.2f} +/- {ead:.2f} A^2")
print(f"  P-P thickness:                                   {t:6.2f} +/- {et:.2f} A")
print(f"  box x, y (A): {np.sqrt(A).mean():.2f}   (errors from 5 time blocks)")
