import sys, warnings
import numpy as np
import MDAnalysis as mda
warnings.filterwarnings("ignore")

gro, xtc = sys.argv[1], sys.argv[2]
t0, t1 = float(sys.argv[3]), float(sys.argv[4])      # window in ns
stride = int(sys.argv[5]) if len(sys.argv) > 5 else 1
tag = sys.argv[6] if len(sys.argv) > 6 else "dens"

u = mda.Universe(gro, xtc.split(","))
lip   = u.select_atoms("resname DPPC CHL1")
P     = u.select_atoms("resname DPPC and name P")
dox   = u.select_atoms("resname LIG")
water = u.select_atoms("resname TIP3 and name OH2")
tails = u.select_atoms("resname DPPC and name C2* C3* and not name C2 C3 C21 C31")
groups = {"phosphate P": P, "DOX": dox, "water O": water, "lipid tails": tails}
if len(u.select_atoms("resname CHL1")):
    groups["cholesterol OH (O3)"] = u.select_atoms("resname CHL1 and name O3")

def cmean(z, L):
    th = 2*np.pi*z/L
    return (np.angle(np.exp(1j*th).mean()) % (2*np.pi))*L/(2*np.pi)
def wrap(d, L): return (d + L/2) % L - L/2

dt = u.trajectory.dt
f0 = int(round(t0*1000/dt)); f1 = min(u.trajectory.n_frames, int(round(t1*1000/dt)))
frames = list(range(f0, f1, stride))
print(f"{len(frames)} frames, {t0}-{t1} ns")

zmax, dz = 60.0, 1.0
edges = np.arange(-zmax, zmax+dz/2, dz); cen = 0.5*(edges[1:]+edges[:-1])
H = {k: np.zeros(len(cen)) for k in groups}
for fr in frames:
    ts = u.trajectory[fr]; L = ts.dimensions[2]; area = ts.dimensions[0]*ts.dimensions[1]
    zc = cmean(lip.positions[:, 2], L)
    for k, g in groups.items():
        z = wrap(g.positions[:, 2] - zc, L)
        H[k] += np.histogram(z, edges)[0]/(area*dz)*1000.0
n = len(frames)
out = np.column_stack([cen] + [H[k]/n for k in groups])
np.savetxt(f"{tag}.csv", out, delimiter=",", header="z_A,"+",".join(k.replace(",", "") for k in groups), comments="")
zA = cen
print("\nWhere each group peaks (A from bilayer centre, folded onto the + side):")
for k in groups:
    prof = (H[k] + H[k][::-1])[zA > 0]; zp = zA[zA > 0]
    if k == "DOX": prof = H[k][zA > 0] + H[k][zA < 0][::-1]
    print(f"  {k:22s} peak at {zp[np.argmax(prof)]:5.1f} A")
d = H["DOX"]; tot = d.sum()
if tot > 0:
    mz = (np.abs(zA)*d).sum()/tot
    print(f"\nDOX mean |z| = {mz:.1f} A   (phosphate peak at ~{zA[zA>0][np.argmax((H['phosphate P']+H['phosphate P'][::-1])[zA>0])]:.1f} A)")
print(f"saved {tag}.csv")
