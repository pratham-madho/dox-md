import sys, warnings
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
warnings.filterwarnings("ignore")

gro, xtc = sys.argv[1], sys.argv[2]
t0, t1 = float(sys.argv[3]), float(sys.argv[4])          # ns
stride = int(sys.argv[5]) if len(sys.argv) > 5 else 10   # every Nth saved frame
win = float(sys.argv[6]) if len(sys.argv) > 6 else 10.0  # ns per window

u = mda.Universe(gro, xtc.split(","))
lipids = u.select_atoms("resname DPPC CHL1")
heavy = u.select_atoms("resname LIG and not name H*")
N = u.select_atoms("resname LIG and name N")
assert len(N) == 1, "expected exactly one atom named N in LIG"

def cmean(z, L):
    th = 2 * np.pi * z / L
    return (np.angle(np.exp(1j * th).mean()) % (2 * np.pi)) * L / (2 * np.pi)
def wrap(d, L):
    return (d + L / 2) % L - L / 2

dt = u.trajectory.dt
f0 = int(round(t0 * 1000 / dt)); f1 = min(u.trajectory.n_frames, int(round(t1 * 1000 / dt)))

rows = []
for fr in range(f0, f1, stride):
    ts = u.trajectory[fr]
    L = ts.dimensions[2]
    zc = cmean(lipids.positions[:, 2], L)                # robust: centre of all lipid atoms
    zcom = cmean(heavy.positions[:, 2], L)
    s = np.sign(wrap(zcom - zc, L))                      # +1: drug on the upper side, -1: lower side
    d = distance_array(N.positions, heavy.positions, box=ts.dimensions)[0]
    j = d.argmax()                                       # heavy atom farthest from the NH3+ nitrogen
    hN = s * wrap(N.positions[0, 2] - zc, L)             # heights are positive on the drug's own water side
    hC = s * wrap(zcom - zc, L)
    hF = s * wrap(heavy.positions[j, 2] - zc, L)
    cos_in = -s * wrap(heavy.positions[j, 2] - N.positions[0, 2], L) / d[j]   # N->far atom vs the inward normal
    rows.append((ts.time / 1000, hN, hC, hF, cos_in, d[j]))
r = np.array(rows)
print(f"{len(r)} frames from {t0} to {t1} ns")
print(" window(ns)   N height   centre height   far-end height   angle of N->far end vs inward normal   N-far length")
for a in np.arange(t0, t1, win):
    m = (r[:, 0] >= a) & (r[:, 0] < a + win)
    if m.sum() < 3: continue
    ang = np.degrees(np.arccos(np.clip(r[m, 4].mean(), -1, 1)))
    print(f"  {a:4.0f}-{a+win:<4.0f}   {r[m,1].mean():6.1f}      {r[m,2].mean():6.1f}         {r[m,3].mean():6.1f}              {ang:6.0f} deg                       {r[m,5].mean()/10:.2f} nm")
