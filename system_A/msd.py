import sys, warnings
import numpy as np
import MDAnalysis as mda
warnings.filterwarnings("ignore")
gro, xtc, segs = sys.argv[1], sys.argv[2], sys.argv[3]
label = sys.argv[4] if len(sys.argv) > 4 else ""
u = mda.Universe(gro, xtc.split(","))
dox = u.select_atoms("resname LIG and not name H*")
P = u.select_atoms("resname DPPC and name P")
dt = u.trajectory.dt
def cmean(z, L):
    th = 2*np.pi*z/L
    return (np.angle(np.exp(1j*th).mean()) % (2*np.pi))*L/(2*np.pi)
def wrap(d, L): return (d + L/2) % L - L/2
def track(f0, f1):
    out, prevd, prevP, acc, drift = [], None, None, None, np.zeros(2)
    for ts in u.trajectory[f0:f1]:
        L2 = ts.dimensions[:2]
        d = np.array([cmean(dox.positions[:, 0], L2[0]), cmean(dox.positions[:, 1], L2[1])])
        pP = P.positions[:, :2].copy()
        if prevd is None: acc = d.copy()
        else:
            acc = acc + wrap(d - prevd, L2)
            drift = drift + wrap(pP - prevP, L2).mean(0)
        prevd, prevP = d, pP
        out.append(acc - drift)
    return np.array(out)
print(f"lateral (xy) diffusion of DOX {label} - MSD = 4 D t, fitted between 0.1 ns and 1/5 of each window")
print("  window (ns)    D (10^-7 cm^2/s, = A^2/ns)   fit range (ns)")
for s in segs.split(","):
    a, b = map(float, s.split("-"))
    f0, f1 = int(round(a*1000/dt)), min(u.trajectory.n_frames, int(round(b*1000/dt)))
    x = track(f0, f1); n = len(x)
    lags = np.unique(np.linspace(1, n // 5, 40).astype(int))
    msd = np.array([((x[l:] - x[:-l])**2).sum(1).mean() for l in lags])
    tns = lags * dt / 1000.0
    m = tns >= 0.1
    if m.sum() < 3: print(f"  {s:12s}  window too short"); continue
    D = np.polyfit(tns[m], msd[m], 1)[0] / 4
    print(f"  {s:12s}   {D:9.2f}                    {tns[m][0]:.2f}-{tns[m][-1]:.2f}")
print("(One window = one trajectory piece, so treat each D as +/- 30-50%.)")
