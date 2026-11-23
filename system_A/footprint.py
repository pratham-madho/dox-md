import sys, warnings
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
warnings.filterwarnings("ignore")
gro, xtc = sys.argv[1], sys.argv[2]
t0, t1 = float(sys.argv[3]), float(sys.argv[4])
stride = int(sys.argv[5]) if len(sys.argv) > 5 else 10
R = float(sys.argv[6]) if len(sys.argv) > 6 else 10.0
FAR = 15.0
u = mda.Universe(gro, xtc.split(","))
dppc = u.select_atoms("resname DPPC").residues
P = u.select_atoms("resname DPPC and name P")
lipids = u.select_atoms("resname DPPC CHL1")
dox = u.select_atoms("resname LIG and not name H*")
chlO = u.select_atoms("resname CHL1 and name O3")
assert len(P) == len(dppc), "one P per DPPC expected"
pairs = []
for n in range(3, 9):
    for c, hs in ((2, (f"H{n}R", f"H{n}S")), (3, (f"H{n}X", f"H{n}Y"))):
        C = u.select_atoms(f"resname DPPC and name C{c}{n}")
        for h in hs:
            H = u.select_atoms(f"resname DPPC and name {h}")
            assert len(C) == len(H) == len(dppc), f"selection problem C{c}{n}/{h}"
            pairs.append((C, H))
def cmean(z, L):
    th = 2*np.pi*z/L
    return (np.angle(np.exp(1j*th).mean()) % (2*np.pi))*L/(2*np.pi)
def wrap(d, L): return (d + L/2) % L - L/2
dt = u.trajectory.dt
f0 = int(round(t0*1000/dt)); f1 = min(u.trajectory.n_frames, int(round(t1*1000/dt)))
rows = []; nskip = 0
for fr in range(f0, f1, stride):
    ts = u.trajectory[fr]; box = ts.dimensions; Lx, Ly, Lz = box[:3]
    if distance_array(dox.positions, P.positions, box=box).min() > 5.0:
        nskip += 1; continue
    zc = cmean(lipids.positions[:, 2], Lz)
    side = np.sign(wrap(cmean(dox.positions[:, 2], Lz) - zc, Lz))
    dx = np.array([cmean(dox.positions[:, 0], Lx), cmean(dox.positions[:, 1], Ly)])
    hP = side * wrap(P.positions[:, 2] - zc, Lz)
    mine = hP > 0
    lat = np.hypot(wrap(P.positions[:, 0] - dx[0], Lx), wrap(P.positions[:, 1] - dx[1], Ly))
    S = np.zeros(len(dppc))
    for C, H in pairs:
        v = H.positions - C.positions; v -= box[:3] * np.round(v / box[:3])
        S += -(3 * v[:, 2]**2 / (v**2).sum(1) - 1) / 2
    S /= len(pairs)
    near, far = mine & (lat < R), mine & (lat > FAR)
    if near.sum() == 0 or far.sum() == 0: nskip += 1; continue
    row = [S[near].mean(), S[far].mean(), hP[near].mean(), hP[far].mean()]
    if len(chlO):
        hc = side * wrap(chlO.positions[:, 2] - zc, Lz)
        latc = np.hypot(wrap(chlO.positions[:, 0] - dx[0], Lx), wrap(chlO.positions[:, 1] - dx[1], Ly))
        cn, cf = ((hc > 0) & (latc < R)).sum(), ((hc > 0) & (latc > FAR)).sum()
        row += [cn / (cn + near.sum()), cf / (cf + far.sum())]
    rows.append(row)
r = np.array(rows)
print(f"{len(r)} frames with the drug touching the membrane, {nskip} skipped ({t0}-{t1} ns, every {stride})")
print(f"'near' = lipids in the drug's leaflet within {R:.0f} A laterally of the drug; 'far' = beyond {FAR:.0f} A\n")
def blk(x):
    b = np.array_split(x, 5); m = np.array([bb.mean() for bb in b if len(bb)])
    return x.mean(), m.std(ddof=1) / np.sqrt(len(m))
labels = ["tail order S_CD (C3-C8)", "phosphate height (A)"] + (["cholesterol fraction"] if r.shape[1] > 4 else [])
for i, lab in enumerate(labels):
    (mn, en), (mf, ef) = blk(r[:, 2*i]), blk(r[:, 2*i+1])
    print(f"  {lab:25s} near {mn:6.3f} +/- {en:.3f}   far {mf:6.3f} +/- {ef:.3f}   difference {mn-mf:+.3f}")
