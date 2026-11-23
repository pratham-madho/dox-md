import sys, warnings
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis as HBA
warnings.filterwarnings("ignore")
gro, xtc = sys.argv[1], sys.argv[2]
t0, t1 = float(sys.argv[3]), float(sys.argv[4])
stride = int(sys.argv[5]) if len(sys.argv) > 5 else 10
u = mda.Universe(gro, xtc.split(","))
dt = u.trajectory.dt
f0 = int(round(t0*1000/dt)); f1 = min(u.trajectory.n_frames, int(round(t1*1000/dt)))
nfr = len(range(f0, f1, stride))
dox_don = "resname LIG and (name N or name O3 O4 O5 O6 O8)"
acc = ("(resname DPPC and name O11 O12 O13 O14 O21 O22 O31 O32)"
       " or (resname CHL1 and name O3)"
       " or (resname TIP3 and name OH2 and around 6 resname LIG)")
runs = {"DOX donates": dict(donors_sel=dox_don, hydrogens_sel="resname LIG and name H*", acceptors_sel=acc)}
if len(u.select_atoms("resname CHL1")):
    runs["cholesterol OH donates to DOX"] = dict(donors_sel="resname CHL1 and name O3",
        hydrogens_sel="resname CHL1 and name H*", acceptors_sel="resname LIG and name O*")
def cls(a):
    if a.resname == "TIP3": return "water"
    if a.resname == "CHL1": return "cholesterol OH"
    if a.resname == "LIG": return "DOX oxygen"
    return "phosphate O" if a.name in ("O11", "O12", "O13", "O14") else "ester O"
print(f"{nfr} frames, {t0}-{t1} ns (every {stride}); criterion D...A < 3.5 A, D-H...A > 150 deg\n")
for title, kw in runs.items():
    h = HBA(u, d_a_cutoff=3.5, d_h_a_angle_cutoff=150, d_h_cutoff=1.2, update_selections=True, **kw)
    h.run(start=f0, stop=f1, step=stride)
    counts = {}
    for row in h.results.hbonds:
        a = u.atoms[int(row[3])]; dn = u.atoms[int(row[1])]
        key = (("NH3+ -> " if dn.name == "N" else "DOX OH -> ") + cls(a)) if title.startswith("DOX") else "DOX oxygen"
        counts[key] = counts.get(key, 0) + 1
    print(f"  {title}:")
    if not counts: print("      none found")
    for k in sorted(counts, key=lambda k: -counts[k]):
        print(f"      {k:32s} {counts[k]/nfr:5.2f} per frame")
