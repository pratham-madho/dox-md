import sys, numpy as np
f = sys.argv[1] if len(sys.argv) > 1 else "dox_track_B.csv"
cut = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
minlen = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
T = float(sys.argv[4]) if len(sys.argv) > 4 else 303.15
d = np.loadtxt(f, delimiter=",", skiprows=1)
t, dmin = d[:, 0], d[:, 3]
dt = np.median(np.diff(t))
raw = dmin < cut
def runs(s):
    out, start = [], 0
    for i in range(1, len(s) + 1):
        if i == len(s) or s[i] != s[start]:
            out.append([bool(s[start]), start, i]); start = i
    return out
s = raw.copy()
while True:
    r = runs(s)
    short = [k for k, (st, a, b) in enumerate(r) if (b - a) * dt < minlen]
    if not short or len(r) == 1: break
    k = min(short, key=lambda k: r[k][2] - r[k][1])
    nbr = [r[j] for j in (k - 1, k + 1) if 0 <= j < len(r)]
    st = max(nbr, key=lambda x: x[2] - x[1])[0]
    s[r[k][1]:r[k][2]] = st
r = runs(s)
pb = s.mean()
print(f"{len(t)} frames, {t[0]:.1f}-{t[-1]:.1f} ns | bound = DOX within {cut} A of a phosphorus, blips under {minlen} ns ignored\n")
print("  state      from (ns)   to (ns)   length (ns)")
for st, a, b in r:
    print(f"  {'BOUND  ' if st else 'unbound'}   {t[a]:8.1f}  {t[b-1]+dt:8.1f}   {(b-a)*dt:8.1f}")
inner = r[1:-1]
on  = [(b-a)*dt for st, a, b in inner if st]
off = [(b-a)*dt for st, a, b in inner if not st]
print(f"\nfraction of time bound: {100*pb:.0f}%")
nb = sum(1 for st, a, b in r if st)
print(f"binding events: {nb}  (complete bound periods: {len(on)}, complete unbound periods: {len(off)})")
if on:  print(f"mean residence time at the surface: {np.mean(on):.1f} ns  (from {len(on)} complete periods)")
if off: print(f"mean time to rebind from water:     {np.mean(off):.1f} ns  (from {len(off)} complete periods)")
if 0 < pb < 1:
    K = pb / (1 - pb); kT = 0.0083145 * T
    print(f"\napparent K = bound/unbound = {K:.2f}  ->  dG_apparent = -kT ln K = {-kT*np.log(K):+.1f} kJ/mol ({-kT*np.log(K)/4.184:+.1f} kcal/mol) at {T} K")
    print("(apparent: valid only for this box size and one drug molecule; few events, so order-of-magnitude)")
