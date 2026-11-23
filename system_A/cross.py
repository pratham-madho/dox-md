import sys, numpy as np
f = sys.argv[1] if len(sys.argv) > 1 else "dox_track.csv"
d = np.loadtxt(f, delimiter=",", skiprows=1)
t, z, pp = d[:, 0], d[:, 1], d[:, 2]
half = pp / 2
a = np.abs(z)
reg = np.where(a > half, "water", np.where(a > 10, "headgroup", np.where(a > 5, "upper tails", "core")))
n = len(z)
print(f"{n} frames, {t[0]:.1f}-{t[-1]:.1f} ns, phosphate planes at about +/-{half.mean():.1f} A\n")
for lab, desc in [("water", "in water (outside phosphate plane)"),
                  ("headgroup", "inside phosphate plane, 10 A or more from centre"),
                  ("upper tails", "5-10 A from centre (tail region)"),
                  ("core", "within 5 A of the bilayer centre")]:
    print(f"  {desc:50s} {100*np.mean(reg == lab):5.1f}%")
print(f"\ndeepest point: |height| {a.min():.1f} A at {t[a.argmin()]:.2f} ns")
inside = a < half
cross = [(t[i], z[i-1], z[i]) for i in range(1, n)
         if inside[i-1] and inside[i] and np.sign(z[i-1]) != np.sign(z[i])]
if cross:
    print(f"\nCROSSED the bilayer centre {len(cross)} time(s):")
    for c in cross[:10]: print(f"  at {c[0]:.2f} ns: {c[1]:+.1f} -> {c[2]:+.1f} A")
else:
    side = "+" if np.median(z[inside]) > 0 else "-" if inside.any() else "none"
    print(f"\nNo crossing: it never passed through the centre while inside the membrane"
          f"{'' if side == 'none' else f' (stayed in the {side} leaflet when inside)'}.")
print("(sign flips while in water are the drug moving through the periodic box, not a crossing)")
