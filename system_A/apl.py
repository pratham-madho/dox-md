import sys
import numpy as np

NLIP = 64  # lipids per leaflet
for f in sys.argv[1:]:
    d = np.loadtxt(f, comments=("#", "@"))
    t, x, y, z = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    apl = x * y * 100 / NLIP          # nm^2 -> A^2, per lipid in one leaflet
    h = len(t) // 2
    print(f"{f}: {t[0]:.0f}-{t[-1]:.0f} ps | area/lipid: all {apl.mean():.1f}, "
          f"last half {apl[h:].mean():.1f} +/- {apl[h:].std():.1f} A^2 | box z {z.mean():.3f} nm")
