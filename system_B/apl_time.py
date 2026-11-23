import sys
import numpy as np

f = sys.argv[1]
block = float(sys.argv[2]) if len(sys.argv) > 2 else 2000.0   # ps per block
d = np.loadtxt(f, comments=("#", "@"))
t, x, y, z = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
apl = x * y * 100 / 64                      # A^2 per lipid (64 lipids per leaflet)
print(" time(ns)    area/lipid (A^2)    box z (nm)")
for t0 in np.arange(t[0], t[-1], block):
    m = (t >= t0) & (t < t0 + block)
    if m.sum() > 5:
        print(f"{t0/1000:5.0f}-{(t0+block)/1000:<5.0f}   {apl[m].mean():8.1f} +/- {apl[m].std():3.1f}     {z[m].mean():.3f}")
