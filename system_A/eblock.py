import sys, numpy as np
f = sys.argv[1]; win = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
d = np.loadtxt(f, comments=("#", "@"))
t = d[:, 0] / 1000.0; cols = d.shape[1] - 1
names = [l.split('"')[1] for l in open(f) if l.startswith("@ s") and "legend" in l][:cols]
print(f"{len(t)} frames, {t[0]:.1f}-{t[-1]:.1f} ns   (kJ/mol, kcal/mol in brackets)")
print("  window(ns) " + "".join(f"{n:>24s}" for n in names))
for a in np.arange(0, t[-1], win):
    m = (t >= a) & (t < a + win)
    if m.sum() < 3: continue
    print(f"  {a:4.0f}-{a+win:<4.0f}" + "".join(f"   {d[m,i+1].mean():9.1f} ({d[m,i+1].mean()/4.184:6.1f})" for i in range(cols)))
