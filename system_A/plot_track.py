import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.loadtxt("dox_track.csv", delimiter=",", skiprows=1)
t, dz, pp, dm = d.T
h = np.abs(dz)                                   # distance from the bilayer centre, either side
w = 100                                          # 100 frames = 1 ns moving average
ma = lambda a: np.convolve(a, np.ones(w) / w, mode="valid")
tm = t[w - 1:]

fig, ax = plt.subplots(2, 1, figsize=(6.5, 5), sharex=True)
ax[0].plot(t, h, lw=0.3, color="tab:blue", alpha=0.5)
ax[0].plot(tm, ma(h), lw=1.6, color="tab:blue", label="DOX+ (1 ns average)")
ax[0].axhline(pp.mean() / 2, color="k", ls="-", lw=1, label="phosphate plane (this run)")
ax[0].axhline(10.5, color="tab:red", ls="--", lw=1, label="paper, DPPC Ld: 10.5 A")
ax[0].axhline(25.5, color="tab:green", ls="--", lw=1, label="paper, DPPC So: 25.5 A")
ax[0].set_ylabel("distance from bilayer centre (A)")
ax[0].set_ylim(0, 45); ax[0].legend(fontsize=6.5, loc="lower left", ncol=2)
ax[1].plot(t, dm, lw=0.3, color="tab:orange", alpha=0.6)
ax[1].plot(tm, ma(dm), lw=1.6, color="tab:orange")
ax[1].axhline(5, color="k", ls=":", lw=1)
ax[1].set_ylabel("closest DOX-P distance (A)"); ax[1].set_xlabel("time (ns)")
plt.tight_layout(); plt.savefig("dox_track.png", dpi=200)
print("saved dox_track.png:", len(t), "frames,", f"{t[-1]:.1f} ns")
