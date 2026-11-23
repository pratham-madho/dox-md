import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.loadtxt("dox_track_B.csv", delimiter=",", skiprows=1)
t, dz, pp, dP, dC, nC = d.T
h = np.abs(dz)
w = max(1, min(100, len(t) // 5))             # 100 frames = 1 ns moving average (shorter if the file is still small)
ma = lambda a: np.convolve(a, np.ones(w) / w, mode="valid")
tm = t[w - 1:]

fig, ax = plt.subplots(3, 1, figsize=(6.5, 7), sharex=True)
ax[0].plot(t, h, lw=0.3, color="tab:blue", alpha=0.5)
ax[0].plot(tm, ma(h), lw=1.6, color="tab:blue", label="DOX+ (1 ns average)")
ax[0].axhline(pp.mean() / 2, color="k", lw=1, label="phosphate plane (this run)")
ax[0].axhline(21.5, color="tab:red", ls="--", lw=1, label="paper, DPPC/Chol (Lo): 21.5 A")
ax[0].set_ylabel("distance from bilayer centre (A)")
ax[0].set_ylim(0, 50); ax[0].legend(fontsize=6.5, loc="lower left")
ax[1].plot(t, dP, lw=0.3, color="tab:orange", alpha=0.6); ax[1].plot(tm, ma(dP), lw=1.6, color="tab:orange")
ax[1].axhline(5, color="k", ls=":", lw=1); ax[1].set_ylabel("closest DOX-P (A)")
ax[2].plot(t, dC, lw=0.3, color="tab:green", alpha=0.6); ax[2].plot(tm, ma(dC), lw=1.6, color="tab:green")
ax[2].axhline(3.5, color="k", ls=":", lw=1); ax[2].set_ylabel("closest DOX O - Chol OH (A)"); ax[2].set_xlabel("time (ns)")
plt.tight_layout(); plt.savefig("dox_track_B.png", dpi=200)
print("saved dox_track_B.png:", len(t), "frames,", f"{t[-1]:.1f} ns")
