import sys

import numpy as np
import matplotlib.pyplot as plt

def load(PATH: str, height: int, width: int) -> np.ndarray:
  return np.fromfile(
    PATH, dtype=np.float32, count=height*width
  ).reshape([height, width], order='C')

def plot_seismogram(
    seismogram: np.ndarray, 
    nt: int, dt: float, 
    offset: int, 
    nrec: int, 
    perc=99
) -> None:
  
  tloc = np.linspace(0, nt - 1, 11, dtype=int)
  tlab = np.around(tloc * dt, decimals=1)

  xloc = np.linspace(0, nrec - 1, 9)
  xlab = np.array(offset * xloc, dtype=int)

  scale_min = np.percentile(seismogram, 100 - perc)
  scale_max = np.percentile(seismogram, perc)

  fig, ax = plt.subplots(figsize=(10, 8))

  img = ax.imshow(seismogram, aspect="auto", cmap="Greys",
                    vmin=scale_min, vmax=scale_max)

  ax.set_yticks(tloc)
  ax.set_yticklabels(tlab)

  ax.set_xticks(xloc)
  ax.set_xticklabels(xlab)

  ax.set_xlabel("Offset (m)", fontsize=13)
  ax.set_ylabel("TWT (s)", fontsize=13)

  plt.show()

t0 = 30

if len(sys.argv) > 1:
  t0 = int(sys.argv[1])

"""
nt, nrec = 2001, 40
dt = 1e-3
offset = 5

size = 101

ref_alpha = 3

alpha_min = 0.1
alpha_max = 5

alphas = np.linspace(alpha_min, alpha_max, size)
idx = np.abs(alphas - ref_alpha).argmin()
alphas[idx] = ref_alpha

H_cor = np.fromfile(f"data/H_cor_size101_t0_{t0}.bin", dtype=np.float32, count=size)

equation = r"$\sum_{\tau=0}^{nt} \; \sum_{n=0}^{nrec} (P(\tau)c(\tau, x))^2$"

plt.figure(figsize=(6, 10))

plt.plot(H_cor, alphas, label=fr"$t0 = {t0}$")

#plt.title(f"Cross-correlation Objetive Function", fontsize=13)
plt.xlabel(r"$\H_cor$", fontsize=13)
plt.ylabel("alphas", fontsize=13)

plt.tight_layout()
plt.grid(True)
plt.legend(fontsize=13)

plt.savefig(f"H_cor_ref_alpha-{ref_alpha}.png", dpi=500)

plt.show()
"""
nt = 200

h_cor = np.fromfile("data/1d/h_cor.bin", dtype=np.float32, count=nt)
l2 = np.fromfile("data/1d/l2.bin", dtype=np.float32, count=nt)
l1 = np.fromfile("data/1d/l1.bin", dtype=np.float32, count=nt)
decon = np.fromfile("data/1d/decon.bin", dtype=np.float32, count=nt)

fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(h_cor, label="Cross-correlation")
ax.plot(l2, label="L2 Norm")
ax.plot(l1, label="L1 Norm")
ax.plot(decon, label="Deconvolution")

ax.grid(True)
ax.legend(loc="upper right")

plt.tight_layout()
plt.show()
