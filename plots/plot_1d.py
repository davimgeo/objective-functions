import numpy as np
import matplotlib.pyplot as plt

nt = 200

h_cor = np.fromfile("data/1d/h_cor.bin", dtype=np.float32, count=nt)
l2 = np.fromfile("data/1d/l2.bin", dtype=np.float32, count=nt)
l1 = np.fromfile("data/1d/l1.bin", dtype=np.float32, count=nt)
decon = np.fromfile("data/1d/decon.bin", dtype=np.float32, count=nt)

fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(h_cor, label="Cross-correlation")
ax.plot(l2, label="L2 Norm")
ax.plot(l1, label="L1 Norm")
#ax.plot(decon, label="Deconvolution")

ax.grid(True)
ax.legend(loc="upper right")

plt.tight_layout()
plt.show()
