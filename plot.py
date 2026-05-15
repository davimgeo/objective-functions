import numpy as np
import matplotlib.pyplot as plt

SIZE = 1001

ricker = np.fromfile("output/ricker.bin", dtype=np.float32)
arr = np.fromfile("output/mag_ricker.bin", dtype=np.float32)
freq = np.fromfile("output/freq.bin", dtype=np.float32)

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10, 8))

ax[0].plot(np.arange(SIZE) * 1e-3, ricker, label=f"dt={1e-3}, fmax=30hz")
ax[0].legend()

ax[1].plot(freq[:(SIZE // 2)], arr[:(SIZE // 2)])

plt.show()
