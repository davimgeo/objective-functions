import numpy as np
import matplotlib.pyplot as plt

SIZE = 1001

recoved_ricker = np.fromfile("data/recovered_wavelet.bin", dtype=np.float32)
ricker = np.fromfile("data/ricker.bin", dtype=np.float32)
arr = np.fromfile("data/mag_ricker.bin", dtype=np.float32)
freq = np.fromfile("data/freq.bin", dtype=np.float32)

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10, 8))

ax[0].plot(np.arange(SIZE) * 1e-3, ricker, label=f"dt={1e-3}, fmax=30hz")
ax[0].legend()

ax[1].plot(freq[:(SIZE // 2)], arr[:(SIZE // 2)])

for axs in ax:
  axs.grid(True)

plt.show()

plt.plot(np.arange(SIZE) * 1e-3, ricker)

plt.plot(np.arange(SIZE) * 1e-3, recoved_ricker)

plt.show()
