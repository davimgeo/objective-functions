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

  img = ax.imshow(seismogram, aspect="auto", cmap="jet",
                    vmin=scale_min, vmax=scale_max)

  ax.set_yticks(tloc)
  ax.set_yticklabels(tlab)

  ax.set_xticks(xloc)
  ax.set_xticklabels(xlab)

  ax.set_xlabel("Offset (m)", fontsize=13)
  ax.set_ylabel("TWT (s)", fontsize=13)

  plt.show()


nt, nrec = 4001, 113
dt = 1e-3
offset = 15

seismogram_mag = load("data/seismogram_mag.bin", nt, nrec)

plot_seismogram(seismogram_mag, nt, dt, offset, nrec)
