import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter

def load(PATH: str, height: int, width: int) -> np.ndarray:
  return np.fromfile(
    PATH, dtype=np.float32, count=height*width
  ).reshape([height, width], order='C')

def plot_objf(objf: np.ndarray, trace_size: int, label=None, save=False) -> None:
  if label == None:
    label = "Objective Funcion"

  fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

  ax[0].imshow(objf)
  line_trace, = ax[0].plot([], [], 'b--')

  line_objf, = ax[1].plot([], [], label=label)
  line_l2, = ax[1].plot([], [], label="L2 Norm")

  ax[1].set_xlim(0, trace_size - 1)
  ax[1].set_ylim(
      min(objf.min(), l2.min()),
      max(objf.max(), l2.max())
  )

  def update(trace):
    line_trace.set_data(
        np.arange(trace_size),
        np.full(trace_size, trace)
    )

    line_objf.set_data(np.arange(trace_size), objf[trace, :])
    line_l2.set_data(np.arange(trace_size), l2[trace, :])

    return line_trace, line_objf, line_l2  

  ani = FuncAnimation(
    fig,
    update,
    frames=trace_size,
    interval=300,
    blit=True,
    repeat=True,
  )

  ax[1].grid(True)
  ax[1].legend(loc="upper right")

  plt.tight_layout()
  if save:
    ani.save("animation.mp4", writer="ffmpeg", fps=20, dpi=200)
  plt.show()

PATH = "data/temp/"

result_size = 200

cross = load("data/temp/cross.bin", result_size, result_size)
l2 = load("data/temp/l2.bin", result_size, result_size)
l1 = load("data/temp/l1.bin", result_size, result_size)
decon = load("data/temp/decon.bin", result_size, result_size)

plot_objf(decon, result_size, label="Deconvution", save=True)


