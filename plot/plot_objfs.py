import numpy as np
import matplotlib.pyplot as plt

def plot_objf(objf: np.ndarray, trace_size: int, trace: int, label=None) -> None:
  if label == None:
    label = "Objective Funcion"

  fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(10, 5))

  ax[0].imshow(objf)
  ax[0].plot(trace*np.ones(trace_size), np.arange(trace_size), 'b--')

  ax[1].plot(objf[:, trace], label=label)
  ax[1].plot(l2[:, trace], label="L2 Norm")
  ax[1].grid(True)

  ax[1].legend(loc="upper right")

  plt.tight_layout()
  plt.show()

PATH = "data/temp/"

result_size = 200

cross = np.zeros((result_size, result_size))
decon = np.zeros((result_size, result_size))
l1 = np.zeros((result_size, result_size))
l2 = np.zeros((result_size, result_size))

for i in range(result_size):

  cross[:, i] = np.fromfile(PATH + f"h_cor_{i}.bin", dtype=np.float32, count=result_size)
  decon[:, i] = np.fromfile(PATH + f"decon_{i}.bin", dtype=np.float32, count=result_size)
  l2[:, i] = np.fromfile(PATH + f"l2_{i}.bin", dtype=np.float32, count=result_size)
  l1[:, i] = np.fromfile(PATH + f"l1_{i}.bin", dtype=np.float32, count=result_size)

plot_objf(cross, trace_size=result_size, trace=182)
