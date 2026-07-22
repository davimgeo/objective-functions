import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

# Parâmetros
nt = 1001
dt = 1e-3
fmax = 30.0

# Wavelet original
wavelet = np.fromfile(
    "data/ricker2.bin",
    dtype=np.float32,
    count=nt
)

time = np.arange(nt) * dt


def get_ricker(nt, fmax, dt, phase):
    t0 = 2.0 * np.pi / fmax

    fc = fmax / (3.0 * np.sqrt(np.pi))

    ricker = np.zeros(nt, dtype=np.float32)

    for i in range(nt):
        t = (i * dt) - t0 - phase

        arg = np.pi * (
            np.pi * np.pi * fc * fc * t * t
        )

        ricker[i] = (
            (1.0 - 2.0 * arg)
            * np.exp(-arg)
        )

    return ricker


# Figura
fig, ax = plt.subplots(figsize=(10, 5))

# Wavelet original fixa
ax.plot(
    time,
    wavelet,
    label="Original Wavelet"
)

# Ricker inicial
ricker_line, = ax.plot(
    time,
    get_ricker(nt, fmax, dt, phase=0.0),
    label="Moving Ricker",
    color="red"
)

ax.set_xlabel("Time [s]", fontsize=13)
ax.set_ylabel("Amplitude", fontsize=13)

ax.grid()
plt.tight_layout()
ax.legend()


# Valores de phase
phases = np.linspace(
    0.0,
    0.7,
    200
)


def update(frame):
    phase = phases[frame]

    ricker = get_ricker(
        nt,
        fmax,
        dt,
        phase
    )

    ricker_line.set_ydata(ricker)

    #ax.set_title(
    #    f"Phase = {phase:.4f} s"
    #)

    return ricker_line,

animation = FuncAnimation(
    fig,
    update,
    frames=len(phases),
    interval=30,
    blit=True
)

writer = PillowWriter(
    fps=30
)

animation.save(
    "ricker_animation.gif",
    writer=writer
)

plt.show()

plt.show()
