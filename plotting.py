import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pathlib import Path

# ----------------------------
# Settings
# ----------------------------

DT = 60.0                  # Simulation timestep in seconds
HISTORY_LENGTH = 30000     # Number of previous points to show in Earth's trail
LIMIT = 1.6e11             # Fixed camera range in metres
FRAME_SKIP = 1250      # Increase this to make the animation faster
FIGURE_DIR = Path("figure")
FIGURE_DIR.mkdir(exist_ok=True)

# ----------------------------
# Load simulation data
# ----------------------------


body_history = pd.read_csv("body_history_rk4.csv")
earth_history = (
    body_history[body_history["body"] == "Earth"]
    .reset_index(drop=True)
)

sun_history = (
    body_history[body_history["body"] == "Sun"]
    .reset_index(drop=True)
)

if earth_history.empty or sun_history.empty:
    raise ValueError("CSV must contain rows for both Earth and Sun.")

if len(earth_history) != len(sun_history):
    raise ValueError("Earth and Sun must have the same number of timesteps.")

simulation_years = (
    (len(earth_history) - 1) * DT / (60 * 60 * 24 * 365.25)
)
time_array= np.arange(0, simulation_years, DT / (60 * 60 * 24 * 365.25))
sun_velocity = np.sqrt(sun_history["vx"]**2 + sun_history["vy"]**2)
max_sun_velocity = np.max(sun_velocity)
sun_speed = np.sqrt(
    sun_history["vx"]**2 +
    sun_history["vy"]**2 +
    sun_history["vz"]**2
)


# ----------------------------
# Figure setup
# ----------------------------

plt.style.use("dark_background")

fig, ax = plt.subplots(figsize=(10, 10))
ax.set_xlim(-LIMIT, LIMIT)
ax.set_ylim(-LIMIT, LIMIT)
ax.set_aspect("equal")
ax.set_xlabel("X position (m)")
ax.set_ylabel("Y position (m)")
ax.set_title("Earth-Sun N-body Simulation")
ax.grid(alpha=0.2)

# ----------------------------
# Initial plot objects
# ----------------------------

# Sun and Earth markers
sun_point, = ax.plot([], [], "o", color="orange", markersize=12, label="Sun")
earth_point, = ax.plot([], [], "o", color="dodgerblue", markersize=7, label="Earth")

# Earth's trail
earth_trail, = ax.plot(
    [], [],
    color="dodgerblue",
    alpha=0.65,
    linewidth=1.5,
)

# Velocity and acceleration arrows
sun_velocity_arrow = None
sun_acceleration_arrow = None
earth_velocity_arrow = None
earth_acceleration_arrow = None

ax.legend(loc="upper right")

# ----------------------------
# Arrow helper
# ----------------------------

def make_arrow(x, y, vx, vy, color, length):
    magnitude = np.hypot(vx, vy)

    if magnitude == 0:
        return None

    dx = (vx / magnitude) * length
    dy = (vy / magnitude) * length

    return ax.quiver(
        x, y, dx, dy,
        angles="xy",
        scale_units="xy",
        scale=1,
        color=color,
        width=0.004
    )

# ----------------------------
# Animation update
# ----------------------------

def update(frame):
    global sun_velocity_arrow
    global sun_acceleration_arrow
    global earth_velocity_arrow
    global earth_acceleration_arrow

    # Get current positions
    earth = earth_history.iloc[frame]
    sun = sun_history.iloc[frame]

    # Update body markers
    earth_point.set_data([earth["x"]], [earth["y"]])
    sun_point.set_data([sun["x"]], [sun["y"]])

    # Update Earth's trail
    start = max(0, frame - HISTORY_LENGTH)

    earth_trail.set_data(
        earth_history["x"].iloc[start:frame + 1],
        earth_history["y"].iloc[start:frame + 1]
    )

    # Remove previous arrows
    for arrow in (
        sun_velocity_arrow,
        sun_acceleration_arrow,
        earth_velocity_arrow,
        earth_acceleration_arrow
    ):
        if arrow is not None:
            arrow.remove()

    # Earth's velocity arrow
    earth_velocity_arrow = make_arrow(
        earth["x"], earth["y"],
        earth["vx"], earth["vy"],
        "green", 1.5e10
    )

    # Earth's acceleration arrow
    earth_acceleration_arrow = make_arrow(
        earth["x"], earth["y"],
        earth["ax"], earth["ay"],
        "red", 1.5e10
    )


    # Display elapsed simulation time
    elapsed_days = frame * DT / 86400

    ax.set_title(
        f"Earth-Sun N-body Simulation — {elapsed_days:.2f} days"
        f"\n({simulation_years:.2f} years, dt={DT:g} s)"
    )

    return (
        earth_point,
        sun_point,
        earth_trail
    )

# ----------------------------
# Run animation
# ----------------------------

frames = range(0, len(earth_history), FRAME_SKIP)

animation = FuncAnimation(
    fig,
    update,
    frames=frames,
    interval=20,
    blit=False,
    repeat=False,
    cache_frame_data=False
)

gif_path = FIGURE_DIR / (
    f"earth_sun_simulation_{simulation_years:.2f}years_dt{DT:g}s.gif"
)
animation.save(gif_path, writer="pillow", fps=10)

