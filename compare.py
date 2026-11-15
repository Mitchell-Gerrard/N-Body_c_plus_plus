# Compare Euler and RK4 position differences over time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

DT = 60
Steps =  5256000
years = Steps * DT / (60 * 60 * 24 * 365.25)  # Convert steps to years
FIGURE_DIR = Path("figure")
FIGURE_DIR.mkdir(exist_ok=True)


def save_fig(name):
    """Save the current figure with its simulation duration and timestep."""
    figure = plt.gcf()
    metadata = f"{years:.2f} years, dt={DT} s"
    figure.axes[0].set_title(f"{figure.axes[0].get_title()}\n({metadata})")
    filename = f"{name}_{years:.2f}years_dt{DT}s.png"
    figure.savefig(FIGURE_DIR / filename, dpi=300, bbox_inches="tight")


print(f"Simulation duration: {years:.2f} years")
body_history_euler = pd.read_csv("body_history_euler.csv")
body_history_rk4 = pd.read_csv("body_history_rk4.csv")
# Extract Earth data from each simulation
earth_euler = body_history_euler[body_history_euler["body"] == "Earth"].reset_index(drop=True)
earth_rk4 = body_history_rk4[body_history_rk4["body"] == "Earth"].reset_index(drop=True)

# Calculate the 3D distance between their positions at each timestep
position_difference = np.sqrt(
    (earth_euler["x"] - earth_rk4["x"])**2 +
    (earth_euler["y"] - earth_rk4["y"])**2 +
    (earth_euler["z"] - earth_rk4["z"])**2
)

# Convert metres to kilometres
position_difference_km = position_difference / 1000

# Plot
plt.figure(figsize=(10, 5))
plt.plot(earth_euler["step"], position_difference_km)

plt.title("Euler vs RK4: Earth Position Difference")
plt.xlabel("Timestep")
plt.ylabel("Position Difference (km)")
plt.grid(True, alpha=0.3)
plt.tight_layout()
save_fig("earth_position_difference")
plt.show()
# Compare orbital radius differences (Earth-Sun distance)

# Extract Sun data
sun_euler = body_history_euler[body_history_euler["body"] == "Sun"].reset_index(drop=True)
sun_rk4 = body_history_rk4[body_history_rk4["body"] == "Sun"].reset_index(drop=True)

# Calculate Earth-Sun distance for Euler
radius_euler = np.sqrt(
    (earth_euler["x"].values - sun_euler["x"].values)**2 +
    (earth_euler["y"].values - sun_euler["y"].values)**2 +
    (earth_euler["z"].values - sun_euler["z"].values)**2
)

# Calculate Earth-Sun distance for RK4
radius_rk4 = np.sqrt(
    (earth_rk4["x"].values - sun_rk4["x"].values)**2 +
    (earth_rk4["y"].values - sun_rk4["y"].values)**2 +
    (earth_rk4["z"].values - sun_rk4["z"].values)**2
)

# Difference in orbital radius (km)
radius_difference = np.abs(radius_euler - radius_rk4) / 1000

# Plot
plt.figure(figsize=(10, 5))
plt.plot(earth_euler["step"], radius_difference)

plt.title("Euler vs RK4: Difference in Orbital Radius")
plt.xlabel("Timestep")
plt.ylabel("Orbital Radius Difference (km)")
plt.grid(True, alpha=0.3)
plt.tight_layout()
save_fig("orbital_radius_difference")
plt.show()
# Compare orbital angle differences (Euler vs RK4)

# Calculate Earth's angle around the Sun for each simulation
angle_euler = np.arctan2(
    earth_euler["y"].values - sun_euler["y"].values,
    earth_euler["x"].values - sun_euler["x"].values
)

angle_rk4 = np.arctan2(
    earth_rk4["y"].values - sun_rk4["y"].values,
    earth_rk4["x"].values - sun_rk4["x"].values
)

# Unwrap angles to prevent jumps at -180 / +180 degrees
angle_euler = np.unwrap(angle_euler)
angle_rk4 = np.unwrap(angle_rk4)

# Calculate angular difference in degrees
angle_difference = np.abs(np.degrees(angle_euler - angle_rk4))

# Plot
plt.figure(figsize=(10, 5))
plt.plot(earth_euler["step"], angle_difference)

plt.title("Euler vs RK4: Orbital Angular Difference")
plt.xlabel("Timestep")
plt.ylabel("Angular Difference (degrees)")
plt.grid(True, alpha=0.3)
plt.tight_layout()
save_fig("orbital_angular_difference")
plt.show()
kinectic_energy_euler_earth = body_history_euler[body_history_euler["body"] == "Earth"]["kinetic_energy"].values
kinetic_energy_euler_sun = body_history_euler[body_history_euler["body"] == "Sun"]["kinetic_energy"].values
kinetic_energy_rk4_earth = body_history_rk4[body_history_rk4["body"] == "Earth"]["kinetic_energy"].values
kinetic_energy_rk4_sun = body_history_rk4[body_history_rk4["body"] == "Sun"]["kinetic_energy"].values
potential_energy_euler = body_history_euler[body_history_euler["body"] == "Earth"]["potential_energy"].values
potential_energy_rk4 = body_history_rk4[body_history_rk4["body"] == "Earth"]["potential_energy"].values
total_energy_euler = kinectic_energy_euler_earth + kinetic_energy_euler_sun - potential_energy_euler
total_energy_rk4 = kinetic_energy_rk4_earth + kinetic_energy_rk4_sun - potential_energy_rk4
plt.figure(figsize=(10, 5))
plt.plot(body_history_euler[body_history_euler["body"] == "Earth"]["step"], (total_energy_euler-total_energy_euler[0])/total_energy_euler[0], label="Euler Total Energy")
plt.plot(body_history_rk4[body_history_rk4["body"] == "Earth"]["step"], (total_energy_rk4-total_energy_rk4[0])/total_energy_euler[0], label="RK4 Total Energy")
plt.title("Total Energy Comparison: Euler vs RK4")
plt.xlabel("Timestep")
plt.ylabel("Total Energy ")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
save_fig("total_energy_comparison")
plt.show()