import matplotlib.pyplot as plt

# Energy levels normalized to ground state (E - Egs)
Z = 82
A = 208
N = A- Z

energy_levels = [
      -0.06,
      -0.24,
      -0.29,
      -0.48,
      -0.73,
      -0.82,
      -0.85,
      -1.19,
      -1.24,
      -1.88,
      -2.25,
      -2.34,
      -2.36,
      -2.43,
      -2.48,
      -4.07,
      -5.96,
      -6.70,
      -7.15,
      -7.33,
      -8.13,
      -8.71,
      -9.40,
      -9.52,
      -9.79,
     -10.28,
     -10.41,
     -10.57,
     -11.01,
     -11.02,
     -13.39,
     -14.15,
     -14.40,
     -14.50,
     -14.78,
     -16.39,
     -17.75,
     -18.02,
     -18.31,
     -18.84,
     -18.89,
     -19.22,
     -19.56,
     -21.77,
     -22.46,
     -22.54,
     -23.29,
     -23.78,
     -23.93,
     -24.25,
     -24.28,
     -24.64,
     -26.24,
     -27.75,
     -28.05,
     -28.19,
     -28.73,
     -28.98,
     -29.05,
     -30.68,
     -30.94,
     -31.84,
     -32.19,
     -32.20,
     -32.29,
     -32.98,
]

print(energy_levels)
# Plot setup
fig, ax = plt.subplots(figsize=(6, 10))
ax.set_title("Energy Level Scheme")
ax.set_ylabel("Energy (MeV above ground state)")
protons = 0
# Draw horizontal lines for each level
for i, energy in enumerate(reversed(energy_levels)):
    protons = protons + 2
    color = 'black' if protons <= Z else 'red'
    ax.hlines(energy, xmin=-0.4, xmax=0.4, color=color)

    ax.text(0.45, energy, f"{protons}", ha='center', fontsize=6)
    # Optionally label the energy value

ax.set_xlim(-1, 1)
plt.tight_layout()
plt.show()