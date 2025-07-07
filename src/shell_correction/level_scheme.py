import matplotlib.pyplot as plt

# Energy levels normalized to ground state (E - Egs)
Z = 82
A = 208
N = A-Z

def readlevelfile(file):
    p = []
    n = []
    with open(file, 'r') as f:
        hdr = f.readline().split(',')
        Z = int(hdr[0])
        A = int(hdr[1])
        N = A-Z
        
        for line in f:
            split = line.split(',')
            p.append(float(split[0]))
            n.append(float(split[1]))
    p.sort()
    n.sort()
    return p,n
# Plot setup
e_p, e_n = readlevelfile('data/out/levels.dat')
fig, ax = plt.subplots(figsize=(6, 10))
ax.set_title("Energy Level Scheme")
ax.set_ylabel("Energy (MeV)")
# Draw proton levels (left)
protons = 0
for energy in (e_p):
    if(energy > 0): continue
    protons += 2
    color = 'red' if protons <= Z else 'grey'
    ax.hlines(energy, xmin=-0.8, xmax=-0.2, color=color)
    t = ax.text(-0.85, energy, f"{protons}", ha='right', va='center', fontsize=6)
    t.set_bbox(dict(facecolor='white',alpha=1,edgecolor='white'))

# Draw neutron levels (right)
neutrons = 0
for energy in (e_n):
    if(energy > 0): continue
    neutrons += 2
    color = 'blue' if neutrons <= N else 'grey'
    ax.hlines(energy, xmin=0.2, xmax=0.8, color=color)
    t = ax.text(0.85, energy, f"{neutrons}", ha='left', va='center', fontsize=6)
    t.set_bbox(dict(facecolor='white',alpha=1,edgecolor='white'))

# Label sides
ax.text(-0.5, 2, "Protons", ha='center', fontsize=10, color='red')
ax.text(0.5, 2, "Neutrons", ha='center', fontsize=10, color='blue')
#ax.text(0,2, r"$^{240}$Pu", fontsize=12)
# Format x-axis
ax.set_xticks([])
ax.set_xlim(-1, 1)
ax.set_ylim(min(min(e_p),min(e_n))-5,0)

plt.tight_layout()
plt.savefig('data/out/levelscheme.png')
plt.show()
