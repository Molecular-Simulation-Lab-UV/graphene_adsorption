import numpy as np
import matplotlib.pyplot as plt
from glob import glob

files = sorted(glob("*.dat"))

data = {x[:-4]:np.loadtxt(x) for x in files}

fig, axs = plt.subplots(figsize=(20, 8), nrows=3, ncols=4)

# ax0 = totals
# ax1 = interactions
# ax2 = internal

names0 = ["totcrf", "totlj", "totpot", "totene"]
names1 = ["E_WW", "E_AW", "E_GA", "E_GW"]
names2 = ["int_LJ", "int_CRF", "int_angles", "int_bonds", "int_cross", "int_dihedral", "int_improper", "int_E"]

def myplot(n, x, y):
    X = data[n][:, 0]
    Y = data[n][:, 1]
    E = data[n][:, 2]
    axs[x, y].errorbar(X, Y, E, label=n)

# ax0 = totals
myplot("totcrf", 0, 0)
myplot("totlj", 0, 1)
myplot("totpot", 0, 2)
myplot("totene", 0, 3)

# ax1 interactions
myplot("E_WW", 1, 0)
myplot("E_GW", 1, 1)
myplot("E_AW", 1, 2)
myplot("E_GA", 1, 3)

# ax2 internal
myplot("int_LJ", 2, 0)
myplot("int_CRF", 2, 3)
myplot("int_angles", 2, 1)
#myplot("int_bonds", 2, 1)
#myplot("int_cross", 2, 2)
myplot("int_dihedral", 2, 2)
#myplot("int_improper", 2, 3)
myplot("int_E", 2, 3)

for ax in axs.flatten():
    ax.legend()
    ax.set_xlim(0, 1.5)
    ax.set_xlabel(r"$\xi [nm]$")

fig.tight_layout()
fig.savefig("energy.png", dpi=300)

print("Done!")
