import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import linregress
from glob import glob

output = "../figures"
imagename = "FigS3"
ylabel = r"PMF$(\xi) (kJ \cdot mol^{-1})$"

files = sorted(glob("../data/pmfs/*.dat"))
files = [file for file in files if "25o" in file]

# Use this to fill in integrated parts
boundaries = (0.9, 1.4)

# https://personal.sron.nl/~pault/#sec:qualitative
colors = {
#    "Apolar":"C0",
    "Apolar":"#4477aa",
#    "Aromatic":"C1",
    "Aromatic":"#66ccee",
#    "Negative":"C2",
    "Negative":"#228833",
#    "Polar":"C3",
    "Polar":"#ccbb44",
#    "Positive":"C4"
    "Positive":"#ee6677"
}



fig, axs = plt.subplots(figsize=(12, 4), ncols=2, nrows=1)

ax1 = axs[0]
ax2 = axs[1]

# common
c = colors["Apolar"]
ls ="solid"
# Ax1: Helix
ax1.set_title("5-Alanine (Helix)")

data = np.loadtxt("../data/5ala/helix.dat")

ax1.errorbar(data[:, 0], data[:, 1], data[:, 2], color=c,
             linewidth=2, elinewidth=2, linestyle=ls)

# Ax2: Stretched
ax2.set_title("5-Alanine (Stretched)")

data = np.loadtxt("../data/5ala/stretched.dat")

ax2.errorbar(data[:, 0], data[:, 1], data[:, 2], color=c,
             linewidth=2, elinewidth=2, linestyle=ls)
# Common:
for ax in axs:
    ax.axvline(boundaries[0], color="#BBCCEE", linewidth=2)
    ax.axvline(boundaries[1], color="#FFCCCC", linewidth=2)

    ax.set_xlim(0, 2)
    #ax.set_xticks([0, .5, 1., 1.5])
    #ax.set_xticklabels([0, 0.5, 1.0, 1.5], fontsize=12.5)
    ax.set_xlabel(r"$\xi (nm)$", fontsize=12.5)
    ax.set_ylim(0, 50)
    #ax.set_yticks([0, 10, 20, 30])
    #ax.set_yticklabels([0, 10, 20, 30], fontsize=12.5)
    ax.set_ylabel(ylabel, fontsize=11)

    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(1.5)


fig.tight_layout()
fig.savefig(f"{output}/{imagename}.pdf", dpi=300)
print("DONE!")
