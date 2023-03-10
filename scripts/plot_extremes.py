import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from seaborn import histplot

fig, axs = plt.subplots(figsize=(12, 12), nrows=3, ncols=2)

def readtser(file):
    data = np.loadtxt(file)
    X = data[:, 0] / 1000
    Y1 = data[:, 1]
    Y2 = data[:, 2]
    return X, Y1, Y2

def readener(file):
    data = np.loadtxt(file)
    X = data[:, 0] / 1000
    Y = data[:, 1]
    return X, Y

def mplot(file_tser, file_ener, column, label, color="C0", cut=0):
    X, Y1, Y2 = readtser(file_tser)
    XE, YE = readener(file_ener)

    lr1 = linregress(YE, Y1[cut:])
    lr2 = linregress(YE, Y2[cut:])

#    axs[0, column].hist(Y1, bins=100, density=True, alpha=0.8,
#                        label=f"${label}, r^2 = {lr1.rvalue**2:.2f}$")
    histplot(Y1, bins=100, ax=axs[0, column], binrange=(.2, .6),
             stat="probability", color=color,
             label=f"{label}. $r^2 = {lr1.rvalue**2:.2f}$")
#    axs[1, column].hist(Y2, bins=100, density=True, alpha=0.8,
#                        label=f"${label}, r^2 = {lr2.rvalue**2:.2f}$")
    histplot(Y2, bins=100, ax=axs[1, column], binrange=(.2, .6),
             stat="probability", color=color,
             label=f"{label}, $r^2 = {lr2.rvalue**2:.2f}$")
    #axs[2, column].hist(YE, bins=100, density=True)
    histplot(YE, bins=100, ax=axs[2, column], binrange=(-225, 25),
             stat="probability", color=color, label=f"{label}")


# Plots
mplot("../data/glutamic/r0450_tser_extremes.out",
      "../data/glutamic/r0450_internal_energy.dat",
      0, "Bound", cut=1000)

mplot("../data/glutamic/r1500_tser_extremes.out",
      "../data/glutamic/r1500_internal_energy.dat",
      0, "Unbound", color="C1", cut=1000)

mplot("../data/glutamic/unrestricted_tser_extremes.out",
      "../data/glutamic/unrestricted_internal_energy.dat",
      1, "Unrestricted")

mplot("../data/glutamic/nograph_tser_extremes.out",
      "../data/glutamic/nograph_internal_energy.dat",
      1, "No Graphene", color="C1")

# Aesthetics
axs[0, 0].set_title("Restrained Distance")
axs[0, 0].set_xlabel("Distance (nm)")
#axs[0, 0].set_ylabel("p()")
axs[0, 0].set_ylabel("CD-N Distribution")
axs[0, 0].legend(loc="upper left")
axs[0, 0].set_xlim(0.2, 0.6)
axs[0, 0].set_ylim(0, .075)

axs[1, 0].set_xlabel("Distance (nm)")
#axs[1, 0].set_ylabel("p()")
axs[1, 0].set_ylabel("CD-C Distribution")
axs[1, 0].legend(loc="upper left")
axs[1, 0].set_xlim(0.2, 0.6)
axs[1, 0].set_ylim(0, .1)

axs[2, 0].set_xlabel("Energy (Kj/mol)")
#axs[2, 0].set_ylabel("p()")
axs[2, 0].set_ylabel("Internal Energy Dist.")
axs[2, 0].set_xlim(-225, 25)
axs[2, 0].legend(loc="upper left")
axs[2, 0].set_ylim(0, .05)

axs[0, 1].set_title("Unrestrained")
axs[0, 1].set_xlabel("Distance (nm)")
#axs[0, 1].set_ylabel("p()")
axs[0, 1].set_ylabel("CD-N Distribution")
axs[0, 1].legend(loc="upper left")
axs[0, 1].set_xlim(0.2, 0.6)
axs[0, 1].set_ylim(0, .075)

axs[1, 1].set_xlabel("Distance (nm)")
#axs[1, 1].set_ylabel("p()")
axs[1, 1].set_ylabel("CD-C Distribution")
axs[1, 1].legend(loc="upper left")
axs[1, 1].set_xlim(0.2, 0.6)
axs[1, 1].set_ylim(0, .1)

axs[2, 1].set_xlabel("Energy (Kj/mol)")
#axs[2, 1].set_ylabel("p()")
axs[2, 1].set_ylabel("Internal Energy Dist.")
axs[2, 1].set_xlim(-225, 25)
axs[2, 1].legend(loc="upper left")
axs[2, 1].set_ylim(0, .05)

fig.tight_layout()
fig.savefig("../figures/tser_extremes.pdf", dpi=300)
print("DONE!")
