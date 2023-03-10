import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import linregress
from glob import glob

table_tot = "/home/mbarria/don-elias/AA_GO0US/tables/2022_04_05_parallel_asterisk/Table1.csv"
table_pairs = "/home/mbarria/don-elias/AA_GO0US/tables/2022_04_05_parallel_asterisk/Table7.csv"

output = "../figures"
imagename = "FigS9"

# How to group amino acids
type_dict = {
    "Alanine":"Apolar",
    "Arginine":"Positive",
    "Arginine_H+":"Positive",
    "Asparagine":"Polar",
    "Aspartic_Acid":"Negative",
    "Aspartic_Acid_H+":"Negative",
    "Cysteine_H+":"Polar",
    "Glutamic_Acid":"Negative",
    "Glutamic_Acid*":"Negative",
    "Glutamic_Acid_H+":"Negative",
    "Glutamine":"Polar",
    "Glycine":"Apolar",
    "Histidine_2H+":"Positive",
    "Histidine_H+_ND1":"Positive",
    "Isoleucine":"Apolar",
    "Leucine":"Apolar",
    "Lysine":"Positive",
    "Lysine_H+":"Positive",
    "Methionine":"Apolar",
    "Phenylalanine":"Aromatic",
    "Proline":"Apolar",
    "Serine":"Polar",
    "Threonine":"Polar",
    "Tryptophan":"Aromatic",
    "Tyrosine":"Aromatic",
    "Valine":"Apolar",
    #"Water":"Polar"
}

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

# Better names
name_dict = {
    "Alanine":"Alanine",
    "Arginine":"Arginine",
    "Arginine_H+":"Arginine (+)",
    "Asparagine":"Asparagine",
    "Aspartic_Acid":"Aspartic Acid (-)",
    "Aspartic_Acid_H+":"Aspartic Acid",
    "Cysteine_H+":"Cysteine",
    "Glutamic_Acid":"Glutamic Acid (-)",
    "Glutamic_Acid*":"Glutamic Acid* (-)",
    "Glutamic_Acid_H+":"Glutamic Acid",
    "Glutamine":"Glutamine",
    "Glycine":"Glycine",
    "Histidine_2H+":"Histidine (+)",
    "Histidine_H+_ND1":"Histidine",
    "Isoleucine":"Isoleucine",
    "Leucine":"Leucine",
    "Lysine":"Lysine",
    "Lysine_H+":"Lysine (+)",
    "Methionine":"Methionine",
    "Phenylalanine":"Phenylalanine",
    "Proline":"Proline",
    "Serine":"Serine",
    "Threonine":"Threonine",
    "Tryptophan":"Tryptophan",
    "Tyrosine":"Tyrosine",
    "Valine":"Valine",
    #"Water":#"Water"
}

# 3 Letter Codes
code_dict = {
    "Alanine":"Ala",
    "Arginine":"Arg",
    "Arginine_H+":"Arg+",
#    "Arginine_H+":"Arg",
    "Asparagine":"Asn",
    "Aspartic_Acid":"Asp-",
#    "Aspartic_Acid":"Asp",
    "Aspartic_Acid_H+":"Asp",
    "Cysteine_H+":"Cys",
    "Glutamic_Acid":"Glu-",
    "Glutamic_Acid*":"Glu-*",
#    "Glutamic_Acid":"Glu",
    "Glutamic_Acid_H+":"Glu",
    "Glutamine":"Gln",
    "Glycine":"Gly",
    "Histidine_2H+":"His+",
#    "Histidine_2H+":"His",
    "Histidine_H+_ND1":"His",
    "Isoleucine":"Ile",
    "Leucine":"Leu",
    "Lysine":"Lys",
    "Lysine_H+":"Lys+",
#    "Lysine_H+":"Lys",
    "Methionine":"Met",
    "Phenylalanine":"Phe",
    "Proline":"Pro",
    "Serine":"Ser",
    "Threonine":"Thr",
    "Tryptophan":"Trp",
    "Tyrosine":"Tyr",
    "Valine":"Val",
    #"Water":"$H_2O$"
}

# 1 Letter Codeds
letter_dict = {
    "Alanine":"A",
    "Arginine":"R",
    "Arginine_H+":"R+",
    "Asparagine":"N",
    "Aspartic_Acid":"D-",
    "Aspartic_Acid_H+":"D",
    "Cysteine_H+":"C",
    "Glutamic_Acid":"E-",
    "Glutamic_Acid*":"E-^*",
    "Glutamic_Acid_H+":"E",
    "Glutamine":"Q",
    "Glycine":"G",
    "Histidine_2H+":"H+",
    "Histidine_H+_ND1":"H",
    "Isoleucine":"I",
    "Leucine":"L",
    "Lysine":"K",
    "Lysine_H+":"K+",
    "Methionine":"M",
    "Phenylalanine":"F",
    "Proline":"P",
    "Serine":"S",
    "Threonine":"T",
    "Tryptophan":"W",
    "Tyrosine":"Y",
    "Valine":"V",
    #"Water":"H_2O"
}


# Markers
marker_dict = {
    "Alanine":"o",
    "Arginine":"o",
    "Arginine_H+":"^",
    "Asparagine":"o",
    "Aspartic_Acid":"^",
    "Aspartic_Acid_H+":"o",
    "Cysteine_H+":"o",
    "Glutamic_Acid":"^",
    "Glutamic_Acid*":"^",
    "Glutamic_Acid_H+":"o",
    "Glutamine":"o",
    "Glycine":"o",
    "Histidine_2H+":"^",
    "Histidine_H+_ND1":"o",
    "Isoleucine":"o",
    "Leucine":"o",
    "Lysine":"o",
    "Lysine_H+":"^",
    "Methionine":"o",
    "Phenylalanine":"o",
    "Proline":"o",
    "Serine":"o",
    "Threonine":"o",
    "Tryptophan":"o",
    "Tyrosine":"o",
    "Valine":"o",
    #"Water":"o"
}

# MWeights
weight_dict = {
    "Alanine": 89.09,
    "Arginine": 174.20,
    "Arginine_H+": 175.21,
    "Asparagine": 132.12,
    "Aspartic_Acid": 132.10,
    "Aspartic_Acid_H+": 133.10,
    "Cysteine_H+": 121.15,
    "Glutamic_Acid": 146.12,
    "Glutamic_Acid*": 146.12,
    "Glutamic_Acid_H+": 147.13,
    "Glutamine": 146.15,
    "Glycine": 75.07,
    "Histidine_2H+": 156.16,
    "Histidine_H+_ND1": 155.16,
    "Isoleucine": 131.18,
    "Leucine": 131.18,
    "Lysine": 146.19,
    "Lysine_H+": 147.20,
    "Methionine": 149.21,
    "Phenylalanine": 165.19,
    "Proline": 116.14,
    "Serine": 105.09,
    "Threonine": 119.12,
    "Tryptophan": 204.23,
    "Tyrosine": 181.19,
    "Valine": 117.15,
    #"Water": 18.02,
}

# \Delta SASA
sasa_dict = {
    "Alanine"          : -1.02,
    "Arginine"         : -1.59,
    "Arginine_H+"      : -1.68,
    "Asparagine"       : -1.20,
    "Aspartic_Acid"    : -1.08,
    "Aspartic_Acid_H+" : -1.17,
    "Cysteine_H+"      : -1.13,
    "Glutamic_Acid"    : -1.14,
    "Glutamic_Acid*"    : -1.14,
    "Glutamic_Acid_H+" : -1.35,
    "Glutamine"        : -1.41,
    "Glycine"          : -1.01,
    "Histidine_2H+"    : -1.38,
    "Histidine_H+_ND1" : -1.36,
    "Isoleucine"       : -1.22,
    "Leucine"          : -1.28,
    "Lysine"           : -1.47,
    "Lysine_H+"        : -1.38,
    "Methionine"       : -1.38,
    "Phenylalanine"    : -1.45,
    "Proline"          : -1.13,
    "Serine"           : -1.09,
    "Threonine"        : -1.11,
    "Tryptophan"       : -1.69,
    "Tyrosine"         : -1.55,
    "Valine"           : -1.12,
    #"Water"            : -0.76
}

# Linear regressions
# | - - - - | - - - - |
# | Energy  | Entropy |
# | - - - - | - - - - |
# | AA-Gra  | AA-H2O  |
# | - - - - | - - - - |
# | H2O-Gra | H2O-H2O |
# | - - - - | - - - - |

fig, axs = plt.subplots(figsize=(14, 16.5), ncols=2, nrows=3)

Data_tot = {}

with open(table_tot, 'r') as file:
    file.readline()
    file.readline()
    for line in file:
        aa, dA, errA, dE, errE, TdS, errS = line.split(',')
        values = [float(x) for x in [dA, errA, dE, errE, TdS, errS]]
        Data_tot[aa] = values

Data_pairs = {}

with open(table_pairs, 'r') as file:
    file.readline()
    file.readline()
    for line in file:
        aa, dGA, errGA, dGW, errGW, dWW, errWW, dAW, errAW  = line.split(',')
        values = [float(x) for x in [dGA, errGA, dGW, errGW, dWW, errWW, dAW, errAW]]
        Data_pairs[aa] = values


# Sort amino acids by category
sort_dict_list = sorted(type_dict.items(), key=lambda x:x[1])
sorted_aa = dict(sort_dict_list)
sorted_codes = [code_dict[aa] for aa in sorted_aa.keys()]  # for xlabels

# To set glycine to 0
correction_E = Data_tot["Glycine"][2]
correction_S = Data_tot["Glycine"][4]
correction_GA = Data_pairs["Glycine"][0]
correction_GW = Data_pairs["Glycine"][2]
correction_WW = Data_pairs["Glycine"][4]
correction_AW = Data_pairs["Glycine"][6]

# Store values for the linear regression
weights = []
sasas = []
ener_tot = []
entropy = []
ener_GA = []
ener_GW = []
ener_WW = []
ener_AW = []

for aa in sorted_aa.keys():
    category = type_dict[aa]
    color = colors[category]
    marker = marker_dict[aa]
    # Marker size goes from 2 to 9, proportional to size
    weight = weight_dict[aa]
    sasa = sasa_dict[aa]
    size = 2 + 5 * (weight - min(weight_dict.values()))/(max(weight_dict.values()) - min(weight_dict.values()))
    code = code_dict[aa]
    # Values
    E      = Data_tot[aa][2] - correction_GA
    Eerr   = Data_tot[aa][3]
    S      = Data_tot[aa][4] - correction_GA
    Serr   = Data_tot[aa][5]
    GA     = Data_pairs[aa][0] - correction_GA
    GAerr  = Data_pairs[aa][1]
    GW     = Data_pairs[aa][2] - correction_GW
    GWerr  = Data_pairs[aa][3]
    WW     = Data_pairs[aa][4] - correction_WW
    WWerr  = Data_pairs[aa][5]
    AW     = Data_pairs[aa][6] - correction_AW
    AWerr  = Data_pairs[aa][7]

    # Weights
    letter = letter_dict[aa]
#    size2 = 10 if len(letter) == 1 else 20

    if  (l := len(letter)) == 1:
        size2 = 10
    elif l == 2:
        size2 = 20
    elif l == 4:
        size2 = 30
    else:
        size2 = 0
    condition = len(letter) == 1 # and aa != "Proline" and aa != "Tyrosine"

    if condition:
        weights.append(weight)
        sasas.append(sasa)
        ener_tot.append(E)
        entropy.append(S)
        ener_GA.append(GA)
        ener_GW.append(GW)
        ener_WW.append(WW)
        ener_AW.append(AW)

    # 0, 0 = Energy
    axs[0, 0].plot([sasa], [E], color=color, marker=f"${letter}$", markersize=size2)
    # 0, 1 = Entropy
    axs[0, 1].plot([sasa], [S], color=color, marker=f"${letter}$", markersize=size2)
    # 1, 0 = GA
    axs[1, 0].plot([sasa], [GA], color=color, marker=f"${letter}$", markersize=size2)
    # 1, 1 = AW
    axs[1, 1].plot([sasa], [AW], color=color, marker=f"${letter}$", markersize=size2)
    # 2, 0 = GW
    axs[2, 0].plot([sasa], [GW], color=color, marker=f"${letter}$", markersize=size2)
    # 2, 1 = WW
    axs[2, 1].plot([sasa], [WW], color=color, marker=f"${letter}$", markersize=size2)

# Data for regressions
sasas = np.asarray(sasas)
ener_tot = np.asarray(ener_tot)
entropy  = np.asarray(entropy)
ener_GA  = np.asarray(ener_GA)
ener_GW  = np.asarray(ener_GW)
ener_WW  = np.asarray(ener_WW)
ener_AW  = np.asarray(ener_AW)
# Regressions
lr_ener_tot = linregress(sasas, ener_tot)
lr_entropy = linregress(sasas, entropy)
lr_ener_GA = linregress(sasas, ener_GA)
lr_ener_GW = linregress(sasas, ener_GW)
lr_ener_WW = linregress(sasas, ener_WW)
lr_ener_AW = linregress(sasas, ener_AW)

lr_list = [lr_ener_tot, lr_entropy, lr_ener_GA, lr_ener_AW, lr_ener_GW, lr_ener_WW]

# y labels
axs[0, 0].set_ylabel("$\Delta E^{ads}_{pot}~[kJ \cdot mol^{-1}]$", fontsize=13)
axs[0, 1].set_ylabel("$T \Delta S^{ads}~[kJ \cdot mol^{-1}]$", fontsize=13)
axs[1, 0].set_ylabel("$\Delta E^{ads}_{AA-Graph}~[kJ \cdot mol^{-1}]$", fontsize=13)
axs[1, 1].set_ylabel("$\Delta E^{ads}_{AA-H_2O}~[kJ \cdot mol^{-1}]$", fontsize=13)
axs[2, 0].set_ylabel("$\Delta E^{ads}_{H_2O-Graph}~[kJ \cdot mol^{-1}]$", fontsize=13)
axs[2, 1].set_ylabel("$\Delta E^{ads}_{H_2O-H_2O}~[kJ \cdot mol^{-1}]$", fontsize=13)

# Custom handles por labels
custom_handles = [Line2D([0], [0], color=colors[key], marker="$X$", markersize=5, fillstyle='none', linewidth=0, label=key) for key in colors.keys()]

w = np.linspace(-1.8, -0.9, 50)
# For all axes
for ind, ax in enumerate(axs.flatten()):
    lr = lr_list[ind]
    # use lr to decide legend locations
    if lr.slope >= 0:
        loc_leg = "lower right"
        loc_r   = "upper left"
    elif ind == 3:
        loc_leg = "upper right"
        loc_r   = "upper left"
    else:
        loc_leg = "lower left"
        loc_r   = "upper right"

    ax.set_xlabel(r"$\Delta$SASA$^{ads}$ $[nm^{2}]$", fontsize=14)
    ax.set_xlim(-1.8, -0.9)
    leg_colors = ax.legend(handles=custom_handles, loc=loc_leg, frameon=False, fontsize=12)

    ax.plot(w, w*lr.slope + lr.intercept, color="k",
             label=f"$r^2_{{neutral}}={lr.rvalue**2:.2f}$",
             linewidth=1, linestyle="--")
    leg_lr = ax.legend(loc=loc_r, frameon=False, fontsize=14)

    ax.add_artist(leg_lr)
    ax.add_artist(leg_colors)

    ax.text(-0.1, 1.05, "abcdef"[ind], transform=ax.transAxes,
            size=24, weight='bold')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)

fig.tight_layout()
fig.savefig(f"{output}/{imagename}.pdf", dpi=600)

print("DONE!")
