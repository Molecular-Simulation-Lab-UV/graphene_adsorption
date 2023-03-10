import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import linregress
# Path to csv file
#table = "../data/tables/Free_Energy_Decomposition_SimpleDiff.csv"
table = "/home/mbarria/don-elias1/AA_GO25US/tables/2022_02_15_stationary_bootstrap_MyWHAM/Table1.csv"
#output
output = "../figures"
imagename = "Fig7"

# How to group amino acids
type_dict = {
#    "Alanine":"Apolar",
#    "Arginine":"Positive",
    "Arginine_H+":"Positive",
    "Asparagine":"Polar",
#    "Aspartic_Acid":"Negative",
#    "Aspartic_Acid_H+":"Negative",
#    "Cysteine_H+":"Polar",
    "Glutamic_Acid":"Negative",
#    "Glutamic_Acid_H+":"Negative",
#    "Glutamine":"Polar",
    "Glycine":"Apolar",
#    "Histidine_2H+":"Positive",
#    "Histidine_H+_ND1":"Positive",
    "Isoleucine":"Apolar",
#    "Leucine":"Apolar",
#    "Lysine":"Positive",
#    "Lysine_H+":"Positive",
#    "Methionine":"Apolar",
    "Phenylalanine":"Aromatic",
#    "Proline":"Apolar",
#    "Serine":"Polar",
#    "Threonine":"Polar",
#    "Tryptophan":"Aromatic",
#    "Tyrosine":"Aromatic",
#    "Valine":"Apolar",
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

Data = {}

with open(table, 'r') as file:
    file.readline()
    file.readline()
    for line in file:
        aa, dA, errA, dE, errE, TdS, errS = line.split(',')
        values = [float(x) for x in [dA, errA, dE, errE, TdS, errS]]
        Data[aa] = values

# Sort amino acids by category
sort_dict_list = sorted(type_dict.items(), key=lambda x:x[1])
#sort_dict_list = sorted(weight_dict.items(), key=lambda x:x[1])
sorted_aa = dict(sort_dict_list)
sorted_codes = [code_dict[aa] for aa in sorted_aa.keys()]  # for xlabels

fig, axs = plt.subplots(figsize=(7, 12), ncols=1, nrows=2)

correction_A = Data["Glycine"][0]

# Store values for the linear regression
weights = []
sasas = []
free_energies = []

for aa in sorted_aa.keys():
    category = type_dict[aa]
    color = colors[category]
    marker = marker_dict[aa]
    # Marker size goes from 2 to 9, proportional to size
    weight = weight_dict[aa]
    sasa = sasa_dict[aa]
    size = 2 + 5 * (weight - min(weight_dict.values()))/(max(weight_dict.values()) - min(weight_dict.values()))
    code = code_dict[aa]
    values = Data[aa]
    # Weights
    letter = letter_dict[aa]
    size2 = 10 if len(letter) == 1 else 20
    condition = len(letter) == 1 # and aa != "Proline" and aa != "Tyrosine"

    if condition:
        weights.append(weight)
        free_energies.append(values[0] - correction_A)
        sasas.append(sasa)

    # Free energy
    # Top by amino_acid
    axs[0].bar([code], [values[0] - correction_A], yerr=[values[1]], color=color, align='center', alpha=0.8, ecolor='black', capsize=1)
    # Bottom by SASA
    axs[1].plot([sasa], [values[0] - correction_A], color=color, marker=f"${letter}$", markersize=size2)

# y axis changes
axs[0].set_ylabel("$\Delta A_{ads}~[kJ \cdot mol^{-1}]$", fontsize=13)
axs[0].set_title("Adsorption Free Energy", fontsize=15)
#axs[0].set_yticklabels(["-25", "-20", "-15", "-10", "-5", "0", "5"], fontsize=13)
#axs[0].set_yticks([-25, -20, -15, -10, -5, 0, 5])
axs[0].set_ylim(-3, 7)

axs[1].set_ylabel("$\Delta A_{ads}~[kJ \cdot mol^{-1}]$", fontsize=13)
#axs[1].set_yticklabels(["-25", "-20", "-15", "-10", "-5", "0", "5"], fontsize=13)
#axs[1].set_yticks([-25, -20, -15, -10, -5, 0, 5])
axs[1].set_ylim(-3, 7)


# For top row
custom_handles1 = [Line2D([0], [0], color=colors[key], linewidth=3, label=key) for key in colors.keys()]
custom_handles2 = [
    Line2D([0], [0], color='k', marker="o", markersize=5, fillstyle='none', linewidth=0, label="Neutral"),
    Line2D([0], [0], color='k', marker="^", markersize=5, fillstyle='none', linewidth=0, label="Charged")
]
# Middle and Bottom rows
custom_handles3 = [Line2D([0], [0], color=colors[key], marker="$X$", markersize=5, fillstyle='none', linewidth=0, label=key) for key in colors.keys()]
custom_handles4 = [Line2D([0], [0], color=colors[key], marker="$X$", markersize=5, fillstyle='none', linewidth=0, label=key) for key in colors.keys()]
# Changes to top ax

axs[0].tick_params(labelrotation=60)
axs[0].set_xlabel("Amino Acid", fontsize=14)
leg1 = axs[0].legend(handles=custom_handles1, loc="lower left", frameon=False, fontsize=12)
axs[0].add_artist(leg1)
axs[0].set_xticklabels(sorted_codes, fontsize=10)

# Bottom ax

axs[1].set_xlabel(r"$\Delta$SASA$^{ads}$ $[nm^{2}]$", fontsize=14)
#axs[2].set_xticklabels(["-1.8", "-1.6", "-1.4", "-1.2", "-1.0", "-0.8"], fontsize=13)
#axs[2].set_xticks([-1.8, -1.6, -1.4, -1.2, -1.0, -0.8])
axs[1].set_xlim(-1.8, -0.9)
leg3_colors = axs[1].legend(handles=custom_handles4, loc="upper left", frameon=False, fontsize=12)

# Linear regression in free energy vs weight
weights = np.asarray(weights)
free_energies = np.asarray(free_energies)
sasas = np.asarray(sasas)

# bottom ax

lr_sasa = linregress(sasas, free_energies)

w = np.linspace(-1.8, -0.7, 50)
axs[1].plot(w, w*lr_sasa.slope + lr_sasa.intercept, color="k",
            #label=f"{lr_sasa.slope:.1f}m + {lr_sasa.intercept:.1f}, $r^2$={lr_sasa.rvalue**2:.2f}",
            label=f"$r^2$={lr_sasa.rvalue**2:.2f}",
            linewidth=1, linestyle="--")
leg3 = axs[1].legend(loc="upper right", frameon=False, fontsize=14)
axs[1].add_artist(leg3)
axs[1].add_artist(leg3_colors)

# Changes to all axes
for ind, ax in enumerate(axs.T.flatten()):
    ax.text(-0.1, 1.05, "abcdef"[ind], transform=ax.transAxes,
            size=24, weight='bold')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)

fig.tight_layout()
fig.savefig(f"{output}/{imagename}.pdf", dpi=600)

print("DONE!")
