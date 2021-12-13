import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import linregress
# Path to csv file
table = "../data/tables/Free_Energy_Decomposition.csv"
#table = "/home/mbarria/donelias/AA_GO0US/tables/2021_08_28/Table1.csv"
#output
output = "../figures"
imagename = "Differences_by_names_and_weight"
# separate figures?
SEPARATE = False
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
#    "Water":"Polar"
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
#    "Water":"Water"
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
#    "Water":"H2O"
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
#    "Water":"H2O"
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
#    "Water":"o"
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
#    "Water": 18.02,
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

# Graph
if SEPARATE:
    fig1, ax1 = plt.subplots(figsize=(6, 5), ncols=1, nrows=1)
    fig2, ax2 = plt.subplots(figsize=(6, 5), ncols=1, nrows=1)
    fig3, ax3 = plt.subplots(figsize=(6, 5), ncols=1, nrows=1)
    fig4, ax4 = plt.subplots(figsize=(6, 5), ncols=1, nrows=1)
    fig5, ax5 = plt.subplots(figsize=(6, 5), ncols=1, nrows=1)
    fig6, ax6 = plt.subplots(figsize=(6, 5), ncols=1, nrows=1)
    axs = np.array([[ax1, ax2, ax3] , [ax4, ax5, ax6]])
else:
    fig, axs = plt.subplots(figsize=(18, 10), ncols=3, nrows=2)
# Set 0 to glycine
correction_A = Data["Glycine"][0]
correction_E = Data["Glycine"][2]
correction_S = Data["Glycine"][4]
# Store values with respect to weight for a linear regression
weights = []
free_energies = []
energies = []
entropies = []
for aa in sorted_aa.keys():
    category = type_dict[aa]
    color = colors[category]
    marker = marker_dict[aa]
    # Marker size goes from 2 to 9, proportional to size
    weight = weight_dict[aa]
    size = 2 + 5 * (weight - min(weight_dict.values()))/(max(weight_dict.values()) - min(weight_dict.values()))
    code = code_dict[aa]
    values = Data[aa]
    # Weights
    letter = letter_dict[aa]
    size2 = 10 if len(letter) == 1 else 20
    condition = len(letter) == 1 # and aa != "Proline" and aa != "Tyrosine"
    #condition = True
    # Free energy
    axs[0, 0].bar([code], [values[0] - correction_A], yerr=[values[1]], color=color, align='center', alpha=0.8, ecolor='black', capsize=1)
    # Bottom by weight
    axs[1, 0].plot([weight], [values[0] - correction_A], color=color, marker=f"${letter}$", markersize=size2)
#    axs[1, 0].plot([weight], [values[0] - correction_A], color='k', marker=f"o", markersize=.5, linewidth=0)
    if condition:
        weights.append(weight)
        free_energies.append(values[0] - correction_A)
    # Energy
    axs[0, 1].bar([code], [values[2] - correction_E], yerr=[values[3]], color=color, align='center', alpha=0.8, ecolor='black', capsize=1)
    # Bottom by w
    axs[1, 1].plot([weight], [values[2] - correction_E], color=color, marker=f"${letter}$", markersize=size2, fillstyle='none')
    if condition:
        energies.append(values[2] - correction_E)
    # Entropy
    # axs[0, 2].errorbar([code], [values[4] - correction_S], [values[5]], color=color, marker=marker, markersize=size, fillstyle='none')
    axs[0, 2].bar([code], [values[4] - correction_S], yerr=[values[5]], color=color, align='center', alpha=0.8, ecolor='black', capsize=1)
    # Bottom by w
    axs[1, 2].plot([weight], [values[4] - correction_S], color=color, marker=f"${letter}$", markersize=size2, fillstyle='none')
    if condition:
        entropies.append(values[4] - correction_S)
    # comment because my editor forces the below "else" to go after the above "if"

# y axis changes
axs[0, 0].set_ylabel("$\Delta A_{ads}~[kJ \cdot mol^{-1}]$", fontsize=13)
axs[0, 0].set_title("Adsorption Free Energy", fontsize=15)
axs[0, 0].set_yticklabels(["-20", "-15", "-10", "-5", "0"], fontsize=13)
axs[0, 0].set_yticks([-20, -15, -10, -5, 0])

axs[1, 0].set_ylabel("$\Delta A_{ads}~[kJ \cdot mol^{-1}]$", fontsize=13)
axs[1, 0].set_yticklabels(["-20", "-15", "-10", "-5", "0"], fontsize=13)
axs[1, 0].set_yticks([-20, -15, -10, -5, 0])

axs[0, 1].set_ylabel("$\Delta E_{ads}~[kJ \cdot mol^{-1}]$", fontsize=13)
axs[0, 1].set_title("Adsorption Energy", fontsize=15)
axs[0, 1].set_yticklabels(["-40", "-30", "-20", "-10", "0", "10"], fontsize=13)
axs[0, 1].set_yticks([-40, -30, -20, -10, 0, 10])

axs[1, 1].set_ylabel("$\Delta E_{ads}~[kJ \cdot mol^{-1}]$", fontsize=13)
axs[1, 1].set_yticklabels(["-40", "-30", "-20", "-10", "0", "10"], fontsize=13)
axs[1, 1].set_yticks([-40, -30, -20, -10, 0, 10])

axs[0, 2].set_ylabel("$T \Delta S_{ads}~[kJ \cdot mol^{-1}]$", fontsize=13)
axs[0, 2].set_title("Adsorption Entropy", fontsize=15)
axs[0, 2].set_yticklabels(["-25", "-20", "-15", "-10", "-5", "0", "5", "10"], fontsize=13)
axs[0, 2].set_yticks([-25, -20, -15, -10, -5, 0, 5, 10])

axs[1, 2].set_ylabel("$T \Delta S_{ads}~[kJ \cdot mol^{-1}]$", fontsize=13)
axs[1, 2].set_yticklabels(["-25", "-20", "-15", "-10", "-5", "0", "5", "10"], fontsize=13)
axs[1, 2].set_yticks([-25, -20, -15, -10, -5, 0, 5, 10])
# Formats
#custom_handles1 = [Line2D([0], [0], color=colors[key], marker="o", markersize=5, fillstyle='none', linewidth=0, label=key) for key in colors.keys()]
# For top row
custom_handles1 = [Line2D([0], [0], color=colors[key], linewidth=3, label=key) for key in colors.keys()]
custom_handles2 = [
    Line2D([0], [0], color='k', marker="o", markersize=5, fillstyle='none', linewidth=0, label="Neutral"),
    Line2D([0], [0], color='k', marker="^", markersize=5, fillstyle='none', linewidth=0, label="Charged")
]
# Bottom Row
custom_handles3 = [Line2D([0], [0], color=colors[key], marker="$X$", markersize=5, fillstyle='none', linewidth=0, label=key) for key in colors.keys()]
# Changes to all top axes
for ax in axs[0, :]:
    ax.tick_params(labelrotation=60)
    ax.set_xlabel("Amino Acid", fontsize=14)
    leg1 = ax.legend(handles=custom_handles1, loc="lower left", frameon=False, fontsize=12)
    #ax.legend(handles=custom_handles2, loc="upper right")
    #ax.add_artist(leg1)
    ax.set_xticklabels(sorted_codes, fontsize=10)


# Linear regression in free energy vs weight
weights = np.asarray(weights)
free_energies = np.asarray(free_energies)
lr = linregress(weights, free_energies)

w = np.linspace(min(weight_dict.values()), max(weight_dict.values()))
axs[1, 0].plot(w, w*lr.slope + lr.intercept, color="k",
               label=f"{lr.slope:.1f}m + {lr.intercept:.1f}, $r^2$={lr.rvalue**2:.2f}",
               linewidth=1, linestyle="--")
leg1 = axs[1, 0].legend(loc="upper right", frameon=False, fontsize=10)

# Test other regressions
lr2 = linregress(weights, energies)
axs[1, 1].plot(w, w*lr2.slope + lr2.intercept, color="k",
               label=f"{lr2.slope:.1f}m + {lr2.intercept:.1f}, $r^2$={lr2.rvalue**2:.2f}",
               linewidth=1, linestyle="--")
leg2 = axs[1, 1].legend(loc="upper right", frameon=False, fontsize=10)

lr3 = linregress(weights, entropies)
axs[1, 2].plot(w, w*lr3.slope + lr3.intercept, color="k",
               label=f"{lr3.slope:.1f}m + {lr3.intercept:.1f}, $r^2$={lr3.rvalue**2:.2f}",
               linewidth=1, linestyle="--")
leg3 = axs[1, 2].legend(loc="upper right", frameon=False, fontsize=10)

# Changes to bottom axes
for ax in axs[1, :]:
    #ax.tick_params(labelrotation=90)
    ax.set_xlabel("Molecular Weight $[g \cdot mol^{-1}]$", fontsize=14)
    ax.legend(handles=custom_handles3, loc="lower left", frameon=False, fontsize=12)
    ax.set_xticklabels(["80", "100", "120", "140", "160", "180", "200"], fontsize=13)
    ax.set_xticks([80, 100, 120, 140, 160, 180, 200])
    #ax.legend(handles=custom_handles2, loc="upper right")

axs[1, 0].add_artist(leg1)
axs[1, 1].add_artist(leg2)
axs[1, 2].add_artist(leg3)

# Changes to all axes
for ind, ax in enumerate(axs.T.flatten()):
    ax.text(-0.1, 1.05, "ABCDEF"[ind], transform=ax.transAxes,
            size=20, weight='bold')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)

if SEPARATE:
    for ind, fig in enumerate([fig1, fig2, fig3, fig4, fig5, fig6]):
        fig.tight_layout()
        fig.savefig(f"{output}/{imagename}_{ind+1}.png", dpi=600)
else:
    fig.tight_layout()
    fig.savefig(f"{output}/{imagename}.png", dpi=600)
