#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import linregress
from glob import glob


table = "../data/tables/Comparison.csv"
output = "../figures"
imagename = "FigS10"

# How to group amino acids
type_dict = {
    "Alanine":"Apolar",
    # "Arginine":"Positive",
    "Arginine_H+":"Positive",
    "Asparagine":"PolarE",
    "Aspartic_Acid":"PomNegative",
    # "Aspartic_Acid_H+":"Negative",
    "Cysteine_H+":"PolarD",
    "Glutamic_Acid":"PomNegative",
    # "Glutamic_Acid*":"Negative",
    # "Glutamic_Acid_H+":"Negative",
    "Glutamine":"PolarF",
    "Glycine":"AApolar",
    # "Histidine_2H+":"Positive",
    "Histidine_H+_ND1":"ZPositive",
    "Isoleucine":"Apolar",
    "Leucine":"ACpolar",
    # "Lysine":"Positive",
    "Lysine_H+":"Positive",
    "Methionine":"Apolar",
    "Phenylalanine":"AromaticB",
    "Proline":"PolarA",
    "Serine":"PolarB",
    "Threonine":"PolarC",
    "Tryptophan":"AromaticA",
    "Tyrosine":"AromaticC",
    "Valine":"ABpolar",
    #"Water":"Polar"
}

# https://personal.sron.nl/~pault/#sec:qualitative
colors = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']

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

# Comparisons
# | - - - - | - - - - - |
# | Abs     |  -Glycine |
# | - - - - | - - - - - |
#

fig, axs = plt.subplots(figsize=(12, 24), ncols=1, nrows=3)

Data = {}

labels = []

with open(table, 'r') as file:
    file.readline()
    headers = file.readline().split(",")
    for head in headers[1::2]:
        labels.append(head)
    for line in file:
        aa, dA, errA, Das1, err1, Das2, err2, Das3, err3, Das4, err4, Hugh, errH = line.split(',')
        values = [float(x) for x in [dA, errA, Das1, err1, Das2, err2, Das3, err3, Das4, err4, Hugh, errH]]
        Data[aa] = values

labels[1] += "$^{[50]}$"
labels[2] += "$^{[50]}$"
labels[3] += "$^{[50]}$"
labels[4] += "$^{[50]}$"
labels[-1] += "$^{[48]}$"
# Sort amino acids by category
sort_dict_list = sorted(type_dict.items(), key=lambda x:x[1])
sorted_aa = dict(sort_dict_list)
sorted_codes = [code_dict[aa] for aa in sorted_aa.keys()]  # for xlabels

# To set Gly as 0
correction_A  = Data["Glycine"][0]
correction_D1 = Data["Glycine"][2]
correction_D2 = Data["Glycine"][4]
correction_D3 = Data["Glycine"][6]
correction_D4 = Data["Glycine"][8]
correction_H  = Data["Glycine"][10]

kb = 0.00813446
corrections_Gly = [correction_A, correction_D1, correction_D2, correction_D3, correction_D4, correction_H]
corrections_Std = [
    -kb * 298.15 * np.log(5.056*5.004*.400/1.661),
    -kb * 300 * np.log(4.26*4.18*1.1/1.661),
    -kb * 300 * np.log(4.26*4.18*1.1/1.661),
    -kb * 300 * np.log(4.26*4.18*1.1/1.661),
    -kb * 300 * np.log(4.26*4.18*1.1/1.661),
    -kb * 300 * np.log(63.9*59.8*42.5/1661)
]

print(corrections_Std)

for ind, source in enumerate(labels):
    color = colors[ind]

    vals = np.array([Data[aa][ind * 2] for aa in sorted_aa.keys()])
    errs = np.array([Data[aa][ind * 2 + 1] for aa in sorted_aa.keys()])

    axs[0].errorbar(sorted_codes, vals, errs, color=color, label=source)
    axs[1].errorbar(sorted_codes, vals + corrections_Std[ind], errs, color=color, label=source)
    axs[2].errorbar(sorted_codes, vals - corrections_Gly[ind], errs, color=color, label=source)

# y axis changes
axs[0].set_ylabel("$\Delta A^{ads}~[kJ \cdot mol^{-1}]$", fontsize=13)
axs[0].tick_params(labelrotation=60)
axs[0].set_xlabel("Amino Acid", fontsize=14)
axs[0].set_xticklabels(sorted_codes, fontsize=10)

axs[1].set_ylabel("$\Delta A^{ads}~[kJ \cdot mol^{-1}]$", fontsize=13)
axs[1].tick_params(labelrotation=60)
axs[1].set_xlabel("Amino Acid", fontsize=14)
axs[1].set_xticklabels(sorted_codes, fontsize=10)

axs[2].set_ylabel("$\Delta A^{ads}~[kJ \cdot mol^{-1}]$", fontsize=13)
axs[2].tick_params(labelrotation=60)
axs[2].set_xlabel("Amino Acid", fontsize=14)
axs[2].set_xticklabels(sorted_codes, fontsize=10)

# Changes to all axes
for ind, ax in enumerate(axs.T.flatten()):
    ax.text(-0.1, 1.05, "abcdef"[ind], transform=ax.transAxes,
            size=24, weight='bold')
    ax.legend()
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)

fig.tight_layout()
fig.savefig(f"{output}/{imagename}.pdf", dpi=600)

print("DONE!")
