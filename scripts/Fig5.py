import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import linregress
# Path to csv file
#table = "../data/tables/Free_Energy_Decomposition.csv"
#table_intra = "/home/mbarria/don-elias/AA_GO0US/tables/2022_02_07_simple_diff/Table3.csv"
table_intra = "/home/mbarria/don-elias/AA_GO0US/tables/2022_04_05_parallel_asterisk/Table3.csv"
#table_tot = "../data/tables/Free_Energy_Decomposition_SimpleDiff.csv"
#table_ener = "/home/mbarria/don-elias/AA_GO0US/tables/2022_02_07_simple_diff/Table7.csv"
table_ener = "/home/mbarria/don-elias/AA_GO0US/tables/2022_04_05_parallel_asterisk/Table7.csv"
#table_ener = "../data/tables/Energy_Decomposition_SimpleDiff.csv"
#output
output = "../figures"
imagename = "Fig5"

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
    "Glutamic_Acid*":"E-*",
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

Data_Int = {}
Data_E = {}

with open(table_ener, 'r') as file:
    file.readline()
    file.readline()
    for line in file:
        aa, dGA, errGA, dGW, errGW, dWW, errWW, dAW, errAW  = line.split(',')
        values = [float(x) for x in [dGA, errGA, dGW, errGW, dWW, errWW, dAW, errAW]]
        Data_E[aa] = values

with open(table_intra, 'r') as file:
    file.readline()
    file.readline()
    for line in file:
        aa, dInt, errInt  = line.split(',')[:3]
        values = [float(x) for x in [dInt, errInt]]
        Data_Int[aa] = values


# Sort amino acids by category
sort_dict_list = sorted(type_dict.items(), key=lambda x:x[1])
#sort_dict_list = sorted(weight_dict.items(), key=lambda x:x[1])
sorted_aa = dict(sort_dict_list)
sorted_codes = [code_dict[aa] for aa in sorted_aa.keys()]  # for xlabels

# One column
# Two rows: 1) AA-graphene, 2) AA-water
fig, axs = plt.subplots(figsize=(8, 12), ncols=1, nrows=2)

correction_GA = Data_E["Glycine"][0]
correction_GW = Data_E["Glycine"][2]
correction_WW = Data_E["Glycine"][4]
correction_AW = Data_E["Glycine"][6]
correction_Int = Data_Int["Glycine"][0]


# Store values for the linear regression
weights = []
sasas = []
energies_GA  = []
energies_GW  = []
energies_WW  = []
energies_AW  = []
energies_Int = []

## Collect data
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
    GA     = Data_E[aa][0] - correction_GA
    GAerr  = Data_E[aa][1]
    GW     = Data_E[aa][2] - correction_GW
    GWerr  = Data_E[aa][3]
    WW     = Data_E[aa][4] - correction_WW
    WWerr  = Data_E[aa][5]
    AW     = Data_E[aa][6] - correction_AW
    AWerr  = Data_E[aa][7]
    Int    = Data_Int[aa][0] - correction_Int
    Interr = Data_Int[aa][1]
    # Weights
    letter = letter_dict[aa]
    size2 = 10 if len(letter) == 1 else 20
    condition = len(letter) == 1 # and aa != "Proline" and aa != "Tyrosine"

    if condition:
        weights.append(weight)
        sasas.append(sasa)
        energies_GA.append(GA)
        energies_GW.append(GW)
        energies_WW.append(WW)
        energies_AW.append(AW)
        energies_Int.append(Int)
    ## Amino acid interactions
    # Top: AA-graph
    axs[0].bar([code], [GA], yerr=[GAerr], color=color, align='center', alpha=0.8, ecolor='black', capsize=1)
    # Bottom: AA-water
    axs[1].bar([code], [AW], yerr=[AWerr], color=color, align='center', alpha=0.8, ecolor='black', capsize=1)

## Linear Regressions

weights = np.asarray(weights)
sasas = np.asarray(sasas)
energies_GA  = np.asarray(energies_GA)
energies_GW  = np.asarray(energies_GW)
energies_WW  = np.asarray(energies_WW)
energies_AW  = np.asarray(energies_AW)
energies_Int = np.asarray(energies_Int)

w = np.linspace(-1.8, -0.7)

# Graphene--Amino-acid
lr = linregress(sasas, energies_GA)
axs[0].text("Ser", 80, f"$r_{{\Delta SASA, neutral}}^2={lr.rvalue**2:.2f}$", size=14)
# Graphene--Water
lr = linregress(sasas, energies_AW)
axs[1].text("Ser", 80, f"$r_{{\Delta SASA, neutral}}^2={lr.rvalue**2:.2f}$", size=14)

## Changes to all axes
for ind, ax in enumerate(axs.T.flatten()):
    ax.text(-0.1, 1.05, "abcdefghijk"[ind], transform=ax.transAxes,
            size=20, weight='bold')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)

## Changes to all upper axes:
custom_handles1 = [Line2D([0], [0], color=colors[key], linewidth=3, label=key) for key in colors.keys()]

for ax in axs:
    # y
    ax.set_ylim(-55, 95)
    # x
    ax.tick_params(axis='x', labelrotation=60)
    ax.set_xlabel("Amino Acid", fontsize=14)
    leg1 = ax.legend(handles=custom_handles1, loc="upper left", frameon=False, fontsize=12)
    #ax.add_artist(leg1)
    ax.set_xticklabels(sorted_codes, fontsize=10)

# Axis specific changes

#axs[1, 0].legend(handles=custom_handles3, loc="lower right", frameon=True, fontsize=12)
#axs[1, 0].add_artist(legGA)
#axs[0, 0].set_title("Graphene--Amino-acid Interactions", fontsize=15)
axs[0].set_ylabel("$\Delta E^{ads}_{AA-Graph}~[kJ \cdot mol^{-1}]$", fontsize=13)


axs[1].set_ylabel("$\Delta E^{ads}_{AA-H_2O}~[kJ \cdot mol^{-1}]$", fontsize=13)


fig.tight_layout()
fig.savefig(f"{output}/{imagename}.pdf", dpi=600)

print("DONE!")
