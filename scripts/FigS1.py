import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import linregress
from glob import glob

output = "../figures"
imagename = "FigS1"
ylabel = r"PMF$(\xi) (kJ \cdot mol^{-1})$"

files = sorted(glob("../data/pmfs/*.dat"))
files = [file for file in files if "25o" not in file]

# Use this to fill in integrated parts
boundaries = (0.6, 1.1)

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
    "Water":"Polar"
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
    "Water":"Water"
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
    "Water":"H2O"
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
    "Water":"H2O"
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
    "Water":"o"
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
    "Water": 18.02,
}

fig, axs = plt.subplots(figsize=(14, 16.5), ncols=3, nrows=7)

counter = 0
move_next = False

sorted_files = sorted(files, key=lambda x: type_dict[x.split("/")[-1][:-4]])


for file in sorted_files:
    aminoacid = file.split("/")[-1][:-4]
    name = name_dict[aminoacid]

    if move_next:
        counter -= 1
        move_next = False
    elif "(-)" in name or "Histidine (+)" in name:
        move_next = True
    elif "(+)" in name:
        counter -= 1

    ax = axs.flatten()[counter]

    if "(+)" in name or "(-)" in name:
        ls = "dashed"
    else:
        ls ="solid"
        ax.set_title(name)

    t = type_dict[aminoacid]
    c = colors[t]
    data = np.loadtxt(file)
    ax.errorbar(data[:, 0], data[:, 1], data[:, 2], label=name, color=c,
                linewidth=2, elinewidth=2, linestyle=ls)
    ax.axvline(boundaries[0], color="#BBCCEE", linewidth=2)
    ax.axvline(boundaries[1], color="#FFCCCC", linewidth=2)
#    ax.legend(loc="lower right", fontsize=11, frameon=True)
    ax.set_xlim(0, 1.5)
    ax.set_xticks([0, .5, 1., 1.5])
    ax.set_xticklabels([0, 0.5, 1.0, 1.5], fontsize=12.5)
    ax.set_xlabel(r"$\xi (nm)$", fontsize=12.5)
    ax.set_ylim(0, 35)
    ax.set_yticks([0, 10, 20, 30])
    ax.set_yticklabels([0, 10, 20, 30], fontsize=12.5)
    ax.set_ylabel(ylabel, fontsize=11)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    counter += 1

while counter < len(axs.flatten()):
    ax = axs.flatten()[counter]
    fig.delaxes(ax)
    counter += 1

fig.tight_layout()
fig.savefig(f"{output}/{imagename}.pdf", dpi=300)
print("DONE!")
