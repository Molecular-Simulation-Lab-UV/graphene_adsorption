"""
Script to copy data files used to make figures and tables
for the paper.
"""
from glob import glob
from shutil import copy
from os import makedirs

def copycsv(original, new):
    # to use with latex we need to remove the first line (comment)
    with open(original, 'r') as io:
        lines = io.readlines()
    with open(new, 'w') as io:
        for line in lines:
            if line[0] != '#' and "Water" not in line:
                newline = line.replace("_", " ")
                newline = newline.replace("Acid,", "Acid (-),")
                newline = newline.replace("Acid H+", "Acid")
                newline = newline.replace("2H+", "(+)")
                newline = newline.replace(" H+ ND1", "")
                newline = newline.replace("Cysteine H+", "Cysteine")
                newline = newline.replace("H+", "(+)")
                io.write(newline)

target = "/home/mbarria/Dropbox/Papers/graphene_adsorption_paper/data"
# 1.- Tables
makedirs(f"{target}/tables", exist_ok=True)
source = "/home/mbarria/don-elias/AA_GO0US/tables/2021_12_10_SASA"

copycsv(f"{source}/Table1.csv", f"{target}/tables/Free_Energy_Decomposition.csv")
copycsv(f"{source}/Table2.csv", f"{target}/tables/Energy_Decomposition.csv")
copycsv(f"{source}/Table3.csv", f"{target}/tables/Internal_and_External_Energy.csv")
copycsv(f"{source}/Table4.csv", f"{target}/tables/Entropy_and_Diffusion.csv")
copycsv(f"{source}/Table5.csv", f"{target}/tables/SASA.csv")

# 25%
source = "/home/mbarria/don-elias1/AA_GO25US/tables/2021_12_09_SASA"

copycsv(f"{source}/Table1.csv", f"{target}/tables/Free_Energy_Decomposition_25o.csv")
copycsv(f"{source}/Table2.csv", f"{target}/tables/Energy_Decomposition_25o.csv")
copycsv(f"{source}/Table3.csv", f"{target}/tables/Internal_and_External_Energy_25o.csv")
copycsv(f"{source}/Table4.csv", f"{target}/tables/Entropy_and_Diffusion_25o.csv")
copycsv(f"{source}/Table5.csv", f"{target}/tables/SASA_25o.csv")

# Helix 5-Ala

# Stretched 5-Ala


# 2.- PMFs
makedirs(f"{target}/pmfs", exist_ok=True)

source = "/home/mbarria/don-elias/AA_GO0US"
aminoacids = [aa.split("/")[-2] for aa in glob(f"{source}/*/*.top") if "Nosé" not in aa]
for aa in aminoacids:
    copy(f"{source}/{aa}/wham/disres_freefile.dat", f"{target}/pmfs/{aa}.dat")

source = "/home/mbarria/don-elias1/AA_GO25US"
aminoacids = [aa.split("/")[-2] for aa in glob(f"{source}/*/*.top") if "Nosé" not in aa]
for aa in aminoacids:
    copy(f"{source}/{aa}/wham/disres_freefile.dat", f"{target}/pmfs/{aa}_25o.dat")
