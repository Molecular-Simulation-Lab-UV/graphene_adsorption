"""
Script to copy data files used to make figures and tables
for the paper.
"""
from glob import glob
from shutil import copy
from os import makedirs

target = "/home/mbarria/Dropbox/Papers/graphene_adsorption_paper/data"
# 1.- Tables
makedirs(f"{target}/tables", exist_ok=True)
source = "/home/mbarria/don-elias/AA_GO0US/tables/2021_12_10_SASA"

copy(f"{source}/Table1.csv", f"{target}/tables/Free_Energy_Decomposition.csv")
copy(f"{source}/Table2.csv", f"{target}/tables/Energy_Decomposition.csv")
copy(f"{source}/Table3.csv", f"{target}/tables/Internal_and_External_Energy.csv")
copy(f"{source}/Table4.csv", f"{target}/tables/Entropy_and_Diffusion.csv")
copy(f"{source}/Table5.csv", f"{target}/tables/SASA.csv")

# 25%
source = "/home/mbarria/don-elias1/AA_GO25US/tables/2021_12_09_SASA"

copy(f"{source}/Table1.csv", f"{target}/tables/Free_Energy_Decomposition_25o.csv")
copy(f"{source}/Table2.csv", f"{target}/tables/Energy_Decomposition_25o.csv")
copy(f"{source}/Table3.csv", f"{target}/tables/Internal_and_External_Energy_25o.csv")
copy(f"{source}/Table4.csv", f"{target}/tables/Entropy_and_Diffusion_25o.csv")
copy(f"{source}/Table5.csv", f"{target}/tables/SASA_25o.csv")

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
