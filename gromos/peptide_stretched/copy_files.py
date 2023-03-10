from glob import glob
from shutil import copy
from os import makedirs, getcwd

source = "../AA_GO0US_Helix/"

#aminoacids = ['Alanine', 'Arginine', 'Arginine_H+', 'Asparagine', 'Aspartic_Acid', 'Aspartic_Acid_H+',
#              'Cysteine_H+', 'Glutamic_Acid', 'Glutamic_Acid_H+', 'Glutamine', 'Glycine', 'Histidine_2H+',
#              'Histidine_H+_ND1', 'Isoleucine', 'Leucine', 'Lysine', 'Lysine_H+', 'Methionine', 'Phenylalanine',
#              'Proline', 'Serine', 'Threonine', 'Tryptophan', 'Tyrosine', 'Valine', 'Water']
aminoacids = ['Alanine']

extensions = [".por", ".rpr", ".cnf", ".top"]

cwd = getcwd()

mk_script = """@sys	US_pep
@bin	/opt/gromos/md++-1.5.0_gcc-8.5.0/bin/md
@version	md++
@dir	/home/mbarria/AA_GO0US_Stretched/AMINOACID/pep/RFOLDER/
@files
  topo	../../AMINOACID_pep_graphene_54a8.top
  input	US_pep.imd
  coord	../../AMINOACID_pep_graphene_h2o_54a8.cnf
  posresspec	../../AMINOACID_pep_graphene_h2o_54a8.por
  refpos	../../AMINOACID_pep_graphene_h2o_54a8.rpr
  disres	../AMINOACID_pep_graphene_h2o_54a8_RDIST.dsr
@template	/home/mbarria/mkscript_mms_leftraru_openMP.lib
@script	1 100"""

for aa  in aminoacids:
    # Copy topology, coordinates, etc
    makedirs(cwd + "/" + aa, exist_ok=True)
    for path in glob(f"{source}{aa}/*"):
        ext = path[-4:]
        if ext in extensions:
            fname = path.split("/")[-1]
            copy(path, f"{cwd}/{aa}/{fname}")
    # Copy distance restrictions
    makedirs(f"{cwd}/{aa}/pep", exist_ok=True)
    for path in glob(f"{source}{aa}/pep/*"):
        ext = path[-4:]
        if ext == ".dsr":
            fname = path.split("/")[-1]
            with open(path, 'r') as file:
                lines = file.readlines()
            with open(f"{cwd}/{aa}/pep/{fname}", 'w') as file:
                for line in lines:
                    newline = line.replace("0.38", "0.379").replace("0.59", "1.449")
                    file.write(newline)
    # Copy mk_script, imd
    for rpath in glob(f"{source}/{aa}/pep/r*"):
        rfolder = rpath.split("/")[-1]
        rdist = rfolder[1:]
        makedirs(f"{cwd}/{aa}/pep/{rfolder}", exist_ok=True)
        #imd
        imd = rpath + "/pep_US.imd"
        with open(imd, 'r') as file:
            lines = file.readlines()
        with open(f"{cwd}/{aa}/pep/{rfolder}/US_pep.imd", 'w') as file:
            for line in lines:
                file.write(line)
        #arg
        with open(f"{cwd}/{aa}/pep/{rfolder}/US_mk_script.arg", 'w') as file:
            file.write(mk_script.replace("RFOLDER", rfolder).replace("RDIST", rdist).replace("AMINOACID", aa))

