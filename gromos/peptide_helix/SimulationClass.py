#!/usr/bin/env python

import os, subprocess, glob
import random
import MDAnalysis
import numpy as np
import math
import re

##Forcefield name and location
ff = '54a8'
ifp = os.getcwd() + '/FF/54a8.ifp'
mtb = os.getcwd() + '/FF/54a8_CCC_OX.mtb'
##Other Variables
seed = '0'
if len(seed) != 6:
    seed = '%i' % random.randint(100000,999999)
genbox = ['4.935500000','4.865000000','5.000000000']
##Directories
prepDir = os.getcwd() + '/Prep/'
frameoutDir = os.getcwd() + '/FRAMEOUTS/'
templateDir = os.getcwd() + '/Templates/'
outputDir = os.getcwd() + '/Simulation/'
leftraruDir = '/home/username/AA_GO0US_Helix/'
gromos = '/home/mbarria/gromos++-1.4.0/x86/'
gromosJAG = '/home/mbarria/gromos++-1.2.0_JAG/x86/'
md_leftraru = '/home/username/opt/md++-1.4.0/bin/md'
md = '/home/mbarria/md++/bin/md'
#Coordinates need to be supplied, I use VMD to create them.
vmd =  '/home/mbarria/Peptidos/vmd_pdb_penta/'
# Ion dictionaries to ensure prober counter ions
ion_file = {'+e':'Cl','-e':'Na'}
ion_arg = {'+e':'CL-','-e':'NA+'}
ion_charge = {'+e':'negative','-e':'positive'}

w0_dict = {0.3: 25,
           0.35: 20,
           0.4: 20,
           0.45: 20,
           0.5: 20,
           0.55: 20,
           0.6: 15,
           0.65: 15,
           0.7: 15,
           0.75: 15,
           0.8: 10,
           #           0.85:10,
           0.9: 10,
           #           0.95:10,
           1.0: 10,
           #           1.05:10,
           1.1: 10,
           #           1.15:10,
           1.2: 10,
           #           1.25:10,
           1.3: 10,
           #           1.35:10,
           1.4: 10,
           #           1.45:10,
           1.5: 10,
           }

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    deltas = nodes - node
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)

###General purpose functions

def run_program(program, arg, output = False, name = 'test',  directory = prepDir, gromos = gromos, folder = 'programs/'):
    #Run gromos programs
    if output:
        command = '%s @f %s > %s' % (gromos + folder + program, arg, output)
    else:
        command = '%s @f %s' % (gromos + folder + program, arg)
    subprocess.run(command, shell=True, cwd = directory + name)

class Simulation():
    def __init__(self, name, kind, sequence, charge, ff = ff, ifp = ifp, mtb = mtb, seed='0', graphene = False, prepDir = prepDir, outputDir = outputDir, templateDir = templateDir, frameoutDir = frameoutDir, leftraruDir = leftraruDir, gromos = gromos, gromosJAG = gromosJAG, vmd = vmd, md = md, md_leftraru = md_leftraru, local_elevation = False, umbrella_sampling = False, restrain_graphene = True, genbox = genbox, aaUser = False, pepUser = False, sidechain = True, w0_dict = w0_dict, stop = 100, LINCS = False, block = 'CCC'):
        ##Definition
        self.name = name
        self.kind = kind
        self.sequence = sequence
        self.charge = charge
        self.ff = ff
        self.ifp = ifp
        self.mtb = mtb
        if len(seed) != 6:
            self.seed = '%i' % random.randint(100000,999999)
        else:
            self.seed = seed
        self.graphene = graphene
        self.prepDir = prepDir
        self.outputDir = outputDir
        self.templateDir = templateDir
        self.frameoutDir = frameoutDir
        self.leftraruDir = leftraruDir
        self.gromos = gromos
        self.gromosJAG = gromosJAG
        self.vmd = vmd
        self.md = md
        self.md_leftraru = md_leftraru
        self.local_elevation = local_elevation
        self.umbrella_sampling = umbrella_sampling
        self.restrain_graphene = restrain_graphene
        self.genbox = genbox
        self.aaUser = aaUser
        self.pepUser = pepUser
        self.sidechain = sidechain
        self.w0_dict = w0_dict
        self.stop = stop
        self.LINCS = LINCS
        self.block = block
        if self.block != 'CCC':
            self.oxidize = True
        else:
            self.oxidize = False

        ##Filenames
        self.top = Simulation.get_top(self.name, self.kind, self.ff)
        self.cnf = Simulation.get_cnf(self.name, self.kind, self.ff)
        self.cnf_gch = 'gch_%s' % self.cnf
        self.cnf_min = self.cnf.replace('.cnf','_min.cnf')
        self.cnf_solvent = Simulation.get_cnf(self.name, self.kind, self.ff, solvent = True)
        self.top_charged = Simulation.get_top(self.name, self.kind, self.ff, self.charge)
        self.cnf_charged = Simulation.get_cnf(self.name, self.kind, self.ff, solvent = False, charge = self.charge, graphene = False)
        self.cnf_charged_solvent = Simulation.get_cnf(self.name, self.kind, self.ff, solvent = True, charge = self.charge, graphene = False)
        self.top_graph = Simulation.get_top(self.name, self.kind, self.ff, False, self.graphene)
        self.cnf_graph = Simulation.get_cnf(self.name, self.kind, self.ff, solvent = False, charge = False, graphene = self.graphene)
        self.cnf_graph_solvent = Simulation.get_cnf(self.name, self.kind, self.ff, solvent = True, charge = False, graphene = self.graphene)
        self.top_graph_charged = Simulation.get_top(self.name, self.kind, self.ff, self.charge, self.graphene)
        self.cnf_graph_charged = Simulation.get_cnf(self.name, self.kind, self.ff, solvent = False, charge = self.charge, graphene = self.graphene)
        self.cnf_graph_charged_solvent = self.cnf_graph_charged = Simulation.get_cnf(self.name, self.kind, self.ff, solvent = True, charge = self.charge, graphene = self.graphene)

        ##Default values for functions

        self.make_folders = lambda dir = self.prepDir, EQ = False, separate = False : Simulation.make_folders(self.name, self.kind, dir, EQ, separate)
        self.frameout = lambda destination, cnf, charge = False, solvent = False, run = True, top = False, graphene = False, gather = False, move = True : Simulation.frameout(self.name, self.kind, destination, cnf, self.prepDir, self.ff, charge, solvent, self.gromos, run, top, graphene, gather, move)
        self.count_atoms = lambda cnf, charge = False, solvent = False, graphene = False, top = False : Simulation.count_atoms(self.name, self.kind, cnf, self.prepDir, self.frameoutDir, self.ff, charge, solvent, self.gromos, graphene, top)
        self.check_omd = lambda directory : Simulation.check_omd(self.name, self.kind, directory)
        self.get_genbox = lambda cnf : Simulation.get_genbox(self.name, cnf, self.prepDir)
        self.make_top = lambda run = True : Simulation.make_top(self.name, self.kind, self.sequence, self.top, self.prepDir, self.ifp, self.mtb, run)
        self.make_top_ion = lambda run = True : Simulation.make_top_ion(self.name, self.charge, self.prepDir, self.ifp, self.mtb, run)
        self.com_top_ion = lambda run = True : Simulation.com_top_ion(self.name, self.kind, self.charge, self.top_graph_charged, self.prepDir, self.ff, run, self.graphene)
        self.fix_pdb = lambda x = 0 : Simulation.fix_pdb(self.name, self.sequence, self.prepDir, self.vmd)
        self.pdb2g96 = lambda run = True :  Simulation.pdb2g96(self.name, self.sequence, self.kind, self.top, self.prepDir, self.ff, self.gromos, run)
        self.gch = lambda run = True : Simulation.gch(self.name, self.kind, self.top, self.prepDir, self.ff, self.gromos, run)
        self.make_imd = lambda cnf, template, destination = self.prepDir, directory = self.prepDir, solvent = False, graphene = False, charge = False, TSTRA = '500', TSENERG = '100', TAU = '0.1', filename = False, top = False, add_restraint = False, sidechain = False, r = False: Simulation.make_imd(self.name, self.kind, cnf, template, destination, self.prepDir,self.frameoutDir, self.ff, self.gromos, sidechain, solvent, graphene, charge,  TSTRA, TSENERG, self.seed, filename, top, TAU, add_restraint, r)
        self.min_run = lambda top, cnf, imd, run = True : Simulation.min_run(self.name, self.kind, top, cnf, imd, self.prepDir, self.ff, self.gromos, self.md, run = True)
        self.sim_box = lambda top, cnf, genbox = False, run = True : Simulation.sim_box(self.name, self.kind, top, cnf, self.prepDir, self.templateDir, self.gromos, genbox, run)
        self.restrain = lambda cnf, destination = self.prepDir, restrain_graphene = False : Simulation.restrain(self.name, self.kind, cnf, destination, self.graphene, restrain_graphene)
        self.solv_min = lambda top, cnf, imd, graphene = False, run = True : Simulation.solv_min(self.name, self.kind, top, cnf, imd, self.prepDir, graphene, run)
        self.add_ion = lambda run = True :  Simulation.add_ion(self.name, self.kind, self.top_graph, self.cnf_graph_solvent, self.charge, self.prepDir, self.ff, self.gromos, run, self.graphene)
        self.make_graphene = lambda x = 0 : Simulation.make_graphene(self.name, self.kind, self.cnf_solvent, self.prepDir, self.templateDir, self.ff, self.gromos, self.genbox, self.block)
        self.graphite = lambda lastgraph, top = 'graphene_separate_%s_%s.top' % (self.ff, self.kind), cnf = 'graphene_%s.cnf' % self.kind : Simulation.graphite(self.name, self.kind, top, cnf, self.prepDir, False, True, self.ff, self.gromos, self.gromosJAG, lastgraph, self.oxidize)
        self.fix_graphene_top = lambda name = self.name, top = 'graphene_%s_%s.top' % (self.ff, self.kind), coord = 'graphene_%s.cnf' % self.kind, directory = self.prepDir : Simulation.fix_graphene_top(name, top, coord, directory)
        self.make_box = lambda cnf = self.cnf_min, below = '1.500000000', SSD0 = '0.300000000' : Simulation.make_box(self.name,self.kind, cnf, self.prepDir, False, self.ff, self.gromos, True, below, SSD0)
        self.graphene_com_top = lambda top = self.top : Simulation.graphene_com_top(self.name, self.kind, top, self.prepDir, False, False, self.frameoutDir, self.ff, self.gromos)
        self.equilibrate = lambda top, cnf, run = True : Simulation.equilibrate(self.name, self.kind, top, cnf, self.prepDir, self.gromos, run, self.graphene, self.restrain_graphene)
        self.copy_files = lambda eq_cnf = self.prepDir + '%s/eq_%s_5.cnf' % (self.name, self.kind) : Simulation.copy_files(self.name, self.kind, eq_cnf, self.cnf_graph_charged_solvent, self.top_graph_charged, self.prepDir, self.outputDir)
        self.leftraru_mk_script = lambda cnf = self.cnf_graph_charged_solvent, top = self.top_graph_charged, job = 'md', r = '' : Simulation.leftraru_mk_script(self.name, self.kind, cnf, top, self.leftraruDir, self.outputDir, self.md_leftraru, self.pepUser, self.aaUser, self.graphene, self.restrain_graphene, self.local_elevation, self.umbrella_sampling, self.stop, job, r)
        self.leftraru_jobs = lambda cores = 10, job = 'md' , r = False: Simulation.leftraru_jobs(self.name, self.kind, self.templateDir, self.leftraruDir, self.outputDir, cores, self.pepUser, self.aaUser, self.graphene, job, r, stop = self.stop)
        self.local_elevation_file = lambda cnf = self.cnf_graph_charged_solvent : Simulation.local_elevation_file(self.name, self.kind, cnf, self.outputDir)
        self.umbrella_sampling_file = lambda r, outpuDir = self.outputDir, cnf = self.cnf_graph_charged_solvent : Simulation.umbrella_sampling_file(self.name,self.kind, cnf, r, outputDir, self.w0_dict)

    ##Files and file operations
    def make_folders(name, kind, directory = prepDir, EQ = False, separate = False):
        """
        Make folders to store preparation files
        """
        folder = directory + name + '/'
        if EQ:
            os.makedirs(folder + 'EQ_%s' % kind,exist_ok=True)
        elif separate:
            os.makedirs(folder + '%s' % kind,exist_ok=True)
        else:
            os.makedirs(folder,exist_ok=True)


    def get_top(name, kind, ff=ff, charge = False, graphene = False):
        """
        Format topology filename
        """
        top = '%s_%s' % (name, kind)
        if charge:
            top += '_%s' % ion_file[charge]
        if graphene:
            top += '_graphene'
        top += '_%s.top' % ff
        return top

    def get_cnf(name, kind, ff=ff, solvent = False, charge = False, graphene = False):
        """
        Format coordinate filename
        """
        cnf = '%s_%s' % (name, kind)
        if charge:
            cnf += '_1%s' % ion_file[charge]
        if graphene:
            cnf += '_graphene'
        if solvent:
            cnf += '_h2o'
        cnf += '_%s.cnf' % ff
        return cnf

    def frameout(name, kind, destination, cnf, directory = prepDir, ff = ff, charge = False, solvent = False, gromos = gromos, run = True, top = False, graphene = False, gather = False, move = True):
        """
        Creates a FRAMEOUT pdb and stores it in destination
        """
        file = 'frameout_%s.arg' % cnf[:-4]
        if not top: top = Simulation.get_top(name, kind, ff, charge, graphene)
        if gather:
    #        boundary = 'r 3 refg %s'  % cnf
            boundary = 'r 2'
        elif solvent:
            boundary = 'r'
        else:
            boundary = 'v'
        with open(directory + name + '/%s' % file,'w') as arg:
            arg.write('@topo \t%s\n' % top)
            arg.write('@pbc \t%s\n' % boundary)
            arg.write('@outformat \tpdb\n')
            arg.write('@notimeblock\n')
            arg.write('@traj \t%s\n' % (cnf))
            arg.write('@include ALL')
        if run:
            run_program('frameout',file, name = name, directory = directory, gromos = gromos)
            if move:
                Simulation.make_folders(name, kind, frameoutDir + destination)
                subprocess.run(['mv',directory + name + '/FRAME_00001.pdb', frameoutDir + destination + name + '/FRAME_%s_%s.pdb' % (name, kind)])

    def count_atoms(name, kind, cnf, directory = prepDir, frameoutDir = frameoutDir, ff = ff, charge = False, solvent = False, gromos = gromos, graphene = False, top = False):
        """
        counts last atoms and molecules (solute, then ion, then solvent) of a given cnf file
        """
        Simulation.frameout(name,kind,'TEMP/',cnf,directory,ff,charge,solvent,gromos, graphene = graphene, top = top)
        with open(frameoutDir + 'TEMP/' + name + '/FRAME_%s_%s.pdb' % (name, kind)) as pdb:
            tempAtom = tempMol = ''
            Atoms = []
            Mols = []
            for line in pdb:
                if 'TER' in line:
                    Atoms.append(tempAtom)
                    Mols.append(tempMol)
                else:
                    try:
                        tempAtom = line.split()[1]
                        tempMol = line.split()[4]
                    except IndexError:
                        pass

        return [Atoms, Mols]

    def check_omd(name, kind, directory):
        """
        Checks all omd files within a directory
        """
        for i in glob.iglob(directory +'**.omd'):
            with open(i,'r') as omd:
                lines = omd.readlines()
                if 'success' in lines[-2] or 'success' in lines[-3]:
                    print(i, 'OK')
                else:
                    print(i, 'FAIL')

    def get_genbox(name, cnf, directory = prepDir):
        with open(directory + name + '/%s' % cnf,'r') as file:
            flag = False
            genbox = ''
            for line in file:
                if 'GENBOX' in line:
                    flag = True
                elif 'END' in line:
                    flag = False
                elif flag:
                    genbox += line
        return genbox

    ###Topology

    ### Turns out the graphite program has serious issues with the direction of the improper dihedrals
    ### These functions aim to fix that
    def fix_graphene_top(name, top, coord, directory = prepDir):
        #fixed = top.replace('54a8','54a8_fixed')
        energies = []
        angles = []
        command = ['check_top','@topo',directory + name + '/' + top,'@coord',directory + name + '/'  + coord,'>',directory + name + '/'  + 'check_top.out']
        subprocess.getoutput(' '.join(command))

        with open(directory + name + '/'  + 'check_top.out') as file:
            dihedrals = False
            read = False
            for line in file:
                if len(line) < 2:
                    continue
                if 'DIHEDRAL ANGLES' in line:
                    break
                if 'IMPROPER DIHEDRALS' in line:
                    dihedrals = True
                    continue
                if dihedrals and '   1' in line:
                    read = True
                if read:
                    energy = float(line.split()[-1])
                    angle = line.strip()[1:25].strip().replace('-',' ')
                    angle = [int(x) for x in angle.split()]
                    energies.append(energy)
                    if energy > 1:
                        #print(line)
                        angles.append(angle)
        #            if energy > 1:
        #                print(angle, energy)
                else:
                    continue

        with open(directory + name + '/'  + top) as original:
            lines = original.readlines()
        with open(directory + name + '/'  + top,'w') as new:
            for line in lines:
                try:
                    angle = [int(x) for x in line.split()[0:4]]
                    if angle in angles:
                        newline = line[:-2] + '4' + '\n'
                        #print(newline)
                        new.write(newline)
                    else:
                        new.write(line)
                except (IndexError, ValueError):
                    new.write(line)

    def make_top(name, kind, sequence, top = False, directory = prepDir, ifp = ifp, mtb = mtb, run = True):
        """
        Runs make_top and stores an .arg file
        """
        file = 'make_top_%s.arg' % kind
        sequence = ' '.join(sequence)
        with open(directory + name + '/' + file, 'w') as arg:
            arg.write('# %s based peptide' % (name))
            arg.write('\n@build \t %s' % (mtb))
            arg.write('\n@param \t %s' % (ifp))
            arg.write('\n@seq \t %s' % (sequence))
            arg.write('\n@solv \t H2O')
        if run:
            if not top:
                top = Simulation.get_top(name, kind, ff)
            run_program('make_top', file, top , name, directory, gromos)

    def make_top_ion(name, charge, directory = prepDir, ifp = ifp, mtb = mtb, run = True):
        #Runs make top for ions when necessary
        file = 'make_top_%s.arg' % ion_file[charge]
        with open(directory + name + '/' + file,'w') as arg:
            arg.write('@build\t%s\n' % mtb)
            arg.write('@param\t%s\n' % ifp)
            arg.write('@seq\t%s\n' % ion_arg[charge])
            arg.write('@solv\tH2O')
        if run:
            run_program('make_top', file, '%s_%s.top' % (ion_file[charge], ff),name , directory, gromos)

    def com_top_ion(name, kind, charge, top = False, directory = prepDir, ff = ff, run = True, graphene = False):
        """
        Combines ion and solute topology
        """
        file = 'com_top_%s_ion.arg' % kind
        with open(directory + name + '/' + file, 'w') as arg:
            arg.write('@topo\t%s %s_%s.top\n' % (Simulation.get_top(name,kind,ff, graphene = graphene), ion_file[charge], ff))
            arg.write('@param\t1\n')
            arg.write('@solv\t1\n')
        if run:
            if not top:
                top = Simulation.get_top(name, kind, ff, charge, graphene)
            run_program('com_top', file, top, name, directory, gromos)


    def fix_pdb(name, sequence, directory = prepDir, vmd = vmd):
        """
        Pdb files by vmd need "Fixing"
        This replaces values as used by vmd when generating coordinates to values as used by gromos
        """
        for i in glob.iglob(vmd + '*/*.pdb') :
            destination = '/'.join(i.split('/')[5:])
            destination = directory + '_'.join(destination.split())
            subprocess.run(['cp',i,destination])
        #Fix aa coordinates
        with open(directory + name + '/aa.pdb', 'r') as file:
            with open(directory + name + '/aa_fixed.pdb', 'w') as pdb:
                pdb.write(file.readline())
                resnum = 0
                lineCounter = 0
                for line in file:
                    try:
                        resnum = 0
                        if line.split()[2] == 'N': #Si falla el error temrinara el archivo
                            resnum += 1
                    #Enumerar oxígenos
                        newline = line.replace('HC1','O2 ')
                        newline = newline.replace('HC','O2')
                        newline = 'ATOM' + newline[4:].replace(' O ',' O1')
                        pdb.write(newline)
                    except IndexError:
                        pdb.write('END')
                        pass
        #Fix peptide coordinates
        coord_file = open(directory + name + '/peptide_coord.pdb', 'r').readlines()
        with open(directory + name + '/peptide_nombres.pdb', 'r') as file:
            with open(directory + name + '/pep_fixed.pdb', 'w') as pdb2:
                pdb2.write(file.readline())
                resnum = 0
                lineCounter = 0
                for line in file:
                    try:
                        #Poner residuo correcto
                        if line.split()[2] == 'N':
                            resnum += 1
                        newline = line.replace('THR',sequence[resnum])
                        newline = newline.replace('X   1','X   ' + str(resnum))
        #                    newline = newline.replace('HN ','H  ')
                        #Enumerar oxígenos del último residuo
                        if resnum == 9:
                            newline = newline.replace('HC1','O2 ')
                            newline = newline.replace('HC','O2')
                            newline = 'ATOM' + newline[4:].replace('O ','O1')
                        #Corregir Coordenadas
                        lineCounter += 1
                        newline = newline[0:27] + coord_file[lineCounter][27:]

                        pdb2.write(newline)
                    except IndexError:
                        pdb2.write('END')
                        pass

    def pdb2g96(name, sequence, kind, top = False, directory = prepDir, ff = ff, gromos = gromos, run = True):
        #Converts VMD pdb files to gromos cnf files
        res = sequence[(len(sequence))//2]
        #Lib File
        with open(directory + name + '/pdb2g96_%s.lib' % kind,'w') as lib:
            lib.write('TITLE')
            lib.write('\n  Library file to specify residue and atom names')
            lib.write('\n  that are to be considered similar')
            lib.write('\nEND')
            lib.write('\nRESIDUENAMELIB')
            lib.write('\n# name in my pdb   name in topology')
            if len(res) == 4:
                lib.write('\n\t' + res[0:3] + '\t' + res)
            lib.write('\nEND')
            lib.write('\nATOMNAMELIB')
            lib.write('\n#\tresidue in topology \tatom in pdb \tatom in topology')
            lib.write('\n\t' + res + '\t O \t O1')
            lib.write('\n\t' + res + '\t HC1 \t O2')
            lib.write('\nEND')
        #arg file
        file = 'pdb2g96_%s.arg' % kind
        with open(directory + name + '/%s' % file,'w') as arg:
            if not top: top = Simulation.get_top(name, kind, ff)
            arg.write('@topo \t%s' % top)
            arg.write('\n@pdb  \t%s_fixed.pdb' % kind)
            arg.write('\n@lib  \tpdb2g96_%s.lib' % kind)
        if run:
            run_program('pdb2g96',file, 'pdb2g96_%s.cnf' % kind,name,  directory, gromos)

    def gch(name, kind, top = False, directory = prepDir, ff = ff, gromos = gromos, run = True):
        #Adds hydrogen to coordinate file
        file = 'gch_%s.arg' % kind
        with open(directory + name +'/%s' % file,'w') as arg:
            if not top: top = Simulation.get_top(name, kind, ff)
            arg.write('@topo \t%s' % top)
            arg.write('\n@pos \tpdb2g96_%s.cnf' % kind)
            arg.write('\n@tol \t0.1')
        if run:
            run_program('gch',file, 'gch_%s_%s_%s.cnf' % (name, kind, ff), name, directory, gromos)


    def make_imd(name, kind, cnf, template, destination, directory = prepDir, frameoutDir = frameoutDir, ff = ff, gromos = gromos, sidechain = False, solvent = False, graphene = False, charge = False, TSTRA = '500', TSENERG = '100', seed = seed, filename = False, top = False, TAU = '0.1', add_restraint = False, r = False):
        #Basesd on atom count and other parameters replaces variables in a template imd and creates one for the corresponding molecule
        NRE = []
        atoms, mols = Simulation.count_atoms(name, kind, cnf, directory, frameoutDir, ff, charge, solvent, gromos, graphene, top = top)
        LASTSOLUTEATOM = atoms[0]
        if graphene:
            GRAPHENE = atoms[1]
        if charge:
            ION = atoms[-2]
        if solvent:
            LASTSOLVENTMOL = mols[-1]
            LASTSOLVENTATOM = atoms[-1]
        else:
            LASTSOLVENTMOL = '0'
            LASTSOLVENTATOM = '0'
        if sidechain and name != 'Glycine':
            if kind == 'aa':
                NRE.append('4')
                NRE.append(str(int(LASTSOLUTEATOM) - 4))
                NRE.append(LASTSOLUTEATOM)
            else:
                NRE.append('33')
                NRE.append(str(int(LASTSOLUTEATOM) - 33))
                NRE.append(LASTSOLUTEATOM)
        else:
            NRE.append(LASTSOLUTEATOM)
        if graphene: NRE.append(GRAPHENE)
        if charge: NRE.append(ION)
        if solvent: NRE.append(LASTSOLVENTATOM)
        if not filename:
            filename = '%s_%s' % (kind, template)
        if not r:
            outputFile = destination + name + '/%s' % filename
        else:
            outputFile = destination + name + '/%s/r%s/%s' % (kind,('%1.3f' % r).replace('.','_'),filename)
        with open(templateDir + template,'r') as tFile:
            with open(outputFile,'w') as imd:
                if graphene:
                    bathnum = 4
                    dofset = 3
                else:
                    bathnum = dofset = 2
                for line in tFile:
                    line = line.replace('SEED',seed)
                    line = line.replace('TSTRA',TSTRA)
                    line = line.replace('TSENERG',TSENERG)
                    line = line.replace('LASTSOLVENTMOL',LASTSOLVENTMOL)
                    line = line.replace('LASTSOLVENTATOM',LASTSOLVENTATOM)
                    line = line.replace('BATHNUM',str(bathnum))
                    line = line.replace('DOFSET',str(dofset))
                    if graphene: line = line.replace('LASTSOLUTEATOM',GRAPHENE)
                    else: line = line.replace('LASTSOLUTEATOM',LASTSOLUTEATOM)
                    if 'NRETOTAL' in line:
                        line = '\t%s\t%s\n' % (len(NRE),'\t'.join(NRE))
                    if '     0.1     ' in line:
                        line = line[0:18].replace('0.1',TAU) * bathnum + '\n'
                    if 'BATHLINE' in line and not graphene:
                        line = '\t%s\t1\t1\t%s\t2\t2\n' % (LASTSOLUTEATOM, LASTSOLVENTATOM)
                    elif 'BATHLINE' in line and graphene:
                        line = '\t%s\t1\t1\t%s\t2\t3\t%s\t4\t4\n' % (LASTSOLUTEATOM, GRAPHENE, LASTSOLVENTATOM)
                    imd.write(line)
            if add_restraint:
                with open(outputFile,'a') as imd:
                    imd.write('\nPOSITIONRES\n')
                    imd.write('#\tNTPOR\tNTPORB\tNTPORS\tCPOR\n')
                    imd.write('\t3\t1\t0\t25000\n')
                    imd.write('END\n')
#            if local_elevation:
#                with open(destination + name + '/%s' % filename,'a') as imd:
#                    imd.write('LOCALELEV\n')
#                    imd.write('#\tNTLES\tNLEPOT\tNTLESA\tNTWLE\n')
#                    imd.write('\t1\t1\t2\t100\n')
#                    imd.write('#\tNLEPID[1..NLEPOT]\n')
#                    imd.write('\t1\n')
#                    imd.write('#\tNLEPFR[1..NLEPOT]\n')
#                    imd.write('\t0\n')
#                    imd.write('END\n')

    def min_run(name, kind, top, cnf, imd, directory = prepDir, ff = ff, gromos = gromos, md = md, run = True):
        #Runs a minimization Simulation for the molecule in a vacuum
        with open(directory + name + '/em_%s.run' % kind,'w') as run:
            run.write('#!/bin/sh\n')
            run.write('GROMOS=%s\n' % md)
            run.write('$GROMOS \\\n')
            run.write('  @topo %s \\\n' % top )
            run.write('  @conf %s \\\n' % cnf )
            run.write('  @fin %s \\\n' % cnf.replace('.cnf','_min.cnf').replace('gch_','') )
            run.write('  @input %s > %s' % (imd, imd.replace('imd','omd')))
        subprocess.run(['chmod','+x',directory + name + '/em_%s.run' % kind], cwd = directory + name)
        if run:
            subprocess.run(['./em_%s.run' % kind], cwd = directory + name)

    def sim_box(name, kind, top, cnf, directory = prepDir, templateDir = templateDir, gromos = gromos, genbox = False, run = True):
        #Puts the solvent in a solute box
        if kind == 'pep':
            minwall = '0.8'
        else:
            minwall = '1.3'
        file = 'sim_box_%s.arg' % kind
        with open(directory + name + '/' + file, 'w') as arg:
            arg.write('@topo\t%s\n' % top)
            arg.write('@pos\t%s\n' % cnf)
            arg.write('@solvent\t%s\n' % (templateDir + 'spc.cnf'))
            if not genbox:
                arg.write('@pbc\tr\n')
                arg.write('@minwall\t%s\n' % minwall)
                arg.write('@rotate\n')
            elif genbox:
                arg.write('@pbc\tr 3 refg %s\n' % cnf)
                arg.write('@boxsize\n')
            arg.write('@thresh\t0.23\n')

        if run:
            run_program('sim_box',file,file.replace('arg','cnf'),name, directory,gromos)

    def restrain(name, kind, cnf, directory = prepDir, graphene = False, restrain_graphene = False):
        #Creates restraining files for Simulations of solvent
        with open(directory + name + '/%s' % cnf,'r') as coord:
            with open(directory + name + '/%s' % cnf.replace('cnf','por'),'w') as por:
                write_flag = False
                por.write('TITLE\nsolute atoms to be positionally restrained\nEND\n')
                if not restrain_graphene:
                    for line in coord:
                        if 'POSITION' in line and 'REF' not in line:
                            write_flag = True
                            por.write('POSRESSPEC\n')
                        elif write_flag == True and graphene and 'SOLV'  in line:
                            write_flag = False
                            por.write('END\n')
                        elif write_flag == True and not graphene and 'CCC'  in line:
                            write_flag = False
                            por.write('END\n')
                        elif write_flag == False:
                            pass
                        elif write_flag == True:
                            por.write(line)
                elif restrain_graphene:
                    for line in coord:
                        if 'POSITION' in line and 'REF' not in line:
                            write_flag = True
                            por.write('POSRESSPEC\n')
                        elif write_flag and 'SOLV'  in line:
                            write_flag = False
                            por.write('END\n')
                        elif not write_flag:
                            pass
                        elif write_flag and 'CCC' in line:
                            por.write(line)
        with open(directory + name + '/%s' % cnf,'r') as coord:
            with open(directory + name + '/%s' % cnf.replace('cnf','rpr'),'w') as rpr:
                write_flag = False
                rpr.write('TITLE\nreference positions for restraining solute atoms\nEND\n')
                for line in coord:
                    if 'POSITION' in line and 'REF' not in line:
                        write_flag = True
                        rpr.write('REFPOSITION\n')
                    elif 'END' in line and write_flag:
                        write_flag = False
                        rpr.write('END\n')
                    elif write_flag == False:
                        pass
                    elif write_flag == True:
                        rpr.write(line)

    def local_elevation_file(name, kind, cnf, directory = outputDir):
        led = cnf.replace(".cnf",".led")
        lud = cnf.replace(".cnf",".lud")
        with open(directory + name + '/%s' % led,'w') as led:
            led.write("TITLE\n")
            led.write("Distance definition file for local elevation.\n")
            led.write('END\n')
            led.write("LOCALELEVSPEC\n")
            led.write("# here we only define one distance from the center of mass of each molecule\n")
            led.write("1\t2\t4\t100\t1\t1\n")
            led.write('END\n')
        with open(directory + name + '/%s' % lud,'w') as lud:
            lud.write("TITLE\n")
            lud.write("Local elevation umbrella definition file\n")
            lud.write('END\n')
            lud.write("LEUSBIAS\n")
            lud.write("#\tNUMB\n")
            lud.write("\t1\n")
            lud.write("#\tNLEPID\tNLEDIM\tCLES\n")
            lud.write("\t1\t1\t0.001\n")
            lud.write("#\tVARTYPE(N)\tNTLEFU(N)\tWLES(N)\tRLES(N)\tNGRID(N)\tGRIDMIN(N)\tGRIDMAX(N)\n")
            lud.write("\t2\t0\t1.0     1.0       60        0.0        3.0\n") #TODO!!
            lud.write("#\tNCONLE\n")
            lud.write("\t0\n")
            lud.write("#\tNVISLE\tICONF\n")
            lud.write("END\n")

    def umbrella_sampling_file(name,kind, cnf, r0, directory = outputDir, w0_dict = w0_dict):
        dsr = cnf.replace(".cnf","_%s.dsr" % ('%1.3f' % r0).replace('.','_'))
        carbons = Simulation.find_alpha_carbon(directory + name + '/' + cnf, kind)
        with open(directory + name + '/%s/%s' % (kind,dsr),'w') as dsr:
            dsr.write("TITLE\n")
            dsr.write("Distance restraint specification for umbrella sampling\n")
            dsr.write('END\n')
            dsr.write("DISTANCERESSPEC\n")
            dsr.write("#\tDISH\tDISC\n")
            dsr.write("\t0.1\t0.153\n") #TODO!
            if kind=="pep":
                dsr.write("#\tNDR\n")
                dsr.write("\t6\n") # Add restraints between all alpha carbons + center alpha to graphene
            dsr.write("#\ti j k l type\ti j k l type\tr0\tw0\trahn")
            if kind == "aa":
                dsr.write("\n %s 0 0 0 0\t %s -2\t %1.3f %i 60\n" % (carbons[0],Simulation.find_ccc_corners(directory + name + '/' + cnf), r0, w0_dict[round(r0,4)]))
            elif kind=="pep":
                dsr.write("\n %s 0 0 0 0\t %s -2\t %1.3f %i 60\n" % (carbons[2],Simulation.find_ccc_corners(directory + name + '/' + cnf), r0, w0_dict[round(r0,4)]))               
                dsr.write(f"\t{carbons[0]} 0 0 0 0\t{carbons[1]} 0 0 0 0\tr0.38\t5\t60\n")
                dsr.write(f"\t{carbons[1]} 0 0 0 0\t{carbons[2]} 0 0 0 0\tr0.38\t5\t60\n")
                dsr.write(f"\t{carbons[2]} 0 0 0 0\t{carbons[3]} 0 0 0 0\tr0.38\t5\t60\n")
                dsr.write(f"\t{carbons[3]} 0 0 0 0\t{carbons[4]} 0 0 0 0\tr0.38\t5\t60\n")
                dsr.write(f"\t{carbons[0]} 0 0 0 0\t{carbons[4]} 0 0 0 0\tr0.59\t5\t60\n")
            dsr.write("END\n")

    def find_ccc_corners(cnf):
        px = []
        py = []
        num = []
        with open(cnf,'r') as file:
            ccc = False
            for line in file:
                if ccc and 'END' in line:
                    break
                elif 'CCC' in line:
                    ccc = True
                    n, x, y = line.split()[3:6]
                    px.append(float(x))
                    py.append(float(y))
                    num.append(int(n))
        xy = np.asarray([[x, py[ind]] for ind, x in enumerate(px)])

        c1 = np.asarray([np.min(px), np.min(py)])
        c2 = np.asarray([np.min(px), np.max(py)])
        c3 = np.asarray([np.max(px), np.min(py)])
        c4 = np.asarray([np.max(px), np.max(py)])
        corners = []
        for i in [c1,c2,c3,c4]:
            corners.append(closest_node(i, xy))
        return ' '.join([str(x) for x in (np.asarray(corners) + num[0])])

    def find_alpha_carbon(cnf, kind = 'aa'):
        carbons = []
        with open(cnf,'r') as file:
            position = False
            if kind == 'aa':
                for line in file:
                    if 'POSITION' in line:
                        position = True
                    elif position and 'CA' in line:
                        carbons.append(line.split()[3])
                        return carbons
            else:
                for line in file:
                    if 'POSITION' in line:
                        position = True
                    elif position and 'CA' in line:
                        carbons.append(line.split()[3])
                    elif 'CCC' in line:
                        return carbons
    def solv_min(name, kind, top, cnf, imd, directory = prepDir, graphene = False, run = True):
        #Runs a minimization Simulation for the solvent
        file = 'em_solvent_%s.run' % kind
        with open(directory + name + '/%s' % file, 'w') as run:
            run.write('#!/bin/sh\n')
            run.write('GROMOS=/home/mbarria/md++/bin/md\n')
            run.write('$GROMOS \\\n')
            run.write('  @topo %s \\\n' % top)
            run.write('  @conf %s \\\n' % cnf)
            run.write('  @fin %s \\\n' % Simulation.get_cnf(name, kind, ff, solvent = True, graphene = graphene))
            run.write('  @refpos %s \\\n' % cnf.replace('cnf','rpr'))
            run.write('  @posresspec %s \\\n' % cnf.replace('cnf','por'))
            run.write('  @input %s > %s' % (imd, imd.replace('imd','omd')))
        subprocess.run(['chmod','+x',file],cwd = directory + name)
        if run:
            subprocess.run(['./%s' % file],cwd = directory + name)

    def add_ion(name, kind, top, cnf, charge, directory = prepDir, ff = ff, gromos = gromos, run = True, graphene = False):
        #Adds a counter-ion to coordinate files
        file = 'ion_%s.arg' % kind

        with open(directory + name + '/%s' % file, 'w') as arg:
            arg.write('@topo\t%s\n' % top )
            arg.write('@pbc\tr\n')
            arg.write('@%s\t1\t%s\n' % (ion_charge[charge], ion_arg[charge]))
            arg.write('@potential\t0.8\n')
            arg.write('@mindist\t0.35\n')
            arg.write('@pos\t%s' % cnf)
        if run:
            run_program('ion', file, Simulation.get_cnf(name, kind, ff, True, charge, graphene), name, directory, gromos)


    def make_graphene(name, kind, cnf, directory = prepDir, templates = templateDir, ff = ff, gromos = gromos, genbox = False, block = 'CCC'):
        #Makes a graphene layer based on solvent box size
        PDB = directory + name + "/graphene_%s_vmd.pdb" % kind
        #Coordinates
        if not genbox:
            genbox = Simulation.get_genbox(name, cnf, directory).split()[1:]
        lx = genbox[0]
        ly = genbox[1]
        lz = genbox[2]

        bond = '0.139'
        ydist = float(bond)
        xdist = 0.1205 #https://rechneronline.de/pi/hexagon.php
        newlx = '%1.9f' % (float(lx) + xdist)
        newly = '%1.9f' % (float(ly) + ydist)
        shape = 'armchair'
        command1 = 'graphene -lx %s -ly %s -type %s -cc %s' % (lx,ly,shape,bond)
        command2 = 'animate write pdb {%s} beg 0 end 0 skip 1 0' % (PDB)
        with open(directory + name + '/%s_graphene.tcl' % kind, 'w') as tcl:
            tcl.write(command1)
            tcl.write('\n')
            tcl.write(command2)
            tcl.write('\n')
            tcl.write('quit')

        subprocess.run(['vmd','-e', '%s_graphene.tcl' % kind],cwd = directory + name)
        #Basado en correo de Yerko pero actualizado a nuevos nombres en MDAnalysis
        u = MDAnalysis.Universe(PDB)
        sele = u.select_atoms("all")
        resnum = len(sele.residues)
        sele.residues.resnames = ['CCC'] * resnum
        for i in sele:
            i.name = "C"+str(i.ix%4 +1)
        r = 0
        for i in sele:
            if i.ix%4 +1 == 1:
                r+=1
            i.residue.resid = r
        sele.write(directory + name + "/graphene_%s.pdb" % kind)
        ##Modify in case of oxidation
        if block == 'CCCOH':
            with open(directory + name + "/graphene_%s.pdb" % kind) as file:
                lines = file.readlines()
            with open(directory + name + "/graphene_%s_notox.pdb" % kind,'w') as file:
                for line in lines:
                    file.write(line)
            with open(directory + name + "/graphene_%s.pdb" % kind,'w') as file:
                for line in lines[0:3]: #Header
                    file.write(line)
                count = 0
                original_count = 0
                for line in lines[3:-1]:
                    count += 1
                    original_count += 1
                    if line.split()[2] != 'C4':
                        newline = line[:11].replace(str(original_count),str(count)) + line[12:]
                        file.write(newline)
                    else:
                        ccc_line = line[:11].replace(str(original_count),str(count)) + line[12:]
                        file.write(ccc_line)

                        count += 1
                        o_line = line[:11].replace(str(original_count),str(count)) + line[12:].replace('C4',' O').replace(' 0.000','13.600')
                        file.write(o_line)

                        count += 1
                        h_line = line[:11].replace(str(original_count),str(count)) + line[12:].replace('C4',' H').replace(' 0.000','14.600')
                        file.write(h_line)
                file.write(lines[-1])
        #Topology
        file = 'make_top_grapehene.arg'
        with open(directory + name + '/%s' % file,'w') as arg:
            arg.write('#Graphene topology\n')
            arg.write('@build\t%s\n' % mtb)
            arg.write('@param\t%s\n' % ifp)
            arg.write('@seq\t%s\n' % block)
            arg.write('@solv\tH2O\n')

        run_program( 'make_top',file,'graphene_%s_1.top' % ff, name, directory, gromos)
        command3 = [gromos + 'programs/com_top','@topo','%i:' % r + 'graphene_%s_1.top' % ff,'@param','1','@solv','1','>', 'graphene_separate_%s_%s.top' % (ff, kind)]
        subprocess.run(' '.join(command3),shell = True, cwd = directory + name)
        command4 = [gromos + 'programs/pdb2g96', '@topo', 'graphene_separate_%s_%s.top' % (ff, kind),'@lib', templates + '/pdb2g96.lib', '@pdb', 'graphene_%s.pdb' % kind, '>', 'graphene_separate_%s.cnf' % kind]
    #    print(' '.join(command4))
        subprocess.run(' '.join(command4), shell = True, cwd = directory + name)
        if block == 'CCCOH':
            file = 'make_top_grapehene_notox.arg'
            with open(directory + name + '/%s' % file,'w') as arg:
                arg.write('#Graphene topology\n')
                arg.write('@build\t%s\n' % mtb)
                arg.write('@param\t%s\n' % ifp)
                arg.write('@seq\t%s\n' % 'CCC')
                arg.write('@solv\tH2O\n')
            run_program( 'make_top',file,'graphene_%s_1_notox.top' % ff, name, directory, gromos)
            subprocess.run(' '.join(command3).replace('.top','_notox.top'),shell = True, cwd = directory + name)
            subprocess.run(' '.join(command4).replace('.pdb','_notox.pdb').replace('.cnf','_notox.cnf').replace('.top','_notox.top'), shell = True, cwd = directory + name)
        last = ''
        with open(directory + name + '/graphene_separate_%s.cnf' % kind, 'r') as ifile:
            with open(directory + name + '/graphene_%s.cnf' % kind,'w') as ofile:
                position = False
                for line in ifile:
                    if 'POSITION' in line and 'REF' not in line:
                        newline = 'POSITION\n'
                        position = True
                    elif 'END' in line:
                        newline = 'END\n'
                        position = False
                    elif position:
                        split = re.split(r'(\s+)', line)
                        split[2] = ' ' * (len(split[2]) - 1) + '1'
                        # Adjustments to account for borders of box
                        split[10] = "%.9f" % (float(split[10]) + xdist*1.5)
                        split[12] = "%.9f" % (float(split[12]) + ydist/2)
                        newline = ''.join(split)
                        last = split[8]
                    elif 'GENBOX' in line:

                            ofile.write('GENBOX\n')
                            ofile.write('\t1\n')
                            ofile.write('\t%s\t%s\t%s\n' % (newlx, newly, lz))
                            ofile.write('\t90.000000000\t90.000000000\t90.000000000\n')
                            ofile.write('\t0.000000000\t0.000000000\t0.000000000\n')
                            ofile.write('\t0.000000000\t0.000000000\t0.000000000\n')
                            ofile.write('END\n')
                            break
                    else:
                        newline = line
                    ofile.write(newline)
        if block != 'CCC':
            with open(directory + name + '/graphene_separate_%s_notox.cnf' % kind, 'r') as ifile:
                with open(directory + name + '/graphene_%s_notox.cnf' % kind,'w') as ofile:
                    position = False
                    for line in ifile:
                        if 'POSITION' in line and 'REF' not in line:
                            newline = 'POSITION\n'
                            position = True
                        elif 'END' in line:
                            newline = 'END\n'
                            position = False
                        elif position:
                            split = re.split(r'(\s+)', line)
                            split[2] = ' ' * (len(split[2]) - 1) + '1'
                            # Adjustments to account for borders of box
                            split[10] = "%.9f" % (float(split[10]) + xdist*1.5)
                            split[12] = "%.9f" % (float(split[12]) + ydist/2)
                            newline = ''.join(split)
                            #last = split[8] PASS THE LAST ATOM OF THE OXIDIZED FILE
                        elif 'GENBOX' in line:

                                ofile.write('GENBOX\n')
                                ofile.write('\t1\n')
                                ofile.write('\t%s\t%s\t%s\n' % (newlx, newly, lz))
                                ofile.write('\t90.000000000\t90.000000000\t90.000000000\n')
                                ofile.write('\t0.000000000\t0.000000000\t0.000000000\n')
                                ofile.write('\t0.000000000\t0.000000000\t0.000000000\n')
                                ofile.write('END\n')
                                break
                        else:
                            newline = line
                        ofile.write(newline)
        return last

    def graphite(name, kind, top, cnf = False, directory = prepDir, charge = False, solvent = True, ff = ff, gromos = gromos, gromosJAG = gromosJAG, lastgraph = '1', oxidize = False):
        #runs graphite program to fix graphene topology

        if cnf == False:
            cnf = Simulation.get_cnf(name, kind, ff, solvent, charge, graphene = True)

        file = 'graphite_%s.arg' % kind

        if not oxidize:
            with open(directory + name + '/' + file,'w') as arg:
                arg.write('@topo\t%s\n' % top)
                arg.write('@pos\t%s\n' % cnf)
                arg.write('@pbc\tr')

            run_program( 'graphite',file,'graphene_graphite_%s_%s.top' % (ff, kind), name, directory, gromosJAG, 'contrib/')
            with open(directory + name + '/graphene_graphite_%s_%s.top' % (ff, kind), 'r') as ifile:
                with open(directory + name + '/graphene_%s_%s.top' % (ff, kind), 'w') as file:
                    temp_flag = False
                    pres_flag = False
                    for line in ifile:
                        if 'TEMPERATUREGROUPS' in line:
                            temp_flag = True
                            file.write('TEMPERATUREGROUPS\n')
                            file.write('# NSTM: number of temperature atom groups\n')
                            file.write('# NST[1...NSTM]: atom sequence number of last atom\n')
                            file.write('#                of the successive temperature atom groups\n')
                            file.write('#      NSTM  NST[1...NSTM]\n')
                            file.write('         1\n')
                            file.write('\t%s\n' % (lastgraph))
                            file.write('END\n')
                        if 'PRESSUREGROUPS' in line:
                            pres_flag = True
                            file.write('PRESSUREGROUPS\n')
                            file.write('# NSVM: number of pressure atom groups\n')
                            file.write('# NSV[1...NSVM]: atom sequence number of last atom\n')
                            file.write('#                of the successive pressure atom groups\n')
                            file.write('#      NSVM  NSV[1...NSVM]\n')
                            file.write('         1\n')
                            file.write('\t%s\n' % (lastgraph))
                            file.write('END\n')
                        if not temp_flag and not pres_flag:
                            file.write(line)
                        elif 'END' in line:
                            temp_flag = False
                            pres_flag = False
        else:
            def convert(x):
                return x + 2 * ((x-1)  // 4)
            def index_convert(y):
                if y % 7 == 5 or y % 6 == 0:
                    x = 4*((y-1)//6 + 1)
                else:
                    x =  y - 2*(y//6)
                return x - 1

            def parse_bond(line):
                i = int(line[:7])
                j = int(line[7:14])
                return line[:7].replace(str(i),str(convert(i))) + line[7:14].replace(str(j),str(convert(j))) + line[14:]
            def parse_angle(line):
                i = int(line[:7])
                j = int(line[7:14])
                k = int(line[14:21])
                return line[:7].replace(str(i),str(convert(i))) + line[7:14].replace(str(j),str(convert(j))) + line[14:21].replace(str(k),str(convert(k)))    + line[21:]
            def parse_dihedral(line):
                i = int(line[:7])
                j = int(line[7:14])
                k = int(line[14:21])
                l = int(line[21:28])
                return line[:7].replace(str(i),str(convert(i))) + line[7:14].replace(str(j),str(convert(j))) + line[14:21].replace(str(k),str(convert(k))) + line[21:28].replace(str(l),str(convert(l)))   + line[28:]
            def parse_exclusions(line):
                start = line[:2]
                end = line[2:]
                for i in end.split():
                    if i != '0':
                        old = ' ' + i + ' '
                        new = ' ' + str(convert(int(i))) + ' '
                        end.replace(old,new)
                return start + end

            with open(directory + name + '/' + file,'w') as arg:
                arg.write('@topo\t%s\n' % top.replace('.top','_notox.top'))
                arg.write('@pos\t%s\n' % cnf.replace('.cnf','_notox.cnf'))
                arg.write('@pbc\tr')
            run_program( 'graphite',file,'graphene_graphite_%s_%s.top' % (ff, kind), name, directory, gromosJAG, 'contrib/')

            with open(directory + name + '/graphene_graphite_%s_%s.top' % (ff, kind), 'r') as graphite:
                atoms_flag = False
                exclusion_flag = False
                bond_flag = False
                bond_angle_flag = False
                imp_dihedral_flag = False
                bonds = []
                angles = []
                dihedrals = []
                exclusions = []
                exclusion = ''
                for line in graphite:
                    if len(line.strip()) == 0:
                        continue
                    elif line.strip()[0] == '#':
                        continue
                    elif line == 'END\n':
                        bond_flag = bond_angle_flag = imp_dihedral_flag = atoms_flag = exclusion_flag = False
                        continue
                    elif atoms_flag and not exclusion_flag:
                        if len(line)>50:
                            exclusion_flag = True
                            exclusion = line[45:]
                        continue
                    elif atoms_flag and exclusion_flag:
                        if line.strip() == '0':
                            exclusion_flag = False
                            exclusion += line
                            exclusions.append(parse_exclusions(exclusion))
                        else:
                            exclusion += line
                        continue
                    elif bond_flag:
                        if len(line.strip()) > 5:
                            bonds.append(parse_bond(line))
                        continue
                    elif bond_angle_flag:
                        if len(line.strip()) > 5:
                            angles.append(parse_angle(line))
                        continue
                    elif imp_dihedral_flag:
                        if len(line.strip()) > 5:
                            dihedrals.append(parse_dihedral(line))
                        continue
                    elif line == 'BOND\n':
                        bond_flag = True
                        continue
                    elif line == 'BONDANGLE\n':
                        bond_angle_flag = True
                        continue
                    elif line == 'IMPDIHEDRAL\n':
                        imp_dihedral_flag = True
                        continue
                    elif line == 'SOLUTEATOM\n':
                        atoms_flag = True
                with open(directory + name + '/%s' % top) as oxide:
                    with open(directory + name + '/graphene_%s_%s.top' % (ff, kind), 'w') as file:
                        temp_flag = False
                        pres_flag = False
                        solute_flag = False
                        pasted_bonds = False
                        atomcount = 0
                        for line in oxide:
                            if len(line.strip()) == 0:
                                file.write(line)
                                continue
                            elif line.strip()[0] == '#':
                                file.write(line)
                                continue
                            elif 'TEMPERATUREGROUPS' in line:
                                temp_flag = True
                                file.write('TEMPERATUREGROUPS\n')
                                file.write('# NSTM: number of temperature atom groups\n')
                                file.write('# NST[1...NSTM]: atom sequence number of last atom\n')
                                file.write('#                of the successive temperature atom groups\n')
                                file.write('#      NSTM  NST[1...NSTM]\n')
                                file.write('         1\n')
                                file.write('\t%s\n' % (lastgraph))
                                file.write('END\n')
                                continue
                            elif 'PRESSUREGROUPS' in line:
                                pres_flag = True
                                file.write('PRESSUREGROUPS\n')
                                file.write('# NSVM: number of pressure atom groups\n')
                                file.write('# NSV[1...NSVM]: atom sequence number of last atom\n')
                                file.write('#                of the successive pressure atom groups\n')
                                file.write('#      NSVM  NSV[1...NSVM]\n')
                                file.write('         1\n')
                                file.write('\t%s\n' % (lastgraph))
                                file.write('END\n')
                                continue
                            elif 'SOLUTEMOLECULES' in line:
                                pres_flag = True
                                file.write('SOLUTEMOLECULES\n')
                                file.write('# NSPM: number of separate molecules in solute block\n')
                                file.write('# NSP[1...NSPM]: atom sequence number of last atom\n')
                                file.write('#                of the successive submolecules\n')
                                file.write('#      NSPM  NSP[1...NSPM]\n')
                                file.write('         1\n')
                                file.write('\t%s\n' % (lastgraph))
                                file.write('END\n')
                                continue

                            elif 'END' in line:
                                temp_flag = False
                                pres_flag = False
                                solute_flag = False
                                bond_flag = bond_angle_flag = imp_dihedral_flag = atoms_flag = exclusion_flag = pasted_bonds = False
                            elif line == 'SOLUTEATOM\n':
                                atoms_flag = True
                            elif atoms_flag:
                                if len(line.strip()) < 10:
                                    continue
                                else:
                                    atomcount += 1
                                    try:
                                        newline = line[:45] + exclusions[index_convert(atomcount)]
                                    except IndexError:
                                        print(line)
                                        print(exclusions)
                                        print(index_convert(atomcount))
                                        exit()
                                    file.write(newline)
                            elif line == 'BOND\n':
                                bond_flag = True
                                continue
                            elif bond_flag:
                                if len(line.strip()) > 5:
                                    if not pasted_bonds:
                                        for bond in bonds:
                                            file.write(bond)
                                        pasted_bonds = True
                                        continue
                                    else:
                                        if line.strip()[-2:] == '13':
                                            file.write(line)
                                            continue
                                        else:
                                            continue
                            elif line == 'BONDANGLE\n':
                                bond_angle_flag = True
                                continue
                            elif bond_angle_flag:
                                if len(line.strip()) > 5:
                                    if not pasted_bonds:
                                        for angle in angles:
                                            file.write(angle)
                                        pasted_bonds = True
                                        continue
                                    else:
                                        if line.strip()[-1] == '1':
                                            file.write(line)
                                            continue
                                        else:
                                            continue
                            elif line == 'IMPDIHEDRAL\n':
                                imp_dihedral_flag = True
                                continue
                            elif imp_dihedral_flag:
                                if len(line.strip()) > 5:
                                    if not pasted_bonds:
                                        for imp in dihedrals:
                                            file.write(imp)
                                        pasted_bonds = True
                                        continue
                                    else:
                                        continue
                            elif not temp_flag and not pres_flag and not solute_flag and not bond_flag and not imp_dihedral_flag and not bond_angle_flag and not atoms_flag:
                                file.write(line)
                                continue


    def make_box(name, kind, cnf = False, directory = prepDir, charge = False, ff = ff, gromos = gromos, rotate = True, below = '1.500000000', SSD0 = '0.300000000'):
        if not cnf:
            cnf = Simulation.get_cnf(name, kind, ff=ff, charge = charge)
        #Creates a box with specified size and places graphene and peptide at specific positions and angles
        graphene_pos = []
        graphene_lattice = []
        graphene_genbox = []
        solute_pos = []
        solute_lattice = []
        solute_timestep = []
        #Read graphene
        with open(directory + name + '/graphene_%s_min.cnf' % kind, 'r') as file:
            pos = False
            lattice = False
            genbox = False
            for line in file:
                if 'POSITION' in line:
                    pos = True
                elif 'LATTICESHIFTS' in line:
                    lattice = True
                elif 'GENBOX' in line:
                    genbox = True
                elif 'END' in line:
                    pos = False
                    lattice = False
                    genbox = False
                elif pos:
                    graphene_pos.append(line)
                elif lattice:
                    graphene_lattice.append(line)
                elif genbox:
                    graphene_genbox.append(line)

        #Read Solute
        with open(directory + name + '/' + cnf,'r') as file:
            pos = False
            lattice = False
            timestep = False
            for line in file:
                if 'POSITION' in line:
                    pos = True
                elif 'LATTICESHIFTS' in line:
                    lattice = True
                elif 'TIMESTEP' in line:
                    timestep = True
                elif 'END' in line:
                    pos = False
                    lattice = False
                    timestep = False
                elif pos:
                    solute_pos.append(line)
                elif lattice:
                    solute_lattice.append(line)
                elif timestep:
                    solute_timestep.append(line)

        #Read pos as vectors
        last_sol = int(solute_pos[-1].split()[3])
        solute_pos_moved = []
        coords = np.empty((last_sol, 3))
        for counter, i in enumerate(solute_pos[1:]): #Take coords as vectors
            x, y, z =i.split()[4:7]
            fx, fy, fz = float(x), float(y), float(z)
            vector = np.array([fx, fy, fz])
            coords[counter] = vector
        #Orientation
        #Calculate max distance
        max_dist = 0
        A = 0
        B = 0
        for a, i in enumerate(coords): #Find most distant atoms
            for b, j in enumerate(coords):
                distance = np.linalg.norm(i - j)
                if distance > max_dist:
                    max_dist = distance
                    A = a
                    B = b
        orientation_vector = coords[A] - coords[B]
        Simulation.frameout(name, kind, 'TEMP/', cnf , move = False)
        PDB = '%s_moved.pdb' % kind
        subprocess.run(['mv','FRAME_00001.pdb',PDB], cwd = directory + name)
        with open(directory + name + '/%s_rotation.tcl' % kind, 'w') as tcl:
            x , y, z = orientation_vector
            tcl.write('set sel [atomselect top all]\n')
            tcl.write('set M [transvecinv {%f %f %f}]\n' % (x, y, z))
            tcl.write('$sel move $M\n')
            tcl.write('set M [transaxis y -180]\n')
            tcl.write('$sel move $M\n')
            tcl.write('animate write pdb {%s} beg 0 end 0 skip 1 0\n' % (PDB))
            tcl.write('quit')
        subprocess.run(['vmd',PDB,'-e', '%s_rotation.tcl' % kind],cwd = directory + name)
        new_coords = np.empty((last_sol, 3))
        pdb_coords = []
        with open(directory + name + '/' + PDB,'r') as pdb:
            for line in pdb:
                if 'ATOM' in line:
                    pdb_coords.append(line[30:56])
        for counter, i in enumerate(pdb_coords):
            x, y , z = i.split()
            x, y, z= float(x)/10, float(y)/10, float(z)/10
            vector = np.array([x,y,z])
            new_coords[counter] = vector
        #Center
        center_x = new_coords[:,0].sum()/last_sol
        center_y = new_coords[:,1].sum()/last_sol
        center_z = new_coords[:,2].sum()/last_sol
        center_solute = np.array((center_x, center_y, center_z))
    #    new_center = np.array([2.5,2.5,float(below) + float(SSD0)])
        new_center = np.array([2.5,2.5,2.5])
        displacement = new_center - center_solute
        new_coords += displacement
        #Find min
    #    minZ = 5
    #    for i in new_coords:
    #        if i[2] < minZ:
    #            minZ = i[2]
    #    displacement = np.array([0,0,float(below) + float(SSD0)] - minZ)
    #    new_coords += displacement


        for counter, i in enumerate(solute_pos[1:]):
            fx, fy, fz = new_coords[counter]
            line = re.split(r'(\s+)', i)
            if len(line[10]) > 11 and fx >0:
                line[9] += ' '
            elif len(line[10]) == 11 and fx < 0:
                line[9] = line[9][:-1]
            if len(line[12]) > 11 and fx >0:
                line[11] += ' '
            elif len(line[12]) == 11 and fx < 0:
                line[11] = line[11][:-1]
            if len(line[14]) > 11 and fx >0:
                line[13] += ' '
            elif len(line[14]) == 11 and fx < 0:
                line[13] = line[13][:-1]
            line[10] = '%1.9f' % fx
            line[12] = '%1.9f' % fy
            line[14] = '%1.9f' % fz
            solute_pos_moved.append(''.join(line))
        #Move and Modify Graphene
        graphene_pos_moved = []
        for i in graphene_pos[1:]:
            line = re.split(r'(\s+)', i)
            line[14] = below
            line[1] = line[1][:-1]
            line[2] = '10'
            new_number = str(int(line[8]) + last_sol)
            len_diff = len(new_number) - len(line[8])
            line[7] = ' ' * (len(line[7]) - len_diff)
            line[8] = new_number
            graphene_pos_moved.append(''.join(line))


        #Write
        with open(directory + name + '/%s' % Simulation.get_cnf(name, kind, ff, False, charge, graphene = True),'w') as file:
            file.write('TITLE\n\tCombined coordinates of solute and graphene\nEND\n')
            file.write('POSITION\n')
            for i in solute_pos_moved:
                file.write(i)
            for i in graphene_pos_moved:
                file.write(i)
            file.write('END\nLATTICESHIFTS\n')
            for i in solute_lattice:
                file.write(i)
            for i in graphene_lattice:
                file.write(i)
            file.write('END\nGENBOX\n')
            for i in graphene_genbox:
                file.write(i)
            file.write('END\n')

    def graphene_com_top(name, kind, top = False, directory = prepDir, charge = False, solvent = False, frameoutDir = frameoutDir, ff = ff, gromos = gromos):
        #Topology
        if not top:
            top = Simulation.get_top(name,kind,ff,charge)
        file = 'com_top_%s_graphene.arg' % kind
        with open(directory + name + '/' + file, 'w') as arg:
            arg.write('@topo\t%s graphene_%s_%s.top\n' % (top, ff, kind))
            arg.write('@param\t1\n')
            arg.write('@solv\t1\n')
        run_program( 'com_top', file, 'graphene_%s_com_top.top' % kind, name, directory, gromos)

    #    graphene_count = count_atoms(name, kind, get_cnf(name, kind, ff ,solvent, charge, True), directory, ff, charge, solvent, gromos, True, top = 'graphene_%s_com_top.top' % kind)[0][1]

        with open(directory + name + '/graphene_%s_com_top.top' % kind,'r') as template:
            with open(directory + name + '/' + Simulation.get_top(name, kind, ff, charge, True),'w') as file:
                temp_flag = False
                pres_flag = False
                atoms, mols = Simulation.count_atoms(name, kind, Simulation.get_cnf(name, kind, ff, solvent, charge, graphene = True),directory, frameoutDir, ff, solvent, gromos, graphene = True, top = 'graphene_%s_com_top.top' % kind)
                lastsol = atoms[0]
                lastgraph = atoms[1]

                for line in template:
                    if 'TEMPERATUREGROUPS' in line:
                        temp_flag = True
                        file.write('TEMPERATUREGROUPS\n')
                        file.write('# NSTM: number of temperature atom groups\n')
                        file.write('# NST[1...NSTM]: atom sequence number of last atom\n')
                        file.write('#                of the successive temperature atom groups\n')
                        file.write('#      NSTM  NST[1...NSTM]\n')
                        file.write('         2\n')
                        file.write('\t%s\t%s\n' % (lastsol, lastgraph))
                        file.write('END\n')
                    if 'PRESSUREGROUPS' in line:
                        pres_flag = True
                        file.write('PRESSUREGROUPS\n')
                        file.write('# NSVM: number of pressure atom groups\n')
                        file.write('# NSV[1...NSVM]: atom sequence number of last atom\n')
                        file.write('#                of the successive pressure atom groups\n')
                        file.write('#      NSVM  NSV[1...NSVM]\n')
                        file.write('         2\n')
                        file.write('\t%s\t%s\n' % (lastsol, lastgraph))
                        file.write('END\n')
                    if not temp_flag and not pres_flag:
                        file.write(line)
                    elif 'END' in line:
                        temp_flag = False
                        pres_flag = False



    ###Equilibration
    def equilibrate(name, kind, top, cnf, directory = prepDir, gromos = gromos,  run = True, graphene = False, restrain_graphene = False):
        #runs equilibration Simulation for the solvent box
    #    #Copy jobs template
    #    if graphene:
    #        origin = templateDir + 'equilibration_graphene.jobs'
    #    else:
    #        origin = templateDir + 'equilibration.jobs'
    #    destination = directory + name + '/equilibration.jobs'
    #    subprocess.run(['cp',origin,destination],cwd=directory)
        #Use mk_script program
        file = 'eq_mk_script_%s.arg' % kind
        with open(directory + name + '/%s' % file,'w') as arg:
            arg.write('@sys\teq_%s\n' % kind)
            arg.write('@bin\t%s\n' % md)
            arg.write('@dir\t%s\n' % (directory + name))
            arg.write('@files\n')
            arg.write('  topo\t%s\n' %  top)
            arg.write('  input\t%s_equilibration.imd\n' % kind)
            arg.write('  coord\t%s\n' %  cnf)
            arg.write('  posresspec\t%s\n' %  cnf.replace('cnf','por'))
            arg.write('  refpos\t%s\n' %  cnf.replace('cnf','rpr'))
            arg.write('@template\t%smk_script.lib\n' % templateDir)
            arg.write('@version\tmd++\n')
            if not graphene:
                arg.write('@joblist\t%sequilibration.jobs' % templateDir)
            if graphene:
                arg.write('@joblist\t%sequilibration_graphene.jobs' % templateDir)
        if run:
            run_program( 'mk_script',file, file.replace('arg','out'), name, directory = directory, gromos = gromos)
            if restrain_graphene:
                    subprocess.run(['mv','eq_%s_5.run'%kind,'eq_%s_5.back'%kind],cwd = directory + name)
                    with open(directory + name + '/eq_%s_5.back'%kind,'r') as back:
                        with open(directory + name + '/eq_%s_5.run'%kind,'w') as run:
                            for line in back:
                                if 'POSRESSPEC=${SIMULDIR}/' in line:
                                    run.write('POSRESSPEC=${SIMULDIR}/eq_%s_4.por' % kind)
                                else:
                                    run.write(line)
                    subprocess.run(['chmod','+x','eq_%s_5.run' % kind],cwd = directory + name)
                    subprocess.run(['mv','eq_%s_4.run'%kind,'eq_%s_4.back'%kind],cwd = directory + name)
                    with open(directory + name + '/eq_%s_4.back'%kind,'r') as back:
                        with open(directory + name + '/eq_%s_4.run'%kind,'w') as run:
                            for line in back:
                                if '# perform last command (usually submit next job)' in line:
                                    break
                                else:
                                    run.write(line)
                    subprocess.run(['chmod','+x','eq_%s_4.run' % kind],cwd = directory + name)
            subprocess.run(['./eq_%s_1.run' % kind],shell = True, cwd = directory + name)
            if restrain_graphene:
                Simulation.restrain(name, kind, 'eq_%s_4.cnf' % kind, directory, graphene = True, restrain_graphene=True)
                subprocess.run(['./eq_%s_5.run' % kind],shell = True, cwd = directory + name)
    #Run Simulation


    ###Prep for Leftraru
    def copy_files(name, kind, eq_cnf, cnf, top, inputDir = prepDir, outputDir = outputDir):
        #Copies cnf and topo files to the output folder
        subprocess.run(['cp',eq_cnf,outputDir + name + '/' + cnf],cwd = inputDir + name)
        subprocess.run(['cp',top,outputDir + name + '/' + top],cwd = inputDir + name)

    def leftraru_mk_script(name, kind, cnf, top, leftraruDir = leftraruDir, outputDir = outputDir, md = md_leftraru, pepUser = False, aaUser = False, graphene = False, restrain_graphene = False, local_elevation = False, umbrella_sampling = False, stop = 100, job = 'md', r = ''):
        if pepUser and kind != 'aa':
            leftraruDir = leftraruDir.replace('username',pepUser)
            md = md.replace('username',pepUser)
        if aaUser and kind == 'aa':
            leftraruDir = leftraruDir.replace('username',aaUser)
            md = md.replace('username',aaUser)
        #Prepares .arg file for
        if not umbrella_sampling:
            with open('%s/%s/%s/%s_mk_script.arg' % (outputDir,name, kind,job) ,'w') as arg:
                    arg.write('@sys\t%s_%s\n' % (job,kind) )
                    arg.write('@bin\t%s\n' % md)
                    arg.write('@version\tmd++\n')
                    arg.write('@dir\t%s%s/%s\n' % (leftraruDir, name, kind))
                    arg.write('@files\n')
                    arg.write('  topo\t../%s\n' % top)
                    arg.write('  input\t%s_%s.imd\n' % (kind,job))
                    arg.write('  coord\t../%s\n' % cnf)
                    if restrain_graphene:
                        arg.write('  posresspec\t../%s\n' %  cnf.replace('cnf','por'))
                        arg.write('  refpos\t../%s\n' %  cnf.replace('cnf','rpr'))
                    if local_elevation:
                        arg.write('  ledih\t../%s\n' % cnf.replace(".cnf",".led"))
                        arg.write('  leumb\t../%s\n' % cnf.replace(".cnf",".lud"))
                    arg.write('@template\t%smkscript_mms_leftraru_openMP.lib\n' % leftraruDir)
                    arg.write('@script\t1 %i' % stop)
        else:
            with open('%s/%s/%s/r%s/%s_mk_script.arg' % (outputDir,name, kind,('%1.3f' % r).replace('.','_'),job) ,'w') as arg:
                    arg.write('@sys\t%s_%s\n' % (job,kind) )
                    arg.write('@bin\t%s\n' % md)
                    arg.write('@version\tmd++\n')
                    arg.write('@dir\t%s%s/%s/r%s\n' % (leftraruDir, name, kind,('%1.3f' % r).replace('.','_')))
                    arg.write('@files\n')
                    arg.write('  topo\t../../%s\n' % top)
                    arg.write('  input\t%s_%s.imd\n' % (kind,job))
                    arg.write('  coord\t../../%s\n' % cnf)
                    if restrain_graphene:
                        arg.write('  posresspec\t../../%s\n' %  cnf.replace('cnf','por'))
                        arg.write('  refpos\t../../%s\n' %  cnf.replace('cnf','rpr'))
                    if local_elevation:
                        arg.write('  ledih\t../../%s\n' % cnf.replace(".cnf",".led"))
                        arg.write('  leumb\t../../%s\n' % cnf.replace(".cnf",".lud"))
                    arg.write('  disres\t../%s\n' % cnf.replace(".cnf","_%s.dsr" % ('%1.3f' % r).replace('.','_')))
                    arg.write('@template\t%smkscript_mms_leftraru_openMP.lib\n' % leftraruDir)
                    arg.write('@script\t1 %i' % stop)

    def leftraru_jobs(name, kind, templateDir = templateDir, leftraruDir = leftraruDir, outputDir = outputDir, cores = 10, pepUser = False, aaUser = False, graphene = False, job = 'md', r = False, stop = 100):
        #Prepares .job files for leftraru
        if pepUser and kind != 'aa':
            leftraruDir = leftraruDir.replace('username',pepUser)
        if aaUser and kind == 'aa':
            leftraruDir = leftraruDir.replace('username',aaUser)

        if not r:
            jobsFile = outputDir + name + '/' + kind + '/nlhpc.jobs'
            submitFile = outputDir + name + '/' + kind + '/submit.txt'
        else:
            jobsFile = outputDir + name + '/' + kind + '/r%s/nlhpc.jobs' % ('%1.3f' % r).replace('.','_')
            submitFile = outputDir + name + '/' + kind + '/r%s/submit.txt' % ('%1.3f' % r).replace('.','_')

        with open(templateDir + 'nlhpc.jobs' ,'r') as ifile:
            template = ifile.readlines()
        with open(jobsFile,'w') as ofile:
            jobname = '%s_%s_h2o_100NS_%ic_%s' % (name, kind, cores, job)
            if graphene:
                jobname.replace('h2o_','h2o_graph_')
            for line in template:
                newline = line.replace('TI_L2', jobname)
                newline = newline.replace('totcores','%i' % cores)
                newline = newline.replace('1-100%1','1-%i%%1' % stop)
                if cores >= 20:
                    newline = newline.replace('cores', '20')
                else:
                    newline = newline.replace('cores', '%i' % cores)
                ofile.write(newline)
        with open(submitFile,'w') as ofile:
            for j in range(1,stop + 1):
                if r:
                    ofile.write('%s%s/%s/r%s/%s_%s_%i.run\n' % (leftraruDir, name,kind,('%1.3f' % r).replace('.','_'),job,kind,j) )
                else:
                    ofile.write('%s%s/%s/%s_%s_%i.run\n' % (leftraruDir, name,kind,job,kind,j) )



    def prepare(self):
        self.make_folders(self.prepDir)
        self.make_folders(self.outputDir, separate = True)
        self.fix_pdb()
        if self.charge: self.make_top_ion()
        self.make_top()
        self.pdb2g96()
        self.gch()
        self.frameout('gch/',self.cnf_gch)
        self.make_imd(self.cnf_gch, 'minimization.imd', self.prepDir,TSENERG = '10')
        self.min_run(self.top, self.cnf_gch, '%s_minimization.imd' % self.kind)
        self.frameout('min/',self.cnf_min)
        if self.graphene:
            #Prepare
            self.lastgraph = self.make_graphene()
            self.graphite(self.lastgraph)
            self.fix_graphene_top()

            self.graphene_top = 'graphene_%s_%s.top' % (self.ff, self.kind)
            self.graphene_cnf = 'graphene_%s.cnf' % (self.kind)
            self.graphene_cnf_min = 'graphene_%s_min.cnf' % (self.kind)
            self.graphite_top = 'graphene_graphite_%s_%s.top' % (self.ff, self.kind)

            self.frameout('graphene/',self.graphene_cnf, top = self.graphite_top, gather = True)
            self.make_imd(self.graphene_cnf, 'minimization_graphene.imd', self.prepDir, TSENERG='10', top = self.graphite_top)
            self.min_run(self.graphene_top, self.graphene_cnf, '%s_minimization_graphene.imd' % self.kind)
            self.frameout('graph_min/', self.graphene_cnf_min, top = self.graphite_top)

            #Solvate
            self.make_box()
            self.graphene_com_top()
            self.frameout('combined/',self.cnf_graph_charged, gather = True, top = self.top_graph_charged)
            self.sim_box(self.top_graph, self.cnf_graph, genbox = self.genbox)
        else:
            self.sim_box(self.top, self.cnf_min)
        self.sim_box_cnf = 'sim_box_%s.cnf' % self.kind
        self.frameout('/simbox',self.sim_box_cnf, gather = True, top = self.top_graph)
        self.restrain(self.sim_box_cnf)
        self.make_imd(self.sim_box_cnf, 'solvent.imd',solvent = True, TSENERG = '10', graphene=self.graphene, top = self.top_graph)
        self.solv_min(self.top_graph, self.sim_box_cnf, '%s_solvent.imd' % self.kind, self.graphene)

        if self.charge:
            self.com_top_ion()
            self.add_ion()

        #Equilibrate Solvent Box
        self.frameout('simbox_min/', self.cnf_graph_charged_solvent, gather = True, top = self.top_graph_charged)
        self.make_imd(self.cnf_graph_charged_solvent, 'equilibration.imd', solvent = True, graphene = self.graphene, charge = self.charge, TSTRA = '100', TSENERG='100',TAU='0.01')
        self.restrain(self.cnf_graph_charged_solvent)
        self.equilibrate(self.top_graph_charged,self.cnf_graph_charged_solvent)
        self.frameout('thermalized/','eq_%s_5.cnf' % self.kind, gather = True, top = self.top_graph_charged)

        #LeftraruPrep
        self.copy_files()
        self.fix_graphene_top(self.name, self.top_graph_charged, self.cnf_graph_charged_solvent, self.outputDir)
        if self.restrain_graphene: self.restrain(self.cnf_graph_charged_solvent, self.outputDir, True)

        if self.local_elevation:
            self.local_elevation_file()
            self.make_imd(self.cnf_graph_charged_solvent, 'LE.imd', self.outputDir, self.outputDir, TAU = '0.01', solvent = True, sidechain = self.sidechain, graphene = self.graphene, charge = self.charge, top = self.top_graph_charged, add_restraint = self.restrain_graphene)
            subprocess.run(['mv',self.outputDir + self.name + '/%s_LE.imd' % self.kind,self.outputDir + self.name + '/%s/LE.imd' % self.kind])
            self.leftraru_mk_script(job = 'LE')
            self.leftraru_jobs(job = 'LE')

        elif self.umbrella_sampling:
            for r in np.arange(0.3,1.50001,0.05):
                subdir = self.outputDir + '%s/%s/r%s/' % (self.name,self.kind,('%1.3f' % r).replace('.','_'))
                os.makedirs(subdir, exist_ok = True)
                if self.LINCS:
                    self.us_imd = 'US_LINCS.imd'
                else:
                    self.us_imd = 'US.imd'
                self.make_imd(self.cnf_graph_charged_solvent, self.us_imd, self.outputDir, self.outputDir, TAU = '0.01', solvent = True, sidechain = self.sidechain, graphene = self.graphene, charge = self.charge, top = self.top_graph_charged, add_restraint = self.restrain_graphene, r = r)
                self.leftraru_mk_script(job = self.us_imd, r = r)
                self.leftraru_jobs(job = self.us_imd, r = r)
                self.umbrella_sampling_file(r)



        else:
            if not self.LINCS:
                self.make_imd( self.cnf_graph_charged_solvent, 'md.imd', self.outputDir, self.outputDir, TAU = '0.01', solvent = True, sidechain = self.sidechain, graphene = self.graphene, charge = self.charge, top = self.top_graph_charged, add_restraint = self.restrain_graphene, filename = '%s/%s_md.imd' % (self.kind,self.kind))
            else:
                self.make_imd(self.cnf_graph_charged_solvent, 'md_graphene(links).imd', self.outputDir, self.outputDir, TAU = '0.01', solvent = True, sidechain = self.sidechain, graphene = self.graphene, charge = self.charge, top = self.top_graph_charged, add_restraint = self.restrain_graphene, filename = '%s/%s_md.imd' % (self.kind,self.kind))
            subprocess.run(['mv',self.outputDir + self.name + '/%s_md.imd' % self.kind,self.outputDir + self.name + '/%s/md.imd' % self.kind])
            self.leftraru_mk_script()
            self.leftraru_jobs()
