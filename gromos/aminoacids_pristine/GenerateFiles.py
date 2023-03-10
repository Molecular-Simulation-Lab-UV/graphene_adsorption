#!/usr/bin/env python 
from SimulationClass import Simulation

import csv, os
import random
import numpy as np

##Forcefield name and location
ff = '54a8'
ifp = os.getcwd() + '/FF/54a8.ifp'
mtb = os.getcwd() + '/FF/54a8_CCC.mtb'
##Other Variables
seed = '893230'
if len(seed) != 6:
    seed = '%i' % random.randint(100000,999999)
genbox = ['4.935500000','4.865000000','5.000000000']
##Directories
prepDir = os.getcwd() + '/Prep/'
frameoutDir = os.getcwd() + '/FRAMEOUTS/'
templateDir = os.getcwd() + '/Templates/'
outputDir = os.getcwd() + '/Simulation/'
leftraruDir = '/home/username/US_aa/'
gromos = '/home/mbarria/gromos++-1.4.0/x86/'
gromosJAG = '/home/mbarria/gromos++-1.2.0_JAG/x86/'
md_leftraru = '/home/username/opt/md++-1.4.0/bin/md'
md = '/home/mbarria/md++/bin/md'
vmd = '/home/mbarria/Peptidos/vmd_pdb/' #Coordinates need to be supplied, I use VMD to create them.
##Ion dictionaries to ensure prober counter ions
ion_file = {'+e':'Cl','-e':'Na'}
ion_arg = {'+e':'CL-','-e':'NA+'}
ion_charge = {'+e':'negative','-e':'positive'}

w0_dict = {0.3:10,
           0.35:10,
           0.4:10,
           0.45:10,
           0.5:10,
           0.55:10,
           0.6:10,
           0.65:10,
           0.7:10,
           0.75:10,
           0.8:5,
           0.85:5,
           0.9:5,
           0.95:5,
           1.0:5,
           1.05:5,
           1.1:5,
           1.15:5,
           1.2:5,
           1.25:5,
           1.3:5,
           1.35:5,
           1.4:5,
           1.45:5,
           1.5:5,
           }

if __name__ == '__main__':
    with open('Lista.csv', newline = '') as csvfile:
        print('Using seed: ',seed)
        csvfile.readline()
        file = csv.reader(csvfile, delimiter=',')        
        for row in file:   
            name = "_".join(row[0].split()[0:-1])
            if name == 'Tryptophan': continue
            print('Starting in %s' % name)
            sequence_dic = {'pep':row[3].split(),'aa':row[5].split()}
            charge = row[1]
            for kind in ['aa']:
                sequence = sequence_dic[kind]
                system = Simulation(name, kind, sequence, charge, ff = ff, ifp = ifp, mtb = mtb, seed=seed, graphene = True, 
					prepDir = prepDir, outputDir = outputDir, templateDir = templateDir, 
					frameoutDir = frameoutDir, leftraruDir = leftraruDir, gromos = gromos, gromosJAG = gromosJAG, vmd = vmd, md = md, 
					md_leftraru = md_leftraru, local_elevation = False, umbrella_sampling = True, restrain_graphene = True, 
					genbox = genbox, aaUser = 'rhonorato', pepUser = 'rhonorato', sidechain = True, w0_dict = w0_dict, stop = 25)
                system.prepare()