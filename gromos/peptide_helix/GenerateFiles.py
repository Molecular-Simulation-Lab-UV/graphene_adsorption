#!/usr/bin/env python
from SimulationClass import Simulation

import csv
import os
import random

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
vmd = '/home/mbarria/Peptidos/vmd_pdb_penta/'
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

if __name__ == '__main__':
    with open('Lista.csv', newline='') as csvfile:
        print('Using seed: ', seed)
        csvfile.readline()
        file = csv.reader(csvfile, delimiter=',')
        user = ''
        for row in file:
            user = 'mbarria' #if user == 'rhonorato' else 'rhonorato'
#            Alternate between users
            name = "_".join(row[0].split()[0:-1])
            if name != 'Alanine':
                continue
            print('Starting in %s' % name)
            sequence_dic = {'pep': row[3].split(), 'aa': row[5].split()}
            charge = row[1]
            for kind in ['pep']:
                sequence = sequence_dic[kind]
                system = Simulation(name, kind, sequence, charge, ff=ff,
                                    ifp=ifp, mtb=mtb, seed=seed,
                                    graphene=True, prepDir=prepDir,
                                    outputDir=outputDir,
                                    templateDir=templateDir,
                                    frameoutDir=frameoutDir,
                                    leftraruDir=leftraruDir, gromos=gromos,
                                    gromosJAG=gromosJAG, vmd=vmd, md=md,
                                    md_leftraru=md_leftraru,
                                    local_elevation=False,
                                    umbrella_sampling=True,
                                    restrain_graphene=True, genbox=genbox,
                                    aaUser=user, pepUser=user,
                                    sidechain=True, w0_dict=w0_dict, stop=25,
                                    LINCS=True, block='CCC')
                system.prepare()
