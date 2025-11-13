import os
import pandas as pd
from calvados.cfg import Config, Job, Components
import subprocess
import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--name',nargs='?',required=True,type=str)
args = parser.parse_args()

cwd = os.getcwd()
sysname = f'{args.name:s}'

# set the side length of the cubic box
L = 20

# set the saving interval (number of integration steps)
N_save = 1000

# set final number of frames to save
N_frames = 1000

residues_file = f'{cwd}/input/residues.csv'

N_aa=97
L = int(np.ceil((N_aa - 1) * 0.38 + 4))


config = Config(
  # GENERAL
  sysname = sysname, # name of simulation system
  box = [L, L, L], # nm
  temp = 298, # K
  ionic = 0.15, # molar
  pH = 7.0, # 7.5
  topol = 'center',
  # RUNTIME SETTINGS
  wfreq = N_save, # dcd writing interval, 1 = 10 fs
  steps = 2*N_frames*N_save, # number of simulation steps
  runtime = 0, # overwrites 'steps' keyword if > 0
  platform = 'CPU', # CPU only for DIH
  restart = 'checkpoint',
  frestart = 'restart.chk',
  verbose = True,
)

# PATH
path = f'{cwd}/{sysname:s}'
subprocess.run(f'mkdir -p {path}',shell=True)

analyses = f"""

from calvados.analysis import save_rg

save_rg("{path:s}","{sysname:s}","{residues_file:s}",".",10)
"""

config.write(path,name='config.yaml',analyses=analyses)

components = Components(
  # Defaults
  molecule_type = 'protein',
  nmol = 1, # number of molecules
  restraint = True, # apply restraints
  charge_termini = 'both', # charge N or C or both  
  
  # INPUT
  fresidues = residues_file, # residue definitions
  ffasta = fastafile, # residue definitions
  pdb_folder = f'{cwd}/input', # directory for pdb and PAE files
  fdomains = f'{cwd}/input/domains.yaml', #only used if harmonic
  
  #DIH restraints
  dihdomains =  f'{cwd}/input/dih_domains.yaml', 
  fresidues_epsds = f'{cwd}/input/epsd.csv',
  use_dih = True, #apply dih restrains
  
  # RESTRAINTS
  restraint_type = 'harmonic', # harmonic or go go uses the pae json file harmonic the fdomain file
  use_com = False, # apply on centers of mass instead of CA
  #colabfold = 0, # PAE format (EBI AF=0, Colabfold=1&2)
  #k_go = 10., # Restraint force constant
)

components.add(name=args.name)

components.write(path,name='components.yaml')

