# CALVADOS integration
Replace the [components.py, sim.py and interactions.py] in your calvados package (located in site-packages in your conda environment) with those provided here  

Then you can prepare a calvados run as normal by typing   

`python prepare.py --name [NameOfRunHERE]`  

Note that the prepare .py file provided here requires the following to be set  

`dihdomains =  f'{cwd}/input/dih_domains.yaml'` this specifies the residues for which you wish to apply the dih restraints  

`fresidues_epsds = f'{cwd}/input/epsd.csv'`  the epsd values used  

`use_dih = True`  # set this to true to use the restraints

An example run for DNAJB6's JD-GF is provided



