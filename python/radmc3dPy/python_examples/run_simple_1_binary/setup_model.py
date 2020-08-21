# Import the radmc3dPy module
import radmc3dPy
import os 

# Write the parameter file with the default parameters
radmc3dPy.analyze.writeDefaultParfile('simple_1')

# Dust model setup with binary input files
radmc3dPy.setup.problemSetupDust('simple_1', binary=True)

# Copy the dust opacity and co data files from the datafiles directory
os.system('cp -v ../datafiles/dustkappa_silicate.inp .')
os.system('cp -v ../datafiles/molecule_co.inp .')
