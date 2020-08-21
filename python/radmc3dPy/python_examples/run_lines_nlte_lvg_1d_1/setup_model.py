# Import the radmc3dPy module
import radmc3dPy
import os

# Write the parameter file with the default parameters
radmc3dPy.analyze.writeDefaultParfile('lines_nlte_lvg_1d_1')

# Dust model setup with ascii input files
radmc3dPy.setup.problemSetupDust('lines_nlte_lvg_1d_1', binary=False, write_dusttemp=True)

# Gas model setup with ascii input files
radmc3dPy.setup.problemSetupGas('lines_nlte_lvg_1d_1', binary=False, write_gastemp=True)

# Copy the dust opacity and co data files from the datafiles directory
os.system('cp -v ../datafiles/dustkappa_silicate.inp .')
os.system('cp -v ../datafiles/molecule_co.inp .')
