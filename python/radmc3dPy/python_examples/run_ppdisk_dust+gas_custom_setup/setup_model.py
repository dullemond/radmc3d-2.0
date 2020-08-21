# Import the radmc3dPy module
import radmc3dPy
import os
import radmc3dPy.natconst as nc

# This model shows how to use a modular/custom setup

# ------------------------------------------------------------------------------------------------------------
# Dust setup
# ------------------------------------------------------------------------------------------------------------
# Get an instance of a radmc3dModel class using the ppdisk model. The 'model' keyword specifies the model
#  to be used to calculate the physical variables (dust/gas density, velocity, etc). The boolean 'binary' keyword 
#  can also be used if we wish the data files to be written in binary instead of formatted ascii format. 
#  One can also pass any parameter in the problem_params.inp file as keyword argument to the radmc3dModel 
#  contructor to update its value. By default if any keyword argument is passed the problem_params.inp file
#  will be overwritten with the new parameters. To prevent the overwriting of the parameterfile use the 
#  parfile_update=False keyword.
model = radmc3dPy.setup.radmc3dModel(model='ppdisk', xbound='[5.0*au,5.05*au, 100.0*au]', nx='[30, 200]', 
        ybound='[0., pi/3., pi/2., 2.*pi/3., pi]', ny='[10,60,60,10]', mdisk='1e-4*ms', binary=False)

# Generate the spatial and frequency grid and write them to file
model.makeGrid(writeToFile=True)

# Generate the radiation sources and write them to file
model.makeRadSources(writeToFile=True)

# Generate the dust opacity and write it to file
model.makeDustOpac()

# Generate the dust density structure and write it to file
model.makeVar(ddens=True, writeToFile=True)

# Write the radmc3d.inp file with all code parameters
model.writeRadmc3dInp()


# Copy the dust opacity 
os.system('cp -v ../datafiles/dustkappa_silicate.inp .')

# ------------------------------------------------------------------------------------------------------------
# Run thermal Monte Carlo simulation to calculate the dust temperature
# ------------------------------------------------------------------------------------------------------------

# Calculate the dust temperature
os.system('radmc3d mctherm')

# ------------------------------------------------------------------------------------------------------------
# Gas setup
# ------------------------------------------------------------------------------------------------------------

# Generate all gas input files (gas number density, microturbulent velocity field, macroscopic gas velocity) and
#  write all of them to files
model.makeVar(gdens=True, gvel=True, vturb=True, writeToFile=True)

# Write lines.inp for the line RT control variables
model.writeLinesInp()

# Write the radmc3d.inp file with all code parameters
model.writeRadmc3dInp()

# Copy the dust opacity and co data files from the datafiles directory
os.system('cp -v ../datafiles/molecule_co.inp .')
