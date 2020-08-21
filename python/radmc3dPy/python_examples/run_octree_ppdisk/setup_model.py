import radmc3dPy
import os
import matplotlib.pyplot as plt
# -----------------------------------------------------------------------------
# Dust model
# -----------------------------------------------------------------------------

# Write the parameter file with the default parameters
radmc3dPy.analyze.writeDefaultParfile('ppdisk_amr')

# Dust model setup with formatted ascii input files
radmc3dPy.setup.problemSetupDust(model='ppdisk_amr', nsample=30, threshold=0.9, binary=False)

exit()
# Copy the dust opacity 
os.system('cp -v ../datafiles/dustkappa_silicate.inp .')

# Calculate the dust temperature
os.system('radmc3d mctherm')

# Calculate a continuum image at 1.3mm
os.system('radmc3d image npix 400 sizeau 250 incl 45. lambda 1300')
os.system('cp image.out image_cont_1300um.out')

#
# Visualisation
#

# Read the dust density to make some plots
d = radmc3dPy.analyze.readData(ddens=True, octree=True, binary=False)

# Plot the density distribution in a vertical slice of the model in the xz plane.
fig = plt.figure()
radmc3dPy.analyze.plotSlice2D(data=d, plane='xz', var='ddens', showgrid=True, linunit='au', gridalpha=0.1,
       xlim=(5., 100.), ylim=(-50., 50.), log=True, vmin=1e-25, nx=100, ny=100, nproc=3)


# Plot the density distribution in a horizontal slice of the model in the xy plane.
fig = plt.figure()
radmc3dPy.analyze.plotSlice2D(data=d, plane='xy', var='ddens', showgrid=True, linunit='au', gridalpha=0.1,
        xlim=(-100., 100.), ylim=(-100., 100.), log=True, vmin=1e-22, nx=100, ny=100, nproc=3)


# Plot the calculated continuum image
im = radmc3dPy.image.readImage()
fig = plt.figure()
radmc3dPy.image.plotImage(im, au=True)
dum = raw_input()
# -----------------------------------------------------------------------------
# Gas model
# -----------------------------------------------------------------------------
# Dust model setup with formatted ascii input files
radmc3dPy.setup.problemSetupGas(model='ppdisk_amr', binary=False)

# Copy the dust opacity and co data files from the datafiles directory
os.system('cp -v ../datafiles/molecule_co.inp .')

# Calculate a CO channel map at the systemic velocity (v=0km/s)
os.system('radmc3d image npix 400 sizeau 250 incl 45. iline 3 vkms 0.')
os.system('cp image.out image_co.out')

# Plot the calculated continuum image
im = radmc3dPy.image.readImage()
fig = plt.figure()
radmc3dPy.image.plotImage(im, au=True)


dum = raw_input()
