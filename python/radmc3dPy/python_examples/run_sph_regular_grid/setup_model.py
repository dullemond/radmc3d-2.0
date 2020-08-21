import sphtool
import numpy as np
import radmc3dPy.natconst as nc
import matplotlib.pyplot as plt
import radmc3dPy
import os

os.system('cp -v ../datafiles/dustkappa_silicate.inp .')

# Filename of the SPH output
sph_fname = 'snap_00000.ascii'
# SPH unit length scale
unit_length = 5.0*nc.au
# SPH unit mass scale
unit_mass   = nc.ms
# radmc3d model name to use to generate the grid and all other input files except for the density
model = 'ppdisk'
# Grid boundaries in the radial direction
xbound ='[1.0*au, 150.*au]'
# Nr of grid points in the radial direction
nx = '[100]'
# Grid boundaries in the polar direction
ybound = '[0., pi]'
# Nr of grid points in the polar direction
ny = '[100]'
# Grid boundaries in the azimuthal direction
zbound = '[0., 2.0*pi]'
# Nr of grid points in the azimuthal direction
nz = '[100]'
# Binary output files for radmc3d variables
binary = True

#
# Generate the RADMC-3D model (except for the dust density!)
#
model = radmc3dPy.setup.radmc3dModel(model=model, xbound=xbound, ybound=ybound, zbound=zbound, nx=nx, ny=ny, nz=nz, binary=binary)
model.setupDust(ddens=False)

#
# Read the SPH data
#
d = np.genfromtxt(sph_fname, skip_header=12)

# Select the SPH particles
ii = (d[:,-1] == 1)
x = d[ii,0]
y = d[ii,1]
z = d[ii,2]
pmass = d[ii,3]
h = d[ii,4]
rho = d[ii,5]
vx = d[ii,6]
vy = d[ii,7]
vz = d[ii,8]

#
# Scale the simulation
#
x *= unit_length
y *= unit_length
z *= unit_length
h *= unit_length
#
# Note that the SPH simulation contains the gas density, thus we need to scale it with the dust-to-gas ratio in order to 
# generate the dust density required by RADMC-3D.
#
rho *= unit_mass * model.par.ppar['dusttogas'] / unit_length**3
pmass *= unit_mass * model.par.ppar['dusttogas']  
vx *= np.sqrt(nc.gg*unit_mass/unit_length)
vy *= np.sqrt(nc.gg*unit_mass/unit_length)
vz *= np.sqrt(nc.gg*unit_mass/unit_length)

vel = np.zeros([x.shape[0], 1, 3], dtype=np.float64)
vel[:,0,0] = vx
vel[:,0,1] = vy 
vel[:,0,2] = vz

#
# Get the gridder tool
#
tool = sphtool.SPHTool(maxSPHParticlesPerCell=10, maxSPHTreeDepth=20, nThreads=4)
#
# Set the sph simulation data
#
tool.setSPHData(x=x, y=y, z=z, h=h, rho=rho, pmass=pmass, vectors=vel)
#
# Initialise the gridding
#
tool.init()
# 
# Now the final thing is to pass the generated read to the gridder tool and re-interpolate the SPH simulation to this grid 
#
tool.regridToRegularGrid(rgrid=model.grid, scalarFloor=[1e-30])

#
# Write the generated dust density to file
#
tool.writeScalar(fname="dust_density.binp", ivar=0, nscal=1, binary=True)   
#
# Run the thermal monte carlo simulation to calculate the dust temperature
#
os.system('radmc3d mctherm')
#
# Calculate an image at 1.3mm
#
radmc3dPy.image.makeImage(npix=400, sizeau=200., incl=45., phi=0., wav=1300.)

#
# Read in the image and display it
#
im = radmc3dPy.image.readImage()
radmc3dPy.image.plotImage(im, au=True, log=True, maxlog=3)

dum = input()


