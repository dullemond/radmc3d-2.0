#
# Import NumPy for array handling
#
import numpy as np
#
# Import the Voronoi grid generator tool
#
from radmc3d_tools.radmc3d_voronoi_grid import *
#
# Import plotting libraries (start Python with ipython --matplotlib)
#
#from mpl_toolkits.mplot3d import axes3d
#from matplotlib import pyplot as plt
#
# Some natural constants
#
au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]
#
# Monte Carlo parameters
#
nphot    = 1000000
#
# Grid parameters
#
npt      = 30000       # Nr of grid points
npttry   = npt*4       # Start with more than necessary, select those that are inside sphere
size     = 10*au       # Radius of the sphere of grid points
#
# Model parameters
#
radius   = 5*au
rho0     = 1e-16
#
# Star parameters
#
mstar    = ms
rstar    = rs
tstar    = ts
pstar    = np.array([0.,0.,0.])
#
# Make the coordinates
#
# Make a sphere of random points using rejection sampling
#
x        = (2*np.random.random(npttry)-1)*size
y        = (2*np.random.random(npttry)-1)*size
z        = (2*np.random.random(npttry)-1)*size
rr       = np.sqrt(x**2+y**2+z**2)
mask     = rr<=size
x        = x[mask]
y        = y[mask]
z        = z[mask]
assert len(x)>=npt, 'Increase npttry'
x        = x[:npt]
y        = y[:npt]
z        = z[:npt]
rr       = np.sqrt(x**2+y**2+z**2)
cellpts  = np.vstack([x,y,z]).T
#
# Call the Voronoi grid generator
#
grid     = Voronoigrid(cellpts)
#
# Make the dust density model
#
rhod     = rho0 * np.exp(-(rr**2/radius**2)/2.0)
#
# Make sure that the density is zero in the "open"
# (unbounded) cells
#
rhod[grid.cell_iopen] = 0.0
#
# Write the wavelength_micron.inp file
#
lam1     = 0.1e0
lam2     = 7.0e0
lam3     = 25.e0
lam4     = 1.0e4
n12      = 20
n23      = 100
n34      = 30
lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
lam      = np.concatenate([lam12,lam23,lam34])
nlam     = lam.size
#
# Write the wavelength file
#
with open('wavelength_micron.inp','w+') as f:
    f.write('%d\n'%(nlam))
    for value in lam:
        f.write('%13.6e\n'%(value))
#
#
# Write the stars.inp file
#
with open('stars.inp','w+') as f:
    f.write('2\n')
    f.write('1 %d\n\n'%(nlam))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
    for value in lam:
        f.write('%13.6e\n'%(value))
    f.write('\n%13.6e\n'%(-tstar))
#
# Write the grid file
#
grid.write_radmc3d_unstr_grid(bin=False)
#
# Write the density file
#
with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(grid.ncells))        # Nr of cells
    f.write('1\n')                       # Nr of dust species
    data = rhod.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
#
# Dust opacity control file
#
with open('dustopac.inp','w+') as f:
    f.write('2               Format number of this file\n')
    f.write('1               Nr of dust species\n')
    f.write('============================================================================\n')
    f.write('1               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('silicate        Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')
#
# Write the radmc3d.inp control file
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot = %d\n'%(nphot))
    f.write('scattering_mode_max = 0\n')   # Put this to 1 for isotropic scattering
    f.write('iranfreqmode = 1\n')

