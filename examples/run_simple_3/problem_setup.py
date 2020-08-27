#
# Import NumPy for array handling
#
import numpy as np
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
nx       = 32
ny       = 32
nz       = 32
sizex    = 10*au
sizey    = 10*au
sizez    = 10*au
#
# Model parameters
#
radius   = 1*au
rho0     = 1e-13
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
xi       = np.linspace(-sizex,sizex,nx+1)
yi       = np.linspace(-sizey,sizey,ny+1)
zi       = np.linspace(-sizez,sizez,nz+1)
xc       = 0.5 * ( xi[0:nx] + xi[1:nx+1] )
yc       = 0.5 * ( yi[0:ny] + yi[1:ny+1] )
zc       = 0.5 * ( zi[0:nz] + zi[1:nz+1] )
#
# Make the grid
#
qq       = np.meshgrid(xc,yc,zc,indexing='ij')
xx       = qq[0]
yy       = qq[1]
zz       = qq[2]
#
# Now put a couple of blobs on the grid
#
pos      = np.array([[3,2,7],[-1,-6,-3],[-7,6,2]])*au
rhod     = np.zeros((nx,ny,nz))
for p in pos:
    rr       = np.sqrt((xx-p[0])**2+(yy-p[1])**2+(zz-p[2])**2)
    rhod     += rho0 * np.exp(-(rr**2/radius**2)/2.0)
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
    np.savetxt(f,lam.T,fmt=['%13.6e'])
#
#
# Write the stars.inp file
#
with open('stars.inp','w+') as f:
    f.write('2\n')
    f.write('1 %d\n\n'%(nlam))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
    np.savetxt(f,lam.T,fmt=['%13.6e'])
    f.write('\n%13.6e\n'%(-tstar))
#
# Write the grid file
#
with open('amr_grid.inp','w+') as f:
    f.write('1\n')                       # iformat
    f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
    f.write('0\n')                       # Coordinate system
    f.write('0\n')                       # gridinfo
    f.write('1 1 1\n')                   # Include x,y,z coordinate
    f.write('%d %d %d\n'%(nx,ny,nz))     # Size of grid
    np.savetxt(f,xi.T,fmt=['%13.6e'])    # X coordinates (cell walls)
    np.savetxt(f,yi.T,fmt=['%13.6e'])    # Y coordinates (cell walls)
    np.savetxt(f,zi.T,fmt=['%13.6e'])    # Z coordinates (cell walls)
#
# Write the density file
#
with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    f.write('1\n')                       # Nr of dust species
    data = rhod.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    np.savetxt(f,data.T,fmt=['%13.6e'])  # The data
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
    f.write('scattering_mode_max = 0\n')
    f.write('iranfreqmode = 1\n')

