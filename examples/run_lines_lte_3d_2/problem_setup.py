#
# Import NumPy for array handling
#
import numpy as np
from radmc3dPy.natconst import *
pi  = 3.1415926535897932385e0
#
# Import plotting libraries (start Python with ipython --matplotlib)
#
#from mpl_toolkits.mplot3d import axes3d
#from matplotlib import pyplot as plt
#
# Monte Carlo parameters
#
nphot    = 1000000
#
# Grid parameters
#
nx       = 64
ny       = 64
nz       = 64
sizex    = 100*au
sizey    = 100*au
sizez    = 100*au
#
# Model parameters
#
omega    = 2*pi/(50*year)
rhogas0  = 1e-16
temp0    = 50.e0
dusttogas= 0.01
vturb0   = 3.*1e5
#
# Star parameters
#
mstar    = ms
rstar    = rs
tstar    = ts
pstar    = np.array([0.,0.,0.])
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
# Make grid with size interpreted as half-width size here
#
xi       = np.linspace(-sizex,sizex,nx+1)
yi       = np.linspace(-sizey,sizey,ny+1)
zi       = np.linspace(-sizez,sizez,nz+1)
xc       = 0.5 * ( xi[0:nx] + xi[1:nx+1] )
yc       = 0.5 * ( yi[0:ny] + yi[1:ny+1] )
zc       = 0.5 * ( zi[0:nz] + zi[1:nz+1] )
qq       = np.meshgrid(xc,yc,zc,indexing='ij')
xx       = qq[0]
yy       = qq[1]
zz       = qq[2]
rrcyl    = np.sqrt(xx**2+yy**2)
#
# Make a simple solid-body rotating gas flow
#
rhogas  = np.zeros((nx,ny,nz)) + rhogas0
rhod    = rhogas*dusttogas
tgas    = np.zeros((nx,ny,nz)) + temp0
vcir    = rrcyl * omega
vx      =  vcir * yy / ( rrcyl + 1e-10*sizex )
vy      = -vcir * xx / ( rrcyl + 1e-10*sizey )
vz      = np.zeros((nx,ny,nz))
vturb   = np.zeros((nx,ny,nz)) + vturb0
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
    for value in xi:
        f.write('%13.6e\n'%(value))      # X coordinates (cell walls)
    for value in yi:
        f.write('%13.6e\n'%(value))      # Y coordinates (cell walls)
    for value in zi:
        f.write('%13.6e\n'%(value))      # Z coordinates (cell walls)
#
# Write the dust density file. Here we use a dust-to-gas ratio of 0.01
#
with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    f.write('1\n')                       # Nr of dust species
    data = rhod.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
#
# Write the molecule number density file. 
#
abunco = 1e-4
factco = abunco/(2.3*mp)
nco    = rhogas*factco
with open('numberdens_co.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    data = nco.ravel(order='F')          # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
#
# Write the gas velocity field
#
with open('gas_velocity.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                f.write('%13.6e %13.6e %13.6e\n'%(vx[ix,iy,iz],vy[ix,iy,iz],vz[ix,iy,iz]))
#
# Write the microturbulence file
#
with open('microturbulence.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    data = vturb.ravel(order='F')          # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
#
# Write the gas temperature
#
with open('gas_temperature.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    data = tgas.ravel(order='F')          # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
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
# Write the lines.inp control file
#
with open('lines.inp','w') as f:
    f.write('1\n')
    f.write('1\n')
    f.write('co    leiden    0    0\n')
#
# Write the radmc3d.inp control file
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot = %d\n'%(nphot))
    f.write('scattering_mode_max = 0\n')   # Put this to 1 for isotropic scattering
    f.write('tgas_eq_tdust   = 1')
