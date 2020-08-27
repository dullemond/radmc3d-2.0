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
# The layer patch
#
laynlev  = 1
laynr    = 1
lix      = nx//4+1
liy      = ny//4+1
liz      = nz//4+1
lnx      = nx//2
lny      = ny//2
lnz      = nz//2
lnnx     = lnx*2
lnny     = lny*2
lnnz     = lnz*2
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
xi       = np.linspace(-sizex,sizex,nx+1)
yi       = np.linspace(-sizey,sizey,ny+1)
zi       = np.linspace(-sizez,sizez,nz+1)
xc       = 0.5 * ( xi[0:nx] + xi[1:nx+1] )
yc       = 0.5 * ( yi[0:ny] + yi[1:ny+1] )
zc       = 0.5 * ( zi[0:nz] + zi[1:nz+1] )
#
# Make the coordinates of the layer patch of refinement
#
xi_l1    = np.linspace(-sizex/2,sizex/2,lnnx+1)
yi_l1    = np.linspace(-sizey/2,sizey/2,lnny+1)
zi_l1    = np.linspace(-sizez/2,sizez/2,lnnz+1)
xc_l1    = 0.5 * ( xi_l1[:-1] + xi_l1[1:] )
yc_l1    = 0.5 * ( yi_l1[:-1] + yi_l1[1:] )
zc_l1    = 0.5 * ( zi_l1[:-1] + zi_l1[1:] )
#
# Make the dust density model
#
qq       = np.meshgrid(xc,yc,zc,indexing='ij')
xx       = qq[0]
yy       = qq[1]
zz       = qq[2]
rr       = np.sqrt(xx**2+yy**2+zz**2)
rhod     = rho0 * np.exp(-(rr**2/radius**2)/2.0)
#
# Make the dust density model in the layer
#
qq_l1    = np.meshgrid(xc_l1,yc_l1,zc_l1,indexing='ij')
xx_l1    = qq_l1[0]
yy_l1    = qq_l1[1]
zz_l1    = qq_l1[2]
rr_l1    = np.sqrt(xx_l1**2+yy_l1**2+zz_l1**2)
rhod_l1  = rho0 * np.exp(-(rr_l1**2/radius**2)/2.0)
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
with open('amr_grid.inp','w+') as f:
    f.write('1\n')                       # iformat
    f.write('10\n')                      # AMR grid style  (10=layer-style AMR)
    f.write('0\n')                       # Coordinate system
    f.write('0\n')                       # gridinfo
    f.write('1 1 1\n')                   # Include x,y,z coordinate
    f.write('%d %d %d\n'%(nx,ny,nz))     # Size of grid
    f.write('%d %d\n'%(laynlev,laynr))   # Layers: nr of levels, nr of layers
    for value in xi:
        f.write('%13.6e\n'%(value))      # X coordinates (cell walls)
    for value in yi:
        f.write('%13.6e\n'%(value))      # Y coordinates (cell walls)
    for value in zi:
        f.write('%13.6e\n'%(value))      # Z coordinates (cell walls)
    f.write('%d %d %d %d %d %d %d\n'%(0,lix,liy,liz,lnx,lny,lnz))     # Info for layer
#
# Write the density file
#
with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz+lnnx*lnny*lnnz)) # Nr of spatial data points (incl redundant ones)
    f.write('1\n')                       # Nr of dust species
    data = rhod.ravel(order='F')         # Create a 1-D view of rhod
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
    data = rhod_l1.ravel(order='F')      # Create a 1-D view of rhod_l1
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

