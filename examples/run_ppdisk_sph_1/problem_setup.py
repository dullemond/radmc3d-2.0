#
# Model setup for mapping an SPH (Smooth Particle Hydrodynamics) model
# into RADMC-3D. NOTE: The example here does not use a real SPH model,
# but creates a dummy SPH model.
#
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import default_rng
from radmc3d_tools.sph_to_sphergrid import sph_to_sphergrid
import matplotlib.pyplot as plt
import os
#
# Some natural constants
#
au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]
ss  = 5.6703e-5      # Stefan-Boltzmann const  [erg/cm^2/K^4/s]
kk  = 1.3807e-16     # Bolzmann's constant     [erg/K]
mp  = 1.6726e-24     # Mass of proton          [g]
GG  = 6.67408e-08    # Gravitational constant  [cm^3/g/s^2]
LS  = 3.8525e33      # Solar luminosity        [erg/s]
RS  = 6.96e10        # Solar radius            [cm]
MS  = 1.98892e33     # Solar mass              [g]
TS  = 5.78e3         # Solar temperature       [K]
pi  = np.pi          # Pi

#
# Star parameters
#
mstar    = 2.4*ms
rstar    = 2.4*rs
tstar    = 1e4
pstar    = np.array([0.,0.,0.])

#
# Parameters of the dummy model
#
mdisk      = 1e-3*MS      # Disk gas mass
rin        = 10*au        # Inner disk radius
rout       = 100*au       # Outer disk radius
r0         = 10*au        # r0 for use with H_p/r
hpr0       = 0.05         # H_p/r at r=r0
plh        = 0.2          # Flaring
dtg        = 0.01         # Dust-to-gas ratio

#
# Parameter of the smoothing kernel
#
hrkernel0  = 0.03
#hrkernel   = hrkernel0*np.array([1,1,1])       # The smoothing kernel (spherical)
hrkernel   = hrkernel0*np.array([2,1,2])      # The smoothing kernel (triaxial)

#
# Dummy model for the SPH particles
#
rng        = default_rng()
nsph       = 1000000
#sph_r      = rin*(rout/rin)**(rng.uniform(size=nsph))
sph_r      = rin+(rout-rin)*(rng.uniform(size=nsph))
hpr        = hpr0*(sph_r/r0)**plh
sph_th     = np.pi/2 - hpr*rng.standard_normal(size=nsph)
sph_phi    = 2*np.pi*rng.uniform(size=nsph)
sph_lr     = np.log(sph_r)
rthphi     = np.vstack((sph_r,sph_th,sph_phi)).T
masses     = dtg*mdisk/nsph

#
# Parameters of the grid
#
grid_nr    = 100
grid_nth   = 100
grid_nph   = 100
grid_hrup  = np.abs(np.pi/2-sph_th).max()+hrkernel[1]

#
# Map the SPH particles onto the grid
#
sphgrid    = sph_to_sphergrid(rthphi=rthphi,masses=masses,nr=grid_nr, \
                              ntheta=grid_nth,nphi=grid_nph,hrup=grid_hrup,hrkernel=hrkernel)
rhod       = sphgrid.rho
ri         = sphgrid.grid_ri
thetai     = sphgrid.grid_thetai
phii       = sphgrid.grid_phii
nr         = sphgrid.grid_nr
ntheta     = sphgrid.grid_ntheta
nphi       = sphgrid.grid_nphi

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
# Monte Carlo parameters
#
#nphot_therm    = 1000000000
nphot_therm    = 100000000
nphot_scat     = 10000000

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
    f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
    f.write('100\n')                     # Coordinate system: spherical
    f.write('0\n')                       # gridinfo
    f.write('1 1 1\n')                   # Include r,theta,phi coordinates
    f.write('%d %d %d\n'%(nr,ntheta,nphi))  # Size of grid
    for value in ri:
        f.write('%13.6e\n'%(value))      # X coordinates (cell walls)
    for value in thetai:
        f.write('%17.10e\n'%(value))     # Y coordinates (cell walls) (use higher precision here)
    for value in phii:
        f.write('%13.6e\n'%(value))      # Z coordinates (cell walls)
#
# Write the density file
#
with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nr*ntheta*nphi))     # Nr of cells
    f.write('1\n')                       # Nr of dust species
    data = rhod.ravel(order='F')   # Create a 1-D view, fortran-style indexing
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
    f.write('nphot_therm = %d\n'%(nphot_therm))
    f.write('nphot_scat = %d\n'%(nphot_scat))
    f.write('scattering_mode_max = 1\n')        # Only isotropic scattering
    f.write('iranfreqmode = 1\n')
#
# Show the mapping in r,theta, where we average over phi. Shown is rho*r^p
# where you can choose p. 
#
p    = 2.
rc   = np.sqrt(ri[1:]*ri[:-1])
thc  = 0.5*(thetai[1:]+thetai[:-1])
lrc  = np.log10(rc/au)
plt.figure()
plt.imshow((rhod.mean(axis=2)*rc[:,None]**p).T,origin='lower',extent=[lrc.min(),lrc.max(),thc.min(),thc.max()])
plt.xlabel(r'$^{10}\log(r/\mathrm{au})$')
plt.ylabel(r'$\pi/2-\theta$')
plt.figure()
plt.imshow((rhod.sum(axis=1)*rc[:,None]**p).T,origin='lower',extent=[lrc.min(),lrc.max(),0,2*np.pi],aspect=0.1)
plt.xlabel(r'$^{10}\log(r/\mathrm{au})$')
plt.ylabel(r'$\phi$')
plt.show()
