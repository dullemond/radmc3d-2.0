import matplotlib.pyplot as plt
import numpy as np
import os
from radmc3d_tools import bplanck
#
# Some natural constants
#
au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
MS  = 1.98892e33     # Solar mass              [g]
TS  = 5.78e3         # Solar temperature       [K]
LS  = 3.8525e33      # Solar luminosity        [erg/s]
RS  = 6.96e10        # Solar radius            [cm]
cc  = 2.9979245800000e10      # Light speed             [cm/s]
pi  = 3.1415926535897932385e0 
#
# Monte Carlo parameters
#
nphot    = 100000
#
# Grid parameters
#
nx       = 1
ny       = 1
nz       = 200
sizex    = 1e90
sizey    = 1e90
sizez    = 10*au
#
# Model parameters
#
rho0     = 1e-16
#
# Illuminating parallel flux beams
#
#   *** Uncomment one of the options below ***
#
#  Test 1: No illuminating beams from the top.
#
inclillum = [0.]
tempillum = [0.]
dilillum  = [0.]
#
#  Test 2: Vertical incident beam of T=100 K. We should see that
#          the temperature drops a bit at the edge, because even
#          though the flux is consistent with a blackbody of 100 K,
#          it is not the same as being in a thermal bath of 100 K.
#
#          You can see this as follows: A grey dust grain of radius
#          a in a parallel beam of flux of F_nu = pi * B_nu(100K) 
#          will absorb int_0^infty pi*B_nu(100K)*pi*a^2 dnu erg/s.
#          It will emit 4*pi*a^2*int_0^infty pi*B_nu(T) dnu erg/s.
#          This leads to T=(1/4)^(1/4)*100K=70K. In our example
#          the flux from this beam only needs to account for half
#          of the radiation (the other half comes from below). But
#          even if you have two of these beams, you still get just
#          84 K. This means that if you put the sigma_SB * T^4 
#          flux in a parallel beam, you will not get the same
#          effect as a thermal boundary. 
# 
#inclillum = [0.]
#tempillum = [100.]
#dilillum  = [1.e0]
#
#  Test 3: According to the Rybicki & Lightman book, the two-stream
#          approximation is a decent approximation of radiative 
#          transfer in plane-parallel geometries. The two rays have
#          angle cos(theta) = 1/sqrt(3) or cos(theta) = -1/sqrt(3).
#          According to this logic, if we make the incident beam
#          under such an angle (which is 54.735 degrees), and if
#          we boost the unprojected flux such that the projected
#          flux is still equal to the thermal flux of 100 K, then
#          this should be as good as a thermal boundary. 
#
#inclillum = [54.735e0]
#tempillum = [100.]
#dilillum  = [1.e0*sqrt(3.)]
#
#  Test 4: Take a much stronger flux that the 100 K flux in the
#          previous test cases. Put the inclination angle 0, i.e.
#          the radiation goes vertically downward into the atmosphere.
#          You see that the temperature on the top of the atmosphere
#          is clearly higher than the base, but at the very top the
#          temperature drops a bit again. That drop is because due to
#          the incidence being vertical, the energy is injected deep
#          into the atmosphere, where, due to random direction
#          re-emission, it gets a bit "stuck", and thus heats the
#          interior more than the surface.
#
#inclillum = [0.]
#tempillum = [400.]
#dilillum  = [1.e0]
#
#  Test 5: Same as Test 4, but now with grazing incidence angle. Here
#          you see the opposite happening: the very surface layer has
#          a temperature increase. This is because of the fact that 
#          dust grains at the surface see the full blast of the flux
#          while dust grains in the interior only "feel" the projected
#          flux.
#
#inclillum = [80.]
#tempillum = [400.]
#dilillum  = [1.e0]
#
#----------------------------------------------------------------------
#
# Count the illumination beams
#
nillum    = len(inclillum)
if tempillum[0] == 0.0: nillum=0
#
# Make the coordinates
#
xi       = np.array([-1e90,1e90])
yi       = np.array([-1e90,1e90])
zi       = np.linspace(0,sizez,nz+1)
zc       = 0.5 * ( zi[0:nz] + zi[1:nz+1] )
#
# Make the dust density model
#
rhod     = rho0 + np.zeros(nz)
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
# Write the grid file
#
with open('amr_grid.inp','w+') as f:
    f.write('1\n')                       # iformat
    f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
    f.write('10\n')                      # Coordinate system
    f.write('0\n')                       # gridinfo
    f.write('0 0 1\n')                   # Include only z coordinate
    f.write('%d %d %d\n'%(1,1,nz))       # Size of grid
    for value in xi:
        f.write('%13.6e\n'%(value))      # X coordinates (cell walls)
    for value in yi:
        f.write('%13.6e\n'%(value))      # Y coordinates (cell walls)
    for value in zi:
        f.write('%13.6e\n'%(value))      # Z coordinates (cell walls)
#
# Write the density file
#
with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nz))                 # Nr of cells
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
# Add an illuminating beam
#
if(nillum>0):
    with open('illum.inp','w') as f:
        f.write('2\n')
        f.write('{0} {1}\n'.format(nillum,nlam))
        for illum in range(nillum):
            f.write('{0} 0.0\n'.format(inclillum[illum]))
        f.write('\n')
        for ilam in range(nlam):
            f.write('{}\n'.format(lam[ilam]))
        for illum in range(nillum):
            f.write('\n')
            for ilam in range(nlam):
                nu  = 1e4*cc/lam[ilam]
                bpl = bplanck(nu,tempillum[illum])
                f.write('{}\n'.format(dilillum[illum]*pi*bpl[0]))
else:
    os.system('rm -f illum.inp')
#
# Write the radmc3d.inp control file
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot = %d\n'%(nphot))
    f.write('scattering_mode_max = 0\n')   # Put this to 1 for isotropic scattering
    f.write('iranfreqmode = 1\n')
    f.write('thermal_boundary_zl = 100.\n')
