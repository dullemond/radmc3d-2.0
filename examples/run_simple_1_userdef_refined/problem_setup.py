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
radius   = 5*au
rho0     = 1e-16
temp0    = 20.e0                # Put this to 0 if you dont want temp to be set
tnoise   = 0.e0                 # Noise in temperature = 0
#tnoise   = 10.e0               # To clearly see the AMR refinement, make noise
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
# Write the radmc3d.inp control file
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot = %d\n'%(nphot))
    f.write('scattering_mode_max = 0\n')
    f.write('iranfreqmode = 1\n')
    f.write('userdef_nx = %d\n'%(nx))
    f.write('userdef_ny = %d\n'%(ny))
    f.write('userdef_nz = %d\n'%(nz))
    f.write('userdef_sizex = %13.6e\n'%(sizex))
    f.write('userdef_sizey = %13.6e\n'%(sizey))
    f.write('userdef_sizez = %13.6e\n'%(sizez))
    f.write('userdef_radius = %13.6e\n'%(radius))
    f.write('userdef_rho0 = %13.6e\n'%(rho0))
    f.write('userdef_levelmax = %d\n'%(10))
    f.write('userdef_nrefinefact = %d\n'%(4.0))
    f.write('userdef_amr_relcellsize = %13.6e\n'%(0.05))
    f.write('userdef_amr_refregion = %13.6e\n'%(0.5))
    f.write('userdef_temp0 = %13.6e\n'%(temp0))
    f.write('userdef_tempnoise = %13.6e\n'%(tnoise))

