import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from radmc3dPy.image import *
from radmc3dPy.analyze import *
from radmc3dPy.natconst import *
import os

#
# Make sure to have done the following beforhand:
#
#  First compile RADMC-3D
#  Then run:
#   python problem_setup.py
#   radmc3d mctherm
#

wav    = 0.2      # Wavelength in microns
nthrd  = 4        # Use 4 OpenMP parallel threads
nphot  = 100000
#nphot  = 10000   # Set to this small nr of photons if you want to get quick-n-dirty results

#
# Make and plot image of full disk at 1 microns: scattered light
#
# ... First isotropic scattering
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot_scat = {}\n'.format(nphot))
    f.write('scattering_mode_max = 1\n')
    f.write('iranfreqmode = 1\n')
makeImage(npix=200,incl=60.,phi=30.,wav=wav,sizeau=200,setthreads=nthrd)   # This calls radmc3d
os.system('cp image.out image_iso.out')

#
# ... Then with the full scattering
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot_scat = {}\n'.format(nphot))
    f.write('scattering_mode_max = 10\n')
    f.write('iranfreqmode = 1\n')
makeImage(npix=200,incl=60.,phi=30.,wav=wav,sizeau=200,setthreads=nthrd)   # This calls radmc3d 
os.system('cp image.out image_full.out')

#
# Now plot both
#
im_iso  = readImage(fname='image_iso.out')
im_full = readImage(fname='image_full.out')

plt.figure()
plotImage(im_iso,au=True,log=True,vmax=-10,vmin=-15,bunit='inu',cmap='hot')
plt.text(-90,80,'Isotropic scattering',color='white')
plt.figure()
plotImage(im_full,au=True,log=True,vmax=-10,vmin=-15,bunit='inu',cmap='hot')
plt.text(-90,80,'Full scattering',color='white')

#
# Now compare the cross-sections of the images along the
# major and minor axis
#
plt.figure()
plt.semilogy(im_iso.x/au,im_iso.image[:,100],label='iso major')
plt.semilogy(im_full.x/au,im_full.image[:,100],label='full major')
plt.semilogy(im_iso.x/au,im_iso.image[100,:],label='iso minor')
plt.semilogy(im_full.x/au,im_full.image[100,:],label='full minor')
plt.xlabel('x [au]')
plt.ylabel(r'$I_\nu\;[CGS]$')
plt.ylim((1e-16,1e-10))
plt.title(r'$\lambda\;=\;0.2\;\mu\mathrm{m}$')
plt.legend()

plt.show()
