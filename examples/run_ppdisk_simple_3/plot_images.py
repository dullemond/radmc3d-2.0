import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from radmc3dPy.image import *
from radmc3dPy.analyze import *
from radmc3dPy.natconst import *

#
# Make sure to have done the following beforhand:
#
#  First compile RADMC-3D
#  Then run:
#   python problem_setup.py
#   radmc3d mctherm
#

dpc = 140.     # Distance in parsec (for conversion to Jy/pixel in 1.3 mm map)

#
# Make and plot image of full disk at 1.3 mm: thermal dust emission
#
makeImage(npix=200,incl=60.,phi=30.,wav=1.3e3,sizeau=200)   # This calls radmc3d 
im_mm = readImage()
plt.figure()
plotImage(im_mm,au=True,log=True,maxlog=3,bunit='jy/pixel',dpc=dpc,cmap='magma')

#
# Make and plot image of full disk at 1 microns: scattered light
# and show the difference with and without scattering
#
makeImage(npix=200,incl=60.,phi=30.,wav=1.,sizeau=200)   # This calls radmc3d 
im_1 = readImage()
makeImage(npix=200,incl=60.,phi=30.,wav=1.,sizeau=200,noscat=True)   # This calls radmc3d 
im_1_noscat = readImage()
plt.figure()
plotImage(im_1_noscat,au=True,log=True,vmax=-10,vmin=-15,bunit='inu',cmap='hot')
plt.text(-90,80,'No scattering',color='white')
plt.figure()
plotImage(im_1,au=True,log=True,vmax=-10,vmin=-15,bunit='inu',cmap='hot')
plt.text(-90,80,'With isotropic scattering',color='white')

#
# Make and plot image of inner rim at 3 microns, this time in linear color scale
#
makeImage(npix=200,incl=60.,phi=30.,wav=3.,sizeau=2)   # This calls radmc3d 
im_rim = readImage()
plt.figure()
plotImage(im_rim,au=True,bunit='inu',vmax=3e-7,cmap='hot')
