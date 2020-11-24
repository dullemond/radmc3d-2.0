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

#setthreads = None
setthreads = 4

dpc = 140.     # Distance in parsec (for conversion to Jy/pixel in 1.3 mm map)

#
# Make and plot image of full disk at 1.3 mm: thermal dust emission
#
makeImage(npix=200,incl=60.,phi=0.,wav=1.3e3,setthreads=setthreads)   # This calls radmc3d 
im_mm = readImage()
plt.figure()
plotImage(im_mm,au=True,log=True,maxlog=3,bunit='jy/pixel',dpc=dpc,cmap='magma')

#
# Make and plot image of full disk at 1 microns: scattered light
#
makeImage(npix=200,incl=60.,phi=0.,wav=1.,setthreads=setthreads)   # This calls radmc3d 
im_1 = readImage()
plt.figure()
plotImage(im_1,au=True,log=True,vmax=-10,vmin=-15,bunit='inu',cmap='hot')
plt.text(-150,140,'With isotropic scattering',color='white')
