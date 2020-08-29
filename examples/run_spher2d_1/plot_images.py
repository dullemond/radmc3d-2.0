from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import math
import os

from radmc3dPy.image import *    
from radmc3dPy.analyze import *  
from radmc3dPy.natconst import * 

#
# First set up the model with
#
#   python problem_setup.py
#
# Then make sure to have run
#
#   radmc3d mctherm
#
#os.system('radmc3d mctherm')

#
# Now make the image
#
#   radmc3d image lambda 10 incl 60
#
os.system('radmc3d image lambda 10 incl 60')
#
# to get the image
#
im   = readImage()
dum  = im.image.shape
nx   = dum[0]
nx2  = nx//2
plt.figure()
plotImage(im,log=True,maxlog=6,au=True,bunit='inu',cmap=cm.hot)
plt.show()

#
# Now make a strongly zoomed-in image
# This time we plot linear intensity, not logarithmic,
# but we saturate the starlight
#
#   radmc3d image lambda 10 incl 60 sizeau 10
#
os.system('radmc3d image lambda 10 incl 60 sizeau 10')
#
# to get the image
#
imz  = readImage()
dum  = imz.image.shape
nx   = dum[0]
nx2  = nx//2
vmax = 5e-10

plt.figure()
plotImage(imz,au=True,bunit='inu',vmax=vmax,cmap=cm.hot)
plt.show()

#
# Now make sure to have run 
#
#   radmc3d image circ lambda 10 incl 60
#
os.system('radmc3d image circ lambda 10 incl 60')
#
# to get the "circular image". Circular images (available only for
# models in spherical coordinates) are pixel arrangements that are
# not as rows and columns, but instead as concentric circles. This 
# can be useful for models with radial coordinates spanning a huge
# range (from very small r to very large r, several orders of 
# magnitude). For circular images we do not need "subpixeling"
# or "zoom ins", because the pixels are arranged in circles, each
# corresponding to a radial coordinate of the model. So this "zoom
# in" is automatic: all scales are represented in one "circular
# image". They are also useful for 1-D models, since 2-D images of
# 1-D spherically symmetric models are overkill: we only need to 
# compute the intensity as a function of impact parameter. So for
# 1-D models the circular images have no subdivision in phi, but
# only intensity as a function of r. That makes it, as a bonus,
# also much faster.
#
# Here we overplot the result of the spherical image with that
# of the (normal) rectangular image.
#
imcirc = readcircimage()
plt.figure()
plt.plot(imcirc.rc,imcirc.image[:,0,0,0],label='spherical image')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('r [cm]')
plt.ylabel(r'$I_\nu [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}\,\mathrm{Hz}^{-1}]$')
plt.plot(im.x[nx2:],im.image[nx2:,nx2,0],'o',label='rectangular image')
plt.legend()
plt.show()

#
# You can also show a "heatmap" of the circular image
#
#   The radial span
#   (but note that the r-grid is not fully log, so this is not 100% correct)
#
lgrin  = np.log10(imcirc.rc[1]/au)
lgrout = np.log10(imcirc.rc[-1]/au)
#
#   Make the "heatmap" figure of the 10log of the circular image
#
plt.figure()
plt.imshow(np.log10(imcirc.image[1:,:,0,0]+1e-23),vmin=-16,aspect='auto',cmap=cm.hot,origin='lower',extent=[0,360,lgrin,lgrout])
plt.title(r'$\lambda = 10\,\mu\mathrm{m}$')
plt.xlabel(r'$\phi\; [deg]$')
plt.ylabel(r'$^{10}\log(r)\; [\mathrm{AU}]$')
cbar=plt.colorbar()
cbar.set_label(r'$^{10}\log(I_\nu)\;[\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}\,\mathrm{Hz}^{-1}\,\mathrm{ster}^{-1}\,]$')
plt.show()
