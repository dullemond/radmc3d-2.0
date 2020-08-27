from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import math
import os

# Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory
from radmc3dPy.image import *    
from radmc3dPy.analyze import *  
from radmc3dPy.natconst import * 

#
# Make sure to have run
#
#   radmc3d mctherm
#
os.system('radmc3d mctherm')
#
# to compute the dust temperature before you run this plotting session.
#
# Now plot the temperature profile
#
a    = readData(dtemp=True,binary=False)
r    = a.grid.x[:]
nt   = a.dusttemp.shape[1]
temp = a.dusttemp[:,nt/2,0,0]
plt.figure(1)
plt.plot(r/au,temp)
plt.xlabel('r [au]')
plt.ylabel('T [K]')
plt.show()

#
# Now make sure to have run 
#
#   radmc3d sed incl 60
#
os.system('radmc3d sed incl 60')
#
# to get the spectral energy distribution
#
s    = readSpectrum()
plt.figure(2)
lammic = s[:,0]
flux   = s[:,1]
nu     = 1e4*cc/lammic
nufnu  = nu*flux
nulnu  = nufnu*4*math.pi*pc*pc
plt.plot(lammic,nulnu/ls)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\lambda$ [$\mu$m]')
plt.ylabel(r'$\nu L_{\nu}$ [$L_{\odot}$]')
plt.axis([1e-1,1e4, 1e-6,1e1])
plt.show()

#
# Now make sure to have run 
#
#   radmc3d image lambda 10 incl 60
#
os.system('radmc3d image lambda 10 incl 60')
#
# to get the image
#
im   = readImage()
plt.figure(3)
plotImage(im,log=True,maxlog=6,au=True,bunit='inu')
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
plt.figure(4)
plt.plot(imcirc.rc,imcirc.image[:,0,0,0],label='spherical image')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('r [cm]')
plt.ylabel(r'$I_\nu [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}\,\mathrm{Hz}^{-1}]$')
plt.plot(im.x[50:],im.image[50:,50,0],'o',label='rectangular image')
plt.legend()
plt.show()

#
# Now compare the SED created using the rectangular images
# (standard method) and the SED created using the circular
# images.
#
# Make sure to have run 
#
#   radmc3d sed circ
#   cp spectrum.out spectrum_circ.out
#   radmc3d sed
#   cp spectrum.out spectrum_rect.out
#
os.system('radmc3d sed circ incl 60')
os.system('cp spectrum.out spectrum_circ.out')
os.system('radmc3d sed incl 60')
os.system('cp spectrum.out spectrum_rect.out')
#
# and notice the difference in speed. The results should be 
# mostly identical, but of course because the pixel arrangements
# are a bit different the results will deviate a bit. Comparison
# between these two spectra can also be used to check the accuracy
# of the spectrum.
#
srect  = readSpectrum(fname='spectrum_rect.out')
scirc  = readSpectrum(fname='spectrum_circ.out')
plt.figure(5)
lammic = srect[:,0]
nu     = 1e4*cc/lammic
fluxr  = srect[:,1]
nufnur = nu*fluxr
nulnur = nufnur*4*math.pi*pc*pc
fluxc  = scirc[:,1]
nufnuc = nu*fluxc
nulnuc = nufnuc*4*math.pi*pc*pc
plt.plot(lammic,nulnur/ls,label='SED (from rect images)')
plt.plot(lammic,nulnuc/ls,'o',label='SED (from circ images)')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\lambda$ [$\mu$m]')
plt.ylabel(r'$\nu L_{\nu}$ [$L_{\odot}$]')
plt.axis([1e-1,1e4, 1e-6,1e1])
plt.legend()
plt.show()
