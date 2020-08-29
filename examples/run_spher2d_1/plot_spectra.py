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
# Now make sure to have run 
#
#   radmc3d sed incl 60
#
os.system('radmc3d sed incl 60')
#
# to get the spectral energy distribution
#
s    = readSpectrum()
plt.figure()
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
plt.figure()
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
