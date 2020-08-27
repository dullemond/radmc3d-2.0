import problem_setup as p
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt
from matplotlib import cm
from radmc3dPy.image import *    # Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory
from radmc3dPy.analyze import *  # Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory
from radmc3dPy.natconst import * # Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory

#
# Make sure to have done the following beforhand:
#
#  First compile RADMC-3D
#  Then run:
#   python problem_setup.py
#   radmc3d mctherm
#

#
# Make and plot an example image at edge-on view (to better see the effect
# of the two dust species)
#
fig2  = plt.figure()
makeImage(npix=100,incl=60.,wav=30,sizeau=300)   # This calls radmc3d 
a=readImage()
plotImage(a,log=True,au=True,maxlog=4,cmap='hot')

#
# Make and plot the SED as seen at 1 pc distance
#
#os.system("radmc3d sed incl 60 phi 30")
#fig3  = plt.figure()
#s     = readSpectrum()
#lam   = s[:,0]
#nu    = 1e4*cc/lam
#fnu   = s[:,1]
#nufnu = nu*fnu
#plt.plot(lam,nufnu)
#plt.xscale('log')
#plt.yscale('log')
#plt.axis([1e-1, 1e4, 1e-10, 1e-4])
#plt.xlabel('$\lambda\; [\mu \mathrm{m}$]')
#plt.ylabel('$\\nu F_\\nu \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
