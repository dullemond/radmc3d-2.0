from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import math
import os

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
# to compute the dust temperature before you run this plotting session.
#
# Now plot the temperature profile
#
a    = readData()
r    = a.grid.x[:]
temp = a.dusttemp[:,-1,0,0]
plt.figure()
plt.plot(r/au,temp)
plt.xlabel('r [au]')
plt.ylabel('T [K]')
plt.show()

#
# The "heatmap" of the 2-D temperature profile
#
lgrin  = np.log10(a.grid.x[0]/au)
lgrout = np.log10(a.grid.x[-1]/au)
plt.figure()
plt.imshow(a.dusttemp[:,:,0,0].T,extent=[lgrin,lgrout,0,np.pi/2-a.grid.y[0]],aspect='auto',cmap=cm.hot)
plt.xlabel(r'$^{10}\log(r)\; [\mathrm{AU}]$')
plt.ylabel(r'$\pi/2-\theta$')
cbar=plt.colorbar()
cbar.set_label(r'$T\;[\mathrm{K}]$')
plt.show()

#
# The "heatmap" of the 2-D density profile
#
lgrin  = np.log10(a.grid.x[0]/au)
lgrout = np.log10(a.grid.x[-1]/au)
plt.figure()
plt.imshow(np.log10(a.rhodust[:,:,0,0].T),extent=[lgrin,lgrout,0,np.pi/2-a.grid.y[0]],aspect='auto',cmap=cm.Blues)
plt.xlabel(r'$^{10}\log(r)\; [\mathrm{AU}]$')
plt.ylabel(r'$\pi/2-\theta$')
cbar=plt.colorbar()
cbar.set_label(r'$^{10}\log(\rho)\;[\mathrm{g}/\mathrm{cm}^3]$')
plt.show()
