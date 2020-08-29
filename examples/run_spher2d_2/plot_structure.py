from matplotlib import cm
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
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
plt.plot(r/au,temp,'.')
plt.xlabel('r [au]')
plt.ylabel('T [K]')
plt.show()

#
# The "heatmap" of the 2-D temperature profile
#
# Due to the non-uniform grid in r, we cannot simply use
# plt.imshow(), but must use NonUniformImage(). 
#
# We zoom in to the inner region: from 5 to 9 au, to show the
# grid refinement near the inner edge. But you can set the
# zoomr to 1. to see the entire grid from 5 to 100 au.
#
cmap   = cm.hot
zoomr  = 5.    # Zoom factor, to see better the inner rim
fig,ax = plt.subplots()
x      = np.log10(a.grid.x/au)
y      = (np.pi/2-a.grid.y)[::-1]
xi     = np.log10(a.grid.xi/au)
yi     = (np.pi/2-a.grid.yi)[::-1]
z      = a.dusttemp[:,::-1,0,0].T
im     = NonUniformImage(ax,interpolation='nearest',cmap=cmap)
im.set_data(x,y,z)
ax.images.append(im)
ax.set_xlim((xi[0],(xi[-1]-xi[0])/zoomr+xi[0]))
ax.set_ylim((yi[0],yi[-1]))
plt.xlabel(r'$^{10}\log(r)\; [\mathrm{AU}]$')
plt.ylabel(r'$\pi/2-\theta$')
norm=colors.Normalize(vmin=z.min(),vmax=z.max())
cbar=fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap), ax=ax)
cbar.set_label(r'$T\;[\mathrm{K}]$')
plt.show()
