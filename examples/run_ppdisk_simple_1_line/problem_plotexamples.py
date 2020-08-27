import problem_setup as p
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt
from matplotlib import cm
from radmc3dPy.image import *    # Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory
from radmc3dPy.analyze import *  # Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory
from radmc3dPy.natconst import * # Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory
import os

#
# Make sure to have done the following beforhand:
#
#  First compile RADMC-3D
#  Then run:
#   python problem_setup.py
#   radmc3d mctherm
#

#
# Now let's make a set of channel maps
#
vkms = 1.5
#command = 'radmc3d image imolspec 1 iline 2 linenlam 20 widthkms 10. theta 45 npix 200'
command  = 'radmc3d image imolspec 1 iline 2 vkms {} theta 45 npix 200'.format(vkms) 
os.system(command)

#
# Read the image
#
im = readImage()

#
# Plot
#
#plotImage(im,ifreq=18,au=True,cmap=cm.hot)
plotImage(im,au=True,cmap=cm.hot)

plt.show()
