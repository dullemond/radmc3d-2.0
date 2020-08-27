from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import os

from radmc3dPy.image import *    # Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory
from radmc3dPy.analyze import *  # Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory
from radmc3dPy.natconst import *  # Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory

#
# Call RADMC-3D to make the spectrum
#
os.system('radmc3d spectrum iline 1 widthkms 1')

#
# Read the spectrum
#
s = readSpectrum()

#
# Read the molecule data
#
mol = readMol('co')

#
# Plot line spectrum with doppler shift as x-axis
#
plotSpectrum(s,mol=mol,ilin=1)


plt.show()
