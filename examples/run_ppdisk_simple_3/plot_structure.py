import problem_setup as p
import numpy as np
from mpl_toolkits.mplot3d import axes3d
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

#
# View a 2-D slice of the 3-D array of the setup
#
rr   = p.rr[:,:,0]
tt   = p.tt[:,:,0]
zzr  = np.pi/2-tt
rhod = p.rhod[:,:,0]
fig1 = plt.figure()
ax   = fig1.gca(projection='3d')
ax.plot_surface(np.log10(rr)/au, zzr, np.log10(rhod), rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=1, antialiased=False)
#ax.set_zlim3d(0, 10)

