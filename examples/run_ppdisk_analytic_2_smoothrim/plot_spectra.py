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
# Make and plot the SED as seen at 1 pc distance
#
os.system("radmc3d sed incl 60 phi 30")
fig3  = plt.figure()
s     = readSpectrum()
star  = readStars()
lam   = s[:,0]
nu    = 1e4*cc/lam
fnu   = s[:,1]
nufnu = nu*fnu
plt.loglog(lam,nufnu,label='Total SED')
plt.loglog(lam,nu*star.fnustar[:,0],label='Star')
plt.axis([1e-1, 1e4, 1e-8, 1e-2])
plt.xlabel('$\lambda\; [\mu \mathrm{m}$]')
plt.ylabel('$\\nu F_\\nu \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
plt.legend()
plt.show()
