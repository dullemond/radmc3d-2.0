import problem_setup as p
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt
from matplotlib import cm
from radmc3dPy.image import *
from radmc3dPy.analyze import *

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
xx   = p.xx[:,:,16]
yy   = p.yy[:,:,16]
data = p.rhod[:,:,16]
fig1  = plt.figure()
ax   = fig1.gca(projection='3d')
#ax.plot_wireframe(xx, yy, data, rstride=1, cstride=1)
ax.plot_surface(xx, yy, data, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)

#
# Make and plot an example image
#
makeImage(npix=50,incl=60.,phi=30.,wav=300.,sizeau=45)   # This calls radmc3d 
fig2  = plt.figure()
a=readImage()
plotImage(a,log=True,au=True,maxlog=6,cmap='hot')

#
# Make the SED, by calling RADMC-3D
#
os.system("radmc3d sed incl 60 phi 30")

#
# Plotting it "by hand", the SED as seen at 1 pc distance
#
fig3  = plt.figure()
s     = readSpectrum()
lam   = s[:,0]
nu    = 1e4*cc/lam
fnu   = s[:,1]
nufnu = nu*fnu
plt.plot(lam,nufnu)
plt.xscale('log')
plt.yscale('log')
plt.axis([1e-1, 1e4, 1e-10, 1e-4])
plt.xlabel('$\lambda\; [\mu \mathrm{m}$]')
plt.ylabel('$\\nu F_\\nu \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')

#
# Use the radmc3dPy.analyze tool set for plotting the SED, 
# this time let's plot nuLnu in units of Lsun
#
fig4  = plt.figure()
plotSpectrum(s,nulnu=True,lsun=True,xlg=True,ylg=False,micron=True)
plt.axis([1e-1,1e4,1e-8,1])

#
# Use the radmc3dPy.analyze tool set to read in the dust temperatures
# calculated by the command-line command "radmc3d mctherm"
#
q     = readData(dtemp=True,binary=False)
fig5  = plt.figure()
ay    = fig5.gca(projection='3d')
qq    = np.meshgrid(q.grid.x,q.grid.y,q.grid.z,indexing='ij')
xx    = qq[0][:,:,16]
yy    = qq[1][:,:,16]
dd    = q.dusttemp[:,:,16,0]  # The extra "0" is because of possible multiple dust species; here only 1 dust species
ay.plot_surface(xx, yy, dd, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=1)


plt.show()
