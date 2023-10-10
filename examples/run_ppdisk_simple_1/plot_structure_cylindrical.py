import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt
from matplotlib import cm
from radmc3dPy.image import *
from radmc3dPy.analyze import *
from radmc3dPy.natconst import *
from scipy.interpolate import RegularGridInterpolator

#
# Make sure to have done the following beforhand:
#
#  First compile RADMC-3D
#  Then run:
#   python problem_setup.py
#   radmc3d mctherm
#

#
# Read the data
#
d     = readData()
rhod  = d.rhodust[:,:,0,0]
temp  = d.dusttemp[:,:,0,0]

#
# Prepare mapping to cylindrical coordinates
#
rc    = d.grid.x                               # Choose 1D cylindrical radius grid equal to 1D spherical radius grid
theta = d.grid.y
zrc   = 1/np.tan(theta)[::-1]                  # zrc=z/rc increasing from midplane upward
rrc,zzrc = np.meshgrid(rc,zrc,indexing='ij')   # Make a regular 2D cylindrical grid from these
rrs   = rrc*np.sqrt(1+zzrc**2)                 # For each of these 2D points, compute spherical radius rrs
tts   = np.arctan(1/zzrc)                      # For each of these 2D points, compute spherical thetas
pts   = np.array((rrs,tts)).transpose([1,2,0]) # Combine these as a 2D array of 2D points

#
# Now do the mapping
#
fill_value = np.nan   # For plotting this is better (fill out-of-bounds cells with NaN)
#fill_value = 0.       # For further data manipulation this is better (fill out-of-bounds cells with 0)
method     = 'linear' # Simplest interpolation
#method     = 'cubic'  # Better interpolation
intrp = RegularGridInterpolator((d.grid.x,d.grid.y),rhod,bounds_error=False,fill_value=fill_value,method=method)
rhodc = intrp(pts)
intrp = RegularGridInterpolator((d.grid.x,d.grid.y),temp,bounds_error=False,fill_value=fill_value,method=method)
tempc = intrp(pts)

#
# View a surface plot of the density structure
#
fig1 = plt.figure()
ax   = fig1.gca(projection='3d')
ax.plot_surface(np.log10(rrc)/au, zzrc, np.log10(rhodc), rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=1, antialiased=False)

#
# Plot the vertical density structure at different radii
#
irr = [0,10,20,30]
plt.figure()
for ir in irr:
    r    = d.grid.x[ir]
    rstr = '{0:4.0f}'.format(r/au)
    rstr = 'r = '+rstr.strip()+' au'
    plt.semilogy(zzrc[ir,:],rhodc[ir,:],label=rstr)
plt.ylim((1e-25,1e-15))
plt.xlabel(r'$z/r_{\mathrm{cyl}}$')
plt.ylabel(r'$\rho_{\mathrm{dust}}\;[\mathrm{g}/\mathrm{cm}^3]$')
plt.legend()

#
# Plot the vertical temperature structure at different radii
#
irr = [0,10,20,30]
plt.figure()
for ir in irr:
    r    = d.grid.x[ir]
    rstr = '{0:4.0f}'.format(r/au)
    rstr = 'r = '+rstr.strip()+' au'
    plt.plot(zzrc[ir,:],tempc[ir,:],label=rstr)
plt.xlabel(r'$z/r_{\mathrm{cyl}}$')
plt.ylabel(r'$T_{\mathrm{dust}}\;[\mathrm{K}]$')
plt.legend()

plt.show()
