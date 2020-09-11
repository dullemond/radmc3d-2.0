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
# Read the data
#
d     = readData()
rr,tt = np.meshgrid(d.grid.x,d.grid.y,indexing='ij')
zzr   = np.pi/2-tt
rhod_smlgr  = d.rhodust[:,:,0,0]
rhod_biggr  = d.rhodust[:,:,0,1]
if(len(d.dusttemp)>0):
    temp_smlgr  = d.dusttemp[:,:,0,0]
    temp_biggr  = d.dusttemp[:,:,0,1]

#
# Set the radii where to make the following plots
#
rpl = np.array([0.5,3,15,90])*au   # Radii where to make the plots
irr = np.array(np.interp(rpl,d.grid.x,np.arange(len(d.grid.x)))+0.5,dtype=int)  # Nearest radial grid point

#
# Plot the vertical density structure at different radii
#
plt.figure()
icol = 0
for ir in irr:
    r    = d.grid.x[ir]
    rstr = '{0:4.0f}'.format(r/au)
    rstr = 'r = '+rstr.strip()+' au'
    rstr_smlgr = rstr+' (small grains)'
    rstr_biggr = rstr+' (big grains)'
    plt.semilogy(zzr[ir,:],rhod_smlgr[ir,:],label=rstr_smlgr,color='C{}'.format(icol))
    plt.semilogy(zzr[ir,:],rhod_biggr[ir,:],'--',label=rstr_biggr,color='C{}'.format(icol))
    icol += 1
plt.ylim((1e-22,1e-12))
plt.xlabel(r'$\pi/2-\theta\simeq z/r$')
plt.ylabel(r'$\rho_{\mathrm{dust}}\;[\mathrm{g}/\mathrm{cm}^3]$')
plt.legend()

#
# Same as above, but now with grid points, and zoom-in to midplane to
# see the effect of the vertical grid refinement
#
plt.figure()
icol = 0
for ir in irr:
    r    = d.grid.x[ir]
    rstr = '{0:4.0f}'.format(r/au)
    rstr = 'r = '+rstr.strip()+' au'
    rstr_smlgr = rstr+' (small grains)'
    rstr_biggr = rstr+' (big grains)'
    plt.semilogy(zzr[ir,:],rhod_smlgr[ir,:],color='C{}'.format(icol))
    plt.semilogy(zzr[ir,:],rhod_biggr[ir,:],'--',color='C{}'.format(icol))
    plt.semilogy(zzr[ir,:],rhod_smlgr[ir,:],'.',label=rstr_smlgr,color='C{}'.format(icol))
    plt.semilogy(zzr[ir,:],rhod_biggr[ir,:],'o',label=rstr_biggr,color='C{}'.format(icol))
    icol += 1
plt.xlim((0,0.1))
plt.ylim((1e-22,1e-10))
plt.xlabel(r'$\pi/2-\theta\simeq z/r$')
plt.ylabel(r'$\rho_{\mathrm{dust}}\;[\mathrm{g}/\mathrm{cm}^3]$')
plt.legend()

#
# Plot the vertical temperature structure at different radii
#
if(len(d.dusttemp)>0):
    plt.figure()
    icol = 0
    for ir in irr:
        r    = d.grid.x[ir]
        rstr = '{0:4.0f}'.format(r/au)
        rstr = 'r = '+rstr.strip()+' au'
        rstr_smlgr = rstr+' (small grains)'
        rstr_biggr = rstr+' (big grains)'
        plt.semilogy(zzr[ir,:],temp_smlgr[ir,:],label=rstr_smlgr,color='C{}'.format(icol))
        plt.semilogy(zzr[ir,:],temp_biggr[ir,:],'--',label=rstr_biggr,color='C{}'.format(icol))
        icol += 1
    plt.xlabel(r'$\pi/2-\theta\simeq z/r$')
    plt.ylabel(r'$T_{\mathrm{dust}}\;[\mathrm{K}]$')
    plt.legend()

plt.show()
