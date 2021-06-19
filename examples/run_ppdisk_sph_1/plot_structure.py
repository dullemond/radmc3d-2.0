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
plt.title('Midplane dust temperature')
plt.show()

#
# The "heatmap" of the 2-D temperature profile
#
iphi   = 10    # Azimuthal index at which the temperature structure is shown
cmap   = cm.hot
zoomr  = 1.    # Zoom factor
fig,ax = plt.subplots()
x      = np.log10(a.grid.x/au)
y      = (np.pi/2-a.grid.y)[::-1]
xi     = np.log10(a.grid.xi/au)
yi     = (np.pi/2-a.grid.yi)[::-1]
z      = a.dusttemp[:,::-1,iphi,0].T
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
plt.title(r'Dust temperature structure at $\phi=$'+'{0:3d}'.format(int(a.grid.z[iphi]*180/np.pi))+r'$^{\circ}$')

#
# Compute the RADIAL tau=1 surface for stellar radiation
# We focus only on the small grain component here
#
osmall  = readOpac(ext='silicate',scatmat=False)
lamstar = 0.45
kappa   = np.interp(lamstar,osmall.wav[0],osmall.kabs[0]+osmall.ksca[0])
alpha   = kappa*a.rhodust[:,:,:,0]
rri,tt,pp = np.meshgrid(a.grid.xi,a.grid.y,a.grid.z,indexing='ij')
drr     = rri[1:,:,:]-rri[:-1,:,:]
taurad  = (alpha*drr).cumsum(axis=0)
tausurf = np.zeros((a.grid.ny,a.grid.nz))
for iz in range(a.grid.nz):
    for iy in range(a.grid.ny):
        tausurf[iy,iz] = np.interp(1.0,taurad[:,iy,iz],rri[:-1,iy,iz])
mask = tausurf>a.grid.xi[-3]
tausurf[mask] = 1.1*a.grid.xi[-1]
assert len(mask[0])>0, "ERROR in setup: Vertical extent of the grid too small (upper layers of disk are cut off)."

plt.plot(np.log10(tausurf[:,0]/au),np.pi/2-a.grid.y,color='blue')

#
# Plot the optical depth
#
sig_d   = []
opacs   = []
with open('dustopac.inp','r') as f:
    str=f.readline()
    str=f.readline()
    str=f.readline()
    for i in range(a.rhodust.shape[-1]):
        a.getSigmaDust(i)
        sig_d.append(a.sigmadust)
        str=f.readline()
        str=f.readline()
        str=f.readline()
        ext=str.split()[0]
        o=readOpac(ext=ext,scatmat=False)
        opacs.append(o)
        str=f.readline()
        
lams    = [0.45,15,1300]
taus    = []
for lam in lams:
    kappas = []
    for o in opacs:
        kap = np.interp(lam,o.wav[0],o.kabs[0]+o.ksca[0])
        kappas.append(kap)
    tau = np.zeros((len(opacs),a.grid.nx,a.grid.nz))
    for i,s in zip(np.arange(len(sig_d)),sig_d):
        tau[i,:,:] = s*kappas[i]
    taus.append(tau)

plt.figure()
for lam,tau in zip(lams,taus):
    plt.loglog(a.grid.x/au,tau.sum(0)[:,0],label=r'$\lambda=$'+'{}'.format(lam)+r' $\mu\mathrm{m}$')
plt.xlabel('r [au]')
plt.ylabel(r'$\tau$ (at $\phi=0$)')
plt.legend()

#
# The "heatmap" of the 2-D density profile
#
cmap   = cm.Blues
zoomr  = 1.    # Zoom factor
fig,ax = plt.subplots()
x      = np.log10(a.grid.x/au)
y      = (np.pi/2-a.grid.y)[::-1]
xi     = np.log10(a.grid.xi/au)
yi     = (np.pi/2-a.grid.yi)[::-1]
z      = a.rhodust[:,::-1,0,:].sum(axis=-1).T.copy()
zmax   = z.max()
zmin   = zmax*1e-10
z[z<zmin]=np.nan
z      = np.log10(z)
im     = NonUniformImage(ax,interpolation='nearest',cmap=cmap)
im.set_data(x,y,z)
ax.images.append(im)
ax.set_xlim((xi[0],(xi[-1]-xi[0])/zoomr+xi[0]))
ax.set_ylim((yi[0],yi[-1]))
plt.xlabel(r'$^{10}\log(r)\; [\mathrm{AU}]$')
plt.ylabel(r'$\pi/2-\theta$')
norm=colors.Normalize(vmin=np.log10(zmin),vmax=np.log10(zmax))
cbar=fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap), ax=ax)
cbar.set_label(r'$^{10}\log(\rho)\;[\mathrm{g/cm}^3]$')
plt.title('Dust density structure')

#
# Compute the VERTICAL tau=1 surface for various wavelengths
#
tauverts  = []
zrsurf    = []
rr,tti,pp = np.meshgrid(a.grid.x,a.grid.yi,a.grid.z,indexing='ij')
dzz       = rr[:,:-1,:]*(tti[:,1:,:]-tti[:,:-1,:])
for lam in lams:
    kappas = []
    for o in opacs:
        kap = np.interp(lam,o.wav[0],o.kabs[0]+o.ksca[0])
        kappas.append(kap)
    tau = np.zeros((len(opacs),a.grid.nx,a.grid.ny,a.grid.nz))
    for i in range(len(opacs)):
        tau[i,:,:,:] = (kappas[i]*a.rhodust[:,:,:,i]*dzz).cumsum(1)
    tauverts.append(tau)
    zr = np.zeros((a.grid.nx,a.grid.nz))
    for ir in range(a.grid.nx):
        for ip in range(a.grid.nz):
            zr[ir,ip] = np.pi/2-np.interp(1.0,np.squeeze(tau[:,ir,:,ip].sum(0)),np.squeeze(tti[ir,1:,ip]))
    zrsurf.append(zr)
    plt.plot(np.log10(a.grid.x/au),zr[:,0],label=r'$\tau_{\mathrm{vert}}=1\;\mathrm{at}\;\lambda=$'+'{}'.format(lam)+r' $\mu\mathrm{m}$')
plt.legend()

plt.show()
