import numpy as np
from natconst import *
from numpy.random import default_rng
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve

def kernel1d(r,h):
    # Spline kernel of Gadget-2 (Springel MNRAS 364, 1105, 2005, Eq. 4)
    x  = r/h
    w0 = 8/(np.pi*h**3)
    w  = np.zeros_like(r)
    i0 = x>=1
    i1 = np.logical_and(x>=0,x<0.5)
    i2 = np.logical_and(x>=0.5,x<1)
    w[i1] = 1-6*x[i1]**2+6*x[i1]**3
    w[i2] = 2*(1-x[i2])**3
    return w

def kernel3d(dlgr,dth,dphi,hr):
    n_r   = (int(2*hr/dlgr)//2)*2+1
    n_th  = (int(2*hr/dth)//2)*2+1
    n_phi = (int(2*hr/dphi)//2)*2+1
    rc    = np.arange(-n_r//2+1,n_r//2+1)*dlgr
    thc   = np.arange(-n_th//2+1,n_th//2+1)*dth
    phic  = np.arange(-n_phi//2+1,n_phi//2+1)*dphi
    r     = np.sqrt(rc[:,None,None]**2+thc[None,:,None]**2+phic[None,None,:]**2)
    w     = kernel1d(r,hr)
    w     = w / w.sum()
    return w

def uniform_spher_loggrid(rin,rout,hrup,nr,nth,nph):
    ri    = rin * (rout/rin)**np.linspace(0,1,nr+1)
    thup  = np.pi/2-hrup
    thlo  = np.pi/2+hrup
    thi   = thup + (thlo-thup)*np.linspace(0,1,nth+1)
    phi   = 2*np.pi*np.linspace(0,1,nph+1)
    return ri,thi,phi

rng     = default_rng()

#
# Parameters of the dummy model
#
mdisk      = 1e-5*MS
hpr        = 0.05
hrkernel   = 0.03

#
# Dummy model for the SPH particles
#
nsph       = 1000000
sph_r      = au*10**(2*rng.uniform(size=nsph))
sph_th     = np.pi/2 - hpr*rng.standard_normal(size=nsph)
sph_phi    = 2*np.pi*rng.uniform(size=nsph)
sph_lr     = np.log(sph_r)

#
# Set up the grid
#
grid_nr    = 100
grid_nth   = 100
grid_nph   = 100
grid_rin   = 1*au
grid_rout  = 100*au
grid_hrup  = 0.3
grid_ri, grid_thi, grid_phi = uniform_spher_loggrid(grid_rin,grid_rout,grid_hrup,grid_nr,grid_nth,grid_nph)
grid_lri   = np.log(grid_ri)

#
# Find particles in the grid
#
sph_ir     = np.array(np.interp(sph_lr,grid_lri,np.arange(len(grid_lri))),dtype=int)
sph_ith    = np.array(np.interp(sph_th,grid_thi,np.arange(len(grid_thi))),dtype=int)
sph_iphi   = np.array(np.interp(sph_phi,grid_phi,np.arange(len(grid_phi))),dtype=int)

#
# Check
#
sph_r_eps  = (sph_r-grid_ri[sph_ir])/(grid_ri[sph_ir+1]-grid_ri[sph_ir])
sph_th_eps = (sph_th-grid_thi[sph_ith])/(grid_thi[sph_ith+1]-grid_thi[sph_ith])
sph_phi_eps= (sph_phi-grid_phi[sph_iphi])/(grid_phi[sph_iphi+1]-grid_phi[sph_iphi])
assert len(np.where(sph_r_eps<0)[0])==0, 'Error in r-grid-finding'
assert len(np.where(sph_r_eps>1)[0])==0, 'Error in r-grid-finding'
assert len(np.where(sph_th_eps<0)[0])==0, 'Error in th-grid-finding'
assert len(np.where(sph_th_eps>1)[0])==0, 'Error in th-grid-finding'
assert len(np.where(sph_phi_eps<0)[0])==0, 'Error in phi-grid-finding'
assert len(np.where(sph_phi_eps>1)[0])==0, 'Error in phi-grid-finding'

#
# Add particles to the grid
#
print('Adding particles to the grid')
grid_count = np.zeros((grid_nr,grid_nth,grid_nph))
for isph in range(nsph):
    ir   = sph_ir[isph]
    ith  = sph_ith[isph]
    iphi = sph_iphi[isph]
    grid_count[ir,ith,iphi] += 1

#
# Create a kernel
#
dlgr = grid_lri[1]-grid_lri[0]
dth  = grid_thi[1]-grid_thi[0]
dphi = grid_phi[1]-grid_phi[0]
kern = kernel3d(dlgr,dth,dphi,hrkernel)

#
# Convolve
#
print('Convolving with SPH kernel using scipy.signal.fftconvolve')
grid_convolve = fftconvolve(grid_count,kern)

#
# Show
#
print('Plotting the two grids')
plt.figure()
plt.imshow(grid_count.sum(axis=2).T,origin='lower')
plt.figure()
plt.imshow(grid_convolve.sum(axis=2).T,origin='lower')
plt.show()
