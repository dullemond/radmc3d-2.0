from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import math

def rebin( a, newshape ):
    '''Rebin an array to a new shape.
    '''
    assert len(a.shape) == len(newshape)
    slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
    coordinates = np.mgrid[slices]
    indices = coordinates.astype('i')   #choose the biggest smaller integer index
    return a[tuple(indices)]

def plotpoldir(x,y,ii,qq,uu,nx=20,ny=20):
    xr   = rebin(x,[nx])
    yr   = rebin(y,[ny])
    mesh = np.meshgrid(xr,yr,indexing='ij')
    xxr  = mesh[0]
    yyr  = mesh[1]
    qqr  = rebin(np.squeeze(qq)/(np.squeeze(ii)+1e-60),[ny,nx])
    uur  = rebin(np.squeeze(uu)/(np.squeeze(ii)+1e-60),[ny,nx])
    lpol = np.sqrt(qqr*qqr+uur*uur)
    qqr  = qqr / (lpol+1e-60)
    uur  = uur / (lpol+1e-60)
    angg = np.arccos(qqr)/2.e0
    angg[uur<0] = math.pi - angg[uur<0]
    vx   = np.cos(angg)
    vy   = np.sin(angg)
    vx[lpol<1e-6]=0.001
    vy[lpol<1e-6]=0.001
    q = plt.quiver(xxr,yyr,vx,vy,
      color='white',pivot='mid',scale=2.*np.max([nx,ny]),
      headwidth=1e-10,headlength=1e-10,headaxislength=1e-10)

