# This module is meant to allow easier quick-and-dirty analysis of
# data without having to cut-and-paste stuff from templates. These
# routines are presumably not really useful for final paper-grade
# plotting. 
#
# Usage example:
#
#   from simpleplot import *
#   f=np.load('mydata.npy')
#   surface(f)
#   image(f)
#   plot(f[:,30])
#
# or more sophisticated:
#
#   from plot import *
#   f=np.load('mydata.npy')
#   x=np.load('x.npy')
#   y=np.load('y.npy')
#   fig,ax = surface(f,x=x,y=y,cstride=4,rstride=4)
#   ax.set_xlabel('X Label')
#   ax.set_ylabel('Y Label')
#   fig,ax = image(f,x=x,y=y,aspect=1,cmap=cm.hot,range=[-4.5,-3.])
#   fig,ax = plot(f[:,30],x=x)
#   oplot(f[:,100],x=x)
#
# C.P. Dullemond (2020)

from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

def plot(x,f=None,oplot=None,xlog=False,ylog=False,xrange=None,yrange=None):
    if oplot is None:
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        if xlog:
            plt.xscale('log')
        if ylog:
            plt.yscale('log')
        if xrange is not None:
            plt.xlim(xrange)
        if yrange is not None:
            plt.ylim(yrange)
    else:
        fig = None
        ax  = None
    if f is None:
        func = x
        nx   = len(x)
        xx   = np.linspace(0,nx-1,nx)
    else:
        func = f
        xx   = x
    plt.plot(xx,func)
    return fig,ax
    

def surface(ff,x=None,y=None,rstride=None,cstride=None,stride=None,
            xlabel=None,ylabel=None,zlabel=None):
    f = np.squeeze(ff)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    nx = f.shape[0]
    ny = f.shape[1]
    if x is None:
        x = np.linspace(0,nx-1,nx)
    if y is None:
        y = np.linspace(0,ny-1,ny)
    if stride is not None:
        rstride = stride
        cstride = stride
    if rstride is None:
        rstride = 1
    if cstride is None:
        cstride = 1
    xx, yy = np.meshgrid(x, y,indexing='ij')
    ax.plot_wireframe(xx, yy, f, rstride=rstride, cstride=cstride)
    if xlabel is not None: ax.set_xlabel(xlabel)
    if ylabel is not None: ax.set_ylabel(ylabel)
    if zlabel is not None: ax.set_zlabel(zlabel)
    return fig,ax

def image(ff,x=None,y=None,aspect='equal',cmap=None,range=None):
    f = np.squeeze(ff)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if x is None:
        nx = f.shape[1]
        x  = np.array((0,nx-1))
    if y is None:
        ny = f.shape[0]
        y  = np.array((0,ny-1))
    extent = (x[0],x[len(x)-1],y[0],y[len(y)-1])
    if range is None:
        ff = f.copy()
    else:
        ff = f.copy()
        mi = range[0]
        ma = range[1]
        if mi > ma:
            ma = range[0]
            mi = range[1]        
        ff[ff>ma] = ma
        ff[ff<mi] = mi
    plt.imshow(ff,origin='lower',aspect=aspect,extent=extent,cmap=cmap)
    return fig,ax
    
