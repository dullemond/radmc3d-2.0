#
# For this to work you have to have the python script "interactive_plot.py" somehow reachable
# by Python (e.g. using PYTHON_PATH or copying it here). You can obtain it from:
#
#   https://github.com/dullemond/interactive_plot
#

import numpy as np
from interactive_plot import *
from matplotlib import cm
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage

from radmc3dPy.analyze import *  
from radmc3dPy.natconst import * 

a        = readData()

iplot    = 'rhodust'
#iplot    = 'dusttemp'

if iplot=='rhodust':
    data     = a.rhodust.copy()
    dmax     = data.max()
    dmin     = dmax*1e-10
    data[data<dmin]=np.nan
    data     = np.log10(data)
    dmin     = np.log10(dmin)
    dmax     = np.log10(dmax)
    label    = r'$\rho_d\;[\mathrm{g/cm}^3]$'
    cmap     = cm.Blues
if iplot=='dusttemp':
    data     = a.dusttemp
    dmin     = data.min()
    dmax     = data.max()
    label    = r'$T\;[\mathrm{K}]$'
    cmap     = cm.hot

x        = a.grid.x
y        = (np.pi/2-a.grid.y)[::-1]
z        = data[:,::-1,:,:]
norm     = colors.Normalize(vmin=dmin,vmax=dmax)
fig,ax   = plt.subplots()
im       = NonUniformImage(ax,interpolation='nearest',cmap=cmap,norm=norm)
im.set_data(x,y,z[:,:,0,0].T)
ax.images.append(im)
ax.set_xlim((x[0]-0.5*(x[1]-x[0]),x[-1]+0.5*(x[-1]-x[-2])))
ax.set_ylim((y[0]-0.5*(y[1]-y[0]),y[-1]+0.5*(y[-1]-y[-2])))
cbar=fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap), ax=ax)
cbar.set_label(label)
def img_func(param,fixedpar={}): return fixedpar['z'][:,:,param[0],0]
params = [np.arange(z.shape[-2])] # Choices of parameter values
fixedpar = {}
fixedpar["z"]=z
paramsalt = [a.grid.z*180/np.pi,np.arange(data.shape[-1])]
interactive_plot(None, None, params, fixedpar=fixedpar,       \
                 img_x=x,img_y=y,img_func=img_func,img_im=im, \
                 fig=fig,ax=ax,paramsalt=paramsalt,           \
                 parnames=['phi [deg] = '],          \
                 altformat='3.0f')
