from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import math
import os

from radmc3dPy.analyze import *  
from radmc3dPy.natconst import * 

def plot_spher_grid_rtheta(grid,rmax,unit=1.,linewidth=1):
    if grid.y.max()>np.pi/2:
        extent = [0,rmax/unit,-rmax/unit,rmax/unit]
        nth    = 181
    else:
        extent = [0,rmax/unit,0,rmax/unit]
        nth    = 91
    th_min     = grid.yi[0]
    th_max     = grid.yi[-1]
    th         = np.linspace(th_min,th_max,nth)
    halfcirc_x = np.sin(th)
    halfcirc_y = np.cos(th)
    rin        = grid.xi[0]
    rout       = grid.xi[-1]
    for r in grid.xi:
        plt.plot((r/unit)*halfcirc_x,(r/unit)*halfcirc_y,color='black',linewidth=linewidth)
    for t in grid.yi:
        plt.plot([(rin/unit)*np.sin(t),(rout/unit)*np.sin(t)],
                 [(rin/unit)*np.cos(t),(rout/unit)*np.cos(t)],color='black',linewidth=linewidth)
    plt.xlim(extent[0:2])
    plt.ylim(extent[2:4])
    plt.show()

#
# Get the grid from the model
#
grid   = readGrid()

#
# Plot the location of the radial grid walls
#
plt.figure()
plt.semilogy(grid.xi/au,label='r')
plt.semilogy((grid.xi[1:]-grid.xi[:-1])/au,label=r'$\Delta r$')
plt.xlabel('radial grid index')
plt.ylabel('r [au]')
plt.legend()
plt.show()

#
# Make a 2-D plot of the gridlines
#
rout   = 10*au       # Zoom-in to the inner edge
# rout   = 100*au    # Uncomment this to see the full grid
fig,ax = plt.subplots()
ax.set_aspect(1)
plot_spher_grid_rtheta(grid,rout,unit=au,linewidth=0.5)
plt.xlabel(r'$r_{\mathrm{cyl}}\;[\mathrm{au}]$')
plt.ylabel(r'$z\;[\mathrm{au}]$')
plt.show()
