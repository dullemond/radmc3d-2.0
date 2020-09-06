import plotly.graph_objects as go
import numpy as np
from radmc3dPy.subbox import *
from natconst import *

#
# Example how to use RADMC-3D subbox method to sample an irregular grid
# in a regular way, making it easier to plot an analyze stuff.
#

n      = 100
nxyz   = [n,n,n]

#sz     = 1200*au  # See the full disk
#floor  = 1e-20

sz     = 150*au  # Zoom in
floor  = 1e-15

box  = np.array([-sz,sz,-sz,sz,-sz,sz])

s=subBox()
s.makeSubbox('rhodust',box,nxyz,phi1=0.,theta=0.,phi2=0.)
s.readSubbox()

values = np.log10(s.data+floor)
X,Y,Z  = np.meshgrid(s.x,s.y,s.z,indexing='ij')

#
# Use plotly library to make a volume rendering of the disk
# https://plotly.com/python/3d-volume-plots/
#
fig = go.Figure(data=go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=values.flatten(),
    isomin=values.min(),
    isomax=values.max(),
    opacity=0.1, # needs to be small to see through all surfaces
    surface_count=17, # needs to be a large number for good volume rendering
    ))
fig.show()
