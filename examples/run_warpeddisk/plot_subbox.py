import numpy as np
import matplotlib.pyplot as plt
from radmc3dPy.subbox import *
from natconst import *

#
# Example how to use RADMC-3D subbox method to sample an irregular grid
# in a regular way, making it easier to plot an analyze stuff.
#

n    = 100
nxyz = [n,n,n]
sz   = 1200*au
box  = np.array([-sz,sz,-sz,sz,-sz,sz])

s=subBox()
s.makeSubbox('rhodust',box,nxyz,phi1=0.,theta=0.,phi2=0.)
s.readSubbox()

floor = 1e-20
plt.figure()
#plt.imshow(np.log10(s.data[:,:,n//2,0].T+floor),extent=box[0:4]/au)
#plt.imshow(np.log10(s.data[n//2,:,:,0].T+floor),extent=box[2:6]/au)
plt.imshow(np.log10(s.data[:,n//2,:,0].T+floor),extent=[box[0]/au,box[1]/au,box[4]/au,box[5]/au])
plt.show()
