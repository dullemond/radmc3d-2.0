import numpy as np
import matplotlib.pyplot as plt
from radmc3dPy.subbox import *
from natconst import *

#
# Same as subbox, but now just 2D slice at some angle
#

n     = 400
nxyz  = [n,n,1]
sz    = 1200*au
box   = np.array([-sz,sz,-sz,sz,-sz/n,sz/n])
phi1  = 0.
theta = 90.  # Vertical plane
phi2  = 0.   # Along the warped direction 
#phi2  = 90.  # Along the unwarped direction

s=subBox()
s.makeSubbox('rhodust',box,nxyz,phi1=phi1,theta=theta,phi2=phi2)
s.readSubbox()

floor = 1e-20
plt.figure()
plt.imshow(np.log10(s.data[:,:,0,0].T+floor),extent=box[0:4]/au)
plt.show()
