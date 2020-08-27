from radmc3d_tools.simpleread import *
from radmc3d_tools.simpleplot import *
import os

# An image (without the star) pole-on, so that the grid structure is visible
os.system('radmc3d image nostar incl 0 phi 0 lambda 1000 npix 200 zoomau -10 10 -10 10')

im = read_image()

image(np.log10(im.image[:,:,0]+1e-15),cmap=cm.rainbow)

plt.show()
