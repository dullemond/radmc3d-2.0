from radmc3dPy.image import *
from matplotlib import cm
import os

os.system('radmc3d image iline 1 vkms 10 incl 60 phi 30')
im_nocatch = readImage()

os.system('radmc3d image iline 1 vkms 10 incl 60 phi 30 doppcatch')
im_catch = readImage()

plt.figure(1)
plotImage(im_nocatch,log=True,cmap=cm.hot,bunit='snu',dpc=140,arcsec=True,vmin=-10.5,vmax=-4.0)

plt.figure(2)
plotImage(im_catch,log=True,cmap=cm.hot,bunit='snu',dpc=140,arcsec=True,vmin=-10.5,vmax=-4.0)
