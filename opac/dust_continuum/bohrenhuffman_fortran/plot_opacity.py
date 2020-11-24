from plot import *
from radmc3dPy.analyze import *
import os

#
# To make the comparison between the three different Fortran-90
# versions, do the following:
#
#   make
#   ./makeopac
#   ipython --matplotlib
#   %run plot_compare_opacity.py
#


o=readOpac(ext='pyrmg70',scatmat=True)

plt.figure()
plt.loglog(o.wav[0],o.kabs[0],label='Absorption')
plt.loglog(o.wav[0],o.ksca[0],':',label='scattering')
plt.legend()

plt.figure()
ilam = 0
plt.plot(o.scatang[0],o.z11[0][ilam,:])
plt.yscale('log')

plt.show()
