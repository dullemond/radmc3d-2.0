from plot import *
from radmc3dPy.analyze import *
import os

#
# To make the comparison between the three different Fortran-90
# versions, do the following:
#
#   python python_version.py   (only for creating the wavelength_micron.inp file)
#   make
#   ./makeopac
#   mv dustkapscatmat_pyrmg70.inp dustkapscatmat_pyrmg70_fortran_single.inp
#   ./makeopac_smoothed
#   mv dustkapscatmat_pyrmg70.inp dustkapscatmat_pyrmg70_fortran_smoothed.inp
#   ./makeopac_distr
#   mv dustkapscatmat_pyrmg70.inp dustkapscatmat_pyrmg70_fortran_distr.inp
#   ipython --matplotlib
#   %run plot_compare_distr.py
#


osi=readOpac(ext='pyrmg70_fortran_single',scatmat=True)

ods=readOpac(ext='pyrmg70_fortran_distr',scatmat=True)

osm=readOpac(ext='pyrmg70_fortran_smoothed',scatmat=True)

plt.figure()
plt.loglog(osi.wav[0],osi.kabs[0],label='Single')
plt.loglog(osm.wav[0],osm.kabs[0],label='Smoothed')
plt.loglog(ods.wav[0],ods.kabs[0],label='Distribution')
plt.loglog(osi.wav[0],osi.ksca[0],':',label='Single')
plt.loglog(osm.wav[0],osm.ksca[0],':',label='Smoothed')
plt.loglog(ods.wav[0],ods.ksca[0],':',label='Distribution')
plt.legend()

plt.figure()
ilam = 0
plt.plot(osi.scatang[0],osi.z11[0][ilam,:],label='Single')
plt.plot(osm.scatang[0],osm.z11[0][ilam,:],label='Smoothed')
plt.plot(ods.scatang[0],ods.z11[0][ilam,:],label='Distribution')
plt.yscale('log')
plt.legend()

plt.show()
