from plot import *                                                                                                                    
from radmc3dPy.analyze import *                                                                                                       
import os                                                                                                                             

#
# To make the comparison between the Fortran-90 and Python versions,
# do the following:
#
#   python python_version.py
#   make
#   ./makeopac_smoothed
#   mv dustkapscatmat_pyrmg70.inp dustkapscatmat_pyrmg70_fortran.inp
#   ipython --matplotlib
#   %run plot_compare_python.py
#


of=readOpac(ext='pyrmg70_fortran',scatmat=True)

op=readOpac(ext='pyrmg70_python',scatmat=True)                                                                                                

plt.figure()
plt.loglog(of.wav[0],of.kabs[0],label='F90')
plt.loglog(op.wav[0],op.kabs[0],label='Python')
plt.loglog(of.wav[0],of.ksca[0],':',label='F90')
plt.loglog(op.wav[0],op.ksca[0],':',label='Python')
plt.legend()

plt.figure()
ilam = 0
plt.plot(of.scatang[0],of.z11[0][ilam,:],label='F90')
plt.plot(op.scatang[0],op.z11[0][ilam,:],label='Python')
plt.yscale('log')
plt.legend()

plt.show()
