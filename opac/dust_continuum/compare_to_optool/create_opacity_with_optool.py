from radmc3dPy.analyze import *
from makedustopacfortran import *
import os

#
# Make sure to first install optool by cloning:
#
#   git clone https://github.com/cdominik/optool.git
#
# And then put this directory into your shell PATH.
# Or specify PATH here

amic      = 10
chopangle = 5.
errortol  = 100

path      = '/Users/cornelisdullemond/science/software/optool/'
command   = path+'optool pyr-mg70 -a {} -s -p 0 -fmax 0 -radmc -chop {}'.format(amic,chopangle)
#command   = path+'optool -c pyrmg70.lnk 1.0 3.01 -a {} -s -p 0 -fmax 0 -radmc -chop {}'.format(amic,chopangle)
os.system(command)
os.system('mv -f dustkapscatmat.inp dustkapscatmat_pyrmg70_optool.inp')

#
# Now create the same opacity with the RADMC-3D BHMie code
#
optool = readOpac(ext='pyrmg70_optool',scatmat=True)
with open('wavelength_micron.inp','w') as f:
    f.write('{}\n'.format(len(optool.wav[0])))
    for w in optool.wav[0]:
        f.write('{0:13.6e}\n'.format(w))
create_dustkapscatmat_file(amic,'pyrmg70',errortol=errortol,chopangle=chopangle,renametosize=False,command="./makeopac")
os.system('mv -f dustkapscatmat_pyrmg70.inp dustkapscatmat_pyrmg70_fortran.inp')

#
# Now compare
#
of=readOpac(ext='pyrmg70_fortran',scatmat=True)

op=readOpac(ext='pyrmg70_optool',scatmat=True)                                                                                                

plt.figure()
plt.loglog(of.wav[0],of.kabs[0],label='F90')
plt.loglog(op.wav[0],op.kabs[0],label='Optool')
plt.loglog(of.wav[0],of.ksca[0],':',label='F90')
plt.loglog(op.wav[0],op.ksca[0],':',label='Optool')
plt.legend()

plt.figure()
plt.loglog(of.wav[0],of.phase_g[0],label='F90')
plt.loglog(op.wav[0],op.phase_g[0],label='Optool')
plt.legend()

ilam = 100
plt.figure()
plt.plot(of.scatang[0],of.z11[0][ilam,:],label='F90')
plt.plot(op.scatang[0],op.z11[0][ilam,:],label='Optool')
plt.yscale('log')
plt.legend()

plt.show()
