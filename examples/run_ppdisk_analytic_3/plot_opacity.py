import numpy as np
from matplotlib import pyplot as plt
from radmc3dPy.analyze import *
from radmc3dPy.natconst import *

#
# Plot the opacity table
#
o1    = readOpac(ext='0.1_micron')
o2    = readOpac(ext='100_micron')
plt.figure()
plt.loglog(o1.wav[0],o1.kabs[0],label=r'$\kappa_\nu^{\mathrm{abs}}$ (absorption)',color='C0')
plt.loglog(o1.wav[0],o1.ksca[0],':',label=r'$\kappa_\nu^{\mathrm{scat}}$ (scattering)',color='C0')
plt.loglog(o2.wav[0],o2.kabs[0],label=r'$\kappa_\nu^{\mathrm{abs}}$ (absorption)',color='C1')
plt.loglog(o2.wav[0],o2.ksca[0],':',label=r'$\kappa_\nu^{\mathrm{scat}}$ (scattering)',color='C1')
plt.ylim((1e-2,1e5))
plt.xlabel(r'$\lambda\;[\mu\mathrm{m}]$')
plt.ylabel(r'$\kappa_\nu\;[\mathrm{cm}^2/\mathrm{g}]$')
plt.title(r'Dust opacity (olivine, $a=0.1$ and $100\,\mu\mathrm{m}$)')
plt.legend()
