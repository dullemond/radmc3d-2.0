import numpy as np
from matplotlib import pyplot as plt
from radmc3dPy.analyze import *
from radmc3dPy.natconst import *

#
# Plot the opacity table
#
o    = readOpac(ext='silicate')
plt.figure()
plt.loglog(o.wav[0],o.kabs[0],label=r'$\kappa_\nu^{\mathrm{abs}}$ (absorption)')
plt.loglog(o.wav[0],o.ksca[0],':',label=r'$\kappa_\nu^{\mathrm{scat}}$ (scattering)')
plt.ylim((1e-2,1e5))
plt.xlabel(r'$\lambda\;[\mu\mathrm{m}]$')
plt.ylabel(r'$\kappa_\nu\;[\mathrm{cm}^2/\mathrm{g}]$')
plt.title(r'Dust opacity (olivine, $a=0.1\,\mu\mathrm{m}$)')
plt.legend()
