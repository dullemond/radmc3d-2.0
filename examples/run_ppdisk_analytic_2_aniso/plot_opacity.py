import numpy as np
from matplotlib import pyplot as plt
from radmc3dPy.analyze import *
from radmc3dPy.natconst import *

#
# Plot the opacity table
#
o    = readOpac(ext='silicate',scatmat=True)
plt.figure()
plt.loglog(o.wav[0],o.kabs[0],label=r'$\kappa_\nu^{\mathrm{abs}}$ (absorption)')
plt.loglog(o.wav[0],o.ksca[0],':',label=r'$\kappa_\nu^{\mathrm{scat}}$ (scattering)')
plt.ylim((1e-2,1e5))
plt.xlabel(r'$\lambda\;[\mu\mathrm{m}]$')
plt.ylabel(r'$\kappa_\nu\;[\mathrm{cm}^2/\mathrm{g}]$')
plt.title(r'Dust opacity (olivine, $a=0.1\,\mu\mathrm{m}$)')
plt.legend()

#
# Plot the phase function at 0.2, 0.5 and 1.0 micron
#
lam  = np.array([0.2,0.5,1.0])   # Wavelengths where to plot the phase function
ilam = np.array(np.interp(lam,o.wav[0],np.arange(len(o.wav[0])))+0.5,dtype=int)  # Nearest radial grid point
lamn = [r'$\lambda=0.2\,\mu\mathrm{m}$',r'$\lambda=0.5\,\mu\mathrm{m}$',r'$\lambda=1\,\mu\mathrm{m}$']
plt.figure()
for i in range(len(ilam)):
    plt.plot(o.scatang[0],4*np.pi*o.z11[0][ilam[i]]/o.ksca[0][ilam[i]],label=lamn[i])
plt.yscale('log')
plt.xlabel(r'$\theta [\mathrm{deg}]$')
plt.ylabel(r'$4\pi\,Z_{11}/\kappa_{s}$')
plt.title('Scattering phase function of olivine (a=0.1 micron)')
plt.legend()
plt.show()
