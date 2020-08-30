#
# In version 0.41 of RADMC-3D the dustkappa_silicate.inp was extrapolated
# a bit differently to wavelengths outside of the lab-measured domain.
# Here we compare the new opacities (as of version 2.0 of RADMC-3D) with
# the old ones. At wavelengths <0.2 micron there are differences, and there
# are also minor differences >100 micron. Otherwise they are identical.
#

from radmc3dPy.analyze import *

o  = readOpac(ext=['olivine_amorph_mg50_jaeger94dorschner95'],scatmat=[True])
oo = readOpac('silicate_v0.41')

plt.figure()
plt.plot(o.wav[0],o.kabs[0],label=r'$\kappa_{a}$ (absorption)')
plt.plot(o.wav[0],o.ksca[0],label=r'$\kappa_{s}$ (scattering)')
plt.plot(oo.wav[0],oo.kabs[0],':',label='orig (absorption)')
plt.plot(oo.wav[0],oo.ksca[0],':',label='orig (scattering)')
plt.xscale('log')
plt.yscale('log')
plt.ylim(ymin=1e-2)
plt.xlabel(r'$\lambda [\mu\mathrm{m}]$')
plt.ylabel(r'$\kappa [\mathrm{cm}^2/\mathrm{g}]$')
plt.legend()
plt.show()
