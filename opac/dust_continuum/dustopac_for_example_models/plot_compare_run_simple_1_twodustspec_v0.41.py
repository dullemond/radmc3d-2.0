#
# The run_simple_1_twodustspec/ model in version 0.41 of RADMC-3D
# had a different amorphous carbon opacity. The amorphous silicate
# opacity is (nearly) the same.
#

from radmc3dPy.analyze import *

o  = readOpac(ext=['olivine_amorph_mg50_jaeger94dorschner95','carbon_amorph_preibisch93'],scatmat=[True,True])
oo = readOpac(ext=['silicate_v0.41','carbon_v0.41'])

plt.figure()
plt.plot(o.wav[0],o.kabs[0],label=r'$\kappa_{a}$ (Silicate)')
plt.plot(oo.wav[0],oo.kabs[0],'--',label='orig (Silicate)')
plt.plot(o.wav[1],o.kabs[1],label=r'$\kappa_{a}$ (Carbon)')
plt.plot(oo.wav[1],oo.kabs[1],'--',label='orig (Carbon)')
plt.xscale('log')
plt.yscale('log')
plt.ylim(ymin=1e-2)
plt.xlabel(r'$\lambda [\mu\mathrm{m}]$')
plt.ylabel(r'$\kappa [\mathrm{cm}^2/\mathrm{g}]$')
plt.legend()
plt.title('Model run_simple_1_twodustspec: New opacities')
plt.show()
