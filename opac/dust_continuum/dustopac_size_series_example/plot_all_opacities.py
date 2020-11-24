import numpy as np
import matplotlib.pyplot as plt
from radmc3dPy.analyze import readOpac

with open("dustspecies.inp","r") as f:
    species = f.readline()

a_grains_string = []
with open("sizes.inp","r") as f:
    for s in f:
        a_grains_string.append(s.strip())
a_grains = np.array(a_grains_string,dtype="float")
a_grains_cm   = a_grains*1e-4

opacs = []
for astr in a_grains_string:
    o = readOpac(ext=astr,scatmat=True)
    opacs.append(o)

#
# Get the material density (same for all sizes)
#
materialdensity = o.materialdensity[0]

#
# Compute the geometric opacities of these grains
#
m_grains_g    = (4*np.pi/3.)*materialdensity*a_grains_cm**3
s_grains_cm2  = np.pi*a_grains_cm**2       # Geometric cross section of the grains
k_geom_grains = s_grains_cm2/m_grains_g    # Geometric opacity of the grains

#
# Plot the absorption opacities
#
plt.figure()
for astr,o in zip(a_grains_string,opacs):
    plt.loglog(o.wav[0],o.kabs[0],label='a = '+astr+' mu')
plt.xlabel(r'$\lambda\;[\mu\mathrm{m}]$')
plt.ylabel(r'$\kappa_\nu\;[\mathrm{cm}^2/\mathrm{g}]$')
plt.ylim((1e-3,1e5))
plt.legend()

#
# Plot the Q values (cross section in units of pi*a^2)
#
plt.figure()
for astr,o,k in zip(a_grains_string,opacs,k_geom_grains):
    plt.loglog(o.wav[0],(o.kabs[0]+o.ksca[0])/k,label='a = '+astr+' mu')
plt.xlabel(r'$\lambda\;[\mu\mathrm{m}]$')
plt.ylabel(r'$Q_\nu$')
plt.ylim((1e-5,10))
plt.legend()
plt.title('Total opacity in units of geometric opacity')

#
# Plot the albedos
#
plt.figure()
for astr,o in zip(a_grains_string,opacs):
    plt.loglog(o.wav[0],o.ksca[0]/(o.kabs[0]+o.ksca[0]),label='a = '+astr+' mu')
plt.xlabel(r'$\lambda\;[\mu\mathrm{m}]$')
plt.ylabel('albedo')
plt.ylim((2e-8,2))
plt.legend()
plt.title('Albedo as a function of wavelength')

#
# Plot the albedos as a function of grain size
#
albedos  = np.zeros((len(o.wav[0]),len(a_grains_string)))
for i in range(len(a_grains_string)):
    albedos[:,i] = opacs[i].ksca[0]/(opacs[i].kabs[0]+opacs[i].ksca[0])
lams     = np.array([0.25,0.5,1.0,2.0,4.0])
ilams    = np.array(np.interp(lams,o.wav[0],np.arange(len(o.wav[0]))),dtype=int)
lams     = o.wav[0][ilams]

plt.figure()
for il,l in zip(ilams,lams):
    lamstr = "{}".format(l)
    plt.loglog(a_grains,albedos[il,:],label=r'$\lambda = '+lamstr+'\;\mu\mathrm{m}$')
plt.xlabel(r'$a\;[\mu\mathrm{m}]$')
plt.ylabel('albebo')
plt.title('Albedo as a function of grain size')
plt.legend()

plt.show()
