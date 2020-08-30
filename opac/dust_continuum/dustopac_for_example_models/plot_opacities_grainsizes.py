from radmc3dPy.analyze import *

species   = ["0.1",
             "100.0"]

names     = ['0.1 micron','100 micron']


opacs     = []

for s in species:
    opacs.append(readOpac(ext=s,scatmat=True))

plt.figure()
for o,n in zip(opacs,names):
    plt.plot(o.wav[0],o.kabs[0],label=r'$\kappa_{a}$ ('+n+')')
    plt.plot(o.wav[0],o.ksca[0],':',label=r'$\kappa_{s}$ ('+n+')')
plt.xscale('log')
plt.yscale('log')
plt.ylim(ymin=1e-2)
plt.xlabel(r'$\lambda [\mu\mathrm{m}]$')
plt.ylabel(r'$\kappa [\mathrm{cm}^2/\mathrm{g}]$')
plt.legend()

plt.show()
