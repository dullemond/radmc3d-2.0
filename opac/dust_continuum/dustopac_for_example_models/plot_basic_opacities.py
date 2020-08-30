from radmc3dPy.analyze import *

species   = ["pyroxene_amorph_mg70_jaeger94dorschner95",
             "olivine_amorph_mg50_jaeger94dorschner95",
             "carbon_amorph_preibisch93"]

names     = ['Pyroxene','Olivine','Carbon']


opacs     = []

for s in species:
    opacs.append(readOpac(ext=s,scatmat=True))

plt.figure()
for o,n in zip(opacs,names):
    plt.plot(o.wav[0],o.kabs[0],label=r'$\kappa_{a}$ ('+n+')')
    plt.plot(o.wav[0],o.ksca[0],'-',label=r'$\kappa_{s}$ ('+n+')')
plt.xscale('log')
plt.yscale('log')
plt.ylim(ymin=1e-2)
plt.xlabel(r'$\lambda [\mu\mathrm{m}]$')
plt.ylabel(r'$\kappa [\mathrm{cm}^2/\mathrm{g}]$')
plt.legend()

i = 1

o    = opacs[i]
n    = names[i]
ilam = [50,74,100]
lamn = [r'$\lambda=1\,\mu\mathrm{m}$',r'$\lambda=3\,\mu\mathrm{m}$',r'$\lambda=10\,\mu\mathrm{m}$']
plt.figure()
for i in range(len(ilam)):
    plt.plot(o.scatang[0],4*np.pi*o.z11[0][ilam[i]]/o.ksca[0][ilam[i]],label=lamn[i])
plt.yscale('log')
plt.xlabel(r'$\theta [\mathrm{deg}]$')
plt.ylabel(r'$4\pi\,Z_{11}/\kappa_{s}$')
plt.title('Phase function for '+n)
plt.legend()

plt.show()
