from makedustopac import *
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

agraincm  = 10 * 1e-4     # Grain size in cm
logawidth = 0.05          # Smear out the grain size by 5% in both directions
na        = 20            # Use 20 grain size samples
chop      = 5.            # Remove forward scattering within an angle of 5 degrees
optconst  = "pyrmg70"     # The optical constants name
matdens   = 3.0           # The material density in gram / cm^3
extrapol  = True          # Extrapolate optical constants beyond its wavelength grid, if necessary
verbose   = True          # If True, then write out status information
ntheta    = 181           # Number of scattering angle sampling points

#
# Set up a wavelength grid upon which we want to compute the opacities
#
lamcm     = 10.0**np.linspace(-1,3,200)*1e-4

#
# Set up an angular grid for which we want to compute the scattering matrix Z
#
theta     = np.linspace(0.,180.,ntheta)

#-----------------------------------------------------------------------------
#
# Now make the opacity with the bhmie code
#
optconstfile= optconst+'.lnk'
#
# ...First with a (narrow) Gaussian grain size distribution to 
#    suppress the resonances
#
opac       = compute_opac_mie(optconstfile,matdens,agraincm,lamcm,theta=theta,
                              extrapolate=extrapol,logawidth=logawidth,na=na,
                              chopforward=chop,verbose=verbose)
#
# ...Then for comparison a single grain size, i.e. with the resonances
#
opacsingle = compute_opac_mie(optconstfile,matdens,agraincm,lamcm,theta=theta,
                              extrapolate=extrapol,
                              chopforward=chop,verbose=verbose)

#-----------------------------------------------------------------------------
#
# Now make some plots
#

fig1 = plt.figure(1)
plt.figure(1)
plt.cla()
plt.plot(opac["lamcm"]*1e4,opac["kabs"],label='Absorption')
plt.plot(opac["lamcm"]*1e4,opac["kscat"],label='Scattering')
plt.plot(opacsingle["lamcm"]*1e4,opacsingle["kscat"],label='Scattering (no smear, yes chop)')
plt.plot(opac["lamcm"]*1e4,opac["kscat_nochop"],label='Scattering (yes smear, no chop)')
plt.xlabel(r'$\lambda\; [\mu\mathrm{m}]$')
plt.ylabel(r'$\kappa\; [\mathrm{cm}^2/\mathrm{g}]$')
plt.xscale("log")
plt.yscale("log")
plt.title('Absorption and scattering opacity')
plt.legend(loc=3)
plt.show()

fig2 = plt.figure(2)
plt.figure(2)
plt.cla()
inu = 50
ll  = lamcm[inu]*1e4
plt.plot(opac["theta"],opac["zscat"][inu,:,0],label='With smearing, yes chop')
plt.plot(opacsingle["theta"],opacsingle["zscat"][inu,:,0],label='No smearing, yes chop')
plt.plot(opacsingle["theta"],opac["zscat_nochop"][inu,:,0],label='With smearing, no chop')
plt.yscale("log")
plt.xlabel(r'$\theta\; [\mathrm{degrees}]$')
plt.ylabel(r'$Z_{11}\; [\mathrm{cm}^2/\mathrm{g.ster}]$')
plt.title(r'Scattering phase function at $\lambda='+'%13.6f'%(ll)+'\,\mu\mathrm{m}$')
plt.legend()
plt.show()

fig2 = plt.figure(3)
plt.figure(3)
plt.cla()
plt.plot(opac["agr"]*1e4,opac["wgt"]/opac["wgt"].max())
plt.xlabel(r'$a_{\mathrm{grain}}\;[\mu\mathrm{m}]$')
plt.ylabel('Contribution (normalized)')
plt.title('Grain size distribution used for smearing')
plt.show()
