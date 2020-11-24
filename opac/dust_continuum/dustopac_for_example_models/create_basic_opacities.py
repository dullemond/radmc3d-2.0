#
# This Python script uses the Mie algorithm of Bohren & Huffmann (see the bhmie.py
# code) to compute the opacities of:
#
#   pure amorphous pyroxene with 70% Mg   (Jaeger et al. 1994; Dorschner et al. 1995)
#   pure amorphous olivine with 50% Mg    (Jaeger et al. 1994; Dorschner et al. 1995)
#   pure amorphous carbon                 (Preibisch et al. 1993)
#
# All these opacities are for a single grain size. To avoid unrealistic resonant
# features, the single grain size can be replaced by a narrow Gaussian distribution
# of grain sizes. Since the optical constant tables from the laboratory have a limited
# wavelength range, they are extrapolated toward shorter wavelength by assuming them
# constant, and toward longer wavelength by assuming them to be a powerlaw.
#

from makedustopac import *
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import copy
import pickle

agraincm  = 0.1 * 1e-4    # Grain size in cm
logawidth = 0.05          # Smear out the grain size by 5% in both directions
na        = 20            # Use 20 grain size samples
#logawidth = None          # Uncomment this if you do _not_ want a smeared out distribution
#na        = None          # Uncomment this if you do _not_ want a smeared out distribution
chop      = 5.            # Remove forward scattering within an angle of 5 degrees
optconst  = ["pyroxene_amorph_mg70_jaeger94dorschner95",
             "olivine_amorph_mg50_jaeger94dorschner95",
             "carbon_amorph_preibisch93"]                # The optical constants names
descript  = ["Amorphous Pyroxene with 70% Mg and 30% Fe",
             "Amorphous Olivine with 50% Mg and 50% Fe",
             "Amorphous Carbon"]
extrapol  = True          # Extrapolate optical constants beyond its wavelength grid, if necessary
verbose   = False         # If True, then write out status information
ntheta    = 181           # Number of scattering angle sampling points

#
# Set up a wavelength grid upon which we want to compute the opacities
#
lamcm     = 10.0**np.linspace(-1,4,1000)*1e-4

#
# Set up an angular grid for which we want to compute the scattering matrix Z
#
theta     = np.linspace(0.,180.,ntheta)

#-----------------------------------------------------------------------------
# Now create the opacities for each species separately, and then mix.
# This is, of course, a crude mixing method, but in the absence of detailed
# and reliable mixing methods, this is a reasonable approximation.
#-----------------------------------------------------------------------------

opaclist = []

for ispec in range(len(optconst)):
    #
    # Now make the opacity with the bhmie code
    #
    optc         = optconst[ispec]
    optconstfile = optc+'.lnk'
    descr        = descript[ispec]
    print("Running the code. Please wait...")
    opac       = compute_opac_mie(optconstfile,None,agraincm,lamcm,theta=theta,
                              extrapolate=extrapol,logawidth=logawidth,na=na,
                              chopforward=chop,verbose=verbose)
    opaclist.append(copy.deepcopy(opac))
    #
    # Now write it out to a RADMC-3D opacity file:
    #
    # ...The full scattering matrix file
    #
    print("Writing the opacity to scatmat file")
    write_radmc3d_scatmat_file(opac,optc,descr=descr)
    #
    # ...Only the opacity file with simple scattering info
    #    (uncomment the next two commands if you wish to use this)
    #
    print("Writing the opacity to kappa file")
    write_radmc3d_kappa_file(opac,optc,descr=descr)
    #
    # Note that RADMC-3D does not like it when both files are there, so
    # you must choose whether you want to do the full scattering or not
    # (advice: yes, do the full scattering, which is safer as it is more
    # accurate).
    #

#
# Save the opacity structures in case we want to mix using
# mix_opacities.py
#
with open('opacities.pickle','wb') as f:
    pickle.dump(opaclist,f)
