#
# This Python script uses the Mie algorithm of Bohren & Huffmann (see the bhmie.py
# code) to compute the opacities of:
#
#   pure amorphous olivine with 50% Mg    (Jaeger et al. 1994; Dorschner et al. 1995)
#
# At multiple grain sizes.
#

from makedustopac import *
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import copy
import pickle

agraincm  = np.array([0.1,100]) * 1e-4    # Grain size in cm
logawidth = 0.05          # Smear out the grain size by 5% in both directions
na        = 20            # Use 20 grain size samples
#logawidth = None          # Uncomment this if you do _not_ want a smeared out distribution
#na        = None          # Uncomment this if you do _not_ want a smeared out distribution
chop      = 5.            # Remove forward scattering within an angle of 5 degrees
optconst  = "olivine_amorph_mg50_jaeger94dorschner95"
descript  = "Amorphous Olivine with 50% Mg and 50% Fe"
extrapol  = True          # Extrapolate optical constants beyond its wavelength grid, if necessary
verbose   = True          # If True, then write out status information
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
# Now create the opacities for each grain size separately
#-----------------------------------------------------------------------------

opaclist = []

for isize in range(len(agraincm)):
    agr          = agraincm[isize]
    #
    # Now make the opacity with the bhmie code
    #
    sz           = '{}'.format(agraincm[isize]/1e-4)
    optc         = optconst
    optconstfile = optc+'.lnk'
    descr        = descript
    print("Running the code. Please wait...")
    opac       = compute_opac_mie(optconstfile,None,agr,lamcm,theta=theta,
                              extrapolate=extrapol,logawidth=logawidth,na=na,
                              chopforward=chop,verbose=verbose)
    opaclist.append(copy.deepcopy(opac))
    #
    # Now write it out to a RADMC-3D opacity file:
    #
    # ...The full scattering matrix file
    #
    print("Writing the opacity to scatmat file")
    write_radmc3d_scatmat_file(opac,sz,descr=descr)
    #
    # ...Only the opacity file with simple scattering info
    #    (uncomment the next two commands if you wish to use this)
    #
    print("Writing the opacity to kappa file")
    write_radmc3d_kappa_file(opac,sz,descr=descr)
    #
    # Note that RADMC-3D does not like it when both files are there, so
    # you must choose whether you want to do the full scattering or not
    # (advice: yes, do the full scattering, which is safer as it is more
    # accurate).
    #
