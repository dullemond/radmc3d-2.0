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
verbose   = False         # If True, then write out status information
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
print "Running the code. Please wait..."
opac       = compute_opac_mie(optconstfile,matdens,agraincm,lamcm,theta=theta,
                              extrapolate=extrapol,logawidth=logawidth,na=na,
                              chopforward=chop,verbose=verbose)
#
# Now write it out to a RADMC-3D opacity file:
#
# ...The full scattering matrix file
#
print "Writing the opacity to scatmat file"
write_radmc3d_scatmat_file(opac,optconst)
#
# ...Only the opacity file with simple scattering info
#    (uncomment the next two commands if you wish to use this)
#
#print "Writing the opacity to kappa file"
#write_radmc3d_kappa_file(opac,optconst)
#
# Now that RADMC-3D does not like it when both files are there, so
# you must choose whether you want to do the full scattering or not
# (advice: yes, do the full scattering, which is safer as it is more
# accurate).
#
