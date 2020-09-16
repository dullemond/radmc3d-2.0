from makedustopac import *
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

agraincm  = 10 * 1e-4     # Grain size in cm
logawidth = 0.05          # Smear out the grain size by 5% in both directions
na        = 20            # Use 20 grain size samples
chop      = 5.            # Remove forward scattering within an angle of 5 degrees
optconst  = "pyrmg70"     # The optical constants name
descr     = "Amorphous Pyroxene with 70% Mg and 30% Fe"
extrapol  = True          # Extrapolate optical constants beyond its wavelength grid, if necessary
verbose   = False         # If True, then write out status information
ntheta    = 181           # Number of scattering angle sampling points

#
# Set up a wavelength grid upon which we want to compute the opacities
#
#data = np.loadtxt(optconst+'.lnk')
#lammic, ncoef, kcoef = data.T
#lamcm     = lammic*1e-4
#
lammin    = 0.1
lammax    = 1e4
nlam      = 200
lammic    = lammin * (lammax/lammin)**np.linspace(0,1,nlam)
lamcm     = lammic*1e-4
with open('wavelength_micron.inp','w') as f:
    f.write('{}\n'.format(nlam))
    for l in lammic:
        f.write('{0:13.7e}\n'.format(l))

#
# Set up an angular grid for which we want to compute the scattering matrix Z
#
theta     = np.linspace(0.,180.,ntheta)

#-----------------------------------------------------------------------------
#
# Now make the opacity with the bhmie code
#
optconstfile= optconst+'.lnk'
print("Running the code. Please wait...")
opac       = compute_opac_mie(optconstfile,None,agraincm,lamcm,theta=theta,
                              extrapolate=extrapol,logawidth=logawidth,na=na,
                              chopforward=chop,verbose=verbose)
#
# Now write it out to a RADMC-3D opacity file:
#
# ...The full scattering matrix file
#
print("Writing the opacity to scatmat file")
write_radmc3d_scatmat_file(opac,optconst+'_python',descr=descr)
#
# ...Only the opacity file with simple scattering info
#    (uncomment the next two commands if you wish to use this)
#
print("Writing the opacity to kappa file")
write_radmc3d_kappa_file(opac,optconst+'_python',descr=descr)
#
# Note that RADMC-3D does not like it when both files are there, so
# you must choose whether you want to do the full scattering or not
# (advice: yes, do the full scattering, which is safer as it is more
# accurate).
#
