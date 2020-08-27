from makedustopac import *
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import copy

agraincm  = 0.1 * 1e-4    # Grain size in cm
logawidth = 0.05          # Smear out the grain size by 5% in both directions
na        = 20            # Use 20 grain size samples
chop      = 5.            # Remove forward scattering within an angle of 5 degrees
optconst  = ["pyroxene_amorph_mg70_jaeger94","carbon_amorph_preibisch93"]   # The optical constants names
matdens   = [3.0,1.8]     # The material density in gram / cm^3
abun      = [0.8,0.2]     # The mixing ratios by mass
extrapol  = True          # Extrapolate optical constants beyond its wavelength grid, if necessary
verbose   = False         # If True, then write out status information
ntheta    = 181           # Number of scattering angle sampling points

#
# Set up a wavelength grid upon which we want to compute the opacities
#
lamcm     = 10.0**np.linspace(-1,3,1000)*1e-4

#
# Set up an angular grid for which we want to compute the scattering matrix Z
#
theta     = np.linspace(0.,180.,ntheta)

#
# Normalize the mixing ratios
#
abun = np.array(abun)
dum = abun.sum()
abun /= dum

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
    dens         = matdens[ispec]
    print "Running the code. Please wait..."
    opac       = compute_opac_mie(optconstfile,dens,agraincm,lamcm,theta=theta,
                              extrapolate=extrapol,logawidth=logawidth,na=na,
                              chopforward=chop,verbose=verbose)
    opaclist.append(copy.deepcopy(opac))
    #
    # Now write it out to a RADMC-3D opacity file:
    #
    # ...The full scattering matrix file
    #
    print "Writing the opacity to scatmat file"
    write_radmc3d_scatmat_file(opac,optc)
    #
    # ...Only the opacity file with simple scattering info
    #    (uncomment the next two commands if you wish to use this)
    #
    #print "Writing the opacity to kappa file"
    #write_radmc3d_kappa_file(opac,optc)
    #
    # Now that RADMC-3D does not like it when both files are there, so
    # you must choose whether you want to do the full scattering or not
    # (advice: yes, do the full scattering, which is safer as it is more
    # accurate).
    #

#
# Now mix
#
opacmix = copy.deepcopy(opac)
opacmix["kabs"][:]  = 0.0
opacmix["kscat"][:] = 0.0
opacmix["gscat"][:] = 0.0
opacmix["zscat"][:] = 0.0
if "zscat_nochop" in opacmix: opac["zscat_nochop"][:] = 0.0
if "kscat_nochop" in opacmix: opac["kscat_nochop"][:] = 0.0
mdav = 0.0
for ispec in range(len(optconst)):
    opacmix["kabs"][:]  += abun[ispec]*opaclist[ispec]["kabs"]
    opacmix["kscat"][:] += abun[ispec]*opaclist[ispec]["kscat"]
    opacmix["gscat"][:] += abun[ispec]*opaclist[ispec]["gscat"]
    opacmix["zscat"][:] += abun[ispec]*opaclist[ispec]["zscat"]
    if "zscat_nochop" in opacmix: opac["zscat_nochop"][:] += abun[ispec]*opaclist[ispec]["zscat_nochop"]
    if "kscat_nochop" in opacmix: opac["kscat_nochop"][:] += abun[ispec]*opaclist[ispec]["kscat_nochop"]
    mdav += abun[ispec]/opaclist[ispec]["matdens"]
mdav = 1./mdav
opacmix["matdens"] = mdav
#
# Now write it out to a RADMC-3D opacity file:
#
optc = "mix"
#
# ...The full scattering matrix file
#
print "Writing the opacity to scatmat file"
write_radmc3d_scatmat_file(opacmix,optc)
#
# ...Only the opacity file with simple scattering info
#    (uncomment the next two commands if you wish to use this)
#
#print "Writing the opacity to kappa file"
#write_radmc3d_kappa_file(opacmix,optc)
#
# Now that RADMC-3D does not like it when both files are there, so
# you must choose whether you want to do the full scattering or not
# (advice: yes, do the full scattering, which is safer as it is more
# accurate).
#
    
