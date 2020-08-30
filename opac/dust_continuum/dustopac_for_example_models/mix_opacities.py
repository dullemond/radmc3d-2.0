from makedustopac import *
import numpy as np
import copy
import pickle

#
# Read the opacity structures created by create_basic_opacities.py
#
with open('opacities.pickle','rb') as f:
    opaclist = pickle.load(f)

#
# Set the mixing ratios
#
abun      = [0.8,0.1,0.1]     # The mixing ratios by mass

#
# Normalize the mixing ratios
#
abun = np.array(abun)
dum = abun.sum()
abun /= dum

#
# Now mix
#
opacmix = copy.deepcopy(opaclist[0])
opacmix["kabs"][:]  = 0.0
opacmix["kscat"][:] = 0.0
opacmix["gscat"][:] = 0.0
opacmix["zscat"][:] = 0.0
if "zscat_nochop" in opacmix: opacmix["zscat_nochop"][:] = 0.0
if "kscat_nochop" in opacmix: opacmix["kscat_nochop"][:] = 0.0
mdav = 0.0
for ispec in range(len(opaclist)):
    opacmix["kabs"][:]  += abun[ispec]*opaclist[ispec]["kabs"]
    opacmix["kscat"][:] += abun[ispec]*opaclist[ispec]["kscat"]
    opacmix["gscat"][:] += abun[ispec]*opaclist[ispec]["gscat"]
    opacmix["zscat"][:] += abun[ispec]*opaclist[ispec]["zscat"]
    if "zscat_nochop" in opacmix: opacmix["zscat_nochop"][:] += abun[ispec]*opaclist[ispec]["zscat_nochop"]
    if "kscat_nochop" in opacmix: opacmix["kscat_nochop"][:] += abun[ispec]*opaclist[ispec]["kscat_nochop"]
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
print("Writing the opacity to scatmat file")
write_radmc3d_scatmat_file(opacmix,optc)
#
# ...Only the opacity file with simple scattering info
#    (uncomment the next two commands if you wish to use this)
#
print("Writing the opacity to kappa file")
write_radmc3d_kappa_file(opacmix,optc)
#
# Now that RADMC-3D does not like it when both files are there, so
# you must choose whether you want to do the full scattering or not
# (advice: yes, do the full scattering, which is safer as it is more
# accurate).
#
    
