import numpy as np
from makedustopacfortran import *

n_grains        = 7          # How many separate grain size opacity files to make
a_gr_mic_min    = 1e-1       # Smallest grain size in micron
a_gr_mic_max    = 1e5        # Largest grain size in micron
sizeformat      = "8.2e"     # Format of the grain size in the filenames
species         = "olivine"  # Name of the optical constants file for olivine: olivine.lnk
errortol        = 1e99       # The errors come from super-forward-scattering. With chopping we fix that.
chopangle       = 5.         # The chopping angle cone.
lammin          = 0.1        # Smallest wavelength in micron
lammax          = 1e4        # Largest wavelength in micron
nlam            = 200        # Nr of wavelength sampling points 

#
# Create the wavelength grid for the opacity files
#
lammic    = lammin * (lammax/lammin)**np.linspace(0,1,nlam)
lamcm     = lammic*1e-4
with open('wavelength_micron.inp','w') as f:
    f.write('{}\n'.format(nlam))
    for l in lammic:
        f.write('{0:13.7e}\n'.format(l))

#
# Creat the grain sizes
#
a_grains_mic    = a_gr_mic_min * (a_gr_mic_max/a_gr_mic_min)**np.linspace(0,1,n_grains)

with open("dustspecies.inp","w") as f:
    f.write(species+"\n")

with open("sizes.inp","w") as f:
    for a in a_grains_mic:
        astr = ("{0:"+sizeformat+"}").format(a)
        print("Creating opacity for grain size "+astr+" micron radius")
        create_dustkapscatmat_file(a,species,sizeformat=sizeformat, \
                                   errortol=errortol,chopangle=chopangle)
        f.write("{}\n".format(astr))
