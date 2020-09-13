import numpy as np
from makedustopacfortran import *

n_grains        = 20
a_gr_mic_min    = 0.1
a_gr_mic_max    = 1e3
sizeformat      = "8.2e"
species         = "olivine"  # Name of the optical constants file for olivine: olivine.lnk
materialdensity = 3.71       # Material density of olivine in g/cm^3
errortol        = 1e7        # The errors come from super-forward-scattering. With chopping we fix that.
chopangle       = 5.         # The chopping angle cone.

a_grains_mic    = a_gr_mic_min * (a_gr_mic_max/a_gr_mic_min)**np.linspace(0,1,n_grains)

with open("dustproperties.inp","w") as f:
    f.write(species+"\n")
    f.write("{}".format(materialdensity)+"\n")

with open("sizes.inp","w") as f:
    for a in a_grains_mic:
        astr = ("{0:"+sizeformat+"}").format(a)
        print("Creating opacity for grain size "+astr+" micron radius")
        create_dustkapscatmat_file(a,species,materialdensity,sizeformat=sizeformat, \
                                   errortol=errortol,chopangle=chopangle)
        f.write("{}\n".format(astr))
