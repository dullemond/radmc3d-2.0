import numpy as np
import os
#
# Function for creating a dust opacity file
#
def create_dustkapscatmat_file(a,species,nsample=20,logawidth=0.05,       \
                               errortol=0.1,chopangle=5.,renametosize=True,sizeformat="", \
                               command="./makeopac_smoothed",scatmat=True,kappa=False):
    """
    Make a dust opacity file. Make sure that you type "make" to compile the fortran code
    makeopac_smoothed in this directory.

    ARGUMENTS:

      a                 Grain radius in micrometer
      species           Dust species name (there should be a file <species>.lnk in the directory)

    OPTIONAL ARGS:

      nsample           The dust opacity is computed not for exactly one size, but
                        for a narrow Gaussian distribution of sizes, so that the 
                        worst resonance wiggles are somewhat smeared out. Nsample is
                        the number of size samplings used.
      logawidth         Is the width of the Gaussian smoothing distribution.
      errtol            The relative error tolerance (for very large grains you
                        may need to set this to rather large values, but then you
                        should use "chopangle" below for consistency).
      chopangle         For very strongly forward-peaked scattering, the most 
                        extreme parts of the phase function are limited. This sets
                        the angle where this limitation is done. If 0, no limiting is
                        done, if 5, limiting is done within a forward-peaked cone of
                        5 degrees, and the larger chopangle, the more limiting is 
                        done (the less accurate, thought). 
      renametosize      If True, then rename the output opacity file to a filename 
                        that identifies the size of the grain (see sizeformat below, too)
      sizeformat        The format string how to format the size into a string. 
                        Default = "". Another option could be e.g. "8.2e" or "9.3e".
      command           This is the command to be used to call the Fortran90 Mie code.
      scatmat           If True and if renametosize is True, then copy the 
                        dustkapscatmat_<species>.inp files to dustkapscatmat_<size>.inp
      kappa             If True and if renametosize is True, then copy the 
                        dustkappat_<species>.inp files to dustkappa_<size>.inp

    
    """
    if(not os.path.isfile(command)):
        msg = "The Fortran90 code "+command+" does not exist. Please compile this code first (type 'make')"
        raise RuntimeError(msg)
    with open("param_smoothed.inp",'w') as f:
        f.write("{}     ; Name of the grain species\n".format(species))
        f.write("{}          ; Nr of grain sizes in this distribution\n".format(nsample))
        f.write("{0:10.3e}   ; Mean grain radius in cm\n".format(a*1e-4))
        f.write("{}        ; logawidth\n".format(logawidth))
        f.write("3.0         ; wfact\n")
        f.write("0.0         ; Set to zero, so that the material density is read from the optical constants file\n")
        f.write("181         ; Nr of angle sampling points between 0 and 180\n")
        f.write("{}          ; Error tolerance for testing kappa_scat integral\n".format(errortol))
        f.write("0           ; Keep this to 0 (only for backward compatibility with old version)\n")
        f.write("{}          ; Chopping angle\n".format(chopangle))
        f.write("wavelength_micron.inp\n")
        f.write("1e4         ; Maximum value of x for Mie; beyond that we use a scaling relation\n")
    os.system(command)
    if renametosize:
        astr      = ("{0:"+sizeformat+"}").format(a)
        if scatmat:
            file_from = "dustkapscatmat_"+species+".inp"
            file_to   = "dustkapscatmat_"+astr+".inp"
            os.system("mv -f "+file_from+" "+file_to)
        if kappa:
            file_from = "dustkapscatmat_"+species+".inp"
            file_to   = "dustkapscatmat_"+astr+".inp"
            os.system("mv -f "+file_from+" "+file_to)
        os.system("rm -f dustkappa_"+species+".inp")
        os.system("rm -f dustkapscatmat_"+species+".inp")
