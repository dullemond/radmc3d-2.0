import numpy as np
import os

#
# Function for creating a dust opacity file, using the native RADMC-3D Mie code
# tools, which is based on the Bohren & Huffman code (from the book from the same
# authors), and the corresponding F77 code written by Bruce Draine. 
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

#
# Function for creating a dust opacity file, using the optool code by Carsten Dominik
# and Michiel Min: https://github.com/cdominik/optool
# See paper Min, M., Hovenier, J. W., & De Koter, A. (2005). Modeling optical properties
# of cosmic dust grains using a distribution of hollow spheres. Astronomy and Astrophysics,
# 432(3), 909–920. http://doi.org/10.1051/0004-6361:20041920.
# And also see the DIANA project: https://dianaproject.wp.st-andrews.ac.uk/papers/
# This tool is more powerful than the Mie code native to RADMC-3D, as it includes
# things like ice mantels, mixing of materials, porosity, DHS method etc.
#
# NOTE: Compared to the create_dustkapscatmat_file() function above, the
#       smoothing is done not with a Gaussian, but with a powerlaw distribution
#       within a narrow wavelength range: Hence the amin and amax, instead of a.
#
def create_dustkapscatmat_file_optool(a,amin,amax,species,mfrac=None,nsample=20, \
                               dhs_fmax=0.,porosity=0.,                          \
                               chopangle=5.,renametosize=True,sizeformat="",     \
                               command="optool",scatmat=True,kappa=False):
    """
    Make a dust opacity file using the optool code. Make sure that you 
    installed the 'optool' code beforehand. You can do that by cloning the 
    optool from github:

      git clone https://github.com/cdominik/optool.git

    and then install the code according to the directions. As usual, you 
    must ensure that the shell knows where to find the code. The easiest is
    to specify the entire path to the code by setting the command keyword
    (see below).

    Note that the opacity model is more general than the simple Mie model.
    It is the Distribution of Hollow Spheres (DHS) model by Min, M., 
    Hovenier, J. W., & De Koter, A. (2005) A&A 432, 909 mode. However,
    by default the simplest version of DHS is chosen, which equals the
    Mie scattering model again. You can switch on DHS by setting e.g.
    dhs_fmax=0.8. 

    ARGUMENTS:

      a                 Grain radius in micrometer. Note that this is only
                        used for when you use renametosize=True (default).
                        The more relevant size parameters are amin,amax.
      amin,amax         Grain radius interval in micrometer. The reason why not 
                        just a is specified is that the code will average over a
                        small interval between amin and amax, to smooth out the
                        resonances a bit. Typically you take amin and amax to be
                        reasonably close: e.g. amin=0.9 and amax=1.1 . 
      species           Dust species name. To know which are available, type
                        optool -c in the bash shell. It can also be a file name
                        of an optical constants file. 
                        You can also make species a list or array of species. 
                        In that case you must also specify mfrac, which also
                        must be a list or array.

    OPTIONAL ARGS:

      mfrac             (only if species is a list or array): The mass fraction 
                        of the various dust species (must have the same nr of
                        elements as species). 
      nsample           The dust opacity is computed not for exactly one size, but
                        for a narrow powerlaw distribution of sizes, so that the 
                        worst resonance wiggles are somewhat smeared out. Nsample is
                        the number of size samplings used.
      dhs_fmax          The maximum volume fraction of the central cavity in the 
                        DHS grain model. Default=0 (which equals Mie scattering).
      porosity          The porosity factor used. Default=0.
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
    assert amin<=a, "Amin cannot be larger than a"
    assert amax>=a, "Amax cannot be smaller than a"
    fullcommand = command+' '
    if(type(species)==list):
        assert type(mfrac)==list, "If species is a list, also mfrac must be a list"
        assert len(mfrac)==len(species), "Number of elements of mfrac must be same as that of species"
        for s,f in zip(species,mfrac):
            fullcommand += s+' '+'{}'.format(f)
    else:
        fullcommand += species
    assert not (scatmat and kappa), 'Please set EITHER scatmat=True OR kappa=True'
    assert scatmat or kappa, 'Please set EITHER scatmat=True OR kappa=True'
    fullcommand += ' -a {0} {1} {2}'.format(amin,amax,nsample)
    fullcommand += ' -p {}'.format(porosity)
    fullcommand += ' -fmax {}'.format(dhs_fmax)
    fullcommand += ' -chop {}'.format(chopangle)
    fullcommand += ' -radmc'
    if scatmat: fullcommand += ' -s'
    os.system(fullcommand)
    if renametosize:
        astr      = ("{0:"+sizeformat+"}").format(a)
        if scatmat:
            file_from = "dustkapscatmat.inp"
            file_to   = "dustkapscatmat_"+astr+".inp"
            os.system("mv -f "+file_from+" "+file_to)
        if kappa:
            file_from = "dustkapscatmat.inp"
            file_to   = "dustkapscatmat_"+astr+".inp"
            os.system("mv -f "+file_from+" "+file_to)
        os.system("rm -f dustkappa.inp")
        os.system("rm -f dustkapscatmat.inp")
