"""
This module contains functions to read and write input/output data for RADMC-3D and
to do some simple analysis/diagnostics of the model.

Originally from the analyze.pro IDL package by Cornelis Dullemond,
translated into Python and subsequently improved by Attila Juhasz 
(between 2011-2018).
"""
from __future__ import absolute_import
from __future__ import print_function
import traceback
from multiprocessing import Pool
from functools import partial

try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use the python module of RADMC-3D you need to install Numpy')
    print(traceback.format_exc())

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None
    print('Warning')
    print('matplotlib.pyplot cannot be imported')
    print('Without matplotlib you can use the python module to set up a model but you will not be able to plot things')
    print('or display images')

from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import matplotlib.lines as ml

from . natconst import *

from . dustopac import *
from . radsources import *
from . params import *
from . data import *
from . octree import *
from . reggrid import *
from . molecule import *


def readData(ddens=False, dtemp=False, gdens=False, gtemp=False, gvel=False, ispec=None, vturb=False, grid=None,
             old=False, binary=None):
    """Reads the physical variables of the model (e.g. density, velocity, temperature).

    Parameters
    ----------

    ddens : bool
            If True dust density will be read (all dust species and grain sizes)

    dtemp : bool
            If True dust temperature will be read (all dust species and grain sizes)

    gdens : bool
            If True gas density will be read (NOTE: the gas density will be number density in 1/cm^3)

    gtemp : bool
            If True gas temperature will be read (all dust species and grain sizes)

    gvel  : bool
            If True the velocity field will be read

    ispec : str
            Name of the molecule in the 'molecule_ispec.inp' filename

    vturb : bool
            If True the microturbulent velocity field will be read

    grid  : radmc3dGrid
            An instance of radmc3dGrid containing the spatial and frequency grid of the model. If the grid
            is passed to the function it will not be read again from file. This can be useful for octree
            models to save time. 

    old   : bool, optional
            If set to True the file format of the previous, 2D version of radmc will be used

    Returns
    -------
    Returns an instance of the radmc3dData class 
    """
    if binary is not None: print('Note: keyword binary is depricated (binary file format is now automatically recognized)')
    if grid is not None:
        res = radmc3dData(grid=grid)
    else:
        res = radmc3dData()

    # By default: read everything
    if not(ddens or dtemp or gvel or gtemp or vturb or gdens):
        ddens=True
        dtemp=True
        gvel=True
        gtemp=True
        vturb=True
        gdens=True
        
    if ddens:
        res.readDustDens(old=old)
    if dtemp:
        res.readDustTemp(old=old)
    if gvel:
        res.readGasVel()
    if gtemp:
        res.readGasTemp()
    if vturb:
        res.readVTurb()
    if gdens:
        if ispec is not None:
            res.readGasDens(ispec=ispec)

    return res


def readStars(fname=''):
    """
    Reads the data (mass, radius, temperature, spectrum) of discrete stellar sources

    Parameters
    ----------
    fname       : str
                  Name of the file to be read (if omitted the default value is stars.inp)

    Returns
    -------
    An instance of radmc3dRadSources containing the stellar data
    """
    if fname == '':
        fname = 'stars.inp'

    res = radmc3dRadSources()
    res.readStarsinp(fname=fname)

    return res


def readOpac(ext=None, idust=None, scatmat=None, old=False):
    """Reads the dust opacity files.
    This function is an interface to radmc3dDustOpac.readOpac()

    Parameters
    ----------
    ext   : list
            Each element of the list is be a string, the file name extension
            (file names should look like 'dustkappa_ext.inp')

    idust : list
            Each element of the list is an integer, the index of the dust species in the master opacity file
            (dustopac.inp')

    scatmat: list
            If specified, its elements should be booleans indicating whether the opacity file
            contains also the full scattering matrix (True) or only dust opacities (False)

    old   : bool, optional
            If set to True the file format of the previous, 2D version of radmc will be used

    Returns
    -------
        Returns an instance of the radmc3dDustOpac class
    """

    res = radmc3dDustOpac()
    res.readOpac(ext=ext, idust=idust, scatmat=scatmat, old=old)

    return res


def readParams():
    """Reads the problem_params.inp file.
    This function is an interface to radmc3dPar.readPar().

    Returns
    -------
    Returns an instance of the radmc3dPar class 

    """

    dum = radmc3dPar()
    dum.readPar()
    return dum


def writeDefaultParfile(model='', fname=''):
    """Writes a parameter file (problem_params.inp) with default parameters for a given model.

    Parameters
    ----------

    model : str
            Name of the model whose parameter should be written to the file

    fname : str, optional
            Name of the parameter file to be written (if omitted problem_params.inp will be used)
    """

    if model == '':
        raise ValueError('Unknown model. \n No model name is given. ')

    dum = radmc3dPar()
    dum.loadDefaults(model=model)
    dum.writeParfile(fname=fname)


def readSpectrum(fname='', old=False):
    """Reads the spectrum / SED


    Parameters
    -----------
    fname : str, optional
            Name of the file to be read

    old   : bool, optional
            If set to True the file format of the previous, 2D version of radmc will be used


    Returns
    -------

        Returns an ndarray with [Nwavelength, 2] dimensions 
        [Nwavelength,0] is the wavelength / velocity and
        [Nwavelength,1] is the flux density

    """

    if not old:
        if fname.strip() == '':
            fname = 'spectrum.out'

        with open(fname, 'r') as rfile:
            # Read the format number
            dum = rfile.readline()
            # Read the number of wavelengths
            nwav = int(rfile.readline())
            # Read a blank line
            dum = rfile.readline()

            res = np.zeros([nwav, 2], dtype=np.float64)
            for iwav in range(nwav):
                dum = rfile.readline().split()
                res[iwav, 0] = float(dum[0])
                res[iwav, 1] = float(dum[1])

    else:
        if fname.strip() == '':
            fname = 'spectrum.dat'

        with open(fname, 'r') as rfile:
            # Read the number of wavelengths
            nwav = int(rfile.readline())
            rfile.readline()

            res = np.zeros([nwav, 2], dtype=float)
            for iwav in range(nwav):
                dum = rfile.readline().split()
                res[iwav, 0] = cc / float(dum[0]) * 1e4
                res[iwav, 1] = float(dum[1])

    return res


def getDensVstruct(data=None, vmean_temp=False, ispec_tgas=0, gsize=None, idust=None, mstar=None, mu=None):
    """Calculates the vertical hydrostatic equilibrium

    Parameters
    ----------
    data        : radmc3dData
                  An instance of the radmc3DData class containing the density structure of the model

    vmean_temp  : bool
                  If True (T(z) = T(-z) = 0.5*(T(z) + T(-z))) if False (T(z)!=T(-z)) 

    idust       : list
                  List of dust indices whose structure must be calculated

    mstar       : float
                  Stellar mass

    ispec_tgas  : int
                  Index of dust species whose temperature is taken to be the gas temperature

    gsize       : ndarray, optional
                  Dust grain sizes - If specified, the gas temperature is calculated as the average temperature
                  of all dust grains in the grid cell weighted by the total surface area of dust grains with given
                  size - NOTE: this approach assumes that all dust grains of a given size have the same bulk density

    mu          : float, optional
                  Mean molecular weight (default: 2.3)
    Returns
    -------
    Returns an ndarray with the dust density
    """

    if data.grid.crd_sys != 'sph':
        msg = 'Vertical hydrostatic equlibrium structure iteration has been implemented for spherical grids only ' \
              '(i.e. no cartesian grid yet).'
        raise RuntimeError(msg)

    if isinstance(data.grid, radmc3dOctree):
        msg = 'Vertical hydrostatic equilibrium structure iteration has been implemented for regular grids only ' \
               '(i.e. no Octree AMR yet)'
        raise RuntimeError(msg)

    if mu is None:
        # Fix the mean molecular weight to 2.3
        mu = 2.3

    if isinstance(gsize, float) | isinstance(gsize, np.float64):
        if data.rhodust.shape[3] == 1:
            gsize = [gsize]
        else:
            msg = 'The input data contains more than one dust species, but only a single grain size is given. ' \
                  'The number of grain sizes in gsize and data should match.'
            raise ValueError(msg)

    # Pre-calculate some constants
    A = mu * nc.mp * nc.gg * mstar / nc.kk
    cost = np.cos(data.grid.y)
    costi = np.cos(data.grid.yi)


    if mstar is None:
        raise ValueError('Unkonwn mstar. \n The stellar mass is required to calculate the '
                         + ' vertical structure of the disk')

    if idust is None:
        print(' Unknown idust. No dust index was given for which the vertical structure should be calculated, '
              + ' so we do for all dust species')
        idust = range(data.rhodust.shape[3])
    else:
        if isinstance(idust, int) | isinstance(idust, float):
            idust = [int(idust)]

    #
    # Calculate the initial surface density
    #

    # Get the cell volumes
    vol = data.grid.getCellVolume()
    # Calculate the surface of each grid facet in the midplane
    surf = np.zeros([data.grid.nx, data.grid.nz], dtype=np.float64)
    diff_r2 = (data.grid.xi[1:] ** 2 - data.grid.xi[:-1] ** 2) * 0.5
    diff_phi = data.grid.zi[1:] - data.grid.zi[:-1]
    for ix in range(data.grid.nx):
        surf[ix, :] = diff_r2[ix] * diff_phi

    mass = np.zeros([data.grid.nx, data.grid.nz, data.rhodust.shape[3]], dtype=np.float64)
    sigma_init = np.zeros([data.grid.nx, data.grid.nz, data.rhodust.shape[3]], dtype=np.float64)
    for i in range(data.rhodust.shape[3]):
        mass[:, :, i] = (vol * data.rhodust[:, :, :, i]).sum(1)
        sigma_init[:, :, i] = mass[:, :, i] / surf

    # mass = np.array([(vol * data.rhodust[:, :, :, i]).sum(1) for i in idust])
    # sigma_init = mass / surf

    # To improve the smoothness of the temperature structure, if the density structure is
    #  symmetric to the disk midplane we use T_new(theta) = T_new(pi-theta) = 0.5 * (T(theta) + T(pi-theta))
    if vmean_temp:
        if abs(data.grid.yi[data.grid.nyi - 1] - np.pi / 2.) < 1e-8:
            raise RuntimeError("Cannot average temperature in the vertical direction if theta mirroring is active")
        else:
            print(' Smoothing the vertical temperature structure by averaging the temperature of the two half \n'
                  ' planes above and below the disk midplane')
            dusttemp_dummy = np.zeros(data.dusttemp.shape, dtype=np.float64)
            for iy in range(int(data.grid.ny / 2)):
                print(iy)
                dusttemp_dummy[:, iy, :, :] = 0.5 * (data.dusttemp[:, iy, :, :]
                                                     + data.dusttemp[:, data.grid.ny - 1 - iy, :, :])
                dusttemp_dummy[:, data.grid.ny - 1 - iy, :, :] = dusttemp_dummy[:, iy, :, :]

    # Take the temperature of the dust component that represents the gas temperature
    dusttemp_dummy = data.dusttemp[:, :, :, ispec_tgas]

    # rho_new = np.zeros(data.rhodust.shape, dtype=np.float64)
    rho_new = np.array(data.rhodust)

    if gsize is not None:
        if len(gsize) != 0:
            dusttemp = np.zeros([data.grid.nx, data.grid.ny, data.grid.nz], dtype=np.float64)
            w = np.zeros(data.rhodust.shape, dtype=np.float64)
            for ispec in idust:
                w[:, :, :, ispec] = gsize[ispec] ** 2 * (data.rhodust[:, :, :, ispec] / gsize[ispec] ** 3)

            wnorm = w.sum(3)
            for ispec in idust:
                w[:, :, :, ispec] = w[:, :, :, ispec] / wnorm

            for ispec in idust:
                dusttemp = dusttemp + data.dusttemp[:, :, :, ispec] * w[:, :, :, ispec]
        else:
            dusttemp = np.array(dusttemp_dummy)
    else:
        dusttemp = np.array(dusttemp_dummy)

    # Loop over all dust species where we should calculate the vertical structure
    for ispec in idust:
        rho_new[:, :, :, ispec] = 0.
        for ir in range(data.grid.nx):
            print(ir, data.grid.nx - 1)
            r = data.grid.x[ir]
            z = r * cost
            zi = r * costi
            dz = z[:-1] - z[1:]
            const = A / r ** 3

            # Do we have theta mirroring active?
            if abs(data.grid.yi[data.grid.nyi - 1] - np.pi / 2.) < 1e-8:
                for ip in range(data.grid.nz):
                    # dlgrho = np.log(data.rhodust[ir, 1:, ip, ispec]) - np.log(data.rhodust[ir, :-1, ip, ispec])
                    temp = dusttemp[ir, :, ip]

                    it = data.grid.ny - 1
                    temp[it] = 0.5 * (temp[it] + temp[it - 1])

                    dlgtemp = np.log(temp[1:]) - np.log(temp[:-1])
                    zpt = z / temp
                    zpt = 0.5 * (zpt[1:] + zpt[:-1])

                    # Calculate the normalized (rho[z=0] = 1.0) density
                    rho_new[ir, data.grid.ny - 1, ip, ispec] = 1.0

                    for it in range(data.grid.ny - 1, 0, -1):
                        rho_new[ir, it, ip, ispec] = rho_new[ir, it + 1, ip, ispec] * np.exp(
                            -(const * zpt[it] + dlgtemp[it] / dz[it]) * dz[it])

                    rho_new = rho_new.clip(1e-90, 1e90)

                    # # Now re-normalize the surface density to the input value
                    # sigma = (data.rhodust[ir, :, ip, ispec] * (zi[1:] - zi[:-1])).sum()
                    # sigma_new = (rho_new[ir, :, ip, ispec] * (zi[1:] - zi[:-1])).sum()
                    #
                    # rho_new[ir, :, ip, ispec] = rho_new[ir, :, ip, ispec] * sigma / sigma_new

            else:
                for ip in range(data.grid.nz):
                    temp = dusttemp[ir, :, ip]
                    dlgtemp = np.log(temp[1:]) - np.log(temp[:-1])
                    zpt = z / temp
                    zpt = 0.5 * (zpt[1:] + zpt[:-1])

                    # Calculate the normalized (rho[z=0] = 1.0) density
                    rho_new[ir, int(data.grid.ny / 2) - 1, ip, ispec] = 1.0
                    rho_new[ir, int(data.grid.ny / 2), ip, ispec] = 1.0

                    #
                    # From the midplane to the north pole
                    #
                    for it in range(int(data.grid.ny / 2), 0, -1):
                        rho_new[ir, it - 1, ip, ispec] = rho_new[ir, it, ip, ispec] \
                                                         * np.exp(-(const * zpt[it - 1] + dlgtemp[it - 1] / dz[it - 1])
                                                                  * dz[it - 1])
                    #
                    # From the midplane to the north pole
                    #
                    for it in range(int(data.grid.ny / 2), data.grid.ny):
                        rho_new[ir, it, ip, ispec] = rho_new[ir, it - 1, ip, ispec] \
                                                     * np.exp((const * zpt[it - 1] + dlgtemp[it - 1]
                                                               / dz[it - 1]) * dz[it - 1])

                    # # Now re-normalize the surface density to the input value
                    # sigma = (data.rhodust[ir, :, ip, ispec] * (zi[1:] - zi[:-1])).sum()
                    # sigma_new = (rho_new[ir, :, ip, ispec] * (zi[1:] - zi[:-1])).sum()
                    #
                    # rho_new[ir, :, ip, ispec] = rho_new[ir, :, ip, ispec] * sigma / sigma_new
                    # rho_new = rho_new.clip(1e-90, 1e90)


        # Renormalize the density
        mass = (vol * rho_new[:, :, :, ispec]).sum(1)
        sigma = mass / surf
        for it in range(data.grid.ny):
            rho_new[:, it, :, ispec] *= (sigma_init[:, :, ispec] / sigma)
        rho_new[:, it, :, ispec].clip(1e-90, 1e+90)

    return rho_new


def readMol(mol=None, fname=None):
    """ Wrapper around the radmc3dMolecule.read() method

       Parameters
       ----------
       mol             : str
                        molecule name (e.g. 'co') if the file name is in the form of 'molecule_<mol>.inp'

       fname           : str
                        full file name
    """

    m = radmc3dMolecule()
    if m.read(mol=mol, fname=fname) is True:
        return m
    else:
        return


def plotSpectrum(a, ev=False, kev=False, micron=False, jy=False, lsun=False,
                 lnu=False, nulnu=False, fnu=False, nufnu=False, dpc=1.e0,
                 oplot=False, xlg=False, ylg=False, obs=False,
                 mol=None, ilin=None, **kwargs):
    """Plot the spectrum / SED 

    Parameters
    ----------
    a               : ndarray
                     A 2D array of size [Nfreq,2] returned by readSpectrum(). 
                     [:,0] - wavelength in micrometer, or for line data the velocity in km/s
                     [:,1] - flux density in erg/s/cm/cm/Hz
    ev              : bool
                     True --> energy in electronvolt (default=Hz)

    kev             : bool 
                     True --> energy in kiloelectronvolt (default=Hz)

    micron          : bool
                     True --> wavelength in micron (default=Hz)

    jy              : bool
                     True --> Flux in Jansky

    lnu             : bool
                     True --> L_nu (default L_nu)

    fnu             : bool
                     True --> F_nu in units of erg/s/cm^2/Hz(default L_nu)

    nufnu           : bool
                     True --> nu*F_nu in units of erg/s/cm^2 (default L_nu)

    nulnu           : bool
                     True --> nu*L_nu (default F_nu)

    lsun            : bool
                     True --> nu*L_nu in units of solar luminosity

    dpc             : bool
                     Distance of observer in units of parsec (Default: 1 pc)

    oplot           : bool
                     True --> Plot without refreshing subplot

    xlg             : bool
                     True --> logarithmic x-axis

    ylg             : bool
                     True --> logarithmic y-axis

    obs             : bool
                     True --> Treat the spectrum as an observation
                               (i.e. do not scale with dpc^(-2))

    mol             : radmc3dMolecule
                     (optional) Molecule data (see radmc3dMolecule class)
                      This is required if you want to plot a line spectrum
                      with on the x-axis the radial velocity in km/s

    ilin            : bool
                     (if set) the index of the line (of mol; starting,
                      as in RADMC-3D, with the index 1) which shall act
                      as the 0 km/s wavelength reference. If ilin is set
                      the x axis will be in km/s (overriding other settings)

    """
    #
    # Basic
    #
    lam = a[:, 0]
    fluxnu = a[:, 1]
    #
    # Calculate frequency in Hz
    #
    freq = 1e4 * nc.cc / lam
    #
    # Default: frequency in Hz
    #
    xcoord = freq
    xtitle = r'$\nu [\mathrm{Hz}]$'
    #
    # If ev: electronvolt
    #
    if ev:
        xcoord = 4.13568842841e-15 * freq
        xtitle = r'$h\nu [\mathrm{eV}]$'
    #
    # If kev: kiloelectronvolt
    #
    if kev:
        xcoord = 4.13568842841e-18 * freq
        xtitle = r'$h\nu [\mathrm{KeV}]$'
    #
    # If micron
    #
    if micron:
        xcoord = lam
        xtitle = r'$\lambda [\mu\mathrm{m}]$'
    #
    # Plot nuFnu or Fnu (same with Lnu)? And what about Fnu vs Lnu?
    #
    # Default:
    sed = True
    ylum = False
    # The flags:
    if jy:
        sed = False
    if fnu:
        sed = False
        ylum = False
    if lnu:
        sed = False
        ylum = True
    if nulnu:
        sed = True
        ylum = True
    if fnu:
        sed = False
        ylum = False
    if nufnu:
        sed = True
        ylum = False
    if jy:
        ylum = False
    if lsun:
        ylum = True
        sed = True
    #
    # If ilin is set, then override the above and use instead the line
    # as a reference and use km/s as x-axis
    #
    if ilin is not None:
        if mol is None:
            raise ValueError("Unknown mol. If ilin is set, the molecular data should also be provided as mol=...")
        else:
            freq0 = mol.freq[ilin - 1]
            xcoord = nc.cc * (freq0 - freq) / freq0 / 1.e5
            xtitle = '$\Delta v [\mathrm{km/s}]$'
    #
    # Which plot to make? Lum or flux?
    #
    if not ylum:
        #
        # Plot spectrum as flux at a certain distance
        #
        if not obs:
            distfact = 1.0 / (dpc ** 2)
        else:
            distfact = 1.0
        #
        # Set the vertical axis name
        #
        if not jy:
            if not sed:
                lumfact = 1.0
                ytitle = r'$F_{\nu}\; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{Hz}^{-1}\, \mathrm{s}^{-1}]$'
            else:
                lumfact = 1.0 * freq
                ytitle = r'$\nu F_{\nu}\; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$'
        else:
            if not sed:
                lumfact = 1e+23
                ytitle = r'$F_{\nu} [Jy]$'
            else:
                lumfact = 1e+23 * freq
                ytitle = r'$\nu F_{\nu} [JyHz]$'
    else:
        #
        # Plot spectrum as luminosity
        #
        if not obs:
            distfact = 1.1965280793e38  # = 4*pi*(1 parsec)^2 = 1.19d38 cm^2
        else:
            distfact = dpc ** 2 * 1.1965280793e38

        if not sed:
            lumfact = 1.e0
            ytitle = r'$L_{\nu}\; [\mathrm{erg}\,\mathrm{Hz}^{-1}\, \mathrm{s}^{-1}]$'
        else:
            if not lsun:
                lumfact = 1.0 * freq
                ytitle = r'$\nu L_{\nu}\; [\mathrm{erg}\, \mathrm{s}^{-1}]$'
            else:
                lumfact = freq * 2.5956986e-34
                ytitle = r'$\nu L_{\nu}\; [L_{\odot}]$'

    #
    # The data on the y axis
    #
    ycoord = distfact * lumfact * fluxnu
    #
    # If not oplot, then reset the subplot and set the axes
    #
    if not oplot:
        plt.cla()
        if xlg:
            plt.xscale('log')
        if ylg:
            plt.yscale('log')
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
    #
    # Now plot
    #
    plt.plot(xcoord, ycoord, **kwargs)


def gmass(x=None, y=None, z=None, dx=None, dy=None, dz=None, model=None, ppar=None, **kwargs):
    """
    Example function to be used as decision function for resolving cells in tree building. It calculates the gas density
    at a random sample of coordinates within a given cell than take the ratio of the max/min density. If it is larger
    than a certain threshold value it will return True (i.e. the cell should be resolved) if the density variation is 
    less than the threshold it returns False (i.e. the cell should not be resolved)

    Parameters
    ----------

    x       : ndarray
              Cell centre coordinates of the cells in the first dimension

    y       : ndarray
              Cell centre coordinates of the cells in the second dimension

    z       : ndarray
              Cell centre coordinates of the cells in the third dimension

    dx      : ndarray
              Half size of the cells in the first dimension

    dy      : ndarray
              Half size of the cells in the second dimension

    dz      : ndarray
              Half size of the cells in the third dimension

    model   : object
              A radmc3dPy model (must contain a getGasDensity() function) 

    ppar    : dictionary
              All parameters of the problem (from the problem_params.inp file). It is not used here, but must be present 
              for compatibility reasons.

    **kwargs: dictionary
              Parameters used to decide whether the cell should be resolved. It should contain the following keywords; 
              'nsample', which sets the number of random points the gas desity is sampled at within the cell and 
              'threshold' that sets the threshold value for max(gasdens)/min(gasdens) above which the cell should 
              be resolved.
    """

    ncell = x.shape[0]
    rho = np.zeros([ncell, kwargs['nsample']], dtype=np.float64)

    for isample in range(int(kwargs['nsample'])):
        xoffset = (np.random.random_sample(ncell) - 0.5) * dx * 4.0
        yoffset = (np.random.random_sample(ncell) - 0.5) * dy * 4.0
        zoffset = (np.random.random_sample(ncell) - 0.5) * dz * 4.0
        rho[:, isample] = model.getGasDensity(x + xoffset, y + yoffset, z + zoffset, ppar=ppar)

    mass = rho.max(1) * dx * dy * dz * 8.0
    jj = (mass > ppar['threshold'])

    decision = np.zeros(ncell, dtype=bool)
    if True in jj:
        decision[jj] = True

    return decision


def gdensMinMax(x=None, y=None, z=None, dx=None, dy=None, dz=None, model=None, ppar=None, **kwargs):
    """
    Example function to be used as decision function for resolving cells in tree building. It calculates the gas density
    at a random sample of coordinates within a given cell than take the ratio of the max/min density. If it is larger
    than a certain threshold value it will return True (i.e. the cell should be resolved) if the density variation is 
    less than the threshold it returns False (i.e. the cell should not be resolved)

    Parameters
    ----------

    x       : ndarray
              Cell centre coordinates of the cells in the first dimension

    y       : ndarray
              Cell centre coordinates of the cells in the second dimension

    z       : ndarray
              Cell centre coordinates of the cells in the third dimension

    dx      : ndarray
              Half size of the cells in the first dimension

    dy      : ndarray
              Half size of the cells in the second dimension

    dz      : ndarray
              Half size of the cells in the third dimension

    model   : object
              A radmc3dPy model (must contain a getGasDensity() function) 

    ppar    : dictionary
              All parameters of the problem (from the problem_params.inp file). It is not used here, but must be present 
              for compatibility reasons.

    **kwargs: dictionary
              Parameters used to decide whether the cell should be resolved. It should the following keywords; 
              'nsample', which sets the number of random points the gas desity is sampled at within the cell and 
              'threshold' that sets the threshold value for max(gasdens)/min(gasdens) above which the cell should 
              be resolved.
    """

    ncell = x.shape[0]
    rho = np.zeros([ncell, kwargs['nsample']], dtype=np.float64)

    for isample in range(kwargs['nsample']):
        xoffset = (np.random.random_sample(ncell) - 0.5) * dx * 4.0
        yoffset = (np.random.random_sample(ncell) - 0.5) * dy * 4.0
        zoffset = (np.random.random_sample(ncell) - 0.5) * dz * 4.0
        rho[:, isample] = model.getGasDensity(x + xoffset, y + yoffset, z + zoffset, ppar=ppar)

    rho_max = rho.max(axis=1)
    rho_min = rho.min(axis=1)
    # jj      = ((rho_max/rho_min)>ppar['threshold'])
    jj = ((rho_max - rho_min) / rho_max > ppar['threshold'])

    decision = np.zeros(ncell, dtype=bool)
    if True in jj:
        decision[jj] = True
    return decision


def findContainerLeafID(cellCRD=None, cellHW=None, xi=None, yi=None, zi=None, childID=None, isLeaf=None, nChild=None,
                        crd=None):
    """
    Function to find the tree index of a leaf cell containing a given point in space, i.e. if the following is true : 
    xcell - dxcell <= xpoint < xcell + dxcell for each dimension. This function is to be used in multiprocessing.

    Parameters
    ----------

    cellCRD         : ndarray
                      Array with dimensions [ncell, 3] containing the cell centre coordiantes of the tree

    cellHW          : ndarray
                      Array with dimensions [ncell, 3] containing the half width of cells in the tree

    xi              : ndarray
                      Array of cell interface indices in the base grid in the first dimension

    yi              : ndarray
                      Array of cell interface indices in the base grid in the second dimension

    zi              : ndarray
                      Array of cell interface indices in the base grid in the third dimension

    childID         : ndarray
                      Child index array

    isLeaf          : ndarray
                      Boolean array containing the node type for each cell (True - leaf, False - branch) 

    nChild          : int
                      Number of children (8,4,2 for 3,2,1 active dimensions)

    crd             : ndarray
                      Array of length 3 containing the coordinates of the point whose container
                      leaf is to be found

    Returns
    -------
    Tree index of the container leaf if it is found. If the point is outside of the base grid -1 is returned.

    """

    leafID = -1

    if (crd[0] < xi[0]) | (crd[0] > xi[-1]):
        return leafID
    if (crd[1] < yi[0]) | (crd[1] > yi[-1]):
        return leafID
    if (crd[2] < zi[0]) | (crd[2] > zi[-1]):
        return leafID

    ix = np.searchsorted(xi, crd[0])
    iy = np.searchsorted(yi, crd[1])
    iz = np.searchsorted(zi, crd[2])

    if xi[ix] != crd[0]:
        ix -= 1
    if yi[iy] != crd[1]:
        iy -= 1
    if zi[iz] != crd[2]:
        iz -= 1

    nxRoot = xi.shape[0] - 1
    nyRoot = yi.shape[0] - 1
    nzRoot = zi.shape[0] - 1

    if crd[0] == xi[-1]:
        ix = nxRoot - 1
    if crd[1] == yi[-1]:
        iy = nyRoot - 1
    if crd[2] == zi[-1]:
        iz = nzRoot - 1

    ind = iz * nyRoot * nxRoot + iy * nxRoot + ix
    dum = findContainerLeafIDRec(cellCRD[:, 0], cellCRD[:, 1], cellCRD[:, 2], cellHW[:, 0], cellHW[:, 1], cellHW[:, 2],
                                 childID, isLeaf, nChild, crd, ind)
    if dum is None:
        leafID = -1
    else:
        leafID = dum

    return leafID


def findContainerLeafIDRec(x=None, y=None, z=None, dx=None, dy=None, dz=None, childID=None, isLeaf=None, nChild=None,
                           crd=(), cellID=None):
    """
    Recursive function to find the leaf cell in the tree that contains a given point in space

    Parameters
    ----------

    x                 : ndarray
                        Tree cell center array in the first dimension

    y                 : ndarray
                        Tree cell center array in the second dimension

    z                 : ndarray
                        Tree cell center array in the tird dimension

    dx                : ndarray
                        Tree cell halfwidth array in the first dimension

    dy                : ndarray
                        Tree cell halfwidth array in the second dimension

    dz                : ndarray
                        Tree cell halfwidth array in the third dimension

    childID           : list
                        List of children indices. Each list element is an ndarray with nChild elements containing 
                        the child indices

    isLeaf            : ndarray
                        Boolean array for the cell type (True - leaf, False - branch)

    nChild            : int
                        Nr of children (i.e. 8, 4, or 2 for 3, 2, 1 active dimensions, respectively)

    crd               : ndarray
                        Three element list/tuple/array containing the point coordinates

    cellID            : int
                        Index of cell to be tested
    """

    xmin = x[cellID] - dx[cellID]
    xmax = x[cellID] + dx[cellID]
    ymin = y[cellID] - dy[cellID]
    ymax = y[cellID] + dy[cellID]
    zmin = z[cellID] - dz[cellID]
    zmax = z[cellID] + dz[cellID]

    if isLeaf[cellID]:
        if (((crd[0] >= xmin) & (crd[0] < xmax)) &
                ((crd[1] >= ymin) & (crd[1] < ymax)) &
                ((crd[2] >= zmin) & (crd[2] < zmax))):
            return cellID
        else:
            return None

    else:
        dum = None
        for i in range(nChild):
            dum = findContainerLeafIDRec(x, y, z, dx, dy, dz, childID, isLeaf, nChild, crd, childID[cellID][i])
            if dum is not None:
                break

        return dum


def interpolateOctree(data=None, x=None, y=None, z=None, var=None, nproc=1):
    """
    Nearest neighbour inteprolation on an octree

    data        : radmc3dData
                  Data container

    x           : ndarray
                  Coordiantes of the point to be interpolated on in the first dimension

    y           : ndarray
                  Coordiantes of the point to be interpolated on in the second dimension

    z           : ndarray
                  Coordiantes of the point to be interpolated on in the third dimension

    var         : list
                  Name of the variables to be interpolated, supported names are:
                  ddens, dtemp, gdens, ndens, gtemp, gvel, vturb

    nproc       : int
                  Number of processes to be used (for parallel computing) 

    Returns:
    --------
    A dictionary with the interpolated fields 

    """

    if var is None:
        var = ['ddens']

    if nproc == 1:
        print("Nearest neighbour interpolation using " + ("%d" % nproc) + ' process')
    else:
        print("Nearest neighbour interpolation using " + ("%d" % nproc) + ' processes')
    nx = x.shape[0]
    ny = y.shape[0]
    nz = z.shape[0]

    npoint = nx * ny * nz

    if nproc > 1:
        cellCRD = np.zeros([data.grid.nCell, 3], dtype=np.float64)
        cellHW = np.zeros([data.grid.nCell, 3], dtype=np.float64)

        cellCRD[:, 0] = data.grid.x
        cellCRD[:, 1] = data.grid.y
        cellCRD[:, 2] = data.grid.z
        cellHW[:, 0] = data.grid.dx
        cellHW[:, 1] = data.grid.dy
        cellHW[:, 2] = data.grid.dz
        childID = np.array(data.grid.childID)

        chunkSize = int(np.ceil(float(npoint) / nproc))

        crdList = np.zeros([npoint, 3], dtype=np.float64)
        for ind in range(npoint):
            ix = int(np.floor(ind / ny / nz))
            iy = int(np.floor((ind - ix * ny * nz) / nz))
            iz = int(ind - ix * ny * nz - iy * nz)

            crdList[ind, 0] = x[ix]
            crdList[ind, 1] = y[iy]
            crdList[ind, 2] = z[iz]

        pool = Pool(processes=nproc)
        target = partial(findContainerLeafID, cellCRD, cellHW, data.grid.xi, data.grid.yi, data.grid.zi, childID,
                         data.grid.isLeaf, data.grid.nChild)
        res = pool.map(target, crdList, chunksize=chunkSize)
        res = np.array(res)
        pool.close()

        idata = {'cellID': None, 'rhodust': None, 'dusttemp': None, 'rhogas': None, 'ndens_mol': None, 'gtemp': None,
                 'vturb': None, 'gvel': None}

        idata['cellID'] = res

        if 'ddens' in var:
            ndust = data.rhodust.shape[1]
            idata['rhodust'] = np.zeros([nx, ny, nz, ndust], dtype=np.float64)
            for ind in range(npoint):
                ix = int(np.floor(ind / ny / nz))
                iy = int(np.floor((ind - ix * ny * nz) / nz))
                iz = int(ind - ix * ny * nz - iy * nz)

                if res[ind] >= 0:
                    idata['rhodust'][ix, iy, iz, :] = data.rhodust[data.grid.leafID[res[ind]], :]
                else:
                    idata['rhodust'][ix, iy, iz, :] = 0

        if 'dtemp' in var:
            ndust = data.dusttemp.shape[1]
            idata['dusttemp'] = np.zeros([nx, ny, nz, ndust], dtype=np.float64)
            for ind in range(npoint):
                ix = int(np.floor(ind / ny / nz))
                iy = int(np.floor((ind - ix * ny * nz) / nz))
                iz = int(ind - ix * ny * nz - iy * nz)

                if res[ind] is not None:
                    idata['dusttemp'][ix, iy, iz, :] = data.dusttemp[data.grid.leafID[res[ind]], :]
                else:
                    idata['dusttemp'][ix, iy, iz, :] = 0

        if 'gdens' in var:
            idata['rhogas'] = np.zeros([nx, ny, nz], dtype=np.float64)
            for ind in range(npoint):
                ix = int(np.floor(ind / ny / nz))
                iy = int(np.floor((ind - ix * ny * nz) / nz))
                iz = int(ind - ix * ny * nz - iy * nz)

                if res[ind] is not None:
                    idata['rhogas'][ix, iy, iz] = data.rhogas[data.grid.leafID[res[ind]]]
                else:
                    idata['rhogas'][ix, iy, iz] = 0

        if 'ndens' in var:
            idata['ndens_mol'] = np.zeros([nx, ny, nz], dtype=np.float64)
            for ind in range(npoint):
                ix = int(np.floor(ind / ny / nz))
                iy = int(np.floor((ind - ix * ny * nz) / nz))
                iz = int(ind - ix * ny * nz - iy * nz)

                if res[ind] is not None:
                    idata['ndens_mol'][ix, iy, iz] = data.ndens_mol[data.grid.leafID[res[ind]]]
                else:
                    idata['ndens_mol'][ix, iy, iz] = 0

        if 'gtemp' in var:
            idata['gastemp'] = np.zeros([nx, ny, nz], dtype=np.float64)
            for ind in range(npoint):
                ix = int(np.floor(ind / ny / nz))
                iy = int(np.floor((ind - ix * ny * nz) / nz))
                iz = int(ind - ix * ny * nz - iy * nz)

                if res[ind] is not None:
                    idata['gastemp'][ix, iy, iz] = data.gastemp[data.grid.leafID[res[ind]]]
                else:
                    idata['gastemp'][ix, iy, iz] = 0

        if 'vturb' in var:
            idata['vturb'] = np.zeros([nx, ny, nz], dtype=np.float64)
            for ind in range(npoint):
                ix = int(np.floor(ind / ny / nz))
                iy = int(np.floor((ind - ix * ny * nz) / nz))
                iz = int(ind - ix * ny * nz - iy * nz)

                if res[ind] is not None:
                    idata['vturb'][ix, iy, iz] = data.vturb[data.grid.leafID[res[ind]]]
                else:
                    idata['vturb'][ix, iy, iz] = 0

        if 'gvel' in var:
            idata['gvel'] = np.zeros([nx, ny, nz, 3], dtype=np.float64)
            for ind in range(npoint):
                ix = int(np.floor(ind / ny / nz))
                iy = int(np.floor((ind - ix * ny * nz) / nz))
                iz = int(ind - ix * ny * nz - iy * nz)

                if res[ind] is not None:
                    idata['gasvel'][ix, iy, iz, :] = data.gasvel[data.grid.leafID[res[ind]], :]
                else:
                    idata['gasvel'][ix, iy, iz, :] = 0
    else:

        idata = {}
        if 'ddens' in var:
            ndust = data.rhodust.shape[1]
            idata['rhodust'] = np.zeros([nx, ny, nz, ndust], dtype=np.float64)
        if 'dtemp' in var:
            ndust = data.dusttemp.shape[1]
            idata['dusttemp'] = np.zeros([nx, ny, nz, ndust], dtype=np.float64)
        if 'gdens' in var:
            idata['rhogas'] = np.zeros([nx, ny, nz], dtype=np.float64)
        if 'ndens' in var:
            idata['ndens_mol'] = np.zeros([nx, ny, nz], dtype=np.float64)
        if 'gtemp' in var:
            idata['gastemp'] = np.zeros([nx, ny, nz], dtype=np.float64)
        if 'vturb' in var:
            idata['vturb'] = np.zeros([nx, ny, nz], dtype=np.float64)
        if 'gvel' in var:
            idata['gvel'] = np.zeros([nx, ny, nz, 3], dtype=np.float64)

        idata['cellID'] = np.zeros(npoint, dtype=np.int_)

        for ind in range(npoint):

            ix = int(np.floor(ind / ny / nz))
            iy = int(np.floor((ind - ix * ny * nz) / nz))
            iz = int(ind - ix * ny * nz - iy * nz)

            cellID = data.grid.getContainerLeafID((x[ix], y[iy], z[iz]))

            if cellID is not None:
                idata['cellID'][ind] = cellID
                if 'ddens' in var:
                    if cellID >= 0:
                        idata['rhodust'][ix, iy, iz, :] = data.rhodust[data.grid.leafID[cellID]]
                    else:
                        idata['rhodust'][ix, iy, iz, :] = 0.0
                if 'dtemp' in var:
                    if cellID >= 0:
                        idata['dusttemp'][ix, iy, iz, :] = data.dusttemp[data.grid.leafID[cellID]]
                    else:
                        idata['dusttemp'][ix, iy, iz, :] = 0.0
                if 'gdens' in var:
                    if cellID >= 0:
                        idata['rhogas'][ix, iy, iz] = data.ndens_mol[data.grid.leafID[cellID]]
                    else:
                        idata['rhogas'][ix, iy, iz] = 0.0
                if 'ndens' in var:
                    if cellID >= 0:
                        idata['ndens_mol'][ix, iy, iz] = data.ndens_mol[data.grid.leafID[cellID]]
                    else:
                        idata['ndens_mol'][ix, iy, iz] = 0.0
                if 'gtemp' in var:
                    if cellID >= 0:
                        idata['gastemp'][ix, iy, iz] = data.gastemp[data.grid.leafID[cellID]]
                    else:
                        idata['gastemp'][ix, iy, iz] = 0.0
                if 'vturb' in var:
                    if cellID >= 0:
                        idata['vturb'][ix, iy, iz] = data.vturb[data.grid.leafID[cellID]]
                    else:
                        idata['vturb'][ix, iy, iz] = 0.0
                if 'gvel' in var:
                    if cellID >= 0:
                        idata['gasvel'][ix, iy, iz, :] = data.gasvel[data.grid.leafID[cellID], :]
                    else:
                        idata['gasvel'][ix, iy, iz, :] = 0.0

    return idata


def plotSlice2D(data=None, var='ddens', plane='xy', crd3=0.0, icrd3=None, ispec=-1, xlim=(), ylim=(), log=False,
                linunit='cm', angunit='rad',
                nx=100, ny=100, showgrid=False, gridcolor='k', gridalpha=1.0, nproc=1,
                contours=False, clev=None, clmin=None, clmax=None, ncl=None, cllog=False, clcol='k',
                cllabel=False, cllabel_fontsize=10, cllabel_fmt="%.1f", clalpha=1.0, ax=None, lattitude=True,
                **kwargs):
    """
    Function to plot an axis-aligned 2D slice of the variables in the model. Any additional keyword
    argument above the listed ones will be passed on to matplotlib.pylab.pcolormesh(). For an octree grid the variables 
    are interpolated onto a regular grid using nearest neighbour interpolation before plotting. 
    The size and resolution of the regular image grid can be set at input. 


    Parameters
    ----------

    data            : radmc3dData
                      Instance of radmc3dData containing the field variable to be displayed

    var             : {'ddens', 'dtemp', 'gdens', 'ndens', 'gtemp', 'vturb', 'vx', 'vy', 'vz', 'taux', 'tauy'}
                      Variable to be displayed

    plane           : {'xy', 'xz', 'yz', 'yx, 'zx', 'yz'}
                      Plane to be displayed           

    crd3            : float
                      Coordinate of the third dimension (i.e. when plotting a slice in the x-y plane, crd3 is the 
                      z-coordinate)

    icrd3           : int
                      Index of the third coordinate in the grid (only for regular grid!)

    ispec           : int
                      Dust species index. If negative dust densities will be summed up and the total cumulative density 
                      will be displayed

    xlim            : tuple
                      Coordinate boundaries in the first dimension of the plot (also the coordinate boundary of the 
                      regular grid data on
                      AMR grids are interpolated to)

    ylim            : tuple
                      Coordinate boundaries in the second dimension of the plot (also the coordinate boundary of the 
                      regular grid data on
                      AMR grids are interpolated to)

    log             : bool
                      If True the contour/image will be displayed on a logarithmic stretch

    linunit         : {'cm', 'au', 'pc', 'rs'}
                      Unit selection for linear image coordinate axes.

    nx              : int
                      Number of horizontal pixels in the interpolated image if the data is defined in an Octree

    ny              : int
                      Number of vertical pixels in the interpolated image if the data is defined in an Octree

    showgrid        : bool
                      If True the spatial grid will be overlayed 

    gridcolor       : str
                      Color of the spatial grid overlay

    gridalpha       : float
                      Opacity of the lines in the spatial grid overlay (0.0 - fully transparent, 1.0 - fully opaque)

    angunit         : {'rad', 'deg'}
                      Unit selection for angular image coordinate axes (only if spherical coordinate system is used).

    nproc           : int
                      Number of parallel processes to be used for interpolation. 

    contours        : bool
                      If True contour lines are plotted, if False a colorscale plot will be created

    clev            : ndarray  
                      A numpy ndarray containing the levels to be displayed with contour lines. If clev is set
                      then clmin, clmax and ncl are omitted

    clmin           : float
                      Min. contour level (for setting auto-contours between clmin and clmax at ncl values)

    clmax           : float
                      Max. contour level (for setting auto-contours between clmin and clmax at ncl values)

    ncl             : float
                      Number of contour levels (for setting auto-contours between clmin and clmax at ncl values)

    cllog           : bool
                      If clmin, clmax and ncl are used to generate the contour levels, then if cllog is True
                      the contours will be log-scaled

    clcol           : str
                      Color-code for the contour lines for single color contours

    cllabel         : bool
                      If True the contour line values will be displayed, if False only the contour lines will be
                      displayed (default = False)

    cllabel_fontsize: int
                      Size of the font used to displaye the contour line values

    cllabel_fmt     : str
                      Format of the contour line labels (default "%.1f")

    clalpha         : float
                      Transparency of the contour lines (1.0 fully opaque, 0.0 fully transparent)

    lattitude       : bool
                      If the coordinate sytem used in RADMC-3D is spherical, then the 2nd coordiante is the 
                      co-lattitude. If lattitude is set to True then the 2nd coordinate in the RADMC-3D grid will be 
                      transformet to true lattitude (i.e. pi/2.-colattitude). If set to false the original co-lattitude
                      will be used. 

    ax              : matplotlib.axes.Axes
                      Matplotlib axis to plot to


    Keyword Arguments :
                      All other keyword arugments will be passed to pcolormesh() or countour()

    """

    octree = False
    #
    # Check the input consistency
    #
    if data is None:
        raise ValueError('Unkonwn data. Data to be plotted is not specified.')

    if data.grid is None:
        raise AttributeError('Missing grid information in data. Plots cannot be made without a spatial grid.')
    else:
        if isinstance(data.grid, radmc3dOctree):
            octree = True

    if not octree:
        if icrd3 is None:
            if crd3 is None:
                raise ValueError('Unknown coordinate for the third dimension (icrd3/crd3)')

    var = var.strip().lower()
    varFound = False
    if var == 'ddens':
        varFound = True
    if var == 'dtemp':
        varFound = True
    if var == 'gdens':
        varFound = True
    if var == 'ndens':
        varFound = True
    if var == 'gtemp':
        varFound = True
    if var == 'vx':
        varFound = True
    if var == 'vy':
        varFound = True
    if var == 'vz':
        varFound = True
    if var == 'vturb':
        varFound = True
    if var == 'taux':
        varFound = True
        if octree:
            raise RuntimeError('Optical depth calculation has not yet been implemented for octrees')
    if var == 'tauy':
        varFound = True
        if octree:
            raise RuntimeError('Optical depth calculation has not yet been implemented for octrees')

    if not varFound:
        raise ValueError('Unknown variable to be plotted : ', var, '\n Allowed variable names are : ddens, dtemp, '
                         + 'gdens, ndens, dtemp, vx, vy, vz, vturb, taux, '
                         + 'tauy')

    #
    # Get the units
    #
    if linunit.strip().lower() == 'cm':
        linunit_label = '[cm]'
        linunit_norm = 1.
    elif linunit.strip().lower() == 'au':
        linunit_label = '[au]'
        linunit_norm = 1. / nc.au
    elif linunit.strip().lower() == 'pc':
        linunit_label = '[pc]'
        linunit_norm = 1. / nc.pc
    elif linunit.strip().lower() == 'rs':
        linunit_label = r'[R$_\odot$]'
        linunit_norm = 1. / nc.rs
    else:
        raise ValueError('Unknown linunit ', linunit, '\nSupported units are : cm, au, pc, rs')

    if angunit.strip().lower() == 'rad':
        angunit_label = '[rad]'
        angunit_norm = 1.0
    elif angunit.strip().lower() == 'deg':
        angunit_label = '[deg]'
        angunit_norm = np.pi / 180.
    else:
        raise ValueError('Unknown angunit ', angunit, '\nSupported units are : deg, rad')

    #
    # Now check which plane to be plotted
    #
    swapDim = False
    xnorm = 1.0
    ynorm = 1.0
    znorm = 1.0
    iplane = -1

    if octree:
        plot_x = xlim[0] + (xlim[1] - xlim[0]) * np.arange(nx, dtype=float) / float(nx - 1)
        plot_y = ylim[0] + (ylim[1] - ylim[0]) * np.arange(ny, dtype=float) / float(ny - 1)
        plot_z = np.array([crd3])
    else:
        plot_x = None
        plot_y = None
        plot_z = None

    if 'x' in plane:

        # xy plane
        if data.grid.x.shape[0] <= 1:
            msg = 'The x dimension is switched off or has only a single grid cell, thus a 2D slice plot cannot ' \
                  'be done in the '+plane+' plane.'
            raise ValueError(msg)
        if 'y' in plane:
            if data.grid.y.shape[0] <= 1:
                msg = 'The y dimension is switched off or has only a single grid cell, thus a 2D slice plot cannot ' \
                      'be done in the '+plane+' plane.'
                raise ValueError(msg)

            iplane = 2
            if not octree:
                plot_x = np.copy(data.grid.x)
                plot_y = np.copy(data.grid.y)
                if lattitude:
                    plot_y = np.pi / 2.0 - data.grid.y

                if icrd3 is None:
                    icrd3 = np.abs(data.grid.z - crd3).argmin()
            if data.grid.crd_sys == 'car':
                xlabel = r'x ' + linunit_label
                ylabel = r'y ' + linunit_label
                xnorm = linunit_norm
                ynorm = linunit_norm
            else:
                xlabel = 'r ' + linunit_label

                if lattitude:
                    ylabel = '$\\pi/2-\\theta$ ' + angunit_label
                else:
                    ylabel = '$\\theta$ ' + angunit_label

                xnorm = linunit_norm
                ynorm = angunit_norm

            if plane == 'yx':
                swapDim = True
                if data.grid.crd_sys == 'car':
                    xlabel = r'y ' + linunit_label
                    ylabel = r'x ' + linunit_label
                    xnorm = linunit_norm
                    ynorm = linunit_norm
                else:
                    xlabel = '$\\theta$ ' + angunit_label
                    ylabel = 'r ' + linunit_label
                    xnorm = angunit_norm
                    ynorm = linunit_norm

            if octree:
                idata = interpolateOctree(data, x=plot_x / xnorm, y=plot_y / ynorm, z=plot_z, var=var, nproc=nproc)
        # xz plane
        elif 'z' in plane:
            if data.grid.z.shape[0] <= 1:
                msg = 'The z dimension is switched off or has only a single grid cell, thus a 2D slice plot cannot ' \
                      'be done in the '+plane+' plane.'
                raise ValueError(msg)

            iplane = 1
            if not octree:
                plot_x = np.copy(data.grid.x)
                plot_y = np.copy(data.grid.z)
                if icrd3 is None:
                    icrd3 = np.abs(data.grid.y - crd3).argmin()
                    if lattitude:
                        icrd3 = np.abs((np.pi / 2. - data.grid.y) - crd3).argmin()

            if data.grid.crd_sys == 'car':
                xlabel = r'x ' + linunit_label
                ylabel = r'z ' + linunit_label
                xnorm = linunit_norm
                ynorm = linunit_norm
            else:
                xlabel = 'r ' + linunit_label
                ylabel = '$\phi$ ' + angunit_label
                xnorm = linunit_norm
                ynorm = angunit_norm
            if plane == 'zx':
                swapDim = True
                if data.grid.crd_sys == 'car':
                    xlabel = r'z ' + linunit_label
                    ylabel = r'x ' + linunit_label
                    xnorm = linunit_norm
                    ynorm = linunit_norm
                else:
                    xlabel = '$\phi$ ' + angunit_label
                    ylabel = 'r ' + linunit_label
                    xnorm = angunit_norm
                    ynorm = linunit_norm

            if octree:
                idata = interpolateOctree(data, x=plot_x / xnorm, y=plot_z, z=plot_y / ynorm, var=var, nproc=nproc)
    # yz plane
    else:
        if data.grid.y.shape[0] <= 1:
            msg = 'The y dimension is switched off or has only a single grid cell, thus a 2D slice plot cannot ' \
                  'be done in the '+plane+' plane.'
            raise ValueError(msg)

        if data.grid.z.shape[0] <= 1:
            msg = 'The z dimension is switched off or has only a single grid cell, thus a 2D slice plot cannot ' \
                  'be done in the '+plane+' plane.'
            raise ValueError(msg)
        iplane = 0
        if not octree:
            plot_x = np.copy(data.grid.y)
            plot_y = np.copy(data.grid.z)
            if lattitude:
                plot_x = (np.pi / 2.0 - data.grid.y)
            if icrd3 is None:
                icrd3 = np.abs(data.grid.x - crd3).argmin()

        if data.grid.crd_sys == 'car':
            xlabel = r'y ' + linunit_label
            ylabel = r'z ' + linunit_label
            xnorm = linunit_norm
            ynorm = linunit_norm

        else:
            if lattitude:
                xlabel = '$\\pi/2-\\theta$ ' + angunit_label
            else:
                xlabel = '$\\theta$ ' + angunit_label
            ylabel = '$\phi$ ' + angunit_label
            xnorm = angunit_norm
            ynorm = angunit_norm
        if plane == 'zy':
            swapDim = True
            if data.grid.crd_sys == 'car':
                xlabel = r'z ' + linunit_label
                ylabel = r'y ' + linunit_label
                xnorm = linunit_norm
                ynorm = linunit_norm
            else:
                xlabel = '$\phi$ ' + angunit_label
                if lattitude:
                    ylabel = '$\\pi/2-\\theta$ ' + angunit_label
                else:
                    ylabel = '$\\theta$ ' + angunit_label
                xnorm = angunit_norm
                ynorm = angunit_norm

        if octree:
            idata = interpolateOctree(data, x=plot_z, y=plot_x / xnorm, z=plot_y / ynorm, var=var, nproc=nproc)

    #
    # Get the variable to be plotted
    #
    if var == 'ddens':
        if isinstance(data.rhodust, int):
            raise ValueError('Dust density is not present in the passed radmc3dData instance')
        else:
            if octree:
                if ispec >= 0:
                    pdata = np.squeeze(idata['rhodust'][:, :, :, ispec])
                else:
                    pdata = np.squeeze(idata['rhodust'].sum(3))

            else:
                if iplane == 0:
                    if ispec >= 0:
                        pdata = data.rhodust[icrd3, :, :, ispec]
                    else:
                        pdata = data.rhodust[icrd3, :, :, :].sum(2)
                elif iplane == 1:
                    if ispec >= 0:
                        pdata = data.rhodust[:, icrd3, :, ispec]
                    else:
                        pdata = data.rhodust[:, icrd3, :, :].sum(2)
                elif iplane == 2:
                    if ispec >= 0:
                        pdata = data.rhodust[:, :, icrd3, ispec]
                    else:
                        pdata = data.rhodust[:, :, icrd3, :].sum(2)

            cblabel = r'$\rho_{\rm dust}$ [g/cm$^3$]'

    elif var == 'dtemp':
        if isinstance(data.dusttemp, int):
            raise ValueError('Dust temperature is not present in the passed radmc3dData instance')
        else:
            if octree:
                if ispec >= 0:
                    pdata = np.squeeze(idata['dusttemp'][:, :, :, ispec])
                else:
                    raise IndexError('Negative dust species index.')

            else:
                if iplane == 0:
                    if ispec >= 0:
                        pdata = data.dusttemp[icrd3, :, :, ispec]
                    else:
                        raise IndexError('Negative dust species index.')
                elif iplane == 1:
                    if ispec >= 0:
                        pdata = data.dusttemp[:, icrd3, :, ispec]
                    else:
                        raise IndexError('Negative dust species index.')
                elif iplane == 2:
                    if ispec >= 0:
                        pdata = data.dusttemp[:, :, icrd3, ispec]
                    else:
                        raise IndexError('Negative dust species index : ', ispec)
            cblabel = r'$T_{\rm dust}$ [K]'

    elif var == 'gdens':
        if isinstance(data.rhogas, int):
            raise ValueError('Gas density is not present in the passed radmc3dData instance')
        else:
            if octree:
                pdata = np.squeeze(idata['rhogas'])
            else:
                if iplane == 0:
                    pdata = data.rhogas[icrd3, :, :]
                elif iplane == 1:
                    pdata = data.rhogas[:, icrd3, :]
                elif iplane == 2:
                    pdata = data.rhogas[:, :, icrd3]
            cblabel = r'$\rho_{\rm gas}$ [g/cm$^3$]'

    elif var == 'ndens':
        if isinstance(data.ndens_mol, int):
            raise ValueError('Gas number density is not present in the passed radmc3dData instance')
        else:
            if octree:
                pdata = np.squeeze(idata['ndens_mol'])
            else:
                if iplane == 0:
                    pdata = data.ndens_mol[icrd3, :, :]
                elif iplane == 1:
                    pdata = data.ndens_mol[:, icrd3, :]
                elif iplane == 2:
                    pdata = data.ndens_mol[:, :, icrd3]
            cblabel = r'$n_{\rm gas}$ [molecule/cm$^3$]'

    elif var == 'gtemp':
        if isinstance(data.gastemp, int):
            raise ValueError('Gas temperature is not present in the passed radmc3dData instance')
        else:
            if octree:
                pdata = np.squeeze(idata['gastemp'])
            else:
                if iplane == 0:
                    pdata = data.gastemp[icrd3, :, :]
                elif iplane == 1:
                    pdata = data.gastemp[:, icrd3, :]
                elif iplane == 2:
                    pdata = data.gastemp[:, :, icrd3]
            cblabel = r'$T_{\rm gas}$ [K]'

    elif var == 'vx':
        if isinstance(data.gasvel, int):
            raise ValueError('Gas velocity is not present in the passed radmc3dData instance')
        else:
            if octree:
                pdata = np.squeeze(idata['gasvel'][:, :, :, 0])
            else:
                if iplane == 0:
                    pdata = data.gasvel[icrd3, :, :, 0]
                elif iplane == 1:
                    pdata = data.gasvel[:, icrd3, :, 0]
                elif iplane == 2:
                    pdata = data.gasvel[:, :, icrd3, 0]
            cblabel = r'$v_{\rm x}$ [cm/s]'

    elif var == 'vy':
        if isinstance(data.gasvel, int):
            raise ValueError('Gas velocity is not present in the passed radmc3dData instance')
        else:
            if octree:
                pdata = np.squeeze(idata['gasvel'][:, :, :, 1])
            else:
                if iplane == 0:
                    pdata = data.gasvel[icrd3, :, :, 1]
                elif iplane == 1:
                    pdata = data.gasvel[:, icrd3, :, 1]
                elif iplane == 2:
                    pdata = data.gasvel[:, :, icrd3, 1]
            cblabel = r'$v_{\rm y}$ [cm/s]'
    elif var == 'vz':
        if isinstance(data.gasvel, int):
            raise ValueError('Gas velocity is not present in the passed radmc3dData instance')
        else:
            if octree:
                pdata = np.squeeze(idata['gasvel'][:, :, :, 2])
            else:
                if iplane == 0:
                    pdata = data.gasvel[icrd3, :, :, 2]
                elif iplane == 1:
                    pdata = data.gasvel[:, icrd3, :, 2]
                elif iplane == 2:
                    pdata = data.gasvel[:, :, icrd3, 2]
            cblabel = r'$v_{\rm z}$ [cm/s]'

    elif var == 'vturb':
        if isinstance(data.vturb, int):
            raise ValueError('Microturbulent velocity is not present in the passed radmc3dData instance')
        else:
            if octree:
                pdata = np.squeeze(idata['vturb'])
            else:
                if iplane == 0:
                    pdata = data.vturb[icrd3, :, :]
                elif iplane == 1:
                    pdata = data.vturb[:, icrd3, :]
                elif iplane == 2:
                    pdata = data.vturb[:, :, icrd3]
            cblabel = r'$v_{\rm turb}$ [cm/s]'

    if var == 'taux':
        if isinstance(data.taux, int):
            raise ValueError('Optical depth is not present in the passed radmc3dData instance')
        else:
            if octree:
                raise RuntimeError('Optical depth calculation has not yet been implemented for octrees.')
            else:
                if data.taux.shape[0] == 0:
                    raise ValueError('Optical depth has not been calculated yet. Run radmc3dData.getTau(wav=X) to '
                                     'calculate the optical depth before calling plotSlice2D')
                if iplane == 0:
                    pdata = data.taux[icrd3, :, :]
                elif iplane == 1:
                    pdata = data.taux[:, icrd3, :]
                elif iplane == 2:
                    pdata = data.taux[:, :, icrd3]
            cblabel = r'$\tau_{\rm r}$'

    if var == 'tauy':
        if isinstance(data.tauy, int):
            raise ValueError('Optical depth is not present in the passed radmc3dData instance')
        else:
            if octree:
                raise RuntimeError('Optical depth calculation has not yet been implemented for octrees.')
            else:
                if iplane == 0:
                    pdata = data.tauy[icrd3, :, :]
                elif iplane == 1:
                    pdata = data.tauy[:, icrd3, :]
                elif iplane == 2:
                    pdata = data.tauy[:, :, icrd3]
            cblabel = r'$\tau_{\rm \theta}$'

    #
    # Apparently there is some inconsistency in the dimensionality of the gas arrays. I.e. when data is read from file
    #  there is always four dimension with the last dimension, for dust used for the dust species index, set to one.
    #  I guess it was meant to be if we have multiple gas species. However, if data is created by the model functions
    #  directly the gas density is usually returned as a three dimensional array. So the quick and dirty fix is to
    #  check for the dimensionality here and decrease it by one stripping the last dimension either by a specified
    #  'ispec' index
    #

    if len(pdata.shape) == 3:
        if pdata.shape[2] == 1:
            pdata = pdata[:, :, 0]
        else:
            raise ValueError('The plotted data has four dimension (3 spatial + 1 species) but the species index is '
                             + ' not set. Specify ispec keyword which of the dimensinos should be plotted')

    #
    # If the dimensions are flipped in the plotted plane (i.e. yx, zy, or xz) flip the array dimensions
    #
    if swapDim:
        pdata = pdata.swapaxes(0, 1)
        dum = np.array(plot_x)
        plot_y = np.array(plot_x)
        plot_x = dum

    if not octree:
        plot_x *= xnorm
        plot_y *= ynorm
    #
    # Get the data min/max values
    #

    if ax is None:
        ax = plt.gca()

    #
    # Now do the colorscale plotting but only if contours is set to False
    #
    if not contours:
        if log:
            if 'vmin' not in kwargs:
                vmin = pdata.min()
            else:
                vmin = kwargs['vmin']

            if 'vmax' not in kwargs:
                vmax = pdata.max()
            else:
                vmax = kwargs['vmax']

            if pdata.min() <= 0:
                pdata = pdata.clip(1e-90)
            p = ax.pcolormesh(plot_x, plot_y, pdata.T, norm=LogNorm(vmin, vmax), **kwargs)
            # p = ax.imshow(pdata.T, origin='lower', extent=(plot_x[0]-dx*0.5, plot_x[-1]+dx*0.5,
            # plot_y[0]-dy*0.5, plot_y[-1]+dy*0.5),
            # norm=LogNorm(vmin, vmax), interpolation='nearest', aspect='auto', **kwargs)
            # Enable rasterization to enable easy save to file
            p.set_rasterized(True)
        else:
            p = ax.pcolormesh(plot_x, plot_y, pdata.T, **kwargs)
            # p = ax.imshow(pdata.T, origin='lower', extent=(plot_x[0]-dx*0.5, plot_x[-1]+dx*0.5,
            # plot_y[0]-dy*0.5, plot_y[-1]+dy*0.5),
            # vmin=vmin, vmax=vmax, interpolation='nearest', aspect='auto', **kwargs)
            # Enable rasterization to enable easy save to file
            p.set_rasterized(True)

        #
        # Generate the colorbar
        #
        cb = plt.colorbar(p)
        cb.set_label(cblabel)

    else:
        if clmin is None:
            clmin = pdata.min()
        if clmax is None:
            clmax = pdata.max()
        if ncl is None:
            ncl = 20

        # Generate the contour levels
        if clev is None:
            if cllog is True:
                clev = clmin * (clmax / clmin) ** (np.arange(ncl, dtype=float) / float(ncl - 1))
            else:
                clev = clmin + (clmax - clmin) * (np.arange(ncl, dtype=float) / float(ncl - 1))
        if (clcol == 'none') | (clcol is None):
            if 'cmap' in kwargs:
                c = ax.contour(plot_x, plot_y, pdata, clev, kwargs['cmap'], alpha=clalpha)
            else:
                c = ax.contour(plot_x, plot_y, pdata, clev, alpha=clalpha)
        else:
            c = ax.contour(plot_x, plot_y, pdata.T, clev, colors=clcol, alpha=clalpha)

        if cllabel:
            plt.clabel(c, inline=1, fontsize=cllabel_fontsize, fmt=cllabel_fmt)

    #
    # Show the grid (currently only for otctree AMR)
    #
    if showgrid:
        if octree:
            ind = 0
            plottedInd = np.zeros(idata['cellID'].shape[0], dtype=np.int_) - 1

            if plane.strip().lower() == 'xy':
                for i in idata['cellID']:
                    if i not in plottedInd:
                        ind += 1
                        bottomleft = ((data.grid.x[i] - data.grid.dx[i]) * xnorm,
                                      (data.grid.y[i] - data.grid.dy[i]) * ynorm)
                        ax.add_patch(patches.Rectangle(bottomleft, data.grid.dx[i] * 2 * xnorm,
                                                       data.grid.dy[i] * 2 * ynorm, fill=False, edgecolor=gridcolor,
                                                       alpha=gridalpha))
                        plottedInd[ind] = i

            elif plane.strip().lower() == 'yx':

                for i in idata['cellID']:
                    if i not in plottedInd:
                        ind += 1
                        bottomleft = ((data.grid.y[i] - data.grid.dy[i]) * xnorm,
                                      (data.grid.x[i] - data.grid.dx[i]) * ynorm)
                        ax.add_patch(patches.Rectangle(bottomleft, data.grid.dy[i] * 2 * xnorm,
                                                       data.grid.dx[i] * 2 * ynorm, fill=False, edgecolor=gridcolor,
                                                       alpha=gridalpha))
                        plottedInd[ind] = i

            elif plane.strip().lower() == 'xz':

                for i in idata['cellID']:
                    if i not in plottedInd:
                        ind += 1
                        bottomleft = ((data.grid.x[i] - data.grid.dx[i]) * xnorm,
                                      (data.grid.z[i] - data.grid.dz[i]) * ynorm)
                        ax.add_patch(patches.Rectangle(bottomleft, data.grid.dx[i] * 2 * xnorm,
                                                       data.grid.dz[i] * 2 * ynorm, fill=False, edgecolor=gridcolor,
                                                       alpha=gridalpha))
                        plottedInd[ind] = i

            elif plane.strip().lower() == 'zx':

                for i in idata['cellID']:
                    if i not in plottedInd:
                        ind += 1
                        bottomleft = ((data.grid.z[i] - data.grid.dz[i]) * xnorm,
                                      (data.grid.x[i] - data.grid.dx[i]) * ynorm)
                        ax.add_patch(patches.Rectangle(bottomleft, data.grid.dz[i] * 2 * xnorm,
                                                       data.grid.dx[i] * 2 * ynorm, fill=False, edgecolor=gridcolor,
                                                       alpha=gridalpha))
                        plottedInd[ind] = i

            elif plane.strip().lower() == 'yz':

                for i in idata['cellID']:
                    if i not in plottedInd:
                        ind += 1
                        bottomleft = ((data.grid.y[i] - data.grid.y[i]) * xnorm,
                                      (data.grid.z[i] - data.grid.dz[i]) * ynorm)
                        ax.add_patch(patches.Rectangle(bottomleft, data.grid.dy[i] * 2 * xnorm,
                                                       data.grid.dz[i] * 2 * ynorm, fill=False, edgecolor=gridcolor,
                                                       alpha=gridalpha))
                        plottedInd[ind] = i

            elif plane.strip().lower() == 'zy':

                for i in idata['cellID']:
                    if i not in plottedInd:
                        ind += 1
                        bottomleft = ((data.grid.z[i] - data.grid.dz[i]) * xnorm,
                                      (data.grid.y[i] - data.grid.dy[i]) * ynorm)
                        ax.add_patch(patches.Rectangle(bottomleft, data.grid.dz[i] * 2 * xnorm,
                                                       data.grid.dy[i] * 2 * ynorm, fill=False, edgecolor=gridcolor,
                                                       alpha=gridalpha))
                        plottedInd[ind] = i

        else:

            if plane.strip().lower() == 'xy':
                px = data.grid.xi * xnorm
                if (data.grid.crd_sys == 'sph') & lattitude:
                    py = (np.pi / 2. - data.grid.yi) * ynorm
                else:
                    py = data.grid.yi * ynorm

            elif plane.strip().lower() == 'yx':
                py = data.grid.xi * xnorm
                if (data.grid.crd_sys == 'sph') & lattitude:
                    px = (np.pi / 2. - data.grid.yi) * ynorm
                else:
                    px = data.grid.yi * ynorm

            elif plane.strip().lower() == 'xz':
                px = data.grid.xi * xnorm
                py = data.grid.zi * znorm

            elif plane.strip().lower() == 'zx':
                py = data.grid.xi * xnorm
                px = data.grid.zi * znorm

            elif plane.strip().lower() == 'yz':
                if (data.grid.crd_sys == 'sph') & lattitude:
                    px = (np.pi / 2. - data.grid.yi) * ynorm
                else:
                    px = data.grid.yi * ynorm
                py = data.grid.zi * znorm

            elif plane.strip().lower() == 'zy':
                if (data.grid.crd_sys == 'sph') & lattitude:
                    py = (np.pi / 2. - data.grid.yi) * ynorm
                else:
                    py = data.grid.yi * ynorm
                px = data.grid.zi * znorm

            for ix in range(data.grid.nxi):
                ax.add_line(ml.Line2D((px[ix], px[ix]), (py[0], py[-1]), color=gridcolor, alpha=gridalpha))
            for iy in range(data.grid.nyi):
                ax.add_line(ml.Line2D((px[0], px[-1]), (py[iy], py[iy]), color=gridcolor, alpha=gridalpha))

    #
    # Set the axis limits
    #
    if not octree:
        if len(xlim) == 0:
            xlim = (plot_x[0], plot_x[-1])
        if len(ylim) == 0:
            ylim = (plot_y[0], plot_y[-1])

    plt.xlim(xlim[0], xlim[1])
    plt.ylim(ylim[0], ylim[1])
    #
    # Set the axis labels
    #
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)


def plotDustOpac(opac=None, var='kabs', idust=0, ax=None, xlabel=None, ylabel=None, fmt='', **kwargs):
    """
    Plots the dust opacity as a function of wavelength

    Parameters
    ----------
    opac            : radmc3dDustOpac
                      Dust opacity container

    var             : {'kabs', 'ksca', 'kext', 'g'}
                      Variable to be plotted

    idust           : int, optional
                      Dust index in opac to be plotted (default=0)

    ax              : Axis, optional
                      Matplotlib axis to plot on. If not set the current axis will be used.

    xlabel          : str, optional
                      Label of the x-axis

    ylabel          : str, optional
                      Label of the y-axis

    fmt             : str, optional
                      Format of the plotted line. The same as the third non-keyword argument of matplotlib.pyplot.plot()

    Keyword Arguments:
                    Any further keyword argument that will be passed to matplotlib.pyplot.plot()

    Returns
    -------
    The returned list by matplotlib.pyplot.plot()

    """
    if ax is None:
        ax = plt.gca()

    x = opac.wav[idust]

    if not isinstance(idust, int):
        idust = int(idust)

    if var.lower().strip() == 'kabs':
        if opac.kabs[idust].size != x.size:
            msg = 'opac.kabs['+("%d" % idust)+'] has a different number of elements than opac.wav['+("%d" % idust)+']'
            raise ValueError(msg)
        y = opac.kabs[idust]
        if ylabel is None:
            ylabel = r'$\kappa_{\rm abs}$ [cm$^2$/g]'

    elif var.lower().strip() == 'ksca':
        if opac.ksca[idust].size != x.size:
            msg = 'opac.ksca['+("%d" % idust)+'] has a different number of elements than opac.wav['+("%d" % idust)+']'
            raise ValueError(msg)
        y = opac.ksca[idust]
        if ylabel is None:
            ylabel = r'$\kappa_{\rm sca}$ [cm$^2$/g]'

    elif var.lower().strip() == 'kext':
        if opac.kabs[idust].size != x.size:
            msg = 'opac.kabs['+("%d" % idust)+'] has a different number of elements than opac.wav['+("%d" % idust)+']'
            raise ValueError(msg)
        if opac.ksca[idust].size != x.size:
            msg = 'opac.ksca['+("%d" % idust)+'] has a different number of elements than opac.wav['+("%d" % idust)+']'
            raise ValueError(msg)
        y = opac.kabs[idust] + opac.ksca[idust]
        if ylabel is None:
            ylabel = r'$\kappa_{\rm ext}$ [cm$^2$/g]'

    elif var.lower().strip() == 'g':
        if opac.phase_g[idust].size != x.size:
            msg = 'opac.g['+("%d" % idust)+'] has a different number of elements than opac.wav['+("%d" % idust)+']'
            raise ValueError(msg)
        y = opac.phase_g[idust]
        if ylabel is None:
            ylabel = r'g$_{\rm HG}$'
    else:
        msg = 'Unknown variable to plot ' + var
        raise ValueError(msg)

    p = ax.plot(x, y, fmt, **kwargs)

    if xlabel is None:
        xlabel = r'$\lambda$ [$\mu$m]'

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return p


def plotScatmat(opac=None, var='z11', idust=0, iwav=None, wav=None, xvar='ang', iang=None, ang=None, ax=None,
                xlabel=None, ylabel=None, title=None, fmt='', **kwargs):
    """
    Plots the scattering matrix elements either as a function of scattering angle at a specific wavelength (default)
    or as a function of wavelength at a specific scattering angle

    Parameters
    ----------

    opac            : radmc3dDustOpac
                      Dust opacity container

    var             : {'kabs', 'ksca', 'kext', 'g'}
                      Variable to be plotted

    idust           : int, optional
                      Dust index in opac to be plotted (default=0)

    iwav            : int, optional
                      Wavelength index (used only if xvar='ang') to be plotted.

    wav             : float, optional
                      Wavelength at which the plot should be made (used only if xvar='ang'). In practice, instead of
                      interpolating to this wavelength, the nearest wavelength in the wavelength grid will be used.

    xvar            : {'ang', 'wav'}
                      Variable for plotting the scattering matrix elements against (default='ang')

    iang            : int, optional
                      Scattering angle index (used only if xvar='wav')

    ang             : int, optional
                      Scattering angle in degrees at which the plot should be made (used only if xvar='wav').
                      In practice, instead of interpolating to this scattering angle, the nearest angle in the angular
                      grid will be used.

    ax              : Axis, optional
                      Matplotlib axis to plot on. If not set the current axis will be used.

    xlabel          : str, optional
                      Label of the x-axis

    ylabel          : str, optional
                      Label of the y-axis

    title           : str, optional
                      Title of the plot. If not set either the wavelength (for xvar='ang') or the scattering angle
                      (for xvar='wav') the plot was made at.

    fmt             : str, optional
                      Format of the plotted line. The same as the third non-keyword argument of matplotlib.pyplot.plot()

    Keyword Arguments:
                    Any further keyword argument that will be passed to matplotlib.pyplot.plot()

    Returns
    -------
    The returned list by matplotlib.pyplot.plot()
    """

    if ax is None:
        ax = plt.gca()

    if xvar.lower().strip() == 'ang':
        plot_type = 1
        x = opac.scatang[idust]
        if xlabel is None:
            xlabel = 'Scattering angle [deg]'
        # Get wavelength
        if iwav is None:
            if wav is None:
                msg = 'Either the iwav or the wav keywords should be set to indicate which wavelength the ' \
                      'plot should be made at'
                raise ValueError(msg)
            else:
                if (wav < opac.wav[idust].min()) | (wav > opac.wav[idust].max()):
                    msg = 'The requested wavelength is outside of the wavelength range contained in opac \n'
                    msg += ' min(wav) [micron]: '+("%e" % opac.wav[idust].min())+'\n'
                    msg += ' max(wav) [micron]: '+("%e" % opac.wav[idust].max())+'\n'
                    raise ValueError(msg)
                else:
                    iwav = abs(opac.wav[idust] - wav).argmin()
        else:
            if not isinstance(iwav, int):
                iang = int(iwav)

    elif xvar.lower().strip() == 'wav':
        plot_type = 2
        x = opac.wav[idust]
        if xlabel is None:
            xlabel = r'$\lambda$ [$\mu$m]'

        # Get scattering angle
        if iang is None:
            if ang is None:
                msg = 'Either the iang or the ang keywords should be set to indicate which scattering angle the ' \
                      'plot should be made at'
                raise ValueError(msg)
            else:
                if (ang < opac.scatang[idust].min()) | (ang > opac.scatang[idust].max()):
                    msg = 'The requested angle is outside of the range of scattering angles contained in opac \n'
                    msg += ' min(ang) [deg]: '+("%e" % opac.scatang[idust].min())+'\n'
                    msg += ' max(ang) [deg]: '+("%e" % opac.scatang[idust].max())+'\n'
                    raise ValueError(msg)
                else:
                    iang = abs(opac.scatang[idust] - ang).argmin()
        else:
            if not isinstance(iang, int):
                iang = int(iang)
    else:
        msg = 'Unknown variable in xvar ' + xvar
        raise ValueError(msg)

    if var.lower().strip() == 'z11':
        y = opac.z11[idust]
        if ylabel is None:
            ylabel = r'z$_{11}$'

    elif var.lower().strip() == 'z12':
        y = opac.z12[idust]
        if ylabel is None:
            ylabel = r'z$_{12}$'

    elif var.lower().strip() == 'z22':
        y = opac.z22[idust]
        if ylabel is None:
            ylabel = r'z$_{22}$'

    elif var.lower().strip() == 'z33':
        y = opac.z33[idust]
        if ylabel is None:
            ylabel = r'z$_{33}$'

    elif var.lower().strip() == 'z34':
        y = opac.z34[idust]
        if ylabel is None:
            ylabel = r'z$_{34}$'

    elif var.lower().strip() == 'z44':
        y = opac.z44[idust]
        if ylabel is None:
            ylabel = r'z$_{44}$'

    elif var.lower().strip() == 'linpol':
        y = -opac.z12[idust] / opac.z11[idust]
        if ylabel is None:
            ylabel = r'-z$_{12}$/z$_{11}$ (Fraction of linear polarisation)'

    else:
        msg = 'Unknown variable to plot ' + var
        raise ValueError(msg)

    if plot_type == 1:
        y = y[iwav, :]
        if title is None:
            title = r'$\lambda$ = ' + ("%.3f" % opac.wav[idust][iwav]) + r'$\mu$m'
    elif plot_type == 2:
        y = y[:, iang]
        if title is None:
            title = r'$\phi$ = ' + ("%.3f" % opac.scatang[idust][iang]) + r'$^\circ$'

    p = ax.plot(x, y, fmt, **kwargs)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    return p
