"""This module contains functions to get various radiation sources in RADMC-3D.
For help on the syntax or functionality of each function see the help of the individual functions
"""
from __future__ import absolute_import
from __future__ import print_function
import traceback
import shutil
import os

try:
    import numpy as np
except ImportError:
    np = None
    print(traceback.format_exc())

from . import natconst as nc


def getAtmModel(teff=0., logg=None, mstar=None, lstar=None, rstar=None, iwav=None, model='kurucz',
                modeldir=None, wmax=7.):
    """
    Interpolates the stellar model atmosphere on a pre-defined wavelength grid
    The model atmospheres are interpolated in logg and Teff and rebinned in wavelength 
    to the input wavelength grid and scaled to have the same luminosity as specified 


    Parameters
    ----------

    teff        : float
                  Effective temperature of the star

    logg        : float
                  Logarithm of the surface gravity of the star

    mstar       : float
                  Mass of the star in gramm

    lstar       : float, optional
                  Luminosity of the star (either lstar or rstar should be specified)

    rstar       : float, optional
                  Radius of the star (either lstar or rstar should be specified)

    iwav        : ndarray
                  Wavelength grid for which the stellar atmosphere model should be calculated

    model       : {'kurucz', 'nextgen', 'ames-dusty'}
                  Name of the model atmosphere family

    modeldir    : str       
                  Path to the atmosphere models
                
    wmax        : float
                  Maximum wavelength until the model atmosphere is used on the interpolated grid
                  Longwards of this wavelength the Rayleigh-Jeans approximation (lambda^(-2)) is used. 

    Returns
    -------

    Returns a dictionary with the following keys:

        * wav    : ndarray
                   Wavelength in micron (same as the input iwav keyword)

        * lnu    : ndarray
                   Monochromatic luminosity of the stellar atmosphere in erg/s/Hz

    """
    #
    # Get the stellar parameters
    #
    if rstar is None:
        rstar = np.sqrt(lstar / (4. * np.pi * nc.ss * teff**4.))
    else:
        lstar = 4. * np.pi * rstar**2 * nc.ss * teff**4

    if logg is None:
        if mstar is None:
            raise ValueError('Neither logg nor mstar is set. The atmosphere models depend on the stellar mass ' +
                             ' thus either logg or mstar should be set.')
        logg = np.log10(nc.gg * mstar / rstar**2)

    #
    # Get the atmosphere model
    #
    if model.strip().lower() == 'kurucz':
        if modeldir is None:
            modeldir = '/disk2/juhasz/Data/kurucz/'
        sp = getSpectrumKurucz(teff=teff, logg=logg, lstar=lstar, modeldir=modeldir)
        ilnu = rebinSpectrum(wav=sp['wav'], fnu=sp['lnu'], iwav=iwav)
    elif model.strip().lower() == 'nextgen':
        if modeldir is None:
            modeldir = '/disk2/juhasz/Data/NextGen/SPECTRA/'
        sp = getSpectrumNextGen(teff=teff, logg=logg, lstar=lstar, modeldir=modeldir)
        ilnu = rebinSpectrum(wav=sp['wav'], fnu=sp['lnu'], iwav=iwav)
    elif model.strip().lower() == 'ames-dusty':
        if modeldir is None:
            modeldir = '/disk2/juhasz/Data/Ames-Dusty/SPECTRA/'
        sp = getSpectrumAmesDusty(teff=teff, logg=logg, lstar=lstar, modeldir=modeldir)
        ilnu = rebinSpectrum(wav=sp['wav'], fnu=sp['lnu'], iwav=iwav)
    else:
        raise ValueError('Unknown model : "'+model+'" \n Allowed (implemented) models are "kurucz" or "nextgen".')

    ilnu = ilnu.clip(1e0, 1e99)
    ii = (iwav > wmax)
    if True in ii:
        f0 = 10.**np.interp(np.log10(wmax), np.log10(iwav), np.log10(ilnu))
        ilnu[ii] = f0 * (iwav[ii] / wmax)**(-2.)

    return {'wav': iwav, 'lnu': ilnu}


def getSpectrumKurucz(teff=None, logg=None, mstar=None, lstar=None, rstar=None, modeldir=None, wav=None):
    """
    Interpolates the Kurucz model atmospheres in logg and Teff 

    Parameters
    ----------

    teff        : float
                  Effective temperature of the star

    logg        : float
                  Logarithm of the surface gravity of the star

    mstar       : float, optional
                  Stellar mass in g (either logg or mstar and rstar should be specified)
                  
    lstar       : float, optional
                  Luminosity of the star (either lstar or rstar should be specified)

    rstar       : float, optional
                  Radius of the star (either lstar or rstar should be specified)
    
    modeldir    : str
                  Path to the Kurucz atmosphere model grid
    
    wav         : ndarray
                  Wavelength grid for which the stellar atmosphere model should be calculated

    Returns
    -------

    Returns a dictionary with the following keys:

        * wav     : ndarray
                    Wavelength in micron (same as the input wav keyword)

        * lnu     : ndarray
                    Monochromatic luminosity of the stellar atmosphere in erg/s/Hz

        * lnucont : ndarray
                    Monochromatic luminosity of the continuum stellar atmosphere in erg/s/Hz
    """

    #
    # Sanity check
    #
    if teff is None:
        raise ValueError('Unknown teff. Kurucz spectrum cannot be calculated without the stellar effective'
                         + 'temperature.')
    if logg is None:
        if mstar is None:
            raise ValueError('Unknown logg and mstar. For Kurucz atmosphere models either logg or mstar and '
                             + ' rstar/lstar should be specified')
        else:
            if rstar is None:
                if lstar is None:
                    raise ValueError('Unknown logg, rstar, lstar thus logg cannot be calculated. For Kurucz '
                    + 'atmosphere models either logg should be given or a combination of mstar, lstar, rstar '
                    + 'should be specified from which logg can be calculated for a given effective temperature. ')
                else:
                    rstar = np.sqrt(lstar / (4. * np.pi * nc.ss * teff**4.))
            else:
                lstar = 4. * np.pi * rstar**2 * nc.ss * teff**4

        logg = np.log10(nc.gg * mstar / rstar**2)
    else:
        if rstar is None:
            if lstar is None:
                raise ValueError('Unknown rstar and lstar thus luminosity cannot be calculated. For Kurucz '
                                 + 'atmosphere models either lstar or rstar should be given')
            else:
                rstar = np.sqrt(lstar / (4. * np.pi * nc.ss * teff**4.))
        else:
            lstar = 4. * np.pi * rstar**2 * nc.ss * teff**4

    mstar = 10.**logg * rstar**2 / nc.gg

    assert (teff >= 0.), 'Negative effective temperature!'
    assert (mstar >= 0.), 'Negative stellar mass!'
    assert (lstar >= 0.), 'Negative stellar luminosity!'
    assert (rstar >= 0.), 'Negative stellar radius!'

    print('-------------------------------------------')
    print('Interpolating in Kurucz model atmospheres')
    print('Stellar parameters: ')
    print('-------------------------------------------')
    print('Teff [K]          : ', teff)
    print('Radius [Rsun]     : ', rstar / nc.rs)
    print('Luminosity [Lsun] : ', lstar / nc.ls)
    print('Mass [Msun]       : ', mstar / nc.ms)
    print('logg              : ', logg)
    print('-------------------------------------------')

    dum = readKuruczGrid(fname=modeldir + '/fp00k2.pck')
    #
    # Bracket in Teff
    #
    teff_grid = np.unique(dum['teff'])

    ii = abs(teff_grid - teff).argmin()
    idt1 = ii
    if teff_grid[ii] > teff:
        idt1 = ii - 1
    else:
        idt1 = ii
    idt2 = idt1 + 1

    #
    # Bracket in Logg
    #
    ii = (dum['teff'] == teff_grid[idt1])
    logg_grid_lower = dum['logg'][ii]
    ii = (dum['teff'] == teff_grid[idt2])
    logg_grid_upper = dum['logg'][ii]

    ii = abs(logg_grid_lower - logg).argmin()
    if logg < logg_grid_lower[0]:
        idg1 = -1
        idg2 = 0
    elif logg > logg_grid_lower[-1]:
        idg1 = logg_grid_lower.shape[0] - 1
        idg2 = -1
    else:
        idg1 = ii
        if logg_grid_lower[ii] > logg:
            idg1 = ii - 1
        else:
            idg1 = ii
        idg2 = idg1 + 1

    idgl1 = idg1
    idgl2 = idg2

    ii = abs(logg_grid_upper - logg).argmin()
    if logg < logg_grid_upper[0]:
        idg1 = -1
        idg2 = 0
    elif logg > logg_grid_upper[-1]:
        idg1 = logg_grid_upper.shape[0] - 1
        idg2 = -1
    else:
        idg1 = ii
        if logg_grid_upper[ii] > logg:
            idg1 = ii - 1
        else:
            idg1 = ii
        idg2 = idg1 + 1

    idgu1 = idg1
    idgu2 = idg2

    #
    # Check if we need to do a 3point bilinear interpolation
    #
    if ((idgl1 < 0) | (idgl2 < 0) | (idgu1 < 0) | (idgu2 < 0)):
        x = []
        y = []
        sp = []
        spc = []
        if idgl1 >= 0:
            x.append(teff_grid[idt1])
            y.append(logg_grid_lower[idgl1])
            ii = ((dum['teff'] == teff_grid[idt1]) & (dum['logg'] == logg_grid_lower[idgl1]))
            sp.append(np.squeeze(dum['inu'][ii, :]))
            spc.append(np.squeeze(dum['inucont'][ii, :]))
        if idgl2 >= 0:
            x.append(teff_grid[idt1])
            y.append(logg_grid_lower[idgl2])
            ii = ((dum['teff'] == teff_grid[idt1]) & (dum['logg'] == logg_grid_lower[idgl2]))
            sp.append(np.squeeze(dum['inu'][ii, :]))
            spc.append(np.squeeze(dum['inucont'][ii, :]))
        if idgu1 >= 0:
            x.append(teff_grid[idt2])
            y.append(logg_grid_upper[idgu1])
            ii = ((dum['teff'] == teff_grid[idt2]) & (dum['logg'] == logg_grid_upper[idgu1]))
            sp.append(np.squeeze(dum['inu'][ii, :]))
            spc.append(np.squeeze(dum['inucont'][ii, :]))
        if idgu2 >= 0:
            x.append(teff_grid[idt2])
            y.append(logg_grid_upper[idgu2])
            ii = ((dum['teff'] == teff_grid[idt2]) & (dum['logg'] == logg_grid_upper[idgu2]))
            sp.append(np.squeeze(dum['inu'][ii, :]))
            spc.append(np.squeeze(dum['inucont'][ii, :]))

        if len(x) != 3:
            msg = 'Something went wrong in the interpolation of the stellar atmosphere models..\n ' \
                  + ' for a 3 point bilinear interpolation ' + ("%d" % len(x))+' indices were found.'
            raise ValueError(msg)

        else:
            print('Bracketed spectrum with Teff : ', teff, ' and logg : ', logg)
            print('Teff grid : ', x)
            print('Logg grid : ', y)

        c1 = ((y[1] - y[2]) * (teff - x[2]) + (x[2] - x[1]) * (logg - y[2])) / (
        (y[1] - y[2]) * (x[0] - x[2]) + (x[2] - x[1]) * (y[0] - y[2]))
        c2 = ((y[2] - y[0]) * (teff - x[2]) + (x[0] - x[2]) * (logg - y[2])) / (
        (y[1] - y[2]) * (x[0] - x[2]) + (x[2] - x[1]) * (y[0] - y[2]))
        c3 = 1. - c1 - c2

        lnu = c1 * sp[0] + c2 * sp[1] + c3 * sp[2]
        lnucont = c1 * spc[0] + c2 * spc[1] + c3 * spc[2]


    else:
        print('Bracketed spectrum with Teff : ', teff, ' and logg : ', logg)
        print('Teff grid : ', teff_grid[idt1], teff_grid[idt2])
        print('Logg grid : ', logg_grid_lower[idgl1], logg_grid_lower[idgl2])
        #
        # Do the standard four point bilinear interpolation
        #
        ii = ((dum['teff'] == teff_grid[idt1]) & (dum['logg'] == logg_grid_lower[idgl1]))
        sp11 = np.squeeze(dum['inu'][ii, :])
        ii = ((dum['teff'] == teff_grid[idt1]) & (dum['logg'] == logg_grid_lower[idgl2]))
        sp12 = np.squeeze(dum['inu'][ii, :])
        ii = ((dum['teff'] == teff_grid[idt2]) & (dum['logg'] == logg_grid_upper[idgu1]))
        sp22 = np.squeeze(dum['inu'][ii, :])
        ii = ((dum['teff'] == teff_grid[idt2]) & (dum['logg'] == logg_grid_upper[idgu2]))
        sp21 = np.squeeze(dum['inu'][ii, :])

        c11 = (teff_grid[idt2] - teff) * (logg_grid_upper[idgu2] - logg)
        c12 = (teff_grid[idt2] - teff) * (logg - logg_grid_upper[idgu1])
        c22 = (teff - teff_grid[idt1]) * (logg - logg_grid_lower[idgl1])
        c21 = (teff - teff_grid[idt1]) * (logg_grid_lower[idgl2] - logg)
        c00 = 1. / ((teff_grid[idt2] - teff_grid[idt1]) * (logg_grid_lower[idgl2] - logg_grid_lower[idgl1]))

        lnu = c00 * (c11 * sp11 + c12 * sp12 + c22 * sp22 + c21 * sp21)
        lnucont = c00 * (c11 * sp11 + c12 * sp12 + c22 * sp22 + c21 * sp21)

    nu = nc.cc / dum['wav'] * 1e4
    lum = (0.5 * abs(nu[1:] - nu[:-1]) * (lnu[1:] + lnu[:-1])).sum()
    lnu *= lstar / lum
    lnucont *= lstar / lum

    return {'wav': dum['wav'], 'lnu': lnu, 'lnucont': lnucont}

def readKuruczGrid(fname=''):
    """
    Reads the Kurucz model atmosphere grid. It reads a whole file from the Kurucz grid that contains
    a 2D grid of atmosphere models for Teff and logg.  

    Parameters
    ----------

    fname       : str
                  File name containing the Kurucz model atmosphere (e.g. fp00k2.pck)


    Returns
    -------

    Returns a dictionary with the following keys:

        * wav     : ndarray
                    Wavelength in micron 

        * nwav    : int
                    Number of wavelength points

        * inu     : list
                    Each element of the list contains an ndarray with the specific intensity of the stellar 
                    model atmosphere in erg/s/cm/cm/Hz/ster

        * inucont : list
                    Each element of the list contains an ndarray with the specific intensity of the stellar 
                    model atmosphere continuum in erg/s/cm/cm/Hz/ster

        * teff    : list
                    Contains the Teff grid of the model grid 

        * logg    : list
                    Contains the logg grid of the model grid

    """

    with open(fname, 'r') as rfile:
        #
        # Skip the program part
        #
        for i in range(22):
            dum = rfile.readline()

        #
        # Read the wavelength grid
        #
        wav = []
        n = 10
        for i in range(153):
            dum = rfile.readline().split()
            for j in range(len(dum)):
                wav.append(float(dum[j]))

        #
        # Convert the wavelength in Angstrom to micron
        #
        wav = np.array(wav) * 1e-3
        #
        # Now read the grid of spectra
        #
        nwav = wav.shape[0]
        tgrid_list = []
        logg_list = []
        inu_list = []
        inucont_list = []

        #
        # Read the first section header
        #
        dum = rfile.readline()
        while dum.strip() != '':
            # print '>>>> ', dum, len(dum.strip())
            sdum = dum.split()
            tgrid_list.append(float(sdum[1]))
            logg_list.append(float(sdum[3]))

            #
            # Read the stellar spectrum
            #
            arr = []
            for i in range(152):
                dum = rfile.readline()
                for j in range(8):
                    arr.append(float(dum[j * n:(j + 1) * n]))
            dum = rfile.readline()
            for j in range(5):
                arr.append(float(dum[j * n:(j + 1) * n]))
            inu_list.append(np.array(arr))
            #
            # Read the continuum spectrum
            #
            arr = []
            for i in range(152):
                dum = rfile.readline()
                for j in range(8):
                    arr.append(float(dum[j * n:(j + 1) * n]))
            dum = rfile.readline()
            for j in range(5):
                arr.append(float(dum[j * n:(j + 1) * n]))
            inucont_list.append(np.array(arr))

            #
            # Read the next section header
            #
            dum = rfile.readline()

    teff_grid = np.array(tgrid_list)
    logg_grid = np.array(logg_list)
    inu = np.array(inu_list)
    inucont = np.array(inucont_list)

    return {'wav': wav, 'inu': inu, 'inucont': inucont, 'teff': teff_grid, 'logg': logg_grid, 'nwav': nwav}


def readAmesDustySpectrum(fname=''):
    """
    Reads the Ames-Dusty model atmosphere

    Parameters
    ----------
    fname       : str
                  File name containing the Ames-Dusty model atmosphere (e.g. lte30-3.5-0.0.AMES-dusty.7)

    Returns
    -------

    Returns a dictionary with the following keys:

        * wav     : ndarray
                    Wavelength in micron

        * nwav    : int
                    Number of wavelength points

        * inu     : ndarray
                    Specific intensity of the stellar model atmosphere in erg/s/cm/cm/Hz/ster

        * bnu     : list
                    Specific intensity of a blackbody stellar model atmosphere with the same luminosity
                    and the same effective temperature in erg/s/cm/cm/Hz/ster

        * teff    : float
                    Effective temperature of the model

        * logg    : float
                    Logarithm of the surface gravity of the model

        * mph     : float
                    Metallicity of the atmosphere model

    """
    print('Reading : ', fname)

    # Get the effective temperature, logg and metallicity from the file name
    ind = fname.find('lte')
    fname_tags = fname[ind+3:ind+13].split('-')
    teff = np.float(fname_tags[0]) * 100.
    logg = np.float(fname_tags[1]) * 100.
    mph = np.float(fname_tags[2]) * 100.

    wav = []
    inu = []
    bnu = []
    with open(fname, 'r') as rfile:
        dum = rfile.readline()
        while dum != '':
            dum = str(dum).replace('D', 'E')
            sdum = dum.split()
            wav.append(np.float(sdum[0]))
            inu.append(np.float(sdum[1]))
            bnu.append(np.float(sdum[2]))
            dum = rfile.readline()

        wav = np.array(wav)
        inu = np.array(inu)
        bnu = np.array(bnu)
        ii = wav.argsort()

        wav = wav[ii]
        inu = inu[ii]
        bnu = bnu[ii]

    # "Decode" the intensity arrays
    inu = 10.**(inu - 8.0) * wav
    bnu = 10.**(bnu - 8.0) * wav

    # Convert the wavelength to micron from Angstrom
    wav /= 1e4
    nwav = wav.shape[0]

    return {'teff': teff, 'logg': logg, 'mph': mph, 'nwav': nwav, 'wav': wav, 'inu': inu, 'bnu': bnu}


def readNextGenSpectrum(fname=''):
    """
    Reads the NextGen model atmosphere.

    Parameters
    ----------

    fname       : str
                  File name containing the NextGen model atmosphere (e.g. nlte98-4.5-2.5.NextGen.spec)


    Returns
    -------

    Returns a dictionary with the following keys:

        * wav     : ndarray
                    Wavelength in micron 

        * nwav    : int
                    Number of wavelength points

        * inu     : ndarray
                    Specific intensity of the stellar model atmosphere in erg/s/cm/cm/Hz/ster

        * bnu     : list
                    Specific intensity of a blackbody stellar model atmosphere with the same luminosity
                    and the same effective temperature in erg/s/cm/cm/Hz/ster

        * teff    : float
                    Effective temperature of the model

        * logg    : float
                    Logarithm of the surface gravity of the model

        * mph     : float
                    Metallicity of the atmosphere model

    """

    print('Reading : ', fname)

    with open(fname, 'r') as rfile:
        dum = rfile.readline()
        sdum = dum.split()
        teff = float(sdum[0])
        logg = float(sdum[1])
        mph = float(sdum[2])
        dum = rfile.readline()
        nwav = float(dum.split()[0])

        bigline = []
        dum = rfile.readline()
        while dum.strip() != '':
            sdum = dum.split()
            for i in range(len(sdum)):
                bigline.append(float(sdum[i]))
            dum = rfile.readline()

    bigline = np.array(bigline)
    # Convert wavelength from angstrom to micron
    wav = bigline[:nwav] / 1e4
    inu = bigline[nwav:2 * nwav]
    bnu = bigline[nwav * 2:nwav * 3]

    ii = wav.argsort()
    wav = wav[ii]
    inu = inu[ii] * 1e-8 * wav * 1e4 / np.pi / (29979245800.0 / wav * 1e4)
    bnu = bnu[ii] * 1e-8 * wav * 1e4 / np.pi / (29979245800.0 / wav * 1e4)

    #
    # The unit is now erg/s/cm/Hz/ster
    #

    return {'teff': teff, 'logg': logg, 'mph': mph, 'nwav': nwav, 'wav': wav, 'inu': inu, 'bnu': bnu}


def rebinSpectrum(wav=None, fnu=None, iwav=None):
    """
    Rebins the spectrum to a coarser wavelength grid

    Parameters
    ----------

    wav     : ndarray
              Wavelength grid of the spectrum to be rebinned

    fnu     : ndarray
              Wavelength dependent spectrum (e.g. specific intensity, monochromatic luminosity etc) to be rebbinned

    iwav    : ndarray
              Wavelength grid onto which the spectrum should be rebinned (it is assumed to be logarithmically spaced)

    Returns
    -------

    Returns an ndarray containing the spectrum rebinned to the input wavelength grid

    """
    if wav is None:
        raise ValueError('Unknown wav in rebin. The spectrum cannot be rebinned without the wavelength grid.')
    if fnu is None:
        raise ValueError('Unknown fnu. The spectrum cennot be rebinned.')
    if fnu is None:
        raise ValueError('Unknown iwav. The spectrum cannot be rebinned without knowing the wavelength grid to '
        + 'rebin to.')
    #
    # Generate the interface grid
    #
    iiwav = np.zeros(iwav.shape[0] + 1, dtype=float)
    iiwav[1:-1] = np.sqrt(iwav[1:] * iwav[:-1])
    iiwav[0] = iwav[0]**2 / iwav[1]
    iiwav[-1] = iwav[-1]**2 / iwav[-2]

    ifnu = np.zeros(iiwav.shape[0] - 1, dtype=float)
    for i in range(iiwav.shape[0] - 1):
        ii = ((wav > iiwav[i]) & (wav <= iiwav[i + 1]))
        if ii.__contains__(True):
            x = wav[ii]
            y = fnu[ii]
            ifnu[i] = (0.5 * (x[1:] - x[:-1]) * (y[1:] + y[:-1])).sum() / (iiwav[i + 1] - iiwav[i])

    return ifnu


def getSpectrumNextGen(teff=None, logg=None, mstar=None, lstar=None, rstar=None, modeldir=None, wav=None):
    """
    Interpolates the NextGen model atmospheres in logg and Teff

    Parameters
    ----------

    teff        : float
                  Effective temperature of the star

    logg        : float
                  Logarithm of the surface gravity of the star

    mstar       : float, optional
                  Stellar mass in g (either logg or mstar and rstar should be specified)

    lstar       : float, optional
                  Luminosity of the star (either lstar or rstar should be specified)

    rstar       : float, optional
                  Radius of the star (either lstar or rstar should be specified)

    modeldir    : str
                  Path to the NextGen atmosphere model grid

    wav         : ndarray
                  Wavelength grid for which the stellar atmosphere model should be calculated

    Returns
    -------

    Returns a dictionary with the following keys:

        * wav     : ndarray
                    Wavelength in micron (same as the input wav keyword)

        * lnu     : ndarray
                    Monochromatic luminosity of the stellar atmosphere in erg/s/Hz

        * bnu     : ndarray
                    Monochromatic luminosity of a blackbody stellar atmosphere with the same luminosity and
                    effective temperature as the stellar model in erg/s/Hz

    """

    #
    # Sanity check
    #
    if teff is None:
        raise ValueError('Unknown teff. NextGen spectrum cannot be calculated without the stellar effective'
                         + 'temperature.')
    if logg is None:
        if mstar is None:
            raise ValueError('Unknown logg and mstar. For NextGen atmosphere models either logg or mstar and '
                             + ' rstar/lstar should be specified')
        else:
            if rstar is None:
                if lstar is None:
                    raise ValueError('Unknown logg, rstar, lstar thus logg cannot be calculated. For NextGen '
                    + 'atmosphere models either logg should be given or a combination of mstar, lstar, rstar '
                    + 'should be specified from which logg can be calculated for a given effective temperature. ')
                else:
                    rstar = np.sqrt(lstar / (4. * np.pi * nc.ss * teff**4.))
            else:
                lstar = 4. * np.pi * rstar**2 * nc.ss * teff**4

        logg = np.log10(nc.gg * mstar / rstar**2)
    else:
        if rstar is None:
            if lstar is None:
                raise ValueError('Unknown rstar and lstar thus luminosity cannot be calculated. For NextGen '
                                 + 'atmosphere models either lstar or rstar should be given')
            else:
                rstar = np.sqrt(lstar / (4. * np.pi * nc.ss * teff**4.))
        else:
            lstar = 4. * np.pi * rstar**2 * nc.ss * teff**4

    mstar = 10.**logg * rstar**2 / nc.gg

    assert (teff >= 0.), 'Negative effective temperature!'
    assert (mstar >= 0.), 'Negative stellar mass!'
    assert (lstar >= 0.), 'Negative stellar luminosity!'
    assert (rstar >= 0.), 'Negative stellar radius!'

    print('-------------------------------------------')
    print('Interpolating in NextGen model atmospheres')
    print('Stellar parameters: ')
    print('-------------------------------------------')
    print('Teff [K]          : ', teff)
    print('Radius [Rsun]     : ', rstar / nc.rs)
    print('Luminosity [Lsun] : ', lstar / nc.ls)
    print('Mass [Msun]       : ', mstar / nc.ms)
    print('logg              : ', logg)
    print('-------------------------------------------')

    teff_grid = np.append((np.arange(31.) + 9.), (np.arange(30.) * 2. + 40.))
    logg_grid = np.arange(6.) * 0.5 + 3.5

    # Bracket the input teff and logg values
    ii = abs(teff_grid - teff / 100.).argmin()
    idt1 = ii
    if teff_grid[ii] > teff / 100.:
        idt1 = ii - 1
    else:
        idt1 = ii
    idt2 = idt1 + 1

    ii = abs(logg_grid - logg).argmin()
    if logg < logg_grid[0]:
        idg1 = 0
        idg2 = 0
    elif logg > logg_grid[-1]:
        idg1 = logg_grid.shape[0] - 1
        idg2 = logg_grid.shape[0] - 1
    else:
        idg1 = ii
        if logg_grid[ii] > logg:
            idg1 = ii - 1
        else:
            idg1 = ii
        idg2 = idg1 + 1

    print('Bracketing spectrum : ', teff, logg)
    print('Teff  : ', teff_grid[idt1], teff_grid[idt2])
    print('log(g): ', logg_grid[idg1], logg_grid[idg2])

    # Generate the spectral file names
    mph = '0.0'
    fname_11 = 'lte' + ("%d" % teff_grid[idt1]) + '-' + ("%.1f" % logg_grid[idg1]) + '-' + mph + '.NextGen.spec.gz'
    fname_12 = 'lte' + ("%d" % teff_grid[idt1]) + '-' + ("%.1f" % logg_grid[idg2]) + '-' + mph + '.NextGen.spec.gz'
    fname_22 = 'lte' + ("%d" % teff_grid[idt2]) + '-' + ("%.1f" % logg_grid[idg2]) + '-' + mph + '.NextGen.spec.gz'
    fname_21 = 'lte' + ("%d" % teff_grid[idt2]) + '-' + ("%.1f" % logg_grid[idg1]) + '-' + mph + '.NextGen.spec.gz'

    # Create a directory
    if os.path.exists("./tmp"):
        shutil.rmtree("./tmp")

    os.system('mkdir tmp')
    os.system('cp -v ' + modeldir + '/' + fname_11 + ' ./tmp')
    os.system('cp -v ' + modeldir + '/' + fname_12 + ' ./tmp')
    os.system('cp -v ' + modeldir + '/' + fname_22 + ' ./tmp')
    os.system('cp -v ' + modeldir + '/' + fname_21 + ' ./tmp')

    # Unzip the files
    os.chdir('./tmp')
    os.system('gunzip ' + fname_11)
    os.system('gunzip ' + fname_12)
    os.system('gunzip ' + fname_22)
    os.system('gunzip ' + fname_21)
    os.chdir('../')

    # Read the spectra
    sp11 = readNextGenSpectrum(fname='./tmp/' + fname_11[:-3])
    sp12 = readNextGenSpectrum(fname='./tmp/' + fname_12[:-3])
    sp22 = readNextGenSpectrum(fname='./tmp/' + fname_22[:-3])
    sp21 = readNextGenSpectrum(fname='./tmp/' + fname_21[:-3])

    # Do the interpolation
    c11 = (teff_grid[idt2] - teff / 100.) * (logg_grid[idg2] - logg)
    c12 = (teff_grid[idt2] - teff / 100.) * (logg - logg_grid[idg1])
    c22 = (teff / 100. - teff_grid[idt1]) * (logg - logg_grid[idg1])
    c21 = (teff / 100. - teff_grid[idt1]) * (logg_grid[idg2] - logg)
    c00 = 1. / ((teff_grid[idt2] - teff_grid[idt1]) * (logg_grid[idg2] - logg_grid[idg1]))

    lnu = c00 * (c11 * sp11['inu'] + c12 * sp12['inu'] + c22 * sp22['inu'] + c21 * sp21['inu'])
    bnu = c00 * (c11 * sp11['bnu'] + c12 * sp12['bnu'] + c22 * sp22['bnu'] + c21 * sp21['bnu'])

    shutil.rmtree('./tmp')

    #
    # Scale the spectrum to give the same luminosity as required
    #
    nu = nc.cc / sp11['wav'] * 1e4
    lum = (0.5 * abs(nu[1:] - nu[:-1]) * (lnu[1:] + lnu[:-1])).sum()
    lnu *= lstar / lum

    lum = (0.5 * abs(nu[1:] - nu[:-1]) * (bnu[1:] + bnu[:-1])).sum()
    bnu *= lstar / lum

    return {'wav': sp11['wav'], 'lnu': lnu, 'bnu': bnu}

def getSpectrumAmesDusty(teff=None, logg=None, mstar=None, lstar=None, rstar=None, modeldir=None, wav=None):
    """
    Interpolates the Ames-Dusty model atmospheres in logg and Teff

    Parameters
    ----------

    teff        : float
                  Effective temperature of the star

    logg        : float
                  Logarithm of the surface gravity of the star

    mstar       : float, optional
                  Stellar mass in g (either logg or mstar and rstar should be specified)

    lstar       : float, optional
                  Luminosity of the star (either lstar or rstar should be specified)

    rstar       : float, optional
                  Radius of the star (either lstar or rstar should be specified)

    modeldir    : str
                  Path to the AmesDusty atmosphere model grid

    wav         : ndarray
                  Wavelength grid for which the stellar atmosphere model should be calculated

    Returns
    -------

    Returns a dictionary with the following keys:

        * wav     : ndarray
                    Wavelength in micron (same as the input wav keyword)

        * lnu     : ndarray
                    Monochromatic luminosity of the stellar atmosphere in erg/s/Hz

        * bnu     : ndarray
                    Monochromatic luminosity of a blackbody stellar atmosphere with the same luminosity and
                    effective temperature as the stellar model in erg/s/Hz
    """

    #
    # Sanity check
    #
    if teff is None:
        raise ValueError('Unknown teff. Ames-Dusty spectrum cannot be calculated without the stellar effective'
                         + 'temperature.')
    if logg is None:
        if mstar is None:
            raise ValueError('Unknown logg and mstar. For Ames-Dusty atmosphere models either logg or mstar and '
                             + ' rstar/lstar should be specified')
        else:
            if rstar is None:
                if lstar is None:
                    raise ValueError('Unknown logg, rstar, lstar thus logg cannot be calculated. For Ames-Dusty '
                                     + 'atmosphere models either logg should be given or a combination of mstar, '
                                     + 'lstar, rstar should be specified from which logg can be calculated for a '
                                     + 'given effective temperature. ')
                else:
                    rstar = np.sqrt(lstar / (4. * np.pi * nc.ss * teff**4.))
            else:
                lstar = 4. * np.pi * rstar**2 * nc.ss * teff**4

        logg = np.log10(nc.gg * mstar / rstar**2)
    else:
        if rstar is None:
            if lstar is None:
                raise ValueError('Unknown rstar and lstar thus luminosity cannot be calculated. For Ames-Dusty '
                                 + 'atmosphere models either lstar or rstar should be given')
            else:
                rstar = np.sqrt(lstar / (4. * np.pi * nc.ss * teff**4.))
        else:
            lstar = 4. * np.pi * rstar**2 * nc.ss * teff**4

    mstar = 10.**logg * rstar**2 / nc.gg

    assert (teff >= 0.), 'Negative effective temperature!'
    assert (mstar >= 0.), 'Negative stellar mass!'
    assert (lstar >= 0.), 'Negative stellar luminosity!'
    assert (rstar >= 0.), 'Negative stellar radius!'

    print('-------------------------------------------')
    print('Interpolating in Ames-Dusty model atmospheres')
    print('Stellar parameters: ')
    print('-------------------------------------------')
    print('Teff [K]          : ', teff)
    print('Radius [Rsun]     : ', rstar / nc.rs)
    print('Luminosity [Lsun] : ', lstar / nc.ls)
    print('Mass [Msun]       : ', mstar / nc.ms)
    print('logg              : ', logg)
    print('-------------------------------------------')

    teff_grid = np.append((np.arange(31.) + 9.), (np.arange(30.) * 2. + 40.))
    logg_grid = np.arange(6.) * 0.5 + 3.5

    # Bracket the input teff and logg values
    ii = abs(teff_grid - teff / 100.).argmin()
    idt1 = ii
    if teff_grid[ii] > teff / 100.:
        idt1 = ii - 1
    else:
        idt1 = ii
    idt2 = idt1 + 1

    ii = abs(logg_grid - logg).argmin()
    if logg < logg_grid[0]:
        idg1 = 0
        idg2 = 0
    elif logg > logg_grid[-1]:
        idg1 = logg_grid.shape[0] - 1
        idg2 = logg_grid.shape[0] - 1
    else:
        idg1 = ii
        if logg_grid[ii] > logg:
            idg1 = ii - 1
        else:
            idg1 = ii
        idg2 = idg1 + 1

    print('Bracketing spectrum : ', teff, logg)
    print('Teff  : ', teff_grid[idt1], teff_grid[idt2])
    print('log(g): ', logg_grid[idg1], logg_grid[idg2])

    # Generate the spectral file names
    mph = '0.0'
    fname_11 = 'lte' + ("%d" % teff_grid[idt1]) + '-' + ("%.1f" % logg_grid[idg1]) + '-' + mph + '.AMES-dusty.7.gz'
    fname_12 = 'lte' + ("%d" % teff_grid[idt1]) + '-' + ("%.1f" % logg_grid[idg2]) + '-' + mph + '.AMES-dusty.7.gz'
    fname_22 = 'lte' + ("%d" % teff_grid[idt2]) + '-' + ("%.1f" % logg_grid[idg2]) + '-' + mph + '.AMES-dusty.7.gz'
    fname_21 = 'lte' + ("%d" % teff_grid[idt2]) + '-' + ("%.1f" % logg_grid[idg1]) + '-' + mph + '.AMES-dusty.7.gz'

    # Create a directory
    if os.path.exists("./tmp"):
        shutil.rmtree("./tmp")

    os.system('mkdir tmp')
    os.system('cp -v ' + modeldir + '/' + fname_11 + ' ./tmp')
    os.system('cp -v ' + modeldir + '/' + fname_12 + ' ./tmp')
    os.system('cp -v ' + modeldir + '/' + fname_22 + ' ./tmp')
    os.system('cp -v ' + modeldir + '/' + fname_21 + ' ./tmp')

    # Unzip the files
    os.chdir('./tmp')
    os.system('gunzip ' + fname_11)
    os.system('gunzip ' + fname_12)
    os.system('gunzip ' + fname_22)
    os.system('gunzip ' + fname_21)
    os.chdir('../')

    # Read the spectra
    sp11 = readAmesDustySpectrum(fname='./tmp/' + fname_11[:-3])
    sp12 = readAmesDustySpectrum(fname='./tmp/' + fname_12[:-3])
    sp22 = readAmesDustySpectrum(fname='./tmp/' + fname_22[:-3])
    sp21 = readAmesDustySpectrum(fname='./tmp/' + fname_21[:-3])

    # Do the interpolation
    c11 = (teff_grid[idt2] - teff / 100.) * (logg_grid[idg2] - logg)
    c12 = (teff_grid[idt2] - teff / 100.) * (logg - logg_grid[idg1])
    c22 = (teff / 100. - teff_grid[idt1]) * (logg - logg_grid[idg1])
    c21 = (teff / 100. - teff_grid[idt1]) * (logg_grid[idg2] - logg)
    c00 = 1. / ((teff_grid[idt2] - teff_grid[idt1]) * (logg_grid[idg2] - logg_grid[idg1]))

    lnu = c00 * (c11 * sp11['inu'] + c12 * sp12['inu'] + c22 * sp22['inu'] + c21 * sp21['inu'])
    bnu = c00 * (c11 * sp11['bnu'] + c12 * sp12['bnu'] + c22 * sp22['bnu'] + c21 * sp21['bnu'])

    shutil.rmtree('./tmp')

    #
    # Scale the spectrum to give the same luminosity as required
    #
    nu = nc.cc / sp11['wav'] * 1e4
    lum = (0.5 * abs(nu[1:] - nu[:-1]) * (lnu[1:] + lnu[:-1])).sum()
    lnu *= lstar / lum

    lum = (0.5 * abs(nu[1:] - nu[:-1]) * (bnu[1:] + bnu[:-1])).sum()
    bnu *= lstar / lum

    return {'wav': sp11['wav'], 'lnu': lnu, 'bnu': bnu}

