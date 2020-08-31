"""This module contains classes for radiation sources
"""

from __future__ import absolute_import
from __future__ import print_function
import traceback
import os

try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use the python module of RADMC-3D you need to install Numpy')
    print(traceback.format_exc())


from . import staratm
from . reggrid import *


class radmc3dRadSources(object):
    """Class of the radiation sources.
    Currently discrete stars and continuous starlike source, the latter only in spherical coordinates.


    Attributes
    ----------

    wav          : ndarray
                    Wavelength for the stellar spectrum

    freq         : ndarray
                    Frequency for the stellar spectrum

    nwav         : int
                    Number of wavelenghts in the stellar spectrum

    nfreq        : int
                    Number of frequencies in the stellar spectrum

    mstar        : list
                    List of stellar masses

    tstar        : list
                    List of stellar effective temperatures

    rstar        : list
                    List of stellar radii

    lstar        : list
                    List of stellar luminosities

    nstar        : int
                    Number of stars

    pstar        : list
                    Each element of the list contains a three element list, the cartesian coordinates of the stars

    fnustar      : ndarray
                    Stellar spectrum (flux@1pc in erg/s/cm/cm/Hz)

    csdens       : ndarray
                    Stellar density for continuous starlike source

    csntemplate  : int
                    Number of stellar templates

    cstemp       : ndarray
                    Stellar template

    cstemptype   : int
                    Stellar template type 1 - Blackbody given by the effective temperature
                    2 - Frequency dependent spectrum

    cststar      : ndarray
                    Stellar effective temperature

    csmstar      : ndarray
                    Stellar mass

    csrstar      : ndarray
                    Stellar radius

    tacc         : ndarray
                    Effective temperature of a viscous accretion disk as a function of cylindrical radius

    accrate      : float
                    Accretion rate of the viscous accretion disk [g/s]

    fnuaccdisk   : ndarray
                    Spatially integrated frequency-dependent flux density of the accretion disk @ 1pc distance

    tspot        : float
                    Temperature of the hot spot / boundary layer on the stellar surface

    starsurffrac : float
                    Fraction of the stellar surface covered by the hot spot / boundary layer

    fnustpot     : ndarray
                    Frequency-dependent flux density of the hot spot / boundary layer @ 1pc distance

    """

    def __init__(self, ppar=None, grid=None):

        # Spatial and frequency grid
        self.grid = grid

        # Discrete stars

        self.mstar = []
        self.tstar = []
        self.rstar = []
        self.lstar = []
        self.nstar = 0
        self.pstar = []
        self.fnustar = []
        self.staremis_type = []

        # Continuous starlike source
        self.csdens = []
        self.csntemplate = 0
        self.cstemp = []
        self.cstemptype = 1
        self.cststar = []
        self.csmstar = []
        self.csrstar = []

        # Viscous accretion disk
        self.incl_accretion = False
        self.tacc = []
        self.accrate = 0.
        self.fnuaccdisk = []

        # Hot spot / boundary layer - to model accretion in YSOs
        self.tspot = 0.
        self.starsurffrac = 0.
        self.fnuspot = []

        if ppar:
            if isinstance(ppar['mstar'], list):
                self.mstar = ppar['mstar']
            else:
                self.mstar = [ppar['mstar']]

            if isinstance(ppar['tstar'], list):
                self.tstar = ppar['tstar']
            else:
                self.tstar = [ppar['tstar']]

            if isinstance(ppar['rstar'], list):
                self.rstar = ppar['rstar']
            else:
                self.rstar = [ppar['rstar']]

            for istar in range(self.nstar):
                self.lstar.append(4. * np.pi * self.rstar[istar]**2. * nc.ss * self.tstar[istar]**4.)
            self.pstar = ppar['pstar']

            if 'incl_cont_stellarsrc' in ppar:
                self.incl_accretion = ppar['incl_cont_stellarsrc']
            else:
                self.incl_accretion = False

            if 'accrate' in ppar:
                self.accrate = ppar['accrate']
            else:
                self.accrate = 0.

    def findPeakStarspec(self):

        """Calculates the peak wavelength of the stellar spectrum.

        Returns
        -------

        The peak wavelength of the stellar spectrum in nu*Fnu for all
            stars as a list
        """

        pwav = []

        for istar in range(self.nstar):
            ii = (self.fnustar[:, istar] * self.grid.freq).argmax()
            pwav.append(self.grid.wav[ii])

        return pwav

    def readStarsinp(self, fname=''):
        """Reads the data of discrete stellar sources from the stars.inp file.

        Parameters
        ----------

        fname : str, optional
                File name of the file that should be read (if omitted stars.inp will be used)
        """

        if fname == '':
            fname = 'stars.inp'

        with open(fname, 'r') as rfile:

            dum = rfile.readline()
            iformat = int(dum)
            if iformat != 2:
                raise ValueError(' Unknown file format. Format number ', iformat, ' is unknown.')

            dum = rfile.readline().split()
            self.nstar = int(dum[0])
            # self.grid  = radmc3dGrid()
            self.grid = radmc3dGrid()
            self.grid.readWavelengthGrid()
            self.grid.nwav = int(dum[1])
            self.grid.nfreq = self.grid.nwav
            self.rstar = []
            self.mstar = []
            self.tstar = []
            for istar in range(self.nstar):
                dum = rfile.readline()
                while dum.strip()=='':
                    dum = rfile.readline()
                dum = dum.split()
                self.rstar.append(float(dum[0]))
                self.mstar.append(float(dum[1]))
                self.pstar.append([float(dum[2]), float(dum[3]), float(dum[4])])

            dum = rfile.readline()
            wav = []
            for ilam in range(self.grid.nwav):
                dum = rfile.readline()
                wav.append(float(dum))

            self.grid.wav = np.array(wav, dtype=float)
            self.grid.freq = nc.cc / self.grid.wav * 1e4

            # Empty line
            dum = rfile.readline()
            self.fnustar = np.zeros([self.grid.nfreq, self.nstar], dtype=float)
            for istar in range(self.nstar):
                # First data line containing the stellar spectrum or the effective temperature (if negative)
                dum = rfile.readline()
                # If we have only the stellar temperature
                if float(dum) < 0:
                    self.tstar.append(-float(dum))
                    #
                    # Now calculates the stellar spectrum
                    #
                    self.fnustar[:, istar] = 2. * nc.hh * self.grid.freq**3. / nc.cc**2 / (
                        np.exp(nc.hh * self.grid.freq / nc.kk / self.tstar[istar]) - 1.0) * np.pi \
                        * self.rstar[istar]**2. / nc.pc**2.

                # If we have the full frequency-dependent spectrum
                else:
                    self.tstar.append(0.0)
                    self.fnustar[0, istar] = float(dum)
                    for ifreq in range(1, self.grid.nfreq):
                        self.fnustar[ifreq, istar] = float(rfile.readline())

    def writeStarsinp(self, ppar=None, wav=None, freq=None, old=False):
        """Writes the input file for discrete stellar sources (stars.inp).

        Parameters
        ----------

        ppar  : dictionary
                Dictionary containing all parameters of the model (only mandatory if accretion is switched on)

        wav   : ndarray, optional
                Wavelength grid for the stellar spectrum

        freq  : ndarray, optional
                Frequency grid for the stellar spectrum (either freq or wav should be set)

        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """
        if not old:
            if freq is not None:
                if wav is not None:
                    raise ValueError('Either wav or freq should be set but not both.')
                else:
                    self.grid.wav = nc.cc / np.array(freq) * 1e4
                    self.grid.freq = np.array(freq)
                    self.grid.nwav = self.grid.wav.shape[0]
                    self.grid.nfreq = self.grid.nwav

            elif wav is not None:
                if freq is not None:
                    raise ValueError('Either wav or freq should be set but not both.')
                else:
                    self.grid.wav = np.array(wav)
                    self.grid.freq = nc.cc / self.grid.wav * 1e4
                    self.grid.nwav = self.grid.wav.shape[0]
                    self.grid.nfreq = self.grid.nwav

            self.nstar = len(self.rstar)
            # self.pstar = ppar['pstar']

            print('Writing stars.inp')
            with open('stars.inp', 'w') as wfile:
                wfile.write('%d\n' % 2)
                wfile.write('%d %d\n' % (self.nstar, self.grid.nwav))

                if self.nstar > 1:
                    for istar in range(self.nstar):
                        wfile.write('%.9e %.9e %.9e %.9e %.9e\n' % (self.rstar[istar], self.mstar[istar],
                                                                    self.pstar[istar][0], self.pstar[istar][1],
                                                                    self.pstar[istar][2]))
                else:
                    wfile.write('%.9e %.9e %.9e %.9e %.9e\n' % (self.rstar[0], self.mstar[0],
                                                                self.pstar[0], self.pstar[1], self.pstar[2]))

                wfile.write('%s\n' % ' ')
                for ilam in range(self.grid.nwav):
                    wfile.write('%.9e\n' % self.grid.wav[ilam])
                wfile.write('%s\n' % ' ')

                # If accretion is active write the sum of the spot and the star emission
                # NOTE, for now the spot emission is only added to the first star in the list
                if self.incl_accretion:
                    # Get the stellar spectrum
                    self.getStarSpectrum(tstar=self.tstar, rstar=self.rstar)
                    # Get the spot spectrum
                    self.getSpotSpectrum(ppar=ppar)

                    # Write out the spectrum of the first star with the additional spot luminosity
                    for ilam in range(self.grid.nwav):
                        wfile.write('%.9e\n' % (self.fnustar[ilam, 0] + self.fnuspot[ilam]))

                    # Now write the spectrum of all other discrete stars without the spot emission
                    for istar in range(1, self.nstar):
                        for ilam in range(self.grid.nwav):
                            wfile.write('%.9e\n' % (self.fnustar[ilam, istar]))

                else:
                    self.getStarSpectrum(ppar=ppar)
                    for istar in range(self.nstar):
                        if self.staremis_type[istar].strip().lower() == "blackbody":
                            wfile.write('%.9e\n' % (-self.tstar[istar]))
                        else:
                            for ilam in range(self.grid.nwav):
                                wfile.write('%.9e\n' % (self.fnustar[ilam, istar]))
        else:
            if freq is not None:
                self.grid.wav = nc.cc / np.array(freq) * 1e4
                self.grid.freq = np.array(freq)
                self.grid.nwav = self.grid.wav.shape[0]
                self.grid.nfreq = self.grid.nwav

            if wav is not None:
                self.grid.wav = np.array(wav)
                self.grid.freq = nc.cc / self.grid.wav * 1e4
                self.grid.nwav = self.grid.wav.shape[0]
                self.grid.nfreq = self.grid.nwav

            else:
                raise ValueError('Either wav or freq should be set but not both.')

            self.nstar = len(self.rstar)
            self.getStarSpectrum(ppar=ppar)

            if self.grid.freq[-1] < self.grid.freq[0]:
                self.grid.freq = self.grid.freq[::-1]
                self.grid.wav = self.grid.wav[::-1]
                for istar in range(self.nstar):
                    self.fnustar[:, istar] = self.fnustar[::-1, istar]

            print('Writing starinfo.inp')
            fname = 'starinfo.inp'
            with open(fname, 'w') as wfile:
                wfile.write("1\n")
                wfile.write("%.7e\n" % ppar['rstar'][0])
                wfile.write("%.7e\n" % ppar['mstar'][0])
                wfile.write("%.7e\n" % ppar['tstar'][0])

            print('Writing starspectrum.inp')
            fname = 'starspectrum.inp'
            with open(fname, 'w') as wfile:
                wfile.write("%d\n" % self.grid.nfreq)
                wfile.write(" \n")
                for i in range(self.grid.nfreq):
                    wfile.write("%.7e %.7e\n" % (self.grid.freq[i], self.fnustar[i, 0]))

    def getStarSpectrum(self, tstar=None, rstar=None, lstar=None, mstar=None, ppar=None, grid=None):
        """Calculates a blackbody stellar spectrum.

        Parameters
        ----------

        tstar : list
                Effective temperature of the stars in [K]

        rstar : list
                Radius of the stars in [cm]

        lstar : list
                Bolometric luminosity of the star [erg/s] (either rstar or lstar should be given)

        mstar : list
                Stellar mass in [g] (only required if an atmosphere model is used to calculate logg)

        ppar  : dictionary
                Dictionary containing all input parameters

        grid  : radmc3dGrid, optional
                An instance of a radmc3dGrid class containing the spatial and wavelength grid
        """
        #
        # Check the input which parameters are set and which should be calculated
        #

        if grid is not None:
            self.grid = grid

        if ppar is not None:
            if tstar is None:
                tstar = ppar['tstar']
            if rstar is None:
                rstar = ppar['rstar']
            if mstar is None:
                self.mstar = ppar['mstar']

        if mstar:
            self.mstar = mstar

        if tstar:
            if not isinstance(tstar, list):
                tstar = [tstar]
            dum1 = len(tstar)
            if lstar and rstar:
                raise ValueError(' Only two of the input variables tstar, rstar, lstar should be set not all three')
            elif lstar:
                if len(lstar) != dum1:
                    raise ValueError('lstar and tstar have different number of elements')
                else:
                    self.tstar = np.array(tstar)
                    self.lstar = np.array(lstar)
                    self.nstar = self.lstar.shape[0]
                    self.rstar = np.sqrt(self.lstar / (4. * np.pi * nc.ss * self.tstar**4.))
        else:
            if lstar and rstar:
                if len(lstar) != len(rstar):
                    raise ValueError('lstar and rstar have different number of elements')
                else:
                    self.lstar = np.array(lstar)
                    self.rstar = np.array(rstar)
                    self.nstar = self.rstar.shape[0]
                    self.tstar = (self.lstar / (4. * np.pi * nc.ss * self.rstar**2.))**0.25

        #
        # If we take blackbody spectrum
        #

        if ppar:
            if 'staremis_type' in ppar:
                self.staremis_type = ppar['staremis_type']
            else:
                self.staremis_type = []
                for i in range(self.nstar):
                    self.staremis_type.append("blackbody")
        else:
            self.staremis_type = []
            for i in range(self.nstar):
                self.staremis_type.append("blackbody")

        self.fnustar = np.zeros([self.grid.nwav, len(self.tstar)], dtype=np.float64)
        for istar in range(self.nstar):
            if self.staremis_type[istar].strip().lower() == "blackbody":
                self.fnustar[:, istar] = 2. * nc.hh * self.grid.freq**3. / nc.cc**2 / (
                                         np.exp(nc.hh * self.grid.freq / nc.kk / self.tstar[istar]) - 1.0) \
                                         * np.pi * self.rstar[istar]**2. / nc.pc**2.
            elif self.staremis_type[istar].strip().lower() == "kurucz":
                dum = staratm.getAtmModel(teff=self.tstar[istar], mstar=self.mstar[istar], rstar=self.rstar[istar],
                                          iwav=self.grid.wav, model="kurucz")
                self.fnustar[:, istar] = dum['lnu'] / (4. * np.pi * nc.pc**2)

            elif self.staremis_type[istar].strip().lower() == "nextgen":
                dum = staratm.getAtmModel(teff=self.tstar[istar], mstar=self.mstar[istar], rstar=self.rstar[istar],
                                          iwav=self.grid.wav, model="nextgen")
                self.fnustar[:, istar] = dum['lnu'] / (4. * np.pi * nc.pc**2)

            elif self.staremis_type[istar].strip().lower() == "ames-dusty":
                dum = staratm.getAtmModel(teff=self.tstar[istar], mstar=self.mstar[istar], rstar=self.rstar[istar],
                                          iwav=self.grid.wav, model="ames-dusty")
                self.fnustar[:, istar] = dum['lnu'] / (4. * np.pi * nc.pc**2)

            else:
                raise ValueError('Unknown stellar atmosphere model : ', self.staremis_type[istar],
                                 ' Only "kurucz" or "nextgen" are supported.')

    def getAccdiskTemperature(self, ppar=None, grid=None):
        """Calculates the effective temperature of a viscous accretion disk.

        Parameters
        ----------

        ppar : dictionary
               Dictionary containing all input parameters keys should include
               * mstar   : stellar mass
               * rstar   : stellar radius
               * accrate : accretion rate

               NOTE, that for the calculation of the effective disk temperature only the first
               star is used if more than one values are given in mstar and rstar.

        grid : radmc3dGrid, optional
               An instance of a radmc3dGrid class containing the spatial and wavelength grid

        """

        if grid is not None:
            self.grid = grid

        self.tacc = ((3.0 * nc.gg * ppar['mstar'][0] * ppar['accrate']) / (8.0 * np.pi * self.grid.x**3 * nc.ss)
                     * (1.0 - (ppar['rstar'][0] / self.grid.x)**0.5))**0.25

    # --------------------------------------------------------------------------------------------------
    def getSpotSpectrum(self, ppar=None, grid=None):
        """Calculates the spectrum of a hot spot / boundary layer on the stellar surface

        Parameters
        ----------

        ppar : dictionary
               Dictionary containing all input parameters keys should include
               * mstar   : stellar mass
               * rstar   : stellar radius
               * accrate : accretion rate

               NOTE, that for the calculation of the effective disk temperature only the first
               star is used if more than one values are given in mstar and rstar.

        grid : radmc3dGrid, optional
               An instance of a radmc3dGrid class containing the spatial and wavelength grid
        """
        #
        # Check the input which parameters are set and which should be calculated
        #

        if grid is not None:
            self.grid = grid

        # Calculate the total spot luminosity assuming boundary layer accretion, i.e.
        # that half of the total accretion luminosity is released from the boundary layer

        if 'accrate' in ppar:
            tot_acclum = 0.5 * nc.gg * ppar['mstar'][0] * ppar['accrate'] / ppar['rstar'][0]
            spotsurf = 4. * np.pi * ppar['rstar'][0]**2 * ppar['starsurffrac']
            self.starsurffrac = ppar['starsurffrac']
            if spotsurf == 0.:
                self.tspot = 0.

                # Now calculate the spot spectrum (i.e. the flux density @ 1pc distance)
                self.fnuspot = np.zeros(self.grid.nfreq, dtype=float)
            else:
                self.tspot = (tot_acclum / spotsurf / nc.ss)**0.25

                # Now calculate the spot spectrum (i.e. the flux density @ 1pc distance)
                self.fnuspot = np.pi * ppar['rstar'][0]**2 * ppar['starsurffrac'] / nc.pc**2 \
                               * 2. * nc.hh * self.grid.freq**3 / nc.cc**2 \
                               / (np.exp(nc.hh * self.grid.freq / nc.kk / self.tspot) - 1.0)
        else:
            self.fnuspot = np.zeros(self.grid.nfreq, dtype=float)

    def getTotalLuminosities(self, readInput=True):
        """Calcultes the frequency integrated luminosities of all radiation sources.


        Parameters
        ----------

        readInput : bool, optional
                    If true the input files of the radiation sources are read and the the total luminosities
                    are calculated from them. If readInput is set to False, the luminosities are calculated
                    by semi-analytic spectra.

        Returns
        -------

        Returns a dictionary with the following keys
            * lnu_star    : Luminosity of the discrete stars
            * lnu_accdisk : Luminosity of the accretion disk
            * lnu_spot    : Luminosity of the hot spot / boundary layer on the stellar surface
        """

        res = {'lnu_star': np.zeros(self.nstar, dtype=float),
               'lnu_accdisk': 0.,
               'lnu_spot': 0.}

        # If readIpnut
        if readInput:
            self.readStarsinp()

            # Note the negative sign in dnu is there because the frequency array is ordered in wavelength not in
            # frequency
            dnu = -(self.grid.freq[1:] - self.grid.freq[:-1])

            for istar in range(self.nstar):
                res['lnu_star'][istar] = 0.5 * (
                    (self.fnustar[1:, istar] + self.fnustar[:-1, istar]) * dnu).sum() * 4. * np.pi * nc.pc**2

            # Calculate the luminosity in the continuous stellar sources
            csDensFound = False
            csTempFound = False
            if os.path.exists('stellarsrc_density.inp'):
                self.readStellarsrcDensity(fname='stellarsrc_density.inp', binary=False)
                csDensFound = True
            if os.path.exists('stellarsrc_density.binp'):
                self.readStellarsrcDensity(fname='stellarsrc_density.binp', binary=True)
                csDensFound = True
            if os.path.exists('stellarsrc_templates.inp'):
                self.readStellarsrcTemplates()
                csTempFound = True
            if csDensFound & csTempFound:
                vol = self.grid.getCellVolume()

                dnu = abs(self.grid.freq[1:] - self.grid.freq[:-1])
                lum = 0.
                for itemp in range(self.csntemplate):
                    for ix in range(self.grid.nx):
                        for iy in range(self.grid.ny):
                            for iz in range(self.grid.nz):
                                if self.csdens[ix, iy, iz, itemp] > 0.:
                                    expterm = (nc.hh * self.grid.freq / nc.kk / (-self.cststar[ix])).clip(-600, 600)
                                    bb = 2. * nc.hh * self.grid.freq**3 / nc.cc**2 / (np.exp(expterm) - 1.)
                                    dum = bb * np.pi * vol[ix, iy, iz] * 4. * np.pi * self.csdens[ix, iy, iz, itemp]
                                    lum = lum + ((dum[1:] + dum[:-1]) * 0.5 * dnu).sum()

                res['lnu_accdisk'] = lum

            else:
                res['lnu_spot'] = 0.
                res['lnu_accdisk'] = 0.
        else:

            for istar in range(self.nstar):
                res['lnu_star'][istar] = 4. * np.pi * self.rstar[istar]**2 * nc.ss * self.tstar[istar]**4

            if self.accrate == 0.:
                print('Viscsous accretion is switched off')
                res['lnu_spot'] = 0.
                res['lnu_accdisk'] = 0.
            if not self.incl_accretion:
                print('Viscsous accretion is switched off')
                res['lnu_spot'] = 0.
                res['lnu_accdisk'] = 0.
            else:
                # Calculate the spot luminosity
                res['lnu_spot'] = 0.5 * nc.gg * self.mstar[0] * self.accrate / self.rstar[0]
                # Calculate the accretion disk luminosity
                res['lnu_accdisk'] = 0.5 * nc.gg * self.mstar[0] * self.accrate / self.rstar[0]

        return res

    def getAccdiskSpectra(self, ppar=None, grid=None, incl=0.):
        """Calculates the emergent spectra of an optically thick accretion disk at face-on orientation (incl=0deg).

        Parameters
        ----------

        ppar : dictionary
               Dictionary containing all input parameters keys should include
               * mstar   : stellar mass
               * rstar   : stellar radius
               * accrate : accretion rate

               NOTE, that for the calculation of the effective disk temperature only the first
               star is used if more than one values are given in mstar and rstar.

        incl : float, optional
               Inclination angle in degrees at which the spectrum be calculated  (default - 0deg)

        grid : radmc3dGrid, optional
               An instance of a radmc3dGrid class containing the spatial and wavelength grid

        """

        if 'accrate' in ppar:
            if ppar['accrate'] > 0.:
                self.getAccdiskTemperature(ppar=ppar, grid=grid)
                fnuaccdisk = np.zeros([self.grid.nx, self.grid.nwav], dtype=float)
                for ix in range(self.grid.nx):
                    dum = nc.hh * self.grid.freq / nc.kk / self.tacc[ix]
                    dum = dum.clip(-600., 600.)
                    bb = 2. * nc.hh * self.grid.freq**3 / nc.cc**2 / (np.exp(np.float64(dum)) - 1.0)
                    fnuaccdisk[ix, :] = bb * np.pi * (self.grid.xi[ix + 1]**2 - self.grid.xi[ix]**2) / nc.pc**2
                self.fnuaccdisk = fnuaccdisk.sum(0) * np.cos(incl/180.*np.pi)
            else:
                self.fnuaccdisk = np.zeros([self.grid.nwav], dtype=float)
        else:
            self.fnuaccdisk = np.zeros([self.grid.nwav], dtype=float)

    def getAccdiskStellarTemplates(self, ppar=None, grid=None):
        """Calculates the stellar template for continuous starlike sources for modeling a viscous accretion disk.


        Parameters
        ----------

        ppar : dictionary
               Dictionary containing all input parameters keys should include:
               * mstar   : stellar mass
               * rstar   : stellar radius
               * accrate : accretion rate

               NOTE, that for the calculation of the effective disk temperature only the first
               star is used if more than one values are given in mstar and rstar.

        grid : radmc3dGrid, optional
               An instance of a radmc3dGrid class containing the spatial and wavelength grid

        """

        if self.incl_accretion:
            self.getAccdiskTemperature(ppar=ppar, grid=grid)
            self.cstemptype = 1
            self.cststar = -self.tacc
            self.csmstar = self.cststar * 0. + 1.
            self.csrstar = self.cststar * 0. + 1.
            self.csntemplate = self.grid.nx
        else:

            self.cstemptype = 1
            self.cststar = np.zeros(self.grid.nx)
            self.csmstar = self.cststar * 0.
            self.csrstar = self.cststar * 0.
            self.csntemplate = self.grid.nx

    def getAccdiskStellarDensity(self, grid=None):
        """Calculates the stellar density for continuous starlike sources for modeling a viscous accretion disk.

        Parameters
        ----------
        grid : radmc3dGrid, optional
               An instance of a radmc3dGrid class containing the spatial and wavelength grid

        """

        if grid is not None:
            self.grid = grid

        self.csdens = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz, self.grid.nx], dtype=float)
        vol = self.grid.getCellVolume()

        if self.incl_accretion:
            if self.grid.crd_sys != 'sph':
                raise RuntimeError('Viscous accretion is currently available only in spherical coordinate system')
            else:
                if abs(self.grid.yi[self.grid.ny] - np.pi) < 1e-8:

                    for ix in range(self.grid.nx):
                        dA = 2.0 * (self.grid.xi[ix + 1]**2 - self.grid.xi[ix]**2) * np.pi * (
                            self.grid.zi[1:] - self.grid.zi[:-1]) / (2. * np.pi)
                        dV = vol[ix, self.grid.ny / 2 - 1, :] + vol[ix, self.grid.ny / 2, :]
                        self.csdens[ix, self.grid.ny / 2 - 1, :, ix] = dA / dV / (4. * np.pi)
                        self.csdens[ix, self.grid.ny / 2, :, ix] = dA / dV / (4. * np.pi)

                elif abs(self.grid.yi[self.grid.ny] - np.pi / 2.) < 1e-8:
                    for ix in range(self.grid.nx):
                        dA = 2.0 * (self.grid.xi[ix + 1]**2 - self.grid.xi[ix]**2) * np.pi * (
                            self.grid.zi[1:] - self.grid.zi[:-1]) / (2. * np.pi)
                        dV = vol[ix, self.grid.ny - 1, :] * 2.
                        self.csdens[ix, self.grid.ny - 1, :, ix] = dA / dV / (4. * np.pi)

        return True

    def readStellarsrcTemplates(self, fname='stellarsrc_templates.inp'):
        """Reads the stellar template of a continuous starlike source.

        Parameters
        ----------

        fname : str, optional
                Name of the file from which the stellar templates will be read. If omitted the default
                'stellarsrc_templates.inp' will be used.
        """

        with open(fname, 'r') as rfile:

            self.grid = radmc3dGrid()
            self.grid.readGrid()

            hdr = np.fromfile(fname, count=3, sep="\n", dtype=int)

            if self.grid.nwav != hdr[2]:
                raise ValueError('Number of grid points in wavelength_micron.inp is not equal to that in ' + fname)
            else:
                self.csntemplate = hdr[1]
                dum = ''

                # Read the header
                for i in range(3):
                    dum = rfile.readline()

                # Read the frequency grid
                for i in range(hdr[2]):
                    dummy = float(rfile.readline())
                    if np.abs((nc.cc / dummy * 1e4 - self.grid.wav[i]) / self.grid.wav[i]) > 1e-4:
                        raise ValueError('The wavelength grid in wavelength_micron.inp is different '
                                         + 'from that in ' + fname)

                dum = rfile.readline()
                if float(dum) > 0.0:
                    self.cstemp = np.zeros([self.csntemplate, self.grid.nwav], dtype=float)
                    self.cststar = []
                    self.csrstar = []
                    self.csmstar = []

                    self.cstemp[0, 0] = float(dum)
                    for inu in range(1, self.grid.nwav):
                        dum = rfile.readline()
                        self.cstemp[0, inu] = float(dum)

                    for itemp in range(1, self.csntemplate):
                        for inu in range(self.grid.nwav):
                            dum = rfile.readline()
                            self.cstemp[itemp, inu] = float(dum)

                else:
                    self.cstemp = []
                    self.cststar = np.zeros(self.csntemplate, dtype=float)
                    self.csrstar = np.zeros(self.csntemplate, dtype=float)
                    self.csmstar = np.zeros(self.csntemplate, dtype=float)

                    self.cststar[0] = float(dum)
                    dum = rfile.readline()
                    self.csrstar[0] = float(dum)
                    dum = rfile.readline()
                    self.csmstar[0] = float(dum)

                    for i in range(1, self.csntemplate):
                        dum = rfile.readline()
                        self.cststar[i] = float(dum)
                        dum = rfile.readline()
                        self.csrstar[i] = float(dum)
                        dum = rfile.readline()
                        self.csmstar[i] = float(dum)

    def writeStellarsrcTemplates(self, fname='stellarsrc_templates.inp'):
        """Writes the stellar template of a continuous starlike source.

        Parameters
        ----------

        fname : str, optional
                Name of the file into which the stellar templates will be written. If omitted the default
                'stellarsrc_templates.inp' will be used.
        """
        # First check if we'd need to write anything at al

        if len(self.cststar) == 0:
            if len(self.cstemp) == 0:
                if os.path.exists('stellarsrc_templates.inp'):
                    print('The continuous starlike source seems to be inactive (zero input luminosity) '
                          + 'still stellarsrc_templates.inp file is present in the current working directory.')
                    dum = input('Can it be deleted (yes/no)')
                    if dum.strip().lower()[0] == 'y':
                        os.system('rm stellarsrc_templates.inp')
                    return
                return
            else:
                if np.abs(self.cstemp).max() == 0.:
                    if os.path.exists('stellarsrc_templates.inp'):
                        print('The continuous starlike source seems to be inactive (zero input luminosity)'
                              + ' still stellarsrc_templates.inp file is present in the current working directory.')
                        dum = input('Can it be deleted (yes/no)')
                        if dum.strip().lower()[0] == 'y':
                            os.system('rm stellarsrc_templates.inp')
                        return
                    return

        else:
            if np.abs(self.cststar).max() == 0.:
                if os.path.exists('stellarsrc_templates.inp'):
                    print('The continuous starlike source seems to be inactive (zero input luminosity) '
                          + ' still stellarsrc_templates.inp file is present in the current working directory.')
                    dum = input('Can it be deleted (yes/no)')
                    if dum.strip().lower()[0] == 'y':
                        os.system('rm stellarsrc_templates.inp')
                    return
                return

        print('Writing ' + fname)
        with open(fname, 'w') as wfile:
            # Format number
            wfile.write("%d\n" % 1)
            # Nr of templates
            wfile.write("%d\n" % self.csntemplate)
            # Nr of wavelengths
            wfile.write("%d\n" % self.grid.nwav)
            # Write the wavelength grid (in micron!)
            for ilam in range(self.grid.nwav):
                wfile.write("%.9e\n" % self.grid.freq[ilam])

            # Now write the templates
            # Similar to the discrete stellar imput if the first number is negative it means that
            #  instead of a full frequency-dependent spectrum only the blackbody temperature is given.
            #  Thus I'd only need to give the temperatures as negative numbers and radmc-3d will take care
            #  of calculating the Planck-function. This will save some harddisk space.

            if self.cstemptype == 1:
                for itemp in range(self.csntemplate):
                    # Effective temperature
                    if self.cststar[itemp] > 0:
                        wfile.write("%.9e\n" % (-self.cststar[itemp]))
                    else:
                        wfile.write("%.9e\n" % self.cststar[itemp])
                    # "Stellar radius"
                    wfile.write("%.9e\n" % self.csrstar[itemp])
                    # "Stellar mass"
                    wfile.write("%.9e\n" % self.csmstar[itemp])
            elif self.cstemptype == 2:
                for itemp in range(self.csntemplate):
                    for inu in range(self.grid.nwav):
                        wfile.write("%.9e\n" % self.cstemp[itemp, inu])
            else:
                raise ValueError('Unknown cstemptype for the continuous starlike source')

    def readStellarsrcDensity(self, fname=None, binary=False):
        """Reads the stellar density of a continuous starlike source.

        Parameters
        ----------

        fname  : str, optional
                 Name of the file from which the stellar templates will be read. If omitted the default
                 'stellarsrc_templates.inp' will be used.

        binary : bool, optional
                 If True the file should contain a C-style binary stream, if False the file should be
                 written as formatted ASCII
        """

        if fname is None:
            if binary:
                fname = 'stellarsrc_density.binp'
            else:
                fname = 'stellarsrc_density.inp'

        self.grid = radmc3dGrid()
        self.grid.readGrid()
        self.csdens = None

        if binary:
            hdr = np.fromfile(fname, count=4, dtype=int)

            if hdr[2] != (self.grid.nx * self.grid.ny * self.grid.nz):
                raise ValueError('Number of grid points in ' + fname + ' is different from that in amr_grid.inp\n',
                                 (self.grid.nx * self.grid.ny * self.grid.nz), hdr[2])

            if hdr[1] == 8:
                data = np.fromfile(fname, count=-1, dtype=np.float64)
            elif hdr[1] == 4:
                data = np.fromfile(fname, count=-1, dtype=float)
            else:
                raise TypeError('Unknown datatype in ' + fname + ' radmc3dPy supports only 4 byte float or '
                                + ' 8 byte double types.')

            data = np.reshape(data[4:], [hdr[3], self.grid.nz, self.grid.ny, self.grid.nx])
            data = np.swapaxes(data, 0, 3)
            data = np.swapaxes(data, 1, 2)
        else:
            hdr = np.fromfile(fname, count=3, sep="\n", dtype=int)

            if (self.grid.nx * self.grid.ny * self.grid.nz) != hdr[1]:
                raise ValueError('Number of grid points in amr_grid.inp is not equal to that in ' + fname,
                                 (self.grid.nx * self.grid.ny * self.grid.nz), hdr[1])
            else:

                data = np.fromfile(fname, count=-1, sep="\n", dtype=np.float64)
                data = np.reshape(data[3:], [hdr[2], self.grid.nz, self.grid.ny, self.grid.nx])
                # We need to change the axis orders as Numpy always reads  in C-order while RADMC-3D
                # uses Fortran-order
                data = np.swapaxes(data, 0, 3)
                data = np.swapaxes(data, 1, 2)

                self.csdens = data

    def writeStellarsrcDensity(self, fname='', binary=False):
        """Writes the stellar density of a continuous starlike source.

        Parameters
        ----------

        fname : str, optional
                Name of the file into which the stellar templates will be written. If omitted the default
                'stellarsrc_templates.inp' will be used.

        binary : bool, optional
                If True the output will be written in a C-style binary stream, if False the output will be
                formatted ASCII
        """

        # First check if we'd need to write anything at al
        # Both the stellar temperature and the stellar template spectrum arrays are empty
        # So the continuous starlike source needs to be deactivated
        # Let's check if there are any files, and if so ask them to be deleted
        if len(self.cststar) == 0:
            if len(self.cstemp) == 0:
                if (os.path.exists('stellarsrc_density.inp')) | (os.path.exists('stellarsrc_density.binp')):
                    print('The continuous starlike source seems to be inactive (zero input luminosity still '
                          + 'stellarsrc_density.inp/stellarsrc_density.binp file is present in the current '
                          + 'working directory.')
                    dum = input('Can it be deleted (yes/no)')
                    if dum.strip().lower()[0] == 'y':
                        os.system('rm stellarsrc_density.inp')
                        os.system('rm stellarsrc_density.binp')
                    return
                return
            else:
                if np.abs(self.cstemp).max() == 0.:
                    if (os.path.exists('stellarsrc_density.inp')) | (os.path.exists('stellarsrc_density.binp')):
                        print('The continuous starlike source seems to be inactive (zero input luminosity)'
                              + ' still stellarsrc_density.inp/stellarsrc_density.binp file is present in the current '
                              + ' working directory.')
                        dum = input('Can it be deleted (yes/no)')
                        if dum.strip().lower()[0] == 'y':
                            os.system('rm stellarsrc_density.inp')
                            os.system('rm stellarsrc_density.binp')
                        return
                    return

        else:
            if np.abs(self.cststar).max() == 0.:
                if (os.path.exists('stellarsrc_density.inp')) | (os.path.exists('stellarsrc_density.binp')):
                    print('The continuous starlike source seems to be inactive (zero input luminosity)'
                          + ' still stellarsrc_density.inp/stellarsrc_density.binp file is present in the current '
                          + ' working directory.')
                    dum = input('Can it be deleted (yes/no)')
                    if dum.strip().lower()[0] == 'y':
                        os.system('rm stellarsrc_density.inp')
                        os.system('rm stellarsrc_density.binp')
                    return
                return

        if binary:
            if fname.strip() == '':
                fname = 'stellarsrc_density.binp'
            print('Writing ' + fname)
            wfile = open(fname, 'w')
            hdr = np.array([1, 8, self.grid.nx * self.grid.ny * self.grid.nz, self.csntemplate], dtype=int)
            hdr.tofile(wfile)
            data = np.array(self.csdens)
            data = np.swapaxes(data, 0, 3)
            data = np.swapaxes(data, 1, 2)
            data.tofile(wfile)

        else:
            if fname.strip() == '':
                fname = 'stellarsrc_density.inp'

            with open(fname, 'w') as wfile:
                hdr = np.array([1, self.grid.nx * self.grid.ny * self.grid.nz, self.csntemplate], dtype=int)
                hdr.tofile(wfile, sep=" ", format="%d\n")
                data = np.array(self.csdens)
                data = np.swapaxes(data, 0, 3)
                data = np.swapaxes(data, 1, 2)
                data.tofile(wfile, sep=" ", format="%.9e\n")
