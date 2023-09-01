"""This module contains classes/functions to create and read images with RADMC-3D and to calculate
interferometric visibilities and write fits files
For help on the syntax or functionality of each function see the help of the individual functions
"""
from __future__ import absolute_import
from __future__ import print_function
import traceback
import copy
import subprocess as sp
import os

try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use the python module of RADMC-3D you need to install Numpy')
    print(traceback.format_exc())

try:
    import scipy.special as spc
except ImportError:
    spc = None
    print('scipy.special cannot be imported ')
    print('This module is required to be able to calculate Airy-PSFs. Now PSF calculation is limited to Gaussian.')
    print(traceback.format_exc())

try:
    from astropy.io import fits as pf
except ImportError:
    print('astropy.io.fits cannot be imported trying pyfits')
    try:
        import pyfits as pf
    except ImportError:
        pf = None
        print('pyfits cannot be imported. Either of these modules is needed to write RADMC-3D images '
              + 'to FITS format. The rest of radmc3dPy can be used but fits output is now disabled.')
        print(traceback.format_exc())

try:
    import matplotlib.pylab as plt
except ImportError:
    plt = None
    print('Warning')
    print('matplotlib.pyplot cannot be imported')
    print('Without matplotlib you can use the python module to set up a model but you will not be able to plot things')
    print('or display images')

from . import natconst as nc

class radmc3dImage(object):
    """
    RADMC-3D image class

    Attributes
    ----------

    image       : ndarray
                  The image as calculated by radmc3d (the values are intensities in erg/s/cm^2/Hz/ster)

    imageJyppix : ndarray
                  The image with pixel units of Jy/pixel

    x           : ndarray
                  x coordinate of the image [cm]

    y           : ndarray
                  y coordinate of the image [cm]

    nx          : int
                  Number of pixels in the horizontal direction

    ny          : int
                  Number of pixels in the vertical direction

    sizepix_x   : float
                  Pixel size in the horizontal direction [cm]

    sizepix_y   : float
                  Pixel size in the vertical direction [cm]

    nfreq       : int
                  Number of frequencies in the image cube

    freq        : ndarray
                  Frequency grid in the image cube

    nwav        : int
                  Number of wavelengths in the image cube (same as nfreq)

    wav         : ndarray
                  Wavelength grid in the image cube

    """

    def __init__(self):
        self.image = None
        self.imageJyppix = None
        self.x = None
        self.y = None
        self.nx = 0
        self.ny = 0
        self.sizepix_x = 0
        self.sizepix_y = 0
        self.nfreq = 0
        self.freq = None
        self.nwav = 0
        self.wav = None
        self.stokes = False
        self.psf = {}
        self.fwhm = []
        self.pa = 0
        self.dpc = 0
        self.filename = 'image.out'

    def getClosurePhase(self, bl=None, pa=None, dpc=None):
        """Calculates clusure phases for a given model image for any arbitrary baseline triplet.

        Parameters
        ----------

        bl  : list/ndarray
              A list or ndrray containing the length of projected baselines in meter.

        pa  : list/ndarray
              A list or Numpy array containing the position angles of projected baselines in degree.

        dpc : distance of the source in parsec


        Returns
        -------
        Returns a dictionary with the following keys:

            * bl     : projected baseline in meter
            * pa     : position angle of the projected baseline in degree
            * nbl    : number of baselines
            * u      : spatial frequency along the x axis of the image
            * v      : spatial frequency along the v axis of the image
            * vis    : complex visibility at points (u,v)
            * amp    : correlation amplitude
            * phase  : Fourier phase
            * cp     : closure phase
            * wav    : wavelength
            * nwav   : number of wavelengths

        Notes
        -----

        bl and pa should either be an array with dimension [N,3] or if they are lists each element of
        the list should be a list of length 3, since closure phases are calculated only for closed triangles
        """

        ntri = len(bl)
        res = {'bl': np.array(bl, dtype=np.float64),
               'pa': np.array(pa, dtype=np.float64),
               'ntri': ntri,
               'nbl': 3,
               'nwav': self.nwav,
               'wav': self.wav,
               'u': np.zeros([ntri, 3, self.nwav], dtype=np.float64),
               'v': np.zeros([ntri, 3, self.nwav], dtype=np.float64),
               'vis': np.zeros([ntri, 3, self.nwav], dtype=np.complex64),
               'amp': np.zeros([ntri, 3, self.nwav], dtype=np.float64),
               'phase': np.zeros([ntri, 3, self.nwav], dtype=np.float64),
               'cp': np.zeros([ntri, self.nwav], dtype=np.float64)}

        # l = self.x / nc.au / dpc / 3600. / 180. * np.pi
        # m = self.y / nc.au / dpc / 3600. / 180. * np.pi
        # dl = l[1] - l[0]
        # dm = m[1] - m[0]

        for itri in range(ntri):
            print('Calculating baseline triangle # : ', itri)

            dum = self.getVisibility(bl=res['bl'][itri, :], pa=res['pa'][itri, :], dpc=dpc)
            res['u'][itri, :, :] = dum['u']
            res['v'][itri, :, :] = dum['v']
            res['vis'][itri, :, :] = dum['vis']
            res['amp'][itri, :, :] = dum['amp']
            res['phase'][itri, :, :] = dum['phase']
            res['cp'][itri, :] = (dum['phase'].sum(0) / np.pi * 180.) % 360.
            ii = res['cp'][itri, :] > 180.
            if (res['cp'][itri, ii]).shape[0] > 0:
                res['cp'][itri, ii] = res['cp'][itri, ii] - 360.

        return res

    # --------------------------------------------------------------------------------------------------
    def getVisibility(self, bl=None, pa=None, dpc=None):
        """Calculates visibilities for a given set of projected baselines and position angles
        with Discrete Fourier Transform.

        Parameters
        ----------

        bl  : list/ndarray
              A list or ndrray containing the length of projected baselines in meter.

        pa  : list/ndarray
              A list or Numpy array containing the position angles of projected baselines in degree.

        dpc : distance of the source in parsec

        Returns
        -------
        Returns a dictionary with the following keys:

            * bl     : projected baseline in meter
            * pa     : position angle of the projected baseline in degree
            * nbl    : number of baselines
            * u      : spatial frequency along the x axis of the image
            * v      : spatial frequency along the v axis of the image
            * vis    : complex visibility at points (u,v)
            * amp    : correlation amplitude
            * phase  : Fourier phase
            * wav    : wavelength
            * nwav   : number of wavelengths
        """

        nbl = len(bl)
        res = {'bl': np.array(bl, dtype=np.float64),
               'pa': np.array(pa, dtype=np.float64),
               'nbl': nbl,
               'nwav': self.nwav,
               'wav': self.wav,
               'u': np.zeros([nbl, self.nwav], dtype=np.float64),
               'v': np.zeros([nbl, self.nwav], dtype=np.float64),
               'vis': np.zeros([nbl, self.nwav], dtype=np.complex64),
               'amp': np.zeros([nbl, self.nwav], dtype=np.float64),
               'phase': np.zeros([nbl, self.nwav], dtype=np.float64)}

        # res = {}
        # res['bl'] = np.array(bl, dtype=float)
        # res['pa'] = np.array(pa, dtype=float)
        # res['nbl'] = res['bl'].shape[0]
        #
        # res['wav'] = np.array(self.wav)
        # res['nwav'] = self.nwav
        #
        # res['u'] = np.zeros([res['nbl'], self.nwav], dtype=float)
        # res['v'] = np.zeros([res['nbl'], self.nwav], dtype=float)
        #
        # res['vis'] = np.zeros([res['nbl'], self.nwav], dtype=np.complex64)
        # res['amp'] = np.zeros([res['nbl'], self.nwav], dtype=np.float64)
        # res['phase'] = np.zeros([res['nbl'], self.nwav], dtype=np.float64)

        l = self.x / nc.au / dpc / 3600. / 180. * np.pi
        m = self.y / nc.au / dpc / 3600. / 180. * np.pi
        dl = l[1] - l[0]
        dm = m[1] - m[0]

        for iwav in range(res['nwav']):

            # Calculate spatial frequencies
            res['u'][:, iwav] = res['bl'] * np.cos(res['pa']) * 1e6 / self.wav[iwav]
            res['v'][:, iwav] = res['bl'] * np.sin(res['pa']) * 1e6 / self.wav[iwav]

            for ibl in range(res['nbl']):
                dum = complex(0.)
                imu = complex(0., 1.)

                for il in range(len(l)):
                    phase = 2. * np.pi * (res['u'][ibl, iwav] * l[il] + res['v'][ibl, iwav] * m)
                    cterm = np.cos(phase)
                    sterm = -np.sin(phase)
                    dum = dum + (self.image[il, :, iwav] * (cterm + imu * sterm)).sum() * dl * dm

                res['vis'][ibl, iwav] = dum
                res['amp'][ibl, iwav] = np.sqrt(abs(dum * np.conj(dum)))
                res['phase'][ibl, iwav] = np.arccos(np.real(dum) / res['amp'][ibl, iwav])
                if np.imag(dum) < 0.:
                    res['phase'][ibl, iwav] = 2. * np.pi - res['phase'][ibl, iwav]

                print('Calculating baseline # : ', ibl, ' wavelength # : ', iwav)

        return res

    # --------------------------------------------------------------------------------------------------
    def writeFits(self, fname='', dpc=1., coord='03h10m05s -10d05m30s', bandwidthmhz=2000.0,
                  casa=False, nu0=0., stokes='I', fitsheadkeys=[], ifreq=None):
        """Writes out a RADMC-3D image data in fits format.

        Parameters
        ----------

        fname        : str
                        File name of the radmc3d output image (if omitted 'image.fits' is used)

        dpc          : float
                        Distance of the source in pc

        coord        : str
                        Image center coordinates

        bandwidthmhz : float
                        Bandwidth of the image in MHz (equivalent of the CDELT keyword in the fits header)

        casa         : bool
                        If set to True a CASA compatible four dimensional image cube will be written

        nu0          : float
                        Rest frequency of the line (for channel maps)

        stokes       : {'I', 'Q', 'U', 'V', 'PI'}
                       Stokes parameter to be written if the image contains Stokes IQUV (possible
                       choices: 'I', 'Q', 'U', 'V', 'PI' -Latter being the polarized intensity)

        fitsheadkeys : dictionary
                        Dictionary containing all (extra) keywords to be added to the fits header. If
                        the keyword is already in the fits header (e.g. CDELT1) it will be updated/changed
                        to the value in fitsheadkeys, if the keyword is not present the keyword is added to
                        the fits header.

        ifreq        : int
                       Frequency index of the image array to write. If set only this frequency of a multi-frequency
                       array will be written to file.
        """
        # --------------------------------------------------------------------------------------------------
        istokes = 0


        if self.stokes:
            if fname == '':
                fname = 'image_stokes_' + stokes.strip().upper() + '.fits'

            if stokes.strip().upper() == 'I':
                istokes = 0
            if stokes.strip().upper() == 'Q':
                istokes = 1
            if stokes.strip().upper() == 'U':
                istokes = 2
            if stokes.strip().upper() == 'V':
                istokes = 3
        else:
            if fname == '':
                fname = 'image.fits'

        # Decode the image center cooridnates
        # Check first whether the format is OK
        dum = coord

        ra = []
        delim = ['h', 'm', 's']
        for i in delim:
            ind = dum.find(i)
            if ind <= 0:
                msg = 'coord keyword has a wrong format. The format should be coord="0h10m05s -10d05m30s"'
                raise ValueError(msg)
            ra.append(float(dum[:ind]))
            dum = dum[ind + 1:]

        dec = []
        delim = ['d', 'm', 's']
        for i in delim:
            ind = dum.find(i)
            if ind <= 0:
                msg = 'coord keyword has a wrong format. The format should be coord="0h10m05s -10d05m30s"'
                raise ValueError(msg)
            dec.append(float(dum[:ind]))
            dum = dum[ind + 1:]

        target_ra = (ra[0] + ra[1] / 60. + ra[2] / 3600.) * 15.
        if dec[0] >= 0:
            target_dec = (dec[0] + dec[1] / 60. + dec[2] / 3600.)
        else:
            target_dec = (dec[0] - dec[1] / 60. - dec[2] / 3600.)

        if len(self.fwhm) == 0:
            # Conversion from erg/s/cm/cm/ster to Jy/pixel
            conv = self.sizepix_x * self.sizepix_y / (dpc * nc.pc)**2. * 1e23
        else:
            # If the image has already been convolved with a gaussian psf then self.image has
            # already the unit of erg/s/cm/cm/beam, so we need to multiply it by 10^23 to get
            # to Jy/beam
            conv = 1e23

        # Create the data to be written
        if casa:
            # Put the stokes axis to the 4th dimension
            # data = np.zeros([1, self.nfreq, self.ny, self.nx], dtype=float)
            data = np.zeros([1, self.nfreq, self.ny, self.nx], dtype=float)
            if self.nfreq == 1:
                data[0, 0, :, :] = self.image[:, :] * conv

            else:
                for inu in range(self.nfreq):
                    data[0, inu, :, :] = self.image[:, :, inu] * conv
        else:
            data = np.zeros([self.nfreq, self.nx, self.ny], dtype=float)
            if self.stokes:
                if stokes.strip().upper() != 'PI':
                    if self.nfreq == 1:
                        data[0, :, :] = self.image[:, :, istokes, 0] * conv

                    else:
                        for inu in range(self.nfreq):
                            data[inu, :, :] = self.image[:, :, istokes, inu] * conv
                else:
                    if self.nfreq == 1:
                        data[0, :, :] = np.sqrt(self.image[:, :, 1, 0]**2 + self.image[:, :, 2, 0]**2) * conv

                    else:
                        for inu in range(self.nfreq):
                            data[inu, :, :] = np.sqrt(
                                self.image[:, :, 1, inu]**2 + self.image[:, :, 2, inu]**2) * conv

            else:
                if self.nfreq == 1:
                    data[0, :, :] = self.image[:, :, 0] * conv

                else:
                    for inu in range(self.nfreq):
                        data[inu, :, :] = self.image[:, :, inu] * conv

        if ifreq is not None:
            if len(data.shape) == 3:
                data = data[ifreq, :, :]

        naxis = len(data.shape)
        hdu = pf.PrimaryHDU(data.swapaxes(naxis - 1, naxis - 2))
        hdulist = pf.HDUList([hdu])

        hdulist[0].header.set('CRPIX1', (self.nx + 1.) / 2., ' ')
        hdulist[0].header.set('CDELT1', -self.sizepix_x / nc.au / dpc / 3600., '')
        # hdulist[0].header.set('CRVAL1', self.sizepix_x/1.496e13/dpc/3600.*0.5+target_ra, '')
        hdulist[0].header.set('CRVAL1', target_ra, '')
        hdulist[0].header.set('CUNIT1', '     DEG', '')
        hdulist[0].header.set('CTYPE1', 'RA---SIN', '')

        hdulist[0].header.set('CRPIX2', (self.ny + 1.) / 2., '')
        hdulist[0].header.set('CDELT2', self.sizepix_y / nc.au / dpc / 3600., '')
        # hdulist[0].header.set('CRVAL2', self.sizepix_y/1.496e13/dpc/3600.*0.5+target_dec, '')
        hdulist[0].header.set('CRVAL2', target_dec, '')
        hdulist[0].header.set('CUNIT2', '     DEG', '')
        hdulist[0].header.set('CTYPE2', 'DEC--SIN', '')

        # For ARTIST compatibility put the stokes axis to the 4th dimension
        if casa:
            hdulist[0].header.set('CRPIX4', 1., '')
            hdulist[0].header.set('CDELT4', 1., '')
            hdulist[0].header.set('CRVAL4', 1., '')
            hdulist[0].header.set('CUNIT4', '        ', '')
            hdulist[0].header.set('CTYPE4', 'STOKES  ', '')

            if self.nwav == 1:
                hdulist[0].header.set('CRPIX3', 1.0, '')
                hdulist[0].header.set('CDELT3', bandwidthmhz * 1e6, '')
                hdulist[0].header.set('CRVAL3', self.freq[0], '')
                hdulist[0].header.set('CUNIT3', '      HZ', '')
                hdulist[0].header.set('CTYPE3', 'FREQ-LSR', '')

            else:
                if ifreq is None:
                    hdulist[0].header.set('CRPIX3', 1.0, '')
                    hdulist[0].header.set('CDELT3', (self.freq[1] - self.freq[0]), '')
                    hdulist[0].header.set('CRVAL3', self.freq[0], '')
                    hdulist[0].header.set('CUNIT3', '      HZ', '')
                    hdulist[0].header.set('CTYPE3', 'FREQ-LSR', '')
                    hdulist[0].header.set('RESTFRQ', self.freq[0])
                else:
                    hdulist[0].header.set('CRPIX3', 1.0, '')
                    hdulist[0].header.set('CDELT3', bandwidthmhz * 1e6, '')
                    hdulist[0].header.set('CRVAL3', self.freq[ifreq], '')
                    hdulist[0].header.set('CUNIT3', '      HZ', '')
                    hdulist[0].header.set('CTYPE3', 'FREQ-LSR', '')

        else:
            if self.nwav == 1:
                hdulist[0].header.set('CRPIX3', 1.0, '')
                hdulist[0].header.set('CDELT3', bandwidthmhz * 1e6, '')
                hdulist[0].header.set('CRVAL3', self.freq[0], '')
                hdulist[0].header.set('CUNIT3', '      HZ', '')
                hdulist[0].header.set('CTYPE3', 'FREQ-LSR', '')
            else:
                if ifreq is None:
                    hdulist[0].header.set('CRPIX3', 1.0, '')
                    hdulist[0].header.set('CDELT3', (self.freq[1] - self.freq[0]), '')
                    hdulist[0].header.set('CRVAL3', self.freq[0], '')
                    hdulist[0].header.set('CUNIT3', '      HZ', '')
                    hdulist[0].header.set('CTYPE3', 'FREQ-LSR', '')
                else:
                    hdulist[0].header.set('CRPIX3', 1.0, '')
                    hdulist[0].header.set('CDELT3', bandwidthmhz * 1e6, '')
                    hdulist[0].header.set('CRVAL3', self.freq[ifreq], '')
                    hdulist[0].header.set('CUNIT3', '      HZ', '')
                    hdulist[0].header.set('CTYPE3', 'FREQ-LSR', '')

        if nu0 > 0:
            hdulist[0].header.set('RESTFRQ', nu0, '')
        else:
            if self.nwav == 1:
                hdulist[0].header.set('RESTFRQ', self.freq[0], '')

        if len(self.fwhm) == 0:
            hdulist[0].header.set('BUNIT', 'JY/PIXEL', '')
        else:
            hdulist[0].header.set('BUNIT', 'JY/BEAM', '')
            hdulist[0].header.set('BMAJ', self.fwhm[0] / 3600., '')
            hdulist[0].header.set('BMIN', self.fwhm[1] / 3600., '')
            hdulist[0].header.set('BPA', -self.pa, '')

        hdulist[0].header.set('BTYPE', 'INTENSITY', '')
        hdulist[0].header.set('BZERO', 0.0, '')
        hdulist[0].header.set('BSCALE', 1.0, '')

        hdulist[0].header.set('EPOCH', 2000.0, '')
        hdulist[0].header.set('LONPOLE', 180.0, '')

        if fitsheadkeys:
            if len(fitsheadkeys.keys()) > 0:
                for ikey in fitsheadkeys.keys():
                    # hdulist[0].header.update(ikey, fitsheadkeys[ikey], '')
                    hdulist[0].header.set(ikey, fitsheadkeys[ikey], '')

        if os.path.exists(fname):
            print(fname + ' already exists')
            dum = input('Do you want to overwrite it (yes/no)?')
            if (dum.strip()[0] == 'y') | (dum.strip()[0] == 'Y'):
                os.remove(fname)
                hdu.writeto(fname)
            else:
                print('No image has been written')
        else:
            hdu.writeto(fname)
            # --------------------------------------------------------------------------------------------------

    def plotMomentMap(self, moment=0, nu0=None, wav0=None, dpc=1., au=False, arcsec=False, cmap=None, vclip=None):
        """Plots moment maps

        Parameters
        ----------

        moment : int
                 Moment of the channel maps to be calculated

        nu0    : float
                 Rest frequency of the line in Hz

        wav0   : float
                 Rest wavelength of the line in micron

        dpc    : float
                 Distance of the source in pc

        au     : bool
                 If True displays the image with au as the spatial axis unit

        arcsec : bool
                 If True displays the image with arcsec as the spatial axis unit (dpc should also be set!)

        cmap   : matplotlib colormap
                 Color map to be used to display the moment map

        vclip  : list/ndarray
                 Two element list / Numpy array containin the lower and upper limits for the values in the moment
                  map to be displayed

        """

        # I/O error handling
        if nu0 is None:
            if wav0 is None:
                msg = 'Unknown nu0 and wav0. Neither the rest frequency (nu0) nor the rest wavelength (wav0)'\
                      + ' of the line is specified.'
                raise ValueError(msg)
            else:
                nu0 = nc.cc / wav0 * 1e4

        if len(self.image.shape) != 3:
            msg = 'Wrong image shape. Channel map calculation requires a three dimensional array with '\
                  + '[Nx,  Ny,  Nnu] dimensions. The current image has the shape of ' + str(len(self.image.shape))
            raise ValueError(msg)
        mmap = self.getMomentMap(moment=moment, nu0=nu0, wav0=wav0)

        if moment > 0:
            mmap0 = self.getMomentMap(moment=0, nu0=nu0, wav0=wav0)
            mmap = mmap / mmap0

        # Select the coordinates of the data
        if au:
            x = self.x / nc.au
            y = self.y / nc.au
            xlab = 'X [au]'
            ylab = 'Y [au]'
        elif arcsec:
            x = self.x / nc.au / dpc
            y = self.y / nc.au / dpc
            xlab = 'RA offset ["]'
            ylab = 'DEC offset ["]'
        else:
            x = self.x
            y = self.y
            xlab = 'X [cm]'
            ylab = 'Y [cm]'

        ext = (x[0], x[self.nx - 1], y[0], y[self.ny - 1])

        cb_label = ''
        if moment == 0:
            mmap = mmap / (dpc * dpc)
            cb_label = 'I' + r'$_\nu$' + ' [erg/s/cm/cm/Hz/ster*km/s]'
        if moment == 1:
            mmap = mmap / (dpc * dpc)
            cb_label = 'v [km/s]'
        if moment > 1:
            mmap = mmap / (dpc * dpc)
            powex = str(moment)
            cb_label = r'v$^' + powex + '$ [(km/s)$^' + powex + '$]'

        if vclip is not None:
            if len(vclip) != 2:
                msg = 'Wrong shape in vclip. vclip should be a two element list with (clipmin, clipmax)'
                raise ValueError(msg)
            else:
                mmap = mmap.clip(vclip[0], vclip[1])

        implot = plt.imshow(mmap, extent=ext, cmap=cmap)
        cbar = plt.colorbar(implot)
        cbar.set_label(cb_label)
        plt.xlabel(xlab)
        plt.ylabel(ylab)

    def getMomentMap(self, moment=0, nu0=None, wav0=None):
        """Calculates moment maps.

        Parameters
        ----------

        moment : int
                 Moment of the channel maps to be calculated

        nu0    : float
                 Rest frequency of the line in Hz

        wav0   : float
                 Rest wavelength of the line in micron

        Returns
        -------
        Ndarray with the same dimension as the individual channel maps
        """

        # I/O error handling
        if nu0 is None:
            if wav0 is None:
                msg = 'Unknown nu0 and wav0. Neither the rest frequency (nu0) nor the rest wavelength (wav0)' \
                      + ' of the line is specified.'
                raise ValueError(msg)
            else:
                nu0 = nc.cc / wav0 * 1e4

        if len(self.image.shape) != 3:
            msg = 'Wrong image shape. Channel map calculation requires a three dimensional array with '\
                  + '[Nx,  Ny,  Nnu] dimensions. The current image has the shape of ' + str(len(self.image.shape))
            raise ValueError(msg)

        # First calculate the velocity field
        v_kms = nc.cc * (nu0 - self.freq) / nu0 / 1e5

        vmap = np.zeros([self.nx, self.ny, self.nfreq], dtype=np.float64)
        for ifreq in range(self.nfreq):
            vmap[:, :, ifreq] = v_kms[ifreq]

        # Now calculate the moment map
        y = self.image * (vmap**moment)

        dum = (vmap[:, :, 1:] - vmap[:, :, :-1]) * (y[:, :, 1:] + y[:, :, :-1]) * 0.5

        return dum.sum(2)

    def readImage(self, fname=None, binary=False, old=False):
        """Reads an image calculated by RADMC-3D

        Parameters
        ----------

        fname   : str, optional
                 File name of the radmc3d output image (if omitted 'image.out' is used)

        old     : bool
                 If set to True it reads old radmc-2d style image

        binary  : bool, optional
                 False - the image format is formatted ASCII if True - C-compliant binary (omitted if old=True)
        """
        if old:
            if fname is None:
                fname = 'image.dat'

            self.filename = fname
            print('Reading ' + fname)


            with open(fname, 'r') as rfile:

                dum = rfile.readline().split()
                self.nx = int(dum[0])
                self.ny = int(dum[1])
                self.nfreq = int(dum[2])
                self.nwav = self.nfreq

                dum = rfile.readline().split()
                self.sizepix_x = float(dum[0])
                self.sizepix_y = float(dum[1])
                self.wav = np.zeros(self.nwav, dtype=float) - 1.
                self.freq = np.zeros(self.nwav, dtype=float) - 1.

                self.stokes = False
                self.image = np.zeros([self.nx, self.ny, self.nwav], dtype=np.float64)
                for iwav in range(self.nwav):
                    dum = rfile.readline()
                    for iy in range(self.ny):
                        for ix in range(self.nx):
                            self.image[ix, iy, iwav] = float(rfile.readline())

        else:
            if binary:
                if fname is None:
                    fname = 'image.bout'

                self.filename = fname

                dum = np.fromfile(fname, count=4, dtype=int)
                iformat = dum[0]
                self.nx = dum[1]
                self.ny = dum[2]
                self.nfreq = dum[3]
                self.nwav = self.nfreq
                dum = np.fromfile(fname, count=-1, dtype=np.float64)

                self.sizepix_x = dum[4]
                self.sizepix_y = dum[5]
                self.wav = dum[6:6 + self.nfreq]
                self.freq = nc.cc / self.wav * 1e4

                if iformat == 1:
                    self.stokes = False
                    self.image = np.reshape(dum[6 + self.nfreq:], [self.nfreq, self.ny, self.nx])
                    self.image = np.swapaxes(self.image, 0, 2)
                elif iformat == 3:
                    self.stokes = True
                    self.image = np.reshape(dum[6 + self.nfreq:], [self.nfreq, 4, self.ny, self.nx])
                    self.image = np.swapaxes(self.image, 0, 3)
                    self.image = np.swapaxes(self.image, 1, 2)

            else:

                # Look for the image file

                if fname is None:
                    fname = 'image.out'

                print('Reading '+ fname)

                self.filename = fname
                with open(fname, 'r') as rfile:

                    dum = ''

                    # Format number
                    iformat = int(rfile.readline())

                    # Nr of pixels
                    dum = rfile.readline()
                    dum = dum.split()
                    self.nx = int(dum[0])
                    self.ny = int(dum[1])
                    # Nr of frequencies
                    self.nfreq = int(rfile.readline())
                    self.nwav = self.nfreq
                    # Pixel sizes
                    dum = rfile.readline()
                    dum = dum.split()
                    self.sizepix_x = float(dum[0])
                    self.sizepix_y = float(dum[1])
                    # Wavelength of the image
                    self.wav = np.zeros(self.nwav, dtype=np.float64)
                    for iwav in range(self.nwav):
                        self.wav[iwav] = float(rfile.readline())
                    self.wav = np.array(self.wav)
                    self.freq = nc.cc / self.wav * 1e4

                    # If we have a normal total intensity image
                    if iformat == 1:
                        self.stokes = False

                        self.image = np.zeros([self.nx, self.ny, self.nwav], dtype=np.float64)
                        for iwav in range(self.nwav):
                            # Blank line
                            dum = rfile.readline()
                            for iy in range(self.ny):
                                for ix in range(self.nx):
                                    self.image[ix, iy, iwav] = float(rfile.readline())

                    # If we have the full stokes image
                    elif iformat == 3:
                        self.stokes = True
                        self.image = np.zeros([self.nx, self.ny, 4, self.nwav], dtype=np.float64)
                        for iwav in range(self.nwav):
                            # Blank line
                            dum = rfile.readline()
                            for iy in range(self.ny):
                                for ix in range(self.nx):
                                    dum = rfile.readline().split()
                                    imstokes = [float(i) for i in dum]
                                    self.image[ix, iy, 0, iwav] = float(dum[0])
                                    self.image[ix, iy, 1, iwav] = float(dum[1])
                                    self.image[ix, iy, 2, iwav] = float(dum[2])
                                    self.image[ix, iy, 3, iwav] = float(dum[3])

        # Conversion from erg/s/cm/cm/Hz/ster to Jy/pixel
        conv = self.sizepix_x * self.sizepix_y / nc.pc**2. * 1e23
        self.imageJyppix = self.image * conv

        self.x = ((np.arange(self.nx, dtype=np.float64) + 0.5) - self.nx / 2) * self.sizepix_x
        self.y = ((np.arange(self.ny, dtype=np.float64) + 0.5) - self.ny / 2) * self.sizepix_y

    def imConv(self, dpc=1., psfType='gauss', fwhm=None, pa=None, tdiam_prim=8.2, tdiam_sec=0.94):
        """Convolves a RADMC-3D image with a two dimensional Gaussian psf. The output images will have the same
        brightness units as the input images.

        Parameters
        ----------
        dpc         : float
                      Distance of the source in pc.

        psfType     : {'gauss', 'airy'}
                      Shape of the PSF. If psfType='gauss', fwhm and pa should also be given. If psfType='airy', the
                      tdiam_prim, tdiam_sec and wav parameters should also be specified.

        fwhm        : list, optional
                      A list of two numbers; the FWHM of the two dimensional psf along the two principal axes.
                      The unit is assumed to be arcsec. (should only be set if psfType='gauss')

        pa          : float, optional
                      Position angle of the psf ellipse (counts from North counterclockwise, should only be set
                      if psfType='gauss')

        tdiam_prim  : float, optional
                      Diameter of the primary aperture of the telescope in meter. (should be set only if psfType='airy')

        tdiam_sec   : float, optional
                      Diameter of the secondary mirror (central obscuration), if there is any, in meter. If no
                      secondary mirror/obscuration is present, this parameter should be set to zero.
                      (should only be set if psfType='airy')

        Returns
        -------

        Returns a radmc3dImage
        """

        dx = self.sizepix_x / nc.au / dpc
        dy = self.sizepix_y / nc.au / dpc
        nfreq = self.nfreq
        psf = None
        cimage = None

        if self.stokes:
            if self.nfreq == 1:
                # Generate the  psf
                dum = getPSF(nx=self.nx, ny=self.ny, pscale=[dx, dy], psfType=psfType, fwhm=fwhm, pa=pa,
                             tdiam_prim=tdiam_prim, tdiam_sec=tdiam_sec, wav=self.wav[0])
                psf = dum['psf']
                f_psf = np.fft.fft2(psf)

                cimage = np.zeros([self.nx, self.ny, 4, 1], dtype=np.float64)
                for istokes in range(4):
                    imag = self.image[:, :, istokes, 0]
                    f_imag = np.fft.fft2(imag)
                    f_cimag = f_psf * f_imag
                    cimage[:, :, istokes, 0] = np.real(np.fft.ifftshift(np.fft.ifft2(f_cimag)))
            else:
                # If we have a simple Gaussian PSF it will be wavelength independent so we can take it out from the
                # frequency loop
                if psfType.lower().strip() == 'gauss':
                    # Generate the wavelength independent gaussian psf
                    dum = getPSF(nx=self.nx, ny=self.ny, pscale=[dx, dy], psfType=psfType, fwhm=fwhm, pa=pa,
                                 tdiam_prim=tdiam_prim, tdiam_sec=tdiam_sec, wav=self.wav[0])
                    psf = dum['psf']
                    f_psf = np.fft.fft2(psf)

                    cimage = np.zeros([self.nx, self.ny, 4, self.nfreq], dtype=np.float64)
                    for ifreq in range(nfreq):
                        for istokes in range(4):
                            imag = self.image[:, :, istokes, ifreq]
                            f_imag = np.fft.fft2(imag)
                            f_cimag = f_psf * f_imag
                            cimage[:, :, istokes, ifreq] = np.real(np.fft.ifftshift(np.fft.ifft2(f_cimag)))

                # If we have an Airy-PSF calculated from the aperture size(s) and wavelenght, the PSF will depend
                # on the frequency so it has to be re-calculated for each wavelength
                elif psfType.lower().strip() == 'airy':
                    cimage = np.zeros([self.nx, self.ny, 4, self.nfreq], dtype=np.float64)
                    for ifreq in range(nfreq):
                        # Generate the wavelength-dependent airy-psf
                        dum = getPSF(nx=self.nx, ny=self.ny, pscale=[dx, dy], psfType=psfType, fwhm=fwhm, pa=pa,
                                     tdiam_prim=tdiam_prim, tdiam_sec=tdiam_sec, wav=self.wav[ifreq])
                        psf = dum['psf']
                        f_psf = np.fft.fft2(psf)

                        for istokes in range(4):
                            imag = self.image[:, :, istokes, ifreq]
                            f_imag = np.fft.fft2(imag)
                            f_cimag = f_psf * f_imag
                            cimage[:, :, istokes, ifreq] = np.real(np.fft.ifftshift(np.fft.ifft2(f_cimag)))

        else:
            # If we have a simple Gaussian PSF it will be wavelength independent so we can take it out from the
            # frequency loop
            if psfType.lower().strip() == 'gauss':
                # Generate the wavelength independent gaussian psf
                dum = getPSF(nx=self.nx, ny=self.ny, pscale=[dx, dy], psfType=psfType, fwhm=fwhm, pa=pa,
                             tdiam_prim=tdiam_prim, tdiam_sec=tdiam_sec, wav=self.wav[0])
                psf = dum['psf']
                f_psf = np.fft.fft2(psf)

                cimage = np.zeros([self.nx, self.ny, self.nfreq], dtype=np.float64)
                for ifreq in range(nfreq):
                    imag = self.image[:, :, ifreq]
                    f_imag = np.fft.fft2(imag)
                    f_cimag = f_psf * f_imag
                    cimage[:, :, ifreq] = np.real(np.fft.ifftshift(np.fft.ifft2(f_cimag)))

            # If we have an Airy-PSF calculated from the aperture size(s) and wavelenght, the PSF will depend
            # on the frequency so it has to be re-calculated for each wavelength
            elif psfType.lower().strip() == 'airy':
                cimage = np.zeros([self.nx, self.ny, self.nfreq], dtype=np.float64)
                for ifreq in range(nfreq):
                    # Generate the wavelength-dependent airy-psf
                    dum = getPSF(nx=self.nx, ny=self.ny, pscale=[dx, dy], psfType=psfType, fwhm=fwhm, pa=pa,
                                 tdiam_prim=tdiam_prim, tdiam_sec=tdiam_sec, wav=self.wav[ifreq])
                    psf = dum['psf']
                    f_psf = np.fft.fft2(psf)

                    imag = self.image[:, :, ifreq]
                    f_imag = np.fft.fft2(imag)
                    f_cimag = f_psf * f_imag
                    cimage[:, :, ifreq] = np.real(np.fft.ifftshift(np.fft.ifft2(f_cimag)))

        # cimage = squeeze(cimage)

        # Return the convolved image (copy the image class and replace the image attribute to the convolved image)
        res = copy.deepcopy(self)
        res.image = cimage
        conv = self.sizepix_x * self.sizepix_y / nc.pc**2. * 1e23
        res.imageJyppix = res.image * conv

        res.psf = psf
        res.fwhm = fwhm
        res.pa = pa
        res.dpc = dpc

        return res

    def compute_brightness_temperature(self,linear=True):
        """
        Compute, from the image intensity, the corresponding brightness temperature in Kelvin.

        OPTIONAL ARGUMENT:

          linear         If set, then just use the linear conversion factor, not
                         taking into account the exponential Wien part of the
                         Planck spectrum. If unset, then compute the true brightness
                         temperature, using the full Planck spectrum.
        """
        tbright = np.zeros((self.nx,self.ny,self.nfreq))
        for inu in range(self.nfreq):
            freq = self.freq[inu]
            if self.nfreq == 1:
                if self.stokes:
                    img = self.image[:,:,0]
                else:
                    img = self.image[:,:]
            else:
                if self.stokes:
                    img = self.image[:,:,inu,0]
                else:
                    img = self.image[:,:,inu]
            img = np.squeeze(img)
            if linear:
                tbright[:,:,inu] = img*(nc.cc**2)/(2*(freq**2)*nc.kk)
            else:
                tbright[:,:,inu] = nc.hh*freq/(nc.kk*np.log(2*nc.hh*freq**3/nc.cc**2*1./img+1.))
        return tbright


def getPSF(nx=None, ny=None, psfType='gauss', pscale=None, fwhm=None, pa=None, tdiam_prim=8.2, tdiam_sec=0.94,
           wav=None):
    """Calculates a two dimensional Gaussian PSF.

    Parameters
    ----------
    nx          : int
                  Image size in the first dimension

    ny          : int
                  Image size in the second dimension

    psfType     : {'gauss', 'airy'}
                  Shape of the PSF. If psfType='gauss', fwhm and pa should also be given. If psfType='airy', the
                  tdiam_prim, tdiam_sec and wav parameters should also be specified.

    pscale      : list
                  Pixelscale of the image, if set fwhm should be in the same unit, if not set unit of fwhm is pixels

    fwhm        : list, optional
                  Full width at half maximum of the psf along the two axis (should be set only if psfType='gauss')

    pa          : float, optional
                  Position angle of the gaussian if the gaussian is not symmetric
                  (should be set only if psfType='gauss')

    tdiam_prim  : float, optional
                  Diameter of the primary aperture of the telescope in meter. (should be set only if psfType='airy')

    tdiam_sec   : float, optional
                  Diameter of the secondary mirror (central obscuration), if there is any, in meter. If no secondary
                  mirror/obscuration is present, this parameter should be set to zero.
                  (should be set only if psfType='airy')

    wav         : float, optional
                  Wavelength of observation in micrometer (should be set only if psfType='airy')


    Returns
    -------

    Returns a dictionary with the following keys:

        * psf : ndarray
                The two dimensional psf
        * x   : ndarray
                The x-axis of the psf
        * y   : ndarray
                The y-axis of the psf
    """
    # --------------------------------------------------------------------------------------------------

    # Create the two axes

    if pscale is not None:
        dx, dy = pscale[0], pscale[1]
    else:
        dx, dy = 1., 1.

    x = (np.arange(nx, dtype=np.float64) - nx / 2) * dx
    y = (np.arange(ny, dtype=np.float64) - ny / 2) * dy

    # Create the Gaussian PSF
    psf = None
    if psfType.strip().lower() == 'gauss':

        # Calculate the standard deviation of the Gaussians
        sigmax = fwhm[0] / (2.0 * np.sqrt(2.0 * np.log(2.)))
        sigmay = fwhm[1] / (2.0 * np.sqrt(2.0 * np.log(2.)))
        norm = (2. * np.pi * sigmax * sigmay) / dx / dy
        # Pre-compute sin and cos angles

        sin_pa = np.sin(pa / 180. * np.pi - np.pi / 2.)
        cos_pa = np.cos(pa / 180. * np.pi - np.pi / 2.)

        # Define the psf
        psf = np.zeros([nx, ny], dtype=np.float64)
        cos_pa_x = cos_pa * x
        cos_pa_y = cos_pa * y
        sin_pa_x = sin_pa * x
        sin_pa_y = sin_pa * y
        for ix in range(nx):
            for iy in range(ny):
                # xx = cos_pa * x[ix] - sin_pa * y[iy]
                # yy = sin_pa * x[ix] + cos_pa * y[iy]
                xx = cos_pa_x[ix] - sin_pa_y[iy]
                yy = sin_pa_x[ix] + cos_pa_y[iy]

                psf[ix, iy] = np.exp(-0.5 * xx * xx / sigmax / sigmax - 0.5 * yy * yy / sigmay / sigmay)

        psf /= norm

    elif psfType.strip().lower() == 'airy':

        # Check whether scipy was successfully imported
        if not spc:
            msg = 'scipy.special was not imported. PSF calculation is limited to Gaussian only.'
            raise ImportError(msg)

        # Unit conversion
        x_rad = x / 3600. / 180. * np.pi
        y_rad = y / 3600. / 180. * np.pi
        x2 = x_rad**2
        y2 = y_rad**2
        wav_m = wav * 1e-6
        psf = np.zeros([nx, ny], dtype=np.float64)
        if tdiam_sec == 0.:
            for ix in range(nx):
                r = np.sqrt(x2[ix] + y2)
                u = np.pi / wav_m * tdiam_prim * r

                if 0. in u:
                    ii = (u == 0.)
                    u[ii] = 1e-5
                psf[ix, :] = (2.0 * spc.j1(u) / u)**2
        else:
            for ix in range(nx):
                r = np.sqrt(x2[ix] + y2)
                u = np.pi / wav_m * tdiam_prim * r
                eps = tdiam_sec / tdiam_prim
                if 0. in u:
                    ii = (u == 0.)
                    u[ii] = 1e-5
                psf[ix, :] = 1.0 / (1.0 - eps**2)**2 * ((2.0 * spc.j1(u) / u)
                                                        - (eps**2 * 2.0 * spc.j1(eps * u) / (eps * u)))**2

        dum = 0.44 * (wav * 1e-6 / tdiam_prim / np.pi * 180. * 3600.) * 2. * np.sqrt(2. * np.log(2.))
        fwhm = [dum, dum]
        norm = fwhm[0] * fwhm[1] * np.pi / (4. * np.log(2.)) / dx / dy
        psf /= norm

    res = {'psf': psf, 'x': x, 'y': y}

    return res


def readImage(fname=None, binary=False, old=False):
    """Reads an image calculated by RADMC-3D.
       This function is an interface to radmc3dImage.readImage().

    Parameters
    ----------
        fname   : str, optional
                 File name of the radmc3d output image (if omitted 'image.out' is used)

        old     : bool
                 If set to True it reads old radmc-2d style image

        binary  : bool, optional
                 False - the image format is formatted ASCII if True - C-compliant binary (omitted if old=True)
    """

    dum = radmc3dImage()
    dum.readImage(fname=fname, binary=binary, old=old)
    return dum


def plotPolDir(image=None, arcsec=False, au=False, dpc=None, ifreq=0, cmask_rad=None, color='w', nx=20, ny=20):
    """
    Function to plot the polarisation direction for full stokes images

    Parameters
    ----------

    image         : radmc3dImage
                    A radmc3dImage class returned by readimage

    arcsec        : bool
                    If True image axis will have the unit arcsec (NOTE: dpc keyword should also be set!)

    au            : bool
                    If True image axis will have the unit au

    dpc           : float
                    Distance to the source in parsec (This keywords should be set if arcsec=True, or bunit!='norm')

    ifreq         : int
                    If the image file/array consists of multiple frequencies/wavelengths ifreq denotes the index
                    of the frequency/wavelength in the image array to be plotted

    cmask_rad     : float
                    Simulates coronographyic mask : sets the image values to zero within this radius of the image center
                    The unit is the same as the image axis (au, arcsec, cm)
                    NOTE: this works only on the plot, the image array is not changed (for that used the cmask()
                    function)

    color         : str
                    Color for the polarisation direction plot

    nx            : int
                    Number of grid points along the horizontal axis at which the direction should be displayed

    ny            : int
                    Number of grid points along the vertical axis at which the direction should be displayed
    """

    #
    # First check if the images is a full stokes image
    #

    if not image.stokes:
        msg = 'The image is not a full stokes image. Polarisation direction can only be displayed if '\
              + 'the full stokes vector is present at every pixel of the image'
        raise ValueError(msg)

    if cmask_rad is not None:
        dum_image = cmask(image, rad=cmask_rad, au=au, arcsec=arcsec, dpc=dpc)
    else:
        dum_image = copy.deepcopy(image)

    # Select the coordinates of the data
    if au:
        x = image.x / nc.au
        y = image.y / nc.au
        xlab = 'X [au]'
        ylab = 'Y [au]'
    elif arcsec:
        x = image.x / nc.au / dpc
        y = image.y / nc.au / dpc
        xlab = 'RA offset ["]'
        ylab = 'DEC offset ["]'
    else:
        x = image.x
        y = image.y
        xlab = 'X [cm]'
        ylab = 'Y [cm]'

    # ext = (x[0], x[image.nx-1], y[0], y[image.ny-1])
    iix = [int(np.floor(i)) for i in np.arange(nx) * float(x.shape[0]) / nx]
    iiy = [int(np.floor(i)) for i in np.arange(ny) * float(x.shape[0]) / ny]
    xr = x[iix]
    yr = y[iiy]
    xxr, yyr = np.meshgrid(xr, yr, indexing='ij')
    qqr = (np.squeeze(dum_image.image[:, :, 1, ifreq])
           / np.squeeze(dum_image.image[:, :, 0, ifreq]).clip(1e-60))[np.ix_(iix, iiy)]
    uur = (np.squeeze(dum_image.image[:, :, 2, ifreq])
           / np.squeeze(dum_image.image[:, :, 0, ifreq]).clip(1e-60))[np.ix_(iix, iiy)]
    lpol = np.sqrt(qqr**2 + uur**2).clip(1e-60)
    qqr /= lpol
    uur /= lpol
    ang = np.arccos(qqr) / 2.0
    ii = (uur < 0)
    if True in ii:
        ang[ii] = np.pi - ang[ii]
    vx = np.cos(ang)
    vy = np.sin(ang)
    ii = (lpol < 1e-6)
    vx[ii] = 0.001
    vy[ii] = 0.001
    plt.quiver(xxr, yyr, vx, vy, color=color, pivot='mid', scale=2. * np.max([nx, ny]), headwidth=1e-10,
               headlength=1e-10, headaxislength=1e-10)

    plt.xlabel(xlab)
    plt.ylabel(ylab)

    return


def plotImage(image=None, arcsec=False, au=False, log=False, dpc=None, maxlog=None, saturate=None, bunit='norm',
              ifreq=0, cmask_rad=None, interpolation='nearest', cmap=plt.cm.gist_gray, stokes='I',
              fig=None, ax=None, projection='polar', deg=True, rmax=None, rlog=True, **kwargs):
    """Plots a radmc3d image.


    Parameters
    ----------
    image         : radmc3dImage
                    A radmc3dImage class returned by readimage

    arcsec        : bool
                    If True image axis will have the unit arcsec (NOTE: dpc keyword should also be set!)

    au            : bool
                    If True image axis will have the unit au

    log           : bool
                    If True image scale will be logarithmic, otherwise linear

    dpc           : float
                    Distance to the source in parsec (This keywords should be set if arcsec=True, or bunit!='norm')

    maxlog        : float
                    Logarithm of the lowest pixel value to be plotted, lower pixel values will be clippde

    saturate      : float
                    Highest pixel values to be plotted in terms of the peak value, higher pixel values will be clipped

    bunit         : {'norm', 'inu', 'snu', 'jy/beam', 'jy/pixel'}
                    Unit of the image, ('norm' - Inu/max(Inu), 'inu' - Inu, 'snu' - Jy/pixel, 'jy/pixel' - Jy/pixel,
                    'jy/beam' - Jy/beam), default is 'norm'. The 'snu' keyword value is kept for backward compatibility
                    as it is fully equivalent with the 'jy/pixel' keyword value.

    ifreq         : int
                    If the image file/array consists of multiple frequencies/wavelengths ifreq denotes the index
                    of the frequency/wavelength in the image array to be plotted

    cmask_rad     : float
                    Simulates coronographyic mask : sets the image values to zero within this radius of the image center
                    The unit is the same as the image axis (au, arcsec, cm)
                    NOTE: this works only on the plot, the image array is not changed (for that used the cmask()
                    function)

    interpolation : str
                    interpolation keyword for imshow (e.g. 'nearest', 'bilinear', 'bicubic')

    cmap          : matplotlib color map

    stokes        : {'I', 'Q', 'U', 'V', 'PI', 'P'}
                   What to plot for full stokes images, Stokes I/Q/U/V,
                   PI  - polarised intensity (PI = sqrt(Q^2 + U^2 + V^2))
                   P   - polarisation fraction (i.e. sqrt(Q^2 + U^2 + V^2) / I)
                   PIL - polarised intensity (PI = sqrt(Q^2 + U^2))
                   PL  - fraction of linear polarisation (i.e. sqrt(Q^2 + U^2) / I)


    fig           : matplotlig.figure.Figure, optional
                   An instance of a matplotlib Figure. If not provided a new Figure will be generated. If provided
                   plotImage will add a single Axes to the Figure, using add_subplots() with the appropriate projection.
                   If the desired plot is to be made for a multi-panel plot, the appropriate Axes instance can be
                   passed to the ax keyword. This keyword is only used for circular images.

    ax            : matplotlib.axes.Axes, optional
                   An instance of a matplotlib Axes to draw the plot on. Note, that the projection of the axes should
                   be the same as the projection keyword passed to plotImage. This keyword is only used for circular
                   images.

    projection    : {'polar', 'cartesian'}
                   Projection of the plot. For cartesian plots a rectangular plot will be drawn, with the horizontal
                   axis being the azimuth angle, and the vertical axis the radius.

    deg           : bool
                   If True the unit of the azimuthal coordinates will degree, if False it will be radian. Used only for
                   circular images and for cartesian projection.

    rmax          : float
                   Maximum value of the radial coordinate for polar projection. Used only for circular images.

    rlog          : bool
                   If True the radial coordiante axis will be set to logarithmic for cartesian projection. Used only
                   for circular images.


    Example
    -------

    result = plotImage(image='image.out', arcsec=True, au=False, log=True, dpc=140, maxlog=-6.,
             saturate=0.1, bunit='Jy')
    """

    if isinstance(image, radmc3dImage):

        # Check whether or not we need to mask the image

        dum_image = copy.deepcopy(image)
        if dum_image.stokes:
            if stokes.strip().upper() == 'I':
                if dum_image.nwav == 1:
                    dum_image.image = image.image[:, :, 0]
                else:
                    dum_image.image = image.image[:, :, 0, :]

            if stokes.strip().upper() == 'Q':
                if dum_image.nwav == 1:
                    dum_image.image = image.image[:, :, 1]
                else:
                    dum_image.image = image.image[:, :, 1, :]

            if stokes.strip().upper() == 'U':
                if dum_image.nwav == 1:
                    dum_image.image = image.image[:, :, 2]
                else:
                    dum_image.image = image.image[:, :, 2, :]

            if stokes.strip().upper() == 'V':
                if dum_image.nwav == 1:
                    dum_image.image = image.image[:, :, 3]
                else:
                    dum_image.image = image.image[:, :, 3, :]

            if stokes.strip().upper() == 'PI':
                if dum_image.nwav == 1:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1]**2 + image.image[:, :, 2]**2 + image.image[:, :, 3]**2)
                else:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1, :]**2 + image.image[:, :, 2, :]**2 + image.image[:, :, 3, :]**2)

            if stokes.strip().upper() == 'P':
                if dum_image.nwav == 1:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1]**2 + image.image[:, :, 2]**2 + image.image[:, :, 3]**2) / \
                                      image.image[:, :, 0]

                else:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1, :]**2 + image.image[:, :, 2, :]**2 + image.image[:, :, 3, :]**2) / \
                                      image.image[:, :, 0, :]

            if stokes.strip().upper() == 'PIL':
                if dum_image.nwav == 1:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1]**2 + image.image[:, :, 2]**2)
                else:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1, :]**2 + image.image[:, :, 2, :]**2)

            if stokes.strip().upper() == 'PL':
                if dum_image.nwav == 1:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1]**2 + image.image[:, :, 2]**2) / image.image[:, :, 0]

                else:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1, :]**2 + image.image[:, :, 2, :]**2) / image.image[:, :, 0, :]

        if cmask_rad is not None:
            dum_image = cmask(dum_image, rad=cmask_rad, au=au, arcsec=arcsec, dpc=dpc)
        else:
            dum_image = dum_image

        if ifreq is None:
            ifreq = 0
        data = np.squeeze(dum_image.image[:, ::-1, ifreq].T)

        norm = data.max()
        if bunit == 'norm':
            data = data / norm

        clipnorm = data.max()
        # Check if the data should be plotted on a log scale
        if log:
            clipmin = np.log10(data[data > 0.].min())
            data = np.log10(data.clip(1e-90))

            # Clipping the data
            if maxlog is not None:
                clipmin = -maxlog + np.log10(clipnorm)
        else:
            clipmin = data.min()

        if saturate is not None:
            if saturate > 1.:
                saturate = 1.0
            if log:
                clipmax = np.log10(saturate) + np.log10(clipnorm)
            else:
                clipmax = clipnorm * saturate
        else:
            clipmax = clipnorm

        data = data.clip(clipmin, clipmax)

        # Select the unit of the data

        if bunit.lower() == 'norm':
            if log:
                cb_label = 'log(I' + r'$_\nu$' + '/max(I' + r'$_\nu$' + '))'
            else:
                cb_label = 'I' + r'$_\nu$' + '/max(I' + r'$_\nu$' + ')'
        elif bunit.lower() == 'inu':
            if log:
                cb_label = 'log(I' + r'$_\nu$' + ' [erg/s/cm/cm/Hz/ster])'
            else:
                cb_label = 'I' + r'$_\nu$' + ' [erg/s/cm/cm/Hz/ster]'

        elif (bunit.lower() == 'snu') | (bunit.lower() == 'jy/pixel'):
            if dpc is None:
                msg = 'Unknown dpc. If Jy/pixel is selected for the image unit the dpc keyword should also be set'
                raise ValueError(msg)
            else:
                if log:
                    data = data + np.log10(image.sizepix_x * image.sizepix_y / (dpc * nc.pc)**2. * 1e23)
                    cb_label = 'log(S' + r'$_\nu$' + '[Jy/pixel])'
                else:
                    data = data * (image.sizepix_x * image.sizepix_y / (dpc * nc.pc)**2. * 1e23)
                    cb_label = 'S' + r'$_\nu$' + ' [Jy/pixel]'

        elif bunit.lower() == 'jy/beam':
            if len(image.fwhm) == 0:
                msg = 'The image does not appear to be convolved with a Gaussain (fwhm data attribute is empty). ' \
                      'The intensity unit can only be converted to Jy/beam if the convolving beam size is known'
                raise ValueError(msg)

            pixel_area = (image.sizepix_x * image.sizepix_y)/(dpc * nc.pc)**2 * (180./np.pi*3600.)**2
            beam_area = image.fwhm[0] * image.fwhm[1] * np.pi / 4. / np.log(2.0)

            if log:
                # Convert data to Jy/pixel
                data += np.log10((image.sizepix_x * image.sizepix_y / (dpc * nc.pc)**2. * 1e23))
                # Convert data to Jy/beam
                data += np.log10(beam_area / pixel_area)

                cb_label = 'log(S' + r'$_\nu$' + '[Jy/beam])'
            else:
                # Convert data to Jy/pixel
                data *= (image.sizepix_x * image.sizepix_y / (dpc * nc.pc)**2. * 1e23)
                # Convert data to Jy/beam
                data *= beam_area / pixel_area
                cb_label = 'S' + r'$_\nu$' + ' [Jy/beam]'

        else:
            msg = 'Unknown bunit: ' + bunit + ' Allowed values are "norm", "inu", "snu"'
            raise ValueError(msg)

        # Select the coordinates of the data
        if au:
            x = image.x / nc.au
            y = image.y / nc.au
            xlab = 'X [au]'
            ylab = 'Y [au]'
        elif arcsec:
            x = image.x / nc.au / dpc
            y = image.y / nc.au / dpc
            xlab = 'RA offset ["]'
            ylab = 'DEC offset ["]'
        else:
            x = image.x
            y = image.y
            xlab = 'X [cm]'
            ylab = 'Y [cm]'

        ext = (x[0], x[image.nx - 1], y[0], y[image.ny - 1])

        # Now finally put everything together and plot the data
        plt.delaxes()
        plt.delaxes()

        implot = plt.imshow(data, extent=ext, cmap=cmap, interpolation=interpolation, **kwargs)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(r'$\lambda$=' + ("%.5g" % image.wav[ifreq]) + r'$\mu$m')
        cbar = plt.colorbar(implot)
        cbar.set_label(cb_label)
        plt.show()

    elif isinstance(image, radmc3dCircimage):

        # Check whether or not we need to mask the image

        dum_image = copy.deepcopy(image)
        if dum_image.stokes:
            if stokes.strip().upper() == 'I':
                if dum_image.nwav == 1:
                    dum_image.image = image.image[:, :, 0]
                else:
                    dum_image.image = image.image[:, :, 0, :]

            if stokes.strip().upper() == 'Q':
                if dum_image.nwav == 1:
                    dum_image.image = image.image[:, :, 1]
                else:
                    dum_image.image = image.image[:, :, 1, :]

            if stokes.strip().upper() == 'U':
                if dum_image.nwav == 1:
                    dum_image.image = image.image[:, :, 2]
                else:
                    dum_image.image = image.image[:, :, 2, :]

            if stokes.strip().upper() == 'V':
                if dum_image.nwav == 1:
                    dum_image.image = image.image[:, :, 3]
                else:
                    dum_image.image = image.image[:, :, 3, :]

            if stokes.strip().upper() == 'PI':
                if dum_image.nwav == 1:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1]**2 + image.image[:, :, 2]**2 + image.image[:, :, 3]**2)
                else:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1, :]**2 + image.image[:, :, 2, :]**2 + image.image[:, :, 3, :]**2)

            if stokes.strip().upper() == 'P':
                if dum_image.nwav == 1:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1]**2 + image.image[:, :, 2]**2 + image.image[:, :, 3]**2) / \
                                      image.image[:, :, 0]

                else:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1, :]**2 + image.image[:, :, 2, :]**2 + image.image[:, :, 3, :]**2) / \
                                      image.image[:, :, 0, :]

            if stokes.strip().upper() == 'PIL':
                if dum_image.nwav == 1:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1]**2 + image.image[:, :, 2]**2)
                else:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1, :]**2 + image.image[:, :, 2, :]**2)

            if stokes.strip().upper() == 'PL':
                if dum_image.nwav == 1:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1]**2 + image.image[:, :, 2]**2) / image.image[:, :, 0]

                else:
                    dum_image.image = np.sqrt(
                        image.image[:, :, 1, :]**2 + image.image[:, :, 2, :]**2) / image.image[:, :, 0, :]
        else:
            dum_image.image = image.image[:, :, 0, :]

        if cmask_rad is not None:
            cmask_rad_cm = cmask_rad
            if au:
                cmask_rad_cm = cmask_rad * 1.496e13
            if arcsec:
                cmask_rad_cm = cmask_rad * dpc * nc.pc

            ii = (dum_image.ri <= cmask_rad_cm)

            if len(dum_image.image.shape) == 3:
                dum_image.image[ii, :, :] = 0.

            elif len(dum_image.image.shape) == 4:
                dum_image.image[ii, :, :, :] = 0.

        if ifreq is None:
            ifreq = 0

        if bunit == 'snu':
            if dpc is None:
                msg = 'Unknown bunit: ' + bunit + ' Allowed values are "norm", "inu", "snu"'
                raise ValueError(msg)
            else:
                psize = image.getPixelSize()
                conv = psize / (dpc * nc.pc**2) * 1e23
                data = dum_image.image[:, :, ifreq] * conv
        else:
            data = dum_image.image[:, :, ifreq]

        norm = data.max()
        if bunit == 'norm':
            data = data / norm

        clipnorm = data.max()
        # Check if the data should be plotted on a log scale
        if log:
            clipmin = np.log10(data[data > 0.].min())
            data = np.log10(data.clip(1e-90))

            # Clipping the data
            if maxlog is not None:
                clipmin = -maxlog + np.log10(clipnorm)
        else:
            clipmin = data.min()

        if saturate is not None:
            if saturate > 1.:
                saturate = 1.0
            if log:
                clipmax = np.log10(saturate) + np.log10(clipnorm)
            else:
                clipmax = clipnorm * saturate
        else:
            clipmax = clipnorm

        data = data.clip(clipmin, clipmax)

        # Select the unit of the data

        if bunit == 'norm':
            if log:
                cb_label = 'log(I' + r'$_\nu$' + '/max(I' + r'$_\nu$' + '))'
            else:
                cb_label = 'I' + r'$_\nu$' + '/max(I' + r'$_\nu$' + ')'
        elif bunit == 'inu':
            if log:
                cb_label = 'log(I' + r'$_\nu$' + ' [erg/s/cm/cm/Hz/ster])'
            else:
                cb_label = 'I' + r'$_\nu$' + ' [erg/s/cm/cm/Hz/ster]'
        elif bunit == 'snu':
            if log:
                cb_label = 'log(S' + r'$_\nu$' + '[Jy/pixel])'
            else:
                cb_label = 'S' + r'$_\nu$' + ' [Jy/pixel]'

        else:
            msg = 'Unknown bunit: ' + bunit + ' Allowed values are "norm", "inu", "snu"'
            raise ValueError(msg)

        # Select the coordinates of the data
        x = image.phii
        xlab = 'Azimuth angle [rad]'
        if projection == 'cartesian':
            if deg:
                x = image.phii / np.pi * 180.
                xlab = 'Azimuth angle [deg]'

        if au:
            y = image.rc / 1.496e13
            ylab = 'R [au]'
        elif arcsec:
            y = image.rc / 1.496e13 / dpc
            ylab = 'R [arcsec]'
        else:
            y = image.rc
            ylab = 'R [cm]'

        if rmax is None:
            rmax = x.max()

        if fig is None:
            fig = plt.figure()

        if projection == 'polar':
            ax = fig.add_subplot(111, projection='polar')
            implot = plt.pcolormesh(x, y, data, **kwargs)
            ax.set_rmax(rmax)

        elif projection == 'cartesian':
            if ax is not None:
                plt.sca(ax)

            implot = plt.pcolormesh(x, y, data, **kwargs)
            plt.xlim(x[0], x[-1])
            plt.ylim(y[0], y[-1])
            if rlog:
                plt.yscale('log')
        else:
            msg = 'Unknown projection. Accepted values for projection keyword are "polar" or "cartesian".'
            raise ValueError(msg)

        plt.xlabel(xlab)
        plt.ylabel(ylab)
        cbar = plt.colorbar()
        cbar.set_label(cb_label)
        plt.show()

    return {'implot': implot, 'cbar': cbar}


def makeImage(npix=None, incl=None, wav=None, sizeau=None, phi=None, posang=None, pointau=None,
              fluxcons=True, nostar=False, noscat=False, tracetau=False,
              widthkms=None, linenlam=None, vkms=None, iline=None,
              lambdarange=None, nlam=None, stokes=False, binary=False, setthreads=None):
    """Calculates a rectangular image with RADMC-3D

    Parameters
    ----------

    npix        : int
                  Number of pixels on the rectangular images

    sizeau      : float
                  Diameter of the image in au

    incl        : float
                  Inclination angle of the source

    wav         : float
                  Wavelength of the image in micron

    phi         : float, optional
                  Azimuthal rotation angle of the source in the model space

    posang      : float, optional
                  Position angle of the source in the image plane

    pointau     : Float, optional
                  Three elements list of the cartesian coordinates of the image center

    widthkms    : float, optional
                  Width of the frequency axis of the channel maps

    linenlam    : int, optional
                  Number of wavelengths to calculate images at

    vkms        : float, optional
                  A single velocity value at which a channel map should be calculated

    iline       : int, optional
                  Line transition index

    lambdarange : list, optional
                  Two element list with the wavelenght boundaries between which
                  multiwavelength images should be calculated

    nlam        : int, optional
                  Number of wavelengths to be calculated in lambdarange

    fluxcons    : bool, optional
                  This should not even be a keyword argument, it ensures flux conservation
                  (adaptive subpixeling) in the rectangular images

    nostar      : bool, optional
                  If True the calculated images will not contain stellar emission

    noscat      : bool, optional
                  If True, scattered emission will be neglected in the source function, however,
                   extinction will contain scattering if kappa_scat is not zero.

    tracetau    : bool, optional
                  If True, then instead of the intensity I_nu, the optical depth along the
                  ray tau_nu is computed, and read in as if it is an image.

    stokes      : bool, optional
                  If True, images in all four stokes parameters (IQUV) will be calculated, if
                  False only the intensity will be calculated

    binary      : bool, optional
                  If True the output image will be written in a C-style binary format, if False
                  the image format will be ASCII

    setthreads  : integer, optional
                  If set, then do the Monte Carlo with OpenMP parallellization, with
                  the integer setthreads being the number of parallel threads.

    Example
    -------

    makeImage(npix=100, incl=60.0, wav=10.0, sizeau=300., phi=0., posang=15.,
        pointau=[0., 0.,0.], fluxcons=True, nostar=False, noscat=False)

    """
    #
    # The basic keywords that should be set
    #
    if npix is None:
        msg = 'Unkonwn npix. Number of pixels must be set.'
        raise ValueError(msg)

    if incl is None:
        msg = 'Unkonwn incl. Inclination angle must be set.'
        raise ValueError(msg)

    if wav is None:
        if (lambdarange is None) & (nlam is None):
            if vkms is None:
                if (widthkms is None) & (linenlam is None):
                    msg = 'Neither wavelength nor velocity is specified at which the image should be calculated'
                    raise ValueError(msg)
                else:
                    if iline is None:
                        msg = 'Unknown iline. widthkms, linenlam keywords are set indicating that a line '\
                              + 'channel map should be calculated, but the iline keyword is not specified'
                        raise ValueError(msg)
            else:
                if iline is None:
                    msg = 'Unknown iline. vkms keyword is set indicating that a line channel map should be'\
                          + 'calculated, but the iline keyword is not specified'
                    raise ValueError(msg)
    else:
        if lambdarange is not None:
            msg = 'Either lambdarange or wav should be set but not both'
            raise ValueError(msg)

    if lambdarange is not None:
        if len(lambdarange) != 2:
            msg = 'lambdarange must have two and only two elements'
            raise ValueError(msg)

    #
    # Kees' fix for the case when a locally compiled radmc3d exists in the current directory
    #
    com = ''
    if os.path.isfile('radmc3d'):
        com = com + './'
    com = com + 'radmc3d image'

    # com = 'radmc3d image'
    com = com + ' npix ' + str(int(npix))
    com = com + ' incl ' + str(incl)
    if sizeau is not None:
        com = com + ' sizeau ' + str(sizeau)

    if wav is not None:
        com = com + ' lambda ' + str(wav)
    elif (lambdarange is not None) & (nlam is not None):
        com = com + ' lambdarange ' + str(lambdarange[0]) + ' ' + str(lambdarange[1]) + ' nlam ' + str(int(nlam))
    elif vkms is not None:
        com = com + ' vkms ' + str(vkms)
    elif (widthkms is not None) & (linenlam is not None):
        com = com + ' widthkms ' + str(widthkms) + ' linenlam ' + str(linenlam)

    #
    # Now add additional optional keywords/arguments
    #
    if phi is not None:
        com = com + ' phi ' + str(phi)

    if posang is not None:
        com = com + ' posang ' + str(posang)

    if pointau is not None:
        if len(pointau) != 3:
            msg = ' pointau should be a list of 3 elements corresponding to the  cartesian coordinates of the ' \
                  + 'image center'
            raise ValueError(msg)
        else:
            com = com + ' pointau ' + str(pointau[0]) + ' ' + str(pointau[1]) + ' ' + str(pointau[2])
    else:
        com = com + ' pointau 0.0  0.0  0.0'

    if fluxcons:
        com = com + ' fluxcons'

    if iline:
        com = com + ' iline ' + ("%d" % iline)

    if stokes:
        com = com + ' stokes'

    if binary:
        com = com + ' imageunform'

    if nostar:
        com = com + ' nostar'

    if noscat:
        com = com + ' noscat'

    if tracetau:
        com = com + ' tracetau'

    if setthreads is not None:
        com = com + ' setthreads {}'.format(setthreads)

    #
    # Print the command
    #
    print("Executing RADMC-3D Command:")
    print(com)

    #
    # Now finally run radmc3d and calculate the image
    #
    # dum = sp.Popen([com], stdout=sp.PIPE, shell=True).wait()
    dum = sp.Popen([com], shell=True).wait()

    return 0


def cmask(im=None, rad=0.0, au=False, arcsec=False, dpc=None):
    """Simulates a coronographic mask.
        Sets the image values to zero within circle of a given radius around the
        image center.

    Parameters
    ----------
    im     : radmc3dImage
            A radmc3dImage class containing the image

    rad    : float
            The raadius of the mask

    au     : bool
            If true the radius is taken to have a unit of au

    arcsec : bool
            If true the radius is taken to have a unit of arcsec (dpc
            should also be set)

    dpc    : float
            Distance of the source (required if arcsec = True)

    NOTE if arcsec=False and au=False rad is taken to have a unit of pixel

    Returns
    -------

    Returns a radmc3dImage class containing the masked image
    """

    if au:
        if arcsec:
            msg = ' Either au or arcsec should be set, but not both of them'
            raise ValueError(msg)

        crad = rad * nc.au
    else:
        if arcsec:
            crad = rad * nc.au * dpc
        else:
            crad = rad * im.sizepix_x

    res = copy.deepcopy(im)
    if im.nfreq != 1:
        for ix in range(im.nx):
            r = np.sqrt(im.y**2 + im.x[ix]**2)
            ii = r <= crad
            res.image[ix, ii, :] = 0.0
    else:
        for ix in range(im.nx):
            r = np.sqrt(im.y**2 + im.x[ix]**2)
            ii = r <= crad
            res.image[ix, ii] = 0.0

    return res


class radmc3dCircimage(object):
    """
    RADMC-3D circular image class

    Attributes
    ----------

    image       : ndarray
                  The image as calculated by radmc3d (the values are intensities in erg/s/cm^2/Hz/ster)

    rc          : ndarray
                  Radial cell center coordinate of the image [cm]

    ri          : ndarray
                  Radial cell interface coordinate of the image [cm]

    phic        : ndarray
                  Azimuthal cell center coordinate of the image [rad]

    phii        : ndarray
                  Azimuthal cell interface coordinate of the image [rad]

    nr          : int
                  Number of pixels in the radial direction

    nphi        : int
                  Number of pixels in the azimuthal direction

    nfreq       : int
                  Number of frequencies in the image cube

    freq        : ndarray
                  Frequency grid in the image cube

    nwav        : int
                  Number of wavelengths in the image cube (same as nfreq)

    wav         : ndarray
                  Wavelength grid in the image cube

    filename    : str
                  Name of the file the image data was read from

    stokes      : bool
                  If True the image data contain the full stokes vector (I,Q,U,V)

    """

    def __init__(self):

        self.ri = np.zeros(0, dtype=np.float64)
        self.phii = np.zeros(0, dtype=np.float64)
        self.rc = np.zeros(0, dtype=np.float64)
        self.phic = np.zeros(0, dtype=np.float64)
        self.freq = np.zeros(0, dtype=np.float64)
        self.wav = np.zeros(0, dtype=np.float64)
        self.image = np.zeros((0, 0), dtype=np.float64)
        self.imageJyppix = np.zeros((0, 0), dtype=np.float64)
        self.filename = 'circimage.out'
        self.nphi = 0
        self.nwav = 0
        self.nfreq = 0
        self.stokes = False
        self.npol = 1

        self.nr = 0
        self.nphi = 0

    def getPixelSize(self):
        """
        Calculates the pixel size

        Returns
        -------
        The pixel size in cm^2
        """

        nx, ny = self.image.shape[:2]
        psize = np.zeros([nx, ny], dtype=np.float64)

        x2 = np.pi * (self.ri[1:]**2 - self.ri[:-1]**2)
        dy = self.phii[1:] - self.phii[:-1]

        for ix in range(nx):
            psize[ix, :] = x2[ix] * dy / (2.0 * np.pi)

        return psize

    def readImage(self, filename='circimage.out', old=False):
        """
        Reads a circular image

        Parameters
        ----------

        filename        : str
                          Name of the file to be read.

        old             : bool
                          If True the image format of the old 2D code (radmc) will be used. If False (default) the
                          RADMC-3D format is used.
        """

        self.filename = filename

        if old:
            self.stokes = False
            self.npol = 1

            with open(filename, 'r') as f:

                self.nfreq = int(f.readline())
                self.nwav = self.nfreq

                s = f.readline()
                self.freq = np.zeros(self.nfreq, dtype=np.float64)
                self.wav = np.zeros(self.nfreq, dtype=np.float64)
                for inu in range(self.nfreq):
                    s = f.readline()
                    self.freq[inu] = float(s)
                    self.wav[inu] = nc.cc / self.freq[inu]

                s = f.readline()
                s = f.readline().split()
                self.nr = int(s[0])
                self.nphi = int(s[1])
                self.nfreq = int(s[2])

                s = f.readline()
                self.rc = np.zeros(self.nr, dtype=np.float64)
                for ir in range(self.nr):
                    s = f.readline()
                    self.rc[ir] = float(s)

                s = f.readline()
                self.ri = np.zeros(self.nr + 1, dtype=np.float64)
                for ir in range(self.nr + 1):
                    s = f.readline()
                    self.ri[ir] = float(s)

                s = f.readline()
                self.phic = np.zeros(self.nphi, dtype=np.float64)
                for ip in range(self.nphi):
                    s = f.readline()
                    self.phic[ip] = float(s)

                s = f.readline()
                self.phii = np.zeros(self.nphi + 1, dtype=np.float64)
                for ip in range(self.nphi + 1):
                    s = f.readline()
                    self.phii[ip] = float(s)

                s = f.readline()
                self.image = np.zeros((self.nr, self.nphi, self.npol, self.nfreq), dtype=np.float64)

                for inu in range(self.nfreq):
                    for ir in range(self.nr):
                        s = f.readline().split()
                        for iphi in range(self.nphi):
                            self.image[ir, :, 0, inu] = float(s[iphi])

                psize = self.getPixelSize()
                conv = psize / (nc.pc**2) * 1e23
                self.imageJyppix = np.zeros((self.nr, self.nphi, self.npol, self.nfreq), dtype=np.float64)
                for inu in range(self.nfreq):
                    self.imageJyppix[:, :, 0, inu] = self.image[:, :, 0, inu] * conv
        else:

            with open(filename, 'r') as f:

                iformat = int(f.readline())
                if iformat == 1:
                    stokes = False
                    self.npol = 1
                elif iformat == 3:
                    stokes = True
                    self.npol = 4

                s = f.readline().split()
                self.nr = int(s[0])
                self.nphi = int(s[1])
                self.nfreq = int(f.readline())

                s = f.readline()
                self.ri = np.zeros(self.nr + 2, dtype=np.float64)
                for ir in range(self.nr + 2):
                    s = f.readline()
                    self.ri[ir] = float(s)

                s = f.readline()
                self.rc = np.zeros(self.nr + 1, dtype=np.float64)
                for ir in range(self.nr + 1):
                    s = f.readline()
                    self.rc[ir] = float(s)

                s = f.readline()
                self.phii = np.zeros(self.nphi + 1, dtype=np.float64)
                for ip in range(self.nphi + 1):
                    s = f.readline()
                    self.phii[ip] = float(s)

                s = f.readline()
                self.phic = np.zeros(self.nphi, dtype=np.float64)
                for ip in range(self.nphi):
                    s = f.readline()
                    self.phic[ip] = float(s)

                s = f.readline()
                self.freq = np.zeros(self.nfreq, dtype=np.float64)
                self.wav = np.zeros(self.nfreq, dtype=np.float64)
                for inu in range(self.nfreq):
                    s = f.readline()
                    print(inu, s)
                    self.freq[inu] = float(s)
                    self.wav[inu] = nc.cc / self.freq[inu]

                s = f.readline()
                self.image = np.zeros((self.nr + 1, self.nphi, self.npol, self.nfreq), dtype=np.float64)

                for inu in range(self.nfreq):
                    for iphi in range(self.nphi):
                        for ir in range(self.nr + 1):
                            s = f.readline().split()
                            self.image[ir, iphi, 0, inu] = float(s[0])
                            if self.npol > 1:
                                self.image[ir, iphi, 1, inu] = float(s[1])
                                self.image[ir, iphi, 2, inu] = float(s[2])
                                self.image[ir, iphi, 3, inu] = float(s[3])

            psize = self.getPixelSize()
            conv = psize / (nc.pc**2) * 1e23
            self.imageJyppix = np.zeros((self.nr + 1, self.nphi, self.npol, self.nfreq), dtype=np.float64)
            for ipol in range(self.npol):
                for inu in range(self.nfreq):
                    self.imageJyppix[:, :, ipol, inu] = self.image[:, :, ipol, inu] * conv


def readcircimage(filename='circimage.out', old=False):
    """
    A convenience function to read circular images

    Parameters
    ----------

    filename        : str
                      Name of the file to be read.

    old             : bool
                      If True the image format of the old 2D code (radmc) will be used. If False (default) the
                      RADMC-3D format is used.
    """
    dum = radmc3dCircimage()
    dum.readImage(filename=filename, old=old)
    return dum
