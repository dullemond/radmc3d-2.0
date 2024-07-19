"""This module contains a class for handling regular wavelength and spatial grids
"""

from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
import traceback

try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use the python module of RADMC-3D you need to install Numpy')
    print(traceback.format_exc())

from . import natconst as nc


class radmc3dGrid(object):
    """ Class for spatial and frequency grid used by RADMC-3D.

    Attributes
    ----------

    act_dim    : ndarray
                A three element vector the i-th element is 1 if the i-th dimension is active,
                otherwize the i-th element is zero

    crd_sys    : {'sph', 'cyl', 'car'}
                coordinate system of the spatial grid

    nx         : int
                Number of grid points in the x (cartesian) / r (cylindrical) / r (spherical) dimension

    ny         : int
                Number of grid points in the y (cartesian) / theta (cylindrical) / theta (spherical) dimension

    nz         : int
                Number of grid points in the z (cartesian) / z (cylindrical) / phi (spherical) dimension

    nxi        : int
                Number of cell interfaces in the x (cartesian) / r (cylindrical) / r (spherical) dimension

    nyi        : int
                Number of cell interfaces in the y (cartesian) / theta (cylindrical) / theta (spherical) dimension

    nzi        : int
                Number of cell interfaces in the z (cartesian) / z (cylindrical) / phi (spherical) dimension

    nwav       : int
                Number of wavelengths in the wavelength grid

    nfreq      : int
                Number of frequencies in the grid (equal to nwav)

    x          : ndarray
                Cell centered x (cartesian) / r (cylindrical) / r (spherical)  grid points

    y          : ndarray
                Cell centered y (cartesian) / theta (cylindrical) / theta (spherical)  grid points

    z          : ndarray
                Cell centered z (cartesian) / z (cylindrical) / phi (spherical)  grid points

    xi         : ndarray
                Cell interfaces in the x (cartesian) / r (cylindrical) / r (spherical)  dimension

    yi         : ndarray
                Cell interfaces in the y (cartesian) / theta (cylindrical) / theta (spherical)  dimension

    zi         : ndarray
                Cell interfaces in the z (cartesian) / z (cylindrical) / phi (spherical)  dimension

    wav        : ndarray
                Wavelengh  grid

    freq       : ndarray
                Frequency  grid


    """

    # --------------------------------------------------------------------------------------------------

    def __init__(self):

        self.crd_sys = 'sph'
        self.act_dim = [1, 1, 1]
        self.grid_style = 0
        self.octree = False
        self.nx = 0
        self.ny = 0
        self.nz = 0
        self.nxi = 0
        self.nyi = 0
        self.nzi = 0
        self.x = np.zeros(0, dtype=np.float64)
        self.y = np.zeros(0, dtype=np.float64)
        self.z = np.zeros(0, dtype=np.float64)
        self.xi = np.zeros(0, dtype=np.float64)
        self.yi = np.zeros(0, dtype=np.float64)
        self.zi = np.zeros(0, dtype=np.float64)

        self.nwav = 0
        self.nfreq = 0
        self.wav = np.zeros(0, dtype=np.float64)
        self.freq = np.zeros(0, dtype=np.float64)

    # --------------------------------------------------------------------------------------------------
    def makeWavelengthGrid(self, wbound=None, nw=None, ppar=None):
        """Creates the wavelength/frequency grid.

        Parameters
        ----------

        wbound : list
                 Contains the wavelength boundaries of the wavelength grid (should contain at least two elements)

        nw     : list
                 Contains len(wbound)-1 elements containing the number of wavelengths between the bounds
                 set by wbound

        ppar   : dictionary, optional
                 Contains all input parameters with the parameter names as keys
        """
        if ppar:
            if wbound is None:
                wbound = ppar['wbound']
            if nw is None:
                nw = ppar['nw']
        else:
            if wbound is None:
                raise ValueError(
                    'Unknown wbound. Without the grid boundaries the wavelength grid cannot be generated')
            if nw is None:
                raise ValueError('Unknown nw. Without the number of grid points the wavelength grid cannot be '
                                 + 'generated')

        self.nwav = nw[0]
        self.wav = wbound[0] * (wbound[1] / wbound[0]) ** (np.arange(nw[0], dtype=np.float64) / nw[0])

        for ipart in range(1, len(nw) - 1):
            dum = wbound[ipart] * (wbound[ipart + 1]
                                   / wbound[ipart]) ** (np.arange(nw[ipart], dtype=np.float64) / nw[ipart])
            self.wav = np.append(self.wav, dum)

        ipart = len(nw) - 1
        dum = wbound[ipart] * (wbound[ipart + 1]
                               / wbound[ipart]) ** (np.arange(nw[ipart], dtype=np.float64) / (nw[ipart] - 1.))
        self.wav = np.append(self.wav, dum)
        self.nwav = self.wav.shape[0]
        self.freq = nc.cc / self.wav * 1e4
        self.nfreq = self.nwav

    # --------------------------------------------------------------------------------------------------
    def writeWavelengthGrid(self, fname='', old=False):
        """Wriites the wavelength grid to a file (e.g. wavelength_micron.inp).

        Parameters
        ----------

        fname  : str, optional
                 File name into which the wavelength grid should be written. If omitted 'wavelength_micron.inp'
                 will be used

        old    : bool, optional
                 If set to True the file format of the previous, 2D version of radmc will be used
        """

        if not old:
            if fname == '':
                fname = 'wavelength_micron.inp'

            print('Writing ' + fname)
            with open(fname, 'w') as wfile:
                wfile.write('%d\n' % self.nwav)
                for ilam in range(self.nwav):
                    wfile.write('%.9e\n' % self.wav[ilam])
        else:
            if fname == '':
                fname = 'frequency.inp'

            with open(fname, 'w') as wfile:

                print('Writing ' + fname)
                wfile.write("%d\n" % self.nfreq)
                wfile.write(" \n")
                #
                # Reverse the order of the frequency grid as it is ordered in frequency in radmc
                #
                freq = self.freq[::-1]
                for i in range(self.nfreq):
                    wfile.write("%.7e\n" % freq[i])

                    # --------------------------------------------------------------------------------------------------

    def makeSpatialGrid(self, crd_sys=None, xbound=None, ybound=None, zbound=None, nxi=None, nyi=None, nzi=None,
                        ppar=None):
        """Calculates the spatial grid.

        Parameters
        ----------

        crd_sys : {'sph','car'}
                    Coordinate system of the spatial grid

        xbound  : list
                    (with at least two elements) of boundaries for the grid along the first dimension

        ybound  : list
                    (with at least two elements) of boundaries for the grid along the second dimension

        zbound  : list
                    (with at least two elements) of boundaries for the grid along the third dimension

        nxi     : int
                    Number of grid points along the first dimension. List with len(xbound)-1 elements with
                    nxi[i] being the number of grid points between xbound[i] and xbound[i+1]

        nyi     : int
                    Same as nxi but for the second dimension

        nzi     : int
                    Same as nxi but for the third dimension

        ppar    : Dictionary containing all input parameters of the model (from the problem_params.inp file)
                   if ppar is set all keyword arguments that are not set will be taken from this dictionary
        """

        self.act_dim = [1, 1, 1]
        if ppar:
            if not crd_sys:
                crd_sys = ppar['crd_sys']
            self.crd_sys = crd_sys

            if not xbound:
                if 'xbound' in ppar:
                    xbound = ppar['xbound']
                else:
                    print(' No boundary for the first dimension is given, first dimension is deactivated.')
                    self.act_dim[0] = 0
            if not nxi:
                if 'nx' in ppar:
                    if not isinstance(ppar['nx'], list):
                        ppar['nx'] = [ppar['nx']]
                    nxi = [i + 1 for i in ppar['nx']]
                    if ppar['nx'][0] == 0:
                        self.act_dim[0] = 0
                else:
                    self.act_dim[0] = 0

            if not ybound:
                if 'ybound' in ppar:
                    ybound = ppar['ybound']
                else:
                    print(' No boundary for the second dimension is given, second dimension is deactivated.')
                    self.act_dim[1] = 0
            if not nyi:
                if 'ny' in ppar:
                    if not isinstance(ppar['ny'], list):
                        nyi = [ppar['ny'] + 1]
                        ppar['ny'] = [ppar['ny']]
                    else:
                        ppar['ny'] = ppar['ny']
                        nyi = [i + 1 for i in ppar['ny']]

                    if ppar['ny'][0] == 0:
                        self.act_dim[1] = 0
                else:
                    self.act_dim[1] = 0

            if not zbound:
                if 'zbound' in ppar:
                    zbound = ppar['zbound']
                else:
                    print(' No boundary for the third dimension is given, third dimension is deactivated.')
                    self.act_dim[2] = 0
            if not nzi:
                if 'nz' in ppar:
                    if not isinstance(ppar['nz'], list):
                        ppar['nz'] = [ppar['nz']]
                    if ppar['nz'][0] > 0:
                        nzi = [i + 1 for i in ppar['nz']]
                    if ppar['nz'][0] == 0:
                        self.act_dim[2] = 0
                else:
                    self.act_dim[2] = 0
                    nzi = [0]
        #
        # Type checking
        #
        if not isinstance(nxi, list):
            nxi = [nxi]
        if not isinstance(nyi, list):
            nyi = [nyi]
        if not isinstance(nzi, list):
            nzi = [nzi]

        if crd_sys == 'car':
            #
            # First check whether the grid boundaries are specified
            #
            if xbound is None:
                raise ValueError('Unknown xbound. Boundaries for the cartesian x-axis are not specified. '
                                 + 'Without the boundaries the grid cannot be generated')

            if ybound is None:
                raise ValueError('Unknown ybound. Boundaries for the cartesian y-axis are not specified. '
                                 + 'Without the boundaries the grid cannot be generated')
            if zbound is None:
                raise ValueError('Unknown zbound. Boundaries for the cartesian z-axis are not specified. '
                                 + 'Without the boundaries the grid cannot be generated')

            if nxi is None:
                raise ValueError('Unknown nxi. Number of grid points for the cartesian x-axis is not specified. '
                                 + 'The grid cannot be generated')

            if nyi is None:
                raise ValueError('Unknown nyi. Number of grid points for the cartesian y-axis is not specified. '
                                 + 'The grid cannot be generated')

            if nzi is None:
                raise ValueError('Unknown nzi. Number of grid points for the cartesian z-axis is not specified. '
                                 + 'The grid cannot be generated')

            #
            # Create the x-axis
            #
            if len(nxi) > 1:
                self.nxi = sum(nxi)
                self.nx = self.nxi - 1
                self.xi = xbound[0] + (xbound[1] - xbound[0]) * (
                    np.arange(nxi[0], dtype=np.float64) / float(nxi[0]))
                for ipart in range(1, len(nxi) - 1):
                    dum = xbound[ipart] + (xbound[ipart + 1] - xbound[ipart]) \
                                          * (np.arange(nxi[ipart], dtype=np.float64) / float(nxi[ipart]))
                    self.xi = np.append(self.xi, dum)

                ipart = len(nxi) - 1
                dum = xbound[ipart] + (xbound[ipart + 1] - xbound[ipart]) * (
                    np.arange(nxi[ipart], dtype=np.float64) / float(nxi[ipart] - 1))
                self.xi = np.append(self.xi, dum)
                self.x = 0.5 * (self.xi[0:self.nx] + self.xi[1:self.nx + 1])
            else:
                if self.act_dim[0] == 1:
                    self.nxi = nxi[0]
                    self.xi = xbound[0] + (xbound[1] - xbound[0]) * (
                        np.arange(self.nxi, dtype=np.float64) / float(self.nxi - 1.))
                    self.nx = self.nxi - 1
                    self.x = 0.5 * (self.xi[0:self.nx] + self.xi[1:self.nx + 1])
                else:
                    self.x = [0.]
                    self.xi = [0., 0., ]
                    self.nx = 1
                    self.nxi = 2

                    #
                    # Create the y-ayis
                    #
            if len(nyi) > 1:
                self.nyi = sum(nyi)
                self.ny = self.nyi - 1
                self.yi = ybound[0] + (ybound[1] - ybound[0]) * (
                    np.arange(nyi[0], dtype=np.float64) / float(nyi[0]))
                for ipart in range(1, len(nyi) - 1):
                    dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                        np.arange(nyi[ipart], dtype=np.float64) / float(nyi[ipart]))
                    self.yi = np.append(self.yi, dum)

                ipart = len(nyi) - 1
                dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                    np.arange(nyi[ipart], dtype=np.float64) / float(nyi[ipart] - 1))
                self.yi = np.append(self.yi, dum)
                self.y = 0.5 * (self.yi[0:self.ny] + self.yi[1:self.ny + 1])
            else:
                if self.act_dim[0] == 1:
                    self.nyi = nyi[0]
                    self.yi = ybound[0] + (ybound[1] - ybound[0]) * (
                        np.arange(self.nyi, dtype=np.float64) / float(self.nyi - 1.))
                    self.ny = self.nyi - 1
                    self.y = 0.5 * (self.yi[0:self.ny] + self.yi[1:self.ny + 1])
                else:
                    self.y = [0.]
                    self.yi = [0., 0., ]
                    self.ny = 1
                    self.nyi = 2

                    #
                    # Create the z-azis
                    #
            if len(nzi) > 1:
                self.nzi = sum(nzi)
                self.nz = self.nzi - 1
                self.zi = zbound[0] + (zbound[1] - zbound[0]) * (
                    np.arange(nzi[0], dtype=np.float64) / float(nzi[0]))
                for ipart in range(1, len(nzi) - 1):
                    dum = zbound[ipart] + (zbound[ipart + 1] - zbound[ipart]) \
                                          * (np.arange(nzi[ipart], dtype=np.float64) / float(nzi[ipart]))
                    self.zi = np.append(self.zi, dum)

                ipart = len(nzi) - 1
                dum = zbound[ipart] + (zbound[ipart + 1] - zbound[ipart]) * (
                    np.arange(nzi[ipart], dtype=np.float64) / float(nzi[ipart] - 1))
                self.zi = np.append(self.zi, dum)
                self.z = 0.5 * (self.zi[0:self.nz] + self.zi[1:self.nz + 1])
            else:
                if self.act_dim[0] == 1:
                    self.nzi = nzi[0]
                    self.zi = zbound[0] + (zbound[1] - zbound[0]) * (
                        np.arange(self.nzi, dtype=np.float64) / float(self.nzi - 1.))
                    self.nz = self.nzi - 1
                    self.z = 0.5 * (self.zi[0:self.nz] + self.zi[1:self.nz + 1])
                else:
                    self.z = [0.]
                    self.zi = [0., 0.0]
                    self.nz = 1
                    self.nzi = 2

        if crd_sys == 'sph':
            #
            # r->x, theta->y, phi-z
            #
            if xbound is None:
                raise ValueError('Unknown xbound. Boundaries for the spherical radial grid are not specified. '
                                 + 'Without the boundaries the grid cannot be generated')

            if ybound is None:
                print('Unknown ybound. Setting the spherical co-lattitude grid to [0, pi]')
                ybound = [0.0, np.pi]

            if zbound is None:
                print('Unknown zbound. Setting the spherical azimuth angle grid to [0, 2*pi]')
                zbound = [0.0, 2.0 * np.pi]

            if nxi is None:
                raise ValueError(
                    'Unknown nxi. Number of grid points for the spherical radial grid is not specified. '
                    + 'The grid cannot be generated')

            if nyi is None:
                raise ValueError('Unknown nyi. Number of grid points for the spherical co-lattitude grid is not '
                                 + ' specified. The grid cannot be generated')

            if nzi is None:
                raise ValueError('Unknown nzi. Number of grid points for the spherical azimuthal angle grid is not '
                                 + 'specified. The grid cannot be generated')

            #
            # Create the x axis
            #
            if len(nxi) > 1:
                self.nxi = sum(nxi)
                self.nx = self.nxi - 1
                self.xi = xbound[0] * (xbound[1] / xbound[0]) ** (
                    np.arange(nxi[0], dtype=np.float64) / float(nxi[0]))
                for ipart in range(1, len(nxi) - 1):
                    dum = xbound[ipart] * (xbound[ipart + 1] / xbound[ipart]) ** (
                        np.arange(nxi[ipart], dtype=np.float64) / float(nxi[ipart]))
                    self.xi = np.append(self.xi, dum)

                ipart = len(nxi) - 1
                dum = xbound[ipart] * (xbound[ipart + 1] / xbound[ipart]) ** (
                    np.arange(nxi[ipart], dtype=np.float64) / float(nxi[ipart] - 1))
                self.xi = np.append(self.xi, dum)
                self.x = np.sqrt(self.xi[0:self.nx] * self.xi[1:self.nx + 1])
            else:
                if self.act_dim[0] == 1:
                    self.nxi = nxi[0]
                    self.xi = xbound[0] * (xbound[1] / xbound[0]) ** (
                        np.arange(self.nxi, dtype=np.float64) / float(self.nxi - 1.))
                    self.nx = self.nxi - 1
                    self.x = np.sqrt(self.xi[0:self.nx] * self.xi[1:self.nx + 1])
                else:
                    self.x = [0.]
                    self.xi = [0., 0., ]
                    self.nx = 1
                    self.nxi = 2

            # Refinement of the inner edge of the grid
            # This has to be done properly
            if 'xres_nlev' in ppar:
                if ppar['xres_nlev'] > 0:
                    ri_ext = np.array([self.xi[0], self.xi[ppar['xres_nspan']]])
                    for i in range(ppar['xres_nlev']):
                        dum_ri = ri_ext[0] + (ri_ext[1] - ri_ext[0]) * np.arange(ppar['xres_nstep'] + 1,
                                                                                 dtype=np.float64) / float(
                            ppar['xres_nstep'])
                        # print ri_ext[0:2]/au
                        # print dum_ri/au
                        ri_ext_old = np.array(ri_ext)
                        ri_ext = np.array(dum_ri)
                        ri_ext = np.append(ri_ext, ri_ext_old[2:])

                    r_ext = (ri_ext[1:] + ri_ext[:-1]) * 0.5

                    self.xi = np.append(ri_ext, self.xi[ppar['xres_nspan'] + 1:])
                    self.x = np.append(r_ext, self.x[ppar['xres_nspan']:])
                    self.nx = self.x.shape[0]
                    self.nxi = self.xi.shape[0]

                    #
                    # Create the y axis
                    #
            if len(nyi) > 1:

                # Check if we go to the full [0,pi] interval or only use the upper half-plane [0, pi/2]

                if ybound[len(ybound) - 1] != np.pi / 2.:
                    self.nyi = sum(nyi) + 1
                    self.ny = self.nyi - 1
                    self.yi = ybound[0] + (ybound[1] - ybound[0]) * (
                        np.arange(nyi[0], dtype=np.float64) / float(nyi[0]))

                    for ipart in range(1, len(nyi) - 1):
                        # Now make sure that pi/2 will be a cell interface
                        #
                        # BUGFIX! 16-05-2012
                        # The grid was not symmetric to pi/2 when the grid contained multiple sections (i.e. len(nyi)>1)
                        # This is now fixed
                        if ybound[ipart] < np.pi / 2.:
                            dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                                np.arange(nyi[ipart], dtype=np.float64) / float(nyi[ipart]))
                        else:
                            if ybound[ipart] == np.pi / 2.:
                                dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                                    (np.arange(nyi[ipart] + 1, dtype=np.float64)) / (float(nyi[ipart])))
                            else:
                                dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                                    (np.arange(nyi[ipart], dtype=np.float64) + 1.) / float(nyi[ipart]))

                        self.yi = np.append(self.yi, dum)

                    ipart = len(nyi) - 1
                    if len(nyi) == 2:
                        dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                            (np.arange(nyi[ipart] + 1, dtype=np.float64)) / (float(nyi[ipart])))
                    else:
                        dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                            (np.arange(nyi[ipart], dtype=np.float64) + 1.) / float(nyi[ipart]))

                else:
                    self.nyi = sum(nyi) + 1
                    self.ny = self.nyi - 1
                    self.yi = ybound[0] + (ybound[1] - ybound[0]) * (
                        np.arange(nyi[0], dtype=np.float64) / float(nyi[0]))
                    for ipart in range(1, len(nyi) - 1):
                        # Now make sure that pi/2 will be a cell interface
                        #
                        # BUGFIX! 16-05-2012
                        # The grid was not symmetric to pi/2 when the grid contained multiple sections (i.e. len(nyi)>1)
                        # This is now fixed
                        if ybound[ipart] < np.pi / 2.:
                            dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                                np.arange(nyi[ipart], dtype=np.float64) / float(nyi[ipart]))
                        else:
                            dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                                (np.arange(nyi[ipart] + 1, dtype=np.float64)) / (float(nyi[ipart])))

                        self.yi = np.append(self.yi, dum)

                    ipart = len(nyi) - 1

                    if len(nyi) == 2:
                        dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                            (np.arange(nyi[ipart] + 1, dtype=np.float64)) / (float(nyi[ipart])))
                    else:
                        dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                            (np.arange(nyi[ipart], dtype=np.float64) + 1.) / float(nyi[ipart]))

                self.yi = np.append(self.yi, dum)
                self.y = 0.5 * (self.yi[0:self.ny] + self.yi[1:self.ny + 1])

            else:
                if self.act_dim[1] == 1:
                    self.nyi = nyi[0]
                    self.yi = ybound[0] + (ybound[1] - ybound[0]) * (
                        np.arange(self.nyi, dtype=np.float64) / float(self.nyi - 1.))
                    self.ny = self.nyi - 1
                    self.y = 0.5 * (self.yi[0:self.ny] + self.yi[1:self.ny + 1])
                else:
                    self.y = [0.]
                    self.yi = [0., 0., ]
                    self.ny = 1
                    self.nyi = 2
                    #
                    # Create the z axis

            if len(nzi) > 1:
                self.nzi = sum(nzi)
                self.nz = self.nzi - 1

                self.zi = zbound[0] + (zbound[1] - zbound[0]) * (
                    np.arange(nzi[0], dtype=np.float64) / float(nzi[0]))
                for ipart in range(1, len(nzi) - 1):
                    dum = zbound[ipart] + (zbound[ipart + 1] - zbound[ipart]) * (
                        np.arange(nzi[ipart], dtype=np.float64) / float(nzi[ipart]))
                    self.zi = np.append(self.zi, dum)
                ipart = len(nzi) - 1
                dum = zbound[ipart] + (zbound[ipart + 1] - zbound[ipart]) * (
                    np.arange(nzi[ipart], dtype=np.float64) / float(nzi[ipart] - 1))
                self.zi = np.append(self.zi, dum)
                self.z = 0.5 * (self.zi[0:self.nz] + self.zi[1:self.nz + 1])
            else:
                if self.act_dim[2] == 1:
                    self.nzi = nzi[0]
                    self.zi = zbound[0] + (zbound[1] - zbound[0]) * (
                        np.arange(self.nzi, dtype=np.float64) / float(self.nzi - 1))
                    self.nz = self.nzi - 1
                    self.z = 0.5 * (self.zi[0:self.nz] + self.zi[1:self.nz + 1])
                else:
                    self.z = np.array([0.])
                    self.zi = np.array([0., np.pi * 2.])
                    self.nz = 1
                    self.nzi = 2

    def writeSpatialGrid(self, fname='', old=False):
        """Writes the wavelength grid to a file (e.g. amr_grid.inp).

        Parameters
        ----------

        fname : str, optional
                File name into which the spatial grid should be written. If omitted 'amr_grid.inp' will be used.

        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """

        #
        # Write the spatial grid for radmc3d
        #
        if not old:
            if fname == '':
                fname = 'amr_grid.inp'

            print('Writing ' + fname)
            with open(fname, 'w') as wfile:
                # Format number
                wfile.write('%d\n' % 1)
                # AMR self.style (0=regular self. NO AMR)
                wfile.write('%d\n' % 0)
                # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
                if self.crd_sys == 'car':
                    wfile.write('%d\n' % 0)
                # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
                if self.crd_sys == 'sph':
                    wfile.write('%d\n' % 100)
                # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
                if self.crd_sys == 'cyl':
                    wfile.write('%d\n' % 200)
                # Gridinfo
                wfile.write('%d\n' % 0)

                # Active dimensions
                wfile.write('%d %d %d \n' % (self.act_dim[0], self.act_dim[1], self.act_dim[2]))
                # Grid size (x,y,z or r,phi,theta, or r,phi,z)
                wfile.write('%d %d %d \n' % (self.nx, self.ny, self.nz))
                for i in range(self.nxi):
                    wfile.write('%.9e\n' % self.xi[i])
                for i in range(self.nyi):
                    wfile.write('%.9e\n' % self.yi[i])
                for i in range(self.nzi):
                    wfile.write('%.9e\n' % self.zi[i])
            wfile.close()
        #
        # Write the spatial grid for radmc
        #
        else:

            fname = 'radius.inp'
            with open(fname, 'w') as wfile:

                print('Writing ' + fname)
                x = np.sqrt(self.xi[1:] * self.xi[:-1])
                wfile.write("%d\n" % self.nx)
                wfile.write(" \n")
                for i in range(self.nx):
                    wfile.write("%.7e\n" % x[i])

            fname = 'theta.inp'
            with open(fname, 'w') as wfile:
                print('Writing ' + fname)
                wfile.write("%d 1\n" % (self.ny / 2))
                wfile.write(" \n")
                for i in range(int(self.ny / 2)):
                    wfile.write("%.7e\n" % self.y[i])

    def readWavelengthGrid(self, fname=None, old=False):
        """Reads the wavelength grid

        Parameters
        ----------

        fname : str, optional
                File name from which the spatial grid should be read. If omitted 'wavelength_micron.inp' will be used.

        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """
        #
        # Read the radmc3d format
        #
        if not old:
            if fname is None:
                fname = 'wavelength_micron.inp'
            #
            # Read the frequency grid
            #
            print('Reading ' + fname)
            data = np.fromfile(fname, count=-1, sep=" ", dtype=np.float64)
            self.nfreq = np.int_(data[0])
            self.nwav = self.nfreq
            self.wav = data[1:]
            self.freq = nc.cc / self.wav * 1e4
        #
        # Read the old radmc format
        #
        else:
            if fname is None:
                fname = 'frequency.inp'

            print('Reading '+fname)
            data = np.fromfile(fname, count=-1, sep=" ", dtype=np.float64)
            self.nfreq = np.int_(data[0])
            self.nwav = self.nwav
            self.freq = data[1:]
            self.wav = nc.cc / self.freq * 1e4

        return

    def readSpatialGrid(self, fname='', old=False):
        """Reads the spatial grid

        Parameters
        ----------

        fname : str, optional
                File name from which the spatial grid should be read. If omitted 'amr_grid.inp' will be used.

        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """
        #
        # Read the radmc3d format
        #
        if not old:
            if fname == '':
                fname = 'amr_grid.inp'
                #
                # Read the spatial grid
                #

            print('Reading '+fname)
            data = np.fromfile(fname, count=-1, sep=" ", dtype=np.float64)
            hdr = np.array(data[:10], dtype=np.int_)
            data = data[10:]

            # Check the file format
            if hdr[0] != 1:
                msg = 'Unkonwn format number in amr_grid.inp'
                raise RuntimeError(msg)

            # Check the coordinate system
            if hdr[2] < 100:
                self.crd_sys = 'car'
            elif (hdr[2] >= 100) & (hdr[2] < 200):
                self.crd_sys = 'sph'
            elif (hdr[2] >= 200) & (hdr[2] < 300):
                self.crd_sys = 'cyl'
            else:
                raise ValueError('Unsupported coordinate system identification in the ' + fname + ' file.')

            # Get active dimensions
            self.act_dim = hdr[4:7]

            # Get the number of cells in each dimensions
            self.nx = hdr[7]
            self.ny = hdr[8]
            self.nz = hdr[9]
            self.nxi, self.nyi, self.nzi = self.nx + 1, self.ny + 1, self.nz + 1

            # Get the cell interfaces
            self.xi = data[:self.nxi]
            data = data[self.nxi:]
            self.yi = data[:self.nyi]
            data = data[self.nyi:]
            self.zi = data[:self.nzi]

            if self.crd_sys == 'car':
                self.x = (self.xi[0:self.nx] + self.xi[1:self.nx + 1]) * 0.5
                self.y = (self.yi[0:self.ny] + self.yi[1:self.ny + 1]) * 0.5
                self.z = (self.zi[0:self.nz] + self.zi[1:self.nz + 1]) * 0.5
            else:
                self.x = np.sqrt(self.xi[0:self.nx] * self.xi[1:self.nx + 1])
                self.y = (self.yi[0:self.ny] + self.yi[1:self.ny + 1]) * 0.5
                self.z = (self.zi[0:self.nz] + self.zi[1:self.nz + 1]) * 0.5

        #
        # Read the old radmc format
        #
        else:
            self.crd_sys = 'sph'
            self.act_dim = [1, 1, 0]

            #
            # Read the radial grid
            #
            data = np.fromfile('radius.inp', count=-1, sep=" ", dtype=np.float64)
            self.nx = np.int_(data[0])
            self.nxi = self.nx + 1
            self.x = data[1:]
            self.xi = np.zeros(self.nxi, dtype=float)
            self.xi[1:-1] = 0.5 * (self.x[1:] + self.x[:-1])
            self.xi[0] = self.x[0] - (self.xi[1] - self.x[0])
            self.xi[-1] = self.x[-1] + (self.x[-1] - self.xi[-2])

            #
            # Read the poloidal angular grid
            #

            data = np.fromfile('theta.inp', count=-1, sep=" ", dtype=np.float64)
            self.ny = np.int_(data[0]) * 2
            self.nyi = self.ny + 1
            self.y = np.zeros(self.ny, dtype=float)
            self.y[:self.ny//2] = data[2:]
            self.y[self.ny//2:] = np.pi - data[2:][::-1]
            self.yi = np.zeros(self.nyi, dtype=float)
            self.yi[1:-1] = 0.5 * (self.y[1:] + self.y[:-1])
            self.yi[self.ny] = np.pi * 0.5

            #
            # Create the azimuthal grid
            #
            self.nz = 1
            self.zi = np.array([0., 2. * np.pi], dtype=float)

        return

    def readGrid(self, old=False):
        """Reads the spatial (amr_grid.inp) and frequency grid (wavelength_micron.inp).

        Parameters
        ----------

        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """

        self.readSpatialGrid(old=old)
        self.readWavelengthGrid(old=old)

        return

    def getCellVolume(self):
        """Calculates the volume of grid cells.

        """
        if self.crd_sys == 'sph':

            if self.act_dim[0] == 0:
                raise ValueError('The first (r) dimension of a shserical grid is switched off')
            elif self.act_dim[1] == 0:
                if self.act_dim[2] == 0:
                    vol = np.zeros([self.nx, self.ny, self.nz], dtype=np.float64)
                    diff_r3 = self.xi[1:] ** 3 - self.xi[:-1] ** 3
                    diff_cost = 2.0
                    diff_phi = 2. * np.pi
                    for ix in range(self.nx):
                        vol[ix, 0, 0] = 1. / 3. * diff_r3[ix] * diff_cost * diff_phi

                else:
                    vol = np.zeros([self.nx, self.ny, self.nz], dtype=np.float64)
                    diff_r3 = self.xi[1:] ** 3 - self.xi[:-1] ** 3
                    diff_cost = 2.0
                    diff_phi = self.zi[1:] - self.zi[:-1]
                    for ix in range(self.nx):
                        for iz in range(self.nz):
                            vol[ix, 0, iz] = 1. / 3. * diff_r3[ix] * diff_cost * diff_phi[iz]

            elif self.act_dim[2] == 0:
                vol = np.zeros([self.nx, self.ny, self.nz], dtype=np.float64)
                diff_r3 = self.xi[1:] ** 3 - self.xi[:-1] ** 3
                diff_cost = np.cos(self.yi[:-1]) - np.cos(self.yi[1:])
                diff_phi = 2. * np.pi
                for ix in range(self.nx):
                    for iy in range(self.ny):
                        vol[ix, iy, :] = 1. / 3. * diff_r3[ix] * diff_cost[iy] * diff_phi

            else:
                vol = np.zeros([self.nx, self.ny, self.nz], dtype=np.float64)
                diff_r3 = self.xi[1:] ** 3 - self.xi[:-1] ** 3
                diff_cost = np.cos(self.yi[:-1]) - np.cos(self.yi[1:])
                diff_phi = self.zi[1:] - self.zi[:-1]
                for ix in range(self.nx):
                    for iy in range(self.ny):
                        vol[ix, iy, :] = 1. / 3. * diff_r3[ix] * diff_cost[iy] * diff_phi
        else:
            raise ValueError('Coordinate system ' + self.crd_sys + ' is not yet supported.')

        return vol
