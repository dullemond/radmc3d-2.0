"""
This module contains a class for handling variable data in radmc-3d
Note that since verion 0.30.2 of radmc3dPy, a lot has been changed
to make default actions as simple and automatic as possible. 
"""
from __future__ import absolute_import
from __future__ import print_function
import traceback
import os
import warnings

try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use the python module of RADMC-3D you need to install Numpy')
    print(traceback.format_exc())

from . import natconst as nc
from . import crd_trans
from . dustopac import *
from . reggrid import *
from . octree import *

def readGrid(sgrid=True, wgrid=True, sgrid_fname=None, wgrid_fname=None, old=False):
    """
    Reads the spatial and frequency grid.
    This function is an interface to radmc3dGrid.readGrid() or radmc3dOctree.readGrid().
    It will automatically determine which.

    Parameters
    ----------

    sgrid       : bool
                  If True the spatial grid will be read

    wgrid       : bool
                  If True the wavelength grid will be read

    sgrid_fname : str
                  File containing the spatial grid (default: amr_grid.inp)

    wgrid_fname : str
                  File containing the wavelength grid (default: wavelength_micron.inp)

    old         : bool
                  If True the format of the old 2D version of radmc3d (radmc) will be used
    Returns
    -------

    Returns an instance of the radmc3dGrid (for regular grid) or radmc3dOctree (for octree AMR) class
    """

    grid = None

    if sgrid:
        #grid = radmc3dGrid()
        #
        # Check the grid type
        #
        if not old:
            if sgrid_fname is None:
                sgrid_fname = 'amr_grid.inp'
            hdr = np.fromfile(sgrid_fname, count=7, sep="\n", dtype=np.int_)
            if hdr[1] == 0:
                grid = radmc3dGrid()
            elif hdr[1] == 1:
                grid = radmc3dOctree()
            else:
                raise ValueError('Unsupported amr_style' + ("%d" % hdr[1]) + '\n '
                                 + 'Only regular (0) or octree-like (1) AMR styles are supported')
            grid.readSpatialGrid(fname=sgrid_fname)
        else: 
            grid.readSpatialGrid(fname=sgrid_fname, old=old)

    if wgrid:
        if grid is None:
            grid = radmc3dGrid()

        if sgrid_fname is None:
            sgrid_fname = 'amr_grid.inp'

        if isinstance(grid, radmc3dOctree):
            grid.readWavelengthGrid(fname=wgrid_fname)
        else:
            grid.readWavelengthGrid(fname=wgrid_fname, old=old)

    return grid


class radmc3dData(object):
    """RADMC-3D data class.
        Reading and writing dust density/temperature, gas density/temperature/velocity,
        generating a legacy vtk file for visualization.


    Attributes
    ----------

    grid      : radmc3dGrid, radmc3dOctree
                Instance of the radmc3dGrid class, contains the spatial and frequency grids

    rhodust   : ndarray
                Dust density in g/cm^3

    dusttemp  : ndarray
                Dust temperature in K

    rhogas    : ndarray
                Gas density in g/cm^3

    ndens_mol : ndarray
                Number density of the molecule [molecule/cm^3]

    ndens_cp  : ndarray
                Number density of the collisional partner [molecule/cm^3]

    gasvel    : ndarray
                Gas velocity in cm/s

    gastemp   : ndarray
                Gas temperature in K

    vturb     : ndarray
                Mictroturbulence in cm/s

    taux      : ndarray
                Optical depth along the x (cartesian) / r (cylindrical) / r (spherical) dimension

    tauy      : ndarray
                Optical depth along the y (cartesian) / theta (cylindrical) / theta (spherical) dimension

    tauz      : ndarray
                Optical depth along the z (cartesian) / z (cylindrical) / phi (spherical) dimension

    sigmadust : ndarray
                Dust surface density in g/cm^2

    sigmagas  : ndarray
                Gas surface density in molecule/cm^2 (or g/cm^2 depending on the dimension of rhogas)
    """

    def __init__(self, grid=None):

        if grid:
            self.grid = grid
        else:
            self.grid   = readGrid()
            
        self.octree = self.grid.octree

        self.rhodust = np.zeros(0, dtype=np.float64)
        self.dusttemp = np.zeros(0, dtype=np.float64)
        self.rhogas = np.zeros(0, dtype=np.float64)
        self.ndens_mol = np.zeros(0, dtype=np.float64)
        self.ndens_cp = np.zeros(0, dtype=np.float64)
        self.gasvel = np.zeros(0, dtype=np.float64)
        self.gastemp = np.zeros(0, dtype=np.float64)
        self.vturb = np.zeros(0, dtype=np.float64)
        self.taux = np.zeros(0, dtype=np.float64)
        self.tauy = np.zeros(0, dtype=np.float64)
        self.tauz = np.zeros(0, dtype=np.float64)
        self.sigmadust = np.zeros(0, dtype=np.float64)
        self.sigmagas = np.zeros(0, dtype=np.float64)

    def _isBinary(self, fname):
        ext = fname[-5:]
        if ext=='.binp' or ext=='.bdat' or ext=='.bout':
            return True
        else:
            return False

    def _findDataFile(self, basename):
        hits  = 0
        exts  = ['.inp','.binp','.out','.bout','.dat','.bdat']
        filename = ''
        for x in exts:
            fname = basename+x
            if(os.path.isfile(fname)):
                hits    += 1
                filename = fname
        if hits>1:
            print('Warning: Multiple files basename '+basename+' found. Taking '+filename)
        return filename
        
    def _scalarfieldWriter(self, data=None, fname='', binary=True):
        """Writes a scalar field to a file.

        Parameters
        ----------

        data   : ndarray
                Scalar variable to be written

        fname  : str
                Name of the file containing a scalar variable

        binary : bool
                If True the file will be in binary format, if False the file format is formatted ASCII text

        """
        octree = self.octree
        with open(fname, 'w') as wfile:
            if binary:
                if octree:
                    if len(data.shape) == 1:
                        hdr = np.array([1, 8, data.shape[0], 1], dtype=np.int64)
                    elif len(data.shape) == 2:
                        hdr = np.array([1, 8, data.shape[0], data.shape[1]], dtype=np.int64)
                    else:
                        raise ValueError('Incorrect shape. Data stored in an Octree-type grid should have 1 or 2 '
                                         + 'dimensions, with the second dimension being the dust species. '
                                         + ' The data to be written has a dimension of ' + ("%d" % len(data.shape))
                                         + '\n No data has been written')
                    hdr.tofile(wfile)
                    data.flatten(order='f').tofile(wfile)
                else:
                    if len(data.shape) == 3:
                        hdr = np.array([1, 8, self.grid.nx * self.grid.ny * self.grid.nz], dtype=np.int64)
                        hdr.tofile(wfile)
                    elif len(data.shape) == 4:
                        hdr = np.array([1, 8, self.grid.nx * self.grid.ny * self.grid.nz, data.shape[3]],
                                       dtype=np.int64)
                        hdr.tofile(wfile)
                    # Now we need to flatten the dust density array since the Ndarray.tofile function writes the
                    # array always in C-order while we need Fortran-order to be written
                    if len(data.shape) == 4:
                        data = np.swapaxes(data, 0, 3)
                        data = np.swapaxes(data, 1, 2)
                        data.tofile(wfile)
                    elif len(data.shape) == 3:
                        data = np.swapaxes(data, 0, 2)
                        data.tofile(wfile)
                    else:
                        raise ValueError('Incorrect shape. Data stored in a regular grid should have 3 or 4 dimensions'
                                         + ' with the fourth dimension being the dust species. '
                                         + ' The data to be written has a dimension of ' + ("%d" % len(data.shape))
                                         + '\n No data has been written')
            else:
                if octree:
                    if len(data.shape) == 1:
                        hdr = np.array([1, data.shape[0], 1], dtype=np.int64)
                    elif len(data.shape) == 2:
                        hdr = np.array([1, data.shape[0], data.shape[1]], dtype=np.int64)
                    else:
                        raise ValueError('Incorrect shape. Data stored in an Octree-type grid should have 1 or 2 '
                                         + 'dimensions, with the second dimension being the dust species. '
                                         + ' The data to be written has a dimension of ' + ("%d" % len(data.shape))
                                         + '\n No data has been written')

                    hdr.tofile(wfile, sep=" ", format="%d\n")
                    data.flatten(order='f').tofile(wfile, sep=" ", format="%.9e\n")
                else:
                    if len(data.shape) == 3:
                        hdr = np.array([1, self.grid.nx * self.grid.ny * self.grid.nz], dtype=int)
                        hdr.tofile(wfile, sep=" ", format="%d\n")
                        # Now we need to flatten the dust density array since the Ndarray.tofile function writes the
                        # array always in C-order while we need Fortran-order to be written
                        data = np.swapaxes(data, 0, 2)
                        data.tofile(wfile, sep=" ", format="%.9e\n")

                    elif len(data.shape) == 4:
                        hdr = np.array([1, self.grid.nx * self.grid.ny * self.grid.nz, data.shape[3]], dtype=int)
                        hdr.tofile(wfile, sep=" ", format="%d\n")
                        # Now we need to flatten the dust density array since the Ndarray.tofile function writes the
                        # array always in C-order while we need Fortran-order to be written
                        data = np.swapaxes(data, 0, 3)
                        data = np.swapaxes(data, 1, 2)
                        data.tofile(wfile, sep=" ", format="%.9e\n")
                    else:
                        raise ValueError('Incorrect shape. Data stored in a regular grid should have 3 or 4 dimensions'
                                         + ' with the fourth dimension being the dust species. '
                                         + ' The data to be written has a dimension of ' + ("%d" % len(data.shape))
                                         + '\n No data has been written')

    def _scalarfieldReader(self, fname='', ndim=3):
        """Reads a scalar field from file.

        Parameters
        ----------

        fname  : str
                Name of the file containing a scalar variable

        ndim   : int
                Number of dimension of the data field (3 for gas variables, 4 for dust)
        Returns
        -------

        Returns a numpy Ndarray with the scalar field
        """
        octree = self.octree
        binary = self._isBinary(fname)
        data = None
        with open(fname, 'r') as rfile:
            if binary:
                # Read the header
                # hdr[0] = format number
                # hdr[1] = data precision (4=single, 8=double)
                # hdr[2] = nr of cells
                # hdr[3] = nr of dust species
                if ndim == 3:
                    hdr = np.fromfile(rfile, count=3, dtype=np.int64)
                else:
                    hdr = np.fromfile(rfile, count=4, dtype=np.int64)
                if octree:
                    if hdr[2] != self.grid.nLeaf:
                        print(hdr[1], self.grid.nLeaf)
                        msg = ('Number of cells in ' + fname + ' is different from that in amr_grid.inp'
                               + ' nr cells in ' + fname + ' : ' + ("%d" % hdr[2]) + '\n '
                               + ' nr of cells in amr_grid.inp : ' + ("%d" % self.grid.nLeaf))
                        raise ValueError(msg)

                    if hdr[1] == 8:
                        data = np.fromfile(rfile, count=-1, dtype=np.float64)
                    elif hdr[1] == 4:
                        data = np.fromfile(rfile, count=-1, dtype=np.float32)
                    else:
                        msg = 'Unknown datatype/precision in ' + fname + '. RADMC-3D binary files store 4 byte ' \
                              + 'floats or 8 byte doubles. The precision in the file header is ' + ("%d" % hdr[1])
                        raise TypeError(msg)

                    if ndim == 3:
                        if data.shape[0] == hdr[2]:
                            data = np.reshape(data, [self.grid.nLeaf, 1], order='f')
                        else:
                            msg = 'Internal inconsistency in data file, number of cell entries is different from ' \
                                  'indicated in the file header'
                            raise ValueError(msg)

                    else:
                        if data.shape[0] == hdr[2] * hdr[3]:
                            data = np.reshape(data, [self.grid.nLeaf, hdr[3]], order='f')
                        else:
                            msg = 'Internal inconsistency in data file, number of cell entries is different from ' \
                                  'indicated in the file header'
                            raise ValueError(msg)

                else:
                    if hdr[2] != (self.grid.nx * self.grid.ny * self.grid.nz):
                        msg = 'Number of grid cells in ' + fname + ' is different from that in amr_grid.inp ' \
                              + ' nr cells in ' + fname + ' : ' + ("%d" % hdr[2]) + '\n ' \
                              + ' nr of cells in amr_grid.inp : ' \
                              + ("%d" % (self.grid.nx * self.grid.ny * self.grid.nz))
                        raise ValueError(msg)

                    if hdr[1] == 8:
                        data = np.fromfile(rfile, count=-1, dtype=np.float64)
                    elif hdr[1] == 4:
                        data = np.fromfile(rfile, count=-1, dtype=np.float32)
                    else:
                        msg = 'Unknown datatype/precision in ' + fname + '. RADMC-3D binary files store 4 byte ' \
                              + 'floats or 8 byte doubles. The precision in the file header is ' \
                              + ("%d" % hdr[1])
                        raise TypeError(msg)

                    if ndim == 3:
                        if data.shape[0] == hdr[2]:
                            data = np.reshape(data, [1, self.grid.nz, self.grid.ny, self.grid.nx])
                        else:
                            msg = 'Internal inconsistency in data file, number of cell entries is different from ' \
                                  'indicated in the file header'
                            raise ValueError(msg)
                    else:
                        if data.shape[0] == hdr[2] * hdr[3]:
                            data = np.reshape(data, [hdr[3], self.grid.nz, self.grid.ny, self.grid.nx])
                        else:
                            msg = 'Internal inconsistency in data file, number of cell entries is different from ' \
                                  'indicated in the file header'
                            raise ValueError(msg)

                    # data = reshape(data, [hdr[3],self.grid.nz,self.grid.ny,self.grid.nx])
                    # We need to change the axis orders as Numpy always writes binaries in C-order while RADMC-3D
                    # uses Fortran-order
                    data = np.swapaxes(data, 0, 3)
                    data = np.swapaxes(data, 1, 2)
            else:
                if ndim == 3:
                    hdr = np.fromfile(rfile, count=2, sep=" ", dtype=np.int64)
                else:
                    hdr = np.fromfile(rfile, count=3, sep=" ", dtype=np.int64)

                if octree:
                    if hdr[1] != self.grid.nLeaf:
                        msg = 'Number of cells in ' + fname + ' is different from that in amr_grid.inp' \
                              + ' nr cells in ' + fname + ' : ' + ("%d" % hdr[1]) + '\n '\
                              + ' nr of cells in amr_grid.inp : ' + ("%d" % self.grid.nLeaf)
                        raise ValueError(msg)

                    data = np.fromfile(rfile, count=-1, sep=" ", dtype=np.float64)

                    if ndim == 3:
                        if data.shape[0] == hdr[1]:
                            data = np.reshape(data, [self.grid.nLeaf, 1], order='f')
                        else:
                            msg = 'Internal inconsistency in data file, number of cell entries is different from ' \
                                  'indicated in the file header'
                            raise ValueError(msg)

                    else:
                        if data.shape[0] == hdr[1] * hdr[2]:
                            data = np.reshape(data, [self.grid.nLeaf, hdr[2]], order='f')
                        else:
                            msg = 'Internal inconsistency in data file, number of cell entries is different from ' \
                                  'indicated in the file header. (' + str(data.shape[0]) + ' vs ' + \
                                  str(hdr2[1]*hdr[2]) + ')'
                            raise ValueError(msg)

                    # if ndim == 3:
                    #
                    #     data = data.reshape([hdr[1], hdr[2]], order='f')

                else:
                    if hdr[1] != (self.grid.nx * self.grid.ny * self.grid.nz):
                        msg = 'Number of grid cells in ' + fname + ' is different from that in amr_grid.inp '\
                              + ' nr cells in ' + fname + ' : ' + ("%d" % hdr[1]) + '\n '\
                              + ' nr of cells in amr_grid.inp : '\
                              + ("%d" % (self.grid.nx * self.grid.ny * self.grid.nz))
                        raise ValueError(msg)

                    data = np.fromfile(rfile, count=-1, sep=" ", dtype=np.float64)
                    if ndim == 3:
                        if data.shape[0] == hdr[1]:
                            data = np.reshape(data, [1, self.grid.nz, self.grid.ny, self.grid.nx])
                        else:
                            msg = 'Internal inconsistency in data file, number of cell entries is different from ' \
                                  'indicated in the file header'
                            raise ValueError(msg)
                    else:
                        if data.shape[0] == hdr[1] * hdr[2]:
                            data = np.reshape(data, [hdr[2], self.grid.nz, self.grid.ny, self.grid.nx])
                        else:
                            msg = 'Internal inconsistency in data file, number of cell entries is different from ' \
                                  'indicated in the file header'
                            raise ValueError(msg)

                    # data = reshape(data, [hdr[3],self.grid.nz,self.grid.ny,self.grid.nx])
                    # We need to change the axis orders as Numpy always writes binaries in C-order while RADMC-3D
                    # uses Fortran-order
                    data = np.swapaxes(data, 0, 3)
                    data = np.swapaxes(data, 1, 2)

        return data

    def getTauOneDust(self, idust=0, axis='', kappa=0.):
        """Calculates the optical depth of a single dust species along any given combination of the axes.

        Parameters
        ----------

        idust : int
                Index of the dust species whose optical depth should be calculated

        axis  : str
                Name of the axis/axes along which the optical depth should be calculated
                (e.g. 'x' for the first dimension or 'xyz' for all three dimensions)

        kappa : float
                Mass extinction coefficients of the dust species at the desired wavelength

        Returns
        -------

        Returns a dictionary with the following keys

            taux  : ndarray
                    optical depth along the first dimension
            tauy  : ndarray
                    optical depth along the second dimension

            (tauz is not yet implemented)
        """

        # Check along which axis should the optical depth be calculated
        do_taux = False
        do_tauy = False

        if 'x' in axis:
            do_taux = True
        if 'y' in axis:
            do_tauy = True

        # Calculate the optical depth along the x-axis (r in spherical coordinates)
        if do_taux:
            taux = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=np.float64)
            diff_x = self.grid.xi[1:] - self.grid.xi[:-1]
            taux[0, :, :] = self.rhodust[0, :, :, idust] * kappa * diff_x[0]
            for ix in range(1, self.grid.nx):
                taux[ix, :, :] = taux[ix - 1, :, :] + self.rhodust[ix, :, :, idust] * kappa * diff_x[ix]
        else:
            taux = [-1.]

        # Calculate the optical depth along the theta in spherical coordinates
        # Warning the formulation below is valid only in spherical coordinate sytem

        dum_x = np.zeros([self.grid.nx, self.grid.nz], dtype=np.float64)
        for iz in range(self.grid.nz):
            dum_x[:, iz] = self.grid.x

        if do_tauy:
            tauy = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=np.float64)
            diff_y = self.grid.yi[1:] - self.grid.yi[:-1]
            tauy[:, 0, :] = self.rhodust[:, 0, :, idust] * kappa * diff_y[0] * dum_x
            for iy in range(1, self.grid.ny):
                tauy[:, iy, :] = tauy[:, iy - 1, :] + self.rhodust[:, iy, :, idust] * kappa * diff_y[iy] * dum_x
        else:
            tauy = [-1.]

        return {'taux': taux, 'tauy': tauy}

    def getTau(self, idust=None, axis='xy', wav=0., kappa=None, old=False):
        """Calculates the optical depth along any given combination of the axes.

        Parameters
        ----------

        idust : list
                List of dust component indices whose optical depth should be calculated
                If multiple indices are set the total optical depth is calculated summing
                over all dust species in idust

        axis  : str
                Name of the axis/axes along which the optical depth should be calculated
                (e.g. 'x' for the first dimension or 'xyz' for all three dimensions)

        wav   : float
                Wavelength at which the optical depth should be calculated

        kappa : list, tuple
                If set it should be a list of mass extinction coefficients at the desired wavelength
                The number of elements in the list should be equal to that in the idust keyword

        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """

        # Check if the input idust indices can be found in rhoudust
        if len(self.rhodust.shape) == 3:
            ndust = 1
        else:
            ndust = self.rhodust.shape[3]

        if idust is None:
            idust = np.arange(ndust)

        if max(idust) > ndust:
            raise ValueError(' There are less number of dust species than some of the indices in idust')

        scatmat = None
        # If the kappa keyword is set it should be used during the optical depth calculation
        if kappa is not None:
            # Safety check
            if len(kappa) != len(idust):
                raise ValueError(' The number of kappa values should be identical to the number of specified '
                                 + 'dust species ')
        else:
            if not old:
                # Read the master opacity file to get the dustkappa file name extensions
                dum = radmc3dDustOpac()
                mo = dum.readMasterOpac()
                dummy_ext = mo['ext']
                scatmat = mo['scatmat']
                if len(dummy_ext) <= max(idust):
                    raise ValueError('There are less dust species specified in dustopac.inp than '
                                     + 'some of the specified idust indices')
                else:
                    ext = [dummy_ext[i] for i in idust]

        if 'x' in axis:
            self.taux = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=np.float64)
        if 'y' in axis:
            self.tauy = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=np.float64)

        for i in idust:
            if kappa is None:
                if old:
                    opac = radmc3dDustOpac()
                    opac.readOpac(ext=[("%d" % (i + 1))], old=True)
                    # opac = readOpac(ext=[("%d" % (i + 1))], old=True)
                else:
                    opac = radmc3dDustOpac()


                    ### BUG IN CODE NEXT LINE: FIX COULD BE ext=dummy_ext[i]    2019.10.07

                    
                    opac.readOpac(ext=ext[i], scatmat=scatmat)
                    #opac = readOpac(ext=ext[i], scatmat=scatmat)

                if not opac.ext:
                    return -1
                else:
                    if opac.wav[0][-1] > opac.wav[0][0]:
                        kabs = 10. ** np.interp(np.log10(np.array(wav)), np.log10(opac.wav[0]), np.log10(opac.kabs[0]))
                    else:
                        kabs = 10. ** np.interp(np.log10(np.array(wav)), np.log10(opac.wav[0][::-1]),
                                                np.log10(opac.kabs[0][::-1]))
                if opac.ksca[0][0] > 0:
                    if opac.wav[0][-1] > opac.wav[0][0]:
                        ksca = 10. ** np.interp(np.log10(np.array(wav)), np.log10(opac.wav[0]), np.log10(opac.ksca[0]))
                    else:
                        ksca = 10. ** np.interp(np.log10(np.array(wav)), np.log10(opac.wav[0][::-1]),
                                                np.log10(opac.ksca[0][::-1]))
                else:
                    ksca = np.array(kabs) * 0.

                print('Opacity at ' + ("%.2f" % wav) + 'um : ', kabs + ksca)
                dum = self.getTauOneDust(i, axis=axis, kappa=kabs + ksca)
            else:
                dum = self.getTauOneDust(i, axis=axis, kappa=kappa[i])

            if 'x' in axis:
                self.taux = self.taux + dum['taux']
            if 'y' in axis:
                self.tauy = self.tauy + dum['tauy']

    def getGasMass(self, mweight=2.3, rhogas=False):
        """Calculates the gas mass in radmc3dData.ndens_mol or radmc3dData.rhogas

        Parameters
        ----------
            mweight   : float
                        Molecular weight [atomic mass unit / molecule, i.e. same unit as mean molecular weight]

            rhogas    : bool, optional
                        If True the gas mass will be calculated from radmc3dData.rhogas, while if set to False
                        the gas mass will be calculated from radmc3dData.ndens_mol. The mweight parameter is only
                        required for the latter.

        Returns
        -------
            A single float being the gas mass in gramm
        """

        vol = self.grid.getCellVolume()
        if not rhogas:
            #
            # I'm not sure if this is the right way of doing it but right now I don't have a better idea
            #
            if isinstance(self.grid, radmc3dOctree):
                return (vol * self.ndens_mol[:, 0] * mweight * nc.mp).sum()
            else:
                return (vol * self.ndens_mol[:, :, :, 0] * mweight * nc.mp).sum()

        else:
            if isinstance(self.grid, radmc3dOctree):
                return (vol * self.rhogas[:, 0]).sum()
            else:
                return (vol * self.rhogas[:, :, :, 0]).sum()


    def getDustMass(self, idust=-1):
        """Calculates the dust mass in radmc3dData.rhodust

        Parameters
        ----------
            idust   : int
                      Dust index whose dust should be calculated. If it is set to -1 (default) the total
                      dust mass is calculated summing up all dust species

        Returns
        -------
            A single float being the dust mass in gramm
        """

        vol = self.grid.getCellVolume()
        if idust > 0:
            #
            # I'm not sure if this is the right way of doing it but right now I don't have a better idea
            #
            if isinstance(self.grid, radmc3dOctree):
                dmass = (vol * self.rhodust[:, idust]).sum()
            else:
                dmass = (vol * self.rhodust[:, :, :, idust]).sum()
        else:
            dmass = 0.
            if isinstance(self.grid, radmc3dOctree):
                for i in range(self.rhodust.shape[1]):
                    dmass += (vol * self.rhodust[:, i]).sum()
            else:
                for i in range(self.rhodust.shape[3]):
                    dmass += (vol * self.rhodust[:, :, :, i]).sum()

        return dmass

    def readDustDens(self, fname=None, old=False):
        """Reads the dust density.

        Parameters
        ----------

        fname   : str, optional
                  Name of the file that contains the dust density. If omitted 'dust_density.inp'
                  or 'dust_density.binp' is used.

        old     : bool, optional
                  If set to True the file format of the previous, 2D version of radmc will be used

        """
        assert self.grid is not None, "radmc3dData object MUST have a .grid subobject."

        octree = self.octree
        
        #
        # Read radmc3d output
        #
        if not old:
            if fname is None:
                fname = self._findDataFile('dust_density')
            if os.path.isfile(fname):
                print('Reading '+fname)
                self.rhodust = self._scalarfieldReader(fname=fname, ndim=4)
        #
        # Read the output of the previous 2d version of the code
        #
        else:
            fname = 'dustdens.inp'
            if os.path.isfile(fname):
                print('Reading '+fname)
                data = np.fromfile(fname, count=-1, sep=" ", dtype=np.float64)
                # 4 element header: Nr of dust species, nr, ntheta, ?
                hdr = np.array(data[:4], dtype=np.int64)
                data = np.reshape(data[4:], [hdr[0], hdr[1], hdr[2]])
                data = np.moveaxis(data, 0, -1)[:, :, np.newaxis, :]

                self.rhodust = np.zeros([hdr[1], hdr[2]*2, 1, hdr[0]], dtype=np.float64)
                self.rhodust[:, :hdr[2], 0, :] = data[:, :, 0, :]
                self.rhodust[:, hdr[2]:, 0, :] = data[:, ::-1, 0, :]

    def readDustTemp(self, fname=None, old=False):
        """Reads the dust temperature.

        Parameters
        ----------

        fname   : str, optional
                  Name of the file that contains the dust temperature. If omitted 'dust_temperature.dat'
                  or 'dust_temperature.bdat' is used.

        old     : bool, optional
                  If True dust temperature will be written in the old RADMC format
        """
        assert self.grid is not None, "radmc3dData object MUST have a .grid subobject."

        octree = self.octree

        if not old:
            if fname is None:
                fname = self._findDataFile('dust_temperature')
            if os.path.isfile(fname):
                print('Reading '+fname)
                self.dusttemp = self._scalarfieldReader(fname=fname, ndim=4)
        else:
            fname = 'dusttemp_final.dat'
            if os.path.isfile(fname):
                print('Reading '+fname)

                data = np.fromfile(fname, count=-1, sep=" ", dtype=np.float64)
                # 4 element header: Nr of dust species, nr, ntheta, ?
                hdr = np.array(data[:4], dtype=np.int64)
                data = data[4:]
                self.dusttemp = np.zeros([hdr[1], hdr[2]*2, 1, hdr[0]], dtype=np.float64)

                ncell = hdr[1] * hdr[2]
                for idust in range(hdr[0]):
                    # We need to get rid of the dust species index that is written inbetween the
                    # temperature arrays
                    data = data[1:]
                    self.dusttemp[:, :hdr[2], 0, idust] = np.reshape(data[:ncell], [hdr[1], hdr[2]])
                self.dusttemp[:, hdr[2]:, 0, :] = self.dusttemp[:, :hdr[2], 0, :][:, ::-1, :]

    def readGasVel(self, fname=None):
        """Reads the gas velocity.

        Parameters
        -----------

        fname : str, optional
                Name of the file that contains the gas velocity
                If omitted 'gas_velocity.inp' or 'gas_velocity.binp' is used.
        """
        assert self.grid is not None, "radmc3dData object MUST have a .grid subobject."
        
        octree = self.octree

        if fname is None:
            fname = self._findDataFile('gas_velocity')
        if os.path.isfile(fname):
            binary = self._isBinary(fname)
            if binary:
                print('Reading '+fname)
                # If we have an octree grid
                with open(fname, 'r') as rfile:
                    if isinstance(self.grid, radmc3dOctree):
                        hdr = np.fromfile(rfile, count=3, dtype=np.int64)
                        if hdr[2] != self.grid.nLeaf:
                            raise ValueError('Number of cells in ' + fname + ' is different from that in amr_grid.inp'
                                             + ' nr cells in ' + fname + ' : ' + ("%d" % hdr[2]) + '\n '
                                             + ' nr of cells in amr_grid.inp : '
                                             + ("%d" % self.grid.nLeaf))

                        if hdr[1] == 8:
                            self.gasvel = np.fromfile(rfile, count=-1, dtype=np.float64)
                        elif hdr[1] == 4:
                            self.gasvel = np.fromfile(rfile, count=-1, dtype=np.float32)
                        else:
                            raise TypeError(
                                'Unknown datatype/precision in ' + fname + '. RADMC-3D binary files store 4 byte '
                                + 'floats or 8 byte doubles. The precision in the file header is '
                                + ("%d" % hdr[1]))
                        self.gasvel = np.reshape(self.gasvel, [self.nLeaf, 3], order='f')

                    else:
                        hdr = np.fromfile(rfile, count=3, dtype=np.int64)
                        if hdr[2] != self.grid.nx * self.grid.ny * self.grid.nz:
                            raise ValueError('Number of grid cells in ' + fname + ' is different from that in amr_grid.inp '
                                             + ' nr cells in ' + fname + ' : ' + ("%d" % hdr[2]) + '\n '
                                             + ' nr of cells in amr_grid.inp : '
                                             + ("%d" % (self.grid.nx * self.grid.ny * self.grid.nz)))

                        if hdr[1] == 8:
                            self.gasvel = np.fromfile(rfile, count=-1, dtype=np.float64)
                        elif hdr[1] == 4:
                            self.gasvel = np.fromfile(rfile, count=-1, dtype=float)
                        else:
                            raise TypeError(
                                'Unknown datatype/precision in ' + fname + '. RADMC-3D binary files store 4 byte '
                                + 'floats or 8 byte doubles. The precision in the file header is '
                                + ("%d" % hdr[1]))
                        self.gasvel = np.reshape(self.gasvel, [self.grid.nz, self.grid.ny, self.grid.nx, 3])
                        self.gasvel = np.swapaxes(self.gasvel, 0, 2)
            else:
                print('Reading ' + fname)
                with open(fname, 'r') as rfile:
                    # Two element header 1 - iformat, 2 - Nr of cells
                    hdr = np.fromfile(rfile, count=2, sep=" ", dtype=np.float64)
                    #hdr = np.array(data[:2], dtype=np.int64)
                    if octree:
                        if self.grid.nLeaf != hdr[1]:
                            msg = 'Number of cells in ' + fname + ' is different from that in amr_grid.inp'\
                                  + ' nr cells in ' + fname + ' : ' + ("%d" % hdr[1]) + '\n '\
                                  + ' nr of cells in amr_grid.inp : ' +  ("%d" % self.grid.nLeaf)
                            raise RuntimeError(msg)
                        data = np.fromfile(rfile, count=-1, sep=" ", dtype=np.float64)
                        self.gasvel = np.reshape(data, [hdr[1], 3])
                    else:
                        if (self.grid.nx * self.grid.ny * self.grid.nz) != hdr[1]:
                            msg = 'Number of grid cells in ' + fname + ' is different from that in amr_grid.inp ' \
                                  + ' nr cells in ' + fname + ' : ' + ("%d" % hdr[1]) + '\n '\
                                  + ' nr of cells in amr_grid.inp : '\
                                  + ("%d" % (self.grid.nx * self.grid.ny * self.grid.nz))
                            raise RuntimeError(msg)
                        data = np.fromfile(rfile, count=-1, sep=" ", dtype=np.float64)
                        self.gasvel = np.reshape(data, [self.grid.nz, self.grid.ny, self.grid.nx, 3])
                        self.gasvel = np.swapaxes(self.gasvel, 0, 2)

    def readVTurb(self, fname=None):
        """Reads the turbulent velocity field.

        Parameters
        ----------

        fname   : str, optional
                  Name of the file that contains the turbulent velocity field
                  If omitted 'microturbulence.inp' or 'microturbulence.binp' is used.
        """
        assert self.grid is not None, "radmc3dData object MUST have a .grid subobject."

        octree = self.octree

        if fname is None:
            fname = self._findDataFile('microturbulence')
        if os.path.isfile(fname):
            print('Reading microturbulence')
            self.vturb = self._scalarfieldReader(fname=fname, ndim=3)
            if octree:
                self.vturb = np.squeeze(self.vturb)

    def readGasDens(self, fname=None, ispec=''):
        """Reads the gas density.

        Parameters
        ----------

        fname   : str,optional
                  Name of the file that contains the gas temperature. If omitted 'numberdens_<ispec>.inp'
                  or 'numberdens_<ispec>.binp' is used.

        ispec   : str
                  File name extension of the 'numberdens_ispec.*'

        """
        assert self.grid is not None, "radmc3dData object MUST have a .grid subobject."

        octree = self.octree

        if fname is None:
            fname = self._findDataFile('numberdens_' + ispec)
        if os.path.isfile(fname):
            print('Reading gas density (' + fname + ')')
            self.ndens_mol = self._scalarfieldReader(fname=fname, ndim=3)
            if octree:
                self.ndens_mol = np.squeeze(self.ndens_mol)

    def readGasTemp(self, fname=None):
        """Reads the gas temperature.

        Parameters
        ----------

        fname   : str,optional
                  Name of the file that contains the gas temperature. If omitted 'gas_temperature.inp'
                  or 'gas_tempearture.binp' is used.
        """
        assert self.grid is not None, "radmc3dData object MUST have a .grid subobject."

        octree = self.octree

        if fname is None:
            fname = self._findDataFile('gas_temperature')
        if os.path.isfile(fname):
            print('Reading gas temperature')
            self.gastemp = self._scalarfieldReader(fname=fname, ndim=3)
            if octree:
                self.gastemp = np.squeeze(self.gastemp)

    def writeDustDens(self, fname='', binary=True, old=False):
        """Writes the dust density.

        Parameters
        ----------

        fname   : str, optional
                  Name of the file into which the dust density should be written. If omitted 'dust_density.inp' is used.

        binary  : bool
                  If true the data will be written in binary format, otherwise the file format is ascii

        old     : bool, optional
                  If set to True the file format of the previous, 2D version of radmc will be used

        """
        octree = self.octree
        #
        # Write dust density for radmc3d
        #
        if not old:
            if fname == '':
                if binary:
                    fname = 'dust_density.binp'
                else:
                    fname = 'dust_density.inp'

            print('Writing ' + fname)

            if octree:
                self._scalarfieldWriter(data=self.grid.convArrTree2Leaf(self.rhodust), fname=fname, binary=binary)
            else:
                self._scalarfieldWriter(data=self.rhodust, fname=fname, binary=binary)

        #
        # Write dust density for the previous 2D version of the code
        #
        else:
            if self.rhodust.shape[2] > 1:
                raise ValueError('Wrong dimensions for dust density. RADMC (the predecessor code to RADMC-3D is '
                                 + 'strictly 2D, yet the dust density to be written in the RADMC format is 3D.')

            if fname == '':
                fname = 'dustdens.inp'

            with open(fname, 'w') as wfile:

                ntheta = int(self.grid.ny / 2)
                wfile.write("%d %d %d 1\n" % (self.rhodust.shape[3], self.grid.nx, ntheta))
                wfile.write(" \n")
                for idust in range(self.rhodust.shape[3]):
                    for ix in range(self.grid.nx):
                        for iy in range(ntheta):
                            wfile.write("%.7e\n" % self.rhodust[ix, iy, 0, idust])

    def writeDustTemp(self, fname='', binary=True):
        """Writes the dust density.

        Parameters
        ----------

        fname   : str, optional
                Name of the file into which the dust density should be written. If omitted 'dust_density.inp' is used.

        binary  : bool
                  If true the data will be written in binary format, otherwise the file format is ascii

        octree  : bool, optional
                  If the data is defined on an octree-like AMR
        """
        octree = self.octree
        if fname == '':
            if binary:
                fname = 'dust_temperature.bdat'
            else:
                fname = 'dust_temperature.dat'

        print('Writing ' + fname)
        if octree:
            self._scalarfieldWriter(data=self.grid.convArrTree2Leaf(self.dusttemp), fname=fname, binary=binary)
        else:
            self._scalarfieldWriter(data=self.dusttemp, fname=fname, binary=binary)

    def writeGasDens(self, fname=None, ispec='', binary=True):
        """Writes the gas density.

        Parameters
        ----------

        fname   : str, optional
                  Name of the file into which the data will be written. If omitted "numberdens_xxx.inp" and
                  "numberdens_xxx.binp" will be used for ascii and binary format, respectively
                  (xxx is the name of the molecule).

        ispec   : str
                  File name extension of the 'numberdens_ispec.inp' (if binary=True 'numberdens_ispec.binp')
                  file into which the gas density should be written

        binary  : bool
                  If true the data will be written in binary format, otherwise the file format is ascii

        octree  : bool, optional
                  If the data is defined on an octree-like AMR
        """
        octree = self.octree
        if fname is None:
            if ispec == '':
                raise ValueError('Unknown ispec. Without knowing the name of the gas species the gas density '
                                 + 'cannot be written, as the filename "numberdens_ispec.inp" cannot be generated.')
            if binary:
                fname = 'numberdens_' + ispec + '.binp'
            else:
                fname = 'numberdens_' + ispec + '.inp'

        print('Writing ' + fname)
        if octree:
            self._scalarfieldWriter(data=self.grid.convArrTree2Leaf(self.ndens_mol), fname=fname, binary=binary)
        else:
            self._scalarfieldWriter(data=self.ndens_mol, fname=fname, binary=binary)

    def writeGasTemp(self, fname='', binary=True):
        """Writes the gas temperature.

        Parameters
        ----------

        fname : str, optional
                Name of the file into which the gas temperature should be written. If omitted
                'gas_temperature.inp' (if binary=True 'gas_tempearture.binp') is used.

        binary : bool
                If true the data will be written in binary format, otherwise the file format is ascii

        """
        octree = self.octree
        if fname == '':
            if binary:
                fname = 'gas_temperature.binp'
            else:
                fname = 'gas_temperature.inp'

        print('Writing ' + fname)

        if octree:
            self._scalarfieldWriter(data=self.grid.convArrTree2Leaf(self.gastemp), fname=fname, binary=binary)
        else:
            self._scalarfieldWriter(data=self.gastemp, fname=fname, binary=binary)

    def writeGasVel(self, fname='', binary=True):
        """Writes the gas velocity.

        Parameters
        ----------

        fname  : str, optional
                Name of the file into which the gas temperature should be written.
                If omitted 'gas_velocity.inp' (if binary=True 'gas_velocity.binp') is used.

        binary : bool
                If true the data will be written in binary format, otherwise the file format is ascii
        """
        octree = self.octree
        if binary:
            if fname == '':
                fname = 'gas_velocity.binp'

            print('Writing ' + fname)
            with open(fname, 'w') as wfile:

                if octree:
                    #
                    # Check if the gas velocity contains the full tree or only the leaf nodes
                    #
                    if self.gasvel.shape[0] == self.grid.nLeaf:
                        hdr = np.array([1, 8, self.gasvel.shape[0]], dtype=np.int64)
                        hdr.tofile(wfile)
                        self.gasvel.flatten(order='f').tofile(wfile)
                    else:
                        hdr = np.array([1, 8, self.grid.nLeaf], dtype=np.int64)
                        hdr.tofile(wfile)
                        dummy = self.grid.convArrTree2Leaf(self.gasvel)
                        dummy.flatten(order='f').tofile(wfile)
                        # self.gasvel.flatten(order='f').tofile(wfile)

                else:
                    hdr = np.array([1, 8, self.grid.nx * self.grid.ny * self.grid.nz], dtype=int)
                    hdr.tofile(wfile)
                    # Now we need to change the axis orders since the Ndarray.tofile function writes the
                    # array always in C-order while we need Fortran-order to be written
                    self.gasvel = np.swapaxes(self.gasvel, 0, 2)
                    self.gasvel.tofile(wfile)

                    # Switch back to the original axis order
                    self.gasvel = np.swapaxes(self.gasvel, 0, 2)
        else:
            if fname == '':
                fname = 'gas_velocity.inp'

            print('Writing ' + fname)

            # If we have an octree grid
            if octree:
                #
                # Check if the gas velocity contains the full tree or only the leaf nodes
                #
                if self.gasvel.shape[0] == self.grid.nLeaf:
                    hdr = '1\n'
                    hdr += ("%d\n" % self.gasvel.shape[0])
                    try:
                        np.savetxt(fname, self.gasvel, fmt="%.9e %.9e %.9e", header=hdr, comments='')
                    except Exception as e:
                        print(e)
                else:
                    hdr = '1\n'
                    hdr += ("%d\n" % self.grid.nLeaf)
                    dummy = self.grid.convArrTree2Leaf(self.gasvel)
                    try:
                        np.savetxt(fname, dummy, fmt="%.9e %.9e %.9e", header=hdr, comments='')
                    except Exception as e:
                        print(e)

            else:
                # hdr = "1\n"
                # hdr += ('%d'%(self.grid.nx*self.grid.ny*self.grid.nz))
                # np.savetxt(fname, self.gasvel, fmt="%.9 %.9 %.9", header=hdr, comments='')

                with open(fname, 'w') as wfile:

                    wfile.write('%d\n' % 1)
                    wfile.write('%d\n' % (self.grid.nx * self.grid.ny * self.grid.nz))

                    for iz in range(self.grid.nz):
                        for iy in range(self.grid.ny):
                            for ix in range(self.grid.nx):
                                wfile.write("%9e %9e %9e\n" % (self.gasvel[ix, iy, iz, 0], self.gasvel[ix, iy, iz, 1],
                                                               self.gasvel[ix, iy, iz, 2]))

    def writeVTurb(self, fname='', binary=True):
        """Writes the microturbulence file.

        Parameters
        ----------

        fname   : str, optional
                  Name of the file into which the turubulent velocity field should be written.
                  If omitted 'microturbulence.inp' (if binary=True 'microturbuulence.binp') is used.

        binary  : bool
                  If true the data will be written in binary format, otherwise the file format is ascii
        """
        octree = self.octree
        if fname == '':
            if binary:
                fname = 'microturbulence.binp'
            else:
                fname = 'microturbulence.inp'

        print('Writing ' + fname)

        if octree:
            self._scalarfieldWriter(data=self.grid.convArrTree2Leaf(self.vturb), fname=fname, binary=binary)
        else:
            self._scalarfieldWriter(data=self.vturb, fname=fname, binary=binary)

    def writeVTK(self, vtk_fname='', ddens=False, dtemp=False, idust=None, gdens=False, gvel=False, gtemp=False):
        """Writes physical variables to a legacy vtk file.

        Parameters
        ----------

        vtk_fname : str
                    Name of the file to be written, if not specified 'radmc3d_data.vtk' will be used

        ddens     : bool
                    If set to True the dust density will be written to the vtk file

        dtemp     : bool
                    If set to True the dust temperature will be written to the vtk file

        idust     : list
                    List of indices that specifies which dust component should be written
                    if not set then the first dust species (zero index) will be used

        gdens     : bool
                    If set to True the gas density will be written to the vtk file

        gtemp     : bool
                    If set to True the gas temperature will be written to the vtk file

        gvel      : bool
                    If set to True the gas velocity will be written to the vtk file
        """

        if isinstance(idust, int):
            idust = [idust]

        if (ddens is True) | (dtemp is True):
            if idust is None:
                msg = 'Unknown dust species. Dust density or temperature is to be written but  ' \
                      + 'it has not been specified for which dust species.'
                raise ValueError(msg)


        if vtk_fname == '':
            vtk_fname = 'radmc3d_data.vtk'
        else:
            vtk_fname = str(vtk_fname)

        #
        # Get the grid
        #

        x = self.grid.xi
        # For the theta axis I leave out the poles
        #  The current cell type is hexahedron and the cells near the pole are
        #    rather tetrahedra than hexahedra and this is not yet implemented
        y = np.array(self.grid.yi[1:self.grid.nyi - 1])
        z = self.grid.zi
        nxi = x.shape[0]
        nyi = y.shape[0]
        nzi = z.shape[0]

        #
        # Gas velocity field (Should be corner centered)
        # TODO
        #  The lines below should be double checked and re-implemented as
        #  the re-mapping of the cell centered velocity field to the cell corners
        #  is physically not correct and very messy...
        #
        if gvel:
            vgas = np.zeros([nxi, nyi, nzi, 3], dtype=np.float64)
            vgas[0:nxi - 1, 0:nyi, 0:nzi - 1, :] = self.gasvel[:, 1:nyi + 1, :, :]
            vgas[nxi - 1, :, :, :] = vgas[nxi - 2, :, :, :]
            vgas[:, nyi - 1, :, :] = vgas[:, nyi - 2, :, :]
            vgas[:, :, nzi - 1, :] = vgas[:, :, nzi - 2, :]

        #
        # Header
        #
        with open(vtk_fname, 'w') as wfile:
            wfile.write('%s\n' % '# vtk DataFile Version 3.0')
            wfile.write('%s\n' % 'RADMC-3D Data')
            wfile.write('%s\n' % 'ASCII')
            wfile.write('%s\n' % 'DATASET UNSTRUCTURED_GRID')

            #
            # Write out the coordinates of the cell corners
            #
            wfile.write('%s\n' % ('POINTS ' + str(nxi * nyi * nzi).strip() + ' double'))
            print('Writing POINTS: ')
            for ix in range(nxi):
                print(ix, nxi)
                for iy in range(nyi):
                    for iz in range(nzi):
                        crd = crd_trans.ctransSph2Cart([x[ix], z[iz], y[iy]])
                        wfile.write('%.9e %9e %9e\n' % (crd[0], crd[1], crd[2]))

                        # ---------------------------------------------------------------------------------------------
                        # Write out the indices of the cell interface mesh that define a
                        # hexahedron (VTK cell type #12)
                        #
                        # The indexing of a hexahedron is as follows
                        #
                        #                  7________6
                        #                 /|      / |
                        #                / |     /  |
                        #               4_------5   |             z ^   ^ y
                        #               |  3____|___2               |  /
                        #               | /     |  /                | /
                        #               |/      | /                 |/
                        #               0-------1                   0-----> x
                        #
                        # ---------------------------------------------------------------------------------------------

            wfile.write('%s %d %d\n' % ('CELLS ', ((nxi - 1) * (nyi - 1) * (nzi - 1)),
                                        ((nxi - 1) * (nyi - 1) * (nzi - 1)) * 9))

            for ix in range(nxi - 1):
                print('Writing CELL COORDINATES: ', ix, self.grid.nxi - 2)
                for iy in range(nyi - 1):
                    for iz in range(nzi - 1):
                        id1 = nzi * nyi * ix + nzi * iy + iz
                        id2 = nzi * nyi * ix + nzi * (iy + 1) + iz
                        id4 = nzi * nyi * ix + nzi * iy + ((iz + 1) % (nzi - 1))
                        id3 = nzi * nyi * ix + nzi * (iy + 1) + ((iz + 1) % (nzi - 1))
                        id5 = nzi * nyi * (ix + 1) + nzi * iy + iz
                        id6 = nzi * nyi * (ix + 1) + nzi * (iy + 1) + iz
                        id7 = nzi * nyi * (ix + 1) + nzi * (iy + 1) + ((iz + 1) % (nzi - 1))
                        id8 = nzi * nyi * (ix + 1) + nzi * iy + ((iz + 1) % (nzi - 1))

                        line = np.array([8, id1, id2, id3, id4, id5, id6, id7, id8])
                        line.tofile(wfile, sep=' ', format='%d')
                        wfile.write('\n')
                        #
                        # Now write out the type of each cell (#12)
                        #
            wfile.write('%s %d\n' % ('CELL_TYPES', ((nxi - 1) * (nyi - 1) * (nzi - 1))))

            for ix in range(nxi - 1):
                for iy in range(nyi - 1):
                    for iz in range(nzi - 1):
                        wfile.write('%d\n' % 12)
                        #
                        # Now write out the corner centered velocities
                        #

            if gvel:
                wfile.write('%s %d\n' % ('POINT_DATA', (nxi * nyi * nzi)))
                wfile.write('%s\n' % 'VECTORS gas_velocity double')
                for ix in range(nxi):
                    print('Writing velocity : ', ix, nxi - 1)
                    for iy in range(nyi):
                        for iz in range(nzi):
                            vsph = np.array([vgas[ix, iy, iz, 0], vgas[ix, iy, iz, 2], vgas[ix, iy, iz, 1]])
                            vxyz = crd_trans.vtransSph2Cart([x[ix], z[iz], y[iy]], vsph)

                            wfile.write('%.9e %.9e %.9e\n' % (vxyz[0], vxyz[1], vxyz[2]))

            #
            # Write out the cell centered scalars
            #
            wfile.write('%s %d\n' % ('CELL_DATA', ((nxi - 1) * (nyi - 1) * (nzi - 1))))

            if ddens:
                for ids in idust:
                    wfile.write('%s\n' % ('SCALARS dust_density_' + str(int(ids)) + ' double'))
                    wfile.write('%s\n' % 'LOOKUP_TABLE default')

                    for ix in range(nxi - 1):
                        print('Writing dust density : ', ix, nxi - 2)
                        for iy in range(nyi - 1):
                            for iz in range(nzi - 1):
                                wfile.write('%.9e\n' % self.rhodust[ix, iy, iz, ids])

            if dtemp:
                for ids in idust:
                    wfile.write('%s\n' % ('SCALARS dust_temperature_' + str(int(ids)) + ' double'))
                    wfile.write('%s\n' % 'LOOKUP_TABLE default')

                    for ix in range(nxi - 1):
                        print('Writing dust temperature : ', ix, nxi - 2)
                        for iy in range(nyi - 1):
                            for iz in range(nzi - 1):
                                wfile.write('%.9e\n' % self.dusttemp[ix, iy, iz, ids])

            if gdens:
                wfile.write('%s\n' % 'SCALARS gas_numberdensity double')
                wfile.write('%s\n' % 'LOOKUP_TABLE default')

                for ix in range(nxi - 1):
                    print('Writing gas density : ', ix, nxi - 2)
                    for iy in range(nyi - 1):
                        for iz in range(nzi - 1):
                            wfile.write('%.9e\n' % self.ndens_mol[ix, iy, iz])

            if gtemp:
                wfile.write('%s\n' % 'SCALARS gas_temperature double')
                wfile.write('%s\n' % 'LOOKUP_TABLE default')

                for ix in range(nxi - 1):
                    print('Writing dust temperature : ', ix, nxi - 2)
                    for iy in range(nyi - 1):
                        for iz in range(nzi - 1):
                            wfile.write('%.9e\n' % self.gastemp[ix, iy, iz])

    def getSigmaDust(self, idust=-1):
        """Calculates the dust surface density.

        Parameters
        ----------

        idust : int, optional
                Index of the dust species for which the surface density should be calculated
                if omitted the calculated surface density will be the sum over all dust species
        """

        # Calculate the volume of each grid cell
        vol = self.grid.getCellVolume()
        # Dustmass in each grid cell
        if len(self.rhodust.shape) > 3:
            if idust >= 0:
                mass = vol * self.rhodust[:, :, :, idust]
            else:
                mass = vol * self.rhodust.sum(3)
        else:
            mass = vol * self.rhodust

        # Calculate the surface of each grid facet in the midplane
        surf = np.zeros([self.grid.nx, self.grid.nz], dtype=np.float64)
        diff_r2 = (self.grid.xi[1:] ** 2 - self.grid.xi[:-1] ** 2) * 0.5
        diff_phi = self.grid.zi[1:] - self.grid.zi[:-1]
        for ix in range(self.grid.nx):
            surf[ix, :] = diff_r2[ix] * diff_phi

        # Now get the surface density
        dum = mass.sum(1)
        self.sigmadust = dum / surf

    def getSigmaGas(self):
        """Calculates the gas surface density.
        This method uses radmc3dData.rhogas to calculate the surface density, thus the
        unit of surface density depends on the unit of radmc3dData.rhogas (g/cm^2 or molecule/cm^2)
        """

        # Calculate the volume of each grid cell
        vol = self.grid.getCellVolume()
        # Total number of molecules / gas mass in each grid cell
        mass = vol * self.rhogas
        # Calculate the surface are of each grid facet in the midplane
        surf = np.zeros([self.grid.nx, self.grid.nz], dtype=np.float64)
        diff_r2 = (self.grid.xi[1:] ** 2 - self.grid.xi[:-1] ** 2) * 0.5
        diff_phi = self.grid.zi[1:] - self.grid.zi[:-1]
        for ix in range(self.grid.nx):
            surf[ix, :] = diff_r2[ix] * diff_phi

        # Now get the surface density
        dum = np.squeeze(mass.sum(1))
        self.sigmagas = dum / np.squeeze(surf)

