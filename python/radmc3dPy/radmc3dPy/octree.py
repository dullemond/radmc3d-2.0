"""This module contains a class for handling Octree mesh
"""
from __future__ import absolute_import
from __future__ import print_function
import traceback
import inspect
import importlib

try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use the python module of RADMC-3D you need to install Numpy')
    print(traceback.format_exc())

from . import natconst as nc

class radmc3dOctree(object):
    """
    Octree-like object with switchable resolution in each dimension

    Attributes
    ----------

    xi                : ndarray
                        Base grid cell interface grid in the first dimension

    yi                : ndarray
                        Base grid cell interface grid in the second dimension

    zi                : ndarray
                        Base grid cell interface grid in the third dimension

    xc                : ndarray
                        Base grid cell center grid in the first dimension

    yc                : ndarray
                        Base grid cell center grid in the second dimension

    zc                : ndarray
                        Base grid cell center grid in the third dimension

    x                 : ndarray
                        Tree cell center array in the first dimension

    y                 : ndarray
                        Tree cell center array in the second dimension

    z                 : ndarray
                        Tree cell center array in the third dimension

    dx                : ndarray
                        Tree cell halfwidth array in the first dimension

    dy                : ndarray
                        Tree cell halfwidth array in the second dimension

    dz                : ndarray
                        Tree cell halfwidth array in the third dimension

    leafID            : ndarray
                        Leaf index array,  mapping between a full tree and an array containing only the leaves

    isLeaf            : ndarray
                        Boolean array for the cell type (True - leaf, False - branch)

    level             : ndarray
                        Level array (base grid is level 0)

    parentID          : ndarray
                        Array containing the index of the parent cell (currently unused, only needed if we go up
                        in the tree)

    childID           : list
                        List of children indices. Each list element is an ndarray with nChild elements containing the
                        child indices

    act_dim           : list
                        A three element array to indicate which dimension is active, i.e. which dimensions are the
                        cells resolved (0 - inactive, 1 - active)

    nCell             : int
                        Nr of cells (both branch and leaf) in the tree

    nxRoot            : int
                        Nr of cells in the base grid in the first dimension

    nyRoot            : int
                        Nr of cells in the base grid in the second dimension

    nzRoot            : int
                        Nr of cells in the base grid in the third dimension

    nLeaf             : int
                        Nr of leaf cells (i.e. true, unresolved grid cells)

    nBranch           : int
                        Nr of branches (i.e. resolved cells)

    nChild            : int
                        Nr of children (i.e. 8, 4, or 2 for 3, 2, 1 active dimensions, respectively)

    levelMax          : int
                        Highest actual level in the tree (Base grid has a level of 0 and the level increases)

    levelMaxLimit     : int
                        Highest allowed level in the tree (only used in tree building)

    crd_sys           : {'car', 'sph'}
                        Coordinate system type cartesian or spherical


    """

    def __init__(self):

        #
        # Arrays for the base grid (redundant information but improves speed and the memory overhead is not that large)
        #
        self.xi = np.zeros(0, dtype=np.float64)
        self.yi = np.zeros(0, dtype=np.float64)
        self.zi = np.zeros(0, dtype=np.float64)
        self.xc = np.zeros(0, dtype=np.float64)
        self.yc = np.zeros(0, dtype=np.float64)
        self.zc = np.zeros(0, dtype=np.float64)
        #
        # Grid cell center arrays
        #
        self.x = np.zeros(0, dtype=np.float64)
        self.y = np.zeros(0, dtype=np.float64)
        self.z = np.zeros(0, dtype=np.float64)
        #
        # Grid cell half(!)widths (this is in principle also redundant information as the only relevant information is
        # the cell width in the base grid, which is a single float for each dimension. The cell width at any level can
        # be simply calculated as cellwidth_base * 2^level. I still need to check how the performance would be affected
        # if I'd calculate the cell width on the fly instead of storing it in an array
        #
        self.dx = np.zeros(0, dtype=np.float64)
        self.dy = np.zeros(0, dtype=np.float64)
        self.dz = np.zeros(0, dtype=np.float64)
        #
        # Leaf index array - mapping between a full array and an array containing only leaves
        #
        self.leafID = np.zeros(0, dtype=np.int_)
        #
        # Boolean array for the cell type (True - leaf, False - branch)
        #
        self.isLeaf = np.zeros(0, dtype=bool)
        #
        # Level array (base grid is level 0)
        #
        self.level = np.zeros(0, dtype=np.int_)
        #
        # Array containing the index of the parent cell (currently unused, only needed if we go up in the tree)
        #
        self.parentID = np.zeros(0, dtype=np.int_)
        #
        # List of children indices
        #
        self.childID = []
        self.cID = np.zeros(0, dtype=int)
        #
        # Number of cells in the whole tree and for the base grid
        #
        self.nCell = 0
        self.nxRoot = 0
        self.nyRoot = 0
        self.nzRoot = 0
        #
        # Highest actual level in the tree (Base grid has a level of 0 and the level increases)
        #
        self.levelMax = 0
        #
        # Highest allowed level in the tree (only used in tree building)
        #
        self.levelMaxLimit = 0
        #
        # Nr of branches (i.e. resolved cells)
        #
        self.nBranch = 0
        #
        # Nr of leaf cells (i.e. true, unresolved grid cells)
        #
        self.nLeaf = 0
        #
        # Coordinate sytem and active (resolvable) dimensions
        #
        self.crd_sys = 'car'
        self.act_dim = [1, 1, 1]
        self.grid_style = 1
        self.octree = True
        self.nChild = 0
        #
        # Stuff used for tree building
        #
        self.model = None
        #
        # Variable meant to be used internally only
        #
        self.cellIDCur = -1

        self.nwav = 0
        self.nfreq = 0
        self.wav = np.zeros(0, dtype=np.int_)
        self.freq = np.zeros(0, dtype=np.int_)

        self.counter = -1

    def getCellVolume(self, fullTree=False):
        """
        Calculates the grid cell volume

        Parameters
        ----------

        fullTree    : bool, optional
                      If True the cell volumes of the full tree (including both branches and leaves) will be
                      calculated, while if set to False (default) the volume of only the leaf cells will be calculated

        Returns
        -------
        An linear ndarray containing the cell volumes
        """
        if self.crd_sys == 'car':
            vol = (2. * self.dx) * (2. * self.dy) * (2. * self.dz)
        else:
            vol = 1. / 3. * (self.x + self.dx)**3 - (self.x - self.dx)**3 \
                                                      * (np.cos(self.y - self.dy) - np.cos(self.y + self.dy)) \
                                                      * self.dz

        if not fullTree:
            ii = (self.leafID >= 0)
            dummy_vol = np.array(vol)
            vol = np.zeros(self.nLeaf, dtype=np.float64)
            vol[self.leafID[ii]] = dummy_vol[ii]

        return vol

    def _getContainerLeafIDRec(self, crd=(), cellID=-1):
        """
        Recursive function to find the tree index of a leaf that contains a given coordinate

        Parameters
        ----------
        crd         : tuple
                      List/tuple/ndarray containing the tree dimensional coordinates of the point

        cellID      : int
                      Cell index
        """

        xmin = self.x[cellID] - self.dx[cellID]
        xmax = self.x[cellID] + self.dx[cellID]
        ymin = self.y[cellID] - self.dy[cellID]
        ymax = self.y[cellID] + self.dy[cellID]
        zmin = self.z[cellID] - self.dz[cellID]
        zmax = self.z[cellID] + self.dz[cellID]

        if self.isLeaf[cellID]:
            if (((crd[0] >= xmin) & (crd[0] < xmax)) &
                    ((crd[1] >= ymin) & (crd[1] < ymax)) &
                    ((crd[2] >= zmin) & (crd[2] < zmax))):
                return cellID
            else:
                return None
        else:
            dum = None
            for i in range(self.nChild):
                dum = self._getContainerLeafIDRec(crd, self.childID[cellID][i])
                if dum is not None:
                    break
            return dum

    def getContainerLeafID(self, crd=()):
        """
        Finds the tree index of a leaf that contains a given coordinate

        Parameters
        ----------
        crd         : tuple
                      List/tuple/ndarray containing the tree dimensional coordinates of the point
        """

        leafID = -1

        if (crd[0] < self.xi[0]) | (crd[0] > self.xi[-1]):
            return leafID
        if (crd[1] < self.yi[0]) | (crd[1] > self.yi[-1]):
            return leafID
        if (crd[2] < self.zi[0]) | (crd[2] > self.zi[-1]):
            return leafID

        ix = np.searchsorted(self.xi, crd[0])
        iy = np.searchsorted(self.yi, crd[1])
        iz = np.searchsorted(self.zi, crd[2])

        if self.xi[ix] != crd[0]:
            ix -= 1
        if self.yi[iy] != crd[1]:
            iy -= 1
        if self.zi[iz] != crd[2]:
            iz -= 1

        if crd[0] == self.xi[-1]:
            ix = self.nxRoot - 1
        if crd[1] == self.yi[-1]:
            iy = self.nyRoot - 1
        if crd[2] == self.zi[-1]:
            iz = self.nzRoot - 1

        ind = iz * self.nyRoot * self.nxRoot + iy * self.nxRoot + ix
        leafID = self._getContainerLeafIDRec(crd, ind)

        return leafID

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

        if wbound is None:
            raise ValueError('Unknown wbound. Wavelength grid bondaries are not specified.')
        if nw is None:
            raise ValueError('Unknown nw. Number of wavelength grid points are not specified.')

        self.nwav = nw[0]
        self.wav = wbound[0] * (wbound[1] / wbound[0])**(np.arange(nw[0], dtype=np.float64) / nw[0])

        for ipart in range(1, len(nw) - 1):
            dum = wbound[ipart] * (wbound[ipart + 1]
                                   / wbound[ipart])**(np.arange(nw[ipart], dtype=np.float64) / nw[ipart])
            self.wav = np.append(self.wav, dum)

        ipart = len(nw) - 1
        dum = wbound[ipart] * (wbound[ipart + 1]
                               / wbound[ipart])**(np.arange(nw[ipart], dtype=np.float64) / (nw[ipart] - 1.))
        self.wav = np.append(self.wav, dum)
        self.nwav = self.wav.shape[0]
        self.freq = nc.cc / self.wav * 1e4
        self.nfreq = self.nwav

    def readGrid(self):
        """
        Reads the spatial and wavelength grids from files
        """
        self.readWavelengthGrid()
        self.readSpatialGrid()

    # --------------------------------------------------------------------------------------------------
    def readWavelengthGrid(self, fname='wavelength_micron.inp'):
        """
        Function to read the wavelength/frequency grid

        Parameters
        ----------

        fname       : str, optional
                      Name of the file to read the wavelength grid from (if not specified wavelenth_micron.inp will
                      be used)

        """
        #
        # Read the radmc3d format
        #
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


    def writeWavelengthGrid(self, fname='', old=False):
        """Wriites the wavelength grid to a file (e.g. wavelength_micron.inp).

        Parameters
        ----------

        fname  : str, optional
                 File name into which the wavelength grid should be written. If omitted 'wavelength_micron.inp' will
                 be used

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

    def putNode(self, crd=(), cellsize=(), level=None, parentID=-1, cellID=None):
        """
        Function to put the data of a single node into the tree. This funcion assumes that all the arrays
        have already been allocated for the tree so input cell indices must refer to already existing array elements.

        Parameters
        ----------
        crd      : tuple
                   Cell center coordinates of the node

        cellsize : tuple
                   Full size of the cell in each dimension

        level    : int
                   Level of the cell in the tree

        parentID : int
                   Tree index of the parent cell

        cellID   : int
                   Tree index of the cell to be added
        """

        #
        # Add the cell centre and cell half width to the arrays
        #
        self.x[cellID] = crd[0]
        self.y[cellID] = crd[1]
        self.z[cellID] = crd[2]

        self.dx[cellID] = cellsize[0] * 0.5
        self.dy[cellID] = cellsize[1] * 0.5
        self.dz[cellID] = cellsize[2] * 0.5

        self.isLeaf[cellID] = True
        self.level[cellID] = level
        self.parentID[cellID] = parentID
        self.childID.append(np.zeros(self.nChild, dtype=np.int_))

        return

    def resolveNodes(self, rsIDs=None):
        """
        Resolve multiple nodes simultaneously and add the children of the resolved node to the tree arrays extending
        the tree array

        Parameters
        ----------
        rsIDs       : list
                      List/tuple/array of indices of the resolvable cell in the tree array
        """

        if isinstance(rsIDs, np.ndarray):
            ncell = rsIDs.shape[0]
        else:
            ncell = len(rsIDs)

        self.nLeaf -= ncell
        self.nBranch += ncell
        self.nLeaf += ncell * self.nChild

        x = np.zeros(ncell * self.nChild, dtype=np.float64)
        y = np.zeros(ncell * self.nChild, dtype=np.float64)
        z = np.zeros(ncell * self.nChild, dtype=np.float64)
        dx = np.zeros(ncell * self.nChild, dtype=np.float64)
        dy = np.zeros(ncell * self.nChild, dtype=np.float64)
        dz = np.zeros(ncell * self.nChild, dtype=np.float64)
        isLeaf = np.zeros(ncell * self.nChild, dtype=bool)
        level = np.zeros(ncell * self.nChild, dtype=np.int_)
        parentID = np.zeros(ncell * self.nChild, dtype=np.int_)
        ind = np.arange(ncell, dtype=np.int_) * self.nChild
        nx = self.x.shape[0]

        xc_offset = None
        yc_offset = None
        zc_offset = None
        #
        # Generate the cell center offsets for a proper octree
        #
        if self.nChild == 8:
            xc_offset = np.array([-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dx[rsIDs][0]
            yc_offset = np.array([-0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dy[rsIDs][0]
            zc_offset = np.array([-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5], dtype=np.float64) * self.dz[rsIDs][0]

        #
        # If we are resolving only two dimensions
        #
        elif self.nChild == 4:

            if self.act_dim[0] == 0:
                xc_offset = np.array([0., 0., 0., 0.], dtype=np.float64)
                yc_offset = np.array([-0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dy[rsIDs][0]
                zc_offset = np.array([-0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dz[rsIDs][0]
            elif self.act_dim[1] == 0:
                xc_offset = np.array([-0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dx[rsIDs][0]
                yc_offset = np.array([0., 0., 0., 0.], dtype=np.float64)
                zc_offset = np.array([-0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dz[rsIDs][0]
            elif self.act_dim[2] == 0:
                xc_offset = np.array([-0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dx[rsIDs][0]
                yc_offset = np.array([-0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dy[rsIDs][0]
                zc_offset = np.array([0., 0., 0., 0.], dtype=np.float64)

        #
        # Generate the cell center offsets if only two dimensions are resolved
        #
        elif self.nChild == 2:

            if self.act_dim[0] == 1:
                xc_offset = np.array([-0.5, 0.5], dtype=np.float64) * self.dx[rsIDs][0]
                yc_offset = np.array([0., 0.], dtype=np.float64)
                zc_offset = np.array([0., 0.], dtype=np.float64)

            if self.act_dim[1] == 1:
                xc_offset = np.array([0., 0.], dtype=np.float64)
                yc_offset = np.array([-0.5, 0.5], dtype=np.float64) * self.dy[rsIDs][0]
                zc_offset = np.array([0., 0.], dtype=np.float64)

            if self.act_dim[2] == 1:
                xc_offset = np.array([0., 0.], dtype=np.float64)
                yc_offset = np.array([0., 0.], dtype=np.float64)
                zc_offset = np.array([-0.5, 0.5], dtype=np.float64) * self.dz[rsIDs][0]

        #
        # Generate the cell center offsets if only one dimension is resolved
        #
        for i in range(self.nChild):
            x[ind + i] = self.x[rsIDs] + xc_offset[i]
            y[ind + i] = self.y[rsIDs] + yc_offset[i]
            z[ind + i] = self.z[rsIDs] + zc_offset[i]
            if self.act_dim[0] == 1:
                dx[ind + i] = np.zeros(ncell, dtype=np.float64) + self.dx[rsIDs][0] * 0.5
            else:
                dx[ind + i] = np.zeros(ncell, dtype=np.float64) + self.dx[rsIDs][0]

            if self.act_dim[1] == 1:
                dy[ind + i] = np.zeros(ncell, dtype=np.float64) + self.dy[rsIDs][0] * 0.5
            else:
                dy[ind + i] = np.zeros(ncell, dtype=np.float64) + self.dy[rsIDs][0]

            if self.act_dim[2] == 1:
                dz[ind + i] = np.zeros(ncell, dtype=np.float64) + self.dz[rsIDs][0] * 0.5
            else:
                dz[ind + i] = np.zeros(ncell, dtype=np.float64) + self.dz[rsIDs][0]

            isLeaf[ind + i] = np.ones(ncell, dtype=bool)
            level[ind + i] = np.zeros(ncell, dtype=np.int_) + self.level[rsIDs][0] + 1
            parentID[ind + i] = rsIDs

        childID = []
        cid = np.arange(self.nChild, dtype=np.int_)
        for i in range(ncell * self.nChild):
            childID.append(cid)

        #
        # Add the child cells to the tree
        #
        self.x = np.append(self.x, x)
        self.y = np.append(self.y, y)
        self.z = np.append(self.z, z)
        self.dx = np.append(self.dx, dx)
        self.dy = np.append(self.dy, dy)
        self.dz = np.append(self.dz, dz)
        self.isLeaf = np.append(self.isLeaf, isLeaf)
        self.level = np.append(self.level, level)
        self.parentID = np.append(self.parentID, parentID)
        self.childID.extend(childID)

        #
        # Now add the child indices to the parents
        #
        for i in range(ncell):
            ind = rsIDs[i]
            self.childID[ind] = cid + nx + i * self.nChild
            self.isLeaf[ind] = False
            #
        # Update the array length variables
        #
        self.nCell = self.x.shape[0]
        self.cID = np.arange(self.nCell, dtype=np.int_)

        return

    def makeSpatialGrid(self, ppar=None, levelMaxLimit=None, dfunc=None, model='', **kwargs):
        """
        Function to create an octree-like AMR grid

        Parameters
        ----------

        ppar             : dictionary
                            Dictionary containing all input parameters of the model (from the problem_params.inp file)

        model            : str
                            Name of the model to be used in the tree building

        dfunc            : function
                            A user defined function that decides whether an AMR grid cell should be refined

        levelMaxLimit    : int, optional
                            Highest allowable level of the tree. This keyword is optional. If not specified at input as
                            a separate keyword, levelMaxLimit should be present in the problem_params.inp file.

        """

        #
        # Set the model
        #
        self.setModel(model)

        if dfunc is None:
            fnamelist = [f[0] for f in inspect.getmembers(self.model) if inspect.isfunction(f[1])]
            if 'decisionFunction' in fnamelist:
                dfunc = self.model.decisionFunction

            else:
                raise NameError('decisionFunction is not defined. \n '
                                + 'It is required to decide when a node / cell should be resolved. \n '
                                + 'Decision function should be given either as a dfunc keyword argument in \n'
                                + 'setup function call or should be implemented in the model as decisionFunction()')

        if levelMaxLimit is None:
            if 'levelMaxLimit' not in ppar.keys():
                raise KeyError('Unknown levelMaxLimit. \n It is required for AMR-style grid generation. levelMaxLimit '
                               + 'should be defined either in problem_params.inp or '
                               + 'in the call of the makeSpatialGrid/problemSetupDust/problemSetupGas functions.')
            else:
                self.levelMaxLimit = ppar['levelMaxLimit']
        else:
            self.levelMaxLimit = levelMaxLimit
        self.crd_sys = ppar['crd_sys']
        self.act_dim = [1, 1, 1]
        if ppar['nx'] == 0:
            self.act_dim[0] = 0
        if ppar['ny'] == 0:
            self.act_dim[1] = 0
        if ppar['nz'] == 0:
            self.act_dim[2] = 0

            #
        # Generate the base grid
        #
        self.xi = ppar['xbound'][0] + np.arange(ppar['nx'][0] + 1, dtype=float) / float(ppar['nx'][0]) \
                  * (ppar['xbound'][1] - ppar['xbound'][0])
        self.yi = ppar['ybound'][0] + np.arange(ppar['ny'][0] + 1, dtype=float) / float(ppar['ny'][0]) \
                  * (ppar['ybound'][1] - ppar['ybound'][0])
        self.zi = ppar['zbound'][0] + np.arange(ppar['nz'][0] + 1, dtype=float) / float(ppar['nz'][0]) \
                  * (ppar['zbound'][1] - ppar['zbound'][0])
        self.xc = 0.5 * (self.xi[1:] + self.xi[:-1])
        self.yc = 0.5 * (self.yi[1:] + self.yi[:-1])
        self.zc = 0.5 * (self.zi[1:] + self.zi[:-1])

        cellsize_x = self.xi[1] - self.xi[0]
        cellsize_y = self.yi[1] - self.yi[0]
        cellsize_z = self.zi[1] - self.zi[0]
        self.nxRoot = ppar['nx'][0]
        self.nyRoot = ppar['ny'][0]
        self.nzRoot = ppar['nz'][0]

        self.levelMax = 0

        ind = 0
        for iz in range(self.nzRoot):
            for iy in range(self.nyRoot):
                for ix in range(self.nxRoot):
                    #
                    # Add nodes to the tree
                    #
                    self.x = np.append(self.x, self.xc[ix])
                    self.y = np.append(self.y, self.yc[iy])
                    self.z = np.append(self.z, self.zc[iz])

                    self.dx = np.append(self.dx, [cellsize_x * 0.5])
                    self.dy = np.append(self.dy, [cellsize_y * 0.5])
                    self.dz = np.append(self.dz, [cellsize_z * 0.5])

                    self.isLeaf = np.append(self.isLeaf, True)
                    self.level = np.append(self.level, 0)
                    self.parentID = np.append(self.parentID, -1)
                    self.childID.append(np.zeros(self.nChild, dtype=np.int_))

                    self.nLeaf += 1
                    ind += 1
        self.nCell = self.x.shape[0]

        #
        # Now build the tree
        #
        if 1 in self.act_dim:
            self.nChild = 2**(np.array(self.act_dim, dtype=int).sum())

        if self.nChild > 0:
            print('Adaptive Mesh Refinement (AMR) is active')
        txt = 'Active dimensions : '
        for i in range(3):
            if self.act_dim[i] == 1:
                txt += ("%d " % i)
        print(txt)

        #
        # Now go level by level and check which cells are to be resolved and resolve what's necessary
        #
        for ilev in range(self.levelMaxLimit):
            print('Resolving level ' + ("%d" % ilev))
            #
            # Select the cells at the current level
            #
            cID = np.arange(self.x.shape[0], dtype=int)
            ii = (self.level == ilev)

            if True in ii:
                #
                # Check which cells to resolve
                #
                resolve = dfunc(self.x[ii], self.y[ii], self.z[ii], self.dx[ii][0], self.dy[ii][0], self.dz[ii][0],
                                model=self.model, ppar=ppar, **kwargs)

                jj = resolve[self.level[ii] < self.levelMaxLimit]
                #
                # If there are some to resolve do so
                #
                if True in jj:
                    ncell2resolve = cID[ii][jj].shape[0]
                    print('Cells to resolve at this level : ', ncell2resolve)
                    self.resolveNodes(rsIDs=cID[ii][jj])
                    self.levelMax += 1
                else:
                    print('No cells to resolve at this level')
            else:
                print('No cells to resolve at this level')

        self.childID = np.array(self.childID)
        #
        # Print out some statistics
        #
        print('Tree building done')
        print('Maximum tree depth : ', self.levelMax)
        print('Nr of branches     : ', self.nBranch)
        print('Nr of leaves       : ', self.nLeaf)
        # ncells_fullgrid = self.nChild**self.levelMax * self.nxRoot * self.nyRoot * self.nzRoot
        # cell_fraction = float(self.nLeaf + self.nBranch) / ncells_fullgrid
        # print 'Using '+("%.3f"%(cell_fraction*100))+'% memory of a regular grid at max resolution'

        self.generateLeafID()
        return

    def setModel(self, model=''):
        """
        Sets the model to be used for tree building

        Parameters
        ----------
        model       : str
                      Name of the model
        """

        try:
            self.model = importlib.import_module(model)
        except ImportError:
            try:
                # self.model = __import__('radmc3dPy.models.' + model, fromlist=[''])
                self.model = importlib.import_module('radmc3dPy.models.' + model)
            except ImportError:
                print(model + '.py could not be imported. \n '
                      + 'The model files should either be in the current working directory \n '
                      + 'or in the radmc3d python module directory')
                print(traceback.format_exc())

    def _selfCheckCounterRec(self, cellID=None):
        """
        Recursive function for consistency check of the tree
        """

        if self.isLeaf[cellID]:
            self.counter[0] += 1
            if self.childID[cellID].max() > self.counter[2]:
                self.counter[2] = self.childID[cellID].max()
        else:
            self.counter[1] += 1
            for i in range(self.nChild):
                self._selfCheckCounterRec(cellID=self.childID[cellID][i])

        return

    def selfCheck(self):
        """
        Performs a self-check of the tree allocation and report it to the screen
        """

        self.counter = np.zeros([3], dtype=np.int_)
        nRoot = self.nxRoot * self.nyRoot * self.nzRoot
        for i in range(nRoot):
            self._selfCheckCounterRec(cellID=i)

        print('Tree consistency check')
        print('Tree depth      : ' + ("%d" % self.levelMax))
        print('Nr of leaves    : ' + ("%d" % self.counter[0]) + " should be " + ("%d" % self.nLeaf))
        print('Nr of branches  : ' + ("%d" % self.counter[1]) + " should be " + ("%d" % self.nBranch))
        print('Nr of cells     : ' + ("%d" % self.nCell) + " should be " + ("%d" % (self.nBranch + self.nLeaf)))
        print('Leaf array      : ' + ("%d" % self.isLeaf.shape[0]))
        print('Level array     : ' + ("%d" % self.level.shape[0]))
        print('ParentID array  : ' + ("%d" % self.parentID.shape[0]))
        print('ChildID list    : ' + ("%d" % len(self.childID)))
        print('Max childID     : ' + ("%d" % self.counter[2]))
        print('x array         : ' + ("%d" % self.x.shape[0]))
        print('y array         : ' + ("%d" % self.y.shape[0]))
        print('z array         : ' + ("%d" % self.z.shape[0]))
        print('dx array        : ' + ("%d" % self.dx.shape[0]))
        print('dy array        : ' + ("%d" % self.dy.shape[0]))
        print('dz array        : ' + ("%d" % self.dz.shape[0]))

        return

    def _generateLeafIDRec(self, cellID=None):
        """
        Recursive function to generate the leaf indices
        """

        if self.isLeaf[cellID]:
            self.cellIDCur += 1
            self.leafID[cellID] = self.cellIDCur
        else:
            for i in range(self.nChild):
                self._generateLeafIDRec(self.childID[cellID][i])

    def generateLeafID(self):
        """
        Function to generate the cell index mapping from arrays containing the full tree and those containing
        only the leaves
        """

        print('Generating leaf indices')
        self.leafID = np.zeros(self.nCell, dtype=np.int_) - 1
        self.cellIDCur = -1
        nRoot = self.nxRoot * self.nyRoot * self.nzRoot
        for i in range(nRoot):
            self._generateLeafIDRec(i)

        print('Done')

    def convArrLeaf2Tree(self, var=None):
        """
        Converts a leaf array to full tree size.

        Parameters
        ----------
        var     : ndarray
                  A one or two dimensional ndarray with the first dimension is the size of the full tree

        Returns
        -------
        A one or two dimensional ndarray with size of of the full tree in the first dimension
        """

        ii = (self.leafID >= 0)
        ndim = len(var.shape)
        treeVar = None

        if ndim == 1:
            treeVar = np.zeros(self.nLeaf + self.nBranch, dtype=var.dtype)
            treeVar[ii] = var[self.leafID[ii]]
        elif ndim == 2:
            treeVar = np.zeros([self.nLeaf + self.nBranch, var.shape[1]], dtype=var.dtype)
            for i in range(var.shape[1]):
                treeVar[ii, i] = var[self.leafID[ii], i]
        else:
            assert (ndim >= 2), 'Incorrect shape for input array. Input array has too many dimensions. ' \
                                + 'Octree AMR only supports one or two dimensional arrays with the first dimension ' \
                                + 'being the cell/ spatial dimension'

            assert (ndim < 1), 'Incorrect shape for input array. Input array has too few dimensions (i.e. < 1). '

        return treeVar

    def convArrTree2Leaf(self, var=None):
        """
        Converts a data array to leaf size. The input is a scalar or vector variable defined at all nodes and the
        returned variable will only represent values at leaf nodes thereby reduced in length compared to the input.

        Parameters
        ----------
        var     : ndarray
                  A one or two dimensional ndarray with the first dimension is the size of the full tree

        Returns
        -------
        A one or two dimensional ndarray with size of nLeaf in the first dimension
        """

        ii = (self.leafID >= 0)
        ndim = len(var.shape)
        leafVar = None

        if ndim == 1:
            leafVar = np.zeros(self.nLeaf, dtype=var.dtype)
            leafVar[self.leafID[ii]] = var[ii]
        elif ndim == 2:
            leafVar = np.zeros([self.nLeaf, var.shape[1]], dtype=var.dtype)
            for i in range(var.shape[1]):
                leafVar[self.leafID[ii], i] = var[ii, i]
        else:
            assert (ndim >= 2), 'Incorrect shape for input array. Input array has too many dimensions. ' \
                                'Octree AMR only supports one or two dimensional arrays with the first dimension being'\
                                'the cell/ spatial dimension'

            assert (ndim < 1), 'Incorrect shape for input array. Input array has too few dimensions (i.e. < 1). '

        return leafVar

    def writeSpatialGrid(self, fname=''):
        """
        Writes the wavelength grid to a file (e.g. amr_grid.inp).

        Parameters
        ----------

        fname : str, optional
                File name into which the spatial grid should be written. If omitted 'amr_grid.inp' will be used.

        """

        if fname == '':
            fname = 'amr_grid.inp'

        print('Writing ' + fname)
        with open(fname, 'w') as wfile:
            wfile.write('%d\n' % 1)  # Format number

            wfile.write('\n')

            wfile.write('%d\n' % 1)  # AMR self.style (0=regular, 1 - Octree, 10 - Layered)
            if self.crd_sys == 'car':
                wfile.write('%d\n' % 0)  # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
            if self.crd_sys == 'sph':
                wfile.write('%d\n' % 100)  # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
            if self.crd_sys == 'cyl':
                wfile.write('%d\n' % 200)  # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
            wfile.write('%d\n' % 0)  # Gridinfo

            wfile.write('\n')

            # Active dimensions
            wfile.write('%d %d %d \n' % (self.act_dim[0], self.act_dim[1], self.act_dim[2]))
            # Grid size (x,y,z or r,phi,theta, or r,phi,z)
            wfile.write(
                '%d %d %d \n' % (self.nxRoot, self.nyRoot, self.nzRoot))

            wfile.write('\n')

            # Max. refinement level, Nr of leaves and Nr of branches
            wfile.write('%d %d %d \n' % (self.levelMax, self.nLeaf, self.nBranch + self.nLeaf))
            wfile.write('\n')
            for i in range(self.nxRoot + 1):
                wfile.write('%.9e\n' % self.xi[i])
            wfile.write('\n')
            for i in range(self.nyRoot + 1):
                wfile.write('%.9e\n' % self.yi[i])
            wfile.write('\n')
            for i in range(self.nzRoot + 1):
                wfile.write('%.9e\n' % self.zi[i])
            wfile.write('\n')
            wfile.write('\n')

            nRoot = self.nxRoot * self.nyRoot * self.nzRoot
            for i in range(nRoot):
                self._writeOcTreeNodeTypeRec(cellID=i, wfile=wfile)

        return

    def _writeOcTreeNodeTypeRec(self, cellID=None, wfile=None):
        """
        Recursive function to write the node type to file

        Parameters
        ----------

        cellID      : int
                      Tree index of the cell to be written

        wfile       : file
                      File object to write to
        """
        if self.isLeaf[cellID]:
            wfile.write("0\n")
        else:
            wfile.write("1\n")
            for i in range(self.childID[cellID].shape[0]):
                self._writeOcTreeNodeTypeRec(cellID=self.childID[cellID][i], wfile=wfile)

        return

    def readSpatialGrid(self, fname=''):
        """
        Reads the spatial grid from amr_grid.inp

        Parameters
        ----------

        fname : str, optional
                File name from which the spatial grid should be read. If omitted 'amr_grid.inp' will be used.

        """

        if fname == '':
            fname = 'amr_grid.inp'

        with open(fname, 'r') as rfile:
            print('Reading ' + fname)

            form = float(rfile.readline())
            dum = rfile.readline()
            grid_style = float(rfile.readline())
            if grid_style != 1:
                raise ValueError('Unsupported AMR style in ' + fname + '. Current only Octree AMR is supported.')

            crd_system = int(rfile.readline())
            if crd_system < 100:
                self.crd_sys = 'car'
                print("Reading cartesian grid")
            elif (crd_system >= 100) & (crd_system < 200):
                self.crd_sys = 'sph'
                print("Reading spherical grid")
            elif (crd_system >= 200) & (crd_system < 300):
                self.crd_sys = 'cyl'
                print("Reading cylindrical grid")
            else:
                raise ValueError('Unsupported coordinate system in ' + fname)

            grid_info = float(rfile.readline())

            dum = rfile.readline()

            dum = rfile.readline().split()
            self.act_dim = [int(dum[i]) for i in range(len(dum))]
            print("Active dimensions : ", self.act_dim[0], self.act_dim[1], self.act_dim[2])
            if 1 in self.act_dim:
                self.nChild = 2**(np.array(self.act_dim, dtype=int).sum())
            dum = rfile.readline().split()
            self.nxRoot, self.nyRoot, self.nzRoot = int(dum[0]), int(dum[1]), int(dum[2])
            print("Base grid size : ", self.nxRoot, self.nyRoot, self.nzRoot)
            dum = rfile.readline()

            dum = rfile.readline().split()
            levelMax, nLeaf, nCell = int(dum[0]), int(dum[1]), int(dum[2])
            nBranch = nCell - nLeaf

            print("Tree depth : ", levelMax)
            print("Nr of leaves : ", nLeaf)
            print("Nr of cells : ", nCell)
            dum = rfile.readline()

            self.levelMax = 0
            self.nLeaf = 0
            self.nBranch = 0

            self.x = np.zeros(nCell, dtype=np.float64)
            self.y = np.zeros(nCell, dtype=np.float64)
            self.z = np.zeros(nCell, dtype=np.float64)
            self.dx = np.zeros(nCell, dtype=np.float64)
            self.dy = np.zeros(nCell, dtype=np.float64)
            self.dz = np.zeros(nCell, dtype=np.float64)

            self.isLeaf = np.ones(nCell, dtype=bool)
            self.level = np.zeros(nCell, dtype=np.int_)
            self.parentID = np.zeros(nCell, dtype=np.int_)
            self.childID = []

            #
            # First of all read the base grid
            #
            self.xi = np.zeros(self.nxRoot + 1, dtype=np.float64)
            self.yi = np.zeros(self.nyRoot + 1, dtype=np.float64)
            self.zi = np.zeros(self.nzRoot + 1, dtype=np.float64)

            for i in range(self.nxRoot + 1):
                self.xi[i] = float(rfile.readline())
            dum = rfile.readline()
            for i in range(self.nyRoot + 1):
                self.yi[i] = float(rfile.readline())
            dum = rfile.readline()
            for i in range(self.nzRoot + 1):
                self.zi[i] = float(rfile.readline())

            dum = rfile.readline()
            dum = rfile.readline()

            self.xc = (self.xi[0:self.nxRoot] + self.xi[1:self.nxRoot + 1]) * 0.5
            self.yc = (self.yi[0:self.nyRoot] + self.yi[1:self.nyRoot + 1]) * 0.5
            self.zc = (self.zi[0:self.nzRoot] + self.zi[1:self.nzRoot + 1]) * 0.5

            dx = self.xi[1] - self.xi[0]
            dy = self.yi[1] - self.yi[0]
            dz = self.zi[1] - self.zi[0]

            self.cellIDCur = 0
            for iz in range(self.nzRoot):
                for iy in range(self.nyRoot):
                    for ix in range(self.nxRoot):
                        self.putNode(crd=(self.xc[ix], self.yc[iy], self.zc[iz]), cellsize=(dx, dy, dz), level=0,
                                     parentID=-1, cellID=self.cellIDCur)
                        self.cellIDCur += 1
                        self.nLeaf += 1

            nRoot = self.nxRoot * self.nyRoot * self.nzRoot

            for i in range(nRoot):
                self._readGridNodeTypeOcTreeRec(cellID=i, rfile=rfile)
        #
        # Now read which of the cells should be resolved
        #
        self.nCell = self.nLeaf + self.nBranch
        self.generateLeafID()

        self.selfCheck()
        return

    def _readGridNodeTypeOcTreeRec(self, cellID=None, rfile=None):
        """
        Recursive function to write the node type to file

        Parameters
        ----------

        cellID      : int
                      Tree index of the cell to be read

        rfile       : file
                      File object to read from

        """
        dum = int(rfile.readline())
        if dum == 1:

            self.nLeaf -= 1
            self.nBranch += 1
            self.nLeaf += self.nChild
            if self.levelMax < self.level[cellID] + 1:
                self.levelMax = self.level[cellID] + 1
                print('Tree depth : ', self.levelMax)

            #
            # Generate the cell center offsets for a proper octree
            #
            if self.nChild == 8:
                xc_offset = np.array([-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dx[cellID]
                yc_offset = np.array([-0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dy[cellID]
                zc_offset = np.array([-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5], dtype=np.float64) * self.dz[cellID]
            #
            # If we are resolving only two dimensions
            #
            elif self.nChild == 4:

                if self.act_dim[0] == 0:
                    xc_offset = np.array([0., 0., 0., 0.], dtype=np.float64)
                    yc_offset = np.array([-0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dy[cellID]
                    zc_offset = np.array([-0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dz[cellID]
                elif self.act_dim[1] == 0:
                    xc_offset = np.array([-0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dx[cellID]
                    yc_offset = np.array([0., 0., 0., 0.], dtype=np.float64)
                    zc_offset = np.array([-0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dz[cellID]
                elif self.act_dim[2] == 0:
                    xc_offset = np.array([-0.5, 0.5, -0.5, 0.5], dtype=np.float64) * self.dx[cellID]
                    yc_offset = np.array([-0.5, -0.5, 0.5, 0.5], dtype=np.float64) * self.dy[cellID]
                    zc_offset = np.array([0., 0., 0., 0.], dtype=np.float64)

            #
            # Generate the cell center offsets if only two dimensions are resolved
            #
            elif self.nChild == 2:

                if self.act_dim[0] == 1:
                    xc_offset = np.array([-0.5, 0.5], dtype=np.float64) * self.dx[cellID]
                    yc_offset = np.array([0., 0.], dtype=np.float64)
                    zc_offset = np.array([0., 0.], dtype=np.float64)

                if self.act_dim[1] == 1:
                    xc_offset = np.array([0., 0.], dtype=np.float64)
                    yc_offset = np.array([-0.5, 0.5], dtype=np.float64) * self.dy[cellID]
                    zc_offset = np.array([0., 0.], dtype=np.float64)

                if self.act_dim[2] == 1:
                    xc_offset = np.array([0., 0.], dtype=np.float64)
                    yc_offset = np.array([0., 0.], dtype=np.float64)
                    zc_offset = np.array([-0.5, 0.5], dtype=np.float64) * self.dz[cellID]

            else:
                raise ValueError('Wrong number of child/leaf cells. The number of leaf cells should be 2,4,8 in '
                                 + '1,2,3 dimensions, respectively. Insted ' + ("%d" % self.nChild) + ' have been '
                                 + 'claimed.')

            self.isLeaf[cellID] = False
            self.childID[cellID] = np.zeros(self.nChild, dtype=np.int_)
            for i in range(self.nChild):
                self.childID[cellID][i] = self.cellIDCur
                dx = self.dx[cellID] * (2.0 - self.act_dim[0])
                dy = self.dy[cellID] * (2.0 - self.act_dim[1])
                dz = self.dz[cellID] * (2.0 - self.act_dim[2])

                self.putNode(
                    crd=(self.x[cellID] + xc_offset[i], self.y[cellID] + yc_offset[i], self.z[cellID] + zc_offset[i]),
                    cellsize=(dx, dy, dz), level=self.level[cellID] + 1,
                    parentID=cellID, cellID=self.cellIDCur)
                self.cellIDCur += 1
                self._readGridNodeTypeOcTreeRec(cellID=self.childID[cellID][i], rfile=rfile)

        return

