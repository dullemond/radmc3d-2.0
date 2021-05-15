"""This module contains functions to set up a RADMC-3D model for dust and/or line simulations.
For help on the syntax or functionality of each function see the help of the individual functions
"""
from __future__ import absolute_import
from __future__ import print_function
import traceback
import inspect
import sys
import warnings
import importlib

try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use the python module of RADMC-3D you need to install Numpy')
    print(traceback.format_exc())

from . import natconst as nc
from . import analyze


class radmc3dModel(object):
    """

    Parameters
    ----------
     model           : str
                      Name of the model that should be used to create the density structure.
                      The file should be in a directory from where it can directly be imported
                      (i.e. the directory should be in the PYTHON_PATH environment variable or
                      it should be in the current working directory)
                      and the file name should be 'model_xxx.py', where xxx stands for the string
                      that should be specified in this variable

    binary          : bool, optional
                      If True input files will be written in binary format, if False input files are
                      written as formatted ascii text.

    old             : bool, optional
                      If set to True the input files for the old 2D version of radmc will be created

    dfunc           : function, optional
                      Decision function for octree-like amr tree building. It should take linear arrays of
                      cell centre coordinates (x,y,z) and cell half-widhts (dx,dy,dz) in all three dimensions,
                      a radmc3d model, a dictionary with all parameters from problem_params.inp and an other
                      keyword argument (**kwargs). It should return a boolean ndarray of the same length as
                      the input coordinates containing True if the cell should be resolved and False if not.
                      An example for the implementation of such decision function can be found in radmc3dPy.analyze
                      module (radmc3dPy.analyze.gdensMinMax()).

    dfpar           : dictionary
                      Dicionary of keyword arguments to be passed on to dfunc. These parameters will not be written
                      to problem_params.inp. Parameters can also be passed to dfunc via normal keyword arguments
                      gathered in **kwargs, however all keyword arguments in **kwargs will be written to
                      problem_params.inp

    parfile_update  : bool
                      If True the parameter file (problem_params.inp) will be updated / overwritten with any parameters
                      that are possibly passed as keyword arguments. For False the problem_params.inp file will not be
                      overwritten irrespectively of the parameter values the data members of the radmc3dModel instance
                      would contain.

    Keyword Arguments  :
                      Any varible name in problem_params.inp can be used as a keyword argument.
                      At first all variables are read from problem_params.in to a dictionary called ppar. Then
                      if there is any keyword argument set in the call of problem_setup_dust the ppar dictionary
                      is searched for this key. If found the value belonging to that key in the ppar dictionary
                      is changed to the value of the keyword argument. If no such key is found then the dictionary
                      is simply extended by the keyword argument. Finally the problem_params.inp file is updated
                      with the new parameter values.

    Attributes
    ----------

    binary          : bool
                     If True, it sets the data file format to binary (overwriting also the rto_style parameter!),
                     if False the file format will be formatted ascii

    data            : radmc3dData
                     Container for physical variables in the mode

    grid            : radmc3dGrid
                     Container for spatial and wavelength grid

    model           : str
                     Name of the model that has callable functions to generate physical variables for the model

    old             : bool
                     If True the model is meant to be for radmc, the pre-decessor code of radmc-3d.

    opac            : radmc3dOpac
                     Container for dust opacities

    par             : radmc3dPar
                     Container for the parameters of the model (i.e. the content of the problem_params.inp file)

    radsources      : radmc3dRadsources
                     Container for the radiation sources in the model
    """

    def __init__(self, model=None, binary=None, old=False, dfunc=None, dfpar=None, parfile_update=True, **kwargs):


        self.par = None

        # If kwargs set update the problem_params.inp file
        if kwargs:
            try:
                self.par = analyze.readParams()
            except:
                analyze.writeDefaultParfile(model=model)
                self.par = analyze.readParams()

            for ikey in kwargs.keys():
                self.par.ppar[ikey] = kwargs[ikey]

                if isinstance(kwargs[ikey], float):
                    self.par.setPar([ikey, ("%.7e" % kwargs[ikey]), '', ''])
                elif isinstance(kwargs[ikey], int):
                    self.par.setPar([ikey, ("%d" % kwargs[ikey]), '', ''])
                elif isinstance(kwargs[ikey], str):
                    self.par.setPar([ikey, kwargs[ikey], '', ''])
                elif isinstance(kwargs[ikey], list):
                    dum = '['
                    for i in range(len(kwargs[ikey])):
                        if type(kwargs[ikey][i]) is float:
                            dum = dum + ("%.7e" % kwargs[ikey][i])
                        elif type(kwargs[ikey][i]) is int:
                            dum = dum + ("%d" % kwargs[ikey][i])
                        elif type(kwargs[ikey][i]) is str:
                            dum = dum + (kwargs[ikey][i])
                        else:
                            raise TypeError('Unknown data type in ' + ikey)

                        if i < len(kwargs[ikey]) - 1:
                            dum = dum + ', '
                    dum = dum + ']'
                    self.par.setPar([ikey, dum, '', ''])

        # If a binary setup is used update also the radmc3d output file format to binary
        if binary is not None:
            if binary:
                if self.par is None:
                    try:
                        self.par = analyze.readParams()
                    except:
                        analyze.writeDefaultParfile(model=model)
                        self.par = analyze.readParams()
                    self.par.setPar(['rto_style', '1', '  Format of outpuf files (1-ascii, 2-unformatted f77, 3-binary',
                                     'Code parameters'])

        if parfile_update:
            if self.par is not None:
                self.par.writeParfile()
            else:
                 warnings.warn("The parameter file (problem_params.inp) has not been read. Nothing to update",
                               RuntimeWarning)

        if model is not None:
            self.model = model
        else:
            self.model = None

        if binary is not None:
            self.binary = binary
        else:
            self.binary = True

        if old is not None:
            self.old = old
        else:
            self.old = False

        self.data = None
        self.grid = None
        self.radsources = None
        self.opac = None
        self.dfunc = dfunc
        self.dfpar = dfpar

    def readParams(self):

        self.par = analyze.readParams()

        # Check the grid style. If it's 1 (i.e. AMR is active) make sure old is set to False
        if (self.par.ppar['grid_style'] == 1) & (self.old == True):
            print('WARNING')
            print('Grid style in problem_params.inp is set to 1, i.e. AMR is active, but the model appears to have')
            print('the "old" keyword set to True, meaning the model is meant for RADMC, the predecessor code of ')
            print('RADMC3D. Thus for now, the old style is switched to False.')

            self.old = False


    def makeVar(self, ddens=False, dtemp=False, gdens=False, gtemp=False, gvel=False, vturb=False,
                     writeToFile=False):
        """
        Generates variables and possibly writes them to file

        Parameters
        ----------
        ddens           : bool, optional
                          If True generates the dust density

        dtemp           : bool, optional
                          If True generates the dust temperature (normally this would be calculated by RADMC-3D, but
                            in some cases, e.g. debugging, it might come handy)

        gdens           : bool, optional
                          If True generates the gas density

        gtemp           : bool, optional
                          If True generates the gas temperature

        gvel            : bool, optional
                          If True generates the gas velocity

        vturb           : bool, optional
                          If True generates the microturbulent velocity field

        writeToFile     : bool, optional
                          If True the generated variable(s) will be written to file

        """

        #
        # Get the model
        #
        try:
            mdl = importlib.import_module(self.model)
        except ImportError:
            try:
                mdl = importlib.import_module('radmc3dPy.models.' + self.model)
            except ImportError:
                print(' ' + self.model + '.py could not be imported')
                print(' The model files should either be in the current working directory or')
                print(' in the radmc3d python module directory')
                return
        # --------------------------------------------------------------------------------------------
        # Create the dust density distribution
        # --------------------------------------------------------------------------------------------
        if ddens:
            if dir(mdl).__contains__('getDustDensity'):
                if callable(getattr(mdl, 'getDustDensity')):
                    if self.data is None:
                        self.data = analyze.radmc3dData(self.grid)

                    self.data.rhodust = mdl.getDustDensity(grid=self.grid, ppar=self.par.ppar)
                else:
                    raise RuntimeError('Dust density cannot be calculated in model ' + self.model)
            else:
                raise RuntimeError(' ' + self.model + '.py does not contain a getDustDensity() function, therefore, '
                                   + ' dust_density.inp cannot be written')

            if writeToFile:
                self.data.octree = (self.par.ppar['grid_style'] == 1)
                self.data.writeDustDens(binary=self.binary, old=self.old)
        # --------------------------------------------------------------------------------------------
        # Create the dust temperature distribution if the model has such function
        # --------------------------------------------------------------------------------------------
        if dtemp:
            if dir(mdl).__contains__('getDustTemperature'):
                if callable(getattr(mdl, 'getDustTemperature')):
                    if self.data is None:
                        self.data = analyze.radmc3dData(self.grid)
                    self.data.dusttemp = mdl.getDustTemperature(grid=self.grid, ppar=self.par.ppar)
                else:
                    raise RuntimeError('Dust temperature cannot be calculated in model ' + self.model + ' but it was '
                                       + ' requested to be written')
            else:
                raise RuntimeError(' ' + self.model + '.py does not contain a getDustTemperature() function, therefore,'
                                   + ' dust_temperature.inp cannot be written')

            if writeToFile:
                self.data.octree = (self.par.ppar['grid_style'] == 1)
                self.data.writeDustTemp(binary=self.binary)

        # --------------------------------------------------------------------------------------------
        # Create the molecular abundance
        # --------------------------------------------------------------------------------------------
        if gdens:
            if dir(mdl).__contains__('getGasDensity') & dir(mdl).__contains__('getGasAbundance'):
                if callable(getattr(mdl, 'getGasDensity')) & callable(getattr(mdl, 'getGasAbundance')):
                    for imol in range(len(self.par.ppar['gasspec_mol_name'])):
                        self.data.rhogas = mdl.getGasDensity(grid=self.grid, ppar=self.par.ppar)
                        gasabun = mdl.getGasAbundance(grid=self.grid, ppar=self.par.ppar,
                                                      ispec=self.par.ppar['gasspec_mol_name'][imol])
                        self.data.ndens_mol = self.data.rhogas / (2.4 * nc.mp) * gasabun

                    if abs(self.par.ppar['lines_mode']) > 2:
                        for icp in range(len(self.par.ppar['gasspec_colpart_name'])):
                            self.data.rhogas = mdl.getGasDensity(grid=self.grid, ppar=self.par.ppar)
                            gasabun = mdl.getGasAbundance(grid=self.grid, ppar=self.par.ppar,
                                                          ispec=self.par.ppar['gasspec_colpart_name'][icp])
                            self.data.ndens_mol = self.data.rhogas / (2.4 * nc.mp) * gasabun

            else:
                raise RuntimeError(' ' + self.model + '.py does not contain a getGasAbundance() function, therefore, '
                                   + ' numberdens_***.inp cannot be written')

            if writeToFile:
                # Write the gas density
                self.data.octree = (self.par.ppar['grid_style'] == 1)
                for imol in range(len(self.par.ppar['gasspec_mol_name'])):
                    self.data.writeGasDens(ispec=self.par.ppar['gasspec_mol_name'][imol], binary=self.binary)

                if abs(self.par.ppar['lines_mode']) > 2:
                    for icp in range(len(self.par.ppar['gasspec_colpart_name'])):
                        self.data.writeGasDens(ispec=self.par.ppar['gasspec_colpart_name'][icp], binary=self.binary)

        # --------------------------------------------------------------------------------------------
        # Get the gas velocity field
        # --------------------------------------------------------------------------------------------
        if gvel:
            if dir(mdl).__contains__('getVelocity'):
                if callable(getattr(mdl, 'getVelocity')):
                    self.data.gasvel = mdl.getVelocity(grid=self.grid, ppar=self.par.ppar)
            else:
                raise RuntimeError(' ' + self.model + '.py does not contain a getVelocity() function, therefore, '
                                   + ' gas_velocity.inp cannot be written')

            if writeToFile:
                # Write the gas velocity
                self.data.octree = (self.par.ppar['grid_style'] == 1)
                self.data.writeGasVel(binary=self.binary)
        # --------------------------------------------------------------------------------------------
        # Get the kinetik gas temperature
        # --------------------------------------------------------------------------------------------
        # Write the gas temperature if specified
        if gtemp:
            if dir(mdl).__contains__('getGasTemperature'):
                if callable(getattr(mdl, 'getGasTemperature')):
                    self.data.gastemp = mdl.getGasTemperature(grid=self.grid, ppar=self.par.ppar)
            else:
                raise RuntimeError(' ' + self.model + '.py does not contain a getGasTemperature() function, therefore,'
                                   + ' gas_temperature.inp cannot be written')

            if writeToFile:
                # Write the gas temperature
                self.data.octree = (self.par.ppar['grid_style'] == 1)
                self.data.writeGasTemp(binary=self.binary)
        # --------------------------------------------------------------------------------------------
        # Get the turbulent velocity field
        # --------------------------------------------------------------------------------------------
        if vturb:
            if dir(mdl).__contains__('getVTurb'):
                if callable(getattr(mdl, 'getVTurb')):
                    self.data.vturb = mdl.getVTurb(grid=self.grid, ppar=self.par.ppar)
            else:
                print(' ' + self.model + '.py does not contain a getVTurb() function, therefore, zero microturbulent '
                      + ' velocity will be assumed everywhere in the model.')
                if self.par.ppar['grid_style'] == 1:
                    self.data.vturb = np.zeros([self.grid.x.nLeaf], dtype=np.float64)
                    self.data.vturb[:] = 0.
                else:
                    self.data.vturb = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=np.float64)
                    self.data.vturb[:, :, :] = 0.

            if writeToFile:
                self.data.octree = (self.par.ppar['grid_style'] == 1)
                self.data.writeVTurb(binary=self.binary)

    def makeGrid(self, sgrid=True, wgrid=True, writeToFile=False, **kwargs):
        """
        Generates a spatial and/or wavelength grid

        Parameters
        ----------
        wgrid           : bool
                          Set to True to generate the wavelength grid

        sgrid           : bool
                          Set to True to generate the spatial grid

        writeToFile     : bool
                          If True the grid will be written to amr_grid.inp and/or wavelength_micron.inp

        **kwargs        : Any varible name in problem_params.inp can be used as a keyword argument.
                          At first all variables are read from problem_params.in to a dictionary called ppar. Then
                          if there is any keyword argument set in the call of problem_setup_dust the ppar dictionary
                          is searched for this key. If found the value belonging to that key in the ppar dictionary
                          is changed to the value of the keyword argument. If no such key is found then the dictionary
                          is simply extended by the keyword argument. Finally the problem_params.inp file is updated
                          with the new parameter values.

        """

        self.grid = None
        #
        # If we generate the spatial grid
        #
        if sgrid:
            #
            # Check if AMR is activated or not
            #
            if self.par.ppar['grid_style'] == 1:
                self.grid = analyze.radmc3dOctree()

                # Pass all parameters from dfpar to ppar
                if self.dfpar is not None:
                    for ikey in self.dfpar.keys():
                        self.par.ppar[ikey] = self.dfpar[ikey]

                # Spatial grid
                self.grid.makeSpatialGrid(ppar=self.par.ppar, dfunc=self.dfunc, model=self.model, **kwargs)
            else:
                self.grid = analyze.radmc3dGrid()
                # Spatial grid
                self.grid.makeSpatialGrid(ppar=self.par.ppar)

        #
        # If we generate the wavelength grid
        #
        if wgrid:
            if self.grid is None:
                if self.par.ppar['grid_style'] == 1:
                    self.grid = analyze.radmc3dOctree()
                else:
                    self.grid = analyze.radmc3dGrid()

            # Wavelength grid
            self.grid.makeWavelengthGrid(ppar=self.par.ppar)

        #
        # Now write out the grid into file
        #
        if writeToFile:
            #
            # For AMR grids there is no 'old' option as radmc the predecessor code could not handle AMR-style grid
            # refinement
            #
            if self.par.ppar['grid_style'] == 1:
                # Frequency grid
                self.grid.writeWavelengthGrid(old=False)
                # Spatial grid
                self.grid.writeSpatialGrid(old=False)
            else:
                # Frequency grid
                self.grid.writeWavelengthGrid(old=self.old)
                # Spatial grid
                self.grid.writeSpatialGrid(old=self.old)

    def makeRadSources(self, writeToFile=True):
        """

        Parameters
        ----------
        writeToFile     : bool
                          If True the radiation will be written to stars.inp and or stellarsrc_density.inp,
                          stellarsrc_templates.inp

        """
        try:
            mdl = importlib.import_module(self.model)
        except ImportError:
            try:
                mdl = importlib.import_module('radmc3dPy.models.' + self.model)
            except ImportError:
                print(' ' + self.model + '.py could not be imported')
                print(' The model files should either be in the current working directory or')
                print(' in the radmc3d python module directory')
                return

        self.radsources = analyze.radmc3dRadSources(ppar=self.par.ppar, grid=self.grid)
        self.radsources.getStarSpectrum(tstar=self.par.ppar['tstar'], rstar=self.par.ppar['rstar'], grid=self.grid)

        # Check if the model has functions to set up continuous starlike sources
        if 'incl_cont_stellarsrc' in self.par.ppar:
            stellarsrcEnabled = self.par.ppar['incl_cont_stellarsrc']

            if stellarsrcEnabled:
                if dir(mdl).__contains__('getStellarsrcDensity'):
                    if callable(getattr(mdl, 'getStellarsrcDensity')):
                        if dir(mdl).__contains__('getStellarsrcTemplates'):
                            if callable(getattr(mdl, 'getStellarsrcTemplates')):
                                stellarsrcEnabled = True
                            else:
                                stellarsrcEnabled = False
                        else:
                            stellarsrcEnabled = False
                    else:
                        stellarsrcEnabled = False
                else:
                    stellarsrcEnabled = False
        else:
            stellarsrcEnabled = False

        if stellarsrcEnabled:
            dum = mdl.getStellarsrcTemplates(grid=self.grid, ppar=self.par.ppar)
            if dum[0, 0] <= 0.:
                self.radsources.csntemplate = dum.shape[0]
                self.radsources.cstemp = []
                self.radsources.cstemptype = 1
                self.radsources.cststar = dum[:, 0]
                self.radsources.csrstar = dum[:, 1]
                self.radsources.csmstar = dum[:, 2]
            else:
                self.radsources.csntemplate = dum.shape[0]
                self.radsources.cstemp = dum
                self.radsources.cstemptype = 2
                self.radsources.cststar = []
                self.radsources.csrstar = []
                self.radsources.csmstar = []

            self.radsources.csdens = mdl.getStellarsrcDensity(grid=self.grid, ppar=self.par.ppar)

        if writeToFile:
            # Input radiation field (discrete stars)
            self.radsources.writeStarsinp(ppar=self.par.ppar, old=self.old, wav=self.grid.wav)

            # Continuous starlike sources
            if stellarsrcEnabled:
                self.radsources.writeStellarsrcTemplates()
                self.radsources.writeStellarsrcDensity(binary=self.binary)

            # totlum = radSources.getTotalLuminosities()
            if not self.old:
                print('-------------------------------------------------------------')
                print('Luminosities of radiation sources in the model :')

                totlum = self.radsources.getTotalLuminosities(readInput=True)
                print('As calculated from the input files :')
                print('Stars : ')
                print("  Star #%d + hotspot        : %.6e" % (0, totlum['lnu_star'][0]))
                for istar in range(1, self.radsources.nstar):
                    print("  Star #%d               : %.6e" % (istar, totlum['lnu_star'][istar]))
                print("Continuous starlike source : %.6e" % totlum['lnu_accdisk'])
                print(' ')
                print('-------------------------------------------------------------')

    def makeDustOpac(self):
        """
        Generates dust opacities and writes them to file
        """

        if 'dustkappa_ext' in self.par.ppar:
            self.opac = analyze.radmc3dDustOpac()
            # Master dust opacity file
            self.opac.writeMasterOpac(ext=self.par.ppar['dustkappa_ext'],
                                      scattering_mode_max=self.par.ppar['scattering_mode_max'], old=self.old)
            if self.old:
                if self.grid is None:
                    raise ValueError('Unknown wavelength grid. For old style dust opacity the wavelength grid '
                                     + 'must be set.')
                # Frequency grid
                self.grid.writeWavelengthGrid(old=self.old)
                # Create the dust opacity files
                self.opac.makeopacRadmc2D(ext=self.par.ppar['dustkappa_ext'])
        else:
            # if old:
            # print 'ERROR'
            # print 'Calculating dust opacities for radmc (2D version) is not yet implemented)'
            self.opac = analyze.radmc3dDustOpac()
            # Calculate the opacities and write the master opacity file
            self.opac.makeOpac(ppar=self.par.ppar, old=self.old)

    def writeRadmc3dInp(self):
        """
        Writes the radmc3d.inp master command file for RADMC-3D
        """
        writeRadmc3dInp(self.par)

    def writeLinesInp(self):
        """
        Writes the lines.inp master command file for line simulation in RADMC-3D
        """

        writeLinesInp(ppar=self.par.ppar)


    def setupDust(self, sgrid=True, wgrid=True, radsources=True, dustopac=True, ddens=True, dtemp=False):
        """
        Set up a model for RADMC-3D.

        Parameters
        ----------
        sgrid           : bool
                          Should the spatial grid be generated? (default=True)

        wgrid           : bool
                          Should the wavelength grid be generated? (default=True)

        radsources      : bool
                          Should the radiation sources be generated? (default=True)

        dustopac        : bool
                          Should dust opacities be generated? (default=True)

        ddens           : bool
                          Should the dust density be generated? (default=True)

        dtemp           : bool
                          Should the dust temperature be generated (default=False)

        Note, that by default all generated data will be written to files immediately. In case you wish to generate
            the files but not write the data to files but only keep them instead in data attributes, please use the
            individual generator functions methods (e.g. makeGrid, makeRadSources, makeVar, etc)
        """

        self.writeRadmc3dInp()

        self.makeGrid(sgrid=sgrid, wgrid=wgrid, writeToFile=True)

        if radsources:
            self.makeRadSources(writeToFile=True)

        if dustopac:
            self.makeDustOpac()

        if ddens:
            self.makeVar(ddens=True, writeToFile=True)

        if dtemp:
            self.makeVar(dtemp=True, writeToFile=True)

    def setupGas(self, sgrid=False, wgrid=False, radsources=False, gdens=True, gtemp=True, gvel=True, vturb=True):
            """
            Set up a model for RADMC-3D.

            Parameters
            ----------
            sgrid           : bool
                              Should the spatial grid be generated? (default=False)

            wgrid           : bool
                              Should the wavelength grid be generated? (default=False)

            radsources      : bool
                              Should the radiation sources be generated? (default=False)

            gdens           : bool
                              Should the gas density be generated? (default=True)

            gtemp           : bool
                              Should the gas temperature be generated (default=False)

            gvel            : bool
                              Should the gas veloctiy be generated (default=False)

            vturb           : bool
                              Should the turbulent velocity field be generated (default=False)

            Note, that by default all generated data will be written to files immediately. In case you wish to generate
                the files but not write the data to files but only keep them instead in data attributes, please use the
                individual generator functions methods (e.g. makeGrid, makeRadSources, makeVar, etc)
            """

            self.writeRadmc3dInp()

            self.makeGrid(sgrid=sgrid, wgrid=wgrid, writeToFile=True)

            if radsources:
                self.makeRadSources(writeToFile=True)

            if gdens:
                self.makeVar(gdens=True, writeToFile=True)

            if gtemp:
                self.makeVar(gtemp=True, writeToFile=True)

            if gvel:
                self.makeVar(gvel=True, writeToFile=True)

            if vturb:
                self.makeVar(vturb=True, writeToFile=True)



def problemSetupDust(model=None, binary=True, writeDustTemp=False, old=False, dfunc=None, dfpar=None, **kwargs):
    """
    Function to set up a dust model for RADMC-3D 

    Parameters
    ----------
    model           : str
                      Name of the model that should be used to create the density structure.
                      The file should be in a directory from where it can directly be imported 
                      (i.e. the directory should be in the PYTHON_PATH environment variable or
                      it should be in the current working directory)
                      and the file name should be 'model_xxx.py', where xxx stands for the string
                      that should be specified in this variable

    binary          : bool, optional
                      If True input files will be written in binary format, if False input files are
                      written as formatted ascii text. 

    writeDustTemp   : bool, optional
                      If True a separate dust_temperature.inp/dust_tempearture.binp file will be
                      written under the condition that the model contains a function getDustTemperature() 

    old             : bool, optional
                      If set to True the input files for the old 2D version of radmc will be created

    dfunc           : function, optional
                      Decision function for octree-like amr tree building. It should take linear arrays of 
                      cell centre coordinates (x,y,z) and cell half-widhts (dx,dy,dz) in all three dimensions,
                      a radmc3d model, a dictionary with all parameters from problem_params.inp and an other 
                      keyword argument (**kwargs). It should return a boolean ndarray of the same length as 
                      the input coordinates containing True if the cell should be resolved and False if not. 
                      An example for the implementation of such decision function can be found in radmc3dPy.analyze
                      module (radmc3dPy.analyze.gdensMinMax()).

    dfpar           : dictionary
                      Dicionary of keyword arguments to be passed on to dfunc. These parameters will not be written
                      to problem_params.inp. Parameters can also be passed to dfunc via normal keyword arguments 
                      gathered in **kwargs, however all keyword arguments in **kwargs will be written to 
                      problem_params.inp

    **kwargs        : Any varible name in problem_params.inp can be used as a keyword argument.
                      At first all variables are read from problem_params.in to a dictionary called ppar. Then 
                      if there is any keyword argument set in the call of problem_setup_dust the ppar dictionary 
                      is searched for this key. If found the value belonging to that key in the ppar dictionary 
                      is changed to the value of the keyword argument. If no such key is found then the dictionary 
                      is simply extended by the keyword argument. Finally the problem_params.inp file is updated
                      with the new parameter values.


    Notes
    -----
    Files written by problemSetupDust() for RADMC-3D

        * dustopac.inp             : Dust opacity master file.

        * wavelength_micron.inp    : Wavelength grid.

        * amr_grid.inp             : Spatial grid.

        * stars.inp                : Input radiation field (discrete stellar sources).

        * stellarsrc_density.inp   : Input radiation field (continuous stellar sources).

        * stellarsrc_templates.inp : Input radiation field (continuous stellar sources).

        * dust_density.inp         : Dust density distribution.

        * radmc3d.inp              : Parameters for RADMC-3D (e.g. Nr of photons to be used, scattering type, etc).

    """

    if model is None:
        raise ValueError('Unknown model. No model name is given')

    # Read the parameters from the problem_params.inp file
    modpar = analyze.readParams()
    # Make a local copy of the ppar dictionary
    ppar = modpar.ppar
    if ppar is None:
        raise RuntimeError('Unknown ppar. For some reason the model parameters could not be read from '
                           'problem_params.inp was not found')

    if ppar['grid_style'] != 0:
        if old:
            raise ValueError('grid_style is set to ', ppar['grid_style'], ' which is is not compatible with the old'
                             + ' 2D version of RADMC-3D as it did not support AMR style grids')

    # --------------------------------------------------------------------------------------------
    # If there is any additional keyword argument (**kwargs) then check
    #   if there is such key in the ppar dictionary and if is change its value that of
    #   the keyword argument. If there is no such key in the ppar dictionary then add the keyword
    #   to the dictionary
    # --------------------------------------------------------------------------------------------
    if binary:
        modpar.setPar(['rto_style', '3', '', ''])

    if kwargs:
        for ikey in kwargs.keys():
            modpar.ppar[ikey] = kwargs[ikey]

            if isinstance(kwargs[ikey], float):
                modpar.setPar([ikey, ("%.7e" % kwargs[ikey]), '', ''])
            elif isinstance(kwargs[ikey], int):
                modpar.setPar([ikey, ("%d" % kwargs[ikey]), '', ''])
            elif isinstance(kwargs[ikey], str):
                modpar.setPar([ikey, kwargs[ikey], '', ''])
            elif isinstance(kwargs[ikey], list):
                dum = '['
                for i in range(len(kwargs[ikey])):
                    if isinstance(kwargs[ikey][i], float):
                        dum = dum + ("%.7e" % kwargs[ikey][i])
                    elif isinstance(kwargs[ikey][i], int):
                        dum = dum + ("%d" % kwargs[ikey][i])
                    elif isinstance(kwargs[ikey][i], str):
                        dum = dum + (kwargs[ikey][i])
                    else:
                        raise TypeError('Unknown data type in ' + ikey)

                    if i < len(kwargs[ikey]) - 1:
                        dum = dum + ', '
                dum = dum + ']'
                modpar.setPar([ikey, dum, '', ''])

        modpar.writeParfile()
        ppar = modpar.ppar

    # --------------------------------------------------------------------------------------------
    # Create the grid
    # --------------------------------------------------------------------------------------------
    #
    # Check if AMR is activated or not
    #
    if ppar['grid_style'] == 1:
        grid = analyze.radmc3dOctree()

        # Pass all parameters from dfpar to ppar
        if dfpar is not None:
            for ikey in dfpar.keys():
                ppar[ikey] = dfpar[ikey]

        # Spatial grid
        grid.makeSpatialGrid(ppar=ppar, dfunc=dfunc, model=model, **kwargs)
    else:
        grid = analyze.radmc3dGrid()
        # Spatial grid
        grid.makeSpatialGrid(ppar=ppar)
    # Wavelength grid
    grid.makeWavelengthGrid(ppar=ppar)

    # --------------------------------------------------------------------------------------------
    # Dust opacity
    # --------------------------------------------------------------------------------------------
    if 'dustkappa_ext' in ppar:
        opac = analyze.radmc3dDustOpac()
        # Master dust opacity file
        opac.writeMasterOpac(ext=ppar['dustkappa_ext'], scattering_mode_max=ppar['scattering_mode_max'], old=old)
        if old:
            # Frequency grid
            grid.writeWavelengthGrid(old=old)
            # Create the dust opacity files
            opac.makeopacRadmc2D(ext=ppar['dustkappa_ext'])
    else:
        # if old:
        # print 'ERROR'
        # print 'Calculating dust opacities for radmc (2D version) is not yet implemented)'
        opac = analyze.radmc3dDustOpac()
        # Calculate the opacities and write the master opacity file
        opac.makeOpac(ppar=ppar, old=old)

    # --------------------------------------------------------------------------------------------
    # Try to get the specified model
    # --------------------------------------------------------------------------------------------
    try:
        mdl = importlib.import_module(model)
    except ImportError:
        try:
            mdl = importlib.import_module('radmc3dPy.models.' + model)
        except ImportError:
            print(' ' + model + '.py could not be imported')
            print(' The model files should either be in the current working directory or')
            print(' in the radmc3d python module directory')
            return
    # --------------------------------------------------------------------------------------------
    # Create the input radiation field (stars at this point)
    # --------------------------------------------------------------------------------------------
    radSources = analyze.radmc3dRadSources(ppar=ppar, grid=grid)
    radSources.getStarSpectrum(tstar=ppar['tstar'], rstar=ppar['rstar'])

    # Check if the model has functions to set up continuous starlike sources
    if 'incl_cont_stellarsrc' in ppar:
        stellarsrcEnabled = ppar['incl_cont_stellarsrc']

        if stellarsrcEnabled:
            if dir(mdl).__contains__('getStellarsrcDensity'):
                if callable(getattr(mdl, 'getStellarsrcDensity')):
                    if dir(mdl).__contains__('getStellarsrcTemplates'):
                        if callable(getattr(mdl, 'getStellarsrcTemplates')):
                            stellarsrcEnabled = True
                        else:
                            stellarsrcEnabled = False
                    else:
                        stellarsrcEnabled = False
                else:
                    stellarsrcEnabled = False
            else:
                stellarsrcEnabled = False
    else:
        stellarsrcEnabled = False

    if stellarsrcEnabled:
        dum = mdl.getStellarsrcTemplates(grid=grid, ppar=ppar)
        if dum[0, 0] <= 0.:
            radSources.csntemplate = dum.shape[0]
            radSources.cstemp = []
            radSources.cstemptype = 1
            radSources.cststar = dum[:, 0]
            radSources.csrstar = dum[:, 1]
            radSources.csmstar = dum[:, 2]
        else:
            radSources.csntemplate = dum.shape[0]
            radSources.cstemp = dum
            radSources.cstemptype = 2
            radSources.cststar = []
            radSources.csrstar = []
            radSources.csmstar = []

        radSources.csdens = mdl.getStellarsrcDensity(grid=grid, ppar=ppar)
    # --------------------------------------------------------------------------------------------
    # Create the dust density distribution
    # --------------------------------------------------------------------------------------------
    data = analyze.radmc3dData(grid)
    if dir(mdl).__contains__('getDustDensity'):
        if callable(getattr(mdl, 'getDustDensity')):
            data.rhodust = mdl.getDustDensity(grid=grid, ppar=ppar)
        else:
            raise RuntimeError('Dust density cannot be calculated in model ' + model)
    else:
        raise RuntimeError(' ' + model + '.py does not contain a getDustDensity() function, therefore, '
                           + ' dust_density.inp cannot be written')
    # --------------------------------------------------------------------------------------------
    # Create the dust temperature distribution if the model has such function
    # --------------------------------------------------------------------------------------------
    if writeDustTemp:
        if dir(mdl).__contains__('getDustTemperature'):
            if callable(getattr(mdl, 'getDustTemperature')):
                data.dusttemp = mdl.getDustTemperature(grid=grid, ppar=ppar)
            else:
                raise RuntimeError('Dust temperature cannot be calculated in model ' + model + ' but it was '
                                   + ' requested to be written')
        else:
            raise RuntimeError(' ' + model + '.py does not contain a getDustTemperature() function, therefore, '
                               + ' dust_temperature.inp cannot be written')

    # --------------------------------------------------------------------------------------------
    # Now write out everything
    # --------------------------------------------------------------------------------------------

    if ppar['grid_style'] == 1:
        # Frequency grid
        grid.writeWavelengthGrid()
        # Spatial grid
        grid.writeSpatialGrid()

    else:
        # Frequency grid
        grid.writeWavelengthGrid(old=old)
        # Spatial grid
        grid.writeSpatialGrid(old=old)

    # Input radiation field (discrete stars)
    radSources.writeStarsinp(ppar=ppar, old=old, wav=grid.wav)

    # Continuous starlike sources
    if stellarsrcEnabled:
        radSources.writeStellarsrcTemplates()
        radSources.writeStellarsrcDensity(binary=binary)

    # totlum = radSources.getTotalLuminosities()
    if not old:
        print('-------------------------------------------------------------')
        print('Luminosities of radiation sources in the model :')

        totlum = radSources.getTotalLuminosities(readInput=True)
        print('As calculated from the input files :')
        print('Stars : ')
        print("  Star #%d + hotspot        : %.6e" % (0, totlum['lnu_star'][0]))
        for istar in range(1, radSources.nstar):
            print("  Star #%d               : %.6e" % (istar, totlum['lnu_star'][istar]))
        print("Continuous starlike source : %.6e" % totlum['lnu_accdisk'])
        print(' ')
        print('-------------------------------------------------------------')

    # Dust density distribution
    data.octree = (ppar['grid_style'] == 1)
    data.writeDustDens(binary=binary, old=old)

    # Dust temperature distribution
    if writeDustTemp:
        data.octree = (ppar['grid_style'] == 1)
        data.writeDustTemp(binary=binary)
    # radmc3d.inp
    if not old:
        writeRadmc3dInp(modpar=modpar)
    else:
        writeRadmcInp(modpar=modpar)

    # --------------------------------------------------------------------------------------------
    # Calculate optical depth for diagnostics purposes
    # --------------------------------------------------------------------------------------------
    # radSources = radmc3dRadSources(ppar=ppar, grid=grid)
    # radSources.readStarsinp()
    # pwav = radSources.findPeakStarspec()[0]
    # data.getTau(wav=pwav, usedkappa=False)
    # print 'Radial optical depth at '+("%.2f"%pwav)+'um : ', data.taux.max()


def problemSetupGas(model=None, fullsetup=False, binary=True, writeGasTemp=False, dfunc=None, dfpar=None, **kwargs):
    """
    Function to set up a gas model for RADMC-3D 

    Parameters
    ----------

    model           : str
                      Name of the model that should be used to create the density structure
                      the file should be in a directory from where it can directly be imported 
                      (i.e. the directory should be in the PYTHON_PATH environment variable, or
                      it should be the current working directory)
                      and the file name should be 'MODELNAME.py', where MODELNAME stands for the string
                      that should be specified in this variable

    fullsetup       : bool, optional
                      If False only the files related to the gas simulation is written out
                      (i.e. no grid, stellar parameter file and radmc3d master command file is written)
                      assuming that these files have already been created for a previous continuum simulation.
                      If True the spatial and wavelength grid as well as the input radiation field
                      and the radmc3d master command file will be (over)written. 

    binary          : bool, optional
                      If True input files will be written in binary format, if False input files are
                      written as formatted ascii text. 

    writeGasTemp    : bool, optional
                      If True a separate gas_temperature.inp/gas_tempearture.binp file will be
                      written under the condition that the model contains a function get_gas_temperature() 

    dfunc           : function, optional
                      Decision function for octree-like amr tree building. It should take linear arrays of 
                      cell centre coordinates (x,y,z) and cell half-widhts (dx,dy,dz) in all three dimensions,
                      a radmc3d model, a dictionary with all parameters from problem_params.inp and an other 
                      keyword argument (**kwargs). It should return a boolean ndarray of the same length as 
                      the input coordinates containing True if the cell should be resolved and False if not. 
                      An example for the implementation of such decision function can be found in radmc3dPy.analyze
                      module (radmc3dPy.analyze.gdensMinMax()). 

    dfpar           : dictionary
                      Dicionary of keyword arguments to be passed on to dfunc. These parameters will not be written
                      to problem_params.inp. Parameters can also be passed to dfunc via normal keyword arguments 
                      gathered in **kwargs, however all keyword arguments in **kwargs will be written to 
                      problem_params.inp

    **kwargs        : Any varible name in problem_params.inp can be used as a keyword argument.
                      At first all variables are read from problem_params.in to a dictionary called ppar. Then 
                      if there is any keyword argument set in the call of problem_setup_gas the ppar dictionary 
                      is searched for such key. If found the value belonging to that key in the ppar dictionary 
                      is changed to the value of the keyword argument. If no such key is found then the dictionary 
                      is simply extended by the keyword argument. Finally the problem_params.inp file is updated
                      with the new parameter values.
                      Any additional keyword argument for the octree AMR mesh generation should also be passed here.


    Notes
    -----

    Files written by problemSetupGas()


        * lines.inp             : Line mode master command file.

        * numberdens_xxx.inp    : Number density of molecule/atomic species 'xxx'

        * gas_velocity.inp      : Gas velocity

        * microturbulence.inp   : The standard deviation of the Gaussian line profile caused by turbulent 
                                broadening.

        * gas_temperature.inp   : Gas temperature (which may be different from the dust temperature). If
                                tgas_eq_tdust is set to zero in radmc3d.inp the gas temperature in this
                                file will be used instead of the dust temperature. 

        If fullsetup is set to True the following additional files will be created

        * amr_grid.inp          : Spatial grid.

        * wavelength_micron.inp : Wavelength grid.

        * stars.inp             : Input radiation field.

        * radmc3d.inp           : Parameters for RADMC-3D (e.g. Nr of photons to be used, scattering type, etc).

    """

    if model is None:
        raise ValueError('Unknown model. No model name is given')

    # Read the parameters from the problem_params.inp file
    modpar = analyze.readParams()
    # Make a local copy of the ppar dictionary
    ppar = modpar.ppar
    if ppar is None:
        raise RuntimeError('Unknown ppar. For some reason the model parameters could not be read from '
                           'problem_params.inp was not found')

    # --------------------------------------------------------------------------------------------
    # If there is any additional keyword argument (**kwargs) then check
    #   if there is such key in the ppar dictionary and if is change its value that of
    #   the keyword argument. If there is no such key in the ppar dictionary then add the keyword
    #   to the dictionary
    # --------------------------------------------------------------------------------------------
    if binary:
        modpar.setPar(['rto_style', '3', '', ''])

    if kwargs:
        for ikey in kwargs.keys():
            modpar.ppar[ikey] = kwargs[ikey]

            if isinstance(kwargs[ikey], float):
                modpar.setPar([ikey, ("%.7e" % kwargs[ikey]), '', ''])
            elif isinstance(kwargs[ikey], int):
                modpar.setPar([ikey, ("%d" % kwargs[ikey]), '', ''])
            elif isinstance(kwargs[ikey], str):
                modpar.setPar([ikey, kwargs[ikey], '', ''])
            elif isinstance(kwargs[ikey], list):
                dum = '['
                for i in range(len(kwargs[ikey])):
                    if type(kwargs[ikey][i]) is float:
                        dum = dum + ("%.7e" % kwargs[ikey][i])
                    elif type(kwargs[ikey][i]) is int:
                        dum = dum + ("%d" % kwargs[ikey][i])
                    elif type(kwargs[ikey][i]) is str:
                        dum = dum + (kwargs[ikey][i])
                    else:
                        raise TypeError('Unknown data type in ' + ikey)

                    if i < len(kwargs[ikey]) - 1:
                        dum = dum + ', '
                dum = dum + ']'
                modpar.setPar([ikey, dum, '', ''])

        modpar.writeParfile()
        ppar = modpar.ppar

    # --------------------------------------------------------------------------------------------
    # Try to get the specified model
    # --------------------------------------------------------------------------------------------
    try:
        mdl = importlib.import_module(model)
    except ImportError:
        try:
            mdl = importlib.import_module('radmc3dPy.models.' + model)
        except ImportError:
            print(' ' + model + '.py could not be imported')
            print(' The model files should either be in the current working directory or')
            print(' in the radmc3d python module directory')
            return

    # --------------------------------------------------------------------------------------------
    # If the current working directory is empty (i.e. no dust setup is present) then
    #   we must make a complete setup and dump the spatial and wavelength grids as well
    #   as the parameters in the radmc3d.inp file
    # --------------------------------------------------------------------------------------------
    if fullsetup:
        # ----------------------------------------------------------------------------------------
        # Create the grid
        # ----------------------------------------------------------------------------------------
        #
        # Check if AMR is activated or not
        #
        if ppar['grid_style'] == 1:
            grid = analyze.radmc3dOctree()
            # Pass all parameters from dfpar to ppar
            if dfpar is not None:
                for ikey in dfpar.keys():
                    ppar[ikey] = dfpar[ikey]
            # Spatial grid
            grid.makeSpatialGrid(ppar=ppar, dfunc=dfunc, model=model, **kwargs)
        else:
            grid = analyze.radmc3dGrid()
            # Spatial grid
            grid.makeSpatialGrid(ppar=ppar)

        # Wavelength grid
        grid.makeWavelengthGrid(ppar=ppar)

        # ----------------------------------------------------------------------------------------
        # Create the input radiation field (stars at this point)
        # ----------------------------------------------------------------------------------------
        radSources = analyze.radmc3dRadSources(ppar=ppar, grid=grid)
        radSources.getStarSpectrum(tstar=ppar['tstar'], rstar=ppar['rstar'])

        # Check if the model has functions to set up continuous starlike sources
        if 'incl_accretion' in ppar:
            stellarsrcEnabled = ppar['incl_accretion']
        else:
            stellarsrcEnabled = False

        if dir(mdl).__contains__('getStellarsrcDensity'):
            if callable(getattr(mdl, 'getStellarsrcDensity')):
                if dir(mdl).__contains__('getStellarsrcTemplates'):
                    if not callable(getattr(mdl, 'getStellarsrcTemplates')):
                        stellarsrcEnabled = False
                else:
                    stellarsrcEnabled = False
            else:
                stellarsrcEnabled = False
        else:
            stellarsrcEnabled = False

        if stellarsrcEnabled:
            dum = mdl.getStellarsrcTemplates(grid=grid, ppar=ppar)
            if dum[0, 0] < 0.:
                radSources.csntemplate = dum.shape[0]
                radSources.cstemp = []
                radSources.cstemptype = 1
                radSources.cststar = dum[:, 0]
                radSources.csrstar = dum[:, 1]
                radSources.csmstar = dum[:, 2]
            else:
                radSources.csntemplate = dum.shape[0]
                radSources.cstemp = dum
                radSources.cstemptype = 2
                radSources.cststar = []
                radSources.csrstar = []
                radSources.csmstar = []

            radSources.csdens = mdl.getStellarsrcDensity(grid=grid, ppar=ppar)
        # --------------------------------------------------------------------------------------------
        # Now write out everything
        # --------------------------------------------------------------------------------------------
        # Frequency grid
        grid.writeWavelengthGrid()
        # Spatial grid
        grid.writeSpatialGrid()
        # Input radiation field
        radSources.writeStarsinp(ppar=ppar)
        # Continuous starlike sources
        if stellarsrcEnabled:
            radSources.writeStellarsrcTemplates()
            radSources.writeStellarsrcDensity(binary=binary)

        print('-------------------------------------------------------------')
        print('Luminosities of radiation sources in the model :')
        totlum = radSources.getTotalLuminosities(readInput=True)
        print('As calculated from the input files :')
        print('Stars : ')
        print("  Star #%d + hotspot        : %.6e" % (0, totlum['lnu_star'][0]))
        for istar in range(1, radSources.nstar):
            print("  Star #%d               : %.6e" % (istar, totlum['lnu_star'][istar]))
        print("Continuous starlike source : %.6e" % totlum['lnu_accdisk'])
        print(' ')
        print('-------------------------------------------------------------')

        writeRadmc3dInp(modpar=ppar)
    # --------------------------------------------------------------------------------------------
    # If the current working directory contains already a dust setup then we can use the
    #   already existing grid files
    # --------------------------------------------------------------------------------------------
    else:
        grid = analyze.readGrid()

    # --------------------------------------------------------------------------------------------
    # Create the gas density distribution
    # --------------------------------------------------------------------------------------------
    # Create the data structure
    data = analyze.radmc3dData(grid)
    # Calculate the gas density and velocity
    # NOTE: the density function in the model sub-modules should provide the gas volume density
    #       in g/cm^3 but RADMC-3D needs the number density in 1/cm^3 so we should convert the
    #       output of the get_density() function to number density using ppar['gasspecMolAbun']
    #       which is the abundance of the gas species with respect to hydrogen divided by the
    #       mean molecular weight
    if dir(mdl).__contains__('getGasDensity'):
        if callable(getattr(mdl, 'getGasDensity')):
            data.rhogas = mdl.getGasDensity(grid=grid, ppar=ppar)
    else:
        raise RuntimeError(' ' + model + '.py does not contain a getGasDensity() function, therefore, '
                           + ' numberdens_***.inp cannot be written')

    # --------------------------------------------------------------------------------------------
    # Create the molecular abundance
    # --------------------------------------------------------------------------------------------
    if dir(mdl).__contains__('getGasAbundance'):
        if callable(getattr(mdl, 'getGasAbundance')):
            for imol in range(len(ppar['gasspec_mol_name'])):
                gasabun = mdl.getGasAbundance(grid=grid, ppar=ppar, ispec=ppar['gasspec_mol_name'][imol])
                data.ndens_mol = data.rhogas / (2.4 * nc.mp) * gasabun

                # Write the gas density
                data.octree = (ppar['grid_style'] == 1)
                data.writeGasDens(ispec=ppar['gasspec_mol_name'][imol], binary=binary)

            if abs(ppar['lines_mode']) > 2:
                for icp in range(len(ppar['gasspec_colpart_name'])):
                    gasabun = mdl.getGasAbundance(grid=grid, ppar=ppar, ispec=ppar['gasspec_colpart_name'][icp])
                    data.ndens_mol = data.rhogas / (2.4 * nc.mp) * gasabun
                    # Write the gas density
                    data.writeGasDens(ispec=ppar['gasspec_colpart_name'][icp], binary=binary)

    else:
        raise RuntimeError(' ' + model + '.py does not contain a getGasAbundance() function, therefore, '
                           + ' numberdens_***.inp cannot be written')
    # --------------------------------------------------------------------------------------------
    # Get the gas velocity field
    # --------------------------------------------------------------------------------------------
    if dir(mdl).__contains__('getVelocity'):
        if callable(getattr(mdl, 'getVelocity')):
            data.gasvel = mdl.getVelocity(grid=grid, ppar=ppar)
            # Write the gas velocity
            data.octree = (ppar['grid_style'] == 1)
            data.writeGasVel(binary=binary)
    else:
        raise RuntimeError(' ' + model + '.py does not contain a getVelocity() function, therefore, '
                           + ' gas_velocity.inp cannot be written')
    # --------------------------------------------------------------------------------------------
    # Get the kinetik gas temperature
    # --------------------------------------------------------------------------------------------
    # Write the gas temperature if specified
    if writeGasTemp:
        if dir(mdl).__contains__('getGasTemperature'):
            if callable(getattr(mdl, 'getGasTemperature')):
                data.gastemp = mdl.getGasTemperature(grid=grid, ppar=ppar)
                # Write the gas temperature
                data.octree = (ppar['grid_style'] == 1)
                data.writeGasTemp(binary=binary)
        else:
            raise RuntimeError(' ' + model + '.py does not contain a getGasTemperature() function, therefore, '
                               + ' gas_temperature.inp cannot be written')

    # --------------------------------------------------------------------------------------------
    # Get the turbulent velocity field
    # --------------------------------------------------------------------------------------------
    if dir(mdl).__contains__('getVTurb'):
        if callable(getattr(mdl, 'getVTurb')):
            data.vturb = mdl.getVTurb(grid=grid, ppar=ppar)
            # Write the turbulent velocity field
            data.octree = (ppar['grid_style'] == 1)
            data.writeVTurb(binary=binary)
    else:
        print(' ' + model + '.py does not contain a getVTurb() function, therefore, zero microturbulent velocity'
              + ' will be assumed everywhere in the model.')
        data.vturb = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
        data.vturb[:, :, :] = 0.
        data.writeVTurb(binary=binary)
    # --------------------------------------------------------------------------------------------
    # Write the line RT control file
    # --------------------------------------------------------------------------------------------
    # Write the lines.inp the main control file for the line RT
    writeLinesInp(ppar=ppar)


def writeRadmcInp(modpar=None, nphot=None):
    """Writes the radmc.inp master command file for the 2D version of radmc

    Parameters
    ----------
    modpar : radmc3dPar
             Contains all parameters of a radmc run.

    nphot  : int
             Number of photons used for the MC simulation
    """

    if nphot is None:
        ppar = modpar.ppar
        nphot = ppar['nphot']

    fname = 'radmc.inp'
    with open(fname, 'w') as wfile:
        wfile.write("nphot       =    %d\n" % nphot)
        wfile.write("iseed       =    -17933201\n")
        wfile.write("imethod     =    2\n")
        wfile.write("ifast       =    0\n")
        wfile.write("enthres     =      0.010000000\n")
        wfile.write("cntdump     =    100000000\n")
        wfile.write("irestart    =        0\n")
        wfile.write("itempdecoup =        1\n")
        wfile.write("iquantum    =        0\n")
        wfile.write("istarsurf   =        1\n")
        wfile.write("nphotdiff   =       15\n")
        wfile.write("errtol      =    1.0000000e-10\n")
        wfile.write("nvstr       =        0\n")
        wfile.write("vserrtol    =     0.0100000\n")
        wfile.write("ivstrt      =        1\n")
        wfile.write("ntemp       =     3000\n")
        wfile.write("temp0       =      0.010000000\n")
        wfile.write("temp1       =        1000000.0\n")
        wfile.close()


def writeRadmc3dInp(modpar=None):
    """Writes the radmc3d.inp master command file for RADMC-3D

    Parameters
    ----------
    modpar   : radmc3dPar
               Contains all parameters of a RADMC-3D run.

    """
    # Select those keywords, whose block name is 'Code parameters'
    ppar = {}

    for ikey in modpar.pblock.keys():
        if modpar.pblock[ikey] == 'Code parameters':
            ppar[ikey] = modpar.ppar[ikey]

    print('Writing radmc3d.inp')

    with open('radmc3d.inp', 'w') as wfile:
        keys = sorted(ppar.keys())
        for key in keys:
            wfile.write('%s %d\n' % (key + ' =', ppar[key]))


def writeLinesInp(ppar=None):
    """Writes the lines.inp master command file for line simulation in RADMC-3D

    Parameters
    ----------

    ppar   : dictionary,
            Contains all parameters of a RADMC-3D run 
    """

    # Do a consistency check
    n1 = len(ppar['gasspec_mol_name'])
    n2 = len(ppar['gasspec_mol_abun'])
    n3 = len(ppar['gasspec_mol_dbase_type'])

    if (n1 != n2) | (n2 != n3):
        raise ValueError('gasspec_mol_name, gasspec_mol_abun and gasspec_mol_dbase_type have '
                         + 'different number of elements')

    if ('gasspec_colpart_name' in ppar) & ('gasspec_colpart_abun' in ppar):
        n4 = len(ppar['gasspec_colpart_name'])
        n5 = len(ppar['gasspec_colpart_abun'])
    else:
        n4 = 0
        n5 = 0

    if n4 != n5:
        raise ValueError('gasspec_colpart_name and gasspec_colpart_abun have different number of elements')

    print('Writing lines.inp')
    with open('lines.inp', 'w') as wfile:
        # File format
        wfile.write("%d\n" % ppar['lines_mode'])
        # Nr of gas species
        wfile.write("%d\n" % n1)
        # Gas species name and database type

        if np.abs(ppar['lines_mode']) <= 2:
            for imol in range(n1):
                wfile.write("%s %s %d %d %d\n" % (ppar['gasspec_mol_name'][imol],
                                                  ppar['gasspec_mol_dbase_type'][imol], 0, 0, 0))
        else:
            for imol in range(n1):
                wfile.write("%s %s %d %d %d\n" % (ppar['gasspec_mol_name'][imol],
                                                  ppar['gasspec_mol_dbase_type'][imol], 0, 0, n4))

            if n4 > 0:
                for icp in range(n4):
                    wfile.write("%s\n" % ppar['gasspec_colpart_name'][icp])
            else:
                raise RuntimeError(' An NLTE line excitation method is selected (lines_mode='
                                   + ("%d" % ppar["lines_mode"]) + '), but no collisional partner is '
                                   + ' given in the parameter file. ')


def validateModel(model='', dustModel=False, gasModel=False, writeDustTemp=False, octree=False):
    """
    Function to validate a model. It checks three things: 1) whether or not the model can be imported,
    2) whether the model has all the function to be used as dust and/or gas model, 3) if it has the right
    number of arguments. The function names tested are getDefaultParams, getDustDensity, getGasDensity, 
    getGasAbundance, getVTurb, getVelocity, getDustTempearture (optional).

    Parameters
    ----------

    model       : str
                  Name of the model to be tested

    dustModel   : bool
                  If True the existence of functions getDustDensity() and getDustTemperature() will be checked.
                  The latter is only checked if writeDustTemp is set to True.

    gasModel    : bool
                  If True the existence of functions getGasDensity(), getGasAbundance(), getVTurb(), getVelocity()
                  will be checked.

    writeDustTemp: bool
                   If True the existence of the function getDustTemperature() will be checked.

    octree      : bool
                  If True the number of argument of the model functions will be checked. For regular grids only two 
                  arguments should be present for the grid instance and for the parameter dictionary (grid, ppar). 
                  For a model to be used with octree AMR three additional arguments for the three spatial coordiantes
                  (x,y,z) should be present. The argument sequence should then be x, y, z, grid, ppar.

    Returns
    -------
    A boolean True if the model is valid and False if it is not. 

    """

    #
    # First check if the model can be imported
    #
    try:
        mdl = importlib.import_module(model)
    except ImportError:
        try:
            mdl = importlib.import_module('radmc3dPy.models.' + model)
        except ImportError:
            print(' ' + model + '.py could not be imported'
                  + ' The model files should either be in the current working directory or'
                  + ' in the radmc3d python module directory')
            return False

    isValid = True
    #
    # Now check the function names in the model
    #
    fnamelist = [f[0] for f in inspect.getmembers(mdl) if inspect.isfunction(f[1])]

    if 'getDefaultParams' not in fnamelist:
        print(model + ' does not contain a function to provide default parameters (getDefaultParams)')
        isValid = False

    if dustModel:
        if 'getDustDensity' not in fnamelist:
            print(model + ' does not contain a function to provide dust density (getDustDensity)')
            isValid = False

        if writeDustTemp:
            if 'getDustTemperature' not in fnamelist:
                print(model + ' does not contain a function to provide dust temperature (getDustTemperature)'
                      + 'yet the setup function has been called with the option to write the dust temperature.')
                isValid = False

    if gasModel:
        if 'getGasDensity' not in fnamelist:
            print(model + ' does not contain a function to provide gas density (getGasDensity)')
            isValid = False

        if 'getGasAbundance' not in fnamelist:
            print(model + ' does not contain a function to provide molecular abundance (getGasAbundance)')
            isValid = False

        if 'getVTurb' not in fnamelist:
            print(model + ' does not contain a function to provide turbulent velocity (getVTurb)')
            isValid = False

        if 'getVelocity' not in fnamelist:
            print(model + ' does not contain a function to provide gas velocity (getVelocity)')
            isValid = False

    #
    # Check the number of arguments
    #
    if dustModel:
        if sys.version_info.major == 2:
            arglist = inspect.getargspec(mdl.getDustDensity).args
        else:
            arglist = inspect.signature(mdl.getDustDensity).parameters.keys()

        argnames = ''
        if len(arglist) > 0:
            argnames = arglist[0]
            for iarg in arglist[1:]:
                argnames += ', ' + iarg
        if octree:
            if len(arglist) < 5:
                print(model + '.getDustDensity() has only ' + ("%d" % len(arglist)) + ' arguments : ', argnames,
                      'To use octree the argument list should be : x=None, y=None, z=None, grid=None, ppar=None)')
                isValid = False
        else:
            if len(arglist) < 2:
                print(model + '.getDustDensity() has only ' + ("%d" % len(arglist)) + ' arguments : ', argnames,
                      ' The minimal argument list of a model function should be : (grid=None, ppar=None)')
                isValid = False

        if writeDustTemp:
            if sys.version_info.major == 2:
                arglist = inspect.getargspec(mdl.getDustTemperature).args
            else:
                arglist = inspect.signature(mdl.getDustTemperature).parameters.keys()

            argnames = ''
            if len(arglist) > 0:
                argnames = arglist[0]
                for iarg in arglist[1:]:
                    argnames += ', ' + iarg
            if octree:
                if len(arglist) < 5:
                    print(model + '.getDustTemperature() has only ' + ("%d" % len(arglist)) + ' arguments : ', argnames,
                          ' To use octree the argument list should be x=None, y=None, z=None, grid=None, ppar=None)')
                    isValid = False
            else:
                if len(arglist) < 2:
                    print(model + '.getDustTemperature() has only ' + ("%d" % len(arglist)) + ' arguments : ', argnames,
                          ' The minimal argument list of a model function should be : (grid=None, ppar=None)')
                    isValid = False

    if gasModel:
        if sys.version_info.major == 2:
            arglist = inspect.getargspec(mdl.getGasDensity).args
        else:
            arglist = inspect.signature(mdl.getGasDensity).parameters.keys()
        argnames = ''
        if len(arglist) > 0:
            argnames = arglist[0]
            for iarg in arglist[1:]:
                argnames += ', ' + iarg
        if octree:
            if len(arglist) < 5:
                print(model + '.getGasDensity() has only ' + ("%d" % len(arglist)) + ' arguments : ', argnames,
                      ' To use octree the argument list should be : x=None, y=None, z=None, grid=None, ppar=None)')
                isValid = False
        else:
            if len(arglist) < 2:
                print(model + '.getGasDensity() has only ' + ("%d" % len(arglist)) + ' arguments : ', argnames,
                      ' The minimal argument list of a model function should be : (grid=None, ppar=None)')
                isValid = False

        if sys.version_info.major == 2:
            arglist = inspect.getargspec(mdl.getGasAbundance).args
        else:
            arglist = inspect.signature(mdl.getGasAbundance).parameters.keys()
        argnames = ''
        if len(arglist) > 0:
            argnames = arglist[0]
            for iarg in arglist[1:]:
                argnames += ', ' + iarg
        if octree:
            if len(arglist) < 5:
                print(model + '.getGasAbundance() has only ' + ("%d" % len(arglist)) + ' arguments : ', argnames,
                      ' To use octree the argument list should be : x=None, y=None, z=None, grid=None, ppar=None)')
                isValid = False
        else:
            if len(arglist) < 2:
                print(model + '.getGasAbundance() has only ' + ("%d" % len(arglist)) + ' arguments : ', argnames,
                      ' The minimal argument list of a model function should be : (grid=None, ppar=None)')
                isValid = False

        if sys.version_info.major == 2:
            arglist = inspect.getargspec(mdl.getVTurb).args
        else:
            arglist = inspect.signature(mdl.getVTurb).parameters.keys()
        argnames = ''
        if len(arglist) > 0:
            argnames = arglist[0]
            for iarg in arglist[1:]:
                argnames += ', ' + iarg
        if octree:
            if len(arglist) < 5:
                print(model + '.getVTurb() has only ' + ("%d" % len(arglist)) + ' arguments : ', argnames,
                      ' To use octree the argument list should be : x=None, y=None, z=None, grid=None, ppar=None)')
                isValid = False
        else:
            if len(arglist) < 2:
                print(model + '.getVTurb() has only ' + ("%d" % len(arglist)) + ' arguments : ', argnames,
                      'The minimal argument list of a model function should be : (grid=None, ppar=None)')
                isValid = False

        if sys.version_info.major == 2:
            arglist = inspect.getargspec(mdl.getVelocity).args
        else:
            arglist = inspect.signature(mdl.getVelocity).parameters.keys()
        argnames = ''
        if len(arglist) > 0:
            argnames = arglist[0]
            for iarg in arglist[1:]:
                argnames += ', ' + iarg
        if octree:
            if len(arglist) < 5:
                print(model + '.getVelocity() has only ' + ("%d" % len(arglist)) + ' arguments : ', argnames,
                      'To use octree the argument list should be : x=None, y=None, z=None, grid=None, ppar=None)')
                isValid = False
        else:
            if len(arglist) < 2:
                print(model + '.getVelocity() has only ' + ("%d" % len(arglist)) + ' arguments : ', argnames,
                      ' The minimal argument list of a model function should be : (grid=None, ppar=None)')
                isValid = False

    #
    # If it passed all tests so far then formally the model should be OK. There is no guarantee, though
    # that it will work properly.
    #
    return isValid


