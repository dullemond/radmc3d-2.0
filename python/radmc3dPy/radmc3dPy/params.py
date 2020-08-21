"""This module contains classes for handling of parameters of the model setup and for file I/O
"""

from __future__ import absolute_import
from __future__ import print_function
import traceback
import importlib

try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use the python module of RADMC-3D you need to install Numpy')
    print(traceback.format_exc())


from . natconst import *

class radmc3dPar(object):
    """Parameter class for a RADMC-3D model.


    Attributes
    ----------

    ppar   : dictionary
            Contains parameter values with parameter names as keys

    pdesc  : dictionary
            Contains parameter description (comments in the parameter file) with parameter names as keys

    pblock : dictionary
            Contains the block names in the parameter file and parameter names as values

    pvalstr: dictionary
            Contains parameter values as strings with parameter names as keys

    """

    def __init__(self):

        self.ppar = {}
        self.pdesc = {}
        self.pblock = {}
        self.pvalstr = {}

    def readPar(self, fname=''):
        """Reads a parameter file.
        The parameters in the files should follow python syntax


        Parameters
        ----------

        fname  : str, optional
                File name to be read (if omitted problem_params.inp is used)

        Returns
        -------
            Returns a dictionary with the parameter names as keys

        """

        if fname == '':
            fname = 'problem_params.inp'

        with open(fname, 'r') as rfile:

            cchar = '#'
            lbchar = ""

            # Add pi to the local variables
            pi = np.pi
            # --------------------------------------------------------------------------------------------------------
            # First read every line that is not commented (i.e. does not begin with a comment character)
            # --------------------------------------------------------------------------------------------------------
            dumlist = []
            dumline = '-'

            dumline = rfile.readline()
            while dumline != '':
                # First check if the line is commented out, in which case ignore it
                comment = False
                if dumline[0] == cchar:
                    if dumline.find('Block') < 0:
                        comment = True

                # OK, the line is not commented out, now check if it contains a '=' sign (if not ignore the line)
                if not comment:
                    # Check if we have an empty line in which case also ignore it
                    if dumline.strip() != '':
                        dumlist.append(dumline)

                # Read the next line
                dumline = rfile.readline()

        # -------------------------------------------------------------------------------------------------------------
        # After every line in the file was read try to decode the lines to
        #  [variable name] = [variable value] # [comment]
        # also try to catch if an expression has been broken into multiple lines
        # -------------------------------------------------------------------------------------------------------------

        varlist = []
        iline = 0
        while iline < len(dumlist):
            # First check if the line contains an '=' sign if not we have a problem
            #  expression broken into multiple lines are should already be dealt with
            ind = dumlist[iline].find('=')
            if ind <= 0:
                if dumlist[iline].find('Block') <= 0:
                    raise SyntaxError(' Invalid expression in line ' + ("%d" % iline))
                else:
                    if dumlist[iline].find(':') <= 0:
                        raise SyntaxError('Invalid block identified. \n'
                                          + ' The syntax of the block name field is : \n'
                                          + ' # Block : BlockName ')
                    else:
                        blockname = dumlist[iline].split(':')[1].strip()

            else:
                # The line contains a '=' sign and a variable name, so let's check if the
                #  value expression is broken into multiple lines
                vlist = dumlist[iline].split('=')
                lbind = vlist[1].find('\\')
                cind = vlist[1].find('#')

                # The line is full not broken
                if lbind == -1:
                    # Check if there is a comment field
                    if cind >= 0:
                        vlist = [vlist[0], vlist[1][:cind], vlist[1][cind + 1:], blockname]
                    else:
                        vlist = [vlist[0], vlist[1][:cind], ' ', blockname]

                    varlist.append(vlist)
                # The value expression is broken into multiple lines; take all lines and join the pieces
                else:
                    # Check if there is any comment in the line
                    inBrokenLine = False
                    if cind >= 0:
                        # Part of the line is commented, now check if the line break is before or after the
                        # comment character
                        if lbind > cind:
                            # The line break is in the comment field so there is no real line break
                            vlist = [vlist[0], vlist[1][:cind], vlist[1][cind + 1:], blockname]
                        else:
                            # The line break is before the comment character
                            inBrokenLine = True
                            expr = vlist[1][:lbind]
                            com = vlist[1][cind + 1:]
                    else:
                        inBrokenLine = True
                        expr = vlist[1][:lbind]
                        com = ' '

                    if inBrokenLine:
                        # Now gather all other pieces of this line
                        iline2 = 0
                        while inBrokenLine:
                            iline2 = iline2 + 1
                            dummy = dumlist[iline + iline2]
                            # Search for comments
                            cind2 = dummy.find('#')
                            # Search for another line break
                            lbind2 = dummy.find('\\')

                            # TODO:
                            # At the moment I neglect the possiblity that the second line in a broken long line begins
                            # with a linebreak or commented out

                            # There is comment
                            if cind2 > 0:

                                # There is line break
                                if lbind2 > 0:
                                    # The line break is commented out
                                    if lbind2 > cind:
                                        expr = expr + dummy[:cind2].strip()
                                        com = com + dummy[cind2 + 1:]
                                        inBrokenLine = False
                                    else:
                                        # The line break is not commented out
                                        expr = expr + dummy[:lbind2].strip()
                                        com = com + dummy[cind2 + 1:]
                                else:
                                    # There is no line break
                                    expr = expr + dummy[:cind2].strip()
                                    com = com + dummy[cind2 + 1:]
                                    inBrokenLine = False

                            # There is no comment
                            else:
                                # There is a line break
                                if lbind2 > 0:
                                    expr = expr + dummy[:lbind2].strip()
                                    com = com + dummy[cind2 + 1:]

                                # There is no line break
                                else:
                                    expr = expr + dummy[:cind2].strip()
                                    com = com + ' '
                                    inBrokenLine = False
                        iline = iline + iline2
                        vlist = [vlist[0], expr, com, blockname]
                        varlist.append(vlist)

            iline = iline + 1
        # -------------------------------------------------------------------------------------------------------
        # Now evaluate the expressions in the value field and make the final dictionary
        # -------------------------------------------------------------------------------------------------------
        self.ppar = {}
        glob = globals()
        glob['pi'] = np.pi
        loc = locals()
        loc['pi'] = np.pi
        for i in range(len(varlist)):
            try:
                val = eval(varlist[i][1], glob)
                glob[varlist[i][0].strip()] = val
            except:
                try:
                    val = eval(varlist[i][1], loc)
                    loc[varlist[i][0].strip()] = val
                except:
                    raise SyntaxError('Unknown expression "' + varlist[i][1] + '"')
            self.ppar[varlist[i][0].strip()] = val
            self.pvalstr[varlist[i][0].strip()] = varlist[i][1].strip()
            self.pdesc[varlist[i][0].strip()] = varlist[i][2].strip()
            self.pblock[varlist[i][0].strip()] = varlist[i][3].strip()
        return

    # --------------------------------------------------------------------------------------------------
    def setPar(self, parlist=None):
        """Sets a parameter in the radmc3DPar class.
        If the paramter is already defined its value will be modified

        Parameters
        ----------

        parlist : list
                  If the parameter is already defined parlist should be a two element
                  list 1) parameter name, 2) parameter expression/value as a string

                  If the parameter is not yet defined parlist should be a four element
                  list 1) parameter name, 2) parameter expression/value as a string
                  3) Parameter description (= comment field in the parameter file)
        """

        if parlist is None:
            raise ValueError('Unknown parlist. \n No parameters to set.')
        else:
            parname = parlist[0].strip()

        #
        # Check whether or not the parameter is already defined
        #
        new_par = False
        if len(parlist) == 2:
            if parname not in self.ppar:
                raise ValueError(' parlist has too few elements.\n'
                                 + ' The argument of radmc3dPar.setPar() should be a four element list if a new \n'
                                 + ' parameter is defined 1) parameter name, 2) parameter expression/value as a '
                                 + ' string\n'
                                 + ' 3) Parameter description (= comment field in the parameter file) \n'
                                 + ' 4) The name of the block in which the parameter must be placed in the \n'
                                 + ' problem_params.inp file')
        else:
            new_par = True

        # Add pi to the local variables
        pi = np.pi
        #
        # Add the parameter to the dictionaries /change its value
        #
        glob = globals()
        glob['pi'] = np.pi
        loc = locals()
        loc['pi'] = np.pi

        try:
            self.ppar[parname] = eval(parlist[1].strip(), glob)
            glob[parname] = self.ppar[parname]
        except:
            try:
                self.ppar[parname] = eval(parlist[1].strip(), loc)
                loc[parname] = self.ppar[parname]
            except:
                raise SyntaxError('Unknown expression ' + parlist[1].strip())

        self.pvalstr[parname] = parlist[1].strip()

        if new_par:
            if parname not in self.pdesc:
                self.pdesc[parname] = parlist[2].strip()
            if len(parlist) == 4:
                if parname not in self.pblock:
                    self.pblock[parname] = parlist[3].strip()

    def loadDefaults(self, model='', ppar=None, reset=True):
        """Sets all parameters to default values.

        Parameters
        ----------
        model : str
                Model name whose paraemters should also be loaded

        ppar  : dictionary
                Contains parameter values as string and parameter names as keys
                Default values will be re-set to the values in this dictionary

        reset : bool
                If True the all class attributes will be re-initialized before
                the default values would be loaded. I.e. it will remove all entries
                from the dictionary that does not conain default values either in this
                function or in the optional ppar keyword argument
        """

        if reset:
            self.ppar = {}
            self.pvalstr = {}
            self.pdesc = {}
            self.pblock = {}

        #
        # Radiation sources
        #
        self.setPar(['incl_disc_stellarsrc', 'True', '# Switches on (True) or off (False) discrete stellar sources)',
                     'Radiation sources'])
        self.setPar(['mstar', '[1.0*ms]', '# Mass of the star(s)', 'Radiation sources'])
        self.setPar(['rstar', '[2.0*rs]', '# Radius of the star(s)', 'Radiation sources'])
        self.setPar(['tstar', '[4000.0]', '# Effective temperature of the star(s) [K]', 'Radiation sources'])
        self.setPar(
            ['pstar', '[0.0, 0.0, 0.0]', '# Position of the star(s) (cartesian coordinates)', 'Radiation sources'])
        self.setPar(['staremis_type', '["blackbody"]', '# Stellar emission type ("blackbody", "kurucz", "nextgen")',
                     'Radiation sources'])
        self.setPar(
            ['incl_cont_stellarsrc', 'False', '# Switches on (True) or off (False) continuous stellar sources )',
             'Radiation sources'])
        #
        # Grid parameters
        #
        self.setPar(['crd_sys', "'sph'", '  Coordinate system used (car/sph)', 'Grid parameters'])
        self.setPar(
            ['grid_style', '0', '  0 - Regular grid, 1 - Octree AMR, 10 - Layered/nested grid (not yet supported)',
             'Grid parameters'])
        self.setPar(
            ['nx', '50', '  Number of grid points in the first dimension (to switch off this dimension set it to 0)',
             'Grid parameters'])
        self.setPar(
            ['ny', '30', '  Number of grid points in the second dimension (to switch off this dimension set it to 0)',
             'Grid parameters'])
        self.setPar(
            ['nz', '36', '  Number of grid points in the third dimension (to switch off this dimension set it to 0)',
             'Grid parameters'])
        self.setPar(['xbound', '[1.0*au, 100.*au]', '  Boundaries for the x grid', 'Grid parameters'])
        self.setPar(['ybound', '[0.0, pi]', '  Boundaries for the y grid', 'Grid parameters'])
        self.setPar(['zbound', '[0.0, 2.0*pi]', '  Boundraries for the z grid', 'Grid parameters'])
        self.setPar(['xres_nlev', '3', 'Number of refinement levels (spherical coordinates only', 'Grid parameters'])
        self.setPar(['xres_nspan', '3', 'Number of the original grid cells to refine (spherical coordinates only)',
                     'Grid parameters'])
        self.setPar(
            ['xres_nstep', '3', 'Number of grid cells to create in a refinement level (spherical coordinates only)',
             'Grid parameters'])
        self.setPar(['wbound', '[0.1, 7.0, 25., 1e4]', '  Boundraries for the wavelength grid', 'Grid parameters'])
        self.setPar(['nw', '[19, 50, 30]', '  Number of points in the wavelength grid', 'Grid parameters'])
        self.setPar(['levelMaxLimit', '5', '  Highest refinement level in octree AMR', 'Grid parameters'])

        #
        # Dust opacity
        #
        self.setPar(['lnk_fname',
                     "['/disk2/juhasz/Data/JPDOC/astrosil/astrosil_WD2001_new.lnk']",
                     ' ', 'Dust opacity'])
        self.setPar(['gdens', '[3.6, 1.8]', ' Bulk density of the materials in g/cm^3', 'Dust opacity'])
        self.setPar(['gsmin', '0.1', ' Minimum grain size', 'Dust opacity'])
        self.setPar(['gsmax', '10.0', ' Maximum grain size', 'Dust opacity'])
        self.setPar(['ngs', '1', ' Number of grain sizes', 'Dust opacity'])
        self.setPar(['gsdist_powex', '-3.5', ' Grain size distribution power exponent', 'Dust opacity'])

        self.setPar(['nscatang', '180', ' Number of scattering angles (only for scattering_mode_max=5)',
                     'Dust opacity'])
        self.setPar(['logawidth', '0', ' If >0 the opacity will be averaged over a small sample around the specified '
                                       'grain size, with logawidth being the variance of the Gaussian distribution. ',
                     'Dust opacity'])
        self.setPar(['wfact', '3.0', ' Grid width of na sampling points in units of logawidth.', 'Dust opacity'])
        self.setPar(['na', '20', ' Number of size sampling points (if logawidth set, default=20', 'Dust opacity'])
        self.setPar(['chopforwardt', '0.0', ' If >0 this gives the angle (in degrees from forward) within which the '
                                            'scattering phase function should be kept constant', 'Dust opacity'])
        self.setPar(['errtol', '0.01', ' Tolerance of the relative difference between kscat and the integral over the '
                                       'zscat Z11 element over angle.', 'Dust opacity'])
        self.setPar(['verbose', 'False', ' If set to True, the code will give some feedback so that one knows what it '
                                         'is doing if it becomes slow.', 'Dust opacity'])
        self.setPar(['extrapolate', 'True', ' If True extrapolates for wavelengths outside the ones covered by the '
                                             'optical constants', 'Dust opacity'])

        self.setPar(['mixabun', '[0.75, 0.25]', ' Mass fractions of the dust componetns to be mixed', 'Dust opacity'])
        self.setPar(['dustkappa_ext', "['silicate']", ' ', 'Dust opacity'])

        #
        # Gas line RT
        #
        self.setPar(
            ['gasspec_mol_name', "['co']", '  Name of the gas species - the extension of the molecule_EXT.inp file',
             'Gas line RT'])
        self.setPar(['gasspec_mol_abun', '[1e-4]', '  Abundance of the molecule', 'Gas line RT'])
        self.setPar(['gasspec_mol_dbase_type', "['leiden']", '  leiden or linelist', 'Gas line RT'])
        self.setPar(
            ['gasspec_colpart_name', "['h2']", '  Name of the gas species - the extension of the molecule_EXT.inp file',
             'Gas line RT'])
        self.setPar(['gasspec_colpart_abun', '[1e0]', '  Abundance of the molecule', 'Gas line RT'])
        self.setPar(['gasspec_vturb', '0.1e5', '  Microturbulence', 'Gas line RT'])
        #
        # Code parameters
        #
        self.setPar(['nphot', 'int(1e5)', '  Nr of photons for the thermal Monte Carlo', 'Code parameters'])
        self.setPar(['nphot_scat', 'int(3e5)', '  Nr of photons for the scattering Monte Carlo (for images)',
                     'Code parameters'])
        self.setPar(['nphot_spec', 'int(1e5)', '  Nr of photons for the scattering Monte Carlo (for spectra)',
                     'Code parameters'])
        self.setPar(
            ['scattering_mode_max', '1', '  0 - no scattering, 1 - isotropic scattering, 2 - anizotropic scattering',
             'Code parameters'])
        self.setPar(['lines_mode', '-1', '  Line raytracing mode', 'Code parameters'])
        self.setPar(['istar_sphere', '0',
                     '  1 - take into account the finite size of the star, 0 - take the star to be point-like',
                     'Code parameters'])
        self.setPar(['itempdecoup', '1', '  Enable for different dust components to have different temperatures',
                     'Code parameters'])
        self.setPar(['tgas_eq_tdust', '1', '  Take the dust temperature to identical to the gas temperature',
                     'Code parameters'])
        self.setPar(
            ['rto_style', '1', '  Format of outpuf files (1-ascii, 2-unformatted f77, 3-binary', 'Code parameters'])
        self.setPar(
            ['modified_random_walk', '0', '  Switched on (1) and off (0) modified random walk', 'Code parameters'])
        #
        # Model parameters
        #
        if model != '':
            try:
                mdl = importlib.import_module(model)
            except:
                try:
                    mdl = importlib.import_module('radmc3dPy.models.' + model)
                except:
                    print(model + '.py could not be imported. \n '
                          + 'The model files should either be in the current working directory \n '
                          + 'or in the radmc3d python module directory')
                    print(traceback.format_exc())

            modpar = mdl.getDefaultParams()
            for i in range(len(modpar)):
                dum = modpar[i]
                if len(dum) == 3:
                    dum.append('Model ' + model)
                self.setPar(dum)

    def printPar(self):
        """Prints the parameters of the current model.

        """

        #
        # First get the unique block names
        #

        blocknames = ['Radiation sources', 'Grid parameters', 'Dust opacity', 'Gas line RT', 'Code parameters']
        for key in self.pblock.keys():
            dum = self.pblock[key]
            if not blocknames.__contains__(dum):
                blocknames.append(dum)

        #
        # Get the parameter block names and distionary keys
        #
        par_keys = list(self.pblock.keys())
        par_block = list(self.pblock.values())

        #
        # Print the parameters by blocks
        #
        for iblock in blocknames:
            print(
                '%s' % '# ------------------------------------------------------------------------------------------'
                       + '-------------------------------')
            txt = '# Block: ' + iblock
            print('%s' % txt)
            print(
                '%s' % '# ------------------------------------------------------------------------------------------'
                       + '-------------------------------')

            keys = []
            for i in range(len(par_block)):
                if par_block[i] == iblock:
                    keys.append(par_keys[i])

            keys.sort()
            for key in keys:
                print(key.ljust(25) + ' = ' + self.pvalstr[key].strip() + '  # ' + self.pdesc[key].strip())
                # --------------------------------------------------------------------------------------------------

    def writeParfile(self, fname=''):
        """Writes a parameter file.

        Parameters
        ----------

        fname  : str, optional
                File name to be read (if omitted problem_params.inp is used)

        """

        if fname == '':
            fname = 'problem_params.inp'

        print('Writing ' + fname)

        #
        # First get the uniq block names
        #

        blocknames = ['Radiation sources', 'Grid parameters', 'Dust opacity', 'Gas line RT', 'Code parameters']
        for key in self.pblock:
            dum = self.pblock[key]
            if not blocknames.__contains__(dum):
                blocknames.append(dum)

        with open(fname, 'w') as wfile:
            #
            # Write header
            #

            wfile.write(
                '%s\n' % ('#########################################################################################'
                         + '##################################'))
            wfile.write('%s\n' % '# RADMC-3D PARAMETER SETUP')
            wfile.write('%s\n' % '# Created by the python module of RADMC-3D')
            wfile.write(
                '%s\n' % ('#########################################################################################'
                         + '##################################'))

            #
            # Get the parameter block names and distionary keys
            #
            par_keys = list(self.pblock.keys())
            par_block = list(self.pblock.values())

            #
            # Write the parameterfile
            #
            for iblock in blocknames:
                wfile.write(
                    '%s\n' % ('# -----------------------------------------------------------------------------------'
                             + '--------------------------------------'))
                txt = '# Block: ' + iblock
                wfile.write('%s\n' % txt)
                wfile.write(
                    '%s\n' % ('# -----------------------------------------------------------------------------------'
                             + '--------------------------------------'))

                keys = []
                for i in range(len(par_block)):
                    if par_block[i] == iblock:
                        keys.append(par_keys[i])

                keys.sort()
                for key in keys:
                    wfile.write(key.ljust(25) + ' = ' + self.pvalstr[key].strip() + '  # '
                                + self.pdesc[key].strip() + '\n')
