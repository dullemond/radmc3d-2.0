"""A 2D envelope with density dependence on the theta coordinate in spherical grid

Original IDL model by Kees Dullemond, Python translation by Attila Juhasz
"""
from __future__ import absolute_import
from __future__ import print_function
import traceback
try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use the python module of RADMC-3D you need to install Numpy')
    print(traceback.format_exc())

from .. natconst import *


def getModelDesc():
    """
    Provides a brief description of the model
    """

    return "Example model: A 2D envelope with density dependence on the theta coordinate in spherical grid"
           

def getDefaultParams():
    """
    Provides default parameter values 

    Returns a list whose elements are also lists with three elements:
    1) parameter name, 2) parameter value, 3) parameter description
    All three elements should be strings. The string of the parameter
    value will be directly written out to the parameter file if requested,
    and the value of the string expression will be evaluated and be put
    to radmc3dData.ppar. The third element contains the description of the
    parameter which will be written in the comment field of the line when
    a parameter file is written. 
    """

    defpar = [['mstar', '1.0*ms', 'Mass of the star(s)'],
              ['pstar', '[0., 0., 0.]', 'Position of the star(s) (cartesian coordinates)'],
              ['rstar', '1.0*rs', 'Radius of the star(s)'],
              ['tstar', '1.0*ts', 'Effective temperature of the star(s)'],
              ['crd_sys', "'sph'", 'Coordinate system used (car/sph)'],
              ['nx', '100', 'Number of grid points in the first dimension'],
              ['ny', '60', 'Number of grid points in the second dimension'],
              ['nz', '0', 'Number of grid points in the third dimension'],
              ['xbound', '[5.0*au, 100.0*au]', 'Boundaries for the x-grid'],
              ['ybound', '[0.0, pi/2.0]', 'Boundaries for the y-grid'],
              ['zbound', '[0.0, 2.0*pi]', 'Boundaries for the z-grid'],
              ['nw', '[20,100,30]', 'Number of points in the wavelength grid'],
              ['wbound', '[0.1, 7., 25., 1e4]', 'Boundaries for the wavelength grid'],
              ['dustkappa_ext', "['silicate']", 'Dust opacity file name extension'],
              ['nphot', '100000', 'Number of photons in the thermal Monte Carlo simulation'],
              ['scattering_mode_max', '0', '0 - no scattering, 1 - isotropic scattering, 2 - anizotropic scattering'],
              ['prho', '-2.0', ' Power exponent of the radial density distribution'],
              ['hpr', '0.1', ' Vertical scale height'],
              ['rho0', '1e-16*10000', 'Central density']]

    return defpar


def getDustDensity(grid=None, ppar=None):
    """Calculates the dust density distribution 
    
    Parameters
    ----------
    grid : radmc3dGrid
            An instance of the radmc3dGrid class containing the spatial and wavelength grid
    
    ppar : dictionary
            Dictionary containing all parameters of the model 
    
    Returns
    -------
    Returns the dust volume density in g/cm^3
    """

    xx, yy = np.meshgrid(grid.x, np.pi/2.-grid.y)
    xx = xx.swapaxes(0, 1)
    yy = yy.swapaxes(0, 1)

    rho = np.zeros([grid.nx, grid.ny, grid.nz, 1], dtype=np.float64)

    for iz in range(grid.nz):
        rho[:, :, iz, 0] = ppar['rho0'] * (xx/au)**ppar['prho'] * np.exp(-(yy/ppar['hpr'])**2)

    return rho
