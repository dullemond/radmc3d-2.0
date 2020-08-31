"""A 1D simple velocity gradient model to calculate lines with the LVG method

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
    """Provides a brief description of the model
    """

    return "Example model: A 1D simple velocity gradient model to calculate lines with the LVG method"
           

def getDefaultParams():
    """Provides default parameter values 

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
              ['crd_sys', "'car'", 'Coordinate system used (car/sph)'],
              ['nx', '10', 'Number of grid points in the first dimension'],
              ['ny', '1', 'Number of grid points in the second dimension'],
              ['nz', '1', 'Number of grid points in the third dimension'],
              ['xbound', '[-1000.0*au, 1000.0*au]', 'Boundaries for the x-grid'],
              ['ybound', '[-1000.0*au/nx, 1000.0*au/nx]', 'Boundaries for the y-grid'],
              ['zbound', '[-1000.0*au/nx, 1000.0*au/nx]', 'Boundaries for the z-grid'],
              ['nw', '[20,100,30]', 'Number of points in the wavelength grid'],
              ['wbound', '[0.1, 7., 25., 1e4]', 'Boundaries for the wavelength grid'],
              ['dustkappa_ext', "['silicate']", 'Dust opacity file name extension'],
              ['nphot', '1000000', 'Number of photons in the thermal Monte Carlo simulation'],
              ['lines_mode', '3', ''],
              ['scattering_mode_max', '1', '0 - no scattering, 1 - isotropic scattering, 2 - anizotropic scattering'],
              ['gasspec_mol_name', "['co']", ''],
              ['gasspec_mol_abun', '[1e-4]', ''],
              ['gasspec_mol_dbase_type', "['leiden']", ''],
              ['gasspec_colpart_name', "['h2']", ''],
              ['gasspec_colpart_abun', '[1e0]', ''],
              ['gasspec_vturb', '1e5', 'Microturbulent linewidth'],
              ['abun_h2', '0.5', ''],
              ['abun_he', '0.1', ''],
              ['nh2', '1e5', ''],
              ['temp0', '30.', ''],
              ['tdust0', '30.', ''],
              ['dusttogas', '1e-2', ''],
              ['dvdau', '1e-2*1e5', '']]

    return defpar


def getGasTemperature(grid=None, ppar=None):
    """Calculates the gas temperature
    
    Parameters
    ----------
    grid : radmc3dGrid
            An instance of the radmc3dGrid class containing the spatial and wavelength grid
    
    ppar : dictionary
            Dictionary containing all parameters of the model 
    
    Returns
    -------
    Returns the gas temperature in K
    """

    tgas = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) + ppar['temp0']
    return tgas


def getDustTemperature(grid=None, ppar=None):
    """Calculates/sets the dust temperature
    
    Parameters
    ----------
    grid : radmc3dGrid
            An instance of the radmc3dGrid class containing the spatial and wavelength grid
    
    ppar : dictionary
            Dictionary containing all parameters of the model 
    
    Returns
    -------
    Returns the dust temperature in K
    
    """

    tdust = np.zeros([grid.nx, grid.ny, grid.nz, 1], dtype=np.float64) + ppar['tdust0']
    return tdust


def getGasAbundance(grid=None, ppar=None, ispec=''):
    """Calculates/sets the molecular abundance of species ispec 
    The number density of a molecule is rhogas * abun 
   
    Parameters
    ----------
    grid  : radmc3dGrid
            An instance of the radmc3dGrid class containing the spatial and wavelength grid

    ppar  : dictionary
            Dictionary containing all parameters of the model 

    ispec : str
            The name of the gas species whose abundance should be calculated

    Returns
    -------
    Returns the abundance as an ndarray
    """
    # Mass of gas per H2-molecule
    # mgas    = mp*(2.0*ppar['abun_h2']+4*ppar['abun_he'])/ppar['abun_h2']

    if ispec in ppar['gasspec_mol_name']:
        gasabun = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) 
        ind = ppar['gasspec_mol_name'].index(ispec)
        gasabun[:, :, :] = ppar['gasspec_mol_abun'][ind]  # /mgas
 
    elif ispec in ppar['gasspec_colpart_name']:
        gasabun = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) 
        ind = ppar['gasspec_colpart_name'].index(ispec)
        gasabun[:, :, :] = ppar['gasspec_colpart_abun'][ind]  # /mgas
    else:
        raise ValueError(' The abundance of "'+ispec+'" is not specified in the parameter file')
   
    return gasabun


def getGasDensity(grid=None, ppar=None):
    """Calculates the total gas density distribution 
    
    Parameters
    ----------
    grid : radmc3dGrid
            An instance of the radmc3dGrid class containing the spatial and wavelength grid
    
    ppar : dictionary
            Dictionary containing all parameters of the model 
    
    Returns
    -------
    Returns the gas volume density in g/cm^3
    """
    # Mass of gas per H2-molecule
    mgas = mp*(2.0*ppar['abun_h2']+4*ppar['abun_he'])/ppar['abun_h2']
    rhogas = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) + ppar['nh2'] * mgas
    return rhogas


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
    rhogas = getGasDensity(grid=grid, ppar=ppar)
    rhodust = np.zeros([grid.nx, grid.ny, grid.nz, 1], dtype=np.float64) 
    rhodust[:, :, :, 0] = rhogas * ppar['dusttogas']
    return rhodust


def getVTurb(grid=None, ppar=None):
    """Calculates/sets the turbulent velocity field
    
    Parameters
    ----------
    grid : radmc3dGrid
            An instance of the radmc3dGrid class containing the spatial and wavelength grid
    
    ppar : dictionary
            Dictionary containing all parameters of the model 
    
    Returns
    -------
    Returns the turbulent velocity in cm/s
    """
    vturb = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) + ppar['gasspec_vturb']
    return vturb


def getVelocity(grid=None, ppar=None):
    """Calculates/sets the gas velocity field
    
    Parameters
    ----------
    grid : radmc3dGrid
            An instance of the radmc3dGrid class containing the spatial and wavelength grid
    
    ppar : dictionary
            Dictionary containing all parameters of the model 
    
    Returns
    -------
    Returns the turbulent velocity in cm/s
    """
    vel = np.zeros([grid.nx, grid.ny, grid.nz, 3], dtype=np.float64)
    for ix in range(grid.nx):
        vel[ix, :, :, 0] = ppar['dvdau']*grid.x[ix]/au
    
    return vel
