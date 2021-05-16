"""Generic protoplanetary disk model with a gap - octree AMR example

The density is given by 

    .. math::
        
        \\rho = \\frac{\\Sigma(r,\\phi)}{H_p\\sqrt{(2\\pi)}} \\exp{\\left(-\\frac{z^2}{2H_p^2}\\right)}


    * :math:`\Sigma` - surface density
    * :math:`H_{\\rm p}` - Pressure scale height

There are two options for the functional form of surface density as a function of radius. For a simple
power-law the surface density is given by

    * :math:`\Sigma_{\\rm c}(r) = \\Sigma_0\\left(\\frac{r}{r_{\\rm out}}\\right)^p`

alternatively the surface density can also have an exponential outer tapering:

    * :math:`\Sigma_{\\rm c}(r) = \\Sigma_0\\left(\\frac{r}{r_{\\rm out}}\\right)^p\\exp{\\left\\{-\\left(\\frac{r}{r_{\\rm out}}\\right)^{2+p}\\right\\}}`

The :math:`{\\rm c}` index refers to the continuous, unperturbed surface density as the final surface density may 
contain perturbations due to the presence of gaps. Arbitrary number of gaps can be placed in the surface density of 
the disk. The gap is implemented as a radial gaussian depletion in the form 

    * :math:`G(r) = \\frac{1}{\\delta-1}\exp{(-\\left((r-r_{\\rm c})^2/\\sigma^2\\right)}`

where :math:`\\delta` is the density recution in the center of the gap, :math:`r_{\\rm c}` is the distance of the gap 
center from the star, and :math:`\\sigma` is the standard deviation of the gaussian in the radial direction. The final 
surface density including the gaps will then be calculated as 

    * :math:`\\Sigma(r) =\\Sigma_{\\rm c}(r) / (1 + \\sum_{i=1}^{N}G_{\\rm i}(r))`

The molecular abundance function takes into account the freeze-out of the molecules
For freeze-out the molecular abundance below a threshold temperature is decreased by a given fractor. 


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
from .. import analyze


def getModelDesc():
    """Returns the brief description of the model.
    """

    return "Generic protoplanetary disk model with a gap - octree AMR example"
           

def getDefaultParams():
    """Function to provide default parameter values of the model.

    Returns a list whose elements are also lists with three elements:
    1) parameter name, 2) parameter value, 3) parameter description
    All three elements should be strings. The string of the parameter
    value will be directly written out to the parameter file if requested,
    and the value of the string expression will be evaluated and be put
    to radmc3dData.ppar. The third element contains the description of the
    parameter which will be written in the comment field of the line when
    a parameter file is written. 
    """

    defpar = [
        ['xres_nlev', '3', 'Number of refinement levels'],
        ['xres_nspan', '3', 'Number of the original grid cells to refine'],
        ['xres_nstep', '3', 'Number of grid cells to create in a refinement level'],
        ['crd_sys', "'car'", 'Coordinate system used (car/cyl)'],
        ['grid_style', '1', '0 - Regular grid, 1 - Octree AMR, 10 - Layered/nested grid (not yet supported)'],
        ['nx', '[6]', 'Number of grid points in the first dimension'],
        ['xbound', '[-100.*au, 100.*au]', 'Number of radial grid points'],
        ['ny', '[6]', 'Number of grid points in the first dimension'],
        ['ybound', '[-100.*au, 100.*au]', 'Number of radial grid points'],
        ['nz', '[6]', 'Number of grid points in the first dimension'],
        ['zbound', '[-100.*au, 100.*au]', 'Number of radial grid points'],
        ['gasspec_mol_name', "['co']", ''],
        ['gasspec_mol_abun', '[1e-4]', ''],
        ['gasspec_mol_dbase_type', "['leiden']", ''],
        ['gasspec_mol_freezeout_temp', '[19.0]', 'Freeze-out temperature of the molecules in Kelvin'],
        ['gasspec_mol_freezeout_dfact', '[1e-3]',
         'Factor by which the molecular abundance should be decreased in the frezze-out zone'],
        ['gasspec_vturb', '0.2e5', 'Microturbulent line width'],
        ['rin', '1.0*au', ' Inner radius of the disk'],
        ['rdisk', '100.0*au', ' Outer radius of the disk'],
        ['hrdisk', '0.1', ' Ratio of the pressure scale height over radius at hrpivot'],
        ['hrpivot', "100.0*au", ' Reference radius at which Hp/R is taken'],
        ['plh', '1./7.', ' Flaring index'],
        ['plsig1', '-1.0', ' Power exponent of the surface density distribution as a function of radius'],
        ['sig0', '0.3', ' Surface density at rdisk'],
        ['mdisk', '0.0', ' Mass of the disk (either sig0 or mdisk should be set to zero or commented out)'],
        ['bgdens', '1e-30', ' Background density (g/cm^3)'],
        ['srim_rout', '0.0', 'Outer boundary of the smoothing in the inner rim in terms of rin'],
        ['srim_plsig', '0.0', 'Power exponent of the density reduction inside of srim_rout*rin'],
        ['prim_rout', '0.0', 'Outer boundary of the puffed-up inner rim in terms of rin'],
        ['hpr_prim_rout', '0.0', 'Pressure scale height at rin'],
        ['gap_rc', '[50e0*au]', ' Radial center of the gap'],
        ['gap_sigma', '[5e0*au]', ' Standard deviation of the gap in the radial direction'],
        ['gap_drfact', '[1e-6]', ' Density reduction factor in the gap'],
        ['sigma_type', '0',
         ' Surface density type (0 - polynomial, 1 - exponential outer edge (viscous self-similar solution)'],
        ['dusttogas', '0.01', ' Dust-to-gas mass ratio']]

    return defpar


def getDustDensity(x=None, y=None, z=None, ppar=None, grid=None):
    """Calculates the dust density distribution in a protoplanetary disk.
   
    Parameters
    ----------
    x    : ndarray
           Coordinate of the cell centers in the first dimension
    
    y    : ndarray
           Coordinate of the cell centers in the second dimension
    
    z    : ndarray
           Coordinate of the cell centers in the third dimension

    grid : radmc3dGrid
           An instance of the radmc3dGrid class containing the spatial and frequency/wavelength grid
    
    ppar : dictionary
           A dictionary containing all parameters of the model
    
    Returns
    -------
    Returns the volume density in g/cm^3
    """
    if grid is not None:
        x = grid.x
        y = grid.y
        z = grid.z

# Get the gas density
    rhogas = getGasDensity(x=x, y=y, z=z, ppar=ppar)
    rho = np.array(rhogas) * ppar['dusttogas']

# Split up the disk density distribution according to the given abundances
    if 'ngs' in ppar:
        if ppar['ngs'] > 1:
            ngs = ppar['ngs']
            #
            # WARNING!!!!!!
            # At the moment I assume that the multiple dust population differ from each other only in 
            # grain size but not in bulk density thus when I calculate the abundances / mass fractions 
            # they are independent of the grains bulk density since abundances/mass fractions are normalized
            # to the total mass. Thus I use 1g/cm^3 for all grain sizes.
            # TODO: Add the possibility to handle multiple dust species with different bulk densities and 
            # with multiple grain sizes.
            #
            gdens = np.zeros(ngs, dtype=float) + 1.0
            gs = ppar['gsmin'] * (ppar['gsmax']/ppar['gsmin'])**(np.arange(ppar['ngs'], dtype=np.float64)
                                                                 / (float(ppar['ngs'])-1.))
            gmass = 4./3.*np.pi*gs**3. * gdens
            gsfact = gmass * gs**(ppar['gsdist_powex']+1)
            gsfact = gsfact / gsfact.sum()
        else:
            gsfact = [1.0]
            ngs = 1
    elif 'mfrac' in ppar:
        ngs = len(ppar['mfrac'])
        gsfact = ppar['mfrac'] / ppar['mfrac'].sum()
    
    else:
        ngs = 1
        gsfact = [1.0]
    
    rho_old = np.array(rho)
    rho = rho_old*gsfact

    return rho


def getGasDensity(x=None, y=None, z=None, ppar=None, grid=None):
    """Calculates the gas density distribution in a protoplanetary disk.
    
    Parameters
    ----------
    x    : ndarray
           Coordinate of the cell centers in the first dimension
    
    y    : ndarray
           Coordinate of the cell centers in the second dimension
    
    z    : ndarray
           Coordinate of the cell centers in the third dimension

    grid : radmc3dGrid
           An instance of the radmc3dGrid class containing the spatial and frequency/wavelength grid
    
    ppar : dictionary
           A dictionary containing all parameters of the model
    
    Returns
    -------
    Returns the volume density in g/cm^3
    """
   
    if grid is not None:
        x = grid.x
        y = grid.y
        z = grid.z

    if ppar['crd_sys'] == 'car':
        rcyl = (x**2 + y**2)**0.5
        zz = z
    elif ppar['crd_sys'] == 'sph':
        rcyl = x * np.sin(y)
        zz = x * np.cos(x)
    else:
        raise ValueError('Unknown crd_sys : ', ppar['crd_sys'], '\n Allowed values are "car", "sph"')

    # Calculate the pressure scale height as a function of r, phi
    hp = ppar['hrdisk'] * (rcyl/ppar['hrpivot'])**ppar['plh'] * rcyl

    if 'prim_rout' in ppar:
        if ppar['prim_rout'] >= 1.:
            dum_hrdisk = ppar['hrdisk'] * (rcyl/ppar['hrpivot'])**ppar['plh'] 
            hpr0 = ppar['hrdisk'] * (ppar['prim_rout'] * ppar['rin']/ppar['hrpivot'])**ppar['plh']
            dummy = np.log10(hpr0 / ppar['hpr_prim_rout']) / np.log10(ppar['prim_rout'])
            dum_prim = ppar['hpr_prim_rout'] * (rcyl/ppar['rin'])**dummy
            hp = (dum_hrdisk**8. + dum_prim**8.)**(1./8.) * rcyl

    # Calculate the surface density 
    sigma = 0.0

    # Calculate sigma from sig0, rdisk and plsig1
    if 'sig0' in ppar:
        if ppar['sig0'] != 0.:
            if 'sigma_type' in ppar:
                if ppar['sigma_type'] == 0:
                    dum1 = ppar['sig0'] * (rcyl/ppar['rdisk'])**ppar['plsig1']
                else:
                    expterm = np.exp(-(rcyl/ppar['rdisk'])**(2.0 - ppar['plsig1']))
                    dum1 = ppar['sig0'] * (rcyl/ppar['rdisk'])**(-ppar['plsig1']) * expterm 

            else:
                dum1 = ppar['sig0'] * (rcyl/ppar['rdisk'])**ppar['plsig1']

            if ('srim_rout' in ppar) & ('srim_plsig' in ppar):
                if ppar['srim_rout'] != 0.:
                    
                    if 'sigma_type' in ppar:
                        if ppar['sigma_type'] == 0:
                            # Adding the smoothed inner rim
                            sig_srim = ppar['sig0'] * (ppar['srim_rout']*ppar['rin'] / ppar['rdisk'])**ppar['plsig1']
                            dum2 = sig_srim * (rcyl / (ppar['srim_rout']*ppar['rin']))**ppar['srim_plsig']
                        else:
                            # sig_srim = 1.0 * (ppar['srim_rout']*ppar['rin'] / ppar['rdisk'])**ppar['plsig1']
                            sig_srim = ppar['sig0'] * (ppar['srim_rout']*ppar['rin']
                                                       / ppar['rdisk'])**(-ppar['plsig1']) \
                                       * np.exp(-(rcyl/ppar['rdisk'])**(2.0 - ppar['plsig1']))
                            dum2 = sig_srim * (rcyl / (ppar['srim_rout']*ppar['rin']))**ppar['srim_plsig']
                    else:
                        # Adding the smoothed inner rim
                        sig_srim = ppar['sig0'] * (ppar['srim_rout']*ppar['rin'] / ppar['rdisk'])**ppar['plsig1']
                        dum2 = sig_srim * (rcyl / (ppar['srim_rout']*ppar['rin']))**ppar['srim_plsig']

                    p = -5.0
                    dum = (dum1**p + dum2**p)**(1./p)
                else:
                    dum = dum1
            else:
                dum = dum1

            sigma = dum
        
        else:
            if 'sigma_type' in ppar:
                if ppar['sigma_type'] == 0:
                    dum1 = 1.0 * (rcyl/ppar['rdisk'])**ppar['plsig1']
                else:
                    dum1 = 1.0 * (rcyl/ppar['rdisk'])**(-ppar['plsig1']) * \
                           np.exp(-(rcyl/ppar['rdisk'])**(2.0 - ppar['plsig1']))
            else:
                dum1 = 1.0 * (rcyl/ppar['rdisk'])**ppar['plsig1']

            if ('srim_rout' in ppar) & ('srim_plsig' in ppar):
                if ppar['srim_rout'] != 0.:

                    if 'sigma_type' in ppar:
                        if ppar['sigma_type'] == 0:
                            # Adding the smoothed inner rim
                            sig_srim = 1.0 * (ppar['srim_rout']*ppar['rin'] / ppar['rdisk'])**ppar['plsig1']
                            dum2 = sig_srim * (rcyl / (ppar['srim_rout']*ppar['rin']))**ppar['srim_plsig']
                        else:
                            # sig_srim = 1.0 * (ppar['srim_rout']*ppar['rin'] / ppar['rdisk'])**ppar['plsig1']
                            sig_srim = 1.0 * (ppar['srim_rout']*ppar['rin']
                                              / ppar['rdisk'])**(-ppar['plsig1']) \
                                       * np.exp(-(rcyl/ppar['rdisk'])**(2.0 - ppar['plsig1']))
                            dum2 = sig_srim * (rcyl / (ppar['srim_rout']*ppar['rin']))**ppar['srim_plsig']
                    else:
                        # Adding the smoothed inner rim
                        sig_srim = 1.0 * (ppar['srim_rout']*ppar['rin'] / ppar['rdisk'])**ppar['plsig1']
                        dum2 = sig_srim * (rcyl / (ppar['srim_rout']*ppar['rin']))**ppar['srim_plsig']

                    # # Adding the smoothed inner rim
                    # sig_srim = 1.0 * (ppar['srim_rout']*ppar['rin'] / ppar['rdisk'])**ppar['plsig1']
                    # dum2 = sig_srim * (rcyl / (ppar['srim_rout']*ppar['rin']))**ppar['srim_plsig']

                    p = -5.0
                    dum = (dum1**p + dum2**p)**(1./p)
                else:
                    dum = dum1
            else:
                dum = dum1

            sigma = dum

        if ppar['sigma_type'] == 0.:
            ii = ((rcyl < ppar['rin']) | (rcyl > ppar['rdisk']))
            if True in ii:
                sigma[ii] = 0.

        else:
            ii = (rcyl < ppar['rin'])
            if True in ii:
                sigma[ii] = 0.0

    #
    # Add the gaps
    #

    gap_fact = np.zeros(rcyl.shape[0], dtype=np.float64)
    for igap in range(len(ppar['gap_rc'])):
        gap_fact += np.exp(-0.5 * (rcyl-ppar['gap_rc'][igap])**2 / ppar['gap_sigma'][igap]**2) \
                    * (1. / ppar['gap_drfact'][igap]-1.0)
    
    gap_fact += 1.0
    sigma /= gap_fact

    z0 = 0.0
    rho = sigma / np.sqrt(2.0*np.pi) / hp * np.exp(-0.5 * (zz-z0)**2 / hp**2) + ppar['bgdens']

    return rho


def getGasAbundance(x=None, y=None, z=None, grid=None, ppar=None, ispec=''):
    """Calculates the molecular abundance. 
    
    The number density of a molecule is rhogas * abun 
   
    Parameters
    ----------
    x    : ndarray
           Coordinate of the cell centers in the first dimension
    
    y    : ndarray
           Coordinate of the cell centers in the second dimension
    
    z    : ndarray
           Coordinate of the cell centers in the third dimension

    grid  : radmc3dGrid
            An instance of the radmc3dGrid class containing the spatial and wavelength grid
    
    ppar  : dictionary
            Dictionary containing all parameters of the model 
    
    ispec : str
            The name of the gas species whose abundance should be calculated

    Returns
    -------
    Returns an ndarray containing the molecular abundance at each grid point
    """
    if grid is not None:
        x = grid.x

    nx = x.shape[0]

    if 'grid_style' not in ppar:
        raise KeyError('"grid_style" is missing from the parameter dictionary. '
                       + 'It is either missing from problem_params.inp or its value cannot be evaluated')

    if 'crd_sys' not in ppar:
        raise KeyError('"crd_sys" is missing from the parameter dictionary. '
                       + 'It is either missing from problem_params.inp or its value cannot be evaluated')

    if (ppar['crd_sys'] != 'car') & (ppar['crd_sys'] == 'sph'):
        raise ValueError('Unknown crd_sys : ', ppar['crd_sys'], '\n Allowed values are "car", "sph"')

    #
    # Regular grids
    #
    if ppar['grid_style'] == 0:
        # Read the dust density and temperature
        try:
            data = analyze.radmc3dData(grid)
            data.readDustDens()
            data.readDustTemp()
        except:
            raise RuntimeError('Gas abundance cannot be calculated as the required dust density and/or temperature '
                               + 'could not be read in binary or in formatted ascii format.')

        nspec = len(ppar['gasspec_mol_name'])
        if ispec in ppar['gasspec_mol_name']:

            sid = ppar['gasspec_mol_name'].index(ispec)
            gasabun = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)  
            
            for spec in range(nspec):
                gasabun[:, :, :] = ppar['gasspec_mol_abun'][sid]
               
            for iz in range(data.grid.nz):
                for iy in range(data.grid.ny):
                    ii = (data.taux[:, iy, iz] < ppar['gasspec_mol_dissoc_taulim'][sid])
                    gasabun[ii, iy, iz] = 1e-90

                    ii = (data.dusttemp[:, iy, iz, 0] < ppar['gasspec_mol_freezeout_temp'][sid])
                    gasabun[ii, iy, iz] = ppar['gasspec_mol_abun'][sid] * ppar['gasspec_mol_freezeout_dfact'][sid]

        else:
            gasabun = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) + 1e-10
            print('WARNING !!!')
            print('Molecule name "'+ispec+'" is not found in gasspec_mol_name')
            print('A default 1e-10 abundance will be used')

    #
    # Octree AMR
    #
    elif ppar['grid_style'] == 1:

        # Read the dust density and temperature
        try:
            data = analyze.radmc3dData(grid)
            data.readDustTemp()
            data.dusttemp = grid.convArrLeaf2Tree(data.dusttemp)
        except:
            raise RuntimeError('Gas abundance cannot be calculated as the required dust density and/or temperature '
                               + 'could not be read in binary or in formatted ascii format.')

        nspec = len(ppar['gasspec_mol_name'])
        if ispec in ppar['gasspec_mol_name']:

            sid = ppar['gasspec_mol_name'].index(ispec)
            gasabun = np.zeros(nx, dtype=np.float64)
            
            for spec in range(nspec):
                gasabun[:] = ppar['gasspec_mol_abun'][sid]

            ii = (data.dusttemp[:, 0] < ppar['gasspec_mol_freezeout_temp'][sid])

            if True in ii:
                gasabun[ii] = ppar['gasspec_mol_abun'][sid] * ppar['gasspec_mol_freezeout_dfact'][sid]

        else:
            gasabun = np.zeros(nx, dtype=np.float64) + 1e-10
            print('WARNING !!!')
            print('Molecule name "'+ispec+'" is not found in gasspec_mol_name')
            print('A default 1e-10 abundance will be used')
    else:
        raise ValueError('Unknown grid style : ', ppar['grid_sytle'])

    return gasabun


def getVTurb(x=None, y=None, z=None, grid=None, ppar=None):
    """Calculates the turbulent velocity field
    
    Parameters
    ----------
    x    : ndarray
           Coordinate of the cell centers in the first dimension
    
    y    : ndarray
           Coordinate of the cell centers in the second dimension
    
    z    : ndarray
           Coordinate of the cell centers in the third dimension

    grid : radmc3dGrid
           An instance of the radmc3dGrid class containing the spatial and frequency/wavelength grid
    
    ppar : dictionary
           A dictionary containing all parameters of the model
    
    Returns
    -------
    Returns an ndarray with the turbulent velocity in cm/s
    """

    if grid is not None:
        x = grid.x
        y = grid.y
        z = grid.z

    nx = x.shape[0]
    ny = y.shape[0]
    nz = z.shape[0]

    if 'grid_style' not in ppar:
        raise KeyError('"grid_style" is missing from the parameter dictionary. '
                       + 'It is either missing from problem_params.inp or its value cannot be evaluated')

#
# Regular grid
#
    if ppar['grid_style'] == 0:
        vturb = np.zeros([nx, ny, nz], dtype=np.float64) + ppar['gasspec_vturb']
#
# Octree AMR
#
    elif ppar['grid_style'] == 1:
        vturb = np.zeros(nx, dtype=np.float64) + ppar['gasspec_vturb']
        
    else:
        raise ValueError('Unknown grid style : ', ppar['grid_sytle'])

    return vturb


def getVelocity(x=None, y=None, z=None, grid=None, ppar=None):
    """Calculates the velocity field in a protoplanetary disk.
    
    Parameters
    ----------
    x    : ndarray
           Coordinate of the cell centers in the first dimension
    
    y    : ndarray
           Coordinate of the cell centers in the second dimension
    
    z    : ndarray
           Coordinate of the cell centers in the third dimension

    grid : radmc3dGrid
           An instance of the radmc3dGrid class containing the spatial and frequency/wavelength grid
    
    ppar : dictionary
           A dictionary containing all parameters of the model
    
    Returns
    -------
    Returns the gas velocity in cm/s
    """

    if grid is not None:
        x = grid.x
        y = grid.y
        z = grid.z

    nx = x.shape[0]
    ny = y.shape[0]
    nz = z.shape[0]

    if 'grid_style' not in ppar:
        raise KeyError('"grid_style" is missing from the parameter dictionary. '
                       + 'It is either missing from problem_params.inp or its value cannot be evaluated')

    if 'crd_sys' not in ppar:
        raise KeyError('"crd_sys" is missing from the parameter dictionary. '
                       + 'It is either missing from problem_params.inp or its value cannot be evaluated')

    #
    # Regular grid
    #
    if ppar['grid_style'] == 0:
        if ppar['crd_sys'] == 'car':
            xx, yy = np.meshgrid(x, y, indexing='ij')
            rcyl = np.sqrt(xx**2 + yy**2)
            cos_phi = xx/rcyl
            sin_phi = yy/rcyl
            
            vel = np.zeros([nx, ny, nz, 3], dtype=np.float64)
            vkep = np.sqrt(gg*ppar['mstar'][0]/rcyl)
            for iz in range(nz):
                vel[:, :, iz, 0] = -vkep*sin_phi
                vel[:, :, iz, 1] = vkep*cos_phi

        elif ppar['crd_sys'] == 'sph':
            nr = nx
            nphi = nz
            nz = ny
            rcyl = x

            vel = np.zeros([nr, nz, nphi, 3], dtype=np.float64)
            vkep = np.sqrt(gg*ppar['mstar'][0]/rcyl)
            for iz in range(nz):
                for ip in range(nphi):
                    vel[:, iz, ip, 2] = vkep
        else:
            raise ValueError('Unknown coordinate system : ', ppar['crd_sys'])

    #
    # Octree AMR
    #
    elif ppar['grid_style'] == 1:
        if ppar['crd_sys'] == 'car':
            rcyl = np.sqrt(x**2 + y**2)
            cos_phi = x/rcyl
            sin_phi = y/rcyl
            vkep = np.sqrt(gg*ppar['mstar'][0]/rcyl)
            vel = np.zeros([nx, 3], dtype=np.float64)
            vel[:, 0] = -vkep*sin_phi
            vel[:, 1] = vkep*cos_phi
        
        elif ppar['crd_sys'] == 'sph':
            rcyl = x * np.sin(y)
            vkep = np.sqrt(gg*ppar['mstar'][0]/rcyl)
            vel = np.zeros([nx, 3], dtype=np.float64)
            vel[:, 2] = vkep

        else:
            raise ValueError('Unknown coordinate system : ', ppar['crd_sys'])

    else:
        raise ValueError('Unknown grid style : ', ppar['grid_sytle'])

    return vel


def decisionFunction(x=None, y=None, z=None, dx=None, dy=None, dz=None, model=None, ppar=None, **kwargs):
    """
    Example function to be used as decision function for resolving cells in tree building. It calculates the gas 
    density at a random sample of coordinates within a given cell than takes the quantity 
    :math:`(\max{\\rho_{\\rm i}} - \min{\\rho_\\rm{i}})/\max{\\rho_{\\rm i}}`. If it is larger than a certain threshold
    value it will return True (i.e. the cell should be resolved) if the density variation is less than the threshold 
    it returns False (i.e. the cell should not be resolved)

    Parameters
    ----------
    x       : ndarray
              Cell centre coordinates of the cells in the first dimension
    
    y       : ndarray
              Cell centre coordinates of the cells in the second dimension
    
    z       : ndarray
              Cell centre coordinates of the cells in the third dimension
    
    dx      : ndarray
              Half size of the cells in the first dimension
    
    dy      : ndarray
              Half size of the cells in the second dimension
    
    dz      : ndarray
              Half size of the cells in the third dimension
    
    model   : object
              A radmc3dPy model (must contain a getGasDensity() function) 

    ppar    : dictionary
              All parameters of the problem (from the problem_params.inp file). It is not used here, but must be present 
              for compatibility reasons.
    
    **kwargs: dictionary
              Parameters used to decide whether the cell should be resolved. It should contain the following keywords; 
              'nsample', which sets the number of random points the gas desity is sampled at within the cell and 
              'threshold' that sets the threshold value for max(gasdens)/min(gasdens) above which the cell should be 
              resolved.
    """

    ncell = x.shape[0]
    rho = np.zeros([ncell, kwargs['nsample']], dtype=np.float64)

    for isample in range(kwargs['nsample']):
        xoffset = (np.random.random_sample(ncell) - 0.5)*dx*4.0
        yoffset = (np.random.random_sample(ncell) - 0.5)*dy*4.0
        zoffset = (np.random.random_sample(ncell) - 0.5)*dz*4.0
        rho[:, isample] = model.getGasDensity(x+xoffset, y+yoffset, z+zoffset, ppar=ppar)
    
    rho_max = rho.max(axis=1)
    rho_min = rho.min(axis=1)
    jj = ((rho_max-rho_min)/rho_max > ppar['threshold'])
    
    decision = np.zeros(ncell, dtype=bool)
    if True in jj:
        decision[jj] = True

    return decision
