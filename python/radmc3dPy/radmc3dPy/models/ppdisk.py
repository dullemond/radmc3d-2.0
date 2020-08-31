"""Generic protoplanetary disk model 

The density is given by 

    .. math::
        
        \\rho = \\frac{\\Sigma(r,\\phi)}{H_p\\sqrt{(2\\pi)}} \\exp{\\left(-\\frac{z^2}{2H_p^2}\\right)}


    * :math:`\Sigma` - surface density
    * :math:`H_{\\rm p}` - Pressure scale height

There are two options for the functional form of surface density as a function of radius. For a simple
power-law the surface density is given by

    * :math:`\Sigma(r) = \\Sigma_0\\left(\\frac{r}{r_{\\rm out}}\\right)^p`

alternatively the surface density can also have an exponential outer tapering:

    * :math:`\Sigma(r) = \\Sigma_0\\left(\\frac{r}{r_{\\rm out}}\\right)^p\\exp{\\left\\{-\\left(\\frac{r}{r_{\\rm out}}\\right)^{2+p}\\right\\}}`


The molecular abundance function takes into account dissociation and freeze-out of the molecules
For photodissociation only the continuum (dust) shielding is taken into account in a way that
whenever the continuum optical depth radially drops below a threshold value the molecular abundance
is dropped to zero. For freeze-out the molecular abundance below a threshold temperature is decreased
by a given fractor. 


"""
from __future__ import absolute_import
from __future__ import print_function
import warnings
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

    return "Generic protoplanetary disk model"
           

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
        ['nx', '[30,50]', 'Number of grid points in the first dimension'],
        ['xbound', '[1.0*au,1.05*au, 100.0*au]', 'Number of radial grid points'],
        ['ny', '[10,30,30,10]', 'Number of grid points in the first dimension'],
        ['ybound', '[0., pi/3., pi/2., 2.*pi/3., pi]', 'Number of radial grid points'],
        ['nz', '30', 'Number of grid points in the first dimension'],
        ['zbound', '[0., 2.0*pi]', 'Number of radial grid points'],
        ['gasspec_mol_name', "['co']", ''],
        ['gasspec_mol_abun', '[1e-4]', ''],
        ['gasspec_mol_dbase_type', "['leiden']", ''],
        ['gasspec_mol_dissoc_taulim', '[1.0]', 'Continuum optical depth limit below which all molecules dissociate'],
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
        ['sig0', '0.0', ' Surface density at rdisk'],
        ['mdisk', '1e-3*ms', ' Mass of the disk (either sig0 or mdisk should be set to zero or commented out)'],
        ['bgdens', '1e-30', ' Background density (g/cm^3)'],
        ['srim_rout', '0.0', 'Outer boundary of the smoothing in the inner rim in terms of rin'],
        ['srim_plsig', '0.0', 'Power exponent of the density reduction inside of srim_rout*rin'],
        ['prim_rout', '0.0', 'Outer boundary of the puffed-up inner rim in terms of rin'],
        ['hpr_prim_rout', '0.0', 'Pressure scale height at rin'],
        ['gap_rin', '[0e0*au]', ' Inner radius of the gap'],
        ['gap_rout', '[0e0*au]', ' Outer radius of the gap'],
        ['gap_drfact', '[0e0]', ' Density reduction factor in the gap'],
        ['sigma_type', '0',
         ' Surface density type (0 - polynomial, 1 - exponential outer edge (viscous self-similar solution)'],
        ['dusttogas', '0.01', ' Dust-to-gas mass ratio']]

    return defpar


def getDustDensity(grid=None, ppar=None):
    """Calculates the dust density distribution in a protoplanetary disk.
   
    Parameters
    ----------
    grid : radmc3dGrid
           An instance of the radmc3dGrid class containing the spatial and frequency/wavelength grid
    
    ppar : dictionary
           A dictionary containing all parameters of the model
    
    Returns
    -------
    Returns the volume density in g/cm^3
    """

# Get the gas density
    rhogas = getGasDensity(grid=grid, ppar=ppar)

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
    
    # if ppar.has_key('dustkappa_ext'):
        # ngs  = len(ppar['dustkappa_ext'])
        # if ppar.has_key('mfrac'):
        #     gsfact = ppar['mfrac'] / ppar['mfrac'].sum()
        # else:
        #     ngs = 1
        #     gsfact = [1.0]

    # else:
        # ngs = ppar['ngs']
        #
        # WARNING!!!!!!
        # At the moment I assume that the multiple dust population differ from each other only in
        # grain size but not in bulk density thus when I calculate the abundances / mass fractions
        # they are independent of the grains bulk density since abundances/mass fractions are normalized
        # to the total mass. Thus I use 1g/cm^3 for all grain sizes.
        # TODO: Add the possibility to handle multiple dust species with different bulk densities and
        # with multiple grain sizes.
        #
        # gdens = zeros(ngs, dtype=float) + 1.0
        # gs = ppar['gsmin'] * (ppar['gsmax']/ppar['gsmin'])**(arange(ppar['ngs'], dtype=float64)
        # / (float(ppar['ngs'])-1.))
        # gmass = 4./3.*np.pi*gs**3. * gdens
        # gsfact = gmass * gs**(ppar['gsdist_powex']+1)
        # gsfact = gsfact / gsfact.sum()

    rho_old = np.array(rho)
    rho = np.zeros([grid.nx, grid.ny, grid.nz, ngs], dtype=np.float64)
    for igs in range(ngs):
        rho[:, :, :, igs] = rho_old[:, :, :] * gsfact[igs]

    return rho


def getGasDensity(grid=None, ppar=None):
    """Calculates the gas density distribution in a protoplanetary disk.
    
    Parameters
    ----------
    grid : radmc3dGrid
           An instance of the radmc3dGrid class containing the spatial and frequency/wavelength grid
    
    ppar : dictionary
           A dictionary containing all parameters of the model
    
    Returns
    -------
    Returns the volume density in g/cm^3
    """
    
    rr, th = np.meshgrid(grid.x, grid.y)
    zz = rr * np.cos(th)
    rcyl = rr * np.sin(th)

    # Calculate the pressure scale height as a function of r, phi
    hp = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
    dum = ppar['hrdisk'] * (rcyl/ppar['hrpivot'])**ppar['plh'] * rcyl

    if 'prim_rout' in ppar:
        if ppar['prim_rout'] >= 1.:
            dum_hrdisk = ppar['hrdisk'] * (rcyl/ppar['hrpivot'])**ppar['plh'] 
            hpr0 = ppar['hrdisk'] * (ppar['prim_rout'] * ppar['rin']/ppar['hrpivot'])**ppar['plh']
            dummy = np.log10(hpr0 / ppar['hpr_prim_rout']) / np.log10(ppar['prim_rout'])
            dum_prim = ppar['hpr_prim_rout'] * (rcyl/ppar['rin'])**dummy
            dum = (dum_hrdisk**8. + dum_prim**8.)**(1./8.) * rcyl

    dum = dum.swapaxes(0, 1)
    for iz in range(grid.nz):
        hp[:, :, iz] = dum

    # Calculate the surface density 
    sigma = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
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

            dum = dum.swapaxes(0, 1)

            for iz in range(grid.nz):
                sigma[:, :, iz] = dum
        else:
            if 'sigma_type' in ppar:
                if ppar['sigma_type'] == 0:
                    dum1 = 1.0 * (rcyl/ppar['rdisk'])**ppar['plsig1']
                else:
                    dum1 = 1.0 * (rcyl/ppar['rdisk'])**(-ppar['plsig1']) \
                           * np.exp(-(rcyl/ppar['rdisk'])**(2.0 - ppar['plsig1']))
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
                            sig_srim = 1.0 * (ppar['srim_rout']*ppar['rin'] / ppar['rdisk'])**(-ppar['plsig1']) \
                                       * np.exp(-(rcyl/ppar['rdisk'])**(2.0 - ppar['plsig1']))
                            dum2 = sig_srim * (rcyl / (ppar['srim_rout']*ppar['rin']))**ppar['srim_plsig']
                    else:
                        # Adding the smoothed inner rim
                        sig_srim = 1.0 * (ppar['srim_rout']*ppar['rin'] / ppar['rdisk'])**ppar['plsig1']
                        dum2 = sig_srim * (rcyl / (ppar['srim_rout']*ppar['rin']))**ppar['srim_plsig']

                    p = -5.0
                    dum = (dum1**p + dum2**p)**(1./p)
                else:
                    dum = dum1
            else:
                dum = dum1

            dum = dum.swapaxes(0, 1)

            for iz in range(grid.nz):
                sigma[:, :, iz] = dum

        if 'sigma_type' in ppar:
            if ppar['sigma_type'] == 0:
                for iy in range(grid.ny):
                    ii = (rcyl[iy, :] < ppar['rin']) | (rcyl[iy, :] > ppar['rdisk'])
                    sigma[ii, iy, :] = 0.0
        else:
            for iy in range(grid.ny):
                ii = (rcyl[iy, :] < ppar['rin']) | (rcyl[iy, :] > ppar['rdisk'])
                sigma[ii, iy, :] = 0.0


    z0 = np.zeros([grid.nx, grid.nz, grid.ny], dtype=np.float64)
    rho = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
    for iz in range(grid.nz):
        for iy in range(grid.ny):
            rho[:, iy, iz] = sigma[:, iy, iz] / (hp[:, iy, iz] * np.sqrt(2.0*np.pi)) * \
                              np.exp(-0.5 * ((zz[iy, :])-z0[:, iz, iy])*((zz[iy, :])-z0[:, iz, iy])
                                     / (hp[:, iy, iz]*hp[:, iy, iz])) + ppar['bgdens']




    # Normalize the disk to mdisk if it is given instead of sig0
    if 'mdisk' in ppar:
        if ppar['mdisk'] != 0.:
            # Calculate the volume of each grid cell
            vol = grid.getCellVolume()
            mass = (rho * vol).sum(0).sum(0).sum(0)
            rho = rho * (ppar['mdisk'] / mass)

            if np.abs(ppar['ybound'][-1]-(np.pi/2.)) < 1e-8:
                rho = rho*0.5

    for igap in range(len(ppar['gap_rout'])):
        for ix in range(grid.nx):
            if (grid.x[ix] >= ppar['gap_rin'][igap]) & (grid.x[ix] <= ppar['gap_rout'][igap]):
                rho[ix, :, :] = rho[ix, :, :] * ppar['gap_drfact'][igap]
    return rho


def getGasAbundance(grid=None, ppar=None, ispec=''):
    """Calculates the molecular abundance. 
    
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
    Returns an ndarray containing the molecular abundance at each grid point
    """

    # Read the dust density and temperature
    try: 
        data = analyze.readData(ddens=True, dtemp=True, binary=True)
    except:
        try: 
            data = analyze.readData(ddens=True, dtemp=True, binary=False)
        except:
            msg = 'Gas abundance cannot be calculated as the required dust density and/or temperature '\
                  + 'could not be read in binary or in formatted ascii format.'
            raise RuntimeError(msg)

    # Calculate continuum optical depth 
    data.getTau(axis='xy', wav=0.55)
    
    nspec = len(ppar['gasspec_mol_name'])
    if ppar['gasspec_mol_name'].__contains__(ispec):

        sid = ppar['gasspec_mol_name'].index(ispec)
        # Check where the radial and vertical optical depth is below unity
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
        txt = 'Molecule name "'+ispec+'" is not found in gasspec_mol_name \n A default 1e-10 abundance will be used'
        warnings.warn(txt, RuntimeWarning)

    # gasabun = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
    # gasabun[:,:,:] = ppar['gasspec_mol_abun'][0] / (2.4*mp)

    return gasabun


def getVTurb(grid=None, ppar=None):
    """Calculates the turbulent velocity field
    
    Parameters
    ----------
    grid : radmc3dGrid
           An instance of the radmc3dGrid class containing the spatial and frequency/wavelength grid
    
    ppar : dictionary
           A dictionary containing all parameters of the model
    
    Returns
    -------
    Returns an ndarray with the turbulent velocity in cm/s
    """

    vturb = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) + ppar['gasspec_vturb']
    return vturb


def getVelocity(grid=None, ppar=None):
    """Calculates the velocity field in a protoplanetary disk.
    
    Parameters
    ----------
    grid : radmc3dGrid
           An instance of the radmc3dGrid class containing the spatial and frequency/wavelength grid
    
    ppar : dictionary
           A dictionary containing all parameters of the model
    
    Returns
    -------
    Returns the gas velocity in cm/s
    """

    nr = grid.nx
    nphi = grid.nz
    nz = grid.ny
    rcyl = grid.x

    vel = np.zeros([nr, nz, nphi, 3], dtype=np.float64)
    vkep = np.sqrt(gg * ppar['mstar'][0]/rcyl)
    for iz in range(nz):
        for ip in range(nphi):
            vel[:, iz, ip, 2] = vkep

    return vel
