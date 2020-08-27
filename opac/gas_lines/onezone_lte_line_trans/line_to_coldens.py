from onezone_lte_line_trans import *
from radmc3dMolecule import *
import numpy as np
import math
import matplotlib.pyplot as plt
#import scipy.optimize.brentq
#
# Some useful physical constants in CGS units
#
mp  = 1.6726e-24         # Mass of proton          [g]
kk  = 1.3807e-16         # Bolzmann's constant     [erg/K]
hh  = 6.6262e-27         # Planck's constant       [erg.s]
cc  = 2.9979245800000e10 # Light speed             [cm/s]
ss  = 5.6703e-5          # Stefan-Boltzmann const  [erg/cm^2/K^4/s]


def intensity_to_kelvin(lammic,intensity):
    nu   = 1e4*cc/lammic
    temp = intensity * cc**2 / ( 2 * nu**2 * kk )
    return temp

def kelvin_to_intensity(lammic,temp):
    nu   = 1e4*cc/lammic
    int  = temp * ( 2 * nu**2 * kk ) / cc**2
    return int

def bplfunc(nu,t):
    cc  = 2.9979245800000e10
    hh  = 6.6262e-27
    kk  = 1.3807e-16
    bpl = 0.e0
    if t > 0.0:
        x = hh*nu/(kk*t)
        if x < 1.e2:
            if x > 1.e-3:
                bpl = 2.0*hh*nu*nu*nu/(cc*cc*(math.exp(x)-1.e0))
            else:
                bpl = 2.0*nu*nu*kk*t/(cc*cc)
    return bpl

def bplanck(nu,t):
    if type(nu)==np.ndarray:
        n = nu.size
        bnu = np.zeros(n)
        for i in range(n):
            bnu[i] = bplfunc(nu[i],t)
    else:
        bnu = bplfunc(nu,t)
    return bnu

class line_to_coldens(object):
    """
    Computes the column density of a molecule using LTE line transfer, given the
    intensity at line center and the temperature of the gas, as well as the 
    turbulent line broadning in km/s. Note that this routine does not include any
    background dust continuum.

    ARGUMENTS:
      molname      The name of the molecule. E.g. 'co'. Uses the Leiden LAMDA data format,
                   i.e. with molname='co' it will read the file molecule_co.inp .
      iline        The line to be modeled (iline=0 is the first, which for CO is J 1-0, 
                   iline=1 second, which for CO is J 2-1, etc).
      temp         Gas temperature in K
      vturbkms     Turbulent velocity in km/s
      intensity    The observed intensity in erg/s/cm^2/Hz/ster. Note that you can convert
                   from other units with the subroutines given above.

    RETURNS:
      self.molcoldens     The column density of the molecule in question in units of cm^{-2}.
      self.tau            The line center optical depth
    """

    def __init__(self,molname,iline,temp,vturbkms,intensity):
        colgas = 2.3*mp
        q = onezone_lte_line_trans(molecule=molname,
                                   dvkmsmax=30.0,nv=3,dusttogas=0.0,
                                   colgas=colgas,molabund=1.0,temp=temp,
                                   vturb=1e5*vturbkms,iline=iline)
        q.doall()
        self.molname   = molname
        self.iline     = iline
        self.temp      = temp
        self.vturbkms  = vturbkms
        self.intensity = intensity
        self.nu        = q.molecule.freq[iline]
        self.lammic    = 1e4*cc/self.nu
        self.bpl       = bplanck(self.nu,temp)
        self.bplrj     = 2*self.nu**2*kk*temp/cc**2
        if(intensity>self.bpl):
            raise ValueError('Intensity > B_nu(T) for the given temperature. No solution possible.')
        if(math.fabs((self.bpl-intensity)/(self.bpl+intensity))<1e-6):
            raise ValueError('Intensity too close to B_nu(T). Degenerate solutions (tau>12).')
        frac = intensity/self.bpl
        self.tau        = -math.log(1.0-frac)
        self.molcoldens = self.tau / q.tau[1]
        self.rt         = q

def ltctest():
    intens_temp = 6.    # 6 Kelvin intensity, assuming Rayleigh-Jeans
    temp        = 30.   # Assumed gas temperature
    iline       = 1     # CO J 2-1
    vturbkms    = 0.1   # Assumed turbulent line broadning
    molname     = 'co'  # Name of the molecule
    molabun     = 1e-4  # Abundance of the molecule compared to number density of gas molecules
    #
    # Convert intensity from Kelvin to erg/cm^2/s/Hz/ster
    #
    molecule    = readMol(mol=molname)
    lammic      = 1e4*cc/molecule.freq[iline]
    intensity   = kelvin_to_intensity(lammic,intens_temp)
    #
    # Computing the column density
    #
    q           = line_to_coldens(molname,iline,temp,vturbkms,intensity)
    q.gascoldens= q.molcoldens/molabun
    q.gassigma  = q.gascoldens*2.3*mp
    print("CO Column density  = {0:13.6e} cm^(-2)".format(q.molcoldens))
    print("CO Line center tau = {0:13.6e}".format(q.tau))
    print("B_nu(T)            = {0:13.6e} erg/cm^2/s/Hz/ster".format(q.bpl))
    print("B_nu(T)/RayJeans   = {0:13.6e}".format(q.bpl/q.bplrj))
    print("Gas Column density = {0:13.6e} cm^(-2)".format(q.gascoldens))
    print("Gas Column density = {0:13.6e} gram/cm^(-2)".format(q.gassigma))
    #
    # For plotting, call the onezone model again but now with many wavelength points
    #
    fullrt      = onezone_lte_line_trans(molecule=molname,
                                   dvkmsmax=1.0,nv=301,dusttogas=0.0,
                                   colgas=q.gassigma,molabund=molabun,temp=temp,
                                   vturb=1e5*vturbkms,iline=iline)
    fullrt.doall()
    q.fullrt    = fullrt
    fig         = plt.figure()
    plt.plot(fullrt.dvkms,fullrt.intensity,label='Result')
    plt.plot(fullrt.dvkms,fullrt.bpl,label='Planck')
    plt.legend()
    #plt.yscale('log')
    plt.xlabel(r'$\Delta v [\mathrm{km/s}]$')
    plt.ylabel(r'$I_\nu [\mathrm{erg}\,\mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{Hz}^{-1}\mathrm{ster}^{-1}]$')
    plt.show()
    return q
