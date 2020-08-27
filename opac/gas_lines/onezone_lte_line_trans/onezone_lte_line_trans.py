import numpy as np
import math
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from radmc3dMolecule import *
from radmc3dDustOpac import *
#
# Some useful physical constants in CGS units
#
mp  = 1.6726e-24         # Mass of proton          [g]
kk  = 1.3807e-16         # Bolzmann's constant     [erg/K]
hh  = 6.6262e-27         # Planck's constant       [erg.s]
cc  = 2.9979245800000e10 # Light speed             [cm/s]
ss  = 5.6703e-5          # Stefan-Boltzmann const  [erg/cm^2/K^4/s]

class onezone_lte_line_trans(object):
    """
    A very simple one-zone radiative transfer model for line transfer.
    """

    def __init__(self,dust=None,molecule=None,dvkmsmax=None,nv=101,
                 colgas=1e-4,temp=30.,vturb=1e4,molabund=1e-4,iline=0,
                 dusttogas=0.01):
        #
        # Model parameters (to be modified by the user)
        # 
        self.colgas     = colgas              # Column density of H2+He+XXX gas in gram/cm^2 (default value)
        self.mugas      = 2.3                 # Mean molecular weight of the gas in units of mp
        self.temp       = temp                # Temperature in K (default value)
        self.dusttogas  = dusttogas           # Dust-to-gas mass ratio (default value)
        self.vturb      = vturb               # Microturbulence in cm/s (default value)
        self.molabund   = molabund            # Abundance of molecule: molecules per gas particle (default value)
        self.iline      = iline               # The line to be modeled (default value)
        self.bgint      = 0.0                 # Background intensity for line transfer
        self.nv         = nv                  # Nr of velocity bins
        #
        # Read the dust opacity
        #
        if(dust is not None):
            self.readdust(dust)
        #
        # Read the molecule data
        #
        if(molecule is not None):
            self.readmol(molecule)
        #
        # Set up the nu grid
        #
        if(dvkmsmax is not None):
            self.dvkmsmax = dvkmsmax
            self.setup_dvgrid(nv,dvkmsmax)
        
    def readmol(self,molname):
        self.molecule = radmc3dMolecule()
        self.molecule.read(mol=molname)
        self.molname=molname

    def readdust(self,dustname):
        self.dust = radmc3dDustOpac()
        self.dust.readOpac(ext=dustname)
        self.dustname=dustname

    def setup_dvgrid(self,nv=101,dvkmsmax=1.):
        self.dvkms    = dvkmsmax * np.linspace(-1.,1.,nv)
        nu0           = self.molecule.freq[self.iline]
        self.nu       = nu0*(1.0+self.dvkms*1e5/cc)
        
    def compute_lte(self):
        ndum          = self.molecule.wgt*np.exp(-self.molecule.energy/(kk*self.temp))
        self.zpart    = ndum.sum()            # The partition sum
        self.pop      = ndum/self.zpart       # The fractional level populations

    def compute_numcol(self):
        self.gasnumcol = self.colgas/(self.mugas*mp)
        self.molnumcol = self.molabund*self.gasnumcol
        
    def compute_lineprofile(self):
        if(hasattr(self,'dvkms')==False):
            raise NameError('You must first set up a grid of velocities by calling setup_dvgrid().')
        atherm   = math.sqrt(2*kk*self.temp/(self.molecule.molweight*mp))
        aturb    = self.vturb
        atot     = math.sqrt(atherm**2+aturb**2)
        nu0      = self.molecule.freq[self.iline]
        dv       = self.dvkms*1e5
        self.phi = (cc/(atot*nu0*math.sqrt(math.pi)))*np.exp(-(dv/atot)**2)

    def compute_einstein_b(self):
        iup               = self.molecule.iup[self.iline]-1
        idown             = self.molecule.ilow[self.iline]-1
        self.molecule.bud = self.molecule.aud * cc**2 / (2*hh*(self.molecule.freq**3))
        self.molecule.bdu = self.molecule.bud * self.molecule.wgt[iup] / self.molecule.wgt[idown]
        
    def compute_bplanck(self):
        if(hasattr(self,'dvkms')==False):
            raise NameError('You must first set up a grid of velocities by calling setup_dvgrid().')
        self.bpl = np.zeros(len(self.nu))
        x        = hh*self.nu/(kk*self.temp)
        for i in range(len(self.nu)):
            if x[i] < 1.e2:
                if x[i] > 1.e-3:
                    self.bpl[i] = 2.0*hh*(self.nu[i]**3)/(cc*cc*(math.exp(x[i])-1.e0))
                else:
                    self.bpl[i] = 2.0*(self.nu[i]**2)*kk*self.temp/(cc*cc)

    def compute_taudust(self):
        if(hasattr(self,'nu')==False):
            raise NameError('You must first set up a grid of velocities by calling setup_dvgrid().')
        if(self.dusttogas>0.0):
            if(hasattr(self,'dust')==False):
                raise NameError('You must first read in the dust opacity by calling readdust().')
            f              = interp1d(self.dust.freq[0],self.dust.kabs[0])   # Only 1 dust species allowed
            self.dust_kabs = f(self.nu)
            self.taudust   = self.dust_kabs*self.dusttogas*self.colgas
        else:
            self.dust_kabs = 0.0
            self.taudust   = 0.0

    def compute_tauline(self):
        if(hasattr(self,'nu')==False):
            raise NameError('You must first set up a grid of velocities by calling setup_dvgrid().')
        if(hasattr(self,'pop')==False):
            raise NameError('You must first call compute_lte().')
        if(hasattr(self.molecule,'bud')==False):
            raise NameError('You must first call compute_einstein_b().')
        if(hasattr(self,'molnumcol')==False):
            raise NameError('You must first call compute_numcol().')
        if(hasattr(self,'phi')==False):
            raise NameError('You must first call compute_lineprofile().')
        iup          = self.molecule.iup[self.iline]-1
        idown        = self.molecule.ilow[self.iline]-1
        bud          = self.molecule.bud[self.iline]
        bdu          = self.molecule.bdu[self.iline]
        nu0          = self.molecule.freq[self.iline]
        nup          = self.pop[iup]*self.molnumcol
        ndown        = self.pop[idown]*self.molnumcol
        self.tauline = (hh*nu0/(4*math.pi))*(ndown*bdu-nup*bud)*self.phi

    def compute_intensity(self):
        if(hasattr(self,'tauline')==False):
            raise NameError('You must first call compute_tauline().')
        if(hasattr(self,'taudust')==False):
            raise NameError('You must first call compute_taudust().')
        if(hasattr(self,'bpl')==False):
            raise NameError('You must first call compute_bplanck().')
        self.tau       = self.taudust + self.tauline
        xp             = np.exp(-self.tau)
        xp1            = 1.0-xp
        idx            = np.where(self.tau<1e-6)
        xp1[idx]       = self.tau[idx]
        self.intensity = self.bgint*xp + self.bpl*xp1

    def doall(self):
        if(hasattr(self,'molecule')==False):
            raise NameError('You must first call readmol() with the molecule name as argument.')
        if(hasattr(self,'dust')==False and self.dusttogas!=0.0):
            raise NameError('You must first call readdust() with the dust name as argument.')
        if(hasattr(self,'nu')==False):
            raise NameError('You must first set up a grid of velocities by calling setup_dvgrid().')
        if(hasattr(self,'dvkmsmax')==False):
            raise NameError('You must first set dvkmsmax.')
        self.compute_einstein_b()
        self.compute_lte()
        self.compute_numcol()
        self.setup_dvgrid(self.nv,self.dvkmsmax)
        self.compute_lineprofile()
        self.compute_bplanck()
        self.compute_taudust()
        self.compute_tauline()
        self.compute_intensity()
        
        
def test():
    """
    This is a simple test routine to test the onezone RT model. Make sure that the following datafiles are present:
      molecule_co.inp
      dustkappa_silicate.inp
    otherwise the model will not work. The model parameters are taken to be the default values.
    This routine simply returns the object q which contains the results.
    """
    q = onezone_lte_line_trans(dust='silicate',molecule='co',dvkmsmax=1.0,nv=301,
                               colgas=1e-3,temp=30.,vturb=1e4,iline=1)
    q.doall()
    plt.plot(q.dvkms,q.intensity,label='Result')
    plt.plot(q.dvkms,q.bpl,label='Planck')
    plt.legend()
    plt.yscale('log')
    plt.xlabel(r'$\Delta v [\mathrm{km/s}]$')
    plt.ylabel(r'$I_\nu [\mathrm{erg}\,\mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{Hz}^{-1}\mathrm{ster}^{-1}]$')
    plt.show()
    return q
