#
# This is actually part of radmc3dPy/analyze.py, but here extracted as stand-alone.
# NOTE: If you download a file from the LAMDA database, you must rename it
#       to molecule_xxx.inp. Example: molecule_co.inp for co.dat.
# 2017.08.19
#
import numpy as np

class radmc3dMolecule(object):
   """
   RADMC-3D molecule class
   Based on the Leiden LAMDA database, but is in principle generic

   NOTE: For now only the levels and lines are included, not the 
         collision rates. 

   Attributes
   ----------
   name            : str
                    The name as listed in the molecule file

   molweight       : float
                    Molecular weight in units of proton mass
   
   nlev            : int
                    Nr of energy levels
   
   nlin            : int
                    Nr of lines
   
   energycminv     : float
                    Energy[ilev] of level ilev in 1/cm
   
   energy          : float
                    Energy[ilev] of level ilev in erg
   
   wgt             : float
                    Statistical weight[ilev] of level ilev
   
   jrot            : float
                    Quantum rotational J[ilev] of level ilev
   
   iup             : int 
                    ilev of upper level of line ilin (starting with 0)
   
   ilow            : int
                    ilev of lower level of line ilin (starting with 0)
   
   aud             : float
                    Einstein A up low of line ilin in 1/second
   
   freq            : float
                    Frequency of line ilin in Hz
   
   lam             : float
                    Wavelength of line ilin in micron

   """

   def __init__(self):
       self.name        = ""
       self.molweight   = 0.0
       self.nlev        = 0
       self.nlin        = 0
       self.energycminv = 0.0
       self.energy      = 0.0
       self.wgt         = 0.0
       self.jrot        = 0.0
       self.iup         = 0
       self.ilow        = 0
       self.aud         = 0.0
       self.freq        = 0.0
       self.lam         = 0.0

   # --------------------------------------------------------------------------------------------------
   def read(self,mol='',fname=''):
       """Read the molecule_<mol>.inp file

       The file format is the format of the Leiden LAMDA molecular database

       Parameters
       ----------
       mol             : str
                        molecule name (e.g. 'co') if the file name is in the form of 'molecule_<mol>.inp'
       
       fname           : str
                        full file name
       """

       if fname != '':
           try:
               f = open(fname, 'r')
           except Exception as e:
               print(e)
               return False

       else:
           fname = 'molecule_'+mol+'.inp'
           try:
               f = open(fname, 'r')
           except Exception as e:
               print(e)
               return False
        
       
       print('Reading '+fname+'...')
       #with open(fname,'r') as f:
       dum             = f.readline()
       dum             = f.readline().split()
       self.name       = dum[0]
       dum             = f.readline()
       self.molweight  = float(f.readline())
       dum             = f.readline()
       self.nlev       = int(f.readline())
       dum             = f.readline()
       self.energycminv= np.zeros(self.nlev)
       self.energy     = np.zeros(self.nlev)
       self.wgt        = np.zeros(self.nlev)
       self.jrot       = np.zeros(self.nlev)
       for i in range(self.nlev):
           dum                 = f.readline().split()
           self.energycminv[i] = float(dum[1])
           self.energy[i]      = float(dum[1])*1.9864847851996e-16  # const=h*c
           self.wgt[i]         = float(dum[2])
           self.jrot[i]        = float(dum[3])
       dum             = f.readline()
       self.nlin       = int(f.readline())
       dum             = f.readline()
       self.iup        = np.zeros(self.nlin,dtype=np.int)
       self.ilow       = np.zeros(self.nlin,dtype=np.int)
       self.aud        = np.zeros(self.nlin)
       self.freq       = np.zeros(self.nlin)
       self.lam        = np.zeros(self.nlin)
       for i in range(self.nlin):
           dum            = f.readline().split()
           self.iup[i]    = int(dum[1])   # Use as index: [iup-1]
           self.ilow[i]   = int(dum[2])   # Use as index: [ilow-1]
           self.aud[i]    = float(dum[3])
           self.freq[i]   = float(dum[4])*1e9
           self.lam[i]    = 2.9979245800000e+14/self.freq[i]
        
       f.close()

       return True
# --------------------------------------------------------------------------------------------------
def readMol(mol='', fname=''):
    """ Wrapper around the radmc3dMolecule.read() method

       Parameters
       ----------
       mol             : str
                        molecule name (e.g. 'co') if the file name is in the form of 'molecule_<mol>.inp'

       fname           : str
                        full file name
    """

    m = radmc3dMolecule()
    if m.read(mol=mol, fname=fname) == True:
        return m
    else:
        return

