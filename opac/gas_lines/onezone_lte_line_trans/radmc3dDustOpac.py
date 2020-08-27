#
# READING ROUTINES FOR RADMC-3D-STYLE OPACITY FILES
#
# This is actually part of radmc3dPy/analyze.py, but here extracted as stand-alone, and stripped to
# the basics. It is meant to read RADMC-3D dust opacity files.
#
# 2017.08.19
#
import numpy as np
cc = 29979245800.

class radmc3dDustOpac(object):
    """
    Class to handle dust opacities.


    Attributes
    ----------

    wav     : list
                Each element of the list contains an ndarray with the wavelength grid
    
    freq    : list
                Each element of the list contains an ndarray with the frequency grid
    
    nwav    : list
                Each element of the list contains an integer with the number of wavelengths
    
    kabs    : list
                Each element of the list contains an ndarray with the absorption coefficient per unit mass 
    
    ksca    : list
                Each element of the list contains an ndarray with the scattering coefficient per unit mass
    
    phase_g : list
                Each element of the list contains an ndarray with the hase function
    
    ext     : list
                Each element of the list contains a string wht the file name extension of the duskappa_ext.Kappa file
    
    therm   : list
                Each element of the list contains a bool, if it is set to False the dust grains are quantum-heated (default: True)
    
    idust   : lisintt
                Each element of the list contains an integer with the index of the dust species in the dust density distribution array

    scatmat : list
                Each element is a boolean indicating whether the dust opacity table includes (True) the full scattering matrix or not (False)

    nang    : list
                Each element is a string, containing the number of scattering angles in the scattering matrix if its given

    scatang : list
                Each element is a numpy ndarray containing the scattering angles in the scattering matrix if its given

    z11     : list
                Each element is a numpy ndarray containing the (1,1) element of the scattering angles in the scattering matrix if its given

    z12     : list
                Each element is a numpy ndarray containing the (1,2) element of the scattering angles in the scattering matrix if its given

    z22     : list
                Each element is a numpy ndarray containing the (2,2) element of the scattering angles in the scattering matrix if its given
    
    z33     : list
                Each element is a numpy ndarray containing the (3,3) element of the scattering angles in the scattering matrix if its given
    
    z34     : list
                Each element is a numpy ndarray containing the (3,4) element of the scattering angles in the scattering matrix if its given
    
    z44     : list
                Each element is a numpy ndarray containing the (4,4) element of the scattering angles in the scattering matrix if its given

    """
# --------------------------------------------------------------------------------------------------
    def __init__(self):

        self.wav      = []
        self.freq     = []
        self.nwav     = []
        self.nfreq    = []
        self.kabs     = []
        self.ksca     = []
        self.phase_g  = []
        self.ext      = []
        self.idust    = []
        self.therm    = []
        self.scatmat  = []
        self.z11      = [] 
        self.z12      = [] 
        self.z22      = [] 
        self.z33      = [] 
        self.z34      = [] 
        self.z44      = [] 
        self.scatang  = []
        self.nang     = []

# --------------------------------------------------------------------------------------------------
    def readOpac(self, ext=[''], idust=None, scatmat=None, old=False):
        """Reads the dust opacity files.

        Parameters
        ----------
        
        ext  : list
                File name extension (file names should look like 'dustkappa_ext.inp')
        
        idust: list
                Indices of the dust species in the master opacity file (dustopac.inp') - starts at 0 
        
        scatmat: list
                If specified, its elements should be booleans indicating whether the opacity file 
                contains also the full scattering matrix (True) or only dust opacities (False)
        
        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """
       
        # Check the input keywords and if single strings are given convert them to lists
        # This assumes, though, that there is a single dust opacity file or dust species, though!!
        if (type(ext).__name__=='str'):  ext = [ext]
        if scatmat!=None:
            if (type(scatmat).__name__=='str'):  scatmat = [scatmat]
        else:
            # If the scatmat keyword is not given (i.e. if it is None) then assume that 
            # it is False for all dust species
            scatmat = []
            if idust!=None:
                for i in range(len(idust)):
                    scatmat.append(False)
            else:
                for i in range(len(ext)):
                    scatmat.append(False)
        
        if idust!=None:
            if (type(idust).__name__=='int'):  idust = [idust]

        if (len(ext)==1)&(ext[0]!=''):
            if idust!=None:
                print('ERROR')
                print('Either idust or ext should be specified, but not both')
                print(idust)
                print(ext)
                return [-1]
        
        # Read the master dust opacity file to get the dust indices and dustkappa file name extensions
        #mopac = self.readMasterOpac()
        therm = []
        for i in range(len(ext)):
            therm.append(True)
        mopac = {'ext':ext, 'therm':therm, 'scatmat':scatmat}

        # Find the file name extensions in the master opacity file if idust is specified instead of ext
        if idust:
            ext = []
            for ispec in idust:
                if (ispec+1)>len(mopac['ext']):    
                    print('ERROR')
                    print('No dust species found at index ', ispec)
                    return [-1]
                else:
                    ext.append(mopac['ext'][ispec])

        # If only the extension is specified look for the master opacity file and find the index of this dust species
        #  or set the index to -1 if no such dust species is present in the master opacity file
        else:
            idust = []
            for iext in ext:
                try:
                    dum2 = mopac['ext'].index(iext)
                except:
                    dum2 = -1
                idust.append(dum2)
        
        # Now read all dust opacities
        for i in range(len(ext)):
            if scatmat[i]:
                try:
                    rfile = open('dustkapscatmat_'+ext[i]+'.inp', 'r')
                except:
                    print('ERROR')
                    print(' No dustkapscatmat_'+ext[i]+'.inp file was found')
                    return -1
                
                print('Reading dustkapscatmat_'+ext[i]+'.inp ....')

                self.ext.append(ext[i])
                
                # Read the header/comment field
                dum = rfile.readline()
                while dum.strip()[0]=='#':
                    dum = rfile.readline()


                #for j in range(6):
                    #dum = rfile.readline()

                # Read the file format
                iformat = int(dum)
                #iformat = int(rfile.readline())
                if iformat!=1:
                    print('ERROR')
                    print('Format number of the file dustkapscatmat_'+ext[i]+'.inp (iformat='+("%d"%iformat)+') is unkown')
                    return [-1]

                # Read the number of wavelengths in the file
                dum = int(rfile.readline())
                self.nwav.append(dum)
                self.nfreq.append(dum)
                self.idust.append(idust[i])
                idu = len(self.nwav)-1
                
                # Read the scattering angular grid
                self.nang.append(int(rfile.readline()))
                wav     = np.zeros(self.nwav[idu], dtype = np.float64)
                kabs    = np.zeros(self.nwav[idu], dtype = np.float64)
                ksca    = np.zeros(self.nwav[idu], dtype = np.float64)
                phase_g = np.zeros(self.nwav[idu], dtype = np.float64)
                scatang = np.zeros(self.nang[idu], dtype=np.float64)
                z11     = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64) 
                z12     = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64) 
                z22     = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64) 
                z33     = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64) 
                z34     = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64) 
                z44     = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64) 
            
                print('Reading the opacities..')
                dum = rfile.readline()
                for ilam in range(self.nwav[idu]):
                    dum      = rfile.readline().split()
                    wav[ilam]  = float(dum[0])
                    kabs[ilam] = float(dum[1])
                    ksca[ilam] = float(dum[2])
                    phase_g[ilam] = float(dum[3])

                print('Reading the angular grid..')
                dum = rfile.readline()
                for iang in range(self.nang[idu]):
                    dum        = rfile.readline()
                    scatang[iang] = float(dum)

                print('Reading the scattering matrix..')
                for ilam in range(self.nwav[idu]):
                    dum = rfile.readline()
                    for iang in range(self.nang[idu]):
                        dum        = rfile.readline().split()
                        z11[ilam,iang] = float(dum[0])
                        z12[ilam,iang] = float(dum[1])
                        z22[ilam,iang] = float(dum[2])
                        z33[ilam,iang] = float(dum[3])
                        z34[ilam,iang] = float(dum[4])
                        z44[ilam,iang] = float(dum[5])
                
                self.wav.append(wav)
                self.freq.append(cc/wav*1e4)
                self.kabs.append(kabs)
                self.ksca.append(ksca)
                self.phase_g.append(phase_g)
                self.scatang.append(scatang)
                self.z11.append(z11)
                self.z12.append(z12)
                self.z22.append(z22)
                self.z33.append(z33)
                self.z34.append(z34)
                self.z44.append(z44)
               
                rfile.close()
            else:
                if not old:
                    try:
                        rfile = open('dustkappa_'+ext[i]+'.inp', 'r')
                    except:
                        print('ERROR')
                        print(' No dustkappa_'+ext[i]+'.inp file was found')
                        return -1

                    self.ext.append(ext[i])

                    # Read the file format
                    iformat = int(rfile.readline())
                    if (iformat<1)|(iformat>3):
                        print('ERROR')
                        print('Unknown file format in the dust opacity file')
                        rfile.close()
                        return -1


                    # Read the number of wavelengths in the file
                    dum = rfile.readline()
                    self.nwav.append(int(dum))
                    self.nfreq.append(int(dum))
                    self.idust.append(idust[i])
                    idu = len(self.nwav)-1

                    # If only the absorption coefficients are specified
                    if iformat==1:
                        wav = np.zeros(self.nwav[idu], dtype=np.float64)
                        kabs = np.zeros(self.nwav[idu], dtype=np.float64)
                        for ilam in range(self.nwav[idu]):
                            dum = rfile.readline().split()
                            wav[ilam] = float(dum[0])
                            kabs[ilam] = float(dum[1])
                        self.wav.append(wav)
                        self.freq.append(cc/wav*1e4)
                        self.kabs.append(kabs)
                        self.ksca.append([-1])
                        self.phase_g.append([-1])
                    # If the absorption and scattering coefficients are specified
                    elif iformat==2:
                        wav = np.zeros(self.nwav[idu], dtype=np.float64)
                        kabs = np.zeros(self.nwav[idu], dtype=np.float64)
                        ksca = np.zeros(self.nwav[idu], dtype=np.float64)
                        for ilam in range(self.nwav[idu]):
                            dum = rfile.readline().split()
                            wav[ilam] = float(dum[0])
                            kabs[ilam] = float(dum[1])
                            ksca[ilam] = float(dum[2]) 
                        self.wav.append(wav)
                        self.freq.append(cc/wav*1e4)
                        self.kabs.append(kabs)
                        self.ksca.append(ksca)
                        self.phase_g.append([-1])
                    
                    # If the absorption and scattering coefficients and also the scattering phase function are specified
                    elif iformat==3:
                        wav = np.zeros(self.nwav[idu], dtype=np.float64)
                        kabs = np.zeros(self.nwav[idu], dtype=np.float64)
                        ksca = np.zeros(self.nwav[idu], dtype=np.float64)
                        phase_g = np.zeros(self.nwav[idu], dtype=np.float64)
                        for ilam in range(self.nwav[idu]):
                            dum = rfile.readline().split()
                            wav[ilam] = float(dum[0])
                            kabs[ilam] = float(dum[1])
                            ksca[ilam] = float(dum[2])
                            phase_g[ilam] = float(dum[3]) 
                        
                        self.wav.append(wav)
                        self.freq.append(cc/wav*1e4)
                        self.kabs.append(kabs)
                        self.ksca.append(ksca)
                        self.phase_g.append(phase_g)
               
                    rfile.close()
                else:
                    try:
                        rfile = open('dustopac_'+ext[i]+'.inp', 'r')
                    except:
                        print('ERROR')
                        print(' No dustopac_'+ext[i]+'.inp file was found')
                        return -1
                 
                    freq = np.fromfile('frequency.inp', count=-1, sep="\n", dtype=float)
                    nfreq = int(freq[0])
                    freq = freq[1:]

                    self.ext.append(ext[i])
                    dum   = rfile.readline().split()
                    if int(dum[0])!=nfreq:
                        print('ERROR')
                        print('dustopac_'+ext[i]+'.inp contains a different number of frequencies than frequency.inp')
                        return

                    wav     = cc/freq*1e4
                    kabs    = np.zeros(nfreq, dtype = float)
                    ksca    = np.zeros(nfreq, dtype = float)

                    dum     = rfile.readline()
                    for ilam in range(nfreq):
                        kabs[ilam] = float(rfile.readline())
                    dum     = rfile.readline()
                    for ilam in range(nfreq):
                        ksca[ilam] = float(rfile.readline())
                    
                    rfile.close()

                    self.wav.append(wav[::-1])
                    self.freq.append(freq[::-1])
                    self.kabs.append(kabs[::-1])
                    self.ksca.append(ksca[::-1])
                    self.phase_g.append([-1])


                    
        return 0 



    
def readOpac(ext=[''], idust=None, scatmat=None, old=False):
    """Reads the dust opacity files.
    This function is an interface to radmc3dDustOpac.readOpac()

    Parameters
    ----------
    ext   : list
            Each element of the list is be a string, the file name extension (file names should look like 'dustkappa_ext.inp')
    
    idust : list
            Each element of the list is an integer, the index of the dust species in the master opacity file (dustopac.inp')

    scatmat: list
            If specified, its elements should be booleans indicating whether the opacity file 
            contains also the full scattering matrix (True) or only dust opacities (False)
        
    old   : bool, optional
            If set to True the file format of the previous, 2D version of radmc will be used
    
    Returns
    -------
        Returns an instance of the radmc3dDustOpac class 
    """


    res = radmc3dDustOpac()
    res.readOpac(ext=ext, idust=idust, scatmat=scatmat, old=old)
    
    return res
    
