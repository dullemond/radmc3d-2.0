import numpy as np
import matplotlib.pyplot as plt
import os

class frame(object):
    def __init__(self,iframe,rho=True,rhodust=False,e=False,vphi=False,vr=False,dir='.'):
        phi      = np.fromfile(dir+'/domain_x.dat',sep='\n')
        r        = np.fromfile(dir+'/domain_y.dat',sep='\n')
        self.phi = 0.5 * ( phi[:-1] + phi[1:] )
        self.r   = (r[3:-4]*r[4:-3])**0.5
        nphi     = len(self.phi)
        nr       = len(self.r)
        sframe = '{}'.format(iframe)
        if rho:  self.rho   = np.fromfile(dir+'/gasdens'+sframe+'.dat').reshape(nr,nphi)
        if e:    self.e     = np.fromfile(dir+'/gasenergy'+sframe+'.dat').reshape(nr,nphi)
        if vphi: self.vphi  = np.fromfile(dir+'/gasvx'+sframe+'.dat').reshape(nr,nphi)
        if vr:   self.vr    = np.fromfile(dir+'/gasvy'+sframe+'.dat').reshape(nr,nphi)
        if rhodust:
            self.rhodust = []
            idust = 1
            idusts= '{}'.format(idust)
            file  = dir+'/dust'+idusts+'dens'+sframe+'.dat'
            while os.path.isfile(file):
                rhodust = np.fromfile(file).reshape(nr,nphi)
                self.rhodust.append(rhodust)
                idust += 1
                idusts= '{}'.format(idust)
                file  = dir+'/dust'+idusts+'dens'+sframe+'.dat'

    def show(self,q=None,log=False,min=None,max=None):
        if q is None:
            q = self.rho
        if log:
            qq = np.log10(q)
        else:
            qq = q.copy()
        if min is not None:
            qq[qq<min]=min
        if max is not None:
            qq[qq>min]=max
        plt.imshow(qq,origin='lower',extent=[0,360,self.r.min(),self.r.max()],aspect='auto')
        plt.xlabel(r'$\phi\; [\mathrm{deg}]$')
        plt.ylabel(r'$r\; [\mathrm{code}\;\mathrm{units}]$')
        plt.show()

    def compute_azimuthal_average(self):
        self.rho_av = self.rho.sum(axis=1)/len(self.phi)

    def convert_to_cgs(self,mstar,r0):
        GG  = 6.67408e-08    # Gravitational constant  [cm^3/g/s^2]
        self.r_cgs = self.r * r0
        if hasattr(self,'rho'):
            self.sigma_gas_cgs  = self.rho * mstar / (r0*r0)  # Surface density in g/cm^2
        if hasattr(self,'rhodust'):
            self.sigma_dust_cgs = []
            for rd in self.rhodust:
                sigdust = rd * mstar / (r0*r0)                # Surface density in g/cm^2
                self.sigma_dust_cgs.append(sigdust)
        if hasattr(self,'vphi'):
            self.vphi_cgs = self.vphi * np.sqrt(GG*mstar/r0)  # Phi-velocity in cm/s
        if hasattr(self,'vr'):
            self.vr_cgs = self.vr * np.sqrt(GG*mstar/r0)      # r-velocity in cm/s
        if hasattr(self,'e'):
            self.gamma  = 1.66666667    # Default value of FARGO3D
            self.cs     = np.sqrt((self.gamma-1.)*self.e/self.rho)
            self.cs_cgs = self.cs * np.sqrt(GG*mstar/r0)      # Isothermal sound speed in cm/s
            
