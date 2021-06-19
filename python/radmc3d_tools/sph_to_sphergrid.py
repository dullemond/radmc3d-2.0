import numpy as np
from natconst import *
from numpy.random import default_rng
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
from tqdm import tqdm   # Progress bar

class sph_to_sphergrid(object):
    def __init__(self,rthphi=None,xyz=None,masses=1,rin=None,rout=None,nr=100,ntheta=100,nphi=100,
                 hrup=0.15,hrkernel=0.03):
        """
        Simple conversion of SPH model of a circumstellar disk to spherical coordinates 
        for RADMC-3D. The method simply places particles on the grid by the binning (counting)
        method, where initially each particle is a point particle. Then we use an FFT convolution
        method to smear out each grid cell with a smoothing kernel. 

        The basic assumptions of this method are:
          (1) We assume that all SPH particles have a smoothing radius proportional to
              distance from the center of the coordinate system (the location of the star).
              In other words: you specify a dimensionless smoothing radius hr such that the
              real smoothing radius is h = r * hr. If your SPH model has constant or variable
              SPH kernel size: No problem; the only thing is that the mapping of the SPH 
              model onto the grid will then not be exactly according to the SPH kernel size.
          (2) The grid in r is logarithmic, in the sense that r[i+1]/r[i] is the same for
              all i. The grid in theta is uniform between pi/2-hrup and pi/2+hrup (where
              pi/2 is the equatorial plane). The phi-grid is uniform between 0 and 2*pi, and
              periodic. This grid is generated internally and you can specify the parameters
              with the keywords explained below.
        Assumptions (1) and (2) are unavoidable because (for speed) we use an FFT-based 
        convolution method, meaning that the kernel must be a global kernel, and cannot
        depend on location. With (2) the grid is self-similar (at least, if you do not 
        worry about the theta-dependence of the cell volume, which is only minor near
        the equatorial plane). 

        KEYWORD ARGUMENTS:

          rthphi         A numpy array of dimension rthphi[nsph,3] where nsph is the 
                         number of SPH particles and rthphi[:,0] = r, rthphi[:,1] = theta,
                         rthphi[:,2] = phi. Note: rthphi and xyz are mutually exclusive. Use
                         either one or the other.
          xyz            A numpy array of dimension xyz[nsph,3] where nsph is the 
                         number of SPH particles and xyz[:,0] = x, xyz[:,1] = y,
                         xyz[:,2] = z. Note: rthphi and xyz are mutually exclusive. Use
                         either one or the other.
          rin            The inner radius of the r grid
          rout           The outer radius of the r grid
          nr             The number of grid cells in r
          ntheta         The number of grid cells in th
          nphi           The number of grid cells in phi
          hrup           The vertical extent (in terms of pi/2-theta) of the theta grid.
          masses         The mass of each SPH particle. If a scalar: all masses are the
                         same. If m[nsph] with nsph the number of SPH particles, then 
                         each SPH particle has its own mass.
          hrkernel       The dimensionless radius of the smoothing kernel: dimensionless
                         as in units of ln(r) or theta or phi. If it is a scalar, then the
                         kernel is spherical. If it is an array of 3 values, then it is
                         triaxial (an ellipsoid with different radius in ln(r), theta and
                         phi).

        COMPUTES / RETURNS:

        The results are the grid and the density array:

          self.grid_ri       The radial grid cell interface radii
          self.grid_rc       The radial grid cell center radii
          self.grid_thetai   The theta grid cell interfaces
          self.grid_thetac   The theta grid cell centers
          self.grid_phii     The phi grid cell interfaces
          self.grid_phic     The phi grid cell centers
          self.grid_vol      The cell volumes (3D array with indices [ir,itheta,iphi])
          self.rho           The density values (3D array with indices [ir,itheta,iphi])
                             Note that the total mass is (self.rho*self.vol).sum()

        For analytic purposes these arrays may also be useful:

          self.cellmass_unconvolved    The mass in each cell, before convolution
          self.cellmass                The mass in each cell, after convolution
                                            (the self.rho = self.cellmass/self.grid_vol)
        """

        # Import information about the smoothing kernel

        if(np.isscalar(hrkernel)):
            self.hrkernel = np.zeros(3) + hrkernel
        else:
            self.hrkernel = hrkernel

        # Import the SPH particles

        self.sph_r = None
        if xyz is not None:

            # Convert from (x,y,z) coordinates to (r,theta,phi) coordinates

            assert rthphi is None, 'ERROR: xyz and rthphi are mutually exclusive'
            self.sph_r      = np.sqrt(xyz[:,0]**2+xyz[:,1]**2+xyz[:,2]**2)
            self.sph_rcyl   = np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
            self.sph_theta  = np.pi/2-np.arctan(xyz[:,2]/(self.sph_rcyl+1e-90))
            self.sph_phi    = np.arctan(xyz[:,1]/(xyz[:,0]+1e-90))
            self.sph_phi[xyz[:,0]<0] += np.pi
            self.sph_phi[self.sph_phi<0] += 2*np.pi
            self.sph_phi[self.sph_phi>2*pi] -= 2*np.pi
            assert len(np.where(self.sph_phi<0)[0])==0, 'ERROR: Phi computation went wrong'
            assert len(np.where(self.sph_phi>2*np.pi)[0])==0, 'ERROR: Phi computation went wrong'
        if rthphi is not None:

            # Import the SPH (r,theta,phi) coordinates

            assert xyz is None, 'ERROR: xyz and rthphi are mutually exclusive'
            self.sph_r      = rthphi[:,0]
            self.sph_theta  = rthphi[:,1]
            self.sph_phi    = rthphi[:,2]
            self.sph_phi[self.sph_phi<0] += 2*np.pi
            self.sph_phi[self.sph_phi>2*pi] -= 2*np.pi
            assert len(np.where(self.sph_phi<0)[0])==0, 'ERROR: Phi computation went wrong'
            assert len(np.where(self.sph_phi>2*np.pi)[0])==0, 'ERROR: Phi computation went wrong'
        assert self.sph_r is not None, 'ERROR: Must specify SPH particles (either rthphi or xyz).'
        self.nsph   = len(self.sph_r)
        self.sph_lr = np.log(self.sph_r)

        # Set up the spatial spherical grid, first the cell interfaces

        if np.isscalar(masses):
            self.masses = np.zeros(self.nsph) + masses
        else:
            self.masses = masses
        if rin is None:
            rin = self.sph_r[self.sph_r>0].min() * np.exp(-self.hrkernel[0])
        if rout is None:
            rout = self.sph_r.max() * np.exp(self.hrkernel[0])
        self.grid_ri,self.grid_thetai,self.grid_phii = self.uniform_spher_loggrid(rin,rout,hrup,nr,ntheta,nphi)
        self.grid_nr     = nr
        self.grid_ntheta = ntheta
        self.grid_nphi   = nphi

        # Cell center grid

        self.grid_rc     = np.sqrt(self.grid_ri[:-1]*self.grid_ri[1:])
        self.grid_thetac = 0.5 * ( self.grid_thetai[:-1] + self.grid_thetai[1:] )
        self.grid_phic   = 0.5 * ( self.grid_phii[:-1] + self.grid_phii[1:] )

        # Log(r) grid

        self.grid_lrc    = np.log(self.grid_rc)
        self.grid_lri    = np.log(self.grid_ri)

        # Cell volumes
        self.grid_vol    = (1./3.)*(self.grid_ri[1:,None,None]**3-self.grid_ri[:-1,None,None]**3) * \
                           np.abs( np.cos(self.grid_thetai[None,:-1,None]) - np.cos(self.grid_thetai[None,1:,None]) ) * \
                           np.abs( self.grid_phii[None,None,1:] - self.grid_phii[None,None,:-1] )

        # Now add the particles

        self.add_particles_to_grid()

        # And convolve

        self.convolve_particles()

        # Finally, compute the densities

        self.compute_density()

        # Give warning if necessary

        if(len(self.sph_excluded)>0):
            print('WARNING: Some SPH particles are outside of the domain (see self.sph_excluded)')
        
    def uniform_spher_loggrid(self,rin,rout,hrup,nr,nth,nph):
        ri    = rin * (rout/rin)**np.linspace(0,1,nr+1)
        thup  = np.pi/2-hrup
        thlo  = np.pi/2+hrup
        thi   = thup + (thlo-thup)*np.linspace(0,1,nth+1)
        phi   = 2*np.pi*np.linspace(0,1,nph+1)
        return ri,thi,phi

    def add_particles_to_grid(self):

        # Get arrays

        sph_lr     = self.sph_lr
        sph_th     = self.sph_theta
        sph_phi    = self.sph_phi
        grid_ri    = self.grid_ri
        grid_lri   = self.grid_lri
        grid_thi   = self.grid_thetai
        grid_phi   = self.grid_phii
        masses     = self.masses
        
        # Find particles in the grid

        sph_ir     = np.array(np.interp(sph_lr,grid_lri,np.arange(len(grid_lri))),dtype=int)
        sph_ith    = np.array(np.interp(sph_th,grid_thi,np.arange(len(grid_thi))),dtype=int)
        sph_iphi   = np.array(np.interp(sph_phi,grid_phi,np.arange(len(grid_phi))),dtype=int)

        # Check

        # sph_r_eps  = (sph_r-grid_ri[sph_ir])/(grid_ri[sph_ir+1]-grid_ri[sph_ir])
        # sph_th_eps = (sph_th-grid_thi[sph_ith])/(grid_thi[sph_ith+1]-grid_thi[sph_ith])
        # sph_phi_eps= (sph_phi-grid_phi[sph_iphi])/(grid_phi[sph_iphi+1]-grid_phi[sph_iphi])
        # assert len(np.where(sph_r_eps<0)[0])==0, 'Error in r-grid-finding'
        # assert len(np.where(sph_r_eps>1)[0])==0, 'Error in r-grid-finding'
        # assert len(np.where(sph_th_eps<0)[0])==0, 'Error in th-grid-finding'
        # assert len(np.where(sph_th_eps>1)[0])==0, 'Error in th-grid-finding'
        # assert len(np.where(sph_phi_eps<0)[0])==0, 'Error in phi-grid-finding'
        # assert len(np.where(sph_phi_eps>1)[0])==0, 'Error in phi-grid-finding'
        
        # Add particles to the grid

        print('Adding particles to the grid (may take a while)')
        self.sph_excluded = []
        self.cellmass_unconvolved = np.zeros((self.grid_nr,self.grid_ntheta,self.grid_nphi))
        for isph in tqdm(range(self.nsph)):
            ir   = sph_ir[isph]
            ith  = sph_ith[isph]
            iphi = sph_iphi[isph]
            if( (ir>=0) and (ir<self.grid_nr) and \
                (ith>=0) and (ith<self.grid_ntheta) and \
                (iphi>=0) and (iphi<self.grid_nphi) ):                
                self.cellmass_unconvolved[ir,ith,iphi] += masses[isph]
            else:
                self.sph_excluded.append(isph)
                
    def convolve_particles(self):

        # Get arrays

        grid_lri   = self.grid_lri
        grid_thi   = self.grid_thetai
        grid_phi   = self.grid_phii

        # Create a kernel

        dlgr = grid_lri[1]-grid_lri[0]
        dth  = grid_thi[1]-grid_thi[0]
        dphi = grid_phi[1]-grid_phi[0]
        self.kern = self.kernel3d(dlgr,dth,dphi,self.hrkernel)
        
        # Convolve

        print('Convolving with SPH kernel using scipy.signal.fftconvolve')
        self.cellmass = fftconvolve(self.cellmass_unconvolved,self.kern)
        
        # Compute the number of ghost (padding) cells

        nrg = (self.kern.shape[0]-1)//2
        ntg = (self.kern.shape[1]-1)//2
        npg = (self.kern.shape[2]-1)//2
        
        # Implement the periodic boundary in phi

        if npg>0:
            self.cellmass[:,:,npg:2*npg]   += self.cellmass[:,:,-npg:]
            self.cellmass[:,:,-2*npg:-npg] += self.cellmass[:,:,:npg]
        
        # Cut off the padding

        if nrg>0:
            self.cellmass = self.cellmass[nrg:-nrg,:,:]
        if ntg>0:
            self.cellmass = self.cellmass[:,ntg:-ntg,:]
        if npg>0:
            self.cellmass = self.cellmass[:,:,npg:-npg]

    def kernel1d(self,r,h):
        
        # Spline kernel of Gadget-2 (Springel MNRAS 364, 1105, 2005, Eq. 4)
        
        x  = r/h
        w0 = 8/(np.pi*h**3)
        w  = np.zeros_like(r)
        i0 = x>=1
        i1 = np.logical_and(x>=0,x<0.5)
        i2 = np.logical_and(x>=0.5,x<1)
        w[i1] = 1-6*x[i1]**2+6*x[i1]**3
        w[i2] = 2*(1-x[i2])**3
        return w
    
    def kernel3d(self,dlgr,dth,dphi,hr):

        # Create the 3D stencil of the kernel
        
        if np.isscalar(hr):
            hr3 = hr*np.array([1,1,1])
        else:
            hr3 = hr
        n_r   = (int(2*hr3[0]/dlgr)//2)*2+1
        n_th  = (int(2*hr3[1]/dth)//2)*2+1
        n_phi = (int(2*hr3[2]/dphi)//2)*2+1
        rc    = np.arange(-n_r//2+1,n_r//2+1)*dlgr/hr3[0]
        thc   = np.arange(-n_th//2+1,n_th//2+1)*dth/hr3[1]
        phic  = np.arange(-n_phi//2+1,n_phi//2+1)*dphi/hr3[2]
        r     = np.sqrt(rc[:,None,None]**2+thc[None,:,None]**2+phic[None,None,:]**2)
        w     = self.kernel1d(r,1)
        w     = w / w.sum()
        return w

    def compute_density(self):
        self.rho = self.cellmass / self.grid_vol

# #
# # UNCOMMENT THE FOLLOWING FOR AN EXAMPLE OF USAGE
# #
# #
# # Creat the random number generator for creating the dummy SPH model
# #
# rng     = default_rng()
# 
# #
# # Parameters of the dummy model
# #
# mdisk      = 1e-5*MS
# hpr        = 0.05
# hrkernel   = 0.03*np.array([2,1,20])
# 
# #
# # Dummy model for the SPH particles
# #
# nsph       = 1000000
# sph_r      = au*10**(2*rng.uniform(size=nsph))
# sph_th     = np.pi/2 - hpr*rng.standard_normal(size=nsph)
# sph_phi    = 2*np.pi*rng.uniform(size=nsph)
# sph_lr     = np.log(sph_r)
# rthphi     = np.vstack((sph_r,sph_th,sph_phi)).T
# masses     = mdisk/nsph
# 
# #
# # Set up the grid
# #
# grid_nr    = 200
# grid_nth   = 100
# grid_nph   = 100
# grid_rin   = 0.9*au
# grid_rout  = 110*au
# grid_hrup  = 0.25
# 
# #
# # Make the object to convert to spherical grid
# #
# sphgrid    = sph_to_sphergrid(rthph=rthphi,masses=masses,rin=grid_rin,rout=grid_rout,nr=grid_nr,
#                               ntheta=grid_nth,nphi=grid_nph,hrup=grid_hrup,hrkernel=hrkernel)
# 
# #
# # Show
# #
# print('Plotting the two grids')
# plt.figure()
# plt.imshow(sphgrid.cellmass_unconvolved.sum(axis=2).T,origin='lower')
# plt.figure()
# plt.imshow(sphgrid.cellmass.sum(axis=2).T,origin='lower')
# plt.figure()
# plt.imshow((sphgrid.rho.mean(axis=2)*sphgrid.grid_rc[:,None]**3).T,origin='lower')
# plt.figure()
# plt.plot(sphgrid.rho[100,50,:])
# plt.show()
