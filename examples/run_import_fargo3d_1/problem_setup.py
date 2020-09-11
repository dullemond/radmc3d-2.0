#
# Model setup for reading a 2-D (r,phi) model from FARGO3D with dust
# into RADMC-3D.
#
import numpy as np
import matplotlib.pyplot as plt
import readfargo as fg
#
# A simple grid refinement function
#
def grid_refine_inner_edge(x_orig,nlev,nspan):
    x     = x_orig.copy()
    rev   = x[0]>x[1]
    for ilev in range(nlev):
        x_new = 0.5 * ( x[1:nspan+1] + x[:nspan] )
        x_ref = np.hstack((x,x_new))
        x_ref.sort()
        x     = x_ref
        if rev:
            x = x[::-1]
    return x
#
# Some natural constants
#
au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]
ss  = 5.6703e-5      # Stefan-Boltzmann const  [erg/cm^2/K^4/s]
kk  = 1.3807e-16     # Bolzmann's constant     [erg/K]
mp  = 1.6726e-24     # Mass of proton          [g]
GG  = 6.67408e-08    # Gravitational constant  [cm^3/g/s^2]
pi  = np.pi          # Pi

#
# Star parameters
#
mstar    = 2.4*ms
rstar    = 2.4*rs
tstar    = 1e4
pstar    = np.array([0.,0.,0.])

#
# Read a FARGO3D frame
#
dir      = '/Users/cornelisdullemond/science/software/fargo3d/outputs/fargo_multifluid'
itime    = 85
r0       = 60*au    # The radius corresponding to '1' in dimensionless units
fargo    = fg.frame(itime,rhodust=True,dir=dir)
fargo.convert_to_cgs(mstar,r0)
#fargo.show(q=fargo.sigma_gas_cgs)
#fargo.show(q=fargo.sigma_dust_cgs[1])

#
# Make the RADMC-3D radial and phi grids
#
ri       = 0.5*(fargo.r_cgs[1:]+fargo.r_cgs[:-1])
ri       = np.hstack([2*ri[0]-ri[1],ri,2*ri[-1]-ri[-2]])
rc       = 0.5 * ( ri[:-1] + ri[1:] )
nr       = len(rc)   # Recompute nr, because of refinement at inner edge
r        = rc        # The radial grid of the analytic disk model (below)
phii     = 0.5*(fargo.phi[1:]+fargo.phi[:-1])
phii     = np.hstack([2*phii[0]-phii[1],phii,2*phii[-1]-phii[-2]])
phic     = 0.5 * ( phii[:-1] + phii[1:] )
r_2d,phi_2d = np.meshgrid(rc,phic,indexing='ij')

#
# Now make a simple analytical disk model roughly along the
# lines of Chiang & Goldreich (1997), but with just a single
# vertical layer and with a constant radiative incidence angle.
#
flang    = 0.05                        # The assumed constant radiative incidence angle
lstar    = 4*pi*rstar**2*ss*tstar**4   # Stellar luminosity
firr     = flang*lstar/(4*pi*r**2)     # Irradiative flux
tmid     = (firr/ss)**0.25             # Estimate of midplane temperature
cs       = np.sqrt(kk*tmid/(2.3*mp))   # Isothermal sound speed at midplane
omk      = np.sqrt(GG*mstar/r**3)      # The Kepler angular frequency
hp       = cs/omk                      # The pressure scale height
hpr      = hp/r                        # The dimensionless hp

#
# Compare the results to the AspectRatio and FlaringIndex of
# the FARGO3D setup. Note that we compute them using a simple
# Chiang & Goldreich like model, but this also has a free
# parameter: flang (the incidence angle). Since we keep this
# angle constant, the FlaringIndex should be 0.25, because
# T ~ r**(-1) --> c_s ~ r**(-0.5), and Omega_K ~ r**(-1.5),
# so H_p/r = c_s / (Omega_K * r) ~ r**(0.25).
#
AspectRatio  = np.interp(r0,r,hpr)     # The H_p/r at r=r0
flidx        = ( np.log(hpr) - np.log(AspectRatio) ) / ( np.log(r) - np.log(r0) )
FlaringIndex = np.interp(r0,r,flidx)   # The Flidx at r=r0
print('Please check that the RADMC-3D setup is roughly consistent with the FARGO3D setup:')
print('  AspectRatio in RADMC-3D model at R0  = {}'.format(AspectRatio))
print('  FlaringIndex in RADMC-3D model at R0 = {}'.format(FlaringIndex))
print('If these values are vastly different from those in the FARGO3D .par file,')
print('you should either modify the parameters in the RADMC-3D model, or the parameters')
print('in the FARGO3D model, to make them both mutually self-consistent.')

#
# Next we have to decide how to vertically distribute the dust.
# To do this self-consistently with the FARGO3D model, we have
# to apply Eq. (19) of Fromang & Nelson (2009) A&A 496, 597.
# The (Omega*tau_s)_mid in that equation is the Stokes number.
#
# In the standard FARGO3D model, you specify the inverse Stokes
# number of the dust grains, not the grain size. However, for
# the radiative transfer, we need the grain size, which we keep
# fixed everywhere, and the Stokes number should then be variable
# dependent on the Sigma_gas and T_gas. So you may need to
# modify FARGO3D to keep a_grain constant instead of Stokes.
# But you can also estimate roughly which grain size belongs
# to which Stokes number. Here, however, we specify the grain
# radius in micron.
#
# NOTE: These have to be consistent with the opacity files
#       dustkappa_xxxx.inp you specify. Also the material
#       density must be consistent with that of the material
#       of the dust in dustkappa_xxxx.inp.
#
a_grains_mic= np.array([1e2,1e3,1e4])    # Grain radii in micron
a_grains    = a_grains_mic*1e-4          # Grain radii in cm
xi_grains   = 3.7                        # Grain material density in g/cm^3
m_grains    = (4*np.pi/3)*xi_grains*a_grains**3  # Grain masses
s_grains    = np.pi*a_grains**2                  # Grain cross sections

#
# We now compute the Stokes number of these grains at the
# midplane, based on the gas density and the temperature
#
omk_2d      = np.sqrt(GG*mstar/r_2d**3)  # The Kepler frequency Omega_K in 1/s
cs_2d,dum   = np.meshgrid(cs,phic,indexing='ij')  # The isothermal sound speed c_s in cm/s
hp_2d       = cs_2d/omk_2d               # The press scale height H_p = c_s / Omega_K in cm
hpr_2d      = hp_2d/r_2d                 # The aspect ratio H_p / r
vth_2d      = np.sqrt(8/np.pi)*cs_2d     # The thermal velocity of H2 gas particles in cm/s
rho_gas_mid_2d = fargo.sigma_gas_cgs/(np.sqrt(2*pi)*hp_2d)  # Midplane gas density in g/cm^3
Stokes_mid_2d = []
for m,s in zip(m_grains,s_grains):
    St = (3./4.)*(m/s)*omk_2d/(rho_gas_mid_2d*vth_2d)  # Stokes number assuming Epstein drag law
    Stokes_mid_2d.append(St)

#
# Compute an estimate of the mean free path of the gas molecules, to be
# able to check if the Epstein regime is justified
#
mugas        = 2.3                       # Mean molecular gas in a protoplanetary disk in mp units
n_gas_mid_2d = rho_gas_mid_2d/(mugas*mp) # Number density of gas particles
siggas       = 2e-15                     # Cross section of H2 molecule (assume that to be for all gas molecules)
lmfp         = 1 / (n_gas_mid_2d*siggas) # Mean free path of the gas molecules
assert (a_grains.max()/lmfp).max()<1, "Warning: Epstein regime not everywhere valid"







"""
#
# The FARGO3D model has three dust species: Let's assume the first one consists
# of 0.1 micron dust grains, which we assume will be vertically well-mixed
# with the gas. Let's assume the
# consisting of 100 micron dust grains which forms a settled dust midplane
# layer that is 'settfact' settled.
#
fracbig  = 0.90              # Fraction of dust that is the big grain dust
settfact = 0.1               # Settling factor of the big grains
hp_biggr = hp*settfact
hpr_biggr= hpr*settfact

#
# Vertical grid parameters (theta-grid in spherical coordinates)
#
ntheta   = 32
zrmax    = 0.5
thetaup  = np.pi*0.5 - zrmax

#
# Make the theta coordinate, and refine near the midplane
# (to resolve the dust layer)
#
thetai   = np.linspace(thetaup,0.5e0*np.pi,ntheta+1)
zr       = thetai[::-1]
nlev_zr  = 5         # Grid refinement at the midplane: nr of cycles
nspan_zr = 3         # Grid refinement at the midplane: nr of cells each cycle
zr       = grid_refine_inner_edge(zr,nlev_zr,nspan_zr)
thetai   = zr[::-1]
thetac   = 0.5 * ( thetai[:-1] + thetai[1:] )
ntheta   = len(thetac)

#
# Make the phi coordinate
#
nphi     = 1
phii     = np.linspace(0.e0,np.pi*2.e0,nphi+1)
phic     = 0.5 * ( phii[0:nphi] + phii[1:nphi+1] )

#
# Make the 2-D grid (actually 3-D but with axisymmetry)
#
qq       = np.meshgrid(rc,thetac,phic,indexing='ij')
rr       = qq[0]
tt       = qq[1]
zr       = np.pi/2.e0 - qq[1]

#
# Expand the 1-D analytic model to 2-D
#
sigmad_3d  = np.meshgrid(sigmad,thetac,phic,indexing='ij')[0]
hh         = np.meshgrid(hp,thetac,phic,indexing='ij')[0]
hhr        = np.meshgrid(hpr,thetac,phic,indexing='ij')[0]
hh_biggr   = np.meshgrid(hp_biggr,thetac,phic,indexing='ij')[0]
hhr_biggr  = np.meshgrid(hpr_biggr,thetac,phic,indexing='ij')[0]

#
# Make the dust density model
#
rhod_smlgr = (1.-fracbig) * ( sigmad_3d / (np.sqrt(2.e0*np.pi)*hh) ) * np.exp(-(zr**2/hhr**2)/2.e0)
rhod_biggr = fracbig *      ( sigmad_3d / (np.sqrt(2.e0*np.pi)*hh_biggr) ) * np.exp(-(zr**2/hhr_biggr**2)/2.e0)

#
# Monte Carlo parameters
#
nphot    = 1000000

#
# Write the wavelength_micron.inp file
#
lam1     = 0.1e0
lam2     = 7.0e0
lam3     = 25.e0
lam4     = 1.0e4
n12      = 20
n23      = 100
n34      = 30
lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
lam      = np.concatenate([lam12,lam23,lam34])
nlam     = lam.size
#
# Write the wavelength file
#
with open('wavelength_micron.inp','w+') as f:
    f.write('%d\n'%(nlam))
    for value in lam:
        f.write('%13.6e\n'%(value))
#
#
# Write the stars.inp file
#
with open('stars.inp','w+') as f:
    f.write('2\n')
    f.write('1 %d\n\n'%(nlam))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
    for value in lam:
        f.write('%13.6e\n'%(value))
    f.write('\n%13.6e\n'%(-tstar))
#
# Write the grid file
#
with open('amr_grid.inp','w+') as f:
    f.write('1\n')                       # iformat
    f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
    f.write('100\n')                     # Coordinate system: spherical
    f.write('0\n')                       # gridinfo
    f.write('1 1 0\n')                   # Include r,theta coordinates
    f.write('%d %d %d\n'%(nr,ntheta,1))  # Size of grid
    for value in ri:
        f.write('%13.6e\n'%(value))      # X coordinates (cell walls)
    for value in thetai:
        f.write('%17.10e\n'%(value))     # Y coordinates (cell walls) (use higher precision here)
    for value in phii:
        f.write('%13.6e\n'%(value))      # Z coordinates (cell walls)
#
# Write the density file
#
with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nr*ntheta*nphi))     # Nr of cells
    f.write('2\n')                       # Nr of dust species
    data = rhod_smlgr.ravel(order='F')   # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
    data = rhod_biggr.ravel(order='F')   # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
#
# Dust opacity control file
#
with open('dustopac.inp','w+') as f:
    f.write('2               Format number of this file\n')
    f.write('2               Nr of dust species\n')
    f.write('============================================================================\n')
    f.write('1               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('0.1_micron      Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')
    f.write('1               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('100_micron      Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')
#
# Write the radmc3d.inp control file
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot = %d\n'%(nphot))
    f.write('scattering_mode_max = 1\n')
    f.write('iranfreqmode = 1\n')


"""
