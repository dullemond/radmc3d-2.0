#
# Import NumPy for array handling
#
import numpy as np
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
# Get angle from (x,y)
#
#def get_angle(x,y):
#    if x!=0:
#        ang  = np.arctan(y/x)
#        if x<0:
#            if y>=0:
#                ang+=np.pi
#            else:
#                ang-=np.pi
#    else:
#        if y>=0:
#            ang  = np.pi/2
#        else:
#            ang  = -np.pi/2
#    return ang
def get_angle(x,y,positive=False):
    if np.isscalar(x):
        xx   = np.array([x])
        yy   = np.array([y])
    else:
        xx   = np.array(x)
        yy   = np.array(y)
    xx  += 1e-90*np.pi
    ang  = np.arctan(yy/xx)
    ang[xx<0] -= np.pi
    if positive:
        ang[ang<0] += 2*np.pi
    else:
        ang[ang<=-np.pi] += 2*np.pi
    if len(ang)==1:
        ang = ang[0]
    return ang
#
# Simple rotation counterclockwise
#
def rotate(x,y,ang):
    cos = np.cos(ang)
    sin = np.sin(ang)
    xp  = cos*x - sin*y
    yp  = sin*x + cos*y
    return xp,yp
#
# For disk warping:
# Coordinate transformation to rotate (theta,phi) --> (theta',phi'),
# for a normal vector as a function of radius l(r) = (l_x,l_y,l_z),
# where the normal vector is in the cartesian coordinate system.
# A flat disk with midplane at the theta=pi/2 location is given by
# l(r) = (0,0,1).
#
def warped_coordinate_transformation(theta_3d,phi_3d,r_1d,lunit_1d):
    assert theta_3d.shape[0]==len(r_1d), 'Error in coordinate transformation: nr of elements of r_1d array does not match.'
    assert len(lunit_1d.shape)==2, 'Error in coordinate transformation: l_1d array must have [nr,3] elements.'
    assert lunit_1d.shape[0]==len(r_1d), 'Error in coordinate transformation: nr of elements of l_1d array does not match.'
    assert lunit_1d.shape[1]==3, 'Error in coordinate transformation: l_1d array must have [nr,3] elements.'
    theta_prime_3d = np.zeros_like(theta_3d)
    phi_prime_3d   = np.zeros_like(phi_3d)
    for ir in range(len(r_1d)):
        l     = lunit_1d[ir,:]
        ll    = np.sqrt(l[0]**2+l[1]**2+l[2]**2)
        l     = l/ll
        if l[0]==0 and l[1]==0:
            theta_prime_3d[ir,:,:] = theta_3d[ir,:,:].copy()
            phi_prime_3d[ir,:,:]   = phi_3d[ir,:,:].copy()
        else:
            angh  = get_angle(l[0],l[1])
            incl  = get_angle(l[2],np.sqrt(l[0]**2+l[1]**2))
            r     = r_1d[ir]
            theta = theta_3d[ir,:,:]
            phi   = phi_3d[ir,:,:]
            x     = r*np.sin(theta)*np.cos(phi)
            y     = r*np.sin(theta)*np.sin(phi)
            z     = r*np.cos(theta)
            xt,yt = rotate(x,y,-angh)
            xi,zi = rotate(xt,z,incl)
            xn,yn = rotate(xi,yt,angh)
            zn    = zi
            #rr    = np.sqrt(xn**2+yn**2+zn**2)
            pp    = get_angle(xn,yn,positive=True)
            tt    = get_angle(zn,np.sqrt(xn**2+yn**2))
            theta_prime_3d[ir,:,:] = tt
            phi_prime_3d[ir,:,:]   = pp
    return theta_prime_3d,phi_prime_3d



#r     = np.array([1,1])
#theta = np.zeros((2,2,2))+np.pi/2
#phi   = np.zeros((2,2,2))+np.pi/4
#lunit = np.array([[0.1,0.1,1],[0.1,0.1,1]])
#tt,pp = warped_coordinate_transformation(theta,phi,r,lunit)
#print(tt)
#print(pp)


#
# Some natural constants
#
au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]
#
# Monte Carlo parameters
#
nphot_therm = 1000000
nphot_scat  = 1000000
#
# Grid parameters
#
nr       = 100
ntheta   = 64
nphi     = 128
rin      = 0.5*au
rout     = 100*au
#thetaup  = np.pi*0.5 - 0.7e0
thetaup  = 0.1       # Theta grid starting point (0=pole, but singular, so choose >0)
nlev_rin = 8
nspan_rin= 3
#
# Disk parameters
#
sigmag0  = 1e1               # Sigma gas at 1 AU
sigmad0  = sigmag0 * 0.01    # Sigma dust at 1 AU
plsig    = -1.0e0            # Powerlaw of the surface density
hr0      = 0.05              # H_p/r at 1 AU
plh      = 0.1               # Powerlaw of flaring
#
# Star parameters
#
mstar    = 2.4*ms
rstar    = 2.4*rs
tstar    = 1e4
pstar    = np.array([0.,0.,0.])
#
# Make the coordinates
#
ri       = np.logspace(np.log10(rin),np.log10(rout),nr+1)
ri       = grid_refine_inner_edge(ri,nlev_rin,nspan_rin)   # Refinement at inner edge
thetai   = np.linspace(thetaup,np.pi-thetaup,ntheta+1)
phii     = np.linspace(0.e0,np.pi*2.e0,nphi+1)
rc       = 0.5 * ( ri[:-1] + ri[1:] )
thetac   = 0.5 * ( thetai[:-1] + thetai[1:] )
phic     = 0.5 * ( phii[:-1] + phii[1:] )
nr       = len(rc)     # Recompute nr, because of refinement at inner edge
#
# Make the grid
#
qq       = np.meshgrid(rc,thetac,phic,indexing='ij')
rr       = qq[0]
tt0      = qq[1]
pp0      = qq[2]
#
# Warp the coordinates
#
ll       = np.zeros((nr,3))
incl     = np.linspace(0,45.,nr)*np.pi/180.
ll[:,2]  = np.cos(incl)
ll[:,1]  = np.sin(incl)
tt,pp    = warped_coordinate_transformation(tt0,pp0,rc,ll)
zr       = np.pi/2.e0 - tt
#
# Make the dust density model
#
sigmad   = sigmad0 * (rr/au)**plsig
hhr      = hr0 * (rr/au)**plh
hh       = hhr * rr
rhod     = ( sigmad / (np.sqrt(2.e0*np.pi)*hh) ) * np.exp(-(zr**2/hhr**2)/2.e0)
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
    f.write('1 1 1\n')                   # Include r,theta coordinates
    f.write('%d %d %d\n'%(nr,ntheta,nphi))  # Size of grid
    for value in ri:
        f.write('%13.6e\n'%(value))      # X coordinates (cell walls)
    for value in thetai:
        f.write('%13.6e\n'%(value))      # Y coordinates (cell walls)
    for value in phii:
        f.write('%13.6e\n'%(value))      # Z coordinates (cell walls)
#
# Write the density file
#
with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nr*ntheta*nphi))     # Nr of cells
    f.write('1\n')                       # Nr of dust species
    data = rhod.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
#
# Dust opacity control file
#
with open('dustopac.inp','w+') as f:
    f.write('2               Format number of this file\n')
    f.write('1               Nr of dust species\n')
    f.write('============================================================================\n')
    f.write('1               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('silicate        Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')
#
# Write the radmc3d.inp control file
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot = %d\n'%(nphot_therm))
    f.write('nphot_scat = %d\n'%(nphot_scat))
    f.write('scattering_mode_max = 1\n')
    f.write('iranfreqmode = 1\n')

