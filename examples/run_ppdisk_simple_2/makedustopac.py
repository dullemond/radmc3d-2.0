from bhmie import *
import numpy as np
from scipy.interpolate import interp1d
import math

def compute_opac_mie(optconst_file,matdens,agraincm,lamcm,
                     theta=None,logawidth=None,wfact=3.0,na=20,
                     chopforward=0.0,errtol=0.01,verbose=False,
                     extrapolate=False):
    """
    Compute dust opacity with Mie theory based on the optical constants
    in the optconst_file. Optionally also the scattering phase function
    in terms of the Mueller matrix elements can be computed. To smear out
    the resonances that appear due to the perfect sphere shape, you can
    optionally smear out the grain size distribution a bit with setting
    the width of a Gaussian grain size distribution. 

    Arguments:
      optconst_file = File name of the optical constants file. This file
                      should contain three columns: first the wavelength
                      in micron, then the n-coefficient and then the 
                      k-coefficient. See Jena optical constants database:
                      http://www.astro.uni-jena.de/Laboratory/Database/databases.html
      matdens       = Material density in g/cm^3
      agraincm      = Grain radius in cm
      lamcm         = Wavelength grid in cm (a numpy array)
      theta         = Optional angular grid (a numpy array) between 0 and 180
                      which are the scattering angle sampling points at 
                      which the scattering phase function is computed.
      logawidth     = Optional: if set, the size agrain will instead be a 
                      sample of sizes around agrain. This helps to smooth out
                      the strong wiggles in the phase function and opacity
                      of spheres at an exact size. Since in Nature it rarely
                      happens that grains all have exactly the same size, this
                      is quite natural. The value of logawidth sets the width
                      of the Gauss in ln(agrain), so for logawidth<<1 this
                      give a real width of logawidth*agraincm. 
      wfact         = (default=3.0) Grid width of na sampling points in units
                      of logawidth. The Gauss distribution of grain sizes is 
                      cut off at agrain * exp(wfact*logawidth) and
                      agrain * exp(-wfact*logawidth).
      na            = (default=20) Number of size sampling points 
                      (if logawidth set).
      chopforward   = If >0 this gives the angle (in degrees from forward)
                      within which the scattering phase function should be
                      kept constant, essentially removing the strongly peaked
                      forward scattering. This is useful for large grains
                      (large ratio 2*pi*agraincm/lamcm) where the forward
                      scattering peak is extremely strong, yet extremely
                      narrow. If we are not interested in very forward-peaked
                      scattering (e.g. only relevant when modeling e.g. the
                      halo around the moon on a cold winter night), this will
                      remove this component and allow a lower angular grid
                      resolution for the theta grid.
      errtol        = Tolerance of the relative difference between kscat
                      and the integral over the zscat Z11 element over angle.
                      If this tolerance is exceeded, a warning is given.
      verbose       = If set to True, the code will give some feedback so
                      that one knows what it is doing if it becomes slow.
      extrapolate   = If set to True, then if the wavelength grid lamcm goes
                      out of the range of the wavelength grid of the 
                      optical constants file, then it will make a suitable
                      extrapolation: keeping the optical constants constant
                      for lamcm < minimum, and extrapolating log-log for
                      lamcm > maximum.

    Returns (all in the form of a single dictionary):
      kabs          = Absorption opacity kappa_abs_nu (a numpy array) in 
                      units of cm^2/gram
      kscat         = Scattering opacity kappa_abs_nu (a numpy array) in 
                      units of cm^2/gram
      gscat         = The <cos(theta)> g-factor of scattering

    Returns also (if theta grid is given):
      theta         = The theta grid itself (just a copy of what was given)
      zscat         = The components of the scattering Mueller matrix
                      Z_ij for each wavelength and each scattering angel.
                      The normalization of Z is such that kscat can be
                      reproduced (as can be checked) by the integral:
                      2*pi*int_{-1}^{+1}Z11(mu)dmu=kappa_scat.
                      For symmetry reasons only 6 elements of the Z
                      matrix are returned: Z11, Z12, Z22, Z33, Z34, Z44.
                      Note that Z21 = Z12 and Z43 = -Z34.
                      The scattering matrix is normalized such that 
                      if a plane wave with Stokes flux 
                         Fin = (Fin_I,Fin_Q,Fin_U,Fin_V)
                      hits a dust grain (which has mass mgrain), then
                      the scattered flux 
                         Fout = (Fout_I,Fout_Q,Fout_U,Fout_V)
                      at distance r from the grain at angle theta 
                      is given by
                         Fout(theta) = (mgrain/r^2) * Zscat . Fin
                      where . is the matrix-vector multiplication.
                      Note that the Stokes components must be such
                      that the horizontal axis in the "image" is 
                      pointing in the scattering plane. This means
                      that radiation with Fin_Q < 0 is scattered well,
                      because it is vertically polarized (along the
                      scattering angle axis), while radiation with 
                      Fin_Q > 0 is scatterd less well because it 
                      is horizontally polarized (along the scattering
                      plane). 
      kscat_from_z11= The kscat computed from the (above mentioned)
                      integral of Z11 over all angles. This should be
                      nearly identical to kscat if the angular grid
                      is sufficiently fine. If there are strong 
                      differences, this is an indication that the
                      angular gridding (the theta grid) is not fine
                      enough. But you should have then automatically
                      gotten a warning message as well (see errtol).

    If extrapolate is set to True, it will also return:
      wavmic        = The original wavelength grid from the optical 
                      constants file, with possibly an added extrapolated
                      value at each end.
      ncoef         = The optical constant n at that grid
      kcoef         = The optical constant k at that grid

    If logawidth is set to some value, then a size distribution is 
    used to smear out some of the strong oscillations in the 
    opacity and phase function due to the resonances. Then 
    it will also return this size distribution:
      agr           = The grain sizes
      wgt           = The averaging weights of these grain (not the masses!)
                      The sum of wgt.sum() must be 1.

    If chopforward>0 it will also return the unchopped versions:
      zscat_nochop  = The zscat before the forward scattering was chopped off
      kscat_nochop  = The kscat originally from the bhmie code
    """
    #
    # Load the optical constants
    #
    data = np.loadtxt(optconst_file)
    wavmic, ncoef, kcoef = data.T
    assert wavmic.size > 1, "Error: Optical constants file must have at least two rows with two different wavelengths."
    assert wavmic[1]!=wavmic[0], "Error: Optical constants file must have at least two rows with two different wavelengths."
    #
    # Check range, and if needed and requested, extrapolate the
    # optical constants to longer or shorter wavelengths
    #
    if extrapolate:
        wmin = np.min(lamcm)*1e4 * 0.999
        wmax = np.max(lamcm)*1e4 * 1.001
        if wmin < np.min(wavmic):
            if wavmic[0]<wavmic[1]:
                ncoef  = np.append([ncoef[0]],ncoef)
                kcoef  = np.append([kcoef[0]],kcoef)
                wavmic = np.append([wmin],wavmic)
            else:
                ncoef  = np.append(ncoef,[ncoef[-1]])
                kcoef  = np.append(kcoef,[kcoef[-1]])
                wavmic = np.append(wavmic,[wmin])
        if wmax > np.max(wavmic):
            if wavmic[0]<wavmic[1]:
                ncoef  = np.append(ncoef,[ncoef[-1]*math.exp((math.log(wmax)-math.log(wavmic[-1]))*
                                                             (math.log(ncoef[-1])-math.log(ncoef[-2]))/
                                                             (math.log(wavmic[-1])-math.log(wavmic[-2])))])
                kcoef  = np.append(kcoef,[kcoef[-1]*math.exp((math.log(wmax)-math.log(wavmic[-1]))*
                                                             (math.log(kcoef[-1])-math.log(kcoef[-2]))/
                                                             (math.log(wavmic[-1])-math.log(wavmic[-2])))])
                wavmic = np.append(wavmic,[wmax])
            else:
                ncoef  = np.append(ncoef,[ncoef[0]*math.exp((math.log(wmax)-math.log(wavmic[0]))*
                                                             (math.log(ncoef[0])-math.log(ncoef[1]))/
                                                             (math.log(wavmic[0])-math.log(wavmic[1])))])
                kcoef  = np.append(kcoef,[kcoef[0]*math.exp((math.log(wmax)-math.log(wavmic[0]))*
                                                             (math.log(kcoef[0])-math.log(kcoef[1]))/
                                                             (math.log(wavmic[0])-math.log(wavmic[1])))])
                wavmic = np.append([wmax],wavmic)
    else:
        assert np.min(lamcm) >= np.min(wavmic*1e-4), "Error: wavelength range out of range of the optical constants file.\n"
        assert np.max(lamcm) <= np.max(wavmic*1e-4), "Error: wavelength range out of range of the optical constants file.\n"
    #
    # Interpolate
    # Note: Must be within range, otherwise stop
    #
    f      = interp1d(np.log(wavmic*1e-4),np.log(ncoef))
    ncoefi = np.exp(f(np.log(lamcm)))
    f      = interp1d(np.log(wavmic*1e-4),np.log(kcoef))
    kcoefi = np.exp(f(np.log(lamcm)))
    #
    # Make the complex index of refraction
    #
    refidx = ncoefi + kcoefi*1j
    #
    # Make a grid of angles, if not already given by theta
    #
    if theta is None:
        angles = np.array([0.,90.,180.])  # Minimalistic angular grid
        assert chopforward==0.0, "Sorry: Chopping only possible if theta grid given."
    else:
        angles = theta
    nang = angles.size
    #
    # Check that the theta array goes from 0 to 180 or 
    # 180 to 0, and store which is 0 and which is 180
    #
    if angles[0]==0.0:
        assert angles[nang-1]==180, "Error: Angle grid must extend from 0 to 180 degrees."
    else:
        assert angles[0]==180, "Error: Angle grid must extend from 0 to 180 degrees."
        assert angles[nang-1]==0, "Error: Angle grid must extend from 0 to 180 degrees."
    #
    # Make a size distribution for the grains
    # If width is not set, then take just one size
    #
    if logawidth is None:
        agr   = np.array([agraincm])
        wgt   = np.array([1.0])
    else:
        agr   = np.exp(np.linspace(math.log(agraincm)-wfact*logawidth,
                                   math.log(agraincm)+wfact*logawidth,na))
        wgt   = np.exp(-0.5*((np.log(agr/agraincm))/logawidth)**2)
        wgt   = wgt/wgt.sum()
    #
    # Get the true number of grain sizes
    #
    nagr = agr.size
    #
    # Compute the geometric cross sections
    #
    siggeom = math.pi*agr*agr
    #
    # Compute the mass of the grain
    # 
    mgrain  = (4*math.pi/3.0)*matdens*agr*agr*agr
    #
    # Now prepare arrays
    #
    nlam  = lamcm.size
    kabs  = np.zeros(nlam)
    kscat = np.zeros(nlam)
    gscat = np.zeros(nlam)
    if theta is not None:
        zscat = np.zeros((nlam,nang,6))
        S11   = np.zeros(nang)
        S12   = np.zeros(nang)
        S33   = np.zeros(nang)
        S34   = np.zeros(nang)
        if chopforward>0:
            zscat_nochop = np.zeros((nlam,nang,6))
            kscat_nochop = np.zeros(nlam)
    #
    # Set error flag to False
    #
    error  = False
    errmax = 0.0
    kscat_from_z11 = np.zeros(nlam)
    #
    # Loop over wavelengths
    #
    for i in range(nlam):
        #
        # Message
        #
        if verbose:
            print("Doing wavelength %13.6e cm"%lamcm[i])
        #
        # Now loop over the grain sizes
        #
        for l in range(nagr):
            #
            # Message
            #
            if verbose and nagr>1:
                print("...Doing grain size %13.6e cm"%agr[l])
            #
            # Compute x
            #
            x = 2*math.pi*agr[l]/lamcm[i]
            #
            # Call the bhmie code
            #
            S1, S2, Qext, Qabs, Qsca, Qback, gsca = bhmie(x,refidx[i],angles)
            #
            # Add results to the averaging over the size distribution
            #
            kabs[i]   += wgt[l] * Qabs*siggeom[l]/mgrain[l]
            kscat[i]  += wgt[l] * Qsca*siggeom[l]/mgrain[l]
            gscat[i]  += wgt[l] * gsca
            #
            # If angles were set, then also compute the Z matrix elements
            #
            if theta is not None:
                #
                # Compute conversion factor from the Sxx matrix elements
                # from the Bohren & Huffman code to the Zxx matrix elements we
                # use (such that 2*pi*int_{-1}^{+1}Z11(mu)dmu=kappa_scat).
                # This includes the factor k^2 (wavenumber squared) to get 
                # the actual cross section in units of cm^2 / ster, and there 
                # is the mass of the grain to get the cross section per gram.
                #
                factor = (lamcm[i]/(2*math.pi))**2/mgrain[l]
                #
                # Compute the scattering Mueller matrix elements at each angle
                #
                S11[:]        = 0.5 * ( np.abs(S2[:])**2 + np.abs(S1[:])**2 )
                S12[:]        = 0.5 * ( np.abs(S2[:])**2 - np.abs(S1[:])**2 )
                S33[:]        = np.real(S2[:]*np.conj(S1[:]))
                S34[:]        = np.imag(S2[:]*np.conj(S1[:]))
                zscat[i,:,0] += wgt[l] * S11[:] * factor
                zscat[i,:,1] += wgt[l] * S12[:] * factor
                zscat[i,:,2] += wgt[l] * S11[:] * factor
                zscat[i,:,3] += wgt[l] * S33[:] * factor
                zscat[i,:,4] += wgt[l] * S34[:] * factor
                zscat[i,:,5] += wgt[l] * S33[:] * factor
        #
        # If possible, do a check if the integral over zscat is consistent 
        # with kscat
        #
        if theta is not None:
            mu  = np.cos(angles*math.pi/180.)
            dmu = np.abs(mu[1:nang]-mu[0:nang-1])
            zav = 0.5 * ( zscat[i,1:nang,0] + zscat[i,0:nang-1,0] )
            dum = 0.5 * zav*dmu
            sum = dum.sum() * 4 * math.pi
            kscat_from_z11[i] = sum
            err = abs(sum/kscat[i]-1.0)
            if err>errtol:
                error = True
                errmax = max(err,errmax)
        #
        # If the chopforward angle is set >0, then we will remove
        # excessive forward scattering from the opacity. The reasoning
        # is that extreme forward scattering is, in most cases, equivalent
        # to no scattering at all.
        #
        if chopforward>0:
            iang  = np.where(angles<chopforward)
            if angles[0]==0.0:
                iiang = np.max(iang)+1
            else:
                iiang = np.min(iang)-1
            zscat_nochop[i,:,:] = zscat[i,:,:]  # Backup
            kscat_nochop[i]     = kscat[i]      # Backup
            zscat[i,iang,0]     = zscat[i,iiang,0]
            zscat[i,iang,1]     = zscat[i,iiang,1]
            zscat[i,iang,2]     = zscat[i,iiang,2]
            zscat[i,iang,3]     = zscat[i,iiang,3]
            zscat[i,iang,4]     = zscat[i,iiang,4]
            zscat[i,iang,5]     = zscat[i,iiang,5]
            mu  = np.cos(angles*math.pi/180.)
            dmu = np.abs(mu[1:nang]-mu[0:nang-1])
            zav = 0.5 * ( zscat[i,1:nang,0] + zscat[i,0:nang-1,0] )
            dum = 0.5 * zav*dmu
            sum = dum.sum() * 4 * math.pi
            kscat[i] = sum
    #
    # If error found, then warn
    #
    if error:
        print("Warning: Angular integral of Z11 is not equal to kscat at all wavelength. ")
        print("Maximum error = %13.6e"%errmax)
        if chopforward>0:
            print("But I am using chopforward to remove strong forward scattering, and then renormalized kapscat.")
    #                
    # Now return what we computed in a dictionary
    #
    package = {"lamcm":lamcm, "kabs":kabs, "kscat":kscat, 
               "gscat":gscat, "matdens":matdens, "agraincm":agraincm}
    if theta is not None:
        package["zscat"] = np.copy(zscat)
        package["theta"] = np.copy(angles)
        package["kscat_from_z11"] = np.copy(kscat_from_z11)
    if extrapolate:
        package["wavmic"] = np.copy(wavmic)
        package["ncoef"] = np.copy(ncoef)
        package["kcoef"] = np.copy(kcoef)
    if nagr>1:
        package["agr"] = np.copy(agr)
        package["wgt"] = np.copy(wgt)
        package["wfact"] = wfact
        package["logawidth"] = logawidth
    if chopforward>0:
        package["zscat_nochop"] = np.copy(zscat_nochop)
        package["kscat_nochop"] = np.copy(kscat_nochop)
    return package


def write_radmc3d_scatmat_file(package,name):
    """
    The RADMC-3D radiative transfer package
      http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/
    can perform dust continuum radiative transfer for diagnostic purposes.
    It is designed for astronomical applications. The code
    needs the opacities in a particular form. This subroutine 
    writes the opacities out in that form. It will write it to
    the file dustkapscatmat_<name>.inp.
    """
    filename = 'dustkapscatmat_'+name+'.inp'
    with open(filename,'w') as f:
        f.write('# Opacity and scattering matrix file for '+name+'\n')
        f.write('# Please do not forget to cite in your publications the original paper of these optical constant measurements\n')
        f.write('# Made with the makedustopac.py code by Cornelis Dullemond\n')
        f.write('# using the bhmie.py Mie code of Bohren and Huffman (python version by Cornelis Dullemond, from original bhmie.f code by Bruce Draine)\n')
        f.write('# Grain size = %13.6e cm\n'%(package['agraincm']))
        f.write('# Material density = %6.3f g/cm^3\n'%(package['matdens']))
        f.write('1\n')  # Format number
        f.write('%d\n'%(package['lamcm'].size))
        f.write('%d\n'%(package['theta'].size))
        f.write('\n')
        for i in range(package['lamcm'].size):
            f.write('%13.6e %13.6e %13.6e %13.6e\n'%(package['lamcm'][i]*1e4,
                                                     package['kabs'][i],
                                                     package['kscat'][i],
                                                     package['gscat'][i]))
        f.write('\n')
        for j in range(package['theta'].size):
            f.write('%13.6e\n'%(package['theta'][j]))
        f.write('\n')
        for i in range(package['lamcm'].size):
            for j in range(package['theta'].size):
                f.write('%13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n'%
                        (package['zscat'][i,j,0],package['zscat'][i,j,1],
                         package['zscat'][i,j,2],package['zscat'][i,j,3],
                         package['zscat'][i,j,4],package['zscat'][i,j,5]))
        f.write('\n')


def write_radmc3d_kappa_file(package,name):
    """
    The RADMC-3D radiative transfer package
      http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/
    can perform dust continuum radiative transfer for diagnostic purposes.
    It is designed for astronomical applications. The code
    needs the opacities in a particular form. This subroutine 
    writes the opacities out in that form. It will write it to
    the file dustkappa_<name>.inp. This is the simpler version of
    the opacity files, containing only kabs, kscat, gscat as a function
    of wavelength.
    """
    filename = 'dustkappa_'+name+'.inp'
    with open(filename,'w') as f:
        f.write('3\n')  # Format number
        f.write('%d\n'%(package['lamcm'].size))
        f.write('\n')
        for i in range(package['lamcm'].size):
            f.write('%13.6e %13.6e %13.6e %13.6e\n'%(package['lamcm'][i]*1e4,
                                                     package['kabs'][i],
                                                     package['kscat'][i],
                                                     package['gscat'][i]))
        f.write('\n')
