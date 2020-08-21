"""This module contains functions for Mie-scattering and to write dust opacity files
"""
from __future__ import absolute_import
from __future__ import print_function

import traceback
try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use the python module of RADMC-3D you need to install Numpy')
    print(traceback.format_exc())

import warnings

from scipy.interpolate import interp1d

try:
    from ._bhmie import bhmie
except ImportError:
    print("Fast (Fortran90) Mie-scattering module could not be imported. Falling back to the slower Python version.")

    def bhmie(x=None, refrel=None, theta=None):
        """
        The famous Bohren and Huffman Mie scattering code.
        This version was ported to Python from the f77 code from Bruce
        Draine, which can be downloaded from:
          https://www.astro.princeton.edu/~draine/scattering.html
        The code originates from the book by Bohren & Huffman (1983) on
        "Absorption and Scattering of Light by Small Particles".
        This python version was created by Cornelis Dullemond,
        February 2017.

        Parameters
        ----------

        x           : ndarray
                      Size parameter (2*pi*radius_grain/lambda)

        refrel      : complex
                      Complex index of refraction (example: 1.5 + 0.01*1j)

        theta       : ndarray
                      Scattering angles between 0 and 180 deg

        Returns
        -------
        A list with the following elements:

            * [0] S1    : ndarray
                          Complex phase function S1 (E perpendicular to scattering plane) as a function of theta

            * [1] S2    : ndarray
                          Complex phase function S1 (E parallel to scattering plane) as a function of theta

            * [2] Qext  : ndarray
                          Efficiency factor for extinction (C_ext/pi*a**2)

            * [3] Qsca  : ndarray
                          Efficiency factor for scatterin(C_sca/pi*a**2)

            * [4] Qback : ndarray
                          Backscattering efficiency ((dC_sca/domega)/pi*a**2 )
                          Note, this is (1/4*pi) smaller than the "radar backscattering efficiency" - see Bohren &
                          Huffman 1983 pp. 120-123]

            * [5] gsca  :  <cos(theta)> for scattering
        """

        if x is None:
            msg = 'Unknown scattering parameter x.'
            raise ValueError(msg)

        if x.size == 0:
            msg = 'Scattering parameter x has zero elements'
            raise ValueError(msg)

        if refrel is None:
            msg = 'Unknown refractive index refrel.'
            raise ValueError(msg)

        if theta is None:
            msg = 'Unknown scattering angle theta.'
            raise ValueError(msg)

        if theta.shape[0] == 0:
            msg = 'Scattering angle theta has zero elements'
            raise ValueError(msg)

        #
        # First check that the theta array goes from 0 to 180 or
        # 180 to 0, and store which is 0 and which is 180
        #

        if theta[0] != 0:
            msg = "First element of scattering angle array is not 0. Scattering angle grid must extend from 0 to 180 " \
                  "degrees."
            raise ValueError(msg)
        if theta[-1] != 180:
            msg = "Last element of scattering angle array is not 180. Scattering angle grid must extend from 0 to 180 " \
                  "degrees."
            raise ValueError(msg)

        nang = len(theta)
        #
        # Allocate the complex phase functions with double precision
        #
        S1 = np.zeros(nang, dtype=np.complex128)
        S2 = np.zeros(nang, dtype=np.complex128)
        #
        # Allocate and initialize arrays for the series expansion iteration
        #
        pi = np.zeros(nang, dtype=np.float64)
        pi0 = np.zeros(nang, dtype=np.float64)
        pi1 = np.zeros(nang, dtype=np.float64) + 1.0
        tau = np.zeros(nang, dtype=np.float64)
        #
        # Compute an alternative to x
        #
        y = x * refrel
        #
        # Determine at which n=nstop to terminate the series expansion
        #
        xstop = x + 4 * x ** 0.3333 + 2.0
        nstop = int(np.floor(xstop))
        #
        # Determine the start of the logarithmic derivatives iteration
        #
        nmx = int(np.floor(np.max([xstop, abs(y)])) + 15)
        #
        # Compute the mu = cos(theta*pi/180.) for all scattering angles
        #
        mu = np.cos(theta * np.pi / 180.)
        #
        # Now calculate the logarithmic derivative dlog by downward recurrence
        # beginning with initial value 0.+0j at nmx-1
        #
        dlog = np.zeros(nmx, dtype=np.complex128)
        for n in range(nmx - 1):
            en = float(nmx - n)
            dlog[nmx - n - 2] = en / y - 1.0 / (dlog[nmx - n - 1] + en / y)
        #
        # In preparation for the series expansion, reset some variables
        #
        psi0 = np.cos(x)
        psi1 = np.sin(x)
        chi0 = -np.sin(x)
        chi1 = np.cos(x)
        xi1 = psi1 - chi1 * 1j
        p = -1.0
        Qsca = 0.0
        gsca = 0.0
        an = 0j
        bn = 0j
        #
        # Riccati-Bessel functions with real argument x
        # calculated by upward recurrence. This is where the
        # series expansion is done
        #
        for n in range(nstop):
            #
            # Basic calculation of the iteration
            #
            en = float(n + 1)
            fn = (2 * en + 1.0) / (en * (en + 1.0))
            psi = (2 * en - 1.0) * psi1 / x - psi0
            chi = (2 * en - 1.0) * chi1 / x - chi0
            xi = psi - chi * 1j
            an1 = an
            bn1 = bn
            dum = dlog[n] / refrel + en / x
            an = (dum * psi - psi1) / (dum * xi - xi1)
            dum = dlog[n] * refrel + en / x
            bn = (dum * psi - psi1) / (dum * xi - xi1)
            #
            # Add contributions to Qsca and gsca
            #
            Qsca += (2 * en + 1.0) * (abs(an) ** 2 + abs(bn) ** 2)
            dum = (2 * en + 1.0) / (en * (en + 1.0))
            gsca += dum * (an.real * bn.real + an.imag * bn.imag)
            dum = (en - 1.0) * (en + 1.0) / en
            gsca += dum * (an1.real * an.real + an1.imag * an.imag +
                           bn1.real * bn.real + bn1.imag * bn.imag)
            #
            # Now contribute to scattering intensity pattern as a function of angle
            #
            pi[:] = pi1[:]
            tau[:] = en * np.abs(mu[:]) * pi[:] - (en + 1.0) * pi0[:]
            #
            # For mu>=0
            #
            idx = mu >= 0
            S1[idx] += fn * (an * pi[idx] + bn * tau[idx])
            S2[idx] += fn * (an * tau[idx] + bn * pi[idx])
            #
            # For mu<0
            #
            p = -p
            idx = mu < 0
            S1[idx] += fn * p * (an * pi[idx] - bn * tau[idx])
            S2[idx] += fn * p * (bn * pi[idx] - an * tau[idx])
            #
            # Now prepare for the next iteration
            #
            psi0 = psi1
            psi1 = psi
            chi0 = chi1
            chi1 = chi
            xi1 = psi1 - chi1 * 1j
            pi1[:] = ((2 * en + 1.0) * np.abs(mu[:]) * pi[:] - (en + 1.0) * pi0[:]) / en
            pi0[:] = pi[:]
        #
        # Now do the final calculations
        #
        gsca = 2 * gsca / Qsca
        Qsca = (2.0 / (x * x)) * Qsca
        Qext = (4.0 / (x * x)) * S1[0].real
        Qback = (abs(S1[-1]) / x)**2 / np.pi
        Qabs = Qext - Qsca
        #
        # Return results
        #
        return [S1, S2, Qext, Qabs, Qsca, Qback, gsca]


def compute_opac_mie(fname='', matdens=None, agraincm=None, lamcm=None,
                     theta=None, logawidth=None, wfact=3.0, na=20,
                     chopforward=0.0, errtol=0.01, verbose=False,
                     extrapolate=False):
    """
    Compute dust opacity with Mie theory based on the optical constants
    in the optconst_file. Optionally also the scattering phase function
    in terms of the Mueller matrix elements can be computed. To smear out
    the resonances that appear due to the perfect sphere shape, you can
    optionally smear out the grain size distribution a bit with setting
    the width of a Gaussian grain size distribution.

    Parameters
    ----------
    fname       : str
                  File name of the optical constants file. This file
                  should contain three columns: first the wavelength
                  in micron, then the n-coefficient and then the
                  k-coefficient. See Jena optical constants database:
                  http://www.astro.uni-jena.de/Laboratory/Database/databases.html

    matdens     : float
                  Material density in g/cm^3

    agraincm    : float
                  Grain radius in cm

    lamcm       : ndarray
                  Wavelength grid in cm

    theta       : ndarray, optional
                  Angular grid (a numpy array) between 0 and 180
                  which are the scattering angle sampling points at
                  which the scattering phase function is computed.

    logawidth   : float, optional
                 If set, the size agrain will instead be a
                 sample of sizes around agrain. This helps to smooth out
                 the strong wiggles in the phase function and opacity
                 of spheres at an exact size. Since in Nature it rarely
                 happens that grains all have exactly the same size, this
                 is quite natural. The value of logawidth sets the width
                 of the Gauss in ln(agrain), so for logawidth<<1 this
                 give a real width of logawidth*agraincm.

    wfact       : float
                  Grid width of na sampling points in units
                  of logawidth. The Gauss distribution of grain sizes is
                  cut off at agrain * exp(wfact*logawidth) and
                  agrain * exp(-wfact*logawidth). Default = 3


    na          : int
                  Number of size sampling points (if logawidth set, default=20)

    chopforward : float
                  If >0 this gives the angle (in degrees from forward)
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


    errtol      : float
                  Tolerance of the relative difference between kscat
                  and the integral over the zscat Z11 element over angle.
                  If this tolerance is exceeded, a warning is given.

    verbose     : bool
                  If set to True, the code will give some feedback so
                  that one knows what it is doing if it becomes slow.

    extrapolate : bool
                  If set to True, then if the wavelength grid lamcm goes
                  out of the range of the wavelength grid of the
                  optical constants file, then it will make a suitable
                  extrapolation: keeping the optical constants constant
                  for lamcm < minimum, and extrapolating log-log for
                  lamcm > maximum.

    Returns
    -------
    A dictionary with the following keys:

        * kabs          : ndarray
                          Absorption opacity kappa_abs_nu (a numpy array) in
                          units of cm^2/gram

        * ksca          : ndarray
                          Scattering opacity kappa_abs_nu (a numpy array) in
                          units of cm^2/gram

        * gsca          : ndarray
                          The <cos(theta)> g-factor of scattering

        * theta         : ndarray (optional, only if theta is given at input)
                          The theta grid itself (just a copy of what was given)

        * zscat         : ndarray (optional, only if theta is given at input)
                          The components of the scattering Mueller matrix
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

        * kscat_from_z11 : ndarray  (optional, only if theta is given at input)
                           The kscat computed from the (above mentioned)
                           integral of Z11 over all angles. This should be
                           nearly identical to kscat if the angular grid
                           is sufficiently fine. If there are strong
                           differences, this is an indication that the
                           angular gridding (the theta grid) is not fine
                           enough. But you should have then automatically
                           gotten a warning message as well (see errtol).

        * wavmic        : ndarray (optional, only if extrapolate is set to True)
                          The original wavelength grid from the optical constants file,
                          with possibly an added extrapolated

        * ncoef         : ndarray (optional, only if extrapolate is set to True)
                          The optical constant n at that grid

        * kcoef         : ndarray (optional, only if extrapolate is set to True)
                          The optical constant k at that grid

        * agr           : ndarray (optional, only if logawidth is not None)
                          Grain sizes

        * wgt           : ndarray (optional, only if logawidth is not None)
                          The averaging weights of these grain (not the masses!)
                          The sum of wgt.sum() must be 1.

        * zscat_nochop  : ndarray (optional, only if chopforward > 0)
                          The zscat before the forward scattering was chopped off

        * kscat_nochop  : ndarray (optional, only if chopforward > 0)
                          The kscat originally from the bhmie code
    """
    #
    # Load the optical constants
    #
    if matdens is None:
        msg = "Unknown material density matdens"
        raise ValueError(msg)

    if agraincm is None:
        msg = "Unknown grain size agraincm"
        raise ValueError(msg)

    if lamcm is None:
        msg = "Unknown wavelength grid lamcm"
        raise ValueError(msg)

    if theta is None:
        angles = np.array([0., 90., 180.])  # Minimalistic angular s
        if chopforward != 0.:
            warnings.warn("Chopping disabled. Chopping is only possible if theta grid is given. ", RuntimeWarning)
    else:
        angles = theta

    #
    # Check that the theta array goes from 0 to 180 or
    # 180 to 0, and store which is 0 and which is 180
    #
    if angles[0] != 0:
        msg = "First element of the angular grid array is not 0. Scattering angle grid must extend from 0 to 180 " \
              "degrees."
        raise ValueError(msg)
    if angles[-1] != 180:
        msg = "Last element of the angular grid array is not 180. Scattering angle grid must extend from 0 to 180 " \
              "degrees."
        raise ValueError(msg)

    nang = angles.shape[0]

    #
    # Load the optical constants
    #
    data = np.loadtxt(fname)
    wavmic, ncoef, kcoef = data.T

    if wavmic.size <= 1:
        msg = "Optical constants file must have at least two rows with two different wavelengths"
        raise ValueError(msg)

    if wavmic[1] == wavmic[0]:
        msg = "Optical constants file must have at least two rows with two different wavelengths"
        raise ValueError(msg)

    #
    # Check range, and if needed and requested, extrapolate the
    # optical constants to longer or shorter wavelengths
    #
    if extrapolate:
        wmin = np.min(lamcm)*1e4 * 0.999
        wmax = np.max(lamcm)*1e4 * 1.001
        if wmin < np.min(wavmic):
            if wavmic[0] < wavmic[1]:
                ncoef = np.append([ncoef[0]], ncoef)
                kcoef = np.append([kcoef[0]], kcoef)
                wavmic = np.append([wmin], wavmic)
            else:
                ncoef = np.append(ncoef, [ncoef[-1]])
                kcoef = np.append(kcoef, [kcoef[-1]])
                wavmic = np.append(wavmic, [wmin])
        if wmax > np.max(wavmic):
            if wavmic[0] < wavmic[1]:
                ncoef = np.append(ncoef, [ncoef[-1] * np.exp((np.log(wmax) - np.log(wavmic[-1])) *
                                                             (np.log(ncoef[-1]) - np.log(ncoef[-2])) /
                                                             (np.log(wavmic[-1]) - np.log(wavmic[-2])))])
                kcoef = np.append(kcoef, [kcoef[-1]*np.exp((np.log(wmax) - np.log(wavmic[-1])) *
                                                             (np.log(kcoef[-1]) - np.log(kcoef[-2])) /
                                                             (np.log(wavmic[-1]) - np.log(wavmic[-2])))])
                wavmic = np.append(wavmic, [wmax])
            else:
                ncoef = np.append(ncoef, [ncoef[0]*np.exp((np.log(wmax)-np.log(wavmic[0])) *
                                                             (np.log(ncoef[0]) - np.log(ncoef[1])) /
                                                             (np.log(wavmic[0]) - np.log(wavmic[1])))])
                kcoef = np.append(kcoef, [kcoef[0]*np.exp((np.log(wmax) - np.log(wavmic[0])) *
                                                             (np.log(kcoef[0]) - np.log(kcoef[1])) /
                                                             (np.log(wavmic[0]) - np.log(wavmic[1])))])
                wavmic = np.append([wmax], wavmic)
    else:
        if lamcm.min() < wavmic.min()*1e-4:
            print(lamcm.min(), wavmic.min()*1e4)
            raise ValueError("Wavelength range out of range of the optical constants file")

        if lamcm.max() > wavmic.max()*1e-4:
            print(lamcm.min(), wavmic.min()*1e4)
            raise ValueError("Wavelength range out of range of the optical constants file")

    # Interpolate
    # Note: Must be within range, otherwise stop
    #
    f = interp1d(np.log(wavmic*1e-4), np.log(ncoef))
    ncoefi = np.exp(f(np.log(lamcm)))
    f = interp1d(np.log(wavmic*1e-4), np.log(kcoef))
    kcoefi = np.exp(f(np.log(lamcm)))
    #
    # Make the complex index of refraction
    #
    refidx = ncoefi + kcoefi*1j
    #
    # Make a size distribution for the grains
    # If width is not set, then take just one size
    #
    if logawidth is None:
        agr = np.array([agraincm])
        wgt = np.array([1.0])
    else:
        if logawidth != 0.0:
            agr = np.exp(np.linspace(np.log(agraincm) - wfact * logawidth, np.log(agraincm) + wfact * logawidth, na))
            wgt = np.exp(-0.5*((np.log(agr / agraincm)) / logawidth)**2)
            wgt = wgt / wgt.sum()
        else:
            agr = np.array([agraincm])
            wgt = np.array([1.0])
    #
    # Get the true number of grain sizes
    #
    nagr = agr.size
    #
    # Compute the geometric cross sections
    #
    siggeom = np.pi*agr*agr
    #
    # Compute the mass of the grain
    #
    mgrain = (4*np.pi/3.0)*matdens*agr*agr*agr
    #
    # Now prepare arrays
    #
    nlam = lamcm.size
    kabs = np.zeros(nlam)
    kscat = np.zeros(nlam)
    gscat = np.zeros(nlam)
    if theta is not None:
        zscat = np.zeros((nlam, nang, 6))
        S11 = np.zeros(nang)
        S12 = np.zeros(nang)
        S33 = np.zeros(nang)
        S34 = np.zeros(nang)
        if chopforward > 0:
            zscat_nochop = np.zeros((nlam, nang, 6))
            kscat_nochop = np.zeros(nlam)
    #
    # Set error flag to False
    #
    error = False
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
            print("Doing wavelength %13.6e cm" % lamcm[i])
        #
        # Now loop over the grain sizes
        #
        for l in range(nagr):
            #
            # Message
            #
            if verbose and nagr > 1:
                print("...Doing grain size %13.6e cm" % agr[l])
            #
            # Compute x
            #
            x = 2*np.pi*agr[l]/lamcm[i]

            #
            # Call the bhmie code
            #
            S1, S2, Qext, Qabs, Qsca, Qback, gsca = bhmie(x, refidx[i], angles)
            #
            # Add results to the averaging over the size distribution
            #
            kabs[i] += wgt[l] * Qabs*siggeom[l] / mgrain[l]
            kscat[i] += wgt[l] * Qsca*siggeom[l] / mgrain[l]
            gscat[i] += wgt[l] * gsca
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
                factor = (lamcm[i]/(2*np.pi))**2/mgrain[l]
                #
                # Compute the scattering Mueller matrix elements at each angle
                #
                S11[:] = 0.5 * (np.abs(S2[:])**2 + np.abs(S1[:])**2)
                S12[:] = 0.5 * (np.abs(S2[:])**2 - np.abs(S1[:])**2)
                S33[:] = np.real(S2[:] * np.conj(S1[:]))
                S34[:] = np.imag(S2[:] * np.conj(S1[:]))
                zscat[i, :, 0] += wgt[l] * S11[:] * factor
                zscat[i, :, 1] += wgt[l] * S12[:] * factor
                zscat[i, :, 2] += wgt[l] * S11[:] * factor
                zscat[i, :, 3] += wgt[l] * S33[:] * factor
                zscat[i, :, 4] += wgt[l] * S34[:] * factor
                zscat[i, :, 5] += wgt[l] * S33[:] * factor
        #
        # If possible, do a check if the integral over zscat is consistent
        # with kscat
        #
        if theta is not None:
            mu = np.cos(angles * np.pi / 180.)
            dmu = np.abs(mu[1:nang] - mu[0:nang-1])
            zav = 0.5 * (zscat[i, 1:nang, 0] + zscat[i, 0:nang-1, 0])
            dum = 0.5 * zav * dmu
            kscat_from_z11[i] = dum.sum() * 4 * np.pi
            err = abs(kscat_from_z11[i]/kscat[i]-1.0)
            if err > errtol:
                error = True
                errmax = max(err, errmax)
        #
        # If the chopforward angle is set >0, then we will remove
        # excessive forward scattering from the opacity. The reasoning
        # is that extreme forward scattering is, in most cases, equivalent
        # to no scattering at all.
        #
        if chopforward > 0:
            iang = np.where(angles < chopforward)
            if angles[0] == 0.0:
                iiang = np.max(iang)+1
            else:
                iiang = np.min(iang)-1
            zscat_nochop[i, :, :] = zscat[i, :, :]  # Backup
            kscat_nochop[i] = kscat[i]      # Backup
            zscat[i, iang, 0] = zscat[i, iiang, 0]
            zscat[i, iang, 1] = zscat[i, iiang, 1]
            zscat[i, iang, 2] = zscat[i, iiang, 2]
            zscat[i, iang, 3] = zscat[i, iiang, 3]
            zscat[i, iang, 4] = zscat[i, iiang, 4]
            zscat[i, iang, 5] = zscat[i, iiang, 5]
            mu = np.cos(angles * np.pi / 180.)
            dmu = np.abs(mu[1:nang] - mu[0:nang-1])
            zav = 0.5 * (zscat[i, 1:nang, 0] + zscat[i, 0:nang-1, 0])
            dum = 0.5 * zav * dmu
            kscat[i] = dum.sum() * 4 * np.pi

            zav = 0.5 * (zscat[i, 1:nang, 0] * mu[1:] + zscat[i, 0:nang-1, 0] * mu[:-1])
            dum = 0.5 * zav * dmu
            gscat[i] = dum.sum() * 4 * np.pi / kscat[i]

    #
    # If error found, then warn (Then shouldn't it be called a warning? If it's a true error
    #  shouldn't we stop the execution and raise an exception?)
    #
    if error:
        msg = " Angular integral of Z11 is not equal to kscat at all wavelength. \n"
        msg += "Maximum error = %13.6e" % errmax
        if chopforward > 0:
            msg += "But I am using chopforward to remove strong forward scattering, and then renormalized kapscat."
        warnings.warn(msg, RuntimeWarning)
    #
    # Now return what we computed in a dictionary
    #
    package = {"lamcm": lamcm, "kabs": kabs, "kscat": kscat,
               "gscat": gscat, "matdens": matdens, "agraincm": agraincm}
    if theta is not None:
        package["zscat"] = np.copy(zscat)
        package["theta"] = np.copy(angles)
        package["kscat_from_z11"] = np.copy(kscat_from_z11)
    if extrapolate:
        package["wavmic"] = np.copy(wavmic)
        package["ncoef"] = np.copy(ncoef)
        package["kcoef"] = np.copy(kcoef)
    if nagr > 1:
        package["agr"] = np.copy(agr)
        package["wgt"] = np.copy(wgt)
        package["wfact"] = wfact
        package["logawidth"] = logawidth
    if chopforward > 0:
        package["zscat_nochop"] = np.copy(zscat_nochop)
        package["kscat_nochop"] = np.copy(kscat_nochop)

    return package


def write_radmc3d_scatmat_file(package=None, name=None, comment=False):
    """
    The RADMC-3D radiative transfer package
      http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/
    can perform dust continuum radiative transfer for diagnostic purposes.
    It is designed for astronomical applications. The code
    needs the opacities in a particular form. This subroutine
    writes the opacities out in that form. It will write it to
    the file dustkapscatmat_<name>.inp.
    """
    if package is None:
        raise ValueError("Unknown package. No data to be written.")

    if name is None:
        raise ValueError("Unkonwn name. Without a file name tag for dustkapscatmat_NAME.inp must be given.")

    filename = 'dustkapscatmat_' + name + '.inp'
    with open(filename, 'w') as f:
        if comment:
            f.write('# Opacity and scattering matrix file for ' + name + '\n')
            f.write('# Please do not forget to cite in your publications the original paper of these optical ' +
                    'constant measurements\n')
            f.write('# Made with the makedustopac.py code by Cornelis Dullemond\n')
            f.write('# using the bhmie.py Mie code of Bohren and Huffman (python version by Cornelis Dullemond, ' +
                    'from original bhmie.f code by Bruce Draine)\n')
            f.write('# Grain size = %13.6e cm\n' % (package['agraincm']))
            f.write('# Material density = %6.3f g/cm^3\n' % (package['matdens']))
        f.write('1\n')  # Format number
        f.write('%d\n' % package['lamcm'].size)
        f.write('%d\n' % package['theta'].size)
        f.write('\n')
        for i in range(package['lamcm'].size):
            f.write('%13.6e %13.6e %13.6e %13.6e\n' % (package['lamcm'][i] * 1e4, package['kabs'][i],
                                                       package['kscat'][i], package['gscat'][i]))
        f.write('\n')
        for j in range(package['theta'].size):
            f.write('%13.6e\n' % (package['theta'][j]))
        f.write('\n')
        for i in range(package['lamcm'].size):
            for j in range(package['theta'].size):
                f.write('%13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n' %
                        (package['zscat'][i, j, 0], package['zscat'][i, j, 1],
                         package['zscat'][i, j, 2], package['zscat'][i, j, 3],
                         package['zscat'][i, j, 4], package['zscat'][i, j, 5]))
        f.write('\n')


def write_radmc3d_kappa_file(package=None, name=None):
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
    if package is None:
        raise ValueError("Unknown package. No data to be written.")

    if name is None:
        raise ValueError("Unkonwn name. Without a file name tag for dustkappa_NAME.inp must be given.")

    filename = 'dustkappa_'+name+'.inp'
    with open(filename, 'w') as f:
        f.write('3\n')  # Format number
        f.write('%d\n' % package['lamcm'].size)
        # f.write('\n')
        for i in range(package['lamcm'].size):
            f.write('%13.6e %13.6e %13.6e %13.6e\n' % (package['lamcm'][i] * 1e4,
                                                     package['kabs'][i],
                                                     package['kscat'][i],
                                                     package['gscat'][i]))
        f.write('\n')


