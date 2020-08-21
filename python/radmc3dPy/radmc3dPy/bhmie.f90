!!
!! The famous Bohren and Huffman Mie scattering code.
!! This version was ported to Fortran90/95 from the Python version by Attila Juhasz.
!! The Python version was ported from the original Fortran 77 version by Cornelis Dullemond.
!! The original Fortran 77 code was written by Bruce Drain and can be downloaded from:
!! https://www.astro.princeton.edu/~draine/scattering.html
!! The code originates from the book by Bohren & Huffman (1983) on
!! "Absorption and Scattering of Light by Small Particles".
!!
SUBROUTINE bhmie(x, refrel, theta, nang, S1, S2, Qext, Qabs, Qsca, Qback, gsca)
    implicit none
    
    real(kind=8), intent(in)        :: x
    complex(kind=8), intent(in)     :: refrel
    integer, intent(in)             :: nang
    real(kind=8), intent(in)        :: theta(nang)

    complex(kind=8), intent(out)    :: S1(nang)
    complex(kind=8), intent(out)    :: S2(nang)
    real(kind=8), intent(out)       :: Qext
    real(kind=8), intent(out)       :: Qabs
    real(kind=8), intent(out)       :: Qsca
    real(kind=8), intent(out)       :: Qback
    real(kind=8), intent(out)       :: gsca
    !! ---------------------------------------
    !! Local variables
    !! ---------------------------------------
    real(kind=8)                    :: math_pi=3.14159265358979323846264338328
    integer                         :: i, iang
    integer                         :: nstop
    real(kind=8)                    :: xstop
    complex(kind=8)                 :: y
    integer                         :: nmx
    real(kind=8)                    :: mu(nang)
    complex(kind=8), allocatable    :: log_der(:)
   
    real(kind=8)                    :: chi0, chi1, chi, p
    complex(kind=8)                 :: xi1, xi
    complex(kind=8)                 :: an, bn, an1, bn1
    real(kind=8)                    :: pi(nang), pi1(nang), pi0(nang), tau(nang)
    real(kind=8)                    :: fn
    
    real(kind=8)                    :: psi0, psi1, psi
    real(kind=8)                    :: en
    !!
    !! Reset all output arrays
    !!
    do i=1, nang
        S1(i)    = (0d0, 0d0)
        S2(i)    = (0d0, 0d0)
        pi(i)    = 0d0
        pi0(i)   = 0d0
        pi1(i)   = 1d0
        tau(i)   = 0d0
    end do
    Qext  = 0e0
    Qabs  = 0e0
    Qsca  = 0e0
    Qback = 0e0
    gsca  = 0e0

    !!
    !! Compute an alternative to x
    !!
    y = x * refrel
    !!
    !! Determine at which n=nstop to terminate the series expansion
    !!
    xstop = x + 4.0 * x**0.3333 + 2.0
    nstop = int(floor(xstop))
    !!
    !! Determine the start of the logarithmic derivatives iteration
    !!
    nmx = int(floor(max(xstop, abs(y))) + 15)
    !!
    !! Compute the mu = cos(theta*pi/180.) for all scattering angles
    !!
    do iang=1, nang
        mu(iang) = cos(theta(iang) * math_pi / 180.0)
        if (abs(theta(iang)/90.0 - 1.0) < 1e-8) then
            mu(iang) = 0.0
        end if
    end do
    mu(1) = 1.0d0
    mu(nang) = -1.0d0
    !!
    !! Now calculate the logarithmic derivative dlog by downward recurrence
    !! beginning with initial value 0.+0j at nmx-1
    !!
    allocate(log_der(nmx))
    log_der(nmx) = (0e0, 0e0)
    do i=1, nmx-1
        en = dble(nmx - i + 1)
        log_der(nmx - i) = (en / y) - 1.0 / (log_der(nmx - i + 1) + en / y) 
    end do
    !!
    !! In preparation for the series expansion, reset some variables
    !!
    psi0 = cos(x)
    psi1 = sin(x)
    chi0 = -sin(x)
    chi1 = cos(x)
    xi1 = cmplx(psi1, -chi1, kind=8)
    p = -1.0
    an = (0e0, 0e0) 
    bn = (0e0, 0e0)
    !!
    !! Riccati-Bessel functions with real argument x
    !! calculated by upward recurrence. This is where the
    !! series expansion is done
    !!
    do i=1, nstop
        !!
        !! Basic calculation of the iteration
        !!
        en = dble(i)
        fn = (2 * en + 1.0) / (en * (en + 1.0))
        psi = (2 * en - 1.0) * psi1 / x - psi0
        chi = (2 * en - 1.0) * chi1 / x - chi0
        xi = cmplx(psi, -chi, kind=8)
        an1 = an
        bn1 = bn
        an = log_der(i) / refrel + en / x
        an = (an * psi - psi1) / (an * xi - xi1)
        bn = log_der(i) * refrel + en / x
        bn = (bn * psi - psi1) / (bn * xi - xi1)

        !!
        !! Add contributions to Qsca and gsca
        !!
        Qsca = Qsca + (2d0 * en + 1d0) * (abs(an)**2 + abs(bn)**2)
        gsca = gsca + (2d0 * en + 1d0) / (en * (en + 1d0)) * &
                (real(an) * real(bn) + aimag(an) * aimag(bn))
        gsca = gsca + (en - 1.0) * (en + 1.0) / en * &
                (real(an1) * real(an) + aimag(an1) * aimag(an) + &
                 real(bn1) * real(bn) + aimag(bn1) * aimag(bn))
        !!
        !! Now contribute to scattering intensity pattern as a function of angle
        !!
        pi = pi1
        tau = en * abs(mu) * pi - (en + 1.0) * pi0
        p = -p
        do iang=1, nang
            !!
            !! For mu>=0
            !!
            if (mu(iang).ge.0) then
                S1(iang) = S1(iang) + fn * (an * pi(iang) + bn * tau(iang))
                S2(iang) = S2(iang) + fn * (an * tau(iang) + bn * pi(iang))
            !!
            !! For mu<0
            !!
            else
                S1(iang) = S1(iang) + fn * p * (an * pi(iang) - bn * tau(iang))
                S2(iang) = S2(iang) + fn * p * (bn * pi(iang) - an * tau(iang))
            endif
        end do
        !!
        !! Now prepare for the next iteration
        !!
        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
        xi1  = cmplx(psi1, -chi1, kind=8)
        pi1  = ((2d0 * en + 1d0) * abs(mu) * pi - (en + 1d0) * pi0) / en
        pi0  = pi
    end do

    !!
    !! Now do the final calculations
    !!
    gsca = 2 * gsca / Qsca
    Qsca = (2d0 / (x * x)) * Qsca
    Qext = (4d0 / (x * x)) * real(S1(1))
    Qback = abs(S1(nang)) * abs(S1(nang)) / (x * x * math_pi)
    Qabs = Qext - Qsca
    deallocate(log_der)
     
END SUBROUTINE bhmie



