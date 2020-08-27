!=======================================================================
!            MAKING OPACITY TABLE USING BOHREN-HUFFMAN PROGRAM 
!                        ADAPTED BY B.T. DRAINE
!           MADE INTO F90 AND INTO CURRENT FORM BY C.P. DULLEMOND
!
!***********************************************************************
! COMMENTS FROM ORIGINAL CODE:
! Program to interactively call Bohren-Huffman Mie theory program
!
! CALLBHMIE will interactively prompt for:
! 1. refractive index of surrounding medium
! 2. either refractive index or dielectric constant of sphere
! 3. radius of sphere
! 4. wavelength (in vacuo)
! 5. number of angles at which to calculate scattering intensities
!
! CALLBHMIE will return:
! 1. Q_ext, Q_abs, Q_sca, g, Q_back
! 2. If NANG>0, then will also return scattering matrix elements
!    S_11, S_33, S_34, and POL
!
! Adapted by B.T.Draine, Princeton Univ. Obs.
!***********************************************************************
!=======================================================================
program bhmakeopac
  integer, parameter :: MXNANG=1000
  integer :: IREADEP,J,NAN,NANG,NANG0
  real :: AJ,ANG,DANG,GSCA,PI,POL
  real :: QABS,QBACK,QEXT,QSCA,RAD,REFMED
  real :: S11,S12,S33,S34,WAVEL,X
  complex :: REFREL,CXEPS,S1(2*MXNANG-1),S2(2*MXNANG-1)
  real, allocatable :: lambda_cm(:),optcnst_n(:),optcnst_k(:)
  real :: kappa_abs,kappa_sca,kappa_g
  real :: zscat(2*MXNANG-1,1:6)
  integer :: nlam,ilam
  real :: agrain_cm,xigrain,dum(1:3),siggeom,mgrain
  character*160 :: filename,material
  logical :: notfinished
  PI=4.E0*ATAN(1.E0)
  !
  ! Defaults
  !
  REFMED = 1.d0
  !
  ! Open parameter file
  !
  open(unit=1,file='param.inp')
  read(1,*) material
  read(1,*) agrain_cm
  read(1,*) xigrain
  read(1,*) nang0
  close(1)
  filename = trim(material)//".lnk"
  !
  ! Open optical constants file
  !
  nlam = 0
  notfinished = .true.
  open(unit=1,file=filename)
  do while(notfinished)
     read(1,*,end=20) dum
     nlam = nlam + 1
  enddo
20 continue
  close(1)
  allocate(lambda_cm(nlam),optcnst_n(nlam),optcnst_k(nlam))
  open(unit=1,file=filename)
  do ilam=1,nlam
     read(1,*) lambda_cm(ilam),optcnst_n(ilam),optcnst_k(ilam)
  enddo
  close(1)
  lambda_cm = lambda_cm * 1e-4
  !
  ! Compute geometric cross section
  !
  siggeom = pi*agrain_cm**2
  !
  ! Compute mass of grain
  !
  mgrain = (4.d0*pi/3.d0)*xigrain*agrain_cm**3
  !
  ! Open the output file
  !
  filename = 'dustkappa_'//trim(material)//'.inp'
  open(unit=1,file=filename)
  write(1,*) 3  ! Format number
  write(1,*) nlam
  !
  ! Now do the loop over wavelengths
  !
  do ilam=1,nlam
     !
     ! Prepare the parameters for BHMie
     !
     ! The complex index of refraction
     !
     refrel = cmplx(optcnst_n(ilam),optcnst_k(ilam))/refmed
     !
     ! Radius of the grain in cm
     !
     rad = agrain_cm
     !
     ! Wavelength in cm
     !
     wavel = lambda_cm(ilam)
     !
     ! Number of angles
     !
     IF(NANG0.GT.MXNANG) STOP'***Error: NANG > MXNANG'
     NANG=NANG0
     IF(NANG0.LT.2)NANG=2
     !
     ! Compute the dimensionless grain size size
     !
     X=2.E0*PI*RAD*REFMED/WAVEL
     !
     ! NANG=number of angles between 0 and 90 degrees (incl. 0 and 90)
     ! Scattering matrix elements are calculated for 2*NANG-1 angles
     ! including 0, 90, and 180 degrees.
     !
     IF(NANG.GT.1)DANG=0.5E0*PI/FLOAT(NANG-1)
     !
     ! Call BHMie
     !
     CALL BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA)
     QABS=QEXT-QSCA
     !
     ! Put results into array
     !
     kappa_abs = qabs*siggeom/mgrain
     kappa_sca = qsca*siggeom/mgrain
     kappa_g   = GSCA
     write(1,*) lambda_cm(ilam)*1e4,kappa_abs,kappa_sca,kappa_g
     !
     ! Also store the Z matrix elements
     !
!     IF(NANG0.GT.1)THEN
!        NAN=2*NANG-1
!        DO J=1,NAN
!           AJ=J
!           S11=0.5E0*CABS(S2(J))*CABS(S2(J))
!           S11=S11+0.5E0*CABS(S1(J))*CABS(S1(J))
!           S12=0.5E0*CABS(S2(J))*CABS(S2(J))
!           S12=S12-0.5E0*CABS(S1(J))*CABS(S1(J))
!           POL=-S12/S11
!           S33=REAL(S2(J)*CONJG(S1(J)))
!           S34=AIMAG(S2(J)*CONJG(S1(J)))
!           ANG=DANG*(AJ-1.E0)*180.E0/PI
!           WRITE(7,6075)ANG,S11,POL,S33,S34
!           WRITE(*,6075)ANG,S11,POL,S33,S34
!        enddo
!     ENDIF
  enddo
  close(1)
  !
  ! Deallocate stuff
  !
  deallocate(lambda_cm,optcnst_n,optcnst_k)
end program bhmakeopac
