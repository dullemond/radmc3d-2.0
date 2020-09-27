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
  implicit none
  integer, parameter :: MXNANG=1000
  integer :: IREADEP,J,NAN,NANG
  doubleprecision :: DANG,PI,sum,error,errmax,summu
  real :: QABS,QBACK,QEXT,QSCA,RAD,REFMED,GSCA,POL
  real :: S11,S12,S33,S34,WAVEL,X,XMAX
  complex :: REFREL,CXEPS,S1(2*MXNANG-1),S2(2*MXNANG-1)
  doubleprecision, allocatable :: lambda_cm(:),optcnst_n(:),optcnst_k(:)
  doubleprecision, allocatable :: lambda_cm_orig(:),optcnst_n_orig(:),optcnst_k_orig(:)
  doubleprecision, allocatable :: kappa_abs(:),kappa_sca(:),kappa_g(:)
  doubleprecision, allocatable :: zscat(:,:,:),angle(:),mu(:),scalefact(:)
  doubleprecision, allocatable :: agrain_cm(:),weight(:),mgrain(:)
  integer :: nlam,ilam,nang180,nagr,ia,leng,nrcomments,icomment,iang,irescalez
  integer :: nlam_orig,ilam_orig,ilamtmp,istr,jstr,ihighest_x
  doubleprecision :: agrain_cm_mean,logawidth,xigrain,dum(1:3)
  doubleprecision :: siggeom,factor,wfact,chopforward,lam
  doubleprecision :: lam_orig_min,lam_orig_max,opt_n_min,opt_n_max,opt_k_min,opt_k_max
  doubleprecision :: slope_n,slope_k,eps,max_x,highest_x
  character*160 :: filename,material,str0,str1,wlfile,ref1,ref2,ref3
  logical :: notfinished,wlfile_exists
  PI=4.D0*ATAN(1.D0)
  !
  ! Defaults
  !
  REFMED = 1.d0
  errmax = 0.01       ! Default maximum allowed relative error
  chopforward = 0.d0  ! By default no chopping
  irescalez = 0       ! Major difference with older versions: Now by default do not rescale Z to match kappa_scat
  wlfile = ""         ! If not "", then this is the file with the list of wavelengths to be used (in micron)
  ref1 = ""
  ref2 = ""
  ref3 = ""
  max_x = 1d99
  !
  ! Open parameter file
  !
  open(unit=1,file='param_smoothed.inp')
  read(1,*) material
  read(1,*) nagr
  read(1,*) agrain_cm_mean
  read(1,*) logawidth
  read(1,*) wfact
  read(1,*) xigrain                ! Set to 0 to use material density from optical constants file header
  read(1,*) nang180                ! Nr of angles between 0 and 180 degrees
  read(1,*,end=209) errmax         ! If encountering errors above this, then stop
  read(1,*,end=209) irescalez      ! Keep this 0. Only for backward compatibility with old version (=1)
  read(1,*,end=209) chopforward    ! Forward scattering with angles less than this: apply chopping method
  read(1,*,end=209) wlfile         ! If not "", then this is the file with the list of wavelengths to be used (in micron)
  read(1,*,end=209) max_x          ! Maximum value of X=2*pi*a/lambda. Beyond that we extrapolate.
209 continue
  close(1)
  filename = trim(material)//".lnk"
  wlfile_exists = .false.
  if(wlfile.ne."") then
     inquire(file=wlfile,exist=wlfile_exists)
  endif
  !
  ! Do a check
  !
  if(nang180.lt.3) then
     write(*,*) 'This code is meant for making the full scattering matrix'
     write(*,*) 'at a set of discrete scattering angles. You must therefore'
     write(*,*) 'put nang0 (the fourth line of param.inp) to at least 3, but'
     write(*,*) 'preferably something like 90.'
     stop
  endif
  if(nang180.lt.60) then
     write(*,*) 'Warning: You have less than 60 angles between 0 and 180 deg.'
     write(*,*) 'Are you sure that this is what you want?'
  endif
  !
  ! NANG=number of angles between 0 and 90 degrees (incl. 0 and 90)
  ! Scattering matrix elements are calculated for 2*NANG-1 angles
  ! including 0, 90, and 180 degrees.
  !
  IF(NANG180.GT.2*MXNANG-1)STOP'***Error: NANG > MXNANG'
  NANG=(NANG180+1)/2
  DANG=0.5d0*PI/(NANG-1.d0)
  NAN=NANG180
  allocate(mu(nan),angle(nan))
  do j=1,nan
     angle(j)=DANG*((1.d0*J)-1.E0)*180.E0/PI
     mu(j)=cos(angle(j)*PI/180.)
  enddo
  mu(1) = 1.0
  mu(nan) = -1.0
  !
  ! Open optical constants file
  !
  open(unit=1,file=filename)
  notfinished = .true.
  nrcomments = 0
  do while(notfinished)
     read(1,'(A)') str0
     if(str0(1:1)=='#') then
        nrcomments = nrcomments + 1
        if(str0(3:3)=='@') then
           if(str0(4:12)=='reference') then
              istr = 13
              do while((str0(istr:istr)==' ').or.(str0(istr:istr)=='='))
                 istr = istr+1
              enddo
              if(ref1=="") then
                 ref1 = str0(istr:len(str0))
              elseif(ref2=="") then
                 ref2 = str0(istr:len(str0))
              elseif(ref3=="") then
                 ref3 = str0(istr:len(str0))
              endif
           elseif(str0(4:10)=='density') then
              istr = 11
              do while((str0(istr:istr).eq.' ').or.(str0(istr:istr).eq.'=').or.(istr.ge.160))
                 istr = istr+1
              enddo
              jstr = istr
              do while((str0(jstr:jstr).ne.' ').and.(jstr.lt.160))
                 jstr = jstr+1
              enddo
              jstr = jstr - 1
              if(xigrain.le.0.d0) then
                 read(str0(istr:jstr),*) xigrain
              endif
           endif
        endif
     else
        notfinished = .false.
     endif
  enddo
  close(1)
  if(xigrain.le.0.d0) then
     write(*,*) 'ERROR: material density is not specified.'
     stop
  endif
  nlam_orig = 0
  notfinished = .true.
  open(unit=1,file=filename)
  do icomment=1,nrcomments
     read(1,*) str0
  enddo
  do while(notfinished)
     read(1,*,end=20) dum
     nlam_orig = nlam_orig + 1
  enddo
20 continue
  close(1)
  allocate(lambda_cm_orig(nlam_orig),optcnst_n_orig(nlam_orig),optcnst_k_orig(nlam_orig))
  open(unit=1,file=filename)
  do icomment=1,nrcomments
     read(1,*) str0
  enddo
  do ilam=1,nlam_orig
     read(1,*) lambda_cm_orig(ilam),optcnst_n_orig(ilam),optcnst_k_orig(ilam)
  enddo
  close(1)
  lambda_cm_orig = lambda_cm_orig * 1d-4
  !
  ! Map these optical constants onto the actual wavelength grid
  !
  if(.not.wlfile_exists) then
     !
     ! No wavelength file given, so use the same wavelength grid as the optical constants file
     !
     nlam = nlam_orig
     allocate(lambda_cm(nlam),optcnst_n(nlam),optcnst_k(nlam))
     do ilam=1,nlam
        lambda_cm(ilam) = lambda_cm_orig(ilam)
        optcnst_n(ilam) = optcnst_n_orig(ilam)
        optcnst_k(ilam) = optcnst_k_orig(ilam)
     enddo
  else
     !
     ! Read the wavelength file
     !
     open(unit=1,file=wlfile)
     read(1,*) nlam
     allocate(lambda_cm(nlam),optcnst_n(nlam),optcnst_k(nlam))
     do ilam=1,nlam
        read(1,*) lambda_cm(ilam)
     enddo
     close(1)
     lambda_cm = lambda_cm * 1d-4
     !
     ! Get some information about the original grid
     !
     if(lambda_cm_orig(1).lt.lambda_cm_orig(nlam_orig)) then
        lam_orig_min = lambda_cm_orig(1)
        lam_orig_max = lambda_cm_orig(nlam_orig)
        opt_n_min    = optcnst_n_orig(1)
        opt_n_max    = optcnst_n_orig(nlam_orig)
        opt_k_min    = optcnst_k_orig(1)
        opt_k_max    = optcnst_k_orig(nlam_orig)
        slope_n      = (log(optcnst_n_orig(nlam_orig))-log(optcnst_n_orig(nlam_orig-1)))/ &
                       (log(lambda_cm_orig(nlam_orig))-log(lambda_cm_orig(nlam_orig-1)))
        slope_k      = (log(optcnst_k_orig(nlam_orig))-log(optcnst_k_orig(nlam_orig-1)))/ &
                       (log(lambda_cm_orig(nlam_orig))-log(lambda_cm_orig(nlam_orig-1)))
     else
        lam_orig_min = lambda_cm_orig(nlam_orig)
        lam_orig_max = lambda_cm_orig(1)
        opt_n_min    = optcnst_n_orig(nlam_orig)
        opt_n_max    = optcnst_n_orig(1)
        opt_k_min    = optcnst_k_orig(nlam_orig)
        opt_k_max    = optcnst_k_orig(1)
        slope_n      = (log(optcnst_n_orig(2))-log(optcnst_n_orig(1)))/ &
                       (log(lambda_cm_orig(2))-log(lambda_cm_orig(1)))
        slope_k      = (log(optcnst_k_orig(2))-log(optcnst_k_orig(1)))/ &
                       (log(lambda_cm_orig(2))-log(lambda_cm_orig(1)))
     endif
     !
     ! Now do the inter/extra-polation
     !
     do ilam=1,nlam
        lam = lambda_cm(ilam)
        if(lam*0.999.lt.lam_orig_min) then
           optcnst_n(ilam) = opt_n_min
           optcnst_k(ilam) = opt_k_min
        elseif(lam*1.001.gt.lam_orig_max) then
           optcnst_n(ilam) = opt_n_max*exp((log(lam)-log(lam_orig_max))*slope_n)
           optcnst_k(ilam) = opt_k_max*exp((log(lam)-log(lam_orig_max))*slope_k)
        else
           ilam_orig = 1
           if(lambda_cm_orig(1).lt.lambda_cm_orig(nlam_orig)) then
              do ilamtmp=2,nlam_orig-1
                 if(lambda_cm_orig(ilamtmp).le.lam) ilam_orig=ilamtmp
              enddo
           else
              do ilamtmp=2,nlam_orig-1
                 if(lambda_cm_orig(ilamtmp).ge.lam) ilam_orig=ilamtmp
              enddo
           endif
           eps = (log(lam)-log(lambda_cm_orig(ilam_orig)))/                           &
                 (log(lambda_cm_orig(ilam_orig+1))-log(lambda_cm_orig(ilam_orig)))
           optcnst_n(ilam) = exp((1.d0-eps)*log(optcnst_n_orig(ilam_orig))+           &
                                        eps*log(optcnst_n_orig(ilam_orig+1)))
           optcnst_k(ilam) = exp((1.d0-eps)*log(optcnst_k_orig(ilam_orig))+           &
                                        eps*log(optcnst_k_orig(ilam_orig+1)))
        endif
     enddo
  endif
  !
  ! Make the grain size grid
  !
  if(nagr.lt.4) then
     write(*,*) 'ERROR: Must have at least 4 grain size sampling points.'
     stop
  endif
  allocate(agrain_cm(nagr),weight(nagr),mgrain(nagr))
  do ia=1,nagr
     agrain_cm(ia) = log(agrain_cm_mean)+wfact*logawidth*(2*(ia-1.d0)/(nagr-1.d0)-1.d0)
     agrain_cm(ia) = exp(agrain_cm(ia))
  enddo
  !
  ! Compute mass of grains
  !
  do ia=1,nagr
     mgrain(ia) = (4.d0*pi/3.d0)*xigrain*agrain_cm(ia)**3
  enddo
  !
  ! Make the grain size distribution according to a Gauss
  ! in log(a), centered on log(a_mean), and with a width
  ! in log(a) space logawidth
  !
  sum = 0.d0
  do ia=1,nagr
     weight(ia) = exp(-0.5*((log(agrain_cm(ia)/agrain_cm_mean))/logawidth)**2)
     sum = sum + weight(ia)
  enddo
  weight(:) = weight(:) / sum
  !
  ! Now allocate the opacity arrays
  !
  allocate(kappa_abs(nlam),kappa_sca(nlam),kappa_g(nlam),scalefact(nlam))
  allocate(zscat(6,2*MXNANG-1,nlam))
  !
  ! Reset things
  !
  zscat(:,:,:) = 0.d0
  kappa_abs(:) = 0.d0
  kappa_sca(:) = 0.d0
  kappa_g(:)   = 0.d0
  !
  ! Start the loop over grain sizes
  ! 
  do ia=1,nagr
     !
     ! Compute geometric cross section
     !
     siggeom = pi*agrain_cm(ia)**2
     !
     ! Now do the loop over wavelengths
     !
     highest_x  = 0.d0
     ihighest_x = -1
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
        rad = agrain_cm(ia)
        !
        ! Wavelength in cm
        !
        wavel = lambda_cm(ilam)
        !
        ! Compute the dimensionless grain size size
        !
        X=2.E0*PI*RAD*REFMED/WAVEL
        !
        ! Check the highest x value
        !
        if((X.gt.highest_x).and.(X.lt.max_x)) then
           highest_x  = X
           ihighest_x = ilam
        endif
        !
        ! If X<max_x then call BHMIE
        !
        if(X.lt.max_x) then
           !
           ! Call BHMie
           !
           CALL BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA)
           QABS=QEXT-QSCA
           !
           ! Put results into array
           !
           ! Note: The averaging of g has to be done multiplied by kappa_scat, 
           !       otherwise the answer is wrong.
           !
           kappa_abs(ilam) = kappa_abs(ilam) + weight(ia)*qabs*siggeom/mgrain(ia)
           kappa_sca(ilam) = kappa_sca(ilam) + weight(ia)*qsca*siggeom/mgrain(ia)
           kappa_g(ilam)   = kappa_g(ilam)   + weight(ia)*GSCA*qsca*siggeom/mgrain(ia)
           !
           ! Compute conversion factor from the Sxx matrix elements
           ! from the Bohren & Huffman code to the Zxx matrix elements we
           ! use (such that 2*pi*int_{-1}^{+1}Z11(mu)dmu=kappa_scat).
           ! This includes the factor k^2 (wavenumber squared) to get 
           ! the actual cross section in units of cm^2 / ster, and there 
           ! is the mass of the grain to get the cross section per gram.
           !
           factor = (lambda_cm(ilam)/(2*PI))**2/mgrain(ia)
        else
           !
           ! Here we do a trick to guestimate the value for very high X values
           !
           ! Call BHMie
           !
           XMAX = max_x
           CALL BHMIE(XMAX,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA)
           QABS=QEXT-QSCA
           !
           ! Put results into array
           !
           kappa_abs(ilam) = kappa_abs(ilam) + weight(ia)*qabs*siggeom/mgrain(ia)
           kappa_sca(ilam) = kappa_sca(ilam) + weight(ia)*qsca*siggeom/mgrain(ia)
           kappa_g(ilam)   = kappa_g(ilam)   + weight(ia)*GSCA*qsca*siggeom/mgrain(ia)
           !
           ! Compute conversion factor from the Sxx matrix elements
           ! from the Bohren & Huffman code to the Zxx matrix elements we
           ! use (such that 2*pi*int_{-1}^{+1}Z11(mu)dmu=kappa_scat).
           ! This includes the factor k^2 (wavenumber squared) to get 
           ! the actual cross section in units of cm^2 / ster, and there 
           ! is the mass of the grain to get the cross section per gram.
           !
           factor = (X/max_x)**2*(lambda_cm(ilam)/(2*PI))**2/mgrain(ia)
        endif
        !
        ! Also store the Z matrix elements
        !
        NAN=2*NANG-1
        DO J=1,NAN
           S11=0.5E0*CABS(S2(J))*CABS(S2(J))
           S11=S11+0.5E0*CABS(S1(J))*CABS(S1(J))
           S12=0.5E0*CABS(S2(J))*CABS(S2(J))
           S12=S12-0.5E0*CABS(S1(J))*CABS(S1(J))
           POL=-S12/S11
           S33=REAL(S2(J)*CONJG(S1(J)))
           S34=AIMAG(S2(J)*CONJG(S1(J)))
           zscat(1,j,ilam) = zscat(1,j,ilam) + weight(ia) * S11 * factor
           zscat(2,j,ilam) = zscat(2,j,ilam) + weight(ia) * S12 * factor
           zscat(3,j,ilam) = zscat(3,j,ilam) + weight(ia) * S11 * factor
           zscat(4,j,ilam) = zscat(4,j,ilam) + weight(ia) * S33 * factor
           zscat(5,j,ilam) = zscat(5,j,ilam) + weight(ia) * S34 * factor
           zscat(6,j,ilam) = zscat(6,j,ilam) + weight(ia) * S33 * factor
        enddo
        !
        ! End loop over all wavelengths
        !
     enddo
     !
     ! End loop over all grain sizes
     ! 
  enddo
  !
  ! Now renormalize the g factor
  ! 
  kappa_g(:) = kappa_g(:) / kappa_sca(:)
  !
  ! Check if the sum of the S11 over all angles is indeed kappa_scat
  !
  do ilam=1,nlam
     sum = 0.d0
     do j=2,nan
        sum = sum + 0.25d0*(zscat(1,j-1,ilam)+zscat(1,j,ilam))* &
              abs(mu(j)-mu(j-1))
     enddo
     sum = sum * 4*PI
     error = abs(sum/kappa_sca(ilam)-1.d0)
     if(error.gt.errmax) then
        write(*,*) 'ERROR: At lambda=',lambda_cm(ilam)*1d4,'micron the error in the',&
             ' scattering integral is ',error,'which is larger than the error limit',errmax
        write(*,*) '    kappa_scat                     = ',kappa_sca(ilam)
        write(*,*) '    2*pi*int_{-1}^{+1} Z_11(mu)dmu = ',sum
        write(*,*) '    Please use a larger number of angle points or take a weaker error limit (5th line in param.inp).'
        write(*,*) '    Note: This problem usually happens for large ratio of a/lambda.'
        close(2)
        stop
     else
        if(irescalez.gt.0) then
           ! ONLY FOR BACKWARD COMPATIBILITY - SHOULD NOT BE USED (SET irescalez=0)
           scalefact(ilam) = kappa_sca(ilam) / sum
           zscat(1:6,1:nan,ilam) = zscat(1:6,1:nan,ilam) * scalefact(ilam)
        endif
     endif
     !
     ! Do the "chopping" of excessive forward scattering
     !
     if(chopforward.gt.0.d0) then
        iang = 0
        do j=1,nan
           if(angle(j).lt.chopforward) iang=j
        enddo
        if(iang.gt.0) then
           iang = iang+1
           do j=1,iang-1
              zscat(1:6,j,ilam) = zscat(1:6,iang,ilam)
           enddo
           sum = 0.d0
           do j=2,nan
              sum = sum + 0.25d0*(zscat(1,j-1,ilam)+zscat(1,j,ilam))* &
                   abs(mu(j)-mu(j-1))
           enddo
           sum = sum * 4*PI
           kappa_sca(ilam) = sum
           summu = 0.d0
           do j=2,nan
              summu = summu + 0.25d0*(zscat(1,j-1,ilam)+zscat(1,j,ilam))* &
                   0.5d0*(mu(j)+mu(j-1))*abs(mu(j)-mu(j-1))
           enddo
           summu = summu * 4*PI
           kappa_g(ilam) = summu/sum
        endif
     endif
  enddo
  !
  ! Write the results
  !
  filename = 'dustkapscatmat_'//trim(material)//'.inp'
  open(unit=2,file=filename)
  leng = len_trim(material)
  if(leng.lt.10) then
     write(str0,'(I1)')
  elseif(leng.lt.100) then
     write(str0,'(I2)')
  elseif(leng.lt.1000) then
     write(str0,'(I3)')
  else
     write(*,*) 'Dust opacity name too long'
     stop 
  endif
  str1 = '(A41,A'//trim(str0)//')'
  write(2,str1) '# Opacity and scattering matrix file for ',trim(material)
  write(2,'(A109)') '# Please do not forget to cite in your publications the original ' &
       //'paper of these optical constant measurements'
  if(ref1.ne.'') call write_ref(ref1)
  if(ref2.ne.'') call write_ref(ref2)
  if(ref3.ne.'') call write_ref(ref3)
  write(2,'(A47)') '# Made with the make_scatmat_smoothed.f90 code,'
  write(2,'(A70)') '# using the bhmie.f Mie code of Bohren and Huffman (version by Draine)'
  write(2,'(A26)') '# Grain size distribution:'
  write(2,'(A23,E13.6,A3)') '#   agrain_mean      = ',agrain_cm_mean,' cm'
  write(2,'(A23,E13.6)') '#   logawidth        = ',logawidth
  write(2,'(A23,F13.6)') '#   wfact            = ',wfact
  write(2,'(A23,I4,A3)') '#   nr of sizes used = ',nagr
  write(2,'(A19)') '# Material density:'
  call write_dens(xigrain)
  write(2,*) 1     ! Format number
  write(2,*) nlam
  write(2,*) nan
  write(2,*) 
  do ilam=1,nlam
     write(2,'(4(E13.6,1X))') lambda_cm(ilam)*1e4,kappa_abs(ilam),kappa_sca(ilam),kappa_g(ilam)
  enddo
  write(2,*) 
  do j=1,nan
     write(2,'(F13.6)') angle(j)
  enddo
  write(2,*)
  do ilam=1,nlam
     do j=1,nan
        write(2,'(6(E13.6,1X))') zscat(1:6,j,ilam)
     enddo
     write(2,*)
  enddo
  close(2)
  !
  ! Write the simpler dustkappa_*.inp file
  !
  filename = 'dustkappa_'//trim(material)//'.inp'
  open(unit=2,file=filename)
  leng = len_trim(material)
  if(leng.lt.10) then
     write(str0,'(I1)')
  elseif(leng.lt.100) then
     write(str0,'(I2)')
  elseif(leng.lt.1000) then
     write(str0,'(I3)')
  else
     write(*,*) 'Dust opacity name too long'
     stop 
  endif
  str1 = '(A19,A'//trim(str0)//')'
  write(2,str1) '# Opacity file for ',trim(material)
  write(2,'(A109)') '# Please do not forget to cite in your publications the original ' &
       //'paper of these optical constant measurements'
  if(ref1.ne.'') call write_ref(ref1)
  if(ref2.ne.'') call write_ref(ref2)
  if(ref3.ne.'') call write_ref(ref3)
  write(2,'(A47)') '# Made with the make_scatmat_smoothed.f90 code,'
  write(2,'(A70)') '# using the bhmie.f Mie code of Bohren and Huffman (version by Draine)'
  write(2,'(A26)') '# Grain size distribution:'
  write(2,'(A23,E13.6,A3)') '#   agrain_mean      = ',agrain_cm_mean,' cm'
  write(2,'(A23,E13.6)') '#   logawidth        = ',logawidth
  write(2,'(A23,F13.6)') '#   wfact            = ',wfact
  write(2,'(A23,I4,A3)') '#   nr of sizes used = ',nagr
  write(2,'(A19)') '# Material density:'
  call write_dens(xigrain)
  write(2,*) 3     ! Format number
  write(2,*) nlam
  write(2,*) 
  do ilam=1,nlam
     write(2,'(4(E13.6,1X))') lambda_cm(ilam)*1e4,kappa_abs(ilam),kappa_sca(ilam),kappa_g(ilam)
  enddo
  close(2)
  !
  ! Just for information to the user: write out the scaling factor used
  !
  if(irescalez.gt.0) then
     open(unit=1,file='scalefactor.out')
     write(1,*) nlam
     do ilam=1,nlam
        write(1,'(2(E13.6,1X))') lambda_cm(ilam)*1e4,scalefact(ilam)
     enddo
     close(1)
  endif
  !
  ! Deallocate stuff
  !
  deallocate(lambda_cm,optcnst_n,optcnst_k)
  deallocate(lambda_cm_orig,optcnst_n_orig,optcnst_k_orig)
  deallocate(kappa_abs,kappa_sca,kappa_g)
  deallocate(zscat,angle,mu,scalefact)
  deallocate(agrain_cm,weight,mgrain)
end program bhmakeopac

subroutine write_ref(str)
  implicit none
  character*160 :: str,lenstr,fmt
  write(lenstr,'(I3.3)') len_trim(str)
  fmt='(A15,A'//lenstr(1:3)//')'
  write(2,fmt) '# @reference = ',str
end subroutine write_ref

subroutine write_dens(xi)
  implicit none
  double precision :: xi
  character*160 :: fmt
  fmt='(A13,F9.6,A7)'
  write(2,fmt) '# @density = ',xi,' g/cm^3'
end subroutine write_dens
