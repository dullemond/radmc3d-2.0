!=======================================================================
!                         DUST OPACITY MODULE
!
! This dust opacity module contains pieces of code for Mie scattering that
! the author got from Alex de Koter (University of Amsterdam) who presumably
! (?)  got the Mie routine from Joop Hovenier (?).
! =======================================================================
module dust_module
  use rtglobal_module
  use polarization_module
  use mathroutines_module
  !
  ! The dust scattering mode, which is typically automatically set, 
  ! depending on the kind of dust opacities found. 
  !
  integer :: scattering_mode=0           ! =0 Dust scattering not included
                                         ! =1 Include isotropic scattering
                                         ! =2 Include scattering with 
                                         !    Henyey-Greenstein function
                                         ! =3 Include scattering with 
                                         !    full phase function
                                         ! =4 Include scattering with 
                                         !    full phase function and
                                         !    polarization (randomly
                                         !    oriented particles), but
                                         !    only in the scattering
                                         !    source function.
                                         ! =5 Include scattering with 
                                         !    full phase function and
                                         !    polarization (randomly
                                         !    oriented particles), 
                                         !    full treatment.
  integer :: scattering_mode_def=0       ! The default to be used 
  integer :: scattering_mode_max=9999    ! User-specified limit to mode
  !
  ! The dust grain alignment mode. Default 0 = no alignment. 
  ! =-1 means alignment only for the thermal emission seen in the images/spectra.
  !     This mode is merely for testing purposes, because it is not really 
  !     useful to include polarized thermal emission in the images but not
  !     as a source for the scattering. 
  ! =1  means alignment only for the thermal emission, both for images/spectra 
  !     and for the source of photons in the scattering monte carlo. This is
  !     the mode to use if you want to treat thermal polarization but do not
  !     want to include the full aligned scattering stuff (which is much
  !     slower (but as of now not yet even built in)).
  !
  integer :: alignment_mode=0
  !
  ! For those opacities that are read from a file dustkappa_*.inp,
  ! i.e. which are on their own independent frequency grid, we 
  ! may want to store these opacities in their original form. This
  ! can be useful when wanting to do high-spectral-resolution 
  ! raytracing, at a lambda resolution higher than the original
  ! global frequency grid.
  !
  ! Note: "type" gives the type of opacity data stored here:
  !  type = 1       : Opacity without polarization. Phase function only 
  !                   isotropic or Henyey-Greenstein.
  !  type = 2       : Opacity with polarization: Case of spherical particles
  !                   or randomly oriented non-spherical (but non-helical)
  !                   particles. The Zscat array is the scattering matrix
  !                   array (see polarization_module.f90 for definitions):
  !                     Zscat(1:nmu,1:nfreq,1)  =  Z11(1:nmu,1:nfreq) 
  !                     Zscat(1:nmu,1:nfreq,2)  =  Z12(1:nmu,1:nfreq) 
  !                     Zscat(1:nmu,1:nfreq,3)  =  Z22(1:nmu,1:nfreq) 
  !                     Zscat(1:nmu,1:nfreq,4)  =  Z33(1:nmu,1:nfreq) 
  !                     Zscat(1:nmu,1:nfreq,5)  =  Z34(1:nmu,1:nfreq) 
  !                     Zscat(1:nmu,1:nfreq,6)  =  Z44(1:nmu,1:nfreq) 
  !                   The 1:nmu refers to the mu(i) = cos(theta(i)) angles
  !                   listed in the mu(1:nmu) array.
  !  type = 3       : Simplified treatment of polarized thermal emission from
  !                   aligned non-spherical, non-helical grains: No scattering,
  !                   but different thermal emission for parallel and orthogonal
  !                   linear polarizations. NOTE: This is an approximation,
  !                   because scattering should play a role, which for such
  !                   aligned particles is complex.
  !  type = 4       : Full treatment of aligned particles: The data cannot
  !                   be stored here; instead the optical constants are used.
  !
  type dust_kappaarray_link
     integer :: nrfreq,nmu,type,namu
     double precision                        :: sweight     ! Material density [g/cm^3]
     double precision, dimension(:), pointer :: freq        ! Frequency grid [Hz]
     double precision, dimension(:), pointer :: kappa_a     ! kappa_abs [cm^2/gram-of-dust]
     double precision, dimension(:), pointer :: kappa_s     ! kappa_scat [cm^2/gram-of-dust]
     double precision, dimension(:), pointer :: gfactor     ! Henyey-Greenstein g factor = <cos(theta)>
     double precision, dimension(:), pointer :: n,k         ! Optical constants: n and k values
     double precision, dimension(:), pointer :: mu          ! mu=cos(theta) grid for the zscat array
     double precision, dimension(:,:,:), pointer :: zscat   ! The Z Mueller matrix elements (see polarization module)
     double precision, dimension(:), pointer :: alignmu     ! The mu=cos(theta) grid for the alignment ratio
     double precision, dimension(:,:), pointer :: alignorth ! The alignment factor for orthogonal (see polarization module)
     double precision, dimension(:,:), pointer :: alignpara ! The alignment factor for parallel (see polarization module)
  end type dust_kappaarray_link
  !
  ! For each dust species a link to the above array, if associated
  !
  type(dust_kappaarray_link), allocatable :: dust_kappa_arrays(:)
  !
  ! Some stuff for quantum-heated grains
  !
  integer,allocatable :: dust_quantum(:),dust_n_catom(:)
  !
  ! The main dust opacity information
  !
  integer :: dust_nr_species = 0
  double precision, allocatable :: dust_kappa_abs(:,:), dust_kappa_scat(:,:)
  double precision, allocatable :: dust_gfactor(:,:)
  logical :: dust_use_stokes=.false.
  logical :: dust_use_oriented_grains=.false.
  !
  ! Some additional information such as grain size and mass
  !
  double precision, allocatable :: dust_agrain(:),dust_mgrain(:)
  !
contains

!-------------------------------------------------------------------
!                       READ THE DUST DATA
!
!     This is the main routine for the input of all relevant dust
!     data, and the preprocessing of this data, if necessary.
!
!-------------------------------------------------------------------
subroutine read_dustdata(action)
  implicit none
  !
  integer ispec,inu,iformat,idum,idum2,ierr,nrspec,action
  doubleprecision dummy,temp0,temp1,dtemp0,dtemp1,dum3,agrain
  character*80 comstring,filename,base,ext,iduststring,algstring,filenamealt
  logical fex
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(dust_kappa_abs)) return
  endif
  !
  ! Start with a fresh cleanup
  !   
  call dust_cleanup()
  !
  ! Check if the main frequency grid is already present
  !
  if((freq_nr.le.0).or.(.not.allocated(freq_nu))) then
     write(stdo,*) 'ERROR: While reading dust data: global frequency grid not set.'
     stop
  endif
  !
  ! Now open the master file for the dust
  !
  inquire(file='dustopac.inp',exist=fex)
  if(.not.fex) then
     write(stdo,*) 'ERROR: Could not find file dustopac.inp...'
     stop
  endif
  !
  ! Message
  !
  write(stdo,*) 'Reading dust data...'
  call flush(stdo)
  !
  ! Open file
  !
  open(unit=3,file='dustopac.inp',status='old',err=701)
  read(3,*) iformat
  !
  ! Find how many dust species are to be read
  !
  read(3,*) nrspec
  !
  ! If the nr of species already has been set by some other module of
  ! this program, then we must compare
  !
  if(dust_nr_species.gt.0) then
     !
     ! The dust grain sizes and masses were already set
     !
     if(nrspec.ne.dust_nr_species) then
        write(stdo,*) 'ERROR in dust_module: The dust_nr_species was '
        write(stdo,*) '      already set, but at a different value than'
        write(stdo,*) '      I now read from the dustopac.inp file.'
        stop
     endif
     if((.not.allocated(dust_agrain)).or.(.not.allocated(dust_mgrain))) then
        write(stdo,*) 'ERROR in dust_module: dust_nr_species was already'
        write(stdo,*) '      set beforehand, but agrain and mgrain not.'
        stop
     endif
  else
     !
     ! The dust grain sizes and masses will be set here
     !
     dust_nr_species = nrspec
     if((allocated(dust_agrain)).or.(allocated(dust_mgrain))) then
        write(stdo,*) 'ERROR in dust_module: dust_nr_species was 0'
        write(stdo,*) '      but agrain and/or mgrain were allocated.'
        stop
     endif
     allocate(dust_agrain(1:dust_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate dust_agrain'
        stop
     endif
     allocate(dust_mgrain(1:dust_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate dust_mgrain'
        stop
     endif
     do ispec=1,dust_nr_species
        dust_agrain(ispec) = 0.d0
        dust_mgrain(ispec) = 0.d0
     enddo
  endif
  !
  ! Allocate all the other arrays
  !
  allocate(dust_kappa_arrays(1:dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate dust_kappa_arrays'
     stop
  endif
  do ispec=1,dust_nr_species 
     nullify(dust_kappa_arrays(ispec)%freq)
     nullify(dust_kappa_arrays(ispec)%kappa_a)
     nullify(dust_kappa_arrays(ispec)%kappa_s)
     nullify(dust_kappa_arrays(ispec)%gfactor)
     nullify(dust_kappa_arrays(ispec)%n)
     nullify(dust_kappa_arrays(ispec)%k)
     nullify(dust_kappa_arrays(ispec)%mu)
     nullify(dust_kappa_arrays(ispec)%zscat)
     nullify(dust_kappa_arrays(ispec)%alignmu)
     nullify(dust_kappa_arrays(ispec)%alignorth)
     nullify(dust_kappa_arrays(ispec)%alignpara)
     dust_kappa_arrays(ispec)%nrfreq = 0
     dust_kappa_arrays(ispec)%nmu = 0
     dust_kappa_arrays(ispec)%type = 0
  enddo
  allocate(dust_quantum(1:dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate dust quantum info'
     stop
  endif
  allocate(dust_kappa_abs(1:freq_nr,1:dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate dust opacities'
     stop
  endif
  allocate(dust_kappa_scat(1:freq_nr,1:dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate dust opacities'
     stop
  endif
  allocate(dust_gfactor(1:freq_nr,1:dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate dust opacities'
     stop
  endif
  dust_gfactor(:,:) = 0.d0
  !
  ! Now read
  !
  if(iformat.eq.3) then
     write(stdo,*) 'ERROR: dustopac.inp format 3 no longer in use.'
     stop
  endif
  read(3,*) comstring      
  !
  ! Now a loop over the dust species
  !
  do ispec=1,dust_nr_species
     read(3,*) idum
     if(iformat.ge.2) then
        read(3,*) idum2
        if(idum2.ne.0) then 
            dust_quantum(ispec) = 1
            if(idum2.eq.2) then
               !
               ! Read the nr of C-atoms for this PAH molecule 
               ! (this is only done if quantum is on for this species)
               !
               read(3,*) dum3
               dust_n_catom(ispec) = dum3
            elseif(idum2.eq.3) then
               ! OLD MODE, DONT USE
               stop
               !!!! read(3,*) dum3
               !!!! dust_n_catom(ispec) = dum3
               !!!! read(3,*) dum3
               !!!! dust_tmax(ispec) = dum3
            elseif(idum2.eq.4) then
               !
               ! This is the real mode for the PAH destruction temperature
               ! This is not part of this code yet (18.04.06).
               ! MODIFIED 18.04.06
               !
               !!!! read(3,*) dum3
               !!!! dust_mgrain(ispec) = dum3*12*mp
               !!!! dust_n_catom(ispec) = dum3
               !!!! read(3,*) dum3
               !!!! pahdes_temp(ispec) = dum3
               stop 
            elseif(idum2.eq.5) then
               !
               ! This is the real mode for the PAH destruction temperature
               ! With Habart's recipe for destruction...
               ! This is not part of this code yet (25.11.06).
               ! MODIFIED 25.11.06
               !
               !!!! read(3,*) dum3
               !!!! dust_n_catom(ispec) = dum3
               !!!! read(3,*) dum3
               !!!! read(3,*) dum3
            else
               dust_n_catom(ispec) = 0
            endif
         else
            dust_quantum(ispec)=0
         endif
      else
         dust_quantum(ispec)=0
      endif
      !
      ! Now find out which dust opacity file dustopac_*.inp we should
      ! read for this dust species
      !
      !!!!read(3,*) idustfile
      !
      read(3,*) iduststring
      !
      ! Check out how to read the file
      !
      if(idum.lt.0) then
         !
         ! Use the dustopac_*.inp files. These files simply list the
         ! opacities for each of the wavelength points listed in 
         ! frequency.inp or wavelength_micron.inp
         !
         base='dustopac_'
         ext ='.inp'
         !
         !!!!call make_indexed_filename(base,idustfile,ext,filename)
         !
         call make_indexed_filename_string(base,iduststring,ext,filename)
         !
         if(idum.eq.-1) then
            !
            ! Simple temperature-independent opacities
            !
            call read_dustopac_file(ispec,filename)
            !
         elseif(idum.eq.-2) then
            !
            ! Temperature-dependent opacities. 
            ! In RADMC-3D this mode is not yet reimplemented.
            !
            !!!! read(3,*) temp0
            !!!! read(3,*) temp1
            !!!! dust_tmin(ispec)  = temp0
            !!!! dust_tmax(ispec)  = temp1
            !!!! dust_dtmin(ispec) = 0.d0
            !!!! dust_dtmax(ispec) = 0.d0
            !!!! call read_dustopac_file(ispec,filename)
            stop
         elseif(idum.eq.-3) then
            !
            ! Temperature-dependent opacities
            ! Not only tmin/tmax, but also smooth switch-off
            ! In RADMC-3D this model is not yet reimplemented
            !
            !!!! read(3,*) temp0
            !!!! read(3,*) dtemp0
            !!!! read(3,*) temp1
            !!!! read(3,*) dtemp1
            !!!! dust_tmin(ispec)  = temp0
            !!!! dust_tmax(ispec)  = temp1
            !!!! dust_dtmin(ispec) = abs(dtemp0)
            !!!! dust_dtmax(ispec) = abs(dtemp1)
            !!!! call read_dustopac_file(ispec,filename)
            stop
         elseif(idum.eq.-4) then
            !
            ! Temperature-dependent opacities: full model
            ! In RADMC-3D this model is not yet reimplemented
            !
            stop
         else   
            write(stdo,*) 'While reading dustopac.inp: other ',         &
                   'input modes not yet implemented'
            stop
         endif
      else
         if(idum.eq.1) then
            !
            ! This is the new style of dust opacity files: instead of having
            ! to match the dust opacity values directly to the frequency
            ! points in the frequency.inp or wavelength_micron.inp file as
            ! with the dustopac_*.inp files, these dustkappa_*.inp files
            ! have their own frequency grid, preferably at very high
            ! resolution. RADMC-3D then automatically interpolates these
            ! onto the grid of frequency.inp / wavelength_micron.inp.
            !
            base='dustkappa_'
            ext ='.inp'
            !
            !!!!call make_indexed_filename(base,idustfile,ext,filename)
            !
            call make_indexed_filename_string(base,iduststring,ext,filename)
            if(.not.allocated(freq_nu)) then
               write(stdo,*) 'ERROR while reading dust opacity file ',filename
               write(stdo,*) '   The frequency array is not yet set.'
               write(stdo,*) '   This must be set before reading ',filename
               stop
            endif
            call read_dustkappa_file(ispec,filename)
            !
            ! To prevent confusion, let's check if perhaps also a file
            ! dustkapscatmat_xxx.inp is present. If so, then warn.
            !
            base='dustkapscatmat_'
            ext ='.inp'
            call make_indexed_filename_string(base,iduststring,ext,filenamealt)
            inquire(file=filenamealt,exist=fex)
            if(fex) then
               write(stdo,*) 'WARNING: I also found a file ',trim(filenamealt),&
                    ' in addition to ',trim(filename),' (is ok, but just so you are aware).'
            endif
            !
         elseif((idum.ge.10).and.(idum.lt.30)) then
            !
            ! Dust opacities with the fully mu-dependence and polarization
            ! for scattering, assuming spherical particles and/or randomly
            ! oriented particles. In other words: we read in the Z matrix
            ! elements (Z11, Z12, Z22, Z33, Z34, Z44) as a function of
            ! scattering angle mu=cos(theta). See polarization module
            ! for details.
            !
            base='dustkapscatmat_'
            ext ='.inp'
            !
            !!!!call make_indexed_filename(base,idustfile,ext,filename)
            !
            call make_indexed_filename_string(base,iduststring,ext,filename)
            if(.not.allocated(freq_nu)) then
               write(stdo,*) 'ERROR while reading dust opacity and scattering matrix file ',filename
               write(stdo,*) '   The frequency array is not yet set.'
               write(stdo,*) '   This must be set before reading ',filename
               stop
            endif
            call read_dustkapscatmat_file(ispec,filename)
            !
            ! To prevent confusion, let's check if perhaps also a file
            ! dustkappa_xxx.inp is present. If so, then warn.
            !
            base='dustkappa_'
            ext ='.inp'
            call make_indexed_filename_string(base,iduststring,ext,filenamealt)
            inquire(file=filenamealt,exist=fex)
            if(fex) then
               write(stdo,*) 'WARNING: I also found a file ',trim(filenamealt),&
                    ' in addition to ',trim(filename),' (is ok, but just so you are aware).'
            endif
            !
            ! If necessary, also read the angular-dependency of the
            ! parallel and orthogonal absorption coefficient of the grains,
            ! assuming they are aligned with some direction. Note that this
            ! is a very simplistic treatment of polarized thermal emission
            ! by aligned grains, in which the emission/absorption is treated
            ! as aligned grains while the scattering is treated as randomly
            ! oriented grains (which the read_dustkapscatmat_file handles).
            ! But this is the easiest method, and is presumably sufficient
            ! for most purposes.
            !
            if(idum.eq.20) then
               base='dustkapalignfact_'
               ext ='.inp'
               call make_indexed_filename_string(base,iduststring,ext,filenamealt)
               call read_dustalign_angdep(ispec,filenamealt)
            endif
            !
         elseif(idum.eq.100) then
            !
            ! Using this input style we simply read in the optical constants
            ! (n and k) and compute the opacities here on-the-fly
            !
            base='dustoptnk_'
            ext ='.inp'
            call make_indexed_filename_string(base,iduststring,ext,filename)
            if(.not.allocated(freq_nu)) then
               write(stdo,*) 'ERROR while reading dust opacity file ',filename
               write(stdo,*) '   The frequency array is not yet set.'
               write(stdo,*) '   This must be set before reading ',filename
               stop
            endif
            read(3,*) agrain            ! Grain radius in cm
            read(3,*) algstring         ! 'MIE' or 'CDE'
            if(dust_agrain(ispec).eq.0.d0) then
               dust_agrain(ispec) = agrain
            else
               if(abs(1.d0-agrain/dust_agrain(ispec)).gt.1d-4) then
                  write(stdo,*) 'ERROR in dust_module: The grain size '
                  write(stdo,*) '      specified previously is not equal'
                  write(stdo,*) '      to the one from dustopac.inp.'
                  write(stdo,*) '      ispec           = ',ispec
                  write(stdo,*) '      agrain previous = ',dust_agrain(ispec)
                  write(stdo,*) '      agrain now      = ',agrain
                  stop
               endif
            endif
            call read_optconst_file(ispec,filename)
            call make_dust_opacity(ispec,agrain,algstring)
            !
         else
            write(stdo,*) 'While reading dustopac.inp: other ',         &
                   'input modes not yet implemented'
            stop
         endif
      endif
      read(3,*) comstring                
  enddo
  !
  ! Set the scattering mode to the default
  !
  scattering_mode = scattering_mode_def
  !
  close(3)
  !
  goto 710
701 continue
  write(stdo,*) 'Could not open file dustopac.inp'
  stop 
710 continue
  return
end subroutine read_dustdata


!-------------------------------------------------------------------
!                       REMOVE AN OPACITY 
!-------------------------------------------------------------------
subroutine dust_remove_opacity(ispec)
  implicit none
  integer :: ispec
  if(.not.allocated(dust_kappa_arrays)) then
     write(stdo,*) 'ERROR in dust module: cannot remove opacity if dust_kappa_arrays not present.'
     stop
  endif
  if(associated(dust_kappa_arrays(ispec)%freq))  &
       deallocate(dust_kappa_arrays(ispec)%freq)
  if(associated(dust_kappa_arrays(ispec)%kappa_a)) &
       deallocate(dust_kappa_arrays(ispec)%kappa_a)
  if(associated(dust_kappa_arrays(ispec)%kappa_s)) &
       deallocate(dust_kappa_arrays(ispec)%kappa_s)
  if(associated(dust_kappa_arrays(ispec)%gfactor)) &
       deallocate(dust_kappa_arrays(ispec)%gfactor)
  if(associated(dust_kappa_arrays(ispec)%n)) &
       deallocate(dust_kappa_arrays(ispec)%n)
  if(associated(dust_kappa_arrays(ispec)%k)) &
       deallocate(dust_kappa_arrays(ispec)%k)
  if(associated(dust_kappa_arrays(ispec)%mu)) &
       deallocate(dust_kappa_arrays(ispec)%mu)
  if(associated(dust_kappa_arrays(ispec)%zscat)) &
       deallocate(dust_kappa_arrays(ispec)%zscat)
  if(associated(dust_kappa_arrays(ispec)%alignmu)) &
       deallocate(dust_kappa_arrays(ispec)%alignmu)
  if(associated(dust_kappa_arrays(ispec)%alignorth)) &
       deallocate(dust_kappa_arrays(ispec)%alignorth)
  if(associated(dust_kappa_arrays(ispec)%alignpara)) &
       deallocate(dust_kappa_arrays(ispec)%alignpara)
end subroutine dust_remove_opacity

!-------------------------------------------------------------------
!                       CLEAN UP DUST MODULE
!-------------------------------------------------------------------
subroutine dust_cleanup
  implicit none
  integer :: ispec
  if(allocated(dust_kappa_abs)) deallocate(dust_kappa_abs)
  if(allocated(dust_kappa_scat)) deallocate(dust_kappa_scat)
  if(allocated(dust_gfactor)) deallocate(dust_gfactor)
  if(allocated(dust_quantum)) deallocate(dust_quantum)
  if(allocated(dust_n_catom)) deallocate(dust_n_catom)
  if(allocated(dust_agrain)) deallocate(dust_agrain)
  if(allocated(dust_mgrain)) deallocate(dust_mgrain)
  if(allocated(dust_kappa_arrays)) then
     do ispec=1,dust_nr_species
        call dust_remove_opacity(ispec)
     enddo
     deallocate(dust_kappa_arrays)
  endif
  dust_nr_species = 0
  scattering_mode_def = 0
  scattering_mode = 0
  dust_use_stokes=.false.
  dust_use_oriented_grains=.false.
end subroutine dust_cleanup



!-------------------------------------------------------------------
!        READ THE DUST OPACITY FILES: FREQUENCY-MATCHED STYLE
!
! This is the classic style of opacity input file. It is still the most
! reliable way, because the dust opacity values must be listed precisely for
! the same frequency points as in the frequency.inp or wavelength_micron.inp
! file. In other words: exactly the opacity values are used as listed
! here. In the dustopac.inp file, if you specify a negative number as the
! first line of each dust opacity specification (typically -1, but see
! read_dustdata() for other options), then this dust opacity file style is
! used and read_dustopac_file() is called.
! -------------------------------------------------------------------
subroutine read_dustopac_file(ispec,filename)
  implicit none
  character*80 filename
  integer ispec
  integer ifr,isize,nsize,idum,ierr
  doubleprecision dummy
  logical readtrange,fex
  !
  ! Check if the file exists
  !
  inquire(file=filename,exist=fex)
  if(.not.fex) then
     write(stdo,*) 'ERROR: Could not find file ',trim(filename)
     write(stdo,*) '       See dustopac.inp where this filename is listed'
     stop
  endif
  !
  ! Remove any pre-existing opacity, if present
  !
  call dust_remove_opacity(ispec)
  !
  ! Open the file 
  !
  open(unit=1,file=filename,status='old',err=701)
  read(1,*) ifr,nsize
  if(ifr.eq.-1) then
     !
     ! ifr=-1 is a signal to say that the dust opacity file consists
     ! of opacities at a series of temperature points. This means:
     ! temperature-dependent opacities.
     ! NOTE: For now in RADMC-3D this is disabled, so we quit here.
     !
     write(stdo,*) 'ERROR: For now in RADMC-3D temperature-dependent opacities'
     write(stdo,*) '       are not implemented.'
     stop
  else
     !
     ! Else, there is only one (fixed) opacity table for the dust,
     ! independent on temperature
     !
  endif
  if(nsize.ne.1) then 
     write(stdo,*) 'ERROR: The dustopacity file ',filename,' uses more than'
     write(stdo,*) '       dust size. Use multiple species instead.'
     stop
  endif
  !
  ! Do some checks
  !
  if((freq_nr.le.0).or.(.not.allocated(freq_nu))) then
     write(stdo,*) 'ERROR: While reading dust data: global frequency grid not set.'
     stop
  endif
  if(ifr.ne.freq_nr) then
     write(stdo,*) 'Number of frequencies in ',filename
     write(stdo,*) 'not equal to the number in frequency.inp/wavelength_micron.inp'
     write(stdo,*) ifr,freq_nr
     stop 13
  endif
  if((ispec.lt.1).or.(ispec.gt.dust_nr_species)) then
     write(stdo,*) 'INTERNAL ERROR in dust reading'
     stop
  endif
  if(.not.allocated(dust_kappa_arrays)) then 
     write(stdo,*) 'ERROR: In version 0.32 and higher the dust_kappa_arrays MUST be allocated.'
     stop
  endif
  !
  ! Now allocate the arrays for storing this data
  !
  allocate(dust_kappa_arrays(ispec)%freq(1:freq_nr),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%kappa_a(1:freq_nr),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%kappa_s(1:freq_nr),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%gfactor(1:freq_nr),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  dust_kappa_arrays(ispec)%nrfreq     = freq_nr
  dust_kappa_arrays(ispec)%freq(:)    = freq_nu(:)
  dust_kappa_arrays(ispec)%gfactor(:) = 0.d0
  !
  ! Read the absorption opacity
  !
  do ifr=1,freq_nr
     read(1,*) dummy
     dust_kappa_abs(ifr,ispec) = dummy
     dust_kappa_arrays(ispec)%kappa_a(ifr) = dummy
     !
     ! Check for negative opacities
     !
     if(dummy.le.0.d0) then 
        write(stdo,*) 'ERROR: Found negative or zero absorption opacity at inu=',ifr
        write(stdo,*) '       File = ',filename
        write(stdo,*) '       Aborting...'
        stop 
     endif
     !
  enddo
  !
  ! Read the scattering opacity
  !
  do ifr=1,freq_nr
     read(1,*) dummy
     dust_kappa_scat(ifr,ispec) = dummy
     dust_kappa_arrays(ispec)%kappa_s(ifr) = dummy
     !
     ! Check for negative opacities
     !
     if(dummy.lt.0.d0) then 
        write(stdo,*) 'ERROR: Found negative scattering opacity at inu=',ifr
        write(stdo,*) '       File = ',filename
        write(stdo,*) '       Aborting...'
        stop 
     endif
     !
     ! If one or more of these opacities is larger than 0, then we should put the
     ! scattering mode at least to 1. 
     !
     if(dummy.gt.0.d0) scattering_mode_def = min(max(1,scattering_mode_def),scattering_mode_max)
     !
  enddo
  !
  close(1)
  !
  goto 710
701 continue
  write(stdo,*) 'Could not open file ',filename
  stop 13
710 continue
  return
end subroutine read_dustopac_file


!-------------------------------------------------------------------
!      READ THE DUST OPACITY FILES: INDEPENDENT FREQ-GRID STYLE
!
! This is the new style of the dustopacity files. These files can
! have their own (preferably high resolution) frequency/wavelength
! grid. This is then remapped onto the internally used frequency
! grid, i.e. the one listed in frequency.inp or wavelength_micron.inp
!
!-------------------------------------------------------------------
subroutine read_dustkappa_file(ispec,filename)
  implicit none
  character*80 :: filename
  integer :: ispec
  character*160 :: string
  character*80 :: filenamecheck
  integer :: iformat,nlam,ilam,ierr,iline,filename_len,i,k
  double precision, allocatable :: lambda(:), kappa_a(:), kappa_s(:)
  double precision, allocatable :: nu(:), gfactor(:)
  double precision :: dummy2(1:2),dummy3(1:3),dummy4(1:4),dum
  double precision :: numin,numax,frmin,frmax,sgn
  logical :: flag,fex
  !
  ! First check if the frequency grid is present
  !
  if((.not.allocated(freq_nu)).or.(freq_nr.le.0)) then
     write(stdo,*) 'ERROR: While reading dust opacity: frequency array not yet set.'
     stop
  endif
  !
  ! Check if the dust_kappa_arrays is present
  !
  if(.not.allocated(dust_kappa_arrays)) then 
     write(stdo,*) 'ERROR: In version 0.32 and higher the dust_kappa_arrays MUST be allocated.'
     stop
  endif
  !
  ! Remove any pre-existing opacity, if present
  !
  call dust_remove_opacity(ispec)
  !
  ! Open the file 
  !
  filename_len = len_trim(filename)
  open(unit=1,file=filename,status='old',err=101)
  !
  ! Skip all lines at the beginning of file that start with one of the
  ! symbols ";", "#" or "!". These are meant as comments. 300 such lines
  ! is the maximum. 
  !
  do iline=1,300
     read(1,'(A158)',end=10) string
     i=1
     do k=1,100
        if(string(i:i).ne.' ') goto 30
        i=i+1
     enddo
     write(stdo,*) 'ERROR While reading ',filename(1:filename_len)
     write(stdo,*) '  Detected line full of spaces...'
30   continue
     if((string(i:i).ne.';').and.(string(i:i).ne.'#').and.(string(i:i).ne.'!')) then
        goto 20
     endif
  enddo
  stop 8901
10 continue
  write(stdo,*) 'ERROR while reading file ',filename(1:filename_len)
  write(stdo,*) 'Did not find data before end of file.'
  stop 8902
20 continue
  !
  ! First read the format number. Do so from the above read string.
  !
  read(string,*) iformat
  !
  ! Read the number of frequencies
  !
  read(1,*) nlam
  !
  ! Allocate arrays
  !
  allocate(lambda(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate lambda array'
     stop
  endif
  allocate(nu(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate nu array'
     stop
  endif
  allocate(kappa_a(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappa_a array'
     stop
  endif
  allocate(kappa_s(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappa_s array'
     stop
  endif
  allocate(gfactor(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate gfactor array'
     stop
  endif
  !
  ! Read the data
  ! Now come the different format styles
  !
  if(iformat.eq.1) then
     do ilam=1,nlam
        read(1,*) dummy2
        lambda(ilam)  = dummy2(1)
        kappa_a(ilam) = dummy2(2)
        kappa_s(ilam) = 0.d0
        gfactor(ilam) = 0.d0
     enddo
  elseif(iformat.eq.2) then
     do ilam=1,nlam
        read(1,*) dummy3
        lambda(ilam)  = dummy3(1)
        kappa_a(ilam) = dummy3(2)
        kappa_s(ilam) = dummy3(3)
        gfactor(ilam) = 0.d0
     enddo
  elseif(iformat.eq.3) then
     do ilam=1,nlam
        read(1,*) dummy4
        lambda(ilam)  = dummy4(1)
        kappa_a(ilam) = dummy4(2)
        kappa_s(ilam) = dummy4(3)
        gfactor(ilam) = dummy4(4)
     enddo
  else
     write(stdo,*) 'ERROR while reading dust opacity file ',filename(1:filename_len)
     write(stdo,*) '   Do not support format ',iformat
     stop
  endif
  !
  ! Close the file
  ! 
  close(1)
  !
  ! Check that the lambda array is monotonically increasing/decreasing
  !
  if(lambda(2).gt.lambda(1)) then
     sgn = 1
  elseif(lambda(2).lt.lambda(1)) then
     sgn = -1
  else
     write(stdo,*) 'ERROR while reading opacity file ',filename(1:filename_len)
     write(stdo,*) '    Cannot have twice the same lambda.'
     stop
  endif
  do ilam=3,nlam
     dum = (lambda(ilam)-lambda(ilam-1))*sgn
     if(dum.eq.0.d0) then
        write(stdo,*) 'ERROR while reading opacity file ',filename(1:filename_len)
        write(stdo,*) '    Cannot have twice the same lambda.'
        stop
     elseif(dum.lt.0.d0) then
        write(stdo,*) 'ERROR while reading opacity file ',filename(1:filename_len)
        write(stdo,*) '    Wavelength array is not monotonic.'
        stop
     endif
  enddo
  !
  ! Check that the opacities are non-zero
  !
  do ilam=1,nlam
     if(kappa_a(ilam).lt.0.d0) then
        write(stdo,*) 'ERROR while reading opacity file ',filename(1:filename_len)
        write(stdo,*) '    Negative kappa_a detected at ilam=',ilam
        stop
     endif
     if(kappa_s(ilam).lt.0.d0) then
        write(stdo,*) 'ERROR while reading opacity file ',filename(1:filename_len)
        write(stdo,*) '    Negative kappa_s detected at ilam=',ilam
        stop
     endif
  enddo
  !
  ! Convert lambda [micron] to frequency [Hz]
  !
  do ilam = 1,nlam
     nu(ilam) = 2.9979d14/lambda(ilam)
  enddo
  !
  ! Give a warning if range of nu does not include range of freq_nu
  !
  if(freq_nu(1).lt.freq_nu(freq_nr)) then
     frmin = freq_nu(1)
     frmax = freq_nu(freq_nr)
  else
     frmin = freq_nu(freq_nr)
     frmax = freq_nu(1)
  endif
  if(nu(1).lt.nu(nlam)) then
     numin = nu(1)
     numax = nu(nlam)
  else
     numin = nu(nlam)
     numax = nu(1)
  endif
  if((frmin.lt.numin).or.(frmax.gt.numax)) then
     write(stdo,*) 'Note: Opacity file ',filename(1:filename_len),' does not cover the '
     write(stdo,*) '      complete wavelength range of the model (i.e. the range'
     write(stdo,*) '      in the file frequency.inp or wavelength_micron.inp).'
     write(stdo,201) '       Range of model:     lambda = [',2.9979d14/frmax, &
                     ',',2.9979d14/frmin,']'
     write(stdo,201) '       Range in this file: lambda = [',2.9979d14/numax, &
                     ',',2.9979d14/numin,']'
201  format(a37,e9.2,a1,e9.2,a1) 
     write(stdo,*) '      Simple extrapolations are used to extend this range.'
  endif
  !
  ! Now remap the kappa_abs, kappa_scat and gfactor onto freq_nu
  !
  call remap_function(nlam,nu,kappa_a,freq_nr,freq_nu,dust_kappa_abs(:,ispec),&
                      0,2,1,ierr)
  call remap_function(nlam,nu,kappa_s,freq_nr,freq_nu,dust_kappa_scat(:,ispec),&
                      0,2,1,ierr)
  call remap_function(nlam,nu,gfactor,freq_nr,freq_nu,dust_gfactor(:,ispec),&
                      0,1,1,ierr)
  !
  ! Check that the opacities are not negative
  !
  flag = .false.
  do ilam=1,freq_nr
     if(dust_kappa_abs(ilam,ispec).lt.1d-90) then
        dust_kappa_abs(ilam,ispec) = 1d-90
        flag = .true.
     endif
     if(dust_kappa_scat(ilam,ispec).lt.0.d0) then
        dust_kappa_scat(ilam,ispec) = 0.d0
        flag = .true.
     endif
  enddo
  if(flag) then
     write(stdo,*) 'WARNING: Found zero absorption opacity or negative scattering opacity. Fixed it.'
  endif
  !
  ! If one or more of these opacities is larger than 0, then we should put the
  ! scattering mode at least to 1. If one of the g factors is non-zero, then we must
  ! in fact do non-isotropic scattering, i.e. scattering_mode_def=2.
  !
  do ilam=1,freq_nr
     if(dust_kappa_scat(ilam,ispec).gt.0.d0) then
        scattering_mode_def = min(max(1,scattering_mode_def),scattering_mode_max)
        if(dust_gfactor(ilam,ispec).ne.0.d0) scattering_mode_def =   &
               min(max(2,scattering_mode_def),scattering_mode_max)
     endif
  enddo
!  !
!  ! Write a ".used" output file
!  !
!  filenamecheck = filename(1:filename_len)//".used"
!  open(unit=1,file=filenamecheck)
!  write(1,*) '# This is the opacity from ',filename(1:filename_len),&
!             ' remapped onto the model wavelength grid.'
!  write(1,*) '# This is not an input file. It is only meant for the user to',&
!             ' check if the remapping went OK. '
!  write(1,*) 3
!  write(1,*) freq_nr
!  do ilam=1,freq_nr
!     write(1,*) 2.9979d14/freq_nu(ilam),dust_kappa_abs(ilam,ispec),&
!                dust_kappa_scat(ilam,ispec),dust_gfactor(ilam,ispec)
!  enddo
!  close(1)
  !
  ! Also store these original opacities
  !
  allocate(dust_kappa_arrays(ispec)%freq(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate original kappas arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%kappa_a(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate original kappas arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%kappa_s(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate original kappas arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%gfactor(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate original kappas arrays'
     stop
  endif
  dust_kappa_arrays(ispec)%nrfreq     = nlam
  dust_kappa_arrays(ispec)%freq(:)    = nu(:)
  dust_kappa_arrays(ispec)%kappa_a(:) = kappa_a(:)
  dust_kappa_arrays(ispec)%kappa_s(:) = kappa_s(:)
  dust_kappa_arrays(ispec)%gfactor(:) = gfactor(:)
  !
  ! Deallocate arrays
  !
  ! NOTE: In the near future I want to retain these arrays for each
  !       dust opacity, so that one can do post-processing at higher
  !       frequency resolution. 
  !
  deallocate(lambda,nu,kappa_a,kappa_s,gfactor)
  !
  ! Return
  !
  return
  !
  ! Error handling
  !
101 continue
  write(stdo,*) 'ERROR: Could not open file ',filename(1:filename_len)
  stop
  !
end subroutine read_dustkappa_file


!-------------------------------------------------------------------
!      READ THE DUST SCATTERING MATRIX AND ABSORPTION DATA
!              FOR RANDOMLY ORIENTED PARTICLES
!-------------------------------------------------------------------
subroutine read_dustkapscatmat_file(ispec,filename)
  implicit none
  character*80 :: filename
  integer :: ispec
  character*160 :: string
  character*80 :: filenamecheck
  integer :: iformat,nlam,nmu,ilam,ierr,iline,filename_len,i,k,imu
  double precision :: dummy4(1:4),dummy6(1:6),dum
  double precision :: numin,numax,frmin,frmax,sgn,error,errormax
  double precision, allocatable :: kappa_scat(:),gfactor(:),angle(:)
  logical :: flag,fex,zpol
  !
  ! First check if the frequency grid is present
  !
  if((.not.allocated(freq_nu)).or.(freq_nr.le.0)) then
     write(stdo,*) 'ERROR: While reading dust opacity: frequency array not yet set.'
     stop
  endif
  !
  ! Check if dust_kappa_arrays is present
  !
  if(.not.allocated(dust_kappa_arrays)) then 
     write(stdo,*) 'ERROR: In version 0.32 and higher the dust_kappa_arrays MUST be allocated.'
     stop
  endif
  !
  ! Remove any pre-existing opacity, if present
  !
  call dust_remove_opacity(ispec)
  !
  ! Open the file 
  !
  filename_len = len_trim(filename)
  open(unit=1,file=filename,status='old',err=101)
  !
  ! Skip all lines at the beginning of file that start with one of the
  ! symbols ";", "#" or "!". These are meant as comments. 300 such lines
  ! is the maximum. 
  !
  do iline=1,300
     read(1,'(A158)',end=10) string
     i=1
     do k=1,100
        if(string(i:i).ne.' ') goto 30
        i=i+1
     enddo
     write(stdo,*) 'ERROR While reading ',filename(1:filename_len)
     write(stdo,*) '  Detected line full of spaces...'
30   continue
     if((string(i:i).ne.';').and.(string(i:i).ne.'#').and.(string(i:i).ne.'!')) then
        goto 20
     endif
  enddo
  stop 8901
10 continue
  write(stdo,*) 'ERROR while reading file ',filename(1:filename_len)
  write(stdo,*) 'Did not find data before end of file.'
  stop 8902
20 continue
  !
  ! First read the format number. Do so from the above read string.
  !
  read(string,*) iformat
  if(iformat.ne.1) then
     write(stdo,*) 'ERROR: Format number ',iformat,' of file ',trim(filename),&
          ' not known.'
     stop
  endif
  !
  ! Read the number of frequencies
  !
  read(1,*) nlam
  dust_kappa_arrays(ispec)%nrfreq = nlam
  !
  ! Read the number of scattering angles
  !
  read(1,*) nmu
  dust_kappa_arrays(ispec)%nmu = nmu
  !
  ! Make a check
  !
  if(nmu.lt.5) then
     write(stdo,*) 'ERROR: RADMC-3D cannot accept scattering matrices for fewer than 5 angles'
     write(stdo,*) '       in file ',filename
     stop
  endif
  !
  ! Allocate arrays
  !
  allocate(dust_kappa_arrays(ispec)%freq(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%kappa_a(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%kappa_s(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%gfactor(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%mu(1:nmu),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%zscat(1:nmu,1:nlam,1:6),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  !
  ! Read the frequency grid and the absorption and scattering opacities
  ! (which, for non-aligned grains, still have meaning in this context).
  ! Also read the gfactor. Both the kappa_s and the gfactor should, in 
  ! principle, also be calculatable from the zscat array. That would 
  ! allow us to do an internal test.
  !
  do ilam=1,nlam
     read(1,*) dummy4
     dust_kappa_arrays(ispec)%freq(ilam)    = 2.9979d14/dummy4(1)
     dust_kappa_arrays(ispec)%kappa_a(ilam) = dummy4(2)
     dust_kappa_arrays(ispec)%kappa_s(ilam) = dummy4(3)
     dust_kappa_arrays(ispec)%gfactor(ilam) = dummy4(4)
  enddo
  !
  ! Read the mu grid
  !
  allocate(angle(1:nmu))
  do imu=1,nmu
     read(1,*) dum
     angle(imu) = dum
     dust_kappa_arrays(ispec)%mu(imu) = cos(angle(imu)*pi/180.d0)
  enddo
  if(abs(dust_kappa_arrays(ispec)%mu(1)-1.d0).gt.1d-12) then
     write(stdo,*) 'ERROR in dustkapscatmat file: mu(1) is not equal to 1.0'
     write(stdo,*) '     Instead angle = ',dum,' and mu = ',dust_kappa_arrays(ispec)%mu(1)
     stop
  else
     dust_kappa_arrays(ispec)%mu(1) = 1.d0
  endif
  if(abs(dust_kappa_arrays(ispec)%mu(nmu)+1.d0).gt.1d-12) then
     write(stdo,*) 'ERROR in dustkapscatmat file: mu(nmu) is not equal to -1.0'
     write(stdo,*) '     Instead angle = ',dum,' and mu = ',dust_kappa_arrays(ispec)%mu(nmu)
     stop
  else
     dust_kappa_arrays(ispec)%mu(nmu) = -1.d0
  endif
  ! !
  ! ! For now we only accept opacity files in which the angle grid
  ! ! is regular. Reason: it makes it easier to choose a globally
  ! ! (i.e. for all dust species) valid scattering angle grid, which
  ! ! we need to define for the Monte Carlo module. At some point in
  ! ! the future maybe (!) we might relax this requirement, but for
  ! ! now we keep it.
  ! !
  ! dum = angle(2) - angle(1)
  ! do imu=3,nmu
  !    if(abs((angle(imu)-angle(imu-1))/dum-1.d0).gt.1d-2) then
  !       write(stdo,*) 'ERROR: RADMC-3D requires a regular grid of scattering angles'
  !       write(stdo,*) '       in file ',filename
  !       stop
  !    endif
  ! enddo
  !
  ! Now read the full scattering matrix. Each line contains:
  !  Z11 Z12 Z22 Z33 Z34 Z44
  !
  do ilam=1,nlam
     do imu=1,nmu
        read(1,*) dummy6
        dust_kappa_arrays(ispec)%zscat(imu,ilam,1) = dummy6(1)
        dust_kappa_arrays(ispec)%zscat(imu,ilam,2) = dummy6(2)
        dust_kappa_arrays(ispec)%zscat(imu,ilam,3) = dummy6(3)
        dust_kappa_arrays(ispec)%zscat(imu,ilam,4) = dummy6(4)
        dust_kappa_arrays(ispec)%zscat(imu,ilam,5) = dummy6(5)
        dust_kappa_arrays(ispec)%zscat(imu,ilam,6) = dummy6(6)
     enddo
  enddo
  !
  ! Close the file
  ! 
  close(1)
  !
  ! Check that the lambda array is monotonically increasing/decreasing
  !
  if(dust_kappa_arrays(ispec)%freq(2).lt.dust_kappa_arrays(ispec)%freq(1)) then
     sgn = -1
  elseif(dust_kappa_arrays(ispec)%freq(2).gt.dust_kappa_arrays(ispec)%freq(1)) then
     sgn = 1
  else
     write(stdo,*) 'ERROR while reading opacity file ',filename(1:filename_len)
     write(stdo,*) '    Cannot have twice the same lambda.'
     stop
  endif
  do ilam=3,nlam
     dum = (dust_kappa_arrays(ispec)%freq(ilam)-dust_kappa_arrays(ispec)%freq(ilam-1))*sgn
     if(dum.eq.0.d0) then
        write(stdo,*) 'ERROR while reading opacity file ',filename(1:filename_len)
        write(stdo,*) '    Cannot have twice the same lambda.'
        stop
     elseif(dum.lt.0.d0) then
        write(stdo,*) 'ERROR while reading opacity file ',filename(1:filename_len)
        write(stdo,*) '    Wavelength array is not monotonic.'
        stop
     endif
  enddo
  !
  ! Check that the opacities are non-zero, and put the scattering mode to
  ! 1 or 2 if the scattering opacity is non-zero 
  !
  do ilam=1,nlam
     if(dust_kappa_arrays(ispec)%kappa_a(ilam).lt.0.d0) then
        write(stdo,*) 'ERROR while reading opacity file ',filename(1:filename_len)
        write(stdo,*) '    Negative kappa_a detected at ilam=',ilam
        stop
     endif
     if(dust_kappa_arrays(ispec)%kappa_s(ilam).lt.0.d0) then
        write(stdo,*) 'ERROR while reading opacity file ',filename(1:filename_len)
        write(stdo,*) '    Negative kappa_s detected at ilam=',ilam
        stop
     endif
     if(dust_kappa_arrays(ispec)%kappa_s(ilam).gt.0.d0) then
        !
        ! If we have scattering, then it is reasonable to use at least 
        ! scattering_mode=3, because we have just read the full phase
        ! function in terms of the Z-matrix elements (see Mishchenko 
        ! book). 
        ! 
        ! But if not only Z11 is non-zero, but also the other Z-elements,
        ! then polarization is included in the scattering matrix. Then
        ! the scattering mode must be 4 (polarized scattering for the
        ! scattering source function) or 5 (polarized scattering in
        ! the Monte Carlo simulation and in the scattering source function).
        ! In both cases for randomly oriented particles. We will pick
        ! 5 as a standard, but if you set scattering_mode_max=4, then 
        ! RADMC-3D will not do polarization in the MC calculation, but
        ! only in the scattering source function.
        !
        zpol = .false.
        do imu=1,nmu
           if((dust_kappa_arrays(ispec)%zscat(imu,ilam,2).ne.0.d0).or. &
              (dust_kappa_arrays(ispec)%zscat(imu,ilam,3).ne.0.d0).or. &
              (dust_kappa_arrays(ispec)%zscat(imu,ilam,4).ne.0.d0).or. &
              (dust_kappa_arrays(ispec)%zscat(imu,ilam,5).ne.0.d0).or. &
              (dust_kappa_arrays(ispec)%zscat(imu,ilam,6).ne.0.d0)) then
              zpol = .true.
           endif
        enddo
        if(zpol) then
           scattering_mode_def = min(max(5,scattering_mode_def),scattering_mode_max)
        else
           scattering_mode_def = min(max(3,scattering_mode_def),scattering_mode_max)
        endif
     endif
  enddo
  !
  ! Give a warning if range of nu does not include range of freq_nu
  !
  if(freq_nu(1).lt.freq_nu(freq_nr)) then
     frmin = freq_nu(1)
     frmax = freq_nu(freq_nr)
  else
     frmin = freq_nu(freq_nr)
     frmax = freq_nu(1)
  endif
  if(dust_kappa_arrays(ispec)%freq(1).lt.dust_kappa_arrays(ispec)%freq(nlam)) then
     numin = dust_kappa_arrays(ispec)%freq(1)
     numax = dust_kappa_arrays(ispec)%freq(nlam)
  else
     numin = dust_kappa_arrays(ispec)%freq(nlam)
     numax = dust_kappa_arrays(ispec)%freq(1)
  endif
  if((frmin.lt.numin).or.(frmax.gt.numax)) then
     write(stdo,*) 'Note: Opacity file ',filename(1:filename_len),' does not cover the '
     write(stdo,*) '      complete wavelength range of the model (i.e. the range'
     write(stdo,*) '      in the file frequency.inp or wavelength_micron.inp).'
     write(stdo,201) '       Range of model:     lambda = [',2.9979d14/frmax, &
                     ',',2.9979d14/frmin,']'
     write(stdo,201) '       Range in this file: lambda = [',2.9979d14/numax, &
                     ',',2.9979d14/numin,']'
201  format(a37,e9.2,a1,e9.2,a1) 
     write(stdo,*) '      Simple extrapolations are used to extend this range.'
  endif
  !
  ! Now remap the kappa_abs, kappa_scat and gfactor onto freq_nu
  !
  call remap_function(nlam,dust_kappa_arrays(ispec)%freq,dust_kappa_arrays(ispec)%kappa_a,freq_nr,freq_nu,dust_kappa_abs(:,ispec),&
                      0,2,1,ierr)
  call remap_function(nlam,dust_kappa_arrays(ispec)%freq,dust_kappa_arrays(ispec)%kappa_s,freq_nr,freq_nu,dust_kappa_scat(:,ispec),&
                      0,2,1,ierr)
  call remap_function(nlam,dust_kappa_arrays(ispec)%freq,dust_kappa_arrays(ispec)%gfactor,freq_nr,freq_nu,dust_gfactor(:,ispec),&
                      0,1,1,ierr)
  !
  ! Check that the opacities are not negative
  !
  flag = .false.
  do ilam=1,freq_nr
     if(dust_kappa_abs(ilam,ispec).lt.1d-90) then
        dust_kappa_abs(ilam,ispec) = 1d-90
        flag = .true.
     endif
     if(dust_kappa_scat(ilam,ispec).lt.0.d0) then
        dust_kappa_scat(ilam,ispec) = 0.d0
        flag = .true.
     endif
  enddo
  if(flag) then
     write(stdo,*) 'WARNING: Found zero absorption opacity or negative scattering opacity. Fixed it.'
  endif
  !
  ! Now integrate the Z11 matrix element over 4*pi and test if it is indeed
  ! equal to kappa_scat. At the same time we compute the cumulative 
  ! arrays of Z and the g=<mu> factor at each frequency. The former
  ! are used for the Monte Carlo module. The latter is used for the 
  ! Modified Random Walk method of Min et al. and Robitaille.
  !
  allocate(kappa_scat(1:nlam),gfactor(1:nlam))
  errormax = 0.d0
  do ilam=1,nlam
     call polarization_total_scattering_opacity(     &
          dust_kappa_arrays(ispec)%nmu,              &
          dust_kappa_arrays(ispec)%mu(:),            &
          dust_kappa_arrays(ispec)%zscat(:,ilam,1),  &
          dust_kappa_arrays(ispec)%zscat(:,ilam,2),  &
          dust_kappa_arrays(ispec)%zscat(:,ilam,3),  &
          dust_kappa_arrays(ispec)%zscat(:,ilam,4),  &
          dust_kappa_arrays(ispec)%zscat(:,ilam,5),  &
          dust_kappa_arrays(ispec)%zscat(:,ilam,6),  &
          kappa_scat(ilam),                          &
          gfactor(ilam))
     !
     ! Test kappa_scat
     !
     if(kappa_scat(ilam).lt.0.d0) then
        write(stdo,*) 'ERROR while reading dustkapscat_xxx.inp file:'
        write(stdo,*) '   Calculated negative scattering opacity'
        write(stdo,*) '   from scattering matrix.'
        stop
     endif
     if(kappa_scat(ilam).gt.0.d0) then
        error = abs(dust_kappa_arrays(ispec)%kappa_s(ilam)/kappa_scat(ilam)-1.d0)
        errormax = max(errormax,error)
        !
        ! If the error is really too large, then stop
        !
        if(error.gt.1d-1) then
           write(stdo,*) 'ERROR: The scattering opacities kappa_s in the file ',trim(filename)
           write(stdo,*) '       are inconsistent with the scattering matrix elements in the'
           write(stdo,*) '       same file. From the integration of the scattering matrix I'
           write(stdo,*) '       get, at lambda = ',2.9979d14/dust_kappa_arrays(ispec)%freq(ilam)
           write(stdo,*) '       kappa_scat = ',kappa_scat(ilam),', while in the file it says'
           write(stdo,*) '       kappa_scat = ',dust_kappa_arrays(ispec)%kappa_s(ilam)
           stop
        endif
     endif
     !
     ! Test the g-factor
     !
     if(kappa_scat(ilam).gt.0.d0) then
        error = abs(dust_kappa_arrays(ispec)%gfactor(ilam)-gfactor(ilam))
        errormax = max(errormax,error)
        !
        ! If the error is really too large, then stop
        !
        if(error.gt.1d-1) then
           write(stdo,*) 'ERROR: The g=<cos(theta)> for scattering in the file ',trim(filename)
           write(stdo,*) '       are inconsistent with the scattering matrix elements in the'
           write(stdo,*) '       same file. From the integration of the scattering matrix I'
           write(stdo,*) '       get, at lambda = ',2.9979d14/dust_kappa_arrays(ispec)%freq(ilam)
           write(stdo,*) '       g = ',gfactor(ilam),', while in the file it says'
           write(stdo,*) '       g = ',dust_kappa_arrays(ispec)%gfactor(ilam)
           stop
        endif
     endif
     !
  enddo
  !
  ! If the error is large, but not too extreme, just warn and correct
  !
  if(errormax.gt.1d-4) then
     write(stdo,*) 'Warning: The scattering opacities kappa_s in the file ',trim(filename)
     write(stdo,*) '         have relative differences of up to ',errormax,' with the '
     write(stdo,*) '         scattering matrix elements in the same file. I will correct kappa_s.'
     do ilam=1,nlam
        dust_kappa_arrays(ispec)%kappa_s(ilam) = kappa_scat(ilam)
     enddo
  endif
  !
  ! Deallocate the temporary array
  !
  deallocate(kappa_scat,gfactor,angle)
  !
  ! To avoid confusion when using this new version of RADMC-3D that does no
  ! longer write out the .used files:
  !
  filenamecheck = filename(1:filename_len)//".used"
  inquire(file=filenamecheck,exist=fex)
  if(fex) then
     write(stdo,*) 'CAREFUL: In version 0.30 and later of RADMC-3D the file ',filenamecheck
     write(stdo,*) '         is no longer written out. To avoid confusion: please remove it.'
     stop
  endif
  !
  ! Return
  !
  return
  !
  ! Error handling
  !
101 continue
  write(stdo,*) 'ERROR: Could not open file ',filename(1:filename_len)
  stop
  !
end subroutine read_dustkapscatmat_file


!-------------------------------------------------------------------
! READ ANGULAR DEPENDENCE OF ABSORPTION OPACITY FOR ALIGNED GRAINS
!-------------------------------------------------------------------
subroutine read_dustalign_angdep(ispec,filename)
  implicit none
  integer :: ispec,filename_len,iformat,nlam,nmu,ierr,imu,ilam,i,iline,k,imuface
  character*80 :: filename,string
  doubleprecision :: dum,dmu,dum2(1:2)
  logical :: fex
  !
  ! Check if dust_kappa_arrays is present
  !
  if(.not.allocated(dust_kappa_arrays)) then 
     write(stdo,*) 'ERROR: In version 0.32 and higher the dust_kappa_arrays MUST be allocated.'
     stop
  endif
  !
  ! Check if the file exists
  !
  inquire(file=filename,exist=fex)
  if(.not.fex) then
     write(stdo,*) 'ERROR: Could not find file ',trim(filename)
     stop
  endif

  !
  ! Open the file 
  !
  filename_len = len_trim(filename)
  open(unit=1,file=filename,status='old')
  !
  ! Skip all lines at the beginning of file that start with one of the
  ! symbols ";", "#" or "!". These are meant as comments. 300 such lines
  ! is the maximum. 
  !
  do iline=1,300
     read(1,'(A158)',end=10) string
     i=1
     do k=1,100
        if(string(i:i).ne.' ') goto 30
        i=i+1
     enddo
     write(stdo,*) 'ERROR While reading ',filename(1:filename_len)
     write(stdo,*) '  Detected line full of spaces...'
30   continue
     if((string(i:i).ne.';').and.(string(i:i).ne.'#').and.(string(i:i).ne.'!')) then
        goto 20
     endif
  enddo
  stop 8901
10 continue
  write(stdo,*) 'ERROR while reading file ',filename(1:filename_len)
  write(stdo,*) 'Did not find data before end of file.'
  stop 8902
20 continue
  !
  ! First read the format number. Do so from the above read string.
  !
  read(string,*) iformat
  if(iformat.ne.1) then
     write(stdo,*) 'ERROR: Format number ',iformat,' of file ',trim(filename),&
          ' not known.'
     stop
  endif
  !
  ! Read the number of frequencies
  !
  read(1,*) nlam
  if(dust_kappa_arrays(ispec)%nrfreq.ne.nlam) then
     write(stdo,*) 'ERROR while reading file',filename(1:filename_len)
     write(stdo,*) 'Nr of frequencies not the same as in the corresponding scattering matrix file.'
     stop
  endif
  !
  ! Read the number of scattering angles
  !
  read(1,*) nmu
  dust_kappa_arrays(ispec)%namu = nmu
  !
  ! Make a check
  !
  if(nmu.lt.2) then
     write(stdo,*) 'ERROR: RADMC-3D cannot accept alignment ratios for fewer than 2 angles'
     write(stdo,*) '       in file ',trim(filename)
     stop
  endif
  !
  ! Allocate arrays
  !
  allocate(dust_kappa_arrays(ispec)%alignmu(1:nmu),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate alignment arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%alignorth(1:nmu,1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate alignment arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%alignpara(1:nmu,1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate alignment arrays'
     stop
  endif
  !
  ! Check the wavelength grid, which must be the same as the one
  ! read in for the scattering matrix
  !
  do ilam=1,nlam
     read(1,*) dum
     dum = 2.9979d14/dum
     if(abs(dust_kappa_arrays(ispec)%freq(ilam)/dum-1.d0).gt.1d-4) then
        write(stdo,*) 'ERROR while reading file ',filename(1:filename_len)
        write(stdo,*) 'Wavelength grid not the same as for the scattering matrix file'
        stop
     endif
  enddo
  !
  ! Read the angular grid 
  !
  do imu=1,nmu
     read(1,*) dum
     dust_kappa_arrays(ispec)%alignmu(imu) = cos(dum*pi/180.d0)
  enddo
  if(dust_kappa_arrays(ispec)%alignmu(1).lt.dust_kappa_arrays(ispec)%alignmu(nmu)) then
     if(dust_kappa_arrays(ispec)%alignmu(1).lt.-1d-8) then
        write(stdo,*) 'ERROR while reading file ',filename(1:filename_len)
        write(stdo,*) '   angular grid does not span cos(theta) in [0,1]'
        stop
     endif
     dust_kappa_arrays(ispec)%alignmu(1)   = 0.d0
     dust_kappa_arrays(ispec)%alignmu(nmu) = 1.d0
     imuface = nmu
  else
     if(dust_kappa_arrays(ispec)%alignmu(nmu).lt.-1d-8) then
        write(stdo,*) 'ERROR while reading file ',filename(1:filename_len)
        write(stdo,*) '   angular grid does not span cos(theta) in [0,1]'
        stop
     endif
     dust_kappa_arrays(ispec)%alignmu(nmu) = 0.d0
     dust_kappa_arrays(ispec)%alignmu(1)   = 1.d0
     imuface = 1
  endif
  !
  ! Read the alpha_abs_orth and alpha_abs_para
  !
  do ilam=1,nlam
     do imu=1,nmu
        read(1,*) dum2
        dust_kappa_arrays(ispec)%alignorth(imu,ilam) = dum2(1)
        dust_kappa_arrays(ispec)%alignpara(imu,ilam) = dum2(2)
        if((dum2(1).lt.0.d0).or.(dum2(2).lt.0.d0)) then
           write(stdo,*) 'ERROR: Orthogonal and parallel opacities cannot be negative!'
           stop
        endif
        if(imu.eq.imuface) then
           !
           ! Check that for face-on view the orth and para are
           ! equal. If not, then the basic assumption of axially
           ! symmetric (or fast-spinning) grains is not fulfilled.
           !
           if(abs((dum2(2)-dum2(1))/(dum2(2)+dum2(1)+1d-60)).gt.1d-8) then
              write(stdo,*) 'ERROR while reading alignment data:'
              write(stdo,*) '   A grain seen along the alignment axis must have zero polarization.'
              write(stdo,*) '   Even if it is a prolate grain, because then it is spinning fast,'
              write(stdo,*) '   so that it always should average out. In the file I now read in,'
              write(stdo,*) '   however, the orthogonal and parallel opacity coefficients are not'
              write(stdo,*) '   equal, even for face-on view (seen long alignment axis).'
              stop
           endif
        endif
     enddo
  enddo
  !
  ! Now renormalize the align ratio such that when we integrate the
  ! ratio over random angles, the result is 1. This ensures that
  ! when we place the grains at random orientations, the resulting
  ! absorption opacity is equal to the one we already read in from
  ! the opacity (kappascatmat) file. 
  !
  do ilam=1,nlam
     dum = 0.d0
     do imu=2,nmu
        dmu = abs( dust_kappa_arrays(ispec)%alignmu(imu) -    &
                   dust_kappa_arrays(ispec)%alignmu(imu-1) )
        dum = dum + dmu * 0.5d0 * (                           &
              dust_kappa_arrays(ispec)%alignorth(imu,ilam) +  &
              dust_kappa_arrays(ispec)%alignorth(imu-1,ilam) )
        dum = dum + dmu * 0.5d0 * (                           &
              dust_kappa_arrays(ispec)%alignpara(imu,ilam) +  &
              dust_kappa_arrays(ispec)%alignpara(imu-1,ilam) )
     enddo
     dum = 0.5d0 * dum   ! Because int 0.5*(kap_h+kap_v) dmu is normalized
     dust_kappa_arrays(ispec)%alignorth(:,ilam) =             &
          dust_kappa_arrays(ispec)%alignorth(:,ilam) / dum
     dust_kappa_arrays(ispec)%alignpara(:,ilam) =             &
          dust_kappa_arrays(ispec)%alignpara(:,ilam) / dum
  enddo
  !
  ! Close the file
  !
  close(1)
  !
end subroutine read_dustalign_angdep


!-------------------------------------------------------------------
!                   READ OPTICAL CONSTANT FILE
!-------------------------------------------------------------------
subroutine read_optconst_file(ispec,filename)
  implicit none
  character*80 :: filename
  integer :: ispec
  character*160 :: string
  character*80 :: filenamecheck
  integer :: iformat,nlam,ilam,ierr,iline,filename_len,i,k
  double precision :: dummy3(1:3),dum
  double precision :: numin,numax,frmin,frmax,sgn
  logical :: flag,fex
  !
  ! First check if the frequency grid is present
  !
  if((.not.allocated(freq_nu)).or.(freq_nr.le.0)) then
     write(stdo,*) 'ERROR: While reading dust opacity: frequency array not yet set.'
     stop
  endif
  !
  ! Remove any pre-existing opacity, if present
  !
  call dust_remove_opacity(ispec)
  !
  ! Open the file containing the n and k values
  !
  filename_len = len_trim(filename)
  open(unit=1,file=filename,status='old',err=101)
  !
  ! Skip all lines at the beginning of file that start with one of the
  ! symbols ";", "#" or "!". These are meant as comments. 300 such lines
  ! is the maximum. 
  !
  do iline=1,300
     read(1,'(A158)',end=10) string
     i=1
     do k=1,100
        if(string(i:i).ne.' ') goto 30
        i=i+1
     enddo
     write(stdo,*) 'ERROR While reading ',filename(1:filename_len)
     write(stdo,*) '  Detected line full of spaces...'
30   continue
     if((string(i:i).ne.';').and.(string(i:i).ne.'#').and.(string(i:i).ne.'!')) then
        goto 20
     endif
  enddo
  stop 8901
10 continue
  write(stdo,*) 'ERROR while reading file ',filename(1:filename_len)
  write(stdo,*) 'Did not find data before end of file.'
  stop 8902
20 continue
  !
  ! First read the format number. Do so from the above read string.
  !
  read(string,*) iformat
  !
  ! Read the specific weight of this material in g/cm^3
  !
  read(1,*) dum
  dust_kappa_arrays(ispec)%sweight = dum
  !
  ! Read the number of frequencies
  !
  read(1,*) nlam
  dust_kappa_arrays(ispec)%nrfreq = nlam
  !
  ! Allocate arrays
  !
  if(.not.allocated(dust_kappa_arrays)) then 
     write(stdo,*) 'ERROR: In version 0.30 and higher the dust_kappa_arrays MUST be allocated.'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%freq(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%n(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%k(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  !
  ! Read the data
  !
  do ilam=1,nlam
     read(1,*) dummy3
     dust_kappa_arrays(ispec)%freq(ilam)  = 2.9979d14/dummy3(1)
     dust_kappa_arrays(ispec)%n(ilam)     = dummy3(2)
     dust_kappa_arrays(ispec)%k(ilam)     = dummy3(3)
  enddo
  !
  ! Close the file
  ! 
  close(1)
  !
  ! Check that the lambda array is monotonically increasing/decreasing
  !
  if(dust_kappa_arrays(ispec)%freq(2).lt.dust_kappa_arrays(ispec)%freq(1)) then
     sgn = -1
  elseif(dust_kappa_arrays(ispec)%freq(2).gt.dust_kappa_arrays(ispec)%freq(1)) then
     sgn = 1
  else
     write(stdo,*) 'ERROR while reading opacity file ',filename(1:filename_len)
     write(stdo,*) '    Cannot have twice the same lambda.'
     stop
  endif
  do ilam=3,nlam
     dum = (dust_kappa_arrays(ispec)%freq(ilam)-dust_kappa_arrays(ispec)%freq(ilam-1))*sgn
     if(dum.eq.0.d0) then
        write(stdo,*) 'ERROR while reading opacity file ',filename(1:filename_len)
        write(stdo,*) '    Cannot have twice the same lambda.'
        stop
     elseif(dum.lt.0.d0) then
        write(stdo,*) 'ERROR while reading opacity file ',filename(1:filename_len)
        write(stdo,*) '    Wavelength array is not monotonic.'
        stop
     endif
  enddo
  !
  ! Give a warning if range of nu does not include range of freq_nu
  !
  if(freq_nu(1).lt.freq_nu(freq_nr)) then
     frmin = freq_nu(1)
     frmax = freq_nu(freq_nr)
  else
     frmin = freq_nu(freq_nr)
     frmax = freq_nu(1)
  endif
  if(dust_kappa_arrays(ispec)%freq(1).lt.dust_kappa_arrays(ispec)%freq(nlam)) then
     numin = dust_kappa_arrays(ispec)%freq(1)
     numax = dust_kappa_arrays(ispec)%freq(nlam)
  else
     numin = dust_kappa_arrays(ispec)%freq(nlam)
     numax = dust_kappa_arrays(ispec)%freq(1)
  endif
  if((frmin.lt.numin).or.(frmax.gt.numax)) then
     write(stdo,*) 'Note: Opacity file ',filename(1:filename_len),' does not cover the '
     write(stdo,*) '      complete wavelength range of the model (i.e. the range'
     write(stdo,*) '      in the file frequency.inp or wavelength_micron.inp).'
     write(stdo,201) '       Range of model:     lambda = [',2.9979d14/frmax, &
                     ',',2.9979d14/frmin,']'
     write(stdo,201) '       Range in this file: lambda = [',2.9979d14/numax, &
                     ',',2.9979d14/numin,']'
201  format(a37,e9.2,a1,e9.2,a1) 
     write(stdo,*) '      Simple extrapolations are used to extend this range.'
  endif
  !
  ! Return
  !
  return
  !
  ! Error handling
  !
101 continue
  write(stdo,*) 'ERROR: Could not open file ',filename(1:filename_len)
  stop
  !
end subroutine read_optconst_file



!-------------------------------------------------------------------
!          COMPUTE DUST OPACITY FROM OPTICAL CONSTANTS
!
! This is another new style of the dustopacity files. Here the 
! lambda, n and k are read from a file and using CDE or MIE the
! dust opacity is computed on-the-fly. These n and k files can
! have their own (preferably high resolution) frequency/wavelength
! grid. This is then remapped onto the internally used frequency
! grid, i.e. the one listed in frequency.inp or wavelength_micron.inp
!
! NOTE: AT THIS MOMENT THE SCATTERING IS STILL ASSUMED TO BE 
!       ISOTROPIC IF YOU USE THIS SUBROUTINE.
!
!-------------------------------------------------------------------
subroutine make_dust_opacity(ispec,agrain,algstring)
  implicit none
  integer :: ispec
  character*160 :: string
  character*80 :: algstring
  character*80 :: filename,base,ext
  integer :: iformat,nlam,ilam,ierr,i,k
  double precision :: dummy3(1:3),dum,geocrosssec,ttt,geokappa,mass
  double precision :: numin,numax,frmin,frmax,sgn,agrain,lambda
  double precision :: qe,qs,qa,ge,ce,cs,ca
  logical :: flag
  !
  ! Check if the optical constants have already been loaded
  !
  if((.not.associated(dust_kappa_arrays(ispec)%n)).or.(.not.associated(dust_kappa_arrays(ispec)%n))) then
     write(stdo,*) 'ERROR in dust_module: Wanting to compute opacity using Mie/CDE, but optical'
     write(stdo,*) '      constants are not loaded yet.'
     stop
  endif
  if(dust_kappa_arrays(ispec)%sweight.le.0.d0) then
     write(stdo,*) 'ERROR in dust_module: Material density is 0...'
     stop
  endif
  !
  ! Set nlam
  !
  nlam = dust_kappa_arrays(ispec)%nrfreq
  !
  ! If necessary, deallocate opacities
  !
  if(associated(dust_kappa_arrays(ispec)%kappa_a)) &
       deallocate(dust_kappa_arrays(ispec)%kappa_a)
  if(associated(dust_kappa_arrays(ispec)%kappa_s)) &
       deallocate(dust_kappa_arrays(ispec)%kappa_s)
  if(associated(dust_kappa_arrays(ispec)%gfactor)) &
       deallocate(dust_kappa_arrays(ispec)%gfactor)
  if(associated(dust_kappa_arrays(ispec)%mu)) &
       deallocate(dust_kappa_arrays(ispec)%mu)
  if(associated(dust_kappa_arrays(ispec)%zscat)) &
       deallocate(dust_kappa_arrays(ispec)%zscat)
  !
  ! Now allocate the opacity arrays
  !
  allocate(dust_kappa_arrays(ispec)%kappa_a(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%kappa_s(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  allocate(dust_kappa_arrays(ispec)%gfactor(1:nlam),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in dust module: Could not allocate kappas arrays'
     stop
  endif
  !
  ! Compute cross section of the grain, assuming it to be a compact sphere
  !
  geocrosssec = pi * agrain**2
  !
  ! Compute the mass of the grain
  !
  mass = (4.d0*pi/3.d0) * dust_kappa_arrays(ispec)%sweight * agrain**3
  !
  ! Store this mass, or check with previously computed grain mass
  !
  if(dust_mgrain(ispec).eq.0.d0) then
     dust_mgrain(ispec) = mass
  else
     if(abs(1.d0-mass/dust_mgrain(ispec)).gt.1d-4) then
        write(stdo,*) 'ERROR in dust_module: The grain mass '
        write(stdo,*) '      specified previously is not equal'
        write(stdo,*) '      to the one from dustopac.inp.'
        write(stdo,*) '      ispec           = ',ispec
        write(stdo,*) '      mgrain previous = ',dust_mgrain(ispec)
        write(stdo,*) '      mgrain now      = ',mass
        write(stdo,*) '      Maybe material weight different?'
        write(stdo,*) '      Maybe porosity assumed?'
        stop
     endif
  endif
  !
  ! Compute the geometric kappa
  !
  geokappa = geocrosssec / mass
  !
  ! Now make the opacity on the lab grid
  ! 
  if(algstring(1:3).eq.'MIE') then
     !
     ! Use Mie Theory
     !
     do ilam=1,nlam
        !
        ! Compute the absorption and scattering opacity [cm^2/gram_of_dust]
        !
        ttt = 0.d0
        lambda = 2.9979d14/dust_kappa_arrays(ispec)%freq(ilam)
        call q_mie(dust_kappa_arrays(ispec)%n(ilam),dust_kappa_arrays(ispec)%k(ilam), &
                   lambda,agrain*1d4,ttt,qe,qs,qa,ge)
        dust_kappa_arrays(ispec)%kappa_a(ilam) = qa * geokappa
        dust_kappa_arrays(ispec)%kappa_a(ilam) = qs * geokappa
        !
        ! Now the Henyey-Greenstein g-factor
        !
        dust_kappa_arrays(ispec)%gfactor(ilam) = abs(ge)
        !
     enddo
  elseif(algstring(1:3).eq.'CDE') then
     !
     ! Use CDE Theory
     !
     do ilam=1,nlam
        !
        ! Compute the absorption and scattering cross section
        !
        lambda = 2.9979d14/dust_kappa_arrays(ispec)%freq(ilam)
        call q_cde(dust_kappa_arrays(ispec)%n(ilam),dust_kappa_arrays(ispec)%k(ilam), &
                   lambda,agrain*1d4,ce,cs,ca,ge)
        !
        ! Convert cross section from micron^2 -> cm^2
        !
        ce = ce * 1.d-8
        cs = cs * 1.d-8
        ca = ca * 1.d-8
        !
        ! Convert to kappa [cm^2/gram_of_dust]
        !
        dust_kappa_arrays(ispec)%kappa_a(ilam) = ca / mass
        dust_kappa_arrays(ispec)%kappa_s(ilam) = cs / mass
        !
        ! Now the Henyey-Greenstein g-factor
        !
        dust_kappa_arrays(ispec)%gfactor(ilam) = abs(ge)
        !
     enddo
  else
     write(stdo,*) 'ERROR: Mode ',algstring(1:len_trim(algstring)),  &
          ' not known as a method for opacity computation.'
     stop
  endif
  !
  ! Check that the opacities are not negative
  !
  flag = .false.
  do ilam=1,nlam
     if(dust_kappa_arrays(ispec)%kappa_a(ilam).lt.1d-90) then
        dust_kappa_arrays(ispec)%kappa_a(ilam) = 1d-90
        flag = .true.
     endif
     if(dust_kappa_arrays(ispec)%kappa_s(ilam).lt.0.d0) then
        dust_kappa_arrays(ispec)%kappa_s(ilam) = 0.d0
        flag = .true.
     endif
  enddo
  if(flag) then
     write(stdo,*) 'WARNING: Found zero absorption opacity or negative scattering opacity. Fixed it.'
  endif
!  !
!  ! Write a ".used" output file
!  !
!  base='dustkappa_'
!  ext ='.used'
!  call make_indexed_filename(base,ispec,ext,filename)
!  open(unit=1,file=filename)
!  write(1,*) '# This is the opacity from ',trim(filename),&
!             ' computed from optical constants.'
!  write(1,*) '# This is not an input file. It is only meant for the user to',&
!             ' check if the Mie/CDE computation of the opacities went OK. '
!  write(1,*) 3
!  write(1,*) nlam
!  do ilam=1,nlam
!     write(1,*) 2.9979d14/dust_kappa_arrays(ispec)%freq(ilam),dust_kappa_arrays(ispec)%kappa_a(ilam),&
!                dust_kappa_arrays(ispec)%kappa_s(ilam),dust_kappa_arrays(ispec)%gfactor(ilam)
!  enddo
!  close(1)
  !
  ! Set the scattering_mode_def to 2
  !
  scattering_mode_def = min(2,scattering_mode_max)
  !
  ! Return
  !
  return
  !
end subroutine make_dust_opacity


!-------------------------------------------------------------------
!              FIND DUST OPACITY BETWEEN FREQ POINTS
!-------------------------------------------------------------------
function find_dust_kappa_interpol(freq,ispec,temp,iabs,iscat,igfact)
  implicit none
  double precision :: temp,find_dust_kappa_interpol,freq,eps
  integer :: inu,ispec,iabs,iscat,igfact,nfr,inumin,inumax
  double precision :: margin
  !
  margin = 1d-4
  find_dust_kappa_interpol = 0.d0
  !
  ! First find out if we have the original high-resolution opacities
  ! in memory
  !
  if(allocated(dust_kappa_arrays)) then
     if(dust_kappa_arrays(ispec)%nrfreq.gt.0) then
        if(.not.associated(dust_kappa_arrays(ispec)%freq)) stop 8800
        if(.not.associated(dust_kappa_arrays(ispec)%kappa_a)) stop 8800
        if(.not.associated(dust_kappa_arrays(ispec)%kappa_s)) stop 8800
        if(.not.associated(dust_kappa_arrays(ispec)%gfactor)) stop 8800
        !
        ! Yes, there is something in memory. Now check if freq is 
        ! within the domain of this opacity table. If not, then we 
        ! will try nevertheless with the global array.
        !
        nfr  = dust_kappa_arrays(ispec)%nrfreq
        if((freq-dust_kappa_arrays(ispec)%freq(nfr))*                   &
           (freq-dust_kappa_arrays(ispec)%freq(1)).lt.0.d0) then
           !
           ! Yes, the freq lies within the boundaries.
           !
           call hunt(dust_kappa_arrays(ispec)%freq,                     &
                     dust_kappa_arrays(ispec)%nrfreq,freq,inu)
           if((inu.lt.1).or.(inu.ge.dust_kappa_arrays(ispec)%nrfreq)) stop 8201
           eps = (freq-dust_kappa_arrays(ispec)%freq(inu)) /            &
                 (dust_kappa_arrays(ispec)%freq(inu+1)-                 &
                  dust_kappa_arrays(ispec)%freq(inu))
           if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 8202
           if(iabs.ne.0) then
              find_dust_kappa_interpol = find_dust_kappa_interpol +     & 
                   (1.d0-eps) * dust_kappa_arrays(ispec)%kappa_a(inu) + &
                          eps * dust_kappa_arrays(ispec)%kappa_a(inu+1)
           endif
           if(iscat.ne.0) then
              find_dust_kappa_interpol = find_dust_kappa_interpol +     & 
                   (1.d0-eps) * dust_kappa_arrays(ispec)%kappa_s(inu) + &
                          eps * dust_kappa_arrays(ispec)%kappa_s(inu+1)
           endif
           if(igfact.ne.0) then
              find_dust_kappa_interpol =                                & 
                   (1.d0-eps) * dust_kappa_arrays(ispec)%gfactor(inu) + &
                          eps * dust_kappa_arrays(ispec)%gfactor(inu+1)
           endif
           !
           ! Return
           !
           return
        endif
     endif
  endif
  !
  ! No original array in memory, or the original opacity was read from
  ! dustopac_*.inp, which is mapped on the same grid as wavelength_micron.inp,
  ! or the freq is out of range of the original grid.
  ! So use the global grid.
  !
  ! First check if we are at, or just over the edge of the grid. For security, because
  ! of potential round-off errors, we smoothly switch the opacity off here instead of
  ! using a hard switch-off. 
  !
  inu = -1
  if(freq_nu(freq_nr).gt.freq_nu(1)) then
     inumax = freq_nr
     inumin = 1
  else
     inumax = 1
     inumin = freq_nr
  endif
  if(freq.le.freq_nu(inumin)) then
     if(freq.le.(1.d0-margin)*freq_nu(inumin)) then
        return
     else
        eps = (freq-(1.d0-margin)*freq_nu(inumin))/(margin*freq_nu(inumin))
        inu = inumin
     endif
  endif
  if(freq.ge.freq_nu(inumax)) then
     if(freq.ge.(1.d0+margin)*freq_nu(inumax)) then
        return
     else
        eps = ((1.d0+margin)*freq_nu(inumax)-freq)/(margin*freq_nu(inumax))
        inu = inumax
     endif
  endif
  if(inu.gt.0) then
     !
     ! Yes, we are in this marginal region. 
     !
     if(eps.lt.0.d0) eps=0.d0     ! f-(1-margin)*f is not perfectly margin*f
     if(eps.gt.1.d0) eps=1.d0     ! f-(1-margin)*f is not perfectly margin*f
     if(iabs.ne.0) then
        find_dust_kappa_interpol = find_dust_kappa_interpol + & 
             eps * dust_kappa_abs(inu,ispec)
     endif
     if(iscat.ne.0) then
        find_dust_kappa_interpol = find_dust_kappa_interpol + &
             eps * dust_kappa_scat(inu,ispec)
     endif
     if(igfact.ne.0) then
        find_dust_kappa_interpol =                            &
             eps * dust_gfactor(inu,ispec)
     endif
     !
     ! Return
     !
     return
  endif
  !
  ! We are still inside the global frequency grid
  !
  call hunt(freq_nu,freq_nr,freq,inu)
  if((inu.lt.1).or.(inu.ge.freq_nr)) stop 8301
  eps = (freq-freq_nu(inu)) / (freq_nu(inu+1)-freq_nu(inu))
  if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 8302
  if(iabs.ne.0) then
     find_dust_kappa_interpol = find_dust_kappa_interpol + & 
          (1.d0-eps) * dust_kappa_abs(inu,ispec) +         &
                 eps * dust_kappa_abs(inu+1,ispec)
  endif
  if(iscat.ne.0) then
     find_dust_kappa_interpol = find_dust_kappa_interpol + &
          (1.d0-eps) * dust_kappa_scat(inu,ispec) +        &
                 eps * dust_kappa_scat(inu+1,ispec)
  endif
  if(igfact.ne.0) then
     find_dust_kappa_interpol =                            &
          (1.d0-eps) * dust_gfactor(inu,ispec) +             &
                 eps * dust_gfactor(inu+1,ispec)
  endif
  return
end function find_dust_kappa_interpol


!-------------------------------------------------------------------
!          FIND THE Z MATRIX (SCATTERING MATRIX) ELEMENTS
!-------------------------------------------------------------------
subroutine find_dust_zmatrix_interpol(freq,mu,ispec,zmat,extrapol)
  implicit none
  integer :: ispec
  double precision :: freq,mu,zmat(1:6)
  double precision :: epsf,epsf1,epsm,epsm1
  double precision :: fact,fprev,fmin,fmax,valmn,valprev
  integer :: inu,imu,iz,imin,imax,iprev
  logical :: extrapol
  !
  ! Do some basic checks
  !  
  if(.not.associated(dust_kappa_arrays(ispec)%freq)) stop 8800
  if(.not.associated(dust_kappa_arrays(ispec)%zscat)) stop 8800
  !
  ! Nullify stuff
  !
  zmat(1:6) = 0.d0
  !
  ! Find the mu in the mu grid of this dust species
  !
  call hunt(dust_kappa_arrays(ispec)%mu,                        &
            dust_kappa_arrays(ispec)%nmu,mu,imu)
  if(imu.lt.1) imu=1
  if(imu.ge.dust_kappa_arrays(ispec)%nmu) imu=dust_kappa_arrays(ispec)%nmu-1
  epsm = (mu-dust_kappa_arrays(ispec)%mu(imu)) /                &
         (dust_kappa_arrays(ispec)%mu(imu+1)-                   &
          dust_kappa_arrays(ispec)%mu(imu))
  if((epsm.lt.0.d0).or.(epsm.gt.1.d0)) stop 8202
  epsm1 = 1.d0-epsm
  !
  ! Find the frequency in the frequency grid for this dust
  ! species
  !
  call hunt(dust_kappa_arrays(ispec)%freq,                      &
            dust_kappa_arrays(ispec)%nrfreq,freq,inu)
  if((inu.ge.1).and.(inu.lt.dust_kappa_arrays(ispec)%nrfreq)) then
     !
     ! We are within the frequency range of the stored opacity table
     !
     epsf = (freq-dust_kappa_arrays(ispec)%freq(inu)) /            &
            (dust_kappa_arrays(ispec)%freq(inu+1)-                 &
             dust_kappa_arrays(ispec)%freq(inu))
     if((epsf.lt.0.d0).or.(epsf.gt.1.d0)) stop 8202
     epsf1 = 1.d0-epsf
     !
     ! Do the interpolation
     !
     do iz=1,6
        zmat(iz) = epsf1 * (                                                  &
                   epsm1 * dust_kappa_arrays(ispec)%zscat(imu,inu,iz) +       &
                    epsm * dust_kappa_arrays(ispec)%zscat(imu+1,inu,iz) ) +   &
                    epsf * (                                                  &
                   epsm1 * dust_kappa_arrays(ispec)%zscat(imu,inu+1,iz) +     &
                    epsm * dust_kappa_arrays(ispec)%zscat(imu+1,inu+1,iz) )
     enddo
  else
     !
     ! We are outside of the frequency range
     !
     if(extrapol) then
        !
        ! We have to do extrapolations
        !
        fmin = dust_kappa_arrays(ispec)%freq(1)
        fmax = dust_kappa_arrays(ispec)%freq(dust_kappa_arrays(ispec)%nrfreq)
        imin = 1
        imax = dust_kappa_arrays(ispec)%nrfreq
        if(fmax.lt.fmin) then
           fmin = dust_kappa_arrays(ispec)%freq(dust_kappa_arrays(ispec)%nrfreq)
           fmax = dust_kappa_arrays(ispec)%freq(1)
           imin = dust_kappa_arrays(ispec)%nrfreq
           imax = 1
        endif
        if(freq.le.fmin) then
           !
           ! Long wavelength end: Double-logarithmic extrapolation
           !
           if(imin.eq.1) then
              iprev = 2
           else
              iprev = imin-1
           endif
           fprev   = dust_kappa_arrays(ispec)%freq(iprev)
           valmn   = epsm1 * dust_kappa_arrays(ispec)%zscat(imu,imin,1) + &
                      epsm * dust_kappa_arrays(ispec)%zscat(imu+1,imin,1)
           valprev = epsm1 * dust_kappa_arrays(ispec)%zscat(imu,iprev,1) + &
                      epsm * dust_kappa_arrays(ispec)%zscat(imu+1,iprev,1)
           if((valmn.gt.0.d0).and.(valprev.gt.0.d0)) then
              fact = (valprev/valmn)**((log(freq)-log(fmin))/   &
                                       (log(fprev)-log(fmin)))
              do iz=1,6
                 zmat(iz) = fact * (                                        &
                    epsm1 * dust_kappa_arrays(ispec)%zscat(imu,imin,iz) +   &
                     epsm * dust_kappa_arrays(ispec)%zscat(imu+1,imin,iz) )
              enddo
           else
              zmat(1:6) = 0.d0
           endif
        elseif(freq.ge.fmax) then
           !
           ! Short wavelength end: Constant extrapolation
           !
           do iz=1,6
              zmat(iz) = epsm1 * dust_kappa_arrays(ispec)%zscat(imu,imax,iz) + &
                          epsm * dust_kappa_arrays(ispec)%zscat(imu+1,imax,iz)
           enddo
        else
           write(stdo,*) 'INTERNAL ERROR: Inside freq domain yet outside... Warn author.'
           stop
        endif
     endif
  endif
end subroutine find_dust_zmatrix_interpol


!-------------------------------------------------------------------
!          FIND DUST ALIGNMENT FACTORS BETWEEN FREQ POINTS
!-------------------------------------------------------------------
subroutine find_dust_alignfact_interpol(freq,mu,ispec,iabs,iscat,  &
                                        extrapol,orth,para)
  implicit none
  double precision :: temp,find_dust_kappa_interpol,freq,orth,para,mu
  double precision :: epsfreq,epsmu
  integer :: inu,ispec,iabs,iscat,nfr,nmu,imin,imax,imu
  double precision :: fmin,fmax
  logical :: extrapol
  !
  if(iscat.ne.0) then
     write(stdo,*) 'ERROR: Aligned scattering not yet implemented.'
     stop
  endif
  if(iabs.ne.1) then
     write(stdo,*) 'ERROR: Nothing to do.'
     stop
  endif
  !
  ! First find out if we have the original high-resolution opacities
  ! in memory
  !
  if(.not.allocated(dust_kappa_arrays)) then
     write(stdo,*) 'ERROR: No original dust arrays in memory...'
     stop
  endif
  !
  ! Check if mu in range
  !
  if((mu.lt.0.d0).or.(mu.gt.1.d0)) stop 3029
  !
  ! Check if this information is available for this dust species
  ! 
  if((dust_kappa_arrays(ispec)%nrfreq.gt.0).and. &
     (dust_kappa_arrays(ispec)%namu.gt.0)) then
     if(.not.associated(dust_kappa_arrays(ispec)%freq)) stop 8800
     if(.not.associated(dust_kappa_arrays(ispec)%alignmu)) stop 8800
     if(.not.associated(dust_kappa_arrays(ispec)%alignorth)) stop 8800
     if(.not.associated(dust_kappa_arrays(ispec)%alignpara)) stop 8800
     nfr  = dust_kappa_arrays(ispec)%nrfreq
     nmu  = dust_kappa_arrays(ispec)%namu
     !
     ! Find where we are in the mu grid
     ! 
     call hunt(dust_kappa_arrays(ispec)%alignmu,nmu,mu,imu)
     if(imu.eq.nmu) then
        if(mu.eq.dust_kappa_arrays(ispec)%alignmu(nmu)) then
           imu = nmu-1
        endif
     endif
     if((imu.lt.1).or.(imu.gt.nmu)) stop 8202
     epsmu = (mu-dust_kappa_arrays(ispec)%alignmu(imu)) /            &
              (dust_kappa_arrays(ispec)%alignmu(imu+1)-              &
               dust_kappa_arrays(ispec)%alignmu(imu))
     if((epsmu.lt.0.d0).or.(epsmu.gt.1.d0)) stop 8203
     !
     ! Check if we are in the frequency array domain
     !
     if((freq-dust_kappa_arrays(ispec)%freq(nfr))*                   &
        (freq-dust_kappa_arrays(ispec)%freq(1)).lt.0.d0) then
        !
        ! Yes, the freq lies within the boundaries.
        !
        call hunt(dust_kappa_arrays(ispec)%freq,nfr,freq,inu)
        if((inu.lt.1).or.(inu.ge.dust_kappa_arrays(ispec)%nrfreq)) stop 8204
        epsfreq = (freq-dust_kappa_arrays(ispec)%freq(inu)) /        &
              (dust_kappa_arrays(ispec)%freq(inu+1)-                 &
               dust_kappa_arrays(ispec)%freq(inu))
        if((epsfreq.lt.0.d0).or.(epsfreq.gt.1.d0)) stop 8205
        if(iabs.ne.0) then
           orth = (1.d0-epsmu) * (                                               &
              (1.d0-epsfreq) * dust_kappa_arrays(ispec)%alignorth(imu,inu) +     &
                     epsfreq * dust_kappa_arrays(ispec)%alignorth(imu,inu+1) ) + &
                     epsmu * (                                                   &
              (1.d0-epsfreq) * dust_kappa_arrays(ispec)%alignorth(imu+1,inu) +   &
                     epsfreq * dust_kappa_arrays(ispec)%alignorth(imu+1,inu+1) )
           para = (1.d0-epsmu) * (                                               &
              (1.d0-epsfreq) * dust_kappa_arrays(ispec)%alignpara(imu,inu) +     &
                     epsfreq * dust_kappa_arrays(ispec)%alignpara(imu,inu+1) ) + &
                     epsmu * (                                                   &
              (1.d0-epsfreq) * dust_kappa_arrays(ispec)%alignpara(imu+1,inu) +   &
                     epsfreq * dust_kappa_arrays(ispec)%alignpara(imu+1,inu+1) )
        endif
     else
        !
        ! Freq lies outside the domain, but it can lie very
        ! narrowly outside the domain, so we may want to extrapolate.
        !
        if(extrapol) then
           !
           ! We have to do extrapolations
           !
           fmin = dust_kappa_arrays(ispec)%freq(1)
           fmax = dust_kappa_arrays(ispec)%freq(dust_kappa_arrays(ispec)%nrfreq)
           imin = 1
           imax = dust_kappa_arrays(ispec)%nrfreq
           if(fmax.lt.fmin) then
              fmin = dust_kappa_arrays(ispec)%freq(dust_kappa_arrays(ispec)%nrfreq)
              fmax = dust_kappa_arrays(ispec)%freq(1)
              imin = dust_kappa_arrays(ispec)%nrfreq
              imax = 1
           endif
           if(freq.le.fmin) then
              !
              ! Long wavelength end: Constant extrapolation
              !
              inu = imin
           elseif(freq.ge.fmax) then
              !
              ! Short wavelength end: Constant extrapolation
              !
              inu = imax
           else
              write(stdo,*) 'INTERNAL ERROR: Inside freq domain yet outside... Warn author.'
              stop
           endif
           !
           ! Interpolate in mu
           !
           orth = (1.d0-epsmu) * dust_kappa_arrays(ispec)%alignorth(imu,inu) + &
                        epsmu  * dust_kappa_arrays(ispec)%alignorth(imu+1,inu)
           para = (1.d0-epsmu) * dust_kappa_arrays(ispec)%alignpara(imu,inu) + &
                        epsmu  * dust_kappa_arrays(ispec)%alignpara(imu+1,inu)
        endif
     endif
  else
     !
     ! No data available, so assume orth=1 and para=1 
     !
     orth = 1.d0
     para = 1.d0
  endif
end subroutine find_dust_alignfact_interpol


!-------------------------------------------------------------------
!                    THE MIE ROUTINE (FROM F77)
!-------------------------------------------------------------------
SUBROUTINE Q_MIE(E1,E2,LAM,RAD,T,QEX,QSC,QAB,G)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  DOUBLE PRECISION LAM,RAD,T,QEX,QSC,QAB,E1,E2,G
  !
  ! MIE THEORY EFFICIENCY FACTORS FOR SPHERICAL PARTICLES OF
  ! RADIUS 'RAD' AT WAVELENGTH 'LAM'.
  ! E=E1 + I*E2 IS THE SQUARE OF THE COMPLEX REFRACTIVE INDEX.
  ! THE REFRACTIVE INDEX IS GIVEN BY SUBROUTINE 'EPS'
  !                         
  COMPLEX*16 E,RM,Y,ZN,ZN1,ZN2,AN(700000),C,A,B,AO,RRAT,A1,ANM1,BNM1
  X=6.2831853*RAD/LAM
  E=DCMPLX(E1,-E2)
  E=E**2.
  IF(X.LT.0.001)THEN
     !
     ! USE SMALL PARTICLE FORMULAE.
     ! CHANGED CRITERIION FROM X < 0.01 TO 0.001 BECAUSE SILICATE
     ! SCATTERING WAS NOT CORRECT.
     !
     C=(E-1.)/(E+2.)
     QSC=(X**4/.375)*DBLE(C**2)
     A=DIMAG(-4.*C)
     B=DIMAG(-C*(E*E+27.*E+38.)/(2.*E+3.)/3.75)
     QAB=X*(A+X*X*B)
     QEX=QAB+QSC
     !
     ! G THE ASYMMETRY PARAMETER IS ALWAYS NEGLIGIBLE FOR SMALL PARTICLES.
     !
     G=0.0
     RETURN
  END IF
  !
  ! FULL MIE THEORY CALCULATION.
  ! RM - COMPLEX REFRACTIVE INDEX
  !
  RM=CDSQRT(E)
  EN1=DBLE(RM)
  EN2=DIMAG(RM)
  Y=X*RM
  ZN2=DCMPLX(DCOS(X),-DSIN(X))
  ZN1=DCMPLX(DSIN(X),DCOS(X))
  RIND=EN1**2+EN2**2     ! Rind = |rm|²
  NTIL=1.5*SQRT(RIND)*X+1
  NTOT=MAX0(20,NTIL)
  !
  if (ntot.le.700000) then    ! go ahead with full Mie theory
     !
     AN(NTOT)=DCMPLX(0,0)
     SUME=0.
     SUMS=0.
     SUMG1=0.
     SUMG2=0.
     PSG1=0.
     PSG2=0.
     NTOTA=NTOT
100  P=DFLOAT(NTOTA)
     AN(NTOTA-1)=P/Y-(1./(P/Y+AN(NTOTA)))
     NTOTA=NTOTA-1
     IF(NTOTA.EQ.1) GOTO 101
     GOTO 100
101  AO1=DSIN(EN1*X)*DCOS(EN1*X)
     EN2P=-EN2
     !      IF(EN2P*X.GE.44.)WRITE(6,*)'EN2P,X,LAM,RAD,E1,E2',EN2P,X,LAM,
     !     >RAD,E1,E2
     if(EN2P*X.GE.350.) then
        AO=dcmplx(0.0,1.0)
     else
        AO2=DSINH(EN2P*X)*DCOSH(EN2P*X)
        AO3=(DSIN(EN1*X))**2+(DSINH(EN2P*X))**2
        AO=DCMPLX(AO1,AO2)
        AO=AO/AO3
     endif
     A1=-1./Y+(1./(1./Y-AO))
     RRAT=A1/AN(1)
     f=2.0/(x*x)
     DO  N=1,NTOT
        AN(N)=AN(N)*RRAT
     ENDDO
     DO N=1,NTOT
        P=DFLOAT(N)
        ZN=DFLOAT(2*N-1)*ZN1/X-ZN2
        C=AN(N)/RM+P/X
        A=C*DBLE(ZN)-DBLE(ZN1)
        A=A/(C*ZN-ZN1)
        C=RM*AN(N)+P/X
        B=C*DBLE(ZN)-DBLE(ZN1)
        B=B/(C*ZN-ZN1)
        !
        ! PP, PPG1, PPG2 ARE CONSTANTS CONTAINING THE N TERMS IN THE 
        ! SUMMATIONS.
        !
        PP=DFLOAT(2*N+1)
        !
        PSS=PP*(A*dCONJG(A)+B*dCONJG(B))
        PSE=PP*DBLE(A+B)
        IF(N.GT.1)THEN
           !
           ! CALCULATE G USING FORMULA ON P.128 OF VAN DE HULST'S BOOK.
           ! HAVE REPLACED N BY (N-1) IN THE FORMULA SO THAT WE CAN USE
           ! PREVIOUS A(N) AND B(N) INSTEAD OF A(N+1) AND B(N+1)
           !
           REN=DFLOAT(N)
           PPG1=(REN-1.)*(REN+1.)/REN
           PPG2=(2.*REN-1.)/((REN-1.)*REN)
           PSG1=PPG1*DBLE(ANM1*dCONJG(A)+BNM1*dCONJG(B))
           PSG2=PPG2*DBLE(ANM1*dCONJG(BNM1))
        END IF
        SUME=SUME+PSE
        SUMS=SUMS+PSS
        SUMG1=SUMG1+PSG1
        SUMG2=SUMG2+PSG2
        D1=ABS(PSE/SUME)
        D2=ABS(PSS/SUMS)
        !        IF(D1.LT.1.E-7.AND.D2.LT.1.E-7) GO TO 5
        PT=ABS(PSS/PP)
        IF(PT.LE.1.E-20) GOTO 5
        !
        ! SAVE PREVIOUS A AND B FOR CALCULATION OF G THE ASYMMETRY PARAMETER
        !
        ANM1=A
        BNM1=B
        ZN2=ZN1
        ZN1=ZN
     ENDDO
5    F=2.0/(X*X)
     QEX=F*SUME
     QSC=F*SUMS
     QAB=F*(SUME-SUMS)
     G=2.0*F*(SUMG1+SUMG2)/QSC
     RETURN
  else
     !
     ! Geometrical optics for big spheres
     !
     call geopt(rm,ans)
     qex =2.0d0
     g=9.23d-01   !approx true for D&L silicate.......
     qsc=ans
  end if
  return
END SUBROUTINE Q_MIE


!-------------------------------------------------------------------
!                  FOR THE MIE ROUTINE
!
! intgrates the reflection coefficient
! trapezium rule integration from 0 to pi/2
!-------------------------------------------------------------------
subroutine geopt(m,ans)
  implicit double precision (a-h,o-z)
  complex*16 m
  a=0.0d0
  b=1.570796327d0
  nstrip = 5000
  tot=0
  h=(b-a)/dfloat(nstrip)   !strip width
  tot=tot+0.5*ref(m,a)*h   !1st term
  do i=1,nstrip-1
     x=a+h*dfloat(i)
     tot=tot+ref(m,x)*h      !middle terms
  end do
  tot=tot+0.5*ref(m,b)*h   !last term
  ans=1.+2.*tot    !ans is Qsca
  return
end subroutine geopt
      
!-------------------------------------------------------------------
!                    FOR THE MIE ROUTINE
!
! Calculates Reflection coeffs
!-------------------------------------------------------------------
function ref(m,thetai)
  implicit double precision (a-h,o-z)
  complex*16 sinTHETAt,cosTHETAt ,m,rpll,rper          
  sinTHETAt=sin(THETAi)/m
  cosTHETAt=cdsqrt(1-(sinTHETAt*sinTHETAt))
  !
  !      r for E parallel to plane
  !
  rpll = (cosTHETAt-m*cos(THETAi)) / (cosTHETAt+m*cos(THETAi))
  !
  !       r for E perp. to plane
  !
  rper = (cos(THETAi)-m*cosTHETAt) / (cos(THETAi)+m*cosTHETAt)
  !
  !       R = ½(|rpll|²+|rper|²)
  !
  R= (abs(rpll)*abs(rpll) + abs(rper)*abs(rper))/2.0
  ref=r*sin(THETAi)*cos(THETAi)
  return                                            
end function ref



!-------------------------------------------------------------------
!                      THE CDE ROUTINE 
!-------------------------------------------------------------------
subroutine q_cde(e1,e2,lam,rad,qex,qsc,qabs,g)
  implicit doubleprecision (a-h,o-z)
  double precision lsol,msol
  double precision e1,e2,lam ,rad
  double precision qex,qsc,qabs,g
  double precision fact
  complex*8 refrel,epsil
  double precision x    
  !
  ! size parameter
  !
  X = 2.0d0*PI*rad/lam
  !
  ! The formula for qabs is valid in the rayleigh limit. The formula is
  ! given in Bohren en Huffman p.355-366. The use here is restricted to
  ! dust in a vacuum. 
  !
  REFREL = cmplx(e1,e2)
  EPSIL  = REFREL*REFREL
  FACT   = 4.0d0*X/3.0d0
  FACT   = FACT*PI*rad*rad    !calculate the absorbtion coefficient
  qabs   = 2.0d0*DIMAG(EPSIL/(EPSIL-1.0d0)*LOG(EPSIL))
  qabs   = FACT*qabs
  qsc    = 0.0d0
  qex    = qabs+qsc
  g      = 1.0d0
  
  return
end subroutine q_cde



end module dust_module
