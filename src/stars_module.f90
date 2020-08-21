module stars_module
use rtglobal_module
use amr_module
use amrray_module
use mathroutines_module
use constants_module

integer :: nstars=0,nillum=0
integer :: stellarsrc_nrtemplates
logical :: star_sphere=.false.
doubleprecision,allocatable :: star_spec(:,:)
doubleprecision,allocatable :: star_spec_cum(:,:)
doubleprecision,allocatable :: star_lum(:)
doubleprecision,allocatable :: star_lumcum(:)
doubleprecision,allocatable :: star_r(:)
doubleprecision,allocatable :: star_m(:)
doubleprecision,allocatable :: star_pos(:,:)
doubleprecision,allocatable :: star_fraclum(:)
integer,allocatable         :: star_cellindex(:)
doubleprecision,allocatable :: stellarsrc_dens(:,:)
doubleprecision,allocatable :: stellarsrc_mass(:)
doubleprecision,allocatable :: stellarsrc_templates(:,:)
doubleprecision,allocatable :: stellarsrc_templateslum(:)
doubleprecision,allocatable :: stellarsrc_cum(:)
doubleprecision,allocatable :: stellarsrc_speccum(:,:)
doubleprecision,allocatable :: stellarsrc_lumcum(:)
doubleprecision,allocatable :: extlum_intens(:)
doubleprecision,allocatable :: extlum_intens_cum(:)
doubleprecision,allocatable :: heatsource(:)
doubleprecision,allocatable :: heatsource_cum(:)
doubleprecision :: starlumtot,stellarsrclumtot,extlumtot,heatsourcelumtot
!!doubleprecision :: extlum_rsphere,extlum_xsphere,extlum_ysphere,extlum_zsphere
doubleprecision :: stellarsrc_masstot,star_maxexterndist
doubleprecision,allocatable :: illum_costheta(:),illum_phi(:)
doubleprecision,allocatable :: illum_flux_unprojected(:,:)
doubleprecision,allocatable :: illum_flux_unprojected_cum(:,:)
doubleprecision,allocatable :: illum_fluxtot(:)
doubleprecision,allocatable :: illum_fluxtot_cum(:)

!$OMP THREADPRIVATE(stellarsrc_cum)
!$OMP THREADPRIVATE(star_fraclum)

contains


!-------------------------------------------------------------------------
!                 MAIN ROUTINE FOR STAR READING
!-------------------------------------------------------------------------
subroutine read_stars()
  use constants_module
  implicit none
  logical :: fex,warn,warn_star_in_grid,warn_star_must_be_sphere,fexill
  integer :: iformat,ierr,inu,istar,inupeak,illum
  double precision :: dummy,frq,lam,rpos
  character*80 :: str1
  !
  ! Check if frequency array is already read in
  !
  if(freq_nr.le.0) then
     write(stdo,*) 'ERROR: Internal error in stars_module'
     stop
  endif
  !
  ! Reset some stuff
  !
  if(allocated(amrray_spheres_r)) deallocate(amrray_spheres_r)
  if(allocated(amrray_spheres_pos)) deallocate(amrray_spheres_pos)
  if(allocated(amrray_spheres_outsidegrid)) deallocate(amrray_spheres_outsidegrid)
  if(allocated(amrray_spheres_sphidx)) deallocate(amrray_spheres_sphidx)
  !
  ! Check existence of the stars.inp file
  !
  inquire(file='stars.inp',exist=fex)
  !
  ! Specially for the cartesian 1-D plane-parallel or 2-D pencil-parallel mode
  ! we have illum.inp
  !
  inquire(file='illum.inp',exist=fexill)
  !
  ! If we have cartesian 1-D plane-parallel or 2-D pencil-parallel mode, then stars
  ! are not allowed, as they have no meaning in that context
  !
  if(fex.and.((igrid_coord.eq.10).or.(igrid_coord.eq.20))) then
     write(stdo,*) 'ERROR in stars module: When using cartesian 1-D plane-parallel or 2-D pencil-parallel mode'
     write(stdo,*) '      you are not allowed to use stars in the model, because that has no physical meaning.'
     write(stdo,*) '      Please remove the stars.inp file.'
     stop
  endif
  !
  ! If we do NOT have cartesian 1-D plane-parallel or 2-D pencil-parallel mode, then 
  ! the illum.inp file is not allowed.
  !
  if(fexill.and.((igrid_coord.ne.10).and.(igrid_coord.ne.20))) then
     write(stdo,*) 'ERROR in stars module: The illum.inp file is only meaningful for cartesian 1-D plane-parallel'
     write(stdo,*) '      or 2-D pencil-parallel mode, which you do not use right now.'
     write(stdo,*) '      Please remove the illum.inp file.'
     stop
  endif
  !
  ! Temporarily switch off the possibility of using illumination fluxes
  ! for the 2-D pencil-parallel mode, because this requires a bit of
  ! subtle implementation to the get total fluxes (luminosities per length)
  ! right. Not sure if I will have the time to implement this.
  !
  if(fexill.and.(igrid_coord.eq.20)) then
     write(stdo,*) 'For the moment I am not sure if the illum.inp is consistent'
     write(stdo,*) '  with the cartesian 2-D pencil-parallel geometry. Check this...'
     stop
  endif
  !
  ! Defaults
  !
  nstars = 0
  nillum = 0
  !
  ! Read
  !
  if(fex) then
     !
     ! Open and read the stars.inp file
     !
     open(unit=1,file='stars.inp')
     read(1,*) iformat
     if((iformat.ne.1).and.(iformat.ne.2)) then
        write(stdo,*) 'ERROR: Format number of stars.inp not implemented.'
        stop
     endif
     read(1,*) nstars,inu
     if(nstars.lt.0) then
        write(stdo,*) 'ERROR: Negative nr of stars unallowed.'
        stop
     endif
     if(nstars.gt.0) then
        !
        ! OK, there are stars, so read them
        !
        ! First another check
        !
        if(inu.ne.freq_nr) then
           write(stdo,*) 'ERROR: Nr of frequency points in stars.inp is ',&
                'unequal to that in frequency.inp'
           stop
        endif
        !
        ! Allocate the arrays
        !
        allocate(star_spec(1:freq_nr,1:nstars),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate space for stellar spectrum'
           stop
        endif
        allocate(star_spec_cum(1:freq_nr+1,1:nstars),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate space for cumulative stellar spectrum'
           stop
        endif
        allocate(star_lum(1:nstars),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate space for stellar luminosity'
           stop
        endif
        allocate(star_lumcum(1:nstars+1),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate space for cumul stellar luminosity'
           stop
        endif
        allocate(star_r(1:nstars),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate space for stellar radius'
           stop
        endif
        allocate(star_m(1:nstars),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate space for stellar mass'
           stop
        endif
        allocate(star_pos(1:3,1:nstars),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate space for stellar position'
           stop
        endif
        allocate(star_cellindex(1:nstars),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate space for stellar cell index'
           stop
        endif
        !$OMP PARALLEL
        allocate(star_fraclum(1:nstars),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate space for fraclum array'
           stop
        endif
        star_fraclum(:) = 1.d0
        !$OMP END PARALLEL
        !
        ! Read the stellar information
        !
        do istar=1,nstars
           read(1,*) star_r(istar),star_m(istar),star_pos(1:3,istar)
        enddo
        !
        ! If we have symmetries in our grid, then we should check that the star
        ! position(s) are not located wrongly.
        !
        if(igrid_type.lt.100) then
           !
           ! AMR-type grid (or regular grid)
           !
           if((amr_coordsystem.ge.100).and.(amr_coordsystem.lt.200)) then
              !
              ! Spherical coordinates
              !
              if((amr_zdim.eq.0).and.(amr_ydim.eq.1)) then
                 do istar=1,nstars
                    if((star_pos(1,istar).ne.0.d0).or.(star_pos(2,istar).ne.0.d0)) then
                       write(stdo,*) 'ERROR in star position: In axial symmetry a star is only'
                       write(stdo,*) '      allowed to be located on the z-axis.'
                       stop
                    endif
                 enddo
              endif
              if((amr_zdim.eq.0).and.(amr_ydim.eq.0)) then
                 if(nstars.gt.1) then
                    write(stdo,*) 'ERROR in star position: In spherical symmetry at most 1 star'
                    write(stdo,*) '      is allowed.'
                    stop
                 endif
                 if((star_pos(1,1).ne.0.d0).or.(star_pos(2,1).ne.0.d0).or. &
                      (star_pos(3,1).ne.0.d0)) then
                    write(stdo,*) 'ERROR in star position: In spherical symmetry a star is only'
                    write(stdo,*) '      allowed to be located at (0,0,0).'
                    stop
                 endif
              endif
              warn_star_in_grid = .false.
              warn_star_must_be_sphere = .false.
              do istar=1,nstars
                 rpos = sqrt(star_pos(1,istar)**2+star_pos(2,istar)**2+star_pos(3,istar)**2)
                 if((rpos.ge.amr_grid_xi(1,1)).and. &
                    (rpos.le.amr_grid_xi(amr_grid_nx+1,1))) then
                    warn_star_in_grid = .true.
                 endif
                 if(rpos.lt.amr_grid_xi(1,1)) then
                    if(star_r(istar).ge.0.1d0*amr_grid_xi(1,1)) then
                       warn_star_must_be_sphere = .true.
                    endif
                 endif
              enddo
              if(warn_star_in_grid) then
                 write(stdo,*) 'Notice: Found a star inside the spherical coordinate grid.'
                 write(stdo,*) '        That is perfectly fine, just be aware of this.'
              endif
              if(.not.star_sphere) then
                 if(warn_star_must_be_sphere) then
                    write(stdo,*) '*****************************************************'
                    write(stdo,*) 'WARNING: Central star size is more than 10% of'
                    write(stdo,*) '         inner grid radius. But you still treat'
                    write(stdo,*) '         stars in point-source mode.'
                    write(stdo,*) '         Recommendation: In radmc3d.inp'
                    write(stdo,*) '         set istar_sphere = 1.'
                    write(stdo,*) '*****************************************************'
                 else
                    write(stdo,*) 'Note: Please be aware that you treat the star(s) as'
                    write(stdo,*) '      point source(s) while using spherical coordinate mode.'
                    write(stdo,*) '      Since R_*<<R_in this is probably OK, but if you want'
                    write(stdo,*) '      to play safe, then set istar_sphere = 1 in radmc3d.inp.'
                 endif
              endif
           endif
        endif
        !
        ! Now read the radiation frequencies and check if they are the same
        ! as those of the frequency.inp file.
        !
        if(iformat.eq.1) then 
           do inu=1,freq_nr
              read(1,*) frq
              if(abs(frq-freq_nu(inu))/(frq+freq_nu(inu)).gt.1.d-3) then
                 write(stdo,*) 'PROBLEM: Frequency grid of stellar ',             &
                      'spectrum unequal to frequency.inp or wavelength_micron.inp'
                 write(stdo,*) frq,freq_nu(inu)
                 stop
              endif
           enddo
        else
           do inu=1,freq_nr
              read(1,*) lam
              frq = 1d4*cc/lam
              if(abs(frq-freq_nu(inu))/(frq+freq_nu(inu)).gt.1.d-3) then
                 write(stdo,*) 'PROBLEM: Frequency grid of stellar spectrum ',             &
                      'unequal to frequency.inp or wavelength_micron.inp'
                 write(stdo,*) frq,freq_nu(inu)
                 stop
              endif
           enddo
        endif
        !
        ! Now read all the spectra
        !
        ! IMPORTANT: If the first number is negative, then instead of reading
        !            the entire spectrum, a blackbody is taken with a temperature
        !            equal to the absolute value of that first number.
        !
        do istar=1,nstars
           !
           ! First read the first number. If positive, then read on, else make
           ! planck function.
           !
           read(1,*) dummy
           if(dummy.lt.0.d0) then
              !
              ! Make a simple blackbody star
              !
              if(nstars.lt.10) then
                 write(stdo,202) ' Note: Star ',istar,&
                      ' is taken to be a blackbody star'
202              format(a12,I1,a32)
                 write(str1,275) abs(dummy)
275              format(F9.0)
                 write(stdo,*) '      at a temperature T = ',trim(adjustl(str1)),' Kelvin'
              endif
              do inu=1,freq_nr
                 star_spec(inu,istar) = bplanck(abs(dummy),freq_nu(inu))
              enddo
           else
              !
              ! Read the full intensity spectrum. Note that the spectrum
              ! is given in the usual flux-at-one-parsec format, so we have
              ! to translate this here to intensity-at-stellar-surface. 
              !
              star_spec(1,istar) = 3.0308410d36 * dummy / star_r(istar)**2
              do inu=2,freq_nr
                 read(1,*) dummy
                 ! NOTE: pc^2/pi = 3.0308410d36
                 star_spec(inu,istar) = 3.0308410d36 * dummy / star_r(istar)**2
              enddo
           endif
        enddo
     endif ! if(nstars.gt.0)
     !
     ! Close file
     !
     close(1)
     !
     ! Now compute the cumulative spectrum
     !
     do istar=1,nstars
        star_spec_cum(1,istar) = 0.d0
        do inu=1,freq_nr
           star_spec_cum(inu+1,istar) = star_spec_cum(inu,istar) + &
                     star_spec(inu,istar) * freq_dnu(inu)
        enddo
        if(star_spec_cum(freq_nr+1,istar).le.0.d0) then
           write(stdo,*) 'ERROR in stars module: Zero or negative stellar luminosity.'
           stop
        endif
        !!!star_lum(istar) = pi * star_spec_cum(freq_nr+1,istar) *    &
        !!!                  fourpi * star_r(istar)**2
        star_lum(istar) = -1d99 ! Should not matter - if it does, there is a bug
        do inu=1,freq_nr+1
           star_spec_cum(inu,istar) = star_spec_cum(inu,istar) /   &
                                      star_spec_cum(freq_nr+1,istar)
        enddo
     enddo
     !
     ! Check if the wavelength range is OK.
     !
     do istar=1,nstars
        !
        ! Determine location of the peak of the spectrum
        !
        dummy = 0.d0
        inupeak = 1
        do inu=1,freq_nr
           if(star_spec(inu,istar).gt.dummy) then
              inupeak=inu
              dummy = star_spec(inu,istar)
           endif
        enddo
        !
        ! Now check if this is far enough from the smallest wavelength
        !
        warn = .false.
        if(freq_nu(1).gt.freq_nu(freq_nr)) then
           if(freq_nu(1).lt.3.d0*freq_nu(inupeak)) warn=.true.
        else
           if(freq_nu(freq_nr).lt.3.d0*freq_nu(inupeak)) warn=.true.
        endif
        if(warn) then
           write(stdo,*) '**************************************************************'
           write(stdo,*) 'WARNING: It looks as though the stellar spectrum'
           write(stdo,*) '         does not have its peak well within the wavelength'
           write(stdo,*) '         domain! This could result in wrong dust temperatures.'
           write(stdo,*) '    Tip: extend the wavelength domain to shorter values.'
           write(stdo,*) '**************************************************************'
        endif
     enddo
     !
     ! If the stars should be treated as spheres, then initialize the
     ! spheres
     !
     if(star_sphere) then
        if(nstars.gt.0) then 
           !
           ! Switch on the star sphere mode
           !
           write(stdo,*) 'Treating stars as finite-size spheres...'
           call flush(stdo)
           allocate(amrray_spheres_r(1:nstars),STAT=ierr)
           if(ierr.ne.0) stop 7349
           allocate(amrray_spheres_pos(1:3,1:nstars),STAT=ierr)
           if(ierr.ne.0) stop 7349
           amrray_spheres_nr = nstars
           do istar=1,nstars
              !
              ! Copy information to Amrray module
              !
              amrray_spheres_r(istar)     = star_r(istar)
              amrray_spheres_pos(1,istar) = star_pos(1,istar)
              amrray_spheres_pos(2,istar) = star_pos(2,istar)
              amrray_spheres_pos(3,istar) = star_pos(3,istar)
              !
              ! Check if the stellar radius is not rediculously small,
              ! which would point to a user error
              !
              if(star_r(istar).lt.1d5) then
                 write(stdo,*) 'ERROR: Stellar radius appears to be smaller than 1 km.'
                 write(stdo,*) '       This must be a user-error. Aborting.'
                 stop
              endif
           enddo
           call amrray_install_spheres()
        endif
     else
        amrray_spheres_nr = 0
     endif
  endif
  !
  ! If there is an illum.inp file, read that
  !
  if(fexill) then
     !
     ! Open and read the illum.inp file
     !
     open(unit=1,file='illum.inp')
     read(1,*) iformat
     if((iformat.ne.1).and.(iformat.ne.2)) then
        write(stdo,*) 'ERROR: Format number of illum.inp not implemented.'
        stop
     endif
     read(1,*) nillum,inu
     if(nillum.lt.0) then
        write(stdo,*) 'ERROR: Negative nr of illumination fluxes not allowed.'
        stop
     endif
     if(nillum.gt.0) then
        !
        ! OK, there are illumination fluxes, so read them
        !
        ! First another check
        !
        if(inu.ne.freq_nr) then
           write(stdo,*) 'ERROR: Nr of frequency points in illum.inp is ',&
                'unequal to that in frequency.inp'
           stop
        endif
        !
        ! Allocate the arrays
        !
        allocate(illum_flux_unprojected(1:freq_nr,1:nillum),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate space for illumination flux spectra'
           stop
        endif
        allocate(illum_flux_unprojected_cum(1:freq_nr+1,1:nillum),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate space for cumulative illumination flux spectra'
           stop
        endif
        allocate(illum_fluxtot(1:nillum),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate space for total illumination fluxes'
           stop
        endif
        allocate(illum_fluxtot_cum(1:nillum+1),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate space for total cumulative illumination fluxes'
           stop
        endif
        allocate(illum_costheta(1:nillum),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate space for angles of incidence of illumination fluxes'
           stop
        endif
        allocate(illum_phi(1:nillum),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate space for angles of incidence of illumination fluxes'
           stop
        endif
        !
        ! Read the angles of incidence of the illumination fluxes
        ! NOTE: For 1-D plane parallel the second angle is irrelevant; can keep it 0.
        !
        do illum=1,nillum
           read(1,*) illum_costheta(illum),illum_phi(illum)
           illum_costheta(illum) = cos(illum_costheta(illum)*pi/180.d0)
           illum_phi(illum)      = illum_phi(illum)*pi/180.d0
           if(illum_costheta(illum).eq.0.d0) then
              write(stdo,*) 'ERROR: Cannot put an illumination flux at incl = 90 degrees.'
              stop
           endif
        enddo
        !
        ! Now read the radiation frequencies and check if they are the same
        ! as those of the frequency.inp file.
        !
        if(iformat.eq.1) then 
           do inu=1,freq_nr
              read(1,*) frq
              if(abs(frq-freq_nu(inu))/(frq+freq_nu(inu)).gt.1.d-3) then
                 write(stdo,*) 'PROBLEM: Frequency grid of illumination ',             &
                      'spectrum unequal to frequency.inp or wavelength_micron.inp'
                 write(stdo,*) frq,freq_nu(inu)
                 stop
              endif
           enddo
        else
           do inu=1,freq_nr
              read(1,*) lam
              frq = 1d4*cc/lam
              if(abs(frq-freq_nu(inu))/(frq+freq_nu(inu)).gt.1.d-3) then
                 write(stdo,*) 'PROBLEM: Frequency grid of illumination spectrum ',    &
                      'unequal to frequency.inp or wavelength_micron.inp'
                 write(stdo,*) frq,freq_nu(inu)
                 stop
              endif
           enddo
        endif
        !
        ! Now read all the spectra 
        !
        do illum=1,nillum
           do inu=1,freq_nr
              read(1,*) dummy
              illum_flux_unprojected(inu,illum) = dummy
           enddo
        enddo
     endif
     !
     ! Close file
     !
     close(1)
     !
     ! Now compute the cumulative spectrum
     !
     do illum=1,nillum
        illum_flux_unprojected_cum(1,illum) = 0.d0
        do inu=1,freq_nr
           illum_flux_unprojected_cum(inu+1,illum) = illum_flux_unprojected_cum(inu,illum) + &
                     illum_flux_unprojected(inu,illum) * freq_dnu(inu)
        enddo
        if(illum_flux_unprojected_cum(freq_nr+1,illum).le.0.d0) then
           write(stdo,*) 'ERROR in stars module: Zero or negative illumination flux.'
           stop
        endif
        do inu=1,freq_nr+1
           illum_flux_unprojected_cum(inu,illum) = illum_flux_unprojected_cum(inu,illum) /   &
                                     illum_flux_unprojected_cum(freq_nr+1,illum)
       enddo
     enddo
     !
     ! Check if the wavelength range is OK.
     !
     do illum=1,nillum
        !
        ! Determine location of the peak of the spectrum
        !
        dummy = 0.d0
        inupeak = 1
        do inu=1,freq_nr
           if(illum_flux_unprojected(inu,illum).gt.dummy) then
              inupeak=inu
              dummy = illum_flux_unprojected(inu,illum)
           endif
        enddo
        !
        ! Now check if this is far enough from the smallest wavelength
        !
        warn = .false.
        if(freq_nu(1).gt.freq_nu(freq_nr)) then
           if(freq_nu(1).lt.3.d0*freq_nu(inupeak)) warn=.true.
        else
           if(freq_nu(freq_nr).lt.3.d0*freq_nu(inupeak)) warn=.true.
        endif
        if(warn) then
           write(stdo,*) '**************************************************************'
           write(stdo,*) 'WARNING: It looks as though the illumination spectrum'
           write(stdo,*) '         does not have its peak well within the wavelength'
           write(stdo,*) '         domain! This could result in wrong dust temperatures.'
           write(stdo,*) '    Tip: extend the wavelength domain to shorter values.'
           write(stdo,*) '**************************************************************'
        endif
     enddo
  endif
  !
  ! If there are no stars, then definitely switch off the star sphere mode
  !
  if(nstars.eq.0) then
     star_sphere = .false.
     amrray_spheres_nr = 0
  endif
  !
end subroutine read_stars


!-------------------------------------------------------------------------
!                   FIND THE STARS THEIR CELLS 
!-------------------------------------------------------------------------
subroutine stars_findcell()
  implicit none
  integer :: istar
  double precision :: r,theta,phi
  type(amr_branch), pointer :: a
  integer :: ix,iy,iz
  !
  if(nstars.ge.1) then
     if(igrid_type.lt.100) then
        !
        ! AMR-style grid
        !
        if(amr_coordsystem.lt.100) then
           !
           ! Cartesian coordinates
           !
           do istar=1,nstars
              if(amr_tree_present) then
                 call amr_findcell(star_pos(1,istar),star_pos(2,istar),star_pos(3,istar),a)
                 if(associated(a)) then
                    star_cellindex(istar) = a%leafindex
                 else
                    star_cellindex(istar) = 0
                 endif
              else
                 call amr_findbasecell(star_pos(1,istar),star_pos(2,istar),star_pos(3,istar),ix,iy,iz)
                 if(ix.gt.0) then
                    star_cellindex(istar) = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
                 else
                    star_cellindex(istar) = 0
                 endif
              endif
           enddo
        elseif((amr_coordsystem.ge.100).and.(amr_coordsystem.lt.200)) then
           !
           ! Spherical coordinates
           !
           do istar=1,nstars
              !
              ! Do some sanity checks
              !
              ! If largest theta == pihalf, then we have, by rule,
              ! equatorial mirroring.
              !
              if((amr_coordsystem.ge.100).and.(amr_coordsystem.lt.200).and. &
                   (amr_grid_xi(amr_grid_ny+1,2).eq.pihalf).and.            &
                   (amr_ydim.eq.1)) then
                 !
                 ! Yes, mirror symmetry in the equatorial plane. This means
                 ! that the non-central stars MUST be EXACTLY on the 
                 ! equatorial plane. 
                 !
                 r = sqrt(star_pos(1,istar)**2+star_pos(2,istar)**2)
                 if(abs(star_pos(3,istar)).gt.1d-8*r) then
                    write(stdo,*) 'ERROR: In spherical coordinates with mirror symmetry'
                    write(stdo,*) '       in the equatorial plane, all stars must have '
                    write(stdo,*) '       EXACTLY z=0 (i.e. they must reside on the '
                    write(stdo,*) '       equatorial plane).'
                    stop
                 endif
                 if(abs(star_pos(3,istar)).gt.0.d0) then
                    star_pos(3,istar) = 0.d0
                 endif
              else
                 !
                 ! No, no mirror symmetry in the equatorial plane. This
                 ! means that non-central stars should NOT be exactly on
                 ! the equatorial plane.
                 !
                 r = sqrt(star_pos(1,istar)**2+star_pos(2,istar)**2)
                 if(abs(star_pos(3,istar)).lt.1d-8*r) then
                    !
                    ! Star exactly on equatorial plane. This is asking for 
                    ! trouble! So raise it a tiny bit.
                    !
                    star_pos(3,istar) = 1d-8 * r
                 endif
              endif
              !
              ! Check that in 3-D we don't put non-(0,0,0) stars exactly
              ! at phi=0 or phi=twopi.
              !
              if((amr_zdim.eq.1).and.(star_pos(1,istar).gt.0.d0)) then
                 if(abs(star_pos(2,istar)).lt.1d-8*star_pos(1,istar)) then
                    star_pos(2,istar) = 1d-8*star_pos(1,istar)
                 endif
              endif
              !
              ! Convert x,y,z to spherical coordinates
              !
              call amr_xyz_to_rthphi(star_pos(1,istar),star_pos(2,istar),&
                                     star_pos(3,istar),r,theta,phi)
              !
              ! Find the cell in which the star resides
              !
              if(theta.eq.pihalf) then
                 theta = pihalf - 1d-8   ! Only for cell searching
              endif
              if(amr_tree_present) then
                 call amr_findcell(r,theta,phi,a)
                 if(associated(a)) then
                    star_cellindex(istar) = a%leafindex
                 else
                    star_cellindex(istar) = 0
                 endif
              else
                 call amr_findbasecell(r,theta,phi,ix,iy,iz)
                 if(ix.gt.0) then
                    star_cellindex(istar) = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
                 else
                    star_cellindex(istar) = 0
                 endif
              endif
           enddo
        else
           write(stdo,*) 'ERROR: Other coordinate not yet implemented'
           stop 6491
        endif
     else
        write(stdo,*) 'ERROR: Other grid types not yet implemented'
        stop 2230
     endif
  endif
end subroutine stars_findcell


!-------------------------------------------------------------------------
!                 READ CONTINUUUM STELLAR SOURCE
!-------------------------------------------------------------------------
subroutine read_stellarsource(action)
  use constants_module
  implicit none
  integer :: action
  logical :: fex1,fex2,fex3
  integer :: nrfr,itempl,ierr,inu,i,index,irec,style
  integer :: ileaf,leafindex,idum,precis,reclenn
  double precision :: frq,lam,mstr,rstr,factor,dummy
  character*80 :: str1
  integer(kind=8) :: iiformat,reclen,nn,ik
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(stellarsrc_dens)) return
  endif
  !
  ! Check
  !
  if(.not.(allocated(freq_nu))) then
     write(stdo,*) 'ERROR: freq_nu not allocated when reading stellar source files'
     stop
  endif
  !
  ! Check if we can use the cellindex list
  !
  if(.not.allocated(cellindex)) then
     write(stdo,*) 'ERROR while computing stellar sources mass: The cellindex list is not yet made.'
     stop
  endif
  !
  ! Check if we have the cell volume allocated
  !
  if(.not.allocated(cellvolume)) then
     write(stdo,*) 'ERROR while computing stellar sources mass: The cellvolume list is not yet made.'
     stop
  endif
  !
  ! Reset
  !
  if(allocated(stellarsrc_dens)) deallocate(stellarsrc_dens)
  if(allocated(stellarsrc_mass)) deallocate(stellarsrc_mass)
  if(allocated(stellarsrc_templates)) deallocate(stellarsrc_templates)
  if(allocated(stellarsrc_templateslum)) deallocate(stellarsrc_templateslum)
  !$OMP PARALLEL
  if(allocated(stellarsrc_cum)) deallocate(stellarsrc_cum)
  !$OMP END PARALLEL
  if(allocated(stellarsrc_speccum)) deallocate(stellarsrc_speccum)
  if(allocated(stellarsrc_lumcum)) deallocate(stellarsrc_lumcum)
  !
  ! Read the stellar templates. These are the discrete stellar types
  ! that you can distribute independently. For instance, for a simple
  ! model of a galaxy you want to have old red stars and blue young stars.
  ! So you would need at least 2 templates.
  !
  inquire(file='stellarsrc_templates.inp',exist=fex1)
  if(fex1) then
     !
     ! Found stellar templates, so now we assume that stellar source
     ! are used. Check this.
     !
     inquire(file='stellarsrc_density.inp',exist=fex1)
     inquire(file='stellarsrc_density.uinp',exist=fex2)
     inquire(file='stellarsrc_density.binp',exist=fex3)
     idum=0
     if(fex1) idum=idum+1
     if(fex2) idum=idum+1
     if(fex3) idum=idum+1
     if(idum.gt.1) then
        write(stdo,*) 'ERROR: Found more than one file stellarsrc_density.*inp'
        stop
     endif
     if(idum.eq.0) then
        write(stdo,*) 'ERROR: Found stellarsrc_templates.inp, but not'
        write(stdo,*) '       stellarsrc_density.inp or stellarsrc_density.*inp'
        stop
     endif
     !
     ! Message
     !
     write(stdo,*) 'Reading continuously distributed stellar sources...'
     call flush(stdo)
     !
     ! Now read the templates
     !
     open(unit=1,file='stellarsrc_templates.inp')
     read(1,*) iiformat
     if((iiformat.lt.1).or.(iiformat.gt.2)) then
        write(stdo,*) 'ERROR: in stellarsrc_templates.inp: iformat not known'
        stop
     endif
     read(1,*) stellarsrc_nrtemplates
     if(stellarsrc_nrtemplates.le.0) then
        write(stdo,*) 'ERROR: nr of templates 0 or negative'
        stop
     endif
     read(1,*) nrfr
     if(nrfr.ne.freq_nr) then
        write(stdo,*) 'ERROR while reading stellarsrc_templates.inp:'
        write(stdo,*) '      The number of wavelengths unequal to that of '
        write(stdo,*) '      the global array of wavelengths_micron.inp'
        stop
     endif
     !
     ! Allocate the arrays
     !
     allocate(stellarsrc_dens(1:stellarsrc_nrtemplates,1:nrcells),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate space for stellarsrc density'
        stop
     endif
     allocate(stellarsrc_mass(1:stellarsrc_nrtemplates),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate space for stellarsrc masses'
        stop
     endif
     allocate(stellarsrc_templates(1:freq_nr,1:stellarsrc_nrtemplates),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate space for stellarsrc templates'
        stop
     endif
     allocate(stellarsrc_templateslum(1:stellarsrc_nrtemplates),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate space for stellarsrc templateslum'
        stop
     endif
     !$OMP PARALLEL
     allocate(stellarsrc_cum(1:stellarsrc_nrtemplates+1),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate space for stellarsrc cum'
        stop
     endif
     !$OMP END PARALLEL
     allocate(stellarsrc_speccum(1:freq_nr+1,1:stellarsrc_nrtemplates),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate space for stellarsrc speccum'
        stop
     endif
     allocate(stellarsrc_lumcum(1:nrcells+1),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate space for stellarsrc lumcum'
        stop
     endif
     !
     ! Now read the radiation frequencies and check if they are the same
     ! as those of the frequency.inp or wavelenth_micron.inp file.
     !
     if(iiformat.eq.1) then 
        do inu=1,freq_nr
           read(1,*) frq
           if(abs(frq-freq_nu(inu))/(frq+freq_nu(inu)).gt.1.d-3) then
              write(stdo,*) 'PROBLEM: Frequency grid of stellar source ',             &
                   'spectrum unequal to frequency.inp or wavelength_micron.inp'
              write(stdo,*) frq,freq_nu(inu)
              stop
           endif
        enddo
     else
        do inu=1,freq_nr
           read(1,*) lam
           frq = 1d4*cc/lam
           if(abs(frq-freq_nu(inu))/(frq+freq_nu(inu)).gt.1.d-3) then
              write(stdo,*) 'PROBLEM: Frequency grid of stellar source spectrum ',             &
                   'unequal to frequency.inp or wavelength_micron.inp'
              write(stdo,*) frq,freq_nu(inu)
              stop
           endif
        enddo
     endif
     !
     ! Now read all the spectra
     !
     ! IMPORTANT: If the first number is negative, then instead of reading
     !            the entire spectrum, a blackbody is taken with a temperature
     !            equal to the absolute value of that first number.
     !
     do itempl=1,stellarsrc_nrtemplates
        !
        ! First read the first number. If positive, then read on, else make
        ! planck function.
        !
        read(1,*) dummy
        if(dummy.lt.0.d0) then
           !
           ! Make a simple blackbody star.
           !
           if(stellarsrc_nrtemplates.lt.10) then
              write(stdo,202) ' Note: Stellar template ',itempl,      &
                   ' is taken to be a blackbody star'
202           format(a24,I1,a32)
              write(str1,275) abs(dummy)
275           format(F9.0)
              write(stdo,*) '      at a temperature T = ',trim(adjustl(str1)),' Kelvin'
           endif
           !
           ! But we need more information: the stellar mass and its radius
           !
           read(1,*) rstr           
           read(1,*) mstr
           !
           ! The template is: stellar surface * pi*B_nu(T_star) / mass of star / ster
           ! The factor converts from the planck function to the luminosity
           ! per gram of star divided by 4pi (because the j_nu is always per ster).
           !
           factor = 4*pi*rstr**2 * pi / ( 4*pi*mstr )
           !
           ! Now make the template
           !
           do inu=1,freq_nr
              stellarsrc_templates(inu,itempl) = factor *             &
                     bplanck(abs(dummy),freq_nu(inu))
           enddo
        else
           !
           ! Read the full spectrum. Here it is given in units of erg / sec
           ! / Hz / gram-of-star . So multiply this by the density of stars
           ! in units of gram-of-star / cm^3, and divide by 4*pi to get the
           ! stellar source function in units of erg / src / Hz / cm^3 /
           ! steradian.
           !
           stellarsrc_templates(1,itempl) = dummy / fourpi
           do inu=2,freq_nr
              read(1,*) dummy
              stellarsrc_templates(inu,itempl) = dummy / fourpi
           enddo
        endif
        !
        ! Calculate the stellarsrc_templateslum(). Here we must multiply by 4 pi to get 
        ! rid of the "per steradian". 
        !
        dummy = 0.d0
        do inu=1,freq_nr
           dummy = dummy + fourpi * freq_dnu(inu) * stellarsrc_templates(inu,itempl)
        enddo
        stellarsrc_templateslum(itempl) = dummy
        !
        ! Calculate the stellarsrc_speccum
        !
        dummy = 0.d0
        do inu=1,freq_nr
           stellarsrc_speccum(inu,itempl) = dummy
           dummy = dummy + fourpi * freq_dnu(inu) * stellarsrc_templates(inu,itempl) / &
                           stellarsrc_templateslum(itempl) 
        enddo
        stellarsrc_speccum(freq_nr+1,itempl) = 1.d0
     enddo
     !
     ! Close file
     !
     close(1)
     !
     ! Now read the stellar source spatial distributions
     !
     write(stdo,*) 'Reading smooth stellar source distribution...'
     call flush(stdo)
     if(fex1) then
        !
        ! Formatted input
        !
        style = 1
        open(unit=1,file='stellarsrc_density.inp',status='old')
        read(1,*) iiformat
        if(iiformat.ne.1) then
           write(stdo,*) 'ERROR: Format number of stellarsrc_density.inp is invalid/unknown.'
           stop
        endif
        read(1,*) nn
        read(1,*) ik
     elseif(fex2) then
        !
        ! F77-Unformatted input
        !
        style = 2
        open(unit=1,file='stellarsrc_density.uinp',status='old',form='unformatted')
        read(1) iiformat,reclen
        if(iiformat.ne.1) then
           write(stdo,*) 'ERROR: Format number of stellarsrc_density.uinp is invalid/unknown.'
           write(stdo,*) 'Format number = ',iiformat
           write(stdo,*) 'Record length = ',reclen
           stop
        endif
        reclenn=reclen
        read(1) nn,ik
     elseif(fex3) then
        !
        ! C-compliant binary
        !
        style = 3
        open(unit=1,file='stellarsrc_density.binp',status='old',access='stream')
        read(1) iiformat
        if(iiformat.ne.1) then
           write(stdo,*) 'ERROR: Format number of stellarsrc_density.uinp is invalid/unknown.'
           write(stdo,*) 'Format number = ',iiformat
           stop
        endif
        read(1) nn
        precis = nn
        read(1) nn,ik
     endif
     !
     ! Do some checks
     !
     if(nn.ne.nrcells) then
        write(stdo,*) 'ERROR: stellarsrc_density.inp does not have same number'
        write(stdo,*) '       of cells as the grid.'
        write(stdo,*) nn,nrcells
        stop
     endif
     if(ik.ne.stellarsrc_nrtemplates) then
        write(stdo,*) 'ERROR: In the stellarsrc_density.inp file the number'
        write(stdo,*) '  of stellar templates is different from that '
        write(stdo,*) '  specified in stellarsrc_templates.inp'
        write(stdo,*) ik,stellarsrc_nrtemplates
        stop
     endif
     !
     ! Read the data
     !
     do itempl=1,stellarsrc_nrtemplates
        call read_scalarfield(1,style,precis,nrcells,          &
             stellarsrc_nrtemplates,1,itempl,1,1d-99,reclenn,  &
             scalar1=stellarsrc_dens)
     enddo
     !
     ! Close the file
     !   
     close(1)
     !
     ! Reset the stellar sources masses
     !
     do itempl=1,stellarsrc_nrtemplates
        stellarsrc_mass(itempl) = 0.d0
     enddo
     !
     ! Now do a loop over the leafs and stellar templates
     !
     do ileaf=1,nrcells
        leafindex = cellindex(ileaf)
        do itempl=1,stellarsrc_nrtemplates
           stellarsrc_mass(itempl) = stellarsrc_mass(itempl) +   &
             cellvolume(leafindex)*stellarsrc_dens(itempl,leafindex)
        enddo
     enddo
     !
     ! Now compute total star mass in the smooth stellar source distributions
     ! 
     stellarsrc_masstot = 0.d0
     do itempl=1,stellarsrc_nrtemplates
        stellarsrc_masstot = stellarsrc_masstot + stellarsrc_mass(itempl)
     enddo
     !
     ! Switch on stellar sources
     !
     incl_stellarsrc = 1
     !
  else
     !
     ! Switch off stellar sources
     !
     incl_stellarsrc = 0
     !
  endif
end subroutine read_stellarsource



!-------------------------------------------------------------------------
!                       READ EXTERNAL SOURCE
!-------------------------------------------------------------------------
subroutine read_externalsource(action)
  use constants_module
  implicit none
  integer :: action
  logical :: fex1
  integer :: iformat,nrfr,itempl,ierr,inu
  double precision :: frq,lam,mstr,rstr,factor,dummy
  character*80 :: str1
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(extlum_intens)) return
  endif
  !
  ! Check
  !
  if(.not.(allocated(freq_nu))) then
     write(stdo,*) 'ERROR: freq_nu not allocated when reading external source files'
     stop
  endif
  !
  ! Reset
  !
  if(allocated(extlum_intens)) deallocate(extlum_intens)
  if(allocated(extlum_intens_cum)) deallocate(extlum_intens_cum)
  !
  ! Read the external source
  !
  inquire(file='external_source.inp',exist=fex1)
  if(fex1) then
     !
     ! Message
     !
     write(stdo,*) 'Reading external luminosity source...'
     call flush(stdo)
     !
     ! Now read the templates
     !
     open(unit=1,file='external_source.inp')
     read(1,*) iformat
     if((iformat.lt.1).or.(iformat.gt.2)) then
        write(stdo,*) 'ERROR: in external_source.inp: iformat not known'
        stop
     endif
     read(1,*) nrfr
     if(nrfr.ne.freq_nr) then
        write(stdo,*) 'ERROR while reading stellarsrc_templates.inp:'
        write(stdo,*) '      The number of wavelengths unequal to that of '
        write(stdo,*) '      the global array of wavelengths_micron.inp'
        stop
     endif
     !
     ! Allocate the arrays
     !
     allocate(extlum_intens(1:freq_nr),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate space for extlum_intens'
        stop
     endif
     allocate(extlum_intens_cum(1:freq_nr+1),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate space for extlum_intens_cum'
        stop
     endif
     !
     ! Now read the radiation frequencies and check if they are the same
     ! as those of the frequency.inp or wavelenth_micron.inp file.
     !
     if(iformat.eq.1) then 
        do inu=1,freq_nr
           read(1,*) frq
           if(abs(frq-freq_nu(inu))/(frq+freq_nu(inu)).gt.1.d-3) then
              write(stdo,*) 'PROBLEM: Frequency grid of external source ',          &
                   'spectrum unequal to frequency.inp or wavelength_micron.inp'
              write(stdo,*) frq,freq_nu(inu)
              stop
           endif
        enddo
     else
        do inu=1,freq_nr
           read(1,*) lam
           frq = 1d4*cc/lam
           if(abs(frq-freq_nu(inu))/(frq+freq_nu(inu)).gt.1.d-3) then
              write(stdo,*) 'PROBLEM: Frequency grid of external source spectrum ', &
                   'unequal to frequency.inp or wavelength_micron.inp'
              write(stdo,*) frq,freq_nu(inu)
              stop
           endif
        enddo
     endif
     !
     ! Now read the spectrum of the external source. This is the intensity
     ! I [erg/s/cm^2/Hz/ster]. 
     !
     do inu=1,freq_nr
        read(1,*) dummy
        extlum_intens(inu) = dummy
     enddo
     !
     ! Close file
     !
     close(1)
     !
     ! Now determine the location and radius of the sphere
     !
     if(grid_contsph_r.lt.0.d0) then
        write(stdo,*) 'INTERNAL ERROR: Somehow the containing sphere is not set up yet...'
        stop 3407
     endif
!!     if(igrid_type.lt.100) then
!!        !
!!        ! AMR grid
!!        !
!!        if(amr_coordsystem.lt.100) then
!!           !
!!           ! Cartesian coordinates
!!           !
!!           grid_contsph_x = 0.5d0 * ( amr_grid_xi(amr_grid_nx+1,1) + amr_grid_xi(1,1) )
!!           grid_contsph_y = 0.5d0 * ( amr_grid_xi(amr_grid_ny+1,2) + amr_grid_xi(1,2) )
!!           grid_contsph_z = 0.5d0 * ( amr_grid_xi(amr_grid_nz+1,3) + amr_grid_xi(1,3) )
!!           grid_contsph_r = 0.5001d0 * sqrt( ( amr_grid_xi(amr_grid_nx+1,1) - amr_grid_xi(1,1) )**2 +  &
!!                            ( amr_grid_xi(amr_grid_ny+1,2) - amr_grid_xi(1,2) )**2 +  &
!!                            ( amr_grid_xi(amr_grid_nz+1,3) - amr_grid_xi(1,3) )**2 )
!!        elseif(amr_coordsystem.lt.200) then
!!           !
!!           ! Spherical coordinates
!!           !
!!           grid_contsph_x = 0.d0
!!           grid_contsph_y = 0.d0
!!           grid_contsph_z = 0.d0
!!           grid_contsph_r = 1.001 * amr_grid_xi(amr_grid_nx+1,1)
!!        else
!!           !
!!           ! These coordiantes not yet active
!!           !
!!           stop 7871
!!        endif
!!     else
!!        !
!!        ! Unstructured grids not yet active
!!        !
!!        write(stdo,*) 'INTERNAL ERROR: Unstructured grids not yet implemented'
!!        stop            
!!     endif
     !
     ! Swith on this source
     !
     incl_extlum = 1
     !
  else
     !
     ! Switch off this source
     !
     incl_extlum = 0
     !
  endif
end subroutine read_externalsource



!-------------------------------------------------------------------------
!                    CLEANUP ALL STAR ARRAYS
!-------------------------------------------------------------------------
subroutine stars_cleanup
  implicit none
  if(allocated(star_spec)) deallocate(star_spec)
  if(allocated(star_spec_cum)) deallocate(star_spec_cum)
  if(allocated(star_lum)) deallocate(star_lum)
  if(allocated(star_lumcum)) deallocate(star_lumcum)
  if(allocated(star_r)) deallocate(star_r)
  if(allocated(star_m)) deallocate(star_m)
  if(allocated(star_pos)) deallocate(star_pos)
  if(allocated(star_cellindex)) deallocate(star_cellindex)
  if(allocated(stellarsrc_dens)) deallocate(stellarsrc_dens)
  if(allocated(stellarsrc_mass)) deallocate(stellarsrc_mass)
  if(allocated(stellarsrc_templates)) deallocate(stellarsrc_templates)
  if(allocated(stellarsrc_templateslum)) deallocate(stellarsrc_templateslum)
  if(allocated(stellarsrc_speccum)) deallocate(stellarsrc_speccum)
  if(allocated(stellarsrc_lumcum)) deallocate(stellarsrc_lumcum)
  if(allocated(extlum_intens)) deallocate(extlum_intens)
  if(allocated(extlum_intens_cum)) deallocate(extlum_intens_cum)
  if(allocated(heatsource)) deallocate(heatsource)
  if(allocated(heatsource_cum)) deallocate(heatsource_cum)
  if(allocated(illum_costheta)) deallocate(illum_costheta)
  if(allocated(illum_phi)) deallocate(illum_phi)
  if(allocated(illum_flux_unprojected)) deallocate(illum_flux_unprojected)
  if(allocated(illum_flux_unprojected_cum)) deallocate(illum_flux_unprojected_cum)
  if(allocated(illum_fluxtot)) deallocate(illum_fluxtot)
  if(allocated(illum_fluxtot_cum)) deallocate(illum_fluxtot_cum)
  !$OMP PARALLEL
  if(allocated(star_fraclum)) deallocate(star_fraclum)
  if(allocated(stellarsrc_cum)) deallocate(stellarsrc_cum)
  !$OMP END PARALLEL
  incl_stellarsrc           = 0
  stellarsrc_nrtemplates    = 0
  nstars                    = 0
end subroutine stars_cleanup


!-------------------------------------------------------------------
!              FIND STAR INTENSITY BETWEEN FREQ POINTS
!-------------------------------------------------------------------
function find_starlight_interpol(freq,istar)
  implicit none
  doubleprecision :: find_starlight_interpol,freq,eps
  integer :: inu,istar
  if(freq_nu(freq_nr).gt.freq_nu(1)) then
     if((freq.le.freq_nu(1)).or.(freq.ge.freq_nu(freq_nr))) then
        find_starlight_interpol = 0.d0
        return
     endif
  else
     if((freq.ge.freq_nu(1)).or.(freq.le.freq_nu(freq_nr))) then
        find_starlight_interpol = 0.d0
        return
     endif
  endif
  call hunt(freq_nu,freq_nr,freq,inu)
  if((inu.lt.1).or.(inu.ge.freq_nr)) stop 7201
  eps = (freq-freq_nu(inu)) / (freq_nu(inu+1)-freq_nu(inu))
  if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 7202
  find_starlight_interpol =  (1.d0-eps) * star_spec(inu,istar) +    &
                                    eps * star_spec(inu+1,istar)
  return
end function find_starlight_interpol


!-------------------------------------------------------------------
!       FIND EXTERNAL LUMINOSITY INTENSITY BETWEEN FREQ POINTS
!-------------------------------------------------------------------
function find_extlumintens_interpol(freq)
  implicit none
  doubleprecision :: find_extlumintens_interpol,freq,eps
  integer :: inu
  if(freq_nu(freq_nr).gt.freq_nu(1)) then
     if((freq.le.freq_nu(1)).or.(freq.ge.freq_nu(freq_nr))) then
        find_extlumintens_interpol = 0.d0
        return
     endif
  else
     if((freq.ge.freq_nu(1)).or.(freq.le.freq_nu(freq_nr))) then
        find_extlumintens_interpol = 0.d0
        return
     endif
  endif
  call hunt(freq_nu,freq_nr,freq,inu)
  if((inu.lt.1).or.(inu.ge.freq_nr)) stop 7201
  eps = (freq-freq_nu(inu)) / (freq_nu(inu+1)-freq_nu(inu))
  if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 7203
  find_extlumintens_interpol =  (1.d0-eps) * extlum_intens(inu) +    &
                                       eps * extlum_intens(inu+1)
  return
end function find_extlumintens_interpol


!-------------------------------------------------------------------------
!                    JITTER THE STARS A TINY BIT
!
! Note: Tiny bug: for spherical coordinates the jitter in theta and
!       phi does not go correctly (because star_pos is always in 
!       cartesian). But so far this has not caused problems, and
!       since I don't want (for now) to break the automatic test suite,
!       I'll keep it for now. 2017.03.21
!-------------------------------------------------------------------------
subroutine jitter_stars(fact)
  implicit none
  integer :: istar
  double precision :: szx,szy,szz,fact
  double precision :: smallnumberx,smallnumbery,smallnumberz
  parameter(smallnumberx=0.48957949203816064943d-12)
  parameter(smallnumbery=0.38160649492048957394d-12)
  parameter(smallnumberz=0.64943484920957938160d-12)
  !
  if(nstars.ge.1) then
     szx = amr_grid_xi(amr_grid_nx+1,1) - amr_grid_xi(1,1)
     szy = amr_grid_xi(amr_grid_ny+1,2) - amr_grid_xi(1,2)
     szz = amr_grid_xi(amr_grid_nz+1,3) - amr_grid_xi(1,3)
     do istar=1,nstars
        star_pos(1,istar) = star_pos(1,istar) + fact*szx*smallnumberx
        star_pos(2,istar) = star_pos(2,istar) + fact*szy*smallnumbery
        star_pos(3,istar) = star_pos(3,istar) + fact*szz*smallnumberz
     enddo
  endif
end subroutine jitter_stars

!-------------------------------------------------------------------------
!                    READ INTERNAL HEAT SOURCE
!-------------------------------------------------------------------------
subroutine read_internal_heatsource(action)
  use constants_module
  implicit none
  integer :: action
  logical :: fex1,fex2,fex3
  integer :: nrfr,itempl,ierr,inu,i,index,irec,style
  integer :: ileaf,leafindex,idum,precis,reclenn
  double precision :: mstr,rstr,factor,dummy
  character*80 :: str1
  integer(kind=8) :: iiformat,reclen,nn,ik
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(heatsource)) return
  endif
  !
  ! Check if we can use the cellindex list
  !
  if(.not.allocated(cellindex)) then
     write(stdo,*) 'ERROR while computing stellar sources mass: The cellindex list is not yet made.'
     stop
  endif
  !
  ! Check if we have the cell volume allocated
  !
  if(.not.allocated(cellvolume)) then
     write(stdo,*) 'ERROR while computing stellar sources mass: The cellvolume list is not yet made.'
     stop
  endif
  !
  ! Reset
  !
  if(allocated(heatsource)) deallocate(heatsource)
  !
  ! Now read the heat source spatial distribution
  !
  write(stdo,*) 'Reading the heat source spatial distribution...'
  call flush(stdo)
  inquire(file='heatsource.inp',exist=fex1)
  inquire(file='heatsource.uinp',exist=fex2)
  inquire(file='heatsource.binp',exist=fex3)
  !
  ! If the file exists, then allocate the array
  !
  if(fex1.or.fex2.or.fex3) then
     allocate(heatsource(1:nrcells),STAT=ierr)
     if(ierr.gt.0) then
        write(stdo,*) 'ERROR: Could not allocate heatsource(:)'
        stop
     endif
  endif
  !
  ! Read
  !
  if(fex1) then
     !
     ! Formatted input
     !
     style = 1
     open(unit=1,file='heatsource.inp',status='old')
     read(1,*) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of heatsource.inp is invalid/unknown.'
        stop
     endif
     read(1,*) nn
  elseif(fex2) then
     !
     ! F77-Unformatted input
     !
     style = 2
     open(unit=1,file='heatsource.uinp',status='old',form='unformatted')
     read(1) iiformat,reclen
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of heatsource.uinp is invalid/unknown.'
        write(stdo,*) 'Format number = ',iiformat
        write(stdo,*) 'Record length = ',reclen
        stop
     endif
     reclenn=reclen
     read(1) nn
  elseif(fex3) then
     !
     ! C-compliant binary
     !
     style = 3
     open(unit=1,file='heatsource.binp',status='old',access='stream')
     read(1) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of heatsource.uinp is invalid/unknown.'
        write(stdo,*) 'Format number = ',iiformat
        stop
     endif
     read(1) nn
     precis = nn
     read(1) nn
  else
     !
     ! No heatsource file
     !
     return
  endif
  !
  ! Do some checks
  !
  if(nn.ne.nrcells) then
     write(stdo,*) 'ERROR: heatsource.inp does not have same number'
     write(stdo,*) '       of cells as the grid.'
     write(stdo,*) nn,nrcells
     stop
  endif
  !
  ! Read the data
  !
  call read_scalarfield(1,style,precis,nrcells,     &
             1,1,1,1,1d-99,reclenn,scalar0=heatsource)
  !
  ! Close the file
  !   
  close(1)
  !
  ! Switch on the heat source
  !
  incl_heatsource = 1
  !
end subroutine read_internal_heatsource

end module stars_module
