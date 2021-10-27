module rtglobal_module
  !$ use omp_lib
  use constants_module
  use amr_module
  !
  ! A structure containing all scalar parameters for the Monte Carlo
  ! calculation along the Bjorkman & Wood scheme as well as the
  ! single-frequency Monte Carlo. 
  !
  type mc_params
     integer*8 :: nphot_therm               ! Number of photons emitted in 
  !                                         ! total for the thermal Monte Carlo. 
     integer*8 :: nphot_scat                ! Number of photons emitted in 
  !                                         ! a single-wavelength scattering
  !                                         ! Monte Carlo simulation.
     integer*8 :: nphot_spec                ! Like nphot_scat, but now for
  !                                         ! spectra: typically you will
  !                                         ! have nphot_spec << nphot_scat
  !                                         ! to speed up the spectrum calc,
  !                                         ! and that is OK, because we will
  !                                         ! anyway integrate over the images.
     integer*8 :: nphot_mono                ! Number of photons emitted in 
  !                                         ! a monochromatic 
  !                                         ! Monte Carlo simulation.
     integer :: ifast                       ! ifast.ne.0 speeds up code, but 
  !                                         ! at the cost of less accuracy.
  !                                         ! See enthres below
     integer :: iranfreqmode                ! 1 = Interpolate in temperature
  !                                         !     for random frequency 
  !                                         ! 0 = Simply take nearest temperature
     double precision :: enthres            ! If ifast set, then only update
  !                                         ! cell temperature if energy gain
  !                                         ! since last update exceeds enthres
     integer :: irestart                    ! If set, then the simulations is
  !                                         ! restarted from a safety backup.
     integer :: cntdump                     ! Make safety backup every cntdump
  !                                         ! photons.
     integer :: ntemp                       ! Number of bins of temperature
  !                                         ! for precalculated dust thermal
  !                                         ! balance equation.
     integer :: itempdecoup                 ! =0 All dust species in the same
  !                                         !    cell have the same temperature
  !                                         ! =1 Each dust species in the same
  !                                         !    cell has own temperature
     double precision :: temp0,temp1        ! Lower and upper temperature for 
  !                                         ! the precalculated dust thermal
  !                                         ! balance equation.
     double precision :: extinct_thres      ! For the scattering MC: When do
  !                                         ! we consider the photon to be
  !                                         ! gone? This is the abs tau for that.
     integer :: niter_vstruct               ! >0 = Nr of iterations for the 
  !                                         !    vertical structure.
  !                                         !    (For now: not yet ready)
     double precision :: vserrtol           ! Error tolerance for the vertical
  !                                         !    structure mode.
     integer :: ivstrt                      ! Which dust species represents the
  !                                         !    gas temperature for the 
  !                                         !    vertical structure iteration?
     integer :: nphotdiff                   ! >0 = Use diffusion equation
  !                                         !    where photon statistics is bad.
  !                                         !    (For now: not yet ready)
     double precision :: errtoldiff         ! Error tolerance for the diffusion
  !                                         !    mode.
  !!!   integer :: save_scat                   ! =0 No file scatsource_3d.dat is
  !!!!                                         !    dumped.
  !!!!                                         ! =1 File scatsource_3d.dat is
  !!!!                                         !    written with isotropic scat
  !!!!                                         !    source function. Only works
  !!!!                                         !    with scattering_mode==1. 
  !!!!                                         ! <0 (future) experimental 
  !!!!                                         !    approximate ways to save
  !!!!                                         !    aniso scattering source.
     integer :: debug_write_stats           ! >0 Write Monte Carlo statistics
  !                                         !    to a file.
     integer :: debug_write_path            ! =1 Write the path of the photon
  !   integer :: debug_write_eventcounts     ! =1 Write for each photon nr of events
     integer :: countwrite                  ! Each countwrite, write to standard output
  !   integer :: incl_scatsrc_mctherm        ! Do we make scattering src while doing therm MC?
     logical :: optimized_motion            ! For high optical deth cells
     double precision :: optim_dtau         ! For high optical deth cells
     logical :: mod_random_walk             ! Using MRW?
     integer :: mrw_count_trigger           ! How often say in same cell to trigger MRW
     integer :: mrw_db_ntemp                ! For the Planck opacity table: nr of temperatures
     double precision :: mrw_db_temp0       ! Lower temperature (note: table is logarithmic)
     double precision :: mrw_db_temp1       ! Upper temperature
     double precision :: mrw_enthres        ! The relative energy threshold for recalc MRW opacities
     double precision :: mrw_tauthres       ! Minimal Rosseland optical depth for MRW
     double precision :: mrw_gamma          ! The minimal Rosseland tau distance from wall for MRW
     double precision :: mrw_taustepback    ! The tiny step back from the wall at the end of MRW
     double precision :: mrw_tempthres      ! If T_dust < tempthres, then do not use mrw_enthres
  end type mc_params
  !
  ! The Monte Carlo parameter file
  !
  type(mc_params) :: rt_mcparams
  !
  ! All the flags, switches, constants etc defined here are considered
  ! to be "global variables" for all the radiative transfer modules.
  !
  ! Main settings and switches:
  !
  integer :: mc_safetymode               ! If >0 then sacrifice speed for 
  !                                      ! accuracy and safety. Default=0.
  integer :: igrid_type                  ! The type of gridding used
  !                                      ! <100 = Regular or AMR grid
  !                                      ! 101 = Delaunay grid
  !                                      ! 201 = Voronoi grid
  integer :: igrid_coord                 ! The coordinate system
  !                                      !   0 -  99   = Cartesian 
  !                                      ! 100 - 199   = Spherical
  integer :: igrid_mirror                ! 1 = Midplane mirror symmetry
  !                                      ! (only for spherical coordinates)
  integer :: incl_stellarsrc=0           ! >0 Include smooth stellar source
  integer :: incl_extlum=0
  integer :: incl_heatsource=0           ! Include internal heat source
  integer :: incl_quantum=0              ! Include quantum-heated grains
  logical :: dust_2daniso=.false.        ! If true, then do the special 2-D axisymmetric full scattering stuff
  integer :: dust_2daniso_nphi=360       ! If dust_2daniso==.true., then this is the number of phi-angles
  logical :: do_writemodel               ! Write user-defined model to file
  logical :: do_writegridfile            ! Write the grid
  logical :: do_writepop                 ! Write the internally calculated level populations
  logical :: do_calcpop                  ! Calculate level populations (automatic when making images)
  logical :: do_respondwhenready         ! Write '1' to STDO when ready
  logical :: do_montecarlo_therm         ! Do Bjorkman & Wood to 
  !                                      ! find the dust temperature
  !                                      ! and the scattering src
  logical :: do_montecarlo_mono          ! Do monochromatic Monte Carlo to
  !                                      ! find the mean intensity in the
  !                                      ! model. Useful for e.g. chemistry
  logical :: do_userdef_action           ! Do the userdef action
  logical :: do_vstruct                  ! Do vertical structure
  logical :: do_raytrace_spectrum        ! Make a spectrum with the camera freq array
  logical :: do_raytrace_image           ! Make an image
  logical :: do_raytrace_movie           ! Make a series of images
  logical :: do_raytrace_tausurf         ! Make the tau=tausurface surface
  logical :: do_write_lineinfo           ! Dump some information about a line
  logical :: do_writeimage               ! In child mode you must request 
  !                                      ! radmc3d to write the last image.
  !                                      ! In child mode, making an image
  !                                      ! does not yet write it to stdo.
  logical :: do_writespectrum            ! In child mode you must request 
  !                                      ! radmc3d to write the last spectrum.
  !                                      ! In child mode, making a spectrum
  !                                      ! does not yet write it to stdo.
  logical :: do_writescatsrc             ! If set, then the scattering source is written to a file
  logical :: writeimage_unformatted      ! If set, then use unformatter image output
  logical :: circular_images             ! If set, images are made using circular pixel arrangement (spher coord only)
  logical :: do_writesample              ! Write values sampled at a list of points
  logical :: do_sample_dustdens          ! Sample dust density
  logical :: do_sample_dusttemp          ! Sample dust temperature
  logical :: do_sample_levelpop          ! Sample level populations
  logical :: do_write_subbox_dustdens    ! Write subbox of dust density
  logical :: do_write_subbox_dusttemp    ! Write subbox of dust temperature
  logical :: do_write_subbox_levelpop    ! Write subbox of level populations
  logical :: do_write_vtk_grid           ! Write VTK data for the grid
  logical :: do_write_vtk_dust_density   ! Write VTK data for the grid
  logical :: do_write_vtk_dust_temperature ! Write VTK data for the grid
  logical :: do_write_vtk_gas_density    ! Write VTK data for the grid
  logical :: do_write_vtk_gas_temperature ! Write VTK data for the grid
  logical :: do_write_vtk_chemspec       ! Write VTK data for the grid
  logical :: do_write_vtk_levelpop       ! Write VTK data for the grid
  logical :: do_write_vtk_velocity       ! Write VTK data for the grid
  integer :: vtk_dust_ispec              ! If writing dust VTK data, which dust species to write
  integer :: vtk_lines_ispec             ! If writing dust VTK data, which line species to write
  integer :: vtk_lines_ilevel            ! If writing dust VTK data, which line level to write
  !
  ! Debugging options
  !
  integer :: debug_check_all             ! =1 Do on-the-fly checking 
  !
  ! Switch for IO: formatted or unformatted
  !
  integer :: rto_style                   ! 0=ascii,1=f77-unf,2=binary
  logical :: rto_single=.false.          ! If on, use single-precision output
  integer :: rto_reclen                  ! F77-style record length used for 
  !                                      ! unformatted output
  integer :: rtio_gridinfo               ! If 1, then add additional info to
  !                                      ! the grid file
  !
  ! Internal IO switches
  !
  logical :: grid_was_read_from_file     ! Memory where we got the grid from
  logical :: grid_was_written_to_file    
  !
  ! Seed for random generator
  !
  ! ****************** CHECK IF FOR LARGE PHOTON NUMBERS WE DO NOT
  ! GET A RECURRENCE OF RANDOM NUMBERS ****************************
  !
  integer :: iseed
  integer :: iseed_start
  logical :: do_resetseed=.false.
  !
  ! Number of spatial grid cells
  !
  integer :: nrcells = 0
  integer :: nrcellsmax = 0
  integer :: nrcellsinp = 0
  !
  ! Frequency array
  !
  integer :: freq_nr = 0
  doubleprecision, allocatable :: freq_nu(:),freq_dnu(:)
  !
  ! For anisotropic scattering with tabulated phase function: mu grid
  ! Note that they are grid points (hence mui) instead of grid cells.
  !
  integer :: scat_munr = 0
  !integer :: scat_phinr = 0
  doubleprecision, allocatable :: scat_mui_grid(:),scat_thetai_grid(:)
  !
  ! For aligned grains: mu grid
  ! Note that they are grid points (hence mui) instead of grid cells.
  !
  integer :: align_munr = 0
  doubleprecision, allocatable :: align_mui_grid(:),align_etai_grid(:)
  !
  ! Thermal boundaries (only for cartesian coordinates!)
  !
  integer :: incl_thermbc = 0
  logical :: thermal_bc_active(1:2,1:3) = .false.
  doubleprecision :: thermal_bc_temp(1:2,1:3) = 0.d0
  !
  ! Dust density and temperature arrays
  !
  double precision, allocatable :: dustdens(:,:),dusttemp(:,:)
  !
  ! Mass of the dust
  !
  double precision, allocatable :: dust_massdust(:)
  double precision :: dust_massdusttot=0.d0
!  !
!  ! For isotropic scattering: the scattering source function. Note that the
!  ! nice thing about isotropic scattering is that it allows all the scattering
!  ! information to be stored in a single file that is not too big. Then the
!  ! images and spectra can be made with ray-tracing after the Monte Carlo run,
!  ! and can be recomputed at different inclinations and zoom-factors without
!  ! having to redo the Monte Carlo part. The disadvantage is that isotropic
!  ! scattering is less accurate.
!  !
!  ! NOTE: Because the isoscat array can be very large for 3-D models, 
!  !       I make it single precision here. This saves a factor of 2 in
!  !       memory.
!  !
!  real, allocatable :: isoscatsrc(:,:)
  !
  ! Gas density and temperature
  !
  double precision, allocatable :: gasdens(:),gastemp(:)
  double precision :: gastmax=0.d0
  !
  ! For the gas continuum module, for thermal free-free extinction/emission
  !
  double precision, allocatable :: electron_numdens(:),ion_numdens(:)
  !
  ! For grain alignment in polarization of dust emission (and in future: scattering)
  ! The alignment direction (a unit vector) and alignment efficiency (between 0 and 1)
  !
  double precision, allocatable :: grainalign_dir(:,:)
  double precision, allocatable :: grainalign_eff(:)
  !
  ! The cell volume array
  !
  double precision, allocatable :: cellvolume(:)
  !
  ! Array of indices
  !
  integer, allocatable :: cellindex(:)
  !
  ! Current position along a ray (i.e. photon position in MC) and direction
  !
  double precision :: ray_cart_x,ray_cart_y,ray_cart_z
  double precision :: ray_cart_dirx,ray_cart_diry,ray_cart_dirz
  !
  ! For polarization: current S-vector 
  !
  double precision :: ray_cart_svec(1:3)
  !
  ! (For some applications) The previous position along a ray
  !
  double precision :: ray_prev_x,ray_prev_y,ray_prev_z
  !
  ! Distance to travel
  !
  double precision :: ray_dsend,ray_ds
  !
  ! Index of the cell
  !
  integer :: ray_index,ray_indexnext
  !
  ! Frequency index
  !
  integer :: ray_inu
  !
  ! Number of elements of a 1-D ray, as well as its maximum (estimation)
  !
  integer :: ray_ns    = 0
  integer :: ray_nsmax = 0
  !
  !    Some arrays for the general ray tracing 
  !
  double precision, allocatable :: camera_ray_jnu(:,:)
  double precision, allocatable :: camera_ray_alpnu(:,:)
  double precision, allocatable :: camera_ray_ds(:)
  !
  !   Which processes to include?
  !
  logical :: rt_incl_dust    = .true.
  logical :: rt_incl_lines   = .false.
  logical :: rt_incl_gascont = .false.
  logical :: rt_incl_gascont_freefree = .true.    ! Works only if rt_incl_gascont.eq..true.
  logical :: rt_incl_userdef_srcalp = .false.
  !
  ! The line transfer mode
  !   =1      LTE populations, computed beforehand and stored in a large array
  !   =-1     LTE populations, computed on-the-fly to save memory
  !   =-2     User-defined populations, computed on-the-fly to save memory
  !   =3      LVG populations, computed beforehand and stored in a large array
  !   =-3     LVG populations, computed on-the-fly to save memory
  !   =4      Opt thin populations, computed beforehand and stored in a large array
  !   =-4     Opt thin populations, computed on-the-fly to save memory
  !   =-10    User-defined populations of user-defined molecule; totally free
  !   =50     Simply read populations from files and store in a large array
  !   =20..49 Non-LTE non-local modes, various methods [FUTURE]
  !
  integer :: lines_mode = 1
  integer :: lines_show_pictograms = 0
  logical :: lines_make_linelist = .false.
  logical :: lines_pfunc_outofrange = .false.
  logical :: lines_maser_warning = .false.
  logical :: lines_slowlvg_as_alternative = .false.
  !
  !   Use dust temperature as gas temperature?
  !
  integer :: tgas_eq_tdust = 0
  !
  ! Subbox parameters (not important for models, but useful for analysis
  ! of the 3-D variables)
  !
  integer :: subbox_nx,subbox_ny,subbox_nz
  double precision :: subbox_x0,subbox_x1
  double precision :: subbox_y0,subbox_y1
  double precision :: subbox_z0,subbox_z1
  double precision :: subbox_phi1,subbox_phi2,subbox_theta
  !
  ! Some derived variables of the grid that can be useful
  !
  double precision :: grid_contsph_x=-1d99
  double precision :: grid_contsph_y=-1d99
  double precision :: grid_contsph_z=-1d99
  double precision :: grid_contsph_r=-1d99
  !
  ! Line transfer model arrays
  !
  !
  ! Now the big arrays with the level populations, only for the subset of
  ! active levels.
  !
  double precision, allocatable :: lines_levelpop(:,:,:)
  !
  ! (For LTE, LVG, non-LTE internal computation): the number density of the
  ! various chemical species. This is not necessary if you read the level
  ! populations directly using read_levelpop().
  !
  double precision, allocatable :: gas_chemspec_numberdens(:,:)
  !
  ! (For LVG, non-LTE internal computation): the number density of the
  ! various collisional partners. Currently, only assumes collisional 1 partner
  ! for all molecules.
  !
  double precision, allocatable :: collpartner_numberdens(:,:)
  !
  ! Another big array for the velocity field in cm/s
  ! NOTE: This is not particularly associated to lines, but since it is 
  !       only relevant for the lines we put this array here in the line module
  !
  double precision, allocatable :: gasvelocity(:,:)
  !
  ! Another big array for the microturbulence
  !
  double precision, allocatable :: lines_microturb(:)
  !
  ! Yet another big array for the escape probability algorithm: the
  ! length scale to be used at each location for the EscP algorithm
  ! (which is built into the LVG method).
  !
  double precision, allocatable :: lines_escprob_lengthscale(:)
  !
  ! A set of arrays for the ray tracing 
  !
  double precision, allocatable :: lines_ray_levpop(:,:,:)
  double precision, allocatable :: lines_ray_nrdens(:,:)
  double precision, allocatable :: lines_ray_temp(:)
  double precision, allocatable :: lines_ray_turb(:)
  double precision, allocatable :: lines_ray_doppler(:)
  double precision, allocatable :: lines_ray_lorentz_delta(:)
  !
  ! Variables for maximum doppler shift
  !
  double precision :: lines_maxveloc=0.d0,lines_maxturbc=0.d0
  double precision :: lines_maxtempc=0.d0,lines_maxshift=0.d0
  !
  ! To find out if the grid is fine enough (or if instead the
  ! doppler catching method has to be used) the following 
  ! variable can be computed using the lines_compute_maxrellineshift()
  ! subroutine.
  !
  double precision :: lines_maxrelshift=0.d0
  !
  ! Warning that due to insufficient spatial resolution of the model
  ! there are locations where the Doppler shift between two cells is
  ! larger than the local line width, which means that a line could
  ! accidently be 'lost' in the ray-tracing. 
  !
  !logical :: lines_warn_lineleap
  !
  !!$ setthreads sets the number of threads to be used by subsequent parallel regions &
  !!$ in the parallel version of RADMC-3D
  !$ integer :: setthreads
  !
  !$OMP THREADPRIVATE(iseed)
  !$OMP THREADPRIVATE(ray_cart_x,ray_cart_y,ray_cart_z)
  !$OMP THREADPRIVATE(ray_cart_dirx,ray_cart_diry,ray_cart_dirz)
  !$OMP THREADPRIVATE(ray_cart_svec,ray_prev_x,ray_prev_y,ray_prev_z)
  !$OMP THREADPRIVATE(ray_dsend,ray_ds,ray_index,ray_indexnext)
  !$OMP THREADPRIVATE(ray_inu,ray_ns,ray_nsmax)
  !
  !$OMP THREADPRIVATE(lines_ray_levpop,lines_ray_nrdens,lines_ray_temp)
  !$OMP THREADPRIVATE(lines_ray_turb,lines_ray_doppler,lines_ray_lorentz_delta)
  !
contains


! !-------------------------------------------------------------------
! !                           RTGLOBAL INIT
! !-------------------------------------------------------------------
! subroutine rtglobal_init()
!   use amr_module
!   implicit none
!   if(igrid_type.lt.100) then
!      !
!      ! Normal regular grid or AMR grid
!      !
!      if(amr_nrleafs.lt.1) then
!         write(stdo,*) 'INTERNAL ERROR: Nr of grid cells .lt. 1'
!         stop
!      endif
!      nrcells    = amr_nrleafs
!      nrcellsmax = amr_nrleafs_max
!   else
!      write(stdo,*) 'ERROR: No other grid types than regular are allowed yet.'
!      stop
!   endif
! end subroutine rtglobal_init

!-------------------------------------------------------------------
!                   READ THE FREQUENCY ARRAY 
!
! This routine will allocate the frequency array and read the
! values from the file frequency.inp or wavelength_micron.inp.
! You can choose freely whether to use frequency.inp or the new
! style wavelength_micron.inp. The first lists the frequencies
! in Herz while the second lists them in microns.
!-------------------------------------------------------------------
subroutine read_frequencies(action)
  implicit none
  integer inu,ierr,fnlen,action
  doubleprecision frnu,dfrnu,numax,numin
  logical :: fex_f,fex_lm
  character*100 :: filename  
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(freq_nu)) return
  elseif(action.eq.2) then
     if(allocated(freq_nu)) deallocate(freq_nu)
     if(allocated(freq_dnu)) deallocate(freq_dnu)
  endif
  !
  ! Check which file is present
  !
  inquire(file='frequency.inp',exist=fex_f)
  inquire(file='wavelength_micron.inp',exist=fex_lm)
  if(fex_f) then
     if(fex_lm) then
        write(stdo,*) 'ERROR: Cannot have both frequency.inp and ',&
             'wavelength_micron.inp as input files. Must have one of these only.'
        stop
     endif
  else
     if(.not.fex_lm) then
        write(stdo,*) 'ERROR: Must have EITHER frequency.inp OR ',&
             'wavelength_micron.inp as input files to specify which'
        write(stdo,*) '   wavelength/frequency grid is used for the model.'
        stop
     endif
  endif
  if(fex_f) then
     filename = 'frequency.inp'
     fnlen    = 13
  else
     filename = 'wavelength_micron.inp'
     fnlen    = 21
  endif
  !
  ! Message
  !
  write(stdo,*) 'Reading global frequencies/wavelengths...'
  call flush(stdo)
  !
  ! Read this file
  !
  open(unit=1,file=filename,status='old',err=701)
  read(1,*) freq_nr
  if(freq_nr.lt.1) then
     write(stdo,*) 'ERROR: ',filename(1:fnlen),' claims to have < 1 points.'
     stop
  endif
  !
  ! Allocate the frequency arrays
  !
  allocate(freq_nu(1:freq_nr),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate freq_nu'
     stop 
  endif
  allocate(freq_dnu(1:freq_nr),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate freq_dnu'
     stop 
  endif
  !
  ! Now read the frequencies or wavelengths-in-microns
  !
  do inu=1,freq_nr
     read(1,*,end=702,err=702) freq_nu(inu)
     if(fex_lm) then
        freq_nu(inu) = 1d4*cc / freq_nu(inu) 
     endif
  enddo
  close(1)
  !
  ! Do some elementary checks
  !
  if(freq_nu(freq_nr).gt.freq_nu(1)) then
     do inu=2,freq_nr
        if(freq_nu(inu).le.freq_nu(inu-1)) then
           write(stdo,*) 'ERROR while reading wavelength data from ',filename(1:fnlen)
           write(stdo,*) '      Not monotonic!'
           stop
        endif
     enddo
     numin = freq_nu(1)
     numax = freq_nu(freq_nr)
  else
     do inu=2,freq_nr
        if(freq_nu(inu).ge.freq_nu(inu-1)) then
           write(stdo,*) 'ERROR while reading wavelength data from ',filename(1:fnlen)
           write(stdo,*) '      Not monotonic!'
           stop
        endif
     enddo
     numin = freq_nu(freq_nr)
     numax = freq_nu(1)
  endif
  if(1d4*cc/numax.gt.1.d0) then
     write(stdo,*) '***************************************************************'
     write(stdo,*) 'WARNING: While reading wavelength data from ',filename(1:fnlen)
     write(stdo,*) '         Wavelength range does not extend to < 1 micron.'
     write(stdo,*) '         In virtually all applications this means that the '
     write(stdo,*) '         starlight is not covered... If you use RADMC-3D '
     write(stdo,*) '         only for making images, that is ok. But if you use '
     write(stdo,*) '         it for e.g. computing the dust temperature with '
     write(stdo,*) '         radmc3d mctherm, then this goes wrong.'
     write(stdo,*) '***************************************************************'
     !stop
  endif
  if(1d4*cc/numin.lt.100.d0) then
     write(stdo,*) '***************************************************************'
     write(stdo,*) 'WARNING: While reading wavelength data from ',filename(1:fnlen)
     write(stdo,*) '         Wavelength range does not extend to > 100 micron.'
     write(stdo,*) '         In many applications this means that the dust may'
     write(stdo,*) '         not be able to cool properly, yielding too hot'
     write(stdo,*) '         dust. If you use RADMC-3D only for making images, '
     write(stdo,*) '         that is ok. But if you use it for e.g. computing '
     write(stdo,*) '         the dust temperature with radmc3d mctherm, '
     write(stdo,*) '         then this goes wrong.'
     write(stdo,*) '***************************************************************'
     !stop
  endif
  !
  ! Compute the freq_dnu arrays
  !
  if(freq_nr.eq.1) then
     freq_dnu(1) = 1.d0
  else
     freq_dnu(1)       = 0.5d0 * abs( freq_nu(2) - freq_nu(1) )
     freq_dnu(freq_nr) = 0.5d0 * abs( freq_nu(freq_nr) - freq_nu(freq_nr-1) )
     do inu=2,freq_nr-1
        freq_dnu(inu) = 0.5d0 * abs( freq_nu(inu+1) - freq_nu(inu-1) )
     enddo
  endif
  !
  goto 710
701 continue
  write(stdo,*) 'Could not open file ',filename(1:fnlen)
  stop 13
702 continue
  write(stdo,*) filename(1:fnlen),': either wrong format or other reading error'
  stop 
710 continue
  return
end subroutine read_frequencies


!-------------------------------------------------------------------
!                   READ THE SCATTERING ANGLE ARRAY 
!
! The theta-grid is used for the treatment of anisotropic scattering
! with tabulated phase functions. 
!-------------------------------------------------------------------
subroutine read_scat_angular_grid(action)
  implicit none
  integer :: imu,ierr,action,iformat
  logical :: fex
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(scat_mui_grid)) return
  elseif(action.eq.2) then
     if(allocated(scat_mui_grid)) deallocate(scat_mui_grid)
     if(allocated(scat_thetai_grid)) deallocate(scat_thetai_grid)
  endif
  !
  ! Check which file is present
  !
  inquire(file='scattering_angular_grid.inp',exist=fex)
  if(.not.fex) then
     write(stdo,*) 'ERROR: When using anisotropic scattering with tabulated'
     write(stdo,*) '       angular phase functions you must define a global'
     write(stdo,*) '       angular grid with an input file called:'
     write(stdo,*) '       scattering_angular_grid.inp'
     write(stdo,*) '       (see manual for its contents).'
     stop
  endif
  !
  ! Message
  !
  write(stdo,*) 'Reading global angular grid for anisotropic scattering...'
  call flush(stdo)
  !
  ! Read this file
  !
  open(unit=1,file='scattering_angular_grid.inp',status='old',err=701)
  read(1,*) iformat
  if(iformat.ne.1) then
     write(stdo,*) 'ERROR: File format of scattering_angular_grid.inp must be 1'
     stop
  endif
  read(1,*) scat_munr
  !
  ! Allocate the mu angular array. Note that we do not need a phi angular
  ! array because we will always assume it to be linearly spaced.
  !
  allocate(scat_mui_grid(1:scat_munr),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in RTGLobal Module: Could not allocate scat_mui_grid'
     stop 
  endif
  allocate(scat_thetai_grid(1:scat_munr),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in RTGLobal Module: Could not allocate scat_thetai_grid'
     stop 
  endif
  !
  ! Now read the angular grid from 0 all the way to 180
  ! Bugfix 2017.03.15 : theta must be converted into radian
  !
  do imu=1,scat_munr
     read(1,*,end=702,err=702) scat_thetai_grid(imu)
     if(imu.gt.1) then
        if(scat_thetai_grid(imu).le.scat_thetai_grid(imu-1)) then
           write(stdo,*) 'ERROR in file scattering_angular_grid.inp:'
           write(stdo,*) '   theta grid must be monotonically increasing'
           stop
        endif
     endif
     scat_thetai_grid(imu) = scat_thetai_grid(imu)*pi/180.d0
     scat_mui_grid(imu)    = cos(scat_thetai_grid(imu))
  enddo
  close(1)
  !
  ! Do some elementary checks
  !
  if(scat_thetai_grid(1).ne.0.d0) then
     write(stdo,*) 'ERROR while reading scattering_angular_grid.inp file:'
     write(stdo,*) '      The first theta value MUST be 0'
     stop
  endif
  if(abs(scat_thetai_grid(scat_munr)-pi).gt.1d-6) then
     write(stdo,*) 'ERROR while reading scattering_angular_grid.inp file:'
     write(stdo,*) '      The last theta value MUST be 180'
     stop
  endif
  scat_mui_grid(1)         =  1.d0
  scat_mui_grid(scat_munr) = -1.d0
  !
  goto 710
701 continue
  write(stdo,*) 'Could not open file scattering_angular_grid.inp'
  stop 13
702 continue
  write(stdo,*) 'scattering_angular_grid.inp: either wrong format or other reading error'
  stop 
710 continue
  return
end subroutine read_scat_angular_grid


!--------------------------------------------------------------------------
!                THE BLACKBODY PLANCK FUNCTION B_nu(T)
!
!     This function computes the Blackbody function 
!
!                    2 h nu^3 / c^2
!        B_nu(T)  = ------------------    [ erg / cm^2 s ster Hz ]
!                   exp(h nu / kT) - 1
!
!     ARGUMENTS:
!        nu    [Hz]            = Frequency
!        temp  [K]             = Temperature
!
! Bugfix 22.03.2016 by Jon Ramsey to avoid float overflows
!--------------------------------------------------------------------------
function bplanck(temp,nu)
  implicit none
  doubleprecision temp
  doubleprecision nu
  doubleprecision bplanck,xx
  doubleprecision, parameter :: explargest = 709.78271d0
  if(temp.eq.0.d0) then 
     bplanck = 0.d0
     return
  endif
  xx      = 4.7989d-11 * nu / temp
  ! prevents floating overflows in the exponential
  if(xx .gt. explargest) then
     bplanck = 0.0d0
     return
  endif
  bplanck = 1.47455d-47 * nu * nu * nu /                          &
            (exp(xx)-1.d0) + 1.d-290
  return
end function bplanck


!--------------------------------------------------------------
!      THE TEMPERATURE DERIVATIVE OF PLANCK FUNCTION 
!
! This function computes the temperature derivative of the
! Blackbody function 
! 
!    dB_nu(T)     2 h^2 nu^4      exp(h nu / kT)        1 
!    --------   = ---------- ------------------------  ---
!       dT          k c^2    [ exp(h nu / kT) - 1 ]^2  T^2
!
! ARGUMENTS:
!    nu    [Hz]            = Frequency
!    temp  [K]             = Temperature
!
! Bugfix 22.03.2016 by Jon Ramsey to avoid float overflows
!--------------------------------------------------------------
function bplanckdt(temp,nu)
  implicit none
  doubleprecision temp,nu
  doubleprecision theexp,bplanckdt,xx,yy
  doubleprecision, parameter :: largest = 1.7976879d308
  doubleprecision, parameter :: explargest = 709.78271d0
  xx     = 4.7989d-11*nu/temp
  ! prevents floating overflows in the exponent
  if(xx .gt. explargest) then
     bplanckdt = 0.0d0
     return
  endif
  theexp = exp(xx)
  if(theexp.lt.1.d33) then
     bplanckdt = 7.07661334104d-58 * nu**4 * theexp /  &
                 ( (theexp-1.d0)**2 * temp**2 ) + 1.d-290
  else
     yy = theexp * temp**2
     ! prevents floating outflows in the denominator
     if(yy .gt. largest) then
        bplanckdt = 0.0d0
        return
     endif
     bplanckdt = 7.07661334104d-58 * nu**4 / &
                 ( yy ) + 1.d-290
  endif
  return
end function bplanckdt


!-------------------------------------------------------------------
!                           RTGLOBAL CLEANUP
!-------------------------------------------------------------------
subroutine rtglobal_cleanup
  implicit none
  if(allocated(freq_nu)) deallocate(freq_nu)
  if(allocated(freq_dnu)) deallocate(freq_dnu)
  freq_nr = 0
  if(allocated(scat_mui_grid)) deallocate(scat_mui_grid)
  if(allocated(scat_thetai_grid)) deallocate(scat_thetai_grid)
  scat_munr = 0
  if(allocated(align_mui_grid)) deallocate(align_mui_grid)
  if(allocated(align_etai_grid)) deallocate(align_etai_grid)
  align_munr = 0
  if(allocated(cellvolume)) deallocate(cellvolume)
  if(allocated(cellindex)) deallocate(cellindex)
  if(allocated(dustdens)) deallocate(dustdens)
  if(allocated(dusttemp)) deallocate(dusttemp)
  if(allocated(dust_massdust)) deallocate(dust_massdust)
  if(allocated(electron_numdens)) deallocate(electron_numdens)
  if(allocated(ion_numdens)) deallocate(ion_numdens)
  !
  ! Cleanup the lines model stuff
  !
  if(allocated(lines_levelpop)) deallocate(lines_levelpop)
  if(allocated(gasvelocity)) deallocate(gasvelocity)
  !$OMP PARALLEL
  if(allocated(lines_ray_levpop)) deallocate(lines_ray_levpop)
  if(allocated(lines_ray_nrdens)) deallocate(lines_ray_nrdens)
  if(allocated(lines_ray_temp)) deallocate(lines_ray_temp)
  if(allocated(lines_ray_turb)) deallocate(lines_ray_turb)
  if(allocated(lines_ray_doppler)) deallocate(lines_ray_doppler)
  if(allocated(lines_ray_lorentz_delta)) deallocate(lines_ray_lorentz_delta)
  !$OMP END PARALLEL
  if(allocated(lines_microturb)) deallocate(lines_microturb)
  if(allocated(lines_escprob_lengthscale)) deallocate(lines_escprob_lengthscale)
  if(allocated(gas_chemspec_numberdens)) deallocate(gas_chemspec_numberdens)
  if(allocated(collpartner_numberdens)) deallocate(collpartner_numberdens)
  lines_maxveloc = 0.d0
  lines_maxturbc = 0.d0
  lines_maxtempc = 0.d0
  lines_maxshift = 0.d0
  !
  ! Cleanup polarization by alignment stuff
  !
  if(allocated(grainalign_dir)) deallocate(grainalign_dir)
  if(allocated(grainalign_eff)) deallocate(grainalign_eff)
  !
  ! Issue warning
  !
  !if(lines_warn_lineleap) then
  !   write(stdo,*) '*****************************************************'
  !   write(stdo,*) '  WARNING from the line module:'
  !   write(stdo,*) '  Insufficient spatial resolution: doppler shift '
  !   write(stdo,*) '  between some adjacent cells larger than the '
  !   write(stdo,*) '  local line width.'
  !   write(stdo,*) '*****************************************************'
  !endif
  !if(allocated(isoscatsrc)) deallocate(isoscatsrc)
  !call lines_cleanup()
  !call dust_cleanup()
end subroutine rtglobal_cleanup


!-------------------------------------------------------------------
!                  MAKE INDEXED FILE NAME
!-------------------------------------------------------------------
subroutine make_indexed_filename(base,index,ext,filename)
  implicit none
  character*80 base
  integer index
  character*80 ext
  character*80 filename
  character*12 ch
  !    
  if((index.lt.0).or.(index.ge.1000)) then
     write(stdo,*) 'ERROR in make_indexed_filename()'
     stop 729
  endif
  if(index.lt.10) then
     write(ch,11) index
11   format(I1)
  elseif(index.lt.100) then
     write(ch,12) index
12   format(I2)
  elseif(index.lt.1000) then
     write(ch,13) index
13   format(I3)
  else
     stop 5902
  endif
  filename = base(1:len_trim(base))//ch(1:len_trim(ch))//ext(1:len_trim(ext))
  !
end subroutine make_indexed_filename

!-------------------------------------------------------------------
!      MAKE INDEXED FILE NAME WITH 001, 011, 111 TYPE NUMBERS
!-------------------------------------------------------------------
subroutine make_indexed_filename_fill(base,index,ext,filename,len)
  implicit none
  character*80 base
  integer index,len,i
  character*80 ext
  character*80 filename
  character*12 ch
  !    
  if((len.lt.1).or.(len.gt.5)) stop 7009
  if((index.lt.0).or.(index.ge.100000)) then
     write(stdo,*) 'ERROR in make_indexed_filename_fill(). Number out of range:',index
     stop 729
  endif
  if(index.lt.10) then
     write(ch,11) index
11   format(I1)
     do i=1,len-1 
        ch = "0"//ch
     enddo
  elseif(index.lt.100) then
     if(len.lt.2) then
        write(stdo,*) 'ERROR while making filename: number too large:',index
        stop
     endif
     write(ch,12) index
12   format(I2)
     do i=1,len-2 
        ch = "0"//ch
     enddo
  elseif(index.lt.1000) then
     if(len.lt.3) then
        write(stdo,*) 'ERROR while making filename: number too large:',index
        stop
     endif
     write(ch,13) index
13   format(I3)
     do i=1,len-3
        ch = "0"//ch
     enddo
  elseif(index.lt.10000) then
     if(len.lt.4) then
        write(stdo,*) 'ERROR while making filename: number too large:',index
        stop
     endif
     write(ch,14) index
14   format(I4)
     do i=1,len-4
        ch = "0"//ch
     enddo
  elseif(index.lt.100000) then
     if(len.lt.5) then
        write(stdo,*) 'ERROR while making filename: number too large:',index
        stop
     endif
     write(ch,15) index
15   format(I5)
  else
     write(stdo,*) 'ERROR while making filename: number too large:',index
     stop
  endif
  filename = base(1:len_trim(base))//ch(1:len)//ext(1:len_trim(ext))
  !
end subroutine make_indexed_filename_fill

!-------------------------------------------------------------------
!           MAKE INDEXED FILE NAME WITH STRING-VALUED INDEX
!-------------------------------------------------------------------
subroutine make_indexed_filename_string(base,indexstring,ext,filename)
  implicit none
  character*80 base
  character*80 indexstring,finalstring
  character*80 ext
  character*80 filename
  integer :: len,ii
  !    
!!!!  len         = len_trim(indexstring)
!!!!  ii          = 1
!!!!  do while(indexstring(ii:ii).eq.' ')
!!!!     ii = ii + 1
!!!!     if(ii.gt.len) then
!!!!        write(*,*) 'ERROR in make_indexed_filename_string():'
!!!!        write(*,*) '      indexstring is empty: ',             &
!!!!                   indexstring(1:len_trim(indexstring))
!!!!        stop
!!!!     endif
!!!!  enddo
!!!!  len = len + 1 - ii
!!!!  if(len.lt.1) stop 7701
!!!!  finalstring(1:len) = indexstring(ii:ii+len-1)
!!!!  indexstring(1:len) = finalstring(1:len)
!!!!  do while(indexstring(len:len).eq.' ')
!!!!     len = len - 1
!!!!     if(len.lt.1) then 
!!!!        write(*,*) 'ERROR in make_indexed_filename_string():'
!!!!        write(*,*) '      indexstring is empty: ',             &
!!!!                   indexstring(1:len_trim(indexstring))
!!!!        stop
!!!!     endif
!!!!  enddo
  !
  len      = len_trim(indexstring)
  filename = base(1:len_trim(base))//indexstring(1:len)//ext(1:len_trim(ext))
  !
end subroutine make_indexed_filename_string

!--------------------------------------------------------------
!                   FUNCTION: TEST VALIDITY OF NUMBER
!
!     0 = Number is okay
!     1 = Number is INF
!     2 = Number is NAN
!--------------------------------------------------------------
function number_invalid(a)
  implicit none
  double precision :: a,b,c
  integer :: number_invalid
  logical :: div,sub
  
  b=a*2.d0
  b=b/2.d0
  c=a-1.d100
  
  div = (b.eq.a)
  sub = (c.lt.a)
  
  if(div.and.sub) then
     number_invalid = 0
  elseif(div) then
     number_invalid = 1
  else
     number_invalid = 2
  endif
  
  return
end function number_invalid




!-------------------------------------------------------------------------
!              GET ARGUMENT EITHER FROM ARGV OR FROM STDI
!-------------------------------------------------------------------------
subroutine ggetarg(iarg,buffer,fromstdi)
  use constants_module
  implicit none
  character*100 :: buffer
  integer :: iarg
  logical :: fromstdi
  if(fromstdi) then
     !
     ! Get the argument from standard input
     !
     read(ffli,*,end=103) buffer
  else
     !
     ! Get the argument from the command line
     !
     call getarg(iarg,buffer)
  endif
  return
103 continue
  buffer='enter'
  return
end subroutine ggetarg

!-------------------------------------------------------------------
!      INTEGER TO STRING, TRIMMED (ONLY FOR POSITIVE INTEGERS)
!-------------------------------------------------------------------
subroutine integer_to_string(int,string)
  implicit none
  character*80 :: string
  integer :: int
  if(int.lt.0) then
     write(string,*) int
  elseif(int.lt.10) then
     write(string,11) int
11   format(I1)
  elseif(int.lt.100) then
     write(string,12) int
12   format(I2)
  elseif(int.lt.1000) then
     write(string,13) int
13   format(I3)
  elseif(int.lt.10000) then
     write(string,14) int
14   format(I4)
  elseif(int.lt.100000) then
     write(string,15) int
15   format(I5)
  elseif(int.lt.1000000) then
     write(string,16) int
16   format(I6)
  elseif(int.lt.10000000) then
     write(string,17) int
17   format(I7)
  elseif(int.lt.100000000) then
     write(string,18) int
18   format(I8)
  else
     write(string,*) int
  endif
end subroutine integer_to_string
 

!-------------------------------------------------------------------------
!                 HELPER SUBROUTINE FOR READING DATA
!                         CASE: SCALAR FIELD
!-------------------------------------------------------------------------
subroutine read_scalarfield(unit,style,precis,nc,nv1,nv2,iv1,iv2, &
                            floor,reclen,scalar0,scalar1,scalar2)
  implicit none
  integer :: style,unit,irec,index,nv1,nv2,nc,iv1,iv2,i,ierr,precis,reclend
  integer(kind=8) :: iiformat,nn,kk
  integer, optional :: reclen
  double precision, optional :: floor
  double precision :: dummy,fl
  double precision, optional :: scalar0(1:nc)
  double precision, optional :: scalar1(1:nv1,1:nc)
  double precision, optional :: scalar2(1:nv1,1:nv2,1:nc)
  double precision, allocatable :: data(:)
  real*4 :: sdummy
  !
  ! Check precision (only relevant for style.eq.3)
  !
  if(style.eq.3) then
     if((precis.ne.4).and.(precis.ne.8)) then
        write(stdo,*) 'ERROR: Precision ',precis,' not known (must be 4 or 8)'
        stop
     endif
  endif
  !
  ! If a floor value is given, then use that
  !
  if(present(floor)) then
     fl=floor
  else
     fl=-1d99
  endif
  !
  ! Now read the scalar field
  !
  if(style.eq.1) then
     !
     ! Formatted (ASCII) input
     !
     call amr_resetcount()
     call amr_nextcell(index)
     do while(index.ge.0) 
        read(unit,*) dummy 
        if(index.gt.0) then
           if(present(scalar0)) scalar0(index) = max(dummy,fl)
           if(present(scalar1)) scalar1(iv1,index) = max(dummy,fl)
           if(present(scalar2)) scalar2(iv1,iv2,index) = max(dummy,fl)
        endif
        call amr_nextcell(index)
     enddo
  elseif(style.eq.2) then
     !
     ! Unformatted input, old fortran style (with records)
     !
     if(.not.present(reclen)) then
        write(stdo,*) 'ERROR in call to read_scalarfield(): '
        write(stdo,*) '      reclen not present in argument list.'
        stop
     endif
     reclend=reclen/8
     if(allocated(data)) deallocate(data)
     allocate(data(1:reclend),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate data array'
        stop
     endif
     call amr_resetcount()
     call amr_nextcell(index)
     do irec=1,(nrcellsinp-1)/reclend+1
        read(unit) data
        i = 1
        do while((index.ge.0).and.(i.le.reclend))
           if(index.gt.0) then
              if(present(scalar0)) scalar0(index) = max(data(i),fl)
              if(present(scalar1)) scalar1(iv1,index) = max(data(i),fl)
              if(present(scalar2)) scalar2(iv1,iv2,index) = max(data(i),fl)
           endif
           i     = i + 1
           call amr_nextcell(index)
        enddo
     enddo
  elseif(style.eq.3) then
     !
     ! Binary input: C-compliant unformatted streaming data
     !
     call amr_resetcount()
     call amr_nextcell(index)
     do while(index.ge.0) 
        if(precis.eq.4) then
           read(unit) sdummy 
           if(index.gt.0) then
              if(present(scalar0)) scalar0(index) = max(sdummy,fl)
              if(present(scalar1)) scalar1(iv1,index) = max(sdummy,fl)
              if(present(scalar2)) scalar2(iv1,iv2,index) = max(sdummy,fl)
           endif
        else
           read(unit) dummy 
           if(index.gt.0) then
              if(present(scalar0)) scalar0(index) = max(dummy,fl)
              if(present(scalar1)) scalar1(iv1,index) = max(dummy,fl)
              if(present(scalar2)) scalar2(iv1,iv2,index) = max(dummy,fl)
           endif
        endif
        call amr_nextcell(index)
     enddo
  else
     write(stdo,*) 'INTERNAL ERROR: Input style ',style,' not known.'
     stop
  endif
  if(allocated(data)) deallocate(data)
end subroutine read_scalarfield


!-------------------------------------------------------------------------
!                 HELPER SUBROUTINE FOR READING DATA
!                         CASE: VECTOR FIELD
!-------------------------------------------------------------------------
subroutine read_vectorfield(unit,style,precis,nv0,nv,nc,nv1,iv1,floor, &
                            reclen,vector0,vector1)
  implicit none
  integer :: style,unit,irec,index,nc,i,ierr,precis,nv,nv0,nv1,iv1,reclend
  integer(kind=8) :: iiformat,nn,kk
  integer, optional :: reclen
  double precision, optional :: floor
  double precision :: dummy,fl
  double precision, optional :: vector0(1:nv0,1:nc)
  double precision, optional :: vector1(1:nv0,1:nv1,1:nc)
  double precision, allocatable :: thedata(:)
  real*4, allocatable :: thesdata(:)
  double precision, allocatable :: data(:,:)
  !
  ! Backward compatibility: reclen
  !
  if(present(reclen)) then
     reclend = reclen/8
  else
     reclend = 0
  endif
  !
  ! Check precision (only relevant for style.eq.3)
  !
  if(style.eq.3) then
     if((precis.ne.4).and.(precis.ne.8)) then
        write(stdo,*) 'ERROR: Precision ',precis,' not known (must be 4 or 8)'
        stop
     endif
  endif
  !
  ! If a floor value is given, then use that
  !
  if(present(floor)) then
     fl=floor
  else
     fl=-1d99
  endif
  !
  ! Allocate array
  !
  allocate(thedata(1:nv),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR while reading a file: could not allocate thedata()'
     stop
  endif
  allocate(thesdata(1:nv),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR while reading a file: could not allocate thesdata()'
     stop
  endif
  !
  ! Now read the scalar field
  !
  if(style.eq.1) then
     !
     ! Formatted (ASCII) input
     !
     call amr_resetcount()
     call amr_nextcell(index)
     do while(index.ge.0) 
        read(unit,*) thedata(1:nv)
        do i=1,nv
           thedata(i) = max(fl,thedata(i))
        enddo
        if(index.gt.0) then
           if(present(vector0)) vector0(1:nv,index) = thedata(1:nv)
           if(present(vector1)) vector1(1:nv,iv1,index) = thedata(1:nv)
        endif
        call amr_nextcell(index)
     enddo
  elseif(style.eq.2) then
     !
     ! F77 style unformatted input (with records)
     !
     if(reclend.eq.0) then
        !
        ! The new way: 1 record = 1 cell
        !
        call amr_resetcount()
        call amr_nextcell(index)
        do while(index.ge.0) 
           read(unit) thedata(1:nv)
           do i=1,nv
              thedata(i) = max(fl,thedata(i))
           enddo
           if(index.gt.0) then
              if(present(vector0)) vector0(1:nv,index) = thedata(1:nv)
              if(present(vector1)) vector1(1:nv,iv1,index) = thedata(1:nv)
           endif
           call amr_nextcell(index)
        enddo
     else
        !
        ! The old way: 1 record has a length reclend
        ! 
        if(allocated(data)) deallocate(data)
        ! if(allocated(datas)) deallocate(datas)
        ! if(fex2) then
        reclend=reclen/(8*nv)
        allocate(data(1:nv,1:reclend),STAT=ierr)
        !else
        !reclend=reclen/(4*nv)
        !allocate(datas(1:nv,1:reclend),STAT=ierr)
        !allocate(data(1:nv,1:reclend),STAT=ierr)
        !endif
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate data array'
           stop
        endif
        !
        ! Now do the big loop
        !
        call amr_resetcount()
        call amr_nextcell(index)
        do irec=1,(nrcellsinp-1)/reclend+1
           !if(fex2) then
           read(1) data
           !else
           !read(1) datas
           !data(:,:) = datas(:,:)
           !endif
           i = 1
           do while((index.ge.0).and.(i.le.reclend))
              if(index.gt.0) then
                 if(present(vector0)) vector0(1:nv,index) = data(1:nv,i)
                 if(present(vector1)) vector1(1:nv,iv1,index) = data(1:nv,i)
              endif
              i     = i + 1
              call amr_nextcell(index)
           enddo
        enddo
     endif
  elseif(style.eq.3) then
     !
     ! Binary input: C-compliant unformatted streaming data
     !
     call amr_resetcount()
     call amr_nextcell(index)
     do while(index.ge.0) 
        if(precis.eq.4) then
           read(unit) thesdata(1:nv)
           thedata(1:nv) = thesdata(1:nv)
        else
           read(unit) thedata(1:nv)
        endif
        do i=1,nv
           thedata(i) = max(fl,thedata(i))
        enddo
        if(index.gt.0) then
           if(present(vector0)) vector0(1:nv,index) = thedata(1:nv)
           if(present(vector1)) vector1(1:nv,iv1,index) = thedata(1:nv)
        endif
        call amr_nextcell(index)
     enddo
  else
     write(stdo,*) 'INTERNAL ERROR: Input style ',style,' not known.'
     stop
  endif
  if(allocated(data)) deallocate(data)
  if(allocated(thedata)) deallocate(thedata)
  if(allocated(thesdata)) deallocate(thesdata)
end subroutine read_vectorfield


!-------------------------------------------------------------------------
!                 HELPER SUBROUTINE FOR WRITING DATA
!                         CASE: SCALAR FIELD
!-------------------------------------------------------------------------
subroutine write_scalarfield(unit,style,precis,nc,nv1,nv2,iv1,iv2,reclen, &
                             scalar0,scalar1,scalar2)
  implicit none
  integer :: style,unit,irec,index,nv,nc,iv,i,ierr,reclend,precis
  integer :: nv1,nv2,iv1,iv2
  integer(kind=8) :: iiformat,nn,kk
  integer, optional :: reclen
  double precision :: dummy
  double precision, allocatable :: data(:)
  double precision, optional :: scalar0(1:nc)
  double precision, optional :: scalar1(1:nv1,1:nc)
  double precision, optional :: scalar2(1:nv1,1:nv2,1:nc)
  real*4 :: sdummy
  !
  ! Check precision (only relevant for style.eq.3)
  !
  if(style.eq.3) then
     if((precis.ne.4).and.(precis.ne.8)) then
        write(stdo,*) 'ERROR: Precision ',precis,' not known (must be 4 or 8)'
        stop
     endif
  endif
  !
  ! Now write the scalar field
  !
  if(style.eq.1) then
     !
     ! Formatted (ASCII) output
     !
     call amr_resetcount()
     call amr_nextcell(index)
     do while(index.ge.0) 
        if(index.gt.0) then
           if(present(scalar0)) dummy = scalar0(index)
           if(present(scalar1)) dummy = scalar1(iv1,index)
           if(present(scalar2)) dummy = scalar2(iv1,iv2,index)
        else
           dummy = 0.d0
        endif
        write(unit,*) dummy 
        call amr_nextcell(index)
     enddo
  elseif(style.eq.2) then
     !
     ! Unformatted output, old fortran style (with records)
     !
     if(.not.present(reclen)) then
        write(stdo,*) 'ERROR in call to write_scalarfield(): '
        write(stdo,*) '      reclen not present in argument list.'
        stop
     endif
     reclend=reclen/8
     if(allocated(data)) deallocate(data)
     allocate(data(1:reclend),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate data array'
        stop
     endif
     call amr_resetcount()
     call amr_nextcell(index)
     do irec=1,(nrcellsinp-1)/reclend+1
        i = 1
        data(:) = 0.d0
        do while((index.ge.0).and.(i.le.reclend))
           if(index.gt.0) then
              if(present(scalar0)) data(i) = scalar0(index)
              if(present(scalar1)) data(i) = scalar1(iv1,index)
              if(present(scalar2)) data(i) = scalar2(iv1,iv2,index)
           else
              data(i) = 0.d0
           endif
           i       = i + 1
           call amr_nextcell(index)
        enddo
        write(1) data
     enddo
  elseif(style.eq.3) then
     !
     ! Binary output: C-compliant unformatted streaming data
     !
     call amr_resetcount()
     call amr_nextcell(index)
     do while(index.ge.0) 
        if(index.gt.0) then
           if(present(scalar0)) dummy = scalar0(index)
           if(present(scalar1)) dummy = scalar1(iv1,index)
           if(present(scalar2)) dummy = scalar2(iv1,iv2,index)
        else
           dummy = 0.d0
        endif
        if(precis.eq.4) then
           sdummy = dummy
           write(unit) sdummy
        else
           write(unit) dummy
        endif
        call amr_nextcell(index)
     enddo
  else
     write(stdo,*) 'INTERNAL ERROR: Output style ',style,' not known.'
     stop
  endif
  if(allocated(data)) deallocate(data)
end subroutine write_scalarfield



!-------------------------------------------------------------------------
!                 HELPER SUBROUTINE FOR WRITING DATA
!                         CASE: VECTOR FIELD
!-------------------------------------------------------------------------
subroutine write_vectorfield(unit,style,precis,nv0,nv,nc,nv1,iv1, &
                            vector0,vector1)
  implicit none
  integer :: style,unit,irec,index,nv,nc,ii,ierr,precis,nv1,iv1,nv0
  integer(kind=8) :: iiformat,nn,kk
  double precision :: dummy,fl
  double precision, optional :: vector0(1:nv0,1:nc)
  double precision, optional :: vector1(1:nv0,1:nv1,1:nc)
  double precision, allocatable :: thedata(:)
  real*4, allocatable :: thesdata(:)
  character*80 :: strnv,form
  !
  ! Check precision (only relevant for style.eq.3)
  !
  if(style.eq.3) then
     if((precis.ne.4).and.(precis.ne.8)) then
        write(stdo,*) 'ERROR: Precision ',precis,' not known (must be 4 or 8)'
        stop
     endif
  endif
  !
  ! Allocate array
  !
  allocate(thedata(1:nv),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR while reading a file: could not allocate thedata()'
     stop
  endif
  allocate(thesdata(1:nv),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR while reading a file: could not allocate thesdata()'
     stop
  endif
  !
  ! Now write the scalar field
  !
  if(style.eq.1) then
     !
     ! Formatted (ASCII) output
     !
     call amr_resetcount()
     call amr_nextcell(index)
     do while(index.ge.0) 
        if(index.gt.0) then
           if(present(vector0)) thedata(1:nv) = vector0(1:nv,index)
           if(present(vector1)) thedata(1:nv) = vector1(1:nv,iv1,index) 
        else
           thedata(1:nv) = 0.d0
        endif
        ! Bugfix 22.02.2016
        do ii=1,nv
           if(abs(thedata(ii)).lt.1d-98) thedata(ii)=0.d0
        enddo
        ! End Bugfix 22.02.2016
        call integer_to_string(nv,strnv)
        form='('//trim(strnv)//'(E13.6,1X))'
        write(unit,form) thedata(1:nv)
        call amr_nextcell(index)
     enddo
  elseif(style.eq.2) then
     !
     ! F77 style unformatted output (with records)
     !
     call amr_resetcount()
     call amr_nextcell(index)
     do while(index.ge.0) 
        if(index.gt.0) then
           if(present(vector0)) thedata(1:nv) = vector0(1:nv,index) 
           if(present(vector1)) thedata(1:nv) = vector1(1:nv,iv1,index) 
        else
           thedata(1:nv) = 0.d0
        endif
        write(unit) thedata(1:nv)
        call amr_nextcell(index)
     enddo
  elseif(style.eq.3) then
     !
     ! Binary output: C-compliant unformatted streaming data
     !
     call amr_resetcount()
     call amr_nextcell(index)
     do while(index.ge.0) 
        if(index.gt.0) then
           if(present(vector0)) thedata(1:nv) = vector0(1:nv,index) 
           if(present(vector1)) thedata(1:nv) = vector1(1:nv,iv1,index) 
        else
           thedata(1:nv) = 0.d0
        endif
        if(precis.eq.4) then
           thesdata(1:nv) = thedata(1:nv)
           write(unit) thesdata(1:nv)
        else
           write(unit) thedata(1:nv)
        endif
        call amr_nextcell(index)
     enddo
  else
     write(stdo,*) 'INTERNAL ERROR: Output style ',style,' not known.'
     stop
  endif
  if(allocated(thedata)) deallocate(thedata)
  if(allocated(thesdata)) deallocate(thesdata)
end subroutine write_vectorfield


end module rtglobal_module

