program radmc3d
  use amr_module
  use rtglobal_module
  use dust_module
  use lines_module
  use ioput_module
  use camera_module
  use constants_module
  use mathroutines_module
  use montecarlo_module
  use userdef_module
  use diffusion_module
  implicit none
  !
  ! Local variables for main program
  !
  integer*8 :: nphot,nphot_scat_original
  integer :: inu,iformat,irestart,ir,it,idiff,ifill
  integer :: ispec,cntdump,ifile,iredo,i,l
  integer :: ntemp,nvstr,ivstrt,iterstr,iwarn,ierror_diff,ierr
  integer :: ierror
  doubleprecision :: temp0,temp1,error,incrraycost
  logical :: fex1,fex2,iqy,gotit,quit,selfcheck,inclines_bk
  integer :: maxtasks,itask,nv
  character*160 :: subboxfilename,samplefilename,species
  character*160 :: vtkfilename,vtktitle
  integer :: npt,imx,imy,iinu,nvm,bc_idir,bc_ilr
  double precision, pointer :: xpt(:),ypt(:),zpt(:)
  ! Attila Juhasz  - vtkFieldName is the name of the scalar / vector field in
  ! the VTK file
  character*80 :: vtkFieldName
  !
  ! Default safety settings
  ! 
  mc_safetymode = 0         ! Use fast Monte Carlo, though still accurate enough
  !                         ! For debugging or checking results, try a run with
  !                         ! mc_safetymode = 1 or higher (the higher the safer;
  !                         ! max=999)
  debug_check_all = 0       ! Debugging mode off by default. But if safetymode
  !                         ! is switched on (below), it will be switched on
  !                         ! as well.
  amr_always_use_tree=.false.
  !
  ! Default Standard Input/Output units
  !
  stdi = 5
  stdo = 6
  ffli = 5
  fflo = 1
  !
  ! Defaults for the actions
  !
  do_montecarlo_therm        = .false.
  do_montecarlo_mono         = .false.
  do_userdef_action          = .false.
  do_vstruct                 = .false.
  do_raytrace_spectrum       = .false.
  do_raytrace_image          = .false.
  do_raytrace_movie          = .false.
  do_raytrace_tausurf        = .false.
  do_writeimage              = .false.
  do_writespectrum           = .false.
  do_writescatsrc            = .false.
  do_resetseed               = .false.
  do_write_lineinfo          = .false.
  do_writemodel              = .false.
  do_writegridfile           = .false.
  do_respondwhenready        = .false.
  do_writepop                = .false.
  do_calcpop                 = .false.
  !
  ! Processes to include
  !
  rt_incl_dust               = .false.    ! This is now different!
  rt_incl_lines              = .false.
  rt_incl_gascont            = .false.
  rt_incl_gascont_freefree   = .true.
  rt_incl_userdef_srcalp     = .false.
  !
  dust_2daniso               = .false.
  dust_2daniso_nphi          = 360 
  !
  ! Defaults for the camera position and orientation. 
  !
  writeimage_unformatted     = .false.
  circular_images            = .false.
  camera_observer_degr_theta = 0.d0          ! View pole-on
  camera_observer_degr_phi   = 0.d0
  camera_observer_distance_pc= 0.d0
  camera_pointing_position(1:3) = 0.d0
  camera_pointing_degr_posang= 0.d0 
  camera_image_halfsize_x    = -101.
  camera_image_halfsize_y    = -101.
  camera_zoomcenter_x        = 0.d0
  camera_zoomcenter_y        = 0.d0
  camera_localobs_zenith     = 0.d0
  camera_image_nx            = 100
  camera_image_ny            = 100
  camera_nrfreq              = 0
  camera_theinu              = 0
  camera_lambdamic           = 0.d0
  camera_lambdamic1          = -1.d0
  camera_range_nlam          = -1
  camera_tracemode           = 1             ! Normal mode: real continuum emis
  camera_loadcolor           = .false.
  camera_loadlambdas         = .false.
  camera_setfreq_global      = .false.
  camera_nrrefine            = 100
  camera_incl_stars          = 1
  camera_refine_criterion    = 1.d0          ! The smaller, the safer (but heavier)
  camera_scatsrc_allfreq     = .false.
  camera_starsphere_nrpix    = 20
  camera_localobserver       = .false.
  camera_observer_position(1:3) = 1.d99
  camera_localobs_projection = 1
  camera_use_aperture_info   = .false.
  camera_spher_cavity_relres = 0.05d0
!  camera_min_aspectratio     = 0.05d0
  camera_max_dangle          = 0.3d0
  camera_min_dangle          = 0.05d0
  camera_min_drr             = 0.003d0
  camera_secondorder         = .false.
  sources_interpol_jnu       = .true.
  camera_catch_doppler_line  = .false.
  camera_catch_doppler_resolution =0.2d0
  camera_diagnostics_subpix  = .false.
  camera_stokesvector        = .false.
  lines_autosubset           = .true.
  camera_lambda_starlight_single_scat_mode    = 0
  camera_maxdphi             = 0.1d0   ! New as of March 2017
  !
  ! Line wavelength grid parameter settings
  ! 
  lines_user_widthkms        = 0.d0
  lines_user_kms0            = 0.d0
  lines_user_ispec           = 0
  lines_user_iline           = 0
  lines_user_nrfreq          = 1
  tgas_eq_tdust              = 0
  lines_show_pictograms      = 0
  lines_maser_warning        = .false.
  lines_slowlvg_as_alternative = .false.
  !
  ! Thermal boundary conditions (only for 
  ! cartesian coordinates)
  !
  thermal_bc_active(1:2,1:3) = .false.
  thermal_bc_temp(1:2,1:3)   = 0.d0
  !
  ! Defaults for I/O parameters
  !
  rto_style                 = 1
  rto_single                = .false.
  rto_reclen                = 1024           ! Standard record length for unformatted output
  rtio_gridinfo             = 0
  maxtasks                  = 10000
  grid_was_read_from_file   = .false.
  grid_was_written_to_file  = .false.
  !
  ! Defaults for the parameters for the grid
  !
  igrid_type                = 0
  igrid_coord               = 1
  igrid_mirror              = 0
  !
  ! Other global defaults
  ! 
  lines_mode                = 1        ! Simple LTE mode is default
  lines_partition_ntempint  = 1000     ! For internally computed part func: nr of temps
  lines_partition_temp0     = 0.1d0    ! For internally computed part func: temp0
  lines_partition_temp1     = 1d5      ! For internally computed part func: temp1
  lines_make_linelist       = .false.  ! Do not write out a line list (default)
  lines_pfunc_outofrange    = .false.  ! A warning flag
  lines_widthmargin         = 12.d0    ! The margin by which the doppler jumping condition is checked
  lines_profile             = 0        ! Gaussian line profile
  lines_nonlte_maxiter      = 100      ! Default for the max nr of iterations for non-LTE
  lines_nonlte_convcrit     = 1d-2     ! Default convergence criterion for non-LTE
  scattering_mode           = 0        ! Just an initial setting for the scattering mode
  scattering_mode_def       = 0        ! Just an initial setting for the scattering mode
  scattering_mode_max       = 9999     ! Per default no limit to scattering mode
  alignment_mode            = 0        ! Per default no alignment of grains
  star_sphere               = .false.  ! Per default, stars are treated as point-sources
  radmc_as_child            = .false.
  !
  !!$ Defaults for the parallel version of RADMC-3D
  !$ setthreads=1
  !
  ! Defaults for the parameters for the Monte Carlo
  !
  iseed_start               = -17933201
  iseed                     = iseed_start
  rt_mcparams%nphot_therm      = 100000
  rt_mcparams%nphot_scat       = 100000
  rt_mcparams%nphot_spec       = 10000
  rt_mcparams%nphot_mono       = 100000
  rt_mcparams%ifast            = 0
  rt_mcparams%iranfreqmode     = 0          ! Default: use no interpol in temp 
  !                                         ! in finding new random freq
  !                                         ! in thermal Monte Carlo. NOTE:
  !                                         ! This is new compared to versions
  !                                         ! 0.10 and earlier. That is why
  !                                         ! in some of the selftest models
  !                                         ! this must be explicitly 
  !                                         ! switched on.
  rt_mcparams%enthres          = 1.d-2
  rt_mcparams%irestart         = 0
  rt_mcparams%cntdump          = 10000000
  rt_mcparams%ntemp            = 1000
  rt_mcparams%itempdecoup      = 1
  rt_mcparams%temp0            = 0.01d0
  rt_mcparams%temp1            = 1d5
  rt_mcparams%extinct_thres    = 10.d0
!!!  rt_mcparams%save_scat        = 1
  rt_mcparams%debug_write_stats= 0
  rt_mcparams%debug_write_path = 0
!!  rt_mcparams%debug_write_eventcounts = 0
  rt_mcparams%countwrite       = 1000
  rt_mcparams%ivstrt           = 1       ! Dust species 1 used for vert struct
  rt_mcparams%vserrtol         = 0.d0
  rt_mcparams%niter_vstruct    = 0
  rt_mcparams%errtoldiff       = 1d-10
  rt_mcparams%nphotdiff        = 0
  rt_mcparams%nphotdiff_type   = 1
  !rt_mcparams%incl_scatsrc_mctherm= 0    ! The default is not to compute the 
  !!                                      ! scat source in Bjork&Wood
  rt_mcparams%optimized_motion = .false. ! By default switch optimized photon
  !                                      ! motion off. This is perhaps slower
  !                                      ! but safer.
  rt_mcparams%optim_dtau       = 2.d0
  !
  ! Parameters for the Modified Random Walk:
  !
  rt_mcparams%mod_random_walk  = .false. ! Default is no MRW
  rt_mcparams%mrw_count_trigger= 400     ! How often remain in same cell before MRW
  rt_mcparams%mrw_db_ntemp     = 1000    ! Nr of temperatures of the planck mean database
  rt_mcparams%mrw_db_temp0     = 0.01d0  ! Lowest temperature
  rt_mcparams%mrw_db_temp1     = 1d5     ! Highest temperature
  rt_mcparams%mrw_enthres      = 0.3     ! Fraction by which energy should increase to trigger new opac
  rt_mcparams%mrw_tauthres     = 20.d0   ! Minimal optical depth a cell has to have for MRW
  rt_mcparams%mrw_gamma        = 4.d0    ! The MRW gamma value (see Min et al. 2009)
  rt_mcparams%mrw_taustepback  = 0.001   ! Amount of optical depth to step back at end of MRW
  rt_mcparams%mrw_tempthres    = 3.d0    ! Below T=3K we do not recalculate the mean opacities for MRW
  mcscat_nrdirs = 0 
  mc_scat_maxtauabs            = 30.d0   ! Default value
  mc_max_nr_scat_events        = -1      ! Default value. If >=0 then limit nr of scattering events treated in images.
  mc_weighted_photons          = .true.
  scat_munr                    = 0       ! The dust opacity files will set the nr of scattering angles mu
  !scat_phinr                   = 180     ! The default nr of scattering angles in phi
  !
  ! Which non-standard sources of energy or photons to include?
  !
  incl_stellarsrc           = 0
  incl_extlum               = 0
  incl_heatsource           = 0
  incl_quantum              = 0
  incl_thermbc              = 0
  !
  ! Subbox parameters
  !
  subbox_nx = 64
  subbox_ny = 64
  subbox_nz = 64
  subbox_x0 = 0.d0
  subbox_y0 = 0.d0
  subbox_z0 = 0.d0
  subbox_x1 = 0.d0
  subbox_y1 = 0.d0
  subbox_z1 = 0.d0
  subbox_phi1 = 0.d0
  subbox_phi2 = 0.d0
  subbox_theta = 0.d0
  do_write_subbox_dustdens = .false.
  do_write_subbox_dusttemp = .false.
  do_write_subbox_levelpop = .false.
  !
  ! Sampling parameters
  !
  do_writesample     = .false.
  do_sample_dustdens = .false.
  do_sample_dusttemp = .false.
  do_sample_levelpop = .false.
  !
  ! VTK output (VTK=Visual Tool Kit: www.vtk.org/)
  !
  do_write_vtk_grid             = .false.
  do_write_vtk_dust_density     = .false.
  do_write_vtk_dust_temperature = .false.
  do_write_vtk_gas_density      = .false.
  do_write_vtk_gas_temperature  = .false.
  do_write_vtk_chemspec         = .false.
  do_write_vtk_levelpop         = .false.
  do_write_vtk_velocity         = .false.
  !
  ! Now check out which processes are likely to be included. This is
  ! done by checking the existence of some files. 
  !
  inquire(file='dustopac.inp',exist=fex1)
  if(fex1) rt_incl_dust = .true.
  inquire(file='lines.inp',exist=fex1)
  if(fex1) rt_incl_lines = .true.
  !
  ! Call the userdef defaults routine
  !
  call userdef_defaults()
  !
  ! Very first action after defaults is: interpret command line options
  ! 
  call interpet_command_line_options(gotit,.false.,quit)
  if(.not.gotit) goto 700
  !
  ! If radmc3d is meant to run as a child process of a parent process,
  ! and the communication to/from the parent goes via the STDIO, then
  ! we must redirect the 'normal' STDIO.
  !
  if(radmc_as_child) then
     stdo = 110
     open(unit=stdo,file='radmc3d.out')
     fflo = 6
  else
     open(unit=1,file='radmc3d.out')
     write(1,*) 'Standard output went to screen...'
     close(1)
  endif
  !
  ! Write a banner:
  !
  call write_banner()
  !
  ! Now check consistency of some of the input files.
  !
  inquire(file='dustopac.inp',exist=fex1)
  if(.not.fex1) then
     inquire(file='dust_density.inp',exist=fex2)
     if(fex2) then
        write(stdo,*) 'WARNING: dust_density.inp found, but no dustopac.inp...'
     endif
     inquire(file='dust_density.uinp',exist=fex2)
     if(fex2) then
        write(stdo,*) 'WARNING: dust_density.uinp found, but no dustopac.inp...'
     endif
     inquire(file='dust_temperature.dat',exist=fex2)
     if(fex2) then
        write(stdo,*) 'WARNING: dust_temperature.dat found, but no dustopac.inp...'
     endif
     inquire(file='dust_temperature.udat',exist=fex2)
     if(fex2) then
        write(stdo,*) 'WARNING: dust_temperature.udat found, but no dustopac.inp...'
     endif
  endif
  inquire(file='lines.inp',exist=fex1)
  if(.not.fex1) then
     inquire(file='gas_velocity.inp',exist=fex2)
     if(fex2) then
        write(stdo,*) 'WARNING: gas_velocity.inp found, but no lines.inp...'
     endif
     inquire(file='gas_velocity.uinp',exist=fex2)
     if(fex2) then
        write(stdo,*) 'WARNING: gas_velocity.uinp found, but no lines.inp...'
     endif
     inquire(file='microturbulence.inp',exist=fex2)
     if(fex2) then
        write(stdo,*) 'WARNING: microturbulence.inp found, but no lines.inp...'
     endif
     inquire(file='microturbulence.uinp',exist=fex2)
     if(fex2) then
        write(stdo,*) 'WARNING: microturbulence.uinp found, but no lines.inp...'
     endif
  endif
  !
  ! Check if the main input file exists. If it exists, then read it.
  ! Otherwise make an emergency stop.
  !
  call read_radmcinp_file()
  !
  ! Check that SED and lines are not on simultaneously
  !
  if(rt_incl_lines.and.do_raytrace_spectrum.and.camera_setfreq_global) then
     write(stdo,*) 'ERROR: Cannot make SED with line mode on.'
     stop
  endif
  !
  ! If OMP, then immediately set the nr of threads here, now that we
  ! know which value it should have: the setthreads variable.
  ! Bugfix by Jon Ramsey 22.03.2016
  !
  !$ call OMP_set_num_threads(setthreads)
  !$ nprocs=OMP_get_num_procs()
  !$ write(*,*) 'Number of processors: ',nprocs
  !$ write(*,*) 'Number of threads in use: ',setthreads
  !
  !$ if(setthreads.gt.nprocs) then
  !$		write(*,*) 'Number of threads exceeds the number of processors which are available.'
  !$		write(*,*) 'Please reduce the value of the "setthreads" argument.'
  !$		STOP 3189
  !$ endif
  !
  !$ if(setthreads.eq.0) then
  !$		write(*,*) 'ERROR: Minimum number of threads must be 1.'
  !$		STOP 3189
  !$ endif
  !
  ! Do some checks
  !
!   if(rt_mcparams%nphotdiff.ne.0) then
!      write(stdo,*) 'ERROR: For now, the diffusion mode is not yet ready'
!      write(stdo,*) '       in this 3-D version of RADMC. Use RADMC instead.'
!      stop
!   endif
  if(incl_quantum.eq.-1) then
     write(stdo,*) 'WARNING: PAH emission not included as '
     write(stdo,*) '         heating source for thermal grains.'
     write(stdo,*) '         Energy is therefore not well conserved.'
     write(stdo,*) '         Put incl_quantum to -2 to include it,'
     write(stdo,*) '         and run RADMC-3D both before and after'
     write(stdo,*) '         the external PAH code...'
  endif
  if(rt_mcparams%ifast.eq.1) then
     write(stdo,*) 'NOTE: Using acceleration efficiency 1'
     write(stdo,*) '      (=approximate method!)'
  elseif(rt_mcparams%ifast.eq.2) then
     write(stdo,*) 'NOTE: Using acceleration efficiency 2'
     write(stdo,*) '      (=approximate method!)'
  endif
  if(rt_mcparams%nphotdiff.ne.0) then
     write(stdo,*) 'NOTE: Using diffusion algorithm for cells with '
     write(stdo,*) '      less than ',rt_mcparams%nphotdiff,' photons visited'
  endif
  if(tgas_eq_tdust.lt.0) then
     write(stdo,*) 'ERROR: Flag tgas_eq_tdust cannot be <0'
     stop
  endif
  if(tgas_eq_tdust.gt.0) then
     write(stdo,*) 'Note: T_gas taken to be equal to T_dust of dust species ',&
          tgas_eq_tdust
  endif
  !
  ! Make sure the iseed is negative 
  !
  if(iseed.gt.0) iseed=-iseed
  iseed_start = iseed
  !
  ! Some post-processing after command-line reading
  ! 
  if(do_vstruct) then
     do_montecarlo_therm=.true.
  else
     rt_mcparams%niter_vstruct = 0
  endif
  if(rt_mcparams%niter_vstruct.ne.0) then
     write(stdo,*) 'ERROR: For now, the vertical structure iteration is not yet ready'
     write(stdo,*) '       in this 3-D version of RADMC. Use RADMC instead.'
     stop
  endif
  if(rt_mcparams%niter_vstruct.ne.0) then
     inquire(file='chopdens.inp',exist=fex1)
     if(fex1) then
        write(stdo,*) 'DANGER: The file chopdens.inp is detected'
        write(stdo,*) '        (making me assume that the chopdens'
        write(stdo,*) '        code has been employed)'
        write(stdo,*) '        while the vertical structure mode'
        write(stdo,*) '        of RADMC is switched on. These'
        write(stdo,*) '        two things are strongly mutually'
        write(stdo,*) '        exclusive. Prove to me that no'
        write(stdo,*) '        chopping has been done by removing'
        write(stdo,*) '        the chopdens.inp file. '
        write(stdo,*) '  NOTE: Maybe in the future this is fixed.'
        stop
     endif
  endif
  if(alignment_mode.ne.0) then
     if(.not.camera_stokesvector) then
        camera_stokesvector = .true.
        write(stdo,*) 'Setting camera_stokesvector=.true. because ', &
             'aligned grains mode requires full Stokes treatment.'
     endif
  endif
  !
  ! Handle safety settings
  !
  if(mc_safetymode.gt.0) then
     !
     ! Set some Monte Carlo settings to safe mode, which may result in
     ! slower calculation, but without the danger that some of the optimizations
     ! may yield wrong results.
     !
     write(stdo,*) '**** USING SAFE MODE FOR MONTE CARLO: SLOWER BUT SAFER! ****'
     write(stdo,*) '**** DID YOU COMPILE WITH -FBOUNDS-CHECK OR -C TO CHECK ****'
     write(stdo,*) '**** FOR READING/WRITING OUTSIDE OF ARRAY BOUNDARIES?   ****'
     write(stdo,*) '**** THAT CAN BE USEFUL AS WELL FOR DEBUGGING/CHECKING. ****'
     !
     ! ...First make sure that we use a linearly interpolated cumulative
     !    table in the pick_randomfreq_db() routine in the montecarlo module,
     !    instead of using the nearest temperature in the database. This
     !    should not make much difference if db_ntemp is large enough, but
     !    let us be careful nevertheless.
     !
     rt_mcparams%iranfreqmode     = 1 
     !
     ! ...Also for safety set the nr of temperature points to a very high
     !    value...
     !
     rt_mcparams%ntemp            = 10000
     !
     ! ...Do not use optimized transport of photons. Just use basic methods.
     !
     rt_mcparams%optimized_motion = .false.
     !
     ! ...Do all checks possible
     !
     debug_check_all = 1
     !
  endif
  !
  ! Some postprocessing from userdef?
  !
  call userdef_main_namelist_postprocessing()
  !
  ! Set the coordinate system to regular
  !
  igrid_coord  = -1    ! Put it for now at invalid value
  igrid_mirror = 0
  !
  ! Maybe the userdef module wants to read/set the frequencies
  ! and the basic grid. Allow it here.
  !
  call userdef_prep_model()
  !
  ! Read the frequency array (frequency = c/lambda)
  ! If it is already set, then don't modify.
  !
  call read_frequencies(1)
  ! 
  ! Read the external (insterstellar) radiation field
  !
  !!!! call read_interstellar_radfield()
  !
  ! Now read the grid file, which automatically determines (from which
  ! *_grid.inp file is present) which grid type is used. If the grid is
  ! already set (by the userdef module) then it won't modify it.
  !
  if(.not.allocated(amr_grid_xi)) then
     !
     write(stdo,*) 'Reading grid file and prepare grid tree...'
     call flush(stdo)
     call read_grid_file()
     !
  endif
  !
  ! Read the star spectrum
  !
  if(.not.allocated(star_spec)) then
     write(stdo,*) 'Reading star data...'
     call flush(stdo)
     call read_stars()
  endif
  !
  ! Now initialize the grid and it makes the linear list of cells needed for
  ! loading the grid-based data. 
  !
  ! NOTE: If the user adapts the grid in the userdef_setup_model() routine,
  !       he/she should not forget to call this once more, this time with
  !       argument renew=.true..  Also amr_compute_list_all() must be
  !       re-called (before the postprocess_grid_for_use()) call.
  !
  call postprocess_grid_for_use()
  !
  ! Do a check: Thermal boundaries are only allowed in cartesian 
  ! coordinates.
  !
  if(igrid_coord.ge.100) then
     do bc_idir=1,3
        do bc_ilr=1,2
           if(thermal_bc_active(bc_ilr,bc_idir)) then
              write(stdo,*) 'ERROR: Thermal boundaries are only allowed in cartesian coordinates.'
              write(stdo,*) '       please remove any statements like thermal_boundary_xl = 30 or'
              write(stdo,*) '       so from the radmc3d.inp file.'
              stop
           endif
        enddo
     enddo
  endif
  !
  ! Now allow the userdef module to set up some own model stuff.
  ! Any arrays that will be set by the userdef_setup_model() routine
  ! here will automatically not be read anymore later. 
  !
  call userdef_setup_model()
  !
  ! Read the continuous stellar sources if present, and if the userdef
  ! module has not yet set it up.
  !
  call read_stellarsource(1)
  if(incl_stellarsrc.ne.0) then
     write(stdo,*) 'Using smooth stellar source distributions...'
     call flush(stdo)
  endif
  !
  ! Read the internal heat source if present, and if the userdef
  ! module has not yet set it up.
  ! 
  call read_internal_heatsource(1)
  if(incl_heatsource.ne.0) then
     write(stdo,*) 'Using internal heat source distribution...'
     call flush(stdo)
  endif
  !
  ! Read the external radiation field (interstellar radiation field) if
  ! preset, and if the userdef module has not yet set it up.
  !
  call read_externalsource(1)
  if(incl_extlum.ne.0) then
     write(stdo,*) 'Using external luminosity (interstellar radiation field)...'
     call flush(stdo)
  endif
  !
  ! Displace all stars a tiny little bit just to avoid the risk of
  ! having a star *exactly* on a grid wall, which would cause trouble.
  !
  call jitter_stars(1.d0)
  !
  ! Find out in which cell each star is. This must be done here, not
  ! earlier, because the grid may have changed in the userdef_setup_model()
  ! call.
  !
  call stars_findcell()
  !
  ! Initialize the grid
  !
  if((igrid_type.ge.0).and.(igrid_type.lt.100)) then
     !
     ! AMR-type grid
     !
     ! Initialize the AmrRay module, so that raytracing can be done on this grid
     !
     !if(debug_check_all.ne.0) then
     !   selfcheck = .true.
     !else
     !   selfcheck = .false.
     !endif
     !call amrray_initialize(selfcheck)
     !
     ! NOTE: It has turned out to be safer to keep the self-checking switched
     !       on, even though it might be a bit slower. But we use the motto:
     !       safety first!
     !
     call amrray_initialize(.true.)
     !
     ! If mirror symmetry is used, then make this known to the rest of the
     ! code
     !
     if(amrray_mirror_equator) then
        igrid_mirror = 1
     endif
  else
     stop 1400
  endif
  !
  ! If lines are included in the model, we now check for the maximum
  ! cell-to-cell line shift compared to the local line width. 
  !
  if(rt_incl_lines.and.allocated(gasvelocity).and. &
       allocated(lines_microturb).and.allocated(gastemp)) then
     lines_maxrelshift = 0.d0
     call lines_compute_maxrellineshift()
     if(lines_maxrelshift.gt.0.7d0) then
        write(stdo,*) 'WARNING: cell-to-cell line doppler shifts are large compared to the local line width.'
        write(stdo,*) '   In this model the largest dv_doppler/dv_linewidth = ',lines_maxrelshift
        if(camera_catch_doppler_line) then
           write(stdo,*) '   But since you are using doppler catching mode, you are reasonably safe...'
        else
           write(stdo,*) '   To get reasonable line images/spectra please use the doppler catching mode.'
        endif
     endif
  endif
  !
  !
  ! ====================================================================
  !                     NOW START THE LOOP OVER ACTIONS
  ! ====================================================================
  !
  ! If we are in a "child process" mode, then here is the point where we
  ! must return after each action. 
  !
  itask = 1
  quit = .false.
  do while(.not.quit) 
     !
     ! If in child process, then receive commands from the parent.
     ! Otherwise re-read the command-line options from the command line, so
     ! that they overwrite the radmc3d.inp file settings or other settings 
     ! if necessary.
     ! 
     if(radmc_as_child) then
        write(stdo,*) 'Waiting for commands via standard input....'
        call flush(stdo)
        call interpet_command_line_options(gotit,.true.,quit)
        write(stdo,*) 'OK, received commands via standard input, now go to work.'
        call flush(stdo)
     else
        call interpet_command_line_options(gotit,.false.,quit)
     endif
     !
     ! Do some preparations and checks for the images or spectra
     !
     if(do_raytrace_image.or.do_raytrace_spectrum.or.do_raytrace_tausurf) then
        !
        ! If one of the camera actions is to be taken, then make sure
        ! that the image size specifications are set to useful values.
        !
        if((camera_image_halfsize_x.le.-100.d0).or.  &
           (camera_image_halfsize_y.le.-100.d0)) then
           !
           ! Which grid type
           !
           if(igrid_type.lt.100) then
              !
              ! AMR grid
              !
              if(igrid_coord.lt.100) then
                 !
                 ! Cartesian coordinates
                 !
                 if(igrid_coord.eq.10) then
                    !
                    ! 1-D plane-parallel mode
                    !
                    camera_image_halfsize_x = 0.5d0
                    camera_image_halfsize_y = 0.5d0
                 elseif(igrid_coord.eq.20) then
                    !
                    ! 2-D pencil-parallel mode
                    !
                    camera_image_halfsize_x = 1.01*sqrt(3.)*                      &
                      0.5*max(abs(amr_grid_xi(1,2)-amr_grid_xi(amr_grid_ny+1,2)), &
                              abs(amr_grid_xi(1,3)-amr_grid_xi(amr_grid_nz+1,3)))
                    camera_image_halfsize_y = camera_image_halfsize_x
                 else
                    !
                    ! Full 3-D mode
                    !
                    camera_image_halfsize_x = 1.01*sqrt(3.)*                      &
                      0.5*max(abs(amr_grid_xi(1,1)-amr_grid_xi(amr_grid_nx+1,1)), &
                              abs(amr_grid_xi(1,2)-amr_grid_xi(amr_grid_ny+1,2)), &
                              abs(amr_grid_xi(1,3)-amr_grid_xi(amr_grid_nz+1,3)))
                    camera_image_halfsize_y = camera_image_halfsize_x
                 endif
              elseif(igrid_coord.lt.200) then
                 !
                 ! Spherical coordinates
                 !
                 camera_image_halfsize_x = 1.1*amr_grid_xi(amr_grid_nx+1,1)
                 camera_image_halfsize_y = camera_image_halfsize_x
              else
                 stop 77
              endif
           else
              stop 78
           endif
        endif
     endif
     !
     ! If vertical structure mode is on, then do some stuff
     ! 
     if(rt_mcparams%niter_vstruct.gt.0) then
        !
        ! Calculate the surface densities of the dust (for the vertical
        ! structure calculation, and for some diagnostics)
        !
        !############## FOR NOW THIS IS DISABLED ################
        !
        !!!!! call calc_sigmadust()
     endif
     !
     ! Read the vertical structure input file here, if necessary
     !
     ! ############### NOT YET READY ################
     !
     ! If quantum-heated grains are present, but incl_quantum is not set,
     ! then make a big warning
     !
     iqy = .false.
     do ispec=1,dust_nr_species
        if(dust_quantum(ispec).ne.0) iqy=.true.
     enddo
     if(iqy.and.(incl_quantum.eq.0)) then
        write(stdo,*) '***********************************************'
        write(stdo,*) '***********************************************'
        write(stdo,*) '  WARNING: Treatment of quantum-heated grains'
        write(stdo,*) '  switched off, but some of the opacities are'
        write(stdo,*) '  marked as quantum-heated grains... They are'
        write(stdo,*) '  now treated as thermal grains...'
        write(stdo,*) '***********************************************'
        write(stdo,*) '***********************************************'
     endif
     !
     !
     !
     !
     !================================================================
     !                       NOW EXECUTE THE COMMANDS
     !================================================================
     !
     !
     !----------------------------------------------------------------
     !          DO THE THERMAL MONTE CARLO (BJORKMAN & WOOD)
     !----------------------------------------------------------------
     !
     if(do_montecarlo_therm) then
        !
        ! Line mode must be temporarily switched off in the thermal
        ! Monte Carlo method. 
        !
        inclines_bk = rt_incl_lines
        if(rt_incl_lines) then
           write(stdo,*) 'NOTE: Lines are included in this model. But for the'
           write(stdo,*) '      thermal Monte Carlo run they are swithed off...'
        endif
        rt_incl_lines  = .false.
        !
        ! If gas continuum opacities are included, then warn that (at least
        ! for now) these are not iterated along, i.e. the temperature remains
        ! fixed.
        !
        if(rt_incl_gascont) then
           write(stdo,*) 'NOTE: Gas continuum is included in this model. But for the'
           write(stdo,*) '      thermal Monte Carlo run this is kept fixed, i.e.'
           write(stdo,*) '      the gas temperature is not adjusted.'
           write(stdo,*) '      In a later version of RADMC-3D we plan to allow'
           write(stdo,*) '      that also the gas temperature can be iterated on.'
        endif
        !
        ! Set the camera_frequencies(:) array, if requested
        !
        !call set_camera_frequencies()
        !
        ! Do iteration over the vertical structure
        !
        do iterstr=0,rt_mcparams%niter_vstruct
           !
           ! Do the MC according to the algorithm of Bjorkman & Wood
           !
           call do_monte_carlo_bjorkmanwood(rt_mcparams,ierr,resetseed=do_resetseed)
           !
           ! Messages
           !
           write(stdo,*) 'Done with Bjorkman & Wood Monte Carlo simulation'
           !$ write(stdo,*) 'Using the parallel version of RADMC-3D'
           !$ write(stdo,*) 'Number of parallel threads used: ', setthreads
           write(stdo,*) 'Average number of abs/scat events per photon package = ', &
                      ieventcounttot/(1.d0*rt_mcparams%nphot_therm)
           write(stdo,*) 'Average nr of times a photon stays in the same cell  = ', &
                      mc_revisitcell/(1.d0*mc_visitcell)
           write(stdo,*) 'Maximum nr of times a photon stayed in the same cell = ', &
                      int(mc_revisitcell_max)
           if(.not.rt_mcparams%mod_random_walk) then
              if(mc_revisitcell/(1.d0*mc_visitcell).gt.400) then
                 write(stdo,*) '   ---> This number is very high, and therefore responsible for slow performance.'
                 write(stdo,*) '        You may want to consider using the Modified Random Walk (MRW) method.'
                 write(stdo,*) '        To do this, add the line modified_random_walk = 1 to the radmc3d.inp file.'
                 write(stdo,*) '        If you wish, you can fine-tune the MRW behavior (see manual).'
              endif
           endif
           !
           !
           ifile = 0
           write(stdo,*) 'Saving results of thermal Monte Carlo'
           call write_results_mc_bjorkmanwood(ifile)
           !
           ! Dump the photon statistics
           !
           if(rt_mcparams%debug_write_stats.ne.0) then
              write(stdo,*) 'Writing photon statistics to file (because debug_write_stats=1)'
              call write_photon_statistics_to_file()
           endif
           !
           ! Do diffusion
            if (rt_mcparams%nphotdiff.gt.0) then
               call smooth_by_diffusion
            endif
           ! Now check if we want to do a vertical structure computation
           !
           ! ############# FOR NOW THIS IS NOT YET IMPLEMENTED ###############
           !         
           if(iterstr.lt.rt_mcparams%niter_vstruct) then 
              write(stdo,*) 'ERROR: Vertical structure iteration for now disabled.'
              stop
              !
              ! Call the vertical structure integrator, of course only
              ! for the species that are flagged to be done. So you
              ! can still have an envelope and not flag it for vertical
              ! structure iteration.
              !
              !!!!!! call vertstruct_integrate(tdust,ivstrt,error,iwarn)
              !
              ! Check for warnings
              !
              !!!!!! if(iwarn.ne.0) then
              !!!!!!    write(stdo,*) 'WARNING: Grid near midplane too coarse'
              !!!!!!    write(stdo,*) '  for proper calculation of vert struct.'
              !!!!!! endif
              !
              ! Write the dust density
              !
              !!!!!! call write_dustdens()
              !
              ! Check for error
              !
              !!!!!! write(stdo,*) 'Error in vertical structure = ',error
              !!!!!! if(error.lt.errtol) then
              !!!!!!    write(stdo,*) ' '
              !!!!!!    write(stdo,*) ' '
              !!!!!!    write(stdo,*) 'Convergence in vertical structure '
              !!!!!!    write(stdo,*) '   (as measured in rho at the equator)'
              !!!!!!    write(stdo,*) '   (Error=',error,', errtol=',errtol,')'
              !!!!!!    write(stdo,*) 'DONE...'
              !!!!!!    return
              !!!!!! endif
           endif
           !
           ! Check what is the error in the vertical structure
           !
           if(rt_mcparams%niter_vstruct.gt.0) then
              write(stdo,*) ' '
              write(stdo,*) '================================================================'
              write(stdo,*) 'Could not get error in vertical structure below ',rt_mcparams%vserrtol
              write(stdo,*) 'Final error in vertical structure = ',error
              write(stdo,*) '================================================================'
!           else
              !write(stdo,*) 'DONE...'
           endif
           !
        enddo
        !
        ! Restore line mode if necessary
        !
        rt_incl_lines  = inclines_bk
     endif
     !
     !----------------------------------------------------------------
     !          DO THE SINGLE-FREQ SCATTERING MONTE CARLO
     !----------------------------------------------------------------
     !
     ! NOTE: This is normally automatically done within the imaging or
     !       spectrum generating. Here we may wish to do this separately
     !       simply to be able to get mean intensities at certain
     !       wavelengths for use in other kinds models
     !       (e.g. photodissociation of molecules or so). 
     !
     if(do_montecarlo_mono) then
        !
        ! A message:
        !
        call write_message_rad_processes()
        !
        ! If the dust emission is included, then make sure the dust data,
        ! density and temperature are read. If yes, do not read again.
        !
        if(rt_incl_dust) then
           call read_dustdata(1)
           call read_dust_density(1)
           call read_dust_temperature(1)
        endif
        !
        ! If line emission is included, then make sure the line data are
        ! read. If yes, then do not read it again.
        !
        if(rt_incl_lines) then
           call read_lines_all(1)
        endif
        !
        ! If gas continuum is included, then make sure the gas continuum
        ! data are read. If yes, then do not read it again.
        !
        if(rt_incl_gascont) then
           call gascont_init(1)
        endif
        !
        ! Set the mc_frequencies(:) array
        !
        call read_mc_frequencies()
        !
        ! Now call the monochromatic Monte Carlo
        !
        call do_monte_carlo_scattering(rt_mcparams,ierror,   &
               resetseed=do_resetseed,meanint=.true.)
        !
        ! Write the mean intensities to a file
        !
        write(stdo,*) 'Writing mean intensity file...'
        call write_meanint_to_file()
     endif
     !
     !----------------------------------------------------------------
     !     DO THE USER-DEFINED ACTION (CALLED BY radmc3d myaction)
     !----------------------------------------------------------------
     !
     if(do_userdef_action) then
        call userdef_action()
     endif
     !
     !----------------------------------------------------------------
     !                      MAKE A SPECTRUM 
     !----------------------------------------------------------------
     !
     if(do_raytrace_spectrum) then
        !
        ! Check
        !
        if((.not.allocated(freq_nu)).or.(freq_nr.le.0)) then
           write(stdo,*) 'ERROR: Cannot make SED or spectrum if global frequency array not set.'
           stop
        endif
        !
        ! Message
        !
        write(stdo,*) 'Starting procedure for rendering images for the spectrum...'
        call write_message_rad_processes()
        !
        ! Set the camera_frequencies(:) array
        !
        call set_camera_frequencies()
        !
        ! Back up the original value of the nphot_scat and replace it with
        ! the "quick-and-dirty" value nphot_spec. The idea is that we can afford
        ! much lower nphot for the scattering monte carlo if we are not interested
        ! in the images, but only in their flux. 
        !
        nphot_scat_original    = rt_mcparams%nphot_scat
        rt_mcparams%nphot_scat = rt_mcparams%nphot_spec
        !
        ! Set warning flag to 0
        !
        lines_maser_warning = .false.
        !
        ! Make the spectrum. The camera frequency grid should already
        ! have been set in the call to set_camera_frequencies() above.
        !
        if(.not.circular_images) then
           !
           ! Normal way of making a spectrum: make rectangular images
           ! and integrate over them
           !
           call camera_make_spectrum()
        else
           !
           ! Special way of making a spectrum (only in spherical coordinates):
           ! make 'circular images' and integrate over them.
           !
           if((igrid_coord.lt.100).or.(igrid_coord.ge.200)) then
              write(stdo,*) 'ERROR: circular images only possible for spherical coordinates'
              stop
           endif
           if((amr_ydim.eq.0).and.(amr_zdim.eq.0)) then
              !
              ! For 1D spherical models: no need for phi-dispered circular images
              !
              camera_circ_nrphiinf = 1
           endif
           !
           ! Call the circular spectrum solver
           !
           call camera_make_circ_spectrum()
        endif
        !
        ! Messages
        !
        if(camera_incl_stars.eq.0) then
           write(stdo,*) '******************************************************'
           write(stdo,*) 'WARNING: Spectrum made without the light of the stars!'
           write(stdo,*) '         This is OK, as long as you are aware of this!'
           write(stdo,*) '******************************************************'
        endif
        if(lines_maser_warning) then
           write(stdo,*) 'WARNING: Detected masering in spectrum integration.'
           write(stdo,*) '         Put opacity to 0 where this occurred.'
        endif
        write(stdo,*) 'Done rendering spectrum...'
        call flush(stdo)
        !
        ! Restore the original nphot_scat value
        !
        rt_mcparams%nphot_scat = nphot_scat_original
        !
        ! Write the spectrum to file
        !
        if(.not.radmc_as_child) then
           write(stdo,*) 'Writing spectrum to file...'
           call flush(stdo)
           call camera_write_spectrum()
        endif
        !
        ! Make statement about the scattering mode used
        !
        if(rt_incl_dust) then
           if(scattering_mode.eq.0) then
              write(stdo,*) 'Used scattering_mode=0, meaning: no scattering.'
           elseif(scattering_mode.eq.1) then
              write(stdo,*) 'Used scattering_mode=1, meaning: isotropic scattering approximation.'
           elseif(scattering_mode.eq.2) then
              write(stdo,*) 'Used scattering_mode=2, meaning: Henyey-Greenstein scattering approximation.'
           elseif(scattering_mode.eq.3) then
              write(stdo,*) 'Used scattering_mode=3, meaning: non-polarized scattering, ' &
                   //'but with tabulated phase function.'
           elseif(scattering_mode.eq.4) then
              write(stdo,*) 'Used scattering_mode=4, meaning: non-polarized scattering, ' &
                   //'but with tabulated phase function in the Monte Carlo,'
              write(stdo,*) '      but full polarization upon last scattering.'
           elseif(scattering_mode.eq.5) then
              write(stdo,*) 'Used scattering_mode=5, meaning: fully polarized scattering.'
           else
              write(stdo,*) 'I do not know scattering_mode=',scattering_mode
              stop
           endif
           if(alignment_mode.eq.0) then
           elseif(alignment_mode.eq.-1) then
              write(stdo,*) 'Used alignment_mode=-1, meaning: polarized thermal emission from ', &
                   'aligned grains (only for images/spectra).'
           elseif(alignment_mode.eq.1) then
              write(stdo,*) 'Used alignment_mode=1, meaning: polarized thermal emission from ', &
                   'aligned grains (for images/spectra and scattering monte carlo).'
           else
              write(stdo,*) 'I do not know alignment_mode=',alignment_mode
              stop
           endif
        endif
     endif
     !
     !----------------------------------------------------------------
     !                     MAKE IMAGE
     ! NOTE: If radmc3d is run as a child of a parent process, the mere
     !       rendering of the image does not produce any output. You
     !       must request the image by entering 'writeimage' in the
     !       standard input.
     !----------------------------------------------------------------
     !
     if(do_raytrace_image) then
        write(stdo,*) 'Starting procedure for rendering image...'
        call write_message_rad_processes()
        !
        ! Set the camera_frequencies(:) array, if requested
        !
        call set_camera_frequencies()
        !
        ! Set warning flag to 0
        !
        lines_maser_warning = .false.
        !
        ! Call the image maker
        !
        if(.not.circular_images) then
           !
           ! Standard rectangular images
           !
           call camera_make_rect_image(0)
        else
           !
           ! Special case: spherically arranged pixels 
           ! (only for spherical coordinates)
           !
           if((igrid_coord.lt.100).or.(igrid_coord.ge.200)) then
              write(stdo,*) 'ERROR: circular images only possible for spherical coordinates'
              stop
           endif
           if((amr_ydim.eq.0).and.(amr_zdim.eq.0)) then
              !
              ! For 1D spherical models: no need for phi-dispered circular images
              !
              camera_circ_nrphiinf = 1
           endif
           !
           ! Call the circular image routine
           !
           call camera_make_circ_image()
        endif
        !
        ! Check for warnings
        !
        if(lines_maser_warning) then
           write(stdo,*) 'WARNING: Detected masering in image integration.'
           write(stdo,*) '         Put opacity to 0 where this occurred.'
        endif
        !
        ! Write output
        !
        if(.not.radmc_as_child) then
           write(stdo,*) 'Writing image to file...'
           call flush(stdo)
           if(.not.circular_images) then
              call camera_write_image(0,0)
           else
              call camera_write_circ_image(.false.)
           endif
        endif
        !
        ! Make statement about the scattering mode used
        !
        if(rt_incl_dust) then
           if(scattering_mode.eq.0) then
              write(stdo,*) 'Used scattering_mode=0, meaning: no scattering.'
           elseif(scattering_mode.eq.1) then
              write(stdo,*) 'Used scattering_mode=1, meaning: isotropic scattering approximation.'
           elseif(scattering_mode.eq.2) then
              write(stdo,*) 'Used scattering_mode=2, meaning: Henyey-Greenstein scattering approximation.'
           elseif(scattering_mode.eq.3) then
              write(stdo,*) 'Used scattering_mode=3, meaning: non-polarized scattering, ' &
                   //'but with tabulated phase function.'
           elseif(scattering_mode.eq.4) then
              write(stdo,*) 'Used scattering_mode=4, meaning: non-polarized scattering, ' &
                   //'but with tabulated phase function in the Monte Carlo,'
              write(stdo,*) '      but full polarization upon last scattering.'
           elseif(scattering_mode.eq.5) then
              write(stdo,*) 'Used scattering_mode=5, meaning: fully polarized scattering.'
           else
              write(stdo,*) 'I do not know scattering_mode=',scattering_mode
              stop
           endif
           if(alignment_mode.eq.0) then
           elseif(alignment_mode.eq.-1) then
              write(stdo,*) 'Used alignment_mode=-1, meaning: polarized thermal emission from ', &
                   'aligned grains (only for images/spectra).'
           elseif(alignment_mode.eq.1) then
              write(stdo,*) 'Used alignment_mode=1, meaning: polarized thermal emission from ', &
                   'aligned grains (for images/spectra and scattering monte carlo).'
           else
              write(stdo,*) 'I do not know alignment_mode=',alignment_mode
              stop
           endif
        endif
        !
        ! Make a warning if the Monte Carlo is done with full Stokes
        ! vector while the image is scalar.
        !
        if((scattering_mode.gt.3).and.(.not.camera_stokesvector)) then
           write(stdo,*) 'Warning: The Monte Carlo is done with full Stokes components,'
           write(stdo,*) '         while the image is not. That is ok, but if you wish'
           write(stdo,*) '         to have a polarized image too, then add the keyword'
           write(stdo,*) '         stokes to the command line call of radmc3d.'
        endif
        !
        ! If sub-pixeling done, then include diagnostics
        !
        if((camera_nrrefine.gt.0).and.(.not.circular_images)) then
           write(stdo,*) 'Diagnostics of flux-conservation (sub-pixeling):'
           write(stdo,*) '    Nr of main pixels (nx*ny)   = ',camera_image_nx*camera_image_ny
           write(stdo,*) '    Nr of (sub)pixels raytraced = ',camera_subpixeling_npixtot
           write(stdo,*) '    Nr of (sub)pixels used      = ',camera_subpixeling_npixfine
           incrraycost = camera_subpixeling_npixtot/(camera_image_nx*camera_image_ny*1.d0)
           write(stdo,*) '    Increase of raytracing cost = ',incrraycost
           if((incrraycost.gt.20).and.(camera_image_nx*camera_image_ny.gt.1000)) then
              write(stdo,*) '    ----> This increase seems pretty large, and may cause'
              write(stdo,*) '          the ray-tracing to become very slow. It may be'
              write(stdo,*) '          worthwhile to figure out how/where/why the sub-pixeling'
              write(stdo,*) '          is so strong. Maybe you can do subpixeling diagnostics'
              write(stdo,*) '          by adding the keyword diag_subpix to the command line.'
              write(stdo,*) '          This dumps a file called subpixeling_diagnostics.out'
              write(stdo,*) '          which has, for each pixel or sub-pixel, a line containing:'
              write(stdo,*) '            x  y  dx  dy'
              write(stdo,*) '          where x,y are the image-plane coordinates of the pixel, and'
              write(stdo,*) '          dx,dy are the image-plane pixel sizes.'
              if((igrid_coord.ge.100).and.(igrid_coord.lt.200)) then
                 write(stdo,*) ' '
                 write(stdo,*) '          ** ANOTHER TIP for spherical coordinates: **'
                 write(stdo,*) '          Maybe this happens because of refinement in R or Theta coordinates.'
                 write(stdo,*) '          RADMC-3D will always try to make sub-pixeling such that it resolves'
                 write(stdo,*) '          all cells as seen in projection. If you have extreme R-refinement this'
                 write(stdo,*) '          will cause extreme sub-pixeling. Also if you use a linearly spaced grid'
                 write(stdo,*) '          in R but you have a large span in R (max(R)/min(R)>>1), then the outer'
                 write(stdo,*) '          R-grid will be very over-resolved (or the inner R-grid under-resolved).'
                 write(stdo,*) '          It is then advisable to use a log-spaced grid.'
                 write(stdo,*) ' '
                 write(stdo,*) '          You can moderate excessive sub-pixeling caused by fine R-grids by '
                 write(stdo,*) '          setting the camera_min_drr parameter in radmc3d.inp to a value larger'
                 write(stdo,*) '          than its default. Similar with camera_min_dangle for fine theta-grids.'
                 write(stdo,*) ' '
                 write(stdo,*) '          If you use (instead) the command-line option "sloppy" then RADMC-3D'
                 write(stdo,*) '          will set the camera_min_drr=0.1 and camera_min_dangle=0.1 as well as'
                 write(stdo,*) '          camera_spher_cavity_relres=0.1, which are rather large values, '
                 write(stdo,*) '          but hopefully still small enough.'
                 write(stdo,*) '          '
                 write(stdo,*) '          All of these tricks to speed things up are, however, AT YOUR OWN RISK.'
                 write(stdo,*) '          So you can use it, but it is wise to occasionally check without tricks.'
              endif
           endif
        endif
     endif
     !
     !----------------------------------------------------------------
     !                FIND THE TAU=TAUSURF SURFACE
     ! NOTE: If radmc3d is run as a child of a parent process, the mere
     !       rendering of the image does not produce any output. You
     !       must request the image by entering 'writeimage' in the
     !       standard input.
     !----------------------------------------------------------------
     !
     if(do_raytrace_tausurf) then
        write(stdo,*) 'Starting procedure for rendering tau-surface...'
        call write_message_rad_processes()
        !
        ! Set the camera_frequencies(:) array, if requested
        !
        call set_camera_frequencies()
        !
        ! Call the image maker
        !
        call camera_make_rect_image(0,tausurf=.true.)
        !
        ! Write projected tau surface
        !
        if(.not.radmc_as_child) then
           write(stdo,*) 'Writing projected tau surface to file...'
           call flush(stdo)
           call camera_write_image(0,0,noclip=.true.)
        endif
        !
        ! Write the 3-D tau surface (ASCII format)
        !
        if(allocated(camera_tausurface_x)) then
           write(stdo,*) 'Writing 3-D tau surface to file...'
           open(unit=1,file='tausurface_3d.out')
           write(1,*) 1  ! Format number
           write(1,*) camera_image_nx,camera_image_ny
           write(1,*) camera_nrfreq
           write(1,'(E21.14)') (1d4*cc/camera_frequencies(iinu),iinu=1,camera_nrfreq)
           write(1,*) ' '
           do iinu=1,camera_nrfreq
              do imy=1,camera_image_ny
                 do imx=1,camera_image_nx
                    write(1,'(3(E13.6,1X))')                &
                         camera_tausurface_x(imx,imy,iinu), &
                         camera_tausurface_y(imx,imy,iinu), &
                         camera_tausurface_z(imx,imy,iinu)
                 enddo
              enddo
           enddo
           close(1)
        endif
     endif
     !
     !----------------------------------------------------------------
     !                       MAKE MOVIE
     !----------------------------------------------------------------
     !
     if(do_raytrace_movie) then
        !
        ! First open the movie.inp file
        !
        write(stdo,*) 'Reading movie.inp file...'
        call flush(stdo)
        inquire(file='movie.inp',exist=fex1)
        if(.not.fex1) then
           write(stdo,*) 'ERROR: If you want to make a movie sequence of images,'
           write(stdo,*) '       you must prepare a movie.inp file.'
           stop
        endif
        open(unit=1,file='movie.inp')
        read(1,*) iformat
        read(1,*) cameras_nr_images
        if(cameras_nr_images.lt.10000) then
           ifill = 4
        elseif(cameras_nr_images.lt.100000) then
           ifill = 5
        else
           write(stdo,*) 'ERROR: Too many movie frames. '
           stop
        endif
        !
        ! If iformat >0 then the observer is at infinity, if iformat<0 then
        ! the observer is local
        !
        if(iformat.gt.0) then
           camera_localobserver = .false.
        else
           camera_localobserver = .true.
        endif
        !
        ! For safetly, deallocate the cameras arrays
        !
        if(allocated(cameras_pt_pos)) deallocate(cameras_pt_pos)
        if(allocated(cameras_img_hs_x)) deallocate(cameras_img_hs_x)
        if(allocated(cameras_img_hs_y)) deallocate(cameras_img_hs_y)
        if(allocated(cameras_pt_degr_pa)) deallocate(cameras_pt_degr_pa)
        if(allocated(cameras_zmc_x)) deallocate(cameras_zmc_x)
        if(allocated(cameras_zmc_y)) deallocate(cameras_zmc_y)
        if(allocated(cameras_obs_pos)) deallocate(cameras_obs_pos)
        if(allocated(cameras_obs_degr_th)) deallocate(cameras_obs_degr_th)
        if(allocated(cameras_obs_degr_ph)) deallocate(cameras_obs_degr_ph)
        !
        ! Allocate the necessary arrays
        !
        allocate(cameras_pt_pos(1:3,1:cameras_nr_images),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate cameras_pt_pos array'
           stop
        endif
        allocate(cameras_img_hs_x(1:cameras_nr_images),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate cameras_img_hs_x array'
           stop
        endif
        allocate(cameras_img_hs_y(1:cameras_nr_images),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate cameras_img_hs_y array'
           stop
        endif
        allocate(cameras_pt_degr_pa(1:cameras_nr_images),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate cameras_pt_degr_pa array'
           stop
        endif
        allocate(cameras_zmc_x(1:cameras_nr_images),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate cameras_zmc_x array'
           stop
        endif
        allocate(cameras_zmc_y(1:cameras_nr_images),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate cameras_zmc_y array'
           stop
        endif
        if(camera_localobserver) then
           allocate(cameras_obs_pos(1:3,1:cameras_nr_images),STAT=ierr)
           if(ierr.ne.0) then
              write(stdo,*) 'ERROR: Could not allocate cameras_obs_pos array'
              stop
           endif
        else
           allocate(cameras_obs_degr_th(1:cameras_nr_images),STAT=ierr)
           if(ierr.ne.0) then
              write(stdo,*) 'ERROR: Could not allocate cameras_obs_degr_th array'
              stop
           endif
           allocate(cameras_obs_degr_ph(1:cameras_nr_images),STAT=ierr)
           if(ierr.ne.0) then
              write(stdo,*) 'ERROR: Could not allocate cameras_obs_degr_ph array'
              stop
           endif
        endif
        !
        ! Now do a loop over all the images
        !
        do i=1,cameras_nr_images
           if(iformat.eq.1) then 
              !
              ! A sequence of images with the observer at infinity
              !
              read(1,*) cameras_pt_pos(1:3,i),cameras_img_hs_x(i),cameras_img_hs_y(i),&
                        cameras_pt_degr_pa(i),cameras_obs_degr_th(i),cameras_obs_degr_ph(i)
              cameras_zmc_x(i) = 0.d0
              cameras_zmc_y(i) = 0.d0
           elseif(iformat.eq.2) then
              !
              ! A sequence of images with the observer at infinity, but with
              ! more control options
              !
              write(stdo,*) 'ERROR: Format 2 for the movie.inp is not yet ready.'
              write(stdo,*) '       It should include the zoom geometry information.'
              stop
!                 read(1,*) cameras_pt_pos(1:3,i),cameras_img_hs_x(i),cameras_img_hs_y(i),&
!                           cameras_pt_degr_pa(i),cameras_obs_degr_th(i),cameras_obs_degr_ph(i),&
!                           cameras_zmc_x(i),cameras_zmc_y(i)
!                      NOTE: Instead of cameras_img_hs_x, and cameras_zmc_x I
!                            thing we should read the zoom geometry like in the
!                            case of a normal image.
           elseif(iformat.eq.-1) then
              !
              ! A sequence of images with a local observer
              !
              read(1,*) cameras_pt_pos(1:3,i),cameras_img_hs_x(i),cameras_img_hs_y(i),&
                        cameras_pt_degr_pa(i),cameras_obs_pos(1:3,i)
              cameras_zmc_x(i) = 0.d0
              cameras_zmc_y(i) = 0.d0
           endif
        enddo
        close(1)
        write(stdo,*) 'Rendering movie...'
        call flush(stdo)
        !
        ! A message:
        !
        call write_message_rad_processes()
        !
        ! Set the camera_frequencies(:) array, if requested
        !
        call set_camera_frequencies()
        !
        ! Set warning flag to 0
        !
        lines_maser_warning = .false.
        !
        ! Now make the movie
        !
        open(unit=1,file='movie.info')
        write(1,*) -1
        close(1)
        do i=1,cameras_nr_images
           write(stdo,*) 'Movie frame ',i
           call flush(stdo)
           call camera_make_rect_image(i)
           call camera_write_image(i,ifill)
           open(unit=1,file='movie.info')
           write(1,*) i
           close(1)
        enddo
        !
        ! Warning
        !
        if(lines_maser_warning) then
           write(stdo,*) 'WARNING: Detected masering in movie frame integration.'
           write(stdo,*) '         Put opacity to 0 where this occurred.'
        endif
     endif
     !
     !----------------------------------------------------------------
     !                 CALCULATE LEVEL POPULATIONS
     ! Note: For LTE and *local* non-LTE modes this is not necessary, 
     !       because they are always automatically calculated 
     !       beforehand. But even for those modes it can come in 
     !       handy if you simply want to check out the populations
     !       without making an image or spectrum. 
     !----------------------------------------------------------------
     !
     if(do_calcpop) then
        !
        ! Check
        !
        if(lines_mode.lt.0) then
           write(stdo,*) 'ERROR: Cannot calculate and store the level'
           write(stdo,*) '       populations and write them to file, if'
           write(stdo,*) '       the levels_mode is negative. These modes'
           write(stdo,*) '       are for on-the-fly computation of the'
           write(stdo,*) '       populations only.'
           stop
        endif
        !
        ! If the dust emission is included, then make sure the dust data,
        ! density and temperature are read. If yes, do not read again.
        !
        if(rt_incl_dust) then
           call read_dustdata(1)
           call read_dust_density(1)
           call read_dust_temperature(1)
        endif
        !
        ! If line emission is included, then make sure the line data are
        ! read. If yes, then do not read it again.
        !
        if(rt_incl_lines) then
           call read_lines_all(1)
        endif
        !
        ! If gas continuum is included, then make sure the gas continuum
        ! data are read. If yes, then do not read it again.
        !
        if(rt_incl_gascont) then
           call gascont_init(1)
        endif
        !
        ! Set warning flag to 0
        !
        lines_maser_warning = .false.
        !
        ! If lines are active, and if the level populations are to be calculated
        ! beforehand and stored in the big array, then compute them now, if not
        ! already done.
        !
        if(rt_incl_lines) then
           if((lines_mode.ge.1).and.(lines_mode.le.9)) then
              call lines_compute_and_store_local_populations(1)
           endif
        endif
        !
        ! Warning
        !
        if(lines_maser_warning) then
           write(stdo,*) 'WARNING: Detected masering in calculating level populations.'
           write(stdo,*) '         Put opacity to 0 where this occurred.'
        endif
     endif
     !
     !----------------------------------------------------------------
     !                 WRITE THE LATEST IMAGE AGAIN 
     !----------------------------------------------------------------
     !
     if(do_writeimage.and.radmc_as_child) then
        !
        ! First check if an image has indeed been made
        !
        if((.not.(allocated(camera_rect_image_iquv))).or.&
           (camera_image_nx.le.0).or.(camera_image_ny.le.0)) then
           !
           ! No image to write, so ignore
           !
           write(stdo,*) 'No image to write... Ignoring this command.'
           call flush(stdo)
        else
           !
           ! Yes, there is an image still in memory. So write.
           !
           write(stdo,*) 'Writing image to standard output...'
           call flush(stdo)
           call camera_write_image(0,0)
        endif
     endif
     !
     !----------------------------------------------------------------
     !                 WRITE THE LATEST SPECTRUM AGAIN 
     !----------------------------------------------------------------
     !
     if(do_writespectrum.and.radmc_as_child) then
        !
        ! First check if a spectrum has indeed been made
        !
        if((.not.(allocated(camera_spectrum_iquv))).or.&
           (camera_image_nx.le.0).or.(camera_image_ny.le.0).or.&
           (camera_nrfreq.le.0)) then
           !
           ! No spectrum to write, so ignore
           !
           write(stdo,*) 'No spectrum to write... Ignoring this command.'
           call flush(stdo)
        else
           !
           ! Yes, there is a spectrum still in memory. So write.
           !
           write(stdo,*) 'Writing spectrum to standard output...'
           call flush(stdo)
           call camera_write_spectrum()
        endif
     endif
     !
     !----------------------------------------------------------------
     !                    WRITE SOME LINE INFO
     !----------------------------------------------------------------
     !
     if(do_write_lineinfo.and.rt_incl_lines) then
        if((lines_user_ispec.gt.0).and.(lines_user_iline.gt.0)) then
           call lines_print_lineinfo(lines_user_iline,lines_user_ispec)
        endif
     endif
     !
     !----------------------------------------------------------------
     !       DOES THE USERDEF MODULE WANT TO DO SOMETHING HERE?
     !----------------------------------------------------------------
     !
     call userdef_dostuff()
     !
     !----------------------------------------------------------------
     !         WRITE THE INTERNALLY PRODUCED MODEL SETUP DATA
     !----------------------------------------------------------------
     !
     if(do_writemodel) then
        call userdef_writemodel()
     endif
     !
     !----------------------------------------------------------------
     !                    WRITE THE GRID FILE ONLY
     !----------------------------------------------------------------
     !
     if(do_writegridfile) then
        call write_grid_file()
     endif
     !
     !----------------------------------------------------------------
     !       WRITE THE INTERNALLY PRODUCED LEVEL POPULATION DATA
     !----------------------------------------------------------------
     !
     if(do_writepop) then
        call write_levelpop()
     endif
     !
     !----------------------------------------------------------------
     !          EXTRACT A REGULARLY GRIDDED SUBBOX FROM MODEL
     !                   Only for user convenience
     !----------------------------------------------------------------
     !
     ! First a default value for the subbox size
     !
     if((subbox_x0.eq.0.d0).and.(subbox_x1.eq.0.d0).and.       &
          (subbox_y0.eq.0.d0).and.(subbox_y1.eq.0.d0).and.     &
          (subbox_z0.eq.0.d0).and.(subbox_z1.eq.0.d0)) then
        if(igrid_type.lt.100) then
           !
           ! Regular or AMR grid
           !
           if(igrid_coord.lt.100) then
              !
              ! Cartesian coordinates
              !
              subbox_x0 = amr_grid_xi(1,1)
              subbox_x1 = amr_grid_xi(amr_grid_nx+1,1)
              subbox_y0 = amr_grid_xi(1,2)
              subbox_y1 = amr_grid_xi(amr_grid_ny+1,2)
              subbox_z0 = amr_grid_xi(1,3)
              subbox_z1 = amr_grid_xi(amr_grid_nz+1,3)
           elseif(igrid_coord.lt.200) then
              !
              ! Spherical coordinates
              !
              subbox_x0 = -1.01d0*amr_grid_xi(amr_grid_nx+1,1)
              subbox_x1 =  1.01d0*amr_grid_xi(amr_grid_nx+1,1)
              subbox_y0 = -1.01d0*amr_grid_xi(amr_grid_nx+1,1)
              subbox_y1 =  1.01d0*amr_grid_xi(amr_grid_nx+1,1)
              subbox_z0 = -1.01d0*amr_grid_xi(amr_grid_nx+1,1)
              subbox_z1 =  1.01d0*amr_grid_xi(amr_grid_nx+1,1)
           else
              stop 6503
           endif
        else
           write(stdo,*) 'ERROR: In computing subbox default size: unstructured grids not yet supported'
           stop
        endif
     endif
     !
     ! Now make the box, if requested
     !
     if(do_write_subbox_dustdens) then
        !
        ! Dust density box
        !
        call read_dustdata(1)       ! Only read if not yet present
        call read_dust_density(1)   ! Only read if not yet present
        if(.not.allocated(dustdens)) then
           write(stdo,*) 'ERROR in writing out subbox: dustdens not allocated.'
           stop
        endif
        subboxfilename = 'dust_density_subbox.out'
        nv = dust_nr_species
        write(stdo,*) 'Computing and writing subbox of dust density...'
        call do_subbox_extract_and_write(nv,nrcells,subboxfilename,        &
             subbox_nx,subbox_ny,subbox_nz,subbox_x0,subbox_x1,subbox_y0,  &
             subbox_y1,subbox_z0,subbox_z1,                                &
             subbox_phi1,subbox_theta,subbox_phi2,dustdens)
     endif
     if(do_write_subbox_dusttemp) then
        !
        ! Dust temperature box
        !
        call read_dustdata(1)           ! Only read if not yet present
        call read_dust_temperature(1)   ! Only read if not yet present
        if(.not.allocated(dusttemp)) then
           write(stdo,*) 'ERROR in writing out subbox: dusttemp not allocated.'
           stop
        endif
        subboxfilename = 'dust_temperature_subbox.out'
        nv = dust_nr_species
        write(stdo,*) 'Computing and writing subbox of dust temperature...'
        call do_subbox_extract_and_write(nv,nrcells,subboxfilename,        &
             subbox_nx,subbox_ny,subbox_nz,subbox_x0,subbox_x1,subbox_y0,  &
             subbox_y1,subbox_z0,subbox_z1,                                &
             subbox_phi1,subbox_theta,subbox_phi2,dusttemp)
     endif
     if(do_write_subbox_levelpop) then
        !
        ! Level population box
        !
        ! Check
        !
        if(lines_mode.lt.0) then
           write(stdo,*) 'ERROR: Cannot calculate and store the level'
           write(stdo,*) '       populations and write them to file, if'
           write(stdo,*) '       the levels_mode is negative. These modes'
           write(stdo,*) '       are for on-the-fly computation of the'
           write(stdo,*) '       populations only.'
           stop
        endif
        !
        ! If the dust emission is included, then make sure the dust data,
        ! density and temperature are read. If yes, do not read again.
        !
        if(rt_incl_dust) then
           call read_dustdata(1)
           call read_dust_density(1)
           call read_dust_temperature(1)
        endif
        !
        ! If line emission is included, then make sure the line data are
        ! read. If yes, then do not read it again.
        !
        if(rt_incl_lines) then
           call read_lines_all(1)
        endif
        !
        ! If gas continuum is included, then make sure the gas continuum
        ! data are read. If yes, then do not read it again.
        !
        if(rt_incl_gascont) then
           call gascont_init(1)
        endif
        !
        ! If lines are active, and if the level populations are to be calculated
        ! beforehand and stored in the big array, then compute them now, if not
        ! already done.
        !
        if((lines_mode.ge.1).and.(lines_mode.le.9)) then
           call lines_compute_and_store_local_populations(1)
        endif
        !
        ! Now make a loop over all molecules
        !
        do ispec=1,lines_nr_species
           if(lines_species_fullmolec(ispec)) then
              !
              ! Create the file name
              !
              species = lines_speciesname(ispec)
              subboxfilename = 'levelpop_'//species(1:len_trim(species))//'_subbox.out'
              !
              ! Now write the box
              !
              nv = lines_nrlevels_subset(ispec)
              nvm = lines_nrlevels_subset_max
              write(stdo,*) 'Computing and writing subbox of levelpopulations of '//species(1:len_trim(species))//'...'
              call do_subbox_levelpop_extract_and_write(nvm,nv,lines_nr_species, &
                   nrcells,ispec,lines_levels_subset(:,ispec),subboxfilename, &
                   subbox_nx,subbox_ny,subbox_nz,subbox_x0,subbox_x1,subbox_y0,  &
                   subbox_y1,subbox_z0,subbox_z1,                                &
                   subbox_phi1,subbox_theta,subbox_phi2,lines_levelpop)
           endif
        enddo
     endif
     !
     ! Do the sampling, if requested
     !
     if(do_writesample) then
        !
        ! First read the array of grid points for the sampling
        !
        call read_sampling_points_3d(npt,xpt,ypt,zpt)
        !
        ! Now do the sampling and write the outputs
        !
        if(do_sample_dustdens) then
           call read_dustdata(1)       ! Only read if not yet present
           call read_dust_density(1)   ! Only read if not yet present
           if(.not.allocated(dustdens)) then
              write(stdo,*) 'ERROR in writing out subbox: dustdens not allocated.'
              stop
           endif
           samplefilename = 'dust_density_sample.out'
           nv = dust_nr_species
           write(stdo,*) 'Computing and writing sample of dust density...'
           call do_sample_extract_and_write(nv,nrcells,samplefilename,     &
                npt,xpt,ypt,zpt,dustdens)
        endif
        if(do_sample_dusttemp) then
           call read_dustdata(1)           ! Only read if not yet present
           call read_dust_temperature(1)   ! Only read if not yet present
           if(.not.allocated(dusttemp)) then
              write(stdo,*) 'ERROR in writing out subbox: dusttemp not allocated.'
              stop
           endif
           samplefilename = 'dust_temperature_sample.out'
           nv = dust_nr_species
           write(stdo,*) 'Computing and writing sample of dust temperature...'
           call do_sample_extract_and_write(nv,nrcells,samplefilename,     &
                npt,xpt,ypt,zpt,dusttemp)
        endif
        if(do_sample_levelpop) then
           if(lines_mode.lt.0) then
              write(stdo,*) 'ERROR: Cannot calculate and store the level'
              write(stdo,*) '       populations and write them to file, if'
              write(stdo,*) '       the levels_mode is negative. These modes'
              write(stdo,*) '       are for on-the-fly computation of the'
              write(stdo,*) '       populations only.'
              stop
           endif
           if(rt_incl_dust) then
              call read_dustdata(1)
              call read_dust_density(1)
              call read_dust_temperature(1)
           endif
           call read_lines_all(1)
           !
           ! If gas continuum is included, then make sure the gas continuum
           ! data are read. If yes, then do not read it again.
           !
           if(rt_incl_gascont) then
              call gascont_init(1)
           endif
           !
           ! If lines are active, and if the level populations are to be calculated
           ! beforehand and stored in the big array, then compute them now, if not
           ! already done.
           !
           if((lines_mode.ge.1).and.(lines_mode.le.9)) then
              call lines_compute_and_store_local_populations(1)
           endif
           !
           ! Now make a loop over all molecules
           !
           do ispec=1,lines_nr_species
              if(lines_species_fullmolec(ispec)) then
                 !
                 ! Create the file name
                 !
                 species = lines_speciesname(ispec)
                 samplefilename = 'levelpop_'//species(1:len_trim(species))//'_sample.out'
                 !
                 ! Now write the box
                 !
                 nv = lines_nrlevels_subset(ispec)
                 nvm = lines_nrlevels_subset_max
                 write(stdo,*) 'Computing and writing sample of levelpopulations of '//species(1:len_trim(species))//'...'
                 call do_sample_levelpop_extract_and_write(nvm,nv,lines_nr_species, &
                      nrcells,ispec,lines_levels_subset(:,ispec),                   &
                      samplefilename,npt,xpt,ypt,zpt,lines_levelpop)
              endif
           enddo
        endif
        !
        ! Destroy the sampling points
        !
        call destroy_sampling_points(xpt,ypt,zpt)
        !
     endif
     !
     ! Do output of the VTK data 
     !
     if(do_write_vtk_grid) then
        !
        ! Check if the coordinates are 3-D
        !
        if((amr_xdim.ne.1).or.(amr_ydim.ne.1).or.(amr_zdim.ne.1)) then
           write(stdo,*) 'ERROR: VTK Output only works in 3-D.'
           stop
        endif
        !
        ! If the dust emission is included, then make sure the dust data,
        ! density and temperature are read. If yes, do not read again.
        !
        if(rt_incl_dust) then
           call read_dustdata(1)
           call read_dust_density(1)
           if(do_write_vtk_dust_temperature) call read_dust_temperature(1)
        endif
        !
        ! If line emission is included, then make sure the line data are
        ! read. If yes, then do not read it again.
        !
        if(rt_incl_lines) then
           call read_lines_all(1)
        endif
        !
        ! If gas continuum is included, then make sure the gas continuum
        ! data are read. If yes, then do not read it again.
        !
        if(rt_incl_gascont) then
           call gascont_init(1)
        endif
        !
        ! If lines are active, and if the level populations are to be calculated
        ! beforehand and stored in the big array, then compute them now, if not
        ! already done.
        !
        if((lines_mode.ge.1).and.(lines_mode.le.9)) then
           call lines_compute_and_store_local_populations(1)
        endif
        !
        ! Now write the VTK data
        !
        vtkfilename = 'model.vtk'
        vtktitle    = 'RADMC-3D Model'
        write(stdo,*) 'Opening and writing grid to VTK file ',trim(vtkfilename)
        if (rto_style .eq. 3) then
          call amr_open_vtk_file_b(1,vtkfilename,vtktitle,.true., igrid_coord)
        else
          call amr_open_vtk_file(1,vtkfilename,vtktitle,.true., igrid_coord)
        endif
        if(do_write_vtk_dust_density) then
           if(allocated(dustdens)) then
              write(stdo,*) 'Writing dust density to VTK file'
              ! Attila Juhasz - pass the vtkFieldName to the VTK writer (it will
              ! appear in the scalar field name)
              vtkFieldName = 'dust_density'
              if (rto_style .eq. 3) then
                call amr_write_vtk_file_scalar_alt_b(1,dust_nr_species,amr_nrleafs,&
                   vtk_dust_ispec,dustdens,.true.,vtkFieldName)
              else
                call amr_write_vtk_file_scalar_alt(1,dust_nr_species,amr_nrleafs,&
                   vtk_dust_ispec,dustdens,.true.,vtkFieldName)
              endif
           else
              write(stdo,*) 'WARNING: Since dust density is not allocated,',& 
                   'no VTK data for dust density was written...'
           endif
        endif
        if(do_write_vtk_dust_temperature) then
           if(allocated(dusttemp)) then
              write(stdo,*) 'Writing dust temperature to VTK file'
              ! Attila Juhasz - pass the vtkFieldName to the VTK writer (it will
              ! appear in the scalar field name)
              vtkFieldName = 'dust_temperature'
              if (rto_style .eq. 3) then
                call amr_write_vtk_file_scalar_alt_b(1,dust_nr_species,amr_nrleafs,&
                   vtk_dust_ispec,dusttemp,.false.,vtkFieldName)
              else
                call amr_write_vtk_file_scalar_alt(1,dust_nr_species,amr_nrleafs,&
                   vtk_dust_ispec,dusttemp,.false.,vtkFieldName)
              endif
           else
              write(stdo,*) 'WARNING: Since dust temperature is not allocated,',&
                   ' no VTK data for dust temperature was written...'
           endif
        endif
        if(do_write_vtk_gas_density) then
           if(allocated(gasdens)) then
              write(stdo,*) 'Writing gas density to VTK file'
              ! Attila Juhasz - pass the vtkFieldName to the VTK writer (it will
              ! appear in the scalar field name)
              vtkFieldName = 'gas_density'
              if (rto_style .eq. 3) then
                call amr_write_vtk_file_scalar_b(1,amr_nrleafs,gasdens,.true.,vtkFieldName)
              else
                call amr_write_vtk_file_scalar(1,amr_nrleafs,gasdens,.true.,vtkFieldName)
              endif
           else
              write(stdo,*) 'WARNING: Since gas density is not allocated, ',&
                   'no VTK data for gas density was written...'
           endif
        endif
        if(do_write_vtk_gas_temperature) then
           if(allocated(gastemp)) then
              write(stdo,*) 'Writing gas temperature to VTK file'
              ! Attila Juhasz - pass the vtkFieldName to the VTK writer (it will
              ! appear in the scalar field name)
              vtkFieldName = 'gas_temperature'
              if (rto_style .eq. 3) then
                call amr_write_vtk_file_scalar_b(1,amr_nrleafs,gastemp,.false.,vtkFieldName)
              else
                call amr_write_vtk_file_scalar(1,amr_nrleafs,gastemp,.false.,vtkFieldName)
              endif
           else
              write(stdo,*) 'WARNING: Since gas temperature is not allocated, ',&
                   'no VTK data for gas temperature was written...'
           endif
        endif
        if(do_write_vtk_chemspec) then
           if(allocated(gas_chemspec_numberdens)) then
              write(stdo,*) 'Writing molecule/atom species density to VTK file'
              ! Attila Juhasz - pass the vtkFieldName to the VTK writer (it will
              ! appear in the scalar field name)
              vtkFieldName = 'molecule_density'
              if (rto_style .eq. 3) then
                call amr_write_vtk_file_scalar_alt_b(1,lines_nr_species,amr_nrleafs, &
                   vtk_lines_ispec,gas_chemspec_numberdens,.true.,vtkFieldName)
              else
                call amr_write_vtk_file_scalar_alt(1,lines_nr_species,amr_nrleafs, &
                   vtk_lines_ispec,gas_chemspec_numberdens,.true.,vtkFieldName)
              endif
           else
              write(stdo,*) 'WARNING: Since gas_chemspec_numberdens is not allocated, ',&
                   'no VTK data for molecule/atom density was written...'
           endif
        endif
        if(do_write_vtk_levelpop) then
           if(allocated(lines_levelpop)) then
              write(stdo,*) 'Writing level populations to VTK file'
              ! Attila Juhasz - pass the vtkFieldName to the VTK writer (it will
              ! appear in the scalar field name)
              vtkFieldName = 'level_population'
              if (rto_style .eq. 3) then
                call amr_write_vtk_file_scalar_alt2_b(1,lines_nrlevels_subset_max, &
                   lines_nr_species,amr_nrleafs,vtk_lines_ispec,vtk_lines_ilevel, &
                   lines_levelpop,.true.,vtkFieldName)
              else
                call amr_write_vtk_file_scalar_alt2(1,lines_nrlevels_subset_max, &
                   lines_nr_species,amr_nrleafs,vtk_lines_ispec,vtk_lines_ilevel, &
                   lines_levelpop,.true.,vtkFieldName)
              endif
           else
              write(stdo,*) 'WARNING: Since level population array is not allocated, ',&
                   'no VTK data for level populations was written...'
           endif
        endif
        if(do_write_vtk_velocity) then
           if(allocated(gasvelocity)) then
              write(stdo,*) 'Writing gas velocity vector field to VTK file'
              ! Attila Juhasz - pass the vtkFieldName to the VTK writer (it will
              ! appear in the vector field name) 
              vtkFieldName = 'gas_velocity'
              if (rto_style .eq. 3) then
                call amr_write_vtk_file_vector_b(1,amr_nrleafs,gasvelocity,igrid_coord,vtkFieldName)
              else
                call amr_write_vtk_file_vector(1,amr_nrleafs,gasvelocity,igrid_coord,vtkFieldName)
              endif
           else
              write(stdo,*) 'WARNING: Since gas velocity is not allocated, ',&
                   'no VTK data for gas velocity was written...'
           endif
        endif
        !
        ! Close the file
        !
        close(1)
     endif
     !
     !----------------------------------------------------------------
     !                   GIVE A SMALL SIGNAL IF REQUESTED
     !----------------------------------------------------------------
     !
     if(do_respondwhenready.and.radmc_as_child) then
        write(stdo,*) 'Writing a response to the parent process that we are ready'
        call flush(stdo)        
        write(6,*) 1
        call flush(6)
     endif
     !
     ! A check
     !
     if(lines_pfunc_outofrange) then
        write(stdo,*) '********** WARNING: **********'
        write(stdo,*) 'During line transfer: temperature of gas dropped below table of partition function.'
        lines_pfunc_outofrange = .false.
     endif
     !
     ! ===================================
     !     Now we've done an action. 
     !   We must reset for next action
     ! ===================================
     !
     if(.not.radmc_as_child) then
        !
        ! If we are not in child process mode, then
        ! we should be done now.
        !
        quit = .true.
     else
        !
        ! We are a child process, so let us reset all actions 
        !
        do_montecarlo_therm        = .false.
        do_montecarlo_mono         = .false.
        do_userdef_action          = .false.
        do_vstruct                 = .false.
        do_raytrace_spectrum       = .false.
        do_raytrace_image          = .false.
        do_raytrace_movie          = .false.
        do_raytrace_tausurf        = .false.
        do_writeimage              = .false.
        do_writespectrum           = .false.
        do_write_lineinfo          = .false.
        do_writemodel              = .false.
        do_writegridfile           = .false.
        do_respondwhenready        = .false.
        do_writepop                = .false.
        do_calcpop                 = .false.
        !
        ! Reset other defaults
        !
        camera_theinu              = 0
        camera_lambdamic           = 0.d0
        camera_lambdamic1          = -1.d0
        camera_range_nlam          = -1
        camera_loadcolor           = .false.
        camera_loadlambdas         = .false.
        camera_setfreq_global      = .false.
        camera_secondorder         = .false.
        camera_catch_doppler_line  = .false.
        camera_diagnostics_subpix  = .false.
        lines_user_widthkms        = 0.d0
        lines_user_kms0            = 0.d0
        lines_user_ispec           = 0
        lines_user_iline           = 0
        lines_user_nrfreq          = 1
        do_writescatsrc            = .false.
        lines_make_linelist        = .false.
        grid_was_written_to_file   = .false.
        scattering_mode            = scattering_mode_def
        !
        ! Reset the subbox commands, but not the settings
        !
        do_write_subbox_dustdens   = .false.
        do_write_subbox_dusttemp   = .false.
        do_write_subbox_levelpop   = .false.
        !
        ! Reset sampling commands
        !
        do_writesample             = .false.
        do_sample_dustdens         = .false.
        do_sample_dusttemp         = .false.
        do_sample_levelpop         = .false.
        !
        ! Reset VTK output
        !
        do_write_vtk_grid             = .false.
        do_write_vtk_dust_density     = .false.
        do_write_vtk_dust_temperature = .false.
        do_write_vtk_gas_density      = .false.
        do_write_vtk_gas_temperature  = .false.
        do_write_vtk_chemspec         = .false.
        do_write_vtk_levelpop         = .false.
        do_write_vtk_velocity         = .false.
        !
        ! Reset perhaps some userdef stuff?
        !
        call userdef_reset_flags()
        !
     endif
     !
     ! Check if not exceeding max nr of tasks, to avoid an infinite hang
     ! in case of some strange input
     !
     itask = itask + 1
     if(itask.gt.maxtasks) then
        quit=.true.
     endif
     !
  enddo
  !
  ! Now clean up everything
  !
  call montecarlo_cleanup()
  call dust_cleanup()
  call lines_cleanup()
  call stars_cleanup()
  call rtglobal_cleanup()
  call quantum_cleanup()
  call amrray_cleanup()
  call amr_cleanup()
  call camera_cleanup()
  call tabulated_exp_cleanup()
  nrcells    = 0
  nrcellsmax = 0
  !
  ! Done
  ! 
  write(stdo,*) 'Done...'
  call flush(stdo)
  !
  ! Close the standard out, if necessary
  !
  if(radmc_as_child) then
     close(stdo)
  endif
  !
700 continue
  !
end program radmc3d



!-------------------------------------------------------------------------
!                      READ THE RADMC3D.INP FILE
!-------------------------------------------------------------------------
subroutine read_radmcinp_file()
  use montecarlo_module
  use camera_module
  use lines_module
  use namelist_module
  use stars_module
  use userdef_module
  implicit none
  character*160 :: string
  logical :: fex1,fex2
  integer :: idum,len
  logical :: ierror
  integer :: imethod     ! (only for backward compatibility with RADMC; not used)
  integer :: ioptmot,wghtphot
  integer :: inc_dust,inc_lines,inc_freefree,inc_star_size,interpoljnu
  integer :: inc_usersrc,iautosubset,bc_idir,bc_ilr
  double precision :: mindang,mindrr
  !
  ! Defaults
  !
  ioptmot      = -100
  inc_dust     = -100
  inc_lines    = -100
  inc_freefree = -100
  inc_usersrc  = -100
  inc_star_size= -100
  wghtphot     = -100
  interpoljnu  = -100
  iautosubset  = -100
  mindang      = -100.d0
  mindrr       = -100.d0
  !
  ! Now read the radmc3d.inp file
  !
  inquire(file='radmc3d.inp',exist=fex1)
  inquire(file='radmc.inp',exist=fex2)
  if(fex1.and.fex2) then
     write(stdo,*) 'Warning: Found radmc.inp AND radmc3d.inp. Taking radmc3d.inp!'
  endif
  if(fex1.or.fex2) then 
     if(fex1) then
        open(unit=1,file='radmc3d.inp')
     else
        open(unit=1,file='radmc.inp')
     endif
     call read_input_file()
     close(1)
  !$ call parse_input_integer('setthreads@                   ',setthreads)
     call parse_input_integer('mc_safetymode@                ',mc_safetymode)
     call parse_input_integer8('nphot@                        ',rt_mcparams%nphot_therm)
     call parse_input_integer8('nphot_therm@                  ',rt_mcparams%nphot_therm)
     call parse_input_integer8('nphot_scat@                   ',rt_mcparams%nphot_scat)
     call parse_input_integer8('nphot_spec@                   ',rt_mcparams%nphot_spec)
     call parse_input_integer8('nphot_mono@                   ',rt_mcparams%nphot_mono)
     call parse_input_integer('iseed@                        ',iseed)
     call parse_input_integer('imethod@                      ',imethod)
     call parse_input_integer('ifast@                        ',rt_mcparams%ifast)
     call parse_input_integer('iranfreqmode@                 ',rt_mcparams%iranfreqmode)
     call parse_input_double ('enthres@                      ',rt_mcparams%enthres)
     call parse_input_integer('cntdump@                      ',rt_mcparams%cntdump)
     call parse_input_integer('irestart@                     ',rt_mcparams%irestart)
     call parse_input_integer('itempdecoup@                  ',rt_mcparams%itempdecoup)
     call parse_input_integer('debug_write_stats@            ',rt_mcparams%debug_write_stats)
     call parse_input_integer('iquantum@                     ',incl_quantum)
     call parse_input_integer('incl_quantum@                 ',incl_quantum)
     call parse_input_integer('istar_sphere@                 ',inc_star_size)
     call parse_input_double ('nphotdiff@                    ',rt_mcparams%nphotdiff)
     if (rt_mcparams%nphotdiff.gt.0) rt_mcparams%debug_write_stats=1
     call parse_input_integer('nphotdiff_type@               ',rt_mcparams%nphotdiff_type)
     call parse_input_double ('errtol@                       ',rt_mcparams%errtoldiff)
     call parse_input_double ('errtoldiff@                   ',rt_mcparams%errtoldiff)
     call parse_input_integer('nvstr@                        ',rt_mcparams%niter_vstruct)
     call parse_input_integer('niter_vstruct@                ',rt_mcparams%niter_vstruct)
     call parse_input_double ('vserrtol@                     ',rt_mcparams%vserrtol)
     call parse_input_integer('ivstrt@                       ',rt_mcparams%ivstrt)
     call parse_input_integer('ntemp@                        ',rt_mcparams%ntemp)
     call parse_input_double ('temp0@                        ',rt_mcparams%temp0)
     call parse_input_double ('temp1@                        ',rt_mcparams%temp1)
     call parse_input_integer('countwrite@                   ',rt_mcparams%countwrite)
     call parse_input_integer('camera_tracemode@             ',camera_tracemode)
     idum=-1
     call parse_input_integer('writeimage_unformatted@       ',idum)
     if(idum.eq.0) then
        writeimage_unformatted=.false.
     elseif(idum.eq.1) then
        writeimage_unformatted=.true.
     endif
     call parse_input_integer('scattering_mode@              ',scattering_mode_max)
     call parse_input_integer('scattering_mode_max@          ',scattering_mode_max)
     call parse_input_integer('alignment_mode@               ',alignment_mode)
     call parse_input_double ('mc_scat_maxtauabs@            ',mc_scat_maxtauabs)
     call parse_input_integer('scat_munr@                    ',scat_munr)
     !call parse_input_integer('scat_phinr@                   ',scat_phinr)
     call parse_input_integer('dust_2daniso_nphi@            ',dust_2daniso_nphi)
     call parse_input_integer('mc_weighted_photons@          ',wghtphot)
     call parse_input_integer('optimized_motion@             ',ioptmot)
     call parse_input_double ('optim_dtau@                   ',rt_mcparams%optim_dtau)
     call parse_input_double ('thermal_boundary_xl@          ',thermal_bc_temp(1,1))
     call parse_input_double ('thermal_boundary_xr@          ',thermal_bc_temp(2,1))
     call parse_input_double ('thermal_boundary_yl@          ',thermal_bc_temp(1,2))
     call parse_input_double ('thermal_boundary_yr@          ',thermal_bc_temp(2,2))
     call parse_input_double ('thermal_boundary_zl@          ',thermal_bc_temp(1,3))
     call parse_input_double ('thermal_boundary_zr@          ',thermal_bc_temp(2,3))
     call parse_input_integer('lines_partition_ntempint@     ',lines_partition_ntempint)
     call parse_input_double ('lines_partition_temp0@        ',lines_partition_temp0)
     call parse_input_double ('lines_partition_temp1@        ',lines_partition_temp1)
     call parse_input_integer('lines_show_pictograms@        ',lines_show_pictograms)
     idum=0
     call parse_input_integer('modified_random_walk@         ',idum)
     if(idum.eq.0) then
        rt_mcparams%mod_random_walk = .false.
     else
        rt_mcparams%mod_random_walk = .true.
     endif
     call parse_input_integer('mrw_count_trigger@            ',rt_mcparams%mrw_count_trigger)
     call parse_input_integer('mrw_db_ntemp@                 ',rt_mcparams%mrw_db_ntemp)
     call parse_input_double ('mrw_db_temp0@                 ',rt_mcparams%mrw_db_temp0)
     call parse_input_double ('mrw_db_temp1@                 ',rt_mcparams%mrw_db_temp1)
     call parse_input_double ('mrw_enthres@                  ',rt_mcparams%mrw_enthres)
     call parse_input_double ('mrw_tauthres@                 ',rt_mcparams%mrw_tauthres)
     call parse_input_double ('mrw_gamma@                    ',rt_mcparams%mrw_gamma)
     call parse_input_double ('mrw_taustepback@              ',rt_mcparams%mrw_taustepback)
     call parse_input_double ('mrw_tempthres@                ',rt_mcparams%mrw_tempthres)
     idum=0
     call parse_input_integer('unformatted@                  ',idum)
     if(idum.eq.0) then
        rto_style=1
     else
        rto_style=2
     endif
     call parse_input_integer('rto_style@                    ',rto_style)
     idum=0
     call parse_input_integer('rto_single@                   ',idum)
!     call parse_input_integer('singleprec_output@            ',idum)
     if(idum.eq.0) then
        rto_single=.false.
     else
        rto_single=.true.
     endif
     call parse_input_integer('camera_nrrefine@              ',camera_nrrefine)
     call parse_input_integer('camera_incl_stars@            ',camera_incl_stars)
     call parse_input_double ('camera_refine_criterion@      ',camera_refine_criterion)
     call parse_input_double ('camera_starsphere_nrpix@      ',camera_starsphere_nrpix)
     call parse_input_integer('camera_localobs_projection@   ',camera_localobs_projection)
     call parse_input_double ('camera_localobs_zenith@       ',camera_localobs_zenith)
     call parse_input_double ('camera_spher_cavity_relres@   ',camera_spher_cavity_relres)
!     call parse_input_double ('camera_min_aspectratio@       ',camera_min_aspectratio)
     call parse_input_double ('camera_min_dangle@            ',mindang)
     call parse_input_double ('camera_max_dangle@            ',camera_max_dangle)
     call parse_input_double ('camera_min_drr@               ',mindrr)
     call parse_input_integer('camera_interpol_jnu@          ',interpoljnu)
     call parse_input_double ('camera_maxdphi@               ',camera_maxdphi)
     idum=0
     call parse_input_integer('camera_diagnostics_subpix@    ',idum)
     if(idum.eq.0) then
        camera_diagnostics_subpix = .false.
     else
        camera_diagnostics_subpix = .true.
     endif
     call parse_input_integer('sources_interpol_jnu@         ',interpoljnu)
!     call parse_input_double ('lines_maxdoppler@             ',lines_maxdoppler)
     call parse_input_integer('lines_mode@                   ',lines_mode)
     call parse_input_integer('lines_autosubset@             ',iautosubset)
     call parse_input_double ('lines_widthmargin@            ',lines_widthmargin)
     idum=0
     call parse_input_integer('lines_slowlvg_as_alternative@ ',idum)
     if(idum.eq.0) then
        lines_slowlvg_as_alternative = .false.
     else
        lines_slowlvg_as_alternative = .true.
     endif
     call parse_input_double ('catch_doppler_resolution@     ',camera_catch_doppler_resolution)
     string = '###'
     len    = 3
     call parse_input_word   ('lines_mode_keyword@           ',string,len)
     if(string(1:3).ne.'###') then
        if((string(1:3).eq.'lte').and.(len.eq.3)) then
           lines_mode = -1
        elseif((string(1:7).eq.'ltefast').and.(len.eq.7)) then
           lines_mode = 1
        else
           write(stdo,*) 'ERROR in parsing radmc3d.inp: lines_mode_keyword'
           write(stdo,*) '      value is unknown.'
           stop
        endif
     endif
     ! (Line profile stuff added by Thomas Peters 2011)
     call parse_input_integer('lines_profile@                ',lines_profile)
     call parse_input_double ('lines_tbg@                    ',lines_tbg)
     call parse_input_integer('lines_nonlte_maxiter@         ',lines_nonlte_maxiter)
     call parse_input_double ('lines_nonlte_convcrit@        ',lines_nonlte_convcrit)
     call flush(stdo)
     call parse_input_integer('tgas_eq_tdust@                ',tgas_eq_tdust)
     call parse_input_integer('subbox_nx@                    ',subbox_nx)
     call parse_input_integer('subbox_ny@                    ',subbox_ny)
     call parse_input_integer('subbox_nz@                    ',subbox_nz)
     call parse_input_double ('subbox_x0@                    ',subbox_x0)
     call parse_input_double ('subbox_x1@                    ',subbox_x1)
     call parse_input_double ('subbox_y0@                    ',subbox_y0)
     call parse_input_double ('subbox_y1@                    ',subbox_y1)
     call parse_input_double ('subbox_z0@                    ',subbox_z0)
     call parse_input_double ('subbox_z1@                    ',subbox_z1)
     call parse_input_double ('subbox_phi1@                  ',subbox_phi1)
     call parse_input_double ('subbox_theta@                 ',subbox_theta)
     call parse_input_double ('subbox_phi2@                  ',subbox_phi2)
     !
     ! Switches for including or excluding certain processes
     !
     call parse_input_integer('incl_dust@                    ',inc_dust)
     call parse_input_integer('incl_lines@                   ',inc_lines)
     call parse_input_integer('incl_freefree@                ',inc_freefree)
     call parse_input_integer('incl_userdef_srcalp@          ',inc_usersrc)
     !
     ! Call the userdef module for any user-defined keywords
     !
     call userdef_parse_main_namelist()
     !
     ! Done
     !
     call check_all_lines_ok(ierror)
     if(ierror) stop
  else
     write(stdo,*) 'ERROR: Could not find file radmc3d.inp nor radmc.inp.'
     write(stdo,*) '       In this version of RADMC-3D this file is mandatory.'
     write(stdo,*) '       You can keep it empty if you like, but it must exist.'
     stop
  endif
  !
  ! Interpret some options
  !
  if(ioptmot.eq.0) then
     rt_mcparams%optimized_motion=.false.
     write(stdo,*) 'Optmized motion switched off...'
  endif
  if(ioptmot.eq.1) then
     rt_mcparams%optimized_motion=.true.
     write(stdo,*) 'Optmized motion switched on...'
  endif
  ! 
  ! Interpret the switching on/off of various processes
  !
  if(inc_dust.eq.0) then
     rt_incl_dust=.false.
  endif
  if(inc_dust.eq.1) then
     rt_incl_dust=.true.
  endif
  !
  if(inc_lines.eq.0) then
     rt_incl_lines=.false.
  endif
  if(inc_lines.eq.1) then
     rt_incl_lines=.true.
  endif
  !
  if(inc_freefree.eq.0) then
     rt_incl_gascont_freefree=.false.
  endif
  if(inc_freefree.eq.1) then
     rt_incl_gascont=.true.
     rt_incl_gascont_freefree=.true.
  endif
  if(inc_usersrc.eq.0) then
     rt_incl_userdef_srcalp=.false.
  endif
  if(inc_usersrc.eq.1) then
     rt_incl_userdef_srcalp=.true.
  endif
  if(inc_star_size.eq.0) then
     star_sphere = .false.
  endif
  if(inc_star_size.eq.1) then
     star_sphere = .true.
  endif
  if(wghtphot.eq.0) then
     mc_weighted_photons = .false.
  endif
  if(wghtphot.eq.1) then
     mc_weighted_photons = .true.
  endif
  if(interpoljnu.eq.0) then
     write(stdo,*) 'Switching sources_interpol_jnu to .false.'
     sources_interpol_jnu = .false.
  endif
  if(interpoljnu.eq.1) then
     write(stdo,*) 'Switching sources_interpol_jnu to .true.'
     sources_interpol_jnu = .true.
  endif
  if(iautosubset.eq.0) then
     lines_autosubset = .false.
  endif
  if(iautosubset.eq.1) then
     lines_autosubset = .true.
  endif
  incl_thermbc = 0
  do bc_idir=1,3
     do bc_ilr=1,2
        if(thermal_bc_temp(bc_ilr,bc_idir).gt.0.d0) then
           thermal_bc_active(bc_ilr,bc_idir) = .true.
           incl_thermbc = 1
        else
           thermal_bc_active(bc_ilr,bc_idir) = .false.
        endif
     enddo
  enddo
  !
  ! Make warnings if necessary
  !
  if(mindrr.gt.0.d0) then
     camera_min_drr = mindrr
     write(stdo,*) 'Warning: You set camera_min_drr in radmc3d.inp to ',mindrr
  endif
  if(mindang.gt.0.d0) then
     camera_min_dangle = mindang
     write(stdo,*) 'Warning: You set camera_min_dangle in radmc3d.inp to ',mindang
  endif
  !
  if(camera_catch_doppler_resolution.gt.1.d0) then
     camera_catch_doppler_resolution = 1.d0/camera_catch_doppler_resolution
  endif
  !
end subroutine read_radmcinp_file


!-------------------------------------------------------------------------
!                    INTERPRET COMMAND-LINE OPTIONS
!-------------------------------------------------------------------------
subroutine interpet_command_line_options(gotit,fromstdi,quit)
  use rtglobal_module
  use camera_module
  use constants_module
  use userdef_module
  implicit none
  character*100 :: buffer
  integer :: iarg,numarg
  double precision :: dum,zoom_x0,zoom_x1,zoom_y0,zoom_y1,px,py
  integer :: idum
  logical :: gotit,flagzoom,truepix,truezoom,quit,fromstdi
  !
  gotit    = .false.
  quit     = .false.
  flagzoom = .false.
  truepix  = .false.
  truezoom = .false.
  if(fromstdi) then
     numarg   = 1000
  else
     numarg   = iargc()
  endif
  iarg     = 1
  do while(iarg.le.numarg) 
     !
     ! Get the next argument 
     !
     call ggetarg(iarg,buffer,fromstdi)
     iarg = iarg+1
     !
     ! Interpret this
     !
     if((buffer(1:5).eq.'child').and.(.not.fromstdi)) then
        radmc_as_child = .true.
        gotit = .true.
     elseif((buffer(1:5).eq.'enter').and.fromstdi) then
        numarg = 0
     elseif((buffer(1:4).eq.'quit').and.fromstdi) then
        gotit = .true.
        numarg = 0
        quit = .true.
     elseif((buffer(1:4).eq.'exit').and.fromstdi) then
        gotit = .true.
        numarg = 0
        quit = .true.
     elseif((buffer(1:10).eq.'writeimage').and.fromstdi) then
        do_writeimage = .true.
        gotit = .true.
        numarg = 0
!     elseif((buffer(1:11).eq.'uwriteimage').and.fromstdi) then
!        do_writeimage = .true.
!        writeimage_unformatted = .true.
!        gotit = .true.
!        numarg = 0
     elseif((buffer(1:9).eq.'writespec').and.fromstdi) then
        do_writespectrum = .true.
        gotit = .true.
        numarg = 0
     elseif(buffer(1:12).eq.'writescatsrc') then
        do_writescatsrc = .true.
        gotit = .true.
!     elseif(buffer(1:7).eq.'iobsdir') then
!        call ggetarg(iarg,buffer,fromstdi)
!        iarg = iarg+1
!        read(buffer,*) mcscat_current_dir
!        gotit = .true.
     elseif(buffer(1:13).eq.'writegridfile') then
        do_writegridfile = .true.
        gotit = .true.
     elseif(buffer(1:10).eq.'writemodel') then
        do_writemodel = .true.
        gotit = .true.
     elseif((buffer(1:7).eq.'sample-').or.(buffer(1:7).eq.'sample_')) then
        do_writesample = .true.
        if((buffer(8:15).eq.'dusttemp')) then
           do_sample_dusttemp = .true.
        elseif((buffer(8:15).eq.'dustdens')) then
           do_sample_dustdens = .true.
        elseif((buffer(8:15).eq.'levelpop')) then
           do_sample_levelpop = .true.
        else
           write(stdo,*) 'Sorry, I do not know what to sample: ',trim(buffer)
           do_writesample = .false.
        endif
        gotit = .true.
     elseif((buffer(1:8).eq.'writepop')) then
        do_writepop = .true.
        gotit = .true.
     elseif((buffer(1:8).eq.'calcpop')) then
        do_calcpop = .true.
        do_writepop = .true.
        gotit = .true.
     elseif(buffer(1:16).eq.'respondwhenready') then
        do_respondwhenready = .true.
        gotit = .true.
     elseif(buffer(1:9).eq.'resetseed') then
        do_resetseed = .true.
        gotit = .true.
     elseif(buffer(1:12).eq.'notresetseed') then
        do_resetseed = .false.
        gotit = .true.
     elseif(buffer(1:7).eq.'usetree') then
        amr_always_use_tree = .true.
        gotit = .true.
     elseif(buffer(1:7).eq.'mctherm') then
        !
        ! Do the Bjorkman & Wood Monte Carlo 
        ! This computes the dust temperatures inside the model
        ! Important to do before computing dust emission spectra and images
        !
        do_montecarlo_therm = .true.
        gotit = .true.
     elseif(buffer(1:6).eq.'mcmono') then
        !
        ! Do the monochromatic Monte Carlo 
        ! This computes the local radiation field inside the model
        ! Useful for other models, e.g. chemistry
        !
        do_montecarlo_mono = .true.
        gotit = .true.
     elseif(buffer(1:8).eq.'myaction') then
        !
        ! Do the userdef action 
        !
        do_userdef_action = .true.
        gotit = .true.
     elseif(buffer(1:11).eq.'imageunform') then
        !
        ! Switch on the unformatted writing of the images
        !
        writeimage_unformatted = .true.
        gotit = .true.
     elseif(buffer(1:14).eq.'imageformatted') then
        !
        ! Switch back to the formatted writing of the images
        !
        writeimage_unformatted = .false.
        gotit = .true.
!     elseif(buffer(1:7).eq.'vstruct') then
!        !
!        ! Do the Bjorkman & Wood Monte Carlo coupled to the equation
!        ! of vertical structure
!        !
!        do_vstruct = .true.
!        gotit = .true.
     elseif(buffer(1:8).eq.'spectrum') then
        !
        ! Produce a spectrum
        !
        do_raytrace_spectrum = .true.
        gotit = .true.
     elseif(buffer(1:3).eq.'sed') then
        !
        ! Produce a continuum spectrum
        ! By default the global frequencies are used
        !
        camera_setfreq_global = .true.
        do_raytrace_spectrum = .true.
!        rt_incl_dust = .true.
!        rt_incl_lines = .false.
        gotit = .true.
     elseif(buffer(1:6).eq.'sloppy') then
        camera_min_dangle          = 0.1d0
        camera_min_drr             = 0.1d0
        camera_spher_cavity_relres = 0.1d0
        gotit = .true.
     elseif(buffer(1:5).eq.'image') then
        !
        ! Make an image
        !
        do_raytrace_image    = .true.
        gotit = .true.
     elseif(buffer(1:5).eq.'movie') then
        !
        ! Make an image
        !
        do_raytrace_movie    = .true.
        gotit = .true.
     elseif(buffer(1:7).eq.'tausurf') then
        !
        ! Make the tau=tausurf surface
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read tausurf argument.'
           write(stdo,*) '      Expecting 1 float after tausurf.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) camera_tausurface
        do_raytrace_tausurf  = .true.
        gotit = .true.
     elseif(buffer(1:4).eq.'circ') then
        !
        ! Make an image
        !
        circular_images    = .true.
        gotit = .true.
     elseif(buffer(1:6).eq.'stokes') then
        !
        ! Include all four Stokes components in the image or spectrum
        !
        camera_stokesvector = .true.
        gotit = .true.
     elseif(buffer(1:8).eq.'nostokes') then
        !
        ! Include only intensity in the image or spectrum (not Q U or V)
        !
        camera_stokesvector = .false.
        gotit = .true.
     elseif(buffer(1:10).eq.'autosubset') then
        !
        ! For lines (with 1 <= lines_mode <= 9): Automatically
        ! select the level subset.
        !
        lines_autosubset = .true.
        gotit = .true.
     elseif(buffer(1:12).eq.'noautosubset') then
        !
        ! For lines (with 1 <= lines_mode <= 9): Automatically
        ! select the level subset.
        !
        lines_autosubset = .false.
        gotit = .true.
     elseif(buffer(1:6).eq.'nodust') then
        !
        ! Do not include dust in images/spectra. Use this last on the command line
        ! to ensure it has effect.
        !
        rt_incl_dust     = .false.
        gotit = .true.
     elseif(buffer(1:8).eq.'incldust') then
        !
        ! Include dust in images/spectra. Use this last on the command line
        ! to ensure it has effect.
        !
        rt_incl_dust     = .true.
        gotit = .true.
     elseif(buffer(1:6).eq.'noscat') then
        !
        ! Do not include dust scattering in images/spectra. 
        !
        scattering_mode_max = 0
        gotit = .true.
     elseif(buffer(1:8).eq.'inclscat') then
        !
        ! Include dust scattering in images/spectra.
        !
        scattering_mode_max = scattering_mode_def
        gotit = .true.
     elseif(buffer(1:6).eq.'noline') then
        !
        ! Do not include lines in images/spectra. Use this last on the command line
        ! to ensure it has effect.
        !
        rt_incl_lines   = .false.
        gotit = .true.
     elseif(buffer(1:8).eq.'inclline') then
        !
        ! Include lines in images/spectra. Use this last on the command line
        ! to ensure it has effect.
        !
        rt_incl_lines   = .true.
        gotit = .true.
     elseif(buffer(1:9).eq.'nogascont') then
        !
        ! Do not include gas continuum in images/spectra. Use this last on the command line
        ! to ensure it has effect.
        !
        rt_incl_gascont  = .false.
        gotit = .true.
     elseif(buffer(1:11).eq.'inclgascont') then
        !
        ! Include gas continuum in images/spectra. Use this last on the command line
        ! to ensure it has effect.
        !
        rt_incl_gascont   = .true.
        gotit = .true.
     elseif(buffer(1:10).eq.'nofreefree') then
        !
        ! Do not include free-free continuum in images/spectra. Use this last on the command line
        ! to ensure it has effect.
        !
        rt_incl_gascont_freefree  = .false.
        gotit = .true.
     elseif(buffer(1:11).eq.'inclfreefree') then
        !
        ! Include gas free-free in images/spectra. Use this last on the command line
        ! to ensure it has effect.
        !
        rt_incl_gascont          = .true.
        rt_incl_gascont_freefree = .true.
        gotit = .true.
     elseif(buffer(1:8).eq.'linelist') then
        !
        ! Write a list of line positions to standard output
        !
        rt_incl_lines       = .true.
        lines_make_linelist = .true.
        gotit = .true.
     elseif((buffer(1:5).eq.'color').or.(buffer(1:9).eq.'loadcolor')) then
        !
        ! Use multiple wavelengths for single image. The wavelengths are
        ! selections from the global wavelength array. Read these from 
        ! the file color_inus.inp
        !
        camera_loadcolor = .true.
!        rt_incl_dust = .true.
!        rt_incl_lines = .false.
        gotit = .true.
     elseif(buffer(1:10).eq.'loadlambda') then
        !
        ! Use multiple wavelengths for single image. The wavelengths are
        ! directly read from camera_frequency.inp or 
        ! camera_wavelength_micron.inp.
        !
        camera_loadlambdas = .true.
        gotit = .true.
     elseif(buffer(1:5).eq.'iline') then
        !
        ! Specify a line around which to put a window
        ! Make also a few defaults
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read sizeau.'
           write(stdo,*) '      Expecting 1 integer after iline.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) lines_user_iline
        ! Defaults:
        if(lines_user_ispec.eq.0)  lines_user_ispec = 1
        rt_incl_lines     = .true.
        gotit = .true.
     elseif(buffer(1:8).eq.'imolspec') then
        !
        ! Specify a line around which to put a window
        ! Make also a few defaults
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read sizeau.'
           write(stdo,*) '      Expecting 1 integer after imolspec.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) lines_user_ispec
        ! Defaults:
        if(lines_user_iline.eq.0)  lines_user_iline = 1
        rt_incl_lines     = .true.
        gotit = .true.
     elseif((buffer(1:8).eq.'linenlam').or.(buffer(1:9).eq.'linenfreq')) then
        !
        ! Specify the number of wavelength points around this line
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read sizeau.'
           write(stdo,*) '      Expecting 1 integer after linenlam.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) lines_user_nrfreq
        if((lines_user_nrfreq.gt.1).and.(lines_user_widthkms.eq.0.d0)) then
           lines_user_widthkms = lines_user_widthkms_default
        endif
        gotit = .true.
     elseif(buffer(1:8).eq.'widthkms') then
        !
        ! Specify the half-width of the wavelength grid 
        ! (from line center to one side) in km/s
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read sizeau.'
           write(stdo,*) '      Expecting 1 floating point number after widthkms.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) lines_user_widthkms
        if(lines_user_nrfreq.eq.1) then
           lines_user_nrfreq = lines_user_nrfreq_default
        endif
        gotit = .true.
     elseif(buffer(1:4).eq.'vkms') then
        !
        ! Specify the zero-point of the wavelength grid around a line.
        ! Default is 0, meaning that the wavelength grid is centered
        ! around the line center.
        !
        ! NOTE: If the user wants to make a single-wavelength image
        !       at exactly vkms, then simply specify imolspec, iline 
        !       and vkms. That will do. 
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read sizeau.'
           write(stdo,*) '      Expecting 1 floating point number after vkms.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) lines_user_kms0
        gotit = .true.
     elseif(buffer(1:8).eq.'lineinfo') then
        !
        ! Write some line information to file
        !
        do_write_lineinfo = .true.
        gotit = .true.
     elseif(buffer(1:8).eq.'useapert') then
        !
        ! If you make a spectrum, then use the aperture information 
        ! from the file aperture_info.inp.
        !
        camera_use_aperture_info = .true.
        gotit = .true.
     elseif(buffer(1:7).eq.'noapert') then
        !
        ! If you make a spectrum, then use the default way: compute the
        ! flux inside a rectangular region specified in the same way
        ! as for images. This is the default way.
        !
        camera_use_aperture_info = .false.
        gotit = .true.
     elseif(buffer(1:3).eq.'dpc') then
        !
        ! Set the distance of the observer (in observer-at-infinity mode).
        ! This is only necessary if the user wants to make spectra with
        ! the use of aperture information.
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read dpc.'
           write(stdo,*) '      Expecting 1 float after dpc.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) camera_observer_distance_pc
        gotit = .true.
     !$ elseif(buffer(1:10).eq.'setthreads') then
     !
     !!$ Set the number of threads used for the parallel version of RADMC-3D
     !
     !$   if(iarg.gt.numarg) then
     !$      write(stdo,*) 'ERROR while reading command line options: cannot read setthreads.'
     !$      write(stdo,*) '      Expecting 1 integer after setthreads.'
     !$      stop
     !$   endif
     !$   call ggetarg(iarg,buffer,fromstdi)
     !$   iarg = iarg+1
     !$   read(buffer,*) setthreads
     !$   gotit = .true.
     elseif(buffer(1:11).eq.'nphot_therm') then
        !
        ! Set the nr of photons for the therm monte carlo
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read nphot_therm.'
           write(stdo,*) ' Expecting 1 integer after nphot_therm.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) rt_mcparams%nphot_therm
        gotit = .true.
     elseif(buffer(1:10).eq.'nphot_scat') then
        !
        ! Set the nr of photons for the scattering monte carlo
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read nphot_scat.'
           write(stdo,*) '      Expecting 1 integer after nphot_scat.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) rt_mcparams%nphot_scat
        gotit = .true.
     elseif(buffer(1:10).eq.'nphot_spec') then
        !
        ! Set the nr of photons for the scattering monte carlo, but for
        ! use when spectra are made: here you do not need so many photon
        ! packages, since we will anyway integrate over the image to get the
        ! spectrum.
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read nphot_spec.'
           write(stdo,*) '      Expecting 1 integer after nphot_spec.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) rt_mcparams%nphot_spec
        gotit = .true.
     elseif(buffer(1:10).eq.'nphot_mono') then
        !
        ! Set the nr of photons for the monochromatic monte carlo
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read nphot_mono.'
           write(stdo,*) '      Expecting 1 integer after nphot_mono.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) rt_mcparams%nphot_mono
        gotit = .true.
     elseif(buffer(1:10).eq.'selectscat') then
        !
        ! Select which scattering to include. For instance, if you want to
        ! see the role of multiple scattering, you can do selectscat 2 1000000.
        ! Or if you are interested only in the second scattering you can do
        ! selectscat 2 2. Or only first scattering: selectscat 1 1.
        ! NOTE: Only for analysis or debugging, not for production runs!
        !
        if(iarg+1.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read selectscat.'
           write(stdo,*) '      Expecting 2 integers after selectscat.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) selectscat_iscat_first
        if(selectscat_iscat_first.lt.1) then
           write(stdo,*) 'ERROR: selectscat_iscat_first must be .ge.1'
           stop 3284
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) selectscat_iscat_last
        if(selectscat_iscat_last.lt.selectscat_iscat_first) then
           write(stdo,*) 'ERROR: selectscat_iscat_last must be .ge.selectscat_iscat_first'
           stop 3285
        endif
        gotit = .true.
     elseif(buffer(1:10).eq.'countwrite') then
        !
        ! Set how often RADMC-3D writes a progress line in a Monte Carlo run (default = 1000)
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read countwrite.'
           write(stdo,*) '      Expecting 1 integer after countwrite.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) rt_mcparams%countwrite
        gotit = .true.
     elseif(buffer(1:6).eq.'sizeau') then
        !
        ! Set the image size
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read sizeau.'
           write(stdo,*) '      Expecting 1 float after sizeau.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        camera_image_halfsize_x = 1.496d13*dum/2
        camera_image_halfsize_y = camera_image_halfsize_x
        if(camera_localobserver) then
           write(stdo,*) 'User specifies sizeau: Therefore switching OFF local observer mode'
           write(stdo,*) '                       (observer now at infinity)!'
        endif
        camera_localobserver       = .false.
     elseif(buffer(1:6).eq.'sizepc') then
        !
        ! Set the image size
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read sizepc.'
           write(stdo,*) '      Expecting 1 float after sizepc.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        camera_image_halfsize_x = 3.08572d18*dum/2
        camera_image_halfsize_y = camera_image_halfsize_x
        if(camera_localobserver) then
           write(stdo,*) 'User specifies sizepc: Therefore switching OFF local observer mode'
           write(stdo,*) '                       (observer now at infinity)!'
        endif
        camera_localobserver       = .false.
     elseif(buffer(1:6).eq.'zoomau') then
        !
        ! Set the image geometry precisely, by specifying (in units of AU)
        ! the left and right as well as lower and upper edges of the image,
        ! computed with respect to the pointing location (see pointau below).
        !
        if(iarg+3.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read zoomau.'
           write(stdo,*) '      Expecting 4 floats after zoomau: the x0,x1,y0,y1 of the image.'
           stop
        endif
        flagzoom = .true.
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        zoom_x0 = 1.496d13*dum
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        zoom_x1 = 1.496d13*dum
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        zoom_y0 = 1.496d13*dum
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        zoom_y1 = 1.496d13*dum
        !
        if(camera_localobserver) then
           write(stdo,*) 'User specifies zoomau: Therefore switching OFF local observer mode'
           write(stdo,*) '                       (observer now at infinity)!'
        endif
        camera_localobserver       = .false.
     elseif(buffer(1:6).eq.'zoompc') then
        !
        ! Set the image geometry precisely, by specifying (in units of pc)
        ! the left and right as well as lower and upper edges of the image,
        ! computed with respect to the pointing location (see pointpc below).
        !
        if(iarg+3.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read zoompc.'
           write(stdo,*) '      Expecting 4 floats after zoompc: the x0,x1,y0,y1 of the image.'
           stop
        endif
        flagzoom = .true.
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        zoom_x0 = 3.08572d18*dum
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        zoom_x1 = 3.08572d18*dum
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        zoom_y0 = 3.08572d18*dum
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        zoom_y1 = 3.08572d18*dum
        !
        if(camera_localobserver) then
           write(stdo,*) 'User specifies zoompc: Therefore switching OFF local observer mode'
           write(stdo,*) '                       (observer now at infinity)!'
        endif
        camera_localobserver       = .false.
     elseif(buffer(1:10).eq.'sizeradian') then
        !
        ! Set the image size in radian (for the local observer mode)
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read sizeradian.'
           write(stdo,*) '      Expecting 1 float after sizeradian.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        camera_image_halfsize_x = dum/2
        camera_image_halfsize_y = camera_image_halfsize_x
        if(.not.camera_localobserver) then
           write(stdo,*) 'User specifies sizeradian. Therefore switching on local observer mode!'
        endif
        camera_localobserver       = .true.
     elseif(buffer(1:10).eq.'zoomradian') then
        !
        ! Set the image geometry precisely, by specifying (in radian, i.e.
        ! for use with the local observer mode only!!)
        ! the left and right as well as lower and upper edges of the image,
        ! computed with respect to the pointing location (see pointau below).
        !
        if(iarg+3.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read zoomradian.'
           write(stdo,*) '      Expecting 4 floats after zoomradian: the x0,x1,y0,y1 of the image.'
           stop
        endif
        flagzoom = .true.
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        zoom_x0 = dum
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        zoom_x1 = dum
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        zoom_y0 = dum
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        zoom_y1 = dum
        !
        if(.not.camera_localobserver) then
           write(stdo,*) 'User specifies sizeradian. Therefore switching on local observer mode!'
        endif
        camera_localobserver       = .true.
     elseif(buffer(1:10).eq.'projection') then
        !
        ! For the local observer mode, specify which projection to use
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read projection.'
           write(stdo,*) '      Expecting 1 integer after projection.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) camera_localobs_projection
        if(camera_localobs_projection.eq.1) then
           write(stdo,*) 'Using flat screen projection for local observer'
        elseif(camera_localobs_projection.eq.2) then
           write(stdo,*) 'Using dome projection for local observer'
        else
           write(stdo,*) 'Sorry, do not know local observer projection type ',&
                camera_localobs_projection
           stop
        endif
        gotit = .true.
     elseif(buffer(1:6).eq.'zenith') then
        !
        ! If you use fisheye (Planetarium Dome, OMNIMAX) local observer
        ! projection, you may need to shift the zenith a bit away from the pointing
        ! direction, so that the viewers do not always have to look straight up
        ! to see the interesting stuff.
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read zenith.'
           write(stdo,*) '      Expecting 1 float after zenith: y-offset.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        camera_localobs_zenith = dum
     elseif(buffer(1:7).eq.'truepix') then
        !
        ! Force RADMC with the zoomau or zoompc option to use the actual
        ! number of pixels specified in each direction. This could mean
        ! that non-square pixels (with width unequal to height) would 
        ! result, but perhaps someone may want this badly for some reason.
        !
        truepix = .true.
        !
     elseif(buffer(1:7).eq.'truezoom') then
        !
        ! Normally, if the zoom option is used (zoomau or zoompc) the
        ! number of pixels in one of the two direction may be reduced
        ! so as to ensure that square pixels result if the image itself
        ! it not perfectly square. However, since the number of pixels
        ! is an integer, RADMC-3D might have to slightly reduce the
        ! image size to ensure that the narrower direction is an integer
        ! number times the pixel size. If the user really wants to ensure
        ! that *exactly* the given zoom area is used, and he/she is OK
        ! that the pixels are then only *approximately* square, then 
        ! the 'truezoom' option must be used.
        !
        truezoom = .true.
        !
     elseif(buffer(1:7).eq.'pointau') then
        !
        ! Set the pointing position as a 3-D location in the model
        !
        do idum=1,3
           if(iarg.gt.numarg) then
              write(stdo,*) 'ERROR while reading command line options: cannot read pointau.'
              write(stdo,*) '      Expecting 3 floats after pointau.'
              stop
           endif
           call ggetarg(iarg,buffer,fromstdi)
           iarg = iarg+1
           read(buffer,*) dum
           camera_pointing_position(idum) = 1.496d13*dum
        enddo
     elseif(buffer(1:7).eq.'pointpc') then
        !
        ! Set the pointing position as a 3-D location in the model
        !
        do idum=1,3
           if(iarg.gt.numarg) then
              write(stdo,*) 'ERROR while reading command line options: cannot read pointpc.'
              write(stdo,*) '      Expecting 3 floats after pointpc.'
              stop
           endif
           call ggetarg(iarg,buffer,fromstdi)
           iarg = iarg+1
           read(buffer,*) dum
           camera_pointing_position(idum) = 3.08572d18*dum
        enddo
     elseif(buffer(1:8).eq.'locobsau') then
        !
        ! Set the local observer position in units of au
        !
        do idum=1,3
           if(iarg.gt.numarg) then
              write(stdo,*) 'ERROR while reading command line options: cannot read locobsau.'
              write(stdo,*) '      Expecting 3 floats after locobsau.'
              stop
           endif
           call ggetarg(iarg,buffer,fromstdi)
           iarg = iarg+1
           read(buffer,*) dum
           camera_observer_position(idum) = 1.496d13*dum
        enddo
        if(.not.camera_localobserver) then
           write(stdo,*) 'User specifies locobsau. Therefore switching on local observer mode!'
        endif
        camera_localobserver       = .true.
     elseif(buffer(1:8).eq.'locobspc') then
        !
        ! Set the local observer position in units of pc
        !
        do idum=1,3
           if(iarg.gt.numarg) then
              write(stdo,*) 'ERROR while reading command line options: cannot read locobspc.'
              write(stdo,*) '      Expecting 3 floats after locobspc.'
              stop
           endif
           call ggetarg(iarg,buffer,fromstdi)
           iarg = iarg+1
           read(buffer,*) dum
           camera_observer_position(idum) = 3.08572d18*dum
        enddo
        if(.not.camera_localobserver) then
           write(stdo,*) 'User specifies locobspc. Therefore switching on local observer mode!'
        endif
        camera_localobserver       = .true.
     elseif(((buffer(1:4).eq.'incl').and.(len_trim(buffer).eq.4)).or.(buffer(1:5).eq.'theta')) then
        !
        ! Set the theta angle toward which the observer is. For a disk
        ! in the x-y plane this would be the inclination, hence the name.
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read incl/theta.'
           write(stdo,*) '      Expecting 1 float after incl or theta.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        camera_observer_degr_theta = dum
        if(camera_localobserver) then
           write(stdo,*) 'User specifies incl:   Therefore switching OFF local observer mode'
           write(stdo,*) '                       (observer now at infinity)!'
        endif
        camera_localobserver       = .false.
     elseif(buffer(1:3).eq.'phi') then
        !
        ! Set the phi angle toward which the observer is. 
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read phi.'
           write(stdo,*) '      Expecting 1 float after phi.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        camera_observer_degr_phi = dum
        if(camera_localobserver) then
           write(stdo,*) 'User specifies phi:    Therefore switching OFF local observer mode'
           write(stdo,*) '                       (observer now at infinity)!'
        endif
        camera_localobserver       = .false.
     elseif(buffer(1:6).eq.'posang') then
        !
        ! Set the position angle of the camera
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read posang.'
           write(stdo,*) '      Expecting 1 float after posang.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) dum
        camera_pointing_degr_posang = dum
     elseif(buffer(1:4).eq.'npix') then
        !
        ! Set the nr of pixels of the image
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read npix.'
           write(stdo,*) '      Expecting 1 integer after npix.'
           stop
        endif
        if(buffer(1:5).eq.'npixx') then
           call ggetarg(iarg,buffer,fromstdi)
           iarg = iarg+1
           read(buffer,*) idum
           camera_image_nx = idum
        elseif(buffer(1:5).eq.'npixy') then
           call ggetarg(iarg,buffer,fromstdi)
           iarg = iarg+1
           read(buffer,*) idum
           camera_image_ny = idum
        else
           call ggetarg(iarg,buffer,fromstdi)
           iarg = iarg+1
           read(buffer,*) idum
           camera_image_nx = idum
           camera_image_ny = idum
        endif
     elseif(buffer(1:11).eq.'lambdarange') then
        !
        ! A simple way to set a range of wavelengths for the camera
        !
        if(camera_theinu.gt.0) then
           write(stdo,*) 'ERROR: You cannot set the inu AND the lambdarange simultaneously'
           stop
        endif
        if(camera_loadcolor) then
           write(stdo,*) 'ERROR: You cannot set the colors AND the lambdarange simultaneously'
           stop
        endif
        if(iarg.gt.numarg-1) then
           write(stdo,*) 'ERROR while reading command line options: cannot read lambdarange.'
           write(stdo,*) '      Expecting 2 reals after lambdarange.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) camera_lambdamic
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) camera_lambdamic1
        !
        ! Set a default number of points
        !
        if(camera_range_nlam.le.0) then
           camera_range_nlam = 100
        endif
        !
     elseif(buffer(1:4).eq.'nlam') then
        !
        ! Used in combination with lambdarange above
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read nlam.'
           write(stdo,*) '      Expecting 1 integer after nlam.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) camera_range_nlam
        !
     elseif((buffer(1:3).eq.'inu').or.(buffer(1:7).eq.'ilambda')) then
        !
        ! Set which wavelength point the image should be taken
        ! NOTE: Setting this value also immediately implies that we
        !       will use the global frequency array.
        !
        if(camera_lambdamic.gt.0) then
           write(stdo,*) 'ERROR: You cannot set the inu AND the nuhz or lambda simultaneously'
           stop
        endif
        if(camera_loadcolor) then
           write(stdo,*) 'ERROR: You cannot set the colors AND the nuhz or lambda simultaneously'
           stop
        endif
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read inu/ilambda.'
           write(stdo,*) '      Expecting 1 integer after inu/ilambda.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) idum
        camera_theinu = idum
     elseif((buffer(1:4).eq.'nuhz').or.(buffer(1:6).eq.'lambda')) then
        !
        ! Set the precise wavelength where the image should be taken
        ! NOTE: Contrary to the old RADMC, where the nearest point in
        !       the global wavelength grid is taken, we here use the
        !       precise wavelength that is here specified. We will
        !       interpolate the opacities (using their original 
        !       wavelength grid, if not the same as the global one)
        !       and use these for the imaging and single-frequency
        !       Monte Carlo simulations.
        !
        if(camera_theinu.gt.0) then
           write(stdo,*) 'ERROR: You cannot set the inu AND the nuhz or lambda simultaneously'
           stop
        endif
        if(camera_loadcolor) then
           write(stdo,*) 'ERROR: You cannot set the colors AND the nuhz or lambda simultaneously'
           stop
        endif
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read nuhz/lambdamic.'
           write(stdo,*) '      Expecting 1 real after nuhz/lambdamic.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) camera_lambdamic
        if(buffer(1:4).eq.'nuhz') then
           camera_lambdamic = 1d4*cc/camera_lambdamic
        endif
     elseif((buffer(1:10).eq.'globalfreq').or.(buffer(1:9).eq.'globallam').or.&
            (buffer(1:5).eq.'allwl')) then
        !
        ! Use the global frequency array as the frequencies for the camera
        !
        camera_setfreq_global = .true.
        !
     elseif(buffer(1:8).eq.'nrrefine') then
        !
        ! Set the max number of refinement levels for images
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read nrrefine.'
           write(stdo,*) '      Expecting 1 integer after inu/ilambda.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) idum
        camera_nrrefine = idum
        if(camera_nrrefine.lt.0) then
           write(stdo,*) 'ERROR while reading command line options: nrrefine < 0'
           stop
        endif
     elseif(buffer(1:8).eq.'fluxcons') then
        !
        ! Make sure that as many refinement levels are included in the 
        ! images as necessary to guarantee flux conservation. More than
        ! 2^100 is not necessary, I guess...
        !
        camera_nrrefine = 100
     elseif(buffer(1:10).eq.'nofluxcons') then
        !
        ! Make images in a very simple way: one ray per pixel and no
        ! refinement. This is very fast, but may not give the correct
        ! flux in each pixel, especially if AMR is used.
        !
        camera_nrrefine = -1
     elseif(buffer(1:8).eq.'norefine') then
        !
        ! Identical to 'nofluxcons'
        !
        camera_nrrefine = -1
     elseif(buffer(1:6).eq.'second') then
        !
        ! Use second order integration of RT equation for images/spectra 
        ! (instead of the default first order). Will reset after each action.
        !
        camera_secondorder = .true.
     elseif(buffer(1:11).eq.'diag_subpix') then
        !
        ! Dump diagnostics for subpixeling
        !
        camera_diagnostics_subpix = .true.
     elseif((buffer(1:9).eq.'doppcatch').or.(buffer(1:6).eq.'dcatch').or.&
            (buffer(1:12).eq.'dopplercatch')) then
        !
        ! Use the "doppler catching" method for ray-tracing of molecular/atomic
        ! lines. The idea is to assure that smaller integration steps along the
        ! ray are carried out in those ray segments where a line may otherwise
        ! "doppler jump" over the wavelength-of-sight (being more than one local
        ! line width to the blue/red of the wavelength-of-sight at the start
        ! of the ray segment and more than one line width to the red/blue at
        ! the end). The "doppler catching" method requires second order
        ! integration (i.e. camera_secondorder=.true.), which means that this
        ! will use the vertex grid (grid cell corner vertices). We will
        ! automatically switch on the camera_secondorder here.
        !Use second order integration of RT equation for images/spectra 
        ! (instead of the default first order). Will reset after each action.
        !
        camera_catch_doppler_line = .true.
        camera_secondorder = .true.
     elseif(buffer(1:8).eq.'inclstar') then
        !
        ! Make sure that the stars are included in the images and SEDs,
        ! irrespective of what the radmc3d.inp says
        !
        camera_incl_stars=1
     elseif(buffer(1:6).eq.'nostar') then
        !
        ! Make sure that the stars are NOT included in the images and SEDs,
        ! irrespective of what the radmc3d.inp says
        !
        camera_incl_stars=0
     elseif(buffer(1:16).eq.'lambdasinglescat') then
        !
        ! Use do_lambda_starlight_single_scattering() routine instead of do_monte_carlo_scattering(),
        ! meaning that the scattering source function for the images is computed
        ! only for the direct (extincted) stellar light. No multiple scattering is
        ! computed. Only point source starlight is allowed as source. No thermal emission
        ! from dust or other diffuse sources can be used. 
        !
        ! Note: This mode is not Monte-Carlo style: it is deterministic.
        !
        camera_lambda_starlight_single_scat_mode = 1
        write(stdo,*) 'Warning: using lambda starlight single scattering mode now!'
     elseif(buffer(1:10).eq.'simplescat') then
        !
        ! Use do_lambda_starlight_single_scattering_simple() routine instead of
        ! do_monte_carlo_scattering(). This mode works only under very special
        ! conditions: spherical coordinates, a single pointlike star at the center,
        ! single scattering. The scattering source function for the images is computed
        ! only for the direct (extincted) stellar light. No multiple scattering is
        ! computed. Only point source starlight is allowed as source. No thermal emission
        ! from dust or other diffuse sources can be used. 
        !
        ! Note: This mode is not Monte-Carlo style: it is deterministic.
        !
        camera_lambda_starlight_single_scat_mode = 2
        write(stdo,*) 'Warning: using simple starlight single scattering mode now!'
     elseif(buffer(1:9).eq.'maxnrscat') then
        !
        ! Set the max number of scattering events included in the 
        ! scattering monte carlo for images/spectra
        !
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) idum
        mc_max_nr_scat_events = idum
     elseif(buffer(1:8).eq.'tracetau') then
        !
        ! Image the optical depth at the respective wavelength, instead of the intensity
        !
        camera_tracemode = -2
     elseif(buffer(1:11).eq.'tracenormal') then
        !
        ! Truly ray-trace the emission at the respective wavelength
        ! This is the normal way of making an image.
        !
        camera_tracemode = 1
     elseif(buffer(1:11).eq.'tracecolumn') then
        !
        ! Image the total density column, instead of the intensity
        !
        camera_tracemode = -1
     elseif(buffer(1:19).eq.'subbox_dust_density') then
        !
        ! Write a regularly spaced box of data: the dustdens array
        !
        do_write_subbox_dustdens = .true.
        gotit = .true.
     elseif(buffer(1:23).eq.'subbox_dust_temperature') then
        !
        ! Write a regularly spaced box of data: the dustdens array
        !
        do_write_subbox_dusttemp = .true.
        gotit = .true.
     elseif(buffer(1:15).eq.'subbox_levelpop') then
        !
        ! Write a regularly spaced box of data: the level populations
        !
        do_write_subbox_levelpop = .true.
        gotit = .true.
     elseif(buffer(1:12).eq.'subbox_xyz01') then
        !
        ! Dimensions of the box
        !
        if(iarg+5.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read '
           write(stdo,*) '      sub box dimensions (subbox_xyz01).'
           write(stdo,*) '      Expecting 6 reals: x0 x1 y0 y1 z0 z1 .'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) subbox_x0
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) subbox_x1
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) subbox_y0
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) subbox_y1
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) subbox_z0
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) subbox_z1
     elseif(buffer(1:11).eq.'subbox_phi1') then
        !
        ! Rotation angle phi1
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read '
           write(stdo,*) '      sub box rotation angle phi1.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) subbox_phi1
     elseif(buffer(1:11).eq.'subbox_phi2') then
        !
        ! Rotation angle phi2
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read '
           write(stdo,*) '      sub box rotation angle phi2.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) subbox_phi2
     elseif(buffer(1:12).eq.'subbox_theta') then
        !
        ! Rotation angle theta
        !
        if(iarg.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read '
           write(stdo,*) '      sub box rotation angle theta.'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) subbox_theta
     elseif(buffer(1:11).eq.'subbox_nxyz') then
        !
        ! Integer dimensions of the box
        !
        if(iarg+2.gt.numarg) then
           write(stdo,*) 'ERROR while reading command line options: cannot read '
           write(stdo,*) '      sub box size (subbox_nxyz).'
           write(stdo,*) '      Expecting 3 integers: nx ny nz .'
           stop
        endif
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) subbox_nx
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) subbox_ny
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) subbox_nz
        if((subbox_nx.lt.1).or.(subbox_ny.lt.1).or.(subbox_nz.lt.1)) then 
           write(stdo,*) 'ERROR in command lines: subbox size cannot be smaller than 1.'
           write(stdo,*) subbox_nx,subbox_ny,subbox_nz
           stop
        endif
     elseif(buffer(1:8).eq.'vtk_grid') then
        !
        ! Write out the vtk grid
        !
        amr_always_use_tree = .true.
        do_write_vtk_grid   = .true.
        gotit = .true.
     elseif(buffer(1:16).eq.'vtk_dust_density') then
        !
        ! Write out dust density in VTK format
        !
        amr_always_use_tree = .true.
        do_write_vtk_grid   = .true.
        do_write_vtk_dust_density = .true.
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) vtk_dust_ispec
        gotit = .true.
     elseif(buffer(1:20).eq.'vtk_dust_temperature') then
        !
        ! Write out dust temperature in VTK format
        !
        amr_always_use_tree = .true.
        do_write_vtk_grid   = .true.
        do_write_vtk_dust_temperature = .true.
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) vtk_dust_ispec
        gotit = .true.
     elseif(buffer(1:15).eq.'vtk_gas_density') then
        !
        ! Write out gas density in VTK format
        !
        amr_always_use_tree = .true.
        do_write_vtk_grid   = .true.
        do_write_vtk_gas_density = .true.
        gotit = .true.
     elseif(buffer(1:19).eq.'vtk_gas_temperature') then
        !
        ! Write out gas temperature in VTK format
        !
        amr_always_use_tree = .true.
        do_write_vtk_grid   = .true.
        do_write_vtk_gas_temperature = .true.
        gotit = .true.
     elseif((buffer(1:11).eq.'vtk_molspec').or.(buffer(1:12).eq.'vtk_chemspec')) then
        !
        ! Write out dust density in VTK format
        !
        amr_always_use_tree = .true.
        do_write_vtk_grid   = .true.
        do_write_vtk_chemspec = .true.
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) vtk_lines_ispec
        gotit = .true.
     elseif(buffer(1:12).eq.'vtk_levelpop') then
        !
        ! Write out level populations in VTK format
        !
        amr_always_use_tree = .true.
        do_write_vtk_grid   = .true.
        do_write_vtk_levelpop = .true.
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) vtk_lines_ispec
        call ggetarg(iarg,buffer,fromstdi)
        iarg = iarg+1
        read(buffer,*) vtk_lines_ilevel
        gotit = .true.
     elseif(buffer(1:12).eq.'vtk_velocity') then
        !
        ! Write out gas temperature in VTK format
        !
        amr_always_use_tree = .true.
        do_write_vtk_grid   = .true.
        do_write_vtk_velocity = .true.
        gotit = .true.
     else
        !
        ! Try if perhaps the userdef module knows what to do with
        ! this command line option
        !
        call userdef_commandline(buffer,numarg,iarg,fromstdi,gotit)
     endif
  end do
  !
  ! Some post-processing
  !
  if(flagzoom) then
     !
     ! Do some checks
     !
     if((zoom_x1.le.zoom_x0).or.(zoom_y1.le.zoom_y0)) then 
        write(stdo,*) 'ERROR: zoom_x1 <= zoom_x0 or zoom_y1 <= zoom_y0'
        stop
     endif
     !
     ! Convert the zoom locations x0,x1,y0,y1 into camera_zoomcenter_x,y and
     ! camera_image_halfsize_x,y, and also recompute the camera_image_nx,ny.
     !
     camera_zoomcenter_x     = 0.5d0 * ( zoom_x1 + zoom_x0 )
     camera_zoomcenter_y     = 0.5d0 * ( zoom_y1 + zoom_y0 )
     camera_image_halfsize_x = 0.5d0 * ( zoom_x1 - zoom_x0 )
     camera_image_halfsize_y = 0.5d0 * ( zoom_y1 - zoom_y0 )
     !
     ! Now adapt pixel geometry to try to keep pixels square, unless
     ! the user specifies with 'truepix' that the true number of 
     ! pixels is to be used. 
     !
     if(.not.truepix) then
        !
        ! Check if we should reduce one of the pixel numbers to get nicely
        ! square pixels. Slightly adapt the horizontal or vertical size to
        ! match the precise number of pixels. However, if the resulting
        ! number of pixels along the short dimension becomes <1, then we use
        ! non-rectangular pixels.
        !
        px = 2.d0 * camera_image_halfsize_x / (1.d0*camera_image_nx)
        py = 2.d0 * camera_image_halfsize_y / (1.d0*camera_image_ny)
        if(py.lt.px) then
           !
           ! Adopt pixel size px for y direction. If the 
           ! truezoom option is not set, then also slightly adapt
           ! the image size.
           !
           idum = floor( 2.d0 * camera_image_halfsize_y / px )
           if(idum.gt.1) then
              camera_image_ny         = idum
              if(.not.truezoom) then
                 camera_image_halfsize_y = 0.5d0 * ( camera_image_ny * px )
              endif
           else
              camera_image_ny         = 1
           endif
        elseif(px.lt.py) then
           !
           ! Adopt pixel size py for x direction. If the 
           ! truezoom option is not set, then also slightly adapt
           ! the image size.
           !
           idum = floor( 2.d0 * camera_image_halfsize_x / py )
           if(idum.gt.1) then
              camera_image_nx         = idum
              if(.not.truezoom) then
                 camera_image_halfsize_x = 0.5d0 * ( camera_image_nx * py )
              endif
           else
              camera_image_nx         = 1
           endif
        endif
     else
        if(truezoom) then
           write(stdo,*) 'WARNING: If truepix is set, truezoom is not used.'
           write(stdo,*) '         Ignoring truezoom.'
        endif
     endif
  else
     if(truepix) then
        write(stdo,*) 'WARNING: The option truepix is only doing something'
        write(stdo,*) '         if zoomau or zoompc has been specified.'
        write(stdo,*) '         So I ignore this option now.'
     endif
     if(truezoom) then
        write(stdo,*) 'WARNING: The option truezoom is only doing something'
        write(stdo,*) '         if zoomau or zoompc has been specified.'
        write(stdo,*) '         So I ignore this option now.'
        endif
  endif
  !
  ! Check for camera 
  !
  if((camera_range_nlam.gt.0).and.((camera_lambdamic.le.0.d0).or. &
       (camera_lambdamic1.le.0.d0))) then
     write(stdo,*) 'ERROR: Cannot specify nlam while not specifying lambdarange.'
     stop
  endif
  !
  ! Some postprocessing for the userdef?
  ! 
  call userdef_commandline_postprocessing()
  !
  ! Check if we received any commands at all
  !
  if(.not.gotit) then
     call write_banner()
     write(stdo,*) 'Nothing to do... Use command line options to generate action:'
     write(stdo,*) '  mctherm        : Do Monte Carlo simul of thermal radiation'
     write(stdo,*) '  mcmono         : Do Monte Carlo simul only for computing mean intensity'
     write(stdo,*) '  spectrum       : Make continuum spectrum'
     write(stdo,*) '  image          : Make continuum image'
     quit = .true.
  endif
end subroutine interpet_command_line_options



!-------------------------------------------------------------------------
!                  WRITE THE BANNER ON THE SCREEN
!-------------------------------------------------------------------------
subroutine write_banner()
  use constants_module
  implicit none
  write(stdo,*) ' '
  write(stdo,*) '================================================================'
  write(stdo,*) '     WELCOME TO RADMC-3D: A 3-D CONTINUUM AND LINE RT SOLVER    '
  write(stdo,*) '                                                                '
  write(stdo,*) '                         VERSION 2.0                            '
  write(stdo,*) '                                                                '
  write(stdo,*) '               (c) 2008-2023 Cornelis Dullemond                 '
  write(stdo,*) '                                                                '
  write(stdo,*) '      Please feel free to ask questions. Also please report     '
  write(stdo,*) '       bugs and/or suspicious behavior without hestitation.     '
  write(stdo,*) '     The reliability of this code depends on your vigilance!    '
  write(stdo,*) '                   dullemond@uni-heidelberg.de                  '
  write(stdo,*) '                                                                '
  write(stdo,*) '  To keep up-to-date with bug-alarms and bugfixes, follow the   '
  write(stdo,*) '                    RADMC-3D code on github                     '
  write(stdo,*) '            https://github.com/dullemond/radmc3d-2.0            '
  write(stdo,*) '                                                                '
  write(stdo,*) '             Please visit the RADMC-3D home page at             '
  write(stdo,*) ' http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/ '
  write(stdo,*) '================================================================'
  write(stdo,*) ' '
  call flush(stdo)
end subroutine write_banner


!-------------------------------------------------------------------------
!               WRITE INFO ABOUT THE INCLUDED PROCESSES
!-------------------------------------------------------------------------
subroutine write_message_rad_processes()
  use rtglobal_module
  use constants_module
  implicit none
  if(rt_incl_dust) then
     write(stdo,*) '  --> Including dust'
  else
     write(stdo,*) '      No dust included...'
  endif
  if(rt_incl_lines) then
     write(stdo,*) '  --> Including lines'
  else
     write(stdo,*) '      No lines included...'
  endif
  if(rt_incl_gascont) then
     write(stdo,*) '  --> Including gas continuum'
  else
     write(stdo,*) '      No gas continuum included...'
  endif
  if(rt_incl_userdef_srcalp) then
     write(stdo,*) '  --> Including user-defined emissivity and extinction coefficients'
  endif
  call flush(stdo)     
end subroutine write_message_rad_processes


