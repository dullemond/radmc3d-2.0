!=======================================================================
!                     MOLECULAR/ATOMIC LINE MODULE
!
! NOTE:
!   The molecular line data file format is that of the Leiden Atomic and
!   Molecular Database (LAMBDA), see
!   http://www.strw.leidenuniv.nl/~moldata/.  The paper belonging to that
!   data is: Schöier, F.L., van der Tak, F.F.S., van Dishoeck E.F., Black,
!   J.H. 2005, A&A 432, 369-379.
!
! The lines_mode tells how the level populations are computed. The 
! sign of the lines_mode has the following meaning:
!
!   lines_mode < 0    means on-the-fly computation (during the 
!                     ray-tracing, i.e. local; can be slow)
!   lines_mode > 0    means computing the populations beforehand, 
!                     storing them in a global array, and only then
!                     doing ray-tracing. This requires substantial
!                     memory.
! 
! The precise value of lines_mode has the following meaning: 
!
!   lines_mode = -1,+1    LTE populations
!   lines_mode = -2,+2    User-defined populations. RADMC-3D call the
!                         subroutine userdef_compute_levelpop() for that.
!                         Rest of line transfer is done by RADMC-3D.
!   lines_mode = -10      User-defined mode II: Line transfer is much
!                         more in the hands of the userdef_module.
!                         Using userdef_general_compute_levelpop().
!   lines_mode = 50       Level populations are read from file
!   lines_mode = -3,3     Large Velocity Gradient method (local non-LTE)
!   lines_mode = -4,4     Optically thin populations (local non-LTE)
!   lines_mode = 20..49   [FUTURE] These will be the non-local non-LTE
!                         modes such as lambda iteration, accelerated
!                         lambda iteration etc.
!
!=======================================================================
module lines_module
!$ use omp_lib
use amr_module
use rtglobal_module
use ioput_module
use mathroutines_module
use constants_module
!
! The maximum nr of levels and lines, used for the array sizes
!
integer :: lines_maxnrlevels=0
integer :: lines_maxnrlines=0
!
! The number of different atomic/molecular species included
!
integer :: lines_nr_species=0
!
! Some parameters that the main program can use for specifying the
! frequency array. If all are 0, then the frequencies array is read
! from a file.
!
integer :: lines_user_nrfreq=0
double precision :: lines_user_widthkms=0.d0
double precision :: lines_user_kms0=0.d0
integer :: lines_user_ispec=0
integer :: lines_user_iline=0
!
! Some defaults for these user defined things
! 
double precision :: lines_user_widthkms_default = 10.d0
integer :: lines_user_nrfreq_default = 40
!
! For each species, the species name. This is also used for finding
! the file name of that species, so please use names that can be used
! as parts of file names. Also use no spaces. 
!
character*80, allocatable :: lines_speciesname(:)
!
! For each species, the style of the input file:
!  leiden  = The input format of the LAMBDA database:
!            http://www.strw.leidenuniv.nl/~moldata/
!
character*80, allocatable :: lines_filestyle(:)
!
! For each species the data could either be the full level+line set
! or merely a line list. Since RADMC-3D has to treat them very 
! differently, here it is stored which one it is. If fullmolec==.true.
! then the full level diagram is included, if not only a line list
! is included.
!
logical, allocatable :: lines_species_fullmolec(:)
!
! For each species, the available number of levels and lines.
!
integer, allocatable :: lines_nrlevels(:)
integer, allocatable :: lines_nrlines(:)
!
! Same as above, but now only the number that are actually used
! (a subset of the total)
!
integer, allocatable :: lines_nrlevels_subset(:)
integer, allocatable :: lines_nrlines_subset(:)
integer :: lines_nrlevels_subset_max=0
integer :: lines_nrlines_subset_max=0
!
! The list of levels and lines actually used 
!
integer, allocatable :: lines_levels_subset(:,:)
integer, allocatable :: lines_lines_subset(:,:)
!
! Collisional Data
!
! Number of collision partners for each species
!
integer, allocatable :: lines_collisions_npartners(:)
integer, allocatable :: lines_collisions_npartners_used(:)
!
! Names of collision partners
!
integer, parameter :: lines_collisions_max_nr_partner_names = 20  ! More than you'll ever need, I think...
integer :: lines_collisions_nr_partner_names = 0
character*80 :: lines_collisions_partner_names(lines_collisions_max_nr_partner_names)
integer, allocatable :: lines_collisions_partners(:,:)
!
! Max number of collision partners per molecule for the collisional data
!
integer :: lines_collisions_npartmax = 0
!
! Max number of temperatures and transitions for the
! collisional data
!
integer :: lines_collisions_ntempmax = 0
integer :: lines_collisions_ntransmax = 0
!
! Number of Collisional Temperatures and Transitions
!
integer, allocatable :: lines_collisions_ntemp(:,:)
integer, allocatable :: lines_collisions_ntrans(:,:)
!
! Collisional Transition Temperatures 
!
double precision, allocatable :: lines_collisions_temperatures(:,:,:)
!
! Each Collision connects an upper and lower level
!
integer, allocatable :: lines_collisions_levelup(:,:,:)
integer, allocatable :: lines_collisions_leveldown(:,:,:)
!
! Collisional Transition Rates, for each temperature
!
double precision, allocatable :: lines_collisions_rates(:,:,:,:)
!
! Arrays for the statistical equilibrium solver
!
integer, allocatable :: lines_col_iup(:,:)
integer, allocatable :: lines_col_idown(:,:)
double precision, allocatable :: lines_col_ratedu(:,:)
double precision, allocatable :: lines_col_rateud(:,:)
!
! For each species, the molecular/atomic weight in units of
! the proton mass
!
double precision, allocatable :: lines_umass(:)
double precision :: lines_umass_min=0.d0
double precision :: lines_umass_max=0.d0
!
! For each level of each species, the energy
!
double precision, allocatable :: lines_level_energy_cminv(:,:)
!
! For each level of each species, the degeneracy parameter
!
double precision, allocatable :: lines_level_gdeg(:,:)
!
! For each level of each species, the name of the level
!
character*256, allocatable :: lines_level_name(:,:)
!
! For each line of each species, the A_{up,down}, B_{up,down} 
! and B_{down,up} Einstein coefficients.
!
double precision, allocatable :: lines_aud(:,:)
!
! For each line of each species, the frequency of the line in Hertz
!
double precision, allocatable :: lines_nu0(:,:)
!
! Each line connects an upper and a lower level. 
!
integer, allocatable :: lines_levelup(:,:)
integer, allocatable :: lines_leveldown(:,:)
!
! For linelist data: up and down energy and statistical weight
!
double precision, allocatable :: lines_linelist_eup(:,:)
double precision, allocatable :: lines_linelist_edown(:,:)
double precision, allocatable :: lines_linelist_gup(:,:)
double precision, allocatable :: lines_linelist_gdown(:,:)
!
! (For LVG: slowing parameter to avoid flipflopping)
! levelpopnew = speed*levelpopsolved+(1-speed)*levelpopold
! speed=1.0 is normal, speed=0.0 is no convergence,
!
double precision :: lines_solving_speed=0.5
!
! (For LVG) A counter for the average number of iterations
!
double precision :: lines_lvgesc_nriter, lines_lvgesc_nrsolves
!
! (For LTE, LVG, non-LTE internal computation): The partition function
! as a function of temperature
!
integer :: lines_partition_ntempint=1000
integer :: lines_partition_ntempmax=1000
integer, allocatable :: lines_partition_ntemp(:)
double precision :: lines_partition_temp0=0.1d0
double precision :: lines_partition_temp1=1d5
double precision, allocatable :: lines_partition_temperature(:,:)
double precision, allocatable :: lines_partition_function(:,:)
!
! Temporary array for partition sum (only for local internal usage)
!
double precision, allocatable :: lines_psum_local(:)
!
! Flags if (potential) line list data files contain the partition function
!
logical, allocatable :: lines_linelist_contains_parfunc(:)
!
! Margin factor for line width, for determining which lines contribute
! to any given frequency
!
double precision :: lines_widthmargin = 12.d0
!
! Maximum fraction of a line width the doppler shift between two neighboring
! cells may be before an alarm is set off.
! !REMOVED FOR NOW!
!
! double precision :: lines_maxdoppler = 0.3d0 
!
! For the line src and alpha
!
integer, allocatable :: active_lines(:,:)
integer, allocatable :: active_nrlines(:)
integer, allocatable :: active_levels(:,:)
integer, allocatable :: active_levels_subsetindex(:,:)
integer, allocatable :: active_nrlevels(:)
integer, allocatable :: active_dummy(:)
!
! The variable lines_profile sets the line profile function
! that should be used. The default (lines_profile = 0) selects
! a Gaussian function. If a Voigt function (lines_profile = 1)
! should be used, then the array lines_ray_lorentz_delta with
! the delta parameter of the Lorentz profile must be initialized.
! This is done in the userdef module in the routine
! userdef_compute_lorentz_delta. The lines_profile variable
! can be modified with the radmc3d.inp file.
! (Added by Thomas Peters 2011)
!
integer :: lines_profile=0
!
! Microwave background radiation, necessary for millimeter 
! line non-LTE level population calculcations.
!
double precision :: lines_tbg = 2.73d0
!
! Maximum nr of iterations for non-LTE radiative transfer
! (which includes LVG method as well)
!
integer :: lines_nonlte_maxiter = 100
double precision :: lines_nonlte_convcrit = 1d-2
!
! If lines_autosubset is set, then before each image is 
! produced the subset of levels is selected by checking which
! lines are in the wavelength domain of camera_frequencies().
! NOTE: This is only done for 1 <= lines_mode <= 9, i.e. when
! the level populations have to be stored globally, and when
! no fully-fledged non-LTE method is used (in which case anyway
! all levels would have to be stored globally). 
!
! NOTE: Actually this is a flag that is used mostly by the
!       camera module, but it must be located here for
!       accessibility to other modules.
!
logical :: lines_autosubset=.true.
!
contains


!-------------------------------------------------------------------
!                 CLEANUP THE LINE LAB DATA
!-------------------------------------------------------------------
subroutine lines_cleanup()
  if(allocated(lines_nrlevels)) deallocate(lines_nrlevels)
  if(allocated(lines_nrlines)) deallocate(lines_nrlines)
  if(allocated(lines_nrlevels_subset)) deallocate(lines_nrlevels_subset)
  if(allocated(lines_nrlines_subset)) deallocate(lines_nrlines_subset)
  if(allocated(lines_levels_subset)) deallocate(lines_levels_subset)
  if(allocated(lines_lines_subset)) deallocate(lines_lines_subset)
  if(allocated(lines_collisions_npartners)) deallocate(lines_collisions_npartners)
  if(allocated(lines_collisions_npartners_used)) deallocate(lines_collisions_npartners_used)
  if(allocated(lines_collisions_partners)) deallocate(lines_collisions_partners)
  if(allocated(lines_collisions_ntemp)) deallocate(lines_collisions_ntemp)
  if(allocated(lines_collisions_ntrans)) deallocate(lines_collisions_ntrans)
  if(allocated(lines_collisions_temperatures)) deallocate(lines_collisions_temperatures)
  if(allocated(lines_collisions_levelup)) deallocate(lines_collisions_levelup)
  if(allocated(lines_collisions_leveldown)) deallocate(lines_collisions_leveldown)
  if(allocated(lines_collisions_rates)) deallocate(lines_collisions_rates)
  if(allocated(lines_col_iup)) deallocate(lines_col_iup)
  if(allocated(lines_col_idown)) deallocate(lines_col_idown)
  if(allocated(lines_col_ratedu)) deallocate(lines_col_ratedu)
  if(allocated(lines_col_rateud)) deallocate(lines_col_rateud)
  if(allocated(lines_umass)) deallocate(lines_umass)
  if(allocated(lines_speciesname)) deallocate(lines_speciesname)
  if(allocated(lines_filestyle)) deallocate(lines_filestyle)
  if(allocated(lines_species_fullmolec)) deallocate(lines_species_fullmolec)
  if(allocated(lines_nu0)) deallocate(lines_nu0)
  if(allocated(lines_level_gdeg)) deallocate(lines_level_gdeg)
  if(allocated(lines_level_name)) deallocate(lines_level_name)
  if(allocated(lines_level_energy_cminv)) deallocate(lines_level_energy_cminv)
  if(allocated(lines_aud)) deallocate(lines_aud)
  if(allocated(lines_levelup)) deallocate(lines_levelup)
  if(allocated(lines_leveldown)) deallocate(lines_leveldown)
  if(allocated(lines_linelist_eup)) deallocate(lines_linelist_eup)
  if(allocated(lines_linelist_edown)) deallocate(lines_linelist_edown)
  if(allocated(lines_linelist_gup)) deallocate(lines_linelist_gup)
  if(allocated(lines_linelist_gdown)) deallocate(lines_linelist_gdown)
  if(allocated(lines_partition_temperature)) deallocate(lines_partition_temperature)
  if(allocated(lines_partition_function)) deallocate(lines_partition_function)
  if(allocated(lines_partition_ntemp)) deallocate(lines_partition_ntemp)
  if(allocated(lines_psum_local)) deallocate(lines_psum_local)
  if(allocated(lines_linelist_contains_parfunc)) deallocate(lines_linelist_contains_parfunc)
  if(allocated(active_lines)) deallocate(active_lines)
  if(allocated(active_nrlines)) deallocate(active_nrlines)
  if(allocated(active_levels)) deallocate(active_levels)
  if(allocated(active_levels_subsetindex)) deallocate(active_levels_subsetindex)
  if(allocated(active_nrlevels)) deallocate(active_nrlevels)
  if(allocated(active_dummy)) deallocate(active_dummy)
  lines_maxnrlevels = 0
  lines_maxnrlines  = 0
  lines_partition_ntempmax=0
  lines_nrlevels_subset_max = 0
  lines_nrlines_subset_max = 0
  lines_collisions_npartmax = 0
  lines_collisions_ntempmax = 0
  lines_collisions_ntransmax = 0
  lines_collisions_nr_partner_names = 0
  lines_maser_warning = .false.
end subroutine lines_cleanup

!----------------------------------------------------------------------------
!                       READ ALL FOR THE LINES
!
! You can call this routine also if no line data are present. It 
! will then simply reset all line stuff to zero. In fact, by checking
! whether lines_nr_species.gt.0 you know whether the line data are
! present or not. 
!----------------------------------------------------------------------------
subroutine read_lines_all(action)
  implicit none
  integer :: action
  integer :: ispec
  !
  ! Warning
  !
  if(lines_mode.lt.0) then
     write(stdo,*) 'Note: You use lines_mode < 0, meaning that RADMC-3D will calculate'
     write(stdo,*) '      the level populations on-the-fly. This saves memory, but can'
     write(stdo,*) '      be slow. Putting lines_mode to the same value but > 0 will'
     write(stdo,*) '      make RADMC-3D precalculate the populations before rendering'
     write(stdo,*) '      which can speed up the code considerably.'
  endif
  !
  ! Message
  !
  select case(lines_profile)
  case(0)
     write(stdo,*) "Using Gaussian line profile"
  case(1)
     write(stdo,*) "Using Voigt line profile"
  case default
     write(stdo,*) "Unknown line profile mode, aborting..."
     stop
  end select
  !
  ! Default
  ! (2016.07.23: Removed; Now done in read_linedata())
  !lines_pfunc_outofrange = .false.
  !lines_maser_warning = .false.
  !
  ! Read the line fundamental data such as energy levels, degeneration
  ! coefficients, Einstein A coefficients etc.
  !
  ! This has already been initialized in userdef_prep_model for userdef lines.
  !
  if (lines_mode.ne.-10) then
    call read_linedata(action)
  endif
  !
  ! If lines_mode.ne.50 then we must have the number densities of
  ! the molecules because we will compute the populations internally,
  ! either via local prescriptions or full non-LTE. In principle, if
  ! the full non-LTE is done and these populations are to be re-read
  ! by RADMC-3D for ray-tracing, in that case it is not necessary to
  ! read these, but it can't hurt either. 
  !
  if((lines_mode.ne.50).and.(lines_mode.ne.-10)) then 
     call read_molecules_numberdensities(action)
  endif
  !
  ! Always read the gas temperature, necessary, not only for the
  ! computation of the level populations, but also for the thermal
  ! line width.
  !
  call read_gas_temperature(action)
  lines_maxtempc = sqrt(2*kk*gastmax/(lines_umass_min))/cc
  ! 
  ! For non-LTE transfer, we need the number densities of the collisional
  ! partner, so read those.
  !
  if(((abs(lines_mode).ge.3).and.(abs(lines_mode).le.9)).or. &
     ((lines_mode.ge.20).and.(lines_mode.le.49))) then 
     call read_collpartner_numberdensities(action)
  endif
  !
  ! If abs(lines_mode).eq.1 then read the partition functions
  !
  if((lines_mode.eq.1).or.(lines_mode.eq.-1)) then
     call read_partition_function(action)
  endif
  !
  ! If we have non-local non-LTE transfer, we are not allowed to 
  ! use subsets of levels, because the storage of levels is then 
  ! part of the non-LTE transfer algorithm. Non-local non-LTE transfer
  ! is given by lines_mode in the 20 to 49 range. 
  !
  if((lines_mode.lt.50).and.(lines_mode.ge.20)) then
     do ispec=1,lines_nr_species
        if(lines_nrlevels(ispec).ne.lines_nrlevels_subset(ispec)) then
           write(stdo,*) 'ERROR: When doing *full* non-LTE line transfer you must'
           write(stdo,*) '       use the full level set of the molecule/atom.'
           write(stdo,*) '       You cannot use subsets. The reason is that the'
           write(stdo,*) '       global storage of the levels is then part'
           write(stdo,*) '       of the algorithm.'
           stop
        endif
     enddo
  endif
  !
  ! If lines_mode.ge.50, then read or compute the level populations if they
  ! are not yet in memory already. 
  !
  if(lines_mode.ge.50) then 
     call read_levelpop(action)
  endif
  !
  ! The velocity field always has to be read if not yet done
  !
  call read_velocityfield(action)
  !
  ! Microturbulence is only read if a corresponding file is available
  !
  call read_microturbulence(action)
  !
  ! Now compute the maximum width of an observed line, i.e. how far
  ! away from some wavelength should I still look for lines that might
  ! contribute to this wavelength... For the stochastic parts (the
  ! microturbulence and the temperature) we include a factor 
  ! lines_widthmargin to see how far in the line wing we want to look.
  !
  lines_maxshift = lines_widthmargin * ( lines_maxtempc + lines_maxturbc ) + &
                   lines_maxveloc
  !
  ! If lines mode is 3 (LVG), then check if a file called escprob_lengthscale.inp
  ! is present. If so, then read it.
  ! 
  if(abs(lines_mode).eq.3) then
     call read_escprob_lengthscale(action)
  endif
  !
  ! Since we know the grid, we also know the largest length of a single
  ! ray. This means that we can already set up the 1-D formal transfer
  ! along a ray stuff.
  !
  ! ############## FUTURE: FOR NVIDIA GRAPHICS CARD #########
  !!!call lines_ray1d_init_raytrace(action)
  !
  call lines_serial_init_raytrace(action)
  !
  ! Make a message
  !
  if(lines_nr_species.gt.0) then
     select case(lines_mode)
     case(1)
        write(stdo,*) 'Line transfer method: LTE (populations precalculated)'
     case(-1)
        write(stdo,*) 'Line transfer method: LTE (on-the-fly calculation)'
     case(2)
        write(stdo,*) 'Line transfer method: User-defined (populations precalculated)'
     case(-2)
        write(stdo,*) 'Line transfer method: User-defined (on-the-fly calculation)'
     case(3)
        write(stdo,*) 'Line transfer method: LVG (Sobolev) (populations precalculated)'
        if(allocated(lines_escprob_lengthscale)) then
           write(stdo,*) '                      (including escape probability length scale)'
        endif
     case(-3)
        write(stdo,*) 'Line transfer method: LVG (Sobolev) (on-the-fly calculation)'
        if(allocated(lines_escprob_lengthscale)) then
           write(stdo,*) '                      (including escape probability length scale)'
        endif
     case(4)
        write(stdo,*) 'Line transfer method: Opt Thin NLTE (populations precalculated)'
     case(-4)
        write(stdo,*) 'Line transfer method: Opt Thin NLTE (on-the-fly calculation)'
     case(-10)
        write(stdo,*) 'Line transfer method: User-defined Plus (on-the-fly calculation)'
     case(50)
        write(stdo,*) 'Line transfer method: Populations read from file'
     case default
        write(stdo,*) 'Unknown line transfer method: ',lines_mode
     end select
  endif
  !
end subroutine read_lines_all



!-------------------------------------------------------------------
!                    READ THE MAIN LINE DATA FILE
!-------------------------------------------------------------------
subroutine read_linedata(action)
  implicit none
  integer :: iformat,ierror,idum,ispec,idum1,idum2,idum3,ierr,ilevel
  integer :: action,ipart,iname
  character*80 :: ilinestring,jlinestring,name,str1,str2,str3
  logical :: fex,thereis_linelist
  integer, allocatable :: subset_mode(:)
  integer, allocatable :: subset_temp_ilevels(:,:)
  logical :: manual_subset_selection
  !
  ! For technical reasons, if you wish to specify a subset of levels
  ! one-by-one (subset_mode==1), you can specify only up to 
  ! TEMP_MAXLEVELS of those levels. I chose TEMP_MAXLEVELS=100, which
  ! I think is more than you'll ever need. But if, for some reason,
  ! you want to select more (and if subset_mode==0 does not work for
  ! your purpose), then you'll have to change this number here and
  ! recompile the code.
  !
  integer :: TEMP_MAXLEVELS
  parameter(TEMP_MAXLEVELS=100)
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(lines_nrlevels)) return
  endif
  !
  ! Clean up any remaining stuff
  !
  call lines_cleanup()
  !
  ! Set some default constants
  !
  !lines_widthmargin = 3.d0   ! Default should be from the main routine.
  !lines_warn_lineleap = .false.
  thereis_linelist = .false.
  manual_subset_selection = .false.
  lines_pfunc_outofrange = .false.
  lines_maser_warning = .false.
  !
  ! Open the lines.inp
  !
  inquire(file='lines.inp',exist=fex)
  if(.not.fex) then
     write(stdo,*) 'ERROR in line module: Could not find file lines.inp'
     stop
  endif
  !
  ! Message
  !
  write(stdo,*) 'Reading line data...'
  call flush(stdo)
  !
  ! Open file
  !
  open(unit=3,file='lines.inp',status='unknown')
  !
  ! Format number
  !
  read(3,*) iformat
  !
  ! Read how many species of molecules/atoms we want to read
  !
  read(3,*) lines_nr_species
  !
  ! Read the names of the input files and the max nr of 
  ! levels to take into account
  !
  ! First check if no strange number was entered
  !
  if((lines_nr_species.lt.1).or.(lines_nr_species.gt.100)) then
     write(stdo,*) 'ERROR in line module: lines_nr_species out of range'
     stop
  endif
  !
  ! Allocate some arrays
  !
  allocate(lines_species_fullmolec(1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_species_fullmolec(:).'
     stop 
  endif
  allocate(lines_nrlevels(1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_nrlevels(:).'
     stop 
  endif
  lines_nrlevels(:) = 0
  allocate(lines_nrlines(1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_nrlines(:).'
     stop 
  endif
  lines_nrlines(:) = 0
  allocate(lines_nrlevels_subset(1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_nrlevels_subset(:).'
     stop 
  endif
  lines_nrlevels_subset(:) = 0
  allocate(lines_nrlines_subset(1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_nrlines_subset(:).'
     stop 
  endif
  lines_nrlines_subset(:) = 0
  allocate(lines_umass(1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_umass(:).'
     stop 
  endif
  lines_umass(:) = 0.d0
  allocate(lines_speciesname(1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_speciesname(:).'
     stop 
  endif
  allocate(lines_filestyle(1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_filestyle(:).'
     stop 
  endif
  allocate(lines_partition_ntemp(1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate the lines_partition_ntemp() array'
     stop
  endif
  lines_partition_ntemp(:) = 0
  allocate(lines_psum_local(1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate the lines_psum_local() array'
     stop
  endif
  allocate(subset_mode(1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      subset_mode(:).'
     stop 
  endif
  subset_mode(:) = 0
  allocate(lines_collisions_npartners(1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_collisions_npartners(:).'
     stop 
  endif
  lines_collisions_npartners(:) = 0
  allocate(lines_collisions_npartners_used(1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_collisions_npartners_used(:).'
     stop 
  endif
  lines_collisions_npartners_used(:) = 0
  allocate(lines_collisions_partners(1:lines_collisions_max_nr_partner_names,1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_collisions_partners(:,:).'
     stop 
  endif
  lines_collisions_partners(:,:) = 0
  allocate(lines_linelist_contains_parfunc(1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_linelist_contains_parfunc(:,:).'
     stop 
  endif
  lines_linelist_contains_parfunc(:) = .false.
  allocate(subset_temp_ilevels(1:TEMP_MAXLEVELS,1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      subset_temp_ilevels(:).'
     stop 
  endif
  subset_temp_ilevels(:,:) = 0
  !
  ! Now read the species names, the input file format, the 
  ! way to select the nr of levels (idum1), the max nr of 
  ! levels (idum2), and how many collision partner names
  ! are to be listed (latter only for iformat.ge.2). 
  ! The idum2 number means 'unspecified' if 0. 
  ! If idum1 is 0 then we select the levels only by their number limit.
  ! Example (iff iformat.ge.2):
  !    co    leiden    0   0   0    ==> Read molecule_co.inp from
  !                                     the Leiden database and
  !                                     use all levels and lines.
  !    co    leiden    0   10  0    ==> Read molecule_co.inp from
  !                                     the Leiden database and
  !                                     use only levels 1 to 10
  !    co    leiden    1   10  0    ==> Read molecule_co.inp from
  !                                     the Leiden database and
  !                                     select 10 levels (which are
  !                                     specified further down).
  !    co    leiden    0   0   2    ==> Read molecule_co.inp from
  !    p-h2                             the Leiden database, and
  !    o-h2                             if we do non-LTE transfer,
  !                                     then use the files 
  !                                     numberdens_p-h2.inp and
  !                                     numberdens_o-h2.inp
  !                                     as the collision partner
  !                                     number densities corresponding
  !                                     to the two rate tables in
  !                                     the molecule_co.inp file.
  !    co    leiden    0  10   2    ==> Same as before, but if we
  !    p-h2                             do NOT do full-non-LTE transfer,
  !    o-h2                             then use only 10 levels
  !                                     (with non-LTE transfer,
  !                                     we always use all levels).
  !    co    linelist  0  0    0        Read linelist_co.inp, which
  !                                     is the RADMC-3D standard
  !                                     line list format. Can 
  !                                     only be used for LTE.
  !           
  do ispec=1,lines_nr_species
     !
     ! Read the info
     !
     if(iformat.eq.1) then
        read(3,*) ilinestring,jlinestring,idum1,idum2
        idum3 = 0
     else
        read(3,*) ilinestring,jlinestring,idum1,idum2,idum3
     endif
     lines_speciesname(ispec)      = ilinestring
     lines_filestyle(ispec)        = jlinestring
     subset_mode(ispec)            = idum1
     lines_nrlevels_subset(ispec)  = idum2
     lines_collisions_npartners_used(ispec) = idum3
     !
     ! If one or more species have manual selection of the 
     ! level subset, put this flag up:
     !
     if((idum1.ne.0).or.(idum2.ne.0)) then
        manual_subset_selection = .true.
     endif
     !
     ! If subset_mode is 1, then read subset levels
     !
     if(subset_mode(ispec).eq.1) then
        if(lines_nrlevels_subset(ispec).gt.TEMP_MAXLEVELS) then
           write(stdo,*) 'ERROR while making subset of levels:'
           write(stdo,*) 'For technical reasons, if you wish to specify a subset of levels'
           write(stdo,*) 'one-by-one (subset_mode==1), you can specify only up to '
           write(stdo,*) 'TEMP_MAXLEVELS of those levels. I chose TEMP_MAXLEVELS=100, which'
           write(stdo,*) 'I think is more than you will ever need. But since you get this'
           write(stdo,*) 'error message, you apparently want to select more.'
           write(stdo,*) 'You will have to edit the lines_module.f90 file and change'
           write(stdo,*) 'the value of TEMP_MAXLEVELS and recompile the code.'
           stop
        endif
        read(3,*) (subset_temp_ilevels(ilevel,ispec), &
             ilevel=1,lines_nrlevels_subset(ispec))
     endif
     !
     ! Check if this file is present, and if so, check the nr of levels,
     ! and in fact if it contains level information at all (or is a linelist
     ! instead)
     !
     call read_molecule(ispec,.false.,ierror)
     if(ierror.ne.0) then
        write(stdo,*) 'ERROR in line module: somehow cannot read species ',ispec
        stop
     endif
     !
     ! If necessary, read the names of the collision partners for this
     ! molecule, and check if we already have then from previous molecules.
     !
     ! NOTE: These files (e.g. numberdens_p-h2.inp) are only read
     !       if non-LTE transfer is done.
     !
     if(lines_species_fullmolec(ispec)) then
        if(lines_collisions_npartners_used(ispec).gt.lines_collisions_npartners(ispec)) then
           call integer_to_string(ispec,str1)
           call integer_to_string(lines_collisions_npartners(ispec),str2)
           call integer_to_string(lines_collisions_npartners_used(ispec),str3)
           write(stdo,*) 'ERROR: In the molecular data file for molecule '//trim(str1)
           write(stdo,*) '       there are collisional transition rates '
           write(stdo,*) '       for '//trim(str2)//' collision partners.'
           write(stdo,*) '       However, in lines.inp for this line there are '//trim(str3)
           write(stdo,*) '       names of the collision partners specified.'
           write(stdo,*) '       This is inconsistent.'
           stop
        endif
        if(lines_collisions_npartners_used(ispec).lt.lines_collisions_npartners(ispec)) then
           if(lines_collisions_npartners_used(ispec).gt.0) then
              call integer_to_string(ispec,str1)
              call integer_to_string(lines_collisions_npartners(ispec),str2)
              call integer_to_string(lines_collisions_npartners_used(ispec),str3)
              write(stdo,*) 'WARNING: In the molecular data file for molecule '//trim(str1)
              write(stdo,*) '         there are collisional transition rates specified'
              write(stdo,*) '         for '//trim(str2)//' collision partners.'
              write(stdo,*) '         However, in lines.inp for this line only '//trim(str3)
              write(stdo,*) '         names of the collision partners are given.'
              write(stdo,*) '         We will therefore use not all possible partners.'
           else
              call integer_to_string(ispec,str1)
              call integer_to_string(lines_collisions_npartners(ispec),str2)
              call integer_to_string(lines_collisions_npartners_used(ispec),str3)
              write(stdo,*) 'NOTE: In lines.inp for molecule '//trim(str1)// &
                   ' no collision partners specified, therefore no non-LTE possible.'
           endif
        endif
        do ipart=1,lines_collisions_npartners_used(ispec)
           !
           ! Read name of collision partner
           ! 
           read(3,*) name
           !
           ! Check if we already have used this name before
           ! (for another molecule)
           !
           do iname=1,lines_collisions_nr_partner_names
              if(trim(name).eq.trim(lines_collisions_partner_names(iname))) then
                 lines_collisions_partners(ipart,ispec) = iname
                 goto 200
              endif
           enddo
           !
           ! This collision partner is yet used by other
           ! molecules. So assign it.
           !
           lines_collisions_nr_partner_names = lines_collisions_nr_partner_names + 1
           if(lines_collisions_nr_partner_names.gt.lines_collisions_max_nr_partner_names) then
              write(stdo,*) 'ERROR: Maximum number of different collision partners'
              write(stdo,*) '       exceeded. Please recompile RADMC-3D with a larger'
              write(stdo,*) '       lines_collisions_max_nr_partner_names in the'
              write(stdo,*) '       lines_module.f90. But first check if you are sure'
              write(stdo,*) '       that so many collision partners are needed. Seems strange.'
              stop
           endif
           lines_collisions_partner_names(lines_collisions_nr_partner_names) = trim(name)
           lines_collisions_partners(ipart,ispec) = lines_collisions_nr_partner_names
200        continue
        enddo
        !
        ! Check if the number of subset levels does not exceed the number of levels
        !
        ! The default "co    leiden    0   0   0 " (in which the subset is in
        ! fact the full set of levels) is implemented here.
        !
        if((lines_nrlevels_subset(ispec).le.0).or.                        &
           (lines_nrlevels_subset(ispec).gt.lines_nrlevels(ispec))) then
           lines_nrlevels_subset(ispec) = lines_nrlevels(ispec)
        endif
     else
        !
        ! This is apparently a molecule for which a line list is specified
        ! instead of the full molecular data. Signal for later that we have to
        ! allocate some more stuff.
        !
        thereis_linelist = .true.
        !
     endif
     !
  enddo
  !
  ! If one or more molecules had manual subset selection, but the
  ! lines_autosubset is set to .true., we must abort, because
  ! autosubset selection is only possible if none of the molecules
  ! have manual subset selection. In other words: if you want to
  ! manual select the subset, then you must put autosubset to false.
  !
  if(manual_subset_selection.and.lines_autosubset) then
     write(stdo,*) 'ERROR: In the lines.inp you apparently select a subset of'
     write(stdo,*) '       levels by hand (e.g. "co  leiden  0 10 0"), but '
     write(stdo,*) '       the lines_autosubset flag is .true., meaning'
     write(stdo,*) '       RADMC-3D plans to automatically select the subset.'
     write(stdo,*) '       These are mutually exclusive. '
     write(stdo,*) '       If you want to manually select levels for the subset'
     write(stdo,*) '       (only for experts!) then you must switch off the'
     write(stdo,*) '       autosubset - either with the "noautosubset" option'
     write(stdo,*) '       on the command line, or by adding the line:'
     write(stdo,*) '         lines_autosubset = 0'
     write(stdo,*) '       to the radmc3d.inp file. If you want to go safe, then'
     write(stdo,*) '       do not select a subset manually. You can do this by'
     write(stdo,*) '       setting the first two integers after "leiden" in the'
     write(stdo,*) '       lines.inp to zero, e.g.:'
     write(stdo,*) '         co    leiden    0   0   0 '
     stop
  endif
  !
  ! Find the maximum number of levels and lines available for each
  ! species
  !
  lines_maxnrlevels = 0
  lines_maxnrlines = 0
  do ispec=1,lines_nr_species
     lines_maxnrlevels = max(lines_maxnrlevels,lines_nrlevels(ispec))
     lines_maxnrlines = max(lines_maxnrlines,lines_nrlines(ispec))
  enddo
  !
  ! Find the maximum number of levels used for each species
  !
  lines_nrlevels_subset_max = 0
  do ispec=1,lines_nr_species
     lines_nrlevels_subset_max = max(lines_nrlevels_subset_max,     &
                                     lines_nrlevels_subset(ispec))
  enddo
  !
  ! Allocate and set the subset arrays. These arrays contain the indices of
  ! the levels/lines that are used. Please keep them in ascending
  ! order!!
  !
  ! Note: if all molecular data files are linelists, then we do not have
  !       to do this. Hence the if(lines_maxnrlevels.gt.0) statement.
  !
  if(lines_maxnrlevels.gt.0) then
     !
     ! Allocate
     !
     allocate(lines_levels_subset(1:lines_maxnrlevels,lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      lines_levels_subset(:,:).'
        stop 
     endif
     allocate(lines_lines_subset(1:lines_maxnrlines,lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      lines_lines_subset(:,:).'
        stop 
     endif
     !
     ! Now set the subsets
     !
     do ispec=1,lines_nr_species
        if(subset_mode(ispec).eq.0) then
           !
           ! Just select the first nrlevels(ispec) levels
           !
           do ilevel=1,lines_nrlevels_subset(ispec)
              lines_levels_subset(ilevel,ispec) = ilevel
           enddo
        elseif(subset_mode(ispec).eq.1) then
           !
           ! Read a list of levels
           ! 
           do ilevel=1,lines_nrlevels_subset(ispec)
              lines_levels_subset(ilevel,ispec) =      & 
                   subset_temp_ilevels(ilevel,ispec)
           enddo
           !
           ! Check if all these are within range
           !
           do ilevel=1,lines_nrlevels_subset(ispec)
              if((lines_levels_subset(ilevel,ispec).lt.1).or.    &
                   (lines_levels_subset(ilevel,ispec).gt.lines_nrlevels(ispec))) then
                 write(stdo,*) 'ERROR in line module: subset indices for levels out of range'
                 stop
              endif
           enddo
           !
           ! Check if they are in ascending order
           !
           do ilevel=2,lines_nrlevels_subset(ispec)
              if(lines_levels_subset(ilevel,ispec).le.lines_levels_subset(ilevel-1,ispec)) then
                 write(stdo,*) 'ERROR: Subset of levels for species ',ispec
                 write(stdo,*) '    is not in ascending order...'
                 stop
              endif
           enddo
        else
           write(stdo,*) 'ERROR in line module: do not know subset_mode = ', &
                subset_mode(ispec)
           stop
        endif
     enddo
  endif
  !
  ! Close this file
  !
  close(3)
  !
  ! Now allocate the rest of the arrays
  !
  ! First the general line-related ones
  !
  allocate(lines_nu0(1:lines_maxnrlines,1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_nu0(:).'
     stop 
  endif
  allocate(lines_aud(1:lines_maxnrlines,1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_aud(:).'
     stop 
  endif
  !
  ! Then line-related stuff that is only used for molecules for which
  ! the full level diagram is given.
  !
  if(lines_maxnrlevels.gt.0) then
     allocate(lines_levelup(1:lines_maxnrlines,1:lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      lines_levelup(:).'
        stop 
     endif
     allocate(lines_leveldown(1:lines_maxnrlines,1:lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      lines_leveldown(:).'
        stop 
     endif
  endif
  !
  ! Then line-related stuff that is only used for molecules for which
  ! linelist data is available
  !
  if(thereis_linelist) then
     allocate(lines_linelist_gup(1:lines_maxnrlines,1:lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      lines_linelist_gup(:).'
        stop 
     endif
     allocate(lines_linelist_gdown(1:lines_maxnrlines,1:lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      lines_linelist_gdown(:).'
        stop 
     endif
     allocate(lines_linelist_eup(1:lines_maxnrlines,1:lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      lines_linelist_eup(:).'
        stop 
     endif
     allocate(lines_linelist_edown(1:lines_maxnrlines,1:lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      lines_linelist_edown(:).'
        stop 
     endif
  endif
  !
  ! Then the level-related ones. These only for the case when there
  ! are some molecular data files that contain this information
  !
  if(lines_maxnrlevels.gt.0) then
     !
     ! The general ones
     !
     allocate(lines_level_energy_cminv(1:lines_maxnrlevels,1:lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      lines_level_energy_cminv(:).'
        stop 
     endif
     allocate(lines_level_gdeg(1:lines_maxnrlevels,1:lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      lines_level_gdeg(:).'
        stop 
     endif
     allocate(lines_level_name(1:lines_maxnrlevels,1:lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      lines_level_name(:).'
        stop 
     endif
     !
     ! The ones for the collisional rates, if they are available
     !
     if(lines_collisions_npartmax.gt.0) then
        allocate(lines_collisions_ntemp(1:lines_collisions_npartmax,1:lines_nr_species),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in lines module: Could not allocate '
           write(stdo,*) '      lines_collisions_ntemp(:).'
           stop 
        endif
        lines_collisions_ntemp(:,:) = 0
        allocate(lines_collisions_ntrans(1:lines_collisions_npartmax,1:lines_nr_species),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in lines module: Could not allocate '
           write(stdo,*) '      lines_collisions_ntrans(:).'
           stop 
        endif
        lines_collisions_ntrans(:,:) = 0
        allocate(lines_collisions_temperatures(1:lines_collisions_ntempmax,&
                 1:lines_collisions_npartmax,1:lines_nr_species),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in lines module: Could not allocate '
           write(stdo,*) '      lines_collisions_temperatures(:).'
           stop 
        endif
        allocate(lines_collisions_levelup(1:lines_collisions_ntransmax,&
                 1:lines_collisions_npartmax,1:lines_nr_species),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in lines module: Could not allocate '
           write(stdo,*) '      lines_collisions_levelup(:).'
           stop 
        endif
        allocate(lines_collisions_leveldown(1:lines_collisions_ntransmax,&
                 1:lines_collisions_npartmax,1:lines_nr_species),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in lines module: Could not allocate '
           write(stdo,*) '      lines_collisions_leveldown(:).'
           stop 
        endif
        allocate(lines_collisions_rates(1:lines_collisions_ntempmax,&
                 1:lines_collisions_ntransmax,&
                 1:lines_collisions_npartmax,1:lines_nr_species),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in lines module: Could not allocate '
           write(stdo,*) '      lines_collisions_rates(:).'
           stop 
        endif
        allocate(lines_col_iup(1:lines_collisions_ntransmax,&
                 1:lines_collisions_npartmax),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in lines module: Could not allocate '
           write(stdo,*) '      lines_col_iup(:).'
           stop 
        endif
        allocate(lines_col_idown(1:lines_collisions_ntransmax,&
                 1:lines_collisions_npartmax),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in lines module: Could not allocate '
           write(stdo,*) '      lines_col_idown(:).'
           stop 
        endif
        allocate(lines_col_ratedu(1:lines_collisions_ntransmax,&
                 1:lines_collisions_npartmax),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in lines module: Could not allocate '
           write(stdo,*) '      lines_col_ratedu(:).'
           stop 
        endif
        allocate(lines_col_rateud(1:lines_collisions_ntransmax,&
                 1:lines_collisions_npartmax),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in lines module: Could not allocate '
           write(stdo,*) '      lines_col_rateud(:).'
           stop 
        endif
     endif
  endif
  !
  ! Now read each molecule/atom data
  !
  do ispec=1,lines_nr_species
     call read_molecule(ispec,.true.,ierror)
  enddo
  !
  ! Compute the total number of subset lines
  !
  lines_nrlines_subset_max=0
  if(lines_maxnrlevels.gt.0) then
     do ispec=1,lines_nr_species
        lines_nrlines_subset_max = max(lines_nrlines_subset_max,    &
                                       lines_nrlines_subset(ispec))
     enddo
  endif
  !
  ! Allocate a local arrays for lines_serial_addto_jnu_alpnu() or
  ! lines_ray1d_addto_jnu_alpnu()
  !
  allocate(active_lines(lines_nrlines_subset_max,lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in line module: array active_lines() not allocatable'
     stop
  endif
  allocate(active_nrlines(lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in line module: array active_nrlines() not allocatable'
     stop
  endif
  !
  ! For the levels only if molecules exist with level information
  !
  if(lines_maxnrlevels.gt.0) then
     allocate(active_levels(lines_nrlevels_subset_max,lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in line module: array active_levels() not allocatable'
        stop
     endif
     allocate(active_levels_subsetindex(lines_nrlevels_subset_max,lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in line module: array active_levels_subsetindex() not allocatable'
        stop
     endif
     allocate(active_nrlevels(lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in line module: array active_nrlevels() not allocatable'
        stop
     endif
     allocate(active_dummy(lines_maxnrlevels),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in line module: array active_dummy() not allocatable'
        stop
     endif
  endif
  !
  ! Pre-compute some stuff
  !
  lines_umass_min = 1.d99
  lines_umass_max = 0.d0
  do ispec=1,lines_nr_species
     if(lines_umass(ispec).lt.lines_umass_min) lines_umass_min = lines_umass(ispec)
     if(lines_umass(ispec).gt.lines_umass_max) lines_umass_max = lines_umass(ispec)
  enddo
  !
  ! Deallocate temporary array
  !
  if(allocated(subset_mode)) deallocate(subset_mode)
  if(allocated(subset_temp_ilevels)) deallocate(subset_temp_ilevels)
  !
end subroutine read_linedata


!-------------------------------------------------------------------
! READ A MOLECULAR/ATOMIC DATA FILE FROM ANY KIND OF SUPPORTED FILE
!-------------------------------------------------------------------
subroutine read_molecule(ispec,readdata,ierror)
  implicit none
  character*80 :: linestring
  integer :: ispec,ierror,iline,ilinsub
  logical :: readdata
  character*160 :: species
  !
  ! Call the reading routine for this specific database type
  !
  linestring = lines_filestyle(ispec)
  if(linestring(1:6).eq.'leiden') then
     lines_species_fullmolec(ispec) = .true.
     call read_molecule_leiden(ispec,readdata,ierror)
     if(ierror.ne.0) then
        write(stdo,*) 'ERROR during line reading: failed. Aborting.'
        stop
     endif
  elseif(linestring(1:8).eq.'linelist') then
     lines_species_fullmolec(ispec) = .false.
     call read_linelist(ispec,readdata,ierror)
  else
     write(stdo,*) 'ERROR: Do not know file style ',linestring
     stop
  endif
  !
  ! Print out the line list, if requested
  !
  if(readdata.and.lines_make_linelist.and.lines_species_fullmolec(ispec)) then
     species  = lines_speciesname(ispec)
     write(stdo,*) 'Line list for molecule ',species(1:len_trim(species)),':'
     ilinsub=1
     do iline=1,lines_nrlines(ispec)
        if(lines_lines_subset(ilinsub,ispec).eq.iline) then
           write(stdo,*) 'Transition ',iline,': lambda = ',  &
                1d4*cc/lines_nu0(iline,ispec),' micron'
           ilinsub = ilinsub + 1
        else
           write(stdo,*) 'Transition ',iline,': lambda = ',  &
                1d4*cc/lines_nu0(iline,ispec),' micron [not included in model]'
        endif
     enddo
  endif
  !
end subroutine read_molecule


!-------------------------------------------------------------------
!    READ A MOLECULAR/ATOMIC DATA FILE FROM THE LEIDEN DATABASE
!-------------------------------------------------------------------
subroutine read_molecule_leiden(ispec,readdata,ierror)
  implicit none
  integer :: ispec,ierror,ierr
  logical :: readdata,fex,flag1,flag2
  character*160 filename,header,species
  character*80 :: str1
  integer :: nlevels,nlines,ilevel,iline,iup,idown,ilinesub,ilevelsub
  integer :: ilevel1,ilevel2,itemp,itrans,ipartner,ncolltransmax,ntempmax
  integer, allocatable :: line_iup(:),line_idown(:)
  double precision, allocatable :: line_aud(:),line_nu0(:)
  double precision, allocatable :: level_energy_cminv(:),level_gdeg(:)
  integer :: col_npartners
  integer, allocatable :: col_ntemp(:),col_ntrans(:)
  integer, allocatable :: col_iup(:,:),col_idown(:,:)
  double precision, allocatable :: col_temperatures(:,:),col_rates(:,:,:)
  character*256, allocatable :: level_name(:)
  character*256 :: name
  double precision :: e,g,dum1,dum2,dum3,aud,nu0,nucheck,error_linefreq 
  !
  ! Some check
  !
  if((.not.allocated(lines_aud)).and.readdata) then
     write(stdo,*) 'ERROR in line module: stuff not allocated...'
     stop
  endif
  if(ispec.gt.lines_nr_species) then
     write(stdo,*) 'ERROR in line module: ispec out of range...'
     stop
  endif
  !
  ! Initial setting
  !
  ierror = 0
  !
  ! Get species name
  !
  species  = lines_speciesname(ispec)
  !
  ! Message
  !
  if(readdata) then
     write(stdo,*) 'Reading line data of molecule ',species(1:len_trim(species)), &
          ' ...'
     call flush(stdo)
  endif
  !
  ! Find the file
  !
  filename = 'molecule_'//species(1:len_trim(species))//'.inp'
  inquire(file=filename,exist=fex)
  if(.not.fex) then
     ierror = 1
     return
  endif
  !
  ! Open the file
  !
  open(unit=1,file=filename)
  !
  ! Read the first dummy strings
  !      
  read(1,*) header
  read(1,*) header     ! The name of the molecule
  read(1,*) header
  !
  ! Read the molecular weight
  !
  read(1,*) lines_umass(ispec)
  lines_umass(ispec) = lines_umass(ispec)*mp
  !
  ! Read header
  !
  read(1,*) header
  !
  ! Read the number of levels
  !
  read(1,*) nlevels
  !
  ! Allocate the level energy and g array
  !
  allocate(level_energy_cminv(1:nlevels),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in line module: could not allocate temporary array level_energy_cminv'
     stop
  endif
  allocate(level_gdeg(1:nlevels),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in line module: could not allocate temporary array level_gdeg'
     stop
  endif
  allocate(level_name(1:nlevels),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in line module: could not allocate temporary array level_name'
     stop
  endif
  !
  ! Read header
  !
  read(1,*) header
  !
  ! Now do a loop over the levels and read their energies in cm^-1 and
  ! the degeneration parameters g.
  !
  ! NOTE: I have verified for at least one file of the Leiden LAMBDA
  !       database that the right light speed is used for the link between 
  !       the energies and the line positions: c=2.99792458d10 cm/s.
  !       Note that this value is *by definition* the light speed (it
  !       defines the length of the meter).
  !
  do ilevel=1,nlevels
     !
     ! Read the energy and degeneration
     !
     read(1,*) dum1,e,g,name
     !
     ! Store
     !
     level_energy_cminv(ilevel) = e
     level_gdeg(ilevel)         = g
     level_name(ilevel)         = name
  enddo
  !
  ! Read header
  !
  read(1,*) header
  !
  ! Now read the number of transitions
  !
  read(1,*) nlines
  !
  ! Allocate the line info
  !
  allocate(line_iup(1:nlines),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in line module: could not allocate temporary array line_iup'
     stop
  endif
  allocate(line_idown(1:nlines),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in line module: could not allocate temporary array line_idown'
     stop
  endif
  allocate(line_aud(1:nlines),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in line module: could not allocate temporary array line_aud'
     stop
  endif
  allocate(line_nu0(1:nlines),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in line module: could not allocate temporary array line_nu0'
     stop
  endif
  !
  ! Read header
  !
  read(1,*) header
  !
  ! Now loop over the lines
  !
  do iline=1,nlines
     !
     ! Read
     !
     read(1,*) dum1,iup,idown,aud,nu0,dum3
     !
     ! Store
     !
     line_iup(iline)   = iup
     line_idown(iline) = idown
     line_aud(iline)   = aud
     line_nu0(iline)   = nu0*1d9   ! The frequency is in GHz
  enddo
  !  
  ! Begin Read in Collisional Rates.
  ! 
  ! First reset stuff
  !
  col_npartners = 0
  !
  ! Read header
  !
  read(1,*,end=200) header
  !
  ! Now read the number of Collision Partners
  !
  read(1,*,end=200) col_npartners
  !
  ! If col_npartners=0, then no collision rates are present
  !
  if(col_npartners.gt.0) then
     !
     ! First do a dummy read to get the dimensions, then redo
     ! the entire file reading to get the data. 
     !
     ncolltransmax = 0
     ntempmax      = 0
     allocate(col_ntemp(col_npartners),col_ntrans(col_npartners))
     col_ntemp(:)  = 0
     col_ntrans(:) = 0
     !
     ! Loop over collision partners
     !
     do ipartner=1,col_npartners
        !
        ! Read headers
        !
        read(1,*) header
        read(1,*) header    ! Here the name of the coll partner is hidden...
        read(1,*) header
        !
        ! Now read the number of Collision Transitions
        !
        read(1,*) col_ntrans(ipartner)
        if(col_ntrans(ipartner).gt.ncolltransmax) then
           ncolltransmax = col_ntrans(ipartner)
        endif
        !
        ! Read header
        !
        read(1,*) header
        !
        ! Now read the number of temperatures
        !
        read(1,*) col_ntemp(ipartner)
        if(col_ntemp(ipartner).gt.ntempmax) then
           ntempmax = col_ntemp(ipartner)
        endif
        !
        ! Read headers
        !
        read(1,*) header
        read(1,*) header   ! Here the temperatures are listed
        read(1,*) header
        !
        ! Scan over the data
        !
        if(ipartner.lt.col_npartners) then
           do itrans=1,col_ntrans(ipartner)
              read(1,*) header   ! Scan over the data           
           enddo
        endif
     enddo
     !
     ! Close and reopen the file
     !
     close(1)
     open(unit=1,file=filename)
     !
     ! Quickly scan to beginning of the collision rate data
     !
     read(1,*) header
     read(1,*) header
     read(1,*) header
     read(1,*) header
     read(1,*) header
     read(1,*) header
     read(1,*) header
     do ilevel=1,nlevels
        read(1,*) header
     enddo
     read(1,*) header
     read(1,*) header
     read(1,*) header
     do iline=1,nlines
        read(1,*) header
     enddo
     read(1,*) header
     read(1,*) header
     !
     ! Now allocate the collision rate data arrays
     !
     allocate(col_iup(1:ncolltransmax,1:col_npartners),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in line module: could not allocate temporary array col_iup'
        stop
     endif
     allocate(col_idown(1:ncolltransmax,1:col_npartners),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in line module: could not allocate temporary array col_idown'
        stop
     endif
     allocate(col_temperatures(1:ntempmax,1:col_npartners),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in line module: could not allocate temporary array col_temperatures'
        stop
     endif
     allocate(col_rates(1:ntempmax,1:ncolltransmax,1:col_npartners),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in line module: could not allocate temporary array collision_rates'
        stop
     endif
     !
     ! Now loop over the collision partners and read the data
     !
     do ipartner=1,col_npartners
        !
        ! Read header and scan over already known information
        !
        read(1,*) header
        read(1,*) header     ! Name of partner
        read(1,*) header
        read(1,*) header     ! Nr of trans
        read(1,*) header
        read(1,*) header     ! Nr of temp
        read(1,*) header
        !
        ! Now read in temperatures
        !
        read(1,*) (col_temperatures(itemp,ipartner), &
                  itemp=1,col_ntemp(ipartner))
        !
        ! Read header
        !
        read(1,*) header
        !
        ! Now loop over the transitions
        !
        do itrans=1,col_ntrans(ipartner)
           !
           ! Read
           !
           read(1,*) dum1,iup,idown,(col_rates(itemp,itrans,ipartner), &
                     itemp=1,col_ntemp(ipartner))
           !
           ! Store
           !
           col_iup(itrans,ipartner)   = iup
           col_idown(itrans,ipartner) = idown
        enddo
     enddo
  endif
200 continue
  !
  ! Close this file
  !
  close(1)
  !
  ! Store the information about the number of
  ! levels and lines
  !
  lines_nrlevels(ispec)      = nlevels
  lines_nrlines(ispec)       = nlines
  !
  ! Store the information about the number of 
  ! collisional data
  !
  lines_collisions_npartners(ispec) = col_npartners
  lines_collisions_npartmax = max(lines_collisions_npartmax,col_npartners)
  if(col_npartners.gt.0) then
     lines_collisions_ntempmax  = max(lines_collisions_ntempmax,ntempmax)
     lines_collisions_ntransmax = max(lines_collisions_ntransmax,ncolltransmax)
  endif
  !
  ! OK, if we truly want to read the data, then do it now, else return
  ! with only the lines_nrlevels(ispec), lines_nrlines(ispec),
  ! lines_collisions_npartners(ispec) and the updated
  ! lines_collisions_ntempmax and lines_collisions_ntransmax
  ! information.
  !
  if(readdata) then
     !
     ! Yes, we really want to read and store the entire data
     !
     ! First store the level info
     !
     do ilevel=1,lines_nrlevels(ispec)
        lines_level_energy_cminv(ilevel,ispec) = level_energy_cminv(ilevel)
        lines_level_gdeg(ilevel,ispec)         = level_gdeg(ilevel)
        lines_level_name(ilevel,ispec)         = level_name(ilevel)
     enddo
     !     ! Then store the line info
     !
     do iline=1,lines_nrlines(ispec)
        lines_levelup(iline,ispec)   = line_iup(iline)
        lines_leveldown(iline,ispec) = line_idown(iline)
        lines_aud(iline,ispec)       = line_aud(iline)
        lines_nu0(iline,ispec)       = line_nu0(iline)
     enddo
     !
     ! Check the frequency of the transitions for consistency
     !     
     error_linefreq = 0.d0
     do iline=1,lines_nrlines(ispec)
        nucheck = cc * (                                                     &
             lines_level_energy_cminv(lines_levelup(iline,ispec),ispec) -    &
             lines_level_energy_cminv(lines_leveldown(iline,ispec),ispec) )
        if(abs(nucheck/line_nu0(iline)-1.d0).gt.error_linefreq) then
           error_linefreq = abs(nucheck/line_nu0(iline)-1.d0)
        endif
     enddo
!     if(error_linefreq.gt.1d-7) then
!        if(error_linefreq.gt.1d-3) then
!           call integer_to_string(ispec,str1)
!           write(stdo,*) 'Warning: Level energies of molecule '//trim(str1)//' are very'
!           write(stdo,*) '   inaccurate. Relative error is: ',error_linefreq
!        else
!           call integer_to_string(ispec,str1)
!           write(stdo,*) 'Warning: Level energies of molecule '//trim(str1)//' are slightly'
!           write(stdo,*) '   inaccurate. Relative error is: ',error_linefreq
!        endif
!     endif
     if(error_linefreq.gt.1d-3) then
        call integer_to_string(ispec,str1)
        write(stdo,*) 'Warning: Level energies (or line frequency) of molecule '//trim(str1)//' are very'
        write(stdo,*) '   inaccurate. Relative error is: ',error_linefreq
     endif
     !
     ! Now create the subset of lines for which both the upper and lower
     ! levels are contained in the subset of levels
     !
     ilinesub = 1
     do iline=1,lines_nrlines(ispec)
        flag1     = .false.
        flag2     = .false.
        ilevel    = lines_levelup(iline,ispec)
        do ilevelsub=1,lines_nrlevels_subset(ispec)
           if(lines_levels_subset(ilevelsub,ispec).eq.ilevel) then
              ilevel1=ilevelsub
              flag1=.true.
           endif
        enddo
        ilevel    = lines_leveldown(iline,ispec)
        do ilevelsub=1,lines_nrlevels_subset(ispec)
           if(lines_levels_subset(ilevelsub,ispec).eq.ilevel) then
              ilevel2=ilevelsub
              flag2=.true.
           endif
        enddo
        if(flag1.and.flag2) then
           lines_lines_subset(ilinesub,ispec)  = iline
           ilinesub = ilinesub+1
        endif
     enddo
     lines_nrlines_subset(ispec) = ilinesub-1
     !
     ! Check
     !
     if(lines_nrlines_subset(ispec).lt.1) then
        call integer_to_string(ispec,str1)
        write(stdo,*) 'ERROR in line module, species '//trim(str1)
        write(stdo,*) '      No lines found for the selected subset of levels'
        stop
     endif
     !
     ! Check if the number of temperatures and transitions
     ! does not exceed maximum (in principle this should
     ! never happen, because this maximum is determined
     ! beforehand, based on the number required)
     !
     if(col_npartners.gt.0) then
        if((ntempmax.gt.lines_collisions_ntempmax).or. &
             (ncolltransmax.gt.lines_collisions_ntransmax)) then
           write(stdo,*) 'ERROR: Nr of collision rate data exceeds array sizes...'
           write(stdo,*) '       Strange, this should not happen. Try again.'
           write(stdo,*) ntempmax,lines_collisions_ntempmax, &
                ncolltransmax,lines_collisions_ntransmax
           stop
        endif
     endif
     !
     ! Store the collision data
     !
     if(col_npartners.gt.0) then
        if(.not.allocated(lines_collisions_rates)) then
           write(stdo,*) 'ERROR: Cannot store collision data because '
           write(stdo,*) '       arrays are not allocated.'
           stop
        endif
        lines_collisions_npartners(ispec) = col_npartners
        do ipartner=1,col_npartners
           lines_collisions_ntemp(ipartner,ispec)  = col_ntemp(ipartner)
           lines_collisions_ntrans(ipartner,ispec) = col_ntrans(ipartner)
           do itemp=1,col_ntemp(ipartner)
              lines_collisions_temperatures(itemp,ipartner,ispec) = &
                   col_temperatures(itemp,ipartner)
           enddo
           do itrans=1,col_ntrans(ipartner)
              lines_collisions_levelup(itrans,ipartner,ispec) = &
                   col_iup(itrans,ipartner)
              lines_collisions_leveldown(itrans,ipartner,ispec) = &
                   col_idown(itrans,ipartner)
           enddo
           do itrans=1,col_ntrans(ipartner)
              do itemp=1,col_ntemp(ipartner)
                 lines_collisions_rates(itemp,itrans,ipartner,ispec) = &
                      col_rates(itemp,itrans,ipartner)                   
              enddo
           enddo
        enddo
     endif
     !
     ! If you like, here's a level diagram in ascii format
     !
     if(lines_show_pictograms.eq.1) then 
        write(*,*) 'Level Diagram for ispec=',ispec,':'
        call lines_pictlevels(ispec)
     endif
  endif
  !
  ! Deallocate temporary arrays
  !
  deallocate(level_energy_cminv)
  deallocate(level_gdeg)
  deallocate(level_name)
  deallocate(line_iup)
  deallocate(line_idown)
  deallocate(line_aud)
  deallocate(line_nu0)
  if(allocated(col_ntemp)) deallocate(col_ntemp,col_ntrans)
  if(allocated(col_iup)) deallocate(col_iup,col_idown,col_temperatures,col_rates)
  !
  ! Done.
  !
end subroutine read_molecule_leiden



!-------------------------------------------------------------------
!                    READ LINELIST FILE
!
! Thanks to Attila Juhasz for various bugfixes in October 2011
!-------------------------------------------------------------------
subroutine read_linelist(ispec,readdata,ierror)
  implicit none
  integer :: ispec,ierror,iformat,incl_part,incl_extradata
  integer :: i,idum,iline
  logical :: readdata,fex
  character*160 filename,header,species
  double precision :: dum1,dum2,lambda,aul,elo,eup,glo,gup
  !
  ! Defaults
  ! 
  ierror = 0
  !
  ! Check if the file is present
  !
  species = lines_speciesname(ispec) 
  filename = 'linelist_'//species(1:len_trim(species))//'.inp'
  inquire(file=filename,exist=fex)
  if(.not.fex) then
     ierror = 1
     return
  endif
  !
  ! Open file and skip over some headers
  !
  open(unit=1,file=filename)
  read(1,*) header
  read(1,*) header
  read(1,*) iformat
  if(iformat.ne.1) then
     write(stdo,*) 'ERROR while reading ',trim(filename)
     write(stdo,*) '      Format number unfamiliar.'
     stop
  endif
  read(1,*) header
  read(1,*) header     ! The name of the molecule
  read(1,*) header
  read(1,*) header
  !
  ! Read the molecular weight
  !
  read(1,*) lines_umass(ispec)
  lines_umass(ispec) = lines_umass(ispec)*mp
  !
  ! Read header
  !
  read(1,*) header
  !
  ! Check if a partition function table is there,
  ! so that we can skip it.
  !
  read(1,*) incl_part
  !
  ! Read header
  !
  read(1,*) header
  !
  ! Is extra machine readable data available? 
  ! (not used for now)
  !
  read(1,*) incl_extradata
  !
  ! If partition table is present, we must skip over it
  ! because this is read elsewhere in the code.
  !
  if(incl_part.gt.0) then
     !
     ! Partition function is listed.
     !
     lines_linelist_contains_parfunc(ispec) = .true.
     !
     ! Skip over it.
     !
     read(1,*) header
     read(1,*) idum
     read(1,*) header
     do i=1,idum
        read(1,*) header
     enddo
  else
     lines_linelist_contains_parfunc(ispec) = .false.
  endif
  !
  ! Read header
  !
  read(1,*) header
  !
  ! Nr of lines
  !
  read(1,*) lines_nrlines(ispec)
  !
  ! Now read the data if requested
  ! 
  if(readdata) then
     !
     ! Check
     !
     if(lines_nrlines(ispec).gt.lines_maxnrlines) then
        write(stdo,*) 'INTERNAL ERROR while reading line list for species ',ispec
        write(stdo,*) '      Nr of lines is too large. Should not happen. Warn author.'
        stop
     endif
     !
     ! Loop over lines
     !
     read(1,*) header
     do iline=1,lines_nrlines(ispec)
        read(1,*) idum,lambda,aul,elo,eup,glo,gup
        lines_nu0(iline,ispec)            = 1d4*cc/lambda
        lines_aud(iline,ispec)            = aul
        lines_linelist_eup(iline,ispec)   = eup
        lines_linelist_edown(iline,ispec) = elo
        lines_linelist_gup(iline,ispec)   = gup
        lines_linelist_gdown(iline,ispec) = glo
     enddo
  endif
  !
  ! Close file
  !
  close(1)
end subroutine read_linelist




!-------------------------------------------------------------------
!        EXTRACT PARTITION FUNCTION DATA FROM LINELIST FILE
!
! Thanks to Attila Juhasz for various bugfixes in October 2011
!-------------------------------------------------------------------
subroutine read_parfunc_from_linelist(ispec,readdata,ierror)
  implicit none
  integer :: ispec,ierror,iformat,incl_part,itemp
  logical :: readdata,fex
  character*160 filename,header,species
  double precision :: dum1,dum2
  !
  ! Defaults
  ! 
  ierror = 0
  !
  ! Check if the file is present
  !
  species  = lines_speciesname(ispec)
  filename = 'linelist_'//species(1:len_trim(species))//'.inp'
  inquire(file=filename,exist=fex)
  if(.not.fex) then
     ierror = 1
     return
  endif
  !
  ! Open file and get the nr of partition function temperature points
  ! (or 0 if no partition function table is present)
  !
  open(unit=1,file=filename)
  read(1,*) header
  read(1,*) header
  read(1,*) iformat
  if(iformat.ne.1) then
     write(stdo,*) 'ERROR while reading ',trim(filename)
     write(stdo,*) '      Format number unfamiliar.'
     stop
  endif
  read(1,*) header
  read(1,*) header     ! The name of the molecule
  read(1,*) header
  read(1,*) header
  read(1,*) header     ! Molecular weight
  read(1,*) header
  read(1,*) incl_part
  if(incl_part.eq.0) then
     !
     ! Partition function is not listed in this file
     !
     lines_partition_ntemp(ispec) = 0
  else
     !
     ! Partition function is listed.
     !
     read(1,*) header
     read(1,*) header     ! Machine readable information
     read(1,*) header
     read(1,*) lines_partition_ntemp(ispec)
     if(lines_partition_ntemp(ispec).le.0) then
        write(stdo,*) 'ERROR while reading partition function from linelist file ',trim(filename)
        write(stdo,*) '      Somehow the nr of temperature points is strange: ',lines_partition_ntemp(ispec)
        stop
     endif
     !
     ! Now read the table, but only if requested
     !
     if(readdata) then
        read(1,*) header
        if(lines_partition_ntemp(ispec).gt.lines_partition_ntempmax) then
           write(stdo,*) 'INTENRAL ERROR while reading partition function from linelist file ',trim(filename)
           write(stdo,*) '      Somehow the nr of temperature points is larger than allowed: ',lines_partition_ntemp(ispec)
           write(stdo,*) '      This hsould not have happened...'
           stop
        endif
        do itemp=1,lines_partition_ntemp(ispec)
           read(1,*) dum1,dum2
           lines_partition_temperature(itemp,ispec) = dum1
           lines_partition_function(itemp,ispec)    = dum2
        enddo
     endif
     !
  endif
  !
  ! Close file
  !
  close(1)
end subroutine read_parfunc_from_linelist




!-------------------------------------------------------------------
!           MAKE GRAPHICAL PICTOGRAM OF LEVELS AND TRANSITIONS
!-------------------------------------------------------------------
subroutine lines_pictlevels(ispec)
  implicit none
  integer :: ispec
  integer :: i,k,nlinmax,istart,iend
  character*83 pict
  !
  istart=14
  nlinmax = min(lines_nrlines_subset(ispec),69)
  do k=lines_nrlevels_subset(ispec),1,-1
     write(pict(1:10),403) lines_level_energy_cminv(lines_levels_subset(k,ispec),ispec)
403  format(1F10.5)
     do i=11,istart
        pict(i:i) = ' '
     enddo
     iend = 1
     do i=1,nlinmax
        if((k.ge.lines_leveldown(lines_lines_subset(i,ispec),ispec)).and. &
           (k.le.lines_levelup(lines_lines_subset(i,ispec),ispec))) then
14         pict(i+istart:i+istart) = '|'
           iend = i+istart
        else
           pict(i+istart:i+istart) = ' '
        endif
     enddo
     do i=nlinmax+istart+1,83
        pict(i:i) = ' ' 
     enddo
     write(*,*) pict (1:iend+1)
  enddo
  !
end subroutine lines_pictlevels




!-------------------------------------------------------------------
!                   READ THE LEVEL POPULATIONS
!
! NOTE: These are the full level populations, i.e. in units of
!       1/cm^3. The molecular number density is already included
!       in these numbers. 
!-------------------------------------------------------------------
subroutine read_levelpop(action)
  implicit none
  character*80 :: filename1,filename2,filename3,filename,species
  integer :: ispec,ierr,ilevel,index,i,irec,n,nlevels,precis
  integer :: ntmp,itemp,action,idum,style
  integer(kind=8) :: iiformat,reclen,reclend,nn,kk,nnlevels
  logical :: fex1,fex2,fex3
  integer(kind=8), allocatable :: subset(:)
  !
  ! Default
  !
  precis = 8
  !
  ! Consistency check
  !
  if(lines_mode.lt.0) stop 9351
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(lines_levelpop)) return
  endif
  !
  ! Message
  !
  write(stdo,*) 'Reading level populations...'
  call flush(stdo)
  !
  ! Checks
  !
  if(lines_nrlevels_subset_max.le.0) then
     write(stdo,*) 'ERROR: level subset not yet specified while reading level populations.'
     stop
  endif
  !
  ! Remove any already existing array
  !
  if(allocated(lines_levelpop)) deallocate(lines_levelpop)
  !
  ! Create the level populations arrays
  !
  allocate(lines_levelpop(1:lines_nrlevels_subset_max,1:lines_nr_species,1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate the lines_levelpop() array'
     stop
  endif
  !
  ! Now loop over all species
  !
  do ispec=1,lines_nr_species
     if(allocated(subset)) deallocate(subset)  ! Bugfix 2017.08.30
     if(lines_species_fullmolec(ispec)) then
        !
        ! Create the file name and find if the input file is present and if
        ! it is in text format (.dat) or unformatted (.udat or .bdat)
        !
        species   = lines_speciesname(ispec)
        filename1 = 'levelpop_'//species(1:len_trim(species))//'.dat'
        filename2 = 'levelpop_'//species(1:len_trim(species))//'.udat'
        filename3 = 'levelpop_'//species(1:len_trim(species))//'.bdat'
        inquire(file=filename1,exist=fex1)
        inquire(file=filename2,exist=fex2)
        inquire(file=filename3,exist=fex3)
        idum=0
        if(fex1) idum=idum+1
        if(fex2) idum=idum+1
        if(fex3) idum=idum+1
        if(idum.gt.1) then
           write(stdo,*) 'ERROR: Found more than one file levelpop_'//species(1:len_trim(species))//'.*dat'
           stop
        endif
        if(idum.eq.0) then
           write(stdo,*) 'ERROR: Could not find any levelpop_'//species(1:len_trim(species))//'.*dat file'
           stop
        endif
        !
        ! Read data
        !
        if(fex1) then
           !
           ! Open formatted
           !
           style = 1
           open(unit=1,file=filename1)
           !
           ! Read format number
           !
           read(1,*) iiformat
           !
           ! Read number of grid points
           !
           read(1,*) nn
           !
           ! Read number of levels 
           !
           read(1,*) nnlevels
           !
           ! Check if this equals the number of levels of the subset
           !
           if(nnlevels.ne.lines_nrlevels_subset(ispec)) then
              write(stdo,*) 'ERROR in line module: nr of levels in the file '
              write(stdo,*) filename1(1:len_trim(filename1)),' unequal to the'
              write(stdo,*) 'number of selected levels of this molecule.'
              stop
           endif
           nlevels = nnlevels
           !
           ! Allocate temporary arrays
           !
           allocate(subset(1:nnlevels),STAT=ierr)
           if(ierr.ne.0) then
              write(stdo,*) 'ERROR in line module: could not allocate subset()'
              stop
           endif
           !
           ! Read the list of levels
           !
           read(1,*) (subset(ilevel),ilevel=1,nlevels)
        elseif(fex2) then
           !
           ! Open f77-unformatted (with records)
           !
           style = 2
           open(unit=1,file=filename2,status='old',form='unformatted')
           read(1) iiformat,reclen
           if(iiformat.ne.1) then
              write(stdo,*) 'ERROR: Format number of '//TRIM(filename)//' is invalid/unknown.'
              write(stdo,*) 'Format number = ',iiformat
              write(stdo,*) 'Record length = ',reclen
              stop
           endif
           !
           ! Nr of cells
           ! 
           read(1) nn
           !
           ! Nr of levels
           !
           read(1) nnlevels
           if(nnlevels.ne.lines_nrlevels_subset(ispec)) then
              write(stdo,*) 'ERROR in line module: nr of levels in the file '
              write(stdo,*) trim(filename),' unequal to the'
              write(stdo,*) 'number of selected levels of this molecule.'
              stop
           endif
           nlevels = nnlevels
           !
           ! Read the level subset
           !
           allocate(subset(1:nlevels),STAT=ierr)
           if(ierr.ne.0) then
              write(stdo,*) 'ERROR in line module: could not allocate subset()'
              stop
           endif
           read(1) (subset(ilevel),ilevel=1,nlevels)
        else
           !
           ! C-compliant binary format
           !
           style = 3
           open(unit=1,file=filename3,status='old',access='stream')
           read(1) iiformat
           if(iiformat.ne.1) then
              write(stdo,*) 'ERROR: Format number of '//TRIM(filename)//' is invalid/unknown.'
              write(stdo,*) 'Format number = ',iiformat
              write(stdo,*) 'Record length = ',reclen
              stop
           endif
           read(1) nn
           precis = nn
           !
           ! Nr of cells
           ! 
           read(1) nn
           !
           ! Nr of levels
           !
           read(1) nnlevels
           if(nnlevels.ne.lines_nrlevels_subset(ispec)) then
              write(stdo,*) 'ERROR in line module: nr of levels in the file '
              write(stdo,*) trim(filename),' unequal to the'
              write(stdo,*) 'number of selected levels of this molecule.'
              stop
           endif
           nlevels = nnlevels
           !
           ! Read the level subset
           !
           allocate(subset(1:nlevels),STAT=ierr)
           if(ierr.ne.0) then
              write(stdo,*) 'ERROR in line module: could not allocate subset()'
              stop
           endif
           read(1) (subset(ilevel),ilevel=1,nlevels)
        endif
        !
        ! Check if they are identical to those that were selected in
        ! the lines.inp file
        !
        do ilevel=1,nlevels
           if(subset(ilevel).ne.lines_levels_subset(ilevel,ispec)) then
              write(stdo,*) 'ERROR in line module: the subset of levels in the file'
              write(stdo,*) filename1(1:len_trim(filename1)),' unequal to those'
              write(stdo,*) 'selected in the lines.inp file.'
              write(stdo,*) 'Here: ',subset(:)
              write(stdo,*) 'Lines.inp: ',lines_levels_subset(:,ispec)
              stop
           endif
        enddo
        !
        ! Check if nr of grid points is equal to that of the grid.inp
        !
        if(nn.ne.nrcellsinp) then
           write(stdo,*) 'ERROR in line module: nr of grid points specified in '
           write(stdo,*) '    ',filename1(1:len_trim(filename1)),' unequal to the'
           write(stdo,*) '     nr specified in the grid file.'
           stop
        endif
        !
        ! Now read the levelpopulations
        !
        call read_vectorfield(1,style,precis,lines_nrlevels_subset_max, &
             nlevels,nrcellsinp,lines_nr_species,ispec,floor=1d-99,     &
             vector1=lines_levelpop)
        !
        ! Close the file
        !
        close(1)
     endif
  enddo ! Loop over species
  !
  ! Deallocate
  !
  if(allocated(subset)) deallocate(subset)
  !
end subroutine read_levelpop


!-------------------------------------------------------------------
!         READ THE PARTITION FUNCTION FOR ALL MOLECULES
!
! Read, or generate internally, the partition functions for each of
! the species.
!
! Note: If lines_mode.ge.1 and there exists a level population file
!       for some species (or for all species), then of course the
!       partition function is not required. Otherwise it is required.
!       We will check for this.
!-------------------------------------------------------------------
subroutine read_partition_function(action)
  implicit none
  integer :: action
  integer :: ispec,ntmp,iformat,itmp,ilevel,itemp,ierr,ierror
  double precision :: dum1,dum2,lpop,gast,dummy,dendivk,const
  character*80 :: filename1,filename2,species
  logical :: fex1,fex2
  parameter(const=hh*cc/kk)
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(lines_partition_temperature)) return
  endif
  !
  ! Message
  !
  write(stdo,*) 'Reading or computing partition functions...'
  call flush(stdo)
  !
  ! Do stupidity checks
  !
  if(lines_nr_species.le.0) stop 8010
  if(.not.allocated(lines_speciesname)) stop 8011
  if(.not.allocated(gastemp)) then
     write(stdo,*) 'ERROR in line module: wanting to compute/read the partition function, but'
     write(stdo,*) '     the gas temperature has not been read into memory yet.'
     stop
  endif
  if(.not.allocated(gas_chemspec_numberdens)) then
     write(stdo,*) 'ERROR in line module: while compute/read the partition function:'
     write(stdo,*) '     The number density of the molecule(s) has not been read into memory yet.'
     write(stdo,*) '     While they are not required for the partition function, '
     write(stdo,*) '     they will be necessary later.'
     stop 
  endif
  !-----------------------------------------------------------------
  ! Attila Juhasz
  ! lines_level_energy_cminv is not allocated in linelist mode
  !-----------------------------------------------------------------
  !if(.not.allocated(lines_level_energy_cminv)) stop 8014
  if(lines_maxnrlevels.gt.0) then
     if(.not.allocated(lines_level_energy_cminv)) stop 8014
  endif
  !-----------------------------------------------------------------
  !
  ! First we need to find out the maximum number of temperature sampling
  ! points for the partition function. 
  !
  lines_partition_ntempmax = 0
  do ispec=1,lines_nr_species
     !
     ! First figure out where to seek the partition function data
     !
     if(lines_linelist_contains_parfunc(ispec)) then
        !
        ! This molecule is apparently a linelist molecule, and 
        ! the partition function data is included in this list.
        !
        call read_parfunc_from_linelist(ispec,.false.,ierror)
        if(ierror.ne.0) then
           write(stdo,*) 'INTERNAL ERROR while reading partition function from line list'
           stop
        endif
        ntmp = lines_partition_ntemp(ispec)
        if(ntmp.gt.lines_partition_ntempmax) then
           lines_partition_ntempmax = ntmp
        endif
        lines_partition_ntemp(ispec) = 0
     else
        !
        ! Partition function should be in a separate file, if present
        !
        species   = lines_speciesname(ispec)
        filename2  = 'partitionfunction_'//species(1:len_trim(species))//'.inp'
        inquire(file=filename2,exist=fex2)
        if(fex2) then
           !
           ! A file for the partition function exists for this species. So
           ! get the number of temperature sampling points from this file.
           !
           open(unit=2,file=filename2)
           read(2,*) iformat
           read(2,*) ntmp
           close(2)
           if(ntmp.gt.lines_partition_ntempmax) then
              lines_partition_ntempmax = ntmp
           endif
        else
           !
           ! No file exists for the partition function for this species.
           ! This automatically means that we are going to compute the
           ! partition function internally. For this we need a number of
           ! sampling points given by lines_partition_ntempint.
           !
           if(lines_partition_ntempint.gt.lines_partition_ntempmax) then
              lines_partition_ntempmax = lines_partition_ntempint
           endif
        endif
     endif
  enddo
  !
  ! Deallocate stuff if necessary
  !
  if(allocated(lines_partition_temperature)) deallocate(lines_partition_temperature)
  if(allocated(lines_partition_function)) deallocate(lines_partition_function)
  !
  ! Now allocate the arrays for the partition function
  !
  if(.not.allocated(lines_partition_ntemp)) then
     write(stdo,*) 'INTERNAL ERROR: lines_partition_ntemp() not allocated.'
     stop
  endif
  allocate(lines_partition_temperature(1:lines_partition_ntempmax,1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate the lines_partition_temperature() array'
     stop
  endif
  allocate(lines_partition_function(1:lines_partition_ntempmax,1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate the lines_partition_function() array'
     stop
  endif
  !
  ! Now get the partition function for all species, either by reading
  ! them from a file or by computing them on-the-fly
  !
  do ispec=1,lines_nr_species
     !
     ! Get the name of the species
     !
     species  = lines_speciesname(ispec)
     !
     ! Check if there is a partition function 
     !
     if(lines_linelist_contains_parfunc(ispec)) then
        !
        ! This molecule is apparently a linelist molecule, and 
        ! the partition function data is included in this list.
        !
        write(stdo,*) 'Reading partition function from linelist for ispec = ',ispec
        call read_parfunc_from_linelist(ispec,.true.,ierror)
        if(ierror.ne.0) then
           write(stdo,*) 'INTERNAL ERROR while reading partition function from line list'
           stop
        endif
     else
        !
        ! Partition function should be in a separate file, if present
        !
        filename1 = 'partitionfunction_'//species(1:len_trim(species))//'.inp'
        inquire(file=filename1,exist=fex1)
        if(fex1) then
           !
           ! Now open the file for the partition function
           !
           write(stdo,*) 'Reading ',filename1
           open(unit=1,file=filename1)
           read(1,*) iformat
           read(1,*) ntmp
           lines_partition_ntemp(ispec) = ntmp
           do itemp=1,ntmp
              read(1,*) dum1,dum2
              lines_partition_temperature(itemp,ispec) = dum1
              lines_partition_function(itemp,ispec)    = dum2
           enddo
           close(1)
           !
           ! When using an externally computed partition function we must
           ! warn the user that RADMC-3D uses ilevel=1 as the ground
           ! state energy
           !
           write(*,*) '   Note: When using an external partition function, make sure that level 1 is the ground state.'
           !
        else
           !
           ! No file available, so we create the partition function internally.
           !
           ! First check: If this molecule is a linelist molecule, then
           ! we cannot compute the partition function internally
           !
           if(.not.lines_species_fullmolec(ispec)) then
              write(stdo,*) 'ERROR: Cannot find the partition function data for molecule ',ispec
              write(stdo,*) '       Since this molecule is a linelist molecule, I cannot compute'
              write(stdo,*) '       the partition function internally either.'
              stop
           endif
           !
           !
           lines_partition_ntemp(ispec) = lines_partition_ntempint
           !
           ! Note: hh * cc / kk = 1.4387519
           ! NOTE: Now we use the constants_module.f90 for this. 
           !
           ! Loop over temperatures
           !
           write(stdo,*) 'Computing partition function internally for molecule ',ispec
           do itemp=1,lines_partition_ntemp(ispec)
              !
              ! Create a temperature point
              !
              lines_partition_temperature(itemp,ispec) =                       &
                   lines_partition_temp0 *                                     &
                   (lines_partition_temp1/lines_partition_temp0)**             &
                   (dble(itemp-1.d0)/dble(lines_partition_ntemp(ispec)-1.d0))
              !
              ! Now compute the partition function based on all the available
              ! levels for this molecule. 
              !
              ! Note: In particular for high temperatures this internal
              ! calculation may fail.
              !
              lpop  = lines_level_gdeg(1,ispec)
              dummy = lpop
              gast  = lines_partition_temperature(itemp,ispec)
              do ilevel=2,lines_nrlevels(ispec)
                 dendivk = const *                                           &
                           ( lines_level_energy_cminv(ilevel,ispec) -        &
                             lines_level_energy_cminv(ilevel-1,ispec) )
                 lpop    = lpop * exp(-dendivk/gast) *                       &
                           lines_level_gdeg(ilevel,ispec) /                  &
                           lines_level_gdeg(ilevel-1,ispec)
                 dummy   = dummy + lpop
              enddo
              lines_partition_function(itemp,ispec) = dummy
           enddo
        endif
     endif
  enddo ! Loop over species
end subroutine read_partition_function


!-------------------------------------------------------------------
!              MAKE LTE LEVEL POPULATIONS AT ONE POINT
!-------------------------------------------------------------------
subroutine lines_compute_ltepop(ispec,nlev,gast,nrdens,levelpop)
  implicit none
  integer :: nlev,ispec,itemp,ilevel
  double precision :: levelpop(1:nlev),gast,eps,pfunc,dendivk,nrdens,dummy
  !
  ! Check
  !
  if(nlev.ne.lines_nrlevels(ispec)) then
     write(stdo,*) 'INTERNAL ERROR in lines_compute_ltepop()'
     stop
  endif
  !
  ! Find the location of the current temperature in the partition function table
  !
  call hunt(lines_partition_temperature(:,ispec),               &
            lines_partition_ntemp(ispec),gast,itemp)
  if(itemp.ge.lines_partition_ntemp(ispec)) then
     write(stdo,*) 'ERROR: Temperature out of range of partition function '
     write(stdo,*) '       for ispec=',ispec,' Temperature = ',gast
     !write(stdo,*) '       Smallest temperature in table = ',   &
     !     lines_partition_temperature(1,ispec)
     write(stdo,*) '       Largest temperature in table = ',   &
          lines_partition_temperature(lines_partition_ntemp(ispec),ispec)
     stop
  endif
  if(itemp.ge.1) then
     eps = (gast-lines_partition_temperature(itemp,ispec)) /            &
          (lines_partition_temperature(itemp+1,ispec)-                  &
          lines_partition_temperature(itemp,ispec))
  else
     itemp = 1
     eps   = 0.d0
     lines_pfunc_outofrange = .true.
  endif
  if((eps.lt.0.).or.(eps.gt.1.)) stop 7729
  !
  ! Construct the partition function by linear interpolation in the table
  !
  pfunc = (1-eps)*lines_partition_function(itemp,ispec) +            &
              eps*lines_partition_function(itemp+1,ispec)
  !
  ! Now construct all the populations according to LTE
  ! Note that c*h/k=1.4387519 in CGS units.
  !
  dummy = nrdens / pfunc
  do ilevel=1,nlev
     dendivk = 1.4387519 *                               &
             ( lines_level_energy_cminv(ilevel,ispec)    &
             - lines_level_energy_cminv(1,ispec) )
     levelpop(ilevel) = dummy * exp(-dendivk/gast) *     &
             lines_level_gdeg(ilevel,ispec)
  enddo
  !
end subroutine lines_compute_ltepop


!-------------------------------------------------------------------
!              MAKE LTE LEVEL POPULATIONS AT ONE POINT
!                  VERSION: COMPUTE ONLY SUBSET
!
! Note: Here you can specify a local subset of levels, which does
!       not necessarily have to be the same as the subset specified
!       in lines_levels_subset(:,:). This is useful for on-the-fly
!       computation of level populations, where you would not want
!       to compute level populations that are not involved in this
!       particular spectrum or image. 
!
! BUGFIX: I used the ilevsubsetidx for addressing the levelpop. 
!         But I decided to keep levelpop as the full set of
!         levels. Only for the storage in the big array (popul at
!         every location in the grid) do we use the subset level
!         index. I now simply removed ilevsubsetidx here.
! Bugfix: The levelpop array is bigger than nlev. Now added nlevtot.
!-------------------------------------------------------------------
subroutine lines_compute_ltepop_subset(nlevtot,nlev,ilevsubset,    &
                                ispec,gast,nrdens,levelpop)
  implicit none
  integer :: nlevtot,nlev,ispec,itemp,ilevel
  integer :: ilevsubset(1:nlev)
  double precision :: levelpop(1:nlevtot),gast,eps,pfunc,dendivk,nrdens,dummy
  !
  ! If no levels are selected, then return 
  !
  if(nlev.eq.0) return
  !
  ! Find the location of the current temperature in the partition function table
  !
  call hunt(lines_partition_temperature(:,ispec),               &
            lines_partition_ntemp(ispec),gast,itemp)
  if(itemp.ge.lines_partition_ntemp(ispec)) then
     write(stdo,*) 'ERROR: Temperature out of range of partition function '
     write(stdo,*) '       for ispec=',ispec,' Temperature = ',gast
     !write(stdo,*) '       Smallest temperature in table = ',   &
     !     lines_partition_temperature(1,ispec)
     write(stdo,*) '       Largest temperature in table = ',   &
          lines_partition_temperature(lines_partition_ntemp(ispec),ispec)
     stop
  endif
  if(itemp.ge.1) then
     eps = (gast-lines_partition_temperature(itemp,ispec)) /            &
          (lines_partition_temperature(itemp+1,ispec)-                  &
          lines_partition_temperature(itemp,ispec))
  else
     itemp = 1
     eps   = 0.d0
     lines_pfunc_outofrange = .true.
  endif
  if((eps.lt.0.).or.(eps.gt.1.)) stop 7729
  !
  ! Construct the partition function by linear interpolation in the table
  !
  pfunc = (1-eps)*lines_partition_function(itemp,ispec) +            &
              eps*lines_partition_function(itemp+1,ispec)
  !
  ! Now construct all the populations according to LTE
  ! Note that c*h/k=1.4387519 in CGS units.
  !
  dummy = nrdens / pfunc
  do ilevel=1,nlev
     dendivk = 1.4387519 *                                           &
             ( lines_level_energy_cminv(ilevsubset(ilevel),ispec)    &
             - lines_level_energy_cminv(1,ispec) )
     levelpop(ilevsubset(ilevel)) = dummy * exp(-dendivk/gast) *  &
             lines_level_gdeg(ilevsubset(ilevel),ispec)
  enddo
  !
end subroutine lines_compute_ltepop_subset


!-------------------------------------------------------------------
!           RECURSIVE ROUTINE FOR VELOCITY AVERAGING
!-------------------------------------------------------------------
recursive subroutine lines_amr_recursive_velo(cell,idir,velo)
  implicit none
  integer :: idir
  double precision :: velo,velochild
  type(amr_branch), pointer :: cell,child
  integer :: ix,iy,iz,indexchild,cnt
  !
  ! Reset velo
  !
  velo = 0.d0
  !
  ! Loop over all children
  !
  cnt = 0
  do iz=1,1+amr_zdim
     do iy=1,1+amr_ydim
        do ix=1,1+amr_xdim
           if(.not.associated(cell%child(ix,iy,iz)%link)) then
              write(stdo,*) 'ERROR in LVG: Children lost...'
              stop
           endif
           child => cell%child(ix,iy,iz)%link
           if(child%leaf) then
              indexchild = child%leafindex
              if(indexchild.eq.0) then
                 write(stdo,*) 'ERROR in lines_amr_recursive_velo: indexchild=0'
                 stop
              endif
              velo = velo + gasvelocity(idir,indexchild) 
           else
              call lines_amr_recursive_velo(child,idir,velochild)
              velo = velo + velochild
           endif
           cnt  = cnt + 1
        enddo
     enddo
  enddo
  ! 
  ! Compute the average velo
  !
  velo = velo / cnt
  !
end subroutine lines_amr_recursive_velo

!-------------------------------------------------------------------
!         COMPUTE VELOCITY GRADIENT FOR LVG (SOBOLEV) 
!
! Written by Rahul Shetty and Cornelis Dullemond.
!
! NOTE: This can be stored in an array, since it does not change.
!       Will be done in a later version.
! 
! If maxdiff==.true. then compute *maximum* velocity *difference*
!-------------------------------------------------------------------
subroutine lines_compute_velgradient(index,velgradient,maxdiff)
  implicit none
  double precision :: velgradient,velgradsum,factor(1:3)
  double precision :: pos0(1:3),pos1,vel0(1:3),vel1,adj_ds
  integer :: index,numvels,level
  type(amr_branch), pointer :: cell,neighbor
  integer :: ixx,iyy,izz,indexn,ixn,iyn,izn
  logical, optional :: maxdiff
  logical :: average
  !
  ! Set the average flag: If average is set, then compute the
  ! average velocity gradient (default). If not, compute the
  ! *maximum* velocity *difference*.
  !
  average = .true.
  if(present(maxdiff)) then
     if(maxdiff) then
        average = .false.
     endif
  endif
  !
  ! Check...
  !
  if(.not.allocated(gasvelocity)) then
     write(*,*) 'ERROR: No gas velocity array...'
     stop
  endif
  !
  ! Reset variables
  !
  velgradsum = 0.d0
  numvels    = 0
  !
  ! Get cell position
  !
  if(amr_tree_present) then
     !
     ! The AMR tree is available, so we use it
     !
     !!! BUG (30.12.2011): cell  => amr_theleafs(index)%link
     cell  => amr_index_to_leaf(index)%link
     pos0(1) = amr_finegrid_xc(cell%ixyzf(1),1,cell%level)
     pos0(2) = amr_finegrid_xc(cell%ixyzf(2),2,cell%level)
     pos0(3) = amr_finegrid_xc(cell%ixyzf(3),3,cell%level)
  else
     !
     ! We have a regular grid, so we compute the position
     ! using ixx,iyy,izz
     !
     call amr_regular_get_ixyz(index,ixx,iyy,izz)
     pos0(1) = amr_finegrid_xc(ixx,1,0)
     pos0(2) = amr_finegrid_xc(iyy,2,0)
     pos0(3) = amr_finegrid_xc(izz,3,0)
  endif
  !
  ! Get factors for coordinate systems
  !
  if(igrid_coord.lt.100) then
     factor(1:3) = 1.d0
  elseif(igrid_coord.lt.200) then
     factor(1) = 1.d0
     factor(2) = pos0(1)
     factor(3) = pos0(1)*sin(pos0(2))
  else
     stop 34
  endif
  !
  ! Get velocities in the current cell
  !
  vel0(1) = gasvelocity(1,index)
  vel0(2) = gasvelocity(2,index)
  vel0(3) = gasvelocity(3,index)
  !
  ! Now find the velocities in the neighboring cells and
  ! for each one found, compute a contribution to the 
  ! overall average velocity gradient.
  !
  if(amr_tree_present) then
     !
     ! The AMR tree is available, so we use that, and we
     ! will take care of the refinement levels. Note,
     ! however, that if a much bigger cell is next to
     ! the present one, the gradient is necessarily 
     ! going to be a bit inaccurate (skewed):
     !
     !   +---+-------+
     !   | * |       |
     !   +---+   *   |
     !       |       |
     !       +-------+
     !
     ! BUGFIX 2016.07.26: 
     !   The way the stuff was previously implemented the ds was
     !   multiplied by 2^leveldifference. However, the difference
     !   in cell center position due to the different cell size 
     !   is already accounted for by the amr_finegrid_xc... So
     !   I think I accidently did this twice. 
     !
     ! Get the AMR level
     !
     level = cell%level
     !
     ! Check X-direction only if active
     !
     if(amr_xdim.eq.1) then
        !
        ! Find velocity of adjacent X location (left)
        !
        neighbor => cell%neighbor(1,1)%link
        if (associated(neighbor)) then 
           pos1 = amr_finegrid_xc(neighbor%ixyzf(1),1,neighbor%level)
           adj_ds     = pos0(1) - pos1
           if(neighbor%leaf) then
              vel1 = gasvelocity(1,neighbor%leafindex) 
              ! BUGFIX 2016.07.26: 
              !if(neighbor%level.lt.level) then
              !adj_ds = adj_ds * (2**(level-neighbor%level))
              !endif
           else
              call lines_amr_recursive_velo(neighbor,1,vel1)
           endif
           if(average) then
              velgradsum = velgradsum + abs((vel0(1) - vel1)/adj_ds)
           else
              velgradsum = max(velgradsum,abs((vel0(1) - vel1)))
           endif
           numvels    = numvels + 1
        endif
        !
        ! Find velocity of adjacent X location (right)
        !
        neighbor => cell%neighbor(2,1)%link
        if (associated(neighbor)) then 
           pos1 = amr_finegrid_xc(neighbor%ixyzf(1),1,neighbor%level)
           adj_ds     = pos1 - pos0(1)
           if(neighbor%leaf) then
              vel1 = gasvelocity(1,neighbor%leafindex) 
              ! BUGFIX 2016.07.26: 
              !if(neighbor%level.lt.level) then
              !adj_ds = adj_ds * (2**(level-neighbor%level))
              !endif
           else
              call lines_amr_recursive_velo(neighbor,1,vel1)
           endif
           if(average) then
              velgradsum = velgradsum + abs((vel1 - vel0(1))/adj_ds)
           else
              velgradsum = max(velgradsum,abs((vel1 - vel0(1))))
           endif
           numvels    = numvels + 1
        endif
     endif
     !
     ! Check Y-direction only if active
     !
     if(amr_ydim.eq.1) then
        !
        ! Find velocity of adjacent Y location (left)
        !
        neighbor => cell%neighbor(1,2)%link
        if (associated(neighbor)) then 
           pos1 = amr_finegrid_xc(neighbor%ixyzf(2),2,neighbor%level)
           adj_ds     = ( pos0(2) - pos1 ) * factor(2)
           if(neighbor%leaf) then
              vel1 = gasvelocity(2,neighbor%leafindex) 
              ! BUGFIX 2016.07.26: 
              !if(neighbor%level.lt.level) then
              !adj_ds = adj_ds * (2**(level-neighbor%level))
              !endif
           else
              call lines_amr_recursive_velo(neighbor,2,vel1)
           endif
           if(average) then
              velgradsum = velgradsum + abs((vel0(2) - vel1)/adj_ds)
           else
              velgradsum = max(velgradsum,abs((vel0(2) - vel1)))
           endif
           numvels    = numvels + 1
        endif
        !
        ! Find velocity of adjacent Y location (right)
        !
        neighbor => cell%neighbor(2,2)%link
        if (associated(neighbor)) then 
           pos1 = amr_finegrid_xc(neighbor%ixyzf(2),2,neighbor%level)
           adj_ds     = ( pos1 - pos0(2) ) * factor(2)
           if(neighbor%leaf) then
              vel1 = gasvelocity(2,neighbor%leafindex) 
              ! BUGFIX 2016.07.26: 
              !if(neighbor%level.lt.level) then
              !adj_ds = adj_ds * (2**(level-neighbor%level))
              !endif
           else
              call lines_amr_recursive_velo(neighbor,2,vel1)
           endif
           if(average) then
              velgradsum = velgradsum + abs((vel1 - vel0(2))/adj_ds)
           else
              velgradsum = max(velgradsum,abs((vel1 - vel0(2))))
           endif
           numvels    = numvels + 1
        endif
     endif
     !
     ! Check Z-direction only if active
     !
     if(amr_zdim.eq.1) then
        !
        ! Find velocity of adjacent Z location (left)
        !
        neighbor => cell%neighbor(1,3)%link
        if (associated(neighbor)) then 
           pos1 = amr_finegrid_xc(neighbor%ixyzf(3),3,neighbor%level)
           adj_ds     = ( pos0(3) - pos1 ) * factor(3)
           if(neighbor%leaf) then
              vel1 = gasvelocity(3,neighbor%leafindex) 
              ! BUGFIX 2016.07.26: 
              !if(neighbor%level.lt.level) then
              !adj_ds = adj_ds * (2**(level-neighbor%level))
              !endif
           else
              call lines_amr_recursive_velo(neighbor,3,vel1)
           endif
           if(average) then
              velgradsum = velgradsum + abs((vel0(3) - vel1)/adj_ds)
           else
              velgradsum = max(velgradsum,abs((vel0(3) - vel1)))
           endif
           numvels    = numvels + 1
        endif
        !
        ! Find velocity of adjacent Z location (right)
        !
        neighbor => cell%neighbor(2,3)%link
        if (associated(neighbor)) then 
           pos1 = amr_finegrid_xc(neighbor%ixyzf(3),3,neighbor%level)
           adj_ds     = ( pos1 - pos0(3) ) * factor(3)
           if(neighbor%leaf) then
              vel1 = gasvelocity(3,neighbor%leafindex)
              ! BUGFIX 2016.07.26: 
              !if(neighbor%level.lt.level) then
              !adj_ds = adj_ds * (2**(level-neighbor%level))
              !endif
           else
              call lines_amr_recursive_velo(neighbor,3,vel1)
           endif
           if(average) then
              velgradsum = velgradsum + abs((vel1 - vel0(3))/adj_ds)
           else
              velgradsum = max(velgradsum,abs((vel1 - vel0(3))))
           endif
           numvels    = numvels + 1
        endif
     endif
  else
     !
     ! We have a regular grid, so things are much easier. But we
     ! cannot use the AMR tree. 
     !
     if(amr_xdim.eq.1) then
        !
        ! Find velocity of adjacent X location (left)
        !
        ixn = ixx-1
        if(ixn.lt.1) then
           if(amr_cyclic_xyz(1)) ixn=amr_grid_nx
        else
           pos1       = amr_finegrid_xc(ixn,1,0)
           adj_ds     = pos0(1) - pos1
           indexn     = index-1
           vel1       = gasvelocity(1,indexn)
           if(average) then
              velgradsum = velgradsum + abs((vel0(1) - vel1)/adj_ds)
           else
              velgradsum = max(velgradsum,abs((vel0(1) - vel1)))
           endif
           numvels    = numvels + 1
        endif
        !
        ! Find velocity of adjacent X location (right)
        !
        ixn = ixx+1
        if(ixn.gt.amr_grid_nx) then
           if(amr_cyclic_xyz(1)) ixn=1
        else
           pos1       = amr_finegrid_xc(ixn,1,0)
           adj_ds     = pos1 - pos0(1)
           indexn     = index+1
           vel1       = gasvelocity(1,indexn)
           if(average) then
              velgradsum = velgradsum + abs((vel1 - vel0(1))/adj_ds)
           else
              velgradsum = max(velgradsum,abs((vel1 - vel0(1))))
           endif
           numvels    = numvels + 1
        endif
     endif
     !
     ! Check Y-direction only if active
     !
     if(amr_ydim.eq.1) then
        !
        ! Find velocity of adjacent Y location (left)
        !
        iyn = iyy-1
        if(iyn.lt.1) then
           if(amr_cyclic_xyz(2)) iyn=amr_grid_ny
        else
           pos1       = amr_finegrid_xc(iyn,2,0)
           adj_ds     = ( pos0(2) - pos1 ) * factor(2)
           indexn     = index-amr_grid_nx
           vel1       = gasvelocity(2,indexn)
           if(average) then
              velgradsum = velgradsum + abs((vel0(2) - vel1)/adj_ds)
           else
              velgradsum = max(velgradsum,abs((vel0(2) - vel1)))
           endif
           numvels    = numvels + 1
        endif
        !
        ! Find velocity of adjacent Y location (right)
        !
        iyn = iyy+1
        if(iyn.gt.amr_grid_ny) then
           if(amr_cyclic_xyz(2)) iyn=1
        else
           pos1       = amr_finegrid_xc(iyn,2,0)
           adj_ds     = ( pos1 - pos0(2) ) * factor(2)
           indexn     = index + amr_grid_nx
           vel1       = gasvelocity(2,indexn)
           if(average) then
              velgradsum = velgradsum + abs((vel1 - vel0(2))/adj_ds)
           else
              velgradsum = max(velgradsum,abs((vel1 - vel0(2))))
           endif
           numvels    = numvels + 1
        endif
     endif
     !
     ! Check Z-direction only if active
     !
     if(amr_zdim.eq.1) then
        !
        ! Find velocity of adjacent Z location (left)
        !
        izn = izz-1
        if(izn.lt.1) then
           if(amr_cyclic_xyz(3)) izn=amr_grid_nz
        else
           pos1       = amr_finegrid_xc(izn,3,0)
           adj_ds     = ( pos0(3) - pos1 ) * factor(3)
           indexn     = index-(amr_grid_nx*amr_grid_ny)
           vel1       = gasvelocity(3,indexn)
           if(average) then
              velgradsum = velgradsum + abs((vel0(3) - vel1)/adj_ds)
           else
              velgradsum = max(velgradsum,abs((vel0(3) - vel1)))
           endif
           numvels    = numvels + 1
        endif
        !
        ! Find velocity of adjacent Z location (right)
        !
        izn = izz+1
        if(izn.gt.amr_grid_nz) then
           if(amr_cyclic_xyz(3)) izn=1
        else
           pos1       = amr_finegrid_xc(izn,3,0)
           adj_ds     = ( pos1 - pos0(3) ) * factor(3)
           indexn     = index+(amr_grid_nx*amr_grid_ny)
           vel1       = gasvelocity(3,indexn)
           if(average) then
              velgradsum = velgradsum + abs((vel1 - vel0(3))/adj_ds)
           else
              velgradsum = max(velgradsum,abs((vel1 - vel0(3))))
           endif
           numvels    = numvels + 1
        endif
     endif
  endif
  !
  ! Check
  !
  if(numvels.eq.0) then
     write(stdo,*) 'Warning for LVG method: Apparently there is only 1 cell, because'
     write(stdo,*) '      I cannot find neighbors of this cell to determine the '
     write(stdo,*) '      velocity gradient. Setting velocity gradient to 0.'
     write(stdo,*) '      This means only the escape probability method part is'
     write(stdo,*) '      active (for which a length scale has to be set in a'
     write(stdo,*) '      file called escprob_lengthscale.inp), otherwise you'
     write(stdo,*) '      will automatically get LTE).'
     velgradient = 0.d0
     return
  endif
  !
  ! Now divide velgradsum by the number of contributions to obtain
  ! the average
  !
  if(average) then
     velgradient = velgradsum / numvels
  else
     velgradient = velgradsum     ! Note, here it is not the gradient but the difference
  endif
  !
end subroutine lines_compute_velgradient


!--------------------------------------------------------------------------------
!                 SOLVE EQUATIONS OF STATISTICAL EQUILIBRIUM
!
! This subroutine solves the equations of statistical equilibrium.
!
!       ___n
!        \  /                                                    \
!        /  | n_i A_ik beta_ik - ( n_k*B_ki - n_i*B_ik ) Jbar_ik | 
!       --- \                                                    /
!       i>k
!       ___n
!        \  /                                                    \
!     -  /  | n_k A_ki beta_ik - ( n_i*B_ik - n_k*B_ki ) Jbar_ik | 
!       --- \                                                    /
!       i<k
!       ___n
!        \  /                     \
!     +  /  | n_i C_ik - n_k C_ki | = 0
!       --- \                     /
!        i
!
! where beta_ik should be set to 1 for the "normal" statistical equilibrium
! equations. For the escape probability method you must set beta_ik to the
! escape probability for that line and you must set Jbar_ik=0 for that line.
! For Accelerated Lambda Iteration with a local operator Lambda_ik^* you must
! set beta_ik = ( 1 - Lambda_ik^* ) and replace Jbar_ik with
! ( Jbar_ik - Jbar_ik^* ) = ( Lambda_ik - Lambda_ik^* ) [S_ik]. 
!
! These equations are not entirely independent. Only the first n-1 are
! linearly independent. The last one is therefore rejected and replaced by
! the constraint equation:
!
!       ___n
!        \       
!        /   n_i = 1
!       ---      
!        i
!
! INPUT:
!   ispec                          The molecular species
!   nlev                           Number of levels of this molecule 
!                                  (must be identical to lines_nrlevels(ispec))
!   nlin                           Number of lines of this molecule
!                                  (must be identical to lines_nrlines(ispec))
!   Jbar(1:nlin)                   Line-averaged mean intensity at each line
!   beta(1:nlin)                   The beta_ik explained above at each line
!   col_npartner                   Number of collision partners
!   col_ntrmax                     Max nr of transitions
!   col_ntrans(1:npart)            Number of collision rate specifications
!   col_iup(1:ntrmax,1:npart)      Upper level of each rate
!   col_ilo(1:ntrmax,1:npart)      Lower level of each rate
!   col_rateud(1:ntrmax,1:npart)   Rate for collisional transition from up->lo
!   col_ratedu(1:ntrmax,1:npart)   Rate for collisional transition from lo->up
!
! RESULT:
!   popul(1:nlev)                  The resulting level populations
!   error                          =0 --> All is fine
!                                  =1 --> Negative populations
!                                  =2 --> Singular matrix
!                                  =3 --> Badly normalized
!--------------------------------------------------------------------------------
subroutine lines_solve_statequil(ispec,nlev,nlin,jbar,   &
           beta,col_npartner,col_ntrmax,col_ntrans,      &
           col_iup,col_ilo,col_rateud,col_ratedu,popul,  &
           error)
  implicit none
  integer :: ispec,nlev,nlin,col_ntrmax,iline,ilevel,klevel
  integer :: col_npartner
  integer :: col_ntrans(1:col_npartner)
  double precision :: jbar(nlin),beta(nlin),popul(nlev)
  double precision :: bud(nlin), bdu(nlin)
  double precision :: matrix(nlev,nlev)
  integer :: indx(nlev)
  double precision :: col_rateud(col_ntrmax,col_npartner)
  double precision :: col_ratedu(col_ntrmax,col_npartner)
  integer :: col_iup(col_ntrmax,col_npartner)
  integer :: col_ilo(col_ntrmax,col_npartner)
  double precision :: const1,const2
  double precision :: gratio,dummy
  integer :: ipartner,itrans,error
  logical :: success
  parameter(const1=hh/(4*pi))
  parameter(const2=cc*cc/(2*hh))
  !
  ! Default
  !
  error = 0
  !
  ! Compute Einstein B cooefficients 
  !
  do iline=1,nlin
     gratio     = lines_level_gdeg(lines_levelup(iline,ispec),ispec) /   &
                  lines_level_gdeg(lines_leveldown(iline,ispec),ispec)
     bud(iline) = lines_aud(iline,ispec) * const2 / lines_nu0(iline,ispec)**3
     bdu(iline) = bud(iline) * gratio
  enddo
  !
  ! First clear the matrix
  !
  do ilevel=1,nlev
     do klevel=1,nlev
        matrix(klevel,ilevel) = 0.d0
     enddo
  enddo
  !
  ! Loop over lines
  !
  do iline = 1,nlin
     !     
     ! This line has:
     !
     !    up   = klevel
     !    down = ilevel
     !
     klevel = lines_levelup(iline,ispec)
     ilevel = lines_leveldown(iline,ispec)
     !
     ! So the level k can be populated from level i by 
     !
     !    n_i B_ik Jbar 
     !
     matrix(klevel,ilevel) = matrix(klevel,ilevel) + bdu(iline) * Jbar(iline)
     !
     ! and depopulated to level i by
     !
     !   - n_k ( A_ki + B_ki Jbar )
     !
     matrix(klevel,klevel) = matrix(klevel,klevel) - bud(iline) * Jbar(iline) &
                             - lines_aud(iline,ispec) * beta(iline)
     !
     ! This line has:
     !
     !         up   = ilevel
     !         down = klevel
     !
     klevel = lines_leveldown(iline,ispec)
     ilevel = lines_levelup(iline,ispec)
     !
     ! So the level k can be populated from level i by 
     !
     !         n_i ( A_ik + B_ik Jbar ) 
     !
     matrix(klevel,ilevel) = matrix(klevel,ilevel) + bud(iline) * Jbar(iline) &
                             + lines_aud(iline,ispec) * beta(iline)
     !
     ! and depopulated to level i by
     !
     !   - n_k ( B_ki Jbar )
     !
     matrix(klevel,klevel) = matrix(klevel,klevel) - bdu(iline) * Jbar(iline)
     !
  enddo ! End loop over lines
  !
  ! Now all collisional transitions
  !
  do ipartner=1,col_npartner
     do itrans=1,col_ntrans(ipartner)
        !     
        ! This transition has:
        !
        !     up   = klevel
        !     down = ilevel
        !
        klevel = col_iup(itrans,ipartner)
        ilevel = col_ilo(itrans,ipartner)
        !
        ! So the level k can be populated from level i by 
        !
        !     n_i C_ik 
        !
        matrix(klevel,ilevel) = matrix(klevel,ilevel) +      &
                                col_ratedu(itrans,ipartner)
        !
        ! and depopulated to level i by
        !
        !   - n_k C_ki 
        !
        matrix(klevel,klevel) = matrix(klevel,klevel) -      &
                                col_rateud(itrans,ipartner) 
        !
        ! This transition has:
        !
        !     up   = ilevel
        !     down = klevel
        !
        klevel = col_ilo(itrans,ipartner)
        ilevel = col_iup(itrans,ipartner)
        !
        ! So the level k can be populated from level i by 
        !
        !     n_i Cik
        !
        matrix(klevel,ilevel) = matrix(klevel,ilevel) +      & 
                                col_rateud(itrans,ipartner)
        !
        ! and depopulated to level i by
        !
        !   - n_k C_ki
        !
        matrix(klevel,klevel) = matrix(klevel,klevel) -      &
                                col_ratedu(itrans,ipartner)
        !
     enddo ! End loop over icoltrans
  enddo ! End loop over collision partners
  !
  ! Now we replace the first equation with the consistency equation (see above)
  ! to make all the level populations sum to 1.
  !
  ! The consist eq. is at ilevel=1
  !
  do ilevel=1,nlev
     matrix(1,ilevel) = 1.d0
  enddo
  !
  ! Now the right hand side of the equation. We use the array popul for this.
  !
  do ilevel=1,nlev
     popul(ilevel) = 0.d0
  enddo
  !
  ! Since the consist eq. is at ilevel=1:
  !
  popul(1) = 1.d0
  !
  ! Now we can solve the matrix equation. 
  !
  dummy = 1.d0
  call ludcmp(matrix,nlev,nlev,indx,dummy,success)
  if(success) then
     call lubksb(matrix,nlev,nlev,indx,popul)
  else
     write(stdo,*) 'SINGULAR MATRIX IN SOLVING STATISTICAL EQUILIBRIUM'
     error = 2
     return
  endif
  !
  ! Verify if everything is indeed positive, and check normalization
  !
  dummy = 0.d0
  do ilevel=1,nlev
     dummy = dummy + popul(ilevel) 
     if(popul(ilevel).lt.0.d0) then
        write(stdo,*) 'WARNING: Negative popul for ilev '
        call flush(stdo)
        popul(ilevel) = 0.d0
        error = 1
     endif
  enddo
  !
  if(abs(1.d0-dummy).gt.1.d-2) then
     write(stdo,*) 'Problem: normalization levels is not accurate'
     error = 3
  endif
  !
end subroutine lines_solve_statequil


!--------------------------------------------------------------------------------
!                COMPUTE THE COLLISION RATES FROM THE BASE RATE ARRAY
!
! We have the basic collision rates stored in the lines_collisions_rates() 
! array. But to perform an actual calculation of the statistical equilibrium
! we have to extract the rates for a particular temperature, and compute the
! inverse rates as well. This is what is done here, for a single molecular
! species, but for all the collision partners of that species. This subroutine
! computes the lines_col_** array components.
!--------------------------------------------------------------------------------
subroutine lines_compute_collrates(ispec,nlev,nlin,npartner,npartnermax,   &
                                   gast,collpartner_numden)
  implicit none
  integer :: nlev,nlin,ispec,iline,itemp,itrans,npartner,npartnermax,ipartner
  double precision :: gast,collpartner_numden(1:npartnermax),eps
  double precision :: gratio,dendivk
  !
  ! Loop over collision partners
  !
  do ipartner=1,lines_collisions_npartners_used(ispec)
     !
     ! Need to interpolate from temperature table in molecule information
     ! file to get collision rates
     !
     call hunt(lines_collisions_temperatures(:,ipartner,ispec),   &
               lines_collisions_ntemp(ipartner,ispec),gast,itemp)
     if(itemp.ge.lines_collisions_ntemp(ipartner,ispec)) then
        write(stdo,*) 'ERROR: Temperature out of range in collisional temperature table '
        write(stdo,*) '       for ispec=',ispec,' ipartner=',ipartner,' temperature = ',gast
        write(stdo,*) '       Largest temperature in table = ',   &
             lines_collisions_temperatures(itemp,ipartner,ispec)
        call flush(stdo)
        stop
     endif
     if(itemp.ge.1) then
        eps = (gast-lines_collisions_temperatures(itemp,ipartner,ispec)) / &
             (lines_collisions_temperatures(itemp+1,ipartner,ispec)-       &
              lines_collisions_temperatures(itemp,ipartner,ispec))
     else
        itemp = 1
        eps   = 0.d0
     endif
     if((eps.lt.0.).or.(eps.gt.1.)) stop 7729
     !
     ! Compute the collisional rates by interpolation
     !
     do itrans=1,lines_collisions_ntrans(ipartner,ispec)
        !
        ! Set the iup and idown
        !
        lines_col_iup(itrans,ipartner) =                                  &
             lines_collisions_levelup(itrans,ipartner,ispec)
        lines_col_idown(itrans,ipartner) =                                &
             lines_collisions_leveldown(itrans,ipartner,ispec)
        !
        ! Set collisional de-excitation rates - use simple linear interpolation
        !
        lines_col_rateud(itrans,ipartner) = collpartner_numden(ipartner) *     &
             ((1.d0-eps)*lines_collisions_rates(itemp,itrans,ipartner,ispec) + &
             eps*lines_collisions_rates(itemp+1,itrans,ipartner,ispec))
        !
        ! Ratio of level degeneracies
        !
        gratio = lines_level_gdeg(lines_col_iup(itrans,ipartner),ispec) /   &
                 lines_level_gdeg(lines_col_idown(itrans,ipartner),ispec)
        !
        ! Energy difference between levels
        ! Note that c*h/k=1.4387519 in CGS units.
        !
        dendivk = 1.4387519 *                                                     &
             ( lines_level_energy_cminv(lines_col_iup(itrans,ipartner),ispec)     &
             - lines_level_energy_cminv(lines_col_idown(itrans,ipartner),ispec) )
        !
        ! Compute the collisional excitation rates, using Einstein relation
        !
        lines_col_ratedu(itrans,ipartner) = gratio *                      &
             lines_col_rateud(itrans,ipartner) * exp(-dendivk/gast)
        !
     enddo ! End loop over transitions
  enddo ! End loop over collision partners
  !
end subroutine lines_compute_collrates


!--------------------------------------------------------------------------------
!       COMPUTE LEVEL POPULATIONS AT ONE POINT USING LVG APPROXIMATION
!
! For a given molecular species (ispec), for given number density of the
! collision partners and the molecule itself, for given tempearture and
! for given velocity gradient (|dv/dx| averaged over all directions), 
! this routine will return the level populations of the molecule.
! 
! The method is very similar to the iterative method described in the
! paper by van der Tak et al. 2007, A&& 468, 627. At some point this will
! be accelerated using the method of Rybicki & Hummer (1991), A&A 245,
! 171. 
!
! ARGUMENTS:
!   ispec                          Molecule
!   nlev                           Number of levels of this molecule 
!                                  (must be identical to lines_nrlevels(ispec))
!   nlin                           Number of lines of this molecule
!                                  (must be identical to lines_nrlines(ispec))
!   npartner                       Number of collision partners for this molec
!   npartnermax                    Max nr of col partners (for array sizes)
!   gast                           Gas temperature [K]
!   nrdens                         Number density of the molecule [cm^{-3}]
!   collpartner_numdens            Number densities of all collision partners
!                                  (including those that are not used by this
!                                  particular molecule, but perhaps by another
!                                  molecule: This subroutine will find out which
!                                  to use. This way you only have to give the
!                                  full array, not make a pre-selection).
!   velgrad                        The velocity gradient, averaged over 4*pi.
!                                  The unit is [s^{-1}], because it is |dv/dx|.
!                                  Note: This must be a positive number.
! 
! OPTIONAL THINGS:
!   lengthscale                    If set, then allow also escape because
!                                  of a finite size of the object (EscProb)
!   turb                           If lengthscale is set, then the micro-
!                                  turbulent line width is required too.
!   dbg                            If .true. then write the convergence
!                                  sequence. Only for debugging. 
!
! ADDITIONAL PARAMETERS:
!   lines_nonlte_maxiter           This integer sets the max number of iterations
!                                  for the non-LTE iteration.
!
! RETURNS:
!   levelpop                       The level populations of this molecule. 
!                                  The sum of these equals nrdens.
!                                  (Note: Inside this subroutine levelpop
!                                         is the relative one, summing to 1;
!                                         Only upon return this is multiplied
!                                         by nrdens)
! 
! Programmed by Rahul Shetty and Cornelis Dullemond.
!--------------------------------------------------------------------------------
recursive subroutine lines_compute_lvgpop(ispec,nlev,nlin,npartner,         &
                                npartnermax,gast,nrdens,collpartner_numden, &
                                velgrad,levelpop,lengthscale,turb,dbg)
  implicit none
  integer :: nlev,nlin,ispec,itemp,itrans, iter, converge
  integer :: iline,iup,idown,ilevel,npartner,npartnermax
  double precision :: levelpop(1:nlev),gast,eps,dendivk,nrdens,dummy,gratio
  double precision :: velgrad
  double precision :: Jbar(1:nlin), beta(1:nlin), tau
  double precision :: levelpopold(nlev),levelpopold2(nlev)
  double precision :: collpartner_numden(1:npartnermax)
  double precision, optional :: lengthscale,turb
  integer :: indx(nlev),error
  double precision :: const1,const2,fgaus, newTex, convfactor, thisconvergefactor
  double precision :: tau0,bud,bdu,phi0,a
  logical :: doslow,maser
  logical, optional :: dbg
  parameter(const1=hh/(4*pi))
  parameter(const2=cc*cc/(2*hh))
  parameter(fgaus =cc*cc*cc/(8.0*pi*1.0645))
  !
  ! First, obtain the collisional rates for this molecule, at this
  ! temperature, for all the collision partners of this molecule
  !
  call lines_compute_collrates(ispec,nlev,nlin,npartner,npartnermax,   &
                               gast,collpartner_numden)
  !
  ! Increase a counter
  !
  lines_lvgesc_nrsolves = lines_lvgesc_nrsolves + 1
  !
  ! Reset populations and other things
  !
  levelpop(:)    = 0.d0
  levelpopold(:) = 0.d0
  error = 0
  !
  ! Start the iteration for the LVG method
  !
  do iter = 1, lines_nonlte_maxiter
     !
     ! Increase a counter 
     !
     lines_lvgesc_nriter = lines_lvgesc_nriter + 1
     !
     ! Set the maser warning flag to .false.
     !
     maser = .false.
     !
     ! Backup the old populations, so that we can check for convergence
     !
     levelpopold2(:) = levelpopold(:)
     levelpopold(:)  = levelpop(:)
     !
     ! Calculate the Jbar of all the lines
     !
     do iline=1,nlin
        !
        ! Get the lower and upper level for this lines
        !
        idown = lines_leveldown(iline,ispec)
        iup   = lines_levelup(iline,ispec)
        !
        ! Make distinction between first and later iterations
        !
        if(iter.eq.1) then 
           !
           ! Starting guess
           !
           beta(iline) = 0.
        else
           !
           ! For vanderTak et al. 2007 definition:
           ! Using velgrad
           !
           if(velgrad.gt.0.d0) then
              tau = fgaus / lines_nu0(iline,ispec)**3 *                     &
                    lines_aud(iline,ispec) * nrdens / ( velgrad + 1d-90 ) * &
                    (levelpop(idown)*lines_level_gdeg(iup,ispec)/           &
                    lines_level_gdeg(idown,ispec) - levelpop(iup))
           else
              tau = 1d90
           endif
           !
           ! If requested, then limit the tau according to the
           ! length scale given. This is a simple implementation
           ! of escape probability
           !
           if(present(lengthscale)) then
              !
              ! We have to compute the total optical depth at line
              ! center, given nrdens, the populations and the 
              ! length scale. Let us do this only for the Gaussian
              ! line profile for now.
              !
              if(lines_profile.ne.0) then
                 write(stdo,*) 'ERROR: Escape probability mode in LVG for the moment'
                 write(stdo,*) '       only possible for gaussian line profile.'
                 stop
              endif
              if(.not.present(turb)) then
                 write(stdo,*) 'ERROR: For escape probability in LVG mode we need a'
                 write(stdo,*) '       a value for the turbulence...'
                 stop
              endif
              !
              ! Now calculate the center of the Gaussian line profile
              !
              a      = sqrt(turb**2 + 2*kk*gast/lines_umass(ispec))
              phi0   = cc / (a*lines_nu0(iline,ispec)*sqrtpi)
              !
              ! Now calculate the Bdu and Bud from Aud
              !
              gratio = lines_level_gdeg(iup,ispec) /   &
                       lines_level_gdeg(idown,ispec)
              bud    = lines_aud(iline,ispec) * const2 / lines_nu0(iline,ispec)**3
              bdu    = bud * gratio
              !
              ! Now calculate the line-center optical depth in case
              ! no velocity gradient is present
              !
              tau0   = (hh*lines_nu0(iline,ispec)/(4*pi)) * phi0 *   &
                       ( levelpop(idown)*bdu - levelpop(iup)*bud ) * &
                       nrdens * lengthscale
              !
              ! Now, if tau exceeds tau0, then we keep tau = tau0, so that
              ! the optical depth can never exceed tau0
              !
              if(tau.gt.tau0) then
                 tau = tau0
              endif
           endif
           !
           ! As of version 0.36 we prevent masering here
           !
           if(tau.lt.0.d0) then
              tau   = 0.d0
              maser = .true.
           endif
           !
           ! Now calculate the beta coefficient
           !
           if(tau.gt.1.0d-5) then
              beta(iline) = (1.d0 - exp(-tau))/tau  ! For LVG
              if ( beta(iline) .gt. 1.d0 ) beta(iline) = 1.d0
           else 
              beta(iline) = 1.d0 - 0.5d0*tau + tau**2 / 6.d0
           endif
        endif
        !
        ! The only part of Jbar that we explicitly have to feed
        ! into the statistical equilibrium equations is the 
        ! external part:
        !
        Jbar(iline) = beta(iline) * bplanck(lines_tbg,lines_nu0(iline,ispec))
        !
     enddo ! End loop over lines
     !
     ! Solve the statistical equilibrium equations
     !
     call lines_solve_statequil(ispec,nlev,nlin,jbar,beta, &
           lines_collisions_npartners_used(ispec),         &
           lines_collisions_ntransmax,                     &
           lines_collisions_ntrans(:,ispec),               &
           lines_col_iup,lines_col_idown,                  &
           lines_col_rateud,lines_col_ratedu,              &
           levelpop,error)
     !
     ! In case of a serious (!) error, stop this iteration
     !
     if(error.eq.2) then
        goto 500
     endif
     !
     ! For debugging
     !
     if(present(dbg)) then
        if(dbg) then
           write(stdo,*) 'Iteration ',iter
           write(stdo,*) levelpop
        endif
     endif
     !
     ! Reset very small (or negative!) levelpop:
     !
     do iline = 1,nlin
        !
        iup   = lines_levelup(iline,ispec)
        idown = lines_leveldown(iline,ispec)
        !
        if (levelpop(iup) < 1.0e-30) levelpop(iup)=1.0e-30
        if (levelpop(idown) < 1.0e-30) levelpop(idown)=1.0e-30
     enddo
     !
     ! Renormalize populations (not necessary, therefore removed)
     !
     !     dummy = 0.d0
     !     do ilevel=1,nlev
     !        dummy = dummy + levelpop(ilevel)
     !     enddo
     !     dummy = 1.d0/dummy
     !     do ilevel=1,nlev
     !        levelpop(ilevel) = levelpop(ilevel) * dummy
     !     enddo
     !
     ! Check if we have flipflopping in any of the levels
     !
     doslow = .false.
     do ilevel=1,nlev
        if((levelpop(ilevel)-levelpopold(ilevel))*                  &
           (levelpopold(ilevel)-levelpopold2(ilevel)).lt.0.d0) then
           doslow = .true.
        endif
     enddo
     !
     ! If we have flipflopping in any of the levels, we do
     ! a slowdown of the convergence to prevent overshoot
     !
     if(doslow) then
        levelpop = lines_solving_speed*levelpop +            &
                   (1.d0-lines_solving_speed)*levelpopold
     endif
     !
     ! Check if solution converged
     !
     converge = 0
     if ( iter > 1 ) then 
        do ilevel=1,nlev
           convfactor = (levelpop(ilevel) - levelpopold(ilevel))/(levelpopold(ilevel)+1d-90)
           thisconvergefactor = abs(convfactor)
           if (thisconvergefactor.gt.lines_nonlte_convcrit) then
              exit
           else
              converge = converge + 1
           endif
        enddo
     endif
     !
     ! Calculate total population levels if convergence is achieved 
     !             ... or indicate why convergence was not achieved
     !
     if (converge >= nlev) then 
        !
        ! Convert the level population to the 1/cm^3 level population
        !
        levelpop = nrdens*levelpop
        !
        ! TODO: Generalize the tau calculation (for user define line, not simply 1-0)
        !         poptau(pos)=fgaus/lines_nu0(1,ispec)**3 * lines_aud(1,ispec) * nrdens/velgrad &
        !              *(levelpop(1)*lines_level_gdeg(2,ispec)/ &
        !              lines_level_gdeg(1,ispec) -levelpop(2))
        !
        ! If masers were detected in this last iteration, then
        ! set the global maser warning flag
        !
        if(maser) then
           lines_maser_warning = .true.
        endif
        !
        ! Exit the loop and the subroutine
        !
        return
     endif
  enddo ! End loop over iter
  !
  ! Apparently no convergence
  !
500 continue
  !
  ! If lines_slowlvg_as_alternative==.true. then instead of giving a warning
  ! we do the old (slow but more stable) style LVG
  !
  if(lines_slowlvg_as_alternative) then
     !
     ! Try the old style
     !
     call lines_compute_lvgpop_li(ispec,nlev,nlin,npartner,npartnermax, &
                                  gast,nrdens,collpartner_numden,       &
                                  velgrad,levelpop)
  else
     !
     ! Give up, but get some diagnostics first
     !
     if(.not.present(dbg)) then
        if(error.eq.2) then
           write(stdo,*) 'LVG FAILURE! Singular matrix after ', iter, 'lvg iterations...'
        else
           write(stdo,*) 'LVG FAILURE! No convergence after', iter, 'lvg iterations...'
           write(stdo,*) ' ... no convergence at level', ilevel, '(out of', &
                nlev,'):',thisconvergefactor, '>', lines_nonlte_convcrit
        endif
        call flush(stdo)
        !
        ! Recursively call this routine again, this time with debugging info
        !
        write(stdo,*) '*** REDOING THE ITERATION WHILE WRITING OUT POPULATIONS: ***'
        write(stdo,*) 'The subroutine lines_compute_lvgpop() has been called'
        write(stdo,*) 'with the following parameters:'
        write(stdo,*) 'ispec        = ',ispec,', name = ',trim(lines_speciesname(ispec))
        write(stdo,*) 'nlev         = ',nlev
        write(stdo,*) 'nlin         = ',nlin
        write(stdo,*) 'npartner     = ',npartner
        write(stdo,*) 'gast         = ',gast
        write(stdo,*) 'nrdens       = ',nrdens
        write(stdo,*) 'partn_nrdens = ',collpartner_numden(1:npartner)
        write(stdo,*) 'velgrad      = ',velgrad
        if(present(lengthscale)) then
           write(stdo,*) 'lengthscl    = ',lengthscale
           write(stdo,*) 'turb         = ',turb
        endif
        write(stdo,*) 'Here are the populations for each iteration:'
        if(present(lengthscale)) then
           call lines_compute_lvgpop(ispec,nlev,nlin,npartner,         &
                           npartnermax,gast,nrdens,collpartner_numden, &
                           velgrad,levelpop,lengthscale=lengthscale,   &
                           turb=turb,dbg=.true.)
        else
           call lines_compute_lvgpop(ispec,nlev,nlin,npartner,         &
                           npartnermax,gast,nrdens,collpartner_numden, &
                           velgrad,levelpop,dbg=.true.)
        endif
     endif
     write(stdo,*) '*** FINISHED REDOING THE ITERATION WHILE WRITING OUT POPULATIONS: ***'
     write(stdo,*) 'Aborting RADMC-3D...'
     write(stdo,*) 'Note: You might try setting lines_slowlvg_as_alternative=1'
     write(stdo,*) '      in radmc3d.inp. This will retry with a slower but'
     write(stdo,*) '      more stable LVG method in case of bad convergence.'
     write(stdo,*) '      It is a clunky method, but for now it is the best solution'
     write(stdo,*) '      until a better fix for the fast LVG is found.'
     stop
  endif
  !
  ! Done...
  !
end subroutine lines_compute_lvgpop


!--------------------------------------------------------------------------------
!       COMPUTE LEVEL POPULATIONS AT ONE POINT USING LVG APPROXIMATION
!                     (OLDER, BUT MORE STABLE VERSION)
! 
! The new LVG method is faster but can sometimes lead to flipflopping. Here
! is the older, much slower but more stable, version of LVG. The difference
! is that this older version uses the Lambda Iteration scheme while the
! newer uses the Accelerated Lambda Iteration scheme.
! 
! Programmed by Rahul Shetty and Cornelis Dullemond.
!--------------------------------------------------------------------------------
subroutine lines_compute_lvgpop_li(ispec,nlev,nlin,npartner,npartnermax, &
                                  gast,nrdens,collpartner_numden,       &
                                  velgrad,levelpop)
  implicit none
  integer :: nlev,nlin,ispec,itemp,itrans, iter, converge
  integer :: iline,iup,idown,ilevel,npartner,npartnermax
  double precision :: levelpop(1:nlev),gast,eps,dendivk,nrdens,gratio
  double precision :: velgrad
  double precision :: Jbar(1:nlin), S_nu, beta, tau
  double precision :: levelpopold(nlev),lines_Tex(nlin)
  double precision :: collpartner_numden(1:npartnermax)
  integer :: indx(nlev),error
  double precision :: dummyarr(1:nlin)
  double precision :: const1,const2,fgaus, newTex, convfactor, thisconvergefactor
  parameter(const1=hh/(4*pi))
  parameter(const2=cc*cc/(2*hh))
  parameter(fgaus =cc*cc*cc/(8.0*pi*1.0645))
  !
  ! First, obtain the collisional rates for this molecule, at this
  ! temperature, for all the collision partners of this molecule
  !
  call lines_compute_collrates(ispec,nlev,nlin,npartner,npartnermax,   &
                               gast,collpartner_numden)
  !
  ! Set dummy array to 1
  !
  dummyarr(:) = 1.d0
  !
  ! Reset populations
  !
  levelpop(:) = 0.d0
  !
  ! Start the iteration for the LVG method
  !
  do iter = 1, lines_nonlte_maxiter
     !
     ! Backup the old populations, so that we can check for convergence
     !
     levelpopold(:) = levelpop(:)
     !
     ! Calculate the Jbar of all the lines
     !
     do iline=1,nlin
        !
        ! Get the lower and upper level for this lines
        !
        idown = lines_leveldown(iline,ispec)
        iup   = lines_levelup(iline,ispec)
        !
        ! Make distinction between first and later iterations
        !
        if(iter.eq.1) then 
           !
           ! Starting guess
           !
           beta = 0.
           lines_Tex(iline) = gast
        else
           !
           ! For vanderTak et al. 2007 definition:
           ! Using velgrad
           !
           tau = fgaus / lines_nu0(iline,ispec)**3 *             &
                 lines_aud(iline,ispec) * nrdens / velgrad *     &
                 (levelpop(idown)*lines_level_gdeg(iup,ispec)/   &
                 lines_level_gdeg(idown,ispec) - levelpop(iup))
           !
           ! Now calculate the beta coefficient
           !
           if ( tau >= 1.0d-15) then
              beta = (1.d0 - exp(-tau))/tau  ! For LVG
              ! beta = (1.d0 - exp(-3.0*tau))/(3.0*tau)  ! For Slab (shock)
              if ( beta .gt. 1.d0 ) beta = 1.d0
           else 
              beta = 1.d0 
           endif
           !
           ! Finally calculate the next excitation temperature from last
           ! iteration's populations
           !
           newTex = -1.0*hh*lines_nu0(iline,ispec)/ kk /   &
                    dlog(levelpop(iup)/levelpop(idown)*    &
                    lines_level_gdeg(idown,ispec)/         &
                    lines_level_gdeg(iup,ispec)) 
           !
           ! Only update Tex if new Tex > 0 
           !
           if ( newTex .gt. 0.0 ) then 
              lines_Tex(iline) = 0.5*( lines_Tex(iline)+newTex )
              ! lines_Tex(iline) = newTex
           endif
           !
        endif
        !
        ! Calculate Snu, which is only depends on the line frequency
        ! Ignore dust emissivity and absorbtion coefficients
        !   Source function: S_nu = n_u A_ud / [n_d B_du - n_u B_ud]  
        !                         = 2 h nu^3 / c^2 / [exp(Enu/KT) - 1]
        !                    Sud  = nu^2/const2 / [exp(hnu/KT) - 1]
        !
        S_nu        = lines_nu0(iline,ispec)**3/const2/                          &
                      (exp(hh*lines_nu0(iline,ispec)/lines_Tex(iline)/kk) - 1.0)
        Jbar(iline) = (1.0-beta) * S_nu + beta * bplanck(lines_tbg,lines_nu0(iline,ispec))
        !
     enddo ! End loop over lines
     !
     ! Solve the statistical equilibrium equations
     !
     call lines_solve_statequil(ispec,nlev,nlin,jbar,dummyarr, &
           lines_collisions_npartners_used(ispec),             &
           lines_collisions_ntransmax,                         &
           lines_collisions_ntrans(:,ispec),                   &
           lines_col_iup,lines_col_idown,                      &
           lines_col_rateud,lines_col_ratedu,                  &
           levelpop,error)
     !
     ! Update excitation temperatures, and reset very small (or negative!) levelpop:
     !
     do iline = 1,nlin
        !
        iup   = lines_levelup(iline,ispec)
        idown = lines_leveldown(iline,ispec)
        !
        if (levelpop(iup) < 1.0e-30) levelpop(iup)=1.0e-30
        if (levelpop(idown) < 1.0e-30) levelpop(idown)=1.0e-30
        !
        ! For vanderTak et al. 2007 definition:
        !
        tau = fgaus/lines_nu0(iline,ispec)**3 * lines_aud(iline,ispec) * nrdens/velgrad &
              *(levelpop(idown)*lines_level_gdeg(iup,ispec)/ &
              lines_level_gdeg(idown,ispec) -levelpop(iup))                     
        ! 
        ! New excitation temperature
        !
        newTex = -1.0*hh*lines_nu0(iline,ispec)/ kk &
                 / dlog(levelpop(iup)/levelpop(idown)* &
                 lines_level_gdeg(idown,ispec)/lines_level_gdeg(iup,ispec)) 
        !
        ! Only update Tex if new Tex > 0 (and/or only for thick lines):
        !   if ( tau .gt. 0.01 .and. newTex .gt. 0.0 ) then 
        !
        if ( newTex .gt. 0.0 ) then 
           lines_Tex(iline) = 0.5*( lines_Tex(iline)+newTex )
           ! lines_Tex(iline) = newTex
        endif
     enddo
     !
     ! Check if solution converged
     !
     converge = 0
     if ( iter > 1 ) then 
        do ilevel=1,nlev
           convfactor = (levelpop(ilevel) - levelpopold(ilevel))/(levelpopold(ilevel)+1d-90)
           thisconvergefactor = abs(convfactor)
           if (thisconvergefactor.gt.lines_nonlte_convcrit) then
              exit
           else
              converge = converge + 1
           endif
        enddo
     endif
     !
     ! Calculate total population levels if convergence is achieved 
     !             ... or indicate why convergence was not achieved
     if (converge >= nlev) then 
        levelpop = nrdens*levelpop
        !
        ! TODO: Generalize the tau calculation (for user define line, not simply 1-0)
        !         poptau(pos)=fgaus/lines_nu0(1,ispec)**3 * lines_aud(1,ispec) * nrdens/velgrad &
        !              *(levelpop(1)*lines_level_gdeg(2,ispec)/ &
        !              lines_level_gdeg(1,ispec) -levelpop(2))
        !
        exit
     else if (iter .ge. lines_nonlte_maxiter-2) then 
        write(stdo,*) 'LVG FAILURE! No convergence after', iter, 'lvg iterations...'
        write(stdo,*) ' ... no convergence at level', ilevel, '(out of', &
             nlev,'):',thisconvergefactor, '>', lines_nonlte_convcrit
        call flush(stdo)
     endif
     !
     !
  enddo ! End loop over iter
  !
  ! Done...
  !
end subroutine lines_compute_lvgpop_li


!--------------------------------------------------------------------------------
!    COMPUTE LEVEL POPULATIONS AT ONE POINT USING OPTICALL THIN APPROXIMATION
!
! For a given molecular species (ispec), for given number density of the
! collision partners and the molecule itself and for given tempearture,
! this routine will return the level populations of the molecule. It assumes
! that the line emission is perfectly optically thin, so that there is no
! reabsorption of any line radiation that is emitted. It is valid, for 
! instance, in very optically thin nebulae. 
! 
! ARGUMENTS:
!   ispec                          Molecule
!   nlev                           Number of levels of this molecule 
!                                  (must be identical to lines_nrlevels(ispec))
!   nlin                           Number of lines of this molecule
!                                  (must be identical to lines_nrlines(ispec))
!   npartner                       Number of collision partners for this molec
!   npartnermax                    Max nr of col partners (for array sizes)
!   gast                           Gas temperature [K]
!   nrdens                         Number density of the molecule [cm^{-3}]
!   collpartner_numdens            Number densities of all collision partners
!                                  (including those that are not used by this
!                                  particular molecule, but perhaps by another
!                                  molecule: This subroutine will find out which
!                                  to use. This way you only have to give the
!                                  full array, not make a pre-selection).
!
! RETURNS:
!   levelpop                       The level populations of this molecule. 
!                                  The sum of these equals nrdens.
! 
! Programmed by Rahul Shetty and Cornelis Dullemond.
!--------------------------------------------------------------------------------
subroutine lines_compute_optthinpop(ispec,nlev,nlin,npartner,npartnermax, &
                                gast,nrdens,collpartner_numden,levelpop)
  implicit none
  integer :: nlev,nlin,ispec,itemp,itrans, iter, converge
  integer :: iline,iup,idown,ilevel,npartner,npartnermax
  double precision :: levelpop(1:nlev),gast,eps,dendivk,nrdens,dummy,gratio
  double precision :: collpartner_numden(1:npartnermax)
  double precision :: jbar(1:nlin),beta(1:nlin)
  integer :: indx(nlev),error
  double precision :: const1,const2
  parameter(const1=hh/(4*pi))
  parameter(const2=cc*cc/(2*hh))
  !
  ! First, obtain the collisional rates for this molecule, at this
  ! temperature, for all the collision partners of this molecule
  !
  call lines_compute_collrates(ispec,nlev,nlin,npartner,npartnermax,   &
                               gast,collpartner_numden)
  !
  ! For the moment we simply put the Jbar to zero. But in the future
  ! it should be possible to put this to some value, for instance the
  ! radiation from a nearby star or the interstellar radiation field.
  ! Exactly how this is to be done is still being considered.
  !
  jbar(:) = 0.d0
  !
  ! We put the beta to 1.d0 everywhere
  !
  beta(:) = 1.d0
  !
  ! Solve the statistical equilibrium equations
  !
  call lines_solve_statequil(ispec,nlev,nlin,jbar,beta, &
       lines_collisions_npartners_used(ispec),          &
       lines_collisions_ntransmax,                      &
       lines_collisions_ntrans(:,ispec),                &
       lines_col_iup,lines_col_idown,                   &
       lines_col_rateud,lines_col_ratedu,               &
       levelpop,error)
  !
  ! Check for error
  !
  if(error.ne.0) then
     write(stdo,*) 'ERROR in solving statistical equilibrium:'
     if(error.eq.1) then
        write(stdo,*) '      Negative populations found'
     elseif(error.eq.2) then
        write(stdo,*) '      Singular matrix found'
     elseif(error.eq.3) then
        write(stdo,*) '      Normalization of populations bad'
     endif
     write(stdo,*) 'Normalized level populations are:'
     write(stdo,*) levelpop
     stop
  endif
  !
  ! Normalize
  !
  levelpop = nrdens*levelpop
  !
  ! Done...
  !
end subroutine lines_compute_optthinpop


!-------------------------------------------------------------------
!                   WRITE THE LEVEL POPULATIONS
!
! NOTE: These are the full level populations, i.e. in units of
!       1/cm^3. The molecular number density is already included
!       in these numbers. 
!
! Thanks to Attila Juhasz for a bugfix in October 2011
!-------------------------------------------------------------------
subroutine write_levelpop()
  implicit none
  character*80 :: filename,species,strnlev,form
  integer :: ispec,ierr,iformat,nlevels,ilevel,i,precis
  integer(kind=8) :: iiformat,reclen,reclend,nn,kk
  integer(kind=8), allocatable :: subset(:)
  !
  ! Check if the level populations are there
  !
  if(.not.allocated(lines_levelpop)) then
     write(stdo,*) 'WARNING: Wanted to write the level populations. But the level '
     write(stdo,*) '         population array is not allocated... Skipping this.'
     if(lines_mode.lt.0) then
        write(stdo,*) '         Hint: The lines_mode is negative, which means that '
        write(stdo,*) '               the level populations are calculated on-the-fly'
        write(stdo,*) '               during image/spectrum-making; they are not stored.'
        write(stdo,*) '               Use positive lines_mode value instead.'
     endif
     call flush(stdo)
     return
  endif
  !
  ! Message
  !
  write(stdo,*) 'Writing level populations...'
  call flush(stdo)
  !
  ! Determine the precision
  !
  if(rto_single) then
     precis = 4
  else
     precis = 8
  endif
  !
  ! Loop over all species
  !
  do ispec=1,lines_nr_species
     if(lines_species_fullmolec(ispec)) then
        if(igrid_type.lt.100) then
           !
           ! Regular (AMR) grid
           ! 
           ! Just make sure that the cell list is complete
           !
           if(amr_tree_present) then
              call amr_compute_list_all()
           endif
           !
           ! Do a stupidity check
           !
           if(nrcells.ne.amr_nrleafs) stop 3209
           !
           ! Open file and write header
           !
           if(rto_style.eq.1) then
              !
              ! Ascii output
              !
              species  = lines_speciesname(ispec)
              filename = 'levelpop_'//species(1:len_trim(species))//'.dat'
              open(unit=1,file=filename)
              !
              ! Write format number
              !
              write(1,*) 1
              !
              ! Write number of grid points
              !
              write(1,*) nrcellsinp
              !
              ! write number of levels 
              !
              nlevels = lines_nrlevels_subset(ispec)
              write(1,*) nlevels
              !
              ! Write the list of subset levels
              !
              write(1,*) (lines_levels_subset(ilevel,ispec),ilevel=1,nlevels)
           elseif(rto_style.eq.2) then
              !
              ! F77-unformatted
              !
              species  = lines_speciesname(ispec)
              filename = 'levelpop_'//species(1:len_trim(species))//'.udat'
              reclend = rto_reclen/8
              open(unit=1,file=filename,form='unformatted')
              !
              ! Header
              !
              nn=1
              kk=rto_reclen
              write(1) nn,kk
              nn=nrcellsinp
              write(1) nn
              !
              ! Write the nr of levels in the  subset 
              !
              nlevels = lines_nrlevels_subset(ispec)
              nn=nlevels
              write(1) nn
              !
              ! Write the subset of levels
              !
              allocate(subset(1:nlevels),STAT=ierr)
              if(ierr.ne.0) then
                 write(stdo,*) 'ERROR in line module: could not allocate subset()'
                 stop
              endif
              subset(1:nlevels) = lines_levels_subset(1:nlevels,ispec)
              write(1) (subset(ilevel),ilevel=1,nlevels)
           elseif(rto_style.eq.3) then
              !
              ! C-compliant binary output
              !
              species  = lines_speciesname(ispec)
              filename = 'levelpop_'//species(1:len_trim(species))//'.bdat'
              open(unit=1,file=filename,status='replace',access='stream')
              !
              ! Header
              !
              nn=1
              kk=precis
              write(1) nn,kk
              nn=nrcellsinp
              write(1) nn
              !
              ! Write the nr of levels in the  subset 
              !
              nlevels = lines_nrlevels_subset(ispec)
              nn=nlevels
              write(1) nn
              !
              ! Write the subset of levels
              !
              allocate(subset(1:nlevels),STAT=ierr)
              if(ierr.ne.0) then
                 write(stdo,*) 'ERROR in line module: could not allocate subset()'
                 stop
              endif
              subset(1:nlevels) = lines_levels_subset(1:nlevels,ispec)
              write(1) (subset(ilevel),ilevel=1,nlevels)
           else
              write(stdo,*) 'ERROR: Do not know rto_style ',rto_style
              stop
           endif
           !
           ! Write the level populations
           !
           call write_vectorfield(1,rto_style,precis,         &
                lines_nrlevels_subset_max,nlevels,nrcellsinp, &
                lines_nr_species,ispec,                       &
                vector1=lines_levelpop)
           !
           ! Close file
           !
           close(1)
        else
           write(stdo,*) 'ERROR: For the moment no other grid type supported as AMR for writing out level populations'
           stop
        endif
     endif
  enddo ! Loop over species
  !
end subroutine write_levelpop


!-------------------------------------------------------------------
!    COMPUTE LTE UPPER AND LOWER POPULATIONS FOR LINELIST LINES
!
! For molecules for which a linelist is given instead of the full
! molecular data, we must compute the populations according to LTE.
! And we only compute the populations for this particular line.
!-------------------------------------------------------------------
subroutine lines_linelist_compute_ltepop(ispec,iline,nrdens, &
                                         gast,psum,ndown,nup)
  implicit none
  integer :: ispec,iline
  double precision :: nrdens,gast,psum,ndown,nup,dummy,dendivk
  !
  dummy   = nrdens / psum
  dendivk = 1.4387519d0 * lines_linelist_eup(iline,ispec) 
  nup     = dummy * exp(-dendivk/gast) * lines_linelist_gup(iline,ispec)
  dendivk = 1.4387519d0 * lines_linelist_edown(iline,ispec) 
  ndown   = dummy * exp(-dendivk/gast) * lines_linelist_gdown(iline,ispec)
  !
end subroutine lines_linelist_compute_ltepop


!-------------------------------------------------------------------
!            READ THE NUMBER DENSITIES OF TE MOLECULE(S)
!-------------------------------------------------------------------
subroutine read_molecules_numberdensities(action)
  implicit none
  integer :: action,ispec,index,i,ierr,precis,idum,reclenn,style
  integer(kind=8) :: iiformat,reclen,nn,kk
  character*80 :: filename1,filename2,filename3,species
  double precision :: dummy
  logical :: fex1,fex2,fex3
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(gas_chemspec_numberdens)) return
  endif
  !
  ! Default
  !
  precis = 8
  !
  ! Message
  !
  write(stdo,*) 'Reading molecular/atomic number densities...'
  call flush(stdo)
  !
  ! Deallocate stuff if required
  !
  if(allocated(gas_chemspec_numberdens)) deallocate(gas_chemspec_numberdens)
  !
  ! Allocate memory for the number density of each of the species.
  !
  ! NOTE: We may in the future save some memory by not allocating this
  !       if none of the species have a number density file. For now
  !       we will not yet optimize to that level, for the sake of simplicity
  !
  allocate(gas_chemspec_numberdens(1:lines_nr_species,1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate the gas_chemspec_numberdens() array'
     stop
  endif
  !
  ! Now loop over all species
  !
  do ispec=1,lines_nr_species
     !
     ! Get the name of the species
     !
     species  = lines_speciesname(ispec)
     !
     ! Check the file names for the number densities of this species
     ! If no such file is available, then there MUST be a file with
     ! the level populations. If a file is available, then we have
     ! the possibility to calculate e.g. LTE populations internally.
     !
     filename1 = 'numberdens_'//species(1:len_trim(species))//'.inp'
     filename2 = 'numberdens_'//species(1:len_trim(species))//'.uinp'
     filename3 = 'numberdens_'//species(1:len_trim(species))//'.binp'
     inquire(file=filename1,exist=fex1)
     inquire(file=filename2,exist=fex2)
     inquire(file=filename3,exist=fex3)
     idum=0
     if(fex1) idum=idum+1
     if(fex2) idum=idum+1
     if(fex3) idum=idum+1
     if(idum.gt.1) then
        write(stdo,*) 'ERROR: Found more than one file numberdens_'//species(1:len_trim(species))//'.*inp'
        stop
     endif
     if(idum.eq.0) then
        write(stdo,*) 'ERROR: Could not find any numberdens_'//species(1:len_trim(species))//'.*inp'
        stop
     endif
     !
     ! Message
     !
     write(stdo,*) 'Reading gas species number densities...'
     call flush(stdo)
     !
     ! Open file and read header
     !
     if(fex1) then
        !
        ! Formatted input for the number density of this species
        !
        style = 1
        open(unit=1,file=filename1,status='old')
        read(1,*) iiformat
        if(iiformat.ne.1) then
           write(stdo,*) 'ERROR: Format number of '//TRIM(filename1)//' is invalid/unknown.'
           stop
        endif
        read(1,*) nn
     elseif(fex2) then
        !
        ! F77-Unformatted input for the number density of this species
        !
        style = 2
        open(unit=1,file=filename2,status='old',form='unformatted')
        read(1) iiformat,reclen
        reclenn=reclen
        if(iiformat.ne.1) then
           write(stdo,*) 'ERROR: Format number of '//TRIM(filename2)//' is invalid/unknown.'
           write(stdo,*) 'Format number = ',iiformat
           write(stdo,*) 'Record length = ',reclen
           stop
        endif
        read(1) nn
     else
        !
        ! C-style binary format
        !
        style = 3
        open(unit=1,file=filename3,status='old',access='stream')
        read(1) iiformat
        if(iiformat.ne.1) then
           write(stdo,*) 'ERROR: Format number of '//TRIM(filename3)//' is invalid/unknown.'
           write(stdo,*) 'Format number = ',iiformat
           stop
        endif
        read(1) nn
        precis = nn
        read(1) nn
     endif
     !
     ! Do some checks
     !
     if(nn.ne.nrcellsinp) then
        write(stdo,*) 'ERROR: '//trim(filename1)//' does not have same number'
        write(stdo,*) '       of cells as the grid.'
        write(stdo,*) nn,nrcellsinp
        stop
     endif
     !
     ! Read the number densities
     !
     call read_scalarfield(1,style,precis,nrcellsinp,lines_nr_species,1, &
                           ispec,1,1d-99,reclenn,                        &
                           scalar1=gas_chemspec_numberdens)
     !
     ! Close file
     !
     close(1)
  enddo
end subroutine read_molecules_numberdensities


!-------------------------------------------------------------------------
!                 WRITE MOLECULE CHEMSPEC NUMBER DENSITY
! Note: this file has .inp instead of .dat, because typically the molecule
! number density is a quantity that is inserted into RADMC-3D rather than 
! written out. The write_molecule_numberdensity() routine is, however, 
! useful for models made internally with the userdef_module.f90 for 
! debugging or model analysis.
!-------------------------------------------------------------------------
subroutine write_molecule_numberdensity(ispec)
  implicit none
  integer(kind=8) :: iiformat,nn,kk
  integer :: index,ierr,i,precis,ispec
  character*160 :: filename
  character*80 :: species
  !
  ! Set the output style and precision
  !
  if(rto_single) then
     precis = 4
  else
     precis = 8
  endif
  !
  ! Get the name of the species
  !
  species  = lines_speciesname(ispec)
  !
  ! Open file and write header
  !
  if(rto_style.eq.1) then
     !
     ! Formatted output
     !
     filename = 'numberdens_'//species(1:len_trim(species))//'.inp'
     open(unit=1,file=filename)
     write(1,*) 1
     write(1,*) nrcellsinp
  elseif(rto_style.eq.2) then
     !
     ! F77-Unformatted output
     !
     filename = 'numberdens_'//species(1:len_trim(species))//'.uinp'
     open(unit=1,file=filename,form='unformatted')
     nn=1
     kk=rto_reclen
     write(1) nn,kk
     nn=nrcellsinp
     write(1) nn
  elseif(rto_style.eq.3) then
     !
     ! C-compliant binary output
     !
     filename = 'numberdens_'//species(1:len_trim(species))//'.binp'
     open(unit=1,file=filename,status='replace',access='stream')
     nn=1
     kk=precis
     write(1) nn,kk
     nn=nrcellsinp
     write(1) nn
  else
     write(stdo,*) 'ERROR: I/O Style ',rto_style,' unknown'
     stop
  endif
  !
  ! Message
  !
  write(stdo,*) 'Writing molecule number densities to ',trim(filename)
  call flush(stdo)
  !
  ! Write the molecule number density
  !
  call write_scalarfield(1,rto_style,precis,nrcellsinp,         &
                         lines_nr_species,1,ispec,1,rto_reclen, &
                         scalar1=gas_chemspec_numberdens)
  !
  ! Close file
  !
  close(1)
  !
end subroutine write_molecule_numberdensity


!-------------------------------------------------------------------
!       READ THE NUMBER DENSITIES OF THE COLLISIONAL PARTNERS
!-------------------------------------------------------------------
subroutine read_collpartner_numberdensities(action)
  implicit none
  integer :: action,index,i,ierr,reclenn,precis,style
  integer :: ipartner,idum
  integer(kind=8) :: iiformat,reclen,nn,kk
  character*80 :: filename1,filename2,filename3,species
  real :: dummy
  logical :: fex1,fex2,fex3
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(collpartner_numberdens)) return
  endif
  !
  ! Default
  !
  precis = 8
  !
  ! Message
  !
  write(stdo,*) 'Reading collisional partner number densities...'
  call flush(stdo)
  !
  ! If no collision partners, return with a warning
  !
  if(lines_collisions_nr_partner_names.eq.0) then
     write(stdo,*) 'WARNING: Wanting to read collision partner number densities'
     write(stdo,*) '         (presumably because you want to do non-LTE line transfer)'
     write(stdo,*) '         But no collision partners have been specified in the'
     write(stdo,*) '         lines.inp file.'
     return
  endif
  !
  ! Deallocate stuff if required
  !
  if(allocated(collpartner_numberdens)) deallocate(collpartner_numberdens)
  !
  ! Allocate memory for the number density of each of the collisional partner.
  ! TODO: Currently only reads in densities for 1 collisional partner
  !
  ! NOTE: We may in the future save some memory by not allocating this
  !       if none of the collisional partners have a number density file. For now
  !       we will not yet optimize to that level, for the sake of simplicity
  !
  allocate(collpartner_numberdens(1:lines_collisions_nr_partner_names,1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate the collpartner_numberdens() array'
     stop
  endif
  !
  ! Loop over collision partners
  !
  do ipartner=1,lines_collisions_nr_partner_names
     !
     ! Get the file name
     !
     species   = lines_collisions_partner_names(ipartner)
     filename1 = 'numberdens_'//species(1:len_trim(species))//'.inp'
     filename2 = 'numberdens_'//species(1:len_trim(species))//'.uinp'
     filename3 = 'numberdens_'//species(1:len_trim(species))//'.binp'
     inquire(file=filename1,exist=fex1)
     inquire(file=filename2,exist=fex2)
     inquire(file=filename3,exist=fex3)
     idum=0
     if(fex1) idum=idum+1
     if(fex2) idum=idum+1
     if(fex3) idum=idum+1
     if(idum.gt.1) then
        write(stdo,*) 'ERROR: Found more than one file numberdens_'//species(1:len_trim(species))//'.*inp'
        stop
     endif
     if(idum.eq.0) then
        write(stdo,*) 'ERROR: Could not find any numberdens_'//species(1:len_trim(species))//'.*inp'
        stop
     endif
     !
     ! Open file and read header
     !
     if(fex1) then
        !
        ! Formatted input for the number density of this species
        !
        style = 1
        open(unit=1,file=filename1,status='old')
        read(1,*) iiformat
        if(iiformat.ne.1) then
           write(stdo,*) 'ERROR: Format number of '//TRIM(filename1)//' is invalid/unknown.'
           stop
        endif
        read(1,*) nn
     elseif(fex2) then
        !
        ! F77-Unformatted input for the number density of this species
        !
        style = 2
        open(unit=1,file=filename2,status='old',form='unformatted')
        read(1) iiformat,reclen
        reclenn=reclen
        if(iiformat.ne.1) then
           write(stdo,*) 'ERROR: Format number of '//TRIM(filename2)//' is invalid/unknown.'
           write(stdo,*) 'Format number = ',iiformat
           write(stdo,*) 'Record length = ',reclen
           stop
        endif
        read(1) nn
     else
        !
        ! C-style binary format
        !
        style = 3
        open(unit=1,file=filename3,status='old',access='stream')
        read(1) iiformat
        if(iiformat.ne.1) then
           write(stdo,*) 'ERROR: Format number of '//TRIM(filename3)//' is invalid/unknown.'
           write(stdo,*) 'Format number = ',iiformat
           stop
        endif
        read(1) nn
        precis = nn
        read(1) nn
     endif
     !
     ! Do some checks
     !
     if(nn.ne.nrcellsinp) then
        write(stdo,*) 'ERROR: '//trim(filename1)//' does not have same number'
        write(stdo,*) '       of cells as the grid.'
        write(stdo,*) nn,nrcellsinp
        stop
     endif
     !
     ! Read the number densities
     !
     call read_scalarfield(1,style,precis,nrcellsinp,           &
                           lines_collisions_nr_partner_names,1, &
                           ipartner,1,1d-99,reclenn,            &
                           scalar1=collpartner_numberdens)
     !
     ! Close file
     !
     close(1)
  enddo
end subroutine read_collpartner_numberdensities



!-------------------------------------------------------------------
!                    READ THE VELOCITY FIELD
!-------------------------------------------------------------------
subroutine read_velocityfield(action)
  implicit none
  character*80 :: filename1,filename2,filename3
  integer :: ispec,ierr,index,action,icell,i,idum,precis,style,reclen
  integer(kind=8) :: iiformat,nn,kk,reclen8
  logical :: fex1,fex2,fex3
  double precision :: vel(1:3)
  double precision :: dummy
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(gasvelocity)) then
        !
        ! The gas velocity is already determined. It might be
        ! that also the lines_maxveloc has been calculated from
        ! this, but we don't know for sure. So let's calculate
        ! it here to make sure lines_maxveloc is up-to-date. It
        ! can be important if the gas velocity is created 
        ! using the userdef module.
        !
        ! Thanks, Rainer Rolffs, for pointing this out!
        ! 02.01.2011
        !
        lines_maxveloc = 0.d0
        do icell=1,nrcells
           index = cellindex(icell)
           dummy = gasvelocity(1,index)**2 + gasvelocity(2,index)**2 + gasvelocity(3,index)**2 
           dummy = sqrt(dummy)/cc
           if(dummy.gt.lines_maxveloc) then
              lines_maxveloc = dummy
           endif
        enddo
        return
     endif
  endif
  !
  ! Default
  !
  precis = 8
  !
  ! Message
  !
  write(stdo,*) 'Reading velocity field...'
  call flush(stdo)
  !
  ! Create the file name and find if the input file is present and if
  ! it is in text format (.inp) or unformatted (.uinp)
  !
  filename1 = 'gas_velocity.inp'
  filename2 = 'gas_velocity.uinp'
  filename3 = 'gas_velocity.binp'
  inquire(file=filename1,exist=fex1)
  inquire(file=filename2,exist=fex2)
  inquire(file=filename3,exist=fex3)
  idum=0
  if(fex1) idum=idum+1
  if(fex2) idum=idum+1
  if(fex3) idum=idum+1
  if(idum.gt.1) then
     write(stdo,*) 'ERROR: Found more than one file gas_velocity.*inp'
     stop
  endif
  !
  ! If no gas velocity file is present, then we assume that the
  ! gas velocity is 0 everywhere.
  !
  if(idum.eq.0) then
     write(stdo,*) 'Could not find any gas_velocity.*inp. Assuming zero gas velocity.'
     if(allocated(gasvelocity)) deallocate(gasvelocity)
     allocate(gasvelocity(1:3,1:nrcells),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate the gasvelocity() array'
        stop
     endif
     gasvelocity(:,:) = 0.d0
     lines_maxveloc = 0.d0
     return
  endif
  !
  ! Switch formatted/unformatted
  !
  if(fex1) then
     !
     ! Open formatted ascii style
     !
     style = 1
     open(unit=1,file=filename1)
     !
     ! Read format number
     !
     read(1,*) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of '//TRIM(filename1)//' is invalid/unknown.'
        write(stdo,*) 'Format number = ',iiformat
        stop
     endif
     !
     ! Read number of grid points
     !
     read(1,*) nn
  elseif(fex2) then
     !
     ! Open f77-style unformatted (with records)
     !
     style = 2
     open(unit=1,file=filename2,form='unformatted')
     !
     ! Read format number
     !
     read(1) iiformat,reclen8
     reclen = reclen8
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of '//TRIM(filename2)//' is invalid/unknown.'
        write(stdo,*) 'Format number = ',iiformat
        stop
     endif
     read(1) nn
  else
     !
     ! Open C-compliant binary
     !
     style = 3
     open(unit=1,file=filename3,status='old',access='stream')
     !
     ! Read format number
     !
     read(1) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of '//TRIM(filename3)//' is invalid/unknown.'
        write(stdo,*) 'Format number = ',iiformat
        stop
     endif
     read(1) nn
     precis = nn
     read(1) nn
  endif
  !
  ! Do some checks
  !
  if(nn.ne.nrcellsinp) then
     write(stdo,*) 'ERROR: gas_velocity.*inp does not have same number'
     write(stdo,*) '       of cells as the grid.'
     write(stdo,*) nn,nrcellsinp
     stop
  endif
  !
  ! Create the gas velocity arrays
  !
  if(allocated(gasvelocity)) deallocate(gasvelocity)
  allocate(gasvelocity(1:3,1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate the gasvelocity() array'
     stop
  endif
  !
  ! Now read the velocity field
  !
  call read_vectorfield(1,style,precis,3,3,nrcellsinp,1,1,reclen=reclen, &
                        vector0=gasvelocity)
  !
  ! Close file
  !
  close(1)
  !
  ! Determine the maximum absolute velocity. This is useful for when we
  ! have to figure out which lines may contribute to a certain wavelength.
  !
  lines_maxveloc = 0.d0
  do icell=1,nrcells
     index = cellindex(icell)
     dummy = gasvelocity(1,index)**2 + gasvelocity(2,index)**2 + gasvelocity(3,index)**2 
     dummy = sqrt(dummy)/cc
     if(dummy.gt.lines_maxveloc) then
        lines_maxveloc = dummy
     endif
  enddo
  !
end subroutine read_velocityfield


!-------------------------------------------------------------------
!       READ THE LOCAL UNRESOLVED MICRO TURBULENCE (OPTIONAL)
!-------------------------------------------------------------------
subroutine read_microturbulence(action)
  implicit none
  character*80 :: filename1,filename2,filename3
  integer(kind=8) :: iiformat,reclen,reclend,nn,kk
  integer :: ierr,icell,index,action,i,idum,precis,reclenn,style
  real, allocatable :: datas(:)
  double precision, allocatable :: data(:)
  logical :: fex1,fex2,fex3
  real :: dummy
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(lines_microturb)) then
        !
        ! The microturbulence is already determined. It might be
        ! that also the lines_maxturbc has been calculated from
        ! this, but we don't know for sure. So let's calculate
        ! it here to make sure lines_maxturbc is up-to-date. It
        ! can be important if the microturbulence is created 
        ! using the userdef module.
        !
        ! Thanks, Rainer Rolffs, for pointing this out!
        ! 02.01.2011
        !
        lines_maxturbc = 0.d0
        do icell=1,nrcells
           index = cellindex(icell)
           if(lines_microturb(index).gt.lines_maxturbc) then
              lines_maxturbc = lines_microturb(index)
           endif
        enddo
        lines_maxturbc = lines_maxturbc / cc
        return
     endif
  endif
  !
  ! Default
  !
  precis = 8
  !
  ! Deallocate any possible present array
  !
  if(allocated(lines_microturb)) deallocate(lines_microturb)
  !
  ! Create the file name and find if the input file is present and if
  ! it is in text format (.inp) or unformatted (.uinp)
  !
  filename1 = 'microturbulence.inp'
  filename2 = 'microturbulence.uinp'
  filename3 = 'microturbulence.binp'
  inquire(file=filename1,exist=fex1)
  inquire(file=filename2,exist=fex2)
  inquire(file=filename3,exist=fex3)
  idum=0
  if(fex1) idum=idum+1
  if(fex2) idum=idum+1
  if(fex3) idum=idum+1
  if(idum.gt.1) then
     write(stdo,*) 'ERROR: Found more than one file microturbulence.*inp'
     stop
  endif
  !
  ! If no microturbulence file is present, then do not use microturbulence
  !
  if((.not.fex1).and.(.not.fex2).and.(.not.fex3)) then
     write(stdo,*) 'No microturbulence input file found. Assuming zero microturbulence...'
     return
  endif
  !
  ! Message
  !
  write(stdo,*) 'Reading microturbulence...'
  call flush(stdo)
  !
  ! Switch formatted/unformatted
  !
  if(fex1) then
     !
     ! Formatted ascii style
     !
     style = 1
     open(unit=1,file=filename1)
     !
     ! Read format number
     !
     read(1,*) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of '//TRIM(filename1)//' is invalid/unknown.'
        write(stdo,*) 'Format number = ',iiformat
        write(stdo,*) 'Record length = ',reclen
        stop
     endif
     !
     ! Read number of grid points
     !
     read(1,*) nn
  elseif(fex2) then
     !
     ! F77-style unformatted (records)
     !
     style = 2
     open(unit=1,file=filename2,status='old',form='unformatted')
     read(1) iiformat,reclen
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of '//TRIM(filename2)//' is invalid/unknown.'
        write(stdo,*) 'Format number = ',iiformat
        write(stdo,*) 'Record length = ',reclen
        stop
     endif
     reclenn=reclen
     read(1) nn
  else
     !
     ! C-compliant binary style
     !
     style = 3
     open(unit=1,file=filename3,status='old',access='stream')
     read(1) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of '//TRIM(filename3)//' is invalid/unknown.'
        write(stdo,*) 'Format number = ',iiformat
        stop
     endif
     read(1) nn
     precis = nn
     read(1) nn
  endif
  !
  ! Do some checks
  !
  if(nn.ne.nrcellsinp) then
     write(stdo,*) 'ERROR: microturbulence.*inp does not have same number'
     write(stdo,*) '       of cells as the grid.'
     write(stdo,*) nn,nrcellsinp
     stop
  endif
  !
  ! Create the microturbulence arrays
  !
  if(allocated(lines_microturb)) deallocate(lines_microturb)
  allocate(lines_microturb(1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate the lines_microturb() array'
     stop
  endif
  !
  ! Read the microturbulence
  !
  call read_scalarfield(1,style,precis,nrcellsinp,1,1,1,1,1d-99,reclenn,&
                        scalar0=lines_microturb)
  !
  ! Close
  !
  close(1)
  !
  ! Determine the maximum absolute velocity. This is useful for when we
  ! have to figure out which lines may contribute to a certain wavelength.
  !
  lines_maxturbc = 0.d0
  do icell=1,nrcells
     index = cellindex(icell)
     if(lines_microturb(index).gt.lines_maxturbc) then
        lines_maxturbc = lines_microturb(index)
     endif
  enddo
  lines_maxturbc = lines_maxturbc / cc
  !
end subroutine read_microturbulence


!-------------------------------------------------------------------
!       READ THE LENGTH SCALE FOR ESCAPE PROBABILITY (OPTIONAL)
!
! This file must be read if escape probability is desired as an
! option to e.g. the LVG populations.
!-------------------------------------------------------------------
subroutine read_escprob_lengthscale(action)
  implicit none
  character*80 :: filename1,filename2,filename3
  integer(kind=8) :: iiformat,reclen,reclend,nn,kk
  integer :: ierr,icell,index,action,i,idum,precis,reclenn,style
  real, allocatable :: datas(:)
  double precision, allocatable :: data(:)
  logical :: fex1,fex2,fex3
  real :: dummy
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(lines_escprob_lengthscale)) then
        return
     endif
  endif
  !
  ! Default
  !
  precis = 8
  !
  ! Deallocate any possible present array
  !
  if(allocated(lines_escprob_lengthscale)) deallocate(lines_escprob_lengthscale)
  !
  ! Create the file name and find if the input file is present and if
  ! it is in text format (.inp) or unformatted (.uinp)
  !
  filename1 = 'escprob_lengthscale.inp'
  filename2 = 'escprob_lengthscale.uinp'
  filename3 = 'escprob_lengthscale.binp'
  inquire(file=filename1,exist=fex1)
  inquire(file=filename2,exist=fex2)
  inquire(file=filename3,exist=fex3)
  idum=0
  if(fex1) idum=idum+1
  if(fex2) idum=idum+1
  if(fex3) idum=idum+1
  if(idum.gt.1) then
     write(stdo,*) 'ERROR: Found more than one file escprob_lengthscale.*inp'
     stop
  endif
  !
  ! If no length scale file is present, then do not use escape probability
  !
  if((.not.fex1).and.(.not.fex2).and.(.not.fex3)) then
     write(stdo,*) 'No escprob_lengthscale input file found.'
     write(stdo,*) '   ---> Assuming infinite length scale (no escape probability)...'
     return
  endif
  !
  ! Message
  !
  write(stdo,*) 'Reading escape probability length scale file...'
  call flush(stdo)
  !
  ! Switch formatted/unformatted
  !
  if(fex1) then
     !
     ! Formatted ascii style
     !
     style = 1
     open(unit=1,file=filename1)
     !
     ! Read format number
     !
     read(1,*) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of '//TRIM(filename1)//' is invalid/unknown.'
        write(stdo,*) 'Format number = ',iiformat
        write(stdo,*) 'Record length = ',reclen
        stop
     endif
     !
     ! Read number of grid points
     !
     read(1,*) nn
  elseif(fex2) then
     !
     ! F77-style unformatted (records)
     !
     style = 2
     open(unit=1,file=filename2,status='old',form='unformatted')
     read(1) iiformat,reclen
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of '//TRIM(filename2)//' is invalid/unknown.'
        write(stdo,*) 'Format number = ',iiformat
        write(stdo,*) 'Record length = ',reclen
        stop
     endif
     reclenn=reclen
     read(1) nn
  else
     !
     ! C-compliant binary style
     !
     style = 3
     open(unit=1,file=filename3,status='old',access='stream')
     read(1) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of '//TRIM(filename3)//' is invalid/unknown.'
        write(stdo,*) 'Format number = ',iiformat
        stop
     endif
     read(1) nn
     precis = nn
     read(1) nn
  endif
  !
  ! Do some checks
  !
  if(nn.ne.nrcellsinp) then
     write(stdo,*) 'ERROR: escprob_lengthscale.*inp does not have same number'
     write(stdo,*) '       of cells as the grid.'
     write(stdo,*) nn,nrcellsinp
     stop
  endif
  !
  ! Create the escape probability length scale arrays
  !
  if(allocated(lines_escprob_lengthscale)) deallocate(lines_escprob_lengthscale)
  allocate(lines_escprob_lengthscale(1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate the lines_escprob_lengthscale() array'
     stop
  endif
  !
  ! Read the escprob_lengthscale
  !
  call read_scalarfield(1,style,precis,nrcellsinp,1,1,1,1,1d-99,reclenn,&
                        scalar0=lines_escprob_lengthscale)
  !
  ! Close
  !
  close(1)
  !
end subroutine read_escprob_lengthscale


!-------------------------------------------------------------------
!                  COMPUTE THE MAXIMUM LINE SHIFT
!
! (New 2016.07.23)
! In the past the calculations for computing the maximum line shift
! were spread over the various routines for reading in the various
! line-relevant data. But the disadvantage of that is that if we 
! want to include lines in userdef_module.f90 then we must do that
! all by hand. Here is a subroutine that does that all, so that 
! it is merely one command. This is only necessary to be called
! if you use lines in the userdef_module.f90, not for the normal
! lines where all data are read via files.
!-------------------------------------------------------------------
subroutine lines_compute_maxlineshift()
  implicit none
  doubleprecision :: dummy
  integer :: index,icell
  !
  gastmax = 0.d0
  if(.not.allocated(gastemp)) then
     write(stdo,*) 'ERROR: Gas temperature must be set for lines'
     stop
  endif
  do icell=1,nrcells
     index = cellindex(icell)
     if(gastemp(index).gt.gastmax) then
        gastmax = gastemp(index)
     endif
  enddo
  lines_maxturbc = 0.d0
  if(.not.allocated(lines_microturb)) then
     write(stdo,*) 'ERROR: lines_microturb must be set for lines'
     stop
  endif
  do icell=1,nrcells
     index = cellindex(icell)
     if(lines_microturb(index).gt.lines_maxturbc) then
        lines_maxturbc = lines_microturb(index)
     endif
  enddo
  lines_maxveloc = 0.d0
  if(.not.allocated(gasvelocity)) then
     write(stdo,*) 'ERROR: gasvelocity must be set for lines'
     stop
  endif
  do icell=1,nrcells
     index = cellindex(icell)
     dummy = gasvelocity(1,index)**2 + gasvelocity(2,index)**2 + gasvelocity(3,index)**2 
     dummy = sqrt(dummy)/cc
     if(dummy.gt.lines_maxveloc) then
        lines_maxveloc = dummy
     endif
  enddo
  lines_maxturbc = lines_maxturbc / cc
  lines_maxtempc = sqrt(2*kk*gastmax/(lines_umass_min))/cc
  lines_maxshift = lines_widthmargin * ( lines_maxtempc + lines_maxturbc ) + &
                   lines_maxveloc
end subroutine lines_compute_maxlineshift



!-------------------------------------------------------------------
!           COMPUTE THE MAXIMUM RATIO OF DOPPLER SHIFT 
!                  BETWEEN CELLS VS LINEWIDTH
!
! (New 2016.07.26)
! It remains often hard to figure out whether the grid resolution
! is fine enough to avoid cell-by-cell doppler shifts that are 
! larger than the local line width. If such large doppler shifts
! occur, the doppler-catching method must be used for making 
! images and spectra.
!-------------------------------------------------------------------
subroutine lines_compute_maxrellineshift()
  implicit none
  doubleprecision :: width_local,relshift,dv
  integer :: index,icell
  !
  ! Reset 
  !
  lines_maxrelshift = 0.d0
  !
  ! Visit all cells
  !
  do icell=1,nrcells
     index = cellindex(icell)
     !
     ! Compute the local line width in cm/s for the most massive molecule
     ! in the current model
     !
     width_local = sqrt(2*kk*gastemp(index)/lines_umass_max+lines_microturb(index)**2)
     !
     ! Compute the maximum velocity difference between this
     ! cell and adjacent cells (maximum of the 2x3 directions)
     !
     call lines_compute_velgradient(index,dv,maxdiff=.true.)
     !
     ! Take the ratio 
     !
     relshift = dv / (width_local+1d-99)
     !
     ! Only find the maximal value
     !
     lines_maxrelshift = max(lines_maxrelshift,relshift)
     !
  enddo
end subroutine lines_compute_maxrellineshift



!============================================================================
!                     THE RADIATIVE TRANSFER ROUTINES
!
!               SOME OF THESE ROUTINES CAN ALSO BE USED FOR
!                OPEN MP SHARED MEMORY PARALLELIZATION AND
!              FOR FUTURE NVIDIA GRAPHICS CARD IMPLEMENTATION
!============================================================================


!----------------------------------------------------------------------------
!                SET THE LINE TRANSFER FREQUENCY WINDOW
!----------------------------------------------------------------------------
subroutine lines_serial_init_raytrace(action)
  implicit none
  integer :: action
  !$OMP PARALLEL
  ray_ns    = 1
  ray_nsmax = 1
  !$OMP END PARALLEL
  call lines_ray1d_init_raytrace(action)
end subroutine lines_serial_init_raytrace



!----------------------------------------------------------------------------
!                SET THE LINE TRANSFER FREQUENCY WINDOW
!----------------------------------------------------------------------------
subroutine lines_ray1d_init_raytrace(action)
  implicit none
  integer :: ierr,action
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(lines_ray_doppler)) return
  endif
  !
  ! Check lines
  !
  if(.not.rt_incl_lines) then
     write(stdo,*) 'ERROR: Trying to initialize lines, but incl_lines=F'
     stop
  endif
  !
  ! First clean up things
  !
  !$OMP PARALLEL
  if(allocated(lines_ray_levpop)) deallocate(lines_ray_levpop)
  if(allocated(lines_ray_nrdens)) deallocate(lines_ray_nrdens)
  if(allocated(lines_ray_temp)) deallocate(lines_ray_temp)
  if(allocated(lines_ray_turb)) deallocate(lines_ray_turb)
  if(allocated(lines_ray_doppler)) deallocate(lines_ray_doppler)
  !$OMP END PARALLEL
  !lines_warn_lineleap = .false.
  !
  ! Check some basic things
  !
  !$OMP PARALLEL
  if(ray_nsmax.lt.1) then
     write(stdo,*) 'ERROR in line module: ray_nsmax not set.'
     stop
  endif
  !$OMP END PARALLEL
  !
  ! Check if the rest is allocated
  !
  if(.not.allocated(lines_nrlines_subset)) then
     write(stdo,*) 'ERROR in line module: must first initalize the line module'
     write(stdo,*) '      before setting up the line ray-tracing.'
     stop
  endif
  !
  ! Allocate the various arrays
  !  
  !$OMP PARALLEL
  allocate(lines_ray_levpop(lines_maxnrlevels,lines_nr_species,ray_nsmax),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_ray_levpop(:,:,:).'
     stop 
  endif
  !$OMP END PARALLEL
  !$OMP PARALLEL
  allocate(lines_ray_nrdens(lines_nr_species,ray_nsmax),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_ray_nrdens(:,:).'
     stop 
  endif
  !$OMP END PARALLEL
  !$OMP PARALLEL
  allocate(lines_ray_temp(ray_nsmax),STAT=ierr) 
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_ray_temp(:).'
     stop 
  endif
  !$OMP END PARALLEL
  !$OMP PARALLEL
  allocate(lines_ray_turb(ray_nsmax),STAT=ierr) 
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_ray_turb(:).'
     stop 
  endif
  !$OMP END PARALLEL
  !$OMP PARALLEL
  allocate(lines_ray_doppler(ray_nsmax),STAT=ierr) 
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_ray_doppler(:).'
     stop 
  endif
  !$OMP END PARALLEL
  !$OMP PARALLEL
  allocate(lines_ray_lorentz_delta(ray_nsmax),STAT=ierr) 
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in lines module: Could not allocate '
     write(stdo,*) '      lines_ray_lorentz_delta(:).'
     stop 
  endif
  !$OMP END PARALLEL
  !
end subroutine lines_ray1d_init_raytrace



!----------------------------------------------------------------------------
!        FIND ACTIVE LINES FOR GIVEN WAVELENGTH GRID (linelist mode)
!
! To speed up image and spectrum ray-tracing for the case when we have 
! many lines it is important to be able to identify which lines are going to
! contribute to the spectrum or image. This is what this subroutine does.
! You can ask it to check the full spectral range given by the frequencies
! array, or a part (from inu0 to inu1). If you want to do the full, then
! put inu0=1 and inu1=nrfreq. 
!
! Note that it will only scan for lines for which the levels are in the
! level subset. The subset is a user-specified subset of levels of a molecule
! for which you wish to include the lines. Choosing a small subset has the
! advantage that, if you compute the levels before the ray-tracing (i.e.
! if lines_mode is .gt. 0), the big array of populations is not too big,
! since only the populations for the subset levels will be stored. By
! default, however, all lines are included in the subset.
!
!----------------------------------------------------------------------------
subroutine lines_find_active_lines_linelist(nrfreq,&
                        frequencies,inu0,inu1,maxshift)
  implicit none
  integer :: nrfreq,inu0,inu1,inu,ilinesub,ispec,iline,ilevel,ilevelsub
  double precision :: frequencies(1:nrfreq),maxshift,dummy
  logical :: notdone
  !
  if(debug_check_all.ne.0) then
     ! -------------------------------------------------------------------------------
     !  This is not needed in linelist mode - Attila Juhasz 
     ! -------------------------------------------------------------------------------
     !     if((.not.allocated(active_levels)).or.&
     !        (.not.allocated(active_levels_subsetindex)).or.&
     !        (.not.allocated(active_lines)).or.&
     !        (.not.allocated(active_dummy))) then
     ! -------------------------------------------------------------------------------
     if(.not.allocated(active_lines)) then
        write(stdo,*) 'ERROR in line module: active lines not allocated.'
        stop
     endif
  endif
  !
  ! Loop over all molecular species
  !
  do ispec=1,lines_nr_species
     ! -------------------------------------------------------------------------------
     !  No need for this in linelist mode - Attila Juhasz
     ! -------------------------------------------------------------------------------
     !     do ilevel=1,lines_maxnrlevels
     !        active_dummy(ilevel) = 0
     !     enddo
     ! -------------------------------------------------------------------------------
     active_nrlines(ispec) = 0
     !
     ! BUGFIX 22.11.09: Swapped do ilinesub and do inu, and added notdone flag
     !
     ! Make a loop over all lines for which both levels are in the level subset
     !
     ! -------------------------------------------------------------------------------
     !  No need for this in linelist mode - Attila Juhasz
     ! lines_nrlines_subset(ispec) is zero in linelist mode! 
     ! -------------------------------------------------------------------------------
     !
     !     do ilinesub=1,lines_nrlines_subset(ispec)
     ! -------------------------------------------------------------------------------
     do ilinesub=1,lines_nrlines(ispec)
        !
        ! Get the actual line ID
        !
        ! -------------------------------------------------------------------------------
        !  No need for this in linelist mode - Attila Juhasz
        ! lines_lines_subset(:,:) is not set (=empty) in linelist mode
        ! -------------------------------------------------------------------------------
        !iline = lines_lines_subset(ilinesub,ispec)
        ! -------------------------------------------------------------------------------
        iline = ilinesub
        !
        ! Scan the frequency grid between inu0 and inu1 to 
        ! check if this level must be included
        !
        notdone = .true.
        do inu=inu0,inu1
           if(notdone) then
              !
              ! Check the distance of this line to the 
              ! current frequency
              !
              dummy = abs(1.d0-frequencies(inu)/lines_nu0(iline,ispec))
              if(dummy.lt.maxshift) then
                 !
                 ! Yes, this line is within range, and has to be
                 ! included.
                 !
                 active_nrlines(ispec) = active_nrlines(ispec) + 1
                 active_lines(active_nrlines(ispec),ispec) = iline
                 ! -------------------------------------------------------------------------------
                 !  No need for this in linelist mode - Attila Juhasz
                 ! -------------------------------------------------------------------------------
                 ! 
                 ! write(stdo,*)  '   In image: including Molecule', ispec, '  Transition',iline
                 !                 ilevel=lines_levelup(iline,ispec)
                 !                 active_dummy(ilevel) = 1
                 !                 ilevel=lines_leveldown(iline,ispec)
                 !                 active_dummy(ilevel) = 1
                 ! -------------------------------------------------------------------------------
                 notdone = .false.
              endif
           endif
        enddo
     enddo
     ! -------------------------------------------------------------------------------
     !  No need for this in linelist mode - Attila Juhasz
     ! -------------------------------------------------------------------------------
     !
     ! Now that we know which lines are active, we must
     ! find out which of the levels belong to these.
     !
     !    active_nrlevels(ispec) = 0
     !    do ilevelsub=1,lines_nrlevels_subset(ispec)
     !       ilevel = lines_levels_subset(ilevelsub,ispec)
     !       if(active_dummy(ilevel).eq.1) then
     !          active_nrlevels(ispec) = active_nrlevels(ispec) + 1
     !          active_levels(active_nrlevels(ispec),ispec) = ilevel
     !          active_levels_subsetindex(active_nrlevels(ispec),ispec) = ilevelsub
     !        endif
     !enddo
     ! -------------------------------------------------------------------------------
  enddo
end subroutine lines_find_active_lines_linelist



!----------------------------------------------------------------------------
!        FIND ACTIVE LINES AND LEVELS FOR GIVEN WAVELENGTH GRID
!
! To speed up image and spectrum ray-tracing for the case when we have 
! many lines it is important to be able to identify which lines are going to
! contribute to the spectrum or image. This is what this subroutine does.
! You can ask it to check the full spectral range given by the frequencies
! array, or a part (from inu0 to inu1). If you want to do the full, then
! put inu0=1 and inu1=nrfreq. 
!
! Note that it will only scan for lines for which the levels are in the
! level subset. The subset is a user-specified subset of levels of a molecule
! for which you wish to include the lines. Choosing a small subset has the
! advantage that, if you compute the levels before the ray-tracing (i.e.
! if lines_mode is .gt. 0), the big array of populations is not too big,
! since only the populations for the subset levels will be stored. By
! default, however, all lines are included in the subset.
!----------------------------------------------------------------------------
subroutine lines_find_active_lines_levels(nrfreq,&
                        frequencies,inu0,inu1,maxshift)
  implicit none
  integer :: nrfreq,inu0,inu1,inu,ilinesub,ispec,iline,ilevel,ilevelsub
  double precision :: frequencies(1:nrfreq),maxshift,dummy
  logical :: notdone
  !
  if(debug_check_all.ne.0) then
     if((.not.allocated(active_levels)).or.&
        (.not.allocated(active_levels_subsetindex)).or.&
        (.not.allocated(active_lines)).or.&
        (.not.allocated(active_dummy))) then
        write(stdo,*) 'ERROR in line module: active lines/levels/dummy not allocated.'
        stop
     endif
  endif
  !
  ! Loop over all molecular species
  !
  do ispec=1,lines_nr_species
     do ilevel=1,lines_maxnrlevels
        active_dummy(ilevel) = 0
     enddo
     active_nrlines(ispec) = 0
     !
     ! Make a loop over all lines 
     !
     ilinesub=1
     do iline=1,lines_nrlines(ispec)
        !
        ! Scan the frequency grid between inu0 and inu1 to 
        ! check if this level must be included
        !
        notdone = .true.
        do inu=inu0,inu1
           if(notdone) then
              !
              ! Check the distance of this line to the 
              ! current frequency
              !
              dummy = abs(1.d0-frequencies(inu)/lines_nu0(iline,ispec))
              if(dummy.lt.maxshift) then
                 !
                 ! Yes, this line is within range, and has to be
                 ! included.
                 !
                 notdone = .false.
                 !
                 ! Find the line in the subset
                 !
                 do while((lines_lines_subset(ilinesub,ispec).lt.iline).and.&
                      (ilinesub.lt.lines_nrlines_subset(ispec)))
                    ilinesub = ilinesub + 1
                 enddo
                 if(lines_lines_subset(ilinesub,ispec).eq.iline) then
                    !
                    ! Yes, this line is in the subset
                    !
                    active_nrlines(ispec) = active_nrlines(ispec) + 1
                    active_lines(active_nrlines(ispec),ispec) = iline
                    ! write(stdo,*)  '   In image: including Molecule', ispec, '  Transition',iline
                    ilevel=lines_levelup(iline,ispec)
                    active_dummy(ilevel) = 1
                    ilevel=lines_leveldown(iline,ispec)
                    active_dummy(ilevel) = 1
                 else
                    !
                    ! No, this line is not in the subset, so we must
                    ! warn, because the line will not appear in the
                    ! spectrum even though it should (physically) appear.
                    !
                    write(stdo,*) 'WARNING: imol=',ispec,' iline=',iline,' in wavelength range, but not in subset!'
                 endif
              endif
           endif
        enddo
     enddo
     !
     ! Now that we know which lines are active, we must
     ! find out which of the levels belong to these.
     !
     active_nrlevels(ispec) = 0
     do ilevelsub=1,lines_nrlevels_subset(ispec)
        ilevel = lines_levels_subset(ilevelsub,ispec)
        if(active_dummy(ilevel).eq.1) then
           active_nrlevels(ispec) = active_nrlevels(ispec) + 1
           active_levels(active_nrlevels(ispec),ispec) = ilevel
           active_levels_subsetindex(active_nrlevels(ispec),ispec) = ilevelsub
        endif
     enddo
  enddo
end subroutine lines_find_active_lines_levels


!----------------------------------------------------------------------------
!           AUTOMATICALLY SELECT THE RELEVANT LEVEL SUBSET
!
! This subroutine is, in its functionality, very similar to the 
! lines_find_active_lines_levels() subroutine. But instead of choosing from
! the subset of levels, it chooses from ALL levels of all molecules, and
! also sets the subset selection equal to that set of levels. The active
! levels is then equal to the subset. This reason for doing this is to 
! allow RADMC-3D to automatically find the best subset of levels, so that
! the user does not have to do this him/herself by hand. You just give the
! wavelength range in which you are interested (inu0,inu1) and this routine
! will automatically select the subset. NOTE: You can always, after calling
! this routine, focus on a smaller frequency range (and thus a smaller
! set of active lines) by calling the lines_find_active_lines_levels() 
! again. The active levels are then simply overwritten. 
!----------------------------------------------------------------------------
subroutine lines_automatic_subset_selection(nrfreq,      &
                        frequencies,inu0,inu1,maxshift,  &
                        redo)
  implicit none
  integer :: nrfreq,inu0,inu1,inu,ilinesub,ispec,iline,ilevel,ilevelsub,ierr
  double precision :: frequencies(1:nrfreq),maxshift,dummy
  logical :: notdone,redo
  integer, allocatable :: dummy_active_lines(:,:)
  integer, allocatable :: dummy_active_nrlines(:)
  integer, allocatable :: dummy_active_levels(:,:)
  integer, allocatable :: dummy_active_levels_subsetindex(:,:)
  integer, allocatable :: dummy_active_nrlevels(:)
  integer, allocatable :: dummy_active_dummy(:)
  !
  ! Set a flag to .false.
  !
  redo = .false.
  !
  ! If the active_lines() array is not allocated we 
  ! have to "redo" everything anyway
  !
  if(.not.allocated(active_lines)) then
     redo = .true.
  endif
  if(.not.allocated(lines_levels_subset)) then
     redo = .true.
  endif
  if(.not.allocated(lines_levelpop)) then
     redo = .true.
  endif
  !
  ! First we find the lines that are contributing to the
  ! frequencies given in the frequencies array. This is
  ! nearly identical to what is done in the 
  ! lines_find_active_lines_levels() subroutine, but
  ! now for all possible lines and levels.
  !
  ! ...First allocate temporary arrays
  !
  allocate(dummy_active_lines(lines_maxnrlines,lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in line module: array dummy_active_lines() not allocatable'
     stop
  endif
  allocate(dummy_active_nrlines(lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in line module: array dummy_active_nrlines() not allocatable'
     stop
  endif
  allocate(dummy_active_levels(lines_maxnrlevels,lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in line module: array dummy_active_levels() not allocatable'
     stop
  endif
  allocate(dummy_active_levels_subsetindex(lines_maxnrlevels,lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in line module: array dummy_active_levels_subsetindex() not allocatable'
     stop
  endif
  allocate(dummy_active_nrlevels(lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in line module: array dummy_active_nrlevels() not allocatable'
     stop
  endif
  allocate(dummy_active_dummy(lines_maxnrlevels),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in line module: array dummy_active_dummy() not allocatable'
     stop
  endif
  !
  ! Now do the subset selection...
  !
  ! Loop over all molecular species
  !
  do ispec=1,lines_nr_species
     !
     ! Do this only if this molecule is fully specified (i.e. is not
     ! given by just a linelist).
     !
     if(lines_species_fullmolec(ispec)) then
        !
        ! Reset some stuff
        !
        do ilevel=1,lines_maxnrlevels
           dummy_active_dummy(ilevel) = 0
        enddo
        dummy_active_nrlines(ispec) = 0
        !
        ! Make a loop over all lines for which both levels are in the level subset
        !
        do iline=1,lines_nrlines(ispec)
           !
           ! Scan the frequency grid between inu0 and inu1 to 
           ! check if this level must be included
           !
           notdone = .true.
           do inu=inu0,inu1
              if(notdone) then
                 !
                 ! Check the distance of this line to the 
                 ! current frequency
                 !
                 dummy = abs(1.d0-frequencies(inu)/lines_nu0(iline,ispec))
                 if(dummy.lt.maxshift) then
                    !
                    ! Yes, this line is within range, and has to be
                    ! included.
                    !
                    dummy_active_nrlines(ispec) = dummy_active_nrlines(ispec) + 1
                    dummy_active_lines(dummy_active_nrlines(ispec),ispec) = iline
                    ilevel=lines_levelup(iline,ispec)
                    dummy_active_dummy(ilevel) = 1
                    ilevel=lines_leveldown(iline,ispec)
                    dummy_active_dummy(ilevel) = 1
                    notdone = .false.
                 endif
              endif
           enddo
        enddo
        !
        ! Now that we know which lines are active, we must
        ! find out which of the levels belong to these.
        ! As the active levels and lines we take all that 
        ! are in the subset.
        !
        dummy_active_nrlevels(ispec) = 0
        do ilevel=1,lines_nrlevels(ispec)
           if(dummy_active_dummy(ilevel).eq.1) then
              dummy_active_nrlevels(ispec) = dummy_active_nrlevels(ispec) + 1
              dummy_active_levels(dummy_active_nrlevels(ispec),ispec) = ilevel
              dummy_active_levels_subsetindex(dummy_active_nrlevels(ispec),ispec) = ilevel
           endif
        enddo
        !
        ! Now check if this selection is identical to the
        ! selection we already had
        !
        if(.not.redo) then
           if(dummy_active_nrlevels(ispec).eq.active_nrlevels(ispec)) then
              do ilevel=1,active_nrlevels(ispec)
                 if(dummy_active_levels(ilevel,ispec).ne.active_levels(ilevel,ispec)) then
                    redo = .true.
                 endif
              enddo
           else
              redo = .true.
           endif
           if(dummy_active_nrlines(ispec).eq.active_nrlines(ispec)) then
              do iline=1,active_nrlines(ispec)
                 if(dummy_active_lines(iline,ispec).ne.active_lines(iline,ispec)) then
                    redo = .true.
                 endif
              enddo
           else
              redo = .true.
           endif
           if(active_nrlines(ispec).ne.lines_nrlines_subset(ispec)) then
              redo = .true.
           endif
           if(active_nrlevels(ispec).ne.lines_nrlevels_subset(ispec)) then
              redo = .true.
           endif
        endif
     endif
  enddo
  !
  ! If by now redo is still .false. then we do not need to 
  ! deallocate and redo everything. Any populations that
  ! may have been computed previously and have been stored in 
  ! lines_levelpop() also remain valid.
  !
  if(redo) then
     !
     ! But if redo is .true., as it is here, then we have to start
     ! all over
     !
     ! Deallocate all subset and active level stuff
     !
     if(allocated(lines_nrlevels_subset)) deallocate(lines_nrlevels_subset)
     if(allocated(lines_nrlines_subset)) deallocate(lines_nrlines_subset)
     if(allocated(lines_levels_subset)) deallocate(lines_levels_subset)
     if(allocated(lines_lines_subset)) deallocate(lines_lines_subset)
     if(allocated(active_lines)) deallocate(active_lines)
     if(allocated(active_nrlines)) deallocate(active_nrlines)
     if(allocated(active_levels)) deallocate(active_levels)
     if(allocated(active_levels_subsetindex)) deallocate(active_levels_subsetindex)
     if(allocated(active_nrlevels)) deallocate(active_nrlevels)
     if(allocated(active_dummy)) deallocate(active_dummy)
     lines_nrlevels_subset_max = 0
     lines_nrlines_subset_max = 0
     !
     ! Also deallocate the big global level population array
     !
     if(allocated(lines_levelpop)) deallocate(lines_levelpop)
     !
     ! Allocate active levels stuff
     !
     allocate(active_lines(lines_maxnrlines,lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in line module: array active_lines() not allocatable'
        stop
     endif
     allocate(active_nrlines(lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in line module: array active_nrlines() not allocatable'
        stop
     endif
     allocate(active_levels(lines_maxnrlevels,lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in line module: array active_levels() not allocatable'
        stop
     endif
     allocate(active_levels_subsetindex(lines_maxnrlevels,lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in line module: array active_levels_subsetindex() not allocatable'
        stop
     endif
     allocate(active_nrlevels(lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in line module: array active_nrlevels() not allocatable'
        stop
     endif
     allocate(active_dummy(lines_maxnrlevels),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in line module: array active_dummy() not allocatable'
        stop
     endif
     !
     ! Allocate the subset stuff
     !
     allocate(lines_nrlevels_subset(1:lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      lines_nrlevels_subset(:).'
        stop 
     endif
     lines_nrlevels_subset(:) = 0
     allocate(lines_nrlines_subset(1:lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      lines_nrlines_subset(:).'
        stop 
     endif
     lines_nrlines_subset(:) = 0
     allocate(lines_levels_subset(1:lines_maxnrlevels,lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      lines_levels_subset(:,:).'
        stop 
     endif
     allocate(lines_lines_subset(1:lines_maxnrlines,lines_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      lines_lines_subset(:,:).'
        stop 
     endif
     !
     ! Since we already did the work, we can just copy
     !
     active_dummy              = dummy_active_dummy
     active_nrlines            = dummy_active_nrlines
     active_lines              = dummy_active_lines
     active_nrlevels           = dummy_active_nrlevels
     active_levels             = dummy_active_levels
     active_levels_subsetindex = dummy_active_levels_subsetindex
     !
     ! Now set the level subset to the same selection
     !
     do ispec=1,lines_nr_species
        lines_nrlevels_subset(ispec) = active_nrlevels(ispec)
        lines_nrlines_subset(ispec)  = active_nrlines(ispec)
        do ilevel=1,lines_nrlevels_subset(ispec)
           lines_levels_subset(ilevel,ispec) = active_levels(ilevel,ispec)
        enddo
        do iline=1,lines_nrlines_subset(ispec)
           lines_lines_subset(iline,ispec) = active_lines(iline,ispec)
        enddo
     enddo
     !
     ! Compute the total number of subset lines
     !
     lines_nrlines_subset_max=0
     do ispec=1,lines_nr_species
        lines_nrlines_subset_max = max(lines_nrlines_subset_max,    &
                                   lines_nrlines_subset(ispec))
     enddo
     lines_nrlevels_subset_max=0
     do ispec=1,lines_nr_species
        lines_nrlevels_subset_max = max(lines_nrlevels_subset_max,    &
                                    lines_nrlevels_subset(ispec))
     enddo
     !
     ! Some messages
     !
     do ispec=1,lines_nr_species
        if(lines_species_fullmolec(ispec)) then
           write(stdo,'(A42,I2,A44)') ' Will store level populations of molecule ',ispec,&
                ' into global array for the following levels:'
           write(stdo,*) (lines_levels_subset(ilevelsub,ispec),ilevelsub=1,lines_nrlevels_subset(ispec))
        endif
     enddo
  endif
  !
  ! Deallocate the temporary arrays
  !
  deallocate(dummy_active_lines)
  deallocate(dummy_active_nrlines)
  deallocate(dummy_active_levels)
  deallocate(dummy_active_levels_subsetindex)
  deallocate(dummy_active_nrlevels)
  deallocate(dummy_active_dummy)
  !
end subroutine lines_automatic_subset_selection



!----------------------------------------------------------------------------
!     COMPUTE THE SOURCE AND OPACITY FROM THE LEVEL POPULATIONS
!                  FOR RAY-TRACING METHOD A (SERIAL)
!
! This routine gives the full j_nu and alpha_nu for the lines. If the ray is
! only 1 element, then this can be considered as a point-evaluation. The
! values will be calculated for a predefined set of frequencies. 
! NOTE: The src and alp of the lines will be ADDED to anything that
!       already is in the src and alp (for instance dust emission!)
!----------------------------------------------------------------------------
subroutine lines_serial_addto_jnu_alpnu(nrfreq,inu0,inu1,freq,src,alp)
  implicit none
  integer :: inu,iline,ilinesub,ierr,ispec,nrfreq,inu0,inu1,itemp
  doubleprecision :: src(1:nrfreq),alp(1:nrfreq),freq(1:nrfreq)
  doubleprecision :: dummy,nup,ndown,eps
  doubleprecision :: bdu,bud,gratio,a,dnu,phi
  doubleprecision :: dum,alpha,x,y
  doubleprecision :: const1,const2
  parameter(const1=hh/(4*pi))
  parameter(const2=cc*cc/(2*hh))
  !
  ! Checks
  !
  if(debug_check_all.ne.0) then
     if((.not.allocated(lines_ray_levpop)).or. &
          (.not.allocated(lines_ray_temp)).or. &
          (.not.allocated(lines_ray_turb)).or. &
          (.not.allocated(lines_ray_doppler))) then
        write(stdo,*) 'ERROR in line module: ray arrays not allocated.'
        stop
     endif
  endif
  !
  ! If there are molecules for which linelists are specified instead
  ! of full molecular data, then we must compute the partition sum
  ! for each of these molecules here.
  !
  if(allocated(lines_linelist_eup)) then
     do ispec=1,lines_nr_species
        if(.not.lines_species_fullmolec(ispec)) then
           !
           ! Find the location of the current temperature in the partition function table
           !
           call hunt(lines_partition_temperature(:,ispec),               &
                     lines_partition_ntemp(ispec),lines_ray_temp(1),itemp)
           if(itemp.ge.lines_partition_ntemp(ispec)) then
              write(stdo,*) 'ERROR: Temperature out of range of partition function '
              write(stdo,*) '       for ispec=',ispec,' Temperature = ',lines_ray_temp(1)
              !write(stdo,*) '       Smallest temperature in table = ',   &
              !     lines_partition_temperature(1,ispec)
              write(stdo,*) '       Largest temperature in table = ',   &
                   lines_partition_temperature(lines_partition_ntemp(ispec),ispec)
              stop
           endif
           if(itemp.ge.1) then
              eps = (lines_ray_temp(1)-lines_partition_temperature(itemp,ispec)) / &
                    (lines_partition_temperature(itemp+1,ispec)-                   &
                     lines_partition_temperature(itemp,ispec))
           else
              itemp = 1
              eps   = 0.d0
              lines_pfunc_outofrange = .true.
           endif
           if((eps.lt.0.).or.(eps.gt.1.)) stop 7729
           !
           ! Construct the partition function by linear interpolation in the table
           !
           lines_psum_local(ispec) =                              &
                (1-eps)*lines_partition_function(itemp,ispec) +   &
                    eps*lines_partition_function(itemp+1,ispec)
        endif
     enddo
  endif
  !
  ! Now loop over frequency
  !
  do inu=inu0,inu1
     !
     ! Loop over actively contributing lines
     !
     do ispec=1,lines_nr_species
        !
        ! Compute Gaussian line width
        !
        a = sqrt(lines_ray_turb(1)**2 + 2*kk*lines_ray_temp(1)/lines_umass(ispec))
        !
        ! Loop over all active lines
        !
        do ilinesub=1,active_nrlines(ispec)
           !
           ! Get line index
           !
           iline = active_lines(ilinesub,ispec)
           !
           ! To get the populations and statistical weights, we must 
           ! distinguish between molecules that have level data and
           ! those that have linelists.
           !
           if(lines_species_fullmolec(ispec)) then
              !
              ! Level data available, and populations have been computed.
              !
              ! Find the upper and lower level population
              !
              nup   = lines_ray_levpop(lines_levelup(iline,ispec),ispec,1)
              ndown = lines_ray_levpop(lines_leveldown(iline,ispec),ispec,1)
              !
              ! Ratio of degenerations
              !
              gratio = lines_level_gdeg(lines_levelup(iline,ispec),ispec) /   &
                       lines_level_gdeg(lines_leveldown(iline,ispec),ispec)
           else
              !
              ! Line list molecule. We must compute the levels on-the-fly
              !
              call lines_linelist_compute_ltepop(ispec,iline,         &
                   lines_ray_nrdens(ispec,1),lines_ray_temp(1),       &
                   lines_psum_local(ispec),ndown,nup)
              !
              ! The ratio of statistical weights:
              !
              gratio = lines_linelist_gup(iline,ispec) /              &
                       lines_linelist_gdown(iline,ispec)
              !
           endif
           !
           ! Compute distance from line center
           !
           ! Note: doppler positive means moving toward observer, hence
           !       the minus sign
           !
           dnu   = freq(inu) - lines_nu0(iline,ispec) - &
                   lines_ray_doppler(1)*lines_nu0(iline,ispec)
           !
           ! Now we have to choose the line profile
           !
           select case(lines_profile)
           case(0)
              !
              ! Gaussian profile
              !
              phi   = exp(-(cc*dnu/(lines_nu0(iline,ispec)*a))**2) * &
                            cc / (a*lines_nu0(iline,ispec)*sqrtpi)
           case(1)
              !
              ! Voigt profile
              ! Note: the lines_ray_lorentz_delta must be set.
              ! (Added by Thomas Peters 2011)
              !
              alpha = cc / a / lines_nu0(iline,ispec)
              x     = alpha * dnu
              y     = alpha * lines_ray_lorentz_delta(1)
              phi   = ( alpha / sqrtpi ) * voigt_humlicek(x,y)
           end select
           !
           ! The j_nu is now:
           !
           !                 h nu_0
           !          j_nu = ------ n_up A_ud phi(nu)
           !                  4 pi
           !
           src(inu) = src(inu) +                                   &
                const1 * lines_nu0(iline,ispec) * phi *            &
                lines_aud(iline,ispec) * nup 
           !
           ! Compute the Bud and Bdu
           !
           bud = lines_aud(iline,ispec) * const2 / lines_nu0(iline,ispec)**3
           bdu = bud * gratio
           !
           ! The alpha_nu is now:
           !
           !                 h nu_0
           !      alpha_nu = ------ ( n_down B_du - n_up B_ud ) phi(nu)
           !                  4 pi
           !
           dum      = const1 * lines_nu0(iline,ispec) * phi *      &
                      ( bdu * ndown - bud * nup )
           !
           ! As of version 0.36 we prevent masering here.
           !
           if(dum.lt.0.d0) then
              dum = 1d-60
              lines_maser_warning = .true.
           endif
           !
           ! Add to alp
           !
           alp(inu) = alp(inu) + dum
           !
        enddo
     enddo
     !
     ! Check that the opacity is non-negative. We do not allow masers!
     !
     if(alp(inu).lt.0.d0) then
        write(stdo,*) 'PROBLEM: Encountered negative opacity in lines!'
        write(stdo,*) '         Aborting...'
        stop 
     endif
  enddo
  !
end subroutine lines_serial_addto_jnu_alpnu


!
!    *** As long as we don't use the below subroutine, I comment it out ***
!
! !----------------------------------------------------------------------------
! !     COMPUTE THE SOURCE AND OPACITY FROM THE LEVEL POPULATIONS
! !            FOR RAY-TRACING METHOD B (1-D PREPARED RAY)
! !
! ! This routine gives the full j_nu and alpha_nu for the lines. Here it is
! ! done along a 1-D prepared ray. This is useful for future parallelization
! ! on a GPU.
! ! NOTE: The src and alp of the lines will be ADDED to anything that
! !       already is in the src and alp (for instance dust emission!)
! !----------------------------------------------------------------------------
! subroutine lines_ray1d_addto_jnu_alpnu(nrfreq,inu0,inu1,ns,freq,src,alp)
!   implicit none
!   integer :: is,inu,iline,ilinesub,ierr,ispec,nrfreq,ns,inu0,inu1,itemp
!   doubleprecision :: src(1:nrfreq,1:ns),alp(1:nrfreq,1:ns),freq(1:nrfreq)
!   doubleprecision :: dummy,nup,ndown,maxshift,tmax,eps
!   doubleprecision :: bdu,bud,gratio,a,dnu,phi
!   doubleprecision :: dum,dum1,dum2
!   doubleprecision :: const1,const2
!   parameter(const1=hh/(4*pi))
!   parameter(const2=cc*cc/(2*hh))
!   !
!   ! Checks
!   !
!   if(debug_check_all.ne.0) then
!      if((.not.allocated(lines_ray_levpop)).or. &
!           (.not.allocated(lines_ray_temp)).or. &
!           (.not.allocated(lines_ray_turb)).or. &
!           (.not.allocated(lines_ray_doppler))) then
!         write(stdo,*) 'ERROR in line module: ray arrays not allocated.'
!         stop
!      endif
!      if(ns.ne.ray_ns) stop 5455
!   endif
!   !
!   ! Determine maximum temperature along the ray
!   !
!   tmax = 0.d0
!   do is=1,ns
!      tmax = max(tmax,lines_ray_temp(is))
!   enddo
!   !
!   ! Check and compute the maximum doppler shift + max thermal shift
!   ! The latter is multiplied by widthmargin to ensure that also the
!   ! thermal line wings are inside.
!   !
!   maxshift = 0.d0
!   do is=1,ns
!      maxshift = max(maxshift,abs(lines_ray_doppler(is)),lines_widthmargin*lines_ray_turb(is)/cc)
!   enddo
!   maxshift = maxshift + lines_widthmargin*sqrt(2*kk*tmax/lines_umass_min)/cc
!   !
!   ! Check if anywhere a line has possibly doppler-leapt farther than the
!   ! intrinsic line width.
!   !
!   do is=2,ns-1
!      dum1 = abs( lines_ray_doppler(is) - lines_ray_doppler(is-1) )
!      dum2 = abs( lines_ray_doppler(is+1) - lines_ray_doppler(is) )
!      dum  = cc * max(dum1,dum2) / &
!!!             ( sqrt(2*kk*lines_ray_temp(is)/lines_umass_max) + lines_ray_turb(is) ) 
!               sqrt(2*kk*lines_ray_temp(is)/lines_umass_max + lines_ray_turb(is)**2) 
!!!!      if(dum.gt.lines_maxdoppler) lines_warn_lineleap=.true.
!   enddo
!   !
!   ! Check which lines are actively contributing to the current wavelength
!   !
!   call lines_find_active_lines_levels(nrfreq,freq,inu0,inu1,maxshift)
!   !
!   ! If there are molecules for which linelists are specified instead
!   ! of full molecular data, then we must compute the partition sum
!   ! for each of these molecules here.
!   !
!   if(allocated(lines_linelist_eup)) then
!      do ispec=1,lines_nr_species
!         if(.not.lines_species_fullmolec(ispec)) then
!            !
!            ! Find the location of the current temperature in the partition function table
!            !
!            call hunt(lines_partition_temperature(:,ispec),               &
!                      lines_partition_ntemp(ispec),lines_ray_temp(1),itemp)
!            if(itemp.ge.lines_partition_ntemp(ispec)) then
!               write(stdo,*) 'ERROR: Temperature out of range of partition function '
!               write(stdo,*) '       for ispec=',ispec,' Temperature = ',lines_ray_temp(1)
!               !write(stdo,*) '       Smallest temperature in table = ',   &
!               !     lines_partition_temperature(1,ispec)
!               write(stdo,*) '       Largest temperature in table = ',   &
!                    lines_partition_temperature(lines_partition_ntemp(ispec),ispec)
!               stop
!            endif
!            if(itemp.ge.1) then
!               eps = (lines_ray_temp(1)-lines_partition_temperature(itemp,ispec)) / &
!                     (lines_partition_temperature(itemp+1,ispec)-                   &
!                      lines_partition_temperature(itemp,ispec))
!            else
!               itemp = 1
!               eps   = 0.d0
!               lines_pfunc_outofrange = .true.
!            endif
!            if((eps.lt.0.).or.(eps.gt.1.)) stop 7729
!            !
!            ! Construct the partition function by linear interpolation in the table
!            !
!            lines_psum_local(ispec) =                              &
!                 (1-eps)*lines_partition_function(itemp,ispec) +   &
!                     eps*lines_partition_function(itemp+1,ispec)
!         endif
!      enddo
!   endif
!   !
!   ! Now loop over frequency
!   !
!   ! ############## FUTURE: FOR NVIDIA GRAPHICS CARD: PARALLELIZE IN NU HERE #########
!   !
!   do inu=inu0,inu1
!      !
!      ! Now go along the ray to define the src and alp
!      !
!      do is=1,ray_ns
!         !
!         ! Loop over actively contributing lines
!         !
!         do ispec=1,lines_nr_species
!            do ilinesub=1,active_nrlines(ispec)
!               !
!               ! Get line index
!               !
!               iline = active_lines(ilinesub,ispec)
!               !
!               ! To get the populations and statistical weights, we must 
!               ! distinguish between molecules that have level data and
!               ! those that have linelists.
!               !
!               if(lines_species_fullmolec(ispec)) then
!                  !
!                  ! Level data available, and populations have been computed.
!                  !
!                  ! Find the upper and lower level population
!                  !
!                  nup   = lines_ray_levpop(lines_levelup(iline,ispec),ispec,is)
!                  ndown = lines_ray_levpop(lines_leveldown(iline,ispec),ispec,is)
!                  !
!                  ! Ratio of degenerations
!                  !
!                  gratio = lines_level_gdeg(lines_levelup(iline,ispec),ispec) /   &
!                           lines_level_gdeg(lines_leveldown(iline,ispec),ispec)
!               else
!                  !
!                  ! Line list molecule. We must compute the levels on-the-fly
!                  !
!                  call lines_linelist_compute_ltepop(ispec,iline,         &
!                       lines_ray_nrdens(ispec,1),lines_ray_temp(1),       &
!                       lines_psum_local(ispec),ndown,nup)
!                  !
!                  ! The ratio of statistical weights:
!                  !
!                  gratio = lines_linelist_gup(iline,ispec) /              &
!                           lines_linelist_gdown(iline,ispec)
!                  !
!               endif
!               !
!               ! Determine the value of the line profile function
!               !
!               ! Note: doppler positive means moving toward observer, hence
!               !       the minus sign
!               !
!               dnu   = freq(inu) - lines_nu0(iline,ispec) - &
!                       lines_ray_doppler(is)*lines_nu0(iline,ispec)
!!               a     = lines_ray_turb(is) + sqrt(2*kk*lines_ray_temp(is)/lines_umass(ispec))
!               a     = sqrt(lines_ray_turb(is)**2 + 2*kk*lines_ray_temp(is)/lines_umass(ispec))
!               phi   = exp(-(cc*dnu/(lines_nu0(iline,ispec)*a))**2) * &
!                       cc / (a*lines_nu0(iline,ispec)*sqrtpi)
!               !
!               ! The j_nu is now:
!               !
!               !                 h nu_0
!               !          j_nu = ------ n_up A_ud phi(nu)
!               !                  4 pi
!               !
!               src(inu,is) = src(inu,is) +                             &
!                    const1 * lines_nu0(iline,ispec) * phi *            &
!                    lines_aud(iline,ispec) * nup 
!               !
!               ! Compute the Bud and Bdu
!               !
!               bud = lines_aud(iline,ispec) * const2 / lines_nu0(iline,ispec)**3
!               bdu = bud * gratio
!               !
!               ! The alpha_nu is now:
!               !
!               !                 h nu_0
!               !      alpha_nu = ------ ( n_down B_du - n_up B_ud ) phi(nu)
!               !                  4 pi
!               !
!               alp(inu,is) = alp(inu,is) +                             &
!                    const1 * lines_nu0(iline,ispec) * phi *            &
!                    ( bdu * ndown - bud * nup )
!               !
!            enddo
!         enddo
!      enddo
!   enddo
!   !
! end subroutine lines_ray1d_addto_jnu_alpnu



!----------------------------------------------------------------------------
!     COMPUTE THE SOURCE AND OPACITY FROM THE LEVEL POPULATIONS
!               FOR THE DOPPLER-CATCHING METHOD
!
! This is the line stuff without the line profile. This is needed for the
! line "doppler catching" method. Basically everything of the line is
! calculated, except the final bit: the multiplication with the line 
! profile. The reason is that we will only know the line profile once we
! know exactly the gas velocity, and in the doppler catching method this
! is precisely the thing that is going to be varied in small sub-steps
! along a ray-element to avoid doppler jumps.
!----------------------------------------------------------------------------
subroutine lines_get_active_nup_ndown(nlines,nup,ndown)
   implicit none
   integer :: iline,ilinesub,ispec,nlines,icnt,itemp
   doubleprecision :: nup(1:nlines),ndown(1:nlines)
   doubleprecision :: nu,nd
   doubleprecision :: gratio,eps
   !
   ! Checks
   !
   if(debug_check_all.ne.0) then
      if((.not.allocated(lines_ray_levpop)).or. &
           (.not.allocated(lines_ray_temp)).or. &
           (.not.allocated(lines_ray_turb)).or. &
           (.not.allocated(lines_ray_doppler))) then
         write(stdo,*) 'ERROR in line module: ray arrays not allocated.'
         stop
      endif
   endif
   !
   ! If there are molecules for which linelists are specified instead
   ! of full molecular data, then we must compute the partition sum
   ! for each of these molecules here.
   !
   if(allocated(lines_linelist_eup)) then
      do ispec=1,lines_nr_species
         if(.not.lines_species_fullmolec(ispec)) then
            !
            ! Find the location of the current temperature in the partition function table
            !
            call hunt(lines_partition_temperature(:,ispec),               &
                      lines_partition_ntemp(ispec),lines_ray_temp(1),itemp)
            if(itemp.ge.lines_partition_ntemp(ispec)) then
               write(stdo,*) 'ERROR: Temperature out of range of partition function '
               write(stdo,*) '       for ispec=',ispec,' Temperature = ',lines_ray_temp(1)
               !write(stdo,*) '       Smallest temperature in table = ',   &
               !     lines_partition_temperature(1,ispec)
               write(stdo,*) '       Largest temperature in table = ',   &
                    lines_partition_temperature(lines_partition_ntemp(ispec),ispec)
               stop
            endif
            if(itemp.ge.1) then
               eps = (lines_ray_temp(1)-lines_partition_temperature(itemp,ispec)) / &
                     (lines_partition_temperature(itemp+1,ispec)-                   &
                      lines_partition_temperature(itemp,ispec))
            else
               itemp = 1
               eps   = 0.d0
               lines_pfunc_outofrange = .true.
            endif
            if((eps.lt.0.).or.(eps.gt.1.)) stop 7729
            !
            ! Construct the partition function by linear interpolation in the table
            !
            lines_psum_local(ispec) =                              &
                 (1-eps)*lines_partition_function(itemp,ispec) +   &
                     eps*lines_partition_function(itemp+1,ispec)
         endif
      enddo
   endif
   !
   ! Loop over actively contributing lines
   !
   icnt = 0
   do ispec=1,lines_nr_species
      do ilinesub=1,active_nrlines(ispec)
         !
         ! Increase counter
         !
         icnt  = icnt + 1
         if(icnt.gt.nlines) then
            write(stdo,*) 'INTERNAL ERROR: In preparing doppler catch: Nr of lines inconsistent.'
            stop
         endif
         !
         ! Get line index
         !
         iline = active_lines(ilinesub,ispec)
         !
         ! To get the populations and statistical weights, we must 
         ! distinguish between molecules that have level data and
         ! those that have linelists.
         !
         if(lines_species_fullmolec(ispec)) then
            !
            ! Level data available, and populations have been computed.
            !
            ! Find the upper and lower level population
            !
            nup(icnt)   = lines_ray_levpop(lines_levelup(iline,ispec),ispec,1)
            ndown(icnt) = lines_ray_levpop(lines_leveldown(iline,ispec),ispec,1)
         else
            !
            ! Line list molecule. We must compute the levels on-the-fly
            !
            call lines_linelist_compute_ltepop(ispec,iline,         &
                 lines_ray_nrdens(ispec,1),lines_ray_temp(1),       &
                 lines_psum_local(ispec),nd,nu)
            ndown(icnt) = nd
            nup(icnt)   = nu
         endif
      enddo
   enddo
   !
 end subroutine lines_get_active_nup_ndown


!----------------------------------------------------------------------------
!                  MAKE ANU AND JNU FROM NUP AND NDOWN
!
! This subroutine is used for the doppler catching method. The idea of the
! doppler catching method of ray-tracing is that during the ray-tracing it
! is continuously checked if a ray element doppler-jumps over a line. This
! can happen if the Delta (v.n), where n is the vector along the ray and v
! is the gas velocity, is larger than the local line width (consisting of
! thermal and microturbulent line width). To solve this, the ray element is
! then chopped into sub-elements, and the integration is done over all these
! sub-elements. The subroutine lines_jnu_anu_from_nup_ndown() computes the source
! terms for all these sub-elements, starting from the nup and ndown
! for each of the lines that may contribute.
! ----------------------------------------------------------------------------
subroutine lines_jnu_anu_from_nup_ndown(nlines,nup,ndown,doppler,&
                                  turb,temp,nu,anu,jnu)
  implicit none
  integer :: icnt,nlines,iline,ispec,ilinesub
  doubleprecision :: nup(1:nlines),ndown(1:nlines)
  double precision :: nu,anu,jnu,doppler,turb,temp,dnu,phi,a,alpha,x,y
  double precision :: gratio,jbase,alpbase,bdu,bud
  doubleprecision :: const1,const2
  parameter(const1=hh/(4*pi))
  parameter(const2=cc*cc/(2*hh))
  icnt = 0
  anu  = 0.d0
  jnu  = 0.d0
  do ispec=1,lines_nr_species
     !
     ! Compute the line width from (unresolved=subgrid) turbulence and 
     ! the gas temperature, assuming that these are the only intrinsic 
     ! broadning processes
     !
     a     = sqrt(turb**2 + 2*kk*temp/lines_umass(ispec))
     !
     ! Only include this molecule if a.gt.0. A zero a means that the molecule
     ! is not included here.
     !
     if(a.gt.0.d0) then
        !
        ! Now for each line of this molecule
        !
        do ilinesub=1,active_nrlines(ispec)
           !
           ! Increase counter
           !
           icnt  = icnt + 1
           if(icnt.gt.nlines) then
              write(stdo,*) 'INTERNAL ERROR: In preparing doppler catch: Nr of lines inconsistent. Warn author.'
              stop
           endif
           !
           ! Get line index
           !
           iline = active_lines(ilinesub,ispec)
           !
           ! Compute distance from line center
           !
           ! Note: doppler positive means moving toward observer, hence
           !       the minus sign
           !
           dnu   = nu - lines_nu0(iline,ispec) - &
                        doppler*lines_nu0(iline,ispec)
           !
           ! Now we have to choose the line profile
           !
           select case(lines_profile)
           case(0)
              !
              ! Gaussian profile
              !
              phi   = exp(-(cc*dnu/(lines_nu0(iline,ispec)*a))**2) * &
                            cc / (a*lines_nu0(iline,ispec)*sqrtpi)
           case(1)
              !
              ! Voigt profile
              ! Note: the lines_ray_lorentz_delta must be set.
              ! (Added by Thomas Peters 2011)
              !
              alpha = cc / a / lines_nu0(iline,ispec)
              x     = alpha * dnu
              y     = alpha * lines_ray_lorentz_delta(1)
              phi   = ( alpha / sqrtpi ) * voigt_humlicek(x,y)
           end select
           !
           ! The ratio of statistical weights:
           !
           if(lines_species_fullmolec(ispec)) then
              !
              ! Full molecule
              !
              gratio = lines_level_gdeg(lines_levelup(iline,ispec),ispec) /   &
                       lines_level_gdeg(lines_leveldown(iline,ispec),ispec)
           else
              !
              ! Line list molecule
              !
              gratio = lines_linelist_gup(iline,ispec) /              &
                       lines_linelist_gdown(iline,ispec)
              !
           endif
           !
           ! Compute the Bud and Bdu
           !
           bud = lines_aud(iline,ispec) * const2 / lines_nu0(iline,ispec)**3
           bdu = bud * gratio
           !
           ! The j_base is now 
           !
           !               h nu_0
           !      j_base = ------ n_up A_ud
           !                4 pi
           !
           jbase = const1 * lines_nu0(iline,ispec) * lines_aud(iline,ispec) * nup(icnt)
           !
           ! The alpha_base is now:
           !
           !                   h nu_0
           !      alpha_base = ------ ( n_down B_du - n_up B_ud ) 
           !                    4 pi
           !
           alpbase = const1 * lines_nu0(iline,ispec) * ( bdu * ndown(icnt) - bud * nup(icnt) )
           !
           ! As of version 0.36 we prevent masering here
           !
           if(alpbase.lt.0.d0) then
              alpbase = 1d-60
              lines_maser_warning = .true.
           endif
           !
           ! Add all line contributions to the anu and jnu
           !
           anu  = anu + alpbase * phi
           jnu  = jnu + jbase * phi
           !
        enddo
     else
        icnt = icnt + active_nrlines(ispec)
     endif
  enddo
  !
  ! Do a stupidity check
  !
  if(icnt.ne.nlines) then
     write(stdo,*) 'INTERNAL ERROR while computing a_nu and j_nu for lines nup, ndown:'
     write(stdo,*) '      nr of active lines mismatch. Warn author.'
     write(stdo,*) '      icnt = ',icnt,' nlines = ',nlines
     stop
  endif
end subroutine lines_jnu_anu_from_nup_ndown


!----------------------------------------------------------------------------
!              FOR CONVENIENCE OF THE USER: WRITE LINE INFO
!----------------------------------------------------------------------------
subroutine lines_print_lineinfo(iline,ispec)
  implicit none
  integer :: iline,ispec
  call read_linedata(1)
  if(ispec.lt.1.or.ispec.gt.lines_nr_species) then
     write(stdo,*) 'ERROR: ispec out of range for lines'
     write(stdo,*) ispec
     stop
  endif
  if(iline.lt.1.or.iline.gt.lines_nrlines(ispec)) then
     write(stdo,*) 'ERROR: iline out of range for lines'
     write(stdo,*) iline
     stop
  endif
  open(unit=1,file='lineinfo.out')
  write(1,*) 1,    '            ; iformat'
  write(1,*) lines_nr_species,'            ; total nr of species'
  write(1,*) ispec,'            ; ispec'
  write(1,*) trim(lines_speciesname(ispec)),'                      ; species name'
  write(1,*) lines_umass(ispec)/mp,'  ; species mass in units of m_p'
  write(1,*) lines_maxnrlevels,'            ; max nr of levels'
  write(1,*) lines_nrlevels(ispec),'            ; total nr of levels for this species'
  write(1,*) lines_maxnrlines,'            ; max nr of lines'
  write(1,*) lines_nrlines(ispec),'            ; total nr of lines for this species'
  write(1,*) iline,'            ; iline'
  write(1,*) trim(lines_level_name(lines_levelup(iline,ispec),ispec)),'                       ; upper level name'
  write(1,*) trim(lines_level_name(lines_leveldown(iline,ispec),ispec)),'                       ; lower level name'
  write(1,*) 1d4*cc/lines_nu0(iline,ispec),'  ; line center wavelength in micron'
  close(1)
end subroutine lines_print_lineinfo

!--------------------------------------------------------------------------
!              VOIGT FUNCTION APPROXIMATION BY HUMLICEK
!
! This function is a modified version of a subroutine taken from the 
! following web site: 
!
!    http://www.op.dlr.de/oe/ir/voigt.html
!
! and written by F. Schreier (DLR - Institute for Optoelectronics
! Oberpfaffenhofen, 82230 Wessling) from the publication
!
!  F. Schreier, J. Quant. Spectros. Radiat. Transfer 48, 743-762 (1992)
!
! which bases the algorithm on the paper by Humlicek:
!
!  J. Humlicek, J. Quant. Spectros. Radiat. Transfer 27, 437 (1982)
!
! It calculates the real part of the complex probability function for
! complex argument Z=X+iY. The real part = voigt function K(x,y).
!
! Parameters:                                                      
!   X      X argument                                           in  
!   Y      Voigt function parameter, ratio of lorentz/doppler   in  
!
! Result:
!   function value = Voigt function
!  
! The stated accuracy is claimed to be 1.0E-04 by the author.  R. H. Norton
! has checked the accuracy by comparing values computed using a program
! written by B. H. Armstrong, and the accuracy claim seems to be warranted.
!
! NOTE 1: The normalization of this function is such that 
!
!  /
!  | voigt(x) dx = sqrt(pi)
!  /
!
! NOTE 2: The x coordinate is defined such that the Gaussian component
!         is always exp(-x^2). The parameter y simply changes the 
!         width of the Lorenz component compared to the Gaussian.
!
! fgs 12/91
! Modified: CPD 2011
!--------------------------------------------------------------------------
function voigt_humlicek(x,y)
  implicit none
  double precision :: x, y, s, voigt_humlicek, ax
  COMPLEX*16 :: t, u, z
  COMPLEX*16 :: APPROX1, APPROX2, APPROX3, APPROX4
  !
  ! Four different approximations, as inline formulae
  !
  APPROX1(T)   = (T * .5641896) / (.5 + (T * T))
  APPROX2(T,U) = (T * (1.410474 + U*.5641896))/ (.75 + (U *(3.+U)))
  APPROX3(T)   = ( 16.4955 + T * (20.20933 + T * (11.96482 + &
                 T * (3.778987 + 0.5642236*T)))) &
                 / ( 16.4955 + T * (38.82363 + T * &
                 (39.27121 + T * (21.69274 + T * (6.699398 + T))))) 
  APPROX4(T,U) = (T * (36183.31 - U * (3321.99 - U * (1540.787 - U &
                 *(219.031 - U *(35.7668 - U *(1.320522 - U * .56419)))))) &
                 / (32066.6 - U * (24322.8 - U * (9022.23 - U * (2186.18 &
                 - U * (364.219 - U * (61.5704 - U * (1.84144 - U))))))))
  !
  ! If y=0, then we have a Gaussian profile, so we are done quickly
  !
  if(y.eq.0.d0) then
     voigt_humlicek = exp(-x**2)
     return
  endif
  !
  ! Voigt Profile...
  !
  ! Distinguish different y-regimes
  !
  if(y.gt.15.d0) then
     !
     ! All points are in region I
     !
     t = cmplx(y,-x)
     z = APPROX1(t)
     !
  elseif(y.ge.5.5d0) then
     !
     ! Points are in region I or region II
     !
     t  = cmplx(y,-x)
     s  = abs(x) + y
     if(s.ge.15.d0) then
        z = APPROX1(t)
     else
        u = t * t
        z = APPROX2(t,u)
     endif
     !
  elseif(y.gt.0.75d0) then
     t  = cmplx(y,-x)
     s  = abs(x) + y
     if(s.ge.15.d0) then
        z = APPROX1(t)
     elseif(s.lt.5.5d0) then
        z = APPROX3(t)
     else
        u = t * t
        z = APPROX2(t,u)
     endif
  else
     t  = cmplx(y,-x)
     ax = abs(x)
     s  = ax + y
     if(s.ge.15.d0) then
        !
        ! region I
        !
        z = APPROX1(t)
     elseif(s.ge.5.5d0) then
        !
        ! region II
        !
        u = t * t
        z = APPROX2(t,u)
     elseif(y.ge.(0.195d0*ax-0.176d0)) then
        !
        ! region III
        !
        z = APPROX3(t)
     else
        !
        ! region IV
        !
        u = t * t
        z = cdexp(u) - APPROX4(t,u)
     endif
  endif
  !
  voigt_humlicek = real(z)
  return
end function voigt_humlicek

end module lines_module
