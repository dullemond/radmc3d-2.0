!-------------------------------------------------------------------------
!               MODULE FOR MONTE CARLO RADIATIVE TRANSFER
!
! This module does radiative transfer in dust continuum. The module
! can do the following things:
!
!  1) Thermal Monte Carlo according to the Bjorkman & Wood method 
!     for computing the dust temperatures
!
!  2) Monochromatic Monte Carlo for computing the scattering source 
!     function or for computing the mean intensity.
!
! The dust scattering can be done to various degrees of realism. This is
! regulated by the scattering_mode integer (see dust_module.f90). The
! various possible values of scattering_mode are:
!
!     =0 Dust scattering not included
!
!     =1 Isotropic scattering approximation
!
!     =2 Include anisotropic scattering with Henyey-Greenstein function
!
!     =3 Include anisotropic scattering with full phase function, based
!        on the full scattering matrix information
!
!     =4 Include scattering with full phase function and polarization 
!        (for randomly oriented particles), but only in the scattering 
!        source function (i.e. the Monte Carlo photon packages remain
!        unpolarized, and the polarization of the scattering source
!        function only takes into account the last scattering before
!        observation)
!
!     =5 Include scattering with full phase function and polarization 
!        (for randomly oriented particles), full treatment.
!
!-------------------------------------------------------------------------
module montecarlo_module
!$ use omp_lib
use mathroutines_module
use rtglobal_module
use dust_module
use polarization_module
use stars_module
use quantum_module
use amr_module
use amrray_module
use ioput_module
use constants_module
!!!INCLUDE "omp_lib.h"
!
!-----------------------------------------------------------------------
!
! Flag for OMP parallel version
!
logical :: mc_openmp_parallel = .false.
!
! Current photon number
! (Must ensure huge range here)
!
integer(kind=8) :: ieventcounttot
double precision :: mc_visitcell,mc_revisitcell,mc_revisitcell_max
!
! Flag saying whether or not to interpolate the temperature emission
! database in temperature
!
integer :: intplt
!
!    For multi-frequency images with scattering, this flag tells whether
!    the Monte Carlo scattering source array is also frequency-dispersed, 
!    or whether it is monochromatic (i.e. one single freq-bin) to save
!    memory.
!
logical :: camera_mcscat_monochromatic=.false.
!
! Some defines 
!
double precision, parameter :: minieps = 1d-12
double precision, parameter :: epsmargin=1d-8
!
! The energy of each photon package
!
double precision :: energy
!
! The photon package structure for the polarization module
! This is not needed for non-polarized radiation
!
type(photon) :: photpkg
!
! When using weighted photon packages, the energy of these packages
!
logical :: mc_weighted_photons=.true.
double precision :: mc_energy_extlum,mc_energy_stellarsrc,mc_energy_bc
double precision :: mc_energy_heatsource,mc_energy_quant,mc_energy_thermal
double precision, allocatable :: mc_energy_stars(:)
double precision, allocatable :: mc_energy_illum(:)
!
! The overall averaged output spectrum (averaged over all directions)
!
double precision,allocatable :: mc_average_spectrum(:)
double precision,allocatable :: mc_integerspec(:)
!
! The partial energy for each dust species
!
double precision,allocatable :: mc_enerpart(:)
!
! The cumulative energy in each cell for each dust species
!
double precision,allocatable :: mc_cumulener(:,:),mc_cumulener_bk(:,:)
!
! The opacities used for this MC module
!
double precision,allocatable :: kappa_a(:,:),kappa_s(:,:),kappa_g(:,:)
double precision,allocatable :: alpha_a(:),alpha_s(:),alphacum(:)
double precision :: alpha_a_tot,alpha_s_tot
double precision,allocatable :: zmatrix(:,:,:,:),zcumul(:,:,:,:)
!
! For the grain alignment we need additional arrays
!
integer :: mc_align_munr
double precision, allocatable :: mc_align_mu(:)
double precision, allocatable :: mc_align_orth(:,:,:),mc_align_para(:,:,:)
double precision, allocatable :: mc_align_opcumul(:,:,:)
!
! For photon statistics 
! (Must ensure huge range here)
!
double precision,allocatable :: mc_ilastphot(:),mc_iphotcount(:)
double precision :: mc_iphotcurr
!!!!!integer(kind=8),allocatable :: mc_ilastphot(:),mc_iphotcount(:)
!!!!!integer(kind=8) :: mc_iphotcurr
! 
! For arrays of integrals as a function of temperature
!
double precision, allocatable :: db_temp(:),db_enertemp(:,:),db_cumul(:)
double precision, allocatable :: db_logenertemp(:,:)
double precision, allocatable :: db_emiss(:,:,:),db_cumulnorm(:,:,:)
integer :: db_ntemp
!
! For internal use in the routines that compute the above arrays
!
double precision, allocatable :: cellalpha(:),fnu_diff(:),diffemis(:)
double precision, allocatable :: enercum(:)
!
! For convenience for choosing which source of photons to choose
!
double precision :: mc_cumlum1,mc_cumlum2,mc_cumlum3,mc_cumlum4,mc_cumlum5
!
! Arrays for the MRW routines
!
double precision, allocatable :: mrw_cumulener_bk(:)
integer :: mrw_db_ntemp
double precision, allocatable :: mrw_db_temp(:)
double precision, allocatable :: mrw_db_kappa_abs_planck(:,:)
double precision, allocatable :: mrw_alpha_tot_ross(:)
logical, allocatable :: mrw_cell_uses_mrw(:)
double precision :: alpha_t_rm,alpha_a_pm_tot
double precision, allocatable :: alpha_a_pm(:),mrw_dcumen(:)
!
! --------------------------------------------------------------
!
! The maximum absorption tau a photon package travels in a scattering
! Monte Carlo simulation before the package is dropped.
!
double precision :: mc_scat_maxtauabs
double precision :: mc_scat_energy_rellimit
!
! The maximum number of discrete scattering events treated.
! The default is 0 = infinite number of scatterings allowed.
! If set to 1, then we are in single-scattering mode, meaning
! that photon packages are followed from their source up to
! (but not beyond) the first discrete scattering event. This
! means that the scattering source function (used by the 
! camera module to produce images) only includes single-
! scattering events. This can be useful to check how big
! the effect of multiple scattering is on the image: Just
! make a normal image (mc_max_nr_scat_events=-1) and then
! make an image with single scattering (mc_max_nr_scat_events=1)
! and check the difference between the two.
!
integer :: mc_max_nr_scat_events = -1
!
! For some Monte Carlo calculation types we need a special local
! set of frequencies. This is, by the way, not true for the Bjorkman
! & Wood (thermal Monte Carlo) transfer. That will be always done on the 
! global frequency grid.
!
double precision, allocatable :: mc_frequencies(:)
integer :: mc_nrfreq=0
!
! The number of directions for which the scattering source function must
! be stored.
!
integer :: mcscat_nrdirs=0, mcscat_current_dir=1
double precision, allocatable :: mcscat_dirs(:,:)
double precision, allocatable :: mcscat_svec(:,:)  ! For polarization
!
! The storage array for the scattering source function at the frequencies
! stored in mc_frequencies. This includes the possibility for isotropic
! scattering (with mcscat_nrdirs=1)
!
double precision, allocatable :: mcscat_scatsrc_iquv(:,:,:,:)
!
! The storage array for the mean intensity function at the frequencies
! stored in mc_frequencies. 
!
double precision, allocatable :: mcscat_meanint(:,:)
!
! The array of the thermal emissivity
!
double precision, allocatable :: mc_cumulthermemis(:)
double precision :: mc_thermemistot
!
! Thermal emission from the boundaries, if these thermal
! boundaries are activated (only possible for cartesian
! coordinates) via the flags thermal_bc_active().
!
double precision :: mc_bc_lum(2,3) = 0.d0
double precision :: mc_bc_lumcum(7) = 0.d0
double precision :: mc_bc_lumall = 0.d0
doubleprecision, allocatable :: mc_bc_cumspec(:,:,:)
!
! For anisotropic scattering mode: the phase function, normalized
! such that it is 1 for isotropic scattering.
!
double precision, allocatable :: mcscat_phasefunc(:,:)
!
! For the local observer mode
!
logical :: mcscat_localobserver=.false.
double precision :: mcscat_localobs_pos(1:3)
!
!    Arrays for templates of stellar source spectra
!
double precision, allocatable :: mc_stellarsrc_templates(:,:)
!
! Global flag for Monte Carlo
!
logical :: mc_photon_destroyed
!
! For analysis or debugging: Allowing to record only e.g. second
! to fourth scattering source or mean intensity (only for monochromatic
! monte carlo). Note that this does not yield any real observables, but
! it can help understand the results.
!
integer :: selectscat_iscat
integer :: selectscat_iscat_first = 1
integer :: selectscat_iscat_last  = 1000000000
!
!----TO-ADD----
!
! Array for skipping very optically thick cells
!
! double precision, allocatable :: mc_scat_src_skip(:)
!----TO-ADD----
!
! OpenMP Parallellization:
! Global flag for each cell which denotes if the cell is reserved for one thread and 
! therefore blocked for all other ones
!
!$ integer(kind=OMP_LOCK_KIND),allocatable::lock(:)
!
! OpenMP Parallellization:
! id gives identity number of the thread
! nthreads gives the number of threads used for the parallelization
! nprocs gives the number of processor cores which are available
!
!$ integer :: id,nthreads,nprocs,conflict_counter
!
! OpenMP Parallellization:
! Global variables used in subroutine calls within the parallel region which are threadprivate
!
!$OMP THREADPRIVATE(id)
!$OMP THREADPRIVATE(photpkg,energy)
!$OMP THREADPRIVATE(mc_enerpart)
!$OMP THREADPRIVATE(alpha_a,alpha_s)
!$OMP THREADPRIVATE(alphacum)
!$OMP THREADPRIVATE(alpha_a_tot, alpha_s_tot)
!$OMP THREADPRIVATE(mc_iphotcurr)
!$OMP THREADPRIVATE(enercum)
!$OMP THREADPRIVATE(mrw_alpha_tot_ross)
!$OMP THREADPRIVATE(alpha_t_rm,alpha_a_pm_tot)
!$OMP THREADPRIVATE(alpha_a_pm,mrw_dcumen)
!$OMP THREADPRIVATE(mc_photon_destroyed)
!$OMP THREADPRIVATE(mcscat_phasefunc)
!$OMP THREADPRIVATE(db_cumul)
!$OMP THREADPRIVATE(selectscat_iscat)

contains

!--------------------------------------------------------------------------
!                    Initialize the Monte Carlo arrays
!
! This routine is automatically called by the do_monte_carlo_bjorkmanwood()
! routine. It allocates some more arrays and does a number of consistency
! checks.
!
! ARGUMENTS:
!  params             The parameter structure
!  ierr               Error code. 0 means no error.
!  mcaction           =1   --> Bjorkman & Wood
!                     =101 --> Scattering MC to compute the scattering source
!                     =102 --> Scattering MC to compute the mean intensity
!--------------------------------------------------------------------------
subroutine montecarlo_init(params,ierr,mcaction,resetseed)
  implicit none
  type(mc_params) :: params
  integer :: ierr,mcaction
  logical, optional :: resetseed
  !
  integer :: ispec,inu,iinu,itempl,istar,imu,iz
  double precision :: temp,eps,dum,kapscat,g,orth,para
  double precision :: dumarr6(1:6)
  logical :: iwarned,fex
  !
  ! Set flag if OpenMP-parallel
  !
  !$ mc_openmp_parallel = .true.
  !
  ! If OpenMP Parallel, make a warning
  !
  !$ write(stdo,*) 'Beware: The OpenMP-parallel acceleration of RADMC-3D has '
  !$ write(stdo,*) '        not yet been tested with all of the modes and'
  !$ write(stdo,*) '        features that RADMC-3D offers. Please check '
  !$ write(stdo,*) '        your parallel results against the serial version'
  !$ write(stdo,*) '        (i.e. compiling without -fopenmp). '
  !
  ! Currently the polarization module is not compatible with mirror
  ! symmetry mode in spherical coordinates. 
  !
  if((scattering_mode.ge.4).and.amrray_mirror_equator) then
     write(stdo,*) 'ERROR: Currently polarization is not compatible with'
     write(stdo,*) '       equatorial mirror symmetry mode. This means that'
     write(stdo,*) '       your theta-grid must not only cover the upper'
     write(stdo,*) '       quadrant (0<theta<=pi/2) but both the upper'
     write(stdo,*) '       and the lower quadrant (0<theta<pi).'
     stop
  endif
  !
  ! Small angle scattering mode (including polarization) for 2-D spherical coordinates
  ! requires a special treatment. Check basic things, just to be sure.
  !
  if((scattering_mode.ge.2).and.(igrid_coord.ge.100).and.(amr_dim.ne.3)) then
     if(amr_dim.eq.1) stop 4096
     if(scattering_mode.lt.5) stop 4097
     if(.not.dust_2daniso) stop 4099
     if(mcscat_nrdirs.ne.dust_2daniso_nphi + 1) stop 4098
  endif
  if(dust_2daniso) then
     if(scattering_mode.lt.5) stop 3056
     if(amr_dim.ne.2) stop 3057
     if(igrid_coord.lt.100) stop 3055
     !
     ! If dust_2daniso_nphi is very large, then give a tip for speed up
     !
     if(dust_2daniso_nphi.gt.60) then
        write(stdo,*) 'Tip for speed-up: You using full polarized scattering in 2-D axisymmetry, '
        write(stdo,*) '   meaning that RADMC-3D has to internally store the scattering source function'
        write(stdo,*) '   in 3-D (i.e. also on a phi-grid). The default is to use 360', &
             ' phi-grid points for this (accurate but slow).'
        if(dust_2daniso_nphi.ne.360) then
           write(stdo,'(A28,I4)') '    The current setting is: ',dust_2daniso_nphi
        endif
        write(stdo,*) '   You can speed this up (though at your own risk) by adding ', &
             'the following line to radmc3d.inp'
        write(stdo,*) '   dust_2daniso_nphi = 60'
     else
        if(dust_2daniso_nphi.gt.16) then
           write(stdo,'(A36,I4,A62)') ' Warning: using dust_2daniso_nphi = ',dust_2daniso_nphi, &
                ' (this is fine, but may make the phase function less accurate)'
        else
           write(stdo,'(A36,I4,A62)') ' Warning: using dust_2daniso_nphi = ',dust_2daniso_nphi, &
                ': This can be dangerous!!'
        endif
     endif
  endif
  !
  ! Polarized scattering and multiple vantage points at the same time 
  ! (mcscat_nrdirs.gt.1) is not possible. This should, however, not be
  ! hard to build in: just a loop over idir=1,mcscat_nrdirs and the
  ! allocation of the scattering source array to a bigger idir dimension. 
  !
  if((scattering_mode.ge.4).and.((.not.dust_2daniso).and.(mcscat_nrdirs.gt.1))) then
     write(stdo,*) 'ERROR: Polarized scattering is currently incompatible'
     write(stdo,*) '       with multiple vantage points in parallel.'
     stop
  endif
  !
  !
  ! First clean up any possible remaining crap from earlier runs
  !
  call montecarlo_partial_cleanup()
  !
  ! If requested, then reset the random seed to the original value. This can be useful
  ! if you want to make movies and you don't want to be annoyed by constantly changing
  ! noise for nearly identical frames.
  !
  if(present(resetseed)) then
     if(resetseed) then
        iseed = iseed_start
     endif
  endif
  !
  ! If we must store a scattering source function, we must have at least
  ! mcscat_nrdirs.ge.1. This is true for the scattering Monte Carlo.
  !
  if(mcaction.eq.101) then
     if(mcscat_nrdirs.le.0) then
        mcscat_nrdirs = 1
     endif
  endif
  !
  ! Do some tests
  !
  if(.not.rt_incl_dust) then
     write(stdo,*) 'ERROR in Montecarlo module: rt_incl_dust is set to false, but the Monte Carlo'
     write(stdo,*) '        simulations are meant for dust radiative transfer...'
     write(stdo,*) '      Tip: Are you sure you have the file dustopac.inp present?'
     stop
  endif
  if(nrcells.le.0) then 
     write(stdo,*) 'ERROR in Montecarlo module: zero number of cells'
     stop
  endif
  if(igrid_mirror.eq.1) then
     if(igrid_type.ge.100) then
        write(stdo,*) 'ERROR: Mirror symmetry only available for 2-D/3-D spherical coordinates in AMR- or regular type grid'
        stop
     endif
     if((igrid_coord.lt.100).or.(igrid_coord.ge.200)) then
        write(stdo,*) 'ERROR: Mirror symmetry only available for 2-D/3-D spherical coordinates'
        stop
     endif
     if(amr_ydim.eq.0) then
        write(stdo,*) 'ERROR: Mirror symmetry only available for 2-D/3-D spherical coordinates, not 1-D spherical coordinates.'
        stop
     endif
  endif
  if(nrcells.gt.nrcellsmax) then
     write(stdo,*) 'ERROR: Nr of cells larger than max'
     stop
  endif
  if(.not.allocated(cellindex)) then
     write(stdo,*) 'INTERNAL ERROR IN MONTECARLO: Cellindex array not allocated'
     stop
  endif
  if(.not.allocated(cellvolume)) then
     write(stdo,*) 'INTERNAL ERROR IN MONTECARLO: Cellvolume array not allocated'
     stop
  endif
  if(.not.allocated(freq_nu)) then
     write(stdo,*) 'INTERNAL ERROR IN MONTECARLO: Frequency grid not yet allocated.'
     stop
  endif
  if(freq_nr.lt.1) then
     write(stdo,*) 'INTERNAL ERROR IN MONTECARLO: Nr of frequencies has wrong value.'
     stop
  endif
  if(star_sphere.and.mc_weighted_photons) then
     iwarned = .false.
     do istar=1,nstars
        dum = sqrt((star_pos(1,istar)-grid_contsph_x)**2+&
                   (star_pos(2,istar)-grid_contsph_y)**2+&
                   (star_pos(3,istar)-grid_contsph_z)**2)
        if(dum.ge.grid_contsph_r) then
           if(.not.iwarned) then
              write(stdo,*) 'WARNING: When using non-point-like stars outside the grid'
              write(stdo,*) '         as well as having weighted photons mode switched on,'
              write(stdo,*) '         please be warned that to simplify the launching of'
              write(stdo,*) '         photons we treat (only for the Monte Carlo run, and'
              write(stdo,*) '         only for the external stars) these stars as point'
              write(stdo,*) '         sources. If the star is located far enough out of'
              write(stdo,*) '         the grid, this yields only minimal errors, so it is fine.'
              iwarned = .true.
           endif
           if(dum.lt.10*star_r(istar)) then
              write(stdo,*) ' '
              write(stdo,*) '         BUT: We found external star that is too much in near field.'
              write(stdo,*) '         This means that our approximations are no longer valid.'
              write(stdo,*) '         Either move the star further out, or switch off the'
              write(stdo,*) '         weighted photons mode, or put the star inside the grid...'
              write(stdo,*) '         Aborting... Sorry!'
              stop
           endif
        endif
     enddo
  endif
  !
  ! If the stars are located outside the containing sphere (outside the
  ! grid), and if the weighted photons mode is switched on, then we 
  ! calculate the fraction of coverage of the grid w.r.t. the star
  !
  if(mc_weighted_photons) then
     do istar=1,nstars
        star_fraclum(istar) = 1.d0
        dum = sqrt((star_pos(1,istar)-grid_contsph_x)**2+&
                   (star_pos(2,istar)-grid_contsph_y)**2+&
                   (star_pos(3,istar)-grid_contsph_z)**2)
        if(dum.ge.grid_contsph_r) then
           if(dum.le.grid_contsph_r+star_r(istar)) then
              write(stdo,*) 'ERROR: Star outside of grid still crosses'
              write(stdo,*) '       the grid (or at least the containing '
              write(stdo,*) '       sphere around the grid). Aborting...'
              stop
           endif
           if(grid_contsph_r.gt.1d-4*dum) then
              star_fraclum(istar) = 0.5d0 * ( 1.d0 - &
                   sqrt(1.d0-(grid_contsph_r/dum)**2) )
           else
              star_fraclum(istar) = 0.25d0 * (grid_contsph_r/dum)**2
           endif
        endif
     enddo
  endif
  !
  ! Make sure that the dust optical data is read, but only read if not yet
  ! read before.
  !
  call read_dustdata(1)
  !
  ! Make sure that the dust density is read, but only read it from file
  ! if it is not yet in memory.
  !
  call read_dust_density(1)
  !
  ! Get some information about the grid
  !
  if(igrid_type.eq.101) then
     !
     ! The Delaunay triangulated grid
     !
     write(stdo,*) 'NOT YET READY: Delaunay unstructured grid'
     stop
  endif
  if(igrid_type.eq.201) then
     !
     ! The Voronoi grid
     !
     write(stdo,*) 'NOT YET READY: Voronoi unstructured grid'
     stop
  endif
  !
  ! For which kind of Monte Carlo simulation are we going to initialize?
  !
  if(mcaction.eq.1) then
     !
     ! Bjorkman & Wood thermal Monte Carlo
     !
     !
     ! Allocate the dust temperature array
     !
     if(allocated(dusttemp)) deallocate(dusttemp)
     allocate(dusttemp(1:dust_nr_species,1:nrcells),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate dusttemp array.'
        stop 
     endif
     dusttemp(:,:) = 0.d0
     !
     ! OpenMP Parallellization: Allocate the lock array
     !
     !$ if(allocated(lock)) deallocate(lock)
     !$ allocate(lock(1:nrcells),STAT=ierr)
     !$ if(ierr.ne.0) then
     !$    write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate lock array.'
     !$    stop
     !$ endif
     !
     ! The Bjorkman & Wood algorithm is always done on the global
     ! frequency array. But for consistency (because some arrays are
     ! being used by the B&W thermal Monte Carlo as well as by the
     ! frequency-by-frequency scattering Monte Carlo) we will allocate
     ! the mc_frequencies(:) array nevertheless, and fill it with a copy
     ! of the freq_nu(:) array.
     !
     allocate(mc_frequencies(1:freq_nr),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate mc_frequencies()'
        stop 
     endif
     mc_frequencies(:) = freq_nu(:)
     mc_nrfreq = freq_nr
     !
     ! Allocate the arrays for the cumulative spectra of
     ! the thermal boundaries
     !
     allocate(mc_bc_cumspec(1:freq_nr+1,2,3),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate mc_bc_cumspec()'
        stop 
     endif
     !
     ! Allocate the average output spectrum array
     !
     allocate(mc_average_spectrum(1:freq_nr),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate mc_average_spectrum()'
        stop 
     endif
     mc_average_spectrum(:) = 0.d0
     allocate(mc_integerspec(1:freq_nr),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate mc_integerspec()'
        stop 
     endif
     mc_integerspec(:) = 0
     !
     ! Allocate the cumulative energy arrays
     !
     allocate(mc_cumulener(1:dust_nr_species,1:nrcellsmax),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate mc_cumulener()'
        stop 
     endif
     mc_cumulener(:,:) = 0.d0
     if(params%ifast.ne.0) then
        allocate(mc_cumulener_bk(1:dust_nr_species,1:nrcellsmax),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate mc_cumulener_bk()'
           stop 
        endif
        mc_cumulener_bk(:,:) = 0.d0
     endif
     !
     ! If the MRW method is activated, then we need a few arrays
     !
     if(params%mod_random_walk) then
        !
        ! Allocate the MRW arrays
        !
        allocate(mrw_cumulener_bk(1:nrcellsmax),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in Monte Carlo Module: MRW arrays could not be allocated.'
           stop
        endif
        mrw_cumulener_bk(:) = 0.d0
        !$OMP PARALLEL
        allocate(mrw_alpha_tot_ross(1:nrcellsmax),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in Monte Carlo Module: MRW arrays could not be allocated.'
           stop
        endif
        mrw_alpha_tot_ross(:) = 0.d0
        !$OMP END PARALLEL
        allocate(mrw_cell_uses_mrw(1:nrcellsmax),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in Monte Carlo Module: MRW arrays could not be allocated.'
           stop
        endif
        mrw_cell_uses_mrw(:) = .false.
        !$OMP PARALLEL
        allocate(alpha_a_pm(1:dust_nr_species),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in Monte Carlo Module: MRW arrays could not be allocated.'
           stop
        endif
        alpha_a_pm(:) = 0.d0
        allocate(mrw_dcumen(1:dust_nr_species),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in Monte Carlo Module: MRW arrays could not be allocated.'
           stop
        endif
        mrw_dcumen(:) = 0.d0
        !$OMP END PARALLEL
        !
        ! Make sure that the tauthres is at least 2xgamma
        !
        if(params%mrw_tauthres.lt.2*params%mrw_gamma) then
           write(stdo,*) 'MRW: Tauthres was less than 2*gamma. Putting tauthres = 2*gamma.'
           params%mrw_tauthres=2*params%mrw_gamma
        endif
     endif
     !
     ! Allocate the partial energy arrays
     !
     !$OMP PARALLEL
     allocate(mc_enerpart(1:dust_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate mc_enerpart()'
        stop 
     endif
     mc_enerpart(:) = 0.d0
     !$OMP END PARALLEL
     !
     ! Check the alignment mode
     !
     if(alignment_mode.gt.0) then
        write(stdo,*) 'MAJOR WARNING: Grain alignment is not yet included in ', &
             'Bjorkman and Wood thermal Monte Carlo! So RADMC-3D will ignore ', &
             'grain alignment in thermal MC even though alignment_mode.gt.0!'
     endif
     !
  elseif((mcaction.eq.101).or.(mcaction.eq.102)) then
     !
     ! The frequency-by-frequency scattering Monte Carlo simulation
     !
     ! OpenMP Parallellization: Allocate the lock array
     !
     !$ if(allocated(lock)) deallocate(lock)
     !$ allocate(lock(1:nrcells),STAT=ierr)
     !$ if(ierr.ne.0) then
     !$    write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate lock array.'
     !$    stop
     !$ endif
     !
     ! Check that the dust temperatures are there
     !
     call read_dust_temperature(1)
     !
     ! Check that the Monte Carlo frequencies array is set by the caller
     !
     if((mc_nrfreq.le.0).or.(.not.allocated(mc_frequencies))) then
        write(stdo,*) 'ERROR: When you want to do Monte Carlo scattering'
        write(stdo,*) '       you must set the mc_frequencies(:).'
        stop
     endif
     !
     ! Allocate the average output spectrum array
     !
     allocate(mc_average_spectrum(1:mc_nrfreq),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate mc_average_spectrum()'
        stop 
     endif
     mc_average_spectrum(:) = 0.d0
     allocate(mc_integerspec(1:mc_nrfreq),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate mc_integerspec()'
        stop 
     endif
     mc_integerspec(:) = 0
     !
     ! If scattering source function is to be calculated, then allocate the
     ! scattering source function array.  In the future we may include MC
     ! modes where this array is not required.
     !
     if(mcaction.eq.101) then
        if(mcscat_nrdirs.le.0) then
           write(stdo,*) 'ERROR: mcscat_nrdirs.le.0 while initiating MC scat...'
           stop
        endif
        allocate(mcscat_scatsrc_iquv(1:mc_nrfreq,1:nrcellsmax,1:4,1:mcscat_nrdirs),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate mcscat_scatsrc_iquv()'
           write(stdo,*) '      As this can be a huge array, perhaps you have memory limitations?'
           write(stdo,*) '      This array requires ',mcscat_nrdirs*mc_nrfreq* &
                (nrcellsmax*4.d0)/(1024.*1024),' Mbytes'
           stop 
        else
           mcscat_scatsrc_iquv(:,:,:,:) = 0.d0
        endif
     else
        allocate(mcscat_meanint(1:mc_nrfreq,1:nrcellsmax),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate mcscat_meanint()'
           write(stdo,*) '      As this can be a huge array, perhaps you have memory limitations?'
           write(stdo,*) '      This array requires ',mc_nrfreq* &
                (nrcellsmax*4.d0)/(1024.*1024),' Mbytes'
           stop 
        else
           mcscat_meanint(:,:) = 0.d0
        endif
     endif
     !
     ! 
     ! Check if the nr of viewing directions is set by the caller
     ! Note that even if isotropic scattering is used, nrdirs must be
     ! 1, because if nrdirs is 0 there will be no scattering source
     ! function. In the future we may include Monte Carlo scattering modes
     ! that do not need scattering source functions. 
     !
     if(mcscat_nrdirs.gt.1) then
        if(.not.(allocated(mcscat_dirs))) then
           write(stdo,*) 'ERROR: nrdirs>1, but no directions array set'
           stop
        endif
        if(scattering_mode.eq.1) then
           write(stdo,*) 'ERROR: Isotropic scattering mode must be with nrdirs=1.'
           stop
        endif
     endif
!----TO-ADD----
!!
!!   Allocate the skipping array for skipping extremely optically thick 
!!   cells in the scattering Monte Carlo.
!!
!    allocate(mc_scat_src_skip(1:nrcellsmax),STAT=ierr)
!    if(ierr.ne.0) then
!       write(stdo,*) 'ERROR in Monte Carlo: Could not allocate skipping array.'
!       stop
!    endif
!----TO-ADD----
     !
     ! Check the alignment mode
     !
     if(alignment_mode.gt.0) then
        write(stdo,*) 'MAJOR WARNING: Grain alignment is not yet included in ', &
             'Scattering Monte Carlo! So RADMC-3D will ignore ', &
             'grain alignment in scattering MC even though alignment_mode.gt.0!'
     endif
     !
  else
     stop 7012
  endif
  !
  ! For anisotropic scattering, allocate phasefunction array
  !
  if(mcscat_nrdirs.ge.1) then
     !$OMP PARALLEL
     allocate(mcscat_phasefunc(1:mcscat_nrdirs,1:dust_nr_species),STAT=ierr)
     !$OMP END PARALLEL
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in camera module: Could not allocate phasefunction.'
        stop
     endif
  endif
  !
  ! Allocate arrays for the stellar source templates, and
  ! map the templates onto these frequencies
  !
  if(stellarsrc_nrtemplates.gt.0) then
     allocate(mc_stellarsrc_templates(1:mc_nrfreq,1:stellarsrc_nrtemplates),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in monte carlo module: Could not allocate stellarsrc templates.'
        stop
     endif
     do itempl=1,stellarsrc_nrtemplates
        do inu=1,mc_nrfreq
           call hunt(freq_nu,freq_nr,mc_frequencies(inu),iinu)
           if((iinu.le.0).or.(iinu.gt.freq_nr)) then
              mc_stellarsrc_templates(inu,itempl) = 0.d0
           else
              if(iinu.eq.freq_nr) then
                 if(mc_frequencies(inu).eq.freq_nu(freq_nr)) then
                    mc_stellarsrc_templates(inu,itempl) = stellarsrc_templates(freq_nr,itempl)
                 else
                    mc_stellarsrc_templates(inu,itempl) = 0.d0
                 endif
              else
                 eps = ( mc_frequencies(inu) - freq_nu(iinu) ) / &
                       ( freq_nu(iinu+1) - freq_nu(iinu) )
                 if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 674
                 mc_stellarsrc_templates(inu,itempl) = &
                      (1.d0-eps)*stellarsrc_templates(iinu,itempl) + &
                             eps*stellarsrc_templates(iinu+1,itempl)
              endif
           endif
        enddo
     enddo
  endif
  !
  ! Allocate the opacities used for the MC code. Note that a more complex
  ! and generalized form of the opacities is loaded in the dust_module.f90,
  ! but we want to use a simplified form here, to speed up the code.
  !
  allocate(kappa_a(1:mc_nrfreq,1:dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate kappa_a()'
     stop 
  endif
  kappa_a(:,:) = 0.d0
  allocate(kappa_s(1:mc_nrfreq,1:dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate kappa_s()'
     stop 
  endif
  kappa_s(:,:) = 0.d0
  allocate(kappa_g(1:mc_nrfreq,1:dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate kappa_g()'
     stop 
  endif
  kappa_g(:,:) = 0.d0
  !$OMP PARALLEL
  allocate(alpha_a(1:dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate alpha_a()'
     stop 
  endif
  alpha_a(:) = 0.d0
  allocate(alpha_s(1:dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate alpha_s()'
     stop 
  endif
  alpha_s(:) = 0.d0
  allocate(alphacum(1:dust_nr_species+1),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate alphacum()'
     stop 
  endif
  alphacum(:) = 0.d0
  !$OMP END PARALLEL
  !
  ! If we have tabulated scattering matrices, then allocate 
  ! 
  if(scattering_mode.ge.3) then
     !
     ! Check if a scattering angle grid file is present. If so,
     ! then read the angle grid. If not, then create one.
     !
     inquire(file='scattering_angular_grid.inp',exist=fex)
     if(fex) then
        !
        ! Read the scattering angle grid from the file
        ! scattering_angular_grid.inp
        !
        call read_scat_angular_grid(1)
        !
     else
        !
        ! Make the scattering angle grid internally
        !
        ! Determine the maximum fineness of the angular grids of the
        ! dust opacity files. If scat_munr is already set by the user
        ! in the radmc3d.inp file, then take that instead.
        !
        if(scat_munr.eq.0) then
           do ispec=1,dust_nr_species
              if(dust_kappa_arrays(ispec)%nmu.gt.scat_munr) then
                 scat_munr = dust_kappa_arrays(ispec)%nmu
              endif
           enddo
           if(scat_munr.eq.0) then
              write(stdo,*) 'ERROR: None of the dust opacities has the scattering matrix specified'
              write(stdo,*) '       So scattering_mode should not be larger than 2.'
              write(stdo,*) '       Currently it is: ',scattering_mode
              stop
           endif
        endif
        if(scat_munr.lt.5) then
           write(stdo,*) 'ERROR: Trying to allocate scattering matrix, but nr of'
           write(stdo,*) '       grid points in mu=cos(theta) is <5.'
           stop
        endif
        write(stdo,'(A69,I4,A38)') ' In Monte Carlo we will use a linear-spaced scattering angle grid of ', &
                       scat_munr,' gridpoints between 0 and 180 degrees.'
        !
        ! Make the angular grid in mu
        !
        allocate(scat_mui_grid(1:scat_munr),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate scat_mui_grid'
           stop 
        endif
        allocate(scat_thetai_grid(1:scat_munr),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate scat_thetai_grid'
           stop 
        endif
        !
        do imu=1,scat_munr
           scat_thetai_grid(imu) = pi * (imu-1.d0) / (scat_munr-1.d0)
           scat_mui_grid(imu)    = cos(scat_thetai_grid(imu))
        enddo
        if(abs(scat_mui_grid(1)-1.d0).gt.1d-12) stop 438
        if(abs(scat_mui_grid(scat_munr)+1.d0).gt.1d-12) stop 439
     endif
     !
     ! Now allocate the global Z scattering matrix arrays
     !
     allocate(zmatrix(1:scat_munr,1:mc_nrfreq,1:6,1:dust_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate zmatrix()'
        stop 
     endif
     zmatrix(:,:,:,:) = 0.d0
     allocate(zcumul(1:scat_munr,1:mc_nrfreq,1:2,1:dust_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate zcumul()'
        stop 
     endif
     zcumul(:,:,:,:) = 0.d0
  endif
  !
  ! If grain alignment mode > 0, then we have to allocate and
  ! fill some additional arrays
  !
  if(alignment_mode.gt.0) then
     mc_align_munr = align_munr
     if(allocated(mc_align_mu)) deallocate(mc_align_mu)
     if(allocated(mc_align_orth)) deallocate(mc_align_orth)
     if(allocated(mc_align_para)) deallocate(mc_align_para)
     if(allocated(mc_align_opcumul)) deallocate(mc_align_opcumul)
     allocate(mc_align_mu(mc_align_munr))
     allocate(mc_align_orth(mc_align_munr,mc_nrfreq,dust_nr_species))
     allocate(mc_align_para(mc_align_munr,mc_nrfreq,dust_nr_species))
     ! Note for below: since these values are on gridpoints, not
     ! grid cells, the cumulative quantity also only has munr values
     ! instead of munr+1.
     allocate(mc_align_opcumul(mc_align_munr,mc_nrfreq,dust_nr_species))
     !
     ! Now fill these arrays
     !
     mc_align_mu(:) = align_mui_grid(:)
     do ispec=1,dust_nr_species
        do inu=1,mc_nrfreq
           do imu=1,mc_align_munr
              call find_dust_alignfact_interpol(mc_frequencies(inu), &
                   mc_align_mu(imu),ispec,1,0,.true.,orth,para)
              mc_align_orth(imu,inu,ispec) = orth
              mc_align_para(imu,inu,ispec) = para
           enddo
           mc_align_opcumul(1,inu,ispec) = 0.d0
           do imu=2,mc_align_munr
              mc_align_opcumul(imu,inu,ispec) = mc_align_opcumul(imu-1,inu,ispec) +    &
                   abs(mc_align_mu(imu)-mc_align_mu(imu-1)) * 0.5d0 *                  &
                   ( mc_align_orth(imu-1,inu,ispec) + mc_align_para(imu-1,inu,ispec) + &
                     mc_align_orth(imu,inu,ispec) + mc_align_para(imu,inu,ispec) )
           enddo
           do imu=1,mc_align_munr-1
              mc_align_opcumul(imu,inu,ispec) =               &
                   mc_align_opcumul(imu,inu,ispec) /          &
                   mc_align_opcumul(mc_align_munr,inu,ispec)
           enddo
           mc_align_opcumul(mc_align_munr,inu,ispec) = 1.d0
        enddo
     enddo
  endif
  !
  ! Get the opacities (assume global opacities for the moment)
  !
  temp    = 100.d0     ! Arbitrary temperature...
  do inu=1,mc_nrfreq
     do ispec=1,dust_nr_species
        kappa_a(inu,ispec) = find_dust_kappa_interpol(mc_frequencies(inu),ispec,temp,1,0,0)
        kappa_s(inu,ispec) = find_dust_kappa_interpol(mc_frequencies(inu),ispec,temp,0,1,0)
        kappa_g(inu,ispec) = find_dust_kappa_interpol(mc_frequencies(inu),ispec,temp,0,0,1)
     enddo
  enddo
  !
  ! If scattering_mode is 0, then put kappa_s to zero
  ! 
  if(scattering_mode.eq.0) then
     write(stdo,*) 'In Monte Carlo initialization: found that scattering_mode==0, so putting kappa_scat=0.'
     do inu=1,mc_nrfreq
        do ispec=1,dust_nr_species
           kappa_s(inu,ispec) = 0.d0
        enddo
     enddo
  endif
  !
  ! If we have the scattering matrix, then get the values
  !
  if(scattering_mode.ge.3) then
     do ispec=1,dust_nr_species
        do inu=1,mc_nrfreq
           !
           ! Get the Z matrix for this frequency and species, for
           ! all angles
           !
           do imu=1,scat_munr
              call find_dust_zmatrix_interpol(mc_frequencies(inu), &
                   scat_mui_grid(imu),ispec,dumarr6,.true.)
              do iz=1,6
                 zmatrix(imu,inu,iz,ispec) = dumarr6(iz)
              enddo
           enddo
           !
           ! Compute the integrals
           !
           call polarization_total_scattering_opacity(     &
                scat_munr,scat_mui_grid,                   &
                zmatrix(:,inu,1,ispec),                    &
                zmatrix(:,inu,2,ispec),                    &
                zmatrix(:,inu,3,ispec),                    &
                zmatrix(:,inu,4,ispec),                    &
                zmatrix(:,inu,5,ispec),                    &
                zmatrix(:,inu,6,ispec),                    &
                kapscat,g,                                 &
                zcumul(:,inu,1,ispec),                     &
                zcumul(:,inu,2,ispec))
           !
           ! So a simple check
           !
           if(kapscat.lt.0.d0) then
              write(stdo,*) 'INTERNAL ERROR: Somehow the kappa_scat that '
              write(stdo,*) '   I calculate after interpolation of the '
              write(stdo,*) '   scattering matrix on the MC freq- and mu-grids'
              write(stdo,*) '   is negative...'
              write(stdo,*) kapscat
           endif
           if(kapscat.gt.0.d0) then
              if(abs(kapscat/(kappa_s(inu,ispec)+1d-99)-1.d0).gt.0.1) then
                 write(stdo,*) 'INTERNAL ERROR: Somehow the kappa_scat that '
                 write(stdo,*) '   I calculate after interpolation of the '
                 write(stdo,*) '   scattering matrix on the MC freq- and mu-grids'
                 write(stdo,*) '   is much different from the one I calculated'
                 write(stdo,*) '   in the dust module... Strange! Warn author.'
                 write(stdo,*) kapscat,kappa_s(inu,ispec),inu,ispec
                 stop
              endif
              if(abs(g-kappa_g(inu,ispec)).gt.0.1) then
                 write(stdo,*) 'INTERNAL ERROR: Somehow the g=<mu> that '
                 write(stdo,*) '   I calculate after interpolation of the '
                 write(stdo,*) '   scattering matrix on the MC freq- and mu-grids'
                 write(stdo,*) '   is much different from the one I calculated'
                 write(stdo,*) '   in the dust module... Strange! Warn author.'
                 write(stdo,*) g,kappa_g(inu,ispec),inu,ispec
                 stop
              endif
           else
              if(kappa_s(inu,ispec).ne.0.d0) then
                 write(stdo,*) 'INTERNAL ERROR: Somehow the kappa_scat that '
                 write(stdo,*) '   I calculate after interpolation of the '
                 write(stdo,*) '   scattering matrix on the MC freq- and mu-grids'
                 write(stdo,*) '   is 0 while the one I calculated'
                 write(stdo,*) '   in the dust module is not... Strange! Warn author.'
              endif
           endif
        enddo
     enddo
  endif
  !
  ! Allocate, if necessary, the mc_energy_stars(:) array
  !
  if(mc_weighted_photons.and.(nstars.gt.0)) then
     allocate(mc_energy_stars(1:nstars),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate mc_energy_stars() array'
        stop
     endif
  endif
  !
  ! Allocate, if necessary, the mc_energy_illum(:) array
  !
  if(mc_weighted_photons.and.(nillum.gt.0)) then
     allocate(mc_energy_illum(1:nillum),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate mc_energy_illum() array'
        stop
     endif
  endif
  !
  ! If photon statistics are requested, then allocate the arrays
  !
  if(params%debug_write_stats.ne.0) then
     allocate(mc_ilastphot(1:nrcellsmax),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate mc_ilastphot()'
        stop 
     endif
     mc_ilastphot(:) = 0
     allocate(mc_iphotcount(1:nrcellsmax),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate mc_iphotcount()'
        stop 
     endif
     mc_iphotcount(:) = 0
     mc_iphotcurr = 0
  endif
  !
  ! If the Modified Random Walk method is used, then precalculate the
  ! Planck mean opacities
  !
  if(params%mod_random_walk) then
     call mrw_make_planckopac_dbase(params%mrw_db_ntemp,params%mrw_db_temp0,params%mrw_db_temp1)
  endif
  !
end subroutine montecarlo_init


!--------------------------------------------------------------------------
!                 Finish the Monte Carlo module
!--------------------------------------------------------------------------
subroutine montecarlo_cleanup()
  implicit none
  if(allocated(mc_frequencies)) deallocate(mc_frequencies)
  if(allocated(mcscat_dirs)) deallocate(mcscat_dirs)
  if(allocated(mcscat_svec)) deallocate(mcscat_svec)
  mc_nrfreq      = 0
  mcscat_nrdirs  = 0
  call montecarlo_partial_cleanup()
end subroutine montecarlo_cleanup


!--------------------------------------------------------------------------
!                 Finish the Monte Carlo module
!--------------------------------------------------------------------------
subroutine montecarlo_partial_cleanup()
  implicit none
  if(allocated(mc_average_spectrum)) deallocate(mc_average_spectrum)
  if(allocated(mc_integerspec)) deallocate(mc_integerspec)
  if(allocated(mc_cumulener)) deallocate(mc_cumulener)
  if(allocated(mc_cumulener_bk)) deallocate(mc_cumulener_bk)
  if(allocated(mc_ilastphot)) deallocate(mc_ilastphot)
  if(allocated(mc_iphotcount)) deallocate(mc_iphotcount)
  if(allocated(scat_thetai_grid)) deallocate(scat_thetai_grid)
  if(allocated(scat_mui_grid)) deallocate(scat_mui_grid)
  if(allocated(kappa_a)) deallocate(kappa_a)
  if(allocated(kappa_s)) deallocate(kappa_s)
  if(allocated(kappa_g)) deallocate(kappa_g)
  if(allocated(zmatrix)) deallocate(zmatrix)
  if(allocated(zcumul)) deallocate(zcumul)
  if(allocated(mcscat_scatsrc_iquv)) deallocate(mcscat_scatsrc_iquv)
  if(allocated(mcscat_meanint)) deallocate(mcscat_meanint)
  if(allocated(mc_cumulthermemis)) deallocate(mc_cumulthermemis)
  if(allocated(mc_stellarsrc_templates)) deallocate(mc_stellarsrc_templates)
  if(allocated(mc_energy_stars)) deallocate(mc_energy_stars)
  if(allocated(mc_energy_illum)) deallocate(mc_energy_illum)
  if(allocated(mc_bc_cumspec)) deallocate(mc_bc_cumspec)
  if(allocated(mrw_cumulener_bk)) deallocate(mrw_cumulener_bk)
  if(allocated(mrw_db_temp)) deallocate(mrw_db_temp)
  if(allocated(mrw_db_kappa_abs_planck)) deallocate(mrw_db_kappa_abs_planck)
  if(allocated(mrw_cell_uses_mrw)) deallocate(mrw_cell_uses_mrw)
  if(allocated(mc_align_mu)) deallocate(mc_align_mu)
  if(allocated(mc_align_orth)) deallocate(mc_align_orth)
  if(allocated(mc_align_para)) deallocate(mc_align_para)
  if(allocated(mc_align_opcumul)) deallocate(mc_align_opcumul)
  !
  !$OMP PARALLEL
  if(allocated(mc_enerpart)) deallocate(mc_enerpart)
  if(allocated(alpha_a)) deallocate(alpha_a)
  if(allocated(alpha_s)) deallocate(alpha_s)
  if(allocated(alphacum)) deallocate(alphacum)
  if(allocated(alpha_a_pm)) deallocate(alpha_a_pm)
  if(allocated(mrw_dcumen)) deallocate(mrw_dcumen)
  if(allocated(mrw_alpha_tot_ross)) deallocate(mrw_alpha_tot_ross)
  if(allocated(mcscat_phasefunc)) deallocate(mcscat_phasefunc)
  !$OMP END PARALLEL
  !
  !$ if(allocated(lock)) deallocate(lock)
  mrw_db_ntemp = 0
!----TO-ADD----
! if(allocated(mc_scat_src_skip)) deallocate(mc_scat_src_skip)
!----TO-ADD----
end subroutine montecarlo_partial_cleanup


!----------------------------------------------------------------------------
!          READ THE FREQUENCIES FOR THE MONOCHROMATIC MONTE CARLO
!----------------------------------------------------------------------------
subroutine read_mc_frequencies()
  implicit none
  logical :: flag
  integer :: inu,ierr,i,iformat
  logical :: fex
  !
  ! Deallocate frequencies
  !
  if(allocated(mc_frequencies)) deallocate(mc_frequencies)
  !
  ! Read frequencies / wavelengths from a file
  !
  inquire(file='mcmono_wavelength_micron.inp',exist=fex)
  if(fex) then
     open(unit=1,file='mcmono_wavelength_micron.inp')
     read(1,*) mc_nrfreq
     allocate(mc_frequencies(mc_nrfreq),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in montecarlo module: Could not allocate '
        write(stdo,*) '      mc_frequencies(:).'
        stop 
     endif
     do inu=1,mc_nrfreq
        read(1,*) mc_frequencies(inu)
        mc_frequencies(inu) = 1d4*cc/mc_frequencies(inu)
     enddo
     close(1)
  else
     write(stdo,*) 'ERROR: Cannot find mcmono_wavelength_micron.inp'
     stop
  endif
end subroutine read_mc_frequencies



!--------------------------------------------------------------------------
!       COMPUTE ALL TOTAL LUMINOSITIES AND CUMULATIVE LUMINOSITIES
!--------------------------------------------------------------------------
subroutine montecarlo_compute_total_luminosities(params,incl_therm)
  implicit none
  type(mc_params) :: params
  logical :: incl_therm
  integer :: istar,inu,icell,index,itemplate,ispec,ierr,bc_ilr,bc_idir,illum
  double precision :: dum,ener,surf
  double precision, parameter :: ssb=5.6703d-5   ! Stefan-Boltzmann const  [erg/cm^2/K^4/s]
  !
  ! The stars
  !
  starlumtot = 0.d0
  !
  ! Determine the luminosity of the star(s) in all these wavelength bins
  !
  if(nstars.gt.0) then
     if((igrid_coord.eq.10).or.(igrid_coord.eq.20)) then
        write(stdo,*) 'ERROR: In 1-D plane-parallel or 2-D pencil-parallel mode stars are not allowed.'
        stop
     endif
     if(.not.allocated(star_lum)) stop 8709
     if(.not.allocated(star_lumcum)) stop 8708
     starlumtot = 0.d0
     do istar=1,nstars
        star_lum(istar) = 0.d0
        do inu=1,freq_nr
           dum = pi * star_spec(inu,istar) * fourpi * star_r(istar)**2
           star_lum(istar) = star_lum(istar) + dum * freq_dnu(inu)
        enddo
        starlumtot = starlumtot + star_lum(istar)
     enddo
     if(starlumtot.gt.0.d0) then
        star_lumcum(1) = 0.d0
        do istar=1,nstars
           star_lumcum(istar+1) = star_lum(istar) / starlumtot
        enddo
        star_lumcum(nstars+1) = 1.d0
     endif
  endif
  !
  ! For 1-D plane-parallel mode, instead of stars, you can have 
  ! illumination beams at certain angles and with certain fluxes.
  ! NOTE: The total flux is put into "starlumtot"
  !
  if(nillum.gt.0) then
     if(igrid_coord.ne.10) then
        write(stdo,*) 'ERROR: Illumination beams are only permitted in 1-D plane-parallel mode.'
        stop
     endif
     if(.not.allocated(illum_flux_unprojected)) stop 8711
     do illum=1,nillum
        illum_fluxtot(illum) = 0.d0
        do inu=1,freq_nr
           illum_fluxtot(illum) = illum_fluxtot(illum) + abs(illum_costheta(illum)) * &
                                  illum_flux_unprojected(inu,illum) * freq_dnu(inu)
        enddo
        starlumtot = starlumtot + illum_fluxtot(illum)
     enddo
     illum_fluxtot_cum(1) = 0.d0
     do illum=1,nillum
        illum_fluxtot_cum(illum+1) = illum_fluxtot_cum(illum) + illum_fluxtot(illum)
     enddo
     do illum=1,nillum+1
        illum_fluxtot_cum(illum) = illum_fluxtot_cum(illum) / starlumtot
     enddo
     if(abs(illum_fluxtot_cum(nillum+1)-1.d0).gt.1d-6) stop 2361
     illum_fluxtot_cum(nillum+1) = 1.d0
  endif
  !
  ! Determine the luminosity of the external radiation field.  This
  ! luminosity is defined as the integral of the incoming part of the
  ! external flux at a sphere encompassing the entire grid.  Of course, if
  ! you take the size of the grid larger (for the same problem), then the
  ! external luminosity increases, because the surface of the encompassing
  ! sphere increases. But most of these photons will then simply escape to
  ! infinity again without even interacting with the matter. Therefore this
  ! scaling is not a major problem. One should keep in mind, though, that it
  ! may strain the photon-statistics when you have too larger outer boundary
  ! (photons then rarely reach the inner zones).
  !
  if(incl_extlum.ne.0) then 
     !
     ! Determine the luminosity of the external radiation field
     !
     extlumtot = 0.d0
     do inu=1,freq_nr
        dum = pi * extlum_intens(inu) * fourpi * grid_contsph_r**2
        extlumtot = extlumtot + dum * freq_dnu(inu)
     enddo
     extlum_intens_cum(1) = 0.d0
     do inu=1,freq_nr
        extlum_intens_cum(inu+1) = extlum_intens_cum(inu) + &
             pi * extlum_intens(inu) * fourpi *             &
             grid_contsph_r**2 * freq_dnu(inu) / extlumtot
     enddo
     if(abs(extlum_intens_cum(freq_nr+1)-1.d0).gt.1d-10) stop 82791
     !
     ! Note: in case of mirror symmetry, here we do NOT need to 
     ! multiply by 2, because the above computation already is 
     ! done for both halves (see BUGFIX 17.09.2016).
     ! 
  else
     extlumtot = 0.d0
  endif
  !
  ! Compute the luminosity of the internally produced energy heatsource
  ! and make the cumulative array for the heatsource. Note that the
  ! cumulative array heatsource_cum is indexed with icell=1,nrcells+1,
  ! while most other arrays are indexed through the indexing array
  ! cellindex in order to allow adding and removal of cells without
  ! having to reshuffle data to fill in "holes" in the arrays when
  ! cells are removed.
  !
  if(incl_heatsource.ne.0) then
     !
     ! Check if the viscous heating array is allocated
     !
     if(.not.allocated(heatsource)) then
        write(stdo,*) 'INTERNAL ERROR: Internal heat source array not allocated'
        stop
     endif
     !
     ! Allocate the heatsource_cum(:) array
     ! 
     if(allocated(heatsource_cum)) deallocate(heatsource_cum)
     allocate(heatsource_cum(1:nrcells+1),STAT=ierr)
     if(ierr.gt.0) then
        write(stdo,*) 'ERROR: Could not allocate heatsource_cum(:)'
        stop
     endif
     !
     ! Now generate the cumulative array and calculate the total
     ! total luminosity of the viscous heating source
     !
     heatsourcelumtot = 0.d0
     do icell=1,nrcells
        index = cellindex(icell)
        heatsourcelumtot = heatsourcelumtot + heatsource(index)*cellvolume(index)
     enddo
     if(heatsourcelumtot.gt.0.d0) then 
        dum = 0.d0
        do icell=1,nrcells
           index = cellindex(icell)
           heatsource_cum(icell) = dum
           dum = dum + heatsource(index)*cellvolume(index)/heatsourcelumtot
        enddo
        heatsource_cum(nrcells+1) = dum
        if(abs(dum-1.d0).gt.1d-12) then
           write(stdo,*) 'INTERNAL ERROR: qplus somehow doesnt add up'
           stop
        endif
     endif
     !
     ! If mirror symmetry in the z=0 plane, then the actual heat source
     ! luminosity is twice this value
     ! BUGFIX 17.09.2016
     !
     if(igrid_mirror.eq.1) then
        heatsourcelumtot = 2.0d0 * heatsourcelumtot
     endif
  else
     heatsourcelumtot = 0.d0
  endif
  !
  ! Compute the luminosity of the locally produced continuous stellar source.
  ! This is used in simulations of entire galaxies, where individual stars
  ! are too many to be modeled as discrete objects. 
  !
  ! The stellar source is specified as density for each stellar template.
  ! That is: you pre-specify a number (1 or larger) of star spectra that
  ! are luminosity L_nu per gram of the star. Then in each cell of the grid
  ! for each of the stellar templates a density of stars in g/cm^3 is specified
  ! and the local emissivity of the stars are then determined according to
  ! the product of the stellar density and the template spectrum. Note that
  ! a template can also be a galaxy template (consisting of a mix of many type
  ! of stars). Only if you wish to have different stellar populations in 
  ! different regions will you have to define multiple templates.
  !
  if(incl_stellarsrc.ne.0) then
     if(stellarsrc_nrtemplates.le.0) then
        write(stdo,*) 'ERROR: No stellar templates available for the distributed'
        write(stdo,*) '       stellar sources (used for galaxy simulations)'
        stop
     endif
     if(.not.allocated(stellarsrc_dens)) then
        write(stdo,*) 'INTERNAL ERROR: stellar source density not allocated'
        stop
     endif
     if(.not.allocated(stellarsrc_lumcum)) then
        write(stdo,*) 'INTERNAL ERROR: stellar source cumulative array not allocated'
        stop
     endif
     stellarsrclumtot = 0.d0
     do icell=1,nrcells
        !
        ! NOTE: All cumulative arrays are indexed by "icell", all others by "index"
        ! Note: The templateslum is already including fourpi
        !
        index = cellindex(icell)
        stellarsrc_lumcum(icell) = stellarsrclumtot
        dum = 0.d0
        do itemplate=1,stellarsrc_nrtemplates
           dum = dum + stellarsrc_dens(itemplate,index) *               &
                       stellarsrc_templateslum(itemplate)
        enddo
        stellarsrclumtot = stellarsrclumtot + dum * cellvolume(index)
     enddo
     if(stellarsrclumtot.gt.0.d0) then 
        do icell=1,nrcells
           stellarsrc_lumcum(icell) = stellarsrc_lumcum(icell) / stellarsrclumtot
        enddo
        stellarsrc_lumcum(nrcells+1) = 1.d0
     endif
     !
     ! If mirror symmetry in the z=0 plane, then the actual stellar source
     ! luminosity is twice this value
     ! BUGFIX 17.09.2016
     !
     if(igrid_mirror.eq.1) then
        stellarsrclumtot = 2.0d0 * stellarsrclumtot
     endif
  else
     stellarsrclumtot = 0.d0
  endif
  !
  ! Compute the luminosity of the quantum heated grains
  !
  ! NOTE: All cumulative arrays are indexed by "icell", all others by "index"
  !
  emisquanttot = 0.d0
  if(incl_quantum.ne.0) then
     if(incl_quantum.gt.0) then
        write(stdo,*) 'HALTING: For now only the externally computed quantum '
        write(stdo,*) '         mode is available (incl_quantum<0).'
        stop
     endif
     if(.not.allocated(emisquant_loccum)) then
        write(stdo,*) 'ERROR: When quantum heating is included,'
        write(stdo,*) '       the array emisquant_loccum must be allocated.'
        stop
     endif
     do icell=1,nrcells
        index = cellindex(icell)
        dum = 0.d0
        do inu=1,freq_nr
           emisquant_loccum(inu,index) = dum
           dum = dum + emisquant(inu,index) * fourpi *                   &
                       cellvolume(index) * freq_dnu(inu)
        enddo
        emisquant_loctot(index)           = dum
        emisquanttot                      = emisquanttot + dum
        if(dum.gt.0.d0) then
           dum = 1.d0 / dum
           do inu=1,freq_nr
              emisquant_loccum(inu,index) = emisquant_loccum(inu,index) * dum
           enddo
           emisquant_loccum(freq_nr+1,index) = 1.d0
        else
           do inu=1,freq_nr+1
              emisquant_loccum(inu,index) = 0.d0
           enddo
        endif
     enddo
     if(emisquanttot.gt.0.d0) then 
        dum = 0.d0
        do icell=1,nrcells
           index = cellindex(icell)
           emisquant_cum(icell) = dum
           dum = dum + emisquant_loctot(index)
        enddo
        emisquant_cum(nrcells+1) = dum
        if(abs(dum-1.d0).gt.1d-12) then
           write(stdo,*) 'INTERNAL ERROR: emisquant somehow doesnt add up'
           stop
        endif
        !
        ! If mirror symmetry in the z=0 plane, then the actual quantum source
        ! luminosity is twice this value
        ! BUGFIX 17.09.2016
        !
        if(igrid_mirror.eq.1) then
           emisquanttot = 2.0d0 * emisquanttot
        endif
        write(stdo,*) 'L_emisquant/L_stars = ',emisquanttot/starlumtot
     else
        do icell=1,nrcells
           emisquant_cum(icell) = 0.d0
        enddo
        emisquant_cum(nrcells+1) = 0.d0
     endif
  endif
  !
  ! Include the thermal emission as well, IF requested
  !
  mc_thermemistot = 0.d0
  if(incl_therm) then
     !
     ! Allocate the mc_cumulthermemis(:) array
     ! 
     if(allocated(mc_cumulthermemis)) deallocate(mc_cumulthermemis)
     allocate(mc_cumulthermemis(1:nrcells+1),STAT=ierr)
     if(ierr.gt.0) then
        write(stdo,*) 'ERROR: Could not allocate mc_cumulthermemis(:)'
        stop
     endif
     !
     ! First compute for each cell the total thermal emissivity in 
     ! erg/s/cm^3 (not per sterradian!)
     !
     mc_cumulthermemis(1) = 0.d0
     do icell=1,nrcells
        index = cellindex(icell)
        ener = 0.d0
        do ispec=1,dust_nr_species
           !
           ! Compute the emissivity and the total emitted energy
           ! per second per cm^3 by this species.
           !
           do inu=1,freq_nr
              cellalpha(inu) = kappa_a(inu,ispec)*dustdens(ispec,index)
           enddo
           ener = ener + absevfunc(dusttemp(ispec,index))
        enddo
        !
        ! Now add the contribution 
        ! NOTE: Not the index but the icell!!
        !
        index = cellindex(icell)
        mc_cumulthermemis(icell+1) = mc_cumulthermemis(icell) +    &
                                     ener * cellvolume(index)
     enddo
     !
     ! Store the total thermal emissivity
     !
     mc_thermemistot = mc_cumulthermemis(nrcells+1)
     !
     ! Now normalize to unity
     !
     mc_cumulthermemis(:) = mc_cumulthermemis(:) / ( mc_thermemistot + 1d-90 )
     !
     ! Now make sure to put the normalization at nrcells+1 to 1
     !
     mc_cumulthermemis(nrcells+1) = 1.d0
     !
     ! If mirror symmetry in the z=0 plane, then the actual therm emis source
     ! luminosity is twice this value
     ! BUGFIX 17.09.2016
     !
     if(igrid_mirror.eq.1) then
        mc_thermemistot = 2.0d0 * mc_thermemistot
     endif
   endif
   !
   ! If thermal boundaries are active, then calculate their
   ! luminosities
   !
   mc_bc_lumall   = 0.d0
   mc_bc_lum(2,3) = 0.d0
   mc_bc_lumcum(1) = 0.d0
   if((incl_thermbc.ne.0).and.(igrid_coord.ge.0).and.(igrid_coord.lt.100)) then
      !
      ! X-boundaries
      !
      bc_idir=1
      if(igrid_coord.eq.10) then
         surf = 0.d0 
      elseif(igrid_coord.eq.20) then
         surf = 0.d0 
      else
         surf = abs(amr_grid_xi(amr_grid_ny+1,2)-amr_grid_xi(1,2)) * &
                abs(amr_grid_xi(amr_grid_nz+1,3)-amr_grid_xi(1,3))
      endif
      do bc_ilr=1,2
         if(thermal_bc_active(bc_ilr,bc_idir)) then
            mc_bc_cumspec(1,bc_ilr,bc_idir) = 0.d0
            do inu=1,freq_nr
               mc_bc_cumspec(inu+1,bc_ilr,bc_idir) =            &
                    mc_bc_cumspec(inu,bc_ilr,bc_idir) +         &
                    pi*bplanck(thermal_bc_temp(bc_ilr,bc_idir), &
                    freq_nu(inu))*freq_dnu(inu)
            enddo
!            mc_bc_lum(bc_ilr,bc_idir) = surf * ssb *     &
!                 thermal_bc_temp(bc_ilr,bc_idir)**4
            mc_bc_lum(bc_ilr,bc_idir) = surf *                  &
                 mc_bc_cumspec(freq_nr+1,bc_ilr,bc_idir)
            dum = 1.d0 / mc_bc_cumspec(freq_nr+1,bc_ilr,bc_idir)
            mc_bc_cumspec(:,bc_ilr,bc_idir) = dum *             &
                 mc_bc_cumspec(:,bc_ilr,bc_idir) 
            mc_bc_cumspec(freq_nr+1,bc_ilr,bc_idir) = 1.d0
            mc_bc_lumall = mc_bc_lumall + mc_bc_lum(bc_ilr,bc_idir)
         endif
      enddo
      mc_bc_lumcum(2) = mc_bc_lumcum(1) + mc_bc_lum(1,bc_idir)
      mc_bc_lumcum(3) = mc_bc_lumcum(2) + mc_bc_lum(2,bc_idir)
      !
      ! Y-boundaries
      !
      bc_idir=2
      if(igrid_coord.eq.10) then
         surf = 0.d0 
      elseif(igrid_coord.eq.20) then
         surf = abs(amr_grid_xi(amr_grid_nz+1,3)-amr_grid_xi(1,3))
      else
         surf = abs(amr_grid_xi(amr_grid_nx+1,1)-amr_grid_xi(1,1)) * &
                abs(amr_grid_xi(amr_grid_nz+1,3)-amr_grid_xi(1,3))
      endif
      do bc_ilr=1,2
         if(thermal_bc_active(bc_ilr,bc_idir)) then
            mc_bc_cumspec(1,bc_ilr,bc_idir) = 0.d0
            do inu=1,freq_nr
               mc_bc_cumspec(inu+1,bc_ilr,bc_idir) =            &
                    mc_bc_cumspec(inu,bc_ilr,bc_idir) +         &
                    pi*bplanck(thermal_bc_temp(bc_ilr,bc_idir), &
                    freq_nu(inu))*freq_dnu(inu)
            enddo
!            mc_bc_lum(bc_ilr,bc_idir) = surf * ssb *     &
!                 thermal_bc_temp(bc_ilr,bc_idir)**4
            mc_bc_lum(bc_ilr,bc_idir) = surf *                  &
                 mc_bc_cumspec(freq_nr+1,bc_ilr,bc_idir)
            dum = 1.d0 / mc_bc_cumspec(freq_nr+1,bc_ilr,bc_idir)
            mc_bc_cumspec(:,bc_ilr,bc_idir) = dum *             &
                 mc_bc_cumspec(:,bc_ilr,bc_idir) 
            mc_bc_cumspec(freq_nr+1,bc_ilr,bc_idir) = 1.d0
            mc_bc_lumall = mc_bc_lumall + mc_bc_lum(bc_ilr,bc_idir)
         endif
      enddo
      mc_bc_lumcum(4) = mc_bc_lumcum(3) + mc_bc_lum(1,bc_idir)
      mc_bc_lumcum(5) = mc_bc_lumcum(4) + mc_bc_lum(2,bc_idir)
      !
      ! Z-boundaries
      !
      bc_idir=3
      if(igrid_coord.eq.10) then
         surf = 1.d0 
      elseif(igrid_coord.eq.20) then
         surf = abs(amr_grid_xi(amr_grid_ny+1,2)-amr_grid_xi(1,2))
      else
         surf = abs(amr_grid_xi(amr_grid_nx+1,1)-amr_grid_xi(1,1)) * &
                abs(amr_grid_xi(amr_grid_ny+1,2)-amr_grid_xi(1,2))
      endif
      do bc_ilr=1,2
         if(thermal_bc_active(bc_ilr,bc_idir)) then
            mc_bc_cumspec(1,bc_ilr,bc_idir) = 0.d0
            do inu=1,freq_nr
               mc_bc_cumspec(inu+1,bc_ilr,bc_idir) =            &
                    mc_bc_cumspec(inu,bc_ilr,bc_idir) +         &
                    pi*bplanck(thermal_bc_temp(bc_ilr,bc_idir), &
                    freq_nu(inu))*freq_dnu(inu)
            enddo
!            mc_bc_lum(bc_ilr,bc_idir) = surf * ssb *     &
!                 thermal_bc_temp(bc_ilr,bc_idir)**4
            mc_bc_lum(bc_ilr,bc_idir) = surf *                  &
                 mc_bc_cumspec(freq_nr+1,bc_ilr,bc_idir)
            dum = 1.d0 / mc_bc_cumspec(freq_nr+1,bc_ilr,bc_idir)
            mc_bc_cumspec(:,bc_ilr,bc_idir) = dum *             &
                 mc_bc_cumspec(:,bc_ilr,bc_idir) 
            mc_bc_cumspec(freq_nr+1,bc_ilr,bc_idir) = 1.d0
            mc_bc_lumall = mc_bc_lumall + mc_bc_lum(bc_ilr,bc_idir)
         endif
      enddo
      mc_bc_lumcum(6) = mc_bc_lumcum(5) + mc_bc_lum(1,bc_idir)
      mc_bc_lumcum(7) = mc_bc_lumcum(6) + mc_bc_lum(2,bc_idir)
      !
      ! If mirror symmetry in the z=0 plane, then the actual boundary
      ! luminosity is twice this value
      ! BUGFIX 17.09.2016
      !
      if(igrid_mirror.eq.1) then
         mc_bc_lumall    = 2.0d0 * mc_bc_lumall
         mc_bc_lumcum(:) = 2.0d0 * mc_bc_lumcum(:)
         mc_bc_lum(:,:)  = 2.0d0 * mc_bc_lum(:,:)
      endif
   endif
   !
end subroutine montecarlo_compute_total_luminosities


!--------------------------------------------------------------------------
!    COMPUTE ALL FREQ-SPECIFIC LUMINOSITIES AND CUMULATIVE LUMINOSITIES
!--------------------------------------------------------------------------
subroutine montecarlo_compute_freqdep_luminosities(params,freq)
  implicit none
  type(mc_params) :: params
  integer :: istar,inu,icell,index,ispec,ierr,itemplate,bc_ilr,bc_idir,illum
  double precision :: ener,freq,eps,temp,dum,surf
  double precision, allocatable :: kabs(:)
!----TO-ADD----
! type(amr_branch), pointer :: b
! double precision :: alp,cellsize,tau
!----TO-ADD----
  !
  ! Do checks
  !
  if(.not.allocated(dusttemp)) stop 4098
  if(.not.allocated(dustdens)) stop 4099
  if(freq.le.0.d0) then
     write(stdo,*) 'ERROR: Find frequency that is .le.0'
     stop
  endif
  !
  ! First find the inu and eps of the freq in the global frequency array
  !
  call hunt(freq_nu,freq_nr,freq,inu)
  if(inu.eq.freq_nr) then
     if(abs(freq_nu(freq_nr)/freq-1.d0).lt.1d-12) then
        inu=freq_nr-1
        eps=1.d0
     endif
  endif
  if(inu.eq.0) then
     if(abs(freq_nu(1)/freq-1.d0).lt.1d-12) then
        inu=1
        eps=0.d0
     endif
  endif
  if((inu.lt.1).or.(inu.ge.freq_nr)) then
     write(stdo,*) 'ERROR in montecarlo module: Wavelength out of range of the global wavelength array'
     write(stdo,*) 'inu = ',inu,', freq_nr = ',freq_nr
     write(stdo,*) 'freq = ',freq,'  nu0,nu1 = ',freq_nu(1),freq_nu(freq_nr)
     stop
  endif
  eps = (freq-freq_nu(inu)) / (freq_nu(inu+1)-freq_nu(inu))
  if((eps.gt.1.d0).or.(eps.lt.0.d0)) stop 7301
  !
  ! Find the dust opacities at this wavelength
  !
  allocate(kabs(1:dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate kabs(:)'
     stop 
  endif
  temp = 100   ! Dummy temperature
  do ispec=1,dust_nr_species
     kabs(ispec) = find_dust_kappa_interpol(freq,ispec,temp,1,0,0)
  enddo
  !
  ! The stars
  !
  starlumtot = 0.d0
  !
  ! Determine the luminosity of the star(s) at this wavelength
  !
  if(nstars.gt.0) then
     if((igrid_coord.eq.10).or.(igrid_coord.eq.20)) then
        write(stdo,*) 'ERROR: In 1-D plane-parallel or 2-D pencil-parallel mode stars are not allowed.'
        stop
     endif
     if(.not.allocated(star_lum)) stop 8709
     if(.not.allocated(star_lumcum)) stop 8708
     do istar=1,nstars
        star_lum(istar) = (1.d0-eps) * pi * star_spec(inu,istar) * fourpi * star_r(istar)**2 + &
                                 eps * pi * star_spec(inu+1,istar) * fourpi * star_r(istar)**2
        starlumtot = starlumtot + star_lum(istar)
     enddo
     if(starlumtot.gt.0.d0) then
        star_lumcum(1) = 0.d0
        do istar=1,nstars
           star_lumcum(istar+1) = star_lum(istar) / starlumtot
        enddo
        star_lumcum(nstars+1) = 1.d0
     endif
  endif
  !
  ! For 1-D plane-parallel mode, instead of stars, you can have 
  ! illumination beams at certain angles and with certain fluxes.
  ! NOTE: The total flux is put into "starlumtot"
  !
  if(nillum.gt.0) then
     if(igrid_coord.ne.10) then
        write(stdo,*) 'ERROR: Illumination beams are only permitted in 1-D plane-parallel mode.'
        stop
     endif
     if(.not.allocated(illum_flux_unprojected)) stop 8711
     do illum=1,nillum
        illum_fluxtot(illum) = abs(illum_costheta(illum)) * (  &
             (1.d0-eps) * illum_flux_unprojected(inu,illum) +  &
                   eps  * illum_flux_unprojected(inu+1,illum) )
        starlumtot = starlumtot + illum_fluxtot(illum)
     enddo
     illum_fluxtot_cum(1) = 0.d0
     do illum=1,nillum
        illum_fluxtot_cum(illum+1) = illum_fluxtot_cum(illum) + illum_fluxtot(illum)
     enddo
     do illum=1,nillum+1
        illum_fluxtot_cum(illum) = illum_fluxtot_cum(illum) / starlumtot
     enddo
     if(abs(illum_fluxtot_cum(nillum+1)-1.d0).gt.1d-6) stop 2361
     illum_fluxtot_cum(nillum+1) = 1.d0
  endif
  !
  ! Determine the luminosity of the external radiation field.  This
  ! luminosity is defined as the integral of the incoming part of the
  ! external flux at a sphere encompassing the entire grid.  Of course, if
  ! you take the size of the grid larger (for the same problem), then the
  ! external luminosity increases, because the surface of the encompassing
  ! sphere increases. But most of these photons will then simply escape to
  ! infinity again without even interacting with the matter. Therefore this
  ! scaling is not a major problem. One should keep in mind, though, that it
  ! may strain the photon-statistics when you have too large outer boundary
  ! (photons then rarely reach the inner zones).
  !
  if(incl_extlum.ne.0) then 
     !
     ! Determine the luminosity of the external radiation field
     !
     extlumtot = (1.d0-eps) * pi * extlum_intens(inu) * fourpi * grid_contsph_r**2 + &
                        eps * pi * extlum_intens(inu+1) * fourpi * grid_contsph_r**2
     !
     ! Note: in case of mirror symmetry, here we do NOT need to 
     ! multiply by 2, because the above computation already is 
     ! done for both halves (see BUGFIX 17.09.2016).
     ! 
  else
     extlumtot = 0.d0
  endif
  !
  ! Compute the luminosity of the locally produced continuous stellar source.
  ! This is used in simulations of entire galaxies, where individual stars
  ! are too many to be modeled as discrete objects. 
  !
  ! The stellar source is specified as density for each stellar template.
  ! That is: you pre-specify a number (1 or larger) of star spectra that
  ! are luminosity L_nu per gram of the star. Then in each cell of the grid
  ! for each of the stellar templates a density of stars in g/cm^3 is specified
  ! and the local emissivity of the stars are then determined according to
  ! the product of the stellar density and the template spectrum. Note that
  ! a template can also be a galaxy template (consisting of a mix of many type
  ! of stars). Only if you wish to have different stellar populations in 
  ! different regions will you have to define multiple templates.
  !
  stellarsrclumtot = 0.d0
  if(incl_stellarsrc.ne.0) then
     if(stellarsrc_nrtemplates.le.0) then
        write(stdo,*) 'ERROR: No stellar templates available for the distributed'
        write(stdo,*) '       stellar sources (used for galaxy simulations)'
        stop
     endif
     if(.not.allocated(stellarsrc_dens)) then
        write(stdo,*) 'INTERNAL ERROR: stellar source density not allocated'
        stop
     endif
     if(.not.allocated(stellarsrc_templateslum)) then
        write(stdo,*) 'INTERNAL ERROR: stellar source templateslum not allocated'
        stop
     endif
     if(.not.allocated(stellarsrc_lumcum)) then
        write(stdo,*) 'INTERNAL ERROR: stellar source cumulative array not allocated'
        stop
     endif
     stellarsrclumtot = 0.d0
     do icell=1,nrcells
        !
        ! NOTE: All cumulative arrays are indexed by "icell", all others by "index"
        !
        index = cellindex(icell)
        stellarsrc_lumcum(icell) = stellarsrclumtot
        dum = 0.d0
        do itemplate=1,stellarsrc_nrtemplates
           dum = dum + stellarsrc_dens(itemplate,index) * fourpi *             &
                 ( (1-eps) * stellarsrc_templates(inu,itemplate) +             &
                       eps * stellarsrc_templates(inu+1,itemplate) )
        enddo
        stellarsrclumtot = stellarsrclumtot + dum * cellvolume(index)
     enddo
     if(stellarsrclumtot.gt.0.d0) then 
        do icell=1,nrcells
           stellarsrc_lumcum(icell) = stellarsrc_lumcum(icell) / stellarsrclumtot
        enddo
        stellarsrc_lumcum(nrcells+1) = 1.d0
     endif
     !
     ! If mirror symmetry in the z=0 plane, then the actual stellar source
     ! luminosity is twice this value
     ! BUGFIX 17.09.2016
     !
     if(igrid_mirror.eq.1) then
        stellarsrclumtot = 2.0d0 * stellarsrclumtot
     endif
  else
     stellarsrclumtot = 0.d0
  endif
  !
  ! Compute the luminosity of the quantum heated grains
  !
  ! NOTE: All cumulative arrays are indexed by "icell", all others by "index"
  !
  emisquanttot = 0.d0
!!
!! ******* NOT YET READY *******
!!
!!  if(incl_quantum.ne.0) then
!!     if(incl_quantum.gt.0) then
!!        write(stdo,*) 'HALTING: For now only the externally computed quantum '
!!        write(stdo,*) '         mode is available (incl_quantum<0).'
!!        stop
!!     endif
!!     do icell=1,nrcells
!!        index = cellindex(icell)
!!        dum = 0.d0
!!        do inu=1,freq_nr
!!           emisquant_loccum(inu,index) = dum
!!           dum = dum + emisquant(inu,index) * fourpi *                   &
!!                       cellvolume(index) * freq_dnu(inu)
!!        enddo
!!        emisquant_loctot(index)           = dum
!!        emisquanttot                      = emisquanttot + dum
!!        if(dum.gt.0.d0) then
!!           dum = 1.d0 / dum
!!           do inu=1,freq_nr
!!              emisquant_loccum(inu,index) = emisquant_loccum(inu,index) * dum
!!           enddo
!!           emisquant_loccum(freq_nr+1,index) = 1.d0
!!        else
!!           do inu=1,freq_nr+1
!!              emisquant_loccum(inu,index) = 0.d0
!!           enddo
!!        endif
!!     enddo
!!     if(emisquanttot.gt.0.d0) then 
!!        dum = 0.d0
!!        do icell=1,nrcells
!!           index = cellindex(icell)
!!           emisquant_cum(icell) = dum
!!           dum = dum + emisquant_loctot(index)
!!        enddo
!!        emisquant_cum(nrcells+1) = dum
!!        if(abs(dum-1.d0).gt.1d-12) then
!!           write(stdo,*) 'INTERNAL ERROR: emisquant somehow doesnt add up'
!!           stop
!!        endif
!!        !
!!        ! If mirror symmetry in the z=0 plane, then the actual quantum source
!!        ! luminosity is twice this value
!!        ! BUGFIX 17.09.2016
!!        !
!!        if(igrid_mirror.eq.1) then
!!           emisquanttot = 2.0d0 * emisquanttot
!!        endif
!!        write(stdo,*) 'L_emisquant/L_stars = ',emisquanttot/starlumtot
!!     else
!!        do icell=1,nrcells
!!           emisquant_cum(icell) = 0.d0
!!        enddo
!!        emisquant_cum(nrcells+1) = 0.d0
!!     endif
!!  endif
  !
  ! Allocate the mc_cumulthermemis(:) array
  ! 
  if(allocated(mc_cumulthermemis)) deallocate(mc_cumulthermemis)
  allocate(mc_cumulthermemis(1:nrcells+1),STAT=ierr)
  if(ierr.gt.0) then
     write(stdo,*) 'ERROR: Could not allocate mc_cumulthermemis(:)'
     stop
  endif
  !
  ! Include the thermal emission as well
  !
  mc_thermemistot = 0.d0
  !
  ! First compute for each cell the total thermal emissivity in 
  ! erg/s/cm^3 (not per sterradian!)
  !
  mc_cumulthermemis(1) = 0.d0
  do icell=1,nrcells
     index = cellindex(icell)
     ener  = 0.d0
!----TO-ADD----
!    alp  = 0.d0
!----TO-ADD----
     do ispec=1,dust_nr_species
        !
        ! Compute the total emitted energy per second per cm^3 by this
        ! species at this frequency
        !
        ener = ener + dustdens(ispec,index) * kabs(ispec) *        &
                      bplanck(dusttemp(ispec,index),freq)
!----TO-ADD----
!       alp  = alp + dustdens(ispec,index) * kabs(ispec)
!----TO-ADD----
     enddo
     ener = ener * fourpi
!----TO-ADD----
!    b => amr_theleafs(icell)
!    call amr_min_cellsize(b,cellsize) 
!    tau = alp * cellsize
!    if(tau.gt.mc_scat_src_maxtau) then
!       mc_scat_src_skip(index) = .true.
!       ener = 0.d0
!    else
!       mc_scat_src_skip(index) = .false.
!    endif
!----TO-ADD----
     !
     ! Now compute from this the cumulative thermal emission from
     ! this cell compared to the previous cells.
     !
     mc_cumulthermemis(icell+1) = mc_cumulthermemis(icell) +    &
                                  ener * cellvolume(index)
  enddo
  !
  ! Store the total thermal emissivity
  !
  mc_thermemistot = mc_cumulthermemis(nrcells+1)
  !
  ! Now normalize to unity
  !
  mc_cumulthermemis(:) = mc_cumulthermemis(:) / ( mc_thermemistot + 1d-90 )
  !
  ! Fix a small bug, reported by Attila Juhasz on Jan 20, 2011: If there is
  ! no thermal emission, then the normalization must be by hand normalized
  ! to 1.d0
  ! 
  mc_cumulthermemis(nrcells+1) = 1.d0
  !
  ! If mirror symmetry in the z=0 plane, then the actual therm emis source
  ! luminosity is twice this value
  ! BUGFIX 17.09.2016
  !
  if(igrid_mirror.eq.1) then
     mc_thermemistot = 2.0d0 * mc_thermemistot
  endif
  !
  ! If thermal boundaries are active, then calculate their
  ! luminosities
  !
  mc_bc_lumall   = 0.d0
  mc_bc_lum(2,3) = 0.d0
  mc_bc_lumcum(1) = 0.d0
  if((incl_thermbc.ne.0).and.(igrid_coord.ge.0).and.(igrid_coord.lt.100)) then
     !
     ! X-boundaries
     !
     bc_idir=1
     if(igrid_coord.eq.10) then
        surf = 0.d0 
     elseif(igrid_coord.eq.20) then
        surf = 0.d0 
     else
        surf = abs(amr_grid_xi(amr_grid_ny+1,2)-amr_grid_xi(1,2)) * &
               abs(amr_grid_xi(amr_grid_nz+1,3)-amr_grid_xi(1,3))
     endif
     do bc_ilr=1,2
        if(thermal_bc_active(bc_ilr,bc_idir)) then
           mc_bc_lum(bc_ilr,bc_idir) = surf * pi *     &
                bplanck(thermal_bc_temp(bc_ilr,bc_idir),freq)
           mc_bc_lumall = mc_bc_lumall + mc_bc_lum(bc_ilr,bc_idir)
        endif
     enddo
     mc_bc_lumcum(2) = mc_bc_lumcum(1) + mc_bc_lum(1,bc_idir)
     mc_bc_lumcum(3) = mc_bc_lumcum(2) + mc_bc_lum(2,bc_idir)
     !
     ! Y-boundaries
     !
     bc_idir=2
     if(igrid_coord.eq.10) then
        surf = 0.d0 
     elseif(igrid_coord.eq.20) then
        surf = abs(amr_grid_xi(amr_grid_nz+1,3)-amr_grid_xi(1,3))
     else
        surf = abs(amr_grid_xi(amr_grid_nx+1,1)-amr_grid_xi(1,1)) * &
               abs(amr_grid_xi(amr_grid_nz+1,3)-amr_grid_xi(1,3))
     endif
     do bc_ilr=1,2
        if(thermal_bc_active(bc_ilr,bc_idir)) then
           mc_bc_lum(bc_ilr,bc_idir) = surf * pi *     &
                bplanck(thermal_bc_temp(bc_ilr,bc_idir),freq)
           mc_bc_lumall = mc_bc_lumall + mc_bc_lum(bc_ilr,bc_idir)
        endif
     enddo
     mc_bc_lumcum(4) = mc_bc_lumcum(3) + mc_bc_lum(1,bc_idir)
     mc_bc_lumcum(5) = mc_bc_lumcum(4) + mc_bc_lum(2,bc_idir)
     !
     ! Z-boundaries
     !
     bc_idir=3
     if(igrid_coord.eq.10) then
        surf = 1.d0 
     elseif(igrid_coord.eq.20) then
        surf = abs(amr_grid_xi(amr_grid_ny+1,2)-amr_grid_xi(1,2))
     else
        surf = abs(amr_grid_xi(amr_grid_nx+1,1)-amr_grid_xi(1,1)) * &
               abs(amr_grid_xi(amr_grid_ny+1,2)-amr_grid_xi(1,2))
     endif
     do bc_ilr=1,2
        if(thermal_bc_active(bc_ilr,bc_idir)) then
           mc_bc_lum(bc_ilr,bc_idir) = surf * pi *     &
                bplanck(thermal_bc_temp(bc_ilr,bc_idir),freq)
           mc_bc_lumall = mc_bc_lumall + mc_bc_lum(bc_ilr,bc_idir)
        endif
     enddo
     mc_bc_lumcum(6) = mc_bc_lumcum(5) + mc_bc_lum(1,bc_idir)
     mc_bc_lumcum(7) = mc_bc_lumcum(6) + mc_bc_lum(2,bc_idir)
     !
     ! If mirror symmetry in the z=0 plane, then the actual boundary
     ! luminosity is twice this value
     ! BUGFIX 17.09.2016
     !
     if(igrid_mirror.eq.1) then
        mc_bc_lumall    = 2.0d0 * mc_bc_lumall
        mc_bc_lumcum(:) = 2.0d0 * mc_bc_lumcum(:)
        mc_bc_lum(:,:)  = 2.0d0 * mc_bc_lum(:,:)
     endif
  endif
  !
  ! Deallocate temporary array
  !
  if(allocated(kabs)) deallocate(kabs)
  !
end subroutine montecarlo_compute_freqdep_luminosities



!--------------------------------------------------------------------------
!                 MAIN ROUTINE FOR MONTE CARLO
!
! This is the main routine for the Monte Carlo code ala Bjorkman & Wood
! (2001) ApJ 554, 615.
!
! This routine follows the motion of nphot photon packages as they emerge
! from one or more stars and move through the medium. These photons scatter
! and absorb/reemit multiple times until they escape the system. Each
! absorbtion/reemission event also means updating the temperature of the
! cell, and a change in frequency of the photon package. The temperature of
! the cell is updated in such a way that all the photon packages
! before+including this one in total have a correct distribution over
! frequency. This is the trick first described by Bjorkman & Wood. This
! subroutine returns the average spectrum of the outcoming photons, averaged
! over all outgoing directions. 
!
! --------------------------------------------------------------------------
subroutine do_monte_carlo_bjorkmanwood(params,ierror,resetseed)
  implicit none
  type(mc_params) :: params
  doubleprecision :: tempav,dum
  doubleprecision :: aa,fact,seconds
  doubleprecision :: ener,lumtotinv,entotal
  logical :: ievenodd
  logical,optional :: resetseed
  integer :: ierror,ierrpriv,countwrite,index,illum
  integer :: inu,ispec,istar,icell,nsrc,nstarsrc
  integer*8 :: iphot,ipstart,nphot,cnt,cntdump
  integer :: iseeddum,isd,itemplate
  logical :: mc_emergency_break
  !$ integer :: i
  !$ integer OMP_get_num_threads
  !$ integer OMP_get_thread_num
  !$ integer OMP_get_num_procs
  !$ conflict_counter = 0
  !
  ! Some consistency checks
  ! 
  if((nstars.eq.0).and.(incl_stellarsrc.eq.0).and. &
     (incl_extlum.eq.0).and.(incl_heatsource.eq.0).and. &
     (incl_thermbc.eq.0)) then
     write(stdo,*) 'ERROR: There are no sources of photons in this setup.'
     stop 
  endif
  if(params%countwrite.lt.1) then
     write(stdo,*) 'ERROR: Countwrite must be > 0'
     stop 
  endif
  if(params%cntdump.lt.params%countwrite) then
     write(stdo,*) 'ERROR: Cntdump must be >= ',params%countwrite
     stop 
  endif
  !
  ! Initialize the Monte Carlo module by allocating the required arrays
  ! and doing a number of standard checks
  !
  call montecarlo_init(params,ierror,1,resetseed)
  !
  ! Message
  !
  if(scattering_mode.ne.0) then
     write(stdo,*) 'Using dust scattering mode ',scattering_mode
  endif
  !
  ! Set some local integers and logicals
  !
  countwrite = params%countwrite
  cntdump    = params%cntdump
  nphot      = params%nphot_therm
  intplt     = params%iranfreqmode
  ievenodd   = .true.
  !
  ! Flag called 'mc_emergency_break' which is switched to .false. before the big loop.
  ! If an unrecoverable error occurs, it is switched to .true..
  !
  mc_emergency_break = .false.
  ieventcounttot = 0
  mc_visitcell = 0
  mc_revisitcell = 0
  mc_revisitcell_max = 0
  !
  ! If write statistics, then open this file
  !
  if(params%debug_write_stats.ne.0) then
     open(unit=4,file='statistics.out',status='unknown')
     write(4,*) freq_nr
     write(4,*)
  endif
  !
  ! Put to zero various arrays 
  !
  if(incl_quantum.ne.0) then
     if((.not.(allocated(emisquant))).or.(.not.(allocated(miquant))).or. &
        (.not.(allocated(emisquant_cum)))) then
        write(stdo,*) 'INTERNAL ERROR: Including quantum, but quantum not initalized.'
        stop
     endif
     miquant(:,:)       = 0.d0
     emisquant(:,:)     = 0.d0
     emisquant_cum(:)   = 0.d0
  endif
  if(.not.(allocated(dusttemp))) then
     write(stdo,*) 'INTERNAL ERROR: dusttemp arrays not initalized.'
     stop
  endif
  dusttemp(:,:) = 0.d0
  !
  ! Now check if all input array data are OK
  !
  call mc_check_inputdata()
  !
  ! Prepare the database for the emissivity. This is the trick from
  ! Michiel Min that speeds up the code so much, because the CPU-intensive
  ! calculations with integrations over the frequency domain are now done
  ! beforehand and stored in a look-up table.
  !
  write(stdo,*) 'Computing emissivity database...'
  call flush(stdo)
  call make_emiss_dbase(params%ntemp,params%temp0,params%temp1)
  !
  ! Now compute all the luminosities and cumulative luminosity arrays
  !
  write(stdo,*) 'Computing total input luminosities...'
  call flush(stdo)
  call montecarlo_compute_total_luminosities(params,.false.)
  !
  ! Compute the energy per photon package
  ! 
  if(.not.mc_weighted_photons) then
     !
     ! Standard mode: a single photon package energy. Different
     ! luminosities of the various sources are represented by different
     ! numbers of emitted packages.
     !
     ! The energy per photon package is Ltot divided by nr of photons.
     !
     !$OMP PARALLEL
     energy = ( heatsourcelumtot + stellarsrclumtot + emisquanttot + &
                starlumtot + extlumtot + mc_bc_lumall ) / nphot
     !$OMP END PARALLEL
     !
     ! Check 
     !
     if(energy.le.0.d0) then
        write(stdo,*) 'ERROR in setup: There exists no source of photons.'
        write(stdo,*) '      Must have either stars, viscous heating, external'
        write(stdo,*) '      heating or any other netto energy source.'
        stop
     endif
     !
     ! If there is mirror symmetry (a mode only possible for spherical 
     ! coordinates) then divide this by 2 as we only follow the upper
     ! half
     !
     !$OMP PARALLEL
     if(igrid_mirror.eq.1) then
        energy = energy * 0.5d0
     endif
     !$OMP END PARALLEL
     !
     ! Compute the cumulative luminosities of the various different kinds of
     ! sources of photons
     !
     lumtotinv = 1.d0 / ( heatsourcelumtot + emisquanttot + extlumtot + &
                          stellarsrclumtot + mc_bc_lumall + starlumtot )
     mc_cumlum1  = ( heatsourcelumtot ) * lumtotinv
     mc_cumlum2  = ( heatsourcelumtot + emisquanttot ) * lumtotinv
     mc_cumlum3  = ( heatsourcelumtot + emisquanttot  + extlumtot ) * lumtotinv
     mc_cumlum4  = ( heatsourcelumtot + emisquanttot  + extlumtot + stellarsrclumtot ) * lumtotinv
     mc_cumlum5  = ( heatsourcelumtot + emisquanttot  + extlumtot + stellarsrclumtot + mc_bc_lumall ) * lumtotinv
  else
     !
     ! Weighted photon package mode. The different luminosities of the various
     ! sources are represented by different energies per package, but each
     ! source (if non-zero) will still have the same number of photon
     ! packages emitted. In addition to this, if a star is outside the grid,
     ! this weighted photon package mode will focus the photons toward the
     ! grid (so as not to waste photon packages) and adjusts their energy
     ! accordingly.
     ! 
     if(.not.allocated(mc_energy_stars).and.(nstars.gt.0)) stop 9210
     if(.not.allocated(mc_energy_illum).and.(nillum.gt.0)) stop 9210
     !
     ! Count the nr of sources
     !
     nsrc = 0
     if(heatsourcelumtot.gt.0.d0) nsrc = nsrc + 1
     if(emisquanttot.gt.0.d0) nsrc = nsrc + 1
     if(extlumtot.gt.0.d0) nsrc = nsrc + 1
     if(stellarsrclumtot.gt.0.d0) nsrc = nsrc + 1
     if(mc_bc_lumall.gt.0.d0) nsrc = nsrc + 1
     nstarsrc = 0
     if(igrid_coord.ne.10) then
        do istar=1,nstars
           if(star_lum(istar).gt.0.d0) then
              nsrc = nsrc + 1
              nstarsrc = nstarsrc + 1
           endif
        enddo
     else
        do illum=1,nillum
           if(illum_fluxtot(illum).gt.0.d0) then
              nsrc = nsrc + 1
              nstarsrc = nstarsrc + 1
           endif
        enddo
     endif
     if(nsrc.eq.0) then
        write(stdo,*) 'ERROR in setup: There exists no source of photons.'
        write(stdo,*) '      Must have either stars, viscous heating, external'
        write(stdo,*) '      heating or any other netto energy source.'
        stop
     endif
     !
     ! Take into account the mirroring, if active
     !
     if(igrid_mirror.eq.1) then
        fact = 0.5d0
     else
        fact = 1.d0
     endif
     !
     ! Compute the energy per package for each source
     !
     mc_energy_heatsource = fact * heatsourcelumtot * nsrc / nphot
     mc_energy_quant      = fact * emisquanttot * nsrc / nphot
     mc_energy_extlum     = fact * extlumtot * nsrc / nphot
     mc_energy_stellarsrc = fact * stellarsrclumtot * nsrc / nphot
     mc_energy_bc         = fact * mc_bc_lumall * nsrc / nphot
     do istar=1,nstars
        mc_energy_stars(istar) = star_fraclum(istar) * fact *       &
                                 star_lum(istar) * nsrc / nphot
     enddo
     do illum=1,nillum
        mc_energy_illum(illum) = illum_fluxtot(illum) * nsrc / nphot
     enddo
     !
     ! Compute the thesholds for the random assignment of packages
     ! to the various sources
     !
     if(heatsourcelumtot.gt.0.d0) then
        mc_cumlum1  = 1.d0 / nsrc
     else
        mc_cumlum1  = 0.d0
     endif
     if(emisquanttot.gt.0.d0) then
        mc_cumlum2  = mc_cumlum1 + 1.d0 / nsrc
     else
        mc_cumlum2  = mc_cumlum1
     endif
     if(extlumtot.gt.0.d0) then
        mc_cumlum3  = mc_cumlum2 + 1.d0 / nsrc
     else
        mc_cumlum3  = mc_cumlum2
     endif
     if(stellarsrclumtot.gt.0.d0) then
        mc_cumlum4  = mc_cumlum3 + 1.d0 / nsrc
     else
        mc_cumlum4  = mc_cumlum3
     endif
     if(mc_bc_lumall.gt.0.d0) then
        mc_cumlum5  = mc_cumlum4 + 1.d0 / nsrc
     else
        mc_cumlum5  = mc_cumlum4
     endif
     if(nstarsrc.gt.0) then
        if(igrid_coord.ne.10) then
           star_lumcum(1) = 0.d0
           if(star_lum(1).gt.0.d0) then
              star_lumcum(2) = 1.d0 / nstarsrc
           else
              star_lumcum(2) = 0.d0
           endif
           do istar=2,nstars
              if(star_lum(istar).gt.0.d0) then
                 star_lumcum(istar+1) = star_lumcum(istar) + 1.d0 / nstarsrc
              else
                 write(stdo,*) 'ERROR: Cannot treat star with zero luminosity'
                 write(stdo,*) '       Give tiny (but non-zero) luminosity instead.'
                 stop
              endif
           enddo
           if(abs(star_lumcum(nstars+1)-1.d0).gt.1d-14) stop 7311
           star_lumcum(nstars+1) = 1.d0
        else
           illum_fluxtot_cum(1) = 0.d0
           if(illum_fluxtot(1).gt.0.d0) then
              illum_fluxtot_cum(2) = 1.d0 / nstarsrc
           else
              illum_fluxtot_cum(2) = 0.d0
           endif
           do illum=2,nillum
              if(illum_fluxtot(illum).gt.0.d0) then
                 illum_fluxtot_cum(illum+1) = illum_fluxtot_cum(illum) + 1.d0 / nstarsrc
              else
                 write(stdo,*) 'ERROR: Cannot treat illumination beam with zero flux'
                 write(stdo,*) '       Give tiny (but non-zero) flux instead.'
                 stop
              endif
           enddo
           if(abs(illum_fluxtot_cum(nillum+1)-1.d0).gt.1d-14) stop 7311
           illum_fluxtot_cum(nillum+1) = 1.d0
        endif
     endif
  endif
  !   
  ! Open the path file
  ! For debugging only!
  !
  if(params%debug_write_path.eq.1) then
     open(unit=5,file='path.dat',status='unknown')
  endif
  !
  ! Count set to 1
  !
  cnt = 1
  if(params%mod_random_walk) then
     write(stdo,*) 'Modified Random Walk mode is switched ON'
  else
     write(stdo,*) 'Modified Random Walk mode is switched OFF'
  endif
  write(stdo,*) ' '
  write(stdo,*) 'Starting the thermal Monte Carlo simulation....'
  call flush(stdo)
  !
  ! Default
  !
  ipstart = 1
  !
  ! Now, if the irestart.ge.1, then read the latest safety dump
  !
  if(params%irestart.ge.1) then
      !     
      ! Read the last safety backup
      !
      write(stdo,*) '****************************************'
      write(stdo,*) '    !   RESTARTING THE SIMULATION   !   '
      call flush(stdo)
!!!      call read_safety_backup_mctherm(ievenodd,params)
      write(stdo,*) 'ERRROR: Safety backup is not yet built in.'
      stop
      write(stdo,821) iphot
821   format('     !   AT PHOTON NR ',I9,'      !   ')
      !
      ! Now reinitialize the random number generator
      !
      if(params%irestart.eq.1) then 
         !$ write(stdo,*) 'ERROR: Restart not implemented in the parallel version of RADMC-3D.'
         !$ stop
         !
         ! Rerandomize using the iseed from the save file
         !
         write(stdo,*) '    ( ran2: randomizing using iseed )'
         write(stdo,*) '    (       from save file          )'
         if(iseeddum.gt.0) then
            isd=-iseeddum
         else
            isd=iseeddum
         endif
         aa=ran2(isd)
         iseed=isd
      elseif(params%irestart.eq.2) then
         !
         ! Rerandomize using iseed from montecarlo.inp
         !          
         write(stdo,*) '    ( ran2: randomizing using iseed )'
         write(stdo,*) '    (       from input file         )'
         if(iseed.gt.0) then
            isd=-iseed
         else
            isd=iseed
         endif
         aa=ran2(isd)
         iseed=isd
      elseif(params%irestart.eq.3) then
         !
         ! Use exactly the same random set as during last save
         ! (is already done, so do not need to do anything here)
         !
         write(stdo,*) '    ( ran2: restoring exact random  )'
         write(stdo,*) '    (       set as during save      )'
         iseed=iseeddum
      else
         write(stdo,*) 'ERROR: Restart option must be 1,2 or 3'
         stop 13
      endif
      write(stdo,*) '****************************************'
      ipstart=iphot
      cnt = 0
   endif
   !
   ! OpenMP Parallellization:
   !
   !$ seconds = omp_get_wtime()
   !
   !$ do i=1,size(lock)
   !$    call omp_init_lock(lock(i))
   !$ enddo
   !
   ! Note: Here the nr of threads was set in a previous version. 
   !       This is now done in the main.f90
   !       Bugfix Jon Ramsey 22.03.2016
   !
   !!$ Global variables from 'montecarlo_module.f90' used in the subroutine 'do_monte_carlo_bjorkmanwood'
   !
   !$OMP PARALLEL &
   !
   !!$ Local variables from 'do_monte_carlo_bjorkmanwood'
   !
   !$OMP PRIVATE(ierrpriv) &
   !
   !!$ Global variables from other modules used in the subroutine 'do_monte_carlo_bjorkmanwood'
   !!$ or in other subroutine calls within the parallel section
   !
   !!$ REDUCTION: When a variable has been declared as SHARED because all threads need to modify its value,
   !!$ it is necessary to ensure that only one thread at a time is writing/updating the memory
   !!$ location of the considered variable.
   !
   !$OMP REDUCTION(+:mc_integerspec,mc_visitcell,mc_revisitcell,ieventcounttot) &
   !$OMP REDUCTION(max:mc_revisitcell_max)
   !
   !$ id=OMP_get_thread_num()
   !$ nthreads=OMP_get_num_threads()
   !$ write(stdo,*) 'Thread Nr',id,'of',nthreads,'threads in total'
   !$ iseed=-abs(iseed_start+id)
   !$OMP DO SCHEDULE(dynamic)
   !
   ! Launch all the photons
   !
   do iphot=ipstart,nphot
      if (.not.mc_emergency_break) then
         !
         ! For the photon count routines:
         !
         mc_iphotcurr = iphot
         !
         ! Message
         !
         !$OMP CRITICAL
         cnt   = cnt + 1
         if(mod(cnt,countwrite).eq.0) then
            !$   write(stdo,*) 'Thread:',id,'Photon nr:',cnt
            if(.not.mc_openmp_parallel) write(stdo,*) 'Photon nr ',iphot
            call flush(stdo)
         endif
         !$OMP END CRITICAL
         !
         ! Safety dump? [DISABLED]
         !
!         if((cnt.ge.cntdump).and.(iphot.lt.nphot)) then
!!!!         call make_safety_backup_mctherm(ievenodd,params)
!            write(stdo,*) 'WARNING: Wanting to do safety backup, but this is'
!            write(stdo,*) '         not yet implemented.'
!            call flush(stdo)
!  	    !$omp atomic
!            cnt=0
!         endif
         !
         ! Call the walk routine
         !
         !!$ Some critical sections are hidden within this subroutine call
         call walk_full_path_bjorkmanwood(params,ierrpriv)
         !       
         ! If ierrpriv.ne.0 then an error occurred, return with non-zero error code
         !
         if(ierrpriv.ne.0) then
            !$OMP CRITICAL
            mc_emergency_break = .true.
            !$OMP END CRITICAL
            write(stdo,*) '!!!!!!!!!!!!!!!!!!! MC_EMERGENCY_BREAK !!!!!!!!!!!!!!!!!!!!!!'
         endif
         !
         ! Add photon to outgoing overall spectrum
         !
         if(.not.mc_emergency_break) then
            if(.not.mc_photon_destroyed) then 
               mc_integerspec(ray_inu) = mc_integerspec(ray_inu) + 1
            endif
            !
            ! Write debugging stuff
            !
            !      if(params%debug_write_eventcounts.ne.0) then
            !         write(stdo,*) 'Nr of events = ',ieventcount
            !      endif
         endif
      endif
   enddo
   !$OMP END DO
   !$OMP END PARALLEL
   !
   ! OpenMP Parallellization: destroy locks
   !
   !$ do i=1,size(lock)
   !$    call omp_destroy_lock(lock(i))  
   !$ enddo
   !
   ! OpenMP Parallellization: Diagnostic messages
   !
   !$ write(stdo,*)"Elapsed time:",omp_get_wtime() - seconds;
   !!$ open(file='speed_Up_Lock_Thermal',unit=200,position='append')
   !!$ write(200,*)nthreads,"'",omp_get_wtime() - seconds;
   !!$ close(200)
   
   if(mc_emergency_break) then
      return
   endif
   !
   ! Close path file
   !
   if(params%debug_write_path.eq.1) then
      close(5)
   endif
   !
   ! If statistics writing active, close
   !
   if(params%debug_write_stats.ne.0) then
      close(4)
   endif
   !
   ! Write to the radmc_save.info that we're done
   !
   !open(unit=1,file='radmc_save.info')
   !write(1,*) '-'
   !close(1)
   !
   ! Compute the actual spectrum from the integer spectrum
   !
   do inu=1,freq_nr
      mc_average_spectrum(inu) = mc_integerspec(inu) * energy /       &
                              ( twopi * parsec**2 * freq_dnu(inu) )
   enddo
   !
   ! Compute the actual tdust
   !
   if(params%itempdecoup.eq.1) then
      !
      ! If temperatures are decoupled...
      !
      do icell=1,nrcells
         index = cellindex(icell)
         do ispec=1,dust_nr_species
            ener    = mc_cumulener(ispec,index) /                      &
                     ( dustdens(ispec,index) * cellvolume(index) )
            if(ener.gt.0.d0) then
               dusttemp(ispec,index) = compute_dusttemp_energy_bd(ener,ispec)
            else
               dusttemp(ispec,index) = 0.d0
            endif
         enddo
      enddo
   else
      !
      ! If temperatures are coupled...
      !
      do icell=1,nrcells
         index  = cellindex(icell)
         entotal = 0.d0
         do ispec=1,dust_nr_species
            entotal  = entotal + mc_cumulener(ispec,index)
         enddo
         tempav = compute_dusttemp_coupled_bd(dust_nr_species,        &
                             entotal,dustdens(:,index),cellvolume(index))
         do ispec=1,dust_nr_species
            dusttemp(ispec,index) = tempav
         enddo
      enddo
   endif
   !
   ! Do a check
   !
   if(debug_check_all.eq.1) then
      do icell=1,nrcells
         index  = cellindex(icell)
         do ispec=1,dust_nr_species
            if(number_invalid(dusttemp(ispec,index)).ne.0) then
               if(number_invalid(dusttemp(ispec,index)).eq.1) then
                  write(stdo,*) 'ERROR: Discovered INF dust temperature at ', &
                       'ispec=',ispec,' index=',index
                  stop
               else
                  write(stdo,*) 'ERROR: Discovered NAN dust temperature at ', &
                       'ispec=',ispec,' index=',index
                  stop
               endif
            endif
         enddo
      enddo
   endif
   !
   ! Free the emission database
   !
   call free_emiss_dbase()
   !
   ! Done...
   !
end subroutine do_monte_carlo_bjorkmanwood


!--------------------------------------------------------------------------
!         MAIN ROUTINE FOR SINGLE-LAMBDA SCATTERING-ONLY MONTE CARLO
!
! This routine does a Monte Carlo simulation for scattering only. The
! absorption is indeed seen as extinction, i.e. photons are taken
! away. They are taken away in a continuous fashion, so that photon noise
! is limited as much as possible. Photon packages can be emitted from
! all fundamental photon sources (stars, external radiation etc), as well
! as from the dust (in contrast to the Bjorkman & Wood Monte Carlo above).
!
! --------------------------------------------------------------------------
subroutine do_monte_carlo_scattering(params,ierror,resetseed,scatsrc,meanint)
  implicit none
  type(mc_params) :: params
  integer :: ierror,ierrpriv,inu
  doubleprecision :: ener,lumtotinv,temp,freq,fact
  logical :: ievenodd
  logical,optional :: resetseed
  logical,optional :: scatsrc,meanint
  logical :: compute_scatsrc,compute_meanint
  integer*8 :: iphot,nphot,cnt,cntdump
  integer :: countwrite,index
  integer :: ispec,istar,icell,illum
  integer :: iseeddum,isd,itemplate,nsrc,nstarsrc
  logical :: mc_emergency_break
  doubleprecision:: seconds
  !$ integer :: ierr,i
  !$ integer OMP_get_num_threads
  !$ integer OMP_get_thread_num
  !$ integer OMP_get_num_procs
  !$ conflict_counter = 0
  !
  ! Some consistency checks
  ! 
  if(present(scatsrc)) then
     compute_scatsrc=scatsrc 
  else 
     compute_scatsrc=.false.
  endif
  if(present(meanint)) then
     compute_meanint=meanint
  else 
     compute_meanint=.false.
  endif
  if((.not.compute_scatsrc).and.(.not.compute_meanint)) then
     write(stdo,*) 'Monte Carlo Scattering: Must set scatsrc or meanint.'
     stop 3765
  endif
  if(compute_scatsrc.and.compute_meanint) then
     write(stdo,*) 'Monte Carlo Scattering: Cannot set BOTH scatsrc and meanint.'
     stop 3766
  endif
  if(nrcells.le.0) then
     write(stdo,*) 'ERROR: Spatial grid is not set.'
     stop
  endif
  if(params%countwrite.lt.1) then
     write(stdo,*) 'ERROR: Countwrite must be > 0'
     stop 
  endif
  if(params%cntdump.lt.params%countwrite) then
     write(stdo,*) 'ERROR: Cntdump must be >= ',params%countwrite
     stop 
  endif
  if(incl_quantum.ne.0) then
     write(stdo,*) 'ERROR: At the moment quantum heated grains are not included yet'
     stop
  endif
  if(.not.(allocated(dusttemp))) then
     write(stdo,*) 'INTERNAL ERROR: dusttemp arrays not initalized.'
     stop
  endif
  if((mc_nrfreq.le.0).or.(.not.allocated(mc_frequencies))) then
     write(stdo,*) 'ERROR: Cannot do scattering Monte Carlo without'
     write(stdo,*) '       MC frequency array set.'
     stop
  endif
  if(mc_scat_maxtauabs.le.0.d0) then
     write(stdo,*) 'ERROR: The mc_scat_maxtauabs is le 0, meaning that'
     write(stdo,*) '       the scattering Monte Carlo simulation (for images'
     write(stdo,*) '       and spectra) will go wrong... A reasonable value'
     write(stdo,*) '       would be 10.0 or 30.0 or so.'
     stop
  endif
  if(mc_scat_maxtauabs.lt.2.d0) then
     write(stdo,*) 'WARNING: The mc_scat_maxtauabs is lt 2. While this is'
     write(stdo,*) '       formally OK, we advise making this at least 10.'
     stop
  endif
  !
  ! For the selectscat analysis stuff (normally not important)
  !
  if((selectscat_iscat_first.ne.1).or.(selectscat_iscat_last.ne.1000000000)) then
     write(stdo,*) 'ANALYSIS MODE: Scattering source and mean intensity only for iscat between ', &
          selectscat_iscat_first,' and ',selectscat_iscat_last
     write(stdo,*) 'THE RESULTS ARE NOT MEANT FOR PRODUCTION RUNS, ONLY FOR ANALYSIS.'
  endif
  !
  ! Initialize the Monte Carlo module by allocating the required arrays
  ! and doing a number of standard checks
  !
  if(compute_scatsrc) then
     call montecarlo_init(params,ierror,101,resetseed)
  elseif(compute_meanint) then
     call montecarlo_init(params,ierror,102,resetseed)
  else
     stop 6658
  endif
  !
  ! If mc_max_nr_scat_events=0, then no scattering should be
  ! included. But then we must warn the user! Also warn the
  ! user for single scattering or limited-number-scattering.
  ! For no scattering: we can return directly.
  !
  if(mc_max_nr_scat_events.eq.0) then
     write(stdo,*) 'WARNING: Scattering is switched off (mc_max_nr_scat_events=0).'
     return
  endif
  if(mc_max_nr_scat_events.eq.1) then
     write(stdo,*) 'Warning: Only single scattering (mc_max_nr_scat_events=1).'
  endif
  if(mc_max_nr_scat_events.gt.1) then
     write(stdo,*) 'Warning: Only a limited number of scattering events will be included ', &
          '(mc_max_nr_scat_events=',mc_max_nr_scat_events,').'
  endif
  !
  ! Message
  !
  if(scattering_mode.ne.0) then
     write(stdo,*) 'Using dust scattering mode ',scattering_mode
  endif
  !
  ! Set some local integers and logicals
  !
  countwrite = params%countwrite
  cntdump    = params%cntdump
  if(compute_scatsrc) then
     nphot      = params%nphot_scat
  else
     nphot      = params%nphot_mono
  endif
  intplt     = params%iranfreqmode
  ieventcounttot = 0
  !
  ! Flag called 'mc_emergency_break' which is switched to .false. before the big loop.
  ! If an unrecoverable error occurs, it is switched to .true..
  !
  mc_emergency_break = .false.
  !
  ! If write statistics, then open this file
  !
  if(params%debug_write_stats.ne.0) then
     open(unit=4,file='statistics.out',status='unknown')
     write(4,*) freq_nr
     write(4,*)
  endif
  !
  ! Get the opacities (assume global opacities for the moment)
  ! ---> COMMENTED OUT 06.03.2017: Opacities are already installed since
  !      montecarlo_init().
  !
  ! temp    = 100.d0     ! Arbitrary temperature...
  ! do inu=1,mc_nrfreq
  !    freq = mc_frequencies(inu)
  !    do ispec=1,dust_nr_species
  !       kappa_a(inu,ispec) = find_dust_kappa_interpol(freq,ispec,temp,1,0,0)
  !       kappa_s(inu,ispec) = find_dust_kappa_interpol(freq,ispec,temp,0,1,0)
  !       kappa_g(inu,ispec) = find_dust_kappa_interpol(freq,ispec,temp,0,0,1)
  !    enddo
  ! enddo
  !   
  ! Open the path file
  ! For debugging only!
  !
  if(params%debug_write_path.eq.1) then
     open(unit=5,file='path.dat',status='unknown')
  endif
  !
  ! Count set to 1
  !
  cnt = 1
  !
  ! Message
  !
  !write(stdo,*) ' '
  !write(stdo,*) 'Starting scattering Monte Carlo simulation....'
  !call flush(stdo)
  !
  ! Loop over all mc_frequencies
  !
  do inu=1,mc_nrfreq
     !
     ! Message
     !
     write(stdo,*) 'Wavelength nr ',inu,' corresponding to lambda=',  &
               1d4*cc/mc_frequencies(inu),' micron'
     call flush(stdo)
     !
     ! Determine the luminosity of the star(s) in all these wavelength bins
     !
     call montecarlo_compute_freqdep_luminosities(params,mc_frequencies(inu))
     !
     ! Compute the energy per photon package
     ! 
     if(.not.mc_weighted_photons) then
        !
        ! Standard mode: a single photon package energy. Different
        ! luminosities of the various sources are represented by different
        ! numbers of emitted packages.
        !
        ! The energy per photon package is Lnu divided by nr of photons.
        !
        !$OMP PARALLEL
        energy = ( mc_thermemistot + stellarsrclumtot + emisquanttot +           &
                   starlumtot + extlumtot + mc_bc_lumall ) / nphot
        !$OMP END PARALLEL
        !
        ! Check 
        !
        if(energy.le.0.d0) then
           write(stdo,*) 'ERROR in setup: There exists no source of photons.'
           write(stdo,*) '      Must have either stars, viscous heating, external'
           write(stdo,*) '      heating or any other netto energy source.'
           write(stdo,*) '        E_stars      = ',starlumtot
           write(stdo,*) '        E_thermemis  = ',mc_thermemistot
           write(stdo,*) '        E_starsrc    = ',stellarsrclumtot
           write(stdo,*) '        E_quant      = ',emisquanttot
           write(stdo,*) '        E_extlum     = ',extlumtot
           write(stdo,*) '        E_boundary   = ',mc_bc_lumall
           stop
        endif
        !
        ! If there is mirror symmetry (a mode only possible for spherical 
        ! coordinates) then divide this by 2 as we only follow the upper
        ! half
        !
        !$OMP PARALLEL
        if(igrid_mirror.eq.1) then
           energy = energy * 0.5d0
        endif
        !$OMP END PARALLEL
        !
        ! Compute the cumulative luminosities of the various different kinds of
        ! sources of photons
        !
        lumtotinv = 1.d0 / ( mc_thermemistot + emisquanttot + extlumtot +    &
                             stellarsrclumtot + mc_bc_lumall + starlumtot )
        mc_cumlum1  = ( mc_thermemistot ) * lumtotinv
        mc_cumlum2  = ( mc_thermemistot + emisquanttot ) * lumtotinv
        mc_cumlum3  = ( mc_thermemistot + emisquanttot  + extlumtot ) * lumtotinv
        mc_cumlum4  = ( mc_thermemistot + emisquanttot  + extlumtot + stellarsrclumtot ) * lumtotinv
        mc_cumlum5  = ( mc_thermemistot + emisquanttot  + extlumtot + stellarsrclumtot + mc_bc_lumall ) * lumtotinv
     else
        !
        ! Weighted photon package mode. The different luminosities of the various
        ! sources are represented by different energies per package, but each
        ! source (if non-zero) will still have the same number of photon
        ! packages emitted. In addition to this, if a star is outside the grid,
        ! this weighted photon package mode will focus the photons toward the
        ! grid (so as not to waste photon packages) and adjusts their energy
        ! accordingly.
        ! 
        if(.not.allocated(mc_energy_stars).and.(nstars.gt.0)) stop 9211
        if(.not.allocated(mc_energy_illum).and.(nillum.gt.0)) stop 9211
        !
        ! Count the nr of sources
        !
        nsrc = 0
        if(mc_thermemistot.gt.0.d0) nsrc = nsrc + 1
        if(emisquanttot.gt.0.d0) nsrc = nsrc + 1
        if(extlumtot.gt.0.d0) nsrc = nsrc + 1
        if(stellarsrclumtot.gt.0.d0) nsrc = nsrc + 1
        if(mc_bc_lumall.gt.0.d0) nsrc = nsrc + 1
        nstarsrc = 0
        if(igrid_coord.ne.10) then
           do istar=1,nstars
              if(star_lum(istar).gt.0.d0) then
                 nsrc = nsrc + 1
                 nstarsrc = nstarsrc + 1
              endif
           enddo
        else
           do illum=1,nillum
              if(illum_fluxtot(illum).gt.0.d0) then
                 nsrc = nsrc + 1
                 nstarsrc = nstarsrc + 1
              endif
           enddo
        endif
        if(nsrc.eq.0) then
           write(stdo,*) 'ERROR in setup: There exists no source of photons.'
           write(stdo,*) '      Must have either stars, thermal emission, external'
           write(stdo,*) '      heating or any other netto energy source.'
           stop
        endif
        !
        ! Take into account the mirroring, if active
        !
        if(igrid_mirror.eq.1) then
           fact = 0.5d0
        else
           fact = 1.d0
        endif
        !
        ! Compute the energy per package for each source
        !
        mc_energy_thermal    = fact * mc_thermemistot * nsrc / nphot
        mc_energy_quant      = fact * emisquanttot * nsrc / nphot
        mc_energy_extlum     = fact * extlumtot * nsrc / nphot
        mc_energy_stellarsrc = fact * stellarsrclumtot * nsrc / nphot
        mc_energy_bc         = fact * mc_bc_lumall * nsrc / nphot
        do istar=1,nstars
           mc_energy_stars(istar) = star_fraclum(istar) * fact *       &
                                    star_lum(istar) * nsrc / nphot
        enddo
        do illum=1,nillum
           mc_energy_illum(illum) = illum_fluxtot(illum) * nsrc / nphot
        enddo
        !
        ! Compute the thesholds for the random assignment of packages
        ! to the various sources
        !
        if(mc_thermemistot.gt.0.d0) then
           mc_cumlum1  = 1.d0 / nsrc
        else
           mc_cumlum1  = 0.d0
        endif
        if(emisquanttot.gt.0.d0) then
           mc_cumlum2  = mc_cumlum1 + 1.d0 / nsrc
        else
           mc_cumlum2  = mc_cumlum1
        endif
        if(extlumtot.gt.0.d0) then
           mc_cumlum3  = mc_cumlum2 + 1.d0 / nsrc
        else
           mc_cumlum3  = mc_cumlum2
        endif
        if(stellarsrclumtot.gt.0.d0) then
           mc_cumlum4  = mc_cumlum3 + 1.d0 / nsrc
        else
           mc_cumlum4  = mc_cumlum3
        endif
        if(mc_bc_lumall.gt.0.d0) then
           mc_cumlum5  = mc_cumlum4 + 1.d0 / nsrc
        else
           mc_cumlum5  = mc_cumlum4
        endif
        if(nstarsrc.gt.0.d0) then
           if(igrid_coord.ne.10) then
              star_lumcum(1) = 0.d0
              if(star_lum(1).gt.0.d0) then
                 star_lumcum(2) = 1.d0 / nstarsrc
              else
                 star_lumcum(2) = 0.d0
              endif
              do istar=2,nstars
                 if(star_lum(istar).gt.0.d0) then
                    star_lumcum(istar+1) = star_lumcum(istar) + 1.d0 / nstarsrc
                 else
                    star_lumcum(istar+1) = star_lumcum(istar)
                 endif
              enddo
              if(abs(star_lumcum(nstars+1)-1.d0).gt.1d-14) stop 7311
              star_lumcum(nstars+1) = 1.d0
           else
              illum_fluxtot_cum(1) = 0.d0
              if(illum_fluxtot(1).gt.0.d0) then
                 illum_fluxtot_cum(2) = 1.d0 / nstarsrc
              else
                 illum_fluxtot_cum(2) = 0.d0
              endif
              do illum=2,nillum
                 if(illum_fluxtot(illum).gt.0.d0) then
                    illum_fluxtot_cum(illum+1) = illum_fluxtot_cum(illum) + 1.d0 / nstarsrc
                 else
                    write(stdo,*) 'ERROR: Cannot treat illumination beam with zero flux'
                    write(stdo,*) '       Give tiny (but non-zero) flux instead.'
                    stop
                 endif
              enddo
              if(abs(illum_fluxtot_cum(nillum+1)-1.d0).gt.1d-14) stop 7311
              illum_fluxtot_cum(nillum+1) = 1.d0
           endif
        endif
     endif
     !
     ! Compute the energy limit for the Monte Carlo, i.e. when shall we
     ! throw away the photon and continue with the next one? 
     !
     mc_scat_energy_rellimit = exp(-mc_scat_maxtauabs)
     !
     ! Set ray_inu to inu
     !
     ray_inu = inu
     !
     ! OpenMP Parallellization:
     !
     !$ seconds = omp_get_wtime()
     !
     !$ do i=1,size(lock)
     !$    call omp_init_lock(lock(i))
     !$ enddo
     !
     ! Note: Here the nr of threads was set in a previous version. 
     !       This is now done in the main.f90
     !       Bugfix Jon Ramsey 22.03.2016
     !
     !$OMP PARALLEL &
     !!!!!$OMP shared(countwrite,stdo,nphot,nthreads,lock)
     !!!!!$OMP shared(mc_emergency_break,cntdump,cnt,mc_openmp_parallel,inu,iseed_start,params)
     !
     !$OMP PRIVATE(ierrpriv) &
     !
     !!$ Global variables from 'montecarlo_module.f90' used in the subroutine 'do_monte_carlo_scattering'
     !
     !!$ Local variables from 'do_monte_carlo_scattering'
     !$OMP REDUCTION(+:mc_integerspec,mc_visitcell,mc_revisitcell,ieventcounttot) &
     !$OMP REDUCTION(max:mc_revisitcell_max)
     !
     !$ id=OMP_get_thread_num()
     !$ nthreads=OMP_get_num_threads()
     !$ write(stdo,*) 'Thread Nr',id,'of',nthreads,'threads in total'
     !$ iseed=-abs(iseed_start+id)
     !$OMP DO SCHEDULE(dynamic)
     !
     ! Launch all the photons for this wavelength
     !
     do iphot=1,nphot
        if (.not.mc_emergency_break) then
        !
        ! For the photon count routines:
        !
        mc_iphotcurr = iphot
        !
        ! Message
        !
        !$OMP CRITICAL
        cnt   = cnt + 1
        if(mod(cnt,countwrite).eq.0) then
           !$   write(stdo,*) 'Thread:',id,'Photon nr:',cnt
           if(.not.mc_openmp_parallel) write(stdo,*) 'Photon nr ',iphot
           call flush(stdo)
        endif
        !$OMP END CRITICAL
        !
        ! Safety dump? [DISABLED]
        !
!        if((cnt.ge.cntdump).and.(iphot.lt.nphot)) then
!           write(stdo,*) 'WARNING: Wanting to do safety backup, but this is'
!           write(stdo,*) '         not yet implemented.'
!           call flush(stdo)
!           !$omp atomic
!           cnt=0
!        endif
        !
        ! Call the walk routine
        !
        call walk_full_path_scat(params,inu,ierrpriv)
        !       
        ! If ierrpriv.ne.0 then an error occurred, return with non-zero error code
        !
        if(ierrpriv.ne.0) then
           !$OMP CRITICAL
           mc_emergency_break = .true.
           !$OMP END CRITICAL
           write(stdo,*) '!!!!!!!!!!!!!!!!!!! MC_EMERGENCY_BREAK !!!!!!!!!!!!!!!!!!!!!!'
        endif
        !
        ! Write debugging stuff
        !
        !      if(params%debug_write_eventcounts.ne.0) then
        !         write(stdo,*) 'Nr of events = ',ieventcount
        !      endif
     end if
  enddo
  !$OMP END DO 
  !$OMP END PARALLEL
  !
  ! OpenMP Parallellization: destroy locks
  !
  !$ do i=1,size(lock)
  !$    call omp_destroy_lock(lock(i))  
  !$ enddo
  !
  ! OpenMP Parallellization: Diagnostic messages
  !
  !$ write(stdo,*)"Elapsed time:",omp_get_wtime() - seconds;
  !!$ write(stdo,*)"number of conflicts",conflict_counter;
  !!$ open(file='speed_Up_Lock',unit=200,position='append')
  !!$ write(200,*)nthreads,"'",omp_get_wtime() - seconds;
  !!$ close(200)

  if(mc_emergency_break) then
      return
  endif
  enddo
  !
  ! If the optically very thick cells are skipped as sources of photons,
  ! then also make their scattering source functions equal to the planck
  ! function.
  !
!----TO-ADD----
! if(mc_scat_src_maxtau.lt.1d90) then
!    do icell=1,nrcells
!       index = cellindex(icell)
!       if(mc_scat_src_skip(index)) then
!          temp = dusttemp(1,index)  ! Take only first dust temp (== all others)
!          do inu=1,mc_nrfreq
!             mcscat_scatsrc(:,inu,index) = 0.d0
!             do ispec=1,dust_nr_species
!                 mcscat_scatsrc(:,inu,index) = mcscat_scatsrc(:,inu,index) + &
!                    dustdens(ispec,index) * kappa_s(inu,ispec)
!             enddo
!             mcscat_scatsrc(:,inu,index) = mcscat_scatsrc(:,inu,index) * &
!                         bplanck(temp,mc_frequencies(inu))
!          enddo
!       endif
!    enddo 
! endif
!----TO-ADD----
  !
  ! Close path file
  !
  if(params%debug_write_path.eq.1) then
     close(5)
  endif
  !
  ! If statistics writing active, close
  !
  if(params%debug_write_stats.ne.0) then
     close(4)
  endif
  !
  ! Write to the radmc_save.info that we're done
  !
  !open(unit=1,file='radmc_save.info')
  !write(1,*) '-'
  !close(1)
  !
  ! Done...
  !
end subroutine do_monte_carlo_scattering




!--------------------------------------------------------------------------
!         A SIMPLE APPROXIMATION OF SINGLE-LAMBDA SCATTERING
!
! Often it is sufficient to use the single-scattering approximation, if
! the source of radiation is/are only star(s). Here a method is implemented
! that loops over all cells, and calculates the radiation field caused by
! the light of all stars (with their extinction). This is then used to 
! compute the scattering source function. The advantage of this method
! is that it is faster and cleaner (not dependent on stochastics). But this
! does not work for non-point-source sources.
! --------------------------------------------------------------------------
subroutine do_lambda_starlight_single_scattering(params,ierror,scatsrc,meanint)
  implicit none
  type(mc_params) :: params
  integer :: ierror,inu
  doubleprecision :: cellsize,tau,costheta,g,dum
  doubleprecision :: flux,scatsrc0,distance,temp,freq
  logical,optional :: scatsrc,meanint
  logical :: compute_scatsrc,compute_meanint
  integer :: index,idir,ispec,istar,icell,iddr,icount
  doubleprecision :: pos(1:3),cellx0(1:3),cellx1(1:3),src4(1:4)
  !
  ! Some consistency checks
  ! 
  if(star_sphere) then
     write(stdo,*) 'ERROR in single scattering mode: finite-size stars not allowed.'
     stop 4096
  endif
  if(present(scatsrc)) then
     compute_scatsrc=scatsrc 
  else 
     compute_scatsrc=.false.
  endif
  if(present(meanint)) then
     compute_meanint=meanint
  else 
     compute_meanint=.false.
  endif
  if((.not.compute_scatsrc).and.(.not.compute_meanint)) then
     write(stdo,*) 'Single Scattering: Must set scatsrc or meanint.'
     stop 3765
  endif
  if(compute_scatsrc.and.compute_meanint) then
     write(stdo,*) 'Single Scattering: Cannot set BOTH scatsrc and meanint.'
     stop 3766
  endif
  if(nrcells.le.0) then
     write(stdo,*) 'ERROR: Spatial grid is not set.'
     stop
  endif
  if(incl_quantum.ne.0) then
     write(stdo,*) 'ERROR: At the moment quantum heated grains are not included yet'
     stop
  endif
  if((mc_nrfreq.le.0).or.(.not.allocated(mc_frequencies))) then
     write(stdo,*) 'ERROR: Cannot do scattering Monte Carlo without'
     write(stdo,*) '       MC frequency array set.'
     stop
  endif
  !
  ! Initialize the Monte Carlo module by allocating the required arrays
  ! and doing a number of standard checks
  !
  if(compute_scatsrc) then
     call montecarlo_init(params,ierror,101,.false.)
  elseif(compute_meanint) then
     call montecarlo_init(params,ierror,102,.false.)
  else
     stop 6658
  endif
  !
  ! Message
  !
  if(scattering_mode.ne.0) then
     write(stdo,*) 'Using dust scattering mode ',scattering_mode
  endif
  !
  ! Check if source of light is only stars
  ! We ignore the thermal emission of the dust
  !
  if(stellarsrclumtot+emisquanttot+extlumtot+mc_bc_lumall.gt.0.d0) then
     write(stdo,*) 'WARNING: The single-scattering approximation only includes stars as sources of radiation.'
  endif
  !
  ! Get the opacities (assume global opacities for the moment)
  ! ---> COMMENTED OUT 06.03.2017: Opacities are already installed since
  !      montecarlo_init().
  !
  ! temp    = 100.d0     ! Arbitrary temperature...
  ! do inu=1,mc_nrfreq
  !    freq = mc_frequencies(inu)
  !    do ispec=1,dust_nr_species
  !       kappa_a(inu,ispec) = find_dust_kappa_interpol(freq,ispec,temp,1,0,0)
  !       kappa_s(inu,ispec) = find_dust_kappa_interpol(freq,ispec,temp,0,1,0)
  !       kappa_g(inu,ispec) = find_dust_kappa_interpol(freq,ispec,temp,0,0,1)
  !    enddo
  ! enddo
  !
  ! Warning
  !
  if(nrcells.gt.1000000) then
     write(stdo,*) 'This method (single scattering) does accurate but inefficient ray-tracing'
     write(stdo,*) '   from every star to every cell. It may take some time. Progress reports'
     write(stdo,*) '   will be made.'
  endif
  !
  ! Loop over all mc_frequencies
  !
  do inu=1,mc_nrfreq
     !
     ! Message
     !
     write(stdo,*) 'Wavelength nr ',inu,' corresponding to lambda=',  &
               1d4*cc/mc_frequencies(inu),' micron'
     call flush(stdo)
     !
     ! Determine the luminosity of the star(s) in all these wavelength bins
     !
     call montecarlo_compute_freqdep_luminosities(params,mc_frequencies(inu))
     !
     ! Set ray_inu to inu
     !
     ray_inu = inu
     !
     ! Now do the loop over all cells
     !
     icount = 0
     do icell=1,nrcells
        icount = icount + 1
        if(modulo(icount,1000000).eq.0) then
           write(stdo,*) '...done so far: ',icount,' cells of ',nrcells,  &
                ' (=',100*(icount*1.d0)/(nrcells*1.d0),'%)'
        endif
        !
        ! Get index
        ! 
        index  = cellindex(icell)
        !
        ! Get the location of this cell (using cell-center)
        !
        if(amr_tree_present) then
           amrray_cell => amr_index_to_leaf(index)%link
           do iddr=1,3
              cellx0(iddr) = amr_finegrid_xi(amrray_cell%ixyzf(iddr),iddr,amrray_cell%level)
              cellx1(iddr) = amr_finegrid_xi(amrray_cell%ixyzf(iddr)+1,iddr,amrray_cell%level)
           enddo
        else
           call amr_regular_get_ixyz(index,amrray_ix_curr,amrray_iy_curr,amrray_iz_curr)
           cellx0(1) = amr_finegrid_xi(amrray_ix_curr,1,0)
           cellx1(1) = amr_finegrid_xi(amrray_ix_curr+1,1,0)
           cellx0(2) = amr_finegrid_xi(amrray_iy_curr,2,0)
           cellx1(2) = amr_finegrid_xi(amrray_iy_curr+1,2,0)
           cellx0(3) = amr_finegrid_xi(amrray_iz_curr,3,0)
           cellx1(3) = amr_finegrid_xi(amrray_iz_curr+1,3,0)
        endif
        pos(1) = 0.5d0 * ( cellx0(1) + cellx1(1) )
        pos(2) = 0.5d0 * ( cellx0(2) + cellx1(2) )
        pos(3) = 0.5d0 * ( cellx0(3) + cellx1(3) )
        !
        ! Estimate the "cell size" (for smoothing the stellar flux in the
        ! cell where the star is)
        ! 
        cellsize = max(cellx1(1)-cellx0(1),cellx1(2)-cellx0(2),cellx1(3)-cellx0(3))
        !
        ! Convert to other coordinates if necessary
        !
        if((igrid_coord.ge.100).and.(igrid_coord.lt.200)) then
           !
           ! We use spherical coordinates
           !
           ray_cart_x = pos(1) * sin(pos(2)) * cos(pos(3))
           ray_cart_y = pos(1) * sin(pos(2)) * sin(pos(3))
           ray_cart_z = pos(1) * cos(pos(2))
           pos(1) = ray_cart_x
           pos(2) = ray_cart_y
           pos(3) = ray_cart_z
        else
           if((igrid_coord.ge.200).or.(igrid_coord.lt.0)) then
              write(*,*) 'ERROR in do_lambda_starlight_single_scattering(): Coordinate system not supported'
              stop 5239
           endif
        endif
        !
        ! Do a loop over all stars
        !
        do istar=1,nstars
           !
           ! Start at the star position
           !
           ray_cart_x = star_pos(1,istar)
           ray_cart_y = star_pos(2,istar)
           ray_cart_z = star_pos(3,istar)
           !
           ! Find the distance to this star
           !
           distance   = sqrt((ray_cart_x-pos(1))**2+(ray_cart_y-pos(2))**2+(ray_cart_z-pos(3))**2)
           !
           ! Find the direction to move
           !
           ray_cart_dirx = ( pos(1) - ray_cart_x ) / distance
           ray_cart_diry = ( pos(2) - ray_cart_y ) / distance
           ray_cart_dirz = ( pos(3) - ray_cart_z ) / distance
           !
           ! Get the initial cell index (which is the cell index of the star)
           !
           ray_index = star_cellindex(istar)
           !
           ! Now associate pointer to cell.
           ! This is grid type dependent.
           !
           if(igrid_type.lt.100) then
              !
              ! AMR type grid
              !
              if(ray_index.gt.0) then
                 !
                 ! We are in a cell
                 !
                 if(amr_tree_present) then
                    amrray_cell => amr_index_to_leaf(ray_index)%link
                 else
                    call amr_regular_get_ixyz(ray_index,amrray_ix_curr,amrray_iy_curr,amrray_iz_curr)
                 endif
              else
                 nullify(amrray_cell)
                 amrray_ix_curr = -1
                 amrray_iy_curr = -1
                 amrray_iz_curr = -1
              endif
           else
              write(stdo,*) 'ERROR: This grid type not yet implemented'
              stop
           endif
           !
           ! Now start the cell walking from the star to the current cell,
           ! and calculate the optical depth between the star and the current cell.
           !
           call walk_cells_optical_depth(params,distance,inu,tau)
           !
           ! Smooth the distance, if the cell is too close to a star or the 
           ! star is inside this cell, to not over-estimate the 1/r^2.
           !
           distance   = max(distance,cellsize)
           !
           ! Compute the stellar flux at this position
           !
           flux = exp(-tau) * star_lum(istar) / ( fourpi * distance**2 )
           !
           ! If we have scattering mode 2 or 3 (anisotropic, but no polarization)
           ! then we can pre-compute the phase functions for all dust species.
           !
           ! Note: if we use local observer, then we cannot pre-compute the 
           !       phase function for non-isotropic scattering.
           !
           if(compute_scatsrc) then
              if(scattering_mode.eq.2) then
                 !
                 ! Anisotropic scattering using the Henyey-Greenstein function
                 !
                 if(.not.mcscat_localobserver) then
                    do idir=1,mcscat_nrdirs
                       costheta = ray_cart_dirx*mcscat_dirs(1,idir) +              &
                                  ray_cart_diry*mcscat_dirs(2,idir) +              &
                                  ray_cart_dirz*mcscat_dirs(3,idir)
                       do ispec=1,dust_nr_species
                          g                            = kappa_g(ray_inu,ispec)
                          mcscat_phasefunc(idir,ispec) = henyeygreenstein_phasefunc(g,costheta)
                       enddo
                    enddo
                 else
                    write(stdo,*) 'ERROR: For now, anisotropic scattering not allowed'
                    write(stdo,*) '       in local observer mode.'
                    stop 453
                 endif
              elseif(scattering_mode.eq.3) then
                 !
                 ! Anisotropic scattering using tabulated phase function (the
                 ! Z11 matrix elements)
                 !
                 if(.not.mcscat_localobserver) then
                    do idir=1,mcscat_nrdirs
                       costheta = ray_cart_dirx*mcscat_dirs(1,idir) +              &
                                  ray_cart_diry*mcscat_dirs(2,idir) +              &
                                  ray_cart_dirz*mcscat_dirs(3,idir)
                       do ispec=1,dust_nr_species
                          mcscat_phasefunc(idir,ispec) = anisoscat_phasefunc(ray_inu,ispec,costheta)
                       enddo
                    enddo
                 else
                    write(stdo,*) 'ERROR: For now, anisotropic scattering not allowed'
                    write(stdo,*) '       in local observer mode.'
                    stop 453
                 endif
              endif
              !
              ! Check
              !
              if(.not.allocated(mcscat_scatsrc_iquv)) then
                 write(stdo,*) 'ERROR: Strange thing in single scattering. Warn author.'
                 stop 7449
              endif
              if(index.ne.ray_index) then
                 write(stdo,*) 'INTERNAL ERROR. Warn author.'
                 stop 2095
              endif
              !
              ! Now add to the scattering source function
              !
              do idir=1,mcscat_nrdirs
                 !
                 ! Compute the scattering source function for the given 
                 ! scattering mode
                 !
                 if(scattering_mode.eq.1) then
                    !
                    ! Isotropic scattering: simply add the source term
                    !
                    dum = 0.d0
                    do ispec=1,dust_nr_species
                       dum = dum + dustdens(ispec,index) * kappa_s(inu,ispec)
                    enddo
                    scatsrc0 = dum * flux / fourpi
                    mcscat_scatsrc_iquv(ray_inu,ray_index,1,1) =                   &
                         mcscat_scatsrc_iquv(ray_inu,ray_index,1,1) + scatsrc0
                 elseif((scattering_mode.eq.2).or.(scattering_mode.eq.3)) then
                    !
                    ! Anisotropic scattering: add only for the given directions
                    !
                    ! BUGFIX 2018.06.06
                    ! The phase function must be done for each species separately
                    !
                    do ispec=1,dust_nr_species
                       mcscat_scatsrc_iquv(ray_inu,ray_index,1,idir) =             &
                            mcscat_scatsrc_iquv(ray_inu,ray_index,1,idir) +        &
                            mcscat_phasefunc(idir,ispec) * flux *                  &
                            dustdens(ispec,index) * kappa_s(inu,ispec) / fourpi
                    enddo
                    if(mcscat_localobserver) then
                       write(stdo,*) 'ABORTING: For now, non-isotropic scattering and '
                       write(stdo,*) '          local observer are not yet compatible.'
                       stop
                    endif
                 elseif((scattering_mode.eq.4).or.(scattering_mode.eq.5)) then
                    !
                    ! Polarized scattering, full phase function
                    !
                    if(mcscat_localobserver) then
                       !
                       ! Polarized scattering and local-observer are mutually
                       ! incompatible. NOTE: Also multiple directions (i.e.
                       ! mcscat_nrdirs.gt.1) is not possible.
                       !
                       write(stdo,*) 'ERROR: Polarized scattering and local-observer are incompatible.'
                       stop
                    endif
                    !
                    ! Make the "phot" array for use with the polarization
                    ! module. We choose an arbitrary S-vector, because the
                    ! incoming light is assumed to be unpolarized.
                    !
                    photpkg%E    = flux
                    photpkg%Q    = 0.d0
                    photpkg%U    = 0.d0
                    photpkg%V    = 0.d0
                    photpkg%n(1) = ray_cart_dirx
                    photpkg%n(2) = ray_cart_diry
                    photpkg%n(3) = ray_cart_dirz
                    call polarization_make_s_vector(photpkg%n,photpkg%s)   ! arb S-vector
                    !
                    ! Now add the scattering contribution of each of
                    ! the dust species. Note that we do not
                    ! need to divide by 4*pi here because we use 
                    ! Z instead of kappa_scat.
                    !
                    do ispec=1,dust_nr_species
                       call polarization_randomorient_scatsource(photpkg,  &
                            mcscat_dirs(:,1),mcscat_svec(:,1),             &
                            scat_munr,scat_mui_grid(:),                    &
                            zmatrix(:,ray_inu,1,ispec),                    &
                            zmatrix(:,ray_inu,2,ispec),                    &
                            zmatrix(:,ray_inu,3,ispec),                    &
                            zmatrix(:,ray_inu,4,ispec),                    &
                            zmatrix(:,ray_inu,5,ispec),                    &
                            zmatrix(:,ray_inu,6,ispec),                    &
                            src4)
                       mcscat_scatsrc_iquv(ray_inu,ray_index,1:4,1) =      &
                            mcscat_scatsrc_iquv(ray_inu,ray_index,1:4,1) + &
                            dustdens(ispec,index) * src4(1:4)
                    enddo
                 else
                    stop 2960
                 endif
              enddo
           elseif(compute_meanint) then
              !
              ! Compute the mean intensity: despite the word "single scattering"
              ! this does not include the single scattering effect in the 
              ! mean intensity! It only includes the extinction of the 
              ! starlight. So it is a bit trivial...
              !
              ! BUGFIX: 2018.06.01: icell --> ray_index
              !
              mcscat_meanint(ray_inu,ray_index) = &
                   mcscat_meanint(ray_inu,ray_index) + flux / fourpi
           else
              write(stdo,*) 'Strange error in single scattering mode. Contact Author.'
              stop 5491
           endif
        enddo
     enddo
  enddo
  !
  ! Close path file
  !
  if(params%debug_write_path.eq.1) then
     close(5)
  endif
  !
  ! If statistics writing active, close
  !
  if(params%debug_write_stats.ne.0) then
     close(4)
  endif
  !
  ! Done...
  !
end subroutine do_lambda_starlight_single_scattering


!--------------------------------------------------------------------------
!         AN EVEN SIMPLER APPROXIMATION OF SINGLE-LAMBDA SCATTERING
!
! As do_lambda_starlight_single_scattering(), but now for the simplest
! possible case: spherical coordinate, 1 pointlike star in the center,
! single scattering.
! --------------------------------------------------------------------------
subroutine do_lambda_starlight_single_scattering_simple(params,ierror,scatsrc,meanint)
  implicit none
  type(mc_params) :: params
  integer :: ierror,inu
  doubleprecision :: dtau,dvol,g,luminosity,albedo,ds,costheta
  doubleprecision :: scatsrc0,alpha_tot,xp1,flux,r,theta,phi
  doubleprecision :: cosdphi,sindphi,deltaphi,xbk,ybk,levent,cosphievent,sinphievent
  logical,optional :: scatsrc,meanint
  logical :: compute_scatsrc,compute_meanint
  integer :: idir,ispec,istar,icell,iddr
  integer :: ir,itheta,iphi
  doubleprecision :: pos(1:3),cellx0(1:3),cellx1(1:3),src4(1:4)
  !
  ! Some consistency checks
  ! 
  if(star_sphere) then
     write(stdo,*) 'ERROR in single scattering mode: finite-size stars not allowed.'
     stop 4096
  endif
  if(present(scatsrc)) then
     compute_scatsrc=scatsrc 
  else 
     compute_scatsrc=.false.
  endif
  if(present(meanint)) then
     compute_meanint=meanint
  else 
     compute_meanint=.false.
  endif
  if((.not.compute_scatsrc).and.(.not.compute_meanint)) then
     write(stdo,*) 'Single Scattering: Must set scatsrc or meanint.'
     stop 3765
  endif
  if(compute_scatsrc.and.compute_meanint) then
     write(stdo,*) 'Single Scattering: Cannot set BOTH scatsrc and meanint.'
     stop 3766
  endif
  if(nrcells.le.0) then
     write(stdo,*) 'ERROR: Spatial grid is not set.'
     stop
  endif
  if(incl_quantum.ne.0) then
     write(stdo,*) 'ERROR: At the moment quantum heated grains are not included yet'
     stop
  endif
  if((mc_nrfreq.le.0).or.(.not.allocated(mc_frequencies))) then
     write(stdo,*) 'ERROR: Cannot do scattering Monte Carlo without'
     write(stdo,*) '       MC frequency array set.'
     stop
  endif
  if((igrid_coord.lt.100).or.(igrid_coord.ge.200)) then
     write(stdo,*) 'ERROR: The simple single scattering mode requires spherical coordinates'
     stop
  endif
  if(amr_style.ne.0) then
     write(stdo,*) 'ERROR: The simple single scattering mode cannot use AMR. Use simple grid instead.'
     stop
  endif
  if(nstars.ne.1) then
     write(stdo,*) 'ERROR: The simple single scattering mode has to have exactly 1 star.'
     stop
  endif
  if((incl_stellarsrc.ne.0).or.(incl_extlum.ne.0).or. &
     (incl_heatsource.ne.0).or.(incl_thermbc.ne.0)) then
     write(stdo,*) 'ERROR: The simple single scattering mode can only handle a single energy source: the central star'
     stop
  endif
  !
  ! Initialize the Monte Carlo module by allocating the required arrays
  ! and doing a number of standard checks
  !
  if(compute_scatsrc) then
     call montecarlo_init(params,ierror,101,.false.)
  elseif(compute_meanint) then
     call montecarlo_init(params,ierror,102,.false.)
  else
     stop 6658
  endif
  !
  ! Message
  !
  if(scattering_mode.ne.0) then
     write(stdo,*) 'Using dust scattering mode ',scattering_mode
  endif
  !
  ! Check if source of light is only stars
  ! We ignore the thermal emission of the dust
  !
  if(stellarsrclumtot+emisquanttot+extlumtot+mc_bc_lumall.gt.0.d0) then
     write(stdo,*) 'WARNING: The single-scattering approximation only includes stars as sources of radiation.'
  endif
  !
  ! Get the opacities (assume global opacities for the moment)
  ! ---> COMMENTED OUT 06.03.2017: Opacities are already installed since
  !      montecarlo_init().
  !
  ! temp    = 100.d0     ! Arbitrary temperature...
  ! do inu=1,mc_nrfreq
  !    freq = mc_frequencies(inu)
  !    do ispec=1,dust_nr_species
  !       kappa_a(inu,ispec) = find_dust_kappa_interpol(freq,ispec,temp,1,0,0)
  !       kappa_s(inu,ispec) = find_dust_kappa_interpol(freq,ispec,temp,0,1,0)
  !       kappa_g(inu,ispec) = find_dust_kappa_interpol(freq,ispec,temp,0,0,1)
  !    enddo
  ! enddo
  !
  ! Loop over all mc_frequencies
  !
  do inu=1,mc_nrfreq
     !
     ! Message
     !
     write(stdo,*) 'Wavelength nr ',inu,' corresponding to lambda=',  &
               1d4*cc/mc_frequencies(inu),' micron'
     call flush(stdo)
     !
     ! Determine the luminosity of the star(s) in all these wavelength bins
     !
     call montecarlo_compute_freqdep_luminosities(params,mc_frequencies(inu))
     !
     ! Set ray_inu to inu
     !
     ray_inu = inu
     !
     ! Now do the loop over theta and phi coordinates
     !
     do itheta=1,amr_grid_ny
        do iphi=1,amr_grid_nz
           !
           ! Determine the direction vector, passing through the center of
           ! the theta-phi cell square
           !
           theta         = 0.5*(amr_grid_xi(itheta+1,2)+amr_grid_xi(itheta,2))
           phi           = 0.5*(amr_grid_xi(iphi+1,3)+amr_grid_xi(iphi,3))
           ray_cart_dirx = sin(theta)*cos(phi)
           ray_cart_diry = sin(theta)*sin(phi)
           ray_cart_dirz = cos(theta)
           !
           ! Initialize the ray-tracing at the star
           !
           luminosity = star_lum(1)
           !
           ! Now trace a ray radially from the start outside
           !
           do ir=1,amr_grid_nx
              !
              ! In which cell are we?
              !
              ray_index = ir+(itheta-1)*amr_grid_nx+(iphi-1)*amr_grid_nx*amr_grid_ny
              !
              ! Compute the ds in radial direction
              !
              ds = amr_grid_xi(ir+1,1) - amr_grid_xi(ir,1)
              !
              ! Determine location
              !
              r          = 0.5*(amr_grid_xi(ir+1,1)+amr_grid_xi(ir,1))
              ray_cart_x = r * ray_cart_dirx
              ray_cart_y = r * ray_cart_diry
              ray_cart_z = r * ray_cart_dirz
              !
              ! Compute the extinction within the cell
              !
              alpha_a_tot = 0.d0
              alpha_s_tot = 0.d0
              do ispec=1,dust_nr_species
                 alpha_a(ispec) = dustdens(ispec,ray_index) * &
                                  kappa_a(ray_inu,ispec) + 1d-99
                 alpha_s(ispec) = dustdens(ispec,ray_index) * &
                                  kappa_s(ray_inu,ispec) + 1d-99
                 alpha_a_tot    = alpha_a_tot + alpha_a(ispec)
                 alpha_s_tot    = alpha_s_tot + alpha_s(ispec)
              enddo
              alpha_tot = alpha_a_tot + alpha_s_tot
              dtau      = alpha_tot * ds
              albedo    = alpha_s_tot / alpha_tot
              !
              ! Compute the "volume" (4*pi/3)*(r_{i+1/2}^3-r_{i-1/2}^3)
              !
              dvol = (fourpi/3.) * ( amr_grid_xi(ir+1,1)**3 - amr_grid_xi(ir,1)**3 )
              !
              ! Compute xp1 = 1-exp(-tau)
              !
              if(dtau.lt.1e-6) then
                 xp1 = dtau
              else
                 xp1 = 1.d0 - exp(-dtau)
              endif
              !
              !!!!!!!!!!!!!!!!!!!!!jscat = albedo * luminosity * xp1 / ( fourpi * dvol )
              !
              ! Compute, for the anisotropic scattering, the flux
              !
              flux  = luminosity * xp1 / ( dvol * alpha_tot )
              !
              ! Update the luminosity 
              !
              luminosity = luminosity * exp(-dtau)
              !
              ! Now proceed depending on what we wish to calculate (scat source
              ! or mean int)
              !
              if(compute_scatsrc) then
                 !
                 ! Prepare stuff in case of anisotropic scattering
                 !
                 if(scattering_mode.eq.2) then
                    !
                    ! Anisotropic scattering using the Henyey-Greenstein function
                    !
                    if(.not.mcscat_localobserver) then
                       do idir=1,mcscat_nrdirs
                          costheta = ray_cart_dirx*mcscat_dirs(1,idir) +              &
                                     ray_cart_diry*mcscat_dirs(2,idir) +              &
                                     ray_cart_dirz*mcscat_dirs(3,idir)
                          do ispec=1,dust_nr_species
                             g                            = kappa_g(ray_inu,ispec)
                             mcscat_phasefunc(idir,ispec) = henyeygreenstein_phasefunc(g,costheta)
                          enddo
                       enddo
                    else
                       write(stdo,*) 'ERROR: For now, anisotropic scattering not allowed'
                       write(stdo,*) '       in local observer mode.'
                       stop 453
                    endif
                 elseif(scattering_mode.eq.3) then
                    !
                    ! Anisotropic scattering using tabulated phase function (the
                    ! Z11 matrix elements)
                    !
                    if(.not.mcscat_localobserver) then
                       do idir=1,mcscat_nrdirs
                          costheta = ray_cart_dirx*mcscat_dirs(1,idir) +              &
                                     ray_cart_diry*mcscat_dirs(2,idir) +              &
                                     ray_cart_dirz*mcscat_dirs(3,idir)
                          do ispec=1,dust_nr_species
                             mcscat_phasefunc(idir,ispec) = anisoscat_phasefunc(ray_inu,ispec,costheta)
                          enddo
                       enddo
                    else
                       write(stdo,*) 'ERROR: For now, anisotropic scattering not allowed'
                       write(stdo,*) '       in local observer mode.'
                       stop 453
                    endif
                 endif
                 !
                 ! Check
                 !
                 if(.not.allocated(mcscat_scatsrc_iquv)) then
                    write(stdo,*) 'ERROR: Strange thing in single scattering. Warn author.'
                    stop 7449
                 endif
                 !
                 ! Now add to the scattering source function
                 !
                 ! Compute the scattering source function for the given 
                 ! scattering mode
                 !
                 if(scattering_mode.eq.1) then
                    !
                    ! Isotropic scattering: simply add the source term
                    !
                    scatsrc0 = flux * albedo * alpha_tot / fourpi
                    mcscat_scatsrc_iquv(ray_inu,ray_index,1,1) =                   &
                         mcscat_scatsrc_iquv(ray_inu,ray_index,1,1) + scatsrc0
                 elseif((scattering_mode.eq.2).or.(scattering_mode.eq.3)) then
                    !
                    ! Anisotropic scattering: add only for the given directions
                    !
                    if(mcscat_localobserver) then
                       write(stdo,*) 'ABORTING: For now, non-isotropic scattering and '
                       write(stdo,*) '          local observer are not yet compatible.'
                       stop
                    endif
                    do idir=1,mcscat_nrdirs
                       do ispec=1,dust_nr_species
                          mcscat_scatsrc_iquv(ray_inu,ray_index,1,idir) =             &
                               mcscat_scatsrc_iquv(ray_inu,ray_index,1,idir) +        &
                               mcscat_phasefunc(idir,ispec) * flux *                  &
                               dustdens(ispec,ray_index) * kappa_s(inu,ispec) / fourpi
                       enddo
                    enddo
                 elseif((scattering_mode.eq.4).or.(scattering_mode.eq.5)) then
                    !
                    ! Polarized scattering, full phase function
                    !
                    if(mcscat_localobserver) then
                       !
                       ! Polarized scattering and local-observer are mutually
                       ! incompatible. NOTE: Also multiple directions (i.e.
                       ! mcscat_nrdirs.gt.1) is not possible.
                       !
                       write(stdo,*) 'ERROR: Polarized scattering and local-observer are incompatible.'
                       stop
                    endif
                    !
                    ! Make the "phot" array for use with the polarization
                    ! module. We choose an arbitrary S-vector, because the
                    ! incoming light is assumed to be unpolarized.
                    !
                    photpkg%E    = flux
                    photpkg%Q    = 0.d0
                    photpkg%U    = 0.d0
                    photpkg%V    = 0.d0
                    photpkg%n(1) = ray_cart_dirx
                    photpkg%n(2) = ray_cart_diry
                    photpkg%n(3) = ray_cart_dirz
                    call polarization_make_s_vector(photpkg%n,photpkg%s)   ! arb S-vector
                    !
                    ! If 2-D axisymmetric anisotropic scattering mode,
                    ! then we must rotate the photpkg%n and photpkg%s
                    ! vectors appropriately
                    !
                    if(dust_2daniso) then
                       !
                       ! Safety check. Note: the last phi point is
                       ! 360 degrees == first one. This makes it easier
                       ! lateron to interpolate, even though it costs
                       ! a tiny bit more memory and is a tiny bit slower.
                       !
                       if(mcscat_nrdirs.ne.dust_2daniso_nphi+1) stop 2067
                       !
                       ! Rotate n and s vector such that the event lies
                       ! in the x-z-plane (positive x)
                       !
                       levent      = sqrt(ray_cart_x**2+ray_cart_y**2)
                       if(levent.gt.0.d0) then
                          cosphievent  = ray_cart_x/levent
                          sinphievent  = ray_cart_y/levent
                          xbk          = photpkg%n(1)
                          ybk          = photpkg%n(2)
                          photpkg%n(1) =  cosphievent * xbk + sinphievent * ybk
                          photpkg%n(2) = -sinphievent * xbk + cosphievent * ybk
                          xbk          = photpkg%s(1)
                          ybk          = photpkg%s(2)
                          photpkg%s(1) =  cosphievent * xbk + sinphievent * ybk
                          photpkg%s(2) = -sinphievent * xbk + cosphievent * ybk
                       endif
                       !
                       ! Pre-compute the cos and sin of the small incremental
                       ! rotations
                       !
                       deltaphi = twopi / dust_2daniso_nphi
                       cosdphi  = cos(deltaphi)
                       sindphi  = sin(deltaphi)
                    endif
                    !
                    ! Now add the scattering contribution of each of
                    ! the dust species. Note that we do not
                    ! need to divide by 4*pi here because we use 
                    ! Z instead of kappa_scat.
                    !
                    do idir=1,mcscat_nrdirs
                       do ispec=1,dust_nr_species
                          call polarization_randomorient_scatsource(photpkg,  &
                               mcscat_dirs(:,idir),mcscat_svec(:,idir),       &
                               scat_munr,scat_mui_grid(:),                    &
                               zmatrix(:,ray_inu,1,ispec),                    &
                               zmatrix(:,ray_inu,2,ispec),                    &
                               zmatrix(:,ray_inu,3,ispec),                    &
                               zmatrix(:,ray_inu,4,ispec),                    &
                               zmatrix(:,ray_inu,5,ispec),                    &
                               zmatrix(:,ray_inu,6,ispec),                    &
                               src4)
                          mcscat_scatsrc_iquv(ray_inu,ray_index,1:4,idir) =       &
                               mcscat_scatsrc_iquv(ray_inu,ray_index,1:4,idir)  + &
                               dustdens(ispec,ray_index) * src4(1:4)
                       enddo
                       !
                       ! If 2-D axisymmetric anisotropic scattering mode,
                       ! then we must rotate the photpkg%n and photpkg%s
                       ! vectors appropriately
                       !
                       if(dust_2daniso) then
                          xbk = photpkg%n(1)
                          ybk = photpkg%n(2)
                          photpkg%n(1) =  cosdphi * xbk - sindphi * ybk
                          photpkg%n(2) =  sindphi * xbk + cosdphi * ybk
                          xbk = photpkg%s(1)
                          ybk = photpkg%s(2)
                          photpkg%s(1) =  cosdphi * xbk - sindphi * ybk
                          photpkg%s(2) =  sindphi * xbk + cosdphi * ybk
                       endif
                    enddo
                 else
                    stop 2960
                 endif
              elseif(compute_meanint) then
                 !
                 ! Compute the mean intensity: despite the word "single scattering"
                 ! this does not include the single scattering effect in the 
                 ! mean intensity! It only includes the extinction of the 
                 ! starlight. So it is a bit trivial...
                 !
                 mcscat_meanint(ray_inu,ray_index) = &
                      mcscat_meanint(ray_inu,ray_index) + flux / fourpi
              else
                 write(stdo,*) 'Strange error in single scattering simple mode. Contact Author.'
                 stop 5491
              endif
           enddo
        enddo
     enddo
  enddo
  !
  ! Close path file
  !
  if(params%debug_write_path.eq.1) then
     close(5)
  endif
  !
  ! If statistics writing active, close
  !
  if(params%debug_write_stats.ne.0) then
     close(4)
  endif
  !
  ! Done...
  !
end subroutine do_lambda_starlight_single_scattering_simple



!--------------------------------------------------------------------------
!         WALK THE ENTIRE RANDOM WALK UNTIL PHOTON LEAVES SYSTEM
!
!     This routine does the scatterings while keeping track of 
!     the energy lost by a package. It only stops after leaving 
!     the grid.
!--------------------------------------------------------------------------
subroutine walk_full_path_bjorkmanwood(params,ierror)
  implicit none
  type(mc_params) :: params
  integer :: ierror
  !
  double precision :: taupath,sint,epsdist
  double precision :: dir(1:3),enold,albedo,dum,g,rn
  double precision :: pdirx,pdiry,pdirz,dirnewx,dirnewy,dirnewz
  double precision :: dir_perp,dir_planex,dir_planey,dummy,size
  double precision :: mrw_dusttemp,mrw_energy,enthres,enerphot
  double precision :: pos(1:3),nrevents,cellx0(1:3),cellx1(1:3)
  integer :: idir_cross,ilr_cross
  integer :: ispec,iqactive,istar,icell,itemplate
  integer :: ibnd,bc_idir,bc_ilr,ix,iy,iz,illum
  integer :: index_prev,count_samecell,iddr,idim,icoord
  logical :: ok,arrived,therm,usesphere,nospheres,incell
  type(amr_branch), pointer :: acell
  !$ logical::continue
  !
  ! Do checks
  !
  if(debug_check_all.eq.1) then
     if((incl_quantum.eq.0).and.(emisquanttot.ne.0.d0)) then
        write(stdo,*) 'ERROR: Emisquant.ne.0 while quantum not switched on...'
        stop
     endif
     if((nstars.eq.0).and.(starlumtot.ne.0.d0)) then
        write(stdo,*) 'ERROR: Starlumtot.ne.0 while nstars==0...'
        stop
     endif
  endif
  !
  ! Reset
  !
  iqactive = 0
  mc_photon_destroyed = .false.
  !
  ! First choose whether to launch the photon from (one of) the star(s) or
  ! from the external radiation field or from inside the grid (due to the
  ! internally released energy) or from the locally emitted spectral
  ! source.. We make this choice on the basis of the total luminosity of all
  ! sources.
  !
  if(mc_cumlum5.gt.0.d0) then
     !
     ! There are more kinds of luminosity sources than just the star(s)
     !
     rn = ran2(iseed)
     !
  else
     !
     ! Only one source of luminosity: the star(s)
     ! By choosing rn=2 we avoid a costly random number generation and 
     ! we immediately go to the stellar emission
     !
     rn        = 2.d0
     !
  endif
  if(rn.ge.mc_cumlum5) then
     !
     ! ----- We have a photon emitted by a star (or, for 1dpp, an illum source) -----
     !
     if(nstars.ge.1) then
        !
        ! We have a star or multiple stars
        !
        ! Switch on the quantum flag if quantum heating is on.
        !
        ! This is because in the quick quantum mode we treat stellar
        ! photons differently from dust emission photons. This quick 
        ! method does not include quantum heating by thermal emission.
        ! For PAHs, however, this is rarely the case, so this approximation
        ! is likely to be fine.
        !
        if(incl_quantum.ne.0) then
           if((incl_quantum.lt.0).and.(incl_quantum.gt.-10)) then
              iqactive = 1
           elseif(incl_quantum.eq.-10) then
              iqactive = -1
           elseif(incl_quantum.eq.-11) then
              iqactive = 2
           else
              iqactive = 0
           endif
        else
           iqactive = 0
        endif
        !
        ! Which star?
        !
        if(nstars.gt.1) then
           rn = ran2(iseed)
           call hunt(star_lumcum,nstars+1,rn,istar)
           if(istar.gt.nstars) stop
           if(istar.lt.1) stop
        else
           istar = 1
           if(nstars.eq.0) then
              write(stdo,*) 'INTERNAL ERROR with 0 stars in Monte Carlo module'
              stop
           endif
        endif
        !
        ! If we have the spherical star mode switched on, then 
        ! switch it on here as well
        !
        if(star_sphere) then
           usesphere = .true.
        else
           usesphere = .false.
        endif
        !
        ! If we have the weighted photons mode, we must do some
        ! extra stuff
        !
        if(mc_weighted_photons) then
           !
           ! Get the current energy of this photon package
           !
           energy = mc_energy_stars(istar)
           !
           ! If we have weighted photons mode, and a star is outside
           ! of the grid (more precise: outside of the containing sphere), 
           ! then we (at least for now) do not allow the treatment of the
           ! star as a sphere. NOTE: If the star lies inside the inner cavity
           ! of spherical coordinates, it is inside of the containing sphere,
           ! and it will be still treated as a sphere if that mode is on. This
           ! switching off of the sphere mode is only if the star is far away.
           !
           if(star_fraclum(istar).ne.1.d0) then
              usesphere = .false.
           endif
        endif
        !
        ! First choose a random theta at which the photon from the
        ! (spherically symmetric) star enters the inner edge of the
        ! coordinate system.
        !
        if(usesphere) then
           !
           ! Treat star as a sphere 
           !
           ! NOTE: With a finite stellar size, we also have to absorb photons
           !       that hit the star. 
           !
           ! NOTE: For now finite size stars are only allowed to be OUTSIDE
           !       the grid. Otherwise it becomes a nightmare.
           !
           ! Find a random point on the sphere
           !
           call montecarlo_randomdir(pdirx,pdiry,pdirz)
           ray_cart_x  = star_pos(1,istar) + pdirx*star_r(istar)*1.000000000001d0
           ray_cart_y  = star_pos(2,istar) + pdiry*star_r(istar)*1.000000000001d0
           ray_cart_z  = star_pos(3,istar) + pdirz*star_r(istar)*1.000000000001d0
           !
           ! If mirror symmetry, then enforce this here
           !
           if(igrid_mirror.eq.1) then
              if(ray_cart_z.lt.0.d0) then
                 ray_cart_z = -ray_cart_z
              endif
           endif
           !
           ! NOTE: For now finite size stars are only allowed to be OUTSIDE
           !       the grid. Otherwise it becomes a nightmare.
           !
           ray_index = 0
           !
           ! Now determine a random direction for the photon. One must,
           ! however, be careful here because the direction is NOT simply a 3-D
           ! random direction with all inward directions rejected.  The thing
           ! is this: the intensity of a blackbody surface is the same at all
           ! inclinations. However, if one would simply use a normal random
           ! direction (using montecarlo_randomdir()) then one would see a
           ! brighter intensity as the surface is viewed under a more grazing
           ! inclination. This is incorrect.  So to compensate for this one
           ! must have the change of a photon being emitted in direction
           ! mu=cos(theta) to be proportional to mu. This is Lambertian emission.
           !
           ! NOTE: I have been confused again, here, with this sqrt(rn). But it
           ! is indeed correct! The example that made me understand again why
           ! this is so is the example of a transport company with fast trucks,
           ! sending out per day 10 trucks, and a firm with slow trucks, also
           ! sending out per day 10 trucks. Both firms have the same flux of
           ! goods, even though firm2 has slow trucks. Firm2 just requires more
           ! trucks to keep the flux going. That is why photons with mu\simeq 0
           ! (slow) would carry the same flux as photons with mu\simeq 1 (fast)
           ! if they would be emitted equally often. Since we know that
           ! mu\simeq 0 photons carry less flux, they must be emitted less
           ! often. Hence P(mu)=mu, i.e. Pcum(mu)=mu^2, or so to say:
           ! mu=sqrt(Pcum)=sqrt(rn).
           ! 
           rn            = ran2(iseed)
           ray_cart_dirx = sqrt(rn)
           dum           = twopi*ran2(iseed)
           rn            = sqrt(1.d0-rn)
           ray_cart_diry = rn*cos(dum)
           ray_cart_dirz = rn*sin(dum)
           call montecarlo_rotatevec(ray_cart_dirx,ray_cart_diry,ray_cart_dirz, &
                                     pdirx,pdiry,pdirz)
           !
        else
           !
           ! The star is assumed to be a point source
           !
           ray_cart_x  = star_pos(1,istar)
           ray_cart_y  = star_pos(2,istar)
           ray_cart_z  = star_pos(3,istar)
           !
           ! If mirror symmetry, then enforce this here
           ! NOTE: This should not be necessary: It should have already 
           !       been taken care of. But for now let's keep it.
           !
           if(igrid_mirror.eq.1) then
!              if(ray_cart_z.lt.0.d0) then
!                 ray_cart_z = -ray_cart_z
!              endif
              if(ray_cart_z.ne.0.d0) then
                 write(stdo,*) 'ERROR: With mirror symmetry, stars must lie on z=0 plane.'
                 stop 732
              endif
           endif
           !
           ! Find random direction
           !
           if(star_fraclum(istar).eq.1.d0) then
              !
              ! Pure random direction in 4*pi
              !
              call montecarlo_randomdir(ray_cart_dirx,ray_cart_diry,  &
                                        ray_cart_dirz)
           else
              !
              ! Focused toward the model grid (the star lies outside the grid)
              !
              rn            = ran2(iseed)
              dum           = star_fraclum(istar)*rn
              ray_cart_dirx = 1.d0-2.d0*dum
              sint          = 2.d0*sqrt(dum-dum**2)
              dum           = twopi*ran2(iseed)
              ray_cart_diry = sint*cos(dum)
              ray_cart_dirz = sint*sin(dum)
              pdirx         = grid_contsph_x-star_pos(1,istar)
              pdiry         = grid_contsph_y-star_pos(2,istar)
              pdirz         = grid_contsph_z-star_pos(3,istar)
              dum           = 1.d0 / sqrt(pdirx**2+pdiry**2+pdirz**2)
              pdirx         = pdirx * dum
              pdiry         = pdiry * dum
              pdirz         = pdirz * dum
              call montecarlo_rotatevec(ray_cart_dirx,ray_cart_diry,ray_cart_dirz, &
                                        pdirx,pdiry,pdirz)
           endif
           !
           ! If we have mirror symmetry, then all stars MUST lie on the
           ! equatorial plane and their photons must point upward.
           !
           if(igrid_mirror.eq.1) then
              ray_cart_dirz = abs(ray_cart_dirz)
           endif
           !
           ! Which cell is this star in, if it is in the grid?
           !
           ray_index = star_cellindex(istar)
           !
        endif
        !
        ! Find a random ray_inu for the photon package, based on the stellar
        ! flux
        !
        rn = ran2(iseed)
        call hunt(star_spec_cum(:,istar),freq_nr+1,rn,ray_inu)
        if(debug_check_all.eq.1) then
           if((ray_inu.lt.1).or.(ray_inu.gt.freq_nr)) then
              write(stdo,*) 'INTERNAL ERROR: ray_inu out of range for star'
              stop
           endif
        endif
     else
        !
        ! It is an illumination beam (only for 1-D plane-parallel)
        !
        ! Switch on the quantum flag if quantum heating is on.
        !
        ! This is because in the quick quantum mode we treat stellar
        ! photons differently from dust emission photons. This quick 
        ! method does not include quantum heating by thermal emission.
        ! For PAHs, however, this is rarely the case, so this approximation
        ! is likely to be fine.
        !
        if(incl_quantum.ne.0) then
           if((incl_quantum.lt.0).and.(incl_quantum.gt.-10)) then
              iqactive = 1
           elseif(incl_quantum.eq.-10) then
              iqactive = -1
           elseif(incl_quantum.eq.-11) then
              iqactive = 2
           else
              iqactive = 0
           endif
        else
           iqactive = 0
        endif
        !
        ! Which illumination source?
        !
        if(nillum.gt.1) then
           rn = ran2(iseed)
           call hunt(illum_fluxtot_cum,nillum+1,rn,illum)
           if(illum.gt.nillum) stop
           if(illum.lt.1) stop
        else
           illum = 1
           if(nillum.eq.0) then
              write(stdo,*) 'INTERNAL ERROR with 0 illumination beams in Monte Carlo module'
              stop
           endif
        endif
        !
        ! If we have the weighted photons mode, we must do some
        ! extra stuff
        !
        if(mc_weighted_photons) then
           !
           ! Get the current energy of this photon package
           !
           energy = mc_energy_illum(illum)
        endif
        !
        ! Now put the photon package well above the grid to start
        !
        ray_cart_x  = 0.d0
        ray_cart_y  = 0.d0
        if(illum_costheta(illum).gt.0.d0) then
           ray_cart_z  = amr_grid_xi(amr_grid_nz+1,3)+abs(amr_grid_xi(amr_grid_nz+1,3)-amr_grid_xi(1,3))
        else
           ray_cart_z  = amr_grid_xi(1,3)-abs(amr_grid_xi(amr_grid_nz+1,3)-amr_grid_xi(1,3))
        endif
        !
        ! Set the direction. The direction is defined such that
        ! if theta = theta_observer and phi = phi_observer, the
        ! illuminating source is the observer. 
        !
        ray_cart_dirz = -illum_costheta(illum)
        dummy         = sqrt(1.d0-illum_costheta(illum)**2)
        ray_cart_diry = dummy * cos(illum_phi(illum))
        ray_cart_dirx = dummy * sin(illum_phi(illum))
        !
        ! Find a random ray_inu for the photon package, based on the stellar
        ! flux
        !
        rn = ran2(iseed)
        call hunt(illum_flux_unprojected_cum(:,illum),freq_nr+1,rn,ray_inu)
        if(debug_check_all.eq.1) then
           if((ray_inu.lt.1).or.(ray_inu.gt.freq_nr)) then
              write(stdo,*) 'INTERNAL ERROR: ray_inu out of range for illumination beams'
              stop
           endif
        endif
     endif
     !
  elseif(rn.ge.mc_cumlum4) then
     !
     ! ----- We have a photon emitted by one of the thermal boundaries -----
     !       NOTE: This should only happen in Cartesian coordinates
     !
     ! Find which boundary
     !
     rn = ran2(iseed)*mc_bc_lumcum(7)
     call hunt(mc_bc_lumcum,7,rn,ibnd)
     if(ibnd.eq.1) then
        bc_idir = 1
        bc_ilr  = 1
     elseif(ibnd.eq.2) then
        bc_idir = 1
        bc_ilr  = 2
     elseif(ibnd.eq.3) then
        bc_idir = 2
        bc_ilr  = 1
     elseif(ibnd.eq.4) then
        bc_idir = 2
        bc_ilr  = 2
     elseif(ibnd.eq.5) then
        bc_idir = 3
        bc_ilr  = 1
     elseif(ibnd.eq.6) then
        bc_idir = 3
        bc_ilr  = 2
     else
        stop 72
     endif
     !
     ! Compute Lambertian random scattering direction
     !
     rn            = ran2(iseed)
     dir_perp      = sqrt(rn)
     dum           = twopi*ran2(iseed)
     rn            = sqrt(1.d0-rn)
     dir_planex    = rn*cos(dum)
     dir_planey    = rn*sin(dum)
     !
     ! Now handle all cases
     !
     if(bc_idir.eq.1) then
        if((igrid_coord.eq.10).or.(igrid_coord.eq.20)) then
           write(stdo,*) 'INTERNAL ERROR in thermal boundaries.'
           stop
        endif
        ray_cart_diry = dir_planex
        ray_cart_dirz = dir_planey
        if(bc_ilr.eq.1) then
           epsdist       = 1d-4*(amr_grid_xi(2,1)-amr_grid_xi(1,1))/ &
                           (2**amr_levelmax)
           ray_cart_x    = amr_grid_xi(1,1)+epsdist
           ray_cart_dirx = dir_perp
        else
           epsdist       = 1d-4*(amr_grid_xi(amr_grid_nx,1)-amr_grid_xi(amr_grid_nx-1,1))/ &
                           (2**amr_levelmax)
           ray_cart_x    = amr_grid_xi(amr_grid_nx+1,1)-epsdist
           ray_cart_dirx = -dir_perp
        endif
        rn = ran2(iseed)
        ray_cart_y = amr_grid_xi(1,2) + rn * &
             ( amr_grid_xi(amr_grid_ny,2) - amr_grid_xi(1,2) )
        rn = ran2(iseed)
        ray_cart_z = amr_grid_xi(1,3) + rn * &
             ( amr_grid_xi(amr_grid_nz,3) - amr_grid_xi(1,3) )
     elseif(bc_idir.eq.2) then
        if(igrid_coord.eq.10) then
           write(stdo,*) 'INTERNAL ERROR in thermal boundaries.'
           stop
        endif
        ray_cart_dirx = dir_planex
        ray_cart_dirz = dir_planey
        if(bc_ilr.eq.1) then
           epsdist       = 1d-4*(amr_grid_xi(2,2)-amr_grid_xi(1,2))/ &
                           (2**amr_levelmax)
           ray_cart_y    = amr_grid_xi(1,2)+epsdist
           ray_cart_diry = dir_perp
        else
           epsdist       = 1d-4*(amr_grid_xi(amr_grid_ny,2)-amr_grid_xi(amr_grid_ny-1,2))/ &
                           (2**amr_levelmax)
           ray_cart_y    = amr_grid_xi(amr_grid_ny+1,2)-epsdist
           ray_cart_diry = -dir_perp
        endif
        rn = ran2(iseed)
        if(igrid_coord.eq.10) then
           ray_cart_x = 0.d0
        else
           ray_cart_x = amr_grid_xi(1,1) + rn * &
                ( amr_grid_xi(amr_grid_nx,1) - amr_grid_xi(1,1) )
        endif
        rn = ran2(iseed)
        ray_cart_z = amr_grid_xi(1,3) + rn * &
             ( amr_grid_xi(amr_grid_nz,3) - amr_grid_xi(1,3) )
     else
        ray_cart_dirx = dir_planex
        ray_cart_diry = dir_planey
        if(bc_ilr.eq.1) then
           epsdist       = 1d-4*(amr_grid_xi(2,3)-amr_grid_xi(1,3))/ &
                           (2**amr_levelmax)
           ray_cart_z    = amr_grid_xi(1,3)+epsdist
           ray_cart_dirz = dir_perp
        else
           epsdist       = 1d-4*(amr_grid_xi(amr_grid_nz,3)-amr_grid_xi(amr_grid_nz-1,3))/ &
                           (2**amr_levelmax)
           ray_cart_z    = amr_grid_xi(amr_grid_nz+1,3)-epsdist
           ray_cart_dirz = -dir_perp
        endif
        rn = ran2(iseed)
        if((igrid_coord.eq.10).or.(igrid_coord.eq.20)) then
           ray_cart_x = 0.d0
        else
           ray_cart_x = amr_grid_xi(1,1) + rn * &
                ( amr_grid_xi(amr_grid_nx,1) - amr_grid_xi(1,1) )
        endif
        rn = ran2(iseed)
        if(igrid_coord.eq.10) then
           ray_cart_y = 0.d0
        else
           ray_cart_y = amr_grid_xi(1,2) + rn * &
             ( amr_grid_xi(amr_grid_ny,2) - amr_grid_xi(1,2) )
        endif
     endif
     !
     ! Now let the AMR module find the cell index
     !
     if(amr_tree_present) then
        call amr_findcell(ray_cart_x,ray_cart_y,ray_cart_z,acell)
        ray_index = acell%leafindex
     else
        call amr_findbasecell(ray_cart_x,ray_cart_y,ray_cart_z,ix,iy,iz)
        ray_index = ix+((iy-1)+(iz-1)*amr_grid_ny)*amr_grid_nx
     endif
     !
     ! Find random frequency
     !
     rn = ran2(iseed) 
     call hunt(mc_bc_cumspec(:,bc_ilr,bc_idir),freq_nr+1,rn,ray_inu)
     if(debug_check_all.eq.1) then
        if((ray_inu.lt.1).or.(ray_inu.gt.freq_nr)) then
           write(stdo,*) 'INTERNAL ERROR: ray_inu out of range for thermal boundary'
           stop
        endif
     endif
     !
     ! If weighted photon packages, then get the energy
     !
     if(mc_weighted_photons) energy = mc_energy_bc
     !
  elseif(rn.ge.mc_cumlum3) then
     !
     ! ----- We have a photon emitted by the continuous stellar source -----
     !       (for simulations of galaxies and such)
     !
     ! First determine which cell to emit from
     !
     rn = ran2(iseed)
     call hunt(stellarsrc_lumcum,nrcells+1,rn,icell)
     ray_index = cellindex(icell)
     !
     ! Find the position within the cell from where to emit
     !
     call montecarlo_random_pos_in_cell(icell)
     !
     ! Find the direction in which to emit
     !
     call montecarlo_randomdir(ray_cart_dirx,ray_cart_diry,ray_cart_dirz)
     !
     ! Now determine from which stellar population we will send the photon
     !
     if(stellarsrc_nrtemplates.eq.1) then
        itemplate = 1
     else
        dum = 0.d0
        do itemplate=1,stellarsrc_nrtemplates
           stellarsrc_cum(itemplate) = dum
           dum = dum + stellarsrc_dens(itemplate,ray_index) *           &
                       stellarsrc_templateslum(itemplate)
        enddo
        dum = 1.d0 / dum
        do itemplate=1,stellarsrc_nrtemplates
           stellarsrc_cum(itemplate) = stellarsrc_cum(itemplate) * dum
        enddo
        stellarsrc_cum(stellarsrc_nrtemplates+1) = 1.d0
        rn = ran2(iseed)
        call hunt(stellarsrc_cum,stellarsrc_nrtemplates+1,rn,itemplate)
        if(debug_check_all.eq.1) then
           if((itemplate.lt.1).or.(itemplate.gt.stellarsrc_nrtemplates)) then
              write(stdo,*) 'INTERNAL ERROR: Itemplate out of range'
              stop
           endif
        endif
     endif
     !
     ! Now determine which frequency we send the photon out
     !
     rn = ran2(iseed)
     call hunt(stellarsrc_speccum(:,itemplate),freq_nr+1,rn,ray_inu)
     if(debug_check_all.eq.1) then
        if((ray_inu.lt.1).or.(ray_inu.gt.freq_nr)) then
           write(stdo,*) 'INTERNAL ERROR: ray_inu out of range for continuous stellar source'
           stop
        endif
     endif
     !
     ! If weighted photon packages, then get the energy
     !
     if(mc_weighted_photons) energy = mc_energy_stellarsrc
     !
  elseif(rn.ge.mc_cumlum2) then
     !
     ! ----- We have a photon emitted by the external environment into -----
     !       our computational domain. This is done from a sphere 
     !       centered around the center of the coordinate system.
     !
     ray_index = 0
     !
     ! Find a random point on the sphere
     !
     call montecarlo_randomdir(pdirx,pdiry,pdirz)
     ray_cart_x  = grid_contsph_x + pdirx*grid_contsph_r
     ray_cart_y  = grid_contsph_y + pdiry*grid_contsph_r
     ray_cart_z  = grid_contsph_z + pdirz*grid_contsph_r
     !
     ! If mirror symmetry, then enforce this here
     !
     if(igrid_mirror.eq.1) then
        if(ray_cart_z.lt.0.d0) then
           ray_cart_z = -ray_cart_z
        endif
     endif
     !
     ! Now determine a random direction for the photon. Here we again
     ! have the angle difficulty we encountered for emission from a
     ! stellar surface. See above by sqrt(rn) is used. Lambertian 
     ! emission.
     ! 
     rn            = ran2(iseed)
     ray_cart_dirx = sqrt(rn)
     dum           = twopi*ran2(iseed)
     rn            = sqrt(1.d0-rn)
     ray_cart_diry = rn*cos(dum)
     ray_cart_dirz = rn*sin(dum)
     call montecarlo_rotatevec(ray_cart_dirx,ray_cart_diry,ray_cart_dirz, &
                               -pdirx,-pdiry,-pdirz)
     !
     ! Now determine which frequency we send the photon out
     !
     rn = ran2(iseed)
     call hunt(extlum_intens_cum(:),freq_nr+1,rn,ray_inu)
     if(debug_check_all.eq.1) then
        if((ray_inu.lt.1).or.(ray_inu.gt.freq_nr)) then
           write(stdo,*) 'INTERNAL ERROR: ray_inu out of range for external radiation'
           stop
        endif
     endif
     !
     ! If weighted photon packages, then get the energy
     !
     if(mc_weighted_photons) energy = mc_energy_extlum
     !
  elseif(rn.ge.mc_cumlum1) then
     !
     ! ----- We have a photon emitted by quantum-heated grains -----
     !
     ! First determine which cell to emit from
     !
     rn = ran2(iseed)
     call hunt(emisquant_cum,nrcells+1,rn,icell)
     ray_index = cellindex(icell)
     !
     ! Find the position within the cell from where to emit
     !
     call montecarlo_random_pos_in_cell(icell)
     !
     ! Find the direction in which to emit
     !
     call montecarlo_randomdir(ray_cart_dirx,ray_cart_diry,ray_cart_dirz)
     !
     ! Find the frequency at which to emit
     !
     rn = ran2(iseed)
     call hunt(emisquant_loccum(:,ray_index),freq_nr+1,rn,ray_inu)
     if(debug_check_all.eq.1) then
        if((ray_inu.lt.1).or.(ray_inu.gt.freq_nr)) then
           write(stdo,*) 'INTERNAL ERROR: ray_inu out of range for quantum source'
           stop
        endif
     endif
     !
     ! If weighted photon packages, then get the energy
     !
     if(mc_weighted_photons) energy = mc_energy_quant
     !
  else
     !
     ! ----- We have a photon emitted by internal heat production (e.g. -----
     !       viscous heating in an accretion disk)
     !
     ! First determine which cell to emit from
     !
     rn = ran2(iseed)
     call hunt(heatsource_cum,nrcells+1,rn,icell)
     ray_index = cellindex(icell)
     !
     ! Find the position within the cell from where to emit
     !
     call montecarlo_random_pos_in_cell(icell)
     !
     ! The photon will be emitted by thermal radiation, and will therefore
     ! also affect the local temperature accordingly.  This is
     ! mathematically identical to a discrete absorption and re-emission
     ! event. So therefore we call the discrete absorption-and-reemission
     ! event routine, even though this is an only-emission event. But that
     ! does not make a difference. Basically, the reason why the temperature
     ! is increased here is to ensure that the cell can indeed launch this 
     ! photon package.
     !
     ! NOTE: In the future it may be important to rethink the division of
     !       this primary energy over the dust species. I now put ray_inu=-1,
     !       which is a signal to do_asorption_event() to divide the energy
     !       by the density. This may not be correct since some species may
     !       be much more capable of reemitting energy than others. For now
     !       I simply divide it per unit dust mass because I have no better
     !       idea.
     !
     ! If weighted photon packages, then get the energy
     ! Bugfix 22.03.2016 by Jon Ramsey
     !
     if(mc_weighted_photons) energy = mc_energy_heatsource
     !
     ! Add this energy immediately to the local energy of the cell.
     ! (Bugfix 16.09.2016: This was omitted in earlier versions, 
     ! which is for optically very thick models not so dramatic, but
     ! for optically thin models it is wrong)
     !
     mc_cumulener(1,ray_index) = mc_cumulener(1,ray_index) + energy
     !
     ! Now call the do_absorption_event() to recalculate, if necessary,
     ! the dust temperature and find a new random direction.
     !
     ray_inu  = -1    
     !write(stdo,*) 'CHECK THIS. STOPPING FOR NOW.'
     !stop
     iqactive = 0
     call do_absorption_event(params,iqactive,ierror)
     !
     ! Find the direction in which to emit
     !
     call montecarlo_randomdir(ray_cart_dirx,ray_cart_diry,ray_cart_dirz)
     !
     ! The locally produced heat is usually radiated at 
     ! long wavelengths and thus will not excite PAHs.
     ! So switch off quantum (in old quantum mode).
     ! For new quantum mode, quantum is always active.
     !
     if((incl_quantum.lt.0).and.(incl_quantum.gt.-10)) then
        iqactive = 0
     elseif(incl_quantum.eq.-10) then
        iqactive = -1
     elseif(incl_quantum.eq.-11) then
        iqactive = 2
     else
        iqactive = 0
     endif
     mc_photon_destroyed = .false.
     !
  endif
  !
  ! Now do some checks...
  !
  if(debug_check_all.eq.1) then
     if((ray_inu.lt.1).or.(ray_inu.gt.freq_nr)) then 
        write(stdo,*) 'INTERNAL ERROR: ray_inu.lt.1 or ray_inu.gt.mc_f: rayinu = ',ray_inu
        stop
     endif
     if((ray_index.lt.1).or.(ray_index.gt.nrcellsmax)) then
        write(stdo,*) 'INTERNAL ERROR: index out of range: ray_index = ',ray_index
        stop
     endif
     if(abs(ray_cart_dirx**2+ray_cart_diry**2+ray_cart_dirz**2-1.d0).gt.1d-14) then
        write(stdo,*) 'INTERNAL ERROR: dir vector not length unity: dir = ', &
               ray_cart_dirx,ray_cart_diry,ray_cart_dirz
        stop
     endif
  endif
  !
  ! Now associate pointer to cell.
  ! This is grid type dependent.
  !
  if(igrid_type.lt.100) then
     !
     ! AMR type grid
     !
     if(ray_index.gt.0) then
        !
        ! We are in a cell
        !
        if(amr_tree_present) then
           amrray_cell => amr_index_to_leaf(ray_index)%link
        else
           call amr_regular_get_ixyz(ray_index,amrray_ix_curr,amrray_iy_curr,amrray_iz_curr)
        endif
     else
        nullify(amrray_cell)
        amrray_ix_curr = -1
        amrray_iy_curr = -1
        amrray_iz_curr = -1
     endif
  else
     write(stdo,*) 'ERROR: This grid type not yet implemented'
     stop
  endif
  !
  ! In case we use polarization, init the photpkg structure 
  !
  if(scattering_mode.ge.5) then
     photpkg%E    = energy
     photpkg%Q    = 0.d0
     photpkg%U    = 0.d0
     photpkg%V    = 0.d0
     photpkg%n(1) = ray_cart_dirx
     photpkg%n(2) = ray_cart_diry
     photpkg%n(3) = ray_cart_dirz
     call polarization_make_s_vector(photpkg%n,photpkg%s)
  endif
  !
  !-----------------------------------------------
  ! Now start the random walk
  !-----------------------------------------------
  !
  ok = .true.
  count_samecell = 0
  do while(ok)
     !
     ! Backup the current index
     !
     index_prev = ray_index
     !
     ! Get a random tau
     !
     taupath = -log(1.d0-ran2(iseed))
     !
     ! Move photon to next scattering/absorption event
     ! 
     call walk_cells_thermal(params,taupath,iqactive,arrived,therm,ispec,ierror)
     !
     ! If we are still in the same cell, then we have to count how
     ! often we have been in the same cell.
     !
     if(arrived.or.mc_photon_destroyed.or.(ray_index.ne.index_prev)) then
        !
        ! Either we have not been in this cell last time (in which case we
        ! must reset the counter), or we end our journey. In either case
        ! we update the global counters for the statistics.
        !
        mc_visitcell   = mc_visitcell + count_samecell + 1
        mc_revisitcell = mc_revisitcell + count_samecell*count_samecell
        mc_revisitcell_max = max(mc_revisitcell_max,1.d0*count_samecell)
        count_samecell = 0
     else
        !
        ! We have been in this cell the last time around. So, increase the
        ! counter
        !
        count_samecell = count_samecell + 1
     endif
     !
     ! If arrived at end-point or escaped to infinity, then return
     !
     if(arrived) then
        return
     endif
     !
     ! If error, then also return
     !
     if(ierror.ne.0) then
        stop 1142
        return
     endif
     !
     ! If photon destroyed (by e.g. quantum heated grains which are
     ! dealt with elsewhere), then also return
     !
     if(mc_photon_destroyed) then
        return
     endif
     !
     ! Was it an absorption event or scattering event?
     !
     if(therm) then
        !
        ! Absorption event
        !
        call do_absorption_event(params,iqactive,ierror)
        call montecarlo_randomdir(ray_cart_dirx,ray_cart_diry,ray_cart_dirz)
        !
        ! Reset photpkg
        !
        if(scattering_mode.ge.5) then
           photpkg%E    = energy
           photpkg%Q    = 0.d0
           photpkg%U    = 0.d0
           photpkg%V    = 0.d0
           photpkg%n(1) = ray_cart_dirx
           photpkg%n(2) = ray_cart_diry
           photpkg%n(3) = ray_cart_dirz
           call polarization_make_s_vector(photpkg%n,photpkg%s)
        endif
        !
        ! If photon destroyed (by e.g. quantum heated grains which are
        ! dealt with elsewhere), then also return
        !
        if(mc_photon_destroyed) then
           return
        endif
        !
        ! Switch off quantum
        !
        ! ONLY for old quantum method
        !
        if(incl_quantum.gt.-10) then
           iqactive = 0
        endif
        !
     else
        !
        ! Scattering event
        !
        if(scattering_mode.eq.1) then
           !
           ! Isotropic scattering
           ! 
           call montecarlo_randomdir(ray_cart_dirx,ray_cart_diry,ray_cart_dirz)
           !
        elseif(scattering_mode.eq.2) then
           !
           ! Anisotropic scattering off grain for dust species ispec,
           ! following the Henyey-Greenstein recipe
           !
           g = kappa_g(ray_inu,ispec)
           call henyeygreenstein_randomdir(g,dirnewx,dirnewy,dirnewz)
           call montecarlo_rotatevec(dirnewx,dirnewy,dirnewz,   &
                   ray_cart_dirx,ray_cart_diry,ray_cart_dirz)
           ray_cart_dirx = dirnewx
           ray_cart_diry = dirnewy
           ray_cart_dirz = dirnewz
           !
        elseif((scattering_mode.eq.3).or.(scattering_mode.eq.4)) then
           !
           ! Anisotropic scattering off grain for dust species ispec,
           ! following the tabulated Z matrix. 
           !
           call anisoscat_randomdir(ray_inu,ispec,dirnewx,dirnewy,dirnewz)
           call montecarlo_rotatevec(dirnewx,dirnewy,dirnewz,   &
                   ray_cart_dirx,ray_cart_diry,ray_cart_dirz)
           ray_cart_dirx = dirnewx
           ray_cart_diry = dirnewy
           ray_cart_dirz = dirnewz
           !
        elseif(scattering_mode.eq.5) then
           !
           ! Full polarized scattering
           !
           ! Note: This subroutine will also modify the polarization state
           ! of photpkg. Also note that this subroutine does not need any
           ! transformation to and from the particle frame as was necessary
           ! in the scattering_mode<5 cases.
           !
           call polarization_randomorient_mc_scatter(photpkg,  &
                scat_munr,scat_mui_grid,scat_thetai_grid,      &
                zmatrix(:,ray_inu,1,ispec),                    &
                zmatrix(:,ray_inu,2,ispec),                    &
                zmatrix(:,ray_inu,3,ispec),                    &
                zmatrix(:,ray_inu,4,ispec),                    &
                zmatrix(:,ray_inu,5,ispec),                    &
                zmatrix(:,ray_inu,6,ispec),                    &
                zcumul(:,ray_inu,1,ispec),                     &
                zcumul(:,ray_inu,2,ispec),                     &
                .false.)
           ray_cart_dirx = photpkg%n(1)
           ray_cart_diry = photpkg%n(2)
           ray_cart_dirz = photpkg%n(3)
           !
        else
           if(scattering_mode.eq.0) then
              write(stdo,*) 'ERROR in Monte Carlo: encountered a scattering event but '
              write(stdo,*) '     the scattering mode is 0.'
              stop
           else
              write(stdo,*) 'ERROR: scattering mode ',scattering_mode,' not known'
              stop
           endif
        endif
     endif
     !
     ! For debugging: Increase event counter
     !
     ieventcounttot = ieventcounttot + 1
     !
     ! If requested, see if we can do a Modified Random Walk (MRW) from
     ! this point onward until we exit the cell again. This is 
     ! useful if the cell has a very high optical depth and therefore
     ! making the photon package to undergo many events, effectively
     ! "trapping" the photon package. The MRW method of Min et al. (2009)
     ! can drastically speed this up.
     !
     ! NOTE: So far only for regular and AMR-type grids, not for any
     !       possible future irregular grids.
     !
     if(params%mod_random_walk.and.(igrid_type.lt.100)) then
        !
        ! Check if the stellar spheres are switched on, because the MRW method
        ! cannot deal with spheres inside grid cells. If a sphere is inside
        ! the present grid cell, the Monte Carlo will be done in the usual
        ! slow way for this particular cell. For all other cells it should
        ! still work normally.
        !
        nospheres = .true.
        if(allocated(amrray_spheres_sphidx)) then
           if(amrray_spheres_sphidx(ray_index).ne.0) then
              nospheres = .false.
           endif
        endif
        if(nospheres) then
           !
           ! OK, no spheres in this cell, so we can proceed.
           !
           ! If we stay in the same cell for a minimum number of times, then we
           ! may want to switch to Modified Random Walk.
           !
           ! However, if we already earmarked the cell as a MRW cell, then we
           ! go straight into the MRW unless the new Rosseland opacity is
           ! too low.
           !
           if((count_samecell.ge.params%mrw_count_trigger).or. &
                mrw_cell_uses_mrw(ray_index)) then
              !
              ! OpenMP: Lock this cell
              !
              !$ continue=.true.
              !$ do while(.NOT. omp_test_lock(lock(ray_index)))
              !$    if(continue)then
              !$omp atomic
              !$       conflict_counter=conflict_counter+1
              !$       continue=.false.
              !$    end if
              !$ end do
              !
              ! Set the dust temperature to -1, so that we can see if
              ! we have already calculated it or not (because we don't
              ! want to calculate it if not necessary)
              !
              mrw_dusttemp  = -1.d0
              !
              ! Get the cumulative energy for this cell
              !
              mrw_energy = 0.d0
              do ispec=1,dust_nr_species
                 if((dust_quantum(ispec).eq.0).or.(iqactive.le.0)) then 
                    mrw_energy = mrw_energy + mc_cumulener(ispec,ray_index)
                 endif
              enddo
              !
              ! Check if we have to recalculate the Rosseland opacities or
              ! if we can re-use old ones
              !
              if((abs((mrw_energy/(mrw_cumulener_bk(ray_index)+1d-90))-1.d0).gt. &
                   params%mrw_enthres).or.(mrw_alpha_tot_ross(ray_index).le.0.d0)) then
                 !
                 ! Yes, we must recalculate the Rosseland mean opacities
                 !
                 ! First calculate the temperature
                 !
                 mrw_dusttemp = compute_dusttemp_coupled_bd(dust_nr_species, &
                      mrw_energy,dustdens(:,ray_index),cellvolume(ray_index))
                 !
                 ! With this temperature, calculate the RM opacities
                 !
                 call mrw_calculate_rossmean(dust_nr_species,mrw_dusttemp, &
                                             dustdens(:,ray_index),        &
                                             alpha_t_rm)
                 !
                 ! Store these
                 !
                 mrw_alpha_tot_ross(ray_index) = alpha_t_rm
                 !
                 ! Update the backup
                 !
                 mrw_cumulener_bk(ray_index) = mrw_energy
              else
                 !
                 ! No, we can re-use the ones stored in the arrays; they
                 ! are still sufficiently up-to-date
                 !
                 alpha_t_rm    = mrw_alpha_tot_ross(ray_index)
              endif
              !
              ! Now create some stuff that the MRW subroutines 
              ! need
              !
              if(amr_tree_present) then
                 if(.not.associated(amrray_cell)) stop 7270
                 do iddr=1,3
                    cellx0(iddr) = amr_finegrid_xi(amrray_cell%ixyzf(iddr),iddr,amrray_cell%level)
                    cellx1(iddr) = amr_finegrid_xi(amrray_cell%ixyzf(iddr)+1,iddr,amrray_cell%level)
                 enddo
              else
                 if(amrray_ix_curr.le.0) stop 7270
                 cellx0(1) = amr_finegrid_xi(amrray_ix_curr,1,0)
                 cellx1(1) = amr_finegrid_xi(amrray_ix_curr+1,1,0)
                 cellx0(2) = amr_finegrid_xi(amrray_iy_curr,2,0)
                 cellx1(2) = amr_finegrid_xi(amrray_iy_curr+1,2,0)
                 cellx0(3) = amr_finegrid_xi(amrray_iz_curr,3,0)
                 cellx1(3) = amr_finegrid_xi(amrray_iz_curr+1,3,0)
              endif
              !
              ! Set idim and icoord
              !
              idim   = 3
              icoord = igrid_coord
              if(icoord.lt.100) then
                 !
                 ! Cartesian
                 !
                 ! For cartesian we always assume 3-D. This may
                 ! therefore not work. Check for this. Presumably
                 ! it is easy to fix this, but for now let's first
                 ! get this stuff working in the first place.
                 !
                 if((icoord.eq.10).or.(icoord.eq.20)) then
                    write(stdo,*) 'ERROR in Modified Random Walk: ', &
                         'For now not compatible with 1-D plane-parallel slab, '
                    write(stdo,*) '      or 2-D pencil-parallel slab geometry.'
                    stop
                 endif
              else
                 !
                 ! Spherical
                 !
                 if(amr_zdim.eq.0) idim=2
                 if(amr_ydim.eq.0) idim=1
              endif
              !
              ! Now we have to check if the cell is indeed optically thick
              ! enough to even start thinking of MRW. Note also that if
              ! we do MRW, we will automatically put the dust temperatures
              ! of all dust species in this cell to the same value. We do
              ! not usually want this to happen unless the cell is very 
              ! optically thick.
              !
              ! So, first compute the minimal size of this cell
              !
              if(icoord.le.99) then
                 !
                 ! Cartesian coordinates
                 !
                 call find_smallest_size_cart_lite(cellx0,cellx1,size)
              else
                 !
                 ! Spherical coordinates
                 !
                 call find_smallest_size_spher_lite(cellx0,cellx1,size,idim)
              endif
              !
              ! Now decide
              !
              if(size*alpha_t_rm.gt.params%mrw_tauthres) then
                 !
                 ! OK, so we now start the MRW until we escape the cell
                 !
                 ! Now earmark this cell as MRW cell
                 !
                 mrw_cell_uses_mrw(ray_index) = .true.
                 !
                 ! If the dust temperature has not yet been calculated, 
                 ! then do it here.
                 !
                 if(mrw_dusttemp.lt.0.d0) then
                    mrw_dusttemp = compute_dusttemp_coupled_bd(dust_nr_species, &
                         mrw_energy,dustdens(:,ray_index),cellvolume(ray_index))
                 endif
                 !
                 ! Get the Planck mean opacities from the database
                 !
                 call mrw_get_planckopac_dbase(mrw_dusttemp,dust_nr_species,alpha_a_pm)
                 alpha_a_pm_tot = 0.d0
                 do ispec=1,dust_nr_species
                    alpha_a_pm(ispec) = alpha_a_pm(ispec) * dustdens(ispec,ray_index)
                    alpha_a_pm_tot    = alpha_a_pm_tot + alpha_a_pm(ispec)
                 enddo
                 !
                 ! Determine the energy threshold, i.e. the factor (1+mrw_enthres) by which
                 ! the energy of the cell can increase until the MRW routine should
                 ! stop to recalculate the mean opacities.
                 !
                 ! However, if we enter the cell for the first (or nearly first) time,
                 ! this could lead to a very low enthres, meaning a very often stopping
                 ! of the MRW routine. That is also not very useful. So if the temperature
                 ! of the cell currently is smaller than a threshold value, then put
                 ! enthres to 1d99 (=no stopping).
                 ! 
                 ! NOTE (16.09.2016): This no stopping policy is not really the most
                 !      correct method, because it means that the first photon entering
                 !      a very optically thick cell will keep using the mean opacities
                 !      for the lowest temperature even when the temperature may increase
                 !      quite much. BUT: If a single photon increases the temperature so
                 !      much and no further photons enter the cell, then the photon
                 !      statistics is anyway so bad that it does not matter. And since
                 !      dust opacities for low temperatures are generally lower than for
                 !      high temperatures, if anything it will slightly underpredict the
                 !      cell temperature, but it cannot lead to a runaway effect. 
                 !      Anyway: at some point in the future it might be more elegant to
                 !      compute the enthres for that case to be the energy corresponding
                 !      to the temperature = mrw_tempthres.
                 !
                 if(mrw_dusttemp.gt.params%mrw_tempthres) then
                    enthres = mrw_cumulener_bk(ray_index) * ( 1.d0 + params%mrw_enthres )
                 else
                    enthres = 1d99
                 endif
                 !
                 ! Perform the MRW now
                 !
                 pos(1)   = ray_cart_x
                 pos(2)   = ray_cart_y
                 pos(3)   = ray_cart_z
                 dir(1)   = ray_cart_dirx
                 dir(2)   = ray_cart_diry
                 dir(3)   = ray_cart_dirz
                 enerphot = energy
                 !
                 ! ...but first check if we are indeed inside the cell (this should
                 !    be always the case, but round-off errors can occasionally 
                 !    cause troubles).
                 !
                 incell = mrw_check_if_in_cell(cellx0,cellx1,pos,icoord,idim)
                 !
                 ! If inside of cell, then do MRW, otherwise jiggle the photon position a bit
                 !
                 if(incell) then
                    !
                    ! ...yes we are inside the cell, so let's go 
                    !
                    call modified_random_walk(cellx0,cellx1,pos,dir,mrw_energy,enerphot, &
                                   enthres,alpha_a_pm_tot,alpha_t_rm,params%mrw_gamma,   &
                                   icoord,idim,idir_cross,ilr_cross,nrevents,            &
                                   taumargin=params%mrw_taustepback)
                    ! 
                    ! Calculate the dust temperature one last time
                    !
                    mrw_dusttemp = compute_dusttemp_coupled_bd(dust_nr_species, &
                         mrw_energy,dustdens(:,ray_index),cellvolume(ray_index))
                    !
                    ! Copy this to all the dust temperatures
                    !
                    dusttemp(:,ray_index) = mrw_dusttemp
                    !
                    ! We redistribute the energy over the dust species
                    !
                    mrw_dcumen(:) = mc_cumulener(:,ray_index)
                    do ispec=1,dust_nr_species
                            !$OMP CRITICAL
                       mc_cumulener(ispec,ray_index) = mrw_energy *  &
                            alpha_a_pm(ispec) / alpha_a_pm_tot
                            !$OMP END CRITICAL
                    enddo
                    !
                    ! We calculate the delta cumulener
                    !
                    mrw_dcumen(:) = mc_cumulener(:,ray_index) - mrw_dcumen(:)
                    !
                    ! Physically it should be impossible to have one of these
                    ! becoming negative, but numerically you never know. In
                    ! that case, put that to zero
                    !
                    do ispec=1,dust_nr_species
                       if(mrw_dcumen(ispec).lt.0.d0) mrw_dcumen(ispec)=0.d0
                    enddo
                    !
                    ! We draw a random frequency based on the current temperature
                    ! of the dust in the cell. We use again the dB_nu(T)/dT as the
                    ! weighting function for this, according to the Bjorkman & Wood
                    ! algorithm. 
                    !
                    call pick_randomfreq_db(dust_nr_species,dusttemp(:,ray_index), &
                                            mrw_dcumen(:),ray_inu)
                 else
                    !
                    ! ...no, we are not inside the cell
                    !
                    write(stdo,*) 'WARNING: Photon outside of cell before calling Modified Random Walk.'
                    write(stdo,*) '         Jiggling photon position and continuing...'
                    !
                    !    Let's jiggle the position a bit
                    !
                    rn         = ran2(iseed)
                    ray_cart_x = ray_cart_x + (rn-0.5d0)*size*1d-6
                    rn         = ran2(iseed)
                    ray_cart_y = ray_cart_y + (rn-0.5d0)*size*1d-6
                    rn         = ran2(iseed)
                    ray_cart_z = ray_cart_z + (rn-0.5d0)*size*1d-6
                    pos(1)     = ray_cart_x
                    pos(2)     = ray_cart_y
                    pos(3)     = ray_cart_z
                    !
                    ! ...Now let the AMR module find the cell index
                    !
                    if(igrid_coord.ge.100) then
                       write(stdo,*) 'ERROR: The Modified Random Walk is called outside of cell.'
                       write(stdo,*) '       Please warn the author of RADMC-3D.'
                       stop
                    endif
                    if(amr_tree_present) then
                       call amr_findcell(ray_cart_x,ray_cart_y,ray_cart_z,acell)
                       ray_index = acell%leafindex
                    else
                       call amr_findbasecell(ray_cart_x,ray_cart_y,ray_cart_z,ix,iy,iz)
                       ray_index = ix+((iy-1)+(iz-1)*amr_grid_ny)*amr_grid_nx
                    endif
                 endif
                 !
                 ! Put the new position and direction back into the global variables
                 !
                 ray_cart_x    = pos(1)
                 ray_cart_y    = pos(2)
                 ray_cart_z    = pos(3)
                 ray_cart_dirx = dir(1)
                 ray_cart_diry = dir(2)
                 ray_cart_dirz = dir(3)
              endif
              !
              ! OpenMP unlock
              !
              !$ call omp_unset_lock(lock(ray_index));
              !
           endif
        endif
     endif
     !
     ! Next event
     !
  end do
  !
end subroutine walk_full_path_bjorkmanwood


!--------------------------------------------------------------------------
!         WALK THE RANDOM WALK UNTIL PHOTON ENERGY IS LOW OR ESCAPED
!
!     This routine does the scatterings while keeping track of 
!     the energy lost by a package. It only stops after leaving 
!     the grid.
!--------------------------------------------------------------------------
subroutine walk_full_path_scat(params,inu,ierror)
  implicit none
  type(mc_params) :: params
  integer :: ierror,inu
  !
  double precision :: taupath,ener,sint,epsdist
  double precision :: dir(1:3),enold,albedo,dum,g,rn
  double precision :: pdirx,pdiry,pdirz,dirnewx,dirnewy,dirnewz
  double precision :: dir_perp,dir_planex,dir_planey,dummy
  integer :: ispec,iqactive,istar,icell,itemplate,iscatevent
  integer :: ibnd,bc_idir,bc_ilr,ix,iy,iz,illum,ierr
  logical :: ok,arrived,usesphere,todo_photpkg
  type(amr_branch), pointer :: acell
  !
  ! Do checks
  !
  if(debug_check_all.eq.1) then
     if((incl_quantum.eq.0).and.(emisquanttot.ne.0.d0)) then
        write(stdo,*) 'ERROR: Emisquant.ne.0 while quantum not switched on...'
        stop
     endif
     if((nstars.eq.0).and.(starlumtot.ne.0.d0)) then
        write(stdo,*) 'ERROR: Starlumtot.ne.0 while nstars==0...'
        stop
     endif
  endif
  !
  ! Reset
  !
  iqactive = 0
  mc_photon_destroyed = .false.
  ray_inu  = inu
  todo_photpkg = .true.
  !
  ! First choose whether to launch the photon from (one of) the star(s) or
  ! from the external radiation field or from inside the grid (due to the
  ! internally released energy) or from the locally emitted spectral
  ! source.. We make this choice on the basis of the total luminosity of all
  ! sources.
  !
  if(mc_cumlum5.gt.0.d0) then
     !
     ! There are more kinds of luminosity sources than just the star(s)
     !
     rn = ran2(iseed)
     !
  else
     !
     ! Only one source of luminosity: the star(s)
     ! By choosing rn=2 we avoid a costly random number generation and 
     ! we immediately go to the stellar emission
     !
     rn        = 2.d0
     !
  endif
  if(rn.ge.mc_cumlum5) then
     !
     ! ----- We have a photon emitted by a star (or, for 1dpp, an illum source) -----
     !
     if(nstars.ge.1) then
        !
        ! We have a star or multiple stars
        !
        ! Which star?
        !
        if(nstars.gt.1) then
           rn = ran2(iseed)
           call hunt(star_lumcum,nstars+1,rn,istar)
           if(istar.gt.nstars) stop
           if(istar.lt.1) stop
        else
           istar = 1
           if(nstars.eq.0) then
              write(stdo,*) 'INTERNAL ERROR with 0 stars in Monte Carlo module'
              stop
           endif
        endif
        !
        ! If we have the spherical star mode switched on, then 
        ! switch it on here as well
        !
        if(star_sphere) then
           usesphere = .true.
        else
           usesphere = .false.
        endif
        !
        ! If we have the weighted photons mode, we must do some
        ! extra stuff
        !
        if(mc_weighted_photons) then
           !
           ! Get the current energy of this photon package
           !
           energy = mc_energy_stars(istar)
           !
           ! If we have weighted photons mode, and a star is outside
           ! of the grid (more precise: outside of the containing sphere), 
           ! then we (at least for now) do not allow the treatment of the
           ! star as a sphere. NOTE: If the star lies inside the inner cavity
           ! of spherical coordinates, it is inside of the containing sphere,
           ! and it will be still treated as a sphere if that mode is on. This
           ! switching off of the sphere mode is only if the star is far away.
           !
           if(star_fraclum(istar).ne.1.d0) then
              usesphere = .false.
           endif
        endif
        !
        ! Choose a random theta at which the photon from the
        ! (spherically symmetric) star enters the inner edge of the
        ! coordinate system.
        !
        if(usesphere) then
           !
           ! Treat star as a sphere 
           !
           ! NOTE: With a finite stellar size, we also have to absorb photons
           !       that hit the star. 
           !
           ! NOTE: For now finite size stars are only allowed to be OUTSIDE
           !       the grid. Otherwise it becomes a nightmare.
           !
           ! Find a random point on the sphere
           !
           call montecarlo_randomdir(pdirx,pdiry,pdirz)
           ray_cart_x  = star_pos(1,istar) + pdirx*star_r(istar)*1.000000000001d0
           ray_cart_y  = star_pos(2,istar) + pdiry*star_r(istar)*1.000000000001d0
           ray_cart_z  = star_pos(3,istar) + pdirz*star_r(istar)*1.000000000001d0
           !
           ! If mirror symmetry, then enforce this here
           !
           if(igrid_mirror.eq.1) then
              if(ray_cart_z.lt.0.d0) then
                 ray_cart_z = -ray_cart_z
              endif
           endif
           !
           ! NOTE: For now finite size stars are only allowed to be OUTSIDE
           !       the grid. Otherwise it becomes a nightmare.
           !
           ray_index = 0
           !
           ! Now determine a random direction for the photon. One must,
           ! however, be careful here because the direction is NOT simply a 3-D
           ! random direction with all inward directions rejected.  The thing
           ! is this: the intensity of a blackbody surface is the same at all
           ! inclinations. However, if one would simply use a normal random
           ! direction (using montecarlo_randomdir()) then one would see a
           ! brighter intensity as the surface is viewed under a more grazing
           ! inclination. This is incorrect.  So to compensate for this one
           ! must have the change of a photon being emitted in direction
           ! my=cos(theta) to be proportional to mu.
           !
           ! NOTE: I have been confused again, here, with this sqrt(rn). But it
           ! is indeed correct, as far as I know!!  The example that made me
           ! understand again why this is so is the example of a transport
           ! company with fast trucks, sending out per day 10 trucks, and a
           ! firm with slow trucks, also sending out per day 10 trucks. Both
           ! firms have the same flux of goods, even though firm2 has slow
           ! trucks. Firm2 just requires more trucks to keep the flux
           ! going. That is why photons with mu\simeq 0 (slow) would carry the
           ! same flux as photons with mu\simeq 1 (fast) if they would be
           ! emitted equally often. Since we know that mu\simeq 0 photons carry
           ! less flux, they must be emitted less often. Hence P(mu)=mu,
           ! i.e. Pcum(mu)=mu^2, or so to say: mu=sqrt(Pcum)=sqrt(rn).
           ! 
           rn            = ran2(iseed)
           ray_cart_dirx = sqrt(rn)
           dum           = twopi*ran2(iseed)
           rn            = sqrt(1.d0-rn)
           ray_cart_diry = rn*cos(dum)
           ray_cart_dirz = rn*sin(dum)
           call montecarlo_rotatevec(ray_cart_dirx,ray_cart_diry,ray_cart_dirz, &
                                     pdirx,pdiry,pdirz)
           !
        else
           !
           ! The star is assumed to be a point source
           !
           ray_cart_x  = star_pos(1,istar)
           ray_cart_y  = star_pos(2,istar)
           ray_cart_z  = star_pos(3,istar)
           !
           ! If mirror symmetry, then enforce this here
           ! NOTE: This should not be necessary: It should have already 
           !       been taken care of. But for now let's keep it.
           !
           if(igrid_mirror.eq.1) then
!              if(ray_cart_z.lt.0.d0) then
!                 ray_cart_z = -ray_cart_z
!              endif
              if(ray_cart_z.ne.0.d0) then
                 write(stdo,*) 'ERROR: With mirror symmetry, stars must lie on z=0 plane.'
                 stop 732
              endif
           endif
           !
           ! Find random direction
           !
           if(star_fraclum(istar).eq.1.d0) then
              !
              ! Pure random direction in 4*pi
              !
              call montecarlo_randomdir(ray_cart_dirx,ray_cart_diry,  &
                                        ray_cart_dirz)
           else
              !
              ! Focused toward the model grid (the star lies outside the grid)
              !
              rn            = ran2(iseed)
              dum           = star_fraclum(istar)*rn
              ray_cart_dirx = 1.d0-2.d0*dum
              sint          = 2.d0*sqrt(dum-dum**2)
              dum           = twopi*ran2(iseed)
              ray_cart_diry = sint*cos(dum)
              ray_cart_dirz = sint*sin(dum)
              pdirx         = grid_contsph_x-star_pos(1,istar)
              pdiry         = grid_contsph_y-star_pos(2,istar)
              pdirz         = grid_contsph_z-star_pos(3,istar)
              dum           = 1.d0 / sqrt(pdirx**2+pdiry**2+pdirz**2)
              pdirx         = pdirx * dum
              pdiry         = pdiry * dum
              pdirz         = pdirz * dum
              call montecarlo_rotatevec(ray_cart_dirx,ray_cart_diry,ray_cart_dirz, &
                                        pdirx,pdiry,pdirz)
           endif
           !
           ! If we have mirror symmetry, then all stars MUST lie on the
           ! equatorial plane and their photons must point upward.
           !
           if(igrid_mirror.eq.1) then
              ray_cart_dirz = abs(ray_cart_dirz)
           endif
           !
           ! Which cell is this star in, if it is in the grid?
           !
           ray_index = star_cellindex(istar)
           !
        endif
     else
        !
        ! It is an illumination beam (only for 1-D plane-parallel)
        !
        ! Which illumination source?
        !
        if(nillum.gt.1) then
           rn = ran2(iseed)
           call hunt(illum_fluxtot_cum,nillum+1,rn,illum)
           if(illum.gt.nillum) stop
           if(illum.lt.1) stop
        else
           illum = 1
           if(nillum.eq.0) then
              write(stdo,*) 'INTERNAL ERROR with 0 illumination beams in Monte Carlo module'
              stop
           endif
        endif
        !
        ! If we have the weighted photons mode, we must do some
        ! extra stuff
        !
        if(mc_weighted_photons) then
           !
           ! Get the current energy of this photon package
           !
           energy = mc_energy_illum(illum)
        endif
        !
        ! Now put the photon package well above the grid to start
        !
        ray_cart_x  = 0.d0
        ray_cart_y  = 0.d0
        if(illum_costheta(illum).gt.0.d0) then
           ray_cart_z  = amr_grid_xi(amr_grid_nz+1,3)+abs(amr_grid_xi(amr_grid_nz+1,3)-amr_grid_xi(1,3))
        else
           ray_cart_z  = amr_grid_xi(1,3)-abs(amr_grid_xi(amr_grid_nz+1,3)-amr_grid_xi(1,3))
        endif
        !
        ! Set the direction. The direction is defined such that
        ! if theta = theta_observer and phi = phi_observer, the
        ! illuminating source is the observer. 
        !
        ray_cart_dirz = -illum_costheta(illum)
        dummy         = sqrt(1.d0-illum_costheta(illum)**2)
        ray_cart_diry = dummy * cos(illum_phi(illum))
        ray_cart_dirx = dummy * sin(illum_phi(illum))
     endif
     !
  elseif(rn.ge.mc_cumlum4) then
     !
     ! ----- We have a photon emitted by one of the thermal boundaries -----
     !       NOTE: This should only happen in Cartesian coordinates
     !
     ! Find which boundary
     !
     rn = ran2(iseed)*mc_bc_lumcum(7)
     call hunt(mc_bc_lumcum,7,rn,ibnd)
     if(ibnd.eq.1) then
        bc_idir = 1
        bc_ilr  = 1
     elseif(ibnd.eq.2) then
        bc_idir = 1
        bc_ilr  = 2
     elseif(ibnd.eq.3) then
        bc_idir = 2
        bc_ilr  = 1
     elseif(ibnd.eq.4) then
        bc_idir = 2
        bc_ilr  = 2
     elseif(ibnd.eq.5) then
        bc_idir = 3
        bc_ilr  = 1
     elseif(ibnd.eq.6) then
        bc_idir = 3
        bc_ilr  = 2
     else
        stop 72
     endif
     !
     ! Compute Lambertian random scattering direction
     !
     rn            = ran2(iseed)
     dir_perp      = sqrt(rn)
     dum           = twopi*ran2(iseed)
     rn            = sqrt(1.d0-rn)
     dir_planex    = rn*cos(dum)
     dir_planey    = rn*sin(dum)
     !
     ! Now handle all cases
     !
     if(bc_idir.eq.1) then
        if((igrid_coord.eq.10).or.(igrid_coord.eq.20)) then
           write(stdo,*) 'INTERNAL ERROR in thermal boundaries.'
           stop
        endif
        ray_cart_diry = dir_planex
        ray_cart_dirz = dir_planey
        if(bc_ilr.eq.1) then
           epsdist       = 1d-4*(amr_grid_xi(2,1)-amr_grid_xi(1,1))/ &
                           (2**amr_levelmax)
           ray_cart_x    = amr_grid_xi(1,1)+epsdist
           ray_cart_dirx = dir_perp
        else
           epsdist       = 1d-4*(amr_grid_xi(amr_grid_nx,1)-amr_grid_xi(amr_grid_nx-1,1))/ &
                           (2**amr_levelmax)
           ray_cart_x    = amr_grid_xi(amr_grid_nx+1,1)-epsdist
           ray_cart_dirx = -dir_perp
        endif
        rn = ran2(iseed)
        ray_cart_y = amr_grid_xi(1,2) + rn * &
             ( amr_grid_xi(amr_grid_ny,2) - amr_grid_xi(1,2) )
        rn = ran2(iseed)
        ray_cart_z = amr_grid_xi(1,3) + rn * &
             ( amr_grid_xi(amr_grid_nz,3) - amr_grid_xi(1,3) )
     elseif(bc_idir.eq.2) then
        if(igrid_coord.eq.10) then
           write(stdo,*) 'INTERNAL ERROR in thermal boundaries.'
           stop
        endif
        ray_cart_dirx = dir_planex
        ray_cart_dirz = dir_planey
        if(bc_ilr.eq.1) then
           epsdist       = 1d-4*(amr_grid_xi(2,2)-amr_grid_xi(1,2))/ &
                           (2**amr_levelmax)
           ray_cart_y    = amr_grid_xi(1,2)+epsdist
           ray_cart_diry = dir_perp
        else
           epsdist       = 1d-4*(amr_grid_xi(amr_grid_ny,2)-amr_grid_xi(amr_grid_ny-1,2))/ &
                           (2**amr_levelmax)
           ray_cart_y    = amr_grid_xi(amr_grid_ny+1,2)-epsdist
           ray_cart_diry = -dir_perp
        endif
        rn = ran2(iseed)
        if(igrid_coord.eq.10) then
           ray_cart_x = 0.d0
        else
           ray_cart_x = amr_grid_xi(1,1) + rn * &
                ( amr_grid_xi(amr_grid_nx,1) - amr_grid_xi(1,1) )
        endif
        rn = ran2(iseed)
        ray_cart_z = amr_grid_xi(1,3) + rn * &
             ( amr_grid_xi(amr_grid_nz,3) - amr_grid_xi(1,3) )
     else
        ray_cart_dirx = dir_planex
        ray_cart_diry = dir_planey
        if(bc_ilr.eq.1) then
           epsdist       = 1d-4*(amr_grid_xi(2,3)-amr_grid_xi(1,3))/ &
                           (2**amr_levelmax)
           ray_cart_z    = amr_grid_xi(1,3)+epsdist
           ray_cart_dirz = dir_perp
        else
           epsdist       = 1d-4*(amr_grid_xi(amr_grid_nz,3)-amr_grid_xi(amr_grid_nz-1,3))/ &
                           (2**amr_levelmax)
           ray_cart_z    = amr_grid_xi(amr_grid_nz+1,3)-epsdist
           ray_cart_dirz = -dir_perp
        endif
        rn = ran2(iseed)
        if((igrid_coord.eq.10).or.(igrid_coord.eq.20)) then
           ray_cart_x = 0.d0
        else
           ray_cart_x = amr_grid_xi(1,1) + rn * &
                ( amr_grid_xi(amr_grid_nx,1) - amr_grid_xi(1,1) )
        endif
        rn = ran2(iseed)
        if(igrid_coord.eq.10) then
           ray_cart_y = 0.d0
        else
           ray_cart_y = amr_grid_xi(1,2) + rn * &
             ( amr_grid_xi(amr_grid_ny,2) - amr_grid_xi(1,2) )
        endif
     endif
     !
     ! Now let the AMR module find the cell index
     !
     if(amr_tree_present) then
        call amr_findcell(ray_cart_x,ray_cart_y,ray_cart_z,acell)
        ray_index = acell%leafindex
     else
        call amr_findbasecell(ray_cart_x,ray_cart_y,ray_cart_z,ix,iy,iz)
        ray_index = ix+((iy-1)+(iz-1)*amr_grid_ny)*amr_grid_nx
     endif
     !
     ! If weighted photon packages, then get the energy
     !
     if(mc_weighted_photons) energy = mc_energy_bc
     !
  elseif(rn.ge.mc_cumlum3) then
     !
     ! ----- We have a photon emitted by the continuous stellar source -----
     !       (for simulations of galaxies and such)
     !
     ! First determine which cell to emit from
     !
     rn = ran2(iseed)
     call hunt(stellarsrc_lumcum,nrcells+1,rn,icell)
     ray_index = cellindex(icell)
     !
     ! Find the position within the cell from where to emit
     !
     call montecarlo_random_pos_in_cell(icell)
     !
     ! Find the direction in which to emit
     !
     call montecarlo_randomdir(ray_cart_dirx,ray_cart_diry,ray_cart_dirz)
     !
     ! Now determine from which stellar population we will send the photon
     !
     if(stellarsrc_nrtemplates.eq.1) then
        itemplate = 1
     else
        dum = 0.d0
        do itemplate=1,stellarsrc_nrtemplates
           stellarsrc_cum(itemplate) = dum
           dum = dum + stellarsrc_dens(itemplate,ray_index) *           &
                       mc_stellarsrc_templates(inu,itemplate)
        enddo
        dum = 1.d0 / dum
        do itemplate=1,stellarsrc_nrtemplates
           stellarsrc_cum(itemplate) = stellarsrc_cum(itemplate) * dum
        enddo
        stellarsrc_cum(stellarsrc_nrtemplates+1) = 1.d0
        rn = ran2(iseed)
        call hunt(stellarsrc_cum,stellarsrc_nrtemplates+1,rn,itemplate)
        if(debug_check_all.eq.1) then
           if((itemplate.lt.1).or.(itemplate.gt.stellarsrc_nrtemplates)) then
              write(stdo,*) 'INTERNAL ERROR: Itemplate out of range'
              stop
           endif
        endif
     endif
     !
     ! If weighted photon packages, then get the energy
     !
     if(mc_weighted_photons) energy = mc_energy_stellarsrc
     !
  elseif(rn.ge.mc_cumlum2) then
     !
     ! ----- We have a photon emitted by the external environment into -----
     !       our computational domain. This is done from a sphere 
     !       centered around the center of the coordinate system.
     !
     ray_index = 0
     !
     ! Find a random point on the sphere
     !
     call montecarlo_randomdir(pdirx,pdiry,pdirz)
     ray_cart_x  = grid_contsph_x + pdirx*grid_contsph_r
     ray_cart_y  = grid_contsph_y + pdiry*grid_contsph_r
     ray_cart_z  = grid_contsph_z + pdirz*grid_contsph_r
     !
     ! If mirror symmetry, then enforce this here
     !
     if(igrid_mirror.eq.1) then
        if(ray_cart_z.lt.0.d0) then
           ray_cart_z = -ray_cart_z
        endif
     endif
     !
     ! Now determine a random direction for the photon. Here we again
     ! have the angle difficulty we encountered for emission from a
     ! stellar surface. See above by sqrt(rn) is used.
     ! 
     rn            = ran2(iseed)
     ray_cart_dirx = sqrt(rn)
     dum           = twopi*ran2(iseed)
     rn            = sqrt(1.d0-rn)
     ray_cart_diry = rn*cos(dum)
     ray_cart_dirz = rn*sin(dum)
     call montecarlo_rotatevec(ray_cart_dirx,ray_cart_diry,ray_cart_dirz, &
                               -pdirx,-pdiry,-pdirz)
     !
     ! If weighted photon packages, then get the energy
     !
     if(mc_weighted_photons) energy = mc_energy_extlum
     !
  elseif(rn.ge.mc_cumlum1) then
     !
     ! ----- We have a photon emitted by quantum-heated grains -----
     !
     ! First determine which cell to emit from
     !
     rn = ran2(iseed)
     call hunt(emisquant_cum,nrcells+1,rn,icell)
     ray_index = cellindex(icell)
     !
     ! Find the position within the cell from where to emit
     !
     call montecarlo_random_pos_in_cell(icell)
     !
     ! Find the direction in which to emit
     !
     call montecarlo_randomdir(ray_cart_dirx,ray_cart_diry,ray_cart_dirz)
     !
     ! If weighted photon packages, then get the energy
     !
     if(mc_weighted_photons) energy = mc_energy_quant
     !
  else
     !
     ! ----- We have a photon emitted by thermal emission -----
     !
     ! First determine which cell to emit from
     !
     rn = ran2(iseed)
     call hunt(mc_cumulthermemis,nrcells+1,rn,icell)
     ray_index = cellindex(icell)
     !
     ! Find the position within the cell from where to emit
     !
     call montecarlo_random_pos_in_cell(icell)
     !
     ! If weighted photon packages, then get the energy
     !
     if(mc_weighted_photons) energy = mc_energy_thermal
     !
     ! Find the direction in which to emit
     !
     if(alignment_mode.le.0) then
        !
        ! Simple case: just random direction
        !
        call montecarlo_randomdir(ray_cart_dirx,ray_cart_diry,ray_cart_dirz)
     else
        !
        ! If aligned grains, then we must choose the direction from the
        ! angular probability distribution function for the thermal
        ! emission from the aligned grain. Also the polarization state 
        ! must be randomly chosen from the appropriate distribution function.
        !
        ener = energy
        call montecarlo_aligned_randomphot(ray_index,inu,ener,photpkg)
        ray_cart_dirx = photpkg%n(1)
        ray_cart_diry = photpkg%n(2)
        ray_cart_dirz = photpkg%n(3)
        !
        ! Signal that we already set the photpkg, so that it is not
        ! overridden.
        !
        todo_photpkg = .false.
        !
     endif
     !
  endif
  !
  ! Set the local energy to the starting value
  !
  ener = energy
  !
  ! Now do some checks...
  !
  if(debug_check_all.eq.1) then
     if((ray_index.lt.1).or.(ray_index.gt.nrcellsmax)) then
        write(stdo,*) 'INTERNAL ERROR: index out of range: ray_index = ',ray_index
        stop
     endif
     if(abs(ray_cart_dirx**2+ray_cart_diry**2+ray_cart_dirz**2-1.d0).gt.1d-14) then
        write(stdo,*) 'INTERNAL ERROR: dir vector not length unity: dir = ', &
               ray_cart_dirx,ray_cart_diry,ray_cart_dirz
        stop
     endif
  endif
  !
  ! Now associate pointer to cell.
  ! This is grid type dependent.
  !
  if(igrid_type.lt.100) then
     !
     ! AMR type grid
     !
     if(ray_index.gt.0) then
        !
        ! We are in a cell
        !
        if(amr_tree_present) then
           amrray_cell => amr_index_to_leaf(ray_index)%link
        else
           call amr_regular_get_ixyz(ray_index,amrray_ix_curr,amrray_iy_curr,amrray_iz_curr)
        endif
     else
        nullify(amrray_cell)
        amrray_ix_curr = -1
        amrray_iy_curr = -1
        amrray_iz_curr = -1
     endif
  else
     write(stdo,*) 'ERROR: This grid type not yet implemented'
     stop
  endif
  !
  ! In case we use polarization, init the photpkg structure.
  ! We assume the photon to start non-polarized.
  ! Exception: the thermal dust emission if the grains 
  ! are aligned. If that is the case, the todo_photpkg
  ! is switched to .false. because then it is initialized
  ! above. 
  !
  if((scattering_mode.ge.5).and.todo_photpkg) then
     photpkg%E    = ener
     photpkg%Q    = 0.d0
     photpkg%U    = 0.d0
     photpkg%V    = 0.d0
     photpkg%n(1) = ray_cart_dirx
     photpkg%n(2) = ray_cart_diry
     photpkg%n(3) = ray_cart_dirz
     call polarization_make_s_vector(photpkg%n,photpkg%s)
  endif
  !
  ! Now start the random walk
  !
  ok = .true.
  iscatevent = 0
  selectscat_iscat = 1
!!  ieventcount = 0         ! For debugging
  do while(ok)
     !
     ! Get a random tau, this time it is the scattering tau only
     !
     taupath = -log(1.d0-ran2(iseed))
     !
     ! Move photon to next scattering event
     ! 
     call walk_cells_scat(params,taupath,ener,inu,arrived,ispec,ierror)
     !
     ! If arrived at end-point or escaped to infinity, then return
     !
     if(arrived) then
        return
     endif
     !
     ! If error, then also return
     !
     if(ierror.ne.0) then
        stop 1142
        return
     endif
     !
     ! If photon destroyed (by e.g. quantum heated grains which are
     ! dealt with elsewhere), then also return. 
     !
     ! NOTE: In this MC mode a photon package is also destroyed if
     !       the package has extincted more than the threshold
     !       given by mc_scat_maxtauabs
     !
     if(mc_photon_destroyed) then
        return
     endif
     !
     ! Increase the scattering event counter. This is only important
     ! for the case when we wish to limit the number of scattering
     ! events treated.
     !
     iscatevent = iscatevent + 1
     !
     ! Now check if we have reached the max nr of scattering events
     ! we allow (normally this is set to mc_max_nr_scat_events=-1,
     ! which means infinite number of scatterings allowed).
     !
     if(iscatevent.eq.mc_max_nr_scat_events) then
        return
     endif
     !
     ! Scattering event
     !
     if(scattering_mode.eq.1) then
        !
        ! Isotropic scattering
        ! 
        call montecarlo_randomdir(ray_cart_dirx,ray_cart_diry,ray_cart_dirz)
        !
     elseif(scattering_mode.eq.2) then
        !
        ! Anisotropic scattering
        !
        ! The walk_cells_scat() routine already determined with which
        ! dust species the scattering took place: ispec.
        !
        g     = kappa_g(inu,ispec)
        call henyeygreenstein_randomdir(g,dirnewx,dirnewy,dirnewz)
        call montecarlo_rotatevec(dirnewx,dirnewy,dirnewz,   &
             ray_cart_dirx,ray_cart_diry,ray_cart_dirz)
        ray_cart_dirx = dirnewx
        ray_cart_diry = dirnewy
        ray_cart_dirz = dirnewz
        !
     elseif((scattering_mode.eq.3).or.(scattering_mode.eq.4)) then
        !
        ! Anisotropic scattering off grain for dust species ispec,
        ! following the tabulated Z matrix. 
        !
        call anisoscat_randomdir(inu,ispec,dirnewx,dirnewy,dirnewz)
        call montecarlo_rotatevec(dirnewx,dirnewy,dirnewz,   &
             ray_cart_dirx,ray_cart_diry,ray_cart_dirz)
        ray_cart_dirx = dirnewx
        ray_cart_diry = dirnewy
        ray_cart_dirz = dirnewz
        !
     elseif(scattering_mode.eq.5) then
        !
        ! Full polarized scattering
        !
        ! Note: This subroutine will also modify the polarization state
        ! of photpkg. Also note that this subroutine does not need any
        ! transformation to and from the particle frame as was necessary
        ! in the scattering_mode<5 cases.
        !
        call polarization_randomorient_mc_scatter(photpkg,  &
             scat_munr,scat_mui_grid,scat_thetai_grid,      &
             zmatrix(:,inu,1,ispec),                        &
             zmatrix(:,inu,2,ispec),                        &
             zmatrix(:,inu,3,ispec),                        &
             zmatrix(:,inu,4,ispec),                        &
             zmatrix(:,inu,5,ispec),                        &
             zmatrix(:,inu,6,ispec),                        &
             zcumul(:,inu,1,ispec),                         &
             zcumul(:,inu,2,ispec),                         &
             .false.)
        ray_cart_dirx = photpkg%n(1)
        ray_cart_diry = photpkg%n(2)
        ray_cart_dirz = photpkg%n(3)
        !
     else
        write(stdo,*) 'ERROR: scattering mode ',scattering_mode,' not known'
        stop
     endif
     !
     ! For debugging: Increase event counter
     !
!!     ieventcount = ieventcount + 1
     !$omp atomic
     ieventcounttot = ieventcounttot + 1
     !
     ! For selectscat: Increase counter
     !
     selectscat_iscat = selectscat_iscat + 1
     !
     ! Next event
     !
  end do
  !
end subroutine walk_full_path_scat


!--------------------------------------------------------------------------
!       FOLLOW A PHOTON ALONG A RAY FOR THE THERMAL MONTE CARLO
!
!     This routine follows a photon until the next scattering
!     event. In the mean time it leaves its energy behind by
!     continuous absorption in the cells it crosses. 
!--------------------------------------------------------------------------
subroutine walk_cells_thermal(params,taupath,iqactive,arrived, & 
                              therm,ispecc,ierror)
  implicit none
  type(mc_params) :: params
  integer :: iqactive,ispec,ierror,ispecc,idir
  doubleprecision :: taupath,fr
  doubleprecision :: tau,dtau,albedo,absorb,dum,scatsrc0,costheta
  doubleprecision :: dexp,dener,ds,dummy,alpha_tot,g,rn,dss,src4(1:4)
  logical :: ok,arrived,therm
  doubleprecision :: prev_x,prev_y,prev_z
  !$ logical::continue
  !
  ! Reset
  !
  tau      = 0.d0
  ierror   = 0
  !
  ! Default
  !
  mc_photon_destroyed = .false.
  arrived  = .false.
  !
  ! In this subroutine we will not make use of the amrray option to
  ! advance only partly within a cell. We will check here if we have
  ! an event within this cell and compute here the location of this
  ! event. Therefore the ray_dsend is set to 1d99.
  !
  ray_dsend = 1d99
  !
  ! Now do the loop
  !
  ok = .true.
  do while(ok)
     !
     ! OpenMP Parallellization: Lock this cell (and only continue when succesfully locked)
     !
     !$ continue=.true.
     !$ if(ray_index .gt. 0)then
     !$    do while(.NOT. omp_test_lock(lock(ray_index)))
     !$       if(continue)then
     !$omp atomic
     !$          conflict_counter=conflict_counter+1
     !$          continue=.false.
     !$       end if
     !$    end do
     !$ end if
     !
     ! For cell photon statistics, increase the counter
     !
     if(params%debug_write_stats.ne.0) then
        if(ray_index.ge.1) then
           if(mc_iphotcurr.gt.mc_ilastphot(ray_index)) then
              mc_iphotcount(ray_index) = mc_iphotcount(ray_index) + 1
              mc_ilastphot(ray_index)  = mc_iphotcurr
           endif
        endif
     endif
     !
     ! Back up the current position
     !
     prev_x     = ray_cart_x
     prev_y     = ray_cart_y
     prev_z     = ray_cart_z
     !
     ! Find new position
     !
     if(igrid_type.lt.100) then 
        !
        ! We use a regular or AMR grid
        !
        if(igrid_coord.lt.100) then
           !
           ! We use cartesian coordinates
           !
           call amrray_find_next_location_cart(ray_dsend,            &
                ray_cart_x,ray_cart_y,ray_cart_z,                    &
                ray_cart_dirx,ray_cart_diry,ray_cart_dirz,           &
                ray_index,ray_indexnext,ray_ds,arrived)
        elseif(igrid_coord.lt.200) then
           !
           ! We use spherical coordinates
           !
           call amrray_find_next_location_spher(ray_dsend,           &
                ray_cart_x,ray_cart_y,ray_cart_z,                    &
                ray_cart_dirx,ray_cart_diry,ray_cart_dirz,           &
                ray_index,ray_indexnext,ray_ds,arrived)
        else
           write(stdo,*) 'ERROR: Cylindrical coordinates not yet implemented'
           stop
        endif
        !
        ! In case stars are treated as spheres, copy this information
        ! Note that later we will check if the photon may have reached
        ! its end point (an absorption/scattering event) before reaching
        ! the sphere or a grid wall. In that case, mc_photon_destroyed
        ! will be reset to .false., see below.
        !
        if(amrray_ispherehit.lt.0) then
           mc_photon_destroyed = .true.
           arrived             = .true.    ! To end this segment and return
        else
           mc_photon_destroyed = .false.
        endif
        !
     else
        write(stdo,*) 'SORRY: Delaunay or Voronoi grids not yet implemented'
        stop
     endif
     !
     ! Path length
     ! 
     ! ######### CHECK: WHY NOT USE ray_ds ??? ########
     !
     ds   = sqrt( (ray_cart_x-prev_x)**2 + (ray_cart_y-prev_y)**2 + (ray_cart_z-prev_z)**2 )
     !
     ! Compute the alpha's, dtau, albedo and absorb
     !
     ! NOTE: For non-Cartesian coordinates it can happen that we have not
     !       yet arrived, but are currently not in a cell either. Hence 
     !       we need to check if ray_index.ge.1
     !
     if(ray_index.ge.1) then
        alpha_a_tot = 0.d0
        alpha_s_tot = 0.d0
        do ispec=1,dust_nr_species
           alpha_a(ispec) = dustdens(ispec,ray_index) * kappa_a(ray_inu,ispec)
           alpha_s(ispec) = dustdens(ispec,ray_index) * kappa_s(ray_inu,ispec)
           alpha_a_tot    = alpha_a_tot + alpha_a(ispec)
           alpha_s_tot    = alpha_s_tot + alpha_s(ispec)
        enddo
        alpha_tot = alpha_a_tot + alpha_s_tot
        albedo    = alpha_s_tot / alpha_tot
        dtau      = alpha_tot * ds
     else
        do ispec=1,dust_nr_species
           alpha_a(ispec) = 0.d0
           alpha_s(ispec) = 0.d0
        enddo
        alpha_a_tot = 0.d0
        alpha_s_tot = 0.d0
        alpha_tot   = 0.d0
        albedo      = 0.d0
        dtau        = 0.d0
     endif
     absorb = 1.d0 - albedo
     !
     ! Check if we already reached the end point
     !
     if(tau+dtau.ge.taupath) then
        !
        ! Yes: reached the end point, i.e. reached new scattering or
        ! absorption event
        !
        ! Check...
        !
        if(debug_check_all.eq.1) then
           if(ray_index.lt.1) stop 94001
           if(dtau.eq.0.d0) stop 94002
           if(alpha_a_tot.eq.0.d0) stop 94003
        endif
        !
        ! Compute the fraction of the segment we have advanced
        !
        fr     = (taupath-tau)/dtau
        dss    = fr * ds
        !
        ! Add energy to cell
        !
        dum = absorb * (taupath-tau) * energy / alpha_a_tot
        if(iqactive.le.0) then
           do ispec=1,dust_nr_species
              mc_cumulener(ispec,ray_index) = mc_cumulener(ispec,ray_index) +  &
                   dum * alpha_a(ispec) 
           enddo
        else
           do ispec=1,dust_nr_species
              if(dust_quantum(ispec).eq.0) then 
                 mc_cumulener(ispec,ray_index) = mc_cumulener(ispec,ray_index) +  &
                      dum * alpha_a(ispec) 
              endif
           enddo
        endif
        !
        ! Add photons to mean intensity of 
        ! primary (quantum-heating) photons
        !
        if(iqactive.ne.0) then
           miquant(ray_inu,ray_index) = miquant(ray_inu,ray_index) +               &
                   (taupath-tau) * energy /                                        &
                   ( cellvolume(ray_index) * freq_dnu(ray_inu) *                   &
                     alpha_tot * 12.566371d0 )
        endif
        !
        ! Compute the location of this point
        !
        ray_cart_x = prev_x + fr * ( ray_cart_x - prev_x )
        ray_cart_y = prev_y + fr * ( ray_cart_y - prev_y )
        ray_cart_z = prev_z + fr * ( ray_cart_z - prev_z )
        !
        ! Here we determine if it was an absorption or a scattering
        ! event. 
        !
        rn    = ran2(iseed)
        therm = (rn.ge.albedo)
        if(therm) then
           !
           ! It was absorption. We do not need to know ispecc. Put it to zero, for safety.
           !
           ispecc = 0
        else
           !
           ! It was a scattering event.
           !
           if(scattering_mode.gt.1) then
              !
              ! For anisotropic scattering we need to know off which dust species
              ! the photon has scattered. Determine this here, and call it ispecc.
              !
              alphacum(1) = 0.d0
              do ispec=1,dust_nr_species
                 alphacum(ispec+1) = alphacum(ispec) + alpha_s(ispec)
              enddo
              if(alphacum(dust_nr_species+1).le.0.d0) stop 2244
              do ispec=1,dust_nr_species
                 alphacum(ispec) = alphacum(ispec) / alphacum(dust_nr_species+1)
              enddo
              alphacum(dust_nr_species+1) = 1.d0
              rn = ran2(iseed)
              call hunt(alphacum,dust_nr_species+1,rn,ispecc)
              if((ispecc.lt.1).or.(ispecc.gt.dust_nr_species)) stop 7266
           else
              rn = ran2(iseed)   ! NOTE: This is not necessary, but is just
                                 !       there to pass the selftest. 
                                 !       See notes 15.11.09 in Radmc_3D_LOG.txt
              ispecc=0    ! Force error (with array bound check) if accidently used
           endif
        endif
        !
        ! Done with this segment of the random walk, so return and 
        ! handle the scattering or absorption event
        !
        arrived = .false.
        mc_photon_destroyed = .false.
        !
        ! OpenMP Parallellization: Release lock on this cell
        !
        !$ if(ray_index .ge. 1 )then
        !$    call omp_unset_lock(lock(ray_index));
        !$ endif
        return
        !
     endif
     !
     ! Endpoint of this segment.
     !
     if(ray_index.ge.1) then
        !
        ! We are in a cell. So add energy to cell
        !
        dum = absorb * dtau * energy / alpha_a_tot
        do ispec=1,dust_nr_species
           if((dust_quantum(ispec).eq.0).or.(iqactive.le.0)) then 
              mc_cumulener(ispec,ray_index) =                          &
                   mc_cumulener(ispec,ray_index) + dum * alpha_a(ispec) 
           endif
        enddo
        !
        ! Add photons to mean intensity of 
        ! primary (quantum-heating) photons
        !             
        if(iqactive.ne.0) then
           miquant(ray_inu,ray_index) = miquant(ray_inu,ray_index) +   &
                   dtau * energy /                                     &
                   ( cellvolume(ray_index) * freq_dnu(ray_inu) *       &
                     alpha_tot * 12.566371d0 )
        endif
     endif
     !
     ! Increase the tau
     !
     tau = tau + dtau
     !
     ! OpenMP Parallellization: Release lock on this cell
     !
     !$ if(ray_index .ge. 1 )then
     !$    call omp_unset_lock(lock(ray_index));
     !$ endif
     !
     ! Now make next cell the current cell
     !
     ray_index = ray_indexnext
     !
     ! Do the same for the pointers.
     ! This is grid type dependent.
     !
     if(igrid_type.lt.100) then
        !
        ! AMR type grid
        !
        if(ray_index.gt.0) then
           !
           ! We are in a cell
           !
           if(amr_tree_present) then
              amrray_cell => amr_index_to_leaf(ray_index)%link
           else
              call amr_regular_get_ixyz(ray_index,amrray_ix_curr,amrray_iy_curr,amrray_iz_curr)
           endif
        else
           nullify(amrray_cell)
           amrray_ix_curr = -1
           amrray_iy_curr = -1
           amrray_iz_curr = -1
        endif
     else
        write(stdo,*) 'ERROR: This grid type not yet implemented'
        stop
     endif
     !
     ! Now, if this new/next segment goes off to infinity, then
     ! we are done.
     !
     if(arrived) then
        !
        ! We escaped to infinity. 
        !
        ! OpenMP Parallellization: Lock is already released.
        ! BUGFIX 23.02.2017
        !
        return
     endif
     !
  enddo
  !
end subroutine walk_cells_thermal


!--------------------------------------------------------------------------
!       FOLLOW A PHOTON ALONG A RAY FOR THE SCATTERING MONTE CARLO
!
!     This routine follows a photon until the next scattering
!     event. In the mean time it leaves its energy behind by
!     continuous absorption in the cells it crosses. 
!
! NOTE: inu in the argument list is the position in the mcscat_inus()
!       array. The 'real' inu is, of course, mcscat_inus(inu) or 
!       equivalently ray_inu.
!--------------------------------------------------------------------------
subroutine walk_cells_scat(params,taupath,ener,inu,arrived,ispecc,ierror)
  implicit none
  type(mc_params) :: params
  integer :: ispec,ierror,inu,ispecc,idir,iddr
  doubleprecision :: taupath,fr,ener,enerav
  doubleprecision :: tau,dtau,dum,dtauabs,dtauscat,xptauabs,xxtauabs
  doubleprecision :: ds,rn,scatsrc0,mnint,dss
  logical :: ok,arrived
  doubleprecision :: prev_x,prev_y,prev_z
  doubleprecision :: costheta,g,phasefunc,dummy,src4(1:4)
  doubleprecision :: axi(1:2,1:3),Ebk,Qbk,Ubk,Vbk
  doubleprecision :: nvec_orig(1:3),svec_orig(1:3),levent
  doubleprecision :: xbk,ybk,cosphievent,sinphievent
  doubleprecision :: deltaphi,cosdphi,sindphi
  integer :: idirs
  !$ logical::continue
  !
  ! Reset
  !
  tau      = 0.d0
  ierror   = 0
  !
  ! Default
  !
  mc_photon_destroyed = .false.
  arrived  = .false.
  ray_inu  = inu
  !
  ! In this subroutine we will not make use of the amrray option to
  ! advance only partly within a cell. We will check here if we have
  ! an event within this cell and compute here the location of this
  ! event. Therefore the ray_dsend is set to 1d99.
  !
  ray_dsend = 1d99
  !
  ! If we have scattering included (sort of logical we have it included),
  ! then we can pre-compute the phase functions for all dust species.
  !
  ! Note: if we use local observer, then we cannot pre-compute the 
  !       phase function for non-isotropic scattering.
  !
  if(scattering_mode.eq.2) then
     !
     ! Anisotropic scattering using the Henyey-Greenstein function
     !
     if(.not.mcscat_localobserver) then
        do idir=1,mcscat_nrdirs
           costheta = ray_cart_dirx*mcscat_dirs(1,idir) +              &
                      ray_cart_diry*mcscat_dirs(2,idir) +              &
                      ray_cart_dirz*mcscat_dirs(3,idir)
           do ispec=1,dust_nr_species
              g                            = kappa_g(ray_inu,ispec)
              mcscat_phasefunc(idir,ispec) = henyeygreenstein_phasefunc(g,costheta)
           enddo
        enddo
     else
        write(stdo,*) 'ERROR: For now, anisotropic scattering not allowed'
        write(stdo,*) '       in local observer mode.'
        stop 453
     endif
  elseif(scattering_mode.eq.3) then
     !
     ! Anisotropic scattering using tabulated phase function (the
     ! Z11 matrix elements)
     !
     if(.not.mcscat_localobserver) then
        do idir=1,mcscat_nrdirs
           costheta = ray_cart_dirx*mcscat_dirs(1,idir) +              &
                      ray_cart_diry*mcscat_dirs(2,idir) +              &
                      ray_cart_dirz*mcscat_dirs(3,idir)
           do ispec=1,dust_nr_species
              mcscat_phasefunc(idir,ispec) = anisoscat_phasefunc(ray_inu,ispec,costheta)
           enddo
        enddo
     else
        write(stdo,*) 'ERROR: For now, anisotropic scattering not allowed'
        write(stdo,*) '       in local observer mode.'
        stop 453
     endif
  endif
  !
  ! Now do the loop
  !
  ok = .true.
  do while(ok)
     !
     ! OpenMP Parallellization: Lock this cell (and only continue when succesfully locked)
     !
     !$ continue=.true.
     !$ if(ray_index .gt. 0)then
     !$    do while(.NOT. omp_test_lock(lock(ray_index)))
     !$       if(continue)then
     !$omp atomic
     !$          conflict_counter=conflict_counter+1
     !$          continue=.false.
     !$       end if
     !$    end do
     !$ end if
     !
     ! Do some safety checks
     !
     if(debug_check_all.eq.1) then
        if(ray_index.ge.1) then
           if(igrid_type.lt.100) then
              !
              ! AMR grid
              !
              ! Get the cell walls
              !
              if(amr_tree_present) then
                 if(.not.associated(amrray_cell)) stop 7270
                 do iddr=1,3
                    axi(1,iddr) = amr_finegrid_xi(amrray_cell%ixyzf(iddr),iddr,amrray_cell%level)
                    axi(2,iddr) = amr_finegrid_xi(amrray_cell%ixyzf(iddr)+1,iddr,amrray_cell%level)
                 enddo
              else
                 if(amrray_ix_curr.le.0) stop 7270
                 axi(1,1) = amr_finegrid_xi(amrray_ix_curr,1,0)
                 axi(2,1) = amr_finegrid_xi(amrray_ix_curr+1,1,0)
                 axi(1,2) = amr_finegrid_xi(amrray_iy_curr,2,0)
                 axi(2,2) = amr_finegrid_xi(amrray_iy_curr+1,2,0)
                 axi(1,3) = amr_finegrid_xi(amrray_iz_curr,3,0)
                 axi(2,3) = amr_finegrid_xi(amrray_iz_curr+1,3,0)
              endif
              !
              if(igrid_coord.lt.100) then
                 !
                 ! Cartesian coordinates
                 !
                 if((ray_cart_x.lt.axi(1,1)).or.(ray_cart_x.gt.axi(2,1)).or. &
                    (ray_cart_y.lt.axi(1,2)).or.(ray_cart_y.gt.axi(2,2)).or. &
                    (ray_cart_z.lt.axi(1,3)).or.(ray_cart_z.gt.axi(2,3))) then
                    write(stdo,*) 'INTERNAL ERROR: Point outside of cell in montecarlo module'
                    write(stdo,*) ray_cart_x,axi(1,1),axi(2,1)
                    write(stdo,*) ray_cart_y,axi(1,2),axi(2,2)
                    write(stdo,*) ray_cart_z,axi(1,3),axi(2,3)
                    stop
                 endif
              else
                 !
                 ! These coordiantes not yet active
                 !
                 write(stdo,*) 'INTERNAL ERROR: Non-Cart coordinates not yet implemented'
                 stop 
              endif
           else
              !
              ! Unstructured grids not yet active
              !
              write(stdo,*) 'INTERNAL ERROR: Unstructured grids not yet implemented'
              stop            
           endif
        endif
     endif
     !
     ! For cell photon statistics, increase the counter
     !
     if(params%debug_write_stats.ne.0) then
        if(ray_index.ge.1) then
           if(mc_iphotcurr.gt.mc_ilastphot(ray_index)) then
              mc_iphotcount(ray_index) = mc_iphotcount(ray_index) + 1
              mc_ilastphot(ray_index)  = mc_iphotcurr
           endif
        endif
     endif
     !
     ! Back up the current position
     !
     prev_x     = ray_cart_x
     prev_y     = ray_cart_y
     prev_z     = ray_cart_z
     !
     ! Find new position
     !
     if(igrid_type.lt.100) then 
        !
        ! We use a regular or AMR grid
        !
        if(igrid_coord.lt.100) then
           !
           ! We use cartesian coordinates
           !
           call amrray_find_next_location_cart(ray_dsend,            &
                ray_cart_x,ray_cart_y,ray_cart_z,                    &
                ray_cart_dirx,ray_cart_diry,ray_cart_dirz,           &
                ray_index,ray_indexnext,ray_ds,arrived)
        elseif(igrid_coord.lt.200) then
           !
           ! We use spherical coordinates
           !
           call amrray_find_next_location_spher(ray_dsend,           &
                ray_cart_x,ray_cart_y,ray_cart_z,                    &
                ray_cart_dirx,ray_cart_diry,ray_cart_dirz,           &
                ray_index,ray_indexnext,ray_ds,arrived)
        else
           write(stdo,*) 'ERROR: Cylindrical coordinates not yet implemented'
           stop
        endif
        !
        ! In case stars are treated as spheres, copy this information
        ! Note that later we will check if the photon may have reached
        ! its end point (an absorption/scattering event) before reaching
        ! the sphere or a grid wall. In that case, mc_photon_destroyed
        ! will be reset to .false., see below.
        !
        if(amrray_ispherehit.lt.0) then
           mc_photon_destroyed = .true.
           arrived             = .true.    ! To end this segment and return
        else
           mc_photon_destroyed = .false.
        endif
     else
        write(stdo,*) 'SORRY: Delaunay or Voronoi grids not yet implemented'
        stop
     endif
     !
     ! Path length
     ! 
     ! ######### CHECK: WHY NOT USE ray_ds ??? ########
     !
     ds   = sqrt( (ray_cart_x-prev_x)**2 + (ray_cart_y-prev_y)**2 + (ray_cart_z-prev_z)**2 )
     !
     ! Compute the alpha's, dtau
     !
     ! NOTE: For non-Cartesian coordinates it can happen that we have not
     !       yet arrived, but are currently not in a cell either. Hence 
     !       we need to check if ray_index.ge.1
     !
     ! NOTE: Since we are here interested in Monte Carlo in terms of 
     !       scattering only (treating absorption simply as a loss term),
     !       we define dtau as only the scattering dtau.
     !
     ! NOTE: Since we must be able to use this subroutine also without 
     !       scattering, e.g. if we use radmc3d mcmono, we add 1d-99 to
     !       all opacities to make sure that they are non-zero.
     !
     if(ray_index.ge.1) then
        alpha_a_tot = 0.d0
        alpha_s_tot = 0.d0
        do ispec=1,dust_nr_species
           alpha_a(ispec) = dustdens(ispec,ray_index) * &
                            kappa_a(ray_inu,ispec) + 1d-99
           alpha_s(ispec) = dustdens(ispec,ray_index) * &
                            kappa_s(ray_inu,ispec) + 1d-99
           alpha_a_tot    = alpha_a_tot + alpha_a(ispec)
           alpha_s_tot    = alpha_s_tot + alpha_s(ispec)
        enddo
        dtau      = alpha_s_tot * ds
     else
        do ispec=1,dust_nr_species
           alpha_a(ispec) = 0.d0
           alpha_s(ispec) = 0.d0
        enddo
        alpha_a_tot = 0.d0
        alpha_s_tot = 0.d0
        dtau        = 0.d0
     endif
     !
     ! Check if we already reached the end point
     !
     if(tau+dtau.ge.taupath) then
        !
        ! Yes: reached the end point, i.e. reached new scattering event
        !
        ! Check...
        !
        if(debug_check_all.eq.1) then
           if(ray_index.lt.1) stop 94001
        endif
        !
        ! Compute dtauabs and dtauscat as well as exp(-dtauabs) and
        ! (1-exp(-dtauabs))/dtauabs.
        !
        dtauscat = taupath-tau
        dtauabs  = dtauscat * alpha_a_tot / alpha_s_tot
        xptauabs = exp(-dtauabs)
        if(dtauabs.gt.1d-6) then
           xxtauabs = (1.d0-xptauabs)/dtauabs
        else
           xxtauabs = 1.d0-0.5d0*dtauabs
        endif
        !
        ! Compute enerav
        !
        enerav  = ener * xxtauabs
        !
        ! For the selectscat analysis stuff (normally not important)
        !
        if((selectscat_iscat.ge.selectscat_iscat_first).and.(selectscat_iscat.le.selectscat_iscat_last)) then
        !
        ! Add photons to scattering source
        !
        ! NOTE: In the original RADMC I did not divide by 4*pi yet; only
        !       in the writing to scatsource.dat I did that. Here I do 
        !       it immediately.
        !
        ! Note: Since we follow only the scattering events, the tau
        !       and taupath are scattering tau only. So no `albedo'
        !       is necessary here.
        !
        ! Note: Since the energy of the photon package is peeled off 
        !       continuously by the absorption, we must use the average
        !       energy over this ray element, enerav.
        !
        if((alpha_s_tot.gt.0.d0).and.(allocated(mcscat_scatsrc_iquv))) then
           if(scattering_mode.eq.1) then
              !
              ! Isotropic scattering: simply add the source term
              !
              scatsrc0 = dtauscat * enerav / ( cellvolume(ray_index) * 12.566371d0 )
              mcscat_scatsrc_iquv(ray_inu,ray_index,1,1) =                            &
                      mcscat_scatsrc_iquv(ray_inu,ray_index,1,1) + scatsrc0
           elseif((scattering_mode.eq.2).or.(scattering_mode.eq.3)) then
              !
              ! 2D scattering mode only for scattering_mode 5
              !
              if(dust_2daniso) then
                 write(stdo,*) 'ERROR: non-isotropic scattering for 2-D axisymmetric models only for scattering_mode.ge.5'
                 stop
              endif
              !
              ! Anisotropic scattering: add only for the given directions
              !
              scatsrc0 = dtauscat * enerav / ( cellvolume(ray_index) * 12.566371d0 )
              if(.not.mcscat_localobserver) then
                 !
                 ! For the observer at infinity this is easy. We already
                 ! have the pre-calculated phase functions.
                 !
                 do idir=1,mcscat_nrdirs
                    dum = 0.d0
                    do ispec=1,dust_nr_species
                       dum       = dum + mcscat_phasefunc(idir,ispec) *           &
                                         alpha_s(ispec) / alpha_s_tot
                    enddo
                    mcscat_scatsrc_iquv(ray_inu,ray_index,1,idir) =      &
                         mcscat_scatsrc_iquv(ray_inu,ray_index,1,idir) + &
                         scatsrc0 * dum
                 enddo
              else
                 !
                 ! For now we do not include local observer + anisotropic
                 ! scattering yet. What needs to be done here is to calculate
                 ! on-the-fly the direction toward the observer, and then
                 ! to calculate the phase function from this. 
                 !
                 write(stdo,*) 'ABORTING: For now, non-isotropic scattering and '
                 write(stdo,*) '          local observer are not yet compatible.'
                 stop
                 !
                 ! Here a preliminary attempt at this. But this is not
                 ! at all yet tested!! 
                 !
                 mcscat_dirs(1,1) = mcscat_localobs_pos(1) - ray_cart_x
                 mcscat_dirs(2,1) = mcscat_localobs_pos(2) - ray_cart_y
                 mcscat_dirs(3,1) = mcscat_localobs_pos(3) - ray_cart_z
                 dummy = 1.d0 / sqrt( mcscat_dirs(1,1)**2 + &
                         mcscat_dirs(2,1)**2 + mcscat_dirs(3,1)**2 )
                 mcscat_dirs(1,1) = mcscat_dirs(1,1) * dummy
                 mcscat_dirs(2,1) = mcscat_dirs(2,1) * dummy
                 mcscat_dirs(3,1) = mcscat_dirs(3,1) * dummy
                 costheta = ray_cart_dirx*mcscat_dirs(1,1) +          &
                            ray_cart_diry*mcscat_dirs(2,1) +          &
                            ray_cart_dirz*mcscat_dirs(3,1)
                 do ispec=1,dust_nr_species
                    g                         = kappa_g(ray_inu,ispec)
                    mcscat_phasefunc(1,ispec) = henyeygreenstein_phasefunc(g,costheta)
                 enddo
                 dum = 0.d0
                 do ispec=1,dust_nr_species
                    dum       = dum + mcscat_phasefunc(1,ispec) *          &
                                      alpha_s(ispec) / alpha_s_tot
                 enddo
                 mcscat_scatsrc_iquv(ray_inu,ray_index,1,1) =      &
                      mcscat_scatsrc_iquv(ray_inu,ray_index,1,1) + &
                      scatsrc0 * dum
              endif
           elseif(scattering_mode.eq.4) then
              !
              ! 2D scattering mode only for scattering_mode 5
              !
              if(dust_2daniso) then
                 write(stdo,*) 'ERROR: non-isotropic scattering for 2-D axisymmetric models only for scattering_mode.ge.5'
                 stop
              endif
              !
              ! Polarized scattering, but only for the scattering 
              ! source function (we assume that the current photon is
              ! unpolarized)
              !
              if(mcscat_localobserver) then
                 write(stdo,*) 'ERROR: Polarized scattering and local-observer are incompatible.'
                 stop
              endif
              !
              ! Calculate the distance travelled to the final point
              !
              dss = dtauscat / alpha_s_tot
              !
              ! Make the "phot" array for use with the polarization
              ! module. We choose an arbitrary S-vector, because the
              ! incoming light is assumed to be unpolarized.
              !
              photpkg%E    = enerav
              photpkg%Q    = 0.d0
              photpkg%U    = 0.d0
              photpkg%V    = 0.d0
              photpkg%n(1) = ray_cart_dirx
              photpkg%n(2) = ray_cart_diry
              photpkg%n(3) = ray_cart_dirz
              call polarization_make_s_vector(photpkg%n,photpkg%s)   ! arb S-vector
              !
              ! Now add the scattering contribution of each of
              ! the dust species. Note that we do not
              ! need to divide by 4*pi here because we use 
              ! Z instead of kappa_scat.
              !
              do ispec=1,dust_nr_species
                 call polarization_randomorient_scatsource(photpkg, &
                      mcscat_dirs(:,1),mcscat_svec(:,1),            &
                      scat_munr,scat_mui_grid(:),                   &
                      zmatrix(:,ray_inu,1,ispec),                   &
                      zmatrix(:,ray_inu,2,ispec),                   &
                      zmatrix(:,ray_inu,3,ispec),                   &
                      zmatrix(:,ray_inu,4,ispec),                   &
                      zmatrix(:,ray_inu,5,ispec),                   &
                      zmatrix(:,ray_inu,6,ispec),                   &
                      src4)
                 mcscat_scatsrc_iquv(ray_inu,ray_index,1:4,1) =      &
                      mcscat_scatsrc_iquv(ray_inu,ray_index,1:4,1) + &
                      dustdens(ispec,ray_index) * src4(1:4) * dss /  &
                      cellvolume(ray_index)
              enddo
              !
           elseif(scattering_mode.eq.5) then
              !
              ! Polarized scattering with full polarized Monte Carlo
              ! photon package
              !
              if(mcscat_localobserver) then
                 write(stdo,*) 'ERROR: Polarized scattering and local-observer are incompatible.'
                 stop
              endif
              !
              ! Calculate the distance travelled to the final point
              !
              dss = dtauscat / alpha_s_tot
              !
              ! Since we will need the average energy over this
              ! bit of path, instead of the inital one, we backup
              ! the initial values of E,Q,U,V and replace them
              ! with the average values
              !
              Ebk = photpkg%E
              Qbk = photpkg%Q
              Ubk = photpkg%U
              Vbk = photpkg%V
              photpkg%E = photpkg%E * xxtauabs
              photpkg%Q = photpkg%Q * xxtauabs
              photpkg%U = photpkg%U * xxtauabs
              photpkg%V = photpkg%V * xxtauabs              
              !
              ! If 2-D axisymmetric anisotropic scattering mode,
              ! then we must rotate the photpkg%n and photpkg%s
              ! vectors appropriately
              !
              if(dust_2daniso) then
                 !
                 ! Safety check. Note: the last phi point is
                 ! 360 degrees == first one. This makes it easier
                 ! lateron to interpolate, even though it costs
                 ! a tiny bit more memory and is a tiny bit slower.
                 !
                 if(mcscat_nrdirs.ne.dust_2daniso_nphi+1) stop 2067
                 !
                 ! Backup original direction and s-vector
                 !
                 nvec_orig(:) = photpkg%n(:)
                 svec_orig(:) = photpkg%s(:)
                 !
                 ! Rotate n and s vector such that the event lies
                 ! in the x-z-plane (positive x)
                 !
                 levent      = sqrt(ray_cart_x**2+ray_cart_y**2)
                 if(levent.gt.0.d0) then
                    cosphievent  = ray_cart_x/levent
                    sinphievent  = ray_cart_y/levent
                    xbk          = photpkg%n(1)
                    ybk          = photpkg%n(2)
                    photpkg%n(1) =  cosphievent * xbk + sinphievent * ybk
                    photpkg%n(2) = -sinphievent * xbk + cosphievent * ybk
                    xbk          = photpkg%s(1)
                    ybk          = photpkg%s(2)
                    photpkg%s(1) =  cosphievent * xbk + sinphievent * ybk
                    photpkg%s(2) = -sinphievent * xbk + cosphievent * ybk
                 endif
                 !
                 ! Pre-compute the cos and sin of the small incremental
                 ! rotations
                 !
                 deltaphi = twopi / dust_2daniso_nphi
                 cosdphi  = cos(deltaphi)
                 sindphi  = sin(deltaphi)
              endif
              !
              ! Now add the scattering contribution of each of
              ! the dust species. Note that we do not
              ! need to divide by 4*pi here because we use 
              ! Z instead of kappa_scat.
              !
              ! Loop over all viewing directions (important
              ! for multi-vantage-point images = movies, and
              ! for 2-D full-phase scattering).
              !
              do idirs=1,mcscat_nrdirs
                 do ispec=1,dust_nr_species
                    call polarization_randomorient_scatsource(photpkg, &
                         mcscat_dirs(:,idirs),mcscat_svec(:,idirs),    &
                         scat_munr,scat_mui_grid(:),                   &
                         zmatrix(:,ray_inu,1,ispec),                   &
                         zmatrix(:,ray_inu,2,ispec),                   &
                         zmatrix(:,ray_inu,3,ispec),                   &
                         zmatrix(:,ray_inu,4,ispec),                   &
                         zmatrix(:,ray_inu,5,ispec),                   &
                         zmatrix(:,ray_inu,6,ispec),                   &
                         src4)
                    mcscat_scatsrc_iquv(ray_inu,ray_index,1:4,idirs) =      &
                         mcscat_scatsrc_iquv(ray_inu,ray_index,1:4,idirs) + &
                         dustdens(ispec,ray_index) * src4(1:4) * dss /      &
                         cellvolume(ray_index)
                 enddo
                 !
                 ! If 2-D axisymmetric anisotropic scattering mode,
                 ! then we must rotate the photpkg%n and photpkg%s
                 ! vectors appropriately
                 !
                 if(dust_2daniso) then
                    xbk = photpkg%n(1)
                    ybk = photpkg%n(2)
                    photpkg%n(1) =  cosdphi * xbk - sindphi * ybk
                    photpkg%n(2) =  sindphi * xbk + cosdphi * ybk
                    xbk = photpkg%s(1)
                    ybk = photpkg%s(2)
                    photpkg%s(1) =  cosdphi * xbk - sindphi * ybk
                    photpkg%s(2) =  sindphi * xbk + cosdphi * ybk
                 endif
                 !
              enddo
              !
              ! Restore the original values
              !
              photpkg%E = Ebk
              photpkg%Q = Qbk
              photpkg%U = Ubk
              photpkg%V = Vbk
              if(dust_2daniso) then
                 !
                 ! Restore original direction and s-vector
                 !
                 photpkg%n(:) = nvec_orig(:)
                 photpkg%s(:) = svec_orig(:)
                 !
              endif
           endif
        endif
        !
        ! Add photons to mean intensity
        !
        ! Note: Since the energy of the photon package is peeled off 
        !       continuously by the absorption, we must use the average
        !       energy over this ray element, enerav.
        !
        ! Note: Since we only stop in the middle of a cell if we have
        !       non-zero scattering opacity, we can compute the mean
        !       intensity via the scattering. 
        !
        if(allocated(mcscat_meanint)) then
           ! BUGFIX 14.02.2012: Should use the scat-independent way
           !!!! scatsrc0 = dtauscat * enerav / ( cellvolume(ray_index) * 12.566371d0 )
           !!!! mcscat_meanint(ray_inu,ray_index) =                            &
           !!!!      mcscat_meanint(ray_inu,ray_index) + scatsrc0 / alpha_s_tot
           ! BUGFIX 03.05.2014 (thanks Seokho Lee for pointing this out)
           !        It seems BUGFIX 14.02.2012 was not a bugfix but a 
           !        buginjection... ds is the distance to the next cellwall
           !        whereas here we must use the distance to the next 
           !        scattering event.
           !        However, for now let's add also some checks
           !!!! mnint = ds * enerav / ( cellvolume(ray_index) * 12.566371d0 )
           !
           if(alpha_s_tot.eq.0.d0) stop 487    ! Trivial test (can be removed)
           if(dtauscat.lt.0.d0) stop 488       ! Trivial test (can be removed)
           dss   = dtauscat / alpha_s_tot
           mnint = dss * enerav / ( cellvolume(ray_index) * 12.566371d0 )
           mcscat_meanint(ray_inu,ray_index) =                            &
                mcscat_meanint(ray_inu,ray_index) + mnint
        endif
        endif ! This endif is for the selectscat analysis stuff (normally not important)
        !
        ! Compute the location of this point
        !
        fr     = (taupath-tau)/dtau
        ray_cart_x = prev_x + fr * ( ray_cart_x - prev_x )
        ray_cart_y = prev_y + fr * ( ray_cart_y - prev_y )
        ray_cart_z = prev_z + fr * ( ray_cart_z - prev_z )
        !
        ! For anisotropic scattering we need to know off which dust species
        ! the photon has scattered. Determine this here, and call it ispecc.
        !
        if(scattering_mode.gt.1) then
           alphacum(1) = 0.d0
           do ispec=1,dust_nr_species
              alphacum(ispec+1) = alphacum(ispec) + alpha_s(ispec)
           enddo
           if(alphacum(dust_nr_species+1).le.0.d0) stop 2244
           do ispec=1,dust_nr_species
              alphacum(ispec) = alphacum(ispec) / alphacum(dust_nr_species+1)
           enddo
           alphacum(dust_nr_species+1) = 1.d0
           rn = ran2(iseed)
           call hunt(alphacum,dust_nr_species+1,rn,ispecc)
           if((ispecc.lt.1).or.(ispecc.gt.dust_nr_species)) stop 7266
        else
           rn = ran2(iseed)   ! NOTE: This is not necessary, but is just
                              !       there to pass the selftest. 
                              !       See notes 15.11.09 in Radmc_3D_LOG.txt
           ispecc=0    ! Force error (with array bound check) if accidently used
        endif
        !
        ! Compute the new ener
        !  
        ener    = ener * xptauabs
        !
        ! Update the photon package for the polarization
        !
        if(scattering_mode.ge.5) then
           photpkg%E = photpkg%E * xptauabs
           photpkg%Q = photpkg%Q * xptauabs
           photpkg%U = photpkg%U * xptauabs
           photpkg%V = photpkg%V * xptauabs
        endif
        !
        ! Check if we drop below the "throw photon away" limit
        !
        ! **** NOTE: I should also set mc_photon_destroyed = .true. ****
        ! ****       but it does not matter right now.              ****
        !
        if(ener.le.mc_scat_energy_rellimit*energy) then
           arrived = .true.
        else
           arrived = .false.
        endif
        !
        ! Done with this segment of the random walk, so return and 
        ! handle the scattering or absorption event
        !
        ok=.false.
        !
        ! OpenMP Parallellization: Release lock on this cell
        !
        !$ if(ray_index .ge. 1 )then
        !$    call omp_unset_lock(lock(ray_index));
        !$ endif
        return
        !
     endif
     !
     ! Endpoint of this segment.
     !
     if(ray_index.ge.1) then
        !
        ! We are in a cell. 
        !
        ! Compute dtauabs and dtauscat as well as exp(-dtauabs) and
        ! (1-exp(-dtauabs))/dtauabs.
        !
        ! MAJOR BUGFIX 18.01.10:
        !
        dtauscat = dtau
        dtauabs  = dtauscat * alpha_a_tot / alpha_s_tot
        xptauabs = exp(-dtauabs)
        if(dtauabs.gt.1d-6) then
           xxtauabs = (1.d0-xptauabs)/dtauabs
        else
           xxtauabs = 1.d0-0.5d0*dtauabs
        endif
        !
        ! Compute enerav
        !
        enerav  = ener * xxtauabs
        !
        ! For the selectscat analysis stuff (normally not important)
        !
        if((selectscat_iscat.ge.selectscat_iscat_first).and.(selectscat_iscat.le.selectscat_iscat_last)) then
        !
        ! Add photons to scattering source
        !
        ! NOTE: In the original RADMC I did not divide by 4*pi yet; only
        !       in the writing to scatsource.dat I did that. Here I do 
        !       it immediately.
        !
        ! Note: Since we follow only the scattering events, the tau
        !       and taupath are scattering tau only. So no `albedo'
        !       is necessary here.
        !
        ! Note: Since the energy of the photon package is peeled off 
        !       continuously by the absorption, we must use the average
        !       energy over this ray element, enerav.
        !
        if((alpha_s_tot.gt.0.d0).and.(allocated(mcscat_scatsrc_iquv))) then
           if(scattering_mode.eq.1) then
              !
              ! Isotropic scattering: simply add the source term
              !
              scatsrc0 = dtau * enerav / ( cellvolume(ray_index) * 12.566371d0 )
              mcscat_scatsrc_iquv(ray_inu,ray_index,1,1) =                            &
                      mcscat_scatsrc_iquv(ray_inu,ray_index,1,1) + scatsrc0
           elseif((scattering_mode.eq.2).or.(scattering_mode.eq.3)) then
              !
              ! Anisotropic scattering: add only for the given directions
              !
              scatsrc0 = dtau * enerav / ( cellvolume(ray_index) * 12.566371d0 )
              if(.not.mcscat_localobserver) then
                 !
                 ! For the observer at infinity this is easy. We already
                 ! have the pre-calculated phase functions.
                 !
                 do idir=1,mcscat_nrdirs
                    dum = 0.d0
                    do ispec=1,dust_nr_species
                       dum       = dum + mcscat_phasefunc(idir,ispec) *           &
                                         alpha_s(ispec) / alpha_s_tot
                    enddo
                    mcscat_scatsrc_iquv(ray_inu,ray_index,1,idir) =      &
                         mcscat_scatsrc_iquv(ray_inu,ray_index,1,idir) + &
                         scatsrc0 *dum
                 enddo
              else
                 !
                 ! For now we do not include local observer + anisotropic
                 ! scattering yet. What needs to be done here is to calculate
                 ! on-the-fly the direction toward the observer, and then
                 ! to calculate the phase function from this. 
                 !
                 write(stdo,*) 'ABORTING: For now, non-isotropic scattering and '
                 write(stdo,*) '          local observer are not yet compatible.'
                 stop
                 !
                 ! Here a preliminary attempt at this. But this is not
                 ! at all yet tested!! 
                 !
                 mcscat_dirs(1,1) = mcscat_localobs_pos(1) - ray_cart_x
                 mcscat_dirs(2,1) = mcscat_localobs_pos(2) - ray_cart_y
                 mcscat_dirs(3,1) = mcscat_localobs_pos(3) - ray_cart_z
                 dummy = 1.d0 / sqrt( mcscat_dirs(1,1)**2 + &
                         mcscat_dirs(2,1)**2 + mcscat_dirs(3,1)**2 )
                 mcscat_dirs(1,1) = mcscat_dirs(1,1) * dummy
                 mcscat_dirs(2,1) = mcscat_dirs(2,1) * dummy
                 mcscat_dirs(3,1) = mcscat_dirs(3,1) * dummy
                 costheta = ray_cart_dirx*mcscat_dirs(1,1) +          &
                            ray_cart_diry*mcscat_dirs(2,1) +          &
                            ray_cart_dirz*mcscat_dirs(3,1)
                 do ispec=1,dust_nr_species
                    g                         = kappa_g(ray_inu,ispec)
                    mcscat_phasefunc(1,ispec) = henyeygreenstein_phasefunc(g,costheta)
                 enddo
                 dum = 0.d0
                 do ispec=1,dust_nr_species
                    dum       = dum + mcscat_phasefunc(1,ispec) *          &
                                      alpha_s(ispec) / alpha_s_tot
                 enddo
                 mcscat_scatsrc_iquv(ray_inu,ray_index,1,1) =      &
                      mcscat_scatsrc_iquv(ray_inu,ray_index,1,1) + &
                      scatsrc0 * dum
              endif
           elseif(scattering_mode.eq.4) then
              !
              ! Polarized scattering, but only for the scattering 
              ! source function (we assume that the current photon is
              ! unpolarized)
              !
              if(mcscat_localobserver) then
                 write(stdo,*) 'ERROR: Polarized scattering and local-observer are incompatible.'
                 stop
              endif
              !
              ! Make the "phot" array for use with the polarization
              ! module. We choose an arbitrary S-vector, because the
              ! incoming light is assumed to be unpolarized.
              !
              photpkg%E    = enerav
              photpkg%Q    = 0.d0
              photpkg%U    = 0.d0
              photpkg%V    = 0.d0
              photpkg%n(1) = ray_cart_dirx
              photpkg%n(2) = ray_cart_diry
              photpkg%n(3) = ray_cart_dirz
              call polarization_make_s_vector(photpkg%n,photpkg%s)   ! arb S-vector
              !
              ! Now add the scattering contribution of each of
              ! the dust species. Note that we do not
              ! need to divide by 4*pi here because we use 
              ! Z instead of kappa_scat.
              !
              do ispec=1,dust_nr_species
                 call polarization_randomorient_scatsource(photpkg, &
                      mcscat_dirs(:,1),mcscat_svec(:,1),            &
                      scat_munr,scat_mui_grid(:),                   &
                      zmatrix(:,ray_inu,1,ispec),                   &
                      zmatrix(:,ray_inu,2,ispec),                   &
                      zmatrix(:,ray_inu,3,ispec),                   &
                      zmatrix(:,ray_inu,4,ispec),                   &
                      zmatrix(:,ray_inu,5,ispec),                   &
                      zmatrix(:,ray_inu,6,ispec),                   &
                      src4)
                 mcscat_scatsrc_iquv(ray_inu,ray_index,1:4,1) =      &
                      mcscat_scatsrc_iquv(ray_inu,ray_index,1:4,1) + &
                      dustdens(ispec,ray_index) * src4(1:4) * ds /   &
                      cellvolume(ray_index)
              enddo
              !
           elseif(scattering_mode.eq.5) then
              !
              ! Polarized scattering with full polarized Monte Carlo
              ! photon package
              !
              if(mcscat_localobserver) then
                 write(stdo,*) 'ERROR: Polarized scattering and local-observer are incompatible.'
                 stop
              endif
              !
              ! Since we will need the average energy over this
              ! bit of path, instead of the inital one, we backup
              ! the initial values of E,Q,U,V and replace them
              ! with the average values
              !
              Ebk = photpkg%E
              Qbk = photpkg%Q
              Ubk = photpkg%U
              Vbk = photpkg%V
              photpkg%E = photpkg%E * xxtauabs
              photpkg%Q = photpkg%Q * xxtauabs
              photpkg%U = photpkg%U * xxtauabs
              photpkg%V = photpkg%V * xxtauabs              
              !
              ! If 2-D axisymmetric anisotropic scattering mode,
              ! then we must rotate the photpkg%n and photpkg%s
              ! vectors appropriately
              !
              if(dust_2daniso) then
                 !
                 ! Safety check. Note: the last phi point is
                 ! 360 degrees == first one. This makes it easier
                 ! lateron to interpolate, even though it costs
                 ! a tiny bit more memory and is a tiny bit slower.
                 !
                 if(mcscat_nrdirs.ne.dust_2daniso_nphi+1) stop 2067
                 !
                 ! Backup original direction and s-vector
                 !
                 nvec_orig(:) = photpkg%n(:)
                 svec_orig(:) = photpkg%s(:)
                 !
                 ! Rotate n and s vector such that the event lies
                 ! in the x-z-plane (positive x)
                 !
                 levent      = sqrt(ray_cart_x**2+ray_cart_y**2)
                 if(levent.gt.0.d0) then
                    cosphievent  = ray_cart_x/levent
                    sinphievent  = ray_cart_y/levent
                    xbk          = photpkg%n(1)
                    ybk          = photpkg%n(2)
                    photpkg%n(1) =  cosphievent * xbk + sinphievent * ybk
                    photpkg%n(2) = -sinphievent * xbk + cosphievent * ybk
                    xbk          = photpkg%s(1)
                    ybk          = photpkg%s(2)
                    photpkg%s(1) =  cosphievent * xbk + sinphievent * ybk
                    photpkg%s(2) = -sinphievent * xbk + cosphievent * ybk
                 endif
                 !
                 ! Pre-compute the cos and sin of the small incremental
                 ! rotations
                 !
                 deltaphi = twopi / dust_2daniso_nphi
                 cosdphi  = cos(deltaphi)
                 sindphi  = sin(deltaphi)
              endif
              !
              ! Now add the scattering contribution of each of
              ! the dust species. Note that since src4 is proportional
              ! to Z11*enerav, no division by 4*pi is necessary, since
              ! <Z11> = kappa_scat/(4*pi).
              !
              ! Loop over all viewing directions (important
              ! for multi-vantage-point images = movies, and
              ! for 2-D full-phase scattering).
              !
              do idirs=1,mcscat_nrdirs
                 do ispec=1,dust_nr_species
                    call polarization_randomorient_scatsource(photpkg, &
                         mcscat_dirs(:,idirs),mcscat_svec(:,idirs),    &
                         scat_munr,scat_mui_grid(:),                   &
                         zmatrix(:,ray_inu,1,ispec),                   &
                         zmatrix(:,ray_inu,2,ispec),                   &
                         zmatrix(:,ray_inu,3,ispec),                   &
                         zmatrix(:,ray_inu,4,ispec),                   &
                         zmatrix(:,ray_inu,5,ispec),                   &
                         zmatrix(:,ray_inu,6,ispec),                   &
                         src4)
                    mcscat_scatsrc_iquv(ray_inu,ray_index,1:4,idirs) =      &
                         mcscat_scatsrc_iquv(ray_inu,ray_index,1:4,idirs) + &
                         dustdens(ispec,ray_index) * src4(1:4) * ds /       &
                         cellvolume(ray_index)
                 enddo
                 !
                 ! If 2-D axisymmetric anisotropic scattering mode,
                 ! then we must rotate the photpkg%n and photpkg%s
                 ! vectors appropriately
                 !
                 if(dust_2daniso) then
                    xbk = photpkg%n(1)
                    ybk = photpkg%n(2)
                    photpkg%n(1) =  cosdphi * xbk - sindphi * ybk
                    photpkg%n(2) =  sindphi * xbk + cosdphi * ybk
                    xbk = photpkg%s(1)
                    ybk = photpkg%s(2)
                    photpkg%s(1) =  cosdphi * xbk - sindphi * ybk
                    photpkg%s(2) =  sindphi * xbk + cosdphi * ybk
                 endif
                 !
              enddo
              !
              ! Restore the original values
              !
              photpkg%E = Ebk
              photpkg%Q = Qbk
              photpkg%U = Ubk
              photpkg%V = Vbk
              if(dust_2daniso) then
                 !
                 ! Restore original direction and s-vector
                 !
                 photpkg%n(:) = nvec_orig(:)
                 photpkg%s(:) = svec_orig(:)
                 !
              endif
           endif
        endif
        !
        ! Add photons to mean intensity
        !
        ! Note: Since the energy of the photon package is peeled off 
        !       continuously by the absorption, we must use the average
        !       energy over this ray element, enerav.
        !
        if(allocated(mcscat_meanint)) then
           mnint = ds * enerav / ( cellvolume(ray_index) * 12.566371d0 )
           mcscat_meanint(ray_inu,ray_index) =                            &
                mcscat_meanint(ray_inu,ray_index) + mnint
        endif
        endif ! This endif is for the selectscat analysis stuff (normally not important)
        !
        ! Compute the new ener and the new tauabs
        !  
        ener = ener * xptauabs
        !
        ! Update the photon package for the polarization
        !
        if(scattering_mode.ge.5) then
           photpkg%E = photpkg%E * xptauabs
           photpkg%Q = photpkg%Q * xptauabs
           photpkg%U = photpkg%U * xptauabs
           photpkg%V = photpkg%V * xptauabs
        endif
        !
        ! Check if we drop below the "throw photon away" limit
        !
        ! **** NOTE: I should also set mc_photon_destroyed = .true. ****
        ! ****       but it does not matter right now.              ****
        !
        if(ener.le.mc_scat_energy_rellimit*energy) then
           arrived = .true.
        endif
        !
     endif
     !
     ! Increase the tau
     !
     tau = tau + dtau
     !
     ! OpenMP Parallellization: Release lock on this cell
     !
     !$ if(ray_index .ge. 1 )then
     !$    call omp_unset_lock(lock(ray_index));
     !$ endif
     !
     ! Now make next cell the current cell
     !
     ray_index = ray_indexnext
     !
     ! Do the same for the pointers.
     ! This is grid type dependent.
     !
     if(igrid_type.lt.100) then
        !
        ! AMR type grid
        !
        if(ray_index.gt.0) then
           !
           ! We are in a cell
           !
           if(amr_tree_present) then
              amrray_cell => amr_index_to_leaf(ray_index)%link
           else
              call amr_regular_get_ixyz(ray_index,amrray_ix_curr,amrray_iy_curr,amrray_iz_curr)
           endif
        else
           nullify(amrray_cell)
           amrray_ix_curr = -1
           amrray_iy_curr = -1
           amrray_iz_curr = -1
        endif
     else
        write(stdo,*) 'ERROR: This grid type not yet implemented'
        stop
     endif
     !
     ! Now, if this new/next segment goes off to infinity, then
     ! we are done.
     !
     if(arrived) then
        !
        ! We escaped to infinity. 
        !
        ok=.false.
        !
        ! OpenMP Parallellization: Lock is already released.
        ! BUGFIX 23.02.2017
        !
        return
     endif
     !
  enddo
  !
end subroutine walk_cells_scat



!--------------------------------------------------------------------------
!               CALCULATE THE OPTICAL DEPTH ALONG A RAY
!
! This routine is used by the single scattering approximation to calculate
! the optical depth between a point-like star and a cell.
!--------------------------------------------------------------------------
subroutine walk_cells_optical_depth(params,lenpath,inu,tau)
  implicit none
  type(mc_params) :: params
  integer :: ispec,inu,idir,iddr,index
  doubleprecision :: lenpath
  doubleprecision :: tau,dtau
  doubleprecision :: ds,rn,scatsrc0,mnint,sprev,stot
  logical :: arrived
  doubleprecision :: prev_x,prev_y,prev_z
  doubleprecision :: axi(1:2,1:3)
  !
  ! Reset
  !
  tau       = 0.d0
  stot      = 0.d0
  sprev     = 0.d0
  !
  ! Default
  !
  arrived   = .false.
  ray_inu   = inu
  !
  ! In this subroutine we will not make use of the amrray option to
  ! advance only partly within a cell. We will check here if we have
  ! an event within this cell and compute here the location of this
  ! event. Therefore the ray_dsend is set to 1d99.
  !
  ray_dsend = 1d99
  !
  ! Now start moving this photon
  !
  do while(.not.arrived)
     !
     ! Do some safety checks
     !
     if(debug_check_all.eq.1) then
        if(ray_index.ge.1) then
           if(igrid_type.lt.100) then
              !
              ! AMR grid
              !
              ! Get the cell walls
              !
              if(amr_tree_present) then
                 if(.not.associated(amrray_cell)) stop 7270
                 do iddr=1,3
                    axi(1,iddr) = amr_finegrid_xi(amrray_cell%ixyzf(iddr),iddr,amrray_cell%level)
                    axi(2,iddr) = amr_finegrid_xi(amrray_cell%ixyzf(iddr)+1,iddr,amrray_cell%level)
                 enddo
              else
                 if(amrray_ix_curr.le.0) stop 7270
                 axi(1,1) = amr_finegrid_xi(amrray_ix_curr,1,0)
                 axi(2,1) = amr_finegrid_xi(amrray_ix_curr+1,1,0)
                 axi(1,2) = amr_finegrid_xi(amrray_iy_curr,2,0)
                 axi(2,2) = amr_finegrid_xi(amrray_iy_curr+1,2,0)
                 axi(1,3) = amr_finegrid_xi(amrray_iz_curr,3,0)
                 axi(2,3) = amr_finegrid_xi(amrray_iz_curr+1,3,0)
              endif
              !
              if(igrid_coord.lt.100) then
                 !
                 ! Cartesian coordinates
                 !
                 if((ray_cart_x.lt.axi(1,1)).or.(ray_cart_x.gt.axi(2,1)).or. &
                    (ray_cart_y.lt.axi(1,2)).or.(ray_cart_y.gt.axi(2,2)).or. &
                    (ray_cart_z.lt.axi(1,3)).or.(ray_cart_z.gt.axi(2,3))) then
                    write(stdo,*) 'INTERNAL ERROR: Point outside of cell in montecarlo module'
                    write(stdo,*) ray_cart_x,axi(1,1),axi(2,1)
                    write(stdo,*) ray_cart_y,axi(1,2),axi(2,2)
                    write(stdo,*) ray_cart_z,axi(1,3),axi(2,3)
                    stop
                 endif
              else
                 !
                 ! These coordiantes not yet active
                 !
                 write(stdo,*) 'INTERNAL ERROR: Non-Cart coordinates not yet implemented'
                 stop 
              endif
           else
              !
              ! Unstructured grids not yet active
              !
              write(stdo,*) 'INTERNAL ERROR: Unstructured grids not yet implemented'
              stop            
           endif
        endif
     endif
     !
     ! Back up the current position
     !
     prev_x     = ray_cart_x
     prev_y     = ray_cart_y
     prev_z     = ray_cart_z
     !
     ! Find new position
     !
     if(igrid_type.lt.100) then 
        !
        ! We use a regular or AMR grid
        !
        if(igrid_coord.lt.100) then
           !
           ! We use cartesian coordinates
           !
           call amrray_find_next_location_cart(ray_dsend,            &
                ray_cart_x,ray_cart_y,ray_cart_z,                    &
                ray_cart_dirx,ray_cart_diry,ray_cart_dirz,           &
                ray_index,ray_indexnext,ray_ds,arrived)
        elseif(igrid_coord.lt.200) then
           !
           ! We use spherical coordinates
           !
           call amrray_find_next_location_spher(ray_dsend,           &
                ray_cart_x,ray_cart_y,ray_cart_z,                    &
                ray_cart_dirx,ray_cart_diry,ray_cart_dirz,           &
                ray_index,ray_indexnext,ray_ds,arrived)
        else
           write(stdo,*) 'ERROR: Cylindrical coordinates not yet implemented'
           stop
        endif
        !
        ! For this mode stars are not allowed to be treated as 
        ! finite-size spheres.
        !
        if(amrray_ispherehit.lt.0) then
           write(stdo,*) 'ERROR: Finite-size (non-pointlike) stars not allowed in the single-scattering mode.'
           stop 162
        endif
     else
        write(stdo,*) 'SORRY: Delaunay or Voronoi grids not yet implemented'
        stop
     endif
     !
     ! Path length
     ! 
     ! ######### CHECK: WHY NOT USE ray_ds ??? ########
     !
     ds   = sqrt( (ray_cart_x-prev_x)**2 + (ray_cart_y-prev_y)**2 + (ray_cart_z-prev_z)**2 )
     !
     ! Check if we already arrived
     !
     if(sprev.gt.lenpath) then
        write(stdo,*) 'ERROR in single scattering mode. Contact author.'
        write(stdo,*) sprev,lenpath
        stop 7620
     endif
     sprev = stot
     stot  = stot + ds
     if(stot.gt.lenpath) then
        ds      = lenpath - sprev
        stot    = lenpath
        arrived = .true.
     endif
     !
     ! Compute the alpha's, dtau
     !
     ! NOTE: For non-Cartesian coordinates it can happen that we have not
     !       yet arrived, but are currently not in a cell either. Hence 
     !       we need to check if ray_index.ge.1
     !
     if(ray_index.ge.1) then
        alpha_a_tot = 0.d0
        alpha_s_tot = 0.d0
        do ispec=1,dust_nr_species
           alpha_a(ispec) = dustdens(ispec,ray_index) * &
                            kappa_a(ray_inu,ispec) + 1d-99
           alpha_s(ispec) = dustdens(ispec,ray_index) * &
                            kappa_s(ray_inu,ispec) + 1d-99
           alpha_a_tot    = alpha_a_tot + alpha_a(ispec)
           alpha_s_tot    = alpha_s_tot + alpha_s(ispec)
        enddo
        dtau      = ( alpha_a_tot + alpha_s_tot ) * ds
     else
        do ispec=1,dust_nr_species
           alpha_a(ispec) = 0.d0
           alpha_s(ispec) = 0.d0
        enddo
        alpha_a_tot = 0.d0
        alpha_s_tot = 0.d0
        dtau        = 0.d0
     endif
     !
     ! Increase the tau
     !
     tau = tau + dtau
     !
     ! Now make next cell the current cell
     !
     ray_index = ray_indexnext
     !
     ! Do the same for the pointers.
     ! This is grid type dependent.
     !
     if(igrid_type.lt.100) then
        !
        ! AMR type grid
        !
        if(ray_index.gt.0) then
           !
           ! We are in a cell
           !
           if(amr_tree_present) then
              amrray_cell => amr_index_to_leaf(ray_index)%link
           else
              call amr_regular_get_ixyz(ray_index,amrray_ix_curr,amrray_iy_curr,amrray_iz_curr)
           endif
        else
           nullify(amrray_cell)
           amrray_ix_curr = -1
           amrray_iy_curr = -1
           amrray_iz_curr = -1
        endif
     else
        write(stdo,*) 'ERROR: This grid type not yet implemented'
        stop
     endif
  enddo
  !
end subroutine walk_cells_optical_depth



!-------------------------------------------------------------------
!                DO OVERALL CHECK ON INPUT DATA
!-------------------------------------------------------------------
subroutine mc_check_inputdata
  
  implicit none
  integer :: index,icell,ispec
  !
  ! Check the dust density
  !
  if(allocated(dustdens)) then
     do ispec=1,dust_nr_species
        do icell=1,nrcells
           index = cellindex(icell)
           if((index.le.0).or.(index.gt.nrcellsmax)) then
              write(stdo,*) 'ERROR: Internal error while checking dust density'
              stop
           endif
           if(number_invalid(dustdens(ispec,index)).ne.0) then
              write(stdo,*) 'ERROR: Dust density is NAN or INF at some point!'
              stop
           endif
           if(dustdens(ispec,index).lt.0.d0) then
              write(stdo,*) 'ERROR: Dust density is negative at some point!'
              stop
           endif
           if(dustdens(ispec,index).lt.1d-99) then
              dustdens(ispec,index) = 1d-99
           endif
        enddo
     enddo
  endif
  !
  ! Check the dust temperature
  !
  if(allocated(dusttemp)) then
     do ispec=1,dust_nr_species
        do icell=1,nrcells
           index = cellindex(icell)
           if((index.le.0).or.(index.gt.nrcellsmax)) then
              write(stdo,*) 'ERROR: Internal error while checking dust temperature'
              stop
           endif
           if(number_invalid(dusttemp(ispec,index)).ne.0) then
              write(stdo,*) 'ERROR: Dust temperature is NAN or INF at some point!'
              stop
           endif
           if(dusttemp(ispec,index).lt.0.d0) then
              write(stdo,*) 'ERROR: Dust temperature is negative at some point!'
              stop
           endif
        enddo
     enddo
  endif
end subroutine mc_check_inputdata



!-------------------------------------------------------------------
!                         DO ABSORPTION EVENT
!
!     This subroutine increases the temperature of a cell according
!     to the mass it contains. And it randomly creates a new inucur.
!-------------------------------------------------------------------
subroutine do_absorption_event(params,iqactive,ierror)
  implicit none
  type(mc_params) :: params
  integer :: iqactive,ierror
  !
  integer :: ispec
  doubleprecision :: cumen,tempav,rn,fracquant,alpha_a_nu
  doubleprecision :: dum,dum1,dum2,rhodusttot,entotal
  logical :: calctemp
  !
  if(freq_dnu(1).eq.0.d0) stop 32222
  !
  ! Default
  !
  mc_photon_destroyed = .false.
  !
  ! First determine how the energy is going to be divided over
  ! the dust species and sizes
  !
  if(ray_inu.ge.1) then 
     !
     ! The normal way: dependent on frequency of incoming photon
     !
     if(dust_nr_species.eq.1) then
        !
        ! If there is only a single dust species: then we do not need to
        ! do any effort
        !
        mc_enerpart(1) = energy
        !
     else
        !
        ! If there are multiple dust species, then we need to figure out
        ! how to divide the absorbed energy over the various dust species.
        ! This depends on the frequency of the photon.
        !
        alpha_a_nu = 0.d0
        do ispec=1,dust_nr_species
           mc_enerpart(ispec) = dustdens(ispec,ray_index) *                &
                                kappa_a(ray_inu,ispec) 
           alpha_a_nu         = alpha_a_nu + mc_enerpart(ispec)
        enddo
        do ispec=1,dust_nr_species
           mc_enerpart(ispec) = energy * mc_enerpart(ispec) / alpha_a_nu
        enddo
        !
     endif
     !
     ! If quantum heating active:
     ! Original method:
     !    If the photon is absorbed by a quantum heated grain, then 
     !    destroy photon (because the quantum heating is done elsewhere)
     ! New method (NEW QUANTUM: 28.04.06):
     !    The photon is re-emitted according to the quantum source
     !
     if(iqactive.eq.1) then
        !
        ! ORIGINAL METHOD
        !
        ! Determine what is the chance that the photon is absorbed
        ! by a quantum-heated grain
        !
        fracquant = 0.d0
        do ispec=1,dust_nr_species
           if(dust_quantum(ispec).ne.0) then
              fracquant = fracquant + mc_enerpart(ispec)
           endif
        enddo
        fracquant = fracquant / energy
        !
        ! SMALL BUGFIX (06.07.05): 1.d0 --> 1.d0+1d-10
        !
        if(fracquant.gt.1.d0+MINIEPS) stop 92091
        !
        ! Now pick a random number and determine what happens.
        ! If absorbed by a quantum heated grain, then destroy the
        ! photon.
        !     
        rn = ran2(iseed)
        if(rn.le.fracquant) then
           mc_photon_destroyed = .true.
           return
        endif
     endif
  elseif(ray_inu.eq.-1) then
     !
     ! This is primary energy (e.g. viscous heating), and ray_inu.eq.-1 means
     ! that this energy is divided over the species by mass... This only
     ! works if dust_nr_size is always 1
     !
     ! NEWSTUFF 06.06.04
     !
     rhodusttot=0.d0
     do ispec=1,dust_nr_species
        rhodusttot = rhodusttot + dustdens(ispec,ray_index)
     enddo
     do ispec=1,dust_nr_species
        mc_enerpart(ispec) = energy * dustdens(ispec,ray_index) / rhodusttot
     enddo
  else
     write(stdo,*) 'INTERNAL ERROR: Dont know what to do with ray_inu = ',ray_inu
     stop
  endif
  !
  ! Stupidity check
  !
  if(debug_check_all.eq.1) then
     dum = 0.d0
     do ispec=1,dust_nr_species
        dum = dum + mc_enerpart(ispec)
     enddo
     if(abs(dum/energy-1.d0).gt.1.d-6) stop 34123
  endif
  !
  ! Determine if the dust temperature(s) and the frequency distribution
  ! for this cell must be recomputed
  !
  if(params%ifast.eq.0) then
     !
     ! For the original method, always recompute temperature and
     ! freq dist
     !
     calctemp = .true.
  else
     !
     ! For the new, fast method, check if the energy step is too high.
     !
     calctemp = .false.
     if(params%itempdecoup.eq.1) then
        !
        ! If temperatures are decoupled...
        !
        do ispec=1,dust_nr_species
           if((dust_quantum(ispec).eq.0).or.(iqactive.le.0)) then 
              if(abs((mc_cumulener(ispec,ray_index)+mc_enerpart(ispec))/      &
                      (mc_cumulener_bk(ispec,ray_index)+1d-99)-1.d0)          &
                      .gt.params%enthres) then
                 calctemp = .true.
              endif
           endif
        enddo
     else
        !
        ! If temperatures are coupled...
        !     
        dum1 = 0.d0
        dum2 = 0.d0
        dum  = 1.d-99
        do ispec=1,dust_nr_species
           if((dust_quantum(ispec).eq.0).or.(iqactive.le.0)) then 
              dum1 = dum1 + mc_cumulener(ispec,ray_index)
              dum2 = dum2 + mc_enerpart(ispec)
              dum  = dum  + mc_cumulener_bk(ispec,ray_index)
           endif
        enddo
        if(abs((dum1+dum2)/dum-1.d0).gt.params%enthres) then
           calctemp = .true.
        endif
     endif
  endif
  !
  ! Determine the temperature from the local energy
  ! (this need not be hugely accurate)
  !
  if(calctemp) then
     !
     ! OpenMP Parallellization: Lock this cell (and only continue when succesfully locked)
     !
     !$ if(ray_index .gt. 0)then
     !!$ continue=.true.
     !$ do while(.NOT. omp_test_lock(lock(ray_index)))
     !!$ if(continue)then
     !!$ conflict_counter=conflict_counter+1
     !!$ continue=.false.
     !!$ end if
     !$ end do
     !$ end if
     !
     if(params%itempdecoup.eq.1) then
        !
        ! If temperatures are decoupled...
        !
        do ispec=1,dust_nr_species
           if((dust_quantum(ispec).eq.0).or.(iqactive.le.0)) then 
              cumen = mc_cumulener(ispec,ray_index) /                      &
                      ( dustdens(ispec,ray_index) * cellvolume(ray_index) )
              dusttemp(ispec,ray_index) = compute_dusttemp_energy_bd(cumen,ispec)
           endif
        enddo
     else
        !
        ! If temperatures are coupled...
        !
        ! NOTE: Quantum-heated grains are included in the
        !       temperature determination below, even though
        !       the present photon has not been absorbed by one
        !       of them. This is because all grain temperatures
        !       must remain coupled.
        !
        entotal = 0.d0
        do ispec=1,dust_nr_species
           entotal  = entotal + mc_cumulener(ispec,ray_index)
        enddo
        tempav = compute_dusttemp_coupled_bd(dust_nr_species,           &
                 entotal,dustdens(1,ray_index),cellvolume(ray_index))
        do ispec=1,dust_nr_species
           dusttemp(ispec,ray_index) = tempav
        enddo
     endif
     !
     ! Do some testing, if required
     !
!     if(debug_check_all.eq.1) then
!        do ispec=1,dust_nr_species
!           if(number_invalid(dusttemp(ispec,ray_index)).ne.0) then
!              if(number_invalid(dusttemp(ispec,ray_index)).eq.1) then
!                 write(stdo,*) 'ERROR: Discovered INF dust temperature at ', &
!                      'ispec=',ispec,' index=',ray_index
!                 stop
!              else
!                 write(stdo,*) 'ERROR: Discovered NAN dust temperature at ', &
!                      'ispec=',ispec,' index=',ray_index
!                 stop
!              endif
!           endif
!        enddo
!     endif
  endif
  !
  ! Now pick a new frequency, according to dB_nu(T)/dT
  !
  ! NOTE: This must be mc_enerpart and not mc_cumulener, because
  !       suppose we have two dust species, with spec 2 having 0 opacity
  !       in the optical. If an optical photon is absorbed, it must
  !       have been by spec 1, and thus also spec 1 must re-emit!
  !
  call pick_randomfreq_db(dust_nr_species,dusttemp(1,ray_index), &
                          mc_enerpart,ray_inu)
  !
  ! If the temperature(s) and the freqdist have been recomputed,
  ! then we must reset the backup cumulative energy to the present
  ! one
  !
  if(calctemp.and.(params%ifast.ne.0)) then
     do ispec=1,dust_nr_species
        mc_cumulener_bk(ispec,ray_index) = mc_cumulener(ispec,ray_index)
     enddo
  endif
  !
  ! OpenMP Parallellization: Release lock on this cell
  !
  !$ if(ray_index .ge. 1 )then
  !$    call omp_unset_lock(lock(ray_index));
  !$ endif
  !
  ! Now reset the iqactive, since from now on all the events for
  ! this photon package are thermal (i.e. no quantum-heating
  ! events anymore). ONLY IN ORIGINAL PAH METHOD. NEW QUANTUM: 28.04.06
  !
  if(iqactive.eq.1) then 
     iqactive = 0
  endif
  !
  ! Done...
  !
end subroutine do_absorption_event



!--------------------------------------------------------------------------
!                          Random direction
!--------------------------------------------------------------------------
subroutine montecarlo_randomdir(dirx,diry,dirz)
  implicit none
  doubleprecision :: dirx,diry,dirz,linv,l2
  !
  ! Do a loop: only accept directions that are within unit circle
  !
  l2 = 2.d0
  do while(l2.gt.1.d0)
     dirx = 2.d0*ran2(iseed)-1.d0
     diry = 2.d0*ran2(iseed)-1.d0
     dirz = 2.d0*ran2(iseed)-1.d0
     l2   = dirx*dirx + diry*diry + dirz*dirz
     if(l2.lt.1.d-4) l2=2.0
  enddo
  linv = 1.d0 / sqrt(l2)
  !
  ! Now normalize to unity
  !
  dirx = dirx * linv
  diry = diry * linv
  dirz = dirz * linv
  !
end subroutine montecarlo_randomdir



!--------------------------------------------------------------------------
!                           Rotate vector 
!
! Rotate a vector such that the unit vector (1,0,0) is rotated in the
! direction of (dirx,diry,dirz). Note that this does not uniquely specify
! the rotation of the vector, as it is still possible to start with a
! rotation around the x-axis. But as long as that rotation is anyway
! randomized, this is not a problem.
! **** I checked this on 26.04.09, but perhaps check it again some time ***
! **** I tested it in an IDL program on 02.01.13: works fine ****
! --------------------------------------------------------------------------
subroutine montecarlo_rotatevec(vecx,vecy,vecz,dirx,diry,dirz)
  implicit none
  doubleprecision :: vecx,vecy,vecz,dirx,diry,dirz,l
  doubleprecision :: vx,vy,vz,dx,dy
  !
  ! First a rotation around the y-axis
  !
  l    = sqrt(dirx**2+diry**2)
  vx   = l    * vecx - dirz * vecz
  vy   = vecy
  vz   = dirz * vecx + l   * vecz
  !
  ! Then a rotation around the z-axis
  !
  if(l.gt.1d-10) then
     dx   = dirx / l
     dy   = diry / l
     vecx = dx * vx  - dy * vy
     vecy = dy * vx  + dx * vy
     vecz = vz
  else
     vecx = vx
     vecy = vy
     vecz = vz
  endif
  !
end subroutine montecarlo_rotatevec



!--------------------------------------------------------------------------
!                    Find random position in a cell
!--------------------------------------------------------------------------
subroutine montecarlo_random_pos_in_cell(icell)
  implicit none
  integer :: icell,iddr
  double precision :: rn, dum,r,theta,phi,cost,sint,cosp,sinp
  doubleprecision :: axi(1:2,1:3)
  integer :: ix,iy,iz
  !
  if(igrid_type.lt.100) then
     !
     ! The rectangular AMR grid
     !
     ! Get the cell walls
     !
     if(amr_tree_present) then
        do iddr=1,3
           axi(1,iddr) = amr_finegrid_xi(amr_theleafs(icell)%link%ixyzf(iddr),iddr,amr_theleafs(icell)%link%level)
           axi(2,iddr) = amr_finegrid_xi(amr_theleafs(icell)%link%ixyzf(iddr)+1,iddr,amr_theleafs(icell)%link%level)
        enddo
     else
        call amr_regular_get_ixyz(icell,ix,iy,iz)
        axi(1,1) = amr_finegrid_xi(ix,1,0)
        axi(2,1) = amr_finegrid_xi(ix+1,1,0)
        axi(1,2) = amr_finegrid_xi(iy,2,0)
        axi(2,2) = amr_finegrid_xi(iy+1,2,0)
        axi(1,3) = amr_finegrid_xi(iz,3,0)
        axi(2,3) = amr_finegrid_xi(iz+1,3,0)
     endif
     !
     ! Now do thing different for different coordinate systems
     !
     if(igrid_coord.lt.100) then
        !
        ! Cartesian coordinates
        !
        rn         = ran2(iseed)
        dum        = axi(2,1) - axi(1,1)
        ray_cart_x = axi(1,1) + dum * rn
        rn         = ran2(iseed)
        dum        = axi(2,2) - axi(1,2)
        ray_cart_y = axi(1,2) + dum * rn
        rn         = ran2(iseed)
        dum        = axi(2,3) - axi(1,3)
        ray_cart_z = axi(1,3) + dum * rn
        ! BUGFIX 1-D plane-parallel and 2-D pencil-parallel modes 19.06.2021
        ! Must set x and y to 0, otherwise the ray segment length cannot be
        ! properly computed, because x and y become of order 1e90 large.
        ! It has no consequences for the normal RADMC-3D functionality, only
        ! for the 1-D plane-parallel and 2-D pencil-parallel modes.
        ! Thanks, Til Birnstiel, for spotting and reporting the error.
        if(igrid_coord.eq.10) then
           ray_cart_x = 0.d0
           ray_cart_y = 0.d0
        endif
        if(igrid_coord.eq.20) then
           ray_cart_x = 0.d0
        endif
     elseif(igrid_coord.lt.200) then
        !
        ! Spherical coordinates
        ! 
        ! NOTE: We are going to simplify this here now. In reality one
        !       should account for the curved geometry of the cell. But
        !       that would be costly to compute and for small enough
        !       cells this is not too important.
        !
        rn         = ran2(iseed)
        dum        = axi(2,1) - axi(1,1)
        r          = axi(1,1) + dum * rn
        rn         = ran2(iseed)
        dum        = axi(2,2) - axi(1,2)
        theta      = axi(1,2) + dum * rn
        rn         = ran2(iseed)
        dum        = axi(2,3) - axi(1,3)
        phi        = axi(1,3) + dum * rn
        cost       = cos(theta)
        sint       = sin(theta)
        cosp       = cos(phi)
        sinp       = sin(phi)
        ray_cart_x = r*sint*cosp
        ray_cart_y = r*sint*sinp
        ray_cart_z = r*cost
     else
        stop 4955
     endif
     !
  elseif(igrid_type.eq.101) then
     !
     ! The Delaunay triangulated grid
     !
     write(stdo,*) 'NOT YET READY: Delaunay unstructured grid'
     stop
  elseif(igrid_type.eq.201) then
     !
     ! The Voronoi grid
     !
     ! Here we assume that the photon is emitted from the center of the
     ! cell (**** APPROXIMATION ****), because for this gridding mode the
     ! cell shape is too complex to easily find a correct random position.
     ! 
     ! **** TO BE IMPROVED ****
     !
     write(stdo,*) 'NOT YET READY: Voronoi unstructured grid'
     stop
  else
     write(stdo,*) 'ERROR: Grid type ',igrid_type,' not known.'
     stop
  endif
end subroutine montecarlo_random_pos_in_cell


!--------------------------------------------------------------------------
!              HENYEY-GREENSTEIN RANDOM SCATTERING
!
! This routine gives the random mu_scat and phi_scat according
! to the small-angle-scattering formula of Henyey-Greenstein.
!
!            1         1 - g^2
!    P(mu) = - -------------------------
!            2 (1 + g^2 - 2 g mu )^(3/2)
!
! We do this by introducing a variable xi
!
!            1 1 - g  /        1 + g                  \
!    xi   := - -----  | ------------------------- - 1 |
!            2   g    \ (1 + g^2 - 2 g mu )^(1/2)     /
!
! Tested on 02.01.2013 with the IDL program test_henyeygreenstein.pro
! and compared to henyeygreenstein_phasefunc() : it works excellently.
!--------------------------------------------------------------------------
subroutine henyeygreenstein_randomdir(g,dirx,diry,dirz)
  implicit none
  doubleprecision g,dirx,diry,dirz,xi,g2,l2,linv
  !
  if(g.ne.0.d0) then
     g2   = g*g
     xi   = ran2(iseed)
     dirx = (0.5d0/g)*(1.+g2-((1.-g2)/(1.-g+2*g*xi))**2)
  else
     dirx = 2.d0*ran2(iseed)-1.d0
  endif
  l2   = 2.d0
  do while(l2.gt.1.d0)
     diry = 2.d0*ran2(iseed)-1.d0
     dirz = 2.d0*ran2(iseed)-1.d0
     l2   = diry*diry + dirz*dirz
     if(l2.lt.1.d-4) l2=2.0
  enddo
  linv = sqrt((1.d0-dirx**2)/l2)
  !
  ! Now normalize to sqrt(1-mu^2)  where mu = dirx
  !
  diry = diry * linv
  dirz = dirz * linv
  !
  ! Safety check
  !
  if(debug_check_all.eq.1) then
     if(abs(dirx**2+diry**2+dirz**2-1.d0).gt.1d-14) then
        write(stdo,*) 'INTERNAL ERROR: Henyey-Greenstein direction vector not OK.'
        stop
     endif
  endif
  !
end subroutine henyeygreenstein_randomdir


!--------------------------------------------------------------------------
!        HENYEY GREENSTEIN PHASE FUNCTION FOR ANISOTROPIC SCATTERING
!
! The phase function Phi(mu) is twice the mu-probability function P(mu),
! because the phase function Phi(mu) must obey:
!
!   1    /                  1 /+1
!  ---   0 Phi(mu) dOmega = - |   Phi(mu) dmu = 1
!  4*pi  /                  2 /-1
!
! while P(mu) is normalized as 
!
!   /+1
!   |   P(mu) dmu = 1
!   /-1
!
! Hence
!
!                               1 - g^2
!   Phi(mu) = 2 P(mu) = -------------------------
!                       (1 + g^2 - 2 g mu )^(3/2)
!
! --------------------------------------------------------------------------
function henyeygreenstein_phasefunc(g,costheta)
  implicit none
  double precision :: g,g2,costheta,henyeygreenstein_phasefunc
  g2 = g*g
  henyeygreenstein_phasefunc = (1.d0-g2)/(1.d0-2*g*costheta+g2)**1.5
  return
end function henyeygreenstein_phasefunc


!--------------------------------------------------------------------------
!              ANISOTROPIC SCATTERING RANDOM DIRECTION
!
! Similar to henyeygreenstein_randomdir(), but now using the tabulated
! scattering matrix phase function.
!
! Tested on 02.01.2013 with the IDL program test_anisoscat.pro
! and compared to henyeygreenstein_phasefunc() : it works excellently.
!--------------------------------------------------------------------------
subroutine anisoscat_randomdir(inu,ispec,dirx,diry,dirz)
  implicit none
  integer :: inu,ispec,imu
  doubleprecision :: dirx,diry,dirz
  doubleprecision :: rnumber,eps,mu,phi,mu1
  !
  ! Get a random number
  !
  rnumber = ran2(iseed)*zcumul(scat_munr,inu,1,ispec)
  !
  ! Find it in the cumulative array --> mu
  !
  call hunt(zcumul(:,inu,1,ispec),scat_munr,rnumber,imu)
  if((imu.lt.1).or.(imu.ge.scat_munr)) stop 4777
  eps = (rnumber-zcumul(imu,inu,1,ispec)) /                  &
        (zcumul(imu+1,inu,1,ispec)-zcumul(imu,inu,1,ispec))
  if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 4778
  mu = (1.d0-eps)*scat_mui_grid(imu)+eps*scat_mui_grid(imu+1)
  mu1 = sqrt(1.d0-mu*mu)
  !
  ! Find a random phi
  !
  phi = twopi * ran2(iseed)
  !
  ! Make the direction vector
  !
  dirx = mu
  diry = mu1 * cos(phi)
  dirz = mu1 * sin(phi)
  !
end subroutine anisoscat_randomdir


!--------------------------------------------------------------------------
!            PHASE FUNCTION FOR ANISOTROPIC SCATTERING
! Similar to henyeygreenstein_phasefunc(), but now using the tabulated
! scattering matrix phase function.
!--------------------------------------------------------------------------
function anisoscat_phasefunc(inu,ispec,mu)
  implicit none
  integer :: inu,ispec
  double precision :: mu,anisoscat_phasefunc
  integer :: imu
  double precision :: eps
  call hunt(scat_mui_grid(:),scat_munr,mu,imu)
  if(imu.lt.1) imu=1
  if(imu.ge.scat_munr) imu=scat_munr-1
  eps = (mu-scat_mui_grid(imu))/(scat_mui_grid(imu+1)-scat_mui_grid(imu))
  if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 3702
  anisoscat_phasefunc = (1.d0-eps)*zmatrix(imu,inu,1,ispec) + &
                               eps*zmatrix(imu+1,inu,1,ispec)
  ! BUGFIX 06.12.2016
  anisoscat_phasefunc = anisoscat_phasefunc * fourpi / kappa_s(inu,ispec)
  return
end function anisoscat_phasefunc


!--------------------------------------------------------------------------
!    FOR ALIGNED GRAINS: RANDOM PHOTON EMISSION DIRECTION AND POL STATE
!--------------------------------------------------------------------------
subroutine montecarlo_aligned_randomphot(index,inu,ener,pkg)
  implicit none
  integer :: index,inu
  double precision :: ener
  type(photon) :: pkg
  integer :: ispec,imu
  double precision :: rn
  !
  ! First find out which grain emitted this photon
  ! 
  alphacum(1) = 0.d0
  do ispec=1,dust_nr_species
     alphacum(ispec+1) = alphacum(ispec) + dustdens(ispec,index) * &
                         kappa_a(inu,ispec)
  enddo
  if(alphacum(dust_nr_species+1).le.0.d0) stop 2244
  do ispec=1,dust_nr_species
     alphacum(ispec) = alphacum(ispec) / alphacum(dust_nr_species+1)
  enddo
  alphacum(dust_nr_species+1) = 1.d0
  rn = ran2(iseed)
  call hunt(alphacum,dust_nr_species+1,rn,ispec)
  if((ispec.lt.1).or.(ispec.gt.dust_nr_species)) stop 7266
  !
  ! Now find out which direction the photon goes and
  ! which polarization state it has
  !
  call polarization_random_aligned_thermemis(pkg,mc_align_munr,mc_align_mu, &
                         mc_align_orth(:,inu,ispec),                        &
                         mc_align_para(:,inu,ispec),                        &
                         mc_align_opcumul(:,inu,ispec),                     &
                         grainalign_dir(:,index),grainalign_eff(index))
  !
  ! Now make the dimensionless polarization state 
  ! dimension-full
  !
  pkg%E    = ener*pkg%E
  pkg%Q    = ener*pkg%Q
  pkg%U    = ener*pkg%U
  pkg%V    = ener*pkg%V
  !
end subroutine montecarlo_aligned_randomphot


!--------------------------------------------------------------------------
!                  WRITE MEAN INTENSITY TO FILE
!--------------------------------------------------------------------------
subroutine write_meanint_to_file()
  implicit none
  integer :: icell,index,inu,i,ierr,precis
  integer(kind=8) :: nn,kk
  logical :: fex
  double precision, allocatable :: data(:)
  !
  ! Determine the precision
  !
  if(rto_single) then
     precis = 4
  else
     precis = 8
  endif
  !
  ! Now write the dust temperature
  !
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
     ! Open file and write the mean intensity to it
     !
     if(rto_style.eq.1) then
        !
        ! Write the mean intensity in ascii form
        !
        ! NOTE: The new format is "2", and includes a list of frequencies
        !
        open(unit=1,file='mean_intensity.out')
        write(1,*) 2                                   ! Format number
        write(1,*) nrcellsinp
        write(1,*) mc_nrfreq
        write(1,*) (mc_frequencies(inu),inu=1,mc_nrfreq)
     elseif(rto_style.eq.2) then
        !
        ! Write the mean intensity in f77-style unformatted form,
        ! using a record length given by rto_reclen
        !
        ! NOTE: The new format is "2", and includes a list of frequencies
        !
        open(unit=1,file='mean_intensity.uout',form='unformatted')
        nn = 2
        kk = rto_reclen
        write(1) nn,kk               ! Format number and record length
        nn = nrcellsinp
        kk = mc_nrfreq
        write(1) nn,kk
        write(1) (mc_frequencies(inu),inu=1,mc_nrfreq)
     elseif(rto_style.eq.3) then
        !
        ! C-compliant binary
        !
        ! NOTE: The new format is "2", and includes a list of frequencies
        !
        open(unit=1,file='mean_intensity.bout',status='replace',access='stream')
        nn = 2
        kk = precis
        write(1) nn,kk               ! Format number and precision
        nn = nrcellsinp
        kk = mc_nrfreq
        write(1) nn,kk
        write(1) (mc_frequencies(inu),inu=1,mc_nrfreq)
     else
        write(stdo,*) 'ERROR: Do not know I/O style ',rto_style
        stop
     endif
     !
     ! Now write the mean intensity one wavelength at a time 
     !
     do inu=1,mc_nrfreq
        call write_scalarfield(1,rto_style,precis,nrcellsinp, &
             mc_nrfreq,1,inu,1,rto_reclen,                    &
             scalar1=mcscat_meanint)
     enddo
     !
     ! Close
     !
     close(1)
  else
     !
     ! Other grids not yet implemented
     !
     write(stdo,*) 'ERROR: Only regular and AMR grids implemented'
     stop
  endif
  !
  ! If the grid is internally made, then we must make sure that
  ! the grid has been written to file, otherwise the output file
  ! created here makes no sense.
  !
  if((.not.grid_was_read_from_file).and.(.not.grid_was_written_to_file)) then
     call write_grid_file()
     grid_was_written_to_file = .true.     ! Avoid multiple writings
  endif
end subroutine write_meanint_to_file



!--------------------------------------------------------------------------
!                  WRITE PHOTON STATISTICS ARRAY TO FILE
!--------------------------------------------------------------------------
subroutine write_photon_statistics_to_file()
  implicit none
  integer :: icell,index,inu,i,ierr,precis
  integer(kind=8) :: nn,kk
  logical :: fex
  double precision, allocatable :: data(:)
  !
  ! Check
  !
  if(.not.allocated(mc_iphotcount)) then
     write(stdo,*) 'ERROR: Attempting to write photon statistics file'
     write(stdo,*) '       but array mc_iphotcount is not allocated.'
     stop
  endif
  !
  ! Determine the precision
  !
  if(rto_single) then
     precis = 4
  else
     precis = 8
  endif
  !
  ! Now write the dust temperature
  !
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
     ! Open file and write the mean intensity to it
     !
     if(rto_style.eq.1) then
        !
        ! Write the photon statistics in ascii form
        !
        ! NOTE: The new format is "2", and includes a list of frequencies
        !
        open(unit=1,file='photon_statistics.out')
        write(1,*) 1                                   ! Format number
        write(1,*) nrcellsinp
     elseif(rto_style.eq.2) then
        !
        ! Write the mean intensity in f77-style unformatted form,
        ! using a record length given by rto_reclen
        !
        ! NOTE: The new format is "2", and includes a list of frequencies
        !
        open(unit=1,file='photon_statistics.uout',form='unformatted')
        nn = 1
        kk = rto_reclen
        write(1) nn,kk               ! Format number and record length
        nn = nrcellsinp
        write(1) nn
     elseif(rto_style.eq.3) then
        !
        ! C-compliant binary
        !
        ! NOTE: The new format is "2", and includes a list of frequencies
        !
        open(unit=1,file='photon_statistics.bout',status='replace',access='stream')
        nn = 1
        kk = precis
        write(1) nn,kk               ! Format number and precision
        nn = nrcellsinp
        write(1) nn
     else
        write(stdo,*) 'ERROR: Do not know I/O style ',rto_style
        stop
     endif
     !
     ! Now write the photon statistics
     !
     call write_scalarfield(1,rto_style,precis,nrcellsinp, &
             1,1,1,1,rto_reclen,scalar0=mc_iphotcount)
     !
     ! Close
     !
     close(1)
  else
     !
     ! Other grids not yet implemented
     !
     write(stdo,*) 'ERROR: Only regular and AMR grids implemented'
     stop
  endif
  !
  ! If the grid is internally made, then we must make sure that
  ! the grid has been written to file, otherwise the output file
  ! created here makes no sense.
  !
  if((.not.grid_was_read_from_file).and.(.not.grid_was_written_to_file)) then
     call write_grid_file()
     grid_was_written_to_file = .true.     ! Avoid multiple writings
  endif
end subroutine write_photon_statistics_to_file




!===========================================================================
!              ROUTINES FOR STORING TEMPERATURE ARRAYS
!===========================================================================


!-------------------------------------------------------------------
!               HELPER FUNCTION FOR THE ABSORPTION EVENT
!
!     This function returns 0 for that temperature where the new
!     spectrum of the cell equals the old one plus the energy of
!     the photonj package.
!-------------------------------------------------------------------
function absevfunc(temp)
  implicit none
  doubleprecision absevfunc,temp
  !
  integer inu
  doubleprecision demis
  !
  demis = 0.d0
  do inu=1,freq_nr
     fnu_diff(inu) = cellalpha(inu) * bplanck(temp,freq_nu(inu))
     demis         = demis + fnu_diff(inu) * freq_dnu(inu)
  enddo
  demis = demis * fourpi 
  absevfunc = demis 
  return
end function absevfunc


!----------------------------------------------------------------------
!             COMPUTE THE ENERGY BELONGING TO A TEMPERATURE
!
!     This routine computes the energy output per gram of dust
!     of a given species and size for an array of temperatures.
!----------------------------------------------------------------------
subroutine make_emiss_dbase(ntemp,temp0,temp1)
  implicit none
  integer :: ntemp
  double precision :: temp0,temp1
  !
  integer :: itemp,ispec,ierr,inu
  !
  ! Check
  !
  if((freq_nr.lt.1).or.(dust_nr_species.lt.1)) then
     write(stdo,*) 'INTERNAL ERROR IN MONTECARLO: Could not make temp grid'
     stop
  endif
  !
  ! Allocate arrays
  !
  allocate(db_temp(ntemp),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate db_temp'
     stop 
  endif
  allocate(db_enertemp(ntemp,dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate db_temp'
     stop 
  endif
  allocate(db_logenertemp(ntemp,dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate db_temp'
     stop 
  endif
  allocate(db_emiss(freq_nr,ntemp,dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate db_emiss'
     stop 
  endif
  allocate(db_cumulnorm(freq_nr+1,ntemp,dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate db_cumulnorm'
     stop 
  endif
  !$OMP PARALLEL
  allocate(db_cumul(freq_nr+1),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate db_cumul'
     stop 
  endif
  !$OMP END PARALLEL
  !
  ! Allocate array used internally for picking random frequency
  !
  !$OMP PARALLEL
  allocate(enercum(dust_nr_species+1),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate enercum'
     stop 
  endif
  !$OMP END PARALLEL
  !
  ! Allocate temporary arrays
  !
  allocate(cellalpha(freq_nr),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate cellalpha'
     stop 
  endif
  allocate(fnu_diff(freq_nr),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate fnu_diff'
     stop 
  endif
  allocate(diffemis(freq_nr),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate diffemis'
     stop 
  endif
  !
  ! Make temperature grid
  !
  db_ntemp = ntemp
  do itemp=1,ntemp
     db_temp(itemp) = temp0 * (temp1/temp0)**((itemp-1.d0)/(ntemp-1.d0))
  enddo
  !     
  ! A loop over all dust species
  !
  do ispec=1,dust_nr_species
     !
     ! Compute the emissivity and the total emitted energy
     ! per second per gram dust of this species/size.
     !
     do inu=1,freq_nr
        cellalpha(inu) = kappa_a(inu,ispec)
     enddo
     do itemp=1,db_ntemp
        db_enertemp(itemp,ispec) = absevfunc(db_temp(itemp))
        db_logenertemp(itemp,ispec) = log(db_enertemp(itemp,ispec))
        do inu=1,freq_nr
           db_emiss(inu,itemp,ispec) = fnu_diff(inu)
        enddo
     enddo
     !
     ! Now compute the derivative of the emissivity
     ! with respect to energy, so that we can obtain
     ! a right distribution function for picking a
     ! new photon frequency. We immediately make
     ! a cumulative version of this for each
     ! species, so that the picking of a new frequency
     ! for a photon goes in two steps:
     !  1) Pick a species/size from the energies
     !  2) Pick, within this selected species, a freq
     !
     ! NOTE!!!! VERY IMPORTANT!!!: 
     !   This is done slightly differently from the
     !   original method in which the exact difference
     !   of fnu and fnu_old was taken. But this is not
     !   too impportant, because 1) the normalization
     !   if irrelevant and 2) for sufficient photon
     !   statistics the derivative is close enough to
     !   the difference that it is OK.
     !
     ! For the first temperature we simply take the
     ! emissivity function itself
     !
     itemp = 1
     do inu=1,freq_nr
        diffemis(inu) = db_emiss(inu,itemp,ispec)
     enddo
     db_cumul(1) = 0.d0
     do inu=1,freq_nr
        db_cumul(inu+1) = db_cumul(inu) + diffemis(inu) * freq_dnu(inu)
     enddo
     do inu=1,freq_nr+1
        db_cumulnorm(inu,itemp,ispec) = db_cumul(inu) / db_cumul(freq_nr+1)
     enddo
     !
     ! For the rest take indeed the difference
     !              
     do itemp=2,db_ntemp
        do inu=1,freq_nr
           diffemis(inu) = db_emiss(inu,itemp,ispec) -                      &
                           db_emiss(inu,itemp-1,ispec) 
        enddo
        db_cumul(1) = 0.d0
        do inu=1,freq_nr
           db_cumul(inu+1) = db_cumul(inu) + diffemis(inu) * freq_dnu(inu)
        enddo
        do inu=1,freq_nr+1
           db_cumulnorm(inu,itemp,ispec) = db_cumul(inu) / db_cumul(freq_nr+1)
        enddo
     enddo
  enddo
  !
  ! Free temporary arrays
  !
  deallocate(cellalpha,fnu_diff,diffemis)
  !
end subroutine make_emiss_dbase


!----------------------------------------------------------------------
!                  FREE TEMPERATURE DATABASE
!----------------------------------------------------------------------
subroutine free_emiss_dbase()
  implicit none
  db_ntemp=0
  if(allocated(db_temp)) deallocate(db_temp)
  if(allocated(db_cumulnorm)) deallocate(db_cumulnorm)
  if(allocated(db_enertemp)) deallocate(db_enertemp)
  if(allocated(db_logenertemp)) deallocate(db_logenertemp)
  if(allocated(db_emiss)) deallocate(db_emiss)
  !$OMP PARALLEL
  if(allocated(db_cumul)) deallocate(db_cumul)
  if(allocated(enercum)) deallocate(enercum)
  !$OMP END PARALLEL
end subroutine free_emiss_dbase



!----------------------------------------------------------------------
!               COMPUTE THE TEMPERATURE FOR A GIVEN ENERGY
!
!     This routine computes the dust temperature for a given 
!     energy-per-gram-of-dust for this species.
!----------------------------------------------------------------------
function compute_dusttemp_energy_bd(ener,ispec)
  implicit none
  !
  integer :: ispec
  doubleprecision :: ener,compute_dusttemp_energy_bd,logener
  !
  integer :: itemp
  doubleprecision :: eps,temp
  !
  ! First find the energy in the grid
  !
  ! Linear method:
  !
  !!!!!! call hunt(db_enertemp(1,ispec),db_ntemp,ener,itemp)
  !
  ! Logarithmic method:
  !
  logener = log(ener)
  call hunt(db_logenertemp(1,ispec),db_ntemp,logener,itemp)
  !
  ! Check if we are in range
  !
  if(itemp.ge.db_ntemp) then
     write(stdo,*) 'ERROR: Too high temperature discovered' 
     stop 76823
  endif
  if(itemp.le.0) then
     !
     ! Temperature presumably below lowest temp in dbase
     !
     eps = ener/db_enertemp(1,ispec)
     if(eps.gt.1.d0) stop 9911
     temp = eps*db_temp(1)
  else
     !
     ! Temperature in the right range: default action
     !
     eps = (ener-db_enertemp(itemp,ispec)) /                            &
           (db_enertemp(itemp+1,ispec) - db_enertemp(itemp,ispec))
     if((eps.gt.1.d0).or.(eps.lt.0.d0)) stop 9912
     temp = (1.d0-eps)*db_temp(itemp) + eps*db_temp(itemp+1)
  endif
  !
  ! Return
  !
  compute_dusttemp_energy_bd = temp
  return
  !
end function compute_dusttemp_energy_bd


!----------------------------------------------------------------------
!               COMPUTE THE TEMPERATURE FOR A GIVEN ENERGY
!
!     This routine computes the dust temperature for a given 
!     energy-per-gram-of-dust for this species.
!----------------------------------------------------------------------
function compute_dusttemp_coupled_bd(nspec,entotal,rho,vol)
  implicit none
  !
  integer :: nspec
  doubleprecision :: vol,entotal
  doubleprecision :: rho(nspec)
  doubleprecision :: compute_dusttemp_coupled_bd
  !
  integer :: itemp,ispec
  doubleprecision :: eps,temp,entot,dummy1,dummy2
  doubleprecision :: en1,en2
  !
  ! Compute the total energy per gram dust
  !
  entot = entotal / vol
  !
  ! First find the energy in the grid
  !
  ! NOTE: It is difficult to use the logarithmic energy search as we
  !       can easily do for the decoupled temperatures, because we do
  !       not yet know the composition of the dust globally. One could
  !       adapt hunt_temp but that would require many evaluations of
  !       a log, and that might kill all the advantage of the log
  !       search. So we will keep this linear.
  !
  call hunt_temp(rho,entot,itemp)
  !
  ! Check if we are in range
  !
  if(itemp.ge.db_ntemp) then
     write(stdo,*) 'ERROR: Too high temperature discovered' 
     stop 76823
  endif
  if(itemp.le.0) then
     !
     ! Temperature presumably below lowest temp in dbase
     !
     en2 = 0.d0
     do ispec=1,dust_nr_species
        en2 = en2 + rho(ispec) * db_enertemp(1,ispec)
     enddo
     eps = entot/en2
     if(eps.gt.1.d0) stop 59911
     temp = eps*db_temp(1)
  else
     !
     ! Temperature in the right range: default action
     !
     en1 = 0.d0
     do ispec=1,dust_nr_species
        en1 = en1 + rho(ispec)*db_enertemp(itemp,ispec)
     enddo
     en2 = 0.d0
     do ispec=1,dust_nr_species
        en2 = en2 + rho(ispec)*db_enertemp(itemp+1,ispec)
     enddo
     eps = (entot-en1) / (en2-en1)
     if((eps.gt.1.d0).or.(eps.lt.0.d0)) stop 99124
     temp = (1.d0-eps)*db_temp(itemp) + eps*db_temp(itemp+1)
  endif
  !
  ! Return
  !
  compute_dusttemp_coupled_bd = temp
  return
  !
end function compute_dusttemp_coupled_bd


!----------------------------------------------------------------------
!               COMPUTE NEW IFREQ FROM LOCAL SITUATION
!----------------------------------------------------------------------
subroutine pick_randomfreq_db(nspec,temp,mc_enerpart,inupick)
  implicit none
  !
  integer :: nspec,inupick
  doubleprecision :: temp(nspec)
  !
  double precision :: temp_local(nspec)
  !
  doubleprecision :: mc_enerpart(nspec)
  !
  doubleprecision rn,eps
  integer inu,ispec,itemp
  !
  ! Figure out which species to take. This can 
  ! be determined from the energies
  !
  enercum(1) = 0.d0
  do ispec=1,dust_nr_species
     enercum(ispec+1) = enercum(ispec) + mc_enerpart(ispec)
  enddo
  if(dust_nr_species.gt.1) then 
     rn = ran2(iseed)*enercum(dust_nr_species+1)
!    BUGFIX by Seokho Lee 24.02.2015:
!     call hunt(enercum,dust_nr_species,rn,ispec)
     call hunt(enercum,dust_nr_species+1,rn,ispec)
     if((ispec.lt.1).or.(ispec.gt.dust_nr_species)) stop 50209
  else
     ispec = 1
  endif
  !
  temp_local(1:nspec) = temp(1:nspec)
  !
  ! Find the temperature in the database  !
  !
  call hunt(db_temp,db_ntemp,temp_local(ispec),itemp)
  !
  ! Check if we are in range
  !
  if(itemp.ge.db_ntemp) then
     write(stdo,*) 'ERROR: Too high temperature discovered' 
     stop 77823
  endif
  if(itemp.le.0) then
     !
     ! Temperature presumably below lowest temp in dbase
     !
     itemp = 1
     eps   = 0.d0
  else
     !
     ! Temperature in the right range: default action
     !
     eps = (temp_local(ispec)-db_temp(itemp)) /                             &
           (db_temp(itemp+1)-db_temp(itemp))
     if((eps.gt.1.d0).or.(eps.lt.0.d0)) stop 94941
  endif
  !
  ! Now, within this species/size, find the frequency
  !
  if(intplt.eq.1) then
     !
     ! Do this properly
     !
     do inu=1,freq_nr+1
        db_cumul(inu) = (1.d0-eps)*db_cumulnorm(inu,itemp,ispec) +   &
                        eps*db_cumulnorm(inu,itemp+1,ispec)
     enddo
     if(abs(db_cumul(freq_nr+1)-1.d0).gt.1d-11) stop 43465
     rn = ran2(iseed)
     call hunt(db_cumul,freq_nr,rn,inupick)
     if((inupick.lt.1).or.(inupick.gt.freq_nr)) stop 8189
  else
     !
     ! Do this approximatively
     !
     if(eps.gt.0.5d0) then
        if(itemp.lt.db_ntemp-1) itemp=itemp+1
     endif
     rn = ran2(iseed)
     call hunt(db_cumulnorm(1,itemp,ispec),freq_nr,rn,inupick)
  endif
  !
end subroutine pick_randomfreq_db




!------------------------------------------------------------------
!              HUNT IN TOTAL ENERGY FOR TEMPERATURE
!     (MODIFIED NRECIP ROUTINE)
!        energy = total energy divided by volume
!------------------------------------------------------------------
subroutine hunt_temp(rho,energy,jlo)
  integer :: jlo,n
  doubleprecision :: energy,xx,x
  doubleprecision :: rho(dust_nr_species)
  integer :: inc,jhi,jm
  logical :: ascnd
  !
  n = db_ntemp
  x = energy
  ascnd=.true.
  if(jlo.le.0.or.jlo.gt.n)then
     jlo=0
     jhi=n+1
     goto 3
  endif
  inc=1
  xx=0.d0
  do ispec=1,dust_nr_species
     xx = xx + rho(ispec) * db_enertemp(jlo,ispec)
  enddo
  if(x.ge.xx.eqv.ascnd)then
1    jhi=jlo+inc
     xx=0.d0
     do ispec=1,dust_nr_species
        xx = xx + rho(ispec) * db_enertemp(jhi,ispec)
     enddo
     if(jhi.gt.n)then
        jhi=n+1
     else if(x.ge.xx.eqv.ascnd)then
        jlo=jhi
        inc=inc+inc
        goto 1
     endif
  else
     jhi=jlo
2    jlo=jhi-inc
     xx=0.d0
     do ispec=1,dust_nr_species
        xx = xx + rho(ispec) * db_enertemp(jlo,ispec)
     enddo
     if(jlo.lt.1)then
        jlo=0
     else if(x.lt.xx.eqv.ascnd)then
        jhi=jlo
        inc=inc+inc
        goto 2
     endif
  endif
3 if(jhi-jlo.eq.1)return
  jm=(jhi+jlo)/2
  xx=0.d0
  do ispec=1,dust_nr_species
     xx = xx + rho(ispec) * db_enertemp(jm,ispec)
  enddo
  if(x.gt.xx.eqv.ascnd)then
     jlo=jm
  else
     jhi=jm
  endif
  goto 3
end subroutine hunt_temp
!  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


!=========================================================================
!    BELOW ARE THE SUBROUTINES FOR THE MODIFIED RANDOM WALK METHOD
!
! I chose to make the MRW subroutines an entirely separate submodule,
! independent of stuff that was already done in RADMC-3D (except for the
! subroutines computing the opacities). These subroutines have been tested
! separately. Note that they are designed for regular grid cells only. So
! if, in the future, RADMC-3D will get unstructured grids, then the routines
! below cannot be used, and a separate version of them will have to be
! built.
! =========================================================================

!-------------------------------------------------------------------------
!                COMPUTE THE ROSSELAND MEAN OPACITIES 
!-------------------------------------------------------------------------
subroutine mrw_calculate_rossmean(nrspec,temp,dens,alpha)
  implicit none
  integer :: nrspec
  double precision :: temp,dens(1:nrspec),alpha
  integer :: ispec,inu
  double precision :: dum1,dum2,anu,dbdt
  !
  ! Trivial test
  !
  if(nrspec.ne.dust_nr_species) stop 4967
  !
  ! Reset
  !
  dum1 = 0.d0
  dum2 = 0.d0
  !
  ! Loop over frequencies
  !
  do inu=1,freq_nr
     !
     ! Compute the total opacity at this frequency
     !
     ! The factor (1-g) before the kappa_s is to account for
     ! the fact that for non-isotropic scattering the 
     ! diffusion opacity caused by scattering is reduced. 
     ! See Min et al. 2009, which took the formula from
     ! Ishimaru (1978).
     !
     ! Note: We use the kappa_a etc arrays which formally
     !       are evaluated on the mc_frequencies array.
     !       But for the Bjorkman & Wood method they 
     !       must be the same!
     !
     anu = 0.d0
     do ispec=1,dust_nr_species
        anu = anu + dens(ispec) * ( kappa_a(inu,ispec) +         &
               (1.d0-kappa_g(inu,ispec)) * kappa_s(inu,ispec) )
     enddo
     !
     ! Compute the dB_nu(T)/dT at this frequency
     !
     dbdt = bplanckdt(temp,freq_nu(inu))
     !
     ! Now add to the numerator and denominator
     !
     dum1 = dum1 + dbdt * freq_dnu(inu)
     dum2 = dum2 + dbdt * freq_dnu(inu) / ( anu + 1d-99 ) 
     !
  enddo
  !
  ! Now construct the Rosseland mean opacity
  !
  alpha = dum1 / dum2
  !
end subroutine mrw_calculate_rossmean

!-------------------------------------------------------------------------
!       COMPUTE THE PLANCK MEAN OPACITIES BELONGING TO A TEMPERATURE
!                     AND STORE THEM IN A DATABASE
!
! This routine pre-computes the Planck-mean opacities for each dust 
! species, for a table of temperatures. Note that we unfortunately cannot
! do this for the Rosseland mean, because that is a non-linear formula.
! But for the Planck mean this can be done beforehand, saving a lot of
! computational time.
!-------------------------------------------------------------------------
subroutine mrw_make_planckopac_dbase(ntemp,temp0,temp1)
  implicit none
  integer :: ntemp
  double precision :: temp0,temp1
  double precision, allocatable :: bnu(:,:)
  !
  integer :: itemp,ispec,ierr,inu
  double precision :: kap,kpl,dum,dnu
  !
  ! Check
  !
  if((freq_nr.lt.1).or.(dust_nr_species.lt.1)) then
     write(stdo,*) 'INTERNAL ERROR IN MRW MONTECARLO: Could not make temp grid'
     stop
  endif
  !
  ! Allocate arrays
  !
  allocate(mrw_db_temp(ntemp),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate mrw_db_temp'
     stop 
  endif
  allocate(mrw_db_kappa_abs_planck(ntemp,dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate mrw_db_kappa_abs_planck'
     stop 
  endif
  !
  ! Make temperature grid
  !
  mrw_db_ntemp = ntemp
  do itemp=1,ntemp
     mrw_db_temp(itemp) = temp0 * (temp1/temp0)**((itemp-1.d0)/(ntemp-1.d0))
  enddo
  !
  ! Prepare the B_nu(T) grid
  !
  allocate(bnu(1:freq_nr,1:ntemp))
  do itemp=1,ntemp
     do inu=1,freq_nr
        bnu(inu,itemp) = bplanck(mrw_db_temp(itemp),freq_nu(inu))
     enddo
  enddo
  !     
  ! A loop over all dust species
  !
  do ispec=1,dust_nr_species
     !
     ! Loop over temperatures
     !
     do itemp=1,ntemp
        !
        ! Reset
        !
        kpl = 0.d0
        dum = 0.d0
        !
        ! Loop over frequency
        !
        do inu=1,freq_nr
           kap = find_dust_kappa_interpol(freq_nu(inu),ispec,mrw_db_temp(itemp),1,0,0)
           kpl = kpl + kap * bnu(inu,itemp) * freq_dnu(inu)
           dum = dum + bnu(inu,itemp) * freq_dnu(inu)
        enddo
        kpl = kpl / dum
        !
        ! Store
        !
        mrw_db_kappa_abs_planck(itemp,ispec) = kpl
     enddo
  enddo
  !
  ! Finish
  !
  if(allocated(bnu)) deallocate(bnu)
end subroutine mrw_make_planckopac_dbase

!-------------------------------------------------------------------------
!            FIND THE PLANCK MEAN OPACITIES FROM THE DATABASE
! 
! Make sure to first make the database using mrw_make_planckopac_dbase()
!-------------------------------------------------------------------------
subroutine mrw_get_planckopac_dbase(temp,nspec,kappa_a_pl)
  implicit none
  integer :: nspec
  double precision :: temp
  double precision :: kappa_a_pl(1:nspec)
  double precision :: eps,eps1
  integer :: itemp,ispec
  !
  ! Basic check
  !
  if(nspec.ne.dust_nr_species) stop 440
  if(mrw_db_ntemp.le.0) stop 441
  !
  ! Find temperature in the table
  !
  call hunt(mrw_db_temp,mrw_db_ntemp,temp,itemp)
  if(itemp.lt.1) then
     itemp = 1
     eps   = 0.d0
  endif
  if(itemp.ge.mrw_db_ntemp) then
     write(stdo,*) 'ERROR in Planck Database for Modified Random Walk: Temperature out of range:'
     write(stdo,*) '  Temperature = ',temp
     write(stdo,*) '  T-max       = ',mrw_db_temp(mrw_db_ntemp)
     stop
  endif
  eps = ( temp - mrw_db_temp(itemp) ) / ( mrw_db_temp(itemp+1) - mrw_db_temp(itemp) )
  if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 274
  eps1 = 1.d0-eps
  !
  ! Return the planck mean opacities
  !
  do ispec=1,nspec
     kappa_a_pl(ispec) = eps1*mrw_db_kappa_abs_planck(itemp,ispec)   +  &
                          eps*mrw_db_kappa_abs_planck(itemp+1,ispec)
  enddo
  !
end subroutine mrw_get_planckopac_dbase

!-------------------------------------------------------------------------
!             FIND THE CROSSING WITH THE CELL WALL
!
! This is the lite-version of amrray_find_next_location_cart(). Here it
! merely determines the ds until the next crossing. It will not search for
! the next cell etc. This subroutine is only used for inside the cell. It
! is therefore much more rudimentary than the amrray subroutine. 
!-------------------------------------------------------------------------
subroutine find_next_location_cart_lite(x0,x1,pos,dir,ds,idir,ilr)
  implicit none
  double precision :: x0(1:3),x1(1:3),pos(1:3),dir(1:3),ds
  integer :: idir,ilr
  double precision :: dsx,dsy,dsz
  !
  ! Find the crossings with all possible cell walls
  !
  if(dir(1).gt.0.d0) then
     dsx = (x1(1)-pos(1))/dir(1)
  elseif(dir(1).lt.0.d0) then
     dsx = (x0(1)-pos(1))/dir(1)
  else
     dsx = 1d99
  endif
  if(dsx.lt.0.d0) stop 707
  if(dir(2).gt.0.d0) then
     dsy = (x1(2)-pos(2))/dir(2)
  elseif(dir(2).lt.0.d0) then
     dsy = (x0(2)-pos(2))/dir(2)
  else
     dsy = 1d99
  endif
  if(dsy.lt.0.d0) stop 708
  if(dir(3).gt.0.d0) then
     dsz = (x1(3)-pos(3))/dir(3)
  elseif(dir(3).lt.0.d0) then
     dsz = (x0(3)-pos(3))/dir(3)
  else
     dsz = 1d99
  endif
  if(dsz.lt.0.d0) stop 10
  !
  ! The smallest one is the actual crossing
  !
  ds = min(dsx,dsy,dsz)
  !
  ! Check which cell wall it escaped through
  !
  if(dsx.le.min(dsy,dsz)) then
     idir = 1
  elseif(dsy.le.dsz) then
     idir = 2
  else
     idir = 3
  endif
  if(dir(idir).gt.0.d0) then
     ilr = 2
  else
     ilr = 1
  endif
end subroutine find_next_location_cart_lite

!-------------------------------------------------------------------------
!             FIND THE CROSSING WITH THE CELL WALL
!
! This is the lite-version of amrray_find_next_location_spher(). Here it
! merely determines the ds until the next crossing. It will not search for
! the next cell etc. This subroutine is only used for inside the cell. It
! is therefore much more rudimentary than the amrray subroutine. 
!-------------------------------------------------------------------------
subroutine find_next_location_spher_lite(x0,x1,pos,dir,cross_ds,idir,ilr,idim,&
                                         sincost12,sincosp12)
  implicit none
  double precision :: x0(1:3),x1(1:3),pos(1:3),dir(1:3)
  double precision, optional :: sincost12(1:2,1:2),sincosp12(1:2,1:2)
  integer :: idir,ilr,idim
  double precision :: sct12(1:2,1:2),scp12(1:2,1:2)
  double precision :: r0dotdir,cross_ds,det,ds_try,sgnz,st2,ct2,r02
  double precision :: pa,pb,pc,eps,sdet,dum1,dum2
  integer :: icross
  logical :: topquadrant,crossequator
  doubleprecision :: eps_thres
  parameter(eps_thres=1d-4)
  doubleprecision :: small
  parameter(small=1d-12)
  !
  ! Some calculations in advance
  !
  r02       = pos(1)**2 + pos(2)**2 + pos(3)**2
  r0dotdir  = pos(1) * dir(1) + pos(2) * dir(2) + pos(3) * dir(3)
  icross    = 0
  cross_ds  = 1d99
  !
  ! Try out a crossing with inner radius cell wall.
  ! Only if ray is moving inward. It is always the smaller of the two.
  !
  if(r0dotdir.lt.0.d0) then
     det       = r0dotdir*r0dotdir + x0(1)**2 - r02
     if(det.gt.0.d0) then
        ds_try = -r0dotdir - sqrt(det)
        if(ds_try.lt.0.d0) ds_try = 0.d0   ! Could happen due to round-off errors
        if(ds_try.lt.cross_ds) then
           icross = 1
           cross_ds = ds_try
        endif
     endif
  endif
  !
  ! Try out crossing with outer radius.
  ! Since we are within the cell the solution must be
  ! the largest of the two.
  !
  det       = r0dotdir*r0dotdir + x1(1)**2 - r02
  ds_try    = -r0dotdir + sqrt(det)
  if(ds_try.lt.0.d0) ds_try = 0.d0    ! Could happen due to round-off errors
  if(ds_try.lt.cross_ds) then
     icross = 2
     cross_ds = ds_try
  endif
  ! 
  ! Check if the theta dimension is active (if not, then we are doing
  ! 1-D spherical radiative transfer)
  !
  if(idim.ge.2) then
     !
     ! First check if we are in the top or bottom quadrant
     !
     if(x0(2).lt.pihalf) then
        topquadrant = .true.
        sgnz        = 1.d0
        if(pos(3).lt.0.d0) stop 3331
     else
        topquadrant = .false.
        sgnz        = -1.d0
        if(pos(3).gt.0.d0) then
           write(stdo,*) 'ERROR in spherical coordinates in MRW'
           write(stdo,*) pos(3),x0(2)
           stop 3332
        endif
     endif
     !
     ! Since sines and cosines are expensive, it is best to precalculate
     ! them and give them here. 
     !
     ! BUGFIX 31.03.15: For the lower quadrant we have to swap the
     !                  theta0 and theta1, for algorithmic reasons.
     !
     if(present(sincost12)) then
        if(topquadrant) then
           sct12(1,1) = sincost12(1,1)
           sct12(2,1) = sincost12(2,1)
           sct12(1,2) = sincost12(1,2)
           sct12(2,2) = sincost12(2,2)
        else
           sct12(1,1) = sincost12(1,2)
           sct12(2,1) = sincost12(2,2)
           sct12(1,2) = sincost12(1,1)
           sct12(2,2) = sincost12(2,1)
        endif
     else
        if(topquadrant) then
           sct12(1,1) = sin(x0(2))
           sct12(2,1) = cos(x0(2))
           sct12(1,2) = sin(x1(2))
           sct12(2,2) = cos(x1(2))
        else
           sct12(1,1) = sin(x1(2))
           sct12(2,1) = cos(x1(2))
           sct12(1,2) = sin(x0(2))
           sct12(2,2) = cos(x0(2))
        endif
     endif
     if(idim.ge.3) then
        if(present(sincosp12)) then
           scp12(:,:) = sincosp12(:,:)
        else
           scp12(1,1) = sin(x0(3))
           scp12(2,1) = cos(x0(3))
           scp12(1,2) = sin(x1(3))
           scp12(2,2) = cos(x1(3))
        endif
     endif
     !
     ! Check if the cell touches the equatorial plane
     !
     !if(abs(x1(2)-pihalf).lt.1d-10) then    ! BUG FOUND AND FIXED BY ATTILA JUHASZ
     if((abs(x1(2)-pihalf).lt.1d-10).or.(abs(x0(2)-pihalf).lt.1d-10)) then
        !
        ! This is simple: just the crossing with the z=0-plane
        !
        if(dir(3).ne.0.d0) then
           ds_try = -pos(3) / dir(3)
           if((ds_try.gt.0.d0).and.(ds_try.lt.cross_ds)) then
              !
              ! Yes, we have a valid crossing
              !
              if(pos(3).lt.0.d0) then
                 !
                 ! We cross from below
                 !
                 icross   = 3
                 cross_ds = ds_try
                 crossequator = .true.
              else
                 !
                 ! We cross from above
                 !
                 icross   = 4
                 cross_ds = ds_try
                 crossequator = .true.
              endif
           endif
        endif
     endif
     !
     ! Now check for crossing with theta=constant cone: the one farthest
     ! away from the midplane. Since we are definitely inside this cell
     ! this means that we are currently outside this cone.
     !
     st2 = sct12(1,1)**2
     ct2 = sct12(2,1)**2
     !
     ! Be careful if st2=0.d0: this means that this theta=const wall of
     ! the cell is in fact the z-axis, which is infinitely thin and 
     ! should therefore never yield a hit. Since I do not want to be
     ! compiler dependent, I do not want this to depend on the precise
     ! rounding-off way of formulae, so I do a real check here. Note 
     ! that I only have to check this for the theta1 wall, because that
     ! is by definition the one closest to the z-axis.
     !
     if(st2.gt.0.d0) then
        !
        ! Compute some of the coefficients
        !
        pa    = ct2*(dir(1)**2+dir(2)**2)-st2*dir(3)**2
        pb    = 2.d0*ct2*(dir(1)*pos(1) + dir(2)*pos(2)) - &
                2.d0*st2*dir(3)*pos(3)
        pc    = ct2*(pos(1)**2+pos(2)**2)-st2*pos(3)**2
        !
        ! Compute eps == 4*pa*pc/pb^2
        !
        eps   = 4.d0*pa*pc/(pb*pb+1d-99)
        ! 
        ! Now check out if there is a solution
        !
        det   = pb**2 - 4.d0*pa*pc
        if(det.gt.0.d0) then
           !
           ! If |eps| < eps_thres and pb<0 use, instead of the normal
           ! formula, a Taylor series to remove near-cancellation of terms.
           !
           if((abs(eps).lt.eps_thres).and.(pb.lt.0.d0)) then
              !
              ! Taylor series: 
              !
              !  1-sqrt(1-x) = (1/2)*x + (1/8)*x^2 + (1/16)*x^3
              !                + (5/128)*x^4 + ...
              !
              ds_try = ( pc / abs(pb) ) * ( 1.d0 + 0.25d0*eps + &
                         0.125d0*eps**2 + (5.d0/64.d0)*eps**3 )
              if((ds_try.gt.0.d0).and.(ds_try.lt.cross_ds)) then
                 if(topquadrant) then
                    icross   = 3
                 else
                    icross   = 4
                 endif
                 cross_ds = ds_try
              endif
              !
           else               
              !
              ! Normal formula:
              !
              ! Two scenarios in one:
              !
              ! pa>0:
              ! Ray will cut through either the top or the bottom cone
              ! twice. Since we are outside the cone, the smallest 
              ! positive solution is the one: ds=-0.5*(pb+sdet)/pa.
              !
              ! pa<0:
              ! Ray will cut through both the top and the bottom cone.
              ! There will be a positive and a negative ds. We need
              ! the positive one. Happens also to be -0.5*(pb+sdet)/pa.
              !
              ! NOTE: In principle we should also check if we cross the
              !       right one of the two cones (top/bottom). But this
              !       is done automatically, because if it is the wrong
              !       cone, the ds_try will be overruled (now or later)
              !       because there will be a smaller one.
              !
              if(pa.ne.0.d0) then
                 sdet    = sqrt(det)
                 ds_try  = - 0.5d0 * ( pb + sdet ) / pa
                 if((ds_try.gt.0.d0).and.(ds_try.lt.cross_ds)) then
                    if(topquadrant) then
                       icross   = 3
                    else
                       icross   = 4
                    endif
                    cross_ds = ds_try
                 endif
              endif
           endif
        endif
     endif
     !
     ! Now check for crossing with theta=constant cone: the one closest
     ! to the midplane. Since we are definitely inside this cell
     ! this means that we are currently inside this cone.
     !
     st2 = sct12(1,2)**2
     ct2 = sct12(2,2)**2
     !
     ! But do this check only if it has not yet been done before for
     ! the equator crossing
     !
     if(ct2.ne.0.d0) then
        !
        ! Compute some of the coefficients
        !
        pa    = ct2*(dir(1)**2+dir(2)**2)-st2*dir(3)**2
        pb    = 2.d0*ct2*(dir(1)*pos(1) + dir(2)*pos(2)) - &
                2.d0*st2*dir(3)*pos(3)
        pc    = ct2*(pos(1)**2+pos(2)**2)-st2*pos(3)**2
        !
        ! Compute eps == 4*pa*pc/pb^2
        !
        eps   = 4.d0*pa*pc/(pb*pb+1d-99)
        ! 
        ! Now check out if there is a solution
        !
        det   = pb**2 - 4.d0*pa*pc
        if(det.gt.0.d0) then
           !
           ! If |eps| < eps_thres and pb>0 use, instead of the normal
           ! formula, a Taylor series to remove near-cancellation of terms.
           !
           if((abs(eps).lt.eps_thres).and.(pb.gt.0.d0)) then
              !
              ! Taylor series: 
              !
              !  1-sqrt(1-x) = (1/2)*x + (1/8)*x^2 + (1/16)*x^3
              !                + (5/128)*x^4 + ...
              !
              ds_try = - ( pc / pb ) * ( 1.d0 + 0.25d0*eps +   &
                         0.125d0*eps**2 + (5.d0/64.d0)*eps**3 )
              if((ds_try.gt.0.d0).and.(ds_try.lt.cross_ds)) then
                 if(topquadrant) then
                    icross   = 4
                 else
                    icross   = 3
                 endif
                 cross_ds = ds_try
              endif
              !
           else               
              !
              ! Normal formula:
              !
              ! Two scenarios in one:
              !
              ! pa>0:
              ! Ray will cut through either the correct cone twice.
              ! Since we are inside the cone, the largest of the two 
              ! de values is the one: ds=-0.5*(pb-sdet)/pa.
              !
              ! pa<0:
              ! Ray will cut through both the top and the bottom cone.
              ! There will be either two positive or two negative ds
              ! values. If there are two negative ones, then the ray
              ! is moving away from the cone to infinity. If the two
              ! are positive, then we should take the smallest of the
              ! two. Happens also to be -0.5*(pb-sdet)/pa.
              !
              if(pa.ne.0.d0) then
                 sdet    = sqrt(det)
                 ds_try  = - 0.5d0 * ( pb - sdet ) / pa
                 if((ds_try.gt.0.d0).and.(ds_try.lt.cross_ds)) then
                    if(topquadrant) then
                       icross   = 4
                    else
                       icross   = 3
                    endif
                    cross_ds = ds_try
                 endif
              endif
           endif
        endif
     endif
  endif
  !
  ! Check if phi dimension is active. If not, then we are doing either
  ! 1-D spherical or 2-D axisymmetric radiative transfer
  !
  if(idim.ge.3) then
     !
     ! Find the distance to the crossing with the smallest phi
     !
     dum1 = pos(1)*scp12(1,1) - pos(2)*scp12(2,1)
     dum2 = dir(2)*scp12(2,1) - dir(1)*scp12(1,1)
     if(dum2.lt.-small) then   ! Thereby selecting also only outgoing rays
        ds_try = dum1 / dum2
     else
        ds_try = -1.d0
     endif
     if((ds_try.ge.0.d0).and.(ds_try.lt.cross_ds)) then
        icross   = 5
        cross_ds = ds_try
     endif
     !
     ! Find the distance to the crossing with the largest phi
     !
     dum1 = pos(1)*scp12(1,2) - pos(2)*scp12(2,2)
     dum2 = dir(2)*scp12(2,2) - dir(1)*scp12(1,2)
     if(dum2.gt.small) then    ! Thereby selecting also only outgoing rays
        ds_try = dum1 / dum2
     else
        ds_try = -1.d0
     endif
     if((ds_try.ge.0.d0).and.(ds_try.lt.cross_ds)) then
        icross   = 6
        cross_ds = ds_try
     endif
  endif
  !
  ! Find the cell to which we move.
  !
  select case(icross) 
  case(1)
     idir = 1
     ilr  = 2
  case(2)
     idir = 1
     ilr  = 1
  case(3)
     idir = 2
     ilr  = 2
  case(4)
     idir = 2
     ilr  = 1
  case(5)
     idir = 3
     ilr  = 2
  case(6)
     idir = 3
     ilr  = 1
  end select
end subroutine find_next_location_spher_lite

!-------------------------------------------------------------------------
!               FIND CLOSEST EDGE IN CARTESIAN COORDINATES
!-------------------------------------------------------------------------
subroutine find_closest_wall_cart_lite(x0,x1,pos,dist,idir,ilr)
  implicit none
  double precision :: x0(1:3),x1(1:3),pos(1:3),dist
  integer :: idir,ilr
  double precision :: dxl(1:3),dxr(1:3)
  dist = 1d99
  idir = 0
  ilr  = 0
  dxl(:) = pos(:) - x0(:)
  dxr(:) = x1(:) - pos(:)
  if(dxl(1).lt.dist) then
     dist = dxl(1)
     idir = 1
     ilr  = 1
  endif
  if(dxl(2).lt.dist) then
     dist = dxl(2)
     idir = 2
     ilr  = 1
  endif
  if(dxl(3).lt.dist) then
     dist = dxl(3)
     idir = 3
     ilr  = 1
  endif
  if(dxr(1).lt.dist) then
     dist = dxr(1)
     idir = 1
     ilr  = 2
  endif
  if(dxr(2).lt.dist) then
     dist = dxr(2)
     idir = 2
     ilr  = 2
  endif
  if(dxr(3).lt.dist) then
     dist = dxr(3)
     idir = 3
     ilr  = 2
  endif
  ! CHANGE 16.09.2016: No longer stopping the code if this happens
  if(dist.lt.0.d0) then
     write(stdo,*) 'Warning in Modified Random Walk: Round-off error in calculating closest wall distance.'
     dist=0.d0
  endif
end subroutine find_closest_wall_cart_lite

!-------------------------------------------------------------------------
!                FIND CLOSEST EDGE IN SPHERICAL COORDINATES
! idim = 1, 2 or 3 depending on whether you have (r) = 1D, (r,theta) = 2D,
! or (r,theta,phi) = 3D.
!-------------------------------------------------------------------------
subroutine find_closest_wall_spher_lite(x0,x1,pos,dist,idir,ilr,idim, &
                                        sincost12,sincosp12)
  implicit none
  double precision :: x0(1:3),x1(1:3),pos(1:3),dist
  integer :: idir,ilr,idim
  double precision, optional :: sincost12(1:2,1:2),sincosp12(1:2,1:2)
  double precision :: sct12(1:2,1:2),scp12(1:2,1:2)
  double precision :: dxl(1:3),dxr(1:3),r0,rc
  !
  ! Do a few checks
  !
  if(x1(1).lt.x0(1)) stop 9475
  if(x1(2).lt.x0(2)) stop 9476
  if(x1(3).lt.x0(3)) stop 9477
  !
  ! Since sines and cosines are expensive, it is best to precalculate
  ! them and give them here. 
  !
  if(idim.ge.2) then
     if(present(sincost12)) then
        sct12(:,:) = sincost12(:,:)
     else
        sct12(1,1) = sin(x0(2))
        sct12(2,1) = cos(x0(2))
        sct12(1,2) = sin(x1(2))
        sct12(2,2) = cos(x1(2))
     endif
     if(idim.ge.3) then
        if(present(sincosp12)) then
           scp12(:,:) = sincosp12(:,:)
        else
           scp12(1,1) = sin(x0(3))
           scp12(2,1) = cos(x0(3))
           scp12(1,2) = sin(x1(3))
           scp12(2,2) = cos(x1(3))
        endif
     endif
  endif
  !
  ! Convert from cartesian to spherical
  !
  r0     = sqrt(pos(1)**2+pos(2)**2+pos(3)**2)
  !
  ! Now compute the distances
  !
  dist = 1d99
  idir = 0
  ilr  = 0
  dxl(1) = r0 - x0(1)
  dxr(1) = x1(1) - r0
  if(dxl(1).lt.dist) then
     dist = dxl(1)
     idir = 1
     ilr  = 1
  endif
  if(dxr(1).lt.dist) then
     dist = dxr(1)
     idir = 1
     ilr  = 2
  endif
  if(idim.ge.2) then
     rc     = sqrt(pos(1)**2+pos(2)**2)
     dxl(2) = sct12(2,1)*rc    -sct12(1,1)*pos(3)
     dxr(2) = sct12(1,2)*pos(3)-sct12(2,2)*rc
     ! BUGFIX 09.04.2015: Test was too tight: added epsmargin
     !if(dxl(2).lt.0.d0) stop 3901
     !if(dxr(2).lt.0.d0) stop 3902
     if(dxl(2).lt.dist) then
        dist = dxl(2)
        idir = 2
        ilr  = 1
     endif
     if(dxr(2).lt.dist) then
        dist = dxr(2)
        idir = 2
        ilr  = 2
     endif
     if(idim.ge.3) then
        if((x0(3).gt.1d-10).or.(twopi-x1(3).gt.1d-10)) then  ! Check cyclic phi
           dxl(3) = scp12(2,1)*pos(2)-scp12(1,1)*pos(1)
           dxr(3) = scp12(1,2)*pos(1)-scp12(2,2)*pos(2)
           ! BUGFIX 09.04.2015: Test was too tight: added epsmargin
           !if(dxl(3).lt.0.d0) stop 3911
           !if(dxr(3).lt.0.d0) stop 3912
           if(dxl(3).lt.dist) then
              dist = dxl(3)
              idir = 3
              ilr  = 1
           endif
           if(dxr(3).lt.dist) then
              dist = dxr(3)
              idir = 3
              ilr  = 2
           endif
        endif
     endif
  endif
  ! BUGFIX AND CHANGE 16.09.2016: No longer stopping the code if this happens, and in fact resetting dist always >=0
  if(dist.lt.0.d0) then
     dist = 0.d0
     if(dist.lt.-epsmargin*x0(1)) then
        ! BUGFIX 09.04.2015: Test was too tight: added epsmargin before making a warning
        write(stdo,*) 'Warning in Modified Random Walk: Round-off error in calculating closest wall distance. Continuing...'
     endif
  endif
end subroutine find_closest_wall_spher_lite

!-------------------------------------------------------------------------
!          FIND THE SMALLEST DIMENSION OF THE CELL - CARTESIAN
!-------------------------------------------------------------------------
subroutine find_smallest_size_cart_lite(x0,x1,size)
  implicit none
  double precision :: x0(1:3),x1(1:3),size
  size = 1d99
  if(x1(1)-x0(1).lt.size) size=x1(1)-x0(1)
  if(x1(2)-x0(2).lt.size) size=x1(2)-x0(2)
  if(x1(3)-x0(3).lt.size) size=x1(3)-x0(3)
  if(size.le.0.d0) stop 5931
end subroutine find_smallest_size_cart_lite

!-------------------------------------------------------------------------
!          FIND THE SMALLEST DIMENSION OF THE CELL - SPHERICAL
!-------------------------------------------------------------------------
subroutine find_smallest_size_spher_lite(x0,x1,size,idim)
  implicit none
  double precision :: x0(1:3),x1(1:3),size
  integer :: idim
  double precision :: dist(1:3)
  size = 1d99
  dist(1)=x1(1)-x0(1)
  if(dist(1).lt.size) size=dist(1)
  if(idim.ge.2) then
     dist(2)=x0(1)*(x1(2)-x0(2))
     if(dist(2).lt.size) size=dist(2)
     if(idim.ge.3) then
        dist(3)=x0(1)*sin(x0(2))*(x1(3)-x0(3))
        if(dist(3).lt.size) size=dist(3)
     endif
  endif
  if(size.le.0.d0) stop 5931
end subroutine find_smallest_size_spher_lite

!-------------------------------------------------------------------------
!                    CONVERT X,Y,Z INTO R,THETA,PHI
!-------------------------------------------------------------------------
subroutine convert_xyz_into_rthetaphi(x,y,z,r,theta,phi)
  implicit none
  double precision :: x,y,z,r,theta,phi
  r     = sqrt(x**2+y**2+z**2)
  theta = acos(z/(r+1d-199))
  if(x.eq.0.d0) then
     if(y.ge.0.d0) then
        phi = pihalf
     else
        phi = 3.d0*pihalf 
     endif
  else
     phi = atan(y/x)
     if(x.gt.0.d0) then
        if(y.lt.0.d0) then
           phi = phi + twopi
        endif
     else
        phi = phi + pi
     endif
  endif
  if((phi.lt.0.d0).or.(phi.gt.twopi)) stop 1309
end subroutine convert_xyz_into_rthetaphi

!-------------------------------------------------------------------------
!  MAKE A STEP ACCORDING TO THE MODIFIED RANDOM WALK OF MIN ET AL. 2009
!
! Note: This is a simplified version of the Min et al. 2009. Here we
!       immediately use the average path length of the random walk
!       in the cell. Also the scattering and the average absorption
!       opacity are not calculated as exactly as in Min et al. 2009.
!-------------------------------------------------------------------------
subroutine mrw_step(r0,pos,energy,enerphot,alpha_abs,lfreeross)
  implicit none
  double precision :: r0,pos(1:3),alpha_abs,lfreeross,energy,enerphot
  double precision :: lpath_av,twicelfr,dir(1:3),dummy
  logical :: outside
  !
  ! Just some calculation:
  ! 
  twicelfr = 2.d0*lfreeross
  !
  ! First check if the MRW is justified (if this is broken, then 
  ! we would have superluminal motion!)
  !
  if(r0.lt.twicelfr) then
     write(stdo,*) 'ERROR: Cannot perform MRW when R_0 < 2 l_freeross ',&
          '(would result in superluminal motion)!'
     stop
  endif
  !
  ! Compute the average length of the diffusive path the photon takes before
  ! hitting the edge of the sphere with radius R_0. Note that in the Min 
  ! et al. paper this was chosen randomly from a distribtion. But since
  ! the average is certainly equally right and produces less noise and is
  ! also cheaper, we'll use the average here.
  !
  ! Here is how this works: In Min et al. 2009 we wrote an expression for the
  ! chance P(t) that a photon package has still not yet crossed the radius R_0
  ! after a time t (Eqs. 7,8). The derivative of that -dP(t)/dt is the probability 
  ! distribution function for the time t when the photon package first reaches
  ! R_0. If we integrate this with t (i.e. int_0^infty -(dP(t)/dt)*t*dt) we
  ! get the expectation value <t>, i.e. the average time it takes for a
  ! photon to first reach R_0. This is:
  !
  !   <t> = int_0^infty -(dP(t)/dt)*t*dt = R_0^2/(6Dc)
  !
  ! with D the MRW diffusion coefficient. Here we use Eq. 12 of that paper
  ! (D=<d^2>/(6<d>)) to approximate this diffusion coefficient, where
  !
  !   <d> = 1/(rho*kappa_R)        and      <d^2> = 2/(rho*kappa_R)^2
  !
  ! which leads to
  !
  !   <t> = (R_0^2/2c) rho kappa_R 
  !
  ! The distance the photon moves along its path is then <l>=c<t>
  !
  !   <l> = (R_0^2/2) rho kappa_R
  !
  ! or, if we define l_ross = 1/(rho*kappa_R) we get
  !
  !   <l> = R_0^2/(2*l_ross)
  !
  ! So we use this expression as the typical path distance moved by the photon.
  !
  lpath_av = r0**2 / twicelfr
  !
  ! This allows us to update the energy in the cell
  !
  energy = energy + lpath_av * alpha_abs * enerphot
  !
  ! Find a random direction
  !
  outside = .true.
  do while(outside) 
     dir(1) = 2*ran2(iseed)-1.d0
     dir(2) = 2*ran2(iseed)-1.d0
     dir(3) = 2*ran2(iseed)-1.d0
     dummy  = dir(1)**2 + dir(2)**2 + dir(3)**2 
     if(dummy.le.1.d0) outside=.false.
  enddo
  dummy  = sqrt(dir(1)**2 + dir(2)**2 + dir(3)**2)
  dir(1) = dir(1) / dummy
  dir(2) = dir(2) / dummy
  dir(3) = dir(3) / dummy
  !
  ! Now make the spatial step
  !
  pos(:) = pos(:) + dir(:)*r0
  !
end subroutine mrw_step


!-------------------------------------------------------------------------
!                        CHECK IF PHOTON IN CELL
!
! If the photon is slightly outside of the cell, the function returns
! false. But if the photon is badly outside of the cell then it stops
! because then presumably a bug is occurring.
!-------------------------------------------------------------------------
function mrw_check_if_in_cell(cellx0,cellx1,pos,icoord,idim)
  implicit none
  double precision :: cellx0(1:3),cellx1(1:3),pos(1:3)
  double precision :: r,theta,phi,size
  integer :: icoord,idim
  logical :: mrw_check_if_in_cell,baderror
  mrw_check_if_in_cell = .true.
  baderror = .false.
  if(icoord.le.99) then
     if((pos(1).lt.cellx0(1)).or.(pos(1).gt.cellx1(1))) then
        mrw_check_if_in_cell = .false.
        size = cellx1(1) - cellx0(1)
        if(pos(1).lt.cellx0(1)) then
           if(cellx0(1)-pos(1).gt.1d-1*size) baderror=.true.
        else
           if(pos(1)-cellx1(1).gt.1d-1*size) baderror=.true.
        endif
        if(baderror) then
           write(stdo,*) 'ERROR in MRW: position out of x range'
           write(stdo,*) '  x = ',pos(1),', range = ',cellx0(1),cellx1(1)
           stop
        endif
     endif
     if((pos(2).lt.cellx0(2)).or.(pos(2).gt.cellx1(2))) then
        mrw_check_if_in_cell = .false.
        size = cellx1(2) - cellx0(2)
        if(pos(2).lt.cellx0(2)) then
           if(cellx0(2)-pos(2).gt.1d-1*size) baderror=.true.
        else
           if(pos(2)-cellx1(2).gt.1d-1*size) baderror=.true.
        endif
        if(baderror) then
           write(stdo,*) 'ERROR in MRW: position out of y range'
           write(stdo,*) '  y = ',pos(2),', range = ',cellx0(2),cellx1(2)
           stop
        endif
     endif
     if((pos(3).lt.cellx0(3)).or.(pos(3).gt.cellx1(3))) then
        mrw_check_if_in_cell = .false.
        size = cellx1(3) - cellx0(3)
        if(pos(3).lt.cellx0(3)) then
           if(cellx0(3)-pos(3).gt.1d-1*size) baderror=.true.
        else
           if(pos(3)-cellx1(3).gt.1d-1*size) baderror=.true.
        endif
        if(baderror) then
           write(stdo,*) 'ERROR in MRW: position out of z range'
           write(stdo,*) '  z = ',pos(3),', range = ',cellx0(3),cellx1(3)
           stop
        endif
     endif
  else
     call convert_xyz_into_rthetaphi(pos(1),pos(2),pos(3), &
          r,theta,phi)
     ! BUGFIX 09.04.2015: Test was too tight: added epsmargin
     if((r.lt.cellx0(1)*(1.d0-epsmargin)).or. &
          (r.gt.cellx1(1)*(1.d0+epsmargin))) then
        mrw_check_if_in_cell = .false.
        size = cellx1(1) - cellx0(1)
        if(r.lt.cellx0(1)*(1.d0-epsmargin)) then
           if(cellx0(1)-r.gt.1d-1*size) baderror=.true.
        else
           if(r-cellx1(1).gt.1d-1*size) baderror=.true.
        endif
        if(baderror) then
           write(stdo,*) 'ERROR in MRW: position out of r range'
           write(stdo,*) '  r = ',r,', range = ',cellx0(1),cellx1(1)
           stop
        endif
     endif
     if(idim.ge.2) then
        if((theta.lt.cellx0(2)-epsmargin).or. &
             (theta.gt.cellx1(2)+epsmargin)) then
           mrw_check_if_in_cell = .false.
           size = cellx1(2) - cellx0(2)    ! Size in theta
           if(theta.lt.cellx0(2)-epsmargin) then
              if(cellx0(2)-theta.gt.1d-1*size) baderror=.true.
           else
              if(theta-cellx1(2).gt.1d-1*size) baderror=.true.
           endif
           if(baderror) then
              write(stdo,*) 'ERROR in MRW: position out of theta range'
              write(stdo,*) '  theta = ',theta,', range = ',cellx0(2),cellx1(2)
              stop
           endif
        endif
        if(idim.ge.3) then
           if((phi.lt.cellx0(3)-epsmargin).or. &
                (phi.gt.cellx1(3)+epsmargin)) then
              mrw_check_if_in_cell = .false.
              size = cellx1(3) - cellx0(3)    ! Size in phi
              if(phi.lt.cellx0(3)-epsmargin) then
                 if(cellx0(3)-phi.gt.1d-1*size) baderror=.true.
              else
                 if(phi-cellx1(3).gt.1d-1*size) baderror=.true.
              endif
              if(baderror) then
                 write(stdo,*) 'ERROR in MRW: position out of phi range'
                 write(stdo,*) '  theta = ',phi,', range = ',cellx0(3),cellx1(3)
                 stop
              endif
           endif
        endif
     endif
  endif
  return
end function mrw_check_if_in_cell


!-------------------------------------------------------------------------
!                     MODIFIED RANDOM WALK IN A CELL
!
! Subroutine for performing an absorption+reemission random walk within a
! single cell, using the Modified Random Walk method of Min, Dullemond,
! Dominik, De Koter & Hovenier (2009) A&A 497, 155. The method is based on
! the method of Fleck & Canfield (1984) J. Compt. Phys. 54, 508.
!
! This subroutine only works under the assumption that all the dust in 
! the cell has the same temperature and a common energy. So if you use
! this subroutine in a code that has independent dust temperatures for
! independent dust species, then (a) add all the energies together and
! (b) compute a single alpha_abs and alpha_ross for them all. 
!
! Difference from the original MRW of Min et al.: Here we do not randomly
! select the path length in each MRW step. Instead we use the average
! expected value for this. This saves one random number and a lookup in a
! table, and it saves a bit of programming.
!
! Another difference from the original MRW of Min et al.: This subroutine
! simply uses the Planck mean absorption opacity and the Rosseland mean
! total opacity as was figured out by Robitaille (2010) A&A 520, A70.
!
! ARGUMENTS:
!  cellx0(1:3)       The left coordinates of the cell (cart or spher)
!  cellx1(1:3)       The right coordinates of the cell (cart or spher)
!  pos(1:3)          The current position of the photon (always cart)
!  dir(1:3)          The current direction of the photon (always cart)
!  energy            The current energy of all the dust in the cell
!  enerphot          The energy of the photon package
!  enthres           If energy exceeds enthres, then return (to recompute opacs)
!  alpha_abs         The Planck mean absorption opacity
!  alpha_ross        The Rosseland mean total opacity
!  mrw_gamma         The "gamma" of the Min et al paper. A
!                    good value is 4.0. It MUST be >= 2.0. If set to
!                    1d99 the MRW is not used and a normal random walk
!                    is used.
!  icoord            0..99    means cartesian coordinates
!                    100..199 means spherical coordinates
!  idim              1, 2 or 3 dimensions (used only for spher coord)
!
! OPTIONAL ARGUMENTS:
!  taumargin         If set, then the photon stops short of the boundary
!                    of the cell by an optical depth taumargin. The idea
!                    is that you take this, say, 1E-3, and you can then 
!                    let the remaining be done using the normal 
!                    ray-tracing routines who also find the next cell.
!
! RETURNS:
!  pos(1:3)          The final position when the photon escapes from
!                    the cell.
!  dir(1:3)          The final direction when the photon escapes from
!                    the cell.
!  energy            The updated energy of the cell.
!  idir_cross and    The identifiers which cell wall the photon left from.
!   ilr_cross           If they are 0, then we reached enthress before
!                       we reached the cellwall.
!  nrevents          The number of MRW "scattering" events (just FYI)
!
! Notes:
!  - For the alpha_ross you must also include the scattering. If you
!    have anisotropic scattering, then this anisotropy will reduce the
!    scattering efficiency and therefore its contribution to the 
!    alpha_ross will be less than for isotropic scattering. How to
!    implement this is written in Min et al. (2009).
!
!-------------------------------------------------------------------------
subroutine modified_random_walk(cellx0,cellx1,pos,dir,energy,enerphot,     &
                                enthres,alpha_abs,alpha_ross,mrw_gamma,    &
                                icoord,idim,idir_cross,ilr_cross,nrevents, &
                                taumargin)
  implicit none
  double precision :: cellx0(1:3),cellx1(1:3),pos(1:3),dir(1:3)
  double precision :: alpha_abs,alpha_ross,energy,enerphot,mrw_gamma,enthres
  integer :: icoord,idim
  integer :: idir_cross,ilr_cross
  double precision :: nrevents
  double precision, optional :: taumargin
  !
  double precision :: sincost12(1:2,1:2),sincosp12(1:2,1:2)
  logical :: notescaped,outside,use_mrw
  double precision :: ds,dscross,dsevent,l_freeross
  double precision :: taupath,dummy,ds_margin
  double precision :: distclose,tauclose
  integer :: idirclose,ilrclose
  double precision :: r,theta,phi
  !
  notescaped = .true.
  !
  ! Make sure that the direction is normalized
  !
  dummy = dir(1)**2 + dir(2)**2 + dir(3)**2
  if(dummy.eq.0.d0) stop 5
  dir(:) = dir(:) / sqrt(dummy)
  !
  ! Default
  !
  l_freeross = 1.d0 / alpha_ross
  idir_cross = 0
  ilr_cross  = 0
  if(present(taumargin)) then
     ds_margin = taumargin * l_freeross
  else
     ds_margin = -1.d0
  endif
  !
  ! Check
  !
  if(mrw_gamma.lt.2.d0) then
     write(stdo,*) 'ERROR in Modified Random Walk: The gamma must be >=2.d0 ',&
          'otherwise light would go faster than c.'
     stop
  endif
  !
  ! Pre-compute the sines and cosines
  !
  sincost12(1,1) = sin(cellx0(2))
  sincost12(2,1) = cos(cellx0(2))
  sincost12(1,2) = sin(cellx1(2))
  sincost12(2,2) = cos(cellx1(2))
  sincosp12(1,1) = sin(cellx0(3))
  sincosp12(2,1) = cos(cellx0(3))
  sincosp12(1,2) = sin(cellx1(3))
  sincosp12(2,2) = cos(cellx1(3))
  !
  ! Now do the multiple Monte Carlo steps until escape
  !
  do while(notescaped)
     !
     ! First check if we should do a normal or a modified random walk step
     !
     if(mrw_gamma.ge.1d20) then
        use_mrw = .false.
     else
        !
        ! Find the closest distance and optical depth
        !
        if(icoord.le.99) then
           !
           ! Cartesian coordinates
           !
           call find_closest_wall_cart_lite(cellx0,cellx1,pos, &
                distclose,idirclose,ilrclose)
        else
           !
           ! Spherical coordinates
           !
           call find_closest_wall_spher_lite(cellx0,cellx1,pos, &
                distclose,idirclose,ilrclose,idim,sincost12,sincosp12)
        endif
        !
        ! Make a slight adjustment for stability
        ! IMPROVEMENT 31.03.15
        !
        distclose = distclose * 0.99
        !
        ! Compute the corresponding optical depth
        !
        tauclose = distclose * alpha_ross
        !
        ! Check if this is large enough for the MRW
        !
        if(tauclose.ge.mrw_gamma) then
           use_mrw = .true.
        else
           use_mrw = .false.
        endif
     endif
     !
     ! Now switch between normal RW and modified RW
     !
     if(use_mrw) then
        !
        ! Use the MRW (Modified Random Walk) method
        !
        call mrw_step(distclose,pos,energy,enerphot,alpha_abs,l_freeross)
        !
     else
        !
        ! No MRW, just the normal Monte Carlo stepping
        !
        ! Which is the next crossing, if no event were to take place?
        !
        if(icoord.le.99) then
           !
           ! Cartesian coordinates
           !
           call find_next_location_cart_lite(cellx0,cellx1,pos,dir,dscross, &
                                          idir_cross,ilr_cross)
        else
           !
           ! Spherical coordinates
           !
           call find_next_location_spher_lite(cellx0,cellx1,pos,dir,dscross, &
                                          idir_cross,ilr_cross,idim,         &
                                          sincost12,sincosp12)
        endif
        !
        ! When is the next event?
        !
        taupath = -log(1.d0-ran2(iseed))
        dsevent = taupath * l_freeross
        !
        ! Go to the next event
        !
        ds = min(dscross,dsevent)
        pos(1) = pos(1) + ds*dir(1)
        pos(2) = pos(2) + ds*dir(2)
        pos(3) = pos(3) + ds*dir(3)
        !
        ! Add energy to the cell
        !
        energy = energy + ds*alpha_abs*enerphot
        !
        ! Check which is first: event or crossing
        !
        if(dsevent.ge.dscross) then
           !
           ! Signal that the photon escapes this cell
           !
           notescaped = .false.
           !
           ! For cartesian coordinates, make sure that the position is
           ! exactly on the edge of the cell. For spherical coordinates
           ! this will anyway not be perfectly possible, so let's hope
           ! it's accurate enough.
           !
           if(icoord.le.99) then
              if(ilr_cross.eq.1) then
                 pos(idir_cross) = cellx0(idir_cross)
              else
                 pos(idir_cross) = cellx1(idir_cross)
              endif
           endif
           !
           ! But if taumargin is set, then backtrack a bit, so that
           ! we stop short of the edge. This can be useful when using
           ! this subroutine inside a Monte Carlo code with its own
           ! wall-crossing and next-cell-finding algorithms.
           !
           if(ds_margin.gt.0.d0) then
              ds_margin = min(ds_margin,dscross)
              !
              ! In the event that `ds_margin` or `dscross` is so small that
              ! it is comparable to the machine precision of `pos`, instead
              ! use `ds_margin = small * pos` to avoid round-off errors 
              ! down the road (e.g., when trying to find a cell wall crossing
              ! in `amrray_find_next_location_spher`).
              ! Bugfix: Jon Ramsey, 14.02.2016
              !
              ds_margin = max(ds_margin, maxval(abs(pos(:)))*small)
              pos(1) = pos(1) - ds_margin*dir(1)
              pos(2) = pos(2) - ds_margin*dir(2)
              pos(3) = pos(3) - ds_margin*dir(3)
           endif
           !
        else
           !
           ! Increase the event counter
           !
           nrevents = nrevents + 1
           !
           ! New random direction
           !
           outside = .true.
           do while(outside) 
              dir(1) = 2*ran2(iseed)-1.d0
              dir(2) = 2*ran2(iseed)-1.d0
              dir(3) = 2*ran2(iseed)-1.d0
              dummy  = dir(1)**2 + dir(2)**2 + dir(3)**2 
              if(dummy.le.1.d0) outside=.false.
           enddo
           dummy  = sqrt(dir(1)**2 + dir(2)**2 + dir(3)**2)
           dir(1) = dir(1) / dummy
           dir(2) = dir(2) / dummy
           dir(3) = dir(3) / dummy
        endif
     endif
     !
     ! Check if energy exceeded enthres. If so, then return, because
     ! we have to recalculate the opacities.
     !
     if(energy.ge.enthres) then
        notescaped = .false.     ! Not really escaped, but stop anyway
     endif
  enddo
end subroutine modified_random_walk

end module montecarlo_module
