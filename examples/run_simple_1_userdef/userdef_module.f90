module userdef_module
use amr_module
use amrray_module
use rtglobal_module
use constants_module
!use camera_module
use dust_module
use lines_module
use stars_module
use quantum_module
use montecarlo_module
use namelist_module

!------------------------------------------------------------------------
! Here you can define your own variables, arrays etc. 
!------------------------------------------------------------------------

integer :: userdef_nx,userdef_ny,userdef_nz
doubleprecision :: userdef_sizex,userdef_sizey,userdef_sizez
doubleprecision :: userdef_radius,userdef_rho0

contains

!------------------------------------------------------------------------
! This subroutine allows you to specify defaults values for your 
! variables.
!
! WARNING: In this subroutine you are not allowed to use write(stdo,*),
!          because the stdo is not yet set. Reason: The defaults are 
!          set before the command-line options are interpreted (so that
!          the defaults can be overwritten by command-line options),
!          but the stdo depends on whether the user calls RADMC-3D as
!          a child or not, which is given by the command-line options.
!------------------------------------------------------------------------
subroutine userdef_defaults()
  implicit none
  !
  userdef_nx = 10
  userdef_ny = 10
  userdef_nz = 10
  userdef_sizex = 100 * au
  userdef_sizey = 100 * au
  userdef_sizez = 100 * au
  userdef_radius = 50 * au
  userdef_rho0 = 1d-16
  !
end subroutine userdef_defaults


!------------------------------------------------------------------------
! Here you can interpret your own command-line options.
!------------------------------------------------------------------------
subroutine userdef_commandline(buffer,numarg,iarg,fromstdi,gotit)
  implicit none
  character*100 :: buffer
  integer :: iarg,numarg
  logical :: gotit,fromstdi
  !
  !
end subroutine userdef_commandline


!------------------------------------------------------------------------
! Here you can do some postprocessing after all command line options
! have been read. No example given here.
!------------------------------------------------------------------------
subroutine userdef_commandline_postprocessing()
  implicit none
end subroutine userdef_commandline_postprocessing


!------------------------------------------------------------------------
! Here you can parse your own keyword entries from the radmc3d.inp file.
!------------------------------------------------------------------------
subroutine userdef_parse_main_namelist()
  implicit none
  !
  ! And example how to include an integer keyword. Use
  ! parse_input_double and parse_input_word for reading doubleprecision
  ! variables or string variables. You can delete the below example.
  ! Note that the keyword name string should always have length 30,
  ! hence the many whitespaces.
  !
  call parse_input_integer('userdef_nx@                   ',userdef_nx)
  call parse_input_integer('userdef_ny@                   ',userdef_ny)
  call parse_input_integer('userdef_nz@                   ',userdef_nz)
  call parse_input_double ('userdef_sizex@                ',userdef_sizex)
  call parse_input_double ('userdef_sizey@                ',userdef_sizey)
  call parse_input_double ('userdef_sizez@                ',userdef_sizez)
  call parse_input_double ('userdef_radius@               ',userdef_radius)
  call parse_input_double ('userdef_rho0@                 ',userdef_rho0)
  !
end subroutine userdef_parse_main_namelist


!------------------------------------------------------------------------
! Here you can do some post-processing after the radmc3d.inp namelist
! reading. 
!------------------------------------------------------------------------
subroutine userdef_main_namelist_postprocessing()
  implicit none
end subroutine userdef_main_namelist_postprocessing


!------------------------------------------------------------------------
! This is the place where you can define your own (base) grid setup,
! read your own frequency grid or set up your own stellar sources.
! No example given here, because it would interfere with basic operations.
!------------------------------------------------------------------------
subroutine userdef_prep_model()
  implicit none
  doubleprecision, allocatable :: xi(:),yi(:),zi(:)
  integer :: i,ntot
  !
  ! Tell the code which grid type is used (AMR)
  !
  igrid_type      = 0             ! Regular grid
  igrid_coord     = 3             ! 3-D cartesian coordinates
  !
  ! Communicate these values to the AMR module
  !
  amr_style       = igrid_type  
  amr_coordsystem = igrid_coord 
  !
  ! Set up the regular grid coordinates in temporary arrays
  !
  allocate(xi(userdef_nx+1),yi(userdef_ny+1),zi(userdef_nz+1))
  do i=1,userdef_nx+1
     xi(i) = 2.d0*userdef_sizex*(i-1.d0)/(userdef_nx*1.d0)-userdef_sizex
  enddo
  do i=1,userdef_ny+1
     yi(i) = 2.d0*userdef_sizey*(i-1.d0)/(userdef_ny*1.d0)-userdef_sizey
  enddo
  do i=1,userdef_nz+1
     zi(i) = 2.d0*userdef_sizez*(i-1.d0)/(userdef_nz*1.d0)-userdef_sizez
  enddo
  !
  ! Compute the total number of cells, so that we can tell the AMR
  ! initializer routine below how much memory it has to reserve. Note
  ! that if you want to refine lateron (while building up the model
  ! for instance) you must make sure here to reserve enough extra
  ! space for all the extra branches and leafs to fit in. For models
  ! with deep refinement levels the nr of branches can be up to (at
  ! most) 15% larger than the nr of leafs (0.125+0.125^2+0.125^3....
  ! = 0.142857).
  !
  ntot = userdef_nx*userdef_ny*userdef_nz
  !
  ! Now set up the grid, for now just a regular grid
  !
  call amr_initialize(.true.,.true.,.true.,            &
       userdef_nx,userdef_ny,userdef_nz,xi,yi,zi,1,1,  &
       ntot,ntot,.false.,.false.,.false.)
  !
  ! Set up the index lists so that we know how to connect each AMR grid
  ! cell to a data memory location.
  !
  call amr_compute_list_all()
  !
  ! Deallocate temporary arrays
  !
  deallocate(xi,yi,zi)
  !
  ! Done
  !
end subroutine userdef_prep_model


!------------------------------------------------------------------------
! This is the place where you can set up your own model. By the time this
! subroutine is called the grid (at least a basic grid) is already set
! up. You can still modify the basic grid by adding more refinement, but
! to tell the AMR module to reserve space for more grid points you need
! to take matters into your own hand and create and init the base grid
! yourself in the userdef_prep_model() routine above. 
! No example given here, because it would interfere with basic operations.
!------------------------------------------------------------------------
subroutine userdef_setup_model()
  implicit none
  integer :: ix,iy,iz,ierr,index
  double precision :: xc,yc,zc,rc,rho
  type(amr_branch), pointer :: b
  !
  ! Make sure that the dust data is read, because we will need
  ! that in this model
  !
  call read_dustdata(1)
  if(dust_nr_species.ne.1) stop 991
  !
  ! Allocate the dust density array
  !
  allocate(dustdens(1:dust_nr_species,1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate dust density array'
     stop
  endif
  allocate(dust_massdust(1:dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate dust mass array'
     stop
  endif
  !
  ! Fill the dust density array
  !
  do iz=1,amr_grid_nz
     do iy=1,amr_grid_ny
        do ix=1,amr_grid_nx
           b     => amr_grid_branch(ix,iy,iz)%link
           index = b%leafindex
           !xc    = b%xc(1)
           !yc    = b%xc(2)
           !zc    = b%xc(3)
           xc    = amr_finegrid_xc(b%ixyzf(1),1,b%level)
           yc    = amr_finegrid_xc(b%ixyzf(2),2,b%level)
           zc    = amr_finegrid_xc(b%ixyzf(3),3,b%level)
           rc    = sqrt(xc**2+yc**2+zc**2)
           rho   = userdef_rho0 * exp(-rc**2/userdef_radius**2/2.d0)
           dustdens(1,index) = max(rho,1d-99)   ! Make sure the density can never be exactly 0
        enddo
     enddo
  enddo
  !
end subroutine userdef_setup_model


!------------------------------------------------------------------------
! If you want to do some calculation for the model after the main 
! calculation, you can do it here. Here you can also write stuff to
! file. 
!------------------------------------------------------------------------
subroutine userdef_dostuff()
  implicit none
end subroutine userdef_dostuff


!------------------------------------------------------------------------
! If you want to do design your own action (like 'mctherm' or 'image' but
! now designed by you entirely, and activated with 'radmc3d myaction')
! then here is your chance!
!------------------------------------------------------------------------
subroutine userdef_action()
  implicit none
end subroutine userdef_action


!------------------------------------------------------------------------
! If you want to use a Voigt line profile instead of a Gaussian line
! profile, you must initialize the variable lines_ray_lorentz_delta
! here.
! (Added by Thomas Peters 2011)
!------------------------------------------------------------------------
subroutine userdef_compute_lorentz_delta(ray_index)
  implicit none
  integer :: ray_index
end subroutine userdef_compute_lorentz_delta


!------------------------------------------------------------------------
! If you have a good idea how to calculate the level populations 
! of a molecule on-the-fly, you can do it here. But to activate it,
! you must put the line mode to -2.
! IMPORTANT NOTE: If you use the method of selecting a subset of the
!                 levels of a molecule, then you must be very careful
!                 in this subroutine to do it right. You must then use
!                 the "active_***" arrays and variables in the line
!                 module to figure out which levels are "active" and 
!                 which are not. 
!------------------------------------------------------------------------
subroutine userdef_compute_levelpop(ispec,nlevels,index,x,y,z, &
                                    numberdens,levelpop)
  implicit none
  integer :: ispec,nlevels,index
  double precision :: x,y,z,numberdens,levelpop(nlevels)
end subroutine userdef_compute_levelpop


!------------------------------------------------------------------------
! This routine lets you calculate the level populations entirely from
! scratch. Use lines_mode = -10 to use this routine.
!------------------------------------------------------------------------

subroutine userdef_general_compute_levelpop(ray_index, levelpop)
  implicit none
  integer :: ray_index
  double precision :: levelpop(1:lines_nrlevels_subset_max,1:lines_nr_species)
end subroutine userdef_general_compute_levelpop


!------------------------------------------------------------------------
! This subroutine allows you to specify exactly according to your own
! recipes/ideas the emissivity coefficient j_nu [erg/s/cm^3/Hz/ster]
! and extinction coefficient alpha_nu [1/cm] at the wavelengths given
! and at the location given. 
!
! ARGUMENTS:
!  index        The array index of the cell. This allows you to find
!               e.g. the gas density gasdens(index) or the gas
!               temperature gastemp(index) or any other quantity,
!               provided it is read in into the code. 
!  nrfreq       Nr of frequencies of the freq(:) array
!  freq(:)      Array of frequencies [Hz]
!  inu0,inu1    Starting/ending index: Calculate only src(inu0:inu1) and
!               alp(inu0:inu1) belonging to freq(inu0:inu1).
!
! RESULTS:
!  src(:)       Emissivity [erg/s/cm^3/Hz/ster]
!  alp(:)       Extinction [1/cm]
!
! Note: To activate this, you must set incl_userdef_srcalp = 1 in the
!       radmc3d.inp input file (in the code this is the logical 
!       rt_incl_userdef_srcalp from rtglobal_module.f90).
!
! Note: By the time RADMC-3D call this code, it has already computed
!       its own src(:) and alp(:). So just ADD your own values by e.g.
!       src(:) = src(:) + yourstuff and alp(:) = alp(:) + yourstuff.
!       In this way you won't delete the standard stuff. But if you
!       want to replace RADMC-3D's own stuff, you can also just write
!       src(:) = yourstuff and alp(:) = yourstuff. This is up to you.
!------------------------------------------------------------------------
subroutine userdef_srcalp(index,nrfreq,inu0,inu1,freq,src,alp)
  implicit none
  integer :: index,nrfreq,inu0,inu1
  double precision :: freq(1:nrfreq),src(1:nrfreq),alp(1:nrfreq)
end subroutine userdef_srcalp


!------------------------------------------------------------------------
! Here you can write model setup arrays (the stuff you have set up here
! in the userdef_module.f90) to standard RADMC-3D-readable files.
! This will only be done if radmc3d receives the 'writemodel' command
! on the command line. 
!------------------------------------------------------------------------
subroutine userdef_writemodel()
  implicit none
  call write_grid_file()
  call write_dust_density()
end subroutine userdef_writemodel


!------------------------------------------------------------------------
! Reset some action flags for next command?
!------------------------------------------------------------------------
subroutine userdef_reset_flags()
  implicit none
end subroutine userdef_reset_flags


end module userdef_module
