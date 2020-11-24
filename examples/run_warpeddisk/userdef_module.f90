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

integer :: userdef_nx,userdef_ny,userdef_nz,userdef_levelmax
doubleprecision :: userdef_sizex,userdef_sizey,userdef_sizez
doubleprecision :: userdef_r0,userdef_sigmadust0,userdef_gamma,userdef_t10
doubleprecision :: userdef_iw,userdef_aw,userdef_rw,userdef_mstar
doubleprecision :: userdef_nrefinefact,userdef_cell_z_ratio,userdef_zr_min
doubleprecision :: userdef_rin,userdef_drr_max,userdef_rhomin

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
!
! Here we make a typical Lynden-Bell & Pringle disk, but with a warp,
! according to the recipe by Katherine Rosenfeld.
! 
!------------------------------------------------------------------------
subroutine userdef_defaults()
  implicit none
  doubleprecision, parameter :: Msun=1.9889200d+33
  !
  ! Parameters for the base grid
  !
  userdef_nx = 10                    ! Base grid nr of grid points
  userdef_ny = 10
  userdef_nz = 10
  userdef_sizex = 1000 * au          ! Grid box size
  userdef_sizey = 1000 * au
  userdef_sizez = 1000 * au
  !
  ! Parameters for Katherine Rosenfeld's warped disk model
  !
  userdef_mstar = Msun               ! Mass of the star
  userdef_rin = 10 * AU              ! Inner radius of the disk
  userdef_r0 = 100 * AU              ! LB&P critical radius
  userdef_sigmadust0 = 1.d0          ! Surface density factor
  userdef_rhomin = 1d-30             ! Lowest allowed density in g/cm^3
  userdef_gamma = 1.d0               ! LB&P gamma exponent for the viscosity
  userdef_t10 = 100.d0               ! The temperature at 10 AU
  userdef_iw = 40.d0                 ! Katherine Rosenfeld's warp inclination
  userdef_aw = 1.d0                  ! Katherine Rosenfeld's warp powerlaw
  userdef_rw = 10 * AU               ! Katherine Rosenfeld's warp critial radius
  !
  ! Fine-tuning parameters for the AMR grid refinement
  !
  userdef_nrefinefact = 300.d0       ! Estimated fractional increase in nr of
  !                                  ! cells due to AMR refinement
  userdef_levelmax = 9               ! Maximum allowed levels of refinement
  userdef_cell_z_ratio = 0.5         ! Criterion on cellsize/height_above_midplane ratio
  userdef_zr_min = 0.15              ! Limit on cellsize/height_above_midplane ratio
  userdef_drr_max = 0.5              ! Criterion on cellsize/distance_to_star ratio
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
  ! An example how to include an integer keyword. Use
  ! parse_input_double and parse_input_word for reading doubleprecision
  ! variables or string variables. You can delete the below example.
  ! Note that the keyword name string should always have length 30,
  ! hence the many whitespaces.
  !
  ! NOTE: This stuff is not really necessary: You can also change the
  !       parameters in the userdef_defaults() routine above, and
  !       simply recompile. The idea of the parsing here is that you
  !       would not need to recompile each time you want to make a
  !       new model: You could simply change parameters in the 
  !       radmc3d.inp file and rerun radmc3d.
  !
  call parse_input_integer('userdef_nx@                   ',userdef_nx)
  call parse_input_integer('userdef_ny@                   ',userdef_ny)
  call parse_input_integer('userdef_nz@                   ',userdef_nz)
  call parse_input_integer('userdef_levelmax@             ',userdef_levelmax)
  call parse_input_double ('userdef_nrefinefact@          ',userdef_nrefinefact)
  call parse_input_double ('userdef_cell_z_ratio@         ',userdef_cell_z_ratio)
  call parse_input_double ('userdef_zr_min@               ',userdef_zr_min)
  call parse_input_double ('userdef_sizex@                ',userdef_sizex)
  call parse_input_double ('userdef_sizey@                ',userdef_sizey)
  call parse_input_double ('userdef_sizez@                ',userdef_sizez)
  call parse_input_double ('userdef_mstar@                ',userdef_mstar)
  call parse_input_double ('userdef_rin@                  ',userdef_rin)
  call parse_input_double ('userdef_r0@                   ',userdef_r0)
  call parse_input_double ('userdef_sigmadust0@           ',userdef_sigmadust0)
  call parse_input_double ('userdef_rhomin@               ',userdef_rhomin)
  call parse_input_double ('userdef_gamma@                ',userdef_gamma)
  call parse_input_double ('userdef_t10@                  ',userdef_t10)
  call parse_input_double ('userdef_iw@                   ',userdef_iw)
  call parse_input_double ('userdef_aw@                   ',userdef_aw)
  call parse_input_double ('userdef_rw@                   ',userdef_rw)
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
  integer :: i,nbase,nbrtot,nlftot
  !
  ! Tell the code which grid type is used (AMR)
  !
  igrid_type          = 1         ! Oct-tree AMR
  igrid_coord         = 3         ! 3-D cartesian coordinates
  amr_always_use_tree = .true.    ! Make absolutely sure that AMR is active
  !                                 even if the grid is regular. This is
  !                                 necessary for the way we set up the
  !                                 model below.
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
  nbase  = userdef_nx*userdef_ny*userdef_nz  ! Nr of base grid cells
  nlftot = nbase * userdef_nrefinefact       ! Estimated max nr of leafs
  nbrtot = nlftot + nlftot / 6               ! Estimated max nr of branches
  !
  ! Now set up the grid, for now just a regular grid
  !
  call amr_initialize(.true.,.true.,.true.,                           &
       userdef_nx,userdef_ny,userdef_nz,xi,yi,zi,userdef_levelmax,1,  &
       nbrtot,nlftot,.false.,.false.,.false.)
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
  integer :: ix,iy,iz,ierr,index,icell,nrleafs,ispec
  double precision :: xc,yc,zc,rn
  type(amr_branch), pointer :: b
  !
  ! Make sure that the dust data is read, because we will need
  ! that in this model
  !
  call read_dustdata(1)
  !
  ! In this model we are only going to set up the density distribution
  ! of a single dust species, so check that indeed only 1 dust species
  ! is installed.
  !
  if(dust_nr_species.ne.1) stop 991
  !
  ! Allocate the dust density array. The amr_nrleafs_max is the
  ! maximum number of leaf-like cells that we allocated by 
  ! calling amr_initialize(). We called it nlftot in the argument
  ! list of amr_initialize() and amr_initialize() put that value
  ! into the amr_nrleafs_max variable.
  !
  allocate(dustdens(1:dust_nr_species,1:amr_nrleafs_max),STAT=ierr)
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
  ! Fill the dust density array for the basic grid
  !
  do iz=1,amr_grid_nz
     do iy=1,amr_grid_ny
        do ix=1,amr_grid_nx
           b     => amr_grid_branch(ix,iy,iz)%link
           index = b%leafindex
           xc    = amr_finegrid_xc(b%ixyzf(1),1,b%level)
           yc    = amr_finegrid_xc(b%ixyzf(2),2,b%level)
           zc    = amr_finegrid_xc(b%ixyzf(3),3,b%level)
           dustdens(1,index) = userdef_dustdens(xc,yc,zc)
        enddo
     enddo
  enddo
  !
  ! Now let's see if we can refine the grid. Since we have already
  ! a complete list of current cells (by virtue of the call to 
  ! postprocess_grid_for_use() in the main.f90 prior to this call),
  ! we can simply loop over this list.
  !
  nrleafs = amr_nrleafs     ! Need to copy, as amr_nrleafs changes
  do icell=1,nrleafs
     !
     ! Get the branch. Note: The amr_theleafs list is a list of 
     ! cells (leafs) as it was *before* the AMR refinement done here.
     !
     b => amr_theleafs(icell)%link
     !
     ! Call the userdefined recursive refinement routine.
     !
     call userdef_check_and_refine(b)
     !
  enddo
  !
  ! Set refresh the index lists so that we know how to connect each AMR grid
  ! cell to a data memory location. This was already done before (and
  ! even used above), but since we have new cells, we need to redo this.
  !
  call amr_compute_list_all()
  !
  ! Now refresh the grid for use for RADMC-3D.
  ! This routine was already called before in the main.f90 routine, but
  ! that was before the grid refinement, so we have to redo it here.
  !
  call postprocess_grid_for_use(renew=.true.)
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
!            COORDINATE TRANSFORMATION TO LOCAL DISK PLANE
!
! To implement the warp, the easiest way is to simply transform the
! coordinates to the warped plane.
!------------------------------------------------------------------------
subroutine userdef_xyz_to_diskplane(xc,yc,zc,xd,yd,zd)
  implicit none
  double precision :: xc,yc,zc,xd,yd,zd
  double precision :: r,angle,cosa,sina
  r     = sqrt(xc**2+yc**2+zc**2)
  angle = userdef_iw*(pi/180.)*(r/(userdef_rw))**(-userdef_aw)
  cosa  = cos(angle)
  sina  = sin(angle)
  xd    = cosa * xc - sina * zc
  yd    = yc
  zd    = sina * xc + cosa * zc
end subroutine userdef_xyz_to_diskplane


!------------------------------------------------------------------------
!                       MY DENSITY FUNCTION
!------------------------------------------------------------------------
function userdef_dustdens(xc,yc,zc)
  implicit none
  doubleprecision :: userdef_dustdens
  doubleprecision :: rsph,rcyl,xc,yc,zc,xd,yd,zd
  doubleprecision :: temp,hp,cs
  doubleprecision,parameter :: GG=6.6730000d-08  ! Gravitational constant
  doubleprecision,parameter :: kk  = 1.3807d-16  ! Bolzmann's constant     [erg/K]
  doubleprecision,parameter :: mp  = 1.6726d-24  ! Mass of proton          [g]
  !
  ! Convert xc,yc,zc into xd,yd,zd, which are the coordinates
  ! transformed to the local disk plane
  !
  call userdef_xyz_to_diskplane(xc,yc,zc,xd,yd,zd)
  !
  ! Spherical radial coordinate
  !
  rsph = sqrt(xc**2+yc**2+zc**2)
  !
  ! Cylindrical radial coordinate
  !
  rcyl = sqrt(xd**2+yd**2)
  !
  ! Temperature, assuming that it goes as T ~ 1/sqrt(r)
  !
  temp  = userdef_t10 * (rsph/(10*AU))**(-0.5) 
  !
  ! Isothermal sound speed
  !
  cs    = sqrt(kk*temp/(2.3*mp))
  !
  ! Calculate the local scale height
  !
  hp    = cs/sqrt(GG*userdef_mstar/(rcyl**3))
  !
  ! Now the density formula of Lynden-Bell & Pringle, combined
  ! with a Gaussian vertical structure
  !
  if(rsph.ge.userdef_rin) then
     userdef_dustdens = userdef_sigmadust0/(sqrt(twopi)*hp) * &
          (rcyl/userdef_r0)**(-userdef_gamma) *               &
          exp(-(rcyl/userdef_r0)**(2.d0-userdef_gamma)) *     &
          exp(-0.5d0*(zd/hp)**2)
     userdef_dustdens = max(userdef_dustdens,userdef_rhomin)
  else
     userdef_dustdens = userdef_rhomin
  endif
  return
end function userdef_dustdens


!------------------------------------------------------------------------
!                 MY RECURSIVE REFINEMENT SUBROUTINE
!------------------------------------------------------------------------
recursive subroutine userdef_check_and_refine(b)
  implicit none
  type(amr_branch), pointer :: b,c
  double precision :: xc,yc,zc,xd,yd,zd,dx,rsph,rcyl
  integer :: index,ix,iy,iz
  logical :: refine
  !
  ! If this cell is not a leaf, then return
  !
  if(.not.b%leaf) return
  !
  ! If we have reached the maximum refinement level, then return
  !
  if(b%level.ge.userdef_levelmax) then
     return
  endif
  !
  ! Get the cell-center coordinates
  !
  xc    = amr_finegrid_xc(b%ixyzf(1),1,b%level)
  yc    = amr_finegrid_xc(b%ixyzf(2),2,b%level)
  zc    = amr_finegrid_xc(b%ixyzf(3),3,b%level)
  !
  ! Get cell size (since we have cubic grids, we only need this in 
  ! one direction)
  !
  dx    = abs(amr_finegrid_xi(b%ixyzf(1)+1,1,b%level) - &
              amr_finegrid_xi(b%ixyzf(1),1,b%level) )
  !
  ! Convert xc,yc,zc into xd,yd,zd, which are the coordinates
  ! transformed to the local disk plane
  !
  call userdef_xyz_to_diskplane(xc,yc,zc,xd,yd,zd)
  !
  ! Spherical radial coordinate
  !
  rsph = sqrt(xc**2+yc**2+zc**2)
  !
  ! Cylindrical radial coordinate
  !
  rcyl = sqrt(xd**2+yd**2)
  !
  ! Now follows a set of if-statements that check if we should
  ! refine this cell or not
  !
  ! Check if this cell is outside of userdef_rin, because inside
  ! of this radius we should not refine. The sqrt(3)*dx is there
  ! to assure that even if just a corner of the cell lies outside
  ! of userdef_rin, the refinement can take place.
  !
  if(rsph+sqrt(3.d0)*dx.ge.userdef_rin) then
     !
     ! We want that the grid cell size dx obey two criteria: One for
     ! the vertical resolution and one for the radial resolution. If
     ! either of the two is not fulfilled, we refine the cell by splitting
     ! the cell into 2x2x2 subcells.
     !
     refine = .false.
     !
     ! The z-criterion
     !
     if(dx.gt.userdef_cell_z_ratio*max(abs(zd),userdef_zr_min*rsph)) then
        refine = .true.
     endif
     !
     ! The r-criterion
     !
     if(dx.gt.rsph*userdef_drr_max) then
        refine = .true.
     endif
     !
     if(refine) then
        !
        ! Yep, we need to refine
        !
        ! Now call the refinement routine of the AMR module
        !
        call amr_branch_refine(b,0)
        !
        ! Now reinstall the density in the new cells, and recursively
        ! call this refinement routine
        !
        do iz=1,2
           do iy=1,2
              do ix=1,2
                 !
                 ! Get the child branch
                 !
                 c => b%child(ix,iy,iz)%link
                 !
                 ! Reinstall the density
                 !
                 index = c%leafindex
                 xc    = amr_finegrid_xc(c%ixyzf(1),1,c%level)
                 yc    = amr_finegrid_xc(c%ixyzf(2),2,c%level)
                 zc    = amr_finegrid_xc(c%ixyzf(3),3,c%level)
                 dustdens(1,index) = userdef_dustdens(xc,yc,zc)
                 !
                 ! Recursively call the refiner
                 !
                 call userdef_check_and_refine(c)
                 !
              enddo
           enddo
        enddo
     endif
  endif
  !
end subroutine userdef_check_and_refine


!------------------------------------------------------------------------
! Reset some action flags for next command?
!------------------------------------------------------------------------
subroutine userdef_reset_flags()
  implicit none
end subroutine userdef_reset_flags


end module userdef_module
