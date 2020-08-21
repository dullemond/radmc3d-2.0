!==========================================================================
!           Module for tracking straight line through AMR grid
!              and for other formal radiative transfer stuff
!==========================================================================
!
! Before using this AMRRay module you must have set up a grid with the AMR
! module. This is done by calling amr_initialize() with the appropriate
! arguments. You can also call (one of the?) reading routines of the AMR
! library so that the initialization is done for you. But it is important to
! set nrbranchesmax to the maximum number of branches you think you need.
! This is important because the AMRRay module will allocate some data 
! arrays for internal use, and it has to know how many branches there will
! be at maximum. If you set up a grid with, say, 1000 branches, but you
! wish to be able to add branches/leafs later by adaptive grid refinement,
! then you must specify nrbranchesmax to the max number you think you need.
! But do not take it TOO large, because this will waste computer memory.
!
! Now make sure that the amr_coordsystem value is set to any value between:
!     0 <= amr_coordsystem < 100    = Cartesian
!   100 <= amr_coordsystem < 200    = Spherical
!   200 <= amr_coordsystem < 300    = Cylindrical [not yet implemented]
! The reason why there is this funny range of values all meaning the same
! coordinate system is because the author of this module also made a 
! hydrodynamic package which has more closely specified kinds of coordinate
! systems, and we may in the future have compatibility. So for now you can
! simply take 0, 100 or 200 for cartesian, spherical and cylindrical. 
! Internally in the amrray module this will be translated into 0, 1 or 2
! in the amrray_icoord variable. 
!
! After setting up the AMR grid, the next thing to do is to call the routine
! amrray_initialize(). This will set up all precalculated stuff required for
! the ray-tracing (especially in spherical coordinates this is necessary). 
! And it will perform a series of checks on the grid setup to make sure
! that it is conform what is expected for the AMRRay module. Also this is
! mainly for the spherical coordinates.
!
! Now you can start using the module. Once you want to end the program,
! you should call amrray_cleanup() to deallocate all the allocated arrays.
!
!
! NOTE: The xc values of each branch MUST be exactly in the middle of
!       the cell: a%xc(idir) = 0.5d0 * ( a%xi(1,idir) + a%xi(2,idir) ).
!       This is default in the AMR module, but not required in the AMR
!       module. But it is strictly required here!
!--------------------------------------------------------------------------
module amrray_module
use amr_module
use constants_module
implicit none
public

logical :: amrray_selfcheck = .true.
logical :: amrray_mirror_equator = .false.
integer :: amrray_initialized = 0
integer :: amrray_icoord = 0
integer :: theta1_sgn = 1
integer :: theta2_sgn = -1

integer :: amrray_icross

type(amr_branch), pointer :: amrray_cell,amrray_nextcell
integer :: amrray_ix_curr,amrray_iy_curr,amrray_iz_curr,amrray_ix_next,amrray_iy_next,amrray_iz_next

doubleprecision :: tiny,small,tol
parameter(tiny=1d-14)
parameter(small=1d-12)
!parameter(tol=1d-12)
parameter(tol=1d-10)   ! CHANGE 10.03.2015

doubleprecision :: eps_thres
parameter(eps_thres=1d-4)

doubleprecision :: oneplus,oneminus
doubleprecision :: onepluss,oneminuss
doubleprecision :: oneplust,oneminust
parameter(oneplus =1.00000000000001d0)
parameter(oneminus=0.99999999999999d0)
parameter(onepluss =1.000000000001d0)
parameter(oneminuss=0.999999999999d0)
parameter(oneplust =1.00000001d0)
parameter(oneminust=0.99999999d0)

integer :: pixmap_nx,pixmap_ny
doubleprecision, allocatable :: pixmap_xi(:),pixmap_yi(:)
doubleprecision, allocatable :: pixmap_xc(:),pixmap_yc(:)
doubleprecision, allocatable :: pixmap_cart_xs(:,:)
doubleprecision, allocatable :: pixmap_cart_ys(:,:)
doubleprecision, allocatable :: pixmap_cart_zs(:,:)
doubleprecision, allocatable :: pixmap_cart_dxs(:,:)
doubleprecision, allocatable :: pixmap_cart_dys(:,:)
doubleprecision, allocatable :: pixmap_cart_dzs(:,:)

doubleprecision, allocatable :: amrray_finegrid_sintsq1(:,:)
doubleprecision, allocatable :: amrray_finegrid_sintsq2(:,:)
doubleprecision, allocatable :: amrray_finegrid_costsq1(:,:)
doubleprecision, allocatable :: amrray_finegrid_costsq2(:,:)
doubleprecision, allocatable :: amrray_finegrid_sinp1(:,:)
doubleprecision, allocatable :: amrray_finegrid_sinp2(:,:)
doubleprecision, allocatable :: amrray_finegrid_cosp1(:,:)
doubleprecision, allocatable :: amrray_finegrid_cosp2(:,:)

doubleprecision :: amrray_sint1,amrray_sint2,amrray_cost1,amrray_cost2
doubleprecision :: amrray_sinp1,amrray_sinp2,amrray_cosp1,amrray_cosp2
logical :: thetagrid_cross_equator

integer :: amrray_spheres_nr=0,amrray_ispherehit=0
doubleprecision, allocatable :: amrray_spheres_r(:),amrray_spheres_pos(:,:)
logical, allocatable :: amrray_spheres_outsidegrid(:)
integer, allocatable :: amrray_spheres_sphidx(:)
logical :: amrray_spheres_flag


!$OMP THREADPRIVATE(amrray_icross)
!$OMP THREADPRIVATE(amrray_cell,amrray_nextcell)
!$OMP THREADPRIVATE(amrray_ix_curr,amrray_iy_curr,amrray_iz_curr,amrray_ix_next,amrray_iy_next,amrray_iz_next)
!$OMP THREADPRIVATE(amrray_sint1,amrray_sint2,amrray_cost1,amrray_cost2)
!$OMP THREADPRIVATE(amrray_sinp1,amrray_sinp2,amrray_cosp1,amrray_cosp2)
!$OMP THREADPRIVATE(thetagrid_cross_equator)
!$OMP THREADPRIVATE(amrray_ispherehit)

!doubleprecision :: pi,pihalf,twopi
!parameter(pi=3.141592653589793238d0)
!parameter(twopi=6.283185307179586232d0)
!parameter(pihalf=1.570796326794896558d0)

contains


!--------------------------------------------------------------------------
!                 AMR Ray Initialization routine
!
! Should be called after the gridding is set with the amr_module.f90.
! This does not do much, but it just checks if the gridding is ok, because
! the AMR Ray module relies on certain things to be alright. Also it
! sets certain settings.
!
! NOTE: The amr_coordsystem tells which coordinate system is to be used.
!       Please make sure to set this properly before calling this 
!       amrray_initialize() routine. Any value between the below limits
!       have the below meaning:
!           0 <= amr_coordsystem < 100    = Cartesian
!         100 <= amr_coordsystem < 200    = Spherical
!         200 <= amr_coordsystem < 300    = Cylindrical
!
! ARGUMENTS:
!   selfcheck       If set, then the AMR Ray library is in safe-mode, i.e.
!                   it will always check if intermediate results are OK.
!                   May be a bit slower, but especially in the first 
!                   phases it may be very important. 
!   eqmirror        If set, then the equatorial plane is a mirror, so
!                   that we only should set the upper half of the grid.
!                   (only for spherical coordinates)
!
!--------------------------------------------------------------------------
subroutine amrray_initialize(selfcheck,eqmirror)
implicit none
logical :: selfcheck,ok
logical,optional :: eqmirror
integer :: i,ib,ierr,nn,ilevel
doubleprecision :: theta1,theta2,phi1,phi2
!
! Translate the complex amr_coordsystem into a simple:
!   0 = Cartesian
!   1 = Spherical
!   2 = Cylindrical
! for internal use inside the amrray module
!
if((amr_coordsystem.ge.0).and.(amr_coordsystem.lt.100)) then
   amrray_icoord = 0
elseif((amr_coordsystem.ge.100).and.(amr_coordsystem.lt.200)) then
   amrray_icoord = 1
elseif((amr_coordsystem.ge.200).and.(amr_coordsystem.lt.300)) then
   amrray_icoord = 2
   write(stdo,*) 'ERROR: For now the cylindrical coordinate system is not yet built in...'
   stop
else
   write(stdo,*) 'ERROR: Do not know amr_coordsystem = ',amr_coordsystem
   stop
endif
!
! Set settings
!
amrray_selfcheck      = selfcheck
amrray_ispherehit     = 0
if(present(eqmirror)) then
   !
   ! Set mirroring by hand
   !
   amrray_mirror_equator = eqmirror
   if(amrray_mirror_equator) then
      write(stdo,*) 'Using mirror symmetry in equatorial plane.'
   endif
else
   !
   ! If largest theta == pihalf, then we switch on the equatorial
   ! mirroring. 
   !
   if((amr_coordsystem.ge.100).and.(amr_coordsystem.lt.200).and.    &
        (amr_grid_xi(amr_grid_ny+1,2).eq.pihalf).and.               &
        (amr_ydim.eq.1)) then
      amrray_mirror_equator = .true.
      write(stdo,*) 'Using mirror symmetry in equatorial plane, because max(theta)==pi/2.'
   else
      amrray_mirror_equator = .false.
   endif
endif
!
! Make sure the AMR linear branch lists are made
!
if(amr_tree_present) then
   if(.not.allocated(amr_thebranches)) then
      call amr_compute_list_all()
   endif
endif
!
! Check some stuff
!
if(amrray_mirror_equator.and.(amrray_icoord.ne.1)) then
   write(stdo,*) 'ERROR in AMRRay Module: equatorial plane mirror'
   write(stdo,*) '    can only be used in spherical coordinates'
   stop 8371
endif
!
! Now a whole lot of stuff for the spherical coordinates 
!
if(amrray_icoord.eq.1) then 
   !
   ! Check that all coordinates are increasing, and that the spacing is
   ! not too small.
   !
   do i=1,amr_grid_nx
      if(amr_grid_xi(i+1,1).le.amr_grid_xi(i,1)) then
         write(stdo,*) 'ERROR: Radial grid must be in ascending order'
         stop
      endif
      if(amr_grid_xi(i+1,1)/amr_grid_xi(i,1)-1.d0.le.4*small) then
         write(stdo,*) 'ERROR: Radial grid is too finely spaced at ir=',i
         stop
      endif
   enddo
   if(amr_ydim.ne.0) then
      do i=1,amr_grid_ny
         if(amr_grid_xi(i+1,2).le.amr_grid_xi(i,2)) then
            write(stdo,*) 'ERROR: Theta grid must be in ascending order'
            stop
         endif
         if(amr_grid_xi(i+1,2)-amr_grid_xi(i,2).le.4*small) then
            write(stdo,*) 'ERROR: Theta grid is too finely spaced at itheta=',i
            stop
         endif
      enddo
   endif
   if(amr_zdim.ne.0) then
      do i=1,amr_grid_nz
         if(amr_grid_xi(i+1,3).le.amr_grid_xi(i,3)) then
            write(stdo,*) 'ERROR: Phi grid must be in ascending order'
            stop
         endif
         if(amr_grid_xi(i+1,3)-amr_grid_xi(i,3).le.4*small) then
            write(stdo,*) 'ERROR: Phi grid is too finely spaced at iphi=',i
            stop
         endif
      enddo
   else
      if(amr_grid_nz.ne.1) then
         stop 7021
      endif
      if((amr_grid_xi(1,3).ne.0.d0).or.(amr_grid_xi(2,3).ne.twopi)) then
         write(stdo,*) 'ERROR in Amrray module: if phi direction is not'
         write(stdo,*) '    active (i.e. 2-D axisymmetric RT), then the'
         write(stdo,*) '    phi-grid must be precisely [0,2*pi].'
         write(stdo,*) '    Use the constant called twopi for this from the'
         write(stdo,*) '    constants_module.f90 module.'
         stop
      endif
   endif
   !
   ! For now we enforce that the smallest theta is always < pi/2.
   !
   if(amr_grid_xi(1,2).ge.pihalf) then
      write(stdo,*) 'ERROR in Amrray module: Lowest theta grid point '
      write(stdo,*) '      must always be < pi/2.'
      stop 3901
   endif
   !
   ! Catch the dangerous situation in which the lower r is zero.
   !
   if(amr_grid_xi(1,1).le.0.d0) then
      write(stdo,*) 'ERROR in Amrray module: smallest r is at 0.d0!'
      write(stdo,*) '      This is not allowed. R must always be positive!'
      stop 5091
   endif
   !
   ! Catch the dangerous situation in which the upper theta is close to,
   ! but not exactly equal to pi/2. 
   !
   if((amr_grid_xi(amr_grid_ny+1,2).ne.pihalf).and.  &
      (abs(amr_grid_xi(amr_grid_ny+1,2)-pihalf).lt.4*small)) then
      write(stdo,*) 'ERROR in Amrray module: largest theta lies very '
      write(stdo,*) '      close to the equatorial plane, but not exactly'
      write(stdo,*) '      on it. If you want it to be at the equator, then'
      write(stdo,*) '      use the value of pihalf from the module'
      write(stdo,*) '      constants_module.f90 to make Amrray recognize it.'
      write(stdo,*) '      Otherwise make sure to put this theta at least '
      write(stdo,*) '      ',5*small,' away from pihalf.'
      stop 5092
   endif
   !
   ! Catch the dangerous situation in which the lower theta is close to,
   ! but not exactly equal to 0.
   !
   if((amr_grid_xi(1,2).ne.0.d0).and.  &
      (abs(amr_grid_xi(1,2)).lt.4*small)) then
      write(stdo,*) 'ERROR in Amrray module: smallest theta lies very '
      write(stdo,*) '      close to the positive z-axis, but not exactly'
      write(stdo,*) '      on it. If you want it to be at theta=0, then'
      write(stdo,*) '      put it exactly to 0.d0.'
      write(stdo,*) '      Otherwise make sure to put this theta at least '
      write(stdo,*) '      ',5*small,' larger than 0.'
      stop 5093
   endif
   !
   ! Catch the dangerous situation in which the upper theta is close to,
   ! but not exactly equal to pi. 
   !
   if((amr_grid_xi(amr_grid_ny+1,2).ne.pi).and.  &
      (abs(amr_grid_xi(amr_grid_ny+1,2)-pi).lt.4*small)) then
      write(stdo,*) 'ERROR in Amrray module: largest theta lies very '
      write(stdo,*) '      close to the negative z-axis, but not exactly'
      write(stdo,*) '      on it. If you want it to be at theta=pi, then'
      write(stdo,*) '      use the value of pi from the module'
      write(stdo,*) '      constants_module.f90 to make Amrray recognize it.'
      write(stdo,*) '      Otherwise make sure to put this theta at least '
      write(stdo,*) '      ',5*small,' smaller than pi.'
      stop 5094
   endif
   !
   ! Catch the dangerous situation in which the lower phi is close to,
   ! but not exactly equal to 0.
   !
   if((amr_grid_xi(1,3).ne.0.d0).and.  &
      (abs(amr_grid_xi(1,3)).lt.4*small)) then
      write(stdo,*) 'ERROR in Amrray module: smallest phi lies very '
      write(stdo,*) '      close to 0, but not exactly'
      write(stdo,*) '      on it. If you want it to be at 0, then'
      write(stdo,*) '      use the precise value 0.d0.'
      write(stdo,*) '      Otherwise make sure to put this theta at least '
      write(stdo,*) '      ',5*small,' larger than 0.d0.'
      stop 5095
   endif
   !
   ! Catch the dangerous situation in which the upper phi is close to,
   ! but not exactly equal to 2*pi. 
   !
   if((amr_grid_xi(amr_grid_nz+1,3).ne.twopi).and.  &
      (abs(amr_grid_xi(amr_grid_nz+1,3)-twopi).lt.4*small)) then
      write(stdo,*) 'ERROR in Amrray module: largest phi lies very '
      write(stdo,*) '      close to 2*pi, but not exactly'
      write(stdo,*) '      on it. If you want it to be at 2*pi, then'
      write(stdo,*) '      use the value of twopi from the module'
      write(stdo,*) '      constants_module.f90 to make Amrray recognize it.'
      write(stdo,*) '      Otherwise make sure to put this theta at least '
      write(stdo,*) '      ',5*small,' smaller than twopi.'
      stop 5096
   endif
   !
   ! Check that in spherical coordinates no cell crosses the equatorial plane.
   ! This means that if the theta grid crosses the midplane, then one of the
   ! theta cell interfaces must be identical to pi/2=pihalf.
   !
   ! If cells would cross the midplane, this would make it difficult to
   ! solve for theta-crossings within each cell. Also check that in case of
   ! mirror symmetry, the last theta is exactly on the equatorial plane.
   !
   if(amr_ydim.eq.1) then
      if(amrray_mirror_equator) then
         !
         ! Mirror symmetry: check that last theta is on midplane
         !
         if(amr_grid_xi(amr_grid_ny+1,2).ne.pihalf) then
            write(stdo,*) 'ERROR in amrray_module:'
            write(stdo,*) '      In spherical coordinates with mirror symmetry'
            write(stdo,*) '      in the equatorial plane, the largest theta value'
            write(stdo,*) '      must be identical to pi/2 (= pihalf from the'
            write(stdo,*) '      constants_module.f90 module).'
            stop
         endif
      else
         !
         ! No mirror symmetry: check that no cell crosses the midplane
         !
         if((amr_grid_xi(1,2).lt.pihalf).and.                        &
              (amr_grid_xi(amr_grid_ny+1,2).gt.pihalf)) then
            !
            ! There should be at least one coordinate that is exactly pihalf
            !
            ok = .false.
            do i=2,amr_grid_ny 
               if(amr_grid_xi(i,2).eq.pihalf) ok = .true.
            enddo
            if(.not.ok) then
               write(stdo,*) 'ERROR in amrray_module:'
               write(stdo,*) '      Using spherical coordinates where the theta'
               write(stdo,*) '      grid crosses the equatorial plane, but no'
               write(stdo,*) '      cell boundary at pi/2 is found.'
               write(stdo,*) '      You MUST make sure that no cell crosses the'
               write(stdo,*) '      equatorial plane. For at least one cell the'
               write(stdo,*) '      largest theta must be pi/2 (= pihalf from the'
               write(stdo,*) '      constants_module.f90 module).'
               stop
            endif
         endif
      endif
   endif
   !
   ! Set the theta1_sgn and theta2_sgn values properly. They should say if
   ! the grid boundaries in theta are above or below the equatorial plane.
   ! Also set the thetagrid_cross_equator flag
   !
   theta1_sgn = 1
   if((amr_grid_xi(amr_grid_ny+1,2).gt.pihalf).or.amrray_mirror_equator) then
      thetagrid_cross_equator = .true.
      theta2_sgn = -1
   else
      thetagrid_cross_equator = .false.
      theta2_sgn = 1
   endif
   !
   ! Do more checks
   !
   if((amr_grid_xi(1,2).lt.0.d0).or.&
      (amr_grid_xi(amr_grid_ny+1,2).gt.pi)) then
      write(stdo,*) 'ERROR in AMRRay Module: Theta grid cannot exceed [0,pi]'
      stop 8373
   endif
   if(amr_grid_xi(1,2).ge.pihalf) then
      write(stdo,*) 'ERROR in AMRRay Module: Lowest Theta must be <pi/2'
      stop 8362
   endif
   if((amr_grid_xi(1,3).lt.0.d0).or.&
      (amr_grid_xi(amr_grid_nz+1,3).gt.twopi)) then
      write(stdo,*) 'ERROR in AMRRay Module: Phi grid cannot exceed [0,2*pi]'
      write(stdo,*) '    The constant called twopi from the'
      write(stdo,*) '    constants_module.f90 module is used for this.'
      stop 8374
   endif
   if(amr_zdim.eq.1) then
      if(amr_cyclic_xyz(3)) then
         if(amr_grid_xi(1,3).ne.0.d0) then
            write(stdo,*) 'ERROR in AMRRay Module: Phi is supposed to be cyclic'
            write(stdo,*) '   (as amr_cyclic_xyz(3) is set), but the lower'
            write(stdo,*) '   phi boundary is not 0.d0.'
            stop 8375
         endif
         if(amr_grid_xi(amr_grid_nz+1,3).ne.twopi) then
            write(stdo,*) 'ERROR in AMRRay Module: Phi is supposed to be cyclic'
            write(stdo,*) '   (as amr_cyclic_xyz(3) is set), but the upper'
            write(stdo,*) '   phi boundary is not 2*pi.'
            write(stdo,*) '   The constant called twopi from the'
            write(stdo,*) '   constants_module.f90 module is used for this.'
            stop 8376
         endif
      else
         if((abs(amr_grid_xi(1,3)).le.small).and.&
            (abs(amr_grid_xi(amr_grid_nz+1,3)-twopi).le.small)) then
            write(stdo,*) 'ERROR in AMRRay Module: Phi is supposed to be NOT cyclic'
            write(stdo,*) '   (as amr_cyclic_xyz(3) is unset), but the lower'
            write(stdo,*) '   phi boundary is 0.d0 and the upper phi boundary'
            write(stdo,*) '   is 2*pi. Please switch on cyclic bc for z in AMR module.'
            stop 8377
         endif
      endif
   else
      if(amr_grid_xi(1,3).ne.0.d0) then
         write(stdo,*) 'ERROR in AMRRay Module: Phi-dimension not used,'
         write(stdo,*) '   but phi lower boundary is not 0.d0.'
         stop 8378
      endif
      if(amr_grid_xi(amr_grid_nz+1,3).ne.twopi) then
         write(stdo,*) 'ERROR in AMRRay Module: Phi-dimension not used,'
         write(stdo,*) '   but phi upper boundary is not 2*pi.'
         write(stdo,*) '   The constant called twopi from the'
         write(stdo,*) '   constants_module.f90 module is used for this.'
         stop 8379
      endif
   endif
   if(amr_grid_xi(amr_grid_nx+1,1)/amr_grid_xi(1,1).gt.1d10) then
      write(stdo,*) 'ERROR in AMRRay Module: Radial dynamic range too large.'
      stop 8380
   endif
   !
   ! Allocate some arrays for the spherical coordinates mode.
   !
   if(.not.amr_use_index) then
      write(stdo,*) 'ERROR in AMRRay Module: In AMR module the indexing system is not used.'
      write(stdo,*) '   You should have set nrbranchesmax in the amr_initialize() routine.'
      stop
   endif
   if(amr_nrbranches_max.lt.amr_nrbranches) then 
      write(stdo,*) 'ERROR in AMRRay Module: The amr_nrbranches_max is not set'
      stop
   endif
   allocate(amrray_finegrid_sintsq1(1:amr_nxyzfmax,0:amr_levelmax),STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in AMRRay Module: Could not allocate amrray_finegrid_sintsq1 array'
      stop 
   endif
   allocate(amrray_finegrid_sintsq2(1:amr_nxyzfmax,0:amr_levelmax),STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in AMRRay Module: Could not allocate amrray_finegrid_sintsq2 array'
      stop 
   endif
   allocate(amrray_finegrid_costsq1(1:amr_nxyzfmax,0:amr_levelmax),STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in AMRRay Module: Could not allocate amrray_finegrid_costsq1 array'
      stop 
   endif
   allocate(amrray_finegrid_costsq2(1:amr_nxyzfmax,0:amr_levelmax),STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in AMRRay Module: Could not allocate amrray_finegrid_costsq2 array'
      stop 
   endif
   allocate(amrray_finegrid_sinp1(1:amr_nxyzfmax,0:amr_levelmax),STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in AMRRay Module: Could not allocate amrray_finegrid_sinp1 array'
      stop 
   endif
   allocate(amrray_finegrid_sinp2(1:amr_nxyzfmax,0:amr_levelmax),STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in AMRRay Module: Could not allocate amrray_finegrid_sinp2 array'
      stop 
   endif
   allocate(amrray_finegrid_cosp1(1:amr_nxyzfmax,0:amr_levelmax),STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in AMRRay Module: Could not allocate amrray_finegrid_cosp1 array'
      stop 
   endif
   allocate(amrray_finegrid_cosp2(1:amr_nxyzfmax,0:amr_levelmax),STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in AMRRay Module: Could not allocate amrray_finegrid_cosp2 array'
      stop 
   endif
   !
   ! Pre-compute sin(theta) and cos(theta) for grid edges.
   !
   if(amrray_mirror_equator) then
      amrray_sint1 = sin(amr_grid_xi(1,2))
      amrray_sint2 = sin(pi-amr_grid_xi(1,2)) 
      amrray_cost1 = cos(amr_grid_xi(1,2))
      amrray_cost2 = cos(pi-amr_grid_xi(1,2)) 
   else
      amrray_sint1 = sin(amr_grid_xi(1,2))
      amrray_sint2 = sin(amr_grid_xi(amr_grid_ny+1,2)) 
      amrray_cost1 = cos(amr_grid_xi(1,2))
      amrray_cost2 = cos(amr_grid_xi(amr_grid_ny+1,2)) 
   endif
   if(amr_zdim.eq.0) then
      amrray_sinp1 = 0.d0
      amrray_sinp2 = 0.d0
      amrray_cosp1 = 1.d0
      amrray_cosp2 = 1.d0
   else
      amrray_sinp1 = sin(amr_grid_xi(1,3))
      amrray_sinp2 = sin(amr_grid_xi(amr_grid_nz+1,3)) 
      amrray_cosp1 = cos(amr_grid_xi(1,3))
      amrray_cosp2 = cos(amr_grid_xi(amr_grid_nz+1,3)) 
   endif
   !
   ! Now pre-compute the sin(theta) and cos(theta) for the whole grid. 
   ! 
   ! ...Loop over levels
   !
   do ilevel=0,amr_levelmax
      !
      ! Do the loop over theta
      !
      if(amr_ydim.eq.1) then
         nn = amr_grid_ny*(2**ilevel)
      else
         nn = 1
      endif
      do i=1,nn
         !
         ! Assure that theta1 is always the one farthest away from the equator
         ! and theta2 the closest to the equator.
         !
         if(amr_finegrid_xi(i,2,ilevel).lt.pihalf) then
            theta1 = amr_finegrid_xi(i,2,ilevel)
            theta2 = amr_finegrid_xi(i+1,2,ilevel)
         else
            theta1 = amr_finegrid_xi(i+1,2,ilevel)
            theta2 = amr_finegrid_xi(i,2,ilevel)
         endif
         if(theta1.eq.pihalf) stop 8993
         !
         ! Now compute the squared sine and cosines
         !
         amrray_finegrid_sintsq1(i,ilevel) = sin(theta1)**2
         amrray_finegrid_sintsq2(i,ilevel) = sin(theta2)**2
         amrray_finegrid_costsq1(i,ilevel) = 1.d0 - amrray_finegrid_sintsq1(i,ilevel)
         amrray_finegrid_costsq2(i,ilevel) = 1.d0 - amrray_finegrid_sintsq2(i,ilevel)
         !
         ! Exception handling
         !
         if(theta2.eq.pihalf) then
            amrray_finegrid_sintsq2(i,ilevel) = 1.d0
            amrray_finegrid_costsq2(i,ilevel) = 0.d0
         endif
         if(theta1.eq.0.d0) then
            amrray_finegrid_sintsq1(i,ilevel) = 0.d0
            amrray_finegrid_costsq1(i,ilevel) = 1.d0
         endif
         if(theta2.eq.pi) then
            amrray_finegrid_sintsq2(i,ilevel) = 0.d0
            amrray_finegrid_costsq2(i,ilevel) = -1.d0
         endif
      enddo
      !
      ! Do the loop over phi
      !
      if(amr_zdim.eq.1) then
         nn = amr_grid_nz*(2**ilevel)
      else
         nn = 1
      endif
      do i=1,nn
         !
         ! Get the phi values
         !
         phi1   = amr_finegrid_xi(i,3,ilevel)
         phi2   = amr_finegrid_xi(i+1,3,ilevel)
         !
         ! Compute the sines and cosines
         !
         amrray_finegrid_sinp1(i,ilevel)   = sin(phi1)
         amrray_finegrid_sinp2(i,ilevel)   = sin(phi2)
         amrray_finegrid_cosp1(i,ilevel)   = cos(phi1)
         amrray_finegrid_cosp2(i,ilevel)   = cos(phi2)
      enddo
   enddo
   !
endif
amrray_initialized = 12345
end subroutine amrray_initialize


!--------------------------------------------------------------------------
!                  INITIALIZE THE USE OF SPHERES
!
! In addition to the crossings with grid coordinates you can also make
! amrray check for crossings with a discrete set of spheres of different
! sizes. This is a way to include the finite sizes of stars into the model.
! Of course, checking also for this can slow down the code. So we have to
! do this cleverly. Firstly, if you have assured that all spheres are
! entirely outside the grid, then the checking only needs to be done
! when the ray is indeed outside of the grid. If one or more spheres are
! inside the grid, of have part of it lie inside the grid, then we must
! introduce the amrray_spheres_sphidx array to speed things up.
! 
! Before calling this subroutine the user must have allocated and filled
! the amrray_spheres_r(:) and amrray_spheres_pos(:,:) arrays and set 
! amrray_spheres_nr to the number of spheres. 
!--------------------------------------------------------------------------
subroutine amrray_install_spheres()
implicit none
integer :: isph,ierr,icell,isph1,iddr
logical :: flag,centralsphere
type(amr_branch), pointer :: leaf
doubleprecision :: r,theta,phi,rmin,rmax,thetamin,thetamax,phimin,phimax
doubleprecision :: sintheta,dist
doubleprecision :: axi(1:2,1:3)
integer :: ix,iy,iz
!
! Do a self-consistency check
!
if((.not.allocated(amrray_spheres_r)).or.    &
   (.not.allocated(amrray_spheres_pos)).or.  &
   (amrray_spheres_nr.le.0)) then
   write(stdo,*) 'ERROR: If you wish to use the spheres in the amrray module'
   write(stdo,*) '       then you must allocate amrray_spheres_r,pos...'
   stop
endif
if(allocated(amrray_spheres_outsidegrid)) deallocate(amrray_spheres_outsidegrid)
allocate(amrray_spheres_outsidegrid(amrray_spheres_nr),STAT=ierr)
if(ierr.ne.0) then
   write(stdo,*) 'ERROR: Could not allocate amrray_spheres_outsidegrid().'
   stop
endif
amrray_spheres_outsidegrid(:) = .false.
!
! Now check if all spheres are sufficiently clearly outside of the grid
! to NOT require the allocation of amrray_spheres_sphidx(:)
!
if(allocated(amrray_spheres_sphidx)) deallocate(amrray_spheres_sphidx)
if(amr_coordsystem.lt.100) then
   !
   ! Cartesian coordinates
   !
   ! Check if one or more spheres may overlap with the grid. Note that
   ! this is not a perfect criterion, and we might raise false alarm.
   ! But better this way around than not raise alarm when a sphere
   ! does intersect with the grid.
   !
   flag = .false.
   do isph=1,amrray_spheres_nr
      if((amrray_spheres_pos(1,isph).ge.                             &
          amr_grid_xi(1,1)-amrray_spheres_r(isph)).and.              &
         (amrray_spheres_pos(1,isph).le.                             &
          amr_grid_xi(amr_grid_nx+1,1)+amrray_spheres_r(isph)).and.  &
         (amrray_spheres_pos(2,isph).ge.                             &
          amr_grid_xi(1,2)-amrray_spheres_r(isph)).and.              &
         (amrray_spheres_pos(2,isph).le.                             &
          amr_grid_xi(amr_grid_ny+1,2)+amrray_spheres_r(isph)).and.  &
         (amrray_spheres_pos(3,isph).ge.                             &
          amr_grid_xi(1,3)-amrray_spheres_r(isph)).and.              &
         (amrray_spheres_pos(3,isph).le.                             &
          amr_grid_xi(amr_grid_nz+1,3)+amrray_spheres_r(isph))) then
         flag = .true.
      endif
      if((amrray_spheres_pos(1,isph).ge.                             &
          amr_grid_xi(1,1)+amrray_spheres_r(isph)).and.              &
         (amrray_spheres_pos(1,isph).le.                             &
          amr_grid_xi(amr_grid_nx+1,1)-amrray_spheres_r(isph)).and.  &
         (amrray_spheres_pos(2,isph).ge.                             &
          amr_grid_xi(1,2)+amrray_spheres_r(isph)).and.              &
         (amrray_spheres_pos(2,isph).le.                             &
          amr_grid_xi(amr_grid_ny+1,2)-amrray_spheres_r(isph)).and.  &
         (amrray_spheres_pos(3,isph).ge.                             &
          amr_grid_xi(1,3)+amrray_spheres_r(isph)).and.              &
         (amrray_spheres_pos(3,isph).le.                             &
          amr_grid_xi(amr_grid_nz+1,3)-amrray_spheres_r(isph))) then
         amrray_spheres_outsidegrid(isph) = .false.
      else
         amrray_spheres_outsidegrid(isph) = .true.
      endif
   enddo
   !
   ! If we have one or more spheres overlapping with the grid, then we must
   ! allocate the amrray_spheres_sphidx(:) array
   !
   if(flag) then
      !
      ! Allocate the array
      !
      allocate(amrray_spheres_sphidx(1:amr_nrleafs_max),STAT=ierr)
      if(ierr.ne.0) then
         write(stdo,*) 'ERROR in amrray: Could not allocate amrray_spheres_sphidx array'
         stop
      endif
      amrray_spheres_sphidx(:) = 0
      !
      ! Now check out cell-by-cell 
      !
      if(amr_tree_present) then
         if(.not.allocated(amr_theleafs)) stop 5673
         if(.not.allocated(amr_theleaf_index)) stop 5674
      endif
      do icell=1,amr_nrleafs
         !
         ! Get pointers to this cell and cell walls
         !
         if(amr_tree_present) then
            !
            ! If AMR tree available, then use that
            !
            !!index = amr_theleaf_index(icell)
            leaf  => amr_theleafs(icell)%link
            do iddr=1,3
               axi(1,iddr) = amr_finegrid_xi(leaf%ixyzf(iddr),iddr,leaf%level)
               axi(2,iddr) = amr_finegrid_xi(leaf%ixyzf(iddr)+1,iddr,leaf%level)
            enddo
         else
            !
            ! If regular grid, then use ix, iy, iz
            !
            call amr_regular_get_ixyz(icell,ix,iy,iz)
            axi(1,1) = amr_finegrid_xi(ix,1,0)
            axi(2,1) = amr_finegrid_xi(ix+1,1,0)
            axi(1,2) = amr_finegrid_xi(iy,2,0)
            axi(2,2) = amr_finegrid_xi(iy+1,2,0)
            axi(1,3) = amr_finegrid_xi(iz,3,0)
            axi(2,3) = amr_finegrid_xi(iz+1,3,0)
         endif
         !
         ! Now check out all the spheres
         !
         do isph=1,amrray_spheres_nr
            if((amrray_spheres_pos(1,isph).ge.                    &
                 axi(1,1)-amrray_spheres_r(isph)).and.            &
                 (amrray_spheres_pos(1,isph).le.                  &
                 axi(2,1)+amrray_spheres_r(isph)).and.            &
                 (amrray_spheres_pos(2,isph).ge.                  &
                 axi(1,2)-amrray_spheres_r(isph)).and.            &
                 (amrray_spheres_pos(2,isph).le.                  &
                 axi(2,2)+amrray_spheres_r(isph)).and.            &
                 (amrray_spheres_pos(3,isph).ge.                  &
                 axi(1,3)-amrray_spheres_r(isph)).and.            &
                 (amrray_spheres_pos(3,isph).le.                  &
                 axi(2,3)+amrray_spheres_r(isph))) then
               !
               ! Yes, this sphere may intersect this grid cell
               !
               if(amrray_spheres_sphidx(icell).eq.0) then
                  !
                  ! No previous sphere associated to this cell, so associate
                  ! this sphere to this cell
                  !
                  amrray_spheres_sphidx(icell) = isph
               else
                  !
                  ! Already a previous sphere associated. So we force the
                  ! code to check out all spheres for this cell. This is
                  ! slow, but hopefully only necessary in rare cases.
                  !
                  amrray_spheres_sphidx(icell) = -1
               endif
            endif
         enddo
      enddo
      !
   endif
   !
elseif(amr_coordsystem.lt.200) then
   !
   ! Spherical coordinates
   !
   ! Check if one or more spheres may overlap with the grid. Note that
   ! this is not a perfect criterion, and we might raise false alarm.
   ! But better this way around than not raise alarm when a sphere
   ! does intersect with the grid.
   !
   flag = .false.
   do isph=1,amrray_spheres_nr
      if(amrray_mirror_equator) then
         if(amrray_spheres_pos(3,isph).le.amrray_spheres_r(isph)) then
            write(stdo,*) 'ERROR in amrray module: if using mirror symmetry,'
            write(stdo,*) '      spheres (stars) must be entirely above the'
            write(stdo,*) '      equatorial plane.'
            stop
         endif
      endif
      call amr_xyz_to_rthphi(amrray_spheres_pos(1,isph), &
                             amrray_spheres_pos(2,isph), &
                             amrray_spheres_pos(3,isph), &
                             r,theta,phi)
      rmin     = r - amrray_spheres_r(isph)
      rmax     = r + amrray_spheres_r(isph)
      if(rmin.gt.0.d0) then
         thetamin = theta - amrray_spheres_r(isph)/rmin
         thetamax = theta + amrray_spheres_r(isph)/rmin
         if(thetamin.lt.0.d0) thetamin=0.d0
         if(thetamin.gt.pi) thetamax=pi
         sintheta = min(sin(thetamin),sin(thetamax))
         if(sintheta.gt.0.d0) then
            phimin   = phi - amrray_spheres_r(isph)/rmin
            phimax   = phi + amrray_spheres_r(isph)/rmin
         else
            phimin   = 0.d0
            phimax   = twopi
         endif
         centralsphere = .false.
      else
         rmin     = 0.d0
         thetamin = 0.d0
         thetamax = pi
         phimin   = 0.d0
         phimax   = twopi
         centralsphere = .true.
      endif
      !
      ! Check whether this sphere is at least partly overlapping 
      ! with the grid
      !
      if((rmax.ge.amr_grid_xi(1,1)).and.                  &
         (rmin.le.amr_grid_xi(amr_grid_nx+1,1)).and.      &
         (thetamax.ge.amr_grid_xi(1,2)).and.              &
         (thetamin.le.amr_grid_xi(amr_grid_ny+1,2)).and.  &
         (phimax.ge.amr_grid_xi(1,3)).and.                &
         (phimin.le.amr_grid_xi(amr_grid_nz+1,3))) then
         flag = .true.
      endif
      !
      ! Check whether this sphere is at least partly outside 
      ! the grid.
      !
      if((rmin.ge.amr_grid_xi(1,1)).and.                  &
         (rmax.le.amr_grid_xi(amr_grid_nx+1,1)).and.      &
         (thetamin.ge.amr_grid_xi(1,2)).and.              &
         (thetamax.le.amr_grid_xi(amr_grid_ny+1,2)).and.  &
         (phimin.ge.amr_grid_xi(1,3)).and.                &
         (phimax.le.amr_grid_xi(amr_grid_nz+1,3))) then
         amrray_spheres_outsidegrid(isph) = .false.
      else
         amrray_spheres_outsidegrid(isph) = .true.
      endif
      !
      ! If this sphere is a "central sphere", in the sense that
      ! it includes the (x,y,z)=(0,0,0) point, then we do not
      ! allow the sphere overlap with (or even touch) the 
      ! inner grid cells. The reason is a bit vague: it is simply
      ! because I am not sure if the algorithm manages it well
      ! if one has a perfectly central star with radius exactly
      ! equal to the inner radius of the grid. In that case
      ! the user must make the star a tiny bit smaller than 
      ! the inner grid radius.
      !
      if(centralsphere) then
         if(rmax.ge.oneminuss*amr_grid_xi(1,1)) then
            write(stdo,*) 'ERROR: Central star touches (of overlaps with)'
            write(stdo,*) '       the inner edge of the grid. If you want'
            write(stdo,*) '       to have the star radius equal to the inner'
            write(stdo,*) '       grid radius, then please make it a tiny bit'
            write(stdo,*) '       smaller (say by a factor of 0.999999d0).'
            write(stdo,*) '       But a real overlapping is not allowed for'
            write(stdo,*) '       stars that contain the coordinate center.'
            stop
         endif
      endif
      !
   enddo
   !
   ! If we have one or more spheres overlapping with the grid, then we must
   ! allocate the amrray_spheres_sphidx(:) array
   !
   if(flag) then
      !
      ! Allocate the array
      !
      allocate(amrray_spheres_sphidx(1:amr_nrleafs_max),STAT=ierr)
      if(ierr.ne.0) then
         write(stdo,*) 'ERROR in amrray: Could not allocate amrray_spheres_sphidx array'
         stop
      endif
      amrray_spheres_sphidx(:) = 0
      !
      ! Now check out cell-by-cell 
      !
      if(amr_tree_present) then
         if(.not.allocated(amr_theleafs)) stop 5673
         if(.not.allocated(amr_theleaf_index)) stop 5674
      endif
      do icell=1,amr_nrleafs
         !
         ! Get pointers to this cell and cell walls
         !
         if(amr_tree_present) then
            !
            ! If AMR tree available, then use that
            !
            !!index = amr_theleaf_index(icell)
            leaf  => amr_theleafs(icell)%link
            do iddr=1,3
               axi(1,iddr) = amr_finegrid_xi(leaf%ixyzf(iddr),iddr,leaf%level)
               axi(2,iddr) = amr_finegrid_xi(leaf%ixyzf(iddr)+1,iddr,leaf%level)
            enddo
         else
            !
            ! If regular grid, then use ix, iy, iz
            !
            call amr_regular_get_ixyz(icell,ix,iy,iz)
            axi(1,1) = amr_finegrid_xi(ix,1,0)
            axi(2,1) = amr_finegrid_xi(ix+1,1,0)
            axi(1,2) = amr_finegrid_xi(iy,2,0)
            axi(2,2) = amr_finegrid_xi(iy+1,2,0)
            axi(1,3) = amr_finegrid_xi(iz,3,0)
            axi(2,3) = amr_finegrid_xi(iz+1,3,0)
         endif
         !
         ! Now check out all the spheres
         !
         do isph=1,amrray_spheres_nr
            !
            ! Retrieve the rminmax, thetaminmax, phiminmax the stupid
            ! way: every time again. This can be done more efficiently,
            ! but since I do not really expect this to be often occurring,
            ! I don't want to waste my own time on this right now.
            !
            call amr_xyz_to_rthphi(amrray_spheres_pos(1,isph), &
                                   amrray_spheres_pos(2,isph), &
                                   amrray_spheres_pos(3,isph), &
                                   r,theta,phi)
            rmin     = r - amrray_spheres_r(isph)
            rmax     = r + amrray_spheres_r(isph)
            if(rmin.gt.0.d0) then
               thetamin = theta - amrray_spheres_r(isph)/rmin
               thetamax = theta + amrray_spheres_r(isph)/rmin
               if(thetamin.lt.0.d0) thetamin=0.d0
               if(thetamin.gt.pi) thetamax=pi
               sintheta = min(sin(thetamin),sin(thetamax))
               if(sintheta.gt.0.d0) then
                  phimin   = phi - amrray_spheres_r(isph)/rmin
                  phimax   = phi + amrray_spheres_r(isph)/rmin
               else
                  phimin   = 0.d0
                  phimax   = twopi
               endif
            else
               rmin     = 0.d0
               thetamin = 0.d0
               thetamax = pi
               phimin   = 0.d0
               phimax   = twopi
            endif
            !
            ! Now check if this sphere could be in this cell
            !
            if((rmax.ge.axi(1,1)).and.      &
               (rmin.le.axi(2,1)).and.      &
               (thetamax.ge.axi(1,2)).and.  &
               (thetamin.le.axi(2,2)).and.  &
               (phimax.ge.axi(1,3)).and.    &
               (phimin.le.axi(2,3))) then
               !
               ! Yes, this sphere may intersect this grid cell
               !
               if(amrray_spheres_sphidx(icell).eq.0) then
                  !
                  ! No previous sphere associated to this cell, so associate
                  ! this sphere to this cell
                  !
                  amrray_spheres_sphidx(icell) = isph
               else
                  !
                  ! Already a previous sphere associated. So we force the
                  ! code to check out all spheres for this cell. This is
                  ! slow, but hopefully only necessary in rare cases.
                  !
                  amrray_spheres_sphidx(icell) = -1
               endif
            endif
         enddo
      enddo
      !
   endif
else
   stop 
endif
!
! Check that none of the spheres overlaps with any of the other
! spheres...
!
do isph=1,amrray_spheres_nr-1
   do isph1=isph+1,amrray_spheres_nr
      dist = sqrt((amrray_spheres_pos(1,isph)-amrray_spheres_pos(1,isph1))**2 + &
                  (amrray_spheres_pos(2,isph)-amrray_spheres_pos(2,isph1))**2 + &
                  (amrray_spheres_pos(3,isph)-amrray_spheres_pos(3,isph1))**2 )
      if(dist.le.onepluss*(amrray_spheres_r(isph)+amrray_spheres_r(isph1))) then
         write(stdo,*) 'ERROR: Spheres/stars ',isph,' and ',isph1,' are'
         write(stdo,*) '       overlapping or touching each other.'
         stop
      endif
   enddo
enddo
!
end subroutine amrray_install_spheres


!--------------------------------------------------------------------------
!         Find leaf-cell containing x,y,z, starting from nodal cell
!
! This is just a helper routine for the amrray_find_next_location()
!--------------------------------------------------------------------------
subroutine amrray_find_subcell(cell,idir,ilr,x,y,z)
implicit none
type(amr_branch), pointer :: cell
integer idir,ilr,ix,iy,iz,iddr
doubleprecision x,y,z,dum,axc(1:3)
!
! Do loop until we arrive at smallest cell 
!
do while(.not.cell%leaf) 
   !
   ! Get the cell centers
   !
   do iddr=1,3
      axc(iddr) = amr_finegrid_xc(cell%ixyzf(iddr),iddr,cell%level)
   enddo
   !
   ! Check which direction to test, and get the indices of the children
   !
   select case(idir)
   case(0)
      write(stdo,*) 'INTERNAL ERROR: This mode not yet ready'
      stop 9103
   case(1)
      ix = ilr
      if(amr_ydim.eq.1) then
         if(y.lt.axc(2)) then
            iy = 1
         elseif(y.gt.axc(2)) then
            iy = 2
         else
            ! Bizarre incident: photon enters branch *exactly* in the middle
            ! Must do special things here.
            y  = y - tiny*(abs(y)+1.59283d-95)
            iy = 1
         endif
      else
         iy = 1
      endif
      if(amr_zdim.eq.1) then
         if(z.lt.axc(3)) then
            iz = 1
         elseif(z.gt.axc(3)) then
            iz = 2
         else
            ! Bizarre incident: photon enters branch *exactly* in the middle
            ! Must do special things here.
            z  = z - tiny*(abs(z)+1.25781d-95)
            iz = 1
         endif
      else
         iz = 1
      endif
   case(2)
      iy = ilr
      if(amr_xdim.eq.1) then
         if(x.lt.axc(1)) then
            ix = 1
         elseif(x.gt.axc(1)) then
            ix = 2
         else
            ! Bizarre incident: photon enters branch *exactly* in the middle
            ! Must do special things here.
            x  = x - tiny*(abs(x)+1.65934d-95)
            ix = 1
         endif
      else
         ix = 1
      endif
      if(amr_zdim.eq.1) then
         if(z.lt.axc(3)) then
            iz = 1
         elseif(z.gt.axc(3)) then
            iz = 2
         else
            ! Bizarre incident: photon enters branch *exactly* in the middle
            ! Must do special things here.
            z  = z - tiny*(abs(z)+1.27404d-95)
            iz = 1
         endif
      else
         iz = 1
      endif
   case(3)
      iz = ilr
      if(amr_xdim.eq.1) then
         if(x.lt.axc(1)) then
            ix = 1
         elseif(x.gt.axc(1)) then
            ix = 2
         else
            ! Bizarre incident: photon enters branch *exactly* in the middle
            ! Must do special things here.
            x  = x - tiny*(abs(x)+1.04713d-95)
            ix = 1
         endif
      else
         ix = 1
      endif
      if(amr_ydim.eq.1) then
         if(y.lt.axc(2)) then
            iy = 1
         elseif(y.gt.axc(2)) then
            iy = 2
         else
            ! Bizarre incident: photon enters branch *exactly* in the middle
            ! Must do special things here.
            y  = y - tiny*(abs(y)+1.75498d-95)
            iy = 1
         endif
      else
         iy = 1
      endif
   end select
   !
   ! Internal consistency check
   !
   if((ix.lt.1).or.(ix.gt.2).or.(iy.lt.1).or.(iy.gt.2).or. &
      (iz.lt.1).or.(iz.gt.2)) then
      write(stdo,*) idir,ilr
      write(stdo,*) x,y,z
      write(stdo,*) axc(1:3)
      write(stdo,*) ix,iy,iz
      stop 7263
   endif
   !
   ! Check existence of child
   !
   if(.not.associated(cell%child)) then
      write(stdo,*) 'INTERNAL ERROR: No child associated'
      stop 7202
   endif
   if(.not.associated(cell%child(ix,iy,iz)%link)) then
      write(stdo,*) 'INTERNAL ERROR: No child associated'
      stop 7203
   endif
   !
   ! Now take the child as current cell
   !
   cell => cell%child(ix,iy,iz)%link
   !
enddo
!
end subroutine amrray_find_subcell


!--------------------------------------------------------------------------
!                   Find next position (Cartesian case)
!
! This routine computes the next position within the cell OR the location
! on the cell interface including the pointer to the next cell.
!
! ARGUMENTS:
!   ray_dsend          If next crossing is farther than dsend, then stop
!                      at dsend (i.e. before next crossing).
!   ray_cart_x,_y,_z   Current position of photon
!   ray_cart_dirx/y/z  The direction unit vector. Must be normalized to 1.d0
!   ray_indexcurr      The cell index of the current cell 
!                      If >0, we are in cell with this index
!                      If =0, we are definitely outside of any cell
!                      If <0, we do not know: the subroutine must find out
!   arrived            Is .true. if ray_ds becomes ray_dsend or if the
!                      end of the grid is reached. In other words: if
!                      arrived=.true. then the current step is either zero
!                      in length or it is the last step through a cell before
!                      arriving either at the final point or escaping to
!                      infinity (if ray_dsend=1d99). 
! 
! RESULT:
!   ray_cart_x,_y,_z   New position of photon
!   ray_ds             The actual distance travelled (will be <= dsend)
!   ray_indexnext      Index of the cell in which the photon will next enter.
!                      If the photon stopped within the cell, then this will
!                      be the index of the current cell.
!   distmin (optional) The minimum distance to a cell wall 
!   levelnext (option) The level of the next cell
!
! RESULT STORED IN MODULE:
!   amrray_ispherehit  Only if allocated(amrray_spheres_r):
!                      If  0:  No sphere it hit,
!                      if -1:  Left sphere
!                      if >0:  Hit sphere nr ispherehit.
!--------------------------------------------------------------------------
subroutine amrray_find_next_location_cart(ray_dsend,                 &
                ray_cart_x,ray_cart_y,ray_cart_z,                    &
                ray_cart_dirx,ray_cart_diry,ray_cart_dirz,           &
                ray_indexcurr,ray_indexnext,ray_ds,arrived,          &
                distmin,levelnext)
implicit none
doubleprecision :: ray_cart_x,ray_cart_y,ray_cart_z
doubleprecision :: ray_cart_dirx,ray_cart_diry,ray_cart_dirz
doubleprecision :: ray_dsend,ray_ds,dsa,dsb
doubleprecision, optional :: distmin
integer :: ray_indexcurr,ray_indexnext
integer,optional :: levelnext

doubleprecision :: fact,xt,yt,zt,dummy,dsx,dsy,dsz,rr,ds_sphere,ds_sphere0
integer :: idx,idy,idz,ix,iy,iz,idir,ilr,isph,iddr
logical :: arrived
integer :: ihit
doubleprecision :: axi(1:2,1:3)
doubleprecision :: xyzdum
!
! Default
!
if(present(distmin)) distmin = 0.d0    ! Means: default = no fast subcell motion 
amrray_ispherehit = 0
ds_sphere0 = 1d99
amrray_icross = 0
!
! If cell index = 0, then go look for the cell, else simply link to the cell
!
if(ray_indexcurr.gt.0) then
   !
   ! We are in a cell
   !
   if(amr_tree_present) then
      !
      ! If the AMR tree is present, then link the amrray_cell
      !
      amrray_cell => amr_index_to_leaf(ray_indexcurr)%link
      if(.not.associated(amrray_cell)) then
         write(stdo,*) 'ERROR in amrray module: Cell ',ray_indexcurr,   &
              ' is not assigned. It does not exist.'
         stop 6023
      endif
   else
      !
      ! If we have a regular grid, then find the amrray_ix_curr, amrray_iy_curr, amrray_iz_curr
      !
      call amr_regular_get_ixyz(ray_indexcurr,amrray_ix_curr,amrray_iy_curr,amrray_iz_curr)
   endif
elseif(ray_indexcurr.eq.0) then
   !
   ! We are outside of a cell
   !
   nullify(amrray_cell)
   amrray_ix_curr = -1
   amrray_iy_curr = -1
   amrray_iz_curr = -1
else
   !
   ! We do not know, so we must find out ourselves
   !
   if(amr_tree_present) then
      !
      ! If the AMR tree is present, we use it
      !
      call amr_findcell(ray_cart_x,ray_cart_y,ray_cart_z,amrray_cell)
      if(associated(amrray_cell)) then
         ray_indexcurr = amrray_cell%leafindex
      else
         ray_indexcurr = 0
      endif
   else
      !
      ! In case of a regular grid without AMR tree, we use amrray_ix_curr, amrray_iy_curr, amrray_iz_curr
      !
      call amr_findbasecell(ray_cart_x,ray_cart_y,ray_cart_z,amrray_ix_curr,amrray_iy_curr,amrray_iz_curr)
      if(amrray_ix_curr.gt.0) then
         ray_indexcurr = amrray_ix_curr+(amrray_iy_curr-1)*amr_grid_nx+(amrray_iz_curr-1)*amr_grid_nx*amr_grid_ny
      else
         ray_indexcurr = 0
      endif
   endif
endif
!
! Check length of direction unit vector
!
if(amrray_selfcheck) then
   dummy = sqrt(ray_cart_dirx**2+ray_cart_diry**2+ray_cart_dirz**2) 
   if(abs(dummy-1.d0).gt.1d-10) then
      write(stdo,*) 'ERROR: Direction vector is not unit vector'
      stop 2913
   endif
endif
!
! If outside of grid, then do special stuff
!
if(ray_indexcurr.le.0) then 
   !
   ! Check that we are indeed outside of the grid
   !
   if(amrray_selfcheck) then
      if((ray_cart_x>=amr_grid_xi(1,1)).and.(ray_cart_x<=amr_grid_xi(amr_grid_nx+1,1)).and. &
           (ray_cart_y>=amr_grid_xi(1,2)).and.(ray_cart_y<=amr_grid_xi(amr_grid_ny+1,2)).and. &
           (ray_cart_z>=amr_grid_xi(1,3)).and.(ray_cart_z<=amr_grid_xi(amr_grid_nz+1,3))) then
         write(stdo,*) 'ERROR: Inconsistency in algorithm'
         write(stdo,*) ray_cart_x,ray_cart_y,ray_cart_z
         stop 1091
      endif
   endif
   !
   ! If we have spheres, in addition to the grid, then we must check
   ! out intersections with those.
   !
   if(allocated(amrray_spheres_r)) then
      do isph=1,amrray_spheres_nr
         if(amrray_spheres_outsidegrid(isph)) then
            !
            ! This sphere is (at least partly) outside of the grid. So
            ! we must test whether we intersect with it.
            !
            call amrray_find_crossing_with_sphere(          &
                 ray_cart_x,ray_cart_y,ray_cart_z,          &
                 ray_cart_dirx,ray_cart_diry,ray_cart_dirz, &
                 isph,ds_sphere0,ds_sphere,ihit,.false.)
            !
            ! If we have an intersection:
            !
            if(ihit.ne.0) then
               amrray_ispherehit = isph*ihit
               ds_sphere0 = ds_sphere
            endif
         endif
      enddo
      !
      ! Check if this is before ray_dsend
      !
      if(ds_sphere0.gt.ray_dsend) then
         amrray_ispherehit = 0
         ds_sphere0 = 1.d99
      endif
   endif
   !
   ! Check various scenarios
   !
   if((ray_cart_x<amr_grid_xi(1,1)).and.(ray_cart_dirx.gt.0.d0)) then
      !
      ! In X-direction, possible crossing with lower boundary
      !
      fact = (amr_grid_xi(1,1)-ray_cart_x)/ray_cart_dirx
      yt   = ray_cart_y + fact * ray_cart_diry
      zt   = ray_cart_z + fact * ray_cart_dirz
      if((yt.ge.amr_grid_xi(1,2)).and.&
         (yt.le.amr_grid_xi(amr_grid_ny+1,2)).and.&
         (zt.ge.amr_grid_xi(1,3)).and.&
         (zt.le.amr_grid_xi(amr_grid_nz+1,3))) then
         !
         ! Check if this is before hitting a/the sphere
         !
         if(amrray_ispherehit.ne.0) then
            if(fact.gt.ds_sphere0) then
               if(fact.le.ray_dsend) then
                  ray_cart_x = ray_cart_x + ds_sphere0 * ray_cart_dirx
                  ray_cart_y = ray_cart_y + ds_sphere0 * ray_cart_diry
                  ray_cart_z = ray_cart_z + ds_sphere0 * ray_cart_dirz
                  nullify(amrray_nextcell)
                  ray_indexnext = 0
                  if(present(levelnext)) levelnext = -1
                  return
               endif
            endif
         endif
         !
         ! Check if this is before ray_dsend. If so, do not enter grid
         !
         if(fact.gt.ray_dsend) then
            ray_cart_x = ray_cart_x + ray_dsend * ray_cart_dirx
            ray_cart_y = ray_cart_y + ray_dsend * ray_cart_diry
            ray_cart_z = ray_cart_z + ray_dsend * ray_cart_dirz
            nullify(amrray_nextcell)
            ray_indexnext = 0
            arrived = .true.
            if(present(levelnext)) levelnext = -1
            amrray_ispherehit = 0
            return
         endif
         !
         ! Entering grid
         ! ...coordinates of entry
         !
         ray_cart_x = amr_grid_xi(1,1)
         ray_cart_y = yt
         ray_cart_z = zt
         !
         ! ...finding which cell
         !
         call amrhunt(amr_grid_xi(:,2),amr_grid_ny+1,ray_cart_y,iy)
         call amrhunt(amr_grid_xi(:,3),amr_grid_nz+1,ray_cart_z,iz)
         ix = 1
         if(amrray_selfcheck) then
            if((iy.lt.1).or.(iy.gt.amr_grid_ny).or.&
                 (iz.lt.1).or.(iz.gt.amr_grid_nz)) then
               write(stdo,*) 'ERROR: Internal error'
               stop 6203
            endif
            !if(amr_tree_present) then
            !   if(.not.associated(amr_grid_branch(ix,iy,iz)%link)) then
            !      write(stdo,*) 'ERROR: Internal error'
            !      stop 6204
            !   endif
            !endif
         endif
         !
         ! ...Store this information
         !
         if(amr_tree_present) then
            !
            ! If the AMR tree is present, then use amrray_nextcell to 
            ! store the information
            !
            amrray_nextcell => amr_grid_branch(ix,iy,iz)%link
            !
            ! Now, if the new cell is in fact a nodal branch with children, then
            ! we must find this sub-cell
            !
            if(.not.amrray_nextcell%leaf) then
               call amrray_find_subcell(amrray_nextcell,1,1,&
                    ray_cart_x,ray_cart_y,ray_cart_z)
            endif
            !
            ! Finally, find index
            !
            ray_indexnext = amrray_nextcell%leafindex
            if(present(levelnext)) levelnext = amrray_nextcell%level
         else
            !
            ! If the grid is regular, then use amrray_ix_next, amrray_iy_next, amrray_iz_next
            ! to store the information
            !
            amrray_ix_next = ix
            amrray_iy_next = iy
            amrray_iz_next = iz
            !
            ! Finally, find index
            !
            ray_indexnext = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
            if(present(levelnext)) levelnext = 0
         endif
         !
         ! Finally, some flags
         !
         arrived = .false.
         amrray_ispherehit = 0
         amrray_icross = -1
         return
         !
      endif
   endif
   if((ray_cart_x>amr_grid_xi(amr_grid_nx+1,1)).and.(ray_cart_dirx.lt.0.d0)) then
      !
      ! In X-direction, possible crossing with upper boundary
      !
      fact = (amr_grid_xi(amr_grid_nx+1,1)-ray_cart_x)/ray_cart_dirx
      yt   = ray_cart_y + fact * ray_cart_diry
      zt   = ray_cart_z + fact * ray_cart_dirz
      if((yt.ge.amr_grid_xi(1,2)).and.&
         (yt.le.amr_grid_xi(amr_grid_ny+1,2)).and.&
         (zt.ge.amr_grid_xi(1,3)).and.&
         (zt.le.amr_grid_xi(amr_grid_nz+1,3))) then
         !
         ! Check if this is before hitting a/the sphere
         !
         if(amrray_ispherehit.ne.0) then
            if(fact.gt.ds_sphere0) then
               if(fact.le.ray_dsend) then
                  ray_cart_x = ray_cart_x + ds_sphere0 * ray_cart_dirx
                  ray_cart_y = ray_cart_y + ds_sphere0 * ray_cart_diry
                  ray_cart_z = ray_cart_z + ds_sphere0 * ray_cart_dirz
                  nullify(amrray_nextcell)
                  ray_indexnext = 0
                  if(present(levelnext)) levelnext = -1
                  return
               endif
            endif
         endif
         !
         ! Check if this is before ray_dsend. If so, do not enter grid
         !
         if(fact.gt.ray_dsend) then
            ray_cart_x = ray_cart_x + ray_dsend * ray_cart_dirx
            ray_cart_y = ray_cart_y + ray_dsend * ray_cart_diry
            ray_cart_z = ray_cart_z + ray_dsend * ray_cart_dirz
            nullify(amrray_nextcell)
            ray_indexnext = 0
            arrived = .true.
            if(present(levelnext)) levelnext = -1
            amrray_ispherehit = 0
            return
         endif
         !
         ! Entering grid
         ! ...coordinates of entry
         !
         ray_cart_x = amr_grid_xi(amr_grid_nx+1,1)
         ray_cart_y = yt
         ray_cart_z = zt
         !
         ! ...finding which cell
         !
         call amrhunt(amr_grid_xi(:,2),amr_grid_ny+1,ray_cart_y,iy)
         call amrhunt(amr_grid_xi(:,3),amr_grid_nz+1,ray_cart_z,iz)
         ix = amr_grid_nx
         if(amrray_selfcheck) then
            if((iy.lt.1).or.(iy.gt.amr_grid_ny).or.&
                 (iz.lt.1).or.(iz.gt.amr_grid_nz)) then
               write(stdo,*) 'ERROR: Internal error'
               stop 6203
            endif
            !if(amr_tree_present) then
            !   if(.not.associated(amr_grid_branch(ix,iy,iz)%link)) then
            !      write(stdo,*) 'ERROR: Internal error'
            !      stop 6204
            !   endif
            !endif
         endif
         !
         ! ...Store this information
         !
         if(amr_tree_present) then
            !
            ! If the AMR tree is present, then use amrray_nextcell to 
            ! store the information
            !
            amrray_nextcell => amr_grid_branch(ix,iy,iz)%link
            !
            ! Now, if the new cell is in fact a nodal branch with children, then
            ! we must find this sub-cell
            !
            if(.not.amrray_nextcell%leaf) then
               call amrray_find_subcell(amrray_nextcell,1,2,&
                    ray_cart_x,ray_cart_y,ray_cart_z)
            endif
            !
            ! Finally, find index
            !
            ray_indexnext = amrray_nextcell%leafindex
            if(present(levelnext)) levelnext = amrray_nextcell%level
         else
            !
            ! If the grid is regular, then use amrray_ix_next, amrray_iy_next, amrray_iz_next
            ! to store the information
            !
            amrray_ix_next = ix
            amrray_iy_next = iy
            amrray_iz_next = iz
            !
            ! Finally, find index
            !
            ray_indexnext = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
            if(present(levelnext)) levelnext = 0
         endif
         !
         ! Finally, some flags
         !
         arrived = .false.
         amrray_ispherehit = 0
         amrray_icross = -2
         return
         !
      endif
   endif
   if((ray_cart_y<amr_grid_xi(1,2)).and.(ray_cart_diry.gt.0.d0)) then
      !
      ! In Y-direction, possible crossing with lower boundary
      !
      fact = (amr_grid_xi(1,2)-ray_cart_y)/ray_cart_diry
      xt   = ray_cart_x + fact * ray_cart_dirx
      zt   = ray_cart_z + fact * ray_cart_dirz
      if((xt.ge.amr_grid_xi(1,1)).and.&
         (xt.le.amr_grid_xi(amr_grid_nx+1,1)).and.&
         (zt.ge.amr_grid_xi(1,3)).and.&
         (zt.le.amr_grid_xi(amr_grid_nz+1,3))) then
         !
         ! Check if this is before hitting a/the sphere
         !
         if(amrray_ispherehit.ne.0) then
            if(fact.gt.ds_sphere0) then
               if(fact.le.ray_dsend) then
                  ray_cart_x = ray_cart_x + ds_sphere0 * ray_cart_dirx
                  ray_cart_y = ray_cart_y + ds_sphere0 * ray_cart_diry
                  ray_cart_z = ray_cart_z + ds_sphere0 * ray_cart_dirz
                  nullify(amrray_nextcell)
                  ray_indexnext = 0
                  if(present(levelnext)) levelnext = -1
                  return
               endif
            endif
         endif
         !
         ! Check if this is before ray_dsend. If so, do not enter grid
         !
         if(fact.gt.ray_dsend) then
            ray_cart_x = ray_cart_x + ray_dsend * ray_cart_dirx
            ray_cart_y = ray_cart_y + ray_dsend * ray_cart_diry
            ray_cart_z = ray_cart_z + ray_dsend * ray_cart_dirz
            nullify(amrray_nextcell)
            ray_indexnext = 0
            arrived = .true.
            if(present(levelnext)) levelnext = -1
            amrray_ispherehit = 0
            return
         endif
         !
         ! Entering grid
         ! ...coordinates of entry
         !
         ray_cart_x = xt
         ray_cart_y = amr_grid_xi(1,2)
         ray_cart_z = zt
         !
         ! ...finding which cell
         !
         call amrhunt(amr_grid_xi(:,1),amr_grid_nx+1,ray_cart_x,ix)
         call amrhunt(amr_grid_xi(:,3),amr_grid_nz+1,ray_cart_z,iz)
         iy = 1
         if(amrray_selfcheck) then
            if((ix.lt.1).or.(ix.gt.amr_grid_nx).or.&
                 (iz.lt.1).or.(iz.gt.amr_grid_nz)) then
               write(stdo,*) 'ERROR: Internal error'
               stop 6203
            endif
            !if(amr_tree_present) then
            !   if(.not.associated(amr_grid_branch(ix,iy,iz)%link)) then
            !      write(stdo,*) 'ERROR: Internal error'
            !      stop 6204
            !   endif
            !endif
         endif
         !
         ! ...Store this information
         !
         if(amr_tree_present) then
            !
            ! If the AMR tree is present, then use amrray_nextcell to 
            ! store the information
            !
            amrray_nextcell => amr_grid_branch(ix,iy,iz)%link
            !
            ! Now, if the new cell is in fact a nodal branch with children, then
            ! we must find this sub-cell
            !
            if(.not.amrray_nextcell%leaf) then
               call amrray_find_subcell(amrray_nextcell,2,1,&
                    ray_cart_x,ray_cart_y,ray_cart_z)
            endif
            !
            ! Finally, find index
            !
            ray_indexnext = amrray_nextcell%leafindex
            if(present(levelnext)) levelnext = amrray_nextcell%level
         else
            !
            ! If the grid is regular, then use amrray_ix_next, amrray_iy_next, amrray_iz_next
            ! to store the information
            !
            amrray_ix_next = ix
            amrray_iy_next = iy
            amrray_iz_next = iz
            !
            ! Finally, find index
            !
            ray_indexnext = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
            if(present(levelnext)) levelnext = 0
         endif
         !
         ! Finally, some flags
         !
         arrived = .false.
         amrray_ispherehit = 0
         amrray_icross = -3
         return
         !
      endif
   endif
   if((ray_cart_y>amr_grid_xi(amr_grid_ny+1,2)).and.(ray_cart_diry.lt.0.d0)) then
      !
      ! In Y-direction, possible crossing with upper boundary
      !
      fact = (amr_grid_xi(amr_grid_ny+1,2)-ray_cart_y)/ray_cart_diry
      xt   = ray_cart_x + fact * ray_cart_dirx
      zt   = ray_cart_z + fact * ray_cart_dirz
      if((xt.ge.amr_grid_xi(1,1)).and.&
         (xt.le.amr_grid_xi(amr_grid_nx+1,1)).and.&
         (zt.ge.amr_grid_xi(1,3)).and.&
         (zt.le.amr_grid_xi(amr_grid_nz+1,3))) then
         !
         ! Check if this is before hitting a/the sphere
         !
         if(amrray_ispherehit.ne.0) then
            if(fact.gt.ds_sphere0) then
               if(fact.le.ray_dsend) then
                  ray_cart_x = ray_cart_x + ds_sphere0 * ray_cart_dirx
                  ray_cart_y = ray_cart_y + ds_sphere0 * ray_cart_diry
                  ray_cart_z = ray_cart_z + ds_sphere0 * ray_cart_dirz
                  nullify(amrray_nextcell)
                  ray_indexnext = 0
                  if(present(levelnext)) levelnext = -1
                  return
               endif
            endif
         endif
         !
         ! Check if this is before ray_dsend. If so, do not enter grid
         !
         if(fact.gt.ray_dsend) then
            ray_cart_x = ray_cart_x + ray_dsend * ray_cart_dirx
            ray_cart_y = ray_cart_y + ray_dsend * ray_cart_diry
            ray_cart_z = ray_cart_z + ray_dsend * ray_cart_dirz
            nullify(amrray_nextcell)
            ray_indexnext = 0
            arrived = .true.
            if(present(levelnext)) levelnext = -1
            amrray_ispherehit = 0
            return
         endif
         !
         ! Entering grid
         ! ...coordinates of entry
         !
         ray_cart_x = xt
         ray_cart_y = amr_grid_xi(amr_grid_ny+1,2)
         ray_cart_z = zt
         !
         ! ...finding which cell
         !
         call amrhunt(amr_grid_xi(:,1),amr_grid_nx+1,ray_cart_x,ix)
         call amrhunt(amr_grid_xi(:,3),amr_grid_nz+1,ray_cart_z,iz)
         iy = amr_grid_ny
         if(amrray_selfcheck) then
            if((ix.lt.1).or.(ix.gt.amr_grid_nx).or.&
                 (iz.lt.1).or.(iz.gt.amr_grid_nz)) then
               write(stdo,*) 'ERROR: Internal error'
               stop 6203
            endif
            !if(amr_tree_present) then
            !   if(.not.associated(amr_grid_branch(ix,iy,iz)%link)) then
            !      write(stdo,*) 'ERROR: Internal error'
            !      stop 6204
            !   endif
            !endif
         endif
         !
         ! ...Store this information
         !
         if(amr_tree_present) then
            !
            ! If the AMR tree is present, then use amrray_nextcell to 
            ! store the information
            !
            amrray_nextcell => amr_grid_branch(ix,iy,iz)%link
            !
            ! Now, if the new cell is in fact a nodal branch with children, then
            ! we must find this sub-cell
            !
            if(.not.amrray_nextcell%leaf) then
               call amrray_find_subcell(amrray_nextcell,2,2,&
                    ray_cart_x,ray_cart_y,ray_cart_z)
            endif
            !
            ! Finally, find index
            !
            ray_indexnext = amrray_nextcell%leafindex
            if(present(levelnext)) levelnext = amrray_nextcell%level
         else
            !
            ! If the grid is regular, then use amrray_ix_next, amrray_iy_next, amrray_iz_next
            ! to store the information
            !
            amrray_ix_next = ix
            amrray_iy_next = iy
            amrray_iz_next = iz
            !
            ! Finally, find index
            !
            ray_indexnext = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
            if(present(levelnext)) levelnext = 0
         endif
         !
         ! Finally, some flags
         !
         arrived = .false.
         amrray_ispherehit = 0
         amrray_icross = -4
         return
         !
      endif
   endif
   if((ray_cart_z<amr_grid_xi(1,3)).and.(ray_cart_dirz.gt.0.d0)) then
      !
      ! In Z-direction, possible crossing with lower boundary
      !
      fact = (amr_grid_xi(1,3)-ray_cart_z)/ray_cart_dirz
      xt   = ray_cart_x + fact * ray_cart_dirx
      yt   = ray_cart_y + fact * ray_cart_diry
      if((yt.ge.amr_grid_xi(1,2)).and.&
         (yt.le.amr_grid_xi(amr_grid_ny+1,2)).and.&
         (xt.ge.amr_grid_xi(1,1)).and.&
         (xt.le.amr_grid_xi(amr_grid_nx+1,1))) then
         !
         ! Check if this is before hitting a/the sphere
         !
         if(amrray_ispherehit.ne.0) then
            if(fact.gt.ds_sphere0) then
               if(fact.le.ray_dsend) then
                  ray_cart_x = ray_cart_x + ds_sphere0 * ray_cart_dirx
                  ray_cart_y = ray_cart_y + ds_sphere0 * ray_cart_diry
                  ray_cart_z = ray_cart_z + ds_sphere0 * ray_cart_dirz
                  nullify(amrray_nextcell)
                  ray_indexnext = 0
                  if(present(levelnext)) levelnext = -1
                  return
               endif
            endif
         endif
         !
         ! Check if this is before ray_dsend. If so, do not enter grid
         !
         if(fact.gt.ray_dsend) then
            ray_cart_x = ray_cart_x + ray_dsend * ray_cart_dirx
            ray_cart_y = ray_cart_y + ray_dsend * ray_cart_diry
            ray_cart_z = ray_cart_z + ray_dsend * ray_cart_dirz
            nullify(amrray_nextcell)
            ray_indexnext = 0
            arrived = .true.
            if(present(levelnext)) levelnext = -1
            amrray_ispherehit = 0
            return
         endif
         !
         ! Entering grid
         ! ...coordinates of entry
         !
         ray_cart_x = xt
         ray_cart_y = yt
         ray_cart_z = amr_grid_xi(1,3)
         !
         ! ...finding which cell
         !
         call amrhunt(amr_grid_xi(:,1),amr_grid_nx+1,ray_cart_x,ix)
         call amrhunt(amr_grid_xi(:,2),amr_grid_ny+1,ray_cart_y,iy)
         iz = 1
         if(amrray_selfcheck) then
            if((iy.lt.1).or.(iy.gt.amr_grid_ny).or.&
                 (ix.lt.1).or.(ix.gt.amr_grid_nx)) then
               write(stdo,*) 'ERROR: Internal error'
               stop 6203
            endif
            !if(amr_tree_present) then
            !   if(.not.associated(amr_grid_branch(ix,iy,iz)%link)) then
            !      write(stdo,*) 'ERROR: Internal error'
            !      stop 6204
            !   endif
            !endif
         endif
         !
         ! ...Store this information
         !
         if(amr_tree_present) then
            !
            ! If the AMR tree is present, then use amrray_nextcell to 
            ! store the information
            !
            amrray_nextcell => amr_grid_branch(ix,iy,iz)%link
            !
            ! Now, if the new cell is in fact a nodal branch with children, then
            ! we must find this sub-cell
            !
            if(.not.amrray_nextcell%leaf) then
               call amrray_find_subcell(amrray_nextcell,3,1,&
                    ray_cart_x,ray_cart_y,ray_cart_z)
            endif
            !
            ! Finally, find index
            !
            ray_indexnext = amrray_nextcell%leafindex
            if(present(levelnext)) levelnext = amrray_nextcell%level
         else
            !
            ! If the grid is regular, then use amrray_ix_next, amrray_iy_next, amrray_iz_next
            ! to store the information
            !
            amrray_ix_next = ix
            amrray_iy_next = iy
            amrray_iz_next = iz
            !
            ! Finally, find index
            !
            ray_indexnext = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
            if(present(levelnext)) levelnext = 0
         endif
         !
         ! Finally, some flags
         !
         arrived = .false.
         amrray_ispherehit = 0
         amrray_icross = -5
         return
         !
      endif
   endif
   if((ray_cart_z>amr_grid_xi(amr_grid_nz+1,3)).and.(ray_cart_dirz.lt.0.d0)) then
      !
      ! In Z-direction, possible crossing with upper boundary
      !
      fact = (amr_grid_xi(amr_grid_nz+1,3)-ray_cart_z)/ray_cart_dirz
      xt   = ray_cart_x + fact * ray_cart_dirx
      yt   = ray_cart_y + fact * ray_cart_diry
      if((yt.ge.amr_grid_xi(1,2)).and.&
         (yt.le.amr_grid_xi(amr_grid_ny+1,2)).and.&
         (xt.ge.amr_grid_xi(1,1)).and.&
         (xt.le.amr_grid_xi(amr_grid_nx+1,1))) then
         !
         ! Check if this is before hitting a/the sphere
         !
         if(amrray_ispherehit.ne.0) then
            if(fact.gt.ds_sphere0) then
               if(fact.le.ray_dsend) then
                  ray_cart_x = ray_cart_x + ds_sphere0 * ray_cart_dirx
                  ray_cart_y = ray_cart_y + ds_sphere0 * ray_cart_diry
                  ray_cart_z = ray_cart_z + ds_sphere0 * ray_cart_dirz
                  nullify(amrray_nextcell)
                  ray_indexnext = 0
                  if(present(levelnext)) levelnext = -1
                  return
               endif
            endif
         endif
         !
         ! Check if this is before ray_dsend. If so, do not enter grid
         !
         if(fact.gt.ray_dsend) then
            ray_cart_x = ray_cart_x + ray_dsend * ray_cart_dirx
            ray_cart_y = ray_cart_y + ray_dsend * ray_cart_diry
            ray_cart_z = ray_cart_z + ray_dsend * ray_cart_dirz
            nullify(amrray_nextcell)
            ray_indexnext = 0
            arrived = .true.
            if(present(levelnext)) levelnext = -1
            amrray_ispherehit = 0
            return
         endif
         !
         ! Entering grid
         ! ...coordinates of entry
         !
         ray_cart_x = xt
         ray_cart_y = yt
         ray_cart_z = amr_grid_xi(amr_grid_nz+1,3)
         !
         ! ...finding which cell
         !
         call amrhunt(amr_grid_xi(:,1),amr_grid_nx+1,ray_cart_x,ix)
         call amrhunt(amr_grid_xi(:,2),amr_grid_ny+1,ray_cart_y,iy)
         iz = amr_grid_nz
         if(amrray_selfcheck) then
            if((iy.lt.1).or.(iy.gt.amr_grid_ny).or.&
                 (ix.lt.1).or.(ix.gt.amr_grid_nx)) then
               write(stdo,*) 'ERROR: Internal error'
               stop 6203
            endif
            !if(amr_tree_present) then
            !   if(.not.associated(amr_grid_branch(ix,iy,iz)%link)) then
            !      write(stdo,*) 'ERROR: Internal error'
            !      stop 6204
            !   endif
            !endif
         endif
         !
         ! ...Store this information
         !
         if(amr_tree_present) then
            !
            ! If the AMR tree is present, then use amrray_nextcell to 
            ! store the information
            !
            amrray_nextcell => amr_grid_branch(ix,iy,iz)%link
            !
            ! Now, if the new cell is in fact a nodal branch with children, then
            ! we must find this sub-cell
            !
            if(.not.amrray_nextcell%leaf) then
               call amrray_find_subcell(amrray_nextcell,3,2,&
                    ray_cart_x,ray_cart_y,ray_cart_z)
            endif
            !
            ! Finally, find index
            !
            ray_indexnext = amrray_nextcell%leafindex
            if(present(levelnext)) levelnext = amrray_nextcell%level
         else
            !
            ! If the grid is regular, then use amrray_ix_next, amrray_iy_next, amrray_iz_next
            ! to store the information
            !
            amrray_ix_next = ix
            amrray_iy_next = iy
            amrray_iz_next = iz
            !
            ! Finally, find index
            !
            ray_indexnext = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
            if(present(levelnext)) levelnext = 0
         endif
         !
         ! Finally, some flags
         !
         arrived = .false.
         amrray_ispherehit = 0
         amrray_icross = -6
         return
         !
      endif
   endif
   !
   ! None of the scenarios worked, so the ray misses the grid completely
   ! but may have hit a sphere, so that is why we keep amrray_ispherehit 
   ! to the value it had.
   !
   nullify(amrray_nextcell)
   ray_indexnext = 0
   amrray_ix_next = -1
   amrray_iy_next = -1
   amrray_iz_next = -1
   arrived = .true.
   if(present(levelnext)) levelnext = -1
   return
endif
!
! We start INSIDE of the grid
!
! First plant an entry point for the case that we end up exactly
! at a cell ribbon or corner, and therefore have to retry...
!
100 continue
!
! Get the cell walls
!
if(amr_tree_present) then
   !
   ! If the AMR tree is present, then use that
   !
   do iddr=1,3
      axi(1,iddr) = amr_finegrid_xi(amrray_cell%ixyzf(iddr),iddr,amrray_cell%level)
      axi(2,iddr) = amr_finegrid_xi(amrray_cell%ixyzf(iddr)+1,iddr,amrray_cell%level)
   enddo
else
   !
   ! If we have a regular grid, then use amrray_ix_curr, amrray_iy_curr, amrray_iz_curr, and
   ! use level=0
   !
   axi(1,1) = amr_finegrid_xi(amrray_ix_curr,1,0)
   axi(2,1) = amr_finegrid_xi(amrray_ix_curr+1,1,0)
   axi(1,2) = amr_finegrid_xi(amrray_iy_curr,2,0)
   axi(2,2) = amr_finegrid_xi(amrray_iy_curr+1,2,0)
   axi(1,3) = amr_finegrid_xi(amrray_iz_curr,3,0)
   axi(2,3) = amr_finegrid_xi(amrray_iz_curr+1,3,0)
endif
!
! Self-check
!
if(amrray_selfcheck) then
   if((ray_cart_x>axi(2,1)).or.(ray_cart_x<axi(1,1)).or. &
      (ray_cart_y>axi(2,2)).or.(ray_cart_y<axi(1,2)).or. &
      (ray_cart_z>axi(2,3)).or.(ray_cart_z<axi(1,3))) then
      write(stdo,*) 'ERROR: Photon outside of cell'
      write(stdo,*) ray_cart_x,ray_cart_y,ray_cart_z
      if(amr_tree_present) write(stdo,*) associated(amrray_cell),amrray_cell%id
      write(stdo,*) axi
      stop 1031
   endif
endif
!
! Find distances to all three coordinates
!
if(present(distmin)) then
   distmin = 1.d99
   dsa = axi(2,1)-ray_cart_x
   dsb = ray_cart_x-axi(1,1)
   distmin = min(distmin,dsa,dsb)
   if(ray_cart_dirx.gt.0.d0) then
      dsx = dsa/ray_cart_dirx
      idx = 2
   elseif(ray_cart_dirx.lt.0.d0) then
      dsx = -dsb/ray_cart_dirx
      idx = 1
   else
      dsx = 1d99
      idx = 1
   endif
   dsa = axi(2,2)-ray_cart_y
   dsb = ray_cart_y-axi(1,2)
   distmin = min(distmin,dsa,dsb)
   if(ray_cart_diry.gt.0.d0) then
      dsy = dsa/ray_cart_diry
      idy = 2
   elseif(ray_cart_diry.lt.0.d0) then
      dsy = -dsb/ray_cart_diry
      idy = 1
   else
      dsy = 1d99
      idy = 1
   endif
   dsa = axi(2,3)-ray_cart_z
   dsb = ray_cart_z-axi(1,3)
   distmin = min(distmin,dsa,dsb)
   if(ray_cart_dirz.gt.0.d0) then
      dsz = dsa/ray_cart_dirz
      idz = 2
   elseif(ray_cart_dirz.lt.0.d0) then
      dsz = -dsb/ray_cart_dirz
      idz = 1
   else
      dsz = 1d99
      idz = 1
   endif
else
   !
   ! If distmin does not need to be evaluated, then things can be done
   ! quicker
   !
   if(ray_cart_dirx.gt.0.d0) then
      dsx = (axi(2,1)-ray_cart_x)/ray_cart_dirx
      idx = 2
   elseif(ray_cart_dirx.lt.0.d0) then
      dsx = (axi(1,1)-ray_cart_x)/ray_cart_dirx
      idx = 1
   else
      dsx = 1d99
      idx = 1
   endif
   if(ray_cart_diry.gt.0.d0) then
      dsy = (axi(2,2)-ray_cart_y)/ray_cart_diry
      idy = 2
   elseif(ray_cart_diry.lt.0.d0) then
      dsy = (axi(1,2)-ray_cart_y)/ray_cart_diry
      idy = 1
   else
      dsy = 1d99
      idy = 1
   endif
   if(ray_cart_dirz.gt.0.d0) then
      dsz = (axi(2,3)-ray_cart_z)/ray_cart_dirz
      idz = 2
   elseif(ray_cart_dirz.lt.0.d0) then
      dsz = (axi(1,3)-ray_cart_z)/ray_cart_dirz
      idz = 1
   else
      dsz = 1d99
      idz = 1
   endif
endif
!
! Self-check
!
if(amrray_selfcheck) then
   if((dsx.lt.0.d0).or.(dsy.lt.0.d0).or.(dsz.lt.0.d0)) then
      write(stdo,*) 'ERROR: Negative dsz, dsy or dsz'
      stop 1040
   endif
endif
!
! Check bizarre incidents of photons exactly in corner or ribbon
!
if(dsx.eq.dsy) then
   if((dsx.ne.1d99).and.(dsy.ne.1d99).and.(dsz.ge.dsx)) then
      ! Bizarre incident: photon EXACTLY in cell corner
      ! Have to shift the current point a tiny bit, and redo things
      if(amr_tree_present) then
         xyzdum = amr_finegrid_xc(amrray_cell%ixyzf(1),1,amrray_cell%level)
      else
         xyzdum = amr_finegrid_xc(amrray_ix_curr,1,0)
      endif
      if(ray_cart_x.lt.xyzdum) then
         ray_cart_x = ray_cart_x + tiny*(abs(ray_cart_x)+2.03105d-95)
      else
         ray_cart_x = ray_cart_x - tiny*(abs(ray_cart_x)+2.03105d-95)
      endif
      goto 100
   endif
endif
if(dsx.eq.dsz) then
   if((dsx.ne.1d99).and.(dsz.ne.1d99).and.(dsy.ge.dsz)) then
      ! Bizarre incident: photon EXACTLY in cell corner
      ! Have to shift the current point a tiny bit, and redo things
      if(amr_tree_present) then
         xyzdum = amr_finegrid_xc(amrray_cell%ixyzf(1),1,amrray_cell%level)
      else
         xyzdum = amr_finegrid_xc(amrray_ix_curr,1,0)
      endif
      if(ray_cart_x.lt.xyzdum) then
         ray_cart_x = ray_cart_x + tiny*(abs(ray_cart_x)+2.25671d-95)
      else
         ray_cart_x = ray_cart_x - tiny*(abs(ray_cart_x)+2.25671d-95)
      endif
      goto 100
   endif
endif
if(dsy.eq.dsz) then
   if((dsy.ne.1d99).and.(dsz.ne.1d99).and.(dsx.ge.dsy)) then
      ! Bizarre incident: photon EXACTLY in cell corner
      ! Have to shift the current point a tiny bit, and redo things
      if(amr_tree_present) then
         xyzdum = amr_finegrid_xc(amrray_cell%ixyzf(2),2,amrray_cell%level)
      else
         xyzdum = amr_finegrid_xc(amrray_iy_curr,2,0)
      endif
      if(ray_cart_y.lt.xyzdum) then
         ray_cart_y = ray_cart_y + tiny*(abs(ray_cart_y)+2.55931d-95)
      else
         ray_cart_y = ray_cart_y - tiny*(abs(ray_cart_y)+2.55931d-95)
      endif
      goto 100
   endif
endif
!
! Find the closest wall
!
if(dsx.lt.dsy) then
   if(dsx.lt.dsz) then
      idir  = 1
      ray_ds = dsx
   else
      idir  = 3
      ray_ds = dsz
   endif
else
   if(dsy.lt.dsz) then
      idir  = 2
      ray_ds = dsy
   else
      idir  = 3
      ray_ds = dsz
   endif
endif
!
! Now check for the spheres
!
if(allocated(amrray_spheres_sphidx)) then
   if(ray_indexcurr.le.0) stop 7390
   if(amrray_spheres_sphidx(ray_indexcurr).gt.0) then
      !
      ! There is a single sphere inside or partly inside this cell
      !
      isph = amrray_spheres_sphidx(ray_indexcurr)
      if(isph.gt.amrray_spheres_nr) stop 7391
      !
      ! Now check the crossing with this sphere
      !
      call amrray_find_crossing_with_sphere(                &
                 ray_cart_x,ray_cart_y,ray_cart_z,          &
                 ray_cart_dirx,ray_cart_diry,ray_cart_dirz, &
                 isph,ray_ds,ds_sphere,ihit,.false.)
      !
      ! If we have an intersection:
      !
      if(ihit.ne.0) then
         amrray_ispherehit = isph*ihit
         ray_ds = ds_sphere
         idir = 0
      endif
      !
   elseif(amrray_spheres_sphidx(ray_indexcurr).lt.0) then
      !
      ! There is more than one sphere in this cell
      ! We must check all spheres (sorry, no cleverer method yet)
      !
      ! Do a loop over spheres
      !
      do isph=1,amrray_spheres_nr
         !
         ! Now check the crossing with this sphere
         !
         call amrray_find_crossing_with_sphere(             &
                 ray_cart_x,ray_cart_y,ray_cart_z,          &
                 ray_cart_dirx,ray_cart_diry,ray_cart_dirz, &
                 isph,ray_ds,ds_sphere,ihit,.false.)
         !
         ! If we have an intersection:
         !
         if(ihit.ne.0) then
            amrray_ispherehit = isph*ihit
            ray_ds = ds_sphere
            idir = 0
         endif
      enddo
      !
   endif
endif
!
! Check the ds. If direction points outward from cell at starting
! point, the photon will get stuck. This is the check for this.
!
!if(ray_ds.eq.0.d0) then
!   write(stdo,*) 'ERROR: Photon got stuck on cell boundary'
!   write(stdo,*) ray_cart_x,ray_cart_y,ray_cart_z
!   write(stdo,*) ray_cart_dirx,ray_cart_diry,ray_cart_dirz
!   write(stdo,*) dsx,dsy,dsz
!   stop 2719
!endif
!
! If next end-point is closer than wall:
!
if(ray_ds.gt.ray_dsend) then 
   ray_ds = ray_dsend
   idir  = 0
   arrived = .true.
endif
!
! Update position, either in middle of cell (idir=0) or on cell wall
!
select case(idir) 
   case(0)
      if(amr_tree_present) then
         !
         ! The AMR tree is used, so connect amrray_nextcell pointer.
         ! If we get outside of the domain, this link is automatically
         ! null.
         !
         amrray_nextcell => amrray_cell
      else
         !
         ! A regular grid is used, so set the ixyz_next indices
         !
         amrray_ix_next = amrray_ix_curr
         amrray_iy_next = amrray_iy_curr
         amrray_iz_next = amrray_iz_curr
      endif
      !
      ! Update position and set flags
      !
      ray_cart_x = ray_cart_x + ray_ds*ray_cart_dirx
      ray_cart_y = ray_cart_y + ray_ds*ray_cart_diry
      ray_cart_z = ray_cart_z + ray_ds*ray_cart_dirz
      ilr = 0
      amrray_icross = 0
   case(1)
      !
      ! X-direction
      !
      if(amr_tree_present) then
         !
         ! The AMR tree is used, so connect amrray_nextcell pointer.
         ! If we get outside of the domain, this link is automatically
         ! null.
         !
         amrray_nextcell => amrray_cell%neighbor(idx,1)%link
         !
         ! Update flags
         !
         if(ray_cart_dirx.gt.0) then
            ilr = 1
            amrray_icross = 2
         else
            ilr = 2
            amrray_icross = 1
         endif
      else
         !
         ! A regular grid is used, so set the ixyz_next indices
         !
         ! Depending on direction...
         !
         if(ray_cart_dirx.gt.0) then
            !
            ! Update flags
            !
            ilr = 1
            amrray_icross = 2
            !
            ! Update position
            !
            amrray_ix_next = amrray_ix_curr + 1
            amrray_iy_next = amrray_iy_curr
            amrray_iz_next = amrray_iz_curr
            !
            ! Check if out of bounds
            !
            if(amrray_ix_next.gt.amr_grid_nx) then
               amrray_ix_next = -1
               amrray_iy_next = -1
               amrray_iz_next = -1
            endif
         else
            !
            ! Update flags
            !
            ilr = 2
            amrray_icross = 1
            !
            ! Update position
            !
            amrray_ix_next = amrray_ix_curr - 1
            amrray_iy_next = amrray_iy_curr
            amrray_iz_next = amrray_iz_curr
            !
            ! Check if out of bounds
            !
            if(amrray_ix_next.lt.1) then
               amrray_ix_next = -1
               amrray_iy_next = -1
               amrray_iz_next = -1
            endif
         endif
      endif
      !
      ! Update position 
      !
      ray_cart_x = axi(idx,1)
      ray_cart_y = ray_cart_y + ray_ds*ray_cart_diry
      ray_cart_z = ray_cart_z + ray_ds*ray_cart_dirz
   case(2)
      !
      ! Y-direction
      !
      if(amr_tree_present) then
         !
         ! The AMR tree is used, so connect amrray_nextcell pointer.
         ! If we get outside of the domain, this link is automatically
         ! null.
         !
         amrray_nextcell => amrray_cell%neighbor(idy,2)%link
         !
         ! Update flags
         !
         if(ray_cart_diry.gt.0) then
            ilr = 1
            amrray_icross = 4
         else
            ilr = 2
            amrray_icross = 3
         endif
      else
         !
         ! A regular grid is used, so set the ixyz_next indices
         !
         ! Depending on direction...
         !
         if(ray_cart_diry.gt.0) then
            !
            ! Update flags
            !
            ilr = 1
            amrray_icross = 4
            !
            ! Update position
            !
            amrray_ix_next = amrray_ix_curr
            amrray_iy_next = amrray_iy_curr + 1
            amrray_iz_next = amrray_iz_curr
            !
            ! Check if out of bounds
            !
            if(amrray_iy_next.gt.amr_grid_ny) then
               amrray_ix_next = -1
               amrray_iy_next = -1
               amrray_iz_next = -1
            endif
         else
            !
            ! Update flags
            !
            ilr = 2
            amrray_icross = 3
            !
            ! Update position
            !
            amrray_ix_next = amrray_ix_curr
            amrray_iy_next = amrray_iy_curr - 1
            amrray_iz_next = amrray_iz_curr
            !
            ! Check if out of bounds
            !
            if(amrray_iy_next.lt.1) then
               amrray_ix_next = -1
               amrray_iy_next = -1
               amrray_iz_next = -1
            endif
         endif
      endif
      !
      ! Update position 
      !
      ray_cart_x = ray_cart_x + ray_ds*ray_cart_dirx
      ray_cart_y = axi(idy,2)
      ray_cart_z = ray_cart_z + ray_ds*ray_cart_dirz
   case(3)
      !
      ! Z-direction
      !
      if(amr_tree_present) then
         !
         ! The AMR tree is used, so connect amrray_nextcell pointer.
         ! If we get outside of the domain, this link is automatically
         ! null.
         !
         amrray_nextcell => amrray_cell%neighbor(idz,3)%link
         !
         ! Update flags
         !
         if(ray_cart_dirz.gt.0) then
            ilr = 1
            amrray_icross = 6
         else
            ilr = 2
            amrray_icross = 5
         endif
      else
         !
         ! A regular grid is used, so set the ixyz_next indices
         !
         ! Depending on direction...
         !
         if(ray_cart_dirz.gt.0) then
            !
            ! Update flags
            !
            ilr = 1
            amrray_icross = 6
            !
            ! Update position
            !
            amrray_ix_next = amrray_ix_curr
            amrray_iy_next = amrray_iy_curr
            amrray_iz_next = amrray_iz_curr + 1
            !
            ! Check if out of bounds
            !
            if(amrray_iz_next.gt.amr_grid_nz) then
               amrray_ix_next = -1
               amrray_iy_next = -1
               amrray_iz_next = -1
            endif
         else
            !
            ! Update flags
            !
            ilr = 2
            amrray_icross = 5
            !
            ! Update position
            !
            amrray_ix_next = amrray_ix_curr
            amrray_iy_next = amrray_iy_curr
            amrray_iz_next = amrray_iz_curr - 1
            !
            ! Check if out of bounds
            !
            if(amrray_iz_next.lt.1) then
               amrray_ix_next = -1
               amrray_iy_next = -1
               amrray_iz_next = -1
            endif
         endif
      endif
      !
      ! Update position 
      !
      ray_cart_x = ray_cart_x + ray_ds*ray_cart_dirx
      ray_cart_y = ray_cart_y + ray_ds*ray_cart_diry
      ray_cart_z = axi(idz,3)
end select
!
! Now do a few things differently for AMR or regular grids
!
if(amr_tree_present) then
   !
   ! The AMR tree is present...
   !
   ! Now, if the new cell is in fact a nodal branch with children, then
   ! we must find this sub-cell
   !
   if(idir.ne.0) then
      if(associated(amrray_nextcell)) then
         if(.not.amrray_nextcell%leaf) then
            !
            ! Find this smaller cell
            !
            call amrray_find_subcell(amrray_nextcell,idir,ilr,&
                    ray_cart_x,ray_cart_y,ray_cart_z)
            !
            ! For the corner-based integration we must re-calculate icross
            !
            amrray_icross = -(2*(idir-1)+ilr)
         endif
      endif
   endif
   !
   ! Find index to next cell
   !
   if(associated(amrray_nextcell)) then
      ray_indexnext = amrray_nextcell%leafindex
   else
      ray_indexnext = 0
   endif
   !
   ! Now return the arrived argument and possibly the level of the next cell
   !
   if(associated(amrray_nextcell)) then
      arrived = .false.
      if(present(levelnext)) levelnext = amrray_nextcell%level
   else
      arrived = .true.
      if(present(levelnext)) levelnext = -1
   endif
else
   !
   ! Regular grid, no AMR tree...
   !
   ! Find index to next cell
   !
   if(amrray_ix_next.gt.0) then
      ray_indexnext = amrray_ix_next+(amrray_iy_next-1)*amr_grid_nx+(amrray_iz_next-1)*amr_grid_nx*amr_grid_ny
   else
      ray_indexnext = 0
   endif
   !
   ! Now return the arrived argument and possibly the level of the next cell
   !
   if(ray_indexnext.gt.0) then
      arrived = .false.
      if(present(levelnext)) levelnext = 0
   else
      arrived = .true.
      if(present(levelnext)) levelnext = -1
   endif
endif
!
! Done
!
return
!
end subroutine amrray_find_next_location_cart


!--------------------------------------------------------------------------
!                   Find next position (Spherical case)
!
! This routine computes the next position within the cell OR the location
! on the cell interface including the pointer to the next cell.
!
! In this version of the routine the x,y,z coordinates of the AMR grid
! are interpreted as r, theta and phi of spherical coordinates. 
!
! MAKE SURE TO SET THE FOLLOWING:
!   amrray_cell        Current cell pointer, will contain final cell pointer
!                      upon return.
!                      **** IS THIS STILL CURRENT? I THINK IT SHOULD BE: ray_indexcurr ****
!   ray_cart_x,_y,_z   Current position of photon
!   ray_cart_dirx/y/z  The direction unit vector. Must be normalized to 1.d0
!   ray_dsend          If next crossing is farther than dsend, then stop
!                      at dsend (i.e. before next crossing).
!   amrray_mirror_equator   If set, then the midplane (z=0) acts as 
!                      a mirror. Make sure that there exists no cells
!                      which are (whole or partly) below z=0 if you use this.
!   arrived            Is .true. if ray_ds becomes ray_dsend or if the
!                      end of the grid is reached. In other words: if
!                      arrived=.true. then the current step is either zero
!                      in length or it is the last step through a cell before
!                      arriving either at the final point or escaping to
!                      infinity (if ray_dsend=1d99). 
! 
! OPTIONAL:
!   maxdeltasina       If set, then the next position will not have an
!                      angle wrt to the current position that is larger
!                      than this value, or to be more accurate, the sinus
!                      of that angle will not be larger. This is useful
!                      when trying to avoid too large doppler shifts 
!                      within a given cell or when using the 2-D axisymmetric
!                      full scattering mode. If the next crossing with
!                      a cell wall would be at a larger delta angle, then
!                      the routine will instead return a location inside the
!                      cell (not at the cell wall) and the next cell index
!                      will be the same as the current cell index.
!
! THIS ROUTINE PRODUCES:
!   ray_ds        The actual distance travelled (cannot be bigger than dsend)
!   amrray_nextcell  Pointer to the cell in which the photon will next enter.
!                  If the photon stopped within the cell, then this will
!                  be equal to amrray_cell
!                  **** IS THIS STILL CURRENT? I THINK IT SHOULD BE: ray_indexnext ****
!   distmin (optional) The minimum distance to a cell wall 
!
!--------------------------------------------------------------------------
subroutine amrray_find_next_location_spher(ray_dsend,                &
                ray_cart_x,ray_cart_y,ray_cart_z,                    &
                ray_cart_dirx,ray_cart_diry,ray_cart_dirz,           &
                ray_indexcurr,ray_indexnext,ray_ds,arrived,          &
                distmin,levelnext,maxdeltasina)
implicit none
doubleprecision :: ray_cart_x,ray_cart_y,ray_cart_z
doubleprecision :: ray_cart_dirx,ray_cart_diry,ray_cart_dirz
doubleprecision :: ray_dsend,ray_ds
doubleprecision, optional :: distmin
integer,optional :: levelnext
doubleprecision, optional :: maxdeltasina
integer :: ray_indexcurr,ray_indexnext
logical :: arrived
integer :: ihit
doubleprecision :: axi(1:2,1:3),bxi(1:2,1:3)
doubleprecision :: angvec(1:3),sinang
!
doubleprecision :: fact,xt,yt,zt,dummy,ds_r1,ds_r2,ds_t1,ds_t2,ds_p1,ds_p2
integer :: idy,idz,ix,iy,iz,iddr,ilr,idir
!integer :: idx
doubleprecision :: r02,r0,rc0,theta0,phi0,cross_ds,ct2,st2,det,r002,r00
doubleprecision :: s00,sdet,sgnz,x00,y00,z00,x00dotdir,x0dotdir
doubleprecision :: ds_try,dum,dum1,dum2,pa,pb,pc,eps
doubleprecision :: r_r1,r_r2,r_t1,r_t2,r_p1,r_p2
doubleprecision :: theta_r1,theta_r2,theta_t1,theta_t2,theta_p1,theta_p2
doubleprecision :: phi_r1,phi_r2,phi_t1,phi_t2,phi_p1,phi_p2
doubleprecision :: r,r2,theta,phi,rcylnew
doubleprecision :: x_r1,x_r2,x_t1,x_t2,x_p1,x_p2
doubleprecision :: y_r1,y_r2,y_t1,y_t2,y_p1,y_p2
doubleprecision :: z_r1,z_r2,z_t1,z_t2,z_p1,z_p2
doubleprecision :: ds_sphere,ds_sphere0
integer :: isph
logical :: val_r1,val_r2,val_t1,val_t2,val_p1,val_p2
logical :: topquadrant,crossequator
logical :: doshift,dummyflag
doubleprecision :: margin,err,cossindum,pabc
logical :: hitmaxphi
!!######################################
!doubleprecision :: xbk,ybk,zbk
!!######################################
parameter(margin=5.d0)
!
! Do some self-checks
!
!!#############################################################
!write(stdo,*) '=- ',ray_cart_x,ray_cart_y,ray_cart_z
!!#############################################################
if(amrray_selfcheck) then
   if((amr_grid_xi(1,3).lt.0.d0).or.&
        (amr_grid_xi(amr_grid_nz+1,3).gt.twopi*(1.d0+tiny))) then
      write(stdo,*) 'ERROR in AMRRay Module: Theta grid exceeds [0,pi]'
      write(stdo,*) '  Theta = [',amr_grid_xi(1,3),',',&
           amr_grid_xi(amr_grid_nz+1,3),']'
      stop 80031
   endif
   if(amrray_mirror_equator) then
      if((amr_grid_xi(1,2).lt.0.d0).or.                              &
           (amr_grid_xi(amr_grid_ny+1,2).gt.pihalf*(1.d0+tiny)))     &
           stop 80032
      if(ray_cart_z.lt.0.d0) then
         write(stdo,*) 'ERROR: Mirroring is active, but z<0...'
         stop 80033
      endif
   endif
   if(amrray_initialized.ne.12345) then
      write(stdo,*) 'ERROR: Did not initalize AMRRay yet.'
      stop
   endif
   dummy = sqrt(ray_cart_dirx**2+ray_cart_diry**2+ray_cart_dirz**2) 
   if(abs(dummy-1.d0).gt.1d-10) then
      write(stdo,*) 'ERROR: Direction vector is not unit vector'
      write(stdo,*) ray_cart_dirx,ray_cart_diry,ray_cart_dirz
      stop 2913
   endif
endif
!
! Check
!
if(present(distmin)) then
   write(stdo,*) 'ERROR: For now the spherical mode does not allow'
   write(stdo,*) '       the computation of the smallest distance to a wall'
   stop
endif
!
! Set flags to initial value
!
arrived  = .false.
crossequator = .false.
ds_sphere0 = 1d99
amrray_ispherehit = 0
amrray_icross = 0
hitmaxphi = .false.
!
! If mirror symmetry required, first check and if necessary, flip
! NOTE: Let us simply not accept ray positions that are on the wrong side.
!
if(amrray_mirror_equator) then
   if(ray_cart_z.le.0.d0) then
      if(ray_cart_z.lt.0.d0) then
         stop 2065
!         ray_cart_z    = -ray_cart_z
!         ray_cart_dirz = -ray_cart_dirz 
      elseif(ray_cart_dirz.lt.0.d0) then
         stop 2066
!         ray_cart_dirz = -ray_cart_dirz
      endif
   endif
endif
!
!!######################################
!xbk = ray_cart_x
!ybk = ray_cart_y
!zbk = ray_cart_z
!!######################################
!
! Determine location in (r,theta,phi)
! Note: phi in [0,2*pi], theta in [0,pi]. Possibly mirror in theta=pi/2.
! Note: most of this is only necessary for self-checking, but if the
! amrray_cell is not associated, then we need all this nevertheless
!
r02    = ray_cart_x*ray_cart_x + ray_cart_y*ray_cart_y + ray_cart_z*ray_cart_z
r0     = sqrt(r02)
theta0 = acos(ray_cart_z/(r0+1d-199))
if(ray_cart_x.eq.0.d0) then
   if(ray_cart_y.ge.0.d0) then
      phi0   = pihalf
   else
      phi0   = 3.d0*pihalf 
   endif
else
   phi0   = atan(ray_cart_y/ray_cart_x)
   if(ray_cart_x.gt.0.d0) then
      if(ray_cart_y.lt.0.d0) then
         phi0   = phi0 + twopi
      endif
   else
      phi0   = phi0 + pi
   endif
endif
!################################
if((phi0.lt.0.d0).or.(phi0.gt.twopi)) stop 1309
!################################
!
! If cell index = 0, then go look for the cell, else simply link to the cell
!
if(ray_indexcurr.gt.0) then
   !
   ! We are in a cell. Get cell pointer (or ix,iy,iz) and
   ! find the cell walls
   !
   if(amr_tree_present) then
      !
      ! If the AMR tree is present, then link the amrray_cell
      !
      amrray_cell => amr_index_to_leaf(ray_indexcurr)%link
      if(.not.associated(amrray_cell)) then
         write(stdo,*) 'ERROR in amrray module: Cell ',ray_indexcurr,   &
              ' is not assigned. It does not exist.'
         stop 6023
      endif
      do iddr=1,3
         axi(1,iddr) = amr_finegrid_xi(amrray_cell%ixyzf(iddr),iddr,amrray_cell%level)
         axi(2,iddr) = amr_finegrid_xi(amrray_cell%ixyzf(iddr)+1,iddr,amrray_cell%level)
      enddo
   else
      !
      ! If we have a regular grid, then find the amrray_ix_curr, amrray_iy_curr, amrray_iz_curr
      !
      call amr_regular_get_ixyz(ray_indexcurr,amrray_ix_curr,amrray_iy_curr,amrray_iz_curr)
      axi(1,1) = amr_finegrid_xi(amrray_ix_curr,1,0)
      axi(2,1) = amr_finegrid_xi(amrray_ix_curr+1,1,0)
      axi(1,2) = amr_finegrid_xi(amrray_iy_curr,2,0)
      axi(2,2) = amr_finegrid_xi(amrray_iy_curr+1,2,0)
      axi(1,3) = amr_finegrid_xi(amrray_iz_curr,3,0)
      axi(2,3) = amr_finegrid_xi(amrray_iz_curr+1,3,0)
   endif
   !
   ! Since we computed phi0 from x,y,z, we may have an uncertainty of
   ! phi0 near the phi0=0,2*pi boundary. Check this.
   !
   if((phi0.le.small).or.(phi0.ge.twopi-small)) then
      !
      ! Yes, we may have this problem. So check out in which regime
      ! we are.
      !
      if(axi(1,3).eq.0.d0) then
         phi0=0.d0
      elseif(axi(2,3).eq.twopi) then
         phi0=twopi
      endif
   endif
elseif(ray_indexcurr.eq.0) then
   !
   ! We are outside of a cell
   !
   nullify(amrray_cell)
   amrray_ix_curr = -1
   amrray_iy_curr = -1
   amrray_iz_curr = -1
else
   !
   ! We do not know, so we must find out ourselves
   !
   ! Since we computed phi0 from x,y,z, we may have an uncertainty of
   ! phi0 near the phi0=0,2*pi boundary. Check this: If we are moving 
   ! toward negative phi, then if we are close to phi0=0,2*pi, then we
   ! must choose 2*pi. Etc.
   !
   if(phi0.le.small) then
      if(ray_cart_diry.lt.0.d0) then
         phi0=twopi
      endif
   elseif(phi0.ge.twopi-small) then
      if(ray_cart_diry.gt.0.d0) then
         phi0=0.d0
      endif
   endif
   !
   ! Find the cell
   !
   if(amr_tree_present) then
      !
      ! If the AMR tree is present, then we use that
      !
      call amr_findcell(r0,theta0,phi0,amrray_cell)
      if(associated(amrray_cell)) then
         ray_indexcurr = amrray_cell%leafindex
         do iddr=1,3
            axi(1,iddr) = amr_finegrid_xi(amrray_cell%ixyzf(iddr),iddr,amrray_cell%level)
            axi(2,iddr) = amr_finegrid_xi(amrray_cell%ixyzf(iddr)+1,iddr,amrray_cell%level)
         enddo
      else
         ray_indexcurr = 0
      endif
   else
      !
      ! If we have a regular grid, then use amrray_ix_curr, amrray_iy_curr, amrray_iz_curr
      !
      call amr_findbasecell(r0,theta0,phi0,amrray_ix_curr,amrray_iy_curr,amrray_iz_curr)
      if(amrray_ix_curr.gt.0) then
         ray_indexcurr = amrray_ix_curr+(amrray_iy_curr-1)*amr_grid_nx+(amrray_iz_curr-1)*amr_grid_nx*amr_grid_ny
         axi(1,1) = amr_finegrid_xi(amrray_ix_curr,1,0)
         axi(2,1) = amr_finegrid_xi(amrray_ix_curr+1,1,0)
         axi(1,2) = amr_finegrid_xi(amrray_iy_curr,2,0)
         axi(2,2) = amr_finegrid_xi(amrray_iy_curr+1,2,0)
         axi(1,3) = amr_finegrid_xi(amrray_iz_curr,3,0)
         axi(2,3) = amr_finegrid_xi(amrray_iz_curr+1,3,0)
      else
         ray_indexcurr = 0
      endif
   endif
endif
!
! If outside of grid, then do special stuff
!
if(ray_indexcurr.le.0) then 
   !
   !========================================================================
   !                  We are outside of the grid
   !========================================================================
   !
   ! Here we assume that the r0,theta0,phi0 have already been calculated
   !
   ! Check that we are indeed outside of the grid. We put a small margin
   ! around it. That is: if we are outside of the grid, then we must be
   ! sufficiently far away from the grid boundary to avoid round-off
   ! errors causing problems. We use "small" for that instead of "tiny",
   ! so that we have a margin of 1D-12 instead of 1D-14. This means that
   ! we should not allow any stars or an observer to be inside this margin.
   !
   if(  (r0.ge.amr_grid_xi(1,1)*(1.d0-small)).and.                 &
        (r0.le.amr_grid_xi(amr_grid_nx+1,1)*(1.d0+small)).and.     &
        (theta0.ge.amr_grid_xi(1,2)-small).and.                    &
        (theta0.le.amr_grid_xi(amr_grid_ny+1,2)+small).and.        &
        (phi0.ge.amr_grid_xi(1,3)-small).and.                      &
        (phi0.le.amr_grid_xi(amr_grid_nz+1,3)+small)) then
      if(  (r0.ge.amr_grid_xi(1,1)).and.                        &
           (r0.le.amr_grid_xi(amr_grid_nx+1,1)).and.            &
           (theta0.ge.amr_grid_xi(1,2)).and.                    &
           (theta0.le.amr_grid_xi(amr_grid_ny+1,2)).and.        &
           (phi0.ge.amr_grid_xi(1,3)).and.                      &
           (phi0.le.amr_grid_xi(amr_grid_nz+1,3))) then
         write(stdo,*) 'ERROR: Cell index set to 0, so expect to be outside'
         write(stdo,*) '       of the grid. But I find myself inside the grid.'
         write(stdo,*) '       Clearly an internal error in the code.'
         write(stdo,*) ray_cart_x,ray_cart_y,ray_cart_z
         write(stdo,*) r0,theta0,phi0
         write(stdo,*) ray_cart_dirx,ray_cart_diry,ray_cart_dirz
         stop 1090
      else
         write(stdo,*) 'PROBLEM: It appears that a star or other origin'
         write(stdo,*) '         of rays lies just outside the grid, but within'
         write(stdo,*) '         the safety margin of 1D-12 (relative) of'
         write(stdo,*) '         the grid edge. For algorithmic reasons this'
         write(stdo,*) '         is not allowed. Please put your star or other'
         write(stdo,*) '         photon source either clearly inside or further'
         write(stdo,*) '         outside the grid.'
         stop 1091
      endif
   endif
   !
   ! Set as default the following flags to true (important!) 
   !
   val_r1 = .true.
   val_r2 = .true.
   val_t1 = .true.
   val_t2 = .true.
   val_p1 = .true.
   val_p2 = .true.
   !
   ! Compute, for later, the x.dx
   !
   x0dotdir   = ray_cart_x * ray_cart_dirx +              &
                ray_cart_y * ray_cart_diry +              &
                ray_cart_z * ray_cart_dirz
   !
   ! Now, for numerical accuracy reasons (related to the fact that for the
   ! formulae below we are solving quadratic equations; see below for more
   ! explanation), if the direction of the ray is pointing inward, we first
   ! move linearly along the ray a distance equal to r0.
   !
   if(x0dotdir.ge.0.d0) then
      !
      ! Ray is outward pointing, so no shift is necessary
      !
      s00  = 0.d0
      x00  = ray_cart_x
      y00  = ray_cart_y
      z00  = ray_cart_z
      !
   else
      !
      ! Ray is inward pointing, so do the shift for ensuring stability of
      ! the quadratic equations in case of zooming into extremely small
      ! regions close to the center
      !
      s00  = r0
      x00  = ray_cart_x + s00 * ray_cart_dirx
      y00  = ray_cart_y + s00 * ray_cart_diry
      z00  = ray_cart_z + s00 * ray_cart_dirz
      !
   endif
   !
   ! Some preliminary calculations
   !
   x00dotdir  = x00 * ray_cart_dirx +              &
                y00 * ray_cart_diry +              &
                z00 * ray_cart_dirz
   r002       = x00**2 + y00**2 + z00**2 
   r00        = sqrt( r002 )
   !
   ! If r00 is very very much smaller than r0, then the remaining non-zeroness
   ! is likely a numerical error, and in fact a perfectly radial inward
   ! ray is in fact meant. In that case do not attempt to compute crossings
   ! with theta and phi coordinates, because the formulae would become
   ! singular.
   !
   ! ***** IS TINY NOT TOO SMALL? *****
   !
   if(r00.lt.tiny*r0) then
      val_t1 = .false.
      val_t2 = .false.
      val_p1 = .false.
      val_p2 = .false.
   endif
   !
   ! Same is true for radially outward moving rays
   !
   ! ***** IS TINY NOT TOO SMALL? *****
   !
   dum = (ray_cart_x - r0 * ray_cart_dirx)**2 +    &
         (ray_cart_y - r0 * ray_cart_diry)**2 +    &
         (ray_cart_z - r0 * ray_cart_dirz)**2
   if(dum.lt.(tiny*r0)**2) then
      val_t1 = .false.
      val_t2 = .false.
      val_p1 = .false.
      val_p2 = .false.
   endif
   !
   ! In case of mirror symmetry, we do not need to check for the crossing
   ! with the theta2 boundary, because this should (!) be exactly at the
   ! equatorial plane.
   !
   if(amrray_mirror_equator) then
      val_t2 = .false.
   endif
   !
   ! If cyclic boundary conditions in phi are chosen, or if phi is not
   ! considered, then we do not need to find crossings with the phi
   ! grid boundaries
   !
   if(amr_cyclic_xyz(3).or.(amr_zdim.eq.0)) then
      val_p1 = .false.
      val_p2 = .false.
   endif
   !
   ! If 1-D spherically symmetric radiative transfer, then do not search
   ! for theta crossings
   !
   if(amr_ydim.eq.0) then
      val_t1 = .false.
      val_t2 = .false.
   endif
   !
   ! -------------
   ! Now solve the equation for the crossing with inner R=const surface
   ! and check if the crossing happens within the theta, phi grid
   !
   !   s = -(x00.dirx) +- sqrt( (x00.dirx)^2 + r^2 - r00^2 ) + s00
   !
   ! The s00 must be added because we do the calculation now from a
   ! different point along the ray (not the true starting point).  Without
   ! adding the s00 we obtain the s as measured from this different
   ! point. By adding s00 we get the actual s, as measured from the true
   ! starting point. The reason for this complicated stuff is that the
   ! quadratic equations may become inaccurate if the ray points very
   ! sharply toward the center and crosses the grid only very close to the
   ! center. If the ratio of the radius where this crossing happens and the
   ! radius where the ray starts is larger than about 1d-8, then the square
   ! of this is 1d-16, meaning that we are in the numerical noise of the
   ! double precision. A ratio of 1d-6 is therefore the most extreme one can
   ! afford UNLESS we do the shift as described above. The shift is a linear
   ! calculation, so with a ratio of 1d-8 we still have 8 digits of accuracy
   ! left in the shift. Once the shift is done, the quadratic equations are
   ! "renormalized" and such problems with accuracy do no longer take
   ! place. So when would we need this? Well, 1 parsec / 1 Rsun = 4.4d7, so
   ! if we model a huge cloud around a star we may already get into this
   ! regime of dynamic range, and the above shift/renormalization of the
   ! equations becomes necessary.
   ! -------------
   !
   r         = amr_grid_xi(1,1)
   r2        = r*r
   det       = x00dotdir*x00dotdir + r2 - r002
   if(det.le.0.d0) then
      !
      ! No solution at all
      !
      ds_r1 = -1d99
   else
      !
      ! Possibly a solution
      ! 
      ! If it is, it must be the largest of the two, because that is,
      ! by definition, the one that crosses the r=r1 from smaller r,
      ! i.e. from outside the grid. 
      !
      ds_r1 = -x00dotdir + sqrt(det) + s00
   endif
   !
   ! If ds_r1 is positive, then we have a true solution, but we do not
   ! yet know if it is within the grid in theta and phi
   !
   if(ds_r1.le.0.d0) then
      !
      ! ds_r1 is negative, so no solution
      !
      val_r1 = .false.
   else
      !
      ! Now create the new x,y,z coordinates
      !
      x_r1   = ray_cart_x + ds_r1 * ray_cart_dirx
      y_r1   = ray_cart_y + ds_r1 * ray_cart_diry
      z_r1   = ray_cart_z + ds_r1 * ray_cart_dirz
      !
      ! Do a consistency check
      ! NOTE: This error could easily be much larger than "tiny" or "small",
      !       because of linear precision loss if the radial grid spans
      !       a large dynamic range. But we will correct for this later.
      !       Here we merely check if things don't go badly wrong.
      !
      if(amrray_selfcheck) then
         r_r1   = sqrt(x_r1*x_r1+y_r1*y_r1+z_r1*z_r1)
         if(abs(r_r1/r-1.d0).gt.1d-5) then
            write(stdo,*) 'ERROR in AMRRay Module: Inaccuracy in solving'
            write(stdo,*) '   crossing with inner grid radius.'
            write(stdo,*) '   Error = ',abs(r_r1/r-1.d0)
            stop 8329
         endif
      endif
      !
      ! Now create the new r,theta,phi coordinates but immediately
      ! check if they are on the grid
      !
      r_r1   = amr_grid_xi(1,1)
      !
      ! Compute theta
      !
      theta_r1   = acos(z_r1/r_r1)
      !
      ! Check if theta_r1 is on grid
      !
      ! We use a tiny extra margin to prevent the ray from accidently
      ! "slipping through a corner"
      !
      if(amrray_mirror_equator) then
         if((theta_r1.le.amr_grid_xi(1,2)-tiny).or.&
              (theta_r1.ge.pi-amr_grid_xi(1,2)+tiny)) then
            val_r1 = .false.
         endif
      else
         if((theta_r1.le.amr_grid_xi(1,2)-tiny).or.&
              (theta_r1.ge.amr_grid_xi(amr_grid_ny+1,2)+tiny)) then
            val_r1 = .false.
         endif
      endif
      !
      ! Now compute phi
      !
      if(x_r1.eq.0.d0) then
         if(y_r1.ge.0.d0) then
            phi_r1  = pihalf
         else
            phi_r1  = 3.d0*pihalf 
         endif
      else
         phi_r1   = atan(y_r1/x_r1)
         if(x_r1.gt.0.d0) then
            if(y_r1.lt.0.d0) then
               phi_r1   = phi_r1 + twopi
            endif
         else
            phi_r1   = phi_r1 + pi
         endif
      endif
      !
      ! Check if phi on grid
      !
      ! We use a tiny extra margin to prevent the ray from accidently
      ! "slipping through a corner"
      !
      if((amr_zdim.ne.0).and.(.not.amr_cyclic_xyz(3))) then
         if(phi_r1.lt.amr_grid_xi(1,3)-tiny) then
            val_r1 = .false.
         endif
         if(phi_r1.gt.amr_grid_xi(amr_grid_nz+1,3)+tiny) then
            val_r1 = .false.
         endif
      endif
   endif
   !
   ! -------------
   ! Now solve the equation for the crossing with outer R=const surface
   ! and check if the crossing happens within the theta, phi grid
   !
   !   s = -(x00.dirx) +- sqrt( (x00.dirx)^2 + r^2 - r00^2 ) + s00
   !
   ! The same issues with stiffness are accounted for here, although
   ! it is not expected to be necessary
   ! -------------
   !
   r         = amr_grid_xi(amr_grid_nx+1,1)
   r2        = r*r
   det       = x00dotdir*x00dotdir + r2 - r002
   ds_r2     = -1d99
   if(det.gt.0.d0) then
      if((x0dotdir.lt.0.d0).and.(r2.le.r02)) then
         ds_r2 = -x00dotdir - sqrt(det) + s00
      endif
   endif
   !
   ! If ds_r2 is positive, then we have a true solution, but we do not
   ! yet know if it is within the grid in theta and phi
   !
   if(ds_r2.lt.0.d0) then
      !
      ! ds_r2 is negative, so no solution
      !
      val_r2 = .false.
   else
      !
      ! Now create the new x,y,z coordinates
      !
      x_r2   = ray_cart_x + ds_r2 * ray_cart_dirx
      y_r2   = ray_cart_y + ds_r2 * ray_cart_diry
      z_r2   = ray_cart_z + ds_r2 * ray_cart_dirz
      !
      ! Do a consistency check
      !
      if(amrray_selfcheck) then
         r_r2   = sqrt(x_r2*x_r2+y_r2*y_r2+z_r2*z_r2)
         if(abs(r_r2/r-1.d0).gt.1d-5) then
            write(stdo,*) 'ERROR in AMRRay Module: Inaccuracy in solving'
            write(stdo,*) '   crossing with outer grid radius.'
            write(stdo,*) '   Error = ',abs(r_r2/r-1.d0)
            stop 8330
         endif
      endif
      !
      ! Now create the new r,theta,phi coordinates but immediately
      ! check if they are on the grid
      !
      r_r2   = amr_grid_xi(amr_grid_nx+1,1)
      !
      ! Compute theta
      !
      theta_r2   = acos(z_r2/r_r2)
      !
      ! Check if theta_r2 is on grid
      !
      ! We use a tiny extra margin to prevent the ray from accidently
      ! "slipping through a corner"
      !
      if(amrray_mirror_equator) then
         if((theta_r2.le.amr_grid_xi(1,2)-tiny).or.&
              (theta_r2.ge.pi-amr_grid_xi(1,2)+tiny)) then
            val_r2 = .false.
         endif
      else
         if((theta_r2.le.amr_grid_xi(1,2)-tiny).or.&
              (theta_r2.ge.amr_grid_xi(amr_grid_ny+1,2)+tiny)) then
            val_r2 = .false.
         endif
      endif
      !
      ! Now compute phi
      !
      if(x_r2.eq.0.d0) then
         if(y_r2.ge.0.d0) then
            phi_r2  = pihalf
         else
            phi_r2  = 3.d0*pihalf 
         endif
      else
         phi_r2   = atan(y_r2/x_r2)
         if(x_r2.gt.0.d0) then
            if(y_r2.lt.0.d0) then
               phi_r2   = phi_r2 + twopi
            endif
         else
            phi_r2   = phi_r2 + pi
         endif
      endif
      !
      ! Check if phi on grid
      !
      ! We use a tiny extra margin to prevent the ray from accidently
      ! "slipping through a corner"
      !
      if((amr_zdim.ne.0).and.(.not.amr_cyclic_xyz(3))) then
         if(phi_r2.lt.amr_grid_xi(1,3)-tiny) then
            val_r2 = .false.
         endif
         if(phi_r2.gt.amr_grid_xi(amr_grid_nz+1,3)+tiny) then
            val_r2 = .false.
         endif
      endif
   endif
   !
   ! -------------
   ! Now let us try to find the crossings with the theta1=const cone.
   ! -------------
   !
   if(val_t1) then
      st2   = amrray_sint1**2
      ct2   = amrray_cost1**2
      pa    = ct2*(ray_cart_dirx**2+ray_cart_diry**2)-st2*ray_cart_dirz**2
      pb    = 2*ct2*(ray_cart_dirx*x00 + ray_cart_diry*y00) - &
              2*st2*ray_cart_dirz*z00
      pc    = ct2*(x00**2+y00**2)-st2*z00**2
      det   = pb**2 - 4*pa*pc
      if(det.le.0.d0) then
         !
         ! No solution at all
         !
         ds_t1 = -1d99
         val_t1 = .false.
      else
         !
         ! Possibly a solution
         !
         sdet   = sqrt(det)
         !
         ! Compute a value to check if we might be encountering a
         ! near-cancellation, and thus loss of precision, in which
         ! case we must switch to a special treatment with a Taylor-series.
         ! BUGFIX: 06.03.2015 (Thanks to Andrea Isella for reporting the bug)
         !
         pabc  = 4*pa*pc/(pb*pb)
         !
         ! Now check which of the two solutions to take
         !
         if(ray_cart_dirx*ray_cart_x+ray_cart_diry*ray_cart_y.le.0.d0) then
            !
            ! If there is a solution, then it is "around the corner", 
            ! in the sense that the ray first moves toward the z-axis,
            ! then away again and THEN hits the theta grid surface.
            ! It might, during its approach to the z-axis, cross the
            ! same theta=const cone, but that is then not a solution.
            !
            if((abs(pabc).gt.1d-4).or.(pb.lt.0.d0)) then
               !
               ! Use the full formula
               !
               ds_t1  = - 0.5d0 * ( pb - sdet ) / pa + s00
            else
               !
               ! Use a Taylor expansion to avoid precision loss due to
               ! near-cancellation of pb and sdet
               ! BUGFIX: 06.03.2015 (Thanks to Andrea Isella for reporting the bug)
               !
               ds_t1  = - 0.5d0 * pb * ( 0.5d0*pabc + (1.d0/8.d0)*pabc*pabc + &
                                         (1.d0/16.d0)*pabc**3 +               &
                                         (5.d0/128.d0)*pabc**4 ) / pa + s00
            endif
            z_t1   = ray_cart_z + ds_t1 * ray_cart_dirz
            if((ds_t1.le.0.d0).or.(z_t1*theta1_sgn.lt.0.d0)) then
               val_t1 = .false.
            endif
         else
            !
            ! Normal case: try out one, if that does not work, do the 
            ! other.
            !
            if((abs(pabc).gt.1d-4).or.(pb.gt.0.d0)) then
               !
               ! Use the full formula
               !
               ds_t1  = - 0.5d0 * ( pb + sdet ) / pa + s00
            else
               !
               ! Use a Taylor expansion to avoid precision loss due to
               ! near-cancellation of pb and sdet
               ! BUGFIX: 06.03.2015 (Thanks to Andrea Isella for reporting the bug)
               !
               ds_t1  = - 0.5d0 * pb * ( 0.5d0*pabc + (1.d0/8.d0)*pabc*pabc + &
                                         (1.d0/16.d0)*pabc**3 +               &
                                         (5.d0/128.d0)*pabc**4 ) / pa + s00
            endif
            z_t1   = ray_cart_z + ds_t1 * ray_cart_dirz
            if((ds_t1.le.0.d0).or.(z_t1*theta1_sgn.lt.0.d0)) then
               if((abs(pabc).gt.1d-4).or.(pb.lt.0.d0)) then
                  !
                  ! Use the full formula
                  !
                  ds_t1  = - 0.5d0 * ( pb - sdet ) / pa + s00
               else
                  !
                  ! Use a Taylor expansion to avoid precision loss due to
                  ! near-cancellation of pb and sdet
                  ! BUGFIX: 06.03.2015 (Thanks to Andrea Isella for reporting the bug)
                  !
                  ds_t1  = - 0.5d0 * pb * ( 0.5d0*pabc + (1.d0/8.d0)*pabc*pabc + &
                                            (1.d0/16.d0)*pabc**3 +               &
                                            (5.d0/128.d0)*pabc**4 ) / pa + s00
               endif
               z_t1   = ray_cart_z + ds_t1 * ray_cart_dirz
               if((ds_t1.le.0.d0).or.(z_t1*theta1_sgn.lt.0.d0)) then
                  val_t1 = .false.
               endif
            endif
         endif
      endif
      if(val_t1) then
         !
         ! Now create the new x,y,z coordinates
         !
         x_t1   = ray_cart_x + ds_t1 * ray_cart_dirx
         y_t1   = ray_cart_y + ds_t1 * ray_cart_diry
         z_t1   = ray_cart_z + ds_t1 * ray_cart_dirz
         !
         ! Now create the new r,theta,phi coordinates but immediately
         ! check if they are on the grid
         !
         ! Create theta
         !
         theta_t1 = amr_grid_xi(1,2)
         !
         ! Create r
         !
         r_t1   = sqrt( x_t1**2 + y_t1**2 + z_t1**2 )
         !
         ! Check if r_t1 is on grid
         !
         ! We use a tiny extra margin to prevent the ray from accidently
         ! "slipping through a corner"
         !
         if((r_t1.lt.amr_grid_xi(1,1)*oneminus).or.                    &
            (r_t1.gt.amr_grid_xi(amr_grid_nx+1,1)*oneplus)) then
            val_t1 = .false.
         endif
         !
         ! Now compute phi
         !
         if(x_t1.eq.0.d0) then
            if(y_t1.ge.0.d0) then
               phi_t1  = pihalf
            else
               phi_t1  = 3.d0*pihalf 
            endif
         else
            phi_t1   = atan(y_t1/x_t1)
            if(x_t1.gt.0.d0) then
               if(y_t1.lt.0.d0) then
                  phi_t1   = phi_t1 + twopi
               endif
            else
               phi_t1   = phi_t1 + pi
            endif
         endif
         !
         ! Check if phi on grid
         !
         ! We use a tiny extra margin to prevent the ray from accidently
         ! "slipping through a corner"
         !
         if((amr_zdim.ne.0).and.(.not.amr_cyclic_xyz(3))) then
            if(phi_t1.lt.amr_grid_xi(1,3)-tiny) then
               val_t1 = .false.
            endif
            if(phi_t1.gt.amr_grid_xi(amr_grid_nz+1,3)+tiny) then
               val_t1 = .false.
            endif
         endif
      endif
   endif
   !
   ! -------------
   ! Now let us try to find the crossings with the theta2=const cone.
   ! -------------
   !
   if(val_t2) then
      st2   = amrray_sint2**2
      ct2   = amrray_cost2**2
      pa    = ct2*(ray_cart_dirx**2+ray_cart_diry**2)-st2*ray_cart_dirz**2
      pb    = 2*ct2*(ray_cart_dirx*x00 + ray_cart_diry*y00) - &
              2*st2*ray_cart_dirz*z00
      pc    = ct2*(x00**2+y00**2)-st2*z00**2
      det   = pb**2 - 4*pa*pc
      if(det.le.0.d0) then
         !
         ! No solution at all
         !
         ds_t2 = -1d99
         val_t2 = .false.
      else
         !
         ! Possibly a solution
         !
         sdet   = sqrt(det)
         !
         ! Compute a value to check if we might be encountering a
         ! near-cancellation, and thus loss of precision, in which
         ! case we must switch to a special treatment with a Taylor-series.
         ! BUGFIX: 06.03.2015 (Thanks to Andrea Isella for reporting the bug)
         !
         pabc  = 4*pa*pc/(pb*pb)
         !
         ! Now check which of the two solutions to take
         !
         if(thetagrid_cross_equator.and.&
            (ray_cart_dirx*ray_cart_x+ray_cart_diry*ray_cart_y.le.0.d0)) then
            !
            ! If there is a solution, then it is "around the corner", 
            ! in the sense that the ray first moves toward the z-axis,
            ! then away again and THEN hits the theta grid surface.
            ! It might, during its approach to the z-axis, cross the
            ! same theta=const cone, but that is then not a solution.
            !
            if((abs(pabc).gt.1d-4).or.(pb.lt.0.d0)) then
               !
               ! Use the full formula
               !
               ds_t2  = - 0.5d0 * ( pb - sdet ) / pa + s00
            else
               !
               ! Use a Taylor expansion to avoid precision loss due to
               ! near-cancellation of pb and sdet
               ! BUGFIX: 06.03.2015 (Thanks to Andrea Isella for reporting the bug)
               !
               ds_t2  = - 0.5d0 * pb * ( 0.5d0*pabc + (1.d0/8.d0)*pabc*pabc + &
                                         (1.d0/16.d0)*pabc**3 +               &
                                         (5.d0/128.d0)*pabc**4 ) / pa + s00
            endif
            z_t2   = ray_cart_z + ds_t2 * ray_cart_dirz
            if((ds_t2.le.0.d0).or.(z_t2*theta2_sgn.lt.0.d0)) then
               val_t2 = .false.
            endif
         else
            !
            ! Normal case: try out one, if that does not work, do the 
            ! other.
            !
            if((abs(pabc).gt.1d-4).or.(pb.gt.0.d0)) then
               !
               ! Use the full formula
               !
               ds_t2  = - 0.5d0 * ( pb + sdet ) / pa + s00
            else
               !
               ! Use a Taylor expansion to avoid precision loss due to
               ! near-cancellation of pb and sdet
               ! BUGFIX: 06.03.2015 (Thanks to Andrea Isella for reporting the bug)
               !
               ds_t2  = - 0.5d0 * pb * ( 0.5d0*pabc + (1.d0/8.d0)*pabc*pabc + &
                                         (1.d0/16.d0)*pabc**3 +               &
                                         (5.d0/128.d0)*pabc**4 ) / pa + s00
            endif
            z_t2   = ray_cart_z + ds_t2 * ray_cart_dirz
            if((ds_t2.le.0.d0).or.(z_t2*theta2_sgn.lt.0.d0)) then
               if((abs(pabc).gt.1d-4).or.(pb.lt.0.d0)) then
                  !
                  ! Use the full formula
                  !
                  ds_t2  = - 0.5d0 * ( pb - sdet ) / pa + s00
               else
                  !
                  ! Use a Taylor expansion to avoid precision loss due to
                  ! near-cancellation of pb and sdet
                  ! BUGFIX: 06.03.2015 (Thanks to Andrea Isella for reporting the bug)
                  !
                  ds_t2  = - 0.5d0 * pb * ( 0.5d0*pabc + (1.d0/8.d0)*pabc*pabc + &
                                            (1.d0/16.d0)*pabc**3 +               &
                                            (5.d0/128.d0)*pabc**4 ) / pa + s00
               endif
               z_t2   = ray_cart_z + ds_t2 * ray_cart_dirz
               if((ds_t2.le.0.d0).or.(z_t2*theta2_sgn.lt.0.d0)) then
                  val_t2 = .false.
               endif
            endif
         endif
      endif
      if(val_t2) then
         !
         ! Now create the new x,y,z coordinates
         !
         x_t2   = ray_cart_x + ds_t2 * ray_cart_dirx
         y_t2   = ray_cart_y + ds_t2 * ray_cart_diry
         z_t2   = ray_cart_z + ds_t2 * ray_cart_dirz
         !
         ! Now create the new r,theta,phi coordinates but immediately
         ! check if they are on the grid
         !
         ! Create theta
         !
         if(amrray_mirror_equator) then
            theta_t2 = pi-amr_grid_xi(1,2)
         else
            theta_t2 = amr_grid_xi(amr_grid_ny+1,2)
         endif
         !
         ! Create r
         !
         r_t2   = sqrt( x_t2**2 + y_t2**2 + z_t2**2 )
         !
         ! Check if r_t2 is on grid
         !
         ! We use a tiny extra margin to prevent the ray from accidently
         ! "slipping through a corner"
         !
         if((r_t2.lt.amr_grid_xi(1,1)*oneminus).or.                    &
            (r_t2.gt.amr_grid_xi(amr_grid_nx+1,1)*oneplus)) then
            val_t2 = .false.
         endif
         !
         ! Now compute phi
         !
         if(x_t2.eq.0.d0) then
            if(y_t2.ge.0.d0) then
               phi_t2  = pihalf
            else
               phi_t2  = 3.d0*pihalf 
            endif
         else
            phi_t2   = atan(y_t2/x_t2)
            if(x_t2.gt.0.d0) then
               if(y_t2.lt.0.d0) then
                  phi_t2   = phi_t2 + twopi
               endif
            else
               phi_t2   = phi_t2 + pi
            endif
         endif
         !
         ! Check if phi on grid
         !
         ! We use a tiny extra margin to prevent the ray from accidently
         ! "slipping through a corner"
         !
         if((amr_zdim.ne.0).and.(.not.amr_cyclic_xyz(3))) then
            if(phi_t2.lt.amr_grid_xi(1,3)-tiny) then
               val_t2 = .false.
            endif
            if(phi_t2.gt.amr_grid_xi(amr_grid_nz+1,3)+tiny) then
               val_t2 = .false.
            endif
         endif
      endif
   endif
   !
   ! -------------
   ! Now let us try to find the crossings with the phi1=const cone.
   ! -------------
   !
   if(val_p1) then
      !
      ! Find the distance to the crossing
      !
      dum1 = ray_cart_x*amrray_sinp1 - ray_cart_y*amrray_cosp1
      dum2 = ray_cart_diry*amrray_cosp1 - ray_cart_dirx*amrray_sinp1
      if(dum2.eq.0.d0) then
         val_p1 = .false.
         ds_p1 = -1d99
      else
         ds_p1 = dum1 / dum2
      endif
      if(ds_p1.le.0.d0) then
         val_p1 = .false.
      endif
      !
      ! Now create the new x,y,z coordinates
      !
      x_p1   = ray_cart_x + ds_p1 * ray_cart_dirx
      y_p1   = ray_cart_y + ds_p1 * ray_cart_diry
      z_p1   = ray_cart_z + ds_p1 * ray_cart_dirz
      !
      ! Find the r, theta, phi
      !
      r_p1     = sqrt( x_p1**2 + y_p1**2 + z_p1**2 )
      theta_p1 = acos(z_p1/r_p1)
      if(x_p1.eq.0.d0) then
         if(y_p1.ge.0.d0) then
            phi_p1 = pihalf
         else
            phi_p1 = 3.d0*pihalf 
         endif
      else
         phi_p1 = atan(y_p1/x_p1)
         if(x_p1.gt.0.d0) then
            if(y_p1.lt.0.d0) then
               phi_p1 = phi_p1 + twopi
            endif
         else
            phi_p1 = phi_p1 + pi
         endif
      endif
      if((phi_p1.lt.0.d0).or.(phi_p1.gt.twopi)) stop 1309
      !
      ! Check if phi_p1 is the right one (it could also be
      ! +pi or -pi off)
      !
      err=abs(phi_p1-amr_grid_xi(1,3))
      if(err.lt.1d-4) then 
         !
         ! Make sure it is exactly the right value
         !
         phi_p1 = amr_grid_xi(1,3)
      else
         ! 
         ! It is likely wrong, but if amr_grid_xi(1,3) is
         ! 0 or 2*pi it could be a 0 vs. 2*pi confusion...
         ! Since amr_grid_xi(1,3) can not be 2*pi we only
         ! check for 0
         !
         err=abs(phi_p1-twopi-amr_grid_xi(1,3))
         if(err.lt.1d-4) then
            !
            ! Make sure it is exactly the right value
            !
            phi_p1 = amr_grid_xi(1,3)
         else
            !
            ! It is the wrong phi, so invalidate
            !
            val_p1 = .false.
         endif
      endif
      !
      ! Check if r_p1 is on the grid
      !
      ! We use a tiny extra margin to prevent the ray from accidently
      ! "slipping through a corner"
      !
      if((r_p1.lt.amr_grid_xi(1,1)*oneminus).or.                    &
         (r_p1.gt.amr_grid_xi(amr_grid_nx+1,1)*oneplus)) then
         val_p1 = .false.
      endif
      !
      ! Check if theta_p1 is on grid
      !
      ! We use a tiny extra margin to prevent the ray from accidently
      ! "slipping through a corner"
      !
      if(amrray_mirror_equator) then
         if((theta_p1.le.amr_grid_xi(1,2)-tiny).or.&
              (theta_p1.ge.pi-amr_grid_xi(1,2)+tiny)) then
            val_p1 = .false.
         endif
      else
         if((theta_p1.le.amr_grid_xi(1,2)-tiny).or.&
              (theta_p1.ge.amr_grid_xi(amr_grid_ny+1,2)+tiny)) then
            val_p1 = .false.
         endif
      endif
   endif
   !
   ! -------------
   ! Now let us try to find the crossings with the phi2=const cone.
   ! -------------
   !
   if(val_p2) then
      !
      ! Find the distance to the crossing
      !
      dum1 = ray_cart_x*amrray_sinp2 - ray_cart_y*amrray_cosp2
      dum2 = ray_cart_diry*amrray_cosp2 - ray_cart_dirx*amrray_sinp2
      if(dum2.eq.0.d0) then
         val_p2 = .false.
         ds_p2 = -1d99
      else
         ds_p2 = dum1 / dum2
      endif
      if(ds_p2.le.0.d0) then
         val_p2 = .false.
      endif
      !
      ! Now create the new x,y,z coordinates
      !
      x_p2   = ray_cart_x + ds_p2 * ray_cart_dirx
      y_p2   = ray_cart_y + ds_p2 * ray_cart_diry
      z_p2   = ray_cart_z + ds_p2 * ray_cart_dirz
      !
      ! Find the r, theta, phi
      !
      r_p2     = sqrt( x_p2**2 + y_p2**2 + z_p2**2 )
      theta_p2 = acos(z_p2/r_p2)
      if(x_p2.eq.0.d0) then
         if(y_p2.ge.0.d0) then
            phi_p2 = pihalf
         else
            phi_p2 = 3.d0*pihalf 
         endif
      else
         phi_p2 = atan(y_p2/x_p2)
         if(x_p2.gt.0.d0) then
            if(y_p2.lt.0.d0) then
               phi_p2 = phi_p2 + twopi
            endif
         else
            phi_p2 = phi_p2 + pi
         endif
      endif
      if((phi_p2.lt.0.d0).or.(phi_p2.gt.twopi)) stop 1309
      !
      ! Check if phi_p2 is the right one (it could also be
      ! +pi or -pi off)
      !
      err=abs(phi_p2-amr_grid_xi(amr_grid_nz+1,3))
      if(err.lt.1d-4) then 
         !
         ! Make sure it is exactly the right value
         !
         phi_p2 = amr_grid_xi(amr_grid_nz+1,3)
      else
         ! 
         ! It is likely wrong, but if amr_grid_xi(1,3) is
         ! 0 or 2*pi it could be a 0 vs. 2*pi confusion...
         ! Since amr_grid_xi(amr_grid_nz+1,3) can not be 0 we only
         ! check for 2*pi
         !
         err=abs(phi_p2+twopi-amr_grid_xi(amr_grid_nz+1,3))
         if(err.lt.1d-4) then
            !
            ! Make sure it is exactly the right value
            !
            phi_p2 = amr_grid_xi(amr_grid_nz+1,3)
         else
            !
            ! It is the wrong phi, so invalidate
            !
            val_p2 = .false.
         endif
      endif
      !
      ! Check if r_p2 is on the grid
      !
      ! We use a tiny extra margin to prevent the ray from accidently
      ! "slipping through a corner"
      !
      if((r_p2.lt.amr_grid_xi(1,1)*oneminus).or.                    &
         (r_p2.gt.amr_grid_xi(amr_grid_nx+1,1)*oneplus)) then
         val_p2 = .false.
      endif
      !
      ! Check if theta_p2 is on grid
      !
      ! We use a tiny extra margin to prevent the ray from accidently
      ! "slipping through a corner"
      !
      if(amrray_mirror_equator) then
         if((theta_p2.le.amr_grid_xi(1,2)-tiny).or.&
            (theta_p2.ge.pi-amr_grid_xi(1,2)+tiny)) then
            val_p2 = .false.
         endif
      else
         if((theta_p2.le.amr_grid_xi(1,2)-tiny).or.&
              (theta_p2.ge.amr_grid_xi(amr_grid_ny+1,2)+tiny)) then
            val_p2 = .false.
         endif
      endif
   endif
   !
   ! -------------------------------------------------------------------
   ! Now determine which, if any, of the crossings is the first one
   ! that will be encountered.
   ! -------------------------------------------------------------------
   !
   ! Determine if/where the ray would cross the grid
   !
   amrray_icross   = 0
   cross_ds = 1d99
   if(val_r1) then
      amrray_icross   = 1
      cross_ds = ds_r1
   endif
   if(val_r2.and.(ds_r2.lt.cross_ds)) then
      amrray_icross   = 2
      cross_ds = ds_r2
   endif
   if(val_t1.and.(ds_t1.lt.cross_ds)) then
      amrray_icross   = 3
      cross_ds = ds_t1
   endif
   if(val_t2.and.(ds_t2.lt.cross_ds)) then
      amrray_icross   = 4
      cross_ds = ds_t2
   endif
   if(val_p1.and.(ds_p1.lt.cross_ds)) then
      amrray_icross   = 5
      cross_ds = ds_p1
   endif
   if(val_p2.and.(ds_p2.lt.cross_ds)) then
      amrray_icross   = 6
      cross_ds = ds_p2
   endif
   !
   ! If we have spheres, in addition to the grid, then we must check
   ! out intersections with those.
   !
   if(allocated(amrray_spheres_r)) then
      ! amrray_ispherehit = 0     ! Already done
      do isph=1,amrray_spheres_nr
         if(amrray_spheres_outsidegrid(isph)) then
            !
            ! This sphere is (at least partly) outside of the grid. So
            ! we must test whether we intersect with it.
            !
            call amrray_find_crossing_with_sphere(          &
                 ray_cart_x,ray_cart_y,ray_cart_z,          &
                 ray_cart_dirx,ray_cart_diry,ray_cart_dirz, &
                 isph,ds_sphere0,ds_sphere,ihit,            &
                 amrray_mirror_equator)
            !
            ! If we have an intersection:
            !
            if(ihit.ne.0) then
!!#############################################################
!               write(stdo,*) '!!! ',ihit,isph,ds_sphere
!!#############################################################
               amrray_ispherehit = isph*ihit
               ds_sphere0 = ds_sphere
            endif
         endif
      enddo
      !
      ! Check if this is before ray_dsend or cross_ds
      !
      if(ds_sphere0.gt.min(cross_ds,ray_dsend)) then
         !
         ! No it is after, so the sphere is not (yet) hit. So reset
         ! things. 
         !
         amrray_ispherehit = 0
         ds_sphere0 = 1.d99
      else
         !
         ! Yes it is before, so we must tell the code that the
         ! next crossing is not a cell
         !
         amrray_icross = 0
      endif
   endif
!!################################
!   write(stdo,*) '&2 ',amrray_ispherehit,ray_cart_z
!!################################
   !
   ! If amrray_icross==0 then we do not enter the grid at all, so we have
   ! arrived
   !
   if(amrray_icross.eq.0) then
      !
      ! We do not enter the grid at all
      !
      ! Check if we got stuck at a sphere
      !
      if(amrray_ispherehit.ne.0) then
         !
         ! Yes, hit a sphere, so we update the position and exit
         !
         ray_cart_x = ray_cart_x + ds_sphere0 * ray_cart_dirx
         ray_cart_y = ray_cart_y + ds_sphere0 * ray_cart_diry
         ray_cart_z = ray_cart_z + ds_sphere0 * ray_cart_dirz
         if(amrray_mirror_equator) then
            if(ray_cart_z.eq.0.d0) then
               ray_cart_dirz = abs(ray_cart_dirz)
            elseif(ray_cart_z.lt.0.d0) then
               ray_cart_z    = -ray_cart_z
               ray_cart_dirz = -ray_cart_dirz
            endif
         endif
      else
         !
         ! No sphere hit, just exit and say that we arrived.
         !
         arrived = .true.
      endif
      !
      ! Now return
      !
      nullify(amrray_nextcell)
      ray_indexnext = 0
      return
   endif
   !
   ! We could enter the grid, but we have to check if the cross_ds is
   ! smaller than ray_dsend 
   !
   if(cross_ds.gt.ray_dsend) then
      !
      ! Nope, we arrived at the end point
      !
      ray_cart_x = ray_cart_x + ray_dsend * ray_cart_dirx
      ray_cart_y = ray_cart_y + ray_dsend * ray_cart_diry
      ray_cart_z = ray_cart_z + ray_dsend * ray_cart_dirz
      if(amrray_mirror_equator) then
         if(ray_cart_z.eq.0.d0) then
            ray_cart_dirz = abs(ray_cart_dirz)
         elseif(ray_cart_z.lt.0.d0) then
            ray_cart_z    = -ray_cart_z
            ray_cart_dirz = -ray_cart_dirz
         endif
      endif
      nullify(amrray_nextcell)
      ray_indexnext = 0
      if(ray_dsend.lt.ds_sphere0) then
         arrived = .true.
      endif
      return
   endif
   !
   ! We are entering the grid. Check out where, and update things.
   !
   if(amrray_icross.eq.1) then
      !
      ! Entering the grid from inner radius. 
      !
      ray_cart_x = x_r1
      ray_cart_y = y_r1
      ray_cart_z = z_r1
      !
      ! If mirror, then correct for this
      !
      if(amrray_mirror_equator) then
         if(ray_cart_z.eq.0.d0) then
            ray_cart_dirz = abs(ray_cart_dirz)
         elseif(ray_cart_z.lt.0.d0) then
            ray_cart_z    = -ray_cart_z
            ray_cart_dirz = -ray_cart_dirz
            theta_r1      = pi - theta_r1
         endif
      endif
      !
      ! The calculation of the crossing with the inner radius may suffer
      ! from linear precision loss if the starting point of the ray
      ! is several factors of 10 further out. While the s00 trick (see
      ! above) prevents quadratic loss of precision, which would easily
      ! lead to catastrophic problems, it cannot avoid linear loss of
      ! precision. While linear loss of precision is not so damaging,
      ! it might in very special circumstances cause problems if you have
      ! very finely space gridding in radius, for instance. So let us 
      ! correct for this. Note that we already checked beforehand that
      ! the error is not too large.
      !
      dummy      = r_r1/sqrt(ray_cart_x**2+ray_cart_y**2+ray_cart_z**2)
      ray_cart_x = ray_cart_x * dummy
      ray_cart_y = ray_cart_y * dummy
      ray_cart_z = ray_cart_z * dummy
      !
      ! Find cell
      !
      call amrhunt(amr_grid_xi(:,2),amr_grid_ny+1,theta_r1,iy)
      call amrhunt(amr_grid_xi(:,3),amr_grid_nz+1,phi_r1,iz)
      ix = 1
      if(iy.gt.amr_grid_ny) iy = amr_grid_ny
      if(iy.lt.1) iy = 1
      if(iz.gt.amr_grid_nz) then
         if(amr_cyclic_xyz(3)) then
            iz = 1
         else
            iz = amr_grid_nz
         endif
      endif
      if(iz.lt.1) then
         if(amr_cyclic_xyz(3)) then
            iz = amr_grid_nz
         else
            iz = 1
         endif
      endif
!      if(amrray_selfcheck) then
!         if(amr_tree_present) then
!            if(.not.associated(amr_grid_branch(ix,iy,iz)%link)) then
!               write(stdo,*) 'ERROR: Internal error'
!               stop 6204
!            endif
!         endif
!      endif
      !
      ! ...Store this information
      !
      if(amr_tree_present) then
         !
         ! If the AMR tree is present, then use amrray_nextcell to 
         ! store the information
         !
         amrray_nextcell => amr_grid_branch(ix,iy,iz)%link
         !
         ! Now, if the new cell is in fact a nodal branch with children, then
         ! we must find this sub-cell
         !
         if(.not.amrray_nextcell%leaf) then
            call amrray_find_subcell(amrray_nextcell,1,1,     &
                                     r_r1,theta_r1,phi_r1)
         endif
         !
         ! Finally, find index
         !
         ray_indexnext = amrray_nextcell%leafindex
         if(present(levelnext)) levelnext = amrray_nextcell%level
      else
         !
         ! If the grid is regular, then use amrray_ix_next, amrray_iy_next, amrray_iz_next
         ! to store the information
         !
         amrray_ix_next = ix
         amrray_iy_next = iy
         amrray_iz_next = iz
         !
         ! Finally, find index
         !
         ray_indexnext = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
         if(present(levelnext)) levelnext = 0
      endif
      arrived = .false.
      amrray_icross = -1
      return
      !
   elseif(amrray_icross.eq.2) then
      !
      ! Entering the grid from outer radius. 
      !
      ray_cart_x = x_r2
      ray_cart_y = y_r2
      ray_cart_z = z_r2
      !
      ! If mirror, then correct for this
      !
      if(amrray_mirror_equator) then
         if(ray_cart_z.eq.0.d0) then
            ray_cart_dirz = abs(ray_cart_dirz)
         elseif(ray_cart_z.lt.0.d0) then
            ray_cart_z    = -ray_cart_z
            ray_cart_dirz = -ray_cart_dirz
            theta_r2      = pi - theta_r2
         endif
      endif
      !
      ! Find cell
      !
      call amrhunt(amr_grid_xi(:,2),amr_grid_ny+1,theta_r2,iy)
      call amrhunt(amr_grid_xi(:,3),amr_grid_nz+1,phi_r2,iz)
      ix = amr_grid_nx
      if(iy.gt.amr_grid_ny) iy = amr_grid_ny
      if(iy.lt.1) iy = 1
      if(iz.gt.amr_grid_nz) then
         if(amr_cyclic_xyz(3)) then
            iz = 1
         else
            iz = amr_grid_nz
         endif
      endif
      if(iz.lt.1) then
         if(amr_cyclic_xyz(3)) then
            iz = amr_grid_nz
         else
            iz = 1
         endif
      endif
!      if(amrray_selfcheck) then
!         if(amr_tree_present) then
!            if(.not.associated(amr_grid_branch(ix,iy,iz)%link)) then
!               write(stdo,*) 'ERROR: Internal error'
!               stop 6204
!            endif
!         endif
!      endif
      !
      ! ...Store this information
      !
      if(amr_tree_present) then
         !
         ! If the AMR tree is present, then use amrray_nextcell to 
         ! store the information
         !
         amrray_nextcell => amr_grid_branch(ix,iy,iz)%link
         !
         ! Now, if the new cell is in fact a nodal branch with children, then
         ! we must find this sub-cell
         !
         if(.not.amrray_nextcell%leaf) then
            call amrray_find_subcell(amrray_nextcell,1,2,     &
                                     r_r2,theta_r2,phi_r2)
         endif
         !
         ! Finally, find index
         !
         ray_indexnext = amrray_nextcell%leafindex
         if(present(levelnext)) levelnext = amrray_nextcell%level
      else
         !
         ! If the grid is regular, then use amrray_ix_next, amrray_iy_next, amrray_iz_next
         ! to store the information
         !
         amrray_ix_next = ix
         amrray_iy_next = iy
         amrray_iz_next = iz
         !
         ! Finally, find index
         !
         ray_indexnext = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
         if(present(levelnext)) levelnext = 0
      endif
      arrived = .false.
      amrray_icross = -2
      return
      !
   elseif(amrray_icross.eq.3) then
      !
      ! Entering the grid from upper theta.
      !
      ray_cart_x = x_t1
      ray_cart_y = y_t1
      ray_cart_z = z_t1
      !
      ! If mirror, then correct for this
      !
!      if(amrray_mirror_equator) then
!         if(ray_cart_z.eq.0.d0) then
!            ray_cart_dirz = abs(ray_cart_dirz)
!         elseif(ray_cart_z.lt.0.d0) then
!            ray_cart_z    = -ray_cart_z
!            ray_cart_dirz = -ray_cart_dirz
!            theta_t1      = pi - theta_t1
!         endif
!      endif
      !
      ! Find cell
      !
      call amrhunt(amr_grid_xi(:,1),amr_grid_nx+1,r_t1,ix)
      call amrhunt(amr_grid_xi(:,3),amr_grid_nz+1,phi_t1,iz)
      iy = 1
      if(ix.lt.1) ix = 1
      if(ix.gt.amr_grid_nx) ix = amr_grid_nx
      if(iz.gt.amr_grid_nz) then
         if(amr_cyclic_xyz(3)) then
            iz = 1
         else
            iz = amr_grid_nz
         endif
      endif
      if(iz.lt.1) then
         if(amr_cyclic_xyz(3)) then
            iz = amr_grid_nz
         else
            iz = 1
         endif
      endif
!      if(amrray_selfcheck) then
!         if(amr_tree_present) then
!            if(.not.associated(amr_grid_branch(ix,iy,iz)%link)) then
!               write(stdo,*) 'ERROR: Internal error'
!               stop 6204
!            endif
!         endif
!      endif
      !
      ! ...Store this information
      !
      if(amr_tree_present) then
         !
         ! If the AMR tree is present, then use amrray_nextcell to 
         ! store the information
         !
         amrray_nextcell => amr_grid_branch(ix,iy,iz)%link
         !
         ! Now, if the new cell is in fact a nodal branch with children, then
         ! we must find this sub-cell
         !
         if(.not.amrray_nextcell%leaf) then
            call amrray_find_subcell(amrray_nextcell,2,1,     &
                                     r_t1,theta_t1,phi_t1)
         endif
         !
         ! Finally, find index
         !
         ray_indexnext = amrray_nextcell%leafindex
         if(present(levelnext)) levelnext = amrray_nextcell%level
      else
         !
         ! If the grid is regular, then use amrray_ix_next, amrray_iy_next, amrray_iz_next
         ! to store the information
         !
         amrray_ix_next = ix
         amrray_iy_next = iy
         amrray_iz_next = iz
         !
         ! Finally, find index
         !
         ray_indexnext = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
         if(present(levelnext)) levelnext = 0
      endif
      arrived = .false.
      amrray_icross = -3
      return
      !
   elseif(amrray_icross.eq.4) then
      !
      ! Entering the grid from lower theta.
      !
      ray_cart_x = x_t2
      ray_cart_y = y_t2
      ray_cart_z = z_t2
      !
      ! If mirror, then correct for this
      !
      if(amrray_mirror_equator) then
         if(ray_cart_z.eq.0.d0) then
            ray_cart_dirz = abs(ray_cart_dirz)
         elseif(ray_cart_z.lt.0.d0) then
            ray_cart_z    = -ray_cart_z
            ray_cart_dirz = -ray_cart_dirz
            theta_t2      = pi - theta_t2
         endif
      endif
      !
      ! Find cell
      !
      call amrhunt(amr_grid_xi(:,1),amr_grid_nx+1,r_t2,ix)
      call amrhunt(amr_grid_xi(:,3),amr_grid_nz+1,phi_t2,iz)
      iy = amr_grid_ny
      if(ix.lt.1) ix = 1
      if(ix.gt.amr_grid_nx) ix = amr_grid_nx
      if(iz.gt.amr_grid_nz) then
         if(amr_cyclic_xyz(3)) then
            iz = 1
         else
            iz = amr_grid_nz
         endif
      endif
      if(iz.lt.1) then
         if(amr_cyclic_xyz(3)) then
            iz = amr_grid_nz
         else
            iz = 1
         endif
      endif
!      if(amrray_selfcheck) then
!         if(amr_tree_present) then
!            if(.not.associated(amr_grid_branch(ix,iy,iz)%link)) then
!               write(stdo,*) 'ERROR: Internal error'
!               stop 6204
!            endif
!         endif
!      endif
      !
      ! ...Store this information
      !
      if(amr_tree_present) then
         !
         ! If the AMR tree is present, then use amrray_nextcell to 
         ! store the information
         !
         amrray_nextcell => amr_grid_branch(ix,iy,iz)%link
         !
         ! Now, if the new cell is in fact a nodal branch with children, then
         ! we must find this sub-cell
         !
         if(.not.amrray_nextcell%leaf) then
            call amrray_find_subcell(amrray_nextcell,2,2,     &
                                     r_t2,theta_t2,phi_t2)
         endif
         !
         ! Finally, find index
         !
         ray_indexnext = amrray_nextcell%leafindex
         if(present(levelnext)) levelnext = amrray_nextcell%level
      else
         !
         ! If the grid is regular, then use amrray_ix_next, amrray_iy_next, amrray_iz_next
         ! to store the information
         !
         amrray_ix_next = ix
         amrray_iy_next = iy
         amrray_iz_next = iz
         !
         ! Finally, find index
         !
         ray_indexnext = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
         if(present(levelnext)) levelnext = 0
      endif
      arrived = .false.
      amrray_icross = -4
      return
      !
   elseif(amrray_icross.eq.5) then
      !
      ! Entering the grid from smallest phi.
      !
      ray_cart_x = x_p1
      ray_cart_y = y_p1
      ray_cart_z = z_p1
      !
      ! If mirror, then correct for this
      !
      if(amrray_mirror_equator) then
         if(ray_cart_z.eq.0.d0) then
            ray_cart_dirz = abs(ray_cart_dirz)
         elseif(ray_cart_z.lt.0.d0) then
            ray_cart_z    = -ray_cart_z
            ray_cart_dirz = -ray_cart_dirz
            theta_p1      = pi - theta_p1
         endif
      endif
      !
      ! Find cell
      !
      call amrhunt(amr_grid_xi(:,1),amr_grid_nx+1,r_p1,ix)
      call amrhunt(amr_grid_xi(:,2),amr_grid_ny+1,theta_p1,iy)
      iz = 1
      if(ix.lt.1) ix = 1
      if(ix.gt.amr_grid_nx) ix = amr_grid_nx
      if(iy.lt.1) iy = 1
      if(iy.gt.amr_grid_ny) iy = amr_grid_ny
!      if(amrray_selfcheck) then
!         if(amr_tree_present) then
!            if(.not.associated(amr_grid_branch(ix,iy,iz)%link)) then
!               write(stdo,*) 'ERROR: Internal error'
!               stop 6204
!            endif
!         endif
!      endif
      !
      ! ...Store this information
      !
      if(amr_tree_present) then
         !
         ! If the AMR tree is present, then use amrray_nextcell to 
         ! store the information
         !
         amrray_nextcell => amr_grid_branch(ix,iy,iz)%link
         !
         ! Now, if the new cell is in fact a nodal branch with children, then
         ! we must find this sub-cell
         !
         if(.not.amrray_nextcell%leaf) then
            call amrray_find_subcell(amrray_nextcell,3,1,     &
                                     r_p1,theta_p1,phi_p1)
         endif
         !
         ! Finally, find index
         !
         ray_indexnext = amrray_nextcell%leafindex
         if(present(levelnext)) levelnext = amrray_nextcell%level
      else
         !
         ! If the grid is regular, then use amrray_ix_next, amrray_iy_next, amrray_iz_next
         ! to store the information
         !
         amrray_ix_next = ix
         amrray_iy_next = iy
         amrray_iz_next = iz
         !
         ! Finally, find index
         !
         ray_indexnext = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
         if(present(levelnext)) levelnext = 0
      endif
      arrived = .false.
      amrray_icross = -5
      return
      !
   elseif(amrray_icross.eq.6) then
      !
      ! Entering the grid from largest phi.
      !
      ray_cart_x = x_p2
      ray_cart_y = y_p2
      ray_cart_z = z_p2
      !
      ! If mirror, then correct for this
      !
      if(amrray_mirror_equator) then
         if(ray_cart_z.eq.0.d0) then
            ray_cart_dirz = abs(ray_cart_dirz)
         elseif(ray_cart_z.lt.0.d0) then
            ray_cart_z    = -ray_cart_z
            ray_cart_dirz = -ray_cart_dirz
            theta_p2      = pi - theta_p2
         endif
      endif
      !
      ! Find cell
      !
      call amrhunt(amr_grid_xi(:,1),amr_grid_nx+1,r_p2,ix)
      call amrhunt(amr_grid_xi(:,2),amr_grid_ny+1,theta_p2,iy)
      iz = amr_grid_nz
      if(ix.lt.1) ix = 1
      if(ix.gt.amr_grid_nx) ix = amr_grid_nx
      if(iy.lt.1) iy = 1
      if(iy.gt.amr_grid_ny) iy = amr_grid_ny
!      if(amrray_selfcheck) then
!         if(amr_tree_present) then
!            if(.not.associated(amr_grid_branch(ix,iy,iz)%link)) then
!               write(stdo,*) 'ERROR: Internal error'
!               stop 6204
!            endif
!         endif
!      endif
      !
      ! ...Store this information
      !
      if(amr_tree_present) then
         !
         ! If the AMR tree is present, then use amrray_nextcell to 
         ! store the information
         !
         amrray_nextcell => amr_grid_branch(ix,iy,iz)%link
         !
         ! Now, if the new cell is in fact a nodal branch with children, then
         ! we must find this sub-cell
         !
         if(.not.amrray_nextcell%leaf) then
            call amrray_find_subcell(amrray_nextcell,3,2,     &
                                     r_p2,theta_p2,phi_p2)
         endif
         !
         ! Finally, find index
         !
         ray_indexnext = amrray_nextcell%leafindex
         if(present(levelnext)) levelnext = amrray_nextcell%level
      else
         !
         ! If the grid is regular, then use amrray_ix_next, amrray_iy_next, amrray_iz_next
         ! to store the information
         !
         amrray_ix_next = ix
         amrray_iy_next = iy
         amrray_iz_next = iz
         !
         ! Finally, find index
         !
         ray_indexnext = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
         if(present(levelnext)) levelnext = 0
      endif
      arrived = .false.
      amrray_icross = -6
      return
      !
   endif
   !
else
   !
   !========================================================================
   ! We are inside the grid, in a grid cell, so we have to find the crossing
   ! with cell edges.
   !========================================================================
   !
   ! Self-check
   !
   if(amrray_selfcheck) then
      if((r0>axi(2,1)*oneplust).or.               &
         (r0<axi(1,1)*oneminust).or.              &
         (theta0>axi(2,2)+tol).or.              &
         (theta0<axi(1,2)-tol).or.              &
         (phi0>axi(2,3)+tol).or.                &
         (phi0<axi(1,3)-tol)) then
         write(stdo,*) 'ERROR: Photon outside of cell'
         write(stdo,*) ray_cart_x,ray_cart_y,ray_cart_z
         write(stdo,*) ray_cart_dirx,ray_cart_diry,ray_cart_dirz
         write(stdo,*) r0,theta0,phi0
         if(amr_tree_present) write(stdo,*) associated(amrray_cell),amrray_cell%id
         write(stdo,*) axi
         stop 1131
      endif
   endif
   !
   ! Now do some preliminary calculations for computing the crossings
   ! with the cell walls
   !
   x0dotdir  = ray_cart_x * ray_cart_dirx +              &
               ray_cart_y * ray_cart_diry +              &
               ray_cart_z * ray_cart_dirz
   !
   ! Get the current branch index
   !
   !idx   = amrray_cell%branchindex
   !
   ! Reset some stuff
   !
   amrray_icross   = 0
   cross_ds = ray_dsend
   !
   ! Now try out crossings with all 6 cell walls. Since we are starting
   ! _within_ a cell the true crossing is the one with the smallest ds.
   ! We do not need to check if the crossing happens within the limits
   ! of the cell wall, as we had to do above.
   ! That is why this is much easier than the determination of the
   ! crossings from outside the grid where this simple selection criterion
   ! could not be used.
   !
   ! Try out a crossing with inner radius cell wall.
   ! Only if ray is moving inward. It is always the smaller of the two.
   !
   if(x0dotdir.lt.0.d0) then
      det       = x0dotdir*x0dotdir + axi(1,1)**2 - r02
      if(det.gt.0.d0) then
         ds_try = -x0dotdir - sqrt(det)
         if(ds_try.lt.0.d0) ds_try = 0.d0   ! Could happen due to round-off errors
         if(ds_try.lt.cross_ds) then
            amrray_icross = 1
            cross_ds = ds_try
         endif
      endif
   endif
   !
   ! Try out crossing with outer radius.
   ! Since we are within the cell the solution must be
   ! the largest of the two.
   !
   det       = x0dotdir*x0dotdir + axi(2,1)**2 - r02
   ds_try    = -x0dotdir + sqrt(det)
   if(ds_try.lt.0.d0) ds_try = 0.d0    ! Could happen due to round-off errors
   if(ds_try.lt.cross_ds) then
      amrray_icross = 2
      cross_ds = ds_try
   endif
   ! 
   ! Check if the theta dimension is active (if not, then we are doing
   ! 1-D spherical radiative transfer)
   !
   if(amr_ydim.ne.0) then
      !
      ! First check if we are in the top or bottom quadrant
      !
      if(axi(1,2).lt.pihalf) then
         topquadrant = .true.
         sgnz        = 1.d0
         if(ray_cart_z.lt.0.d0) stop 3331
      else
         topquadrant = .false.
         sgnz        = -1.d0
         if(ray_cart_z.gt.0.d0) then
            write(stdo,*) 'ERROR in spherical coordinates'
            write(stdo,*) ray_cart_z,axi(1,2),r0,theta0,phi0
            stop 3332
         endif
      endif
      !
      ! Check if the cell touches the equatorial plane
      !
      if(amr_tree_present) then
         !
         ! If the AMR tree is present, then we use that
         !
         cossindum = amrray_finegrid_costsq2(amrray_cell%ixyzf(2),amrray_cell%level)
      else
         !
         ! If the grid is regular, we use amrray_ix_curr, amrray_iy_curr, amrray_iz_curr
         !
         cossindum = amrray_finegrid_costsq2(amrray_iy_curr,0)
      endif
      if(cossindum.eq.0.d0) then
         !
         ! This is simple: just the crossing with the z=0-plane
         !
         if(ray_cart_dirz.ne.0.d0) then
            ds_try = -ray_cart_z / ray_cart_dirz
            if((ds_try.gt.0.d0).and.(ds_try.lt.cross_ds)) then
               !
               ! Yes, we have a valid crossing
               !
               if(ray_cart_z.lt.0.d0) then
                  !
                  ! We cross from below
                  !
                  amrray_icross   = 3
                  cross_ds = ds_try
                  crossequator = .true.
               else
                  !
                  ! We cross from above
                  !
                  amrray_icross   = 4
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
      if(amr_tree_present) then
         !
         ! If the AMR tree is present, then use that
         !
         st2   = amrray_finegrid_sintsq1(amrray_cell%ixyzf(2),amrray_cell%level)
         ct2   = amrray_finegrid_costsq1(amrray_cell%ixyzf(2),amrray_cell%level)
      else
         !
         ! If the grid is regular, then use amrray_ix_curr, amrray_iy_curr, amrray_iz_curr
         !
         st2   = amrray_finegrid_sintsq1(amrray_iy_curr,0)
         ct2   = amrray_finegrid_costsq1(amrray_iy_curr,0)
      endif
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
         pa    = ct2*(ray_cart_dirx**2+ray_cart_diry**2)-st2*ray_cart_dirz**2
         pb    = 2.d0*ct2*(ray_cart_dirx*ray_cart_x + ray_cart_diry*ray_cart_y) - &
                 2.d0*st2*ray_cart_dirz*ray_cart_z
         pc    = ct2*(ray_cart_x**2+ray_cart_y**2)-st2*ray_cart_z**2
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
               !  1-sqrt(1-eps) = (1/2)*x + (1/8)*x^2 + (1/16)*x^3
               !                  + (5/128)*x^4 + ...
               !
               ds_try = ( pc / abs(pb) ) * ( 1.d0 + 0.25d0*eps + &
                          0.125d0*eps**2 + (5.d0/64.d0)*eps**3 )
               if((ds_try.gt.0.d0).and.(ds_try.lt.cross_ds)) then
                  if(topquadrant) then
                     amrray_icross   = 3
                  else
                     amrray_icross   = 4
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
                        amrray_icross   = 3
                     else
                        amrray_icross   = 4
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
      if(amr_tree_present) then
         !
         ! If the AMR tree is present, then use that
         !
         st2   = amrray_finegrid_sintsq2(amrray_cell%ixyzf(2),amrray_cell%level)
         ct2   = amrray_finegrid_costsq2(amrray_cell%ixyzf(2),amrray_cell%level)
      else
         !
         ! If the grid is regular, then use amrray_ix_curr, amrray_iy_curr, amrray_iz_curr
         !
         st2   = amrray_finegrid_sintsq2(amrray_iy_curr,0)
         ct2   = amrray_finegrid_costsq2(amrray_iy_curr,0)
      endif
      !
      ! But do this check only if it has not yet been done before for
      ! the equator crossing
      !
      if(ct2.ne.0.d0) then
         !
         ! Compute some of the coefficients
         !
         pa    = ct2*(ray_cart_dirx**2+ray_cart_diry**2)-st2*ray_cart_dirz**2
         pb    = 2.d0*ct2*(ray_cart_dirx*ray_cart_x + ray_cart_diry*ray_cart_y) - &
                 2.d0*st2*ray_cart_dirz*ray_cart_z
         pc    = ct2*(ray_cart_x**2+ray_cart_y**2)-st2*ray_cart_z**2
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
               !  1-sqrt(1-eps) = (1/2)*x + (1/8)*x^2 + (1/16)*x^3
               !                  + (5/128)*x^4 + ...
               !
               ds_try = - ( pc / pb ) * ( 1.d0 + 0.25d0*eps +   &
                         0.125d0*eps**2 + (5.d0/64.d0)*eps**3 )
               if((ds_try.gt.0.d0).and.(ds_try.lt.cross_ds)) then
                  if(topquadrant) then
                     amrray_icross   = 4
                  else
                     amrray_icross   = 3
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
                        amrray_icross   = 4
                     else
                        amrray_icross   = 3
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
   if(amr_zdim.ne.0) then
      !
      ! Find the distance to the crossing with the smallest phi
      !
      if(amr_tree_present) then
         !
         ! If the AMR tree is present, then use that
         !
         dum1 = ray_cart_x*amrray_finegrid_sinp1(amrray_cell%ixyzf(3),amrray_cell%level) - &
                ray_cart_y*amrray_finegrid_cosp1(amrray_cell%ixyzf(3),amrray_cell%level)
         dum2 = ray_cart_diry*amrray_finegrid_cosp1(amrray_cell%ixyzf(3),amrray_cell%level) - &
                ray_cart_dirx*amrray_finegrid_sinp1(amrray_cell%ixyzf(3),amrray_cell%level)
      else
         !
         ! If the grid is regular, then use amrray_ix_curr, amrray_iy_curr, amrray_iz_curr
         !
         dum1 = ray_cart_x*amrray_finegrid_sinp1(amrray_iz_curr,0) - &
                ray_cart_y*amrray_finegrid_cosp1(amrray_iz_curr,0)
         dum2 = ray_cart_diry*amrray_finegrid_cosp1(amrray_iz_curr,0) - &
                ray_cart_dirx*amrray_finegrid_sinp1(amrray_iz_curr,0)
      endif
      if(dum2.lt.-small) then   ! Thereby selecting also only outgoing rays
         ds_try = dum1 / dum2
      else
         ds_try = -1.d0
      endif
      if((ds_try.ge.0.d0).and.(ds_try.lt.cross_ds)) then
         amrray_icross   = 5
         cross_ds = ds_try
      endif
      !
      ! Find the distance to the crossing with the largest phi
      !
      if(amr_tree_present) then
         !
         ! If the AMR tree is present, then use that
         !
         dum1 = ray_cart_x*amrray_finegrid_sinp2(amrray_cell%ixyzf(3),amrray_cell%level) - &
                ray_cart_y*amrray_finegrid_cosp2(amrray_cell%ixyzf(3),amrray_cell%level)
         dum2 = ray_cart_diry*amrray_finegrid_cosp2(amrray_cell%ixyzf(3),amrray_cell%level) - &
                ray_cart_dirx*amrray_finegrid_sinp2(amrray_cell%ixyzf(3),amrray_cell%level)
      else
         !
         ! If the grid is regular, then use amrray_ix_curr, amrray_iy_curr, amrray_iz_curr
         !
         dum1 = ray_cart_x*amrray_finegrid_sinp2(amrray_iz_curr,0) - &
                ray_cart_y*amrray_finegrid_cosp2(amrray_iz_curr,0)
         dum2 = ray_cart_diry*amrray_finegrid_cosp2(amrray_iz_curr,0) - &
                ray_cart_dirx*amrray_finegrid_sinp2(amrray_iz_curr,0)
      endif
      if(dum2.gt.small) then    ! Thereby selecting also only outgoing rays
         ds_try = dum1 / dum2
      else
         ds_try = -1.d0
      endif
      if((ds_try.ge.0.d0).and.(ds_try.lt.cross_ds)) then
         amrray_icross   = 6
         cross_ds = ds_try
      endif
   endif
   !
   ! If a maximum delta angle is set, then we must assure that the
   ! next point does not have a larger angle (wrt to the origin) with
   ! the start of this segment than this maximum delta angle. This is
   ! important for e.g. line transfer where we must assure that the
   ! line does not doppler shift too much within a single cell, and
   ! for the 2-D scattering mode.
   !
   if(present(maxdeltasina)) then
      !
      ! Find the cross product between vec x and vec x + cross_dr * vec dir
      ! Since the cross product of vec x with itself is zero, we can
      ! also calculate the cross product between vec x and the 
      ! vector cross_dr * vec dir. 
      !
      angvec(1) = ray_cart_y*ray_cart_dirz - ray_cart_z*ray_cart_diry
      angvec(2) = ray_cart_z*ray_cart_dirx - ray_cart_x*ray_cart_dirz
      angvec(3) = ray_cart_x*ray_cart_diry - ray_cart_y*ray_cart_dirx
      !
      ! Normalize.
      !
      ! Note: actually it should be divided by |r|*|r+dir*cross_ds|
      ! but this is accurate enough for the present purpose.
      ! It remains an estimation, though
      !
      angvec(:) = angvec(:) * cross_ds / ( r0 * ( r0 + cross_ds ) )
      !
      ! Calculate the sin(angle)
      !
      sinang = sqrt(angvec(1)**2+angvec(2)**2+angvec(3)**2)
      !
      ! Now check
      !
      if(sinang.gt.maxdeltasina) then
         !
         ! So indeed the angle threatens to become too large, so
         ! decrease the step
         !
         cross_ds = cross_ds * maxdeltasina / sinang
         !
         ! And tell the rest of the routine that we are not 
         ! crossing a cell wall after all
         !
         amrray_icross   = 0
         hitmaxphi       = .true.
      endif
   endif
   !
   ! Now check for the spheres
   !
   if(allocated(amrray_spheres_sphidx)) then
      ! amrray_ispherehit = 0     ! Already done
      if(ray_indexcurr.le.0) stop 7390
      if(amrray_spheres_sphidx(ray_indexcurr).gt.0) then
         !
         ! There is a single sphere inside or partly inside this cell
         !
         isph = amrray_spheres_sphidx(ray_indexcurr)
         if(isph.gt.amrray_spheres_nr) stop 7391
         !
         ! Now check the crossing with this sphere
         !
         call amrray_find_crossing_with_sphere(             &
                 ray_cart_x,ray_cart_y,ray_cart_z,          &
                 ray_cart_dirx,ray_cart_diry,ray_cart_dirz, &
                 isph,cross_ds,ds_sphere,ihit,.false.)
         !
         ! If we have an intersection:
         !
         if(ihit.ne.0) then
            amrray_ispherehit = isph*ihit
            cross_ds = ds_sphere
            amrray_icross = 0
         endif
         !
      elseif(amrray_spheres_sphidx(ray_indexcurr).lt.0) then
         !
         ! There is more than one sphere in this cell
         ! We must check all spheres (sorry, no cleverer method yet)
         !
         ! Do a loop over spheres
         !
         do isph=1,amrray_spheres_nr
            !
            ! Now check the crossing with this sphere
            !
            call amrray_find_crossing_with_sphere(          &
                 ray_cart_x,ray_cart_y,ray_cart_z,          &
                 ray_cart_dirx,ray_cart_diry,ray_cart_dirz, &
                 isph,cross_ds,ds_sphere,ihit,.false.)
            !
            ! If we have an intersection:
            !
            if(ihit.ne.0) then
               amrray_ispherehit = isph*ihit
               cross_ds = ds_sphere
               amrray_icross = 0
            endif
         enddo
         !
      endif
   endif
   !
   ! OK, now we know which crossing takes place.
   !
   ! Compute the new position
   !
   ray_cart_x = ray_cart_x + cross_ds * ray_cart_dirx
   ray_cart_y = ray_cart_y + cross_ds * ray_cart_diry
   ray_cart_z = ray_cart_z + cross_ds * ray_cart_dirz
   !
   ! Some stuff for when the theta-direction is active and equatorial
   ! plane is crossed or approached
   !
   if(amr_ydim.ne.0) then
      !
      ! Since we know that there is no cell that crosses the equatorial
      ! plane, we must make sure that this does not happen to ray_cart_z
      ! due to round-off errors.
      !
      if(topquadrant) then
         if(ray_cart_z.lt.0.d0) ray_cart_z=0.d0
      else
         if(ray_cart_z.gt.0.d0) ray_cart_z=0.d0
      endif
      !
      ! Check if equator crossed
      !
      if(crossequator) then
         if(topquadrant) then
            if(amrray_icross.ne.4) crossequator=.false.
         else
            if(amrray_icross.ne.3) crossequator=.false.
         endif
         !
         ! If we know that we cross the equator, we put ray_cart_z to 0.d0
         ! to avoid a not-exactly-zero z value.
         !
         if(crossequator) then
            ray_cart_z=0.d0
            !
            ! If we crossed the equator with mirroring on, then we 
            ! must also flip the ray_cart_dirz.
            !
            if(amrray_mirror_equator) ray_cart_dirz=abs(ray_cart_dirz)
         endif
      endif
   endif
   !
   ! Find the cell to which we move.
   !
   select case(amrray_icross) 
   case(0)
      !
      ! We did not reach a cell wall before reaching the end of our
      ! journey
      !
      ! **** QUESTION: Why is ray_indexnext set here? It is anyway set lateron ****
      !
      if(amr_tree_present) then
         !
         ! The AMR tree is used, so connect amrray_nextcell pointer
         !
         amrray_nextcell => amrray_cell
         ray_indexnext = amrray_nextcell%leafindex
      else
         !
         ! A regular grid is used, so set the ixyz_next indices
         !
         amrray_ix_next = amrray_ix_curr
         amrray_iy_next = amrray_iy_curr
         amrray_iz_next = amrray_iz_curr
         ray_indexnext = amrray_ix_next+((amrray_iy_next-1)+(amrray_iz_next-1)*amr_grid_ny)*amr_grid_nx
      endif
      !
      ! If we have amrray_icross==0 without amrray_ispherehit>0, then we have
      ! arrived
      !
      if((amrray_ispherehit.eq.0).and.(.not.hitmaxphi)) arrived=.true.
      return
   case(1)
      idir = 1
      ilr  = 2
      if(amr_tree_present) then
         !
         ! The AMR tree is used, so connect amrray_nextcell pointer.
         ! If we get outside of the domain, this link is automatically
         ! null.
         !
         amrray_nextcell => amrray_cell%neighbor(1,1)%link
      else
         !
         ! A regular grid is used, so set the ixyz_next indices
         !
         amrray_ix_next = amrray_ix_curr-1
         amrray_iy_next = amrray_iy_curr
         amrray_iz_next = amrray_iz_curr
         !
         ! Check if we are outside of the domain
         !
         if(amrray_ix_next.lt.1) then
            amrray_ix_next = -1
            amrray_iy_next = -1
            amrray_iz_next = -1
         endif
      endif
   case(2)
      idir = 1
      ilr  = 1
      if(amr_tree_present) then
         !
         ! The AMR tree is used, so connect amrray_nextcell pointer.
         ! If we get outside of the domain, this link is automatically
         ! null.
         !
         amrray_nextcell => amrray_cell%neighbor(2,1)%link
      else
         !
         ! A regular grid is used, so set the ixyz_next indices
         !
         amrray_ix_next = amrray_ix_curr+1
         amrray_iy_next = amrray_iy_curr
         amrray_iz_next = amrray_iz_curr
         !
         ! Check if we are outside of the domain
         !
         if(amrray_ix_next.gt.amr_grid_nx) then
            amrray_ix_next = -1
            amrray_iy_next = -1
            amrray_iz_next = -1
         endif
      endif
   case(3)
      idir = 2
      ilr  = 2
      if(amr_tree_present) then
         !
         ! The AMR tree is used, so connect amrray_nextcell pointer.
         ! If we get outside of the domain, this link is automatically
         ! null.
         !
         amrray_nextcell  => amrray_cell%neighbor(1,2)%link
      else
         !
         ! A regular grid is used, so set the ixyz_next indices
         !
         amrray_ix_next = amrray_ix_curr
         amrray_iy_next = amrray_iy_curr-1
         amrray_iz_next = amrray_iz_curr
         !
         ! Check if we are outside of the domain
         !
         if(amrray_iy_next.lt.1) then
            amrray_ix_next = -1
            amrray_iy_next = -1
            amrray_iz_next = -1
         endif
      endif
   case(4)
      idir = 2
      ilr  = 1
      if(crossequator.and.amrray_mirror_equator) then
         !
         ! The equator is crossed, and we have mirror symmetry.
         ! So we end up in the same cell, but now move back
         ! upward.
         !
         if(amr_tree_present) then
            !
            ! The AMR tree is used, so connect amrray_nextcell pointer.
            !
            amrray_nextcell  => amrray_cell
         else
            !
            ! A regular grid is used, so set the ixyz_next indices
            !
            amrray_ix_next = amrray_ix_curr
            amrray_iy_next = amrray_iy_curr
            amrray_iz_next = amrray_iz_curr
         endif
      else
         !
         ! Normal case: move toward larger theta
         !
         if(amr_tree_present) then
            !
            ! The AMR tree is used, so connect amrray_nextcell pointer.
            ! If we get outside of the domain, this link is automatically
            ! null.
            !
            amrray_nextcell => amrray_cell%neighbor(2,2)%link
         else
            !
            ! A regular grid is used, so set the ixyz_next indices
            !
            amrray_ix_next = amrray_ix_curr
            amrray_iy_next = amrray_iy_curr+1
            amrray_iz_next = amrray_iz_curr
            !
            ! Check if we are outside of the domain
            !
            if(amrray_iy_next.gt.amr_grid_ny) then
               amrray_ix_next = -1
               amrray_iy_next = -1
               amrray_iz_next = -1
            endif
         endif
      endif
   case(5)
      idir = 3
      ilr  = 2
      if(amr_tree_present) then
         !
         ! The AMR tree is used, so connect amrray_nextcell pointer.
         ! If we get outside of the domain, this link is automatically
         ! null.
         !
         amrray_nextcell => amrray_cell%neighbor(1,3)%link
      else
         !
         ! A regular grid is used, so set the ixyz_next indices
         !
         amrray_ix_next = amrray_ix_curr
         amrray_iy_next = amrray_iy_curr
         amrray_iz_next = amrray_iz_curr-1
         !
         ! Check if we are outside of the domain
         !
         if(amrray_iz_next.lt.1) then
            !
            ! If we have cyclic boundaries, then connect to
            ! the other side. Otherwise invalidate.
            !
            if(amr_cyclic_xyz(3)) then
               amrray_iz_next = amr_grid_nz
            else
               amrray_ix_next = -1
               amrray_iy_next = -1
               amrray_iz_next = -1
            endif
         endif
      endif
   case(6)
      idir = 3
      ilr  = 1
      if(amr_tree_present) then
         !
         ! The AMR tree is used, so connect amrray_nextcell pointer.
         ! If we get outside of the domain, this link is automatically
         ! null.
         !
         amrray_nextcell => amrray_cell%neighbor(2,3)%link
      else
         !
         ! A regular grid is used, so set the ixyz_next indices
         !
         amrray_ix_next = amrray_ix_curr
         amrray_iy_next = amrray_iy_curr
         amrray_iz_next = amrray_iz_curr+1
         !
         ! Check if we are outside of the domain
         !
         if(amrray_iz_next.gt.amr_grid_nz) then
            !
            ! If we have cyclic boundaries, then connect to
            ! the other side. Otherwise invalidate.
            !
            if(amr_cyclic_xyz(3)) then
               amrray_iz_next = 1
            else
               amrray_ix_next = -1
               amrray_iy_next = -1
               amrray_iz_next = -1
            endif
         endif
      endif
   end select
   !
   ! If we enter vacuum then we need to shift a bit away from the grid in
   ! order not to confuse the algorithm in the next round.
   !
   if(amr_tree_present) then
      if(.not.associated(amrray_nextcell)) then
         doshift = .true.
      else
         doshift = .false.
      endif
   else
      if(amrray_ix_next.le.0) then
         doshift = .true.
      else
         doshift = .false.
      endif
   endif
   if(doshift) then
      !
      ! Entering vacuum
      !
      select case(amrray_icross) 
      case(1)
         !
         ! Shift a tiny bit inward
         !
         dummy      = 1.d0 - margin*small
         ray_cart_x = ray_cart_x * dummy 
         ray_cart_y = ray_cart_y * dummy
         ray_cart_z = ray_cart_z * dummy
      case(2)
         !
         ! Shift a tiny bit outward
         !
         dummy      = 1.d0 + margin*small
         ray_cart_x = ray_cart_x * dummy 
         ray_cart_y = ray_cart_y * dummy
         ray_cart_z = ray_cart_z * dummy
      case(3)
         !
         ! Shift a tiny bit by rotation toward positive z-axis
         ! Note: we cannot have rcylnew=0.d0
         !
         rcylnew  = sqrt(ray_cart_x**2+ray_cart_y**2)
         dummy    = margin*small*ray_cart_z / rcylnew
         ray_cart_x = ray_cart_x - ray_cart_x * dummy
         ray_cart_y = ray_cart_y - ray_cart_y * dummy
         ray_cart_z = ray_cart_z + margin*small*rcylnew
      case(4)
         !
         ! Shift a tiny bit by rotation toward negative z-axis
         ! Note: we cannot have rcylnew=0.d0
         !
         rcylnew  = sqrt(ray_cart_x**2+ray_cart_y**2)
         dummy    = margin*small*ray_cart_z / rcylnew
         ray_cart_x = ray_cart_x + ray_cart_x * dummy
         ray_cart_y = ray_cart_y + ray_cart_y * dummy
         ray_cart_z = ray_cart_z - margin*small*rcylnew
      case(5)
         !
         ! Shift a tiny bit clockwise
         !
         ray_cart_x = ray_cart_x + margin*small*ray_cart_y
         ray_cart_y = ray_cart_y - margin*small*ray_cart_x
      case(6)
         !
         ! Shift a tiny bit counter-clockwise
         !
         ray_cart_x = ray_cart_x - margin*small*ray_cart_y
         ray_cart_y = ray_cart_y + margin*small*ray_cart_x
      end select
   endif
   !
   ! Do a check if all is OK. Only for debugging... 
   !
   !###########################################################
   if(amrray_selfcheck) then
      if(amr_tree_present) then
         if(associated(amrray_nextcell)) then
            dummyflag = .true.
         else
            dummyflag = .false.
         endif
      else
         if(amrray_ix_next.gt.0) then
            dummyflag = .true.
         else
            dummyflag = .false.
         endif
      endif
      if(dummyflag) then
         if(amr_tree_present) then
            do iddr=1,3
               bxi(1,iddr) = amr_finegrid_xi(amrray_nextcell%ixyzf(iddr),iddr,amrray_nextcell%level)
               bxi(2,iddr) = amr_finegrid_xi(amrray_nextcell%ixyzf(iddr)+1,iddr,amrray_nextcell%level)
            enddo
         else
            bxi(1,1) = amr_finegrid_xi(amrray_ix_next,1,0)
            bxi(2,1) = amr_finegrid_xi(amrray_ix_next+1,1,0)
            bxi(1,2) = amr_finegrid_xi(amrray_iy_next,2,0)
            bxi(2,2) = amr_finegrid_xi(amrray_iy_next+1,2,0)
            bxi(1,3) = amr_finegrid_xi(amrray_iz_next,3,0)
            bxi(2,3) = amr_finegrid_xi(amrray_iz_next+1,3,0)
         endif
         r02    = ray_cart_x*ray_cart_x + ray_cart_y*ray_cart_y + ray_cart_z*ray_cart_z
         r0     = sqrt(r02)
         theta0 = acos(ray_cart_z/(r0+1d-199))
         if(ray_cart_x.eq.0.d0) then
            if(ray_cart_y.ge.0.d0) then
               phi0   = pihalf
            else
               phi0   = 3.d0*pihalf 
            endif
         else
            phi0   = atan(ray_cart_y/ray_cart_x)
            if(ray_cart_x.gt.0.d0) then
               if(ray_cart_y.lt.0.d0) then
                  phi0   = phi0 + twopi
               endif
            else
               phi0   = phi0 + pi
            endif
         endif
         if(phi0.le.small) then
            if(bxi(2,3).eq.twopi) then
               phi0=twopi
            endif
         elseif(phi0.ge.twopi-small) then
            if(bxi(1,3).eq.0.d0) then
               phi0=0.d0
            endif
         endif
         if((r0>bxi(2,1)*oneplust).or.                 &
              (r0<bxi(1,1)*oneminust).or.              &
              (theta0>bxi(2,2)+tol).or.              &
              (theta0<bxi(1,2)-tol).or.              &
              (phi0>bxi(2,3)+tol).or.                &
              (phi0<bxi(1,3)-tol)) then
            write(stdo,*) 'ERROR: Photon outside of NEXT cell'
            write(stdo,*) ray_cart_x,ray_cart_y,ray_cart_z
            write(stdo,*) ray_cart_dirx,ray_cart_diry,ray_cart_dirz,  &
                 sqrt(ray_cart_dirx**2+ray_cart_diry**2+ray_cart_dirz**2)
            write(stdo,*) r0,theta0,phi0
            if(amr_tree_present) write(stdo,*) associated(amrray_nextcell),amrray_nextcell%id,crossequator
            write(stdo,*) bxi
            write(stdo,*) 'amrray_icross = ',amrray_icross,', cross_ds = ',cross_ds,', idir = ',idir
            if(amr_tree_present) then
               if(.not.associated(amrray_cell)) then 
                  write(stdo,*) 'Started outside of cell'
               else
                  write(stdo,*) 'Started inside of cell'
                  write(stdo,*) axi
               endif
            else
               if(amrray_ix_curr.le.0) then 
                  write(stdo,*) 'Started outside of cell'
               else
                  write(stdo,*) 'Started inside of cell'
                  write(stdo,*) axi
               endif
            endif
            write(stdo,*) 'sin2 cos2 = ',st2,ct2
            if(amr_tree_present) then
               st2   = amrray_finegrid_sintsq1(amrray_cell%ixyzf(2),amrray_cell%level)
               ct2   = amrray_finegrid_costsq1(amrray_cell%ixyzf(2),amrray_cell%level)
            else
               st2   = amrray_finegrid_sintsq1(amrray_iy_curr,0)
               ct2   = amrray_finegrid_costsq1(amrray_iy_curr,0)
            endif
            pa    = ct2*(ray_cart_dirx**2+ray_cart_diry**2)-st2*ray_cart_dirz**2
            pb    = 2.d0*ct2*(ray_cart_dirx*ray_cart_x + ray_cart_diry*ray_cart_y) - &
                 2.d0*st2*ray_cart_dirz*ray_cart_z
            pc    = ct2*(ray_cart_x**2+ray_cart_y**2)-st2*ray_cart_z**2
            write(stdo,*) pa,pb/r0,pc/r0**2
            write(stdo,*) 4*pa*pc/pb**2
            stop 1132
         endif
      endif
   endif
   !###########################################################
endif
!
! Now, wrap things up
!
if(amr_tree_present) then
   !
   ! If the AMR tree exists, then we must still check if we 
   ! have to dig into subcells
   !
   if(associated(amrray_nextcell)) then
      !
      ! We enter a new cell
      !
      if(.not.amrray_nextcell%leaf) then
         !
         ! The new cell is in fact a nodal branch with children, then
         ! we must find this sub-cell. 
         !
         if(amr_xyzdim(idir).eq.0) ilr=1
         call amr_xyz_to_rthphi(ray_cart_x,ray_cart_y,ray_cart_z,r,theta,phi)
         call amrray_find_subcell(amrray_nextcell,idir,ilr,r,theta,phi)
         !
         ! For the corner-based integration we must re-calculate icross
         !
         amrray_icross = -(2*(idir-1)+ilr)
      endif
      ray_indexnext = amrray_nextcell%leafindex
   else
      !
      ! We enter (or are still in) the outer region (outside the grid)
      !
      ray_indexnext = 0
      !
   endif
else
   !
   ! For a regular grid we just calculate the index of the next cell
   !
   if(amrray_ix_next.gt.0) then
      !
      ! We are in a cell
      !
      ray_indexnext = amrray_ix_next+((amrray_iy_next-1)+(amrray_iz_next-1)*amr_grid_ny)*amr_grid_nx
   else
      !
      ! We enter (or are still in) the outer region (outside the grid)
      !
      ray_indexnext = 0
   endif
endif
!
end subroutine amrray_find_next_location_spher



!--------------------------------------------------------------------------
!                    FIND CROSSING WITH A SPHERE
!
! If in addition to the AMR grid, we also have spheres (i.e. "stars")
! then here is a routine to calculate the crossing of the ray with such
! a sphere. It is written in such a way that it would also work if the
! starting point of the ray is very far away compared to the radius of
! the sphere. This is important if you have very small stars in a very
! large simulation volume: the precision limitations could otherwise 
! cause trouble. This is now taken care of.
!--------------------------------------------------------------------------
subroutine amrray_find_crossing_with_sphere(    &
     ray_cart_x,ray_cart_y,ray_cart_z,          &
     ray_cart_dirx,ray_cart_diry,ray_cart_dirz, &
     isph,ds_end,ds,ihit,mirror)
implicit none
integer :: isph,ihit
logical :: mirror
double precision :: ray_cart_x,ray_cart_y,ray_cart_z
double precision :: ray_cart_dirx,ray_cart_diry,ray_cart_dirz
double precision :: ds,ds_end
double precision :: dx,dy,dz,xx,yy,zz,r2,rr,det,xdd
!
ihit   = 0
ds     = 1.d99
!
! Find the difference vector
!
dx   = ray_cart_x - amrray_spheres_pos(1,isph)
dy   = ray_cart_y - amrray_spheres_pos(2,isph)
dz   = ray_cart_z - amrray_spheres_pos(3,isph)
!
! Calculate the distance between the current position and the
! center of the sphere
!
r2   = dx**2 + dy**2 + dz**2
rr   = sqrt( r2 )
!
! If this is large, then we do a special treatment, otherwise
! not
!
if(rr.gt.1d3*amrray_spheres_r(isph)) then
   !
   ! Special treatment
   !
   ! Move along the direction by that distance. This may either
   ! move it further away from the sphere by at maximum a factor
   ! of two (which is not problematic) or it moves it closer
   ! to the sphere, potentially extremely much closer. If that
   ! happens, we will get far more accurate results for the 
   ! later calculations, hence the shift.
   !
   xx   = dx + rr*ray_cart_dirx
   yy   = dy + rr*ray_cart_diry
   zz   = dz + rr*ray_cart_dirz
   !
   ! Now compute the determinant of the quadratic equation
   !
   r2   = xx**2 + yy**2 + zz**2
   xdd  = xx * ray_cart_dirx + yy * ray_cart_diry + zz * ray_cart_dirz
   det  = xdd**2 + amrray_spheres_r(isph)**2 - r2
   if(det.gt.0.d0) then
      !
      ! There are solutions. Since we are definitely outside
      ! of this sphere, we clearly must take the smallest of
      ! the two solutions.
      ! 
      ds = -xdd - sqrt(det) + rr
      if((ds.gt.0.d0).and.(ds.lt.ds_end)) then
         ihit = -1  ! Signal with -1 that we are *entering* the sphere
      endif
   endif
else
   ! 
   ! Normal treatment
   !
   ! Compute the determinant of the quadratic equation
   !
   xdd  = dx * ray_cart_dirx + dy * ray_cart_diry + dz * ray_cart_dirz
   det  = xdd**2 + amrray_spheres_r(isph)**2 - r2
   if(det.gt.0.d0) then
      !
      ! There are solutions. Simply check if we are moving
      ! inward or outward. Use a small margin, so that we
      ! rather classify as the ray moving outward than inward
      ! in case of uncertainty.
      !
      if(xdd.lt.-tiny*rr) then
         !
         ! Moving inward. This means that there IS a real
         ! crossing with the sphere. Check if we are entering
         ! or leaving the sphere. Use a margin so that we 
         ! rather leave the sphere than enter, in case it
         ! is unclear.
         !
         if(rr.lt.onepluss*amrray_spheres_r(isph)) then
            !
            ! Take the largest of the two solutions: we will
            ! exit the sphere.
            !
            ds = -xdd + sqrt(det)
            if((ds.gt.0.d0).and.(ds.lt.ds_end)) then
               ihit = 1  ! Signal with +1 that we are inside and exiting
!!#######################################
!               write(stdo,*) '** YES! Isph = ',isph
!!#######################################
            endif
         else
            !
            ! Take the smallest of the two solutions: we just
            ! enter the sphere.
            !
            ds = -xdd - sqrt(det)
            if((ds.gt.0.d0).and.(ds.lt.ds_end)) then
               ihit = -1  ! Signal with -1 that we are entering the sphere
            endif
         endif
      else
         !
         ! Moving outward. Then there is only a real solution if
         ! we are already inside a sphere. 
         !
         if(rr.lt.oneminuss*amrray_spheres_r(isph)) then
            ds = -xdd + sqrt(det)
            if((ds.gt.0.d0).and.(ds.lt.ds_end)) then
               ihit = 1  ! Signal with +1 that we are inside and exiting
            endif
         endif
         !
      endif
   endif
endif
!
! If we must check mirror symmetry sphere
!
if(mirror) then
   !
   ! Find the difference vector
   !
   dx   = ray_cart_x - amrray_spheres_pos(1,isph)
   dy   = ray_cart_y - amrray_spheres_pos(2,isph)
   dz   = ray_cart_z + amrray_spheres_pos(3,isph)
   !
   ! Calculate the distance between the current position and the
   ! center of the sphere
   !
   r2   = dx**2 + dy**2 + dz**2
   rr   = sqrt( r2 )
   !
   ! If this is large, then we do a special treatment, otherwise
   ! not
   !
   if(rr.gt.1d3*amrray_spheres_r(isph)) then
      !
      ! Special treatment
      !
      ! Move along the direction by that distance. This may either
      ! move it further away from the sphere by at maximum a factor
      ! of two (which is not problematic) or it moves it closer
      ! to the sphere, potentially extremely much closer. If that
      ! happens, we will get far more accurate results for the 
      ! later calculations, hence the shift.
      !
      xx   = dx + rr*ray_cart_dirx
      yy   = dy + rr*ray_cart_diry
      zz   = dz + rr*ray_cart_dirz
      !
      ! Now compute the determinant of the quadratic equation
      !
      r2   = xx**2 + yy**2 + zz**2
      xdd  = xx * ray_cart_dirx + yy * ray_cart_diry + zz * ray_cart_dirz
      det  = xdd**2 + amrray_spheres_r(isph)**2 - r2
      if(det.gt.0.d0) then
         !
         ! There are solutions. Since we are definitely outside
         ! of this sphere, we clearly must take the smallest of
         ! the two solutions.
         ! 
         ds = -xdd - sqrt(det) + rr
         if((ds.gt.0.d0).and.(ds.lt.ds_end)) then
            ihit = -1  ! Signal with -1 that we are *entering* the sphere
         endif
      endif
   else
      ! 
      ! Normal treatment
      !
      ! Compute the determinant of the quadratic equation
      !
      xdd  = dx * ray_cart_dirx + dy * ray_cart_diry + dz * ray_cart_dirz
      det  = xdd**2 + amrray_spheres_r(isph)**2 - r2
      if(det.gt.0.d0) then
         !
         ! There are solutions. Simply check if we are moving
         ! inward or outward. Use a small margin, so that we
         ! rather classify as the ray moving outward than inward
         ! in case of uncertainty.
         !
         if(xdd.lt.-tiny*rr) then
            !
            ! Moving inward. This means that there IS a real
            ! crossing with the sphere. Check if we are entering
            ! or leaving the sphere. Use a margin so that we 
            ! rather leave the sphere than enter, in case it
            ! is unclear.
            !
            if(rr.lt.onepluss*amrray_spheres_r(isph)) then
               !
               ! Take the largest of the two solutions: we will
               ! exit the sphere.
               !
               ds = -xdd + sqrt(det)
               if((ds.gt.0.d0).and.(ds.lt.ds_end)) then
                  ihit = 1  ! Signal with +1 that we are inside and exiting
               endif
            else
               !
               ! Take the smallest of the two solutions: we just
               ! enter the sphere.
               !
               ds = -xdd - sqrt(det)
               if((ds.gt.0.d0).and.(ds.lt.ds_end)) then
                  ihit = -1  ! Signal with -1 that we are entering the sphere
               endif
            endif
         else
            !
            ! Moving outward. Then there is only a real solution if
            ! we are already inside a sphere. 
            !
            if(rr.lt.oneminuss*amrray_spheres_r(isph)) then
               ds = -xdd + sqrt(det)
               if((ds.gt.0.d0).and.(ds.lt.ds_end)) then
                  ihit = 1  ! Signal with +1 that we are inside and exiting
               endif
            endif
            !
         endif
      endif
   endif
endif
!
end subroutine amrray_find_crossing_with_sphere


!--------------------------------------------------------------------------
!                        DESTRUCTOR FOR AMRRAY
!
! If you use the pixmap stuff, then to deallocate all arrays one must
! call this routine.
!--------------------------------------------------------------------------
subroutine amrray_cleanup()
implicit none
integer :: ierr
if(allocated(pixmap_xi)) then
   deallocate(pixmap_xi,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap grid'
      stop 201
   endif
endif
if(allocated(pixmap_xc)) then
   deallocate(pixmap_xc,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap grid'
      stop 201
   endif
endif
if(allocated(pixmap_yi)) then
   deallocate(pixmap_yi,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap grid'
      stop 201
   endif
endif
if(allocated(pixmap_yc)) then
   deallocate(pixmap_yc,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap grid'
      stop 201
   endif
endif
if(allocated(pixmap_cart_xs)) then
   deallocate(pixmap_cart_xs,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap startpos'
      stop 201
   endif
endif
if(allocated(pixmap_cart_ys)) then
   deallocate(pixmap_cart_ys,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap startpos'
      stop 201
   endif
endif
if(allocated(pixmap_cart_zs)) then
   deallocate(pixmap_cart_zs,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap startpos'
      stop 201
   endif
endif
if(allocated(pixmap_cart_dxs)) then
   deallocate(pixmap_cart_dxs,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap startdir'
      stop 201
   endif
endif
if(allocated(pixmap_cart_dys)) then
   deallocate(pixmap_cart_dys,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap startdir'
      stop 201
   endif
endif
if(allocated(pixmap_cart_dzs)) then
   deallocate(pixmap_cart_dzs,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap startdir'
      stop 201
   endif
endif
if(allocated(amrray_finegrid_sintsq1)) then
   deallocate(amrray_finegrid_sintsq1,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate amrray_finegrid_sintsq1'
      stop 201
   endif
endif
if(allocated(amrray_finegrid_sintsq2)) then
   deallocate(amrray_finegrid_sintsq2,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate amrray_finegrid_sintsq2'
      stop 201
   endif
endif
if(allocated(amrray_finegrid_costsq1)) then
   deallocate(amrray_finegrid_costsq1,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate amrray_finegrid_costsq1'
      stop 201
   endif
endif
if(allocated(amrray_finegrid_costsq2)) then
   deallocate(amrray_finegrid_costsq2,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate amrray_finegrid_costsq2'
      stop 201
   endif
endif
if(allocated(amrray_finegrid_sinp1)) deallocate(amrray_finegrid_sinp1)
if(allocated(amrray_finegrid_sinp2)) deallocate(amrray_finegrid_sinp2)
if(allocated(amrray_finegrid_cosp1)) deallocate(amrray_finegrid_cosp1)
if(allocated(amrray_finegrid_cosp2)) deallocate(amrray_finegrid_cosp2)
if(allocated(amrray_spheres_r)) deallocate(amrray_spheres_r)
if(allocated(amrray_spheres_pos)) deallocate(amrray_spheres_pos)
if(allocated(amrray_spheres_outsidegrid)) deallocate(amrray_spheres_outsidegrid)
if(allocated(amrray_spheres_sphidx)) deallocate(amrray_spheres_sphidx)
!
end subroutine amrray_cleanup


!--------------------------------------------------------------------------
!        MAKE SET OF STARTING POSITIONS FOR RECTANGULAR PIXEL MAP
!
! If you want to make a rectangular image, then you need a set of 
! starting positions for the numerical ray-tracing. This is what this
! subroutine does. Since there could be more variants of this subroutine,
! this version is called amrray_rect_pixmap_startpos_method1().
!
! ARGUMENTS:
!  nx,ny           Number of pixels horizontally (nx) and vertically (ny)
!  xl,xr,yl,yr     The boundaries of the images in x and y. Note that
!                  for incl=0 the x direction corresponds truly to the
!                  x-direction and the same for y. In other words: incl=0
!                  means we view from the top. incl=90 means we view 
!                  from the side.
!  incl            Inclination, see above.
!  phi             x-y rotation in 3-D.
!  psi             x-y rotation in image.
!
!--------------------------------------------------------------------------
subroutine amrray_rect_pixmap_startpos_method1(nx,ny,xl,xr,yl,yr,& 
                incl,phi,psi)
implicit none
integer :: ierr
integer :: nx,ny,ix,iy
doubleprecision :: xl,xr,yl,yr,incl,phi,psi
doubleprecision :: dmax,dx,dy,cosincl,sinincl,cosphi,sinphi,cospsi,sinpsi
doubleprecision :: dirx,diry,dirz,pi,dirxbk,dirybk,dirzbk,xbk,ybk,zbk
parameter(pi=3.1415926535897932385d0)
!
if(.not.allocated(amr_grid_xi)) then
   write(stdo,*) 'ERROR in AMRRAY Module: Spatial grid not defined.'
   stop 8203
endif
dmax = max(abs(amr_grid_xi(1,1)),abs(amr_grid_xi(amr_grid_nx+1,1)),  &
           abs(amr_grid_xi(1,2)),abs(amr_grid_xi(amr_grid_ny+1,2)),  &
           abs(amr_grid_xi(1,3)),abs(amr_grid_xi(amr_grid_nz+1,3)))  &
            * 1.1
!
! Allocate the pixel coordinates
!
if(allocated(pixmap_xi)) then
   deallocate(pixmap_xi,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap grid'
      stop 201
   endif
endif
allocate(pixmap_xi(nx+1),STAT=ierr)
if(ierr.ne.0) then 
   write(stdo,*) 'ERROR in AMRRAY Module: Could not allocate pixmap grid'
   stop 201
endif
if(allocated(pixmap_yi)) then
   deallocate(pixmap_yi,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap grid'
      stop 201
   endif
endif
allocate(pixmap_yi(ny+1),STAT=ierr)
if(ierr.ne.0) then 
   write(stdo,*) 'ERROR in AMRRAY Module: Could not allocate pixmap grid'
   stop 201
endif
if(allocated(pixmap_xc)) then
   deallocate(pixmap_xc,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap grid'
      stop 201
   endif
endif
allocate(pixmap_xc(nx),STAT=ierr)
if(ierr.ne.0) then 
   write(stdo,*) 'ERROR in AMRRAY Module: Could not allocate pixmap grid'
   stop 201
endif
if(allocated(pixmap_yc)) then
   deallocate(pixmap_yc,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap grid'
      stop 201
   endif
endif
allocate(pixmap_yc(ny),STAT=ierr)
if(ierr.ne.0) then 
   write(stdo,*) 'ERROR in AMRRAY Module: Could not allocate pixmap grid'
   stop 201
endif
!
! Allocate the pixmap starting positions
!
if(allocated(pixmap_cart_xs)) then
   deallocate(pixmap_cart_xs,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap startpos'
      stop 201
   endif
endif
allocate(pixmap_cart_xs(1:nx,1:ny),STAT=ierr)
if(ierr.ne.0) then 
   write(stdo,*) 'ERROR in AMRRAY Module: Could not allocate pixmap startpos'
   stop 201
endif
if(allocated(pixmap_cart_ys)) then
   deallocate(pixmap_cart_ys,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap startpos'
      stop 201
   endif
endif
allocate(pixmap_cart_ys(1:nx,1:ny),STAT=ierr)
if(ierr.ne.0) then 
   write(stdo,*) 'ERROR in AMRRAY Module: Could not allocate pixmap startpos'
   stop 201
endif
if(allocated(pixmap_cart_zs)) then
   deallocate(pixmap_cart_zs,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap startpos'
      stop 201
   endif
endif
allocate(pixmap_cart_zs(1:nx,1:ny),STAT=ierr)
if(ierr.ne.0) then 
   write(stdo,*) 'ERROR in AMRRAY Module: Could not allocate pixmap startpos'
   stop 201
endif
!
! Allocate the pixmap starting directions
!
if(allocated(pixmap_cart_dxs)) then
   deallocate(pixmap_cart_dxs,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap startdir'
      stop 201
   endif
endif
allocate(pixmap_cart_dxs(1:nx,1:ny),STAT=ierr)
if(ierr.ne.0) then 
   write(stdo,*) 'ERROR in AMRRAY Module: Could not allocate pixmap startdir'
   stop 201
endif
if(allocated(pixmap_cart_dys)) then
   deallocate(pixmap_cart_dys,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap startdir'
      stop 201
   endif
endif
allocate(pixmap_cart_dys(1:nx,1:ny),STAT=ierr)
if(ierr.ne.0) then 
   write(stdo,*) 'ERROR in AMRRAY Module: Could not allocate pixmap startdir'
   stop 201
endif
if(allocated(pixmap_cart_dzs)) then
   deallocate(pixmap_cart_dzs,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMRRAY Module: Could not dallocate pixmap startdir'
      stop 201
   endif
endif
allocate(pixmap_cart_dzs(1:nx,1:ny),STAT=ierr)
if(ierr.ne.0) then 
   write(stdo,*) 'ERROR in AMRRAY Module: Could not allocate pixmap startdir'
   stop 201
endif
!
! Now set up the grid
!
dx = (xr-xl)/nx
dy = (yr-yl)/ny
do ix=1,nx+1
   pixmap_xi(ix) = xl+(ix-1)*dx
enddo
do iy=1,ny+1
   pixmap_yi(iy) = yl+(iy-1)*dy
enddo
do ix=1,nx
   pixmap_xc(ix) = 0.5d0 * ( pixmap_xi(ix) + pixmap_xi(ix+1) )
enddo
do iy=1,ny
   pixmap_yc(iy) = 0.5d0 * ( pixmap_yi(iy) + pixmap_yi(iy+1) )
enddo
!
! The angles
!
cosincl = cos(incl*pi/180.)
sinincl = sin(incl*pi/180.)
cosphi  = cos(phi*pi/180.)
sinphi  = sin(phi*pi/180.)
cospsi  = cos(psi*pi/180.)
sinpsi  = sin(psi*pi/180.)
!
! All the starting directions are the same
!
dirx = 0.d0
diry = 0.d0
dirz = 1.d0
dirxbk = dirx
dirybk = diry
dirx = cospsi*dirxbk-sinpsi*dirybk
diry = sinpsi*dirxbk+cospsi*dirybk
dirybk = diry
dirzbk = dirz
diry = cosincl*dirybk-sinincl*dirzbk
dirz = sinincl*dirybk+cosincl*dirzbk
dirxbk = dirx
dirybk = diry
dirx = cosphi*dirxbk-sinphi*dirybk
diry = sinphi*dirxbk+cosphi*dirybk
!
! Now make the starting positions for integration
!
do ix=1,nx
   do iy=1,ny
      !
      ! Direction is copied
      !
      pixmap_cart_dxs(ix,iy) = dirx
      pixmap_cart_dys(ix,iy) = diry
      pixmap_cart_dzs(ix,iy) = dirz
      !
      ! First without rotation
      !
      pixmap_cart_xs(ix,iy) = pixmap_xc(ix)
      pixmap_cart_ys(ix,iy) = pixmap_yc(iy)
      pixmap_cart_zs(ix,iy) = -dmax
      !
      ! Now rotate in image plane
      !
      xbk                   = pixmap_cart_xs(ix,iy)
      ybk                   = pixmap_cart_ys(ix,iy)
      pixmap_cart_xs(ix,iy) = cospsi*xbk - sinpsi*ybk
      pixmap_cart_ys(ix,iy) = sinpsi*xbk + cospsi*ybk
      !
      ! Now tilt to the right inclination
      !
      ybk                   = pixmap_cart_ys(ix,iy)
      zbk                   = pixmap_cart_zs(ix,iy)
      pixmap_cart_ys(ix,iy) = cosincl*ybk - sinincl*zbk
      pixmap_cart_zs(ix,iy) = sinincl*ybk + cosincl*zbk
      !
      ! Now rotate in 3-D along z-axis
      !
      xbk                   = pixmap_cart_xs(ix,iy)
      ybk                   = pixmap_cart_ys(ix,iy)
      pixmap_cart_xs(ix,iy) = cosphi*xbk - sinphi*ybk
      pixmap_cart_ys(ix,iy) = sinphi*xbk + cosphi*ybk
   enddo
enddo
!
end subroutine amrray_rect_pixmap_startpos_method1



end module amrray_module
