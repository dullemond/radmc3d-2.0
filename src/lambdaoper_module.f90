!========================================================================
!                       LAMBDA OPERATOR MODULE
!
! For the Lambda Iteration scheme (and its sister, the Accelerated
! Lambda Iteration scheme) a Lambda Operator is required.
!
! A Lambda operator is simply an operator (subroutine) that computes
! the mean intensity at some wavelength using the integrals over a large
! number of rays. In contrast to the mean intensity subroutine in the
! montecarlo_module, for the Lambda operator we do not "follow a photon
! along its way" but we instead "integrate the formal transfer equation
! along a ray". If no an-isotropic dust scattering is included in the
! problem, and the scattering source function is known, this gives a 
! more accurate answer. Also, depending on which method of Lambda
! operator you use, you get less noise. Finally, one can use the
! Accelerated Lambda Iteration scheme to speed up convergence of the
! Lambda Iteration.
!========================================================================
module lambdaoper_module
use mathroutines_module
use rtglobal_module
use montecarlo_module
use amr_module
use amrray_module
use constants_module

double precision :: lambdaoper_shift = 1.d90
integer :: lambdaoper_nfreq,lambdaoper_istar
double precision, allocatable :: lambdaoper_intensity(:),lambdaoper_bgint(:)
double precision, allocatable :: lambdaoper_local_lamdiag(:)
integer :: lambdaoper_order=1
logical :: lambdaoper_addtocell_meanint=.false.
logical :: lambdaoper_vertexbased=.false.
integer :: lambdaoper_random_ray_method=0

doubleprecision :: lambdaoper_snu_prev,lambdaoper_anu_prev,lambdaoper_jnu_prev
doubleprecision :: lambdaoper_snu_curr,lambdaoper_anu_curr,lambdaoper_jnu_curr

double precision, allocatable :: lambdaoper_meanint_cell(:)
double precision, allocatable :: lambdaoper_lambda_diag_cell(:)
double precision, allocatable :: lambdaoper_dstot_cell(:)
double precision, allocatable :: lambdaoper_meanint_vertex(:)
double precision, allocatable :: lambdaoper_lambda_diag_vertex(:)
double precision, allocatable :: lambdaoper_dstot_vertex(:)

contains


!------------------------------------------------------------------------
!                         INIT THIS MODULE 
!------------------------------------------------------------------------
subroutine lambdaoper_init(nf,freq,vertexbased,bgint)
  implicit none
  integer :: nf
  double precision :: freq(1:nf)
  double precision, optional :: bgint(1:nf)
  logical :: vertexbased
  !
  ! First clean up
  !
  call lambdaoper_cleanup()
  !
  ! Then copy some stuff
  !
  lambdaoper_vertexbased = vertexbased
  lambdaoper_nfreq = nf
  !
  ! Allocate and init some arrays
  !
  allocate(lambdaoper_intensity(1:nf),lambdaoper_bgint(1:nf),lambdaoper_local_lamdiag(1:nf))
  lambdaoper_intensity(:) = 0.d0
  lambdaoper_local_lamdiag(:) = 0.d0
  if(present(bgint)) then
     lambdaoper_bgint(:) = bgint(:)
  else
     lambdaoper_bgint(:) = 0.d0
  endif
  !
  ! Since we want to compute the mean intensity, here is the array for that
  !
  if(vertexbased) then
     allocate(lambdaoper_meanint_vertex(amr_nr_vertices_max))
     allocate(lambdaoper_lambda_diag_vertex(amr_nr_vertices_max))
     allocate(lambdaoper_dstot_vertex(amr_nr_vertices_max))
  else
     allocate(lambdaoper_meanint_cell(amr_nrleafs_max))
     allocate(lambdaoper_lambda_diag_cell(amr_nrleafs_max))
     allocate(lambdaoper_dstot_cell(amr_nrleafs_max))
  endif
end subroutine lambdaoper_init


!------------------------------------------------------------------------
!                        CLOSE THIS MODULE
!------------------------------------------------------------------------
subroutine lambdaoper_cleanup()
  implicit none
  if(allocated(lambdaoper_intensity)) deallocate(lambdaoper_intensity)
  if(allocated(lambdaoper_bgint)) deallocate(lambdaoper_bgint)
  if(allocated(lambdaoper_local_lamdiag)) deallocate(lambdaoper_local_lamdiag)
  if(allocated(lambdaoper_meanint_cell)) deallocate(lambdaoper_meanint_cell)
  if(allocated(lambdaoper_meanint_vertex)) deallocate(lambdaoper_meanint_vertex)
  if(allocated(lambdaoper_dstot_cell)) deallocate(lambdaoper_dstot_cell)
  if(allocated(lambdaoper_dstot_vertex)) deallocate(lambdaoper_dstot_vertex)
  if(allocated(lambdaoper_lambda_diag_cell)) deallocate(lambdaoper_lambda_diag_cell)
  if(allocated(lambdaoper_lambda_diag_vertex)) deallocate(lambdaoper_lambda_diag_vertex)
  lambdaoper_nfreq = 0
end subroutine lambdaoper_cleanup


!------------------------------------------------------------------------
!              HELPER ROUTINE FOR RT INTEGRATION STEP
! Make sure that ray_index is properly set before calling this routine.
!------------------------------------------------------------------------
subroutine lambdaoper_integrate_element(inu0,inu1,nf)
  implicit none
  integer :: inu
  double precision :: dtau1(1:nf),theomax(1:nf),e0(1:nf),e1(1:nf)
  double precision :: ax(1:nf),bx(1:nf),xp(1:nf),qdr(1:nf)
  double precision :: src(1:nf),alp(1:nf)
  double precision :: r,theta,phi
  integer :: inu0,inu1,nf
  !
  ! Check if we are/were in a star
  !
  if(lambdaoper_istar.le.0) then
     !
     ! Normal situation: we are not inside a star (istar=0), or we will
     ! hit the star but are not inside the star yet (istar<0).
     !
     ! Decide whether we do 1st order, 2nd order or 3rd order integration
     !
     select case(lambdaoper_order)
     case(-1)
        !
        ! Integrate only the optical depth (useful for escape probability)
        !
        ! Get the src and alp values
        !
        call sources_get_src_alp(inu0,inu1,nf,src,alp,.false.)
        !
        ! Now do the RT, including everything
        !
        lambdaoper_intensity(inu0:inu1) = lambdaoper_intensity(inu0:inu1) + &
             alp(inu0:inu1) * ray_ds
        !
     case(1)
        !
        ! First order integration. Cell-based scheme.
        !
        ! Get the src and alp values
        !
        call sources_get_src_alp(inu0,inu1,nf,src,alp,.false.)
        !
        ! Compute tau, exp(-tau), 1-exp(-tau) and the diagonal of 
        ! (this contribution to) the local approximate lambda operator.
        !
        dtau1(inu0:inu1) = alp(inu0:inu1) * ray_ds
        xp(inu0:inu1)    = exp(-dtau1(inu0:inu1))
        do inu=inu0,inu1
           if(dtau1(inu).lt.0.d0) stop 7329
           if(dtau1(inu).gt.1d-4) then
              bx(inu)  = 1.d0 - xp(inu)
              lambdaoper_local_lamdiag(inu) = ( dtau1(inu) - 1.d0 + xp(inu) ) / dtau1(inu)
           else
              bx(inu)  = dtau1(inu)
              lambdaoper_local_lamdiag(inu) = 0.5d0 * dtau1(inu)
           endif
        enddo
        !
        ! Now do the RT, including everything
        !
        lambdaoper_intensity(inu0:inu1) = xp(inu0:inu1) * lambdaoper_intensity(inu0:inu1) + &
             bx(inu0:inu1) * src(inu0:inu1)
        !
     case(2)
        !
        ! Second order integration. Corner-based (=vertex-based) scheme.
        ! Here all sources have all been pre-calculated on
        ! the cell corners.

################### HERE I STILL HAVE TO ADD THE (INU0:INU1) AND THE CORRECT APPROX LAMBDA OPERATUR ###############
        
        !
        ! Backup snu and anu at previous crossing
        !
        lambdaoper_snu_prev = lambdaoper_snu_curr
        lambdaoper_anu_prev = lambdaoper_anu_curr
        lambdaoper_jnu_prev = lambdaoper_jnu_curr
        !
        ! Find snu and anu at new crossing
        !
        ! NOTE: If lambdaoper_interpol_jnu is set, then snu is actually
        !       j_nu (the emissivity). It will be translated below.
        !
        if(igrid_type.lt.100) then 
           if(igrid_coord.lt.100) then
              !
              ! Cartesian
              !
              call lambdaoper_find_srcalp_interpol(       &
                   ray_cart_x,ray_cart_y,ray_cart_z,  &
                   lambdaoper_snu_curr,lambdaoper_anu_curr)
           elseif(igrid_coord.lt.200) then
              !
              ! Spherical 
              !
              call amr_xyz_to_rthphi(ray_cart_x,ray_cart_y,ray_cart_z, &
                                           r,theta,phi)
              call lambdaoper_find_srcalp_interpol(    &
                   r,theta,phi,lambdaoper_snu_curr,lambdaoper_anu_curr)
           else
              stop 498
           endif
        else
           stop 499
        endif
        !
        ! If lambdaoper_interpol_jnu is set, then what we get back
        ! from the camera_find_srcalp_interpol() is in fact the
        ! emissivity j_nu instead of the source function S_nu.
        ! So translate either snu to jnu or vice versa.
        !
        if(lambdaoper_interpol_jnu) then
           jnu_curr = snu_curr
           snu_curr = jnu_curr / (anu_curr+1d-99)
        else
           jnu_curr = snu_curr * (anu_curr+1d-99)
        endif
        !
        ! Now do the integration step
        !
        if(ray_index.gt.0) then
           dtau1     = 0.5d0 * ( lambdaoper_anu_prev + lambdaoper_anu_curr ) * ray_ds
           theomax   = 0.5d0 * ( lambdaoper_jnu_prev + lambdaoper_jnu_curr ) * ray_ds
           if(dtau1.gt.1.d-6) then
              xp = exp(-dtau1)
              e0 = 1.d0 - xp
              e1 = dtau1 - e0
              bx = e1 / dtau1
              ax = e0 - bx
           else
              ax = 0.5d0 * dtau1
              bx = 0.5d0 * dtau1
              xp = 1.d0 - dtau1 
           endif
           if(dtau1.gt.1e-9) then
              qdr = ax * snu_prev + bx * snu_curr
           else
              qdr = theomax
           endif
           qdr = min(qdr,theomax)
           lambdaoper_intensity(lambdaoper_inu) = xp * lambdaoper_intensity(lambdaoper_inu) + qdr
        endif
     case(3) 
        !
        ! Third order integration
        !
        write(*,*) 'ERROR: Third-order integration not yet implemented'
        stop
     end select
  elseif(camera_incl_stars.ne.0) then
     !
     ! Special situation: we are inside a star.
     ! So set the intensity to that of the star.
     !
     lambdaoper_intensity(lambdaoper_inu) = find_starlight_interpol(sources_frequencies(lambdaoper_inu), &
                                                    lambdaoper_istar)
  endif
end subroutine lambdaoper_integrate_element


!------------------------------------------------------------------------
!            INTEGRATE FORMAL TRANSFER EQUATION ALONG A RAY
!------------------------------------------------------------------------
subroutine lambdaoper_integrate_ray()
  implicit none
  double precision :: s
  doubleprecision :: snu_prev,anu_prev,jnu_prev,snu_curr,anu_curr,jnu_curr
  !
  ! Check
  !
  if(lambdaoper_shift.gt.0.d0) then
     write(stdo,*) 'ERROR in Lambda operator: Shift>0.'
     stop
  endif
  if(.not.allocated(lambdaoper_intensity)) then
     write(stdo,*) 'ERROR in Lambda operator: lambdaoper_intensity not allocated.'
     stop
  endif
  !
  ! Start the intensity
  !
  lambdaoper_intensity(inu0:inu1) = lambdaoper_bgint(inu0:inu1)
  !
  ! Set defaults
  !
  snu_prev = 0.d0
  anu_prev = 0.d0
  jnu_prev = 0.d0
  snu_curr = 0.d0
  anu_curr = 0.d0
  jnu_curr = 0.d0
  !
  ! Check the dust scattering mode, because this only works under
  ! certain conditions.
  !
  if(scattering_mode.ge.2) then
     write(stdo,*) 'ERROR in Lambda module: The Lambda Operator only'
     write(stdo,*) '      works with dust scattering if the scattering'
     write(stdo,*) '      is isotropic.'
     stop
  endif
  if(scattering_mode.ge.1) then
     if(.not.allocated(mcscat_scatsrc_iquv)) stop 3005
  endif
  !
  ! If mirror symmetry is switched on, then we must assure that 
  ! z is above the midplane
  !
  if(amrray_mirror_equator) then
     if(ray_cart_z.eq.0.d0) then
        ray_cart_dirz = abs(ray_cart_dirz)
     elseif(ray_cart_z.lt.0.d0) then
        ray_cart_z    = -ray_cart_z
        ray_cart_dirz = -ray_cart_dirz
     endif
  endif
  !
  ! We start outside the domain, so all indices must be 0
  !
  ray_index     = 0
  ray_indexnext = 0
  !
  ! Now trace
  !
  arrived = .false.
  do while(.not.arrived)
     !
     ! Back up the current position
     !
     ray_prev_x     = ray_cart_x
     ray_prev_y     = ray_cart_y
     ray_prev_z     = ray_cart_z
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
                ray_index,ray_indexnext,ray_ds,arrived,              &
                levelnext=levelnext)
           !
           ! If we went through one of the boundaries of the
           ! domain, we might need to implement the thermal
           ! boundary, if present
           !
           if(incl_thermbc.ne.0) then
              if((ray_index.eq.0).and.(ray_indexnext.ne.0)) then
                 !
                 ! First determine which boundary we just went
                 ! through 
                 !
                 if(amrray_icross.ge.0) stop 3902  ! Self-consistency check
                 if(amrray_icross.eq.-1) then
                    bc_idir = 1
                    bc_ilr  = 1
                 elseif(amrray_icross.eq.-2) then
                    bc_idir = 1
                    bc_ilr  = 2
                 elseif(amrray_icross.eq.-3) then
                    bc_idir = 2
                    bc_ilr  = 1
                 elseif(amrray_icross.eq.-4) then
                    bc_idir = 2
                    bc_ilr  = 2
                 elseif(amrray_icross.eq.-5) then
                    bc_idir = 3
                    bc_ilr  = 1
                 elseif(amrray_icross.eq.-6) then
                    bc_idir = 3
                    bc_ilr  = 2
                 else
                    stop 4944
                 endif
                 !
                 ! Now set the intensity to the thermal
                 ! emission of the boundary, if the
                 ! thermal boundary is switched on
                 !
                 if(thermal_bc_active(bc_ilr,bc_idir)) then
                    do inu=1,lambdaoper_nfreq
                       lambdaoper_intensity(inu) = bplanck(thermal_bc_temp(bc_ilr,bc_idir), &
                                                lambdaoper_frequencies(inu))
                    enddo
                 endif
              endif
           endif
           !
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
        ! Set the lambdaoper_istar, which is used for treating finite-size
        ! stars (instead of point source stars)
        !
        lambdaoper_istar = amrray_ispherehit
        !
     else
        write(stdo,*) 'SORRY: Delaunay or Voronoi grids not yet implemented'
        stop
     endif
     !
     ! Path length
     !
     ray_ds    = sqrt( (ray_cart_x-ray_prev_x)**2 + (ray_cart_y-ray_prev_y)**2 + (ray_cart_z-ray_prev_z)**2 )
     !
     ! Now do the integration along this path element
     !
     call lambdaoper_integrate_element()
     !
     ! Add information to the cell
     !
     ! NOTE: We add the final intensity to the cell, even though perhaps
     !       the average would be better. But the Lambda_diag approximate
     !       operator is calculated on the basis of the final intensity.
     !
     if(lambdaoper_addtocell_meanint) then
        lambdaoper_meanint_cell(ray_index) = lambdaoper_meanint_cell(ray_index) + &
             ray_ds * lambdaoper_intensity(lambdaoper_inu)
        lambdaoper_lambda_diag_cell(ray_index) = lambdaoper_lambda_diag_cell(ray_index) + &
             ray_ds * lambdaoper_local_lamdiag
        lambdaoper_dstot_cell(ray_index) = lambdaoper_dstot_cell(ray_index) + ray_ds 
     endif
     !
     ! Increase s
     !
     s = s + ray_ds
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
  enddo
  !
end subroutine lambdaoper_integrate_ray


!------------------------------------------------------------------------
!                        FIND A RANDOM RAY
!------------------------------------------------------------------------
subroutine lambdaoper_draw_random_ray(doshift)
  implicit none
  integer :: icell
  logical :: doshift
  !
  ! We have different ways to choose random rays
  !
  select case(lambdaoper_random_ray_method)
  case(0)
     !
     ! Simplest method: Each cell has the same weight
     !
     ! Find cell
     !
     icell = floor(ran2(iseed)*nrcells)+1
     if(icell.gt.nrcells) stop 9231
     !
     ! Random location in cell
     !
     call montecarlo_random_pos_in_cell(icell)
     !
     ! Find the direction in which to emit
     !
     call montecarlo_randomdir(ray_cart_dirx,ray_cart_diry,ray_cart_dirz)
     !
!
! Note: When we use weighted cells (rays emanating from some cells are
!       more likely than from others), then we must include an extra
!       counter-weighting of the rays to prevent improper weighting in
!       the calculation of the mean intensity.
  case default
     write(stdo,*) 'ERROR: Do not know lambda random ray method ',lambdaoper_random_ray_method
     stop
  end select
  !
  ! If "doshift" is set, then we shift the starting point of the ray
  !
  if(doshift) then
     !
     ! Compute the base shift
     !
     if(grid_contsph_r.eq.0.d0) then
        write(stdo,*) 'ERROR: The grid_contsph_r is not set. Did you forget to'
        write(stdo,*) '       call the subroutine postprocess_grid_for_use() after'
        write(stdo,*) '       setting up the grid?'
        stop
     endif
     lambdaoper_shift = 2.01d0*grid_contsph_r
     !
     ! If star_maxexterndist>0
     !
     if(star_maxexterndist.gt.0.d0) then
        lambdaoper_shift = lambdaoper_shift + 2*star_maxexterndist
     endif
     !
     ! Now shift
     !
     ray_cart_x = ray_cart_x - lambdaoper_shift * ray_cart_dirx
     ray_cart_y = ray_cart_y - lambdaoper_shift * ray_cart_diry
     ray_cart_z = ray_cart_z - lambdaoper_shift * ray_cart_dirz
     !
     ! The end s-value is the shift
     !
     ray_dsend = lambdaoper_shift
  else
     ray_dsend = 1d99
  endif
end subroutine lambdaoper_draw_random_ray


!------------------------------------------------------------------------
!                      RANDOM RAY LAMBDA OPERATOR
!------------------------------------------------------------------------
subroutine lambdaoper_operator_cellbased_random_rays(nrays)
  implicit none
  integer :: iray,nrays
  !
  ! Check
  !
  if(.not.allocated(lambdaoper_meanint_cell)) then
     write(*,*) 'ERROR: Cannot use Random-Ray Lambda Operator without'
     write(*,*) '       first calling lambdaoper_init() with vertexbased=.false..'
     stop
  endif
  !
  ! Clear the mean intensity array
  !
  lambdaoper_meanint_cell(:) = 0.d0
  lambdaoper_dstot_cell(:) = 0.d0
  !
  ! Do a loop over all rays
  !
  do iray=1,nrays
     !
     ! Find a random ray: This sets the ray_cart_x,y,z and ray_cart_dirx,y,z
     ! according to some recipe
     !
     call lambdaoper_draw_random_ray(.true.)
     !
     ! Now integrate along the ray
     !
     call lambdaoper_integrate_ray()
  enddo
  !
  ! Compute the mean intensity 
  ! 
  lambdaoper_meanint_cell(:) = lambdaoper_meanint_cell(:) / ( lambdaoper_dstot_cell(:) + 1d-99 )
  !
  ! Compute the diagonal elements of the Lambda Operator
  !
  lambdaoper_lambda_diag_cell(:) = lambdaoper_lambda_diag_cell(:) / ( lambdaoper_dstot_cell(:) + 1d-99 )
  !
end subroutine lambdaoper_operator_cellbased_random_rays


!========================================================================
!               LINE-TRANSFER LAMBDA OPERATOR STUFF
!========================================================================

end module lambdaoper_module
