module camera_module
  !$ use omp_lib
  use rtglobal_module
  use amrray_module
  use dust_module
  use lines_module
  use gascontinuum_module
  use stars_module
  use montecarlo_module
  use ioput_module
  use userdef_module
  use mathroutines_module
  use sources_module
  use constants_module
  !
  ! This module computes images using ray-tracing. There are two modes of
  ! viewing: one with the observer at infinity (good for most astrophysical
  ! purposes) and one with a local observer (good for public relations movies
  ! for instance). 
  !
  ! ----------------------------------------------------------------------------
  !
  ! Main settings
  !
  !   What do we truly transfer?
  !
  !       1 = The actual radiative transfer [default]
  !      -1 = The total column depth is 'traced' along each ray. 
  !      -2 = The total optical depth is 'traced' along each ray.
  !      -3 = The total optical depth is 'traced' along each ray, and
  !           the point is found where tau_tot-tau=1, i.e. where the
  !           optical depth toward the user is 1 at that frequency.
  !
  integer :: camera_tracemode=1
  !
  !   Do we do first-order integration of the transfer equation (=cell
  !   based, i.e. the default) or second order (=corner based, i.e. much
  !   smoother, and for line radiative transfer allows capturing 
  !   Doppler crossings).
  !
  logical :: camera_secondorder=.false.
  !
  !   Do we include the full Stokes vector in the ray-tracing?
  !
  logical :: camera_stokesvector=.false.
  !
  !   If we have very narrow lines (i.e. intrinsic line width smaller
  !   than cell-to-cell doppler shifts), then the ray-tracing may go wrong
  !   if no sub-stepping is done in the transfer equation when the line
  !   doppler-shifts through the frequency-of-sight. With this flag we can
  !   switch on this careful kind of integration of the transfer equation.
  !   This assures that the result of the ray-tracing is correct, even 
  !   in the presence of strong velocity gradients and narrow lines. 
  !
  !   NOTE: This does NOT assure that the pixel size of the image is 
  !         sufficiently small to catch the regions of the image where 
  !         the line is visible. In extreme cases (very narrow line, very
  !         strong velocity gradients) these regions can be very narrow,
  !         and thus easily missed with coarse pixel resolution. An example
  !         of such a case would be a keplerian disk around a massive object,
  !         in which T_disk <<< T_virial. See Pontoppidan et al. 2009. 
  !
  !   NOTE: We use the variable lines_widthmargin from the lines_module
  !         to determine the window around the line center within which
  !         we insist on doing sub-stepping should it become necessary.
  !
  logical :: camera_catch_doppler_line=.false.
  !
  !   For the doppler catching algorithm (see above): How fine velocity
  !   resolution do you wish for the sub-stepping? Also in same units.
  !   Currently set to 0.2, which is more than fine enough. Also here
  !   you may experiment with courser values (try 0.4, for instance).
  !   The default is the safe value.
  !
  double precision :: camera_catch_doppler_resolution=0.2d0
  !
  !   This is the main switch for the perspective of the viewer: if this
  !   variable is false, the observer-at-infinity mode is used. If this is
  !   true, the local-observer mode is used.
  !
  logical :: camera_localobserver=.false.
  integer :: camera_localobs_projection=1
  !
  !   If this integer is set to 1, then the stars are included in the images.
  !   In this version they are only included as point sources.
  !
  integer :: camera_incl_stars
  !
  !   Dump subpixeling diagnostics?
  !
  logical :: camera_diagnostics_subpix=.false.
  !
  ! For camera_tracemode=-3 we must set the camera_taustop=1d99 before
  ! the first trace, and then at the second trace it will stop at
  ! optical depth = camera_taustop, where taustop = tautot-tausurface 
  ! so that we can figure out where that point is.
  !
  double precision :: camera_tausurface=1.d0
  double precision, allocatable :: camera_taustop(:)
  double precision, allocatable :: camera_xstop(:),camera_ystop(:),camera_zstop(:)
  double precision, allocatable :: camera_dstop(:)
  double precision, allocatable :: camera_tausurface_x(:,:,:)
  double precision, allocatable :: camera_tausurface_y(:,:,:)
  double precision, allocatable :: camera_tausurface_z(:,:,:)
  !
  ! Variables and arrays
  !
  !   The camera-local frequency array
  !
  double precision, allocatable :: camera_frequencies(:)
  integer :: camera_nrfreq=0
  !
  !   The point to which the camera points is a 3-D point in space. This
  !   is used in both viewing modes. For the observer-at-infinity mode this
  !   is merely responsible for a shift of the camera pointing such that the
  !   center of the image plane is at the pointing position (but see 
  !   camera_zoomcenter_x/y below), while for the local-observer mode this
  !   position is used relative to the camera_observer_position (see below)
  !   to determine the camera viewing direction.
  !
  double precision :: camera_pointing_position(1:3)
  !
  !   The size of the image, from image center to image boundary. I.e. the
  !   total size (from edge to edge) is twice this value. 
  !
  double precision :: camera_image_halfsize_x
  double precision :: camera_image_halfsize_y
  !
  !   In the image plane you can specify where you center your 'ccd camera'.
  !   For the observer-at-infinity mode you can get a similar effect by
  !   adapting the camera_pointing_position. It can then nevertheless be 
  !   useful for zooming-in, not requiring any 3-D transformations to 
  !   convert a sub-field in an existing image into a new camera_pointing_position.
  !   For the local-observer mode the zoomcenter is in fact not reproducable
  !   with a different pointing position because of spherical aberrations. The
  !   zoomcenter stuff is done in the image plane, i.e. the location of the
  !   'ccd camera' in the focal plane of the 'telescope'.
  !
  double precision :: camera_zoomcenter_x
  double precision :: camera_zoomcenter_y
  !
  !   For local observer mode, if a fisheye projection (for planetarium
  !   dome, OMNIMAX) is chosen, his is the shift between the pointing
  !   direction and the zenith of the projection.
  !
  double precision :: camera_localobs_zenith
  !
  !   The orientation of the camera. The center of the image stays put, but
  !   the pixels of the camera are rotated around this point.
  !
  double precision :: camera_pointing_degr_posang
  double precision :: camera_pointing_cos_posang
  double precision :: camera_pointing_sin_posang
  !
  !   For the local-observer mode: this is the position of the observer.
  !   This has no meaning in the observer-at-infinity mode.
  !
  double precision :: camera_observer_position(1:3)
  !
  !   For the flux-conserving mode of imaging, this is the criterion for how
  !   much each pixel should be refined to assure flux conservation. This is
  !   typically of order 1, with smaller values leading to more accurate images
  !   (i.e. each pixel representing the true flux better), but also requiring
  !   more computational power. This does not have meaning for (future??) 
  !   Monte-Carlo style images.
  !
  double precision :: camera_refine_criterion
  double precision :: camera_spher_cavity_relres=0.05d0
!  double precision :: camera_min_aspectratio=0.05d0
  double precision :: camera_min_drr=0.003d0
  double precision :: camera_max_dangle=0.3d0
  double precision :: camera_min_dangle=0.05d0
  integer :: camera_subpixeling_npixfine,camera_subpixeling_npixtot
  !
  !   The number of pixels in x and y image-plane directions.
  !
  integer :: camera_image_nx,camera_image_ny
  !
  !   The number of pixels in r and phi image-plane directions (for circular images)
  !
  integer :: camera_image_nr,camera_image_nphi
  !
  !   Flags and variables for the setting of the camera_frequencies array
  !
  integer :: camera_theinu,camera_range_nlam
  double precision :: camera_lambdamic,camera_lambdamic1
  logical :: camera_loadcolor
  logical :: camera_loadlambdas
  logical :: camera_setfreq_global
  !
  !   For the flux-conserving capabilities of the ray-traced images this number
  !   gives the maximum number of times a pixel is allowed to recursively split
  !   into 2x2 sub-pixels. If this number is set to 0, then the ray-traced images
  !   are one-ray-per-pixel: this means that details in the model that are on
  !   scales smaller than the pixel can resolve may be missed entrely, and may
  !   lead to under-estimation of the flux, or if a very bright un-resolved 
  !   detail is hit by coincidence it may lead to over-estimation of the flux
  !   in this pixel. The flux-conserving scheme to avoid this is to recursively
  !   refine the pixel into sub-pixels and once the required resolution is
  !   achieved, add all contributions up (weighted by their sub-pixel size)
  !   to get the 'true' flux of the pixel which is then converted into an
  !   average intensity of the pixel.
  !
  integer :: camera_nrrefine 
  !
  !   For the observer-at-infinity mode these represent the angular position
  !   of the observer at infinity in units of degrees. If theta=0, then the
  !   observer views the stuff from the top (x=y=0, z=infinity). If theta=
  !   pi the observer views from the bottom (x=y=0, z=-infinity). If theta=
  !   pi/2 then the viewer is in the z=0 plane. Other values are obvious from
  !   these specifications. For values of 0<theta<pi the value of phi determines
  !   the location in the remaining degree of freedom. Example: if theta=pi/2
  !   and phi=0 then the observer is located at y=-infinity, x=0, z=0. If theta=pi/2
  !   and phi=pi/2 the observer is x=-infinity, y=0, z=0. 
  !
  double precision :: camera_observer_degr_theta=0.d0
  double precision :: camera_observer_degr_phi=0.d0
  !
  !   Precalculated cosine and sine of the above variables.
  !
  double precision :: camera_observer_cos_theta=0.d0
  double precision :: camera_observer_cos_phi=0.d0
  double precision :: camera_observer_sin_theta=0.d0
  double precision :: camera_observer_sin_phi=0.d0
  !
  !   The spectrum 
  !
  double precision, allocatable :: camera_spectrum_iquv(:,:)
  !
  !   Flag for use or no use of the aperture option for spectra
  ! 
  logical :: camera_use_aperture_info=.false.
  double precision :: camera_observer_distance_pc=0.d0
  integer :: camera_aperture_nf=0
  double precision, allocatable :: camera_aperture_freq(:)
  double precision, allocatable :: camera_aperture_radius_collectarea(:)
  !
  !   For multiple image vantages (optional). This is useful for making
  !   movies where the observer is changing position. Once can then simply
  !   load these positions into these arrays and do all the image rendering
  !   at once.
  !
  double precision, allocatable :: cameras_pt_pos(:,:)
  double precision, allocatable :: cameras_img_hs_x(:),cameras_img_hs_y(:)
  double precision, allocatable :: cameras_zmc_x(:),cameras_zmc_y(:)
  double precision, allocatable :: cameras_pt_degr_pa(:)
  double precision, allocatable :: cameras_obs_pos(:,:)
  double precision, allocatable :: cameras_obs_degr_th(:)
  double precision, allocatable :: cameras_obs_degr_ph(:)
  integer :: cameras_nr_images=0
  !
  !    Some arrays for the images and for rendering
  !
  double precision, allocatable :: camera_intensity_iquv(:,:)
  double precision, allocatable :: camera_rect_image_iquv(:,:,:,:)
  double precision, allocatable :: camera_circ_image_iquv(:,:,:,:)
  logical :: camera_warn_resolution
  !
  !    Flag to force RADMC-3D to precompute the source function for all frequencies
  !    of the camera_frequencies array at once.
  !
  logical :: camera_scatsrc_allfreq=.false.
  !
  !    Integer specifying if we are/were in a star sphere, and if yes,
  !    which one.
  !
  integer :: camera_istar=0
  !
  !    With how many pixels (in each direction) do we need to resolve the
  !    star(s) if we treat them as spheres?
  !
  double precision :: camera_starsphere_nrpix = 20
  !
  ! The data for the circular images
  !
  integer :: cim_nr,cim_np
  double precision, allocatable :: cim_rc(:),cim_ri(:),cim_pc(:),cim_pi(:)
  double precision, allocatable :: cim_surfarea(:,:)
  integer :: camera_circ_nrphiinf = 128
  integer :: camera_circ_nrext    = 10
  integer :: camera_circ_dbdr     = 1
  integer :: camera_circ_imethod  = 1
  integer :: camera_circ_nrref    = 1
  !
  !    A new mode for scattering images: the lambda single scattering mode. Default
  !    is 0 (i.e. normal Monte Carlo scattering). If set to 1, then use the
  !    do_lambda_starlight_single_scattering() routine instead of do_monte_carlo_scattering().
  !    
  integer :: camera_lambda_starlight_single_scat_mode = 0
  !
  ! For convenience (in particular when using aligned grains): the 
  ! line-of-sight direction vector and the polarization orientation vector
  !
  double precision :: camera_dir(1:3),camera_svec(1:3)
  !
  ! For 2-D axisymmetric models it might be important to limit the 
  ! length of a ray segment to avoid too large changes of the phi
  ! (azimuthal) angle of the ray between two cellwall crossings.
  !
  double precision :: camera_maxdphi = 0.d0
  !
  ! OpenMP Parallellization:
  ! Global variables used in subroutine calls within the parallel region which are threadprivate
  !
  !!!!!!$OMP THREADPRIVATE(camera_nrrefine)
  !$OMP THREADPRIVATE(camera_intensity_iquv)
contains


!-------------------------------------------------------------------------
!                   INITIALIZE THE CAMERA MODULE
!-------------------------------------------------------------------------
subroutine camera_init()
  use amr_module
  implicit none
  integer :: ierr,inu,ispec
  double precision :: temp
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
  ! Currently the aligned grains mode is not compatible with mirror
  ! symmetry mode in spherical coordinates. 
  !
  if((alignment_mode.ne.0).and.amrray_mirror_equator) then
     write(stdo,*) 'ERROR: Currently polarization is not compatible with'
     write(stdo,*) '       equatorial mirror symmetry mode. This means that'
     write(stdo,*) '       your theta-grid must not only cover the upper'
     write(stdo,*) '       quadrant (0<theta<=pi/2) but both the upper'
     write(stdo,*) '       and the lower quadrant (0<theta<pi).'
     stop
  endif
  !
  ! Currently the small angle scattering mode (including polarization) is
  ! not compatible with 2-D spherical coordinates.
  !
  if((scattering_mode.ge.2).and.(igrid_coord.ge.100).and. &
     (amr_dim.eq.1)) then
     write(stdo,*) 'ERROR: Non-isotropic scattering is incompatible with'
     write(stdo,*) '       1-D spherical coordinates.'
     stop
  endif
  !
  ! If you wish to use the doppler-catching of lines, then we must be in second-order mode.
  !
  if(camera_catch_doppler_line) then
     camera_secondorder = .true.
  endif
  !
  ! If you wish to use second order with lines in spherical coordinates, then 
  ! you must have full 3-D spherical coordinates
  !
  if((rt_incl_lines).and.(camera_secondorder).and. &
       (igrid_coord.ge.100).and.(amr_dim.ne.3)) then
     write(stdo,*) 'ERROR: Line radiative transfer with second order'
     write(stdo,*) '       integration is incompatible with'
     write(stdo,*) '       2-D spherical coordinates. Use 3-D spherical '
     write(stdo,*) '       coordinates instead (i.e. include, say, 64 or'
     write(stdo,*) '       more cells in phi-direction).'
     stop
  endif
  !
  ! Aligned grains are incompatible with second order integration 
  !
  if((alignment_mode.ne.0).and.(camera_secondorder)) then
     write(stdo,*) 'ERROR: Grain alignment is incompatible with second order integration.'
     stop
  endif
  !
  ! Aligned grains are incompatible with gas lines
  !
  if((alignment_mode.ne.0).and.(rt_incl_lines)) then
     write(stdo,*) 'ERROR: Grain alignment is incompatible with gas lines.'
     stop
  endif
  !
  ! Aligned grains require full Stokes
  !
  if((alignment_mode.ne.0).and.(.not.camera_stokesvector)) then
     write(stdo,*) 'ERROR: Grain alignment requires full Stokes vector mode.'
     stop
  endif
  !
  ! Aligned grains (for now) only guaranteed to work if 
  ! scattering_mode.ge.4, which automatically includes
  ! the full Stokes vectors in the Monte Carlo module.
  ! Maybe I will later do some checks if it can also
  ! work if MC does not use the full Stokes, but that
  ! is not so urgent.
  !
  if((alignment_mode.ne.0).and.(scattering_mode.lt.4)) then
     write(stdo,*) 'ERROR: Grain alignment requires (at least for now) scattering_mode >= 4.'
     stop
  endif
  !
  ! For now, alignment does not yet work with 2-D axisymmetric
  ! models
  !
  if((alignment_mode.ne.0).and.(igrid_coord.ge.100).and.(amr_dim.ne.3)) then
     write(stdo,*) 'ERROR: Grain alignment does not yet work in 2-D axisymmetry.'
     stop
  endif
  !
  ! Check that the camera frequency array or choice of frequencies from the
  ! global frequency array is made.
  !
  if((camera_nrfreq.le.0).or.(.not.allocated(camera_frequencies))) then
     write(stdo,*) 'ERROR in camera module: cannot init camera if '
     write(stdo,*) '      the camera_frequencies(:) array is not set.'
     stop
  endif
  !
  ! If already initialized before, deallocate some stuff
  !
  call camera_partial_cleanup()
  !
  ! Check some stuff
  !
  if((camera_image_nx.le.0).or.(camera_image_ny.le.0)) then 
     if((camera_image_nr.le.0).or.(camera_image_nphi.le.0)) then
        write(stdo,*) 'ERROR in camera module: Zero number of pixels...'
        stop
     endif
  endif
  if(rt_incl_dust) then
     if(.not.allocated(dust_kappa_abs)) then
        write(stdo,*) 'ERROR in camera module: dust opacities not yet set.'
        stop
     endif
     if(.not.allocated(dustdens)) then
        write(stdo,*) 'ERROR in camera module: dust density array not yet allocated.'
        stop
     endif
     if(.not.allocated(dusttemp)) then
        write(stdo,*) 'ERROR in camera module: dust temperature array not yet allocated.'
        stop
     endif
  else
     if(scattering_mode.ne.0) then
        write(stdo,*) 'ERROR in camera module: if dust not included, cannot have scattering_mode.ne.0.'
        stop
     endif
  endif
  if(rt_incl_lines) then
     if(.not.allocated(lines_aud)) then
        write(stdo,*) 'ERROR in camera module: lines included, but no line data read.'
        stop
     endif
     !
     ! ************ DO MORE CHECKS **************
     !
  endif
  !
  ! Call the sources module initialization
  ! 
  call sources_init(camera_nrfreq,camera_frequencies,camera_secondorder,camera_catch_doppler_line)
  !
  ! Now allocate further arrays
  !
  allocate(camera_spectrum_iquv(1:camera_nrfreq,1:4),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in camera module: Could not allocate spectrum array.'
     stop
  endif
  !$OMP PARALLEL
  allocate(camera_intensity_iquv(1:camera_nrfreq,1:4),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in camera module: Could not allocate camera_intensity_iquv() array'
     stop
  endif
  !$OMP END PARALLEL
  !
  ! Now allocate the image array for the rectangular images
  !
  if((camera_image_nx.ge.1).and.(camera_image_ny.ge.1)) then
     allocate(camera_rect_image_iquv(1:camera_image_nx,1:camera_image_ny,1:camera_nrfreq,1:4),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in camera module: Could not allocate camera_rect_image_iquv() array'
        stop
     endif
     camera_rect_image_iquv(:,:,:,:) = 0.d0
  endif
  !
  ! Now allocate the image array for the circular images
  !
  if((camera_image_nr.ge.1).and.(camera_image_nphi.ge.1)) then
     allocate(camera_circ_image_iquv(0:camera_image_nr,1:camera_image_nphi,1:camera_nrfreq,1:4),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in camera module: Could not allocate camera_circ_image_iquv() array'
        stop
     endif
     camera_circ_image_iquv(:,:,:,:) = 0.d0
  endif
  !
  ! Some further comments:
  !
  ! If using second order integration in spherical coordinates, then
  ! better set the camera_maxdphi to a reasonably small value.
  !
  if((igrid_coord.ge.100).and.(igrid_coord.lt.200).and.camera_secondorder) then
     if(camera_maxdphi.eq.0d0) then
        write(stdo,*) 'WARNING: camera_maxdphi set to 0. Could be dangerous. Better set it to e.g. 0.1.'
     elseif(camera_maxdphi.le.0.01d0) then
        write(stdo,*) 'WARNING: camera_maxdphi set to very small value >0. ',&
             'This is OK, but it may make the code slow.'
     endif
  endif
  !
end subroutine camera_init



!-------------------------------------------------------------------------
!                     FINALIZE THE CAMERA MODULE
!-------------------------------------------------------------------------
subroutine camera_cleanup()
  implicit none
  !
  ! Clean up the private arrays
  !
  call camera_partial_cleanup()
  !
  ! Clean up all the rest. These are the arrays that may have 
  ! to be set by the user or elsewhere in the program.
  !
  ! ...The wavelength specification arrays
  !
  if(allocated(camera_frequencies)) deallocate(camera_frequencies)
  camera_nrfreq = 0
  !
  ! ...The movie-related arrays
  ! 
  if(allocated(cameras_pt_pos)) deallocate(cameras_pt_pos)
  if(allocated(cameras_img_hs_x)) deallocate(cameras_img_hs_x)
  if(allocated(cameras_img_hs_y)) deallocate(cameras_img_hs_y)
  if(allocated(cameras_zmc_x)) deallocate(cameras_zmc_x)
  if(allocated(cameras_zmc_y)) deallocate(cameras_zmc_y)
  if(allocated(cameras_pt_degr_pa)) deallocate(cameras_pt_degr_pa)
  if(allocated(cameras_obs_pos)) deallocate(cameras_obs_pos)
  if(allocated(cameras_obs_degr_th)) deallocate(cameras_obs_degr_th)
  if(allocated(cameras_obs_degr_ph)) deallocate(cameras_obs_degr_ph)
  !
  ! ...The aperture-related arrays
  !
  if(allocated(camera_aperture_freq)) deallocate(camera_aperture_freq)
  if(allocated(camera_aperture_radius_collectarea)) deallocate(camera_aperture_radius_collectarea)
  !
  ! ...Any circular image stuff
  !
  if(allocated(cim_rc)) deallocate(cim_rc)
  if(allocated(cim_ri)) deallocate(cim_ri)
  if(allocated(cim_pc)) deallocate(cim_pc)
  if(allocated(cim_pi)) deallocate(cim_pi)
  if(allocated(cim_surfarea)) deallocate(cim_surfarea)
  !
end subroutine camera_cleanup


!-------------------------------------------------------------------------
!               RESET THE CAMERA MODULE FOR NEXT ACTION
!
! Here we only clean up the arrays that contain temporary information
! and the arrays that contain results. All these arrays can be deleted
! at the start of a new action without loss of important information.
!-------------------------------------------------------------------------
subroutine camera_partial_cleanup()
  implicit none
  call sources_partial_cleanup()
  if(allocated(camera_rect_image_iquv)) deallocate(camera_rect_image_iquv)
  if(allocated(camera_circ_image_iquv)) deallocate(camera_circ_image_iquv)
  if(allocated(camera_spectrum_iquv)) deallocate(camera_spectrum_iquv)
  !$OMP PARALLEL
  if(allocated(camera_intensity_iquv)) deallocate(camera_intensity_iquv)
  !$OMP END PARALLEL
  if(allocated(camera_xstop)) deallocate(camera_xstop)
  if(allocated(camera_ystop)) deallocate(camera_ystop)
  if(allocated(camera_zstop)) deallocate(camera_zstop)
  if(allocated(camera_dstop)) deallocate(camera_dstop)
  if(allocated(camera_taustop)) deallocate(camera_taustop)
  if(allocated(camera_tausurface_x)) deallocate(camera_tausurface_x)
  if(allocated(camera_tausurface_y)) deallocate(camera_tausurface_y)
  if(allocated(camera_tausurface_z)) deallocate(camera_tausurface_z)
end subroutine camera_partial_cleanup


!-------------------------------------------------------------------------
!            SET THE STARTING POINT AND DIRECTION OF RAY IN IMAGE
!
! Note: this subroutine will also set the S-vector for polarization,
!       which is a variable in the rtglobal_module.f90 module.
!-------------------------------------------------------------------------
subroutine camera_set_ray(px,py,x,y,z,dirx,diry,dirz,distance)
  implicit none
  double precision :: px,py
  double precision :: x,y,z
  double precision :: xbk,ybk,zbk
  double precision :: dirx,diry,dirz
  double precision :: distance,dum,shift,dumy,dumz,cs,sn
  !
  ! First find out how far we should push the starting points of the ray
  ! back to be absolutely sure that the start behind the object.
  !
  if(camera_localobserver) then
     shift = sqrt( camera_observer_position(1)**2 + &
                   camera_observer_position(2)**2 + &
                   camera_observer_position(3)**2 ) 
  else
     shift = sqrt( camera_pointing_position(1)**2 + &
                   camera_pointing_position(2)**2 + &
                   camera_pointing_position(3)**2 ) 
  endif
  !
  ! If star_maxexterndist>0
  !
  if(star_maxexterndist.gt.0.d0) then
     shift = shift + 2*star_maxexterndist
  endif
  !
  ! Grid dependent stuff
  !
  if(igrid_type.lt.100) then
     !
     ! Regular and AMR grids
     !
     if(igrid_coord.lt.100) then
        !
        ! Cartesian coordinates
        !
        if(igrid_coord.eq.10) then
           if(camera_observer_cos_theta.eq.0.d0) then
              write(stdo,*) 'ERROR: In plane-parallel grid you cannot choose'
              write(stdo,*) '       cos(theta)=0 for the camera.'
              stop
           endif
           shift = shift + sqrt(3.d0)*max(abs(amr_grid_xi(1,3)),&
                            abs(amr_grid_xi(amr_grid_nz+1,3)))/ &
                            abs(camera_observer_cos_theta)
        elseif(igrid_coord.eq.20) then
           if(camera_observer_cos_theta.eq.0.d0) then
              write(stdo,*) 'ERROR: In pencil-parallel grid you cannot choose'
              write(stdo,*) '       cos(theta)=0 for the camera.'
              stop
           endif
           shift = shift + sqrt(3.d0)*max(abs(amr_grid_xi(1,2)),&
                            abs(amr_grid_xi(amr_grid_ny+1,2)),  &
                            abs(amr_grid_xi(1,3)),              &
                            abs(amr_grid_xi(amr_grid_nz+1,3)))/ &
                            abs(camera_observer_cos_theta)
        else
           shift = shift + sqrt(3.d0)*max(abs(amr_grid_xi(1,1)),&
                                       abs(amr_grid_xi(amr_grid_nx+1,1)),&
                                       abs(amr_grid_xi(1,2)),&
                                       abs(amr_grid_xi(amr_grid_ny+1,2)),&
                                       abs(amr_grid_xi(1,3)),&
                                       abs(amr_grid_xi(amr_grid_nz+1,3)))
        endif
     elseif(igrid_coord.lt.200) then
        !
        ! Spherical coordinates
        !
        shift = shift + amr_grid_xi(amr_grid_nx+1,1)
     else
        stop 612
     endif
  else
     write(stdo,*) 'ERROR in camera module: non-regular grids not yet ready'
     stop
  endif
  shift = shift*3        ! For safety, put it even further back
  !
  ! Switch between perspective view from a local observer (if localobserver==.true.),
  ! or from an observer 'infinitely far' away (if localobserver==.false.)
  !
  if(camera_localobserver) then
     !
     ! Local perspective: an observer close to (or inside) the object.
     ! Use the camera_observer_position() array to find the rotation
     ! and other things. Moreover, px and py are now to be interpreted
     ! as angular coordinates (in radian).
     !
     ! Check
     !
     if(camera_stokesvector) then
        write(stdo,*) 'ERROR: Cannot use full Stokes vector in local observer mode...'
        stop
     endif
     !
     ! First make the position, which is at the observer at the moment
     ! 
     x   = camera_observer_position(1)
     y   = camera_observer_position(2)
     z   = camera_observer_position(3)
     !
     ! Now make the 3-D direction vector, depending on which
     ! projection type we use.
     !
     if(camera_localobs_projection.eq.1) then
        !
        ! Simply project onto a flat screen
        !
        dirx = -px
        diry = -py
        dum  = dirx**2 + diry**2
        if(dum.gt.10.d0) then
           write(stdo,*) 'ERROR in camera module: angular coordinates of image'
           write(stdo,*) '      are out of range.'
           stop
        endif
        dum  = 1.d0/sqrt(dum + 1.d0)
        dirx = dirx*dum
        diry = diry*dum
        dirz = dum
     elseif(camera_localobs_projection.eq.2) then
        !
        ! Project onto a sphere: the center of the image is one point on the
        ! sphere. The px and py now point the way along the sphere to the
        ! direction we want to have. 
        !
        dum  = sqrt(px**2+py**2)
        if(dum.gt.0.d0) then
           dirx = -px/dum
           diry = -py/dum
           dirz = cos(dum)
           dum  = sin(dum)
           dirx = dirx*dum
           diry = diry*dum
        else
           dirx = 0.d0
           diry = 0.d0
           dirz = 1.d0
        endif
        !
     else
        write(stdo,*) 'ERROR in camera module: projection type ',camera_localobs_projection,' not known.'
        stop
     endif
     !
     ! Apply a rotation to allow the zenith of the projection to be above 
     ! the pointing direction (so that the viewer in a dome does not have
     ! to always look straight up).
     !
     cs   = cos(camera_localobs_zenith/57.295780d0)
     sn   = sin(camera_localobs_zenith/57.295780d0)
     dumy = diry
     dumz = dirz
     diry =  cs * dumy - sn * dumz
     dirz =  sn * dumy + cs * dumz
     !
     ! The 'position angle' is such that positive z is 'up' and negative z
     ! is 'down'. With that definition the camera is always oriented
     ! 'horizontally' as one would hold a handheld camera. Only with 
     ! the camera_pointing_degr_posang one can change this.
     !
     ! The exact distance for the ray to travel is precisely 'shift'
     !
     distance = shift
     !
  else
     !
     ! Observer at infinite distance (usually OK for astronomical observations).
     ! Use the angles theta and phi to rotate.
     !
     ! First make the starting position in the x-y plane, with z=0
     !
     x = px
     y = py
     z = 0.d0
     !
     ! Now make the direction vector
     !
     dirx = 0.d0
     diry = 0.d0
     dirz = 1.d0   ! Note: direction of light propatation, i.e. toward observer
     !
     ! Set distance to observer to infinity (=1d99)
     !
     distance = 1d99
     !
     ! For the camera at infinity we must, in addition to the direction
     ! vectors (see below) also rotate the x,y,z positions, while for the
     ! perspective view we only rotate the direction vector. 
     !
     ! NOTE: We use the right-hand-rule for the orientation of the 
     !       x,y,z axes. So if you look at the x-y plane from the top,
     !       with the z-axis pointing toward you, then the x and y 
     !       axes are pointing in the usual way: positive x toward the
     !       right, positive y toward the top. Note that this is also
     !       the convention for the images.
     !
     ! First a rotation of the camera itself. The camera is rotated
     ! in clockwise direction. Therefore, on the camera CCD any image is 
     ! rotated in counter-clockwise direction. 
     !
     xbk = x
     ybk = y
     x   = camera_pointing_cos_posang * xbk + camera_pointing_sin_posang * ybk
     y   =-camera_pointing_sin_posang * xbk + camera_pointing_cos_posang * ybk
     !
     ! Then a rotation in the y,z plane, i.e. going from pole-on
     ! to more face-on (if a disk is assumed to be present in the x-y plane)
     !
     ybk = y
     zbk = z
     y   = camera_observer_cos_theta * ybk - camera_observer_sin_theta * zbk
     z   = camera_observer_sin_theta * ybk + camera_observer_cos_theta * zbk
     !
     ! Then a rotation in the x,y plane, rotating horizontally around
     ! the object. If the camera looks toward the object, then the
     ! camera now moves to the left (clockwise around the object).
     !
     xbk = x
     ybk = y
     x   = camera_observer_cos_phi * xbk + camera_observer_sin_phi * ybk
     y   =-camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
     !
     ! Now shift the entire stuff toward the point of focus
     !
     x   = x + camera_pointing_position(1)
     y   = y + camera_pointing_position(2)
     z   = z + camera_pointing_position(3)
     !
  endif
  !
  ! Now perform the rotation of the direction vectors
  ! 
  ! NOTE: We use the right-hand-rule for the orientation of the 
  !       x,y,z axes. So if you look at the x-y plane from the top,
  !       with the z-axis pointing toward you, then the x and y 
  !       axes are pointing in the usual way: positive x toward the
  !       right, positive y toward the top. Note that this is also
  !       the convention for the images.
  !
  ! First a rotation of the camera itself. The camera is rotated
  ! in clockwise direction. Therefore, on the camera CCD any image is 
  ! rotated in counter-clockwise direction. Note: this is only
  ! necessary for the direction vectors in the perspective mode, hence
  ! the if-statement.
  !
  if(camera_localobserver) then
     xbk  = dirx
     ybk  = diry
     dirx = camera_pointing_cos_posang * xbk + camera_pointing_sin_posang * ybk
     diry =-camera_pointing_sin_posang * xbk + camera_pointing_cos_posang * ybk
  endif
  !
  ! Then a rotation in the y,z plane, i.e. going from pole-on
  ! to more face-on (if a disk is assumed to be present in the x-y plane)
  !
  ybk  = diry
  zbk  = dirz
  diry = camera_observer_cos_theta * ybk - camera_observer_sin_theta * zbk
  dirz = camera_observer_sin_theta * ybk + camera_observer_cos_theta * zbk
  !
  ! Then a rotation in the x,y plane, rotating horizontally around
  ! the object. If the camera looks toward the object, then the
  ! camera now moves to the left (clockwise around the object).
  !
  xbk  = dirx
  ybk  = diry
  dirx = camera_observer_cos_phi * xbk + camera_observer_sin_phi * ybk
  diry =-camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
  if(abs(dirx**2+diry**2+dirz**2-1.d0).gt.1d-6) then
     write(stdo,*) 'ERROR in camera module: direction vector not OK.'
     write(stdo,*) dirx,diry,dirz
     stop
  endif
  !
  ! If polarization, then set the S-vector
  !
  if(camera_stokesvector) then
     !
     ! Set the initial vector
     !
     ray_cart_svec(1) = 0.d0
     ray_cart_svec(2) = 1.d0    ! The S-vector points in the positive y-direction in the image
     ray_cart_svec(3) = 0.d0
     !
     ! Rotate camera along its axis
     !
     xbk  = ray_cart_svec(1)
     ybk  = ray_cart_svec(2)
     ray_cart_svec(1) = camera_pointing_cos_posang * xbk + camera_pointing_sin_posang * ybk
     ray_cart_svec(2) =-camera_pointing_sin_posang * xbk + camera_pointing_cos_posang * ybk
     !
     ! Then a rotation in the y,z plane, i.e. going from pole-on
     ! to more face-on (if a disk is assumed to be present in the x-y plane)
     !
     ybk  = ray_cart_svec(2)
     zbk  = ray_cart_svec(3)
     ray_cart_svec(2) = camera_observer_cos_theta * ybk - camera_observer_sin_theta * zbk
     ray_cart_svec(3) = camera_observer_sin_theta * ybk + camera_observer_cos_theta * zbk
     !
     ! Then a rotation in the x,y plane, rotating horizontally around
     ! the object. If the camera looks toward the object, then the
     ! camera now moves to the left (clockwise around the object).
     !
     xbk  = ray_cart_svec(1)
     ybk  = ray_cart_svec(2)
     ray_cart_svec(1) = camera_observer_cos_phi * xbk + camera_observer_sin_phi * ybk
     ray_cart_svec(2) =-camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
     if(abs(ray_cart_svec(1)**2+ray_cart_svec(2)**2+ray_cart_svec(3)**2-1.d0).gt.1d-6) then
        write(stdo,*) 'ERROR in camera module: S-vector not OK.'
        write(stdo,*) ray_cart_svec(1:3)
        stop
     endif
  endif
  !
  ! Now shift
  !
  x = x - shift * dirx
  y = y - shift * diry
  z = z - shift * dirz
  !
end subroutine camera_set_ray


!-------------------------------------------------------------------------
!            SET THE STARTING POINT AND DIRECTION FOR THE STARS
!
! Note: this subroutine will also set the S-vector for polarization,
!       which is a variable in the rtglobal_module.f90 module.
!-------------------------------------------------------------------------
subroutine camera_set_ray_stars_pntsrc(istar,x,y,z,dirx,diry,dirz,px,py,distance)
  implicit none
  integer :: istar
  double precision :: x,y,z
  double precision :: xbk,ybk,zbk
  double precision :: dirx,diry,dirz
  double precision :: px,py
  double precision :: distance
  !
  ! Check
  !
  if(.not.allocated(star_pos)) then
     write(stdo,*) 'ERROR in camera module: star positions not allocated'
     stop
  endif
  if((istar.lt.1).or.(istar.gt.nstars)) then
     write(stdo,*) 'ERROR in camera module: istar out of range'
     stop
  endif
  !
  ! The starting position is clearly the star position
  !
  ! NOTE: The stars are treated as point sources!!!
  !
  x   = star_pos(1,istar)
  y   = star_pos(2,istar)
  z   = star_pos(3,istar)
  !
  ! Check which perspective we use
  !
  if(camera_localobserver) then
     !
     ! If the local perspective is used
     !
     ! Check
     !
     if(camera_stokesvector) then
        write(stdo,*) 'ERROR: Cannot use Stokes vector in local observer mode...'
        stop
     endif
     !
     ! Make the direction vector, pointing toward the observer
     ! At the same time also compute the distance
     !
     dirx  = camera_observer_position(1) - star_pos(1,istar)
     diry  = camera_observer_position(2) - star_pos(2,istar)
     dirz  = camera_observer_position(3) - star_pos(3,istar)
     distance = sqrt( dirx**2 + diry**2 + dirz**2 )
     dirx  = dirx / distance
     diry  = diry / distance
     dirz  = dirz / distance
     !
     ! Find the location on the image
     !
     ! A rotation in the x,y plane, rotating horizontally around
     ! the object. 
     !
     xbk   = -dirx
     ybk   = -diry
     px    = camera_observer_cos_phi * xbk - camera_observer_sin_phi * ybk
     py    = camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
     !
     ! Then a rotation in the y,z plane
     !
     ybk   = py
     zbk   = -dirz
     py    = camera_observer_cos_theta * ybk + camera_observer_sin_theta * zbk
     !
  else
     !
     ! If observer-at-infinity is used
     !
     ! Make the direction vector
     !
     dirx = 0.d0
     diry = 0.d0
     dirz = 1.d0
     ybk  = diry
     zbk  = dirz
     diry = camera_observer_cos_theta * ybk - camera_observer_sin_theta * zbk
     dirz = camera_observer_sin_theta * ybk + camera_observer_cos_theta * zbk
     xbk  = dirx
     ybk  = diry
     dirx = camera_observer_cos_phi * xbk + camera_observer_sin_phi * ybk
     diry =-camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
     if(abs(dirx**2+diry**2+dirz**2-1.d0).gt.1d-6) then
        write(stdo,*) 'ERROR in camera module: direction vector not OK.'
        write(stdo,*) dirx,diry,dirz
        stop
     endif
     !
     ! If polarization, then set the S-vector
     !
     if(camera_stokesvector) then
        !
        ! Set the initial vector
        !
        ray_cart_svec(1) = 0.d0
        ray_cart_svec(2) = 1.d0    ! The S-vector points in the positive y-direction in the image
        ray_cart_svec(3) = 0.d0
        !
        ! Rotate camera along its axis
        !
        xbk  = ray_cart_svec(1)
        ybk  = ray_cart_svec(2)
        ray_cart_svec(1) = camera_pointing_cos_posang * xbk + camera_pointing_sin_posang * ybk
        ray_cart_svec(2) =-camera_pointing_sin_posang * xbk + camera_pointing_cos_posang * ybk
        !
        ! Then a rotation in the y,z plane, i.e. going from pole-on
        ! to more face-on (if a disk is assumed to be present in the x-y plane)
        !
        ybk  = ray_cart_svec(2)
        zbk  = ray_cart_svec(3)
        ray_cart_svec(2) = camera_observer_cos_theta * ybk - camera_observer_sin_theta * zbk
        ray_cart_svec(3) = camera_observer_sin_theta * ybk + camera_observer_cos_theta * zbk
        !
        ! Then a rotation in the x,y plane, rotating horizontally around
        ! the object. If the camera looks toward the object, then the
        ! camera now moves to the left (clockwise around the object).
        !
        xbk  = ray_cart_svec(1)
        ybk  = ray_cart_svec(2)
        ray_cart_svec(1) = camera_observer_cos_phi * xbk + camera_observer_sin_phi * ybk
        ray_cart_svec(2) =-camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
        if(abs(ray_cart_svec(1)**2+ray_cart_svec(2)**2+ray_cart_svec(3)**2-1.d0).gt.1d-6) then
           write(stdo,*) 'ERROR in camera module: S-vector not OK.'
           write(stdo,*) ray_cart_svec(1:3)
           stop
        endif
     endif
     !
     ! Find the location on the image
     !
     ! A rotation in the x,y plane, rotating horizontally around
     ! the object. 
     !
     xbk   = x - camera_pointing_position(1)
     ybk   = y - camera_pointing_position(2)
     px    = camera_observer_cos_phi * xbk - camera_observer_sin_phi * ybk
     py    = camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
     !
     ! Then a rotation in the y,z plane
     !
     ybk   = py
     zbk   = z - camera_pointing_position(3)
     py    = camera_observer_cos_theta * ybk + camera_observer_sin_theta * zbk
     !
     ! The distance is infinite
     ! 
     distance = 1d99
  endif
  !
end subroutine camera_set_ray_stars_pntsrc


!-------------------------------------------------------------------------
!                    RAY TRACE THROUGH THE MODEL
!
! This is the normal ray-tracing routine (see also the subroutine
! camera_ray1d_raytrace which is an alternative). 
!-------------------------------------------------------------------------
subroutine camera_serial_raytrace(nrfreq,inu0,inu1,x,y,z,dx,dy,dz,distance,   &
                           celldxmin,intensity)
  implicit none
  double precision :: x,y,z,dx,dy,dz,dum,distance,celldxmin
  double precision :: r,theta,phi
  double precision :: intensity(1:nrfreq,1:4)
  integer :: ispec,nrfreq,deepestlevel
  integer :: inu,levelnext,itemplate,inu0,inu1,iddr
  logical :: arrived
  double precision :: xp,xp1,dummy,dummy2,epstau
  double precision :: dumm,dumr,dumt,dump,dumx,dumy,dumz,thscale
  double precision :: src(1:nrfreq,1:4),alp(1:nrfreq)
  type(amr_branch), pointer :: a,b
  doubleprecision :: bxi(1:2,1:3),bxc(1:3)
  doubleprecision :: anu_prev,snu_prev(1:4),jnu_prev(1:4)
  doubleprecision :: anu_curr,snu_curr(1:4),jnu_curr(1:4)
  doubleprecision :: ln_anu_prev,ln_jnu_prev,ln_anu_curr,ln_jnu_curr
  doubleprecision :: dtau1,e0,e1,ax,bx,proj
  doubleprecision :: nstp,nu,nu0,nu0_prev,nu0_curr
  doubleprecision :: ap,jp(1:4),sp(1:4),ac,jc(1:4),sc(1:4),theomax(1:4),qdr(1:4)
  doubleprecision :: duc,eps,eps1,eps0,temp,resol,margin
  integer :: nsteps,iline,ilactive,istep,idir
  integer :: ixx,iyy,izz,bc_idir,bc_ilr
  double precision :: xbk,ybk,zbk,spx,spy,spz,rx1,ry1,rz1,rx0,ry0,rz0
  !
  ! Do checks
  !
  if(nrfreq.ne.camera_nrfreq) stop 9301
  if(debug_check_all.eq.1) then
     dum = sqrt(dx**2 + dy**2 + dz**2)
     if(abs(dum-1.d0).gt.1d-6) stop 3487
     if(.not.rt_incl_dust) then
        write(stdo,*) 'ERROR in camera_serial_raytrace: dust is not activated...'
        stop
     endif
!     if(.not.allocated(camera_ray_jnu)) stop 4560
!     if(.not.allocated(camera_ray_alpnu)) stop 4561
  endif
  if(camera_secondorder.and.(inu0.ne.inu1)) then
     write(stdo,*) 'ERROR: With second order integration it is impossible to '
     write(stdo,*) 'raytrace multiple frequencies at the same time'
     stop
  endif
  !
  ! Set defaults
  !
  anu_prev      = 0.d0
  snu_prev(1:4) = 0.d0
  jnu_prev(1:4) = 0.d0
  anu_curr      = 0.d0
  snu_curr(1:4) = 0.d0
  jnu_curr(1:4) = 0.d0
  if(camera_catch_doppler_line) then
     ln_anu_prev = 0.d0
     ln_jnu_prev = 0.d0
     ln_anu_curr = 0.d0
     ln_jnu_curr = 0.d0
     sources_local_doppler_prev = 0.d0
     sources_local_turb_prev = 0.d0
     sources_local_temp_prev = 0.d0
     sources_local_line_nup_prev(:) = 0.d0
     sources_local_line_ndown_prev(:) = 0.d0
     sources_local_doppler_curr = 0.d0
     sources_local_turb_curr = 0.d0
     sources_local_temp_curr = 0.d0
     sources_local_line_nup_curr(:) = 0.d0
     sources_local_line_ndown_curr(:) = 0.d0
     sources_local_doppler_end = 0.d0
     sources_local_turb_end = 0.d0
     sources_local_temp_end = 0.d0
     sources_local_line_nup_end(:) = 0.d0
     sources_local_line_ndown_end(:) = 0.d0
  endif
  !
  ! Check the dust scattering
  !
  if(scattering_mode.ge.1) then
     if(scattering_mode.eq.1) then
        !
        ! Isotropic scattering mode
        !
        if(mcscat_nrdirs.ne.1) then
           if(mcscat_nrdirs.eq.0) then
              write(stdo,*) 'ERROR: For isotropic scattering mode we must have'
              write(stdo,*) '       the scattering source function computed.'
              stop
           else
              write(stdo,*) 'ERROR: For isotropic scattering mode we need not'
              write(stdo,*) '       (and should not) use multiple directions.'
              stop
           endif
        endif
     elseif(scattering_mode.ge.2) then
        !
        ! Non-isotropic scattering mode. 
        !
        ! Check if this mode is allowed...
        !
        if(igrid_mirror.ne.0) then
           write(stdo,*) 'ERROR: If mirror symmetry is switched on (which is'
           write(stdo,*) '       automatically done if the theta coordinates'
           write(stdo,*) '       are all <=pi/2 and the last theta =pi/2)'
           write(stdo,*) '       then non-isotropic scattering is not (yet)'
           write(stdo,*) '       available.'
           stop
        endif
        if(amr_dim.eq.1) then
           if((igrid_coord.ge.100).and.(igrid_coord.lt.200)) then
              write(stdo,*) 'ERROR: Non-isotropic scattering in spherical coordinates'
              write(stdo,*) '       is (for now) only available in 2-D or 3-D.'
              stop
           endif
        endif
        !
        ! Check if the directions are
        ! set properly and that our current direction is equal to that
        ! of the monte carlo pre-computed source function.
        ! Note: Use mcscat_current_dir to select among the pre-computed
        !       discrete directions. 
        !
        if(mcscat_current_dir.gt.mcscat_nrdirs) stop 3012
        if(.not.allocated(mcscat_dirs)) stop 3013
        if((abs(mcscat_dirs(1,mcscat_current_dir)-dx).gt.1d-6).or.     &
           (abs(mcscat_dirs(2,mcscat_current_dir)-dy).gt.1d-6).or.     &
           (abs(mcscat_dirs(3,mcscat_current_dir)-dz).gt.1d-6)) then
           write(stdo,*) 'ERROR in camera module: attempting to use'
           write(stdo,*) '      scattering source for a direction that'
           write(stdo,*) '      is not the direction of the selected '
           write(stdo,*) '      precomputed scattering source function.'
           stop
        endif
     endif
     !
     ! General checks for scattering
     !
     if(.not.allocated(mcscat_scatsrc_iquv)) stop 3005
     if(camera_mcscat_monochromatic) then
        !
        ! Scattering is done one frequency at a time, so the scattering
        ! source function has only one bin in frequency space.
        !
        ! Check that inu1.eq.inu0, because if not, then this method cannot
        ! work properly.
        !
        if(inu0.ne.inu1) then
           write(stdo,*) 'INTERNAL ERROR in ray-tracing with scattering:'
           write(stdo,*) '      If you use monochromatic scattering source '
           write(stdo,*) '      function, then you can ray trace only one'
           write(stdo,*) '      frequency at a time. However, I obtained'
           write(stdo,*) '      inu0=',inu0,' and inu1=',inu1
           stop
        endif
        !
        ! Check that mc_nrfreq.eq.1
        !
        if(mc_nrfreq.ne.1) then
           write(stdo,*) 'INTERNAL ERROR in ray-tracing with scattering:'
           write(stdo,*) '      mc_nrfreq.ne.1 while the monochromatic'
           write(stdo,*) '      style of scattering Monte Carlo is used.'
           stop
        endif
        !
     else
        !
        ! Scattering is done spectrally dispersed, i.e. the scattering 
        ! source function has been computed for all camera frequencies
        ! beforehand and stored for all these frequencies in presumably
        ! a very large array. This method is presumably not so useful
        ! for large models, because it would require too much memory
        ! space. But it may be faster than the frequency-by-frequency
        ! method (i.e. the monochromatic method, see above).
        !
        if(mc_nrfreq.ne.camera_nrfreq) then
           write(stdo,*) 'ERROR in camera module when using scattering: '
           write(stdo,*) '      the frequency grid of the camera is not'
           write(stdo,*) '      equal to that of the montecarlo...' 
           stop
        endif
        do inu=1,camera_nrfreq
           if(camera_frequencies(inu).ne.mc_frequencies(inu)) then
              write(stdo,*) 'ERROR in camera module when using scattering: '
              write(stdo,*) '      the frequency grid of the camera is not'
              write(stdo,*) '      equal to that of the montecarlo...' 
              stop
           endif
        enddo
     endif
  endif
  !
  ! Set starting value of celldxmin, i.e. the size of the smallest cell
  ! we will pass through. Note that for cartesian coordinates, since cell 
  ! sizes will then always be cubes of the base size times (1/2)^level
  ! we can use just the deepest level to compute celldxmin.
  !
  ! NOTE: In the AMR module we now changed the level-numbering: Now
  !       the base grid level is level=0 (used to be 1). Hence we must
  !       now set deepestlevel = 0 here, and the multiplication factor
  !       is now (1/2)^level instead of (1/2)^(level-1).
  !
  celldxmin     = 1d99
  deepestlevel  = 0
  !
  ! Set the starting point and the direction
  !
  ray_cart_x    = x
  ray_cart_y    = y
  ray_cart_z    = z
  ray_cart_dirx = dx
  ray_cart_diry = dy
  ray_cart_dirz = dz
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
  ! Put the distance from the starting point to the observer here
  ! If observer at infinity, simply put distance=1d99
  !
  ray_dsend     = distance 
  !
  ! Check if we start inside or outside the domain
  !
  if(igrid_type.lt.100) then 
     !
     ! We use a regular or AMR grid
     !
     if(igrid_coord.lt.100) then
        !
        ! We use cartesian coordinates
        !
        if(amr_tree_present) then
           !
           ! We use the AMR tree
           !
           call amr_findcell(x,y,z,a)
           if(associated(a)) then
              !
              ! We start inside the domain
              !
              ray_index     = a%leafindex
              ray_indexnext = 0
           else
              !
              ! We start outside the domain, so all indices must be 0
              !
              ray_index     = 0
              ray_indexnext = 0
           endif
        else
           !
           ! We have a regular grid
           !
           call amr_findbasecell(x,y,z,ixx,iyy,izz)
           if(ixx.gt.0) then
              ray_index     = ixx+((iyy-1)+(izz-1)*amr_grid_ny)*amr_grid_nx
              ray_indexnext = 0
           else
              ray_index     = 0
              ray_indexnext = 0
           endif
        endif
        !
     elseif(igrid_coord.lt.200) then
        !
        ! We use spherical coordinates
        !
        call amr_xyz_to_rthphi(ray_cart_x,ray_cart_y,ray_cart_z,r,theta,phi)
        if(amr_tree_present) then
           !
           ! We use the AMR tree
           !
           call amr_findcell(r,theta,phi,a)
           if(associated(a)) then
              !
              ! We start inside the domain
              !
              ray_index     = a%leafindex
              ray_indexnext = 0
           else
              !
              ! We start outside the domain, so all indices must be 0
              !
              ray_index     = 0
              ray_indexnext = 0
           endif
        else
           !
           ! We have a regular grid
           !
           call amr_findbasecell(r,theta,phi,ixx,iyy,izz)
           if(ixx.gt.0) then
              ray_index     = ixx+((iyy-1)+(izz-1)*amr_grid_ny)*amr_grid_nx
              ray_indexnext = 0
           else
              ray_index     = 0
              ray_indexnext = 0
           endif
        endif
!        write(stdo,*) 'FINDING CELL NOT YET READY FOR SPHERICAL COORDINATES'
!        stop
        !
     else
        write(stdo,*) 'ERROR: Cylindrical coordinates not yet implemented'
        stop
     endif
  else
     write(stdo,*) 'SORRY: Delaunay or Voronoi grids not yet implemented'
     stop
  endif
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
           ! Check if this cell is the smallest so far
           !
           deepestlevel = max(deepestlevel,levelnext)
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
                    do inu=1,camera_nrfreq
                       intensity(inu,1) = bplanck(thermal_bc_temp(bc_ilr,bc_idir), &
                                                camera_frequencies(inu))
                    enddo
                    intensity(:,2) = 0.d0
                    intensity(:,3) = 0.d0
                    intensity(:,4) = 0.d0
                 endif
              endif
           endif
           !
        elseif(igrid_coord.lt.200) then
           !
           ! We use spherical coordinates
           !
           if((camera_maxdphi.le.0.d0).or.((.not.camera_secondorder).and. &
                (.not.dust_2daniso))) then
              !
              ! If camera_maxdphi is not set to a positive value, then
              ! do not cut segments (always go from one cell wall to the
              ! next in one segment). 
              !
              ! Same is true if first order integration is used, and the 2-D
              ! anisotropic scattering mode is not used. 
              !
              call amrray_find_next_location_spher(ray_dsend,           &
                   ray_cart_x,ray_cart_y,ray_cart_z,                    &
                   ray_cart_dirx,ray_cart_diry,ray_cart_dirz,           &
                   ray_index,ray_indexnext,ray_ds,arrived)
           else
              !
              ! The case when we want to prevent too strong jumps
              ! in the angle of the ray wrt the origin
              !
              call amrray_find_next_location_spher(ray_dsend,           &
                   ray_cart_x,ray_cart_y,ray_cart_z,                    &
                   ray_cart_dirx,ray_cart_diry,ray_cart_dirz,           &
                   ray_index,ray_indexnext,ray_ds,arrived,              &
                   maxdeltasina=camera_maxdphi)
           endif
           !
           ! Check if this cell is the smallest so far. However, this is far
           ! more subtle than for the Cartesian case where all cells are
           ! forced to be precisely cubic. Here, for spherical coordinates,
           ! we can have very flattened cells. For instance, near the inner
           ! rim of the coordinate system we may have very much refinement
           ! in R (apart from AMR!), so the cells are nearly like sheats
           ! perpendicular to radially outward rays. Or near the midplane
           ! we could have strong refinement for some other reason.
           ! So we must take projection into account.
           !
           ! Note: We approximate sin(theta)=theta for theta<pi/2
           !       and pi-theta for theta>pi/2
           !
           if(ray_indexnext.gt.0) then
              if(amr_tree_present) then
                 !
                 ! We get our information about the cell dimensions 
                 ! from the AMR tree
                 !
                 b => amr_index_to_leaf(ray_indexnext)%link
                 if(.not.associated(b)) stop 7130
                 !
                 ! Get the cell walls and centers
                 !
                 do iddr=1,3
                    bxc(iddr)   = amr_finegrid_xc(b%ixyzf(iddr),iddr,b%level)
                    bxi(1,iddr) = amr_finegrid_xi(b%ixyzf(iddr),iddr,b%level)
                    bxi(2,iddr) = amr_finegrid_xi(b%ixyzf(iddr)+1,iddr,b%level)
                 enddo
              else
                 !
                 ! We get our information about the cell dimensions
                 ! from the regular grid 
                 !
                 call amr_regular_get_ixyz(ray_indexnext,amrray_ix_next,amrray_iy_next,amrray_iz_next)
                 if(amrray_ix_next.le.0) stop 7130
                 bxc(1)   = amr_finegrid_xc(amrray_ix_next,1,0)
                 bxc(2)   = amr_finegrid_xc(amrray_iy_next,2,0)
                 bxc(3)   = amr_finegrid_xc(amrray_iz_next,3,0)
                 bxi(1,1) = amr_finegrid_xi(amrray_ix_next,1,0)
                 bxi(2,1) = amr_finegrid_xi(amrray_ix_next+1,1,0)
                 bxi(1,2) = amr_finegrid_xi(amrray_iy_next,2,0)
                 bxi(2,2) = amr_finegrid_xi(amrray_iy_next+1,2,0)
                 bxi(1,3) = amr_finegrid_xi(amrray_iz_next,3,0)
                 bxi(2,3) = amr_finegrid_xi(amrray_iz_next+1,3,0)
              endif
              !
              ! The smallest size in R direction
              !
              dumr = bxi(2,1) - bxi(1,1)
              dumr = max(dumr,bxc(1)*camera_min_drr)
              !
              ! The smallest size in Theta direction
              !
              dumt = (bxi(2,2) - bxi(1,2))
              dumt = bxc(1)*max(min(dumt,camera_max_dangle),camera_min_dangle)
              !
              ! The smallest size in Phi direction
              !
              if(bxc(2).lt.1.d0) then
                 dumm = bxc(2)
              elseif(pi-bxc(2).lt.1.d0) then
                 dumm = pi-bxc(2)
              else
                 dumm = 1.d0
              endif
              dump = dumm*abs(bxi(2,2) - bxi(1,2))
              dump = bxc(1)*max(min(dump,camera_max_dangle),camera_min_dangle)
              !
              ! Compute the projections of the direction on the coordinate
              ! vector
              !
              r     = sqrt(ray_cart_x**2+ray_cart_y**2+ray_cart_z**2)
              proj  = abs( ray_cart_dirx*ray_cart_x +  &
                           ray_cart_diry*ray_cart_y +  &
                           ray_cart_dirz*ray_cart_z ) / r
!              proj = max(proj,camera_min_aspectratio)
              dummy = abs(ray_cart_dirz)
              !
              ! Compute estimates (!) of the projection of cell theta-width
              ! and cell phi-width onto image plane. This is just a simple
              ! estimation based on the z-component of the direction vector.
              ! If we view things face-on (ray_cart_dirz=1 or -1) then 
              ! the theta-size will be subject to projection, while
              ! if we view things edge-on (ray_cart_dirz=0) then 
              ! the phi-size will be subject to projection
              !
              celldxmin = min(celldxmin,max(dumr, &
                   min(dumt,dump,(dummy*dumt+(1.d0-dummy)*dump)*proj)))
           else
              !
              ! Check if the ray passes within the central hole of
              ! the coordinate system. 
              !
              dummy = ray_cart_dirx*ray_cart_x +  &
                      ray_cart_diry*ray_cart_y +  &
                      ray_cart_dirz*ray_cart_z 
              dumx  = ray_cart_x - dummy*ray_cart_dirx
              dumy  = ray_cart_y - dummy*ray_cart_diry
              dumz  = ray_cart_z - dummy*ray_cart_dirz
              dumr  = sqrt(dumx**2+dumy**2+dumz**2)
              if(dumr.le.amr_grid_xi(1,1)) then
                 celldxmin = min(celldxmin,                            &
                      camera_spher_cavity_relres * amr_grid_xi(1,1) )
              endif
           endif
           !
        else
           write(stdo,*) 'ERROR: Cylindrical coordinates not yet implemented'
           stop
        endif
        !
        ! Set the camera_istar, which is used for treating finite-size
        ! stars (instead of point source stars)
        !
        camera_istar = amrray_ispherehit
        !   
     else
        write(stdo,*) 'SORRY: Delaunay or Voronoi grids not yet implemented'
        stop
     endif
     !
     ! Path length
     ! 
     ! ######### CHECK: WHY NOT USE ray_ds from the above routines ??? ########
     !
     ray_ds    = sqrt( (ray_cart_x-ray_prev_x)**2 + (ray_cart_y-ray_prev_y)**2 + (ray_cart_z-ray_prev_z)**2 )
     !
     ! Do a simple safet check
     !
     if(ray_index.gt.nrcellsmax) stop 6271
     !
     ! Check if we are/were in a star
     !
     if(camera_istar.le.0) then
        !
        ! Normal situation: we are not inside a star (istar=0), or we will
        ! hit the star but are not inside the star yet (istar<0).
        !
        ! Decide whether we do 1st order or 2nd order integration
        !
        if(camera_secondorder) then
           !
           ! Second order integration. Corner-based (=vertex-based) scheme.
           ! Here the continuum sources have all been pre-calculated on
           ! the cell corners. For the lines: they are also pre-calculated
           ! unless the camera_catch_doppler_line==.true., because in
           ! that case the line stuff may have to be done locally at smaller
           ! steps. 
           !
           ! Switch between integration of radiation field and 
           ! integration of the opacity.
           !
           if((camera_tracemode.eq.1).or.(camera_tracemode.le.-2)) then
              !
              ! Do formal radiative transfer
              !
              inu = inu0
              !
              ! Backup snu and anu at previous crossing
              !
              anu_prev = anu_curr
              if(camera_stokesvector) then
                 snu_prev(1:4) = snu_curr(1:4)
                 jnu_prev(1:4) = jnu_curr(1:4)
              else
                 snu_prev(1) = snu_curr(1)
                 jnu_prev(1) = jnu_curr(1)
              endif
              !
              ! Also backup stuff specifically for lines, if the
              ! doppler catching method is switched on
              !
              if(camera_catch_doppler_line) then
                 ln_anu_prev = ln_anu_curr
                 ln_jnu_prev = ln_jnu_curr
                 sources_local_doppler_prev = sources_local_doppler_curr
                 sources_local_turb_prev = sources_local_turb_curr
                 sources_local_temp_prev = sources_local_temp_curr
                 sources_local_line_nup_prev(:) = sources_local_line_nup_curr(:)
                 sources_local_line_ndown_prev(:) = sources_local_line_ndown_curr(:)
              endif
              !
              ! Find snu and anu at new crossing
              !
              ! NOTE: If sources_interpol_jnu is set, then snu is actually
              !       j_nu (the emissivity). It will be translated below.
              !
              if(igrid_type.lt.100) then 
                 if(igrid_coord.lt.100) then
                    !
                    ! Cartesian
                    !
                    call sources_find_srcalp_interpol(      &
                         ray_cart_x,ray_cart_y,ray_cart_z,  &
                         snu_curr,anu_curr,camera_stokesvector)
                 elseif(igrid_coord.lt.200) then
                    !
                    ! Spherical 
                    !
                    call amr_xyz_to_rthphi(ray_cart_x,ray_cart_y,ray_cart_z, &
                                           r,theta,phi)
                    call sources_find_srcalp_interpol(    &
                         r,theta,phi,snu_curr,anu_curr,camera_stokesvector)
                 else
                    stop 498
                 endif
              else
                 stop 499
              endif
              !
              ! If sources_interpol_jnu is set, then what we get back
              ! from the sources_find_srcalp_interpol() is in fact the
              ! emissivity j_nu instead of the source function S_nu.
              ! So translate either snu to jnu or vice versa.
              !
              if(sources_interpol_jnu) then
                 jnu_curr = snu_curr
                 snu_curr = jnu_curr / (anu_curr+1d-99)
              else
                 jnu_curr = snu_curr * (anu_curr+1d-99)
              endif
              !
              if(ray_index.gt.0) then
                 !
                 ! Decide if we do "normal" second order integration, or
                 ! if we take special care of Doppler shifts, i.e. catching
                 ! doppler jumps.
                 !
                 if(camera_catch_doppler_line) then
                    !
                    ! Special kind of integration: With doppler-catching.
                    ! This means that we are going to check at each integration
                    ! segment along the ray if the Doppler shift between the
                    ! start and finish of this ray segment is larger than 
                    ! some fraction of the local thermal+microturbulent line
                    ! width. If so, we will do substeps along this ray segment.
                    !
                    nu = camera_frequencies(inu)
                    !
                    ! Call the line module to make the final computation
                    ! from source function and nup and ndown to j_nu and 
                    ! alpha_nu for the endpoint of the current ray segment.
                    ! This is for the line contributions only. The continuum
                    ! stuff is already done before.
                    !
                    call lines_jnu_anu_from_nup_ndown(sources_vertex_lines_nractivetot, &
                         sources_local_line_nup_curr(:),&
                         sources_local_line_ndown_curr(:),&
                         sources_local_doppler_curr,&
                         sources_local_turb_curr,&
                         sources_local_temp_curr,&
                         camera_frequencies(inu),ln_anu_curr,ln_jnu_curr)
                    !
                    ! Now check if there is any line that Doppler-jumps
                    ! over the current wavelength-of-sight...
                    !
                    ! Do a loop over all active lines and do that check.
                    !
                    ! **** NOTE: This can be done more efficiently! ****
                    !
                    nsteps = 0
                    do ispec=1,lines_nr_species
                       !
                       ! Determine the local line width for this molecule,
                       ! times some safety margin. To be on the safe side
                       ! we take the smallest possible value.
                       !
                       temp  = min(sources_local_temp_curr, &
                                   sources_local_temp_prev )
                       duc   = ( min(sources_local_turb_curr, &
                                 sources_local_turb_prev ) + &
                                 sqrt(2*kk*temp/lines_umass(ispec)) ) / cc
                       margin = duc * lines_widthmargin
                       resol  = duc * camera_catch_doppler_resolution
                       !
                       ! Now loop over all active lines of this molecule
                       !
                       do ilactive=1,active_nrlines(ispec)
                          !
                          ! Which lines is this?
                          !
                          iline = active_lines(ilactive,ispec)
                          !
                          ! What is the rest-frame frequency?
                          !
                          nu0   = lines_nu0(iline,ispec)
                          !
                          ! What is the doppler-shifted frequency, i.e. the
                          ! frequency of line center in the lab frame (prev
                          ! and curr)? Note: If we use double-precision, we
                          ! are fine here.
                          !
                          nu0_prev = nu0 * ( 1.d0 + sources_local_doppler_prev )
                          nu0_curr = nu0 * ( 1.d0 + sources_local_doppler_curr )
                          !
                          ! Does a jump happen (distinguish two scenarios)?
                          !
                          if(nu0_curr.gt.nu0_prev) then
                             if((nu0_curr-nu0_prev.gt.resol*nu).and.  &
                                (nu0_curr.gt.nu*(1.d0-margin)).and.   &
                                (nu0_prev.lt.nu*(1.d0+margin))) then
                                !
                                ! Yes, we have a potentially dangerous jump!
                                ! Compute the required refinement, where 
                                ! for simplicity we simply refine the entire
                                ! ray segment (i.e. a slightly simpler
                                ! procedure than in Pontoppidan et al.)
                                !
                                nstp = (nu0_curr-nu0_prev)/(resol*nu+1.d-99)+1.d0
                                if(nstp.gt.10000) then
                                   write(stdo,*) 'ERROR during Doppler-catching of lines:'
                                   write(stdo,*) '      Local line width is so small that '
                                   write(stdo,*) '      I require more than 10000 refinement'
                                   write(stdo,*) '      steps in a single cell! This is'
                                   write(stdo,*) '      suspicious, so I stop.'
                                   write(stdo,*) '      nu       = ',nu
                                   write(stdo,*) '      nu0_line = ',nu0
                                   write(stdo,*) '      du/c     = ',duc
                                   write(stdo,*) '      nstp     = ',nstp
                                   stop
                                endif
                                if(nstp.gt.nsteps) nsteps=nstp
                             endif
                          else
                             if((nu0_prev-nu0_curr.gt.resol*nu).and.  &
                                (nu0_prev.gt.nu*(1.d0-margin)).and.   &
                                (nu0_curr.lt.nu*(1.d0+margin))) then
                                !
                                ! Yes, we have a potentially dangerous jump!
                                ! Compute the required refinement, where 
                                ! for simplicity we simply refine the entire
                                ! ray segment (i.e. a slightly simpler
                                ! procedure than in Pontoppidan et al.)
                                !
                                nstp = (nu0_prev-nu0_curr)/(resol*nu+1.d-99)+1.d0
                                if(nstp.gt.10000) then
                                   write(stdo,*) 'ERROR during Doppler-catching of lines:'
                                   write(stdo,*) '      Local line width is so small that '
                                   write(stdo,*) '      I require more than 10000 refinement'
                                   write(stdo,*) '      steps in a single cell! This is'
                                   write(stdo,*) '      suspicious, so I stop.'
                                   write(stdo,*) '      nu       = ',nu
                                   write(stdo,*) '      nu0_line = ',nu0
                                   write(stdo,*) '      du/c     = ',duc
                                   write(stdo,*) '      nstp     = ',nstp
                                   stop
                                endif
                                if(nstp.gt.nsteps) nsteps=nstp
                             endif
                          endif
                       enddo
                    enddo
                    !
                    ! Check if we have to do a substepping or
                    ! if we can integrate this ray segment in 
                    ! one step.
                    !
                    if(nsteps.eq.0) then
                       !
                       ! No substeps required (no doppler jump occurring).
                       ! Normal second order integration step of dust+lines
                       ! over the entire ray segment.
                       !
                       ! Create the total alpha_nu and j_nu
                       !
                       ap    = anu_prev + ln_anu_prev
                       ac    = anu_curr + ln_anu_curr
                       jp(1) = jnu_prev(1) + ln_jnu_prev
                       jc(1) = jnu_curr(1) + ln_jnu_curr
                       if(camera_stokesvector) then
                          jp(2:4) = jnu_prev(2:4)
                          jc(2:4) = jnu_curr(2:4)
                          sp(1:4) = jp(1:4)/(ap+1.d-99)
                          sc(1:4) = jc(1:4)/(ac+1.d-99)
                       else
                          sp(1) = jp(1)/(ap+1.d-99)
                          sc(1) = jc(1)/(ac+1.d-99)
                       endif
                       !
                       if(camera_tracemode.le.-2) then
                          !
                          ! Integration of the opacity and/or finding of the
                          ! tau=1 surface
                          !
                          dtau1     = 0.5d0 * ( ap + ac ) * ray_ds
                          if(camera_tracemode.eq.-2) then
                             !
                             ! Integration of the optacity to obtain the total optical depth
                             ! along the ray
                             !
                             intensity(inu,1) = intensity(inu,1) + dtau1
                          else
                             !
                             ! Find the tau=camera_taustop point in order to create 
                             ! a "tau=1" surface (or better: "tau=camera_tausurface").
                             !
                             dummy2 = intensity(inu,1) + dtau1
                             if(dummy2.ge.camera_taustop(inu)) then
                                !
                                ! We arrived at the taustop point
                                !
                                ! Self-consistency check
                                !
                                if(camera_nrrefine.gt.0) stop 9245
                                !
                                ! If we arrived here because camera_taustop(inu).le.0
                                ! then do a special treatment, otherwise continue as
                                ! normal
                                !
                                if(camera_taustop(inu).le.0.d0) then
                                   camera_xstop(inu) = -1d91
                                   camera_ystop(inu) = -1d91
                                   camera_zstop(inu) = -1d91
                                   camera_dstop(inu) = -1d91
                                else
                                   !
                                   ! Do a linear interpolation between tau and tau+dtau1
                                   !
                                   epstau = (camera_taustop(inu)-intensity(inu,1))/dtau1
                                   if((epstau.gt.1.d0).or.(epstau.lt.0.d0)) then
                                      write(stdo,*) epstau,dtau1,camera_taustop(inu)
                                      stop 9223
                                   endif
                                   !
                                   ! Find the 3-D position of the point where the linear
                                   ! interpolation put the tau=1 point
                                   !
                                   camera_xstop(inu) = ray_prev_x + epstau * (ray_cart_x - ray_prev_x)
                                   camera_ystop(inu) = ray_prev_y + epstau * (ray_cart_y - ray_prev_y)
                                   camera_zstop(inu) = ray_prev_z + epstau * (ray_cart_z - ray_prev_z)
                                   !
                                   ! Now perform the usual rotation into the image plane, so that
                                   ! the "image" is now this distace to the tau=1 plane.
                                   !
                                   if(.not.camera_localobserver) then
                                      !
                                      ! Observer at infinity: the "distance" is measured toward the 
                                      ! plane going through (0,0,0) perpendicular to the line of sight, 
                                      ! with positive value being toward the observer (more positive = 
                                      ! closer to the observer).
                                      !
                                      ! A rotation in the x,y plane, rotating horizontally around
                                      ! the pointing position.
                                      !
                                      xbk = camera_xstop(inu) - camera_pointing_position(1)
                                      ybk = camera_ystop(inu) - camera_pointing_position(2)
                                      !spx = camera_observer_cos_phi * xbk - camera_observer_sin_phi * ybk
                                      spy = camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
                                      !
                                      ! Then a rotation in the y,z plane
                                      !
                                      ybk = spy
                                      zbk = camera_zstop(inu) - camera_pointing_position(3)
                                      !spy = camera_observer_cos_theta * ybk + camera_observer_sin_theta * zbk
                                      spz =-camera_observer_sin_theta * ybk + camera_observer_cos_theta * zbk
                                      !
                                      ! Take the z-value and put this into camera_dstop(inu)
                                      !
                                      camera_dstop(inu) = spz
                                   else 
                                      !
                                      ! Local observer: the "distance" is the true distance to the 
                                      ! observer (more positive = farther away from the observer).
                                      ! It can be useful for analysis purposes to find the tau=1 surface 
                                      ! to some specific point in the model
                                      !
                                      camera_dstop(inu) =                                            &
                                           sqrt((camera_observer_position(1)-camera_xstop(inu))**2 + &
                                           (camera_observer_position(2)-camera_ystop(inu))**2 +      &
                                           (camera_observer_position(3)-camera_zstop(inu))**2)
                                   endif
                                endif
                                !
                                ! Reset the taustop to infinity
                                !
                                camera_taustop(inu) = 1d91
                                !
                             endif
                             !
                             ! Continue integrating the optical depth
                             !
                             intensity(inu,1) = dummy2
                          endif
                       else
                          !
                          ! Integration of the radiative transfer equation.
                          ! Second order integration (i.e. integrating
                          ! a first order approximation of the opacity and source
                          ! functions).
                          !
                          dtau1     = 0.5d0 * ( ap + ac ) * ray_ds
                          if(camera_stokesvector) then
                             theomax(1:4)   = 0.5d0 * ( jp(1:4) + jc(1:4) ) * ray_ds
                          else
                             theomax(1)   = 0.5d0 * ( jp(1) + jc(1) ) * ray_ds
                          endif
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
                          if(camera_stokesvector) then
                             if(dtau1.gt.1e-9) then
                                qdr(1:4) = ax * sp(1:4) + bx * sc(1:4)
                             else
                                qdr(1:4) = theomax(1:4)
                             endif
                             if(qdr(1).gt.theomax(1)) then
                                thscale = theomax(1)/qdr(1)
                                qdr(1:4) = qdr(1:4) * thscale
                             endif
                             qdr(2) = min(qdr(2),theomax(2))
                             qdr(3) = min(qdr(3),theomax(3))
                             qdr(4) = min(qdr(4),theomax(4))
                             intensity(inu,1:4) = xp * intensity(inu,1:4) + qdr(1:4)
                          else
                             if(dtau1.gt.1e-9) then
                                qdr(1) = ax * sp(1) + bx * sc(1)
                             else
                                qdr(1) = theomax(1)
                             endif
                             qdr(1) = min(qdr(1),theomax(1))
                             intensity(inu,1) = xp * intensity(inu,1) + qdr(1)
                          endif
                       endif
                    else
                       !
                       ! Yes, we have a doppler jump here: A line shifts
                       ! across our current wavelength by more than
                       ! camera_catch_doppler_resolution times the local
                       ! (thermal+turbulent) line width.
                       !
                       ! So we have to do sub-stepping to resolve this
                       ! line doppler shift.
                       !
                       ap    = anu_prev    + ln_anu_prev
                       jp(1) = jnu_prev(1) + ln_jnu_prev
                       if(camera_stokesvector) then
                          jp(2:4) = jnu_prev(2:4)
                          sp(1:4) = jp(1:4)/(ap+1.d-99)
                       else
                          sp(1)   = jp(1)/(ap+1.d-99)
                       endif
                       !
                       ! Memorize the "current" quantities at the end of
                       ! the present ray segment
                       !
                       sources_local_line_nup_end(:)     = sources_local_line_nup_curr(:)
                       sources_local_line_ndown_end(:)   = sources_local_line_ndown_curr(:)
                       sources_local_doppler_end         = sources_local_doppler_curr
                       sources_local_turb_end            = sources_local_turb_curr
                       sources_local_temp_end            = sources_local_temp_curr
                       !
                       ! Do loop over substepping
                       !
                       do istep=1,nsteps
                          !
                          ! Fraction of the ray segment that we are now
                          !
                          eps  = (istep*1.d0)/nsteps
                          eps1 = 1.d0-eps
                          !
                          ! Find new "current" quantities, by linear interpolation
                          ! between begin and end of this segment ("prev" and "end")
                          !
                          sources_local_line_nup_curr(:)     = eps1 * sources_local_line_nup_prev(:)  + &
                                                               eps * sources_local_line_nup_end(:)
                          sources_local_line_ndown_curr(:)   = eps1 * sources_local_line_ndown_prev(:) + &
                                                               eps * sources_local_line_ndown_end(:)
                          sources_local_doppler_curr         = eps1 * sources_local_doppler_prev         + &
                                                               eps * sources_local_doppler_end
                          sources_local_turb_curr            = eps1 * sources_local_turb_prev            + &
                                                               eps * sources_local_turb_end
                          sources_local_temp_curr            = eps1 * sources_local_temp_prev            + &
                                                               eps * sources_local_temp_end
                          !
                          ! Now calculate the alpha_nu and j_nu again
                          !
                          call lines_jnu_anu_from_nup_ndown(sources_vertex_lines_nractivetot, &
                               sources_local_line_nup_curr(:),&
                               sources_local_line_ndown_curr(:),&
                               sources_local_doppler_curr,&
                               sources_local_turb_curr,&
                               sources_local_temp_curr,&
                               camera_frequencies(inu),ln_anu_curr,ln_jnu_curr)
                          ac    = eps1*anu_prev    + eps*anu_curr    + ln_anu_curr
                          jc(1) = eps1*jnu_prev(1) + eps*jnu_curr(1) + ln_jnu_curr    ! BUGFIX 02.09.2012
                          if(camera_stokesvector) then
                             jc(2:4) = eps1*jnu_prev(2:4) + eps*jnu_curr(2:4)
                             sc(1:4) = jc(1:4)/(ac+1.d-99)
                          else
                             sc(1) = jc(1)/(ac+1.d-99)
                          endif
                          !
                          ! Now the integration step
                          !
                          if(camera_tracemode.le.-2) then
                             !
                             ! Integration of the optacity and/or finding of the
                             ! tau=1 surface
                             !
                             dtau1 = 0.5d0 * ( ap + ac ) * ray_ds / nsteps
                             if(camera_tracemode.eq.-2) then
                                !
                                ! Integration of the optacity to obtain the optical depth
                                ! (useful for debugging etc)
                                !
                                intensity(inu,1) = intensity(inu,1) + dtau1
                             else
                                !
                                ! Find the tau=camera_taustop point in order to create 
                                ! a "tau=1" surface (or better: "tau=camera_tausurface").
                                !
                                dummy2 = intensity(inu,1) + dtau1
                                if(dummy2.ge.camera_taustop(inu)) then
                                   !
                                   ! We arrived at the taustop point
                                   !
                                   ! Self-consistency check
                                   !
                                   if(camera_nrrefine.gt.0) stop 9245
                                   !
                                   ! If we arrived here because camera_taustop(inu).le.0
                                   ! then do a special treatment, otherwise continue as
                                   ! normal
                                   !
                                   if(camera_taustop(inu).le.0.d0) then
                                      camera_xstop(inu) = -1d91
                                      camera_ystop(inu) = -1d91
                                      camera_zstop(inu) = -1d91
                                      camera_dstop(inu) = -1d91
                                   else
                                      !
                                      ! Do a linear interpolation between tau and tau+dtau1
                                      !
                                      epstau = (camera_taustop(inu)-intensity(inu,1))/dtau1
                                      if((epstau.gt.1.d0).or.(epstau.lt.0.d0)) then
                                         write(stdo,*) epstau,dtau1,camera_taustop(inu)
                                         stop 9223
                                      endif
                                      !
                                      ! Find the 3-D position of the point where the linear
                                      ! interpolation put the tau=1 point. Note that we must
                                      ! here take into account the substepping (doppler catching)
                                      !
                                      eps0 = (istep*1.d0-1.d0)/nsteps
                                      rx0  = ray_prev_x + eps0 * (ray_cart_x - ray_prev_x)
                                      ry0  = ray_prev_y + eps0 * (ray_cart_y - ray_prev_y)
                                      rz0  = ray_prev_z + eps0 * (ray_cart_z - ray_prev_z)
                                      rx1  = ray_prev_x + eps * (ray_cart_x - ray_prev_x)
                                      ry1  = ray_prev_y + eps * (ray_cart_y - ray_prev_y)
                                      rz1  = ray_prev_z + eps * (ray_cart_z - ray_prev_z)
                                      camera_xstop(inu) = rx0 + epstau * (rx1 - rx0)
                                      camera_ystop(inu) = ry0 + epstau * (ry1 - ry0)
                                      camera_zstop(inu) = rz0 + epstau * (rz1 - rz0)
                                      !
                                      ! Now perform the usual rotation into the image plane, so that
                                      ! the "image" is now this distace to the tau=1 plane.
                                      !
                                      if(.not.camera_localobserver) then
                                         !
                                         ! Observer at infinity: the "distance" is measured toward the 
                                         ! plane going through (0,0,0) perpendicular to the line of sight, 
                                         ! with positive value being toward the observer (more positive = 
                                         ! closer to the observer).
                                         !
                                         ! A rotation in the x,y plane, rotating horizontally around
                                         ! the pointing position.
                                         !
                                         xbk = camera_xstop(inu) - camera_pointing_position(1)
                                         ybk = camera_ystop(inu) - camera_pointing_position(2)
                                         !spx = camera_observer_cos_phi * xbk - camera_observer_sin_phi * ybk
                                         spy = camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
                                         !
                                         ! Then a rotation in the y,z plane
                                         !
                                         ybk = spy
                                         zbk = camera_zstop(inu) - camera_pointing_position(3)
                                         !spy = camera_observer_cos_theta * ybk + camera_observer_sin_theta * zbk
                                         spz =-camera_observer_sin_theta * ybk + camera_observer_cos_theta * zbk
                                         !
                                         ! Take the z-value and put this into "intensity(inu,1)"
                                         !
                                         camera_dstop(inu) = spz
                                      else 
                                         !
                                         ! Local observer: the "distance" is the true distance to the 
                                         ! observer (more positive = farther away from the observer).
                                         ! It can be useful for analysis purposes to find the tau=1 surface 
                                         ! to some specific point in the model
                                         !
                                         camera_dstop(inu) =                                            &
                                              sqrt((camera_observer_position(1)-camera_xstop(inu))**2 + &
                                              (camera_observer_position(2)-camera_ystop(inu))**2 +      &
                                              (camera_observer_position(3)-camera_zstop(inu))**2)
                                      endif
                                      !
                                      ! Reset the taustop to infinity
                                      !
                                      camera_taustop(inu) = 1d91
                                      !
                                   endif
                                endif
                                !
                                ! Continue integrating the optical depth
                                !
                                intensity(inu,1) = dummy2
                             endif
                          else
                             !
                             ! Integration of the radiative transfer equation.
                             ! Second order integration (i.e. integrating
                             ! a first order approximation of the opacity and source
                             ! functions).
                             !
                             dtau1     = 0.5d0 * ( ap + ac ) * ray_ds / nsteps
                             if(camera_stokesvector) then
                                theomax(1:4) = 0.5d0 * ( jp(1:4) + jc(1:4) ) * ray_ds / nsteps
                             else
                                theomax(1)   = 0.5d0 * ( jp(1) + jc(1) ) * ray_ds / nsteps
                             endif
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
                             if(camera_stokesvector) then
                                if(dtau1.gt.1e-9) then
                                   qdr(1:4) = ax * sp(1:4) + bx * sc(1:4)
                                else
                                   qdr(1:4) = theomax(1:4)
                                endif
                                if(qdr(1).gt.theomax(1)) then
                                   thscale = theomax(1)/qdr(1)
                                   qdr(1:4) = qdr(1:4) * thscale
                                endif
                                qdr(2) = min(qdr(2),theomax(2))
                                qdr(3) = min(qdr(3),theomax(3))
                                qdr(4) = min(qdr(4),theomax(4))
                                intensity(inu,1:4) = xp * intensity(inu,1:4) + qdr(1:4)
                             else
                                if(dtau1.gt.1e-9) then
                                   qdr(1) = ax * sp(1) + bx * sc(1)
                                else
                                   qdr(1) = theomax(1)
                                endif
                                qdr(1) = min(qdr(1),theomax(1))
                                intensity(inu,1) = xp * intensity(inu,1) + qdr(1)
                             endif
                          endif
                          !
                          ! Now shift curr to prev
                          !
                          ap = ac
                          if(camera_stokesvector) then
                             jp(1:4) = jc(1:4)
                             sp(1:4) = sc(1:4)
                          else
                             jp(1) = jc(1)
                             sp(1) = sc(1)
                          endif
                          !
                          ! Now next substep
                          !
                       enddo
                    endif
                 else
                    !
                    ! Normal kind of integration (no sub-stepping and thus
                    ! no doppler catching). Here the anu and jnu already 
                    ! contain the line contributions.
                    !
                    if(camera_tracemode.le.-2) then
                       !
                       ! Integration of the optacity and/or finding of the
                       ! tau=1 surface
                       !
                       dtau1 = 0.5d0 * ( anu_prev + anu_curr ) * ray_ds
                       if(camera_tracemode.eq.-2) then
                          !
                          ! Integration of the optacity to obtain the optical depth
                          ! (useful for debugging etc)
                          !
                          intensity(inu,1) = intensity(inu,1) + dtau1
                       else
                          !
                          ! Find the tau=camera_taustop point in order to create 
                          ! a "tau=1" surface (or better: "tau=camera_tausurface").
                          !
                          dummy2 = intensity(inu,1) + dtau1
                          if(dummy2.ge.camera_taustop(inu)) then
                             !
                             ! We arrived at the taustop point
                             !
                             ! Self-consistency check
                             !
                             if(camera_nrrefine.gt.0) stop 9245
                             !
                             ! If we arrived here because camera_taustop(inu).le.0
                             ! then do a special treatment, otherwise continue as
                             ! normal
                             !
                             if(camera_taustop(inu).le.0.d0) then
                                camera_xstop(inu) = -1d91
                                camera_ystop(inu) = -1d91
                                camera_zstop(inu) = -1d91
                                camera_dstop(inu) = -1d91
                             else
                                !
                                ! Do a linear interpolation between tau and tau+dtau1
                                !
                                epstau = (camera_taustop(inu)-intensity(inu,1))/dtau1
                                if((epstau.gt.1.d0).or.(epstau.lt.0.d0)) then
                                   write(stdo,*) epstau,dtau1,camera_taustop(inu)
                                   stop 9223
                                endif
                                !
                                ! Find the 3-D position of the point where the linear
                                ! interpolation put the tau=1 point
                                !
                                camera_xstop(inu) = ray_prev_x + epstau * (ray_cart_x - ray_prev_x)
                                camera_ystop(inu) = ray_prev_y + epstau * (ray_cart_y - ray_prev_y)
                                camera_zstop(inu) = ray_prev_z + epstau * (ray_cart_z - ray_prev_z)
                                !
                                ! Now perform the usual rotation into the image plane, so that
                                ! the "image" is now this distace to the tau=1 plane.
                                !
                                if(.not.camera_localobserver) then
                                   !
                                   ! Observer at infinity: the "distance" is measured toward the 
                                   ! plane going through (0,0,0) perpendicular to the line of sight, 
                                   ! with positive value being toward the observer (more positive = 
                                   ! closer to the observer).
                                   !
                                   ! A rotation in the x,y plane, rotating horizontally around
                                   ! the pointing position.
                                   !
                                   xbk = camera_xstop(inu) - camera_pointing_position(1)
                                   ybk = camera_ystop(inu) - camera_pointing_position(2)
                                   !spx = camera_observer_cos_phi * xbk - camera_observer_sin_phi * ybk
                                   spy = camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
                                   !
                                   ! Then a rotation in the y,z plane
                                   !
                                   ybk = spy
                                   zbk = camera_zstop(inu) - camera_pointing_position(3)
                                   !spy = camera_observer_cos_theta * ybk + camera_observer_sin_theta * zbk
                                   spz =-camera_observer_sin_theta * ybk + camera_observer_cos_theta * zbk
                                   !
                                   ! Take the z-value and put this into "intensity(inu,1)"
                                   !
                                   camera_dstop(inu) = spz
                                else 
                                   !
                                   ! Local observer: the "distance" is the true distance to the 
                                   ! observer (more positive = farther away from the observer).
                                   ! It can be useful for analysis purposes to find the tau=1 surface 
                                   ! to some specific point in the model
                                   !
                                   camera_dstop(inu) =                                            &
                                        sqrt((camera_observer_position(1)-camera_xstop(inu))**2 + &
                                        (camera_observer_position(2)-camera_ystop(inu))**2 +      &
                                        (camera_observer_position(3)-camera_zstop(inu))**2)
                                endif
                                !
                                ! Reset the taustop to infinity
                                !
                                camera_taustop(inu) = 1d91
                                !
                             endif
                          endif
                          !
                          ! Continue integrating the optical depth
                          !
                          intensity(inu,1) = dummy2
                       endif
                    else
                       !
                       ! Integration of the radiative transfer equation.
                       ! Second order integration (i.e. integrating
                       ! a first order approximation of the opacity and source
                       ! functions).
                       !
                       dtau1     = 0.5d0 * ( anu_prev + anu_curr ) * ray_ds
                       if(camera_stokesvector) then
                          theomax(1:4) = 0.5d0 * ( jnu_prev(1:4) + jnu_curr(1:4) ) * ray_ds
                       else
                          theomax(1)   = 0.5d0 * ( jnu_prev(1) + jnu_curr(1) ) * ray_ds
                       endif
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
                       if(camera_stokesvector) then
                          if(dtau1.gt.1e-9) then
                             qdr(1:4) = ax * snu_prev(1:4) + bx * snu_curr(1:4)
                          else
                             qdr(1:4) = theomax(1:4)
                          endif
                          if(qdr(1).gt.theomax(1)) then
                             thscale = theomax(1)/qdr(1)
                             qdr(1:4) = qdr(1:4) * thscale
                          endif
                          qdr(2) = sign(min(abs(qdr(2)),abs(theomax(2))),theomax(2))
                          qdr(3) = sign(min(abs(qdr(3)),abs(theomax(3))),theomax(3))
                          qdr(4) = sign(min(abs(qdr(4)),abs(theomax(4))),theomax(4))
                          intensity(inu,1:4) = xp * intensity(inu,1:4) + qdr(1:4)
                       else
                          if(dtau1.gt.1e-9) then
                             qdr(1) = ax * snu_prev(1) + bx * snu_curr(1)
                          else
                             qdr(1) = theomax(1)
                          endif
                          qdr(1) = min(qdr(1),theomax(1))
                          intensity(inu,1) = xp * intensity(inu,1) + qdr(1)
                       endif
                    endif
                 endif
              endif
              !
           elseif(camera_tracemode.eq.-1) then
              !
              ! Simply integrate the dust density
              !
              ! Make 'intensity' the integrated density (column density)
              !
              do ispec=1,dust_nr_species
                 intensity(1,1) = intensity(1,1) + sources_dustdens(ispec) * ray_ds
              enddo
              !
           else
              write(stdo,*) 'ERROR in camera module:'
              write(stdo,*) 'Do not know tracemode = ',camera_tracemode
              stop
           endif
           !
        else
           !
           ! First order integration. Cell-based scheme.
           !
           if(alignment_mode.eq.0) then
              !
              ! Get the src and alp values
              !
              call sources_get_src_alp(inu0,inu1,nrfreq,src,alp,camera_stokesvector)
              !
              ! Now do the integration of the formal transfer equation
              !
              if(camera_tracemode.eq.1) then
                 !
                 ! Do formal radiative transfer
                 !
                 do inu=inu0,inu1
                    !
                    ! Compute the exp(-tau)
                    !
                    xp = alp(inu) * ray_ds
                    if(xp.lt.0.d0) stop 7329
                    if(xp.gt.1d-4) then
                       xp  = exp(-xp)
                       xp1 = 1.d0 - xp
                    else
                       xp1 = xp
                       xp  = exp(-xp)
                    endif
                    !
                    ! Now do the RT, including everything
                    !
                    if(camera_stokesvector) then
                       intensity(inu,1:4) = xp * intensity(inu,1:4) + xp1 * src(inu,1:4)
                    else
                       intensity(inu,1)   = xp * intensity(inu,1)   + xp1 * src(inu,1)
                    endif
                 enddo
                 !
              elseif(camera_tracemode.eq.-1) then
                 !
                 ! Simply integrate the dust density
                 !
                 ! Make 'intensity' the integrated density (column density)
                 !
                 do ispec=1,dust_nr_species
                    intensity(inu0,1) = intensity(inu0,1) + sources_dustdens(ispec) * ray_ds
                 enddo
                 !
              elseif(camera_tracemode.le.-2) then
                 !
                 ! Integration of the optacity and/or finding of the
                 ! tau=1 surface
                 !
                 if(camera_tracemode.eq.-2) then
                    !
                    ! Integration of the optacity to obtain the total optical depth
                    ! along the ray
                    !
                    do inu=inu0,inu1
                       intensity(inu,1) = intensity(inu,1) + alp(inu) * ray_ds
                    enddo
                 else
                    !
                    ! Find the tau=camera_taustop point in order to create 
                    ! a "tau=1" surface (or better: "tau=camera_tausurface").
                    !
                    do inu=inu0,inu1
                       !
                       ! Compute the delta tau at this frequency
                       !
                       dtau1 = alp(inu) * ray_ds
                       !
                       ! Now check if we reached the stopping point yet
                       !
                       dummy2 = intensity(inu,1) + dtau1
                       if(dummy2.ge.camera_taustop(inu)) then
                          !
                          ! We arrived at the taustop point
                          !
                          ! Self-consistency check
                          !
                          if(camera_nrrefine.gt.0) stop 9245
                          !
                          ! If we arrived here because camera_taustop(inu).le.0
                          ! then do a special treatment, otherwise continue as
                          ! normal
                          !
                          if(camera_taustop(inu).le.0.d0) then
                             camera_xstop(inu) = -1d91
                             camera_ystop(inu) = -1d91
                             camera_zstop(inu) = -1d91
                             camera_dstop(inu) = -1d91
                          else
                             !
                             ! Do a linear interpolation between tau and tau+dtau1
                             !
                             epstau = (camera_taustop(inu)-intensity(inu,1))/dtau1
                             if((epstau.gt.1.d0).or.(epstau.lt.0.d0)) then
                                write(stdo,*) epstau,dtau1,camera_taustop(inu)
                                stop 9223
                             endif
                             !
                             ! Find the 3-D position of the point where the linear
                             ! interpolation put the tau=1 point
                             !
                             camera_xstop(inu) = ray_prev_x + epstau * (ray_cart_x - ray_prev_x)
                             camera_ystop(inu) = ray_prev_y + epstau * (ray_cart_y - ray_prev_y)
                             camera_zstop(inu) = ray_prev_z + epstau * (ray_cart_z - ray_prev_z)
                             !
                             ! Now perform the usual rotation into the image plane, so that
                             ! the "image" is now this distace to the tau=1 plane.
                             !
                             if(.not.camera_localobserver) then
                                !
                                ! Observer at infinity: the "distance" is measured toward the 
                                ! plane going through (0,0,0) perpendicular to the line of sight, 
                                ! with positive value being toward the observer (more positive = 
                                ! closer to the observer).
                                !
                                ! A rotation in the x,y plane, rotating horizontally around
                                ! the pointing position.
                                !
                                xbk = camera_xstop(inu) - camera_pointing_position(1)
                                ybk = camera_ystop(inu) - camera_pointing_position(2)
                                !spx = camera_observer_cos_phi * xbk - camera_observer_sin_phi * ybk
                                spy = camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
                                !
                                ! Then a rotation in the y,z plane
                                !
                                ybk = spy
                                zbk = camera_zstop(inu) - camera_pointing_position(3)
                                !spy = camera_observer_cos_theta * ybk + camera_observer_sin_theta * zbk
                                spz =-camera_observer_sin_theta * ybk + camera_observer_cos_theta * zbk
                                !
                                ! Take the z-value and put this into camera_dstop(inu)
                                !
                                camera_dstop(inu) = spz
                             else 
                                !
                                ! Local observer: the "distance" is the true distance to the 
                                ! observer (more positive = farther away from the observer).
                                ! It can be useful for analysis purposes to find the tau=1 surface 
                                ! to some specific point in the model
                                !
                                camera_dstop(inu) =                                            &
                                     sqrt((camera_observer_position(1)-camera_xstop(inu))**2 + &
                                     (camera_observer_position(2)-camera_ystop(inu))**2 +      &
                                     (camera_observer_position(3)-camera_zstop(inu))**2)
                             endif
                             !
                             ! Reset the taustop to infinity
                             !
                             camera_taustop(inu) = 1d91
                             !
                          endif
                       endif
                       !
                       ! Continue integrating the optical depth
                       !
                       intensity(inu,1) = dummy2
                    enddo
                 endif
              else
                 write(stdo,*) 'ERROR in camera module:'
                 write(stdo,*) 'Do not know tracemode = ',camera_tracemode
                 stop
              endif
           else
              !
              ! Aligned grains: thermal polarized emission. First order integration
              ! using a special-purpose routine from the polarization_module.f90
              !
              if(camera_tracemode.ne.1) then
                 write(stdo,*) 'ERROR: When using aligned grains you cannot use ', &
                      'any other camera_tracemode than 1. Sorry...'
                 stop
              endif
              !
              ! Only do RT if iray_index.gt.0 (i.e. if we are inside a cell)
              !
              if(ray_index.gt.0) then
                 !
                 ! Get the src and alp values. Note that sources_get_src_alp() 
                 ! knows itself that thermal dust emission as well as 
                 ! all dust opacities have to be skipped because it checks the
                 ! alignment_mode flag.
                 !
                 call sources_get_src_alp(inu0,inu1,nrfreq,src,alp,camera_stokesvector)
                 !
                 ! Now call the first order integration routine that takes
                 ! care of the alignment effects
                 !
                 call pol_integrate_rt_aligned(intensity,camera_dir,camera_svec,  &
                                     grainalign_dir(:,ray_index),                 &
                                     grainalign_eff(ray_index),inu0,inu1,src,alp, &
                                     dustdens(:,ray_index),dusttemp(:,ray_index), &
                                     ray_ds)
              endif
           endif
        endif
     elseif(camera_incl_stars.ne.0) then
        !
        ! Special situation: we are inside a star.
        ! So set the intensity to that of the star.
        !
        do inu=inu0,inu1
           intensity(inu,1) = find_starlight_interpol(camera_frequencies(inu), &
                                                      camera_istar)
        enddo
        if(camera_stokesvector) then
           do inu=inu0,inu1
              intensity(inu,2:4) = 0.d0
           enddo
        endif
        !
     endif
     !
     ! Decrease ray_dsend
     !
     ray_dsend = ray_dsend - ray_ds
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
  ! Return the smallest cell size we encountered
  !
  if(igrid_type.lt.100) then 
     if(igrid_coord.lt.100) then
        celldxmin = (amr_grid_xi(2,1)-amr_grid_xi(1,1)) / (2**deepestlevel)
     endif
     ! For the other coordinates celldxmin is determined on-the-fly
  else
     write(stdo,*) 'SORRY: Delaunay or Voronoi grids not yet implemented'
     stop
  endif
  !
end subroutine camera_serial_raytrace


!-------------------------------------------------------------------------
!             COMPUTE AVERAGE INTENSITY IN A SINGLE PIXEL, 
!                WITH POSSIBLE RECURSIVE SUB-PIXELING
!-------------------------------------------------------------------------
recursive subroutine camera_compute_one_pixel(nrfreq,inu0,inu1,px,py,pdx,pdy,  &
                                              nrrefine,intensity,istar)
  implicit none
  integer :: nrfreq,istar
  double precision :: intensity(nrfreq,1:4)
  double precision :: intensdum(nrfreq,1:4)
  double precision :: intensdum11(nrfreq,1:4)
  double precision :: intensdum12(nrfreq,1:4)
  double precision :: intensdum13(nrfreq,1:4)
  double precision :: intensdum14(nrfreq,1:4)
  double precision :: intensdum2(nrfreq,1:4)
  double precision :: px,py,pdx,pdy
  double precision :: x1,y1,dx1,dy1
  double precision :: x,y,z,dirx,diry,dirz,distance
  double precision :: celldxmin,dum1,factor,rmin
  integer :: nrrefine,idum,inu0,inu1,inu,ns,istar1,is
  logical :: flag,todo,donerefine
  integer :: id
  !$ integer OMP_get_thread_num
  !
  ! Check
  !
  if(abs(pdx-pdy).gt.1.d-6*(pdx+pdy)) then
     write(stdo,*) 'ERROR in camera module: Non-square pixels not allowed.'
     write(stdo,*) '  pdx = ',pdx,', pdy = ' ,pdy
     stop
  endif
  if(nrfreq.ne.camera_nrfreq) stop 3345
  !
  ! For diagnostics only: write out (if requested) the position of this pixel
  !
  if(camera_diagnostics_subpix) then
     write(10,30) px,py,pdx,pdy
30   format(4(E17.10,1X))
  endif
  !
  ! Find the starting point and direction of this ray
  !
  call camera_set_ray(px,py,x,y,z,dirx,diry,dirz,distance)
  !
  ! Reset intensity
  ! 
  if(incl_extlum.eq.0) then
     intensity(inu0:inu1,1:4) = 0.d0
  else
     do inu=inu0,inu1
        intensity(inu,1)   = find_extlumintens_interpol(camera_frequencies(inu))
        intensity(inu,2:4) = 0.d0
     enddo
  endif
  !
  ! Reset flag
  !
  donerefine = .false.
  !
  ! Call the ray-tracer 
  !
  call camera_serial_raytrace(nrfreq,inu0,inu1,x,y,z,dirx,diry,dirz,distance,celldxmin,intensity)
  !
  ! Increase the counter
  !
  ! Need critical to get the sub-pixeling numbers correct, but slows down the
  ! parallelization a lot. I would suggest not using. - Patrick Sheehan
  !!!!!!$OMP CRITICAL
  camera_subpixeling_npixtot = camera_subpixeling_npixtot + 1
  !!!!!!$OMP END CRITICAL
  !
  ! Check if we need to refine our pixel
  !
  if(max(pdx,pdy).gt.celldxmin*camera_refine_criterion) then
     !
     ! Are we still allowed to refine?
     !
     if(nrrefine.gt.0) then
        !
        ! OK, we drop the above computed intensity and instead create
        ! four recursive new calls to camera_compute_one_pixel
        !
        idum = nrrefine-1
        dx1  = 0.5d0*pdx
        dy1  = 0.5d0*pdy
        intensity(inu0:inu1,1:4) = 0.d0
        !
        ! OPENMP PARALLELIZATION HERE TO SPEED THINGS UP.
        !
        !$OMP TASK PRIVATE(x1,y1) SHARED(intensdum11) &
        !$OMP FIRSTPRIVATE(nrfreq,inu0,inu1,dx1,dy1,idum,istar)
        x1   = px-0.25d0*pdx
        y1   = py-0.25d0*pdy
        call camera_compute_one_pixel(nrfreq,inu0,inu1,x1,y1,dx1,dy1,idum,intensdum11,istar)
        !intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum11(inu0:inu1,1:4)
        !$OMP END TASK
        !
        !$OMP TASK PRIVATE(x1,y1) SHARED(intensdum12) &
        !$OMP FIRSTPRIVATE(nrfreq,inu0,inu1,dx1,dy1,idum,istar)
        x1   = px+0.25d0*pdx
        y1   = py-0.25d0*pdy
        call camera_compute_one_pixel(nrfreq,inu0,inu1,x1,y1,dx1,dy1,idum,intensdum12,istar)
        !intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum12(inu0:inu1,1:4)
        !$OMP END TASK
        !
        !$OMP TASK PRIVATE(x1,y1) SHARED(intensdum13) &
        !$OMP FIRSTPRIVATE(nrfreq,inu0,inu1,dx1,dy1,idum,istar)
        x1   = px-0.25d0*pdx
        y1   = py+0.25d0*pdy
        call camera_compute_one_pixel(nrfreq,inu0,inu1,x1,y1,dx1,dy1,idum,intensdum13,istar)
        !intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum13(inu0:inu1,1:4)
        !$OMP END TASK
        !
        !$OMP TASK PRIVATE(x1,y1) SHARED(intensdum14) &
        !$OMP FIRSTPRIVATE(nrfreq,inu0,inu1,dx1,dy1,idum,istar)
        x1   = px+0.25d0*pdx
        y1   = py+0.25d0*pdy
        call camera_compute_one_pixel(nrfreq,inu0,inu1,x1,y1,dx1,dy1,idum,intensdum14,istar)
        !intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum14(inu0:inu1,1:4)
        !$OMP END TASK
        !
        !$OMP TASKWAIT
        !
        intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum11(inu0:inu1,1:4)
        intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum12(inu0:inu1,1:4)
        intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum13(inu0:inu1,1:4)
        intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum14(inu0:inu1,1:4)
        intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) * 0.25d0
        !
        donerefine = .true.
     else
        !
        ! Need to refine, but not allowed. We must warn that our resolving
        ! power is simply not sufficient.
        !
        camera_warn_resolution = .true.
     endif
  else
     !
     ! No refinement was done, so this pixel also counts as a "fine" pixel
     ! So increase that counter (this is just for diagnostics; it's non-essential)
     !
     !!!!!!$OMP CRITICAL
     camera_subpixeling_npixfine = camera_subpixeling_npixfine + 1
     !!!!!!$OMP END CRITICAL
  endif
  !
  ! Include stellar spheres
  !
  if(star_sphere.and.(.not.donerefine)) then
     !
     ! If the stars are spheres, we have to do a check if we should
     ! continue to refine
     !
     if((camera_incl_stars.eq.1).and.(istar.ne.0)) then
        if(istar.gt.0) then
           !
           ! There is precisely 1 star that contributes to this pixel
           !
           todo = .true.
           !
           ! Check if the current pixel size is already small enough.
           ! If yes, then don't continue
           !
           if(max(pdx,pdy).le.2.d0*star_r(istar)/camera_starsphere_nrpix) then
              todo = .false.
           endif
           !
           ! Check if this star is entirely inside this pixel, which
           ! would make the calculation easier. But only do this if
           ! fluxcons is not on.
           !
           if(todo.and.(nrrefine.lt.0)) then
              !
              ! Check if sphere entirely inside pixel
              !
              flag = camera_check_starsphere_in_pixel(istar,px,py,pdx,pdy,.true.)
              if(flag) then
                 !
                 ! Yes! The star lies entirely inside this pixel, and 
                 ! flux conservation is not strict, so we can now treat
                 ! this star as a point source.
                 !
                 ! Set the starting position and direction of the ray, as well as
                 ! the location in the image
                 !
                 call camera_set_ray_stars_pntsrc(istar,x,y,z,dirx,diry,dirz,px,py,distance)
                 !
                 ! Set the intensity to the stellar spectrum
                 ! 
                 do inu=inu0,inu1
                    intensdum(inu,1) =                                          &
                         find_starlight_interpol(camera_frequencies(inu),istar)
                    intensdum(inu,2:4) = 0.d0
                 enddo
                 !
                 ! Do ray tracing
                 !
                 call camera_serial_raytrace(camera_nrfreq,inu0,inu1,        &
                                             x,y,z,dirx,diry,dirz,distance,  &
                                             dum1,intensdum)
                 !
                 ! Compute the ratio starsurface / pixelsurface
                 !
                 if(camera_localobserver) then
                    stop 7390
                    !                factor = pi*(star_r(istar)/distance)**2 / (pdx*pdy)
                 else
                    factor = pi*star_r(istar)**2 / (pdx*pdy)
                 endif
                 !
                 ! Now modify the intensity of the pixel to include the star
                 !
                 intensity(inu0:inu1,1:4) =                          &
                      (1.d0-factor) * intensity(inu0:inu1,1:4) +     &  
                      factor * intensdum(inu0:inu1,1:4)
                 !
                 ! Signal that we are done
                 !
                 todo = .false.
              endif
           endif
           !
           ! If we are not yet done, we need to find which of the 2x2 
           ! sub-pixels have to be refined
           !
           if(todo) then
              !
              ! Reset some stuff
              !
              factor = 1.d0
              intensdum2(inu0:inu1,1:4) = intensity(inu0:inu1,1:4)
              dx1  = 0.5d0*pdx
              dy1  = 0.5d0*pdy
              intensity(inu0:inu1,1:4) = 0.d0
              !
              ! Pixel 1,1
              !
              x1   = px-0.25d0*pdx
              y1   = py-0.25d0*pdy
              flag = camera_check_starsphere_in_pixel(istar,x1,y1,dx1,dy1,.false.)
              if(flag) then
                 call camera_compute_one_pixel(nrfreq,inu0,inu1,x1,y1,dx1,dy1,idum,intensdum,istar)
                 intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum(inu0:inu1,1:4)
                 factor = factor - 0.25d0
              endif
              !
              ! Pixel 2,1
              !
              x1   = px+0.25d0*pdx
              y1   = py-0.25d0*pdy
              flag = camera_check_starsphere_in_pixel(istar,x1,y1,dx1,dy1,.false.)
              if(flag) then
                 call camera_compute_one_pixel(nrfreq,inu0,inu1,x1,y1,dx1,dy1,idum,intensdum,istar)
                 intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum(inu0:inu1,1:4)
                 factor = factor - 0.25d0
              endif
              !
              ! Pixel 1,2
              !
              x1   = px-0.25d0*pdx
              y1   = py+0.25d0*pdy
              flag = camera_check_starsphere_in_pixel(istar,x1,y1,dx1,dy1,.false.)
              if(flag) then
                 call camera_compute_one_pixel(nrfreq,inu0,inu1,x1,y1,dx1,dy1,idum,intensdum,istar)
                 intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum(inu0:inu1,1:4)
                 factor = factor - 0.25d0
              endif
              !
              ! Pixel 2,2
              !
              x1   = px+0.25d0*pdx
              y1   = py+0.25d0*pdy
              flag = camera_check_starsphere_in_pixel(istar,x1,y1,dx1,dy1,.false.)
              if(flag) then
                 call camera_compute_one_pixel(nrfreq,inu0,inu1,x1,y1,dx1,dy1,idum,intensdum,istar)
                 intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum(inu0:inu1,1:4)
                 factor = factor - 0.25d0
              endif
              !
              intensity(inu0:inu1,1:4) = factor * intensdum2(inu0:inu1,1:4) +    &
                                     intensity(inu0:inu1,1:4) * 0.25d0
              !
           endif
        elseif(istar.le.-1) then
           !
           ! We don't know how many stars. We have to count.
           !
           ns = 0
           istar1 = 0
           rmin = 1d99
           do is=1,nstars
              flag = camera_check_starsphere_in_pixel(is,px,py,pdx,pdy,.false.)
              if(flag) then
                 ns = ns + 1
                 if(star_r(is).lt.rmin) then
                    istar1 = is
                    rmin = star_r(is)
                 endif
              endif
           enddo
           !
           ! Check if we still have to refine even if we have one or multiple
           ! stars
           !
           if(max(pdx,pdy).gt.2.d0*rmin/camera_starsphere_nrpix) then
              !
              ! Yes, we have to refine
              !
              if(ns.eq.1) then
                 !
                 ! If there is one star, then redo this pixel with istar=istar1
                 !
                 idum = nrrefine
                 call camera_compute_one_pixel(nrfreq,inu0,inu1,px,py,pdx,pdy,idum,intensdum,istar1)
                 intensity(inu0:inu1,1:4) = intensdum(inu0:inu1,1:4)
              elseif(ns.gt.1) then
                 !
                 ! More than one star. Refine.
                 !
                 ! Reset some stuff
                 !
                 dx1  = 0.5d0*pdx
                 dy1  = 0.5d0*pdy
                 intensity(inu0:inu1,1:4) = 0.d0
                 !
                 ! Pixel 1,1
                 !
                 x1   = px-0.25d0*pdx
                 y1   = py-0.25d0*pdy
                 istar1 = -1
                 call camera_compute_one_pixel(nrfreq,inu0,inu1,x1,y1,dx1,dy1,idum,intensdum,istar1)
                 intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum(inu0:inu1,1:4)
                 !
                 ! Pixel 2,1
                 !
                 x1   = px+0.25d0*pdx
                 y1   = py-0.25d0*pdy
                 istar1 = -1
                 call camera_compute_one_pixel(nrfreq,inu0,inu1,x1,y1,dx1,dy1,idum,intensdum,istar1)
                 intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum(inu0:inu1,1:4)
                 !
                 ! Pixel 1,2
                 !
                 x1   = px-0.25d0*pdx
                 y1   = py+0.25d0*pdy
                 istar1 = -1
                 call camera_compute_one_pixel(nrfreq,inu0,inu1,x1,y1,dx1,dy1,idum,intensdum,istar1)
                 intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum(inu0:inu1,1:4)
                 !
                 ! Pixel 2,2
                 !
                 x1   = px+0.25d0*pdx
                 y1   = py+0.25d0*pdy
                 istar1 = -1
                 call camera_compute_one_pixel(nrfreq,inu0,inu1,x1,y1,dx1,dy1,idum,intensdum,istar1)
                 intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) + intensdum(inu0:inu1,1:4)
                 !
                 ! Finalize
                 !
                 intensity(inu0:inu1,1:4) = intensity(inu0:inu1,1:4) * 0.25d0
                 !
              endif
           endif
        endif
     endif
  endif
  !
end subroutine camera_compute_one_pixel


!-------------------------------------------------------------------------
!                   CHECK IF A STAR SPHERE IS INSIDE PIXEL
!-------------------------------------------------------------------------
function camera_check_starsphere_in_pixel(istar,px,py,dx,dy,inside)
  implicit none
  logical :: camera_check_starsphere_in_pixel,inside
  integer :: istar
  double precision :: x,y,z,dirx,diry,dirz,px,py,distance,pxs,pys,dx,dy,rs
  rs = star_r(istar)
  call camera_set_ray_stars_pntsrc(istar,x,y,z,dirx,diry,dirz,pxs,pys,distance)
  if(inside) then
     !
     ! Check if stellar disc is entirely inside this pixel (ignore the
     ! roundness)
     !
     if((pxs.gt.px-0.5d0*dx+rs).and.(pxs.lt.px+0.5d0*dx-rs).and. &
          (pys.gt.py-0.5d0*dy+rs).and.(pys.lt.py+0.5d0*dy-rs)) then
        camera_check_starsphere_in_pixel = .true.
     else
        camera_check_starsphere_in_pixel = .false.
     endif
  else
     !
     ! Check if stellar disc can contribute to this pixel (ignore the
     ! roundness)
     !
     if((pxs.gt.px-0.5d0*dx-rs).and.(pxs.lt.px+0.5d0*dx+rs).and. &
          (pys.gt.py-0.5d0*dy-rs).and.(pys.lt.py+0.5d0*dy+rs)) then
        camera_check_starsphere_in_pixel = .true.
     else
        camera_check_starsphere_in_pixel = .false.
     endif
  endif
  return
end function camera_check_starsphere_in_pixel



!-------------------------------------------------------------------------
!                         MAKE A RECTANGULAR IMAGE
!
! If img=0, then do just one image. If img>0, then use vantage point
!    settings and pointing settings from the multiple-images settings
!    array, entry number img. This allows you to make an entire movie
!    of viewing the object under different vantage points all at once.
!    Though for each of these images (each may have multiple wavelengths,
!    se above), you must call camera_make_rect_image separately, and
!    save the image separately. 
!
!-------------------------------------------------------------------------
subroutine camera_make_rect_image(img,tausurf)
  use amr_module
  implicit none
  double precision :: x,y,z,dxx,dyy,dzz,px,py
  integer :: inu,img,ix,iy,idxx,idxy,istar,ierr,inu0,inu1,ierror,ispec
  double precision :: pdx,pdy,dz,d2,d3,factor,celldxmin,quvsq
  double precision :: dirx,diry,dirz,distance,xbk,ybk,zbk,svcx,svcy,svcz
  logical :: domc,dotausurf
  logical, optional :: tausurf
  character*80 :: strint
  integer :: iact,icnt,ilinesub
  logical :: redo
  double precision :: seconds
  !
  ! If "tausurf" is set, then the purpose of this subroutine
  ! changes from being an imager to being a "tau=1 surface finder".
  ! Default is .false.
  !
  dotausurf = .false.
  if(present(tausurf)) then
     dotausurf = tausurf
  endif
  !
  ! If camera_diagnostics_subpix.eq..true. then we will make
  ! a dump of the subpixeling diagnostics, so that the user can see
  ! from the file "subpixeling_diagnostics.out" where the (sub-)pixels
  ! are put.
  !
  if(camera_diagnostics_subpix) then
     open(unit=10,file='subpixeling_diagnostics.out')
  endif
  !
  ! Local observer mode is incompatible with 1-D plane-parallel
  ! or 2-D pencil-parallel modes
  !
  if(camera_localobserver) then
     if(igrid_coord.eq.10) then 
        write(stdo,*) 'ERROR: In 1-D plane-parallel mode, no local observer is allowed.'
        stop
     endif
     if(igrid_coord.eq.20) then
        write(stdo,*) 'ERROR: In 2-D pencil-parallel mode, no local observer is allowed.'
        stop
     endif
  endif
  !
  ! If we have 1-D plane-parallel or 2-D pencil-parallel modes, we
  ! switch off the sub-pixeling. For 1-D this is not necessary anyway.
  ! For 2-D it might still be necessary, but we will leave it to the
  ! user to do it by hand - we'll warn.
  !
  if(igrid_coord.eq.10) then
     camera_nrrefine = -1
  endif
  if(igrid_coord.eq.20) then
     if(camera_nrrefine.gt.0) then
        write(stdo,*) 'WARNING: When using 2-D pencil-parallel coordinates,'
        write(stdo,*) '         RADMC-3D does not provide an automatic'
        write(stdo,*) '         subpixeling method (yet). In order to '
        write(stdo,*) '         ensure flux conservation you must take'
        write(stdo,*) '         care yourself to have enough resolution'
        write(stdo,*) '         to resolve all scales in the y-z plane.'
        write(stdo,*) '         In other words, we now switch to nofluxcons.'
     endif
     camera_nrrefine = -1
  endif
  !
  ! If we have 1-D plane-parallel mode, then we always do just a 1x1
  ! pixel image (it's useless to do more pixels)
  !
  if(igrid_coord.eq.10) then
     camera_image_nx = 1
     camera_image_ny = 1
  endif
  !
  ! If we have 2-D pencil-parallel mode, then we always do just a 1xN
  ! pixel image (it's useless to do more pixels)
  !
  if(igrid_coord.eq.20) then
     camera_image_nx = 1
  endif
  !
  ! Since we have two very different modes of observation (one with the 
  ! observer at infinity, in which the image sizes are specified in cm,
  ! and one with the observer local, in which the image sizes are specified
  ! in radians), we must do some checks for self-consistency, because
  ! the user may easily have forgotten to specify something.
  !
  if(img.eq.0) then
     if(camera_localobserver) then
        if(camera_image_halfsize_x.gt.10.d0) then
           write(stdo,*) 'ERROR in camera module: You have chosen the local observer perspective,'
           write(stdo,*) '      but the image size is still specified in cm instead of radian.'
           stop
        endif
        if((abs(camera_observer_position(1)).ge.1d90).or. &
             (abs(camera_observer_position(2)).ge.1d90).or. &
             (abs(camera_observer_position(3)).ge.1d90)) then
           write(stdo,*) 'ERROR in camera module: You have chosen the local observer perspective,'
           write(stdo,*) '      but you have not yet specified the 3-D position of the observer.'
           stop
        endif
     else
        if((igrid_coord.ne.10).and.(igrid_coord.ne.20)) then
           if(camera_image_halfsize_x.lt.10.d0) then
              write(stdo,*) 'ERROR in camera module: You have chosen the observer at infinity perspective,'
              write(stdo,*) '      but the image size is still specified in radian instead of cm (or AU or pc).'
              stop
           endif
        endif
     endif
  endif
  !
  ! Also let the Monte Carlo module know what mode we have (local observer
  ! or not)
  !
  mcscat_localobserver = camera_localobserver
  if(mcscat_localobserver) then
     mcscat_localobs_pos(1:3) = camera_observer_position(1:3)
  endif
  !
  ! If the dust emission is included, then make sure the dust data,
  ! density and temperature are read. If yes, do not read again.
  !
  if(rt_incl_dust) then
     call read_dustdata(1)
     call read_dust_density(1)
     call read_dust_temperature(1)
     if(alignment_mode.ne.0) then
        call aligned_grains_init(1)
     endif
  endif
  !
  ! If line emission is included, then make sure the line data are
  ! read. If yes, then do not read it again.
  !
  if(rt_incl_lines) then
     call read_lines_all(1)
  endif
  !
  ! Do a check
  !
  if(rt_incl_lines.and.(lines_maxshift.le.0.d0)) then
     write(stdo,*) 'INTERNAL ERROR in lines: lines_maxshift variable not set'
     stop
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
  ! Note: it only stores those levels which are selected in the "subset"
  ! (see lines_module.f90). The idea of the subset is that to save memory
  ! you may not want to store all level populations.  Example: A molecule
  ! may have 30 relevant levels. If you need to store all populations
  ! globally, then you need 30 x 8 bytes x nrcells of memory. For large
  ! grids that could be very much. If, however, you wish to only model one
  ! of the lines, then only 2 level populations have to be stored at each
  ! point (the upper and the lower level belonging to that line). If you
  ! select these 2 levels as your "subset" then only these 2 levels will be
  ! stored in the global lines_levelpop() array, saving a lot of memory. For
  ! the LVG method all levels are needed *locally* to compute the
  ! populations, but once these populations are computed, only the 2
  ! relevant ones are then stored in the *global* array lines_levelpop().
  !
  ! You can select this subset manually (in the lines.inp file) or
  ! automatically (by calling the subroutine
  ! lines_automatic_subset_selection() in the lines_module.f90).  The latter
  ! is done just above here.
  !
  if(rt_incl_lines) then
     if((lines_mode.ge.1).and.(lines_mode.le.9)) then
        !
        ! Level subset selection
        !
        if(lines_autosubset) then
           !
           ! If requested, do an automatic subset selection of the
           ! molecular levels
           ! 
           call lines_automatic_subset_selection(camera_nrfreq,     &
                camera_frequencies,1,camera_nrfreq,lines_maxshift,  &
                redo)
           !
           ! Now force a recomputation of the populations, unless
           ! "redo" is .false., meaning that we have the same
           ! levels as before (and thus the population as computed
           ! before must be still correct). Note that if there
           ! exists no lines_levelpop() array, then redo will also
           ! be .true.
           !
           if(redo) then
              iact=2
           else
              iact=1
           endif
           call lines_compute_and_store_local_populations(iact)
        else
           !
           ! Subset is not automatically selected. So either it has
           ! been selected manually or not at all (in which case the
           ! "subset" is the full set of levels).
           !
           write(stdo,*) 'Will store level populations for all levels or ',   &
                'for a manually (in lines.inp) selected subset of levels.'
           !
           ! Now compute the populations only if they have not yet been
           ! computed. If the lines_popul() array is present, we know that
           ! it contains the correct populations, because the subset
           ! is fixed by the user (or is the complete set). 
           !
           call lines_compute_and_store_local_populations(1)
           !
        endif
     endif
  endif
  !
  ! Do the initializing of the camera
  !
  call camera_init()
  !
  ! If tau surface mode, then allocate the following arrays, too.
  ! These give back the 3-D positions of the points on the tau surface.
  ! The projected position is put into the camera_rect_image_iquv() array.
  !
  if(dotausurf) then
     if(allocated(camera_tausurface_z)) deallocate(camera_tausurface_z)
     if(allocated(camera_tausurface_y)) deallocate(camera_tausurface_y)
     if(allocated(camera_tausurface_x)) deallocate(camera_tausurface_x)
     allocate(camera_tausurface_x(1:camera_image_nx,1:camera_image_ny,1:camera_nrfreq))
     allocate(camera_tausurface_y(1:camera_image_nx,1:camera_image_ny,1:camera_nrfreq))
     allocate(camera_tausurface_z(1:camera_image_nx,1:camera_image_ny,1:camera_nrfreq))
  endif
  !
  ! Warnings
  !
  if(rt_incl_lines.and.(lines_mode.lt.0).and.(scattering_mode.ne.0)) then
     write(stdo,*) 'WARNING: Using dust scattering AND line transfer with '
     write(stdo,*) '         on-the-fly level population determination can'
     write(stdo,*) '         make RADMC-3D very slow, because dust scattering'
     write(stdo,*) '         means that the transfer must be done freq-by-freq,'
     write(stdo,*) '         and thus the populations must be recalculated'
     write(stdo,*) '         at each freq, which is slowing down the code.'
     write(stdo,*) '         Various possible solutions:'
     write(stdo,*) '          1. If possible, switch off dust scattering'
     write(stdo,*) '             (for instance for far-IR or mm wavelengths)'
     write(stdo,*) '             by setting scattering_mode_max=0 in radmc3d.inp'
     write(stdo,*) '          2. Use not-on-the-fly populations (lines_mode>0)'
     write(stdo,*) '             [at the moment this mode is still in prep]'
  endif
  !
  ! Check
  !
  if(camera_nrfreq.lt.1) then
     write(stdo,*) 'ERROR in camera module: camera frequency array'
     write(stdo,*) '      not set when calling camera_make_rect_image()'
     stop
  endif
  if((camera_image_nx.le.0).or.(camera_image_ny.le.0)) then
     write(stdo,*) 'ERROR in camera module: must have >=1 values for'
     write(stdo,*) '    camera_image_nx and camera_image_ny'
     stop
  endif
  if((camera_refine_criterion.le.0.d0).and.(camera_nrrefine.gt.0)) then
     write(stdo,*) 'ERROR in camera module: must set the camera_refine_criterion'
     write(stdo,*) '      when making flux-conserving images.'
     stop
  endif
  if(img.gt.0) then
     if(.not.allocated(cameras_pt_pos)) then
        write(stdo,*) 'ERROR in camera module: multiple image arrays not'
        write(stdo,*) '      yet allocated.'
        stop
     endif
     if(img.gt.cameras_nr_images) then
        write(stdo,*) 'ERROR in camera module: asked to make image nr ',img,&
             ' but that is beyond cameras_nr_images.'
        stop
     endif
  endif
  !
  ! If using multiple images, then we copy all the info to the relevant
  ! arrays here.
  !
  if(img.gt.0) then
     !
     ! Get the data from the movie arrays, and perform some basic
     ! checks.
     !
     camera_pointing_position(1:3)   = cameras_pt_pos(1:3,img)    
     camera_image_halfsize_x         = cameras_img_hs_x(img)    
     camera_image_halfsize_y         = cameras_img_hs_y(img)    
     camera_zoomcenter_x             = cameras_zmc_x(img)
     camera_zoomcenter_y             = cameras_zmc_y(img)
     camera_pointing_degr_posang     = cameras_pt_degr_pa(img)  
     if(camera_localobserver) then
        camera_observer_position(1:3)   = cameras_obs_pos(1:3,img)   
        if((camera_image_halfsize_x.gt.10.d0).or. &
           (camera_image_halfsize_y.gt.10.d0)) then
           write(stdo,*) 'ERROR in camera module: You have chosen the local observer perspective,'
           write(stdo,*) '      but the image size is still specified in cm instead of radian.'
           stop
        endif
     else
        camera_observer_degr_theta      = cameras_obs_degr_th(img) 
        camera_observer_degr_phi        = cameras_obs_degr_ph(img) 
        if((camera_image_halfsize_x.lt.10.d0).or. &
           (camera_image_halfsize_y.lt.10.d0)) then
           write(stdo,*) 'ERROR in camera module: You have chosen the observer at infinity perspective,'
           write(stdo,*) '      but the image size is still specified in radian instead of cm (or AU or pc).'
           stop
        endif
     endif
  endif
  !
  ! Check
  !
  if((camera_image_halfsize_x.le.0.d0).or.(camera_image_halfsize_y.le.0.d0)) then
     write(stdo,*) 'ERROR in camera module: zero image scale'
     stop
  endif
  !
  ! For now the flux-conserving ray-tracing is not yet implemented for
  ! local-observer perspective
  !
  if(camera_localobserver) then
     if(camera_nrrefine.gt.0) then
        write(stdo,*) 'ERROR: For the local observer viewing mode the '
        write(stdo,*) '       flux-conserving ray-tracing is not yet available.'
        stop
     endif
  endif
  !
  ! Pre-compute cosines and sines
  !
  if(camera_localobserver) then
     !
     ! For the local perspective view the angles must first be calculated
     !
     camera_pointing_cos_posang  = cos(camera_pointing_degr_posang*pi/180.)
     camera_pointing_sin_posang  = sin(camera_pointing_degr_posang*pi/180.)
     dxx = camera_pointing_position(1) - camera_observer_position(1)
     dyy = camera_pointing_position(2) - camera_observer_position(2)
     dzz = camera_pointing_position(3) - camera_observer_position(3)
     d3 = sqrt( dxx**2 + dyy**2 + dzz**2 )
     d2 = sqrt( dxx**2 + dyy**2 )
     if(d3.eq.0.d0) then
        write(stdo,*) 'ERROR in camera module: observer and pointing positions'
        write(stdo,*) '      are identical. Cannot computing pointing direction.'
        stop
     endif
     camera_observer_sin_theta   = d2/d3
     camera_observer_cos_theta   = sqrt(1.d0-camera_observer_sin_theta**2)
     if(camera_observer_position(3).lt.camera_pointing_position(3)) then
        camera_observer_cos_theta   = -camera_observer_cos_theta
     endif
     if(d2.gt.0.d0) then
        camera_observer_sin_phi     = dxx/d2
        camera_observer_cos_phi     = dyy/d2
     else
        camera_observer_sin_phi     = 0.d0
        camera_observer_cos_phi     = 1.d0
     endif
  else
     !
     ! For the camera-at-infinity view the angles are given and the
     ! cos and sin can be calculated easily
     !
     camera_pointing_cos_posang  = cos(camera_pointing_degr_posang*pi/180.)
     camera_pointing_sin_posang  = sin(camera_pointing_degr_posang*pi/180.)
     camera_observer_cos_theta   = cos(camera_observer_degr_theta*pi/180.)
     camera_observer_sin_theta   = sin(camera_observer_degr_theta*pi/180.)
     camera_observer_cos_phi     = cos(camera_observer_degr_phi*pi/180.)
     camera_observer_sin_phi     = sin(camera_observer_degr_phi*pi/180.)
  endif
  !
  ! For now we do not allow local observer if scattering is switched on
  !
  if((scattering_mode.gt.0).and.camera_localobserver) then
     write(stdo,*) 'ERROR: Monte Carlo scattering not yet implemented with local observer...'
     stop
  endif
  !
  ! For now we do not allow multiple vantage points
  !
  if(mcscat_nrdirs.gt.1) then
     write(stdo,*) 'ERROR: Multiple vantage points is not yet allowed in this version!'
     stop
  endif
  !
  ! Check if we need to (re)do the one-frequency scattering Monte Carlo simulation
  !
  domc = .false.
  if(scattering_mode.ge.1) then
     !
     ! For now scattering is not allowed in local observer mode
     !
     if(camera_localobserver) then
        write(stdo,*) 'ERROR: At the moment scattering is not'
        write(stdo,*) '       allowed for local observer mode.'
        stop
     endif
     !
     ! Compute the direction vector of the observer. Start with (0,0,1).
     ! Actually this is only necessary for scattering_mode.ge.2, but alas.
     !
     dirx = 0.d0
     diry = 0.d0
     dirz = 1.d0
     !
     ! A rotation in the y,z plane, i.e. going from pole-on
     ! to more face-on (if a disk is assumed to be present in the x-y plane)
     !
     ybk  = diry
     zbk  = dirz
     diry = camera_observer_cos_theta * ybk - camera_observer_sin_theta * zbk
     dirz = camera_observer_sin_theta * ybk + camera_observer_cos_theta * zbk
     !
     ! Then a rotation in the x,y plane, rotating horizontally around
     ! the object. If the camera looks toward the object, then the
     ! camera now moves to the left (clockwise around the object).
     !
     xbk  = dirx
     ybk  = diry
     dirx = camera_observer_cos_phi * xbk + camera_observer_sin_phi * ybk
     diry =-camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
     !
     ! Checks
     !
     if(abs(dirx**2+diry**2+dirz**2-1.d0).gt.1d-6) then
        write(stdo,*) 'ERROR in camera module: direction vector not OK.'
        write(stdo,*) dirx,diry,dirz
        stop
     endif
     !
     ! Compute the S-vector of observer for polarization. Start with (0,1,0).
     ! Actually this is only necessary for scattering_mode.ge.4, but it
     ! never hurts.
     !
     svcx = 0.d0
     svcy = 1.d0
     svcz = 0.d0
     !
     ! Rotate camera along its axis 
     !
     ! Note that the camera is rotated in clockwise direction, so any
     ! image is rotated in counter-clockwise direction on the CCD
     !
     xbk  = svcx
     ybk  = svcy
     svcx = camera_pointing_cos_posang * xbk + camera_pointing_sin_posang * ybk
     svcy =-camera_pointing_sin_posang * xbk + camera_pointing_cos_posang * ybk
     !
     ! Then a rotation in the y,z plane, i.e. going from pole-on
     ! to more face-on (if a disk is assumed to be present in the x-y plane)
     !
     ybk  = svcy
     zbk  = svcz
     svcy = camera_observer_cos_theta * ybk - camera_observer_sin_theta * zbk
     svcz = camera_observer_sin_theta * ybk + camera_observer_cos_theta * zbk
     !
     ! Then a rotation in the x,y plane, rotating horizontally around
     ! the object. If the camera looks toward the object, then the
     ! camera now moves to the left (clockwise around the object).
     !
     xbk  = svcx
     ybk  = svcy
     svcx = camera_observer_cos_phi * xbk + camera_observer_sin_phi * ybk
     svcy =-camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
     !
     ! Checks
     !
     if(abs(svcx**2+svcy**2+svcz**2-1.d0).gt.1d-6) then
        write(stdo,*) 'ERROR in camera module: S-vector not OK.'
        write(stdo,*) svcx,svcy,svcz
        stop
     endif
     if(abs(svcx*dirx+svcy*diry+svcz*dirz).gt.1d-6) then
        write(stdo,*) 'INTERNAL ERROR: Somehow the S-vector is not '
        write(stdo,*) '  perpendicular to the direction vector... Warn author.'
        stop
     endif
     !
     ! Some simple tests to see if a new MC calculation is necessary
     !
     if((.not.allocated(mcscat_dirs)).and.(scattering_mode.gt.1)) domc=.true.
     if(allocated(mcscat_dirs).and.(scattering_mode.eq.1)) then
        write(stdo,*) 'ERROR: if isotropic scattering is set, then do not set the mcscat_dirs(:,:) array...'
        stop
     endif
     if(mcscat_nrdirs.ne.1) domc=.true.
     if(.not.camera_scatsrc_allfreq) domc=.true.
     !
     ! Some more subtle tests to see if a new MC calculation is necessary
     !
     if(camera_scatsrc_allfreq) then
        if(mc_nrfreq.ne.camera_nrfreq) then
           domc=.true.
        else
           do inu=1,camera_nrfreq
              if(abs((camera_frequencies(inu)-mc_frequencies(inu))/   &
                     (camera_frequencies(inu)+mc_frequencies(inu))).gt.1d-10) domc=.true.
           enddo
        endif
        if(allocated(mcscat_dirs)) then
           mcscat_current_dir = 1
           if((abs(mcscat_dirs(1,mcscat_current_dir)-dirx).gt.1d-6).or.   &
              (abs(mcscat_dirs(2,mcscat_current_dir)-diry).gt.1d-6).or.   &
              (abs(mcscat_dirs(3,mcscat_current_dir)-dirz).gt.1d-6)) then
              domc=.true.
           endif
        endif
     endif
  endif
  !
  ! Prepare the direction for the scattering source function
  ! For now allow only one single direction (i.e. one single vantage point)
  !
  if(domc.and.(scattering_mode.gt.1)) then
     if((igrid_coord.ge.100).and.(amr_dim.ne.3)) then
        !
        ! Special case: 2D axisymmetric model in spherical coordinate
        ! but with full scattering mode (only possible for scattering_mode.ge.5).
        !
        if(amr_dim.eq.1) then
           write(stdo,*) 'ERROR: scattering_mode.ge.2 is incompatible with'
           write(stdo,*) '       1-D spherical coordinates.'
           stop
        endif
        if(scattering_mode.lt.5) then
           write(stdo,*) 'ERROR: scattering_mode.lt.5 is incompatible with'
           write(stdo,*) '       2-D spherical coordinates.'
           stop
        endif
        if(camera_secondorder) then
           write(stdo,*) 'ERROR: At the moment the 2-D axisymmetric full-scattering mode is not yet'
           write(stdo,*) '       compatible with second order ray-tracing... :-('
           stop
        endif
        !
        ! Switch on the special treatment
        !
        dust_2daniso = .true.
        !
        ! Extend the scattering source array to dust_2daniso_nphi + 1 phi-angle
        ! points starting with phi=0 and ending with phi=360 degrees
        !
        mcscat_nrdirs = dust_2daniso_nphi + 1
        !
        ! Make sure that the ray tracing cannot make larger steps
        ! than a maximum angle wrt the origin. Reason: for near edge-on
        ! views a ray could otherwise pass through a cell (=annulus) 
        ! almost along the annulus tube, and change phi angle too much,
        ! thereby skipping intermediate angles. That is bad for the 
        ! interpolation of the scattering source function.
        !
        if(camera_maxdphi.eq.0.d0) then
           write(stdo,*) 'WARNING: 2-D anisotropic scattering without camera_maxdphi set... Dangerous.'
        endif
        if(camera_maxdphi.gt.twopi/dust_2daniso_nphi) then
           camera_maxdphi = twopi / dust_2daniso_nphi
        endif
        if(camera_maxdphi.ge.0.5d0) camera_maxdphi=0.5d0
        !
        ! Message
        !
        write(stdo,*) 'Note: Using 2-D full-phase scattering mode. This requires a bit of extra memory.'
     else
        !
        ! Normal case (3-D)
        !
        mcscat_nrdirs = 1
     endif
     !
     ! Make the scattering direction the same as the direction of
     ! viewing. 
     ! For now allow only one single direction (i.e. one single vantage point)
     !
     ! NOTE: For 2-D axisymmetric models in spherical coordinates, we
     !       do things in a special way: we reserve mcscat_nrdirs 
     !       "directions", which are in fact all the same, but for
     !       the scattering source function we will set the scattering
     !       event at different positions. 
     !
     if(allocated(mcscat_dirs)) deallocate(mcscat_dirs)
     allocate(mcscat_dirs(1:3,1:mcscat_nrdirs),STAT=ierr)
     if(ierr.gt.0) then
        write(stdo,*) 'ERROR: Could not allocate mcscat_dirs()'
        stop
     endif
     mcscat_current_dir = 1
     mcscat_dirs(1,1:mcscat_nrdirs) = dirx
     mcscat_dirs(2,1:mcscat_nrdirs) = diry
     mcscat_dirs(3,1:mcscat_nrdirs) = dirz
     !
     ! For convenience, store this also here in the camera module
     ! (works only for a single vantage point)
     !
     camera_dir(1) = dirx
     camera_dir(2) = diry
     camera_dir(3) = dirz
     !
     ! In case you want to include polarization, we set also the
     ! S-vectors, which will be perpendicular to the direction 
     ! vector, and pointing vertically upward in the image plane.
     !
     if(allocated(mcscat_svec)) deallocate(mcscat_svec)
     allocate(mcscat_svec(1:3,1:mcscat_nrdirs),STAT=ierr)
     if(ierr.gt.0) then
        write(stdo,*) 'ERROR: Could not allocate mcscat_svec()'
        stop
     endif
     mcscat_svec(1,1:mcscat_nrdirs) = svcx
     mcscat_svec(2,1:mcscat_nrdirs) = svcy
     mcscat_svec(3,1:mcscat_nrdirs) = svcz
     !
     ! For convenience, store this also here in the camera module
     ! (works only for a single vantage point)
     !
     camera_svec(1) = svcx
     camera_svec(2) = svcy
     camera_svec(3) = svcz
     !
     ! But mirror symmetry is not allowed for anisotropic scattering
     !
     if((igrid_coord.ge.100).and.(igrid_coord.le.199)) then
        if(igrid_mirror.ne.0) then
           write(stdo,*) 'ERROR: Mirror symmetry not compatible with anisotropic scattering.'
           stop
        endif
     endif
  endif
  !
  ! Compute pixel size.
  !
  ! NOTE: For the observer-at-infinity mode this is in centimeters.
  !       For the local-observer mode this is in radian.
  !
  pdx = 2*camera_image_halfsize_x / (1.d0*camera_image_nx)
  pdy = 2*camera_image_halfsize_y / (1.d0*camera_image_ny)
  !
  ! If the tausurface mode is on, then allocate some arrays
  !
  if(dotausurf) then
     if(allocated(camera_dstop)) deallocate(camera_dstop)
     if(allocated(camera_zstop)) deallocate(camera_zstop)
     if(allocated(camera_ystop)) deallocate(camera_ystop)
     if(allocated(camera_xstop)) deallocate(camera_xstop)
     if(allocated(camera_taustop)) deallocate(camera_taustop)
     allocate(camera_taustop(camera_nrfreq))
     allocate(camera_xstop(camera_nrfreq))
     allocate(camera_ystop(camera_nrfreq))
     allocate(camera_zstop(camera_nrfreq))
     allocate(camera_dstop(camera_nrfreq))
  endif
  !
  ! Reset flag
  !
  camera_warn_resolution = .false.
  !
  ! Message
  !
  write(stdo,*) 'Rendering image(s)...'
  !
  ! Now make the image. We can do so in two different ways. The first is to
  ! do all frequencies (wavelengths) at once. If scattering is included this
  ! requires the precomputing of the scattering source function at all
  ! wavelengths.  This can have a too large memory requirement.  The second
  ! is to do it frequency-by-frequency, and doing a MC calculation for each
  ! of them separately. That saves a lot of computer memory, so it is the
  ! preferred method for large simulations when making spectra. But it may
  ! be slower. The way to choose one or the other is with the switch called
  ! camera_scatsrc_allfreq.
  !
  if(((scattering_mode.eq.0).or.camera_scatsrc_allfreq).and.   &
      (.not.camera_secondorder)) then
     !
     ! Multi-wavelength Method 1:
     !
     ! Do all frequencies at once for each pixel. This is the preferred
     ! method for cases without scattering. With scattering this requires a
     ! lot of memory to store the scattering source function for all
     ! frequencies, so it might not be the best for that.
     !
     ! Note that at present this multi-frequency mode is not available
     ! for second order integration of the transfer equation (i.e. with
     ! the corner points). 
     !
     inu0 = 1
     inu1 = camera_nrfreq
     !
     ! Tell the ray tracer that if scattering is done, then it is done
     ! frequency dispersed
     !
     camera_mcscat_monochromatic = .false.
     !
     ! If necessary, then do the scattering source functions at all
     ! frequencies beforehand.  WARNING: This can be a very large array!
     !
     !$ seconds = omp_get_wtime()
     if(domc) then
        if(allocated(mc_frequencies)) deallocate(mc_frequencies)
        mc_nrfreq=camera_nrfreq
        allocate(mc_frequencies(1:mc_nrfreq),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate mc_frequencies(:) array'
           stop
        endif
        mc_frequencies(:) = camera_frequencies(:)
        !
        ! Message
        !
        write(stdo,*) 'Doing scattering Monte Carlo simulation...'
        call flush(stdo)
        !
        ! Do Monte Carlo simulation
        !
        if(camera_lambda_starlight_single_scat_mode.eq.0) then
           call do_monte_carlo_scattering(rt_mcparams,ierror,do_resetseed,&
                                          scatsrc=.true.)
           write(stdo,*) 'Average number of scattering events per photon package = ', &
                      ieventcounttot/(1.d0*rt_mcparams%nphot_scat)
           !
           ! If the Monte Carlo settings are very conservative, then give a warning
           ! that you may want to change this (but at your own risk). 
           !
           if(mc_scat_maxtauabs.gt.5.d0) then
              write(stdo,*) 'Tip for speed-up: By default the settings of RADMC-3D are conservative (i.e. safe but slow).'
              write(stdo,*) '   A photon package in monochromatic Monte Carlo is only destroyed after tau_abs = ',mc_scat_maxtauabs
              write(stdo,*) '   In most cases, however, an optical depth limit of 5 is enough.'
              write(stdo,*) '   You can (though at your own risk) speed this up by adding the following line to radmc3d.inp:'
              write(stdo,*) '   mc_scat_maxtauabs = 5.d0'
           endif
        elseif(camera_lambda_starlight_single_scat_mode.eq.1) then
           call do_lambda_starlight_single_scattering(rt_mcparams,ierror,scatsrc=.true.)
        elseif(camera_lambda_starlight_single_scat_mode.eq.2) then
           call do_lambda_starlight_single_scattering_simple(rt_mcparams,ierror,scatsrc=.true.)
        else
           write(stdo,*) 'Lambda single scattering mode cannot be other than 0 or 1 or 2 for now.'
           stop 8762
        endif
     endif
     !$ write(stdo,*)"Total elapsed time:",omp_get_wtime() - seconds;
     !
     ! Pre-compute which lines and which levels for line transfer may
     ! contribute to these wavelengths. Note that this only has to be
     ! pre-computed if the serial ray tracing is used, but we do it 
     ! nevertheless, as it will not hurt.
     !
     if(rt_incl_lines) then
        ! ------------------------------------------------------------------
        ! Attila Juhasz
        ! lines_find_active_lines_leves -> if all energy levels are known - leiden
        ! format for line data
        ! lines_find_active_lines_linelist -> for linelist mode
        ! ------------------------------------------------------------------
        !   call lines_find_active_lines_levels(camera_nrfreq,             &
        !                   camera_frequencies,inu0,inu1,lines_maxshift)
        !
        if(lines_maxnrlevels.gt.0) then
           call lines_find_active_lines_levels(camera_nrfreq,             &
                camera_frequencies,inu0,inu1,lines_maxshift)
        else
           call lines_find_active_lines_linelist(camera_nrfreq,             &
                camera_frequencies,inu0,inu1,lines_maxshift)
        endif
        ! ------------------------------------------------------------------
     endif
     !
     ! Message
     !
     if(inu0.eq.inu1) then
        write(stdo,*) 'Ray-tracing image for lambda = ', &
                1d4*cc/camera_frequencies(inu0),' micron...'
     else
        call integer_to_string(abs(inu1-inu0)+1,strint)
        write(stdo,*) 'Ray-tracing images: all '//trim(strint)//' wavelength at once...'
     endif
     !
     ! Now make the image at all wavelengths simultaneously
     !
     call camera_make_rect_image_sub(inu0,inu1)
     !
  else
     !
     ! Multi-wavelength Method 2:
     !
     ! If scattering is included, we need the scattering source function. Since this 
     ! array is easily too big to be stored for all wavelengths: In this method we
     ! compute the scattering source function for each frequency separately, and then
     ! make the corresponding image, and then go to the next frequency.
     !
     ! Tell the ray tracer that if scattering is done, then it is done
     ! one frequency at a time, i.e. the scattering source function is
     ! just one frequency bin.
     !
     camera_mcscat_monochromatic = .true.
     !
     ! Since we recompute the scattering source function for each wavelength, the domc
     ! is by default true. Note that if you make an image at a single wavelength you 
     ! may want to be able to make use of a previously computed source function, if for
     ! instance you make your next image at the same vantage point, the same wavelength
     ! but at a different zoom factor or so. In that case, please set the flag
     ! camera_scatsrc_allfreq to .true.. 
     !
     if(domc) then
        if(allocated(mc_frequencies)) deallocate(mc_frequencies)
        mc_nrfreq=1
        allocate(mc_frequencies(1:1),STAT=ierr)
     endif
     !
     ! Do a loop over all camera frequencies
     !
     do inu0=1,camera_nrfreq
        !
        ! Set this by default
        !
        inu1=inu0
        !
        ! If we must do Monte Carlo, then do this
        !
        if(domc) then
           !
           ! Set the wavelength for the Monte Carlo scattering simulation
           !
           mc_frequencies(1) = camera_frequencies(inu0)
           !
           ! Message
           !
           write(stdo,*) 'Doing scattering Monte Carlo simulation for lambda = ', &
                1d4*cc/mc_frequencies(1),' micron...'
           call flush(stdo)
           !
           ! Call the single wavelength Monte Carlo module
           !
           if(camera_lambda_starlight_single_scat_mode.eq.0) then
              call do_monte_carlo_scattering(rt_mcparams,ierror,do_resetseed,&
                                             scatsrc=.true.)
              write(stdo,*) 'Average number of scattering events per photon package = ', &
                      ieventcounttot/(1.d0*rt_mcparams%nphot_scat)
              !
              ! If the Monte Carlo settings are very conservative, then give a warning
              ! that you may want to change this (but at your own risk). 
              !
              if(mc_scat_maxtauabs.gt.5.d0) then
                 write(stdo,*) 'Tip for speed-up: By default the settings of RADMC-3D are conservative ', &
                      '(i.e. safe but slow).'
                 write(stdo,'(A68,A16,F6.2)') '   A photon package in monochromatic Monte Carlo is only destroyed ', &
                      'after tau_abs = ',mc_scat_maxtauabs
                 write(stdo,*) '   In most cases, however, an optical depth limit of 5 is enough.'
                 write(stdo,*) '   You can (though at your own risk) speed this up by adding the following ', &
                      'line to radmc3d.inp:'
                 write(stdo,*) '   mc_scat_maxtauabs = 5.d0'
              else
                 if(mc_scat_maxtauabs.gt.2.d0) then
                    write(stdo,'(A36,F6.2,A36)') ' Warning: Using mc_scat_maxtauabs = ',mc_scat_maxtauabs, &
                         ' (this is fine, but be aware of it).'
                 else
                    write(stdo,'(A36,F6.2,A34)') ' ERROR: Using mc_scat_maxtauabs = ',mc_scat_maxtauabs, &
                         ': This is too low...'
                 endif
              endif
           elseif(camera_lambda_starlight_single_scat_mode.eq.1) then
              call do_lambda_starlight_single_scattering(rt_mcparams,ierror,scatsrc=.true.)
           elseif(camera_lambda_starlight_single_scat_mode.eq.2) then
              call do_lambda_starlight_single_scattering_simple(rt_mcparams,ierror,scatsrc=.true.)
           else
              write(stdo,*) 'Lambda single scattering mode cannot be other than 0 or 1 or 2 for now.'
              stop 8762
           endif
        endif
        !
        ! Pre-compute which lines and which levels for line transfer may
        ! contribute to these wavelengths. Note that this only has to be
        ! pre-computed if the serial ray tracing is used, but we do it 
        ! nevertheless, as it will not hurt.
        !
        if(rt_incl_lines) then
           !
           ! Check which lines to include
           !
           ! ------------------------------------------------------------------
           ! Attila Juhasz
           ! lines_find_active_lines_leves -> if all energy levels are known - leiden
           ! format for line data
           ! lines_find_active_lines_linelist -> for linelist mode
           ! ------------------------------------------------------------------
           !   call lines_find_active_lines_levels(camera_nrfreq,             &
           !                   camera_frequencies,inu0,inu1,lines_maxshift)
           !
            if(lines_maxnrlevels.gt.0) then
                call lines_find_active_lines_levels(camera_nrfreq,             &
                     camera_frequencies,inu0,inu1,lines_maxshift)
            else
                call lines_find_active_lines_linelist(camera_nrfreq,             &
                     camera_frequencies,inu0,inu1,lines_maxshift)
            endif

           ! ------------------------------------------------------------------
           !
           ! If camera_catch_doppler_line.eq..true., then allocate the
           ! corner-based line quantities
           !
           if(camera_catch_doppler_line) then
              if(allocated(sources_vertex_line_nup)) then
                 deallocate(sources_local_line_ndown_curr)
                 deallocate(sources_local_line_nup_curr)
                 deallocate(sources_local_line_ndown_prev)
                 deallocate(sources_local_line_nup_prev)
                 deallocate(sources_local_line_ndown_end)
                 deallocate(sources_local_line_nup_end)
                 deallocate(sources_cell_line_ndown)
                 deallocate(sources_cell_line_nup)
                 deallocate(sources_vertex_line_ndown)
                 deallocate(sources_vertex_line_nup)
              endif
              sources_vertex_lines_nractivetot = 0
              do ispec=1,lines_nr_species
                 sources_vertex_lines_nractivetot = sources_vertex_lines_nractivetot + &
                      active_nrlines(ispec)
              enddo
              allocate(sources_vertex_line_nup(1:sources_vertex_lines_nractivetot,1:amr_nr_vertices_max))
              allocate(sources_vertex_line_ndown(1:sources_vertex_lines_nractivetot,1:amr_nr_vertices_max))
              allocate(sources_cell_line_nup(1:sources_vertex_lines_nractivetot,1:amr_nr_vertices_max))
              allocate(sources_cell_line_ndown(1:sources_vertex_lines_nractivetot,1:amr_nr_vertices_max))
              allocate(sources_local_line_nup_curr(1:sources_vertex_lines_nractivetot))
              allocate(sources_local_line_ndown_curr(1:sources_vertex_lines_nractivetot))
              allocate(sources_local_line_nup_prev(1:sources_vertex_lines_nractivetot))
              allocate(sources_local_line_ndown_prev(1:sources_vertex_lines_nractivetot))
              allocate(sources_local_line_nup_end(1:sources_vertex_lines_nractivetot))
              allocate(sources_local_line_ndown_end(1:sources_vertex_lines_nractivetot))
           endif
        endif
        !
        ! Message
        !
        write(stdo,*) 'Ray-tracing image for lambda = ', &
                1d4*cc/camera_frequencies(inu0),' micron...'
        !
        ! If we do second order integration of the transfer equation,
        ! we must compute the emissivities at the corner points of the
        ! cells (the vertex grid). 
        !
        if(camera_secondorder) then
           !
           ! Set some flags and values
           !
           sources_localobserver = camera_localobserver
           !
           ! Do some preparations 
           !
           if(.not.camera_localobserver) then
              !
              ! Set the direction vector, if observer is at infinity
              ! 
              ! Note: this is necessary only for the inclusion of the scattering
              !       source function from the Monte Carlo module. This is
              !       anyway unavailable for local observer mode. 
              !
              dirx = 0.d0
              diry = 0.d0
              dirz = 1.d0
              ybk  = diry
              zbk  = dirz
              diry = camera_observer_cos_theta * ybk - camera_observer_sin_theta * zbk
              dirz = camera_observer_sin_theta * ybk + camera_observer_cos_theta * zbk
              xbk  = dirx
              ybk  = diry
              dirx = camera_observer_cos_phi * xbk + camera_observer_sin_phi * ybk
              diry =-camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
              if(abs(dirx**2+diry**2+dirz**2-1.d0).gt.1d-6) then
                 write(stdo,*) 'ERROR in camera module: direction vector not OK.'
                 write(stdo,*) dirx,diry,dirz
                 stop
              endif
              ray_cart_dirx = dirx
              ray_cart_diry = diry
              ray_cart_dirz = dirz
           else
              !
              ! Local observer mode, so set the observer position
              !
              sources_observer_position(:) = camera_observer_position(:)
           endif
           !
           ! Then call the subroutine
           !
           call sources_compute_snualphanu_at_vertices(inu0,camera_stokesvector)
        endif
        !
        ! Now make the image for this wavelength only
        !
        call camera_make_rect_image_sub(inu0,inu1)
        !
     enddo
     !
  endif
  !
  ! Check if refinement was sufficient
  !
  if(camera_warn_resolution) then
     write(stdo,*) 'WARNING: The currently produced image is not '
     write(stdo,*) '         guaranteed to have the correct flux '
     write(stdo,*) '         because it does not use the flux'
     write(stdo,*) '         conserving recursive pixeling.'
     write(stdo,*) '    Tip: Use fluxcons argument in command line.'
  endif
  !
  ! Close (if necessary) the "subpixeling_diagnostics.out" file
  ! and switch this diagnostics off, so that it won't do this each
  ! image in case you render a spectrum or SED.
  !
  if(camera_diagnostics_subpix) then
     close(10)
     camera_diagnostics_subpix = .false.
  endif
  !
  ! If the tausurface mode is on, then deallocate some arrays
  !
  if(dotausurf) then
     if(allocated(camera_dstop)) deallocate(camera_dstop)
     if(allocated(camera_zstop)) deallocate(camera_zstop)
     if(allocated(camera_ystop)) deallocate(camera_ystop)
     if(allocated(camera_xstop)) deallocate(camera_xstop)
     if(allocated(camera_taustop)) deallocate(camera_taustop)
  endif
  !
  !--------------------------------------------------------------
  !        A sub-subroutine for making the image
  !--------------------------------------------------------------
  !
  contains
  subroutine camera_make_rect_image_sub(inu00,inu11)
    implicit none
    integer :: inu00,inu11
    integer :: inuu
    integer :: backup_nrrefine,backup_tracemode
    logical :: warn_tausurf_problem,flag_quv_too_big
    integer :: id,nthreads
    double precision :: seconds
    integer :: pixel_count = 0
    !$ integer OMP_get_num_threads
    !$ integer OMP_get_thread_num
    !$ integer OMP_get_num_procs
    !
    ! Reset some non-essential counters
    !
    camera_subpixeling_npixfine = 0
    camera_subpixeling_npixtot  = 0
    !
    ! Here we decide whether to make a "normal" image or
    ! whether we find the tau=1 surface 
    !
    if(.not.dotausurf) then
       !
       ! Make a "normal" image, i.e. compute the intensity of all pixels
       !
       ! *** NEAR FUTURE: PUT OPENMP DIRECTIVES HERE (START) ***
       !
       !$ seconds = omp_get_wtime()
       !
       !$OMP PARALLEL &
       !
       !!$ Local variables from this function.
       !
       !$OMP PRIVATE(px,py,id,nthreads,pixel_count)
       !
       !$ pixel_count = 0
       !
       !$ id=OMP_get_thread_num()
       !$ nthreads=OMP_get_num_threads()
       !$ write(stdo,*) 'Thread Nr',id,'of',nthreads,'threads in total'
       flag_quv_too_big = .false.
       !
       !$OMP DO COLLAPSE(2) SCHEDULE(dynamic)
       !
       do iy=1,camera_image_ny
          do ix=1,camera_image_nx
             !
             ! Set the ray variables
             !
             px = camera_zoomcenter_x + ((ix-1)-(0.5d0*(camera_image_nx-1)))*pdx
             py = camera_zoomcenter_y + ((iy-1)-(0.5d0*(camera_image_ny-1)))*pdy
             !
             ! Now compute the intensity of this pixel
             !
             ! Note that if nrrefine>0 this routine will call itself
             ! recursively until a desired spatial resolution is acquired
             ! so as to guarantee that all flux is picked up, but this
             ! recursion is limited to nrrefine depth levels.
             !
             call camera_compute_one_pixel(camera_nrfreq,inu00,inu11,px,py,pdx,pdy,    &
                                           camera_nrrefine,camera_intensity_iquv,-1)
             !
             ! Put the result into the image array
             !
             camera_rect_image_iquv(ix,iy,inu00:inu11,1) = camera_intensity_iquv(inu00:inu11,1)
             if(camera_stokesvector) then
                !
                ! Copy also the other Stokes components
                !
                camera_rect_image_iquv(ix,iy,inu00:inu11,2) = camera_intensity_iquv(inu00:inu11,2)
                camera_rect_image_iquv(ix,iy,inu00:inu11,3) = camera_intensity_iquv(inu00:inu11,3)
                camera_rect_image_iquv(ix,iy,inu00:inu11,4) = camera_intensity_iquv(inu00:inu11,4)
                !
                ! Self-consistency check
                !
                do inuu=inu00,inu11
                   quvsq = camera_intensity_iquv(inuu,2)**2 + &
                           camera_intensity_iquv(inuu,3)**2 + &
                           camera_intensity_iquv(inuu,4)**2
                   if(quvsq.gt.0.d0) then
                      if(camera_intensity_iquv(inuu,1).eq.0.d0) then
                         write(stdo,*) 'INTERNAL ERROR: Q^2+U^2+V^2>0 but I=0...'
                         write(stdo,*) '    Warn author.'
                         stop
                      endif
                      quvsq = quvsq / camera_intensity_iquv(inuu,1)**2
                      if(quvsq.gt.1.000001d0) then
                         flag_quv_too_big = .true.
                      endif
                   endif
                enddo
             endif
             !
             !$ pixel_count = pixel_count + 1
          enddo
       enddo
       !
       !$OMP END DO
       !
       if(flag_quv_too_big) then
          write(stdo,*) 'WARNING: While making an image, I found an instance of Q^2+U^2+V^2>I^2...'
       endif
       !
       ! *** NEAR FUTURE: PUT OPENMP DIRECTIVES HERE (FINISH) ***
       !
       !$   write(stdo,*) 'Thread:',id,'raytraced:',pixel_count,'pixels'
       !
       !$OMP END PARALLEL
       !
       !$ write(stdo,*)"Elapsed time:",omp_get_wtime() - seconds;
       !
    else
       !
       ! Find the tau=1 surface (or any tau=tausurf surface)
       !
       ! Set the nrrefine to 0, and set the camera_tracemode to -3
       !
       backup_nrrefine  = camera_nrrefine
       backup_tracemode = camera_tracemode
       camera_nrrefine  = 0
       camera_tracemode = -3
       !
       warn_tausurf_problem = .false.
       !
       do iy=1,camera_image_ny
          do ix=1,camera_image_nx
             !
             ! Set the ray variables
             !
             px = camera_zoomcenter_x + ((ix-1)-(0.5d0*(camera_image_nx-1)))*pdx
             py = camera_zoomcenter_y + ((iy-1)-(0.5d0*(camera_image_ny-1)))*pdy
             !
             ! Reset the camera_xyzsstop
             !
             camera_taustop(inu00:inu11) = 1d91
             camera_xstop(inu00:inu11)   = -1d91
             camera_ystop(inu00:inu11)   = -1d91
             camera_zstop(inu00:inu11)   = -1d91
             camera_dstop(inu00:inu11)   = -1d91
             !
             ! Now do a first ray trace, to find the total optical depth
             !
             call camera_compute_one_pixel(camera_nrfreq,inu00,inu11,px,py,pdx,pdy,    &
                                           camera_nrrefine,camera_intensity_iquv,-1)
             !
             ! Check if anywhere the optical depth goes beyond 1d14*camera_tausurface
             !
             do inuu=inu00,inu11
                if(camera_intensity_iquv(inuu,1).gt.1d14*camera_tausurface) then
                   warn_tausurf_problem = .true.
                endif
             enddo
             !
             ! Compute the taustop
             !
             camera_taustop(inu00:inu11) = camera_intensity_iquv(inu00:inu11,1) - camera_tausurface
             !
             ! Reset xyzdstop (should not be necessary, but just for safety)
             !
             camera_xstop(inu00:inu11)   = -1d91
             camera_ystop(inu00:inu11)   = -1d91
             camera_zstop(inu00:inu11)   = -1d91
             camera_dstop(inu00:inu11)   = -1d91
             !
             ! Now do the second ray trace, to find the tau=tausurface surface
             !
             call camera_compute_one_pixel(camera_nrfreq,inu00,inu11,px,py,pdx,pdy,    &
                                           camera_nrrefine,camera_intensity_iquv,-1)
             !
             ! Put the resulting dstop into the image array
             !
             camera_rect_image_iquv(ix,iy,inu00:inu11,1) = camera_dstop(inu00:inu11)
             !
             ! If the big arrays for the xyz stop are available, then also store
             ! these
             !
             if(allocated(camera_tausurface_x)) then
                camera_tausurface_x(ix,iy,inu00:inu11) = camera_xstop(inu00:inu11)
                camera_tausurface_y(ix,iy,inu00:inu11) = camera_ystop(inu00:inu11)
                camera_tausurface_z(ix,iy,inu00:inu11) = camera_zstop(inu00:inu11)
             endif
             !
          enddo
       enddo
       !
       ! Print a warning if the tau exceeded the limit somewhere
       !
       if(warn_tausurf_problem) then
          write(stdo,*) 'WARNING: The optical depth exceeded 1d14*tausurface, so the tau surface determination may go wrong.'
       endif
       !
       ! Reset the original values of the tracemode and nrrefine
       !
       camera_nrrefine  = backup_nrrefine
       camera_tracemode = backup_tracemode
    endif
    !
    ! Add the discrete star sources. 
    !
    ! NOTE: Only for the mode in which stars are treated as point-sources.
    !
    if((camera_incl_stars.ne.0).and.(.not.star_sphere).and.(.not.dotausurf)) then
       !
       ! Do a check
       !
       if(.not.allocated(star_spec).and.(nstars.gt.0)) then
          write(stdo,*) 'WARNING in camera module: Stars not allocated.'
          stop
       endif
       !
       ! Now add all stars
       !
       do istar=1,nstars
          !
          ! Set the starting position and direction of the ray, as well as
          ! the location in the image
          !
          call camera_set_ray_stars_pntsrc(istar,x,y,z,dirx,diry,dirz,px,py,distance)
          !
          ! Check if this star is located in the image
          !
          if((px.ge.camera_zoomcenter_x-camera_image_halfsize_x).and.  &
             (px.le.camera_zoomcenter_x+camera_image_halfsize_x).and.  &
             (py.ge.camera_zoomcenter_y-camera_image_halfsize_y).and.  &
             (py.le.camera_zoomcenter_y+camera_image_halfsize_y)) then
             !
             ! Find the indices of the pixel in which the star is
             !
             idxx = floor((px-(camera_zoomcenter_x-camera_image_halfsize_x))/pdx)+1
             idxy = floor((py-(camera_zoomcenter_y-camera_image_halfsize_y))/pdy)+1
             if((idxx.lt.1).or.(idxx.gt.camera_image_nx).or.&
                (idxy.lt.1).or.(idxy.gt.camera_image_ny)) stop 8309
             !
             ! Set the intensity to the stellar spectrum
             ! 
             do inu=inu00,inu11
                camera_intensity_iquv(inu,1) =                               &
                     find_starlight_interpol(camera_frequencies(inu),istar)
                camera_intensity_iquv(inu,2) = 0.
                camera_intensity_iquv(inu,3) = 0.
                camera_intensity_iquv(inu,4) = 0.
             enddo
             !
             ! Do ray tracing
             !
             call camera_serial_raytrace(camera_nrfreq,inu00,inu11,x,y,z,dirx,diry,dirz,distance,  &
                                  celldxmin,camera_intensity_iquv)
             !
             ! Compute the ratio starsurface / pixelsurface
             !
             if(camera_localobserver) then
                factor = pi*(star_r(istar)/distance)**2 / (pdx*pdy)
             else
                factor = pi*star_r(istar)**2 / (pdx*pdy)
             endif
             !
             ! Check if the point source assumption is violated. 
             !
!             if(factor.gt.1.d0) then
!                write(stdo,*) 'ERROR in camera module: The image resolution is so'
!                write(stdo,*) '      high that the stellar surface is no longer'
!                write(stdo,*) '      small enough to be considered a point source.'
!                write(stdo,*) '  NOTE: In the (hopefully near) future a mode will'
!                write(stdo,*) '        be implemented to include non-point stars.'
!                stop
!             endif
             !
             ! Now modify the intensity of the pixel to include the star
             ! Note: We always assume that the starlight is unpolarized
             !
             camera_rect_image_iquv(idxx,idxy,inu00:inu11,1) =                        &
                   (1.d0-factor) * camera_rect_image_iquv(idxx,idxy,inu00:inu11,1) +  &  
                   factor * camera_intensity_iquv(inu00:inu11,1)
             if(camera_stokesvector) then
                camera_rect_image_iquv(idxx,idxy,inu00:inu11,2) =                     &
                   (1.d0-factor) * camera_rect_image_iquv(idxx,idxy,inu00:inu11,2)
                camera_rect_image_iquv(idxx,idxy,inu00:inu11,3) =                     &
                   (1.d0-factor) * camera_rect_image_iquv(idxx,idxy,inu00:inu11,3)
                camera_rect_image_iquv(idxx,idxy,inu00:inu11,4) =                     &
                   (1.d0-factor) * camera_rect_image_iquv(idxx,idxy,inu00:inu11,4)
             endif
          endif
       enddo
    endif
  end subroutine camera_make_rect_image_sub
  !
end subroutine camera_make_rect_image


!-------------------------------------------------------------------------
!                              WRITE THE IMAGE
!-------------------------------------------------------------------------
subroutine camera_write_image(img,ifill,noclip)
  implicit none
  character*80 :: filename,base,ext
  integer :: ix,iy,img,iinu,is,ns
  integer :: ifill
  logical, optional :: noclip
  logical :: donoclip
  double precision :: dummy(1:4)
  integer(kind=8) :: iformat, nn, kk
  !
  ! Interpret the noclip
  !
  donoclip = .false.
  if(present(noclip)) then
     donoclip = noclip
  endif
  !
  ! Open file
  !
  if(stdo.eq.6) then
     !
     ! Current radmc run is not a child process 
     ! So write output to a file  (otherwise to standard out)
     !
     if(img.le.0) then
        !
        ! A single image
        !
        if(writeimage_unformatted) then
           if(rto_style.eq.2) then
              open(unit=fflo,file='image.uout',form='unformatted')
           else
              open(unit=fflo,file='image.bout',status='replace',access='stream')
           endif
        else
           open(unit=fflo,file='image.out')
        endif
     else
        !
        ! One of a set of images
        !
        base='image_'
        if(writeimage_unformatted) then
           if(rto_style.eq.2) then
              ext ='.uout'
              if(ifill.le.0) then
                 call make_indexed_filename(base,img,ext,filename)
              else
                 call make_indexed_filename_fill(base,img,ext,filename,ifill)
              endif
              open(unit=fflo,file=filename,form='unformatted')
           else
              ext ='.bout'
              if(ifill.le.0) then
                 call make_indexed_filename(base,img,ext,filename)
              else
                 call make_indexed_filename_fill(base,img,ext,filename,ifill)
              endif
              open(unit=fflo,file=filename,status='replace',access='stream')
           endif
        else
           ext ='.out'
           if(ifill.le.0) then
              call make_indexed_filename(base,img,ext,filename)
           else
              call make_indexed_filename_fill(base,img,ext,filename,ifill)
           endif
           open(unit=fflo,file=filename)
        endif
     endif
  else
     !
     ! Current radmc run is a child process, so the image must
     ! be written to standard out. Check if fflo is indeed 6.
     !
     if(fflo.ne.6) then
        write(stdo,*) 'ERROR: fflo.ne.6 while stdo.ne.6'
        stop
     endif
  endif
  !
  ! Decide how many of the Stokes components we have to write
  !
  if(camera_stokesvector) then
     ns=4
  else
     ns=1
  endif
  !
  ! Now decide whether unformatted or formatted
  !
  if(writeimage_unformatted) then
     if (rto_style.eq.2) then
         !
         ! Unformatted output (F77-Unformatted output, faster, more compact)
         !
         if(camera_stokesvector) then
            if(camera_localobserver) then
               write(fflo) 4           ! Format number: size units in radian (angular) + Stokes
            else
               write(fflo) 3           ! Format number: size units in cm (spatial size) + Stokes
            endif
         else
            if(camera_localobserver) then
               write(fflo) 2           ! Format number: size units in radian (angular)
            else
               write(fflo) 1           ! Format number: size units in cm (spatial size)
            endif
         endif
         write(fflo) camera_image_nx,camera_image_ny
         write(fflo) camera_nrfreq
         write(fflo) 2.d0*camera_image_halfsize_x/(1.d0*camera_image_nx), &
              2.d0*camera_image_halfsize_y/(1.d0*camera_image_ny)
         write(fflo) (1d4*cc/camera_frequencies(iinu),iinu=1,camera_nrfreq)
         do iinu=1,camera_nrfreq
            do is=1,ns
               write(fflo) ((camera_rect_image_iquv(ix,iy,iinu,is),  &
                         ix=1,camera_image_nx),iy=1,camera_image_ny)
            enddo
         enddo
      else
         !
         ! DEFAULT: Unformatted output (C-compliant binary output, faster, more compact)
         !
         if(camera_stokesvector) then
            if(camera_localobserver) then
                iformat = 4           ! Format number: size units in radian (angular) + Stokes
            else
                iformat = 3           ! Format number: size units in cm (spatial size) + Stokes
            endif
         else
            if(camera_localobserver) then
                iformat = 2           ! Format number: size units in radian (angular)
            else
                iformat = 1           ! Format number: size units in cm (spatial size)
            endif
         endif
         write(fflo) iformat
         nn = camera_image_nx
         kk = camera_image_ny
         write(fflo) nn, kk
         nn = camera_nrfreq
         write(fflo) nn
         write(fflo) 2.d0*camera_image_halfsize_x/(1.d0*camera_image_nx), &
              2.d0*camera_image_halfsize_y/(1.d0*camera_image_ny)
         write(fflo) (1d4*cc/camera_frequencies(iinu),iinu=1,camera_nrfreq)
         do iinu=1,camera_nrfreq
            do is=1,ns
               write(fflo) ((camera_rect_image_iquv(ix,iy,iinu,is),  &
                         ix=1,camera_image_nx),iy=1,camera_image_ny)
            enddo
         enddo
      endif
  else
     !
     ! Standard: write it formatted way
     !
     if(camera_stokesvector) then
        if(camera_localobserver) then
           write(fflo,*) 4           ! Format number: size units in radian (angular)
        else
           write(fflo,*) 3           ! Format number: size units in cm (spatial size)
        endif
     else
        if(camera_localobserver) then
           write(fflo,*) 2           ! Format number: size units in radian (angular)
        else
           write(fflo,*) 1           ! Format number: size units in cm (spatial size)
        endif
     endif
     write(fflo,*) camera_image_nx,camera_image_ny
     write(fflo,*) camera_nrfreq
     write(fflo,*) 2.d0*camera_image_halfsize_x/(1.d0*camera_image_nx), &
          2.d0*camera_image_halfsize_y/(1.d0*camera_image_ny)
     write(fflo,320) (1d4*cc/camera_frequencies(iinu),iinu=1,camera_nrfreq)
320  format(E21.14)
     write(fflo,*) ' '
     do iinu=1,camera_nrfreq
        do iy=1,camera_image_ny
           do ix=1,camera_image_nx
              if(ns.eq.1) then
                 if((camera_rect_image_iquv(ix,iy,iinu,1).gt.1d-97).or.donoclip) then
                    write(fflo,333) camera_rect_image_iquv(ix,iy,iinu,1)
333                 format(E21.14)
                 else
                    write(fflo,*) 0.d0
                 endif
              else
                 dummy(1:4) = camera_rect_image_iquv(ix,iy,iinu,1:4)
                 do is=1,ns
                    if((abs(dummy(is)).le.1d-97).and.(.not.donoclip)) then
                       dummy(is) = 0.d0
                    endif
                 enddo
                 write(fflo,334) dummy(1:4)
334              format(4(E21.14,1X))
              endif
           enddo
        enddo
        write(fflo,*) ' '
     enddo
  endif
  if(stdo.eq.6) then
     !
     ! Output went to file, so close that file
     !
     close(1)
  else
     !
     ! Output went to standard output, so flush this
     !
     call flush(fflo)
  endif
  !
end subroutine camera_write_image


!-------------------------------------------------------------------------
!                  MAKE A SPECTRUM USING THE IMAGES
!
! This routine calls the make image routine at all frequencies and
! integrates these to make the spectrum. It can do it the simple
! way (just making rectangular images and integrating over them), or
! with aperture information (approximated as a circular mask with
! a wavelength-dependent size given in the aperture_info.inp file).
! Note, however, that the use of the aperture may make the routine
! slower in some cases.
!
! Be sure to have set the following variables to the appropriate values:
!  camera_incl_stars                Put to 1 for full SED with starlight
!  camera_pointing_position(1:3)    As with images
!  camera_image_halfsize_x,y        This is the size of the image 
!                                   in which the fluxes are computed for
!                                   the spectrum
!  camera_zoomcenter_x,y            As with images
!  camera_pointing_degr_*           As with images
!  camera_refine_criterion          As with images
!  camera_observer_degr_theta       As with images
!  camera_observer_degr_phi         As with images
!  camera_use_aperture_info         =.false. --> Integrates square image
!                                   =.true.  --> Uses aperture
!  
!-------------------------------------------------------------------------
subroutine camera_make_spectrum()
  use constants_module
  implicit none
  integer :: inu,backup_nrrefine,ierr,backup_nrfreq
  integer :: backup_lines_mode
  double precision :: pdx,pdy,factor,r_colarea,r,eps,hsx_bk,hsy_bk
  double precision, allocatable :: backup_freq(:),backup_spectrum_iquv(:,:)
  double precision, allocatable :: mask(:,:)
  integer :: ix,iy,apinu
  integer :: iact
  logical :: redo
  !
  ! Check
  !
  if(rt_incl_lines.and.(lines_maxshift.le.0.d0)) then
     write(stdo,*) 'INTERNAL ERROR in lines: lines_maxshift variable not set'
     stop
  endif
  if(igrid_coord.eq.10) then
     write(stdo,*) 'ERROR: In 1-D plane-parallel mode the spectrum mode is not'
     write(stdo,*) '       available, because for an infinitely large slab'
     write(stdo,*) '       it anyway does not have meaning. Instead you should'
     write(stdo,*) '       make multi-wavelength 1x1-pixel images to get the'
     write(stdo,*) '       frequency-dependent intensity.'
     stop
  endif
  if(igrid_coord.eq.20) then
     write(stdo,*) 'ERROR: In 2-D pencil-parallel mode the spectrum mode is not'
     write(stdo,*) '       available, because for an infinitely long slab'
     write(stdo,*) '       it anyway does not have meaning. Instead you should'
     write(stdo,*) '       make multi-wavelength 1xN-pixel images to get the'
     write(stdo,*) '       frequency-dependent intensity.'
     stop
  endif
  if((camera_nrfreq.le.0).or.(.not.allocated(camera_frequencies))) then
     write(stdo,*) 'ERROR: Somehow the frequencies for the spectrum are not set.'
     stop
  endif
  if(camera_localobserver) then
     write(stdo,*) 'ABORTING: For the moment we do not allow the making of a spectrum'
     write(stdo,*) '          in the local observer mode. Use observer at infinity mode,'
     write(stdo,*) '          for instance by specifying the incl and phi on the '
     write(stdo,*) '          command line and setting the image size with sizeau or sizepc.'
     stop
  endif
  if((((scattering_mode.ne.0).and.(.not.camera_scatsrc_allfreq)).or.   &
      camera_secondorder).and.(lines_mode.lt.-1).and.rt_incl_lines) then
     write(stdo,*)    '===> Warning: Combination of modes that may make RADMC-3D very slow... Here is why:'
     if((scattering_mode.ne.0).and.(.not.camera_scatsrc_allfreq)) then
        write(stdo,*) '     You are including dust scattering (doing it freq-by-freq), which means that'
     else
        if(camera_catch_doppler_line) then
           write(stdo,*) '     You are using the doppler catching mode, which requires second order '
           write(stdo,*) '     integration of the RT equation, which (for memory reasons) means that'
        else
           write(stdo,*) '     You are using second order integration of the RT equation,'
           write(stdo,*) '     which (for memory reasons) means that:'
        endif
     endif
     write(stdo,*) '     RADMC-3D must make one image for each frequency at a time to compute each '
     write(stdo,*) '     point of the spectrum. At the same time you use on-the-fly computation of'
     write(stdo,*) '     the line level populations, which are not stored. If you have 100 frequencies'
     write(stdo,*) '     in your spectrum, this means that 100 times an image is made, and a 100 times'
     write(stdo,*) '     the populations are RECALCULATED. If the population calculation is heavy, this'
     write(stdo,*) '     means that you waste a huge amount of time. You may want to calculate these'
     write(stdo,*) '     populations ONCE and store them, which can be done using a positive lines_mode'
     write(stdo,*) '     number (see manual). This may require a lot of memory if you have a molecule'
     write(stdo,*) '     with many levels. To overcome this, you can tell RADMC-3D to store only a subset'
     write(stdo,*) '     of the computed level populations: Only those related to the line you wish to'
     write(stdo,*) '     model. See manual.'
  endif
  !
  ! Make a backup of the camera_nrrefine, and set camera_nrrefine to infinity
  !
  backup_nrrefine = camera_nrrefine 
  camera_nrrefine = 10000
  !
  ! Set some other values
  !
  camera_localobserver = .false.      ! Spectra only for observer at infinity for now
  camera_tracemode = 1
  camera_nrrefine = 100
  !
  ! Check if the aperture file is present and read if needed
  !
  if(camera_use_aperture_info) then
     !
     ! Read the file aperture_info.inp
     ! 
     call camera_read_aperture_info()
     if(.not.allocated(camera_aperture_freq)) stop 8341
     if(.not.allocated(camera_aperture_radius_collectarea)) stop 8342
     !
     ! Check if the current wavelengths are all within range
     !
     if((max(camera_frequencies(1),camera_frequencies(camera_nrfreq)).gt.           &
         max(camera_aperture_freq(1),camera_aperture_freq(camera_aperture_nf))).or. &
        (min(camera_frequencies(1),camera_frequencies(camera_nrfreq)).lt.           &
         min(camera_aperture_freq(1),camera_aperture_freq(camera_aperture_nf)))) then
        write(stdo,*) 'ERROR while making spectrum with aperture information:'
        write(stdo,*) '      The aperture information does not cover all wavelengths'
        write(stdo,*) '      that are to be used for the spectrum.'
        stop
     endif
  endif
  !
  ! If the dust emission is included, then make sure the dust data,
  ! density and temperature are read. If yes, do not read again.
  !
  if(rt_incl_dust) then
     call read_dustdata(1)
     call read_dust_density(1)
     call read_dust_temperature(1)
     if(alignment_mode.ne.0) then
        call aligned_grains_init(1)
     endif
  endif
  !
  ! If line emission is included, then make sure the line data are
  ! read. If yes, then do not read it again.
  !
  if(rt_incl_lines) then
     call read_lines_all(1)
  endif
  !
  ! Now switch between different methods of making the spectrum
  !
  if(.not.camera_use_aperture_info) then
     !
     ! Simplest form of spectrum
     ! 
     ! It has turned out that under some circumstances this may
     ! lead to bad results because the sub-pixeling may accidently skip
     ! certain important parts of the image. Instead, I now decided
     ! to keep the camera_image_nx and camera_image_ny unchanged.
     ! 15.09.2016
     !
     ! Init or re-init the camera
     !
     call camera_init()
     !
     ! Set some values
     !
     pdx = 2*camera_image_halfsize_x / (1.d0*camera_image_nx)
     pdy = 2*camera_image_halfsize_y / (1.d0*camera_image_ny)
     !
     ! Make the 2x2 pixel images
     !
     call camera_make_rect_image(0)
     !
     ! Now make the spectrum at 1 parsec distance
     !
     factor = pdx*pdy/(parsec**2)
     do inu=1,camera_nrfreq
        camera_spectrum_iquv(inu,1:4) = 0.d0
        do iy=1,camera_image_ny
           do ix=1,camera_image_nx
              camera_spectrum_iquv(inu,1:4) = camera_spectrum_iquv(inu,1:4) + &
                   camera_rect_image_iquv(ix,iy,inu,1:4) * factor
           enddo
        enddo
     enddo
     !
     ! Restore the camera_nrrefine and the image nx,ny to their original value
     ! 
     camera_nrrefine = backup_nrrefine 
  else
     !
     ! Use the aperture information to make the spectrum
     !
     ! Use the nx,ny as in the images, but check that they are not too small.
     !
     if((camera_image_nx.lt.10).or.(camera_image_ny.lt.10)) then
        write(stdo,*) 'ERROR in camera module while making spectrum:'
        write(stdo,*) '      Taking npix less than 10 in either direction '
        write(stdo,*) '      is bound to yield wrong results.'
        write(stdo,*) '      nx,ny = ',camera_image_nx,camera_image_ny
        stop
     endif
     !
     ! Now make a separate image for each wavelength, where the size of the image is
     ! tuned to the size of the circular collecting area. 
     !
     ! We must first "fool" the imaging routine in thinking that we have only 1 frequency.
     ! So we must make a backup of the camera_frequencies() array and also allocate
     ! memory for the spectrum.
     !
     ! For 1 <= lines_mode <= 9 we must also tell the imaging routine that the
     ! level populations have already been calculated, and calculate them here.
     ! Otherwise the imaging routine will recalculate the populations every time
     ! again, which is very time-consuming for e.g. LVG or opt. thin. populations.
     !
     allocate(backup_freq(1:camera_nrfreq),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in spectrum: could not allocate backup_freq() array'
        stop
     endif
     allocate(backup_spectrum_iquv(1:camera_nrfreq,1:4),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in spectrum: could not allocate backup_spectrum_iquv() array'
        stop
     endif
     backup_spectrum_iquv(:,:) = -1d33     ! Purposely put to unphysical value, for check
     backup_freq(:) = camera_frequencies(:)
     backup_nrfreq  = camera_nrfreq
     backup_lines_mode = lines_mode
     !
     ! So if 1 <= lines_mode <= 9 we calculate the populations beforehand and
     ! put lines_mode to lines_mode=50, which means that the imaging routine
     ! will simply use the level populations that are present.
     !
     if(rt_incl_lines.and.(lines_mode.ge.1).and.(lines_mode.le.9)) then
        !
        ! Level subset selection and computation of the populations
        !
        if(lines_autosubset) then
           !
           ! If requested, do an automatic subset selection of the
           ! molecular levels
           ! 
           call lines_automatic_subset_selection(camera_nrfreq,     &
                camera_frequencies,1,camera_nrfreq,lines_maxshift,  &
                redo)
           !
           ! Now force a recomputation of the populations, unless
           ! "redo" is .false., meaning that we have the same
           ! levels as before (and thus the population as computed
           ! before must be still correct). Note that if there
           ! exists no lines_levelpop() array, then redo will also
           ! be .true.
           !
           if(redo) then
              iact=2
           else
              iact=1
           endif
           call lines_compute_and_store_local_populations(iact)
        else
           !
           ! Subset is not automatically selected. So either it has
           ! been selected manually or not at all (in which case the
           ! "subset" is the full set of levels).
           !
           write(stdo,*) 'Will store level populations for all levels or ',   &
                'for a manually (in lines.inp) selected subset of levels.'
           !
           ! Compute the populations if they have not yet been
           ! computed beforehand
           !
           call lines_compute_and_store_local_populations(1)
           !
        endif
        !
        ! Set lines_mode to 50 (a backup has been made, so it will be
        ! restored)
        !
        lines_mode = 50
     endif
     !
     ! From now on we will do things 1 frequency at a time
     !
     camera_nrfreq  = 1
     !
     ! Make a circular mask
     !
     allocate(mask(camera_image_nx,camera_image_ny),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in spectrum: could not allocate mask() array'
        stop
     endif
     do iy=1,camera_image_ny
        do ix=1,camera_image_nx
           r = sqrt( ((ix-0.5d0-0.5d0*camera_image_nx)/(0.5d0*camera_image_nx))**2 + &
                     ((iy-0.5d0-0.5d0*camera_image_ny)/(0.5d0*camera_image_ny))**2 )
           if(r.le.1.d0) then
              mask(ix,iy) = 1.d0
           else
              mask(ix,iy) = 0.d0
           endif
        enddo
     enddo
     !
     ! Re-init the camera, now making it think that we have just 1 frequency
     !
     call camera_init()
     !
     ! Perhaps not necessary, but let us backup the image half-sizes
     !
     hsx_bk = camera_image_halfsize_x
     hsy_bk = camera_image_halfsize_y
     !
     ! Now a loop over all frequencies
     !
     do inu=1,backup_nrfreq
        !
        ! Set frequency for this image
        !
        camera_frequencies(1) = backup_freq(inu)
        !
        ! Now interpolate the aperture size from the tables
        ! Use logarithmic interpolation, which preserves
        ! the powerlaw character of the beam size as a
        ! function of wavelength.
        !
        call hunt(camera_aperture_freq,camera_aperture_nf, &
                  camera_frequencies(1),apinu)
        if((apinu.lt.1).or.(apinu.ge.camera_aperture_nf)) then
           write(stdo,*) 'INTERNAL ERROR in spectrum: aperture index out of range...'
           stop 6398
        endif
        eps = (log(camera_frequencies(1))-log(camera_aperture_freq(apinu))) / &
              (log(camera_aperture_freq(apinu+1))-log(camera_aperture_freq(apinu)))
        if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 6399
        r_colarea = exp((1.d0-eps) * log(camera_aperture_radius_collectarea(apinu)) + &
                    eps * log(camera_aperture_radius_collectarea(apinu+1)))
        !
        ! Convert this radius from arcsec into cm for this image, using the
        ! distance specified by the observer
        !
        r_colarea = r_colarea * AU * camera_observer_distance_pc
        !
        ! Set the size of the image
        !
        camera_image_halfsize_x = r_colarea
        camera_image_halfsize_y = r_colarea
        !
        ! Set the pixel sizes
        !
        pdx = 2*r_colarea / (1.d0*camera_image_nx)
        pdy = 2*r_colarea / (1.d0*camera_image_ny)
        !
        ! Make the image
        !
        call camera_make_rect_image(0)
        !
        ! Now compute the flux at 1 parsec distance, using the mask
        !
        factor = pdx*pdy/(parsec**2)
        backup_spectrum_iquv(inu,1) = 0.d0
        do iy=1,camera_image_ny
           do ix=1,camera_image_nx
              if(camera_stokesvector) then
                 backup_spectrum_iquv(inu,1:4) = backup_spectrum_iquv(inu,1:4) + &
                      mask(ix,iy)*camera_rect_image_iquv(ix,iy,1,1:4) * factor
              else
                 backup_spectrum_iquv(inu,1) = backup_spectrum_iquv(inu,1) + &
                      mask(ix,iy)*camera_rect_image_iquv(ix,iy,1,1) * factor
              endif
           enddo
        enddo
     enddo
     !
     ! Reset the camera frequencies
     !
     camera_nrfreq = backup_nrfreq
     camera_frequencies(:) = backup_freq(:)
     !
     ! Restore lines_mode
     !
     lines_mode = backup_lines_mode
     !
     ! Re-init camera, to get the right nr of frequencies 
     ! in the arrays back again.
     !
     call camera_init()
     !
     ! Now reinstall the frequencies and copy the spectrum
     !
     camera_frequencies(:) = backup_freq(:)
     camera_spectrum_iquv(:,:) = backup_spectrum_iquv(:,:)
     !
     ! Perhaps not necessary, but let us restore the image half-sizes
     !
     camera_image_halfsize_x = hsx_bk
     camera_image_halfsize_y = hsy_bk
     !
  endif
  !
  ! Deallocate stuff
  !
  if(allocated(backup_freq)) deallocate(backup_freq)
  if(allocated(backup_spectrum_iquv)) deallocate(backup_spectrum_iquv)
  if(allocated(mask)) deallocate(mask)
  !
end subroutine camera_make_spectrum



!--------------------------------------------------------------------------
!                         WRITE THE SPECTRUM 
!--------------------------------------------------------------------------
subroutine camera_write_spectrum()
  use constants_module
  implicit none
  integer :: inu
  double precision :: lambda
  !
  ! Open file
  !
  if(stdo.eq.6) then
     !
     ! Current radmc run is not a child process 
     ! So write output to a file  (otherwise to standard out)
     !
     open(unit=fflo,file='spectrum.out')
  else
     !
     ! Current radmc run is a child process, so the spectrum must
     ! be written to standard out. Check if fflo is indeed 6.
     !
     if(fflo.ne.6) then
        write(stdo,*) 'ERROR: fflo.ne.6 while stdo.ne.6'
        stop
     endif
  endif
  !
  ! Now write
  !
  if(camera_stokesvector) then
     write(fflo,*) 3                             ! Format number
  else
     write(fflo,*) 1                             ! Format number
  endif
  write(fflo,*) camera_nrfreq
  write(fflo,*) ' '
  do inu=1,camera_nrfreq
     lambda = 1d4*cc/camera_frequencies(inu)
     if(camera_stokesvector) then
        if(camera_spectrum_iquv(inu,1).gt.1.d-97) then
           write(fflo,11) lambda,camera_spectrum_iquv(inu,1),&
                camera_spectrum_iquv(inu,2),&
                camera_spectrum_iquv(inu,3),&
                camera_spectrum_iquv(inu,4)
        else
           write(fflo,11) lambda,0.d0,0.d0,0.d0,0.d0
        endif
     else
        if(camera_spectrum_iquv(inu,1).gt.1.d-97) then
           write(fflo,10) lambda,camera_spectrum_iquv(inu,1)
        else
           write(fflo,10) lambda,0.d0
        endif
     endif
  enddo
10 format(E21.13,1X,E21.13)
11 format(E21.13,1X,E21.13,1X,E21.13,1X,E21.13,1X,E21.13)
  !
  ! Now close file or flush output
  !
  if(stdo.eq.6) then
     !
     ! Output went to file, so close that file
     !
     close(1)
  else
     !
     ! Output went to standard output, so flush this
     !
     call flush(fflo)
  endif
  !
end subroutine camera_write_spectrum



!--------------------------------------------------------------------------
!                      READ THE APERTURE INFORMATION
!--------------------------------------------------------------------------
subroutine camera_read_aperture_info()
  implicit none
  logical :: fex
  integer :: iformat,ierr,inu
  !
  ! Check that the distance is set
  !
  if(camera_observer_distance_pc.eq.0.d0) then
     write(stdo,*) 'ERROR: When using the aperture information, please'
     write(stdo,*) '       also specify the distance with dpc on the command line.'
     stop
  endif
  !
  ! First check if the file is present
  !
  inquire(file='aperture_info.inp',exist=fex)
  if(.not.fex) then
     write(stdo,*) 'ERROR: No file aperture_info.inp found.'
     stop
  endif
  !
  ! If some aperture stuff is already allocated, then deallocate
  !
  if(allocated(camera_aperture_freq)) deallocate(camera_aperture_freq)
  if(allocated(camera_aperture_radius_collectarea)) deallocate(camera_aperture_radius_collectarea)
  !
  ! Now open the file
  !
  open(unit=1,file='aperture_info.inp')
  read(1,*) iformat
  if(iformat.eq.1) then
     !
     ! With iformat==1 the radius of the collecting area as a function of
     ! wavelength is given.
     !
     read(1,*) camera_aperture_nf
     allocate(camera_aperture_freq(1:camera_aperture_nf),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in aperture information reading: could not allocate array'
        stop
     endif
     allocate(camera_aperture_radius_collectarea(1:camera_aperture_nf),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in aperture information reading: could not allocate array'
        stop
     endif
     do inu=1,camera_aperture_nf
        read(1,*) camera_aperture_freq(inu),camera_aperture_radius_collectarea(inu)
        camera_aperture_freq(inu) = 1d4*cc/camera_aperture_freq(inu) ! Convert to Hz
     enddo
     !
     ! Check that the frequencies are monotonic
     !
     if(camera_aperture_freq(camera_aperture_nf).gt.camera_aperture_freq(1)) then
        !
        ! Check monotonicity
        !
        do inu=2,camera_aperture_nf
           if(camera_aperture_freq(inu).le.camera_aperture_freq(inu-1)) then
              write(stdo,*) 'ERROR in aperture_info.inp file:'
              write(stdo,*) '      wavelengths not monotonic.'
              stop
           endif
        enddo
        !
        ! Just for safety, expand the range a tiny bit
        !
        camera_aperture_freq(1)                   = 0.999999d0 * camera_aperture_freq(1)
        camera_aperture_freq(camera_aperture_nf)  = 1.000001d0 * camera_aperture_freq(camera_aperture_nf)
     else
        !
        ! Check monotonicity
        !
        do inu=2,camera_aperture_nf
           if(camera_aperture_freq(inu).ge.camera_aperture_freq(inu-1)) then
              write(stdo,*) 'ERROR in aperture_info.inp file:'
              write(stdo,*) '      wavelengths not monotonic.'
              stop
           endif
        enddo
        !
        ! Just for safety, expand the range a tiny bit
        !
        camera_aperture_freq(1)                   = 1.000001d0 * camera_aperture_freq(1)
        camera_aperture_freq(camera_aperture_nf)  = 0.999999d0 * camera_aperture_freq(camera_aperture_nf)
     endif
  else
     write(stdo,*) 'ERROR: For now ',iformat,' is not known in aperture_info.inp'
     stop
  endif
  close(1)
  !
end subroutine camera_read_aperture_info



!----------------------------------------------------------------------------
!                  RAY TRACING IN GPU-PARALLELIZABLE WAY
!
! This routine is meant to emulate the way a GPU can parallelize the ray-
! tracing in frequency. 
!----------------------------------------------------------------------------
subroutine camera_ray1d_raytrace(nrfreq,x,y,z,dx,dy,dz,distance,celldxmin,&
                                   intensity)
  implicit none
  double precision :: x,y,z,dx,dy,dz,dum,distance,celldxmin
  double precision :: freq,alpnu,jnu,alpha_a,dummy,expt,temp,src
  integer :: ispec,nrfreq,deepestlevel
  integer :: inu,levelnext,is
  logical :: arrived
  double precision :: intensity(nrfreq)
  double precision :: rlen,cosp,sinp,cost,sint,vx,vy,vz
  double precision :: xav,yav,zav
  type(amr_branch), pointer :: a
  integer :: ixx,iyy,izz,bc_idir,bc_ilr
  !
  ! Do checks
  !
  write(stdo,*) 'ERROR: GPU emulator not yet ready'
  stop
  if(nrfreq.ne.camera_nrfreq) stop 3345
  if(camera_tracemode.ne.1) then
     write(stdo,*) 'ERROR in camera module: for lines only tracemode=1 allowed.'
     stop
  endif
  if(debug_check_all.eq.1) then
     if(rt_incl_lines) then
        if(.not.allocated(gastemp)) then
           write(stdo,*) 'ERROR in line ray tracing: gas temperature not set.'
           stop
        endif
        if(.not.allocated(lines_ray_levpop)) then
           write(stdo,*) 'ERROR in camera_ray1d_raytrace: lines module not ready:'
           write(stdo,*) '  lines_ray_levpop not allocated.'
           stop
        endif
        if(lines_mode.ge.1) then
           if(.not.allocated(lines_levelpop)) then
              write(stdo,*) 'ERROR in camera_ray1d_raytrace: lines module not ready:'
              write(stdo,*) '  lines_levelpop not allocated.'
              stop
           endif
        endif
     endif
     dum = sqrt(dx**2 + dy**2 + dz**2)
     if(abs(dum-1.d0).gt.1d-6) stop 3487
     if(.not.allocated(camera_frequencies)) then
        write(stdo,*) 'ERROR: Wanting to ray-trace using general method, but'
        write(stdo,*) '       camera_frequencies is not allocated.'
        stop
     endif
  endif
  !
  ! Set starting value of celldxmin, i.e. the size of the smallest cell
  ! we will pass through. Note that for cartesian coordinates, since cell 
  ! sizes will then always be cubes of the base size times (1/2)^level
  ! we can use just the deepest level to compute celldxmin.
  !
  ! NOTE: In the AMR module we now changed the level-numbering: Now
  !       the base grid level is level=0 (used to be 1). Hence we must
  !       now set deepestlevel = 0 here, and the multiplication factor
  !       is now (1/2)^level instead of (1/2)^(level-1).
  !
  celldxmin     = 1d99
  deepestlevel  = 0
  !
  ! Set the starting point and the direction
  !
  ray_cart_x    = x
  ray_cart_y    = y
  ray_cart_z    = z
  ray_cart_dirx = dx
  ray_cart_diry = dy
  ray_cart_dirz = dz
  !
  ! Put the distance from the starting point to the observer here
  ! If observer at infinity, simply put distance=1d99
  !
  ray_dsend     = distance 
  !
  ! Check if we start inside or outside the domain
  !
  if(igrid_type.lt.100) then 
     !
     ! We use a regular or AMR grid
     !
     if(igrid_coord.lt.100) then
        !
        ! We use cartesian coordinates
        !
        if(amr_tree_present) then
           !
           ! We use the AMR tree
           !
           call amr_findcell(x,y,z,a)
           if(associated(a)) then
              !
              ! We start inside the domain
              !
              ray_index     = a%leafindex
              ray_indexnext = 0
           else
              !
              ! We start outside the domain, so all indices must be 0
              !
              ray_index     = 0
              ray_indexnext = 0
           endif
        else
           !
           ! We have a regular grid
           !
           call amr_findbasecell(x,y,z,ixx,iyy,izz)
           if(ixx.gt.0) then
              ray_index     = ixx+((iyy-1)+(izz-1)*amr_grid_ny)*amr_grid_nx
              ray_indexnext = 0
           else
              ray_index     = 0
              ray_indexnext = 0
           endif
        endif
        !
     elseif(igrid_coord.lt.200) then
        !
        ! We use spherical coordinates
        !
        write(stdo,*) 'FINDING CELL NOT YET READY FOR SPHERICAL COORDINATES'
        stop
        !
     else
        write(stdo,*) 'ERROR: Cylindrical coordinates not yet implemented'
        stop
     endif
  else
     write(stdo,*) 'SORRY: Delaunay or Voronoi grids not yet implemented'
     stop
  endif
  !
  ! Reset arrays for the 1-D line tracing
  !
  camera_ray_jnu(:,:)   = 0.d0
  camera_ray_alpnu(:,:) = 0.d0
  ray_ns                = 0
  is = 1
  !
  ! Now trace
  !
  arrived = .false.
  do while(.not.arrived)
     !
     ! Check if the ray is not too long
     !
     if(is.gt.ray_nsmax) then
        write(stdo,*) 'ERROR in camera/lines: ray length exceeds maximum.'
        stop
     endif
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
           ! Check if this cell is the smallest so far
           !
           deepestlevel = max(deepestlevel,levelnext)
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
                    do inu=1,camera_nrfreq
                       intensity(inu) = bplanck(thermal_bc_temp(bc_ilr,bc_idir), &
                                                camera_frequencies(inu))
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
           !
           ! Check if this cell is the smallest so far
           !
           ! **** NOT YET READY. Use the camera_amr_zoomfact array
           ! **** multiplied by the local radius to get the correct
           ! **** celldxmin
           write(stdo,*) 'FINDING CELLDXMIN NOT YET READY FOR SPHERICAL COORDINATES'
           stop
           !
        else
           write(stdo,*) 'ERROR: Cylindrical coordinates not yet implemented'
           stop
        endif
        !
        ! Set the camera_istar, which is used for treating finite-size
        ! stars (instead of point source stars)
        !
        camera_istar = amrray_ispherehit
        !   
     else
        write(stdo,*) 'SORRY: Delaunay or Voronoi grids not yet implemented'
        stop
     endif
     !
     ! Check if we are in a cell or not
     !
     if(debug_check_all.eq.1) then
        if(ray_index.gt.nrcellsmax) stop 6271
     endif
     if(ray_index.ge.1) then
        !
        ! Set the length 
        !
        ray_ds            = sqrt( (ray_cart_x-ray_prev_x)**2 + (ray_cart_y-ray_prev_y)**2 + (ray_cart_z-ray_prev_z)**2 )
        camera_ray_ds(is) = ray_ds
        !
        ! Set the dust quantities
        !
        if(rt_incl_dust) then
           if(allocated(dustdens)) then
              sources_dustdens(:) = dustdens(:,ray_index)
              sources_dusttemp(:) = dusttemp(:,ray_index)
           else
              sources_dustdens(:) = 0.d0
              sources_dusttemp(:) = 0.d0
           endif
        endif
        !
        ! Find the local dust j_nu and alpha_nu
        !
        !  ############# FUTURE: ADD SCATTERING ALPHA AND SOURCE AS WELL ###########
        !
        if(rt_incl_dust) then
           temp    = 100.d0     ! Arbitrary temperature...
           do inu=1,camera_nrfreq
              freq  = camera_frequencies(inu)
              alpnu = 0.d0
              jnu   = 0.d0
              do ispec=1,dust_nr_species
                 alpha_a = sources_dustdens(ispec) * find_dust_kappa_interpol(freq,ispec,temp,1,0,0)
                 !      alpha_s = sources_dustdens(ispec) * find_dust_kappa_interpol(freq,ispec,temp,0,1,0)
                 alpnu   = alpnu + alpha_a ! + alpha_s
                 jnu     = jnu + alpha_a * bplanck(sources_dusttemp(ispec),camera_frequencies(inu))
              enddo
              camera_ray_jnu(inu,is)   = jnu
              camera_ray_alpnu(inu,is) = alpnu
           enddo
        else
           camera_ray_jnu(:,is)     = 0.d0
           camera_ray_alpnu(:,is) = 0.d0
        endif
        !
        ! Find the local line j_nu and alpha_nu
        !
        if(rt_incl_lines) then
           !
           ! Set the doppler shift
           !
           if(igrid_coord.lt.100) then
              !
              ! Cartesian coordinates:
              !
              dummy = gasvelocity(1,ray_index) * ray_cart_dirx +   &
                      gasvelocity(2,ray_index) * ray_cart_diry +   &
                      gasvelocity(3,ray_index) * ray_cart_dirz
           elseif(igrid_coord.lt.200) then
              !
              ! Spherical coordinates: Here component 1 is the v_r component,
              ! component 2 is the v_theta component (in the same unit as
              ! v_r, i.e. cm/s) and component 3 is the v_phi component (also
              ! in the same unit: cm/s).
              !
              xav   = 0.5d0 * ( ray_cart_x + ray_prev_x )
              yav   = 0.5d0 * ( ray_cart_y + ray_prev_y )
              zav   = 0.5d0 * ( ray_cart_z + ray_prev_z )
              dummy = xav**2+yav**2
              rlen  = sqrt(dummy+zav**2)+1d-90
              dummy = sqrt(dummy)+1d-90
              cosp  = xav/dummy
              sinp  = yav/dummy
              cost  = zav/rlen
              sint  = dummy/rlen
              dummy = sint * gasvelocity(1,ray_index) +   &
                      cost * gasvelocity(2,ray_index) 
              vx    = cosp * dummy - sinp * gasvelocity(3,ray_index)
              vy    = sinp * dummy + cosp * gasvelocity(3,ray_index)
              vz    = cost * gasvelocity(1,ray_index) -   &
                      sint * gasvelocity(2,ray_index)
              dummy = vx * ray_cart_dirx +   &
                      vy * ray_cart_diry +   &
                      vz * ray_cart_dirz
           else
              stop 501
           endif
           lines_ray_doppler(is) = dummy / cc
           !
           ! Set the gas temperature
           !
           lines_ray_temp(is)    = gastemp(ray_index)
           !
           ! Set the microturbulence
           !
           if(allocated(lines_microturb)) then
              lines_ray_turb(is) = lines_microturb(ray_index)
           else
              lines_ray_turb(is) = 0.d0
           endif
           !
           ! Set the level populations
           !
           ! NOTE: THE STUFF BELOW IS STILL NOT UP TO DATE
           !
           if(lines_mode.ge.1) then
              !
              ! We use the level populations from the global array
              ! **** THE LINE BELOW IS WRONG ****
              lines_ray_levpop(:,:,is) = lines_levelpop(:,:,ray_index)
           else
              !
              ! We compute the level populations on-the-fly
              !
              do ispec=1,lines_nr_species
                 call lines_compute_ltepop_subset(lines_nrlevels(ispec),    &
                        active_nrlevels(ispec),                             &
                        active_levels(:,ispec),                             &
                        ispec,lines_ray_temp(is),                           &
                        gas_chemspec_numberdens(ispec,ray_index),           &
                        lines_ray_levpop(:,ispec,is))
              enddo
           endif
           !
        endif
     else
        !
        ! Reset all to zero
        !
        ray_ds                = 0.d0
        camera_ray_ds(is)     = 0.d0
        if(rt_incl_lines) then
           lines_ray_doppler(is) = 0.d0
           lines_ray_temp(is)    = 0.d0
           lines_ray_turb(is)    = 0.d0
           lines_ray_levpop(:,:,is) = 0.d0
        endif
        !
     endif
     !
     ! Increase is
     !
     is = is + 1
     !
     ! Decrease ray_dsend
     !
     ray_dsend = ray_dsend - ray_ds
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
  ! Set the length of the 1-D ray
  !
  ray_ns = is - 1
  !
  ! Return the smallest cell size we encountered
  !
  if(igrid_type.lt.100) then 
     if(igrid_coord.lt.100) then
        celldxmin = (amr_grid_xi(2,1)-amr_grid_xi(1,1)) / (2**deepestlevel)
     elseif(igrid_coord.lt.200) then
        write(stdo,*) 'For spherical coordinates, the celldxmin not yet ready'
        stop
     else
        write(stdo,*) 'ERROR: Cylindrical coordinates not yet implemented'
        stop
     endif
  else
     write(stdo,*) 'SORRY: Delaunay or Voronoi grids not yet implemented'
     stop
  endif
  !
  ! Now install the line sources in the 1-D ray
  !
  if(rt_incl_lines) then
!!!     call lines_addto_jnu_alpnu()
     stop 3333
  endif
  !
  ! Now call the ray tracer of the lines module
  !
  intensity(:) = 0.d0
  !
  ! Do a loop over frequencies
  !
  ! ############## FUTURE: FOR NVIDIA GRAPHICS CARD: PARALLELIZE IN NU HERE #########
  !
  do inu=1,camera_nrfreq
     !
     ! Loop along the ray
     !
     do is=1,ray_ns
        !
        ! Each segment is assumed to be at constant j_nu and alpha_nu.
        !
        ! ############ FUTURE: HERE WE MUST CATCH THIN-LINE DOPPLER 'QUANTUM LEAPS' ############
        !
        expt   = exp(-camera_ray_ds(is)*camera_ray_alpnu(inu,is))
        src    = camera_ray_jnu(inu,is)/camera_ray_alpnu(inu,is)
        intensity(inu) = intensity(inu) * expt + (1.d0-expt) * src
     enddo
  enddo
  !
  ! Done
  !
end subroutine camera_ray1d_raytrace



!----------------------------------------------------------------------------
!                SET THE FREQUENCIES FOR THE CAMERA
!----------------------------------------------------------------------------
subroutine set_camera_frequencies(nfr,freq1,freq2)
  implicit none
  logical :: flag
  integer :: inu,ierr,i,iformat
  integer, intent(in), optional :: nfr
  double precision, intent(in), optional :: freq1,freq2
  logical :: fex1,fex2
  !
  ! Do some checks
  !
  if((lines_user_widthkms.gt.0.d0).or.(lines_user_iline.gt.0)) then
     !
     ! Apparently the user wants to make a spectral line spectrum
     !
     ! Check that at least both widthkms and iline are set
     !
     if((lines_user_widthkms.le.0).and.(lines_user_nrfreq.gt.1)) then
        write(stdo,*) 'ERROR: If you specify iline, you must also specify widthkms'
        write(stdo,*) '       (unless you make a single wavelength image).'
        write(stdo,*) '       Currently nr of freq = ',lines_user_nrfreq
        stop
     endif
     if(lines_user_iline.le.0) then
        write(stdo,*) 'ERROR: If you specify widthkms, you must also specify iline'
        stop
     endif
     !
     ! Some things can be defaults, if they are not yet set
     ! 
     if(lines_user_ispec.le.0) then
        lines_user_ispec = 1
        if(lines_nr_species.gt.1) then
           write(stdo,*) 'For line spectra we, by default, choose molecular species 1.'
           write(stdo,*) '   Specify with imolspec on the command line if you want another molecule.'
        endif
     endif
     if(lines_user_nrfreq.le.0) then
        lines_user_nrfreq = 100
        write(stdo,*) 'For line spectra we, by default, take nr of frequencies to 100.'
        write(stdo,*) '   Specify with linenlam or linenfreq if you want another nr.'
     endif
  endif
  !
  ! Deallocate frequencies
  !
  if(allocated(camera_frequencies)) deallocate(camera_frequencies)
  !
  ! Check whether to use the above nfr, freq1 and freq2 to make the
  ! frequencies
  !
  if(present(nfr)) then
     !
     ! A regularly spaced grid between freq1 and freq2
     !
     if((.not.present(freq1)).or.(.not.present(freq2))) then
        write(stdo,*) 'ERROR in set_camera_frequencies: either set all or no optional args.'
        stop
     endif
     camera_nrfreq = nfr
     allocate(camera_frequencies(camera_nrfreq),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in lines module: Could not allocate '
        write(stdo,*) '      camera_frequencies(:).'
        stop 
     endif
     if(camera_nrfreq.gt.1) then
        do inu=1,camera_nrfreq
           camera_frequencies(inu) = freq1 * (freq2/freq1)**   &
                ((inu-1.d0)/(camera_nrfreq-1.d0))
        enddo
     else
        camera_frequencies(1) = freq1
     endif
  else
     !
     ! Now make decisions where to get our frequency grid from
     !
     if((camera_lambdamic.gt.0.d0).and.(camera_lambdamic1.gt.0.d0).and. &
          (camera_range_nlam.gt.1)) then
        !
        ! Make a simple regular grid between lambdamic and lambdamic1
        ! using logarithmic spacing
        !
        camera_nrfreq = camera_range_nlam
        allocate(camera_frequencies(camera_nrfreq),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in camera module: Could not allocate '
           write(stdo,*) '      camera_frequencies(:).'
           stop 
        endif
        do inu=1,camera_nrfreq
           camera_frequencies(inu) = 1d4*cc / ( camera_lambdamic *   &
                (camera_lambdamic1/camera_lambdamic)**               &
                ((inu-1.d0)/(camera_nrfreq-1.d0)) )
        enddo
        !
     elseif(camera_lambdamic.gt.0.d0) then
        !
        ! Use exactly this wavelength, given on the command line
        !
        camera_nrfreq = 1
        allocate(camera_frequencies(camera_nrfreq),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in camera module: Could not allocate '
           write(stdo,*) '      camera_frequencies(:).'
           stop 
        endif
        camera_frequencies(1) = 1d4*cc/camera_lambdamic
        !
     elseif(camera_theinu.gt.0) then
        !
        ! Use one of the global frequency array frequencies, selected
        ! on the command line
        !
        if(.not.allocated(freq_nu)) stop
        camera_nrfreq = 1
        allocate(camera_frequencies(camera_nrfreq),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in camera module: Could not allocate '
           write(stdo,*) '      camera_frequencies(:).'
           stop 
        endif
        camera_frequencies(1) = freq_nu(camera_theinu)
        !
     elseif((lines_user_nrfreq.ge.1).and. &
        (lines_user_ispec.gt.0).and.(lines_user_iline.gt.0).and.        &
        rt_incl_lines) then
        !
        ! Use the line user-parameters to set the frequency grid 
        !
        ! Make sure the line data is read, otherwise we cannot specify
        ! the frequencies based on the line positions. 
        !
        call read_lines_all(1)
        !
        ! Check stuff
        !
        if(.not.allocated(lines_nu0)) then
           write(stdo,*) 'ERROR in camera module: while creating line frequencies array'
           write(stdo,*) 'line module not yet initialized.'
           stop
        endif
        if(lines_user_ispec.gt.lines_nr_species) then
           write(stdo,*) 'ERROR in line module: lines_user_ispec out of range'
           stop
        endif
        if(lines_user_iline.gt.lines_nrlines(lines_user_ispec)) then
           write(stdo,*) 'ERROR in line module: lines_user_iline out of range'
           stop
        endif
        !
        ! Allocate the camera frequencies array
        !
        camera_nrfreq = lines_user_nrfreq
        allocate(camera_frequencies(camera_nrfreq),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR in lines module: Could not allocate '
           write(stdo,*) '      camera_frequencies(:).'
           stop 
        endif
        !
        ! Make a window around a particular line
        !
        ! * CONVENTION * : Positive velocity = reddening (moving away from obs)
        !
        if(camera_nrfreq.gt.1) then
           do inu=1,camera_nrfreq
              camera_frequencies(inu) = lines_nu0(lines_user_iline,lines_user_ispec) *  &
                   (1.d0 - lines_user_kms0*1d5/cc -                                     &
                    (2*(inu-1.d0)/(lines_user_nrfreq-1.d0)-1.d0)*                       &
                    lines_user_widthkms*1d5/cc)
           enddo
        else
           camera_frequencies(1) = lines_nu0(lines_user_iline,lines_user_ispec) * &
                (1.d0-lines_user_kms0*1d5/cc)
        endif
     elseif(camera_loadlambdas) then
        !
        ! Read frequencies / wavelengths from a file
        !
        inquire(file='camera_frequency.inp',exist=fex1)
        inquire(file='camera_wavelength_micron.inp',exist=fex2)
        if(fex1.and.fex2) then
           write(stdo,*) 'ERROR in camera module: found camera_frequency.inp AND camera_wavelength_micron.inp'
           stop
        endif
        if(fex1) then
           open(unit=1,file='camera_frequency.inp')
           read(1,*) camera_nrfreq
           allocate(camera_frequencies(camera_nrfreq),STAT=ierr)
           if(ierr.ne.0) then
              write(stdo,*) 'ERROR in camera module: Could not allocate '
              write(stdo,*) '      camera_frequencies(:).'
              stop 
           endif
           do inu=1,camera_nrfreq
              read(1,*) camera_frequencies(inu)
           enddo
           close(1)
        elseif(fex2) then
           open(unit=1,file='camera_wavelength_micron.inp')
           read(1,*) camera_nrfreq
           allocate(camera_frequencies(camera_nrfreq),STAT=ierr)
           if(ierr.ne.0) then
              write(stdo,*) 'ERROR in camera module: Could not allocate '
              write(stdo,*) '      camera_frequencies(:).'
              stop 
           endif
           do inu=1,camera_nrfreq
              read(1,*) camera_frequencies(inu)
              camera_frequencies(inu) = 1d4*cc/camera_frequencies(inu)
           enddo
           close(1)
        else
           write(stdo,*) 'ERROR: Cannot find camera_wavelength_micron.inp or'
           write(stdo,*) '       camera_frequency.inp file.'
           stop
        endif
     elseif(camera_loadcolor) then
        !
        ! Use a subset of wavelengths from the global frequency array
        !
        if(.not.allocated(freq_nu)) stop
        if(rt_incl_lines) then
           write(stdo,*) 'ERROR: Cannot use SED colors and do lines at the same time'
           stop
        endif
        if(camera_theinu.ge.1) then
           write(stdo,*) 'ERROR: Cannot use SED colors and specify a single frequency at the same time'
           stop
        endif
        inquire(file='color_inus.inp',exist=fex1)
        if(.not.fex1) then
           write(stdo,*) 'ERROR: Could not find file color_inus.inp'
           stop
        endif
        write(stdo,*) 'Reading color_inus.inp file...'
        call flush(stdo)
        open(unit=1,file='color_inus.inp')
        read(1,*) iformat
        read(1,*) camera_nrfreq
        if(camera_nrfreq.le.0) then
           write(stdo,*) 'ERROR: File color_inus.inp corrupt.'
           stop
        endif
        allocate(camera_frequencies(1:camera_nrfreq),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate camera_frequencies array'
           stop
        endif
        do i=1,camera_nrfreq
           read(1,*) inu
           if((inu.lt.1).or.(inu.gt.freq_nr)) then
              write(stdo,*) 'ERROR: in color_inus.inp inu out of range.'
              stop
           endif
           camera_frequencies(i) = freq_nu(inu)
        enddo
        close(1)
     elseif(camera_setfreq_global) then
        !
        ! Use the global frequency array
        !
        if(.not.allocated(freq_nu)) stop
        camera_nrfreq = freq_nr
        allocate(camera_frequencies(1:camera_nrfreq),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate camera_frequencies array'
           stop
        endif
        do inu=1,camera_nrfreq
           camera_frequencies(inu) = freq_nu(inu)
        enddo
     else
        write(stdo,*) 'ERROR: No method specified for determining the wavelength grid'
        write(stdo,*) '       of the camera module. Use, on the command line, one of the following:'
        write(stdo,*) '          inu or ilambda followed by an integer'
        write(stdo,*) '          nuhz or lambda followed by a number'
        write(stdo,*) '          color or loadcolor to let RADMC-3D read color_inus.inp'
        write(stdo,*) '          loadlambda to let RADMC-3D read camera_wavelength_micron.inp'
        write(stdo,*) '          globalfreq or globallam to use the global wavelength array.'
        stop
     endif
  endif
end subroutine set_camera_frequencies


!==========================================================================
!             SOME BACKWARD COMPATIBILITY ROUTINES WITH RADMC
!==========================================================================

!--------------------------------------------------------------------------
!                  SET UP PIXEL POSITIONS FOR CIRCULAR IMAGES
!
! Since we now have better methods than the circular images, we do not
! really need these. But since RADMC uses them, and since we want to 
! be able to reproduce old RADMC results with RADMC-3D, it might be useful.
! Note that it only works for spherical coordinates.
!--------------------------------------------------------------------------
subroutine setup_pixels_circular(irmin,irmax,nrphiinf,nrext,dbdr,imethod,nrref)
  use amr_module
  implicit none
  integer :: irmin,irmax,nrphiinf,nrext,dbdr
  integer :: imethod,nrref
  integer :: ix,ir,iins,i,iphi
  integer :: iradius,nrrextra
  doubleprecision :: epsrrr,refdum1,refdum2,telesc_eps,dr
  parameter(telesc_eps=1d-10)
  parameter(epsrrr=1.0d3*telesc_eps)
  !
  ! Check if we are indeed in spherical coordinates
  !
  if(amrray_icoord.ne.1) then 
     write(stdo,*) 'ERROR: You can only use the circular images if you'
     write(stdo,*) '       use spherical coordinates.'
     stop
  endif
  if(.not.allocated(star_r)) then
     write(stdo,*) 'ERROR in circular image generator: For circular images at '
     write(stdo,*) '      a central star is required (even if its luminosity is 0).'
     stop
  endif
  !
  ! Check or fix the central star
  !
  if(nstars.ge.1) then
     !
     ! Set the position of the central star 
     !
     if((abs(star_pos(1,1)).gt.1.d-6*amr_grid_xi(1,1)).or.&
        (abs(star_pos(2,1)).gt.1.d-6*amr_grid_xi(1,1)).or.&
        (abs(star_pos(3,1)).gt.1.d-6*amr_grid_xi(1,1))) then
        write(stdo,*) 'Note: Setting position of first star to (0,0,0).'
     endif
     star_pos(:,1) = 0.d0
  endif
  !
  ! Default
  !
  nrrextra = nrext
  if(nrrextra.lt.0) nrrextra = -nrrextra
  !
  ! Decide on which way the radial rays are arranged
  !     
  if(imethod.eq.0) then
     !
     ! The original method from RADICAL
     !
     if(nrext.gt.0) then
        nrrextra = nrext
        imethod = -1
     elseif(nrext.lt.0) then
        nrrextra = -nrext
        imethod = -2
     else
        write(*,*) 'ERROR: Must have non-zero nrrextra'
        stop
     endif
  endif
  !
  ! Allocate arrays
  !
  cim_nr  = (irmax-irmin+1)*dbdr+nrrextra+nrref
  cim_np  = nrphiinf
  if(allocated(cim_rc)) deallocate(cim_rc)
  if(allocated(cim_ri)) deallocate(cim_ri)
  if(allocated(cim_pc)) deallocate(cim_pc)
  if(allocated(cim_pi)) deallocate(cim_pi)
  if(allocated(cim_surfarea)) deallocate(cim_surfarea)
  allocate(cim_rc(0:cim_nr),cim_ri(0:cim_nr+1))  ! Note: cell 0 = star
  allocate(cim_pc(1:cim_np),cim_pi(1:cim_np+1))
  allocate(cim_surfarea(0:cim_nr,1:cim_np))
  !
  ! Do some checks
  !
  if(irmax.le.irmin) then
     write(stdo,*) 'ERROR Telescope: irmax.le.irmin'
     stop
  endif
  if((nrrextra.gt.0).and.(irmin.ne.1)) then
     write(stdo,*) 'ERROR: In setup_pixels_circular():'
     write(stdo,*) '       In You request extra rays at the ', &
          'inner edge while you dont have irmin=1.'
     write(stdo,*) '       This is contradictory. Aborting...'
     stop
  endif
  if(nstars.ne.1) then
     write(stdo,*) 'ERROR: In circular images: must have nstar==1.'
     stop
  endif
  !
  ! Now let us make the set of rays for the image at
  ! infinity... This consists of a series of concentric rings
  ! of pixels around the center of the object. Each pixel of 
  ! the image corresponds to one ray defined here. This
  ! circular image has coordinates r (in [cm], is the radius
  ! of each ringin the circular image) and phi (in [rad], is
  ! the angle along each concentric ring). 
  !
  ! First make the phi grid, which is the easiest
  !
  do i=1,cim_np
     cim_pc(i) = twopi * (i-0.5d0)/cim_np
  enddo
  do i=1,cim_np+1
     cim_pi(i) = twopi * (i-1.d0)/cim_np
  enddo
  !
  ! Now make the r-grid, which must be adapted to the r-grid of the
  ! model itself. NOTE: It will only adapt itself to the REGULAR radial
  ! grid, not to the AMR grid. 
  !
  ! The 0th pixel is the star
  !
  ir              = 1
  iradius         = 0
  cim_ri(0)       = 0.d0
  cim_rc(0)       = 0.d0
  if(nstars.ge.1) then
     cim_ri(1)    = star_r(1)
  else
     cim_ri(1)    = 0.d0
  endif
  !     
  ! Next make the extra rays, which are the rays with impact parametes
  ! smaller than the radius of the inner edge of the (R,Theta) spatial
  ! grid of the object itself. 
  !
  ir              = 2
  iradius         = 1
  !
  ! These are the original RADICAL-style ways to arrange these
  ! rays. But be careful that there are drawbacks for near-face-on
  ! spectra for very optically thick disks with vertical inner walls.
  !
  if(imethod.eq.-2) then 
     !
     ! The old RADICAL way
     !      
     do ix=1,nrrextra 
        if(iradius.gt.cim_nr) then
           write(stdo,*) 'EROR in circular images: exceeded cim_nr'
           stop
        endif
        if(amr_grid_xi(1,1).le.cim_ri(1)) then
           write(stdo,*) 'ERROR in circular images: Grid inner edge is not'
           write(stdo,*) '      larger than stellar radius.'
           stop
        endif
        cim_rc(iradius) = ( ix * ( amr_grid_xi(1,1) - &
                          cim_ri(1) ) / (nrrextra+1.d0) ) + cim_ri(1)
        iradius = iradius + 1
     enddo
  elseif(imethod.eq.1) then
     !
     ! The newer (RAYTRACE) methods for arranging the rays.
     ! The idea here is to make sure that the inner edge is better
     ! sampled. Also, at some point, I want to introduce here the
     ! method that allows the star to be sampled in a better way,
     ! in case the inner edge of the disk/envelope is not so far 
     ! away from the stellar surface.
     !
     if(nrref.le.0) then
        write(stdo,*) 'ERROR: If new method for ray-arrangement'
        write(stdo,*) '  is chosen, then nrref must be set>0, too'
        stop
     endif
     refdum1 = 0.5d0
     refdum2 = 0.d0
     do ix=1,nrrextra+nrref 
        if(iradius.gt.cim_nr) then
           write(stdo,*) 'BUG in setup_pixels_circular()'
           write(stdo,*) '    iradius exceeds max.'
           stop 155
        endif
        !
        ! A method by which a minimal grid fineness is guaranteed,
        ! but also refinement near the inner edge is done.
        ! 
        if(amr_grid_xi(1,1).le.cim_ri(1)) then
           write(stdo,*) 'ERROR in circular images: Grid inner edge is not'
           write(stdo,*) '      larger than stellar radius.'
           stop
        endif
        if(ix.le.nrrextra) then
           cim_rc(iradius) = ( ix * ( amr_grid_xi(1,1) - &
                             cim_ri(1) ) / (nrrextra+1.d0) ) + cim_ri(1)
        else
           refdum2 = refdum2 + refdum1
           refdum1 = refdum1/2
           cim_rc(iradius) = ( ( nrrextra + refdum2 )            &
                             * ( amr_grid_xi(1,1) - cim_ri(1) ) /    &
                             (nrrextra+1.d0) ) + cim_ri(1)
        endif
        iradius = iradius + 1
     enddo
  else
     stop 33720
  endif
  !
  ! Now make the usual rays, which are the rays with impact parameters
  ! larger than the inner edge, but smaller than the outer edge of the
  ! spatial object. These rays will contain most of the information of 
  ! the image, once the telescope imaging routine has done it's job.
  !
  ! irmax = max(irmax,amr_grid_nx)   ! Bugfix 20.03.2017
  irmax = min(irmax,amr_grid_nx)
  if(irmax.eq.0) irmax=amr_grid_nx   ! 0 means: till the end of the grid
  do ix=irmin,irmax-1
     if(iradius.gt.cim_nr) then
        write(stdo,*) 'BUG setup_pixels_circular(): exceed limits'
        stop 155
     endif
     cim_rc(iradius) = amr_grid_xi(ix,1) * (1.d0 + epsrrr )
     iradius = iradius + 1
     !
     ! If I am not at the last point, and if the dbdr
     ! gt 1 then I will add a few points in between this
     ! R_i and R_i+1.... And I will simply count the amount
     ! of added points....
     !
     if((ix.lt.irmax).and.(dbdr.gt.1)) then
        dr = ( amr_grid_xi(ix+1,1) - amr_grid_xi(ix,1) ) / &
             ( 1.d0*dbdr )
        do iins = 1,dbdr-1
           if(iradius.gt.cim_nr) then
              write(stdo,*) 'BUG in circular images: exceeded limit.'
              stop 155
           endif
           cim_rc(iradius) = amr_grid_xi(ix,1) + iins*dr
           iradius = iradius + 1
        enddo
     endif
  enddo
  if(iradius-1.gt.cim_nr) stop 4
  cim_nr  = iradius - 1
  !
  ! Compute the pixel cell walls 
  !
  do ix=2,cim_nr
     cim_ri(ix) = 0.5d0*(cim_rc(ix)+cim_rc(ix-1))
  enddo
  if(cim_ri(1).ge.cim_rc(1)) stop 4091
  cim_ri(cim_nr+1) = 2*cim_ri(cim_nr) - cim_ri(cim_nr-1) 
  !
  ! Set some auxiliary stuff
  !
  camera_image_halfsize_x = cim_ri(cim_nr+1)
  camera_image_halfsize_y = cim_ri(cim_nr+1)
  camera_image_nx = 0
  camera_image_ny = 0
  camera_image_nr   = cim_nr
  camera_image_nphi = cim_np
  !
  ! Calculate the surface areas of these "pixels"
  !
  do ir=0,camera_image_nr
     do iphi=1,camera_image_nphi
        cim_surfarea(ir,iphi) = 0.5d0 * ( cim_ri(ir+1)**2 - cim_ri(ir)**2 ) * &
             ( cim_pi(iphi+1) - cim_pi(iphi) )
     enddo
  enddo
  !
end subroutine setup_pixels_circular

!-------------------------------------------------------------------------
!                         MAKE A CIRCULAR IMAGE
! 
! In parcicular for models in spherical coordinates it can be useful to
! arrange the pixels not in x and y, but instead in r and phi. This allows
! one to automatically adapt to the refining grid toward the origin. It
! also makes it easier to analyze the results of 1-D models, because you
! will then merely get intensity as a function of radius instead of the
! overkill of a full x,y image.
!
! The method used for setting up the pixels is the "tangent ray method",
! with several additions. 
!
! DEPENDENCE ON SETTINGS:
!    camera_circ_nrphiinf  The number of pixels arranged in each circle (i.e.
!                          the number of phi-pixels). For 1-D spherically symmetric
!                          models this can be set to 1, because we do not expect 
!                          any differences in the image along the phi-direction 
!                          for spherically symmetric models. Recommended value
!                          for 2-D axisymmetric and 3-D models: something of the
!                          order of 128. 
!    camera_circ_nrext     The number of extra rays inside of the inner edge
!                          between the stellar surface and the grid inner edge.
!                          Recommended value: something like 10.
!    camera_circ_dbdr      For the original tangent-ray method this should be 
!                          set to 1. If set to 2 or higher, then for each radial
!                          coordinate in the spherical coordinates, extra rays
!                          are inserted. For dbdr=2 there will be 2 radial 
!                          pixel circles for each radial coordinate grid point. 
!                          More accurate but substantially slower. Reccomended
!                          value is 1.
!    camera_circ_imethod   Set this to 1 (the other methods are for backward
!                          compatibility with older radiative transfer programs).
!    camera_circ_nrref     Further refinement near inner grid edge. Recommended
!                          value is 1 (no further refinement).
!-------------------------------------------------------------------------
subroutine camera_make_circ_image()
  use amr_module
  implicit none
  double precision :: r,phi,px,py,x,y,z
  integer :: inu,ir,iphi,ierr,inu0,inu1,ierror,ispec
  double precision :: celldxmin,quvsq
  double precision :: dirx,diry,dirz,distance,xbk,ybk,zbk,svcx,svcy,svcz
  logical :: domc
  character*80 :: strint
  integer :: iact,icnt,ilinesub
  logical :: redo
  integer :: irmin,irmax
  !
  ! Set some defaults
  !
  irmin = 1
  irmax = amr_grid_nx
  !
  ! Circular images are only possible if you use spherical coordinates
  !
  if((igrid_coord.lt.100).or.(igrid_coord.ge.200)) then
     write(stdo,*) 'ERROR: Circular images only allowed for spherical coordinates.'
     stop
  endif
  !
  ! Subpixeling not allowed in spherical images
  !
  if(camera_diagnostics_subpix) then
     write(stdo,*) 'ERROR: Subpixeling is not allowed in circular images.'
     stop
  endif
  !
  ! Local observer not allowed in circular images
  !
  if(camera_localobserver) then
     write(stdo,*) 'ERROR: local observer mode not allowed in circular images.'
     stop
  endif
  !
  ! For now we do not allow more than 1 star in the spherical images
  !
  if(nstars.gt.1) then
     write(stdo,*) 'ERROR: In circular images at the moment only one star allowed.'
     stop
  endif
  !
  ! Set up the circular/radial pixel arrangement for the circular image
  ! 
  call setup_pixels_circular(irmin,irmax,                            &
                             camera_circ_nrphiinf,camera_circ_nrext, &
                             camera_circ_dbdr,camera_circ_imethod,   &
                             camera_circ_nrref)
  !
  ! If the dust emission is included, then make sure the dust data,
  ! density and temperature are read. If yes, do not read again.
  !
  if(rt_incl_dust) then
     call read_dustdata(1)
     call read_dust_density(1)
     call read_dust_temperature(1)
     if(alignment_mode.ne.0) then
        call aligned_grains_init(1)
     endif
  endif
  !
  ! If line emission is included, then make sure the line data are
  ! read. If yes, then do not read it again.
  !
  if(rt_incl_lines) then
     call read_lines_all(1)
  endif
  !
  ! Do a check
  !
  if(rt_incl_lines.and.(lines_maxshift.le.0.d0)) then
     write(stdo,*) 'INTERNAL ERROR in lines: lines_maxshift variable not set'
     stop
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
  ! Note: it only stores those levels which are selected in the "subset"
  ! (see lines_module.f90). The idea of the subset is that to save memory
  ! you may not want to store all level populations.  Example: A molecule
  ! may have 30 relevant levels. If you need to store all populations
  ! globally, then you need 30 x 8 bytes x nrcells of memory. For large
  ! grids that could be very much. If, however, you wish to only model one
  ! of the lines, then only 2 level populations have to be stored at each
  ! point (the upper and the lower level belonging to that line). If you
  ! select these 2 levels as your "subset" then only these 2 levels will be
  ! stored in the global lines_levelpop() array, saving a lot of memory. For
  ! the LVG method all levels are needed *locally* to compute the
  ! populations, but once these populations are computed, only the 2
  ! relevant ones are then stored in the *global* array lines_levelpop().
  !
  ! You can select this subset manually (in the lines.inp file) or
  ! automatically (by calling the subroutine
  ! lines_automatic_subset_selection() in the lines_module.f90).  The latter
  ! is done just above here.
  !
  if(rt_incl_lines) then
     if((lines_mode.ge.1).and.(lines_mode.le.9)) then
        !
        ! Level subset selection
        !
        if(lines_autosubset) then
           !
           ! If requested, do an automatic subset selection of the
           ! molecular levels
           ! 
           call lines_automatic_subset_selection(camera_nrfreq,     &
                camera_frequencies,1,camera_nrfreq,lines_maxshift,  &
                redo)
           !
           ! Now force a recomputation of the populations, unless
           ! "redo" is .false., meaning that we have the same
           ! levels as before (and thus the population as computed
           ! before must be still correct). Note that if there
           ! exists no lines_levelpop() array, then redo will also
           ! be .true.
           !
           if(redo) then
              iact=2
           else
              iact=1
           endif
           call lines_compute_and_store_local_populations(iact)
        else
           !
           ! Subset is not automatically selected. So either it has
           ! been selected manually or not at all (in which case the
           ! "subset" is the full set of levels).
           !
           write(stdo,*) 'Will store level populations for all levels or ',   &
                'for a manually (in lines.inp) selected subset of levels.'
           !
           ! Now compute the populations only if they have not yet been
           ! computed. If the lines_popul() array is present, we know that
           ! it contains the correct populations, because the subset
           ! is fixed by the user (or is the complete set). 
           !
           call lines_compute_and_store_local_populations(1)
           !
        endif
     endif
  endif
  !
  ! Do the initializing of the camera
  !
  call camera_init()
  !
  ! Warnings
  !
  if(rt_incl_lines.and.(lines_mode.lt.0).and.(scattering_mode.ne.0)) then
     write(stdo,*) 'WARNING: Using dust scattering AND line transfer with '
     write(stdo,*) '         on-the-fly level population determination can'
     write(stdo,*) '         make RADMC-3D very slow, because dust scattering'
     write(stdo,*) '         means that the transfer must be done freq-by-freq,'
     write(stdo,*) '         and thus the populations must be recalculated'
     write(stdo,*) '         at each freq, which is slowing down the code.'
     write(stdo,*) '         Various possible solutions:'
     write(stdo,*) '          1. If possible, switch off dust scattering'
     write(stdo,*) '             (for instance for far-IR or mm wavelengths)'
     write(stdo,*) '             by setting scattering_mode_max=0 in radmc3d.inp'
     write(stdo,*) '          2. Use not-on-the-fly populations (lines_mode>0)'
     write(stdo,*) '             [at the moment this mode is still in prep]'
  endif
  !
  ! Check
  !
  if(camera_nrfreq.lt.1) then
     write(stdo,*) 'ERROR in camera module: camera frequency array'
     write(stdo,*) '      not set when calling camera_make_circ_image()'
     stop
  endif
  !
  ! Pre-compute cosines and sines
  ! For the camera-at-infinity view the angles are given and the
  ! cos and sin can be calculated easily
  !
  camera_pointing_cos_posang  = cos(camera_pointing_degr_posang*pi/180.)
  camera_pointing_sin_posang  = sin(camera_pointing_degr_posang*pi/180.)
  camera_observer_cos_theta   = cos(camera_observer_degr_theta*pi/180.)
  camera_observer_sin_theta   = sin(camera_observer_degr_theta*pi/180.)
  camera_observer_cos_phi     = cos(camera_observer_degr_phi*pi/180.)
  camera_observer_sin_phi     = sin(camera_observer_degr_phi*pi/180.)
  !
  ! For now we do not allow multiple vantage points
  !
  if(mcscat_nrdirs.gt.1) then
     write(stdo,*) 'ERROR: Multiple vantage points is not yet allowed in this version!'
     stop
  endif
  !
  ! Check if we need to (re)do the one-frequency scattering Monte Carlo simulation
  !
  domc = .false.
  if(scattering_mode.ge.1) then
     !
     ! For now scattering is not allowed in local observer mode
     !
     if(camera_localobserver) then
        write(stdo,*) 'ERROR: At the moment scattering is not'
        write(stdo,*) '       allowed for local observer mode.'
        stop
     endif
     !
     ! Compute the direction vector of the observer. Start with (0,0,1).
     ! Actually this is only necessary for scattering_mode.ge.2, but alas.
     !
     dirx = 0.d0
     diry = 0.d0
     dirz = 1.d0
     !
     ! A rotation in the y,z plane, i.e. going from pole-on
     ! to more face-on (if a disk is assumed to be present in the x-y plane)
     !
     ybk  = diry
     zbk  = dirz
     diry = camera_observer_cos_theta * ybk - camera_observer_sin_theta * zbk
     dirz = camera_observer_sin_theta * ybk + camera_observer_cos_theta * zbk
     !
     ! Then a rotation in the x,y plane, rotating horizontally around
     ! the object. If the camera looks toward the object, then the
     ! camera now moves to the left (clockwise around the object).
     !
     xbk  = dirx
     ybk  = diry
     dirx = camera_observer_cos_phi * xbk + camera_observer_sin_phi * ybk
     diry =-camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
     !
     ! Checks
     !
     if(abs(dirx**2+diry**2+dirz**2-1.d0).gt.1d-6) then
        write(stdo,*) 'ERROR in camera module: direction vector not OK.'
        write(stdo,*) dirx,diry,dirz
        stop
     endif
     !
     ! Compute the S-vector of observer for polarization. Start with (0,1,0).
     ! Actually this is only necessary for scattering_mode.ge.4, but it
     ! never hurts.
     !
     svcx = 0.d0
     svcy = 1.d0
     svcz = 0.d0
     !
     ! Rotate camera along its axis 
     !
     ! Note that the camera is rotated in clockwise direction, so any
     ! image is rotated in counter-clockwise direction on the CCD
     !
     xbk  = svcx
     ybk  = svcy
     svcx = camera_pointing_cos_posang * xbk + camera_pointing_sin_posang * ybk
     svcy =-camera_pointing_sin_posang * xbk + camera_pointing_cos_posang * ybk
     !
     ! Then a rotation in the y,z plane, i.e. going from pole-on
     ! to more face-on (if a disk is assumed to be present in the x-y plane)
     !
     ybk  = svcy
     zbk  = svcz
     svcy = camera_observer_cos_theta * ybk - camera_observer_sin_theta * zbk
     svcz = camera_observer_sin_theta * ybk + camera_observer_cos_theta * zbk
     !
     ! Then a rotation in the x,y plane, rotating horizontally around
     ! the object. If the camera looks toward the object, then the
     ! camera now moves to the left (clockwise around the object).
     !
     xbk  = svcx
     ybk  = svcy
     svcx = camera_observer_cos_phi * xbk + camera_observer_sin_phi * ybk
     svcy =-camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
     !
     ! Checks
     !
     if(abs(svcx**2+svcy**2+svcz**2-1.d0).gt.1d-6) then
        write(stdo,*) 'ERROR in camera module: S-vector not OK.'
        write(stdo,*) svcx,svcy,svcz
        stop
     endif
     if(abs(svcx*dirx+svcy*diry+svcz*dirz).gt.1d-6) then
        write(stdo,*) 'INTERNAL ERROR: Somehow the S-vector is not '
        write(stdo,*) '  perpendicular to the direction vector... Warn author.'
        stop
     endif
     !
     ! Some simple tests to see if a new MC calculation is necessary
     !
     if((.not.allocated(mcscat_dirs)).and.(scattering_mode.gt.1)) domc=.true.
     if(allocated(mcscat_dirs).and.(scattering_mode.eq.1)) then
        write(stdo,*) 'ERROR: if isotropic scattering is set, then do not set the mcscat_dirs(:,:) array...'
        stop
     endif
     if(mcscat_nrdirs.ne.1) domc=.true.
     if(.not.camera_scatsrc_allfreq) domc=.true.
     !
     ! Some more subtle tests to see if a new MC calculation is necessary
     !
     if(camera_scatsrc_allfreq) then
        if(mc_nrfreq.ne.camera_nrfreq) then
           domc=.true.
        else
           do inu=1,camera_nrfreq
              if(abs((camera_frequencies(inu)-mc_frequencies(inu))/   &
                     (camera_frequencies(inu)+mc_frequencies(inu))).gt.1d-10) domc=.true.
           enddo
        endif
        if(allocated(mcscat_dirs)) then
           mcscat_current_dir = 1
           if((abs(mcscat_dirs(1,mcscat_current_dir)-dirx).gt.1d-6).or.   &
              (abs(mcscat_dirs(2,mcscat_current_dir)-diry).gt.1d-6).or.   &
              (abs(mcscat_dirs(3,mcscat_current_dir)-dirz).gt.1d-6)) then
              domc=.true.
           endif
        endif
     endif
  endif
  !
  ! Prepare the direction for the scattering source function
  ! For now allow only one single direction (i.e. one single vantage point)
  !
  if(domc.and.(scattering_mode.gt.1)) then
     if((igrid_coord.ge.100).and.(amr_dim.ne.3)) then
        !
        ! Special case: 2D axisymmetric model in spherical coordinate
        ! but with full scattering mode (only possible for scattering_mode.ge.5).
        !
        if(amr_dim.eq.1) then
           write(stdo,*) 'ERROR: scattering_mode.ge.2 is incompatible with'
           write(stdo,*) '       1-D spherical coordinates.'
           stop
        endif
        if(scattering_mode.lt.5) then
           write(stdo,*) 'ERROR: scattering_mode.lt.5 is incompatible with'
           write(stdo,*) '       2-D spherical coordinates.'
           stop
        endif
        if(camera_secondorder) then
           write(stdo,*) 'ERROR: At the moment the 2-D axisymmetric full-scattering mode is not yet'
           write(stdo,*) '       compatible with second order ray-tracing... :-('
           stop
        endif
        !
        ! Switch on the special treatment
        !
        dust_2daniso = .true.
        !
        ! Extend the scattering source array to dust_2daniso_nphi + 1 phi-angle
        ! points starting with phi=0 and ending with phi=360 degrees
        !
        mcscat_nrdirs = dust_2daniso_nphi + 1
        !
        ! Make sure that the ray tracing cannot make larger steps
        ! than a maximum angle wrt the origin. Reason: for near edge-on
        ! views a ray could otherwise pass through a cell (=annulus) 
        ! almost along the annulus tube, and change phi angle too much,
        ! thereby skipping intermediate angles. That is bad for the 
        ! interpolation of the scattering source function.
        !
        if(camera_maxdphi.eq.0.d0) then
           write(stdo,*) 'WARNING: 2-D anisotropic scattering without camera_maxdphi set... Dangerous.'
        endif
        if(camera_maxdphi.gt.twopi/dust_2daniso_nphi) then
           camera_maxdphi = twopi / dust_2daniso_nphi
        endif
        if(camera_maxdphi.ge.0.5d0) camera_maxdphi=0.5d0
        !
        ! Message
        !
        write(stdo,*) 'Note: Using 2-D full-phase scattering mode. This requires a bit of extra memory.'
     else
        !
        ! Normal case (3-D)
        !
        mcscat_nrdirs = 1
     endif
     !
     ! Make the scattering direction the same as the direction of
     ! viewing. 
     ! For now allow only one single direction (i.e. one single vantage point)
     !
     ! NOTE: For 2-D axisymmetric models in spherical coordinates, we
     !       do things in a special way: we reserve mcscat_nrdirs 
     !       "directions", which are in fact all the same, but for
     !       the scattering source function we will set the scattering
     !       event at different positions. 
     !
     if(allocated(mcscat_dirs)) deallocate(mcscat_dirs)
     allocate(mcscat_dirs(1:3,1:mcscat_nrdirs),STAT=ierr)
     if(ierr.gt.0) then
        write(stdo,*) 'ERROR: Could not allocate mcscat_dirs()'
        stop
     endif
     mcscat_current_dir = 1
     mcscat_dirs(1,1:mcscat_nrdirs) = dirx
     mcscat_dirs(2,1:mcscat_nrdirs) = diry
     mcscat_dirs(3,1:mcscat_nrdirs) = dirz
     !
     ! For convenience, store this also here in the camera module
     ! (works only for a single vantage point)
     !
     camera_dir(1) = dirx
     camera_dir(2) = diry
     camera_dir(3) = dirz
     !
     ! In case you want to include polarization, we set also the
     ! S-vectors, which will be perpendicular to the direction 
     ! vector, and pointing vertically upward in the image plane.
     !
     if(allocated(mcscat_svec)) deallocate(mcscat_svec)
     allocate(mcscat_svec(1:3,1:mcscat_nrdirs),STAT=ierr)
     if(ierr.gt.0) then
        write(stdo,*) 'ERROR: Could not allocate mcscat_svec()'
        stop
     endif
     mcscat_svec(1,1:mcscat_nrdirs) = svcx
     mcscat_svec(2,1:mcscat_nrdirs) = svcy
     mcscat_svec(3,1:mcscat_nrdirs) = svcz
     !
     ! For convenience, store this also here in the camera module
     ! (works only for a single vantage point)
     !
     camera_svec(1) = svcx
     camera_svec(2) = svcy
     camera_svec(3) = svcz
     !
     ! But mirror symmetry is not allowed for anisotropic scattering
     !
     if((igrid_coord.ge.100).and.(igrid_coord.le.199)) then
        if(igrid_mirror.ne.0) then
           write(stdo,*) 'ERROR: Mirror symmetry not compatible with anisotropic scattering.'
           stop
        endif
     endif
  endif
  !
  ! Reset flag
  !
  camera_warn_resolution = .false.
  !
  ! Message
  !
  write(stdo,*) 'Rendering circular image(s)...'
  !
  ! Now make the image. We can do so in two different ways. The first is to
  ! do all frequencies (wavelengths) at once. If scattering is included this
  ! requires the precomputing of the scattering source function at all
  ! wavelengths.  This can have a too large memory requirement.  The second
  ! is to do it frequency-by-frequency, and doing a MC calculation for each
  ! of them separately. That saves a lot of computer memory, so it is the
  ! preferred method for large simulations when making spectra. But it may
  ! be slower. The way to choose one or the other is with the switch called
  ! camera_scatsrc_allfreq.
  !
  if(((scattering_mode.eq.0).or.camera_scatsrc_allfreq).and.   &
      (.not.camera_secondorder)) then
     !
     ! Multi-wavelength Method 1:
     !
     ! Do all frequencies at once for each pixel. This is the preferred
     ! method for cases without scattering. With scattering this requires a
     ! lot of memory to store the scattering source function for all
     ! frequencies, so it might not be the best for that.
     !
     ! Note that at present this multi-frequency mode is not available
     ! for second order integration of the transfer equation (i.e. with
     ! the corner points). 
     !
     inu0 = 1
     inu1 = camera_nrfreq
     !
     ! Tell the ray tracer that if scattering is done, then it is done
     ! frequency dispersed
     !
     camera_mcscat_monochromatic = .false.
     !
     ! If necessary, then do the scattering source functions at all
     ! frequencies beforehand.  WARNING: This can be a very large array!
     !
     if(domc) then
        if(allocated(mc_frequencies)) deallocate(mc_frequencies)
        mc_nrfreq=camera_nrfreq
        allocate(mc_frequencies(1:mc_nrfreq),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate mc_frequencies(:) array'
           stop
        endif
        mc_frequencies(:) = camera_frequencies(:)
        !
        ! Message
        !
        write(stdo,*) 'Doing scattering Monte Carlo simulation...'
        call flush(stdo)
        !
        ! Do Monte Carlo simulation
        !
        if(camera_lambda_starlight_single_scat_mode.eq.0) then
           call do_monte_carlo_scattering(rt_mcparams,ierror,do_resetseed,&
                                          scatsrc=.true.)
           write(stdo,*) 'Average number of scattering events per photon package = ', &
                      ieventcounttot/(1.d0*rt_mcparams%nphot_scat)
           !
           ! If the Monte Carlo settings are very conservative, then give a warning
           ! that you may want to change this (but at your own risk). 
           !
           if(mc_scat_maxtauabs.gt.5.d0) then
              write(stdo,*) 'Tip for speed-up: By default the settings of RADMC-3D are conservative (i.e. safe but slow).'
              write(stdo,*) '   A photon package in monochromatic Monte Carlo is only destroyed after tau_abs = ',mc_scat_maxtauabs
              write(stdo,*) '   In most cases, however, an optical depth limit of 5 is enough.'
              write(stdo,*) '   You can (though at your own risk) speed this up by adding the following line to radmc3d.inp:'
              write(stdo,*) '   mc_scat_maxtauabs = 5.d0'
           endif
        elseif(camera_lambda_starlight_single_scat_mode.eq.1) then
           call do_lambda_starlight_single_scattering(rt_mcparams,ierror,scatsrc=.true.)
        elseif(camera_lambda_starlight_single_scat_mode.eq.2) then
           call do_lambda_starlight_single_scattering_simple(rt_mcparams,ierror,scatsrc=.true.)
        else
           write(stdo,*) 'Lambda single scattering mode cannot be other than 0 or 1 or 2 for now.'
           stop 8762
        endif
     endif
     !
     ! Pre-compute which lines and which levels for line transfer may
     ! contribute to these wavelengths. Note that this only has to be
     ! pre-computed if the serial ray tracing is used, but we do it 
     ! nevertheless, as it will not hurt.
     !
     if(rt_incl_lines) then
        ! ------------------------------------------------------------------
        ! Attila Juhasz
        ! lines_find_active_lines_leves -> if all energy levels are known - leiden
        ! format for line data
        ! lines_find_active_lines_linelist -> for linelist mode
        ! ------------------------------------------------------------------
        !   call lines_find_active_lines_levels(camera_nrfreq,             &
        !                   camera_frequencies,inu0,inu1,lines_maxshift)
        !
        if(lines_maxnrlevels.gt.0) then
           call lines_find_active_lines_levels(camera_nrfreq,             &
                camera_frequencies,inu0,inu1,lines_maxshift)
        else
           call lines_find_active_lines_linelist(camera_nrfreq,             &
                camera_frequencies,inu0,inu1,lines_maxshift)
        endif
        ! ------------------------------------------------------------------
     endif
     !
     ! Message
     !
     if(inu0.eq.inu1) then
        write(stdo,*) 'Ray-tracing image for lambda = ', &
                1d4*cc/camera_frequencies(inu0),' micron...'
     else
        call integer_to_string(abs(inu1-inu0)+1,strint)
        write(stdo,*) 'Ray-tracing images: all '//trim(strint)//' wavelength at once...'
     endif
     !
     ! Now make the image at all wavelengths simultaneously
     !
     call camera_make_circ_image_sub(inu0,inu1)
     !
  else
     !
     ! Multi-wavelength Method 2:
     !
     ! If scattering is included, we need the scattering source function. Since this 
     ! array is easily too big to be stored for all wavelengths: In this method we
     ! compute the scattering source function for each frequency separately, and then
     ! make the corresponding image, and then go to the next frequency.
     !
     ! Tell the ray tracer that if scattering is done, then it is done
     ! one frequency at a time, i.e. the scattering source function is
     ! just one frequency bin.
     !
     camera_mcscat_monochromatic = .true.
     !
     ! Since we recompute the scattering source function for each wavelength, the domc
     ! is by default true. Note that if you make an image at a single wavelength you 
     ! may want to be able to make use of a previously computed source function, if for
     ! instance you make your next image at the same vantage point, the same wavelength
     ! but at a different zoom factor or so. In that case, please set the flag
     ! camera_scatsrc_allfreq to .true.. 
     !
     if(domc) then
        if(allocated(mc_frequencies)) deallocate(mc_frequencies)
        mc_nrfreq=1
        allocate(mc_frequencies(1:1),STAT=ierr)
     endif
     !
     ! Do a loop over all camera frequencies
     !
     do inu0=1,camera_nrfreq
        !
        ! Set this by default
        !
        inu1=inu0
        !
        ! If we must do Monte Carlo, then do this
        !
        if(domc) then
           !
           ! Set the wavelength for the Monte Carlo scattering simulation
           !
           mc_frequencies(1) = camera_frequencies(inu0)
           !
           ! Message
           !
           write(stdo,*) 'Doing scattering Monte Carlo simulation for lambda = ', &
                1d4*cc/mc_frequencies(1),' micron...'
           call flush(stdo)
           !
           ! Call the single wavelength Monte Carlo module
           !
           if(camera_lambda_starlight_single_scat_mode.eq.0) then
              call do_monte_carlo_scattering(rt_mcparams,ierror,do_resetseed,&
                                             scatsrc=.true.)
              write(stdo,*) 'Average number of scattering events per photon package = ', &
                      ieventcounttot/(1.d0*rt_mcparams%nphot_scat)
              !
              ! If the Monte Carlo settings are very conservative, then give a warning
              ! that you may want to change this (but at your own risk). 
              !
              if(mc_scat_maxtauabs.gt.5.d0) then
                 write(stdo,*) 'Tip for speed-up: By default the settings of RADMC-3D are conservative ', &
                      '(i.e. safe but slow).'
                 write(stdo,'(A68,A16,F6.2)') '   A photon package in monochromatic Monte Carlo is only destroyed ', &
                      'after tau_abs = ',mc_scat_maxtauabs
                 write(stdo,*) '   In most cases, however, an optical depth limit of 5 is enough.'
                 write(stdo,*) '   You can (though at your own risk) speed this up by adding the following ', &
                      'line to radmc3d.inp:'
                 write(stdo,*) '   mc_scat_maxtauabs = 5.d0'
              else
                 if(mc_scat_maxtauabs.gt.2.d0) then
                    write(stdo,'(A36,F6.2,A36)') ' Warning: Using mc_scat_maxtauabs = ',mc_scat_maxtauabs, &
                         ' (this is fine, but be aware of it).'
                 else
                    write(stdo,'(A36,F6.2,A34)') ' ERROR: Using mc_scat_maxtauabs = ',mc_scat_maxtauabs, &
                         ': This is too low...'
                 endif
              endif
           elseif(camera_lambda_starlight_single_scat_mode.eq.1) then
              call do_lambda_starlight_single_scattering(rt_mcparams,ierror,scatsrc=.true.)
           elseif(camera_lambda_starlight_single_scat_mode.eq.2) then
              call do_lambda_starlight_single_scattering_simple(rt_mcparams,ierror,scatsrc=.true.)
           else
              write(stdo,*) 'Lambda single scattering mode cannot be other than 0 or 1 or 2 for now.'
              stop 8762
           endif
        endif
        !
        ! Pre-compute which lines and which levels for line transfer may
        ! contribute to these wavelengths. Note that this only has to be
        ! pre-computed if the serial ray tracing is used, but we do it 
        ! nevertheless, as it will not hurt.
        !
        if(rt_incl_lines) then
           !
           ! Check which lines to include
           !
           ! ------------------------------------------------------------------
           ! Attila Juhasz
           ! lines_find_active_lines_leves -> if all energy levels are known - leiden
           ! format for line data
           ! lines_find_active_lines_linelist -> for linelist mode
           ! ------------------------------------------------------------------
           !   call lines_find_active_lines_levels(camera_nrfreq,             &
           !                   camera_frequencies,inu0,inu1,lines_maxshift)
           !
            if(lines_maxnrlevels.gt.0) then
                call lines_find_active_lines_levels(camera_nrfreq,             &
                     camera_frequencies,inu0,inu1,lines_maxshift)
            else
                call lines_find_active_lines_linelist(camera_nrfreq,             &
                     camera_frequencies,inu0,inu1,lines_maxshift)
            endif

           ! ------------------------------------------------------------------
           !
           ! If camera_catch_doppler_line.eq..true., then allocate the
           ! corner-based line quantities
           !
           if(camera_catch_doppler_line) then
              if(allocated(sources_vertex_line_nup)) then
                 deallocate(sources_local_line_ndown_curr)
                 deallocate(sources_local_line_nup_curr)
                 deallocate(sources_local_line_ndown_prev)
                 deallocate(sources_local_line_nup_prev)
                 deallocate(sources_local_line_ndown_end)
                 deallocate(sources_local_line_nup_end)
                 deallocate(sources_cell_line_ndown)
                 deallocate(sources_cell_line_nup)
                 deallocate(sources_vertex_line_ndown)
                 deallocate(sources_vertex_line_nup)
              endif
              sources_vertex_lines_nractivetot = 0
              do ispec=1,lines_nr_species
                 sources_vertex_lines_nractivetot = sources_vertex_lines_nractivetot + &
                      active_nrlines(ispec)
              enddo
              allocate(sources_vertex_line_nup(1:sources_vertex_lines_nractivetot,1:amr_nr_vertices_max))
              allocate(sources_vertex_line_ndown(1:sources_vertex_lines_nractivetot,1:amr_nr_vertices_max))
              allocate(sources_cell_line_nup(1:sources_vertex_lines_nractivetot,1:amr_nr_vertices_max))
              allocate(sources_cell_line_ndown(1:sources_vertex_lines_nractivetot,1:amr_nr_vertices_max))
              allocate(sources_local_line_nup_curr(1:sources_vertex_lines_nractivetot))
              allocate(sources_local_line_ndown_curr(1:sources_vertex_lines_nractivetot))
              allocate(sources_local_line_nup_prev(1:sources_vertex_lines_nractivetot))
              allocate(sources_local_line_ndown_prev(1:sources_vertex_lines_nractivetot))
              allocate(sources_local_line_nup_end(1:sources_vertex_lines_nractivetot))
              allocate(sources_local_line_ndown_end(1:sources_vertex_lines_nractivetot))
           endif
        endif
        !
        ! Message
        !
        write(stdo,*) 'Ray-tracing image for lambda = ', &
                1d4*cc/camera_frequencies(inu0),' micron...'
        !
        ! If we do second order integration of the transfer equation,
        ! we must compute the emissivities at the corner points of the
        ! cells (the vertex grid). 
        !
        if(camera_secondorder) then
           !
           ! Set some flags and values
           !
           sources_localobserver = camera_localobserver
           !
           ! Do some preparations 
           !
           if(.not.camera_localobserver) then
              !
              ! Set the direction vector, if observer is at infinity
              ! 
              ! Note: this is necessary only for the inclusion of the scattering
              !       source function from the Monte Carlo module. This is
              !       anyway unavailable for local observer mode. 
              !
              dirx = 0.d0
              diry = 0.d0
              dirz = 1.d0
              ybk  = diry
              zbk  = dirz
              diry = camera_observer_cos_theta * ybk - camera_observer_sin_theta * zbk
              dirz = camera_observer_sin_theta * ybk + camera_observer_cos_theta * zbk
              xbk  = dirx
              ybk  = diry
              dirx = camera_observer_cos_phi * xbk + camera_observer_sin_phi * ybk
              diry =-camera_observer_sin_phi * xbk + camera_observer_cos_phi * ybk
              if(abs(dirx**2+diry**2+dirz**2-1.d0).gt.1d-6) then
                 write(stdo,*) 'ERROR in camera module: direction vector not OK.'
                 write(stdo,*) dirx,diry,dirz
                 stop
              endif
              ray_cart_dirx = dirx
              ray_cart_diry = diry
              ray_cart_dirz = dirz
           else
              !
              ! Local observer mode, so set the observer position
              !
              sources_observer_position(:) = camera_observer_position(:)
           endif
           !
           ! Then call the subroutine
           !
           call sources_compute_snualphanu_at_vertices(inu0,camera_stokesvector)
        endif
        !
        ! Now make the image for this wavelength only
        !
        call camera_make_circ_image_sub(inu0,inu1)
        !
     enddo
     !
  endif
  !
  !--------------------------------------------------------------
  !        A sub-subroutine for making the image
  !--------------------------------------------------------------
  !
  contains
  subroutine camera_make_circ_image_sub(inu00,inu11)
    implicit none
    integer :: inu00,inu11
    integer :: inuu
    integer :: backup_nrrefine,backup_tracemode
    logical :: warn_tausurf_problem,flag_quv_too_big
    double precision :: r,phi
    !
    ! *** NEAR FUTURE: PUT OPENMP DIRECTIVES HERE (START) ***
    !
    flag_quv_too_big = .false.
    do ir=1,camera_image_nr
       r = cim_rc(ir)
       do iphi=1,camera_image_nphi
          phi = cim_pc(iphi)
          !
          ! Set the ray variables
          !
          px = r * cos(phi)
          py = r * sin(phi)
          !
          ! Find the starting point and direction of this ray
          !
          call camera_set_ray(px,py,x,y,z,dirx,diry,dirz,distance)
          !
          ! Reset intensity
          ! 
          if(incl_extlum.eq.0) then
             camera_intensity_iquv(inu00:inu11,1:4) = 0.d0
          else
             do inu=inu00,inu11
                camera_intensity_iquv(inu,1)   = find_extlumintens_interpol(camera_frequencies(inu))
                camera_intensity_iquv(inu,2:4) = 0.d0
             enddo
          endif
          !
          ! Call the ray-tracer 
          !
          call camera_serial_raytrace(camera_nrfreq,inu00,inu11,x,y,z,dirx,diry,dirz, &
                                      distance,celldxmin,camera_intensity_iquv)
          !
          ! Put the result into the image array
          !
          camera_circ_image_iquv(ir,iphi,inu00:inu11,1) = camera_intensity_iquv(inu00:inu11,1)
          if(camera_stokesvector) then
             !
             ! Copy also the other Stokes components
             !
             camera_circ_image_iquv(ir,iphi,inu00:inu11,2) = camera_intensity_iquv(inu00:inu11,2)
             camera_circ_image_iquv(ir,iphi,inu00:inu11,3) = camera_intensity_iquv(inu00:inu11,3)
             camera_circ_image_iquv(ir,iphi,inu00:inu11,4) = camera_intensity_iquv(inu00:inu11,4)
             !
             ! Self-consistency check
             !
             do inuu=inu00,inu11
                quvsq = camera_intensity_iquv(inuu,2)**2 + &
                        camera_intensity_iquv(inuu,3)**2 + &
                        camera_intensity_iquv(inuu,4)**2
                if(quvsq.gt.0.d0) then
                   if(camera_intensity_iquv(inuu,1).eq.0.d0) then
                      write(stdo,*) 'INTERNAL ERROR: Q^2+U^2+V^2>0 but I=0...'
                      write(stdo,*) '    Warn author.'
                      stop
                   endif
                   quvsq = quvsq / camera_intensity_iquv(inuu,1)**2
                   if(quvsq.gt.1.000001d0) then
                      flag_quv_too_big = .true.
                   endif
                endif
             enddo
          endif
          !
       enddo
    enddo
    !
    ! *** NEAR FUTURE: PUT OPENMP DIRECTIVES HERE (FINISH) ***
    !
    if(flag_quv_too_big) then
       write(stdo,*) 'WARNING: While making an image, I found an instance of Q^2+U^2+V^2>I^2...'
    endif
    !
    ! Add the central star (star 1)
    ! 
    if((camera_incl_stars.ne.0).and.(nstars.ge.1)) then
       !
       ! The starting point of this ray is at the stellar surface
       ! Assuming that dirx, diry and dirz are already computed above
       !
       px = 0.d0
       py = 0.d0
       x  = dirx * star_r(1)
       y  = diry * star_r(1)
       z  = dirz * star_r(1)
       !
       ! Reset intensity to the stellar intensity
       ! 
       do inu=inu00,inu11
          camera_intensity_iquv(inu,1)   =                               &
               find_starlight_interpol(camera_frequencies(inu),1)
          camera_intensity_iquv(inu,2:4) = 0.d0
       enddo
       !
       ! Call the ray-tracer 
       !
       call camera_serial_raytrace(camera_nrfreq,inu00,inu11,x,y,z,dirx,diry,dirz, &
                                   distance,celldxmin,camera_intensity_iquv)
       !
       ! Put the result into the image array
       !
       ! Note: since the star has no phi-dependence we put this into all phi
       !
       do iphi=1,camera_image_nphi
          camera_circ_image_iquv(0,iphi,inu00:inu11,1) = camera_intensity_iquv(inu00:inu11,1)
          if(camera_stokesvector) then
             !
             ! Copy also the other Stokes components
             !
             camera_circ_image_iquv(0,iphi,inu00:inu11,2) = camera_intensity_iquv(inu00:inu11,2)
             camera_circ_image_iquv(0,iphi,inu00:inu11,3) = camera_intensity_iquv(inu00:inu11,3)
             camera_circ_image_iquv(0,iphi,inu00:inu11,4) = camera_intensity_iquv(inu00:inu11,4)
          endif
       enddo
    else
       !
       ! If there is no star, then put this to 0. 
       !
       camera_circ_image_iquv(0,1:camera_image_nphi,inu00:inu11,:) = 0.d0
    endif
    !
    ! For now we do not allow more than 1 star in the spherical images.
    ! If I add this later, it will be here.
    !
  end subroutine camera_make_circ_image_sub
  !
end subroutine camera_make_circ_image


!-------------------------------------------------------------------------
!                       WRITE CIRCULAR IMAGE
!-------------------------------------------------------------------------
subroutine camera_write_circ_image(noclip)
  implicit none
  integer :: ns,iinu,ir,iphi,is,nn,kk,iformat
  logical, optional :: noclip
  logical :: donoclip
  double precision :: dummy(1:4)
  !
  ! Interpret the noclip
  !
  donoclip = .false.
  if(present(noclip)) then
     donoclip = noclip
  endif
  !
  ! Open the file
  !
  if(writeimage_unformatted) then
     if(rto_style.eq.2) then
        open(unit=fflo,file='circimage.uout',form='unformatted')
     else
        open(unit=fflo,file='circimage.bout',status='replace',access='stream')
     endif
  else
     open(unit=fflo,file='circimage.out')
  endif
  !
  ! Decide how many of the Stokes components we have to write
  !
  if(camera_stokesvector) then
     ns=4
  else
     ns=1
  endif
  !
  ! Now decide whether unformatted or formatted
  !
  if(writeimage_unformatted) then
     if (rto_style.eq.2) then
        !
        ! Unformatted output (F77-Unformatted output, faster, more compact)
        !
        if(camera_stokesvector) then
           !if(camera_localobserver) then
           !   write(fflo) 4           ! Format number: size units in radian (angular) + Stokes
           !else
              write(fflo) 3           ! Format number: size units in cm (spatial size) + Stokes
           !endif
        else
           !if(camera_localobserver) then
           !   write(fflo) 2           ! Format number: size units in radian (angular)
           !else
              write(fflo) 1           ! Format number: size units in cm (spatial size)
           !endif
        endif
        write(fflo) camera_image_nr,camera_image_nphi
        write(fflo) camera_nrfreq
        write(fflo) (cim_ri(ir),ir=1,cim_nr+1)
        write(fflo) (cim_rc(ir),ir=1,cim_nr)
        write(fflo) (cim_pi(iphi),iphi=1,cim_np+1)
        write(fflo) (cim_pc(iphi),iphi=1,cim_np)
        write(fflo) (1d4*cc/camera_frequencies(iinu),iinu=1,camera_nrfreq)
        do iinu=1,camera_nrfreq
           do is=1,ns
              write(fflo) ((camera_circ_image_iquv(ir,iphi,iinu,is),  &
                          ir=0,camera_image_nr),iphi=1,camera_image_nphi)
           enddo
        enddo
     else
        !
        ! Unformatted output (C-compliant binary output, faster, more compact)
        !
        if(camera_stokesvector) then
           if(camera_localobserver) then
              iformat = 4           ! Format number: size units in radian (angular) + Stokes
           else
              iformat = 3           ! Format number: size units in cm (spatial size) + Stokes
           endif
        else
           if(camera_localobserver) then
              iformat = 2           ! Format number: size units in radian (angular)
           else
              iformat = 1           ! Format number: size units in cm (spatial size)
           endif
        endif
        write(fflo) iformat
        nn = camera_image_nr
        kk = camera_image_nphi
        write(fflo) nn, kk
        nn = camera_nrfreq
        write(fflo) nn
        write(fflo) (cim_ri(ir),ir=1,cim_nr+1)
        write(fflo) (cim_rc(ir),ir=1,cim_nr)
        write(fflo) (cim_pi(iphi),iphi=1,cim_np+1)
        write(fflo) (cim_pc(iphi),iphi=1,cim_np)
        write(fflo) (1d4*cc/camera_frequencies(iinu),iinu=1,camera_nrfreq)
        do iinu=1,camera_nrfreq
           do is=1,ns
              write(fflo) ((camera_circ_image_iquv(ir,iphi,iinu,is),  &
                          ir=0,camera_image_nr),iphi=1,camera_image_nphi)
           enddo
        enddo
     endif
  else
     !
     ! Standard: write it formatted way
     !
     if(camera_stokesvector) then
        if(camera_localobserver) then
           write(fflo,*) 4           ! Format number: size units in radian (angular)
        else
           write(fflo,*) 3           ! Format number: size units in cm (spatial size)
        endif
     else
        if(camera_localobserver) then
           write(fflo,*) 2           ! Format number: size units in radian (angular)
        else
           write(fflo,*) 1           ! Format number: size units in cm (spatial size)
        endif
     endif
     write(fflo,*) camera_image_nr,camera_image_nphi
     write(fflo,*) camera_nrfreq
     write(fflo,*) ' '
     write(fflo,320) (cim_ri(ir),ir=0,cim_nr+1)
     write(fflo,*) ' '
     write(fflo,320) (cim_rc(ir),ir=0,cim_nr)
     write(fflo,*) ' '
     write(fflo,320) (cim_pi(iphi),iphi=1,cim_np+1)
     write(fflo,*) ' '
     write(fflo,320) (cim_pc(iphi),iphi=1,cim_np)
     write(fflo,*) ' '
     write(fflo,320) (1d4*cc/camera_frequencies(iinu),iinu=1,camera_nrfreq)
320  format(E21.14)
     write(fflo,*) ' '
     do iinu=1,camera_nrfreq
        do iphi=1,camera_image_nphi
           do ir=0,camera_image_nr
              if(ns.eq.1) then
                 if((camera_circ_image_iquv(ir,iphi,iinu,1).gt.1d-97).or.donoclip) then
                    write(fflo,333) camera_circ_image_iquv(ir,iphi,iinu,1)
333                 format(E21.14)
                 else
                    write(fflo,*) 0.d0
                 endif
              else
                 dummy(1:4) = camera_circ_image_iquv(ir,iphi,iinu,1:4)
                 do is=1,ns
                    if((abs(dummy(is)).le.1d-97).and.(.not.donoclip)) then
                       dummy(is) = 0.d0
                    endif
                 enddo
                 write(fflo,334) dummy(1:4)
334              format(4(E21.14,1X))
              endif
           enddo
        enddo
        write(fflo,*) ' '
     enddo
  endif
  !
  ! Close file
  !
  close(1)
  !
end subroutine camera_write_circ_image


!-------------------------------------------------------------------------
!              MAKE A SPECTRUM USING THE CIRCULAR IMAGES
!-------------------------------------------------------------------------
subroutine camera_make_circ_spectrum()
  use constants_module
  implicit none
  integer :: inu,ierr,backup_nrfreq
  integer :: backup_lines_mode
  double precision :: pdx,pdy,factor,r_colarea,r,eps,hsx_bk,hsy_bk
  double precision :: flux(1:4)
  integer :: apinu,irout,ir,iphi
  integer :: iact
  logical :: redo
  !
  ! Check
  !
  if((igrid_coord.lt.100).or.(igrid_coord.ge.200)) then
     write(stdo,*) 'ERROR: Spectrum made with circular images only possible ', &
          'for spherical coordinates'
     stop
  endif
  if(rt_incl_lines.and.(lines_maxshift.le.0.d0)) then
     write(stdo,*) 'INTERNAL ERROR in lines: lines_maxshift variable not set'
     stop
  endif
  if((camera_nrfreq.le.0).or.(.not.allocated(camera_frequencies))) then
     write(stdo,*) 'ERROR: Somehow the frequencies for the spectrum are not set.'
     stop
  endif
  if(camera_localobserver) then
     write(stdo,*) 'ABORTING: For the moment we do not allow the making of a spectrum'
     write(stdo,*) '          in the local observer mode. Use observer at infinity mode,'
     write(stdo,*) '          for instance by specifying the incl and phi on the '
     write(stdo,*) '          command line and setting the image size with sizeau or sizepc.'
     stop
  endif
  if((((scattering_mode.ne.0).and.(.not.camera_scatsrc_allfreq)).or.   &
      camera_secondorder).and.(lines_mode.lt.-1).and.rt_incl_lines) then
     write(stdo,*)    '===> Warning: Combination of modes that may make RADMC-3D very slow... Here is why:'
     if((scattering_mode.ne.0).and.(.not.camera_scatsrc_allfreq)) then
        write(stdo,*) '     You are including dust scattering (doing it freq-by-freq), which means that'
     else
        if(camera_catch_doppler_line) then
           write(stdo,*) '     You are using the doppler catching mode, which requires second order '
           write(stdo,*) '     integration of the RT equation, which (for memory reasons) means that'
        else
           write(stdo,*) '     You are using second order integration of the RT equation,'
           write(stdo,*) '     which (for memory reasons) means that:'
        endif
     endif
     write(stdo,*) '     RADMC-3D must make one image for each frequency at a time to compute each '
     write(stdo,*) '     point of the spectrum. At the same time you use on-the-fly computation of'
     write(stdo,*) '     the line level populations, which are not stored. If you have 100 frequencies'
     write(stdo,*) '     in your spectrum, this means that 100 times an image is made, and a 100 times'
     write(stdo,*) '     the populations are RECALCULATED. If the population calculation is heavy, this'
     write(stdo,*) '     means that you waste a huge amount of time. You may want to calculate these'
     write(stdo,*) '     populations ONCE and store them, which can be done using a positive lines_mode'
     write(stdo,*) '     number (see manual). This may require a lot of memory if you have a molecule'
     write(stdo,*) '     with many levels. To overcome this, you can tell RADMC-3D to store only a subset'
     write(stdo,*) '     of the computed level populations: Only those related to the line you wish to'
     write(stdo,*) '     model. See manual.'
  endif
  !
  ! Set some other values
  !
  camera_localobserver = .false.      ! Spectra only for observer at infinity for now
  camera_tracemode = 1
  !
  ! Check if the aperture file is present and read if needed
  !
  if(camera_use_aperture_info) then
     !
     ! Read the file aperture_info.inp
     ! 
     call camera_read_aperture_info()
     if(.not.allocated(camera_aperture_freq)) stop 8341
     if(.not.allocated(camera_aperture_radius_collectarea)) stop 8342
     !
     ! Check if the current wavelengths are all within range
     !
     if((max(camera_frequencies(1),camera_frequencies(camera_nrfreq)).gt.           &
         max(camera_aperture_freq(1),camera_aperture_freq(camera_aperture_nf))).or. &
        (min(camera_frequencies(1),camera_frequencies(camera_nrfreq)).lt.           &
         min(camera_aperture_freq(1),camera_aperture_freq(camera_aperture_nf)))) then
        write(stdo,*) 'ERROR while making spectrum with aperture information:'
        write(stdo,*) '      The aperture information does not cover all wavelengths'
        write(stdo,*) '      that are to be used for the spectrum.'
        stop
     endif
  endif
  !
  ! If the dust emission is included, then make sure the dust data,
  ! density and temperature are read. If yes, do not read again.
  !
  if(rt_incl_dust) then
     call read_dustdata(1)
     call read_dust_density(1)
     call read_dust_temperature(1)
     if(alignment_mode.ne.0) then
        call aligned_grains_init(1)
     endif
  endif
  !
  ! If line emission is included, then make sure the line data are
  ! read. If yes, then do not read it again.
  !
  if(rt_incl_lines) then
     call read_lines_all(1)
  endif
  !
  ! Now make all the circular images
  !
  call camera_make_circ_image()
  !
  ! Now calculate the flux for each of these images
  !
  do inu=1,camera_nrfreq
     !
     ! First find the maximum ir up to which to integrate (aperture)
     !
     if(camera_use_aperture_info) then
        !
        ! Include only the part of the image that lies within the
        ! aperture radius
        !
        call hunt(camera_aperture_freq,camera_aperture_nf, &
                  camera_frequencies(inu),apinu)
        if((apinu.lt.1).or.(apinu.ge.camera_aperture_nf)) then
           write(stdo,*) 'INTERNAL ERROR in spectrum: aperture index out of range...'
           stop 6398
        endif
        eps = (log(camera_frequencies(1))-log(camera_aperture_freq(apinu))) / &
              (log(camera_aperture_freq(apinu+1))-log(camera_aperture_freq(apinu)))
        if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 6399
        r_colarea = exp((1.d0-eps) * log(camera_aperture_radius_collectarea(apinu)) + &
                    eps * log(camera_aperture_radius_collectarea(apinu+1)))
        !
        ! Convert this radius from arcsec into cm for this image, using the
        ! distance specified by the observer
        !
        r_colarea = r_colarea * AU * camera_observer_distance_pc
        call hunt(cim_ri,cim_nr,r_colarea,irout)
        if(irout.lt.1) irout=0
        if(irout.gt.cim_nr) irout=cim_nr
     else
        !
        ! Include the entire image
        !
        irout = cim_nr
     endif
     !
     ! Now calculate the flux
     !
     flux(1:4) = 0.d0
     do ir=0,irout
        do iphi=1,cim_np
           flux(1:4) = flux(1:4) + cim_surfarea(ir,iphi)*camera_circ_image_iquv(ir,iphi,inu,1:4)
        enddo
     enddo
     flux(:) = flux(:) / ( parsec**2 )
     !
     ! Put this into the spectrum
     !
     camera_spectrum_iquv(inu,:)   = flux(:)
  enddo
  !
end subroutine camera_make_circ_spectrum



!==========================================================================
!                   ALIGNED GRAINS POLARIZATION STUFF
!==========================================================================


!-----------------------------------------------------------------------
!      FIRST ORDER INTEGRATION OF RT EQUATION WITH ALIGNED GRAINS
!-----------------------------------------------------------------------
subroutine pol_integrate_rt_aligned(int_iquv,dir,svec,aligndir,aligneff,inu0,inu1,  &
                                    src0,alp0,dustdens,dusttemp,ds)
  implicit none
  integer :: iop,inu0,inu1,inu,ispec
  double precision :: int_iquv(1:sources_nrfreq,1:4)
  double precision :: src0(1:sources_nrfreq,1:4),alp0(1:sources_nrfreq)
  double precision :: dustdens(1:dust_nr_species),dusttemp(1:dust_nr_species)
  double precision :: dir(1:3),svec(1:3),aligndir(1:3),aligneff,aligneff1,ds,freq
  double precision :: salign(1:3),dum
  double precision :: aligned_iquv(1:4),aligned_opuv(1:4)
  double precision :: aligned_jnu_iquv(1:4),aligned_jnu_opuv(1:4)
  double precision :: cosa,sina,cos2a,sin2a,cosb
  double precision :: alpabs,para_alpabs,orth_alpabs,epspo,bpl
  double precision :: para_alpabs_eff,orth_alpabs_eff
  double precision :: aligned_alpha_opuv(1:4),xp(1:4)
  double precision :: aligned_snu_opuv(1:4),exptau_opuv(1:4),exptau1_opuv(1:4)
  !
  !################################
  ! Test normalization; can be removed after testing
  dum = aligndir(1)*aligndir(1)+aligndir(2)*aligndir(2)+aligndir(3)*aligndir(3)
  if(abs(dum-1.d0).gt.1d-6) stop 30511
  dum = svec(1)*svec(1)+svec(2)*svec(2)+svec(3)*svec(3)
  if(abs(dum-1.d0).gt.1d-6) stop 30512
  dum = dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3)
  if(abs(dum-1.d0).gt.1d-6) stop 30513
  !################################
  !
  ! Calculate the cos(eta) between the line-of-sight direction and the
  ! alignment direction.
  !
  cosb      = aligndir(1)*dir(1) + aligndir(2)*dir(2) + aligndir(3)*dir(3)
  !
  ! Compute the temporary s-vector that is aligned with the alignment
  ! direction, but still perpendicular to the line-of-sight direction.
  ! In other words: the projected alignment direction.
  !
  salign(1) = aligndir(1) - cosb*dir(1)
  salign(2) = aligndir(2) - cosb*dir(2)
  salign(3) = aligndir(3) - cosb*dir(3)
  dum       = sqrt(salign(1)*salign(1) + salign(2)*salign(2) + salign(3)*salign(3))
  if(dum.lt.1d-8) then
     salign(1) = svec(1)
     salign(2) = svec(2)
     salign(3) = svec(3)
  else
     salign(1) = salign(1) / dum
     salign(2) = salign(2) / dum
     salign(3) = salign(3) / dum
  endif
  !
  !################################
  ! Test if new vector is indeed perpendicular to n
  dum = salign(1)*dir(1) + salign(2)*dir(2) + salign(3)*dir(3)
  if(abs(dum).gt.1d-4) stop 30514
  dum = salign(1)*salign(1) + salign(2)*salign(2) + salign(3)*salign(3)
  if(abs(dum-1.d0).gt.1d-4) stop 30515
  !################################
  !
  ! Determine the cos(ang) between the line-of-sight and the 
  ! aligndir. This is necessary to compute the ratio of the
  ! absorption opacity in parallel and orthogonal directions
  ! with respect to the salign vector. We are only interested
  ! in the absolute value, because we assume that the grains
  ! are ellipsoidal without top/bottom asymmetry. If cosb=1
  ! then the alignment vector was already in the plane of the
  ! sky, which means that we have (assuming the grains are 
  ! oblate with the minor axis aligned with the alignment
  ! vector) the strongest ratio of parallel vs orthogonal.
  ! If cosb=0, then parallel and orthogonal are the same.
  ! Note: we already computed cosb. Now only abs().
  !
  cosb = abs(cosb)
  !################################
  ! Stupidity test. Can be removed after testing.
  if(cosb.gt.1.0001d0) stop 30517
  !################################
  if(cosb.ge.1.d0) cosb = 0.999999d0
  !
  ! Interpolate the angular grid
  !
  call hunt(sources_align_mu,sources_align_munr,cosb,iop)
  if(iop.ge.sources_align_munr) then
     if(cosb.eq.sources_align_mu(sources_align_munr)) then
        iop = sources_align_munr-1
     else
        write(stdo,*) sources_align_mu(:),cosb
        stop 30518
     endif
  endif
  epspo = (cosb-sources_align_mu(iop)) /                  &
          (sources_align_mu(iop+1)-sources_align_mu(iop)) 
  !################################
  ! Stupidity test. Can be removed after testing.
  if(epspo.gt.1.d0) stop 30519
  if(epspo.lt.0.d0) stop 30520
  !################################
  !
  ! Determine the cos(ang) between this salign and the S-vector 
  ! of the incoming light along the line-of-sight (i.e. the 
  ! global S-vector of the image, i.e. svec(:))
  !
  cosa = salign(1)*svec(1) + salign(2)*svec(2) + salign(3)*svec(3)
  !
  ! Now the cross- and inner product formula for finding sin(a). The sign
  ! convention is such that if the alignment S-vector is
  ! counter-clockwise from the incoming light S-vector when the photon
  ! propagation direction is pointing toward the observer, then the angle
  ! is positive.
  !
  sina = ( svec(2)*salign(3) - svec(3)*salign(2) ) * dir(1) + &
         ( svec(3)*salign(1) - svec(1)*salign(3) ) * dir(2) + &
         ( svec(1)*salign(2) - svec(2)*salign(1) ) * dir(3)
  !
  ! Since we need cos(2*ang) and sin(2*ang) for the rotation of the
  ! Stokes vector, we compute them here.
  !
  cos2a = cosa**2 - sina**2
  sin2a = 2d0*sina*cosa
  !
  !################################
  ! Test if sin2a^2 + cos2a^2 = 1
  dum = cos2a*cos2a + sin2a*sin2a
  if(abs(dum-1.0).gt.1d-6) stop 30516
  !################################
  !
  ! Now loop over frequency
  !
  do inu=inu0,inu1
     !
     ! Rotate the Stokes vector of the incoming photon to the new S-vector
     !
     aligned_iquv(1) =  int_iquv(inu,1)
     aligned_iquv(2) =  int_iquv(inu,2)*cos2a + int_iquv(inu,3)*sin2a
     aligned_iquv(3) = -int_iquv(inu,2)*sin2a + int_iquv(inu,3)*cos2a
     aligned_iquv(4) =  int_iquv(inu,4)
     !
     ! Now change from Stokes IQUV to OPUV with O being orthogonal
     ! to the alignment direction salign, and P being parallel
     ! to the alignment direction salign.
     !
     aligned_opuv(1) =  0.5d0 * ( aligned_iquv(1) + aligned_iquv(2) )
     aligned_opuv(2) =  0.5d0 * ( aligned_iquv(1) - aligned_iquv(2) )
     aligned_opuv(3) =  aligned_iquv(3) 
     aligned_opuv(4) =  aligned_iquv(4) 
     !
     ! Now rotate the Stokes vector of the scattering and non-dust 
     ! source to the new S-vector
     !
     aligned_jnu_iquv(1) =  src0(inu,1)
     aligned_jnu_iquv(2) =  src0(inu,2)*cos2a + src0(inu,3)*sin2a
     aligned_jnu_iquv(3) = -src0(inu,2)*sin2a + src0(inu,3)*cos2a
     aligned_jnu_iquv(4) =  src0(inu,4)
     !
     ! Note, however, that src0 is the source function, not the 
     ! emissivity. So multiply by alpha.
     !
     aligned_jnu_iquv(:)  = aligned_jnu_iquv(:) * alp0(inu)
     !
     ! Now change from Stokes IQUV to OPUV with O being orthogonal
     ! to the alignment direction salign, and P being parallel
     ! to the alignment direction salign.
     !
     aligned_jnu_opuv(1)  =  0.5d0 * ( aligned_jnu_iquv(1) + aligned_jnu_iquv(2) )
     aligned_jnu_opuv(2)  =  0.5d0 * ( aligned_jnu_iquv(1) - aligned_jnu_iquv(2) )
     aligned_jnu_opuv(3)  =  aligned_jnu_iquv(3) 
     aligned_jnu_opuv(4)  =  aligned_jnu_iquv(4) 
     !
     ! Now add to the absorption coefficient the non-dust isotropic absorption
     ! and the angle-averaged scattering opacity
     !
     ! NOTE: Here we still assume that the scattering extinction is
     !       angle-independent! In the future, when also scattering 
     !       off aligned grains is included, this must also be
     !       split into orth and para. That would mean that we would
     !       have to subtract again the angle-averaged scattering
     !       opacity and add the orientation-dependent version. That
     !       is rather "ugly" (it would have been cleaner if the
     !       alp0 would not contain the scattering opacity at all),
     !       but since we decided (ages ago!) that sources_get_src_alp()
     !       returns the source function src=S_nu instead of the 
     !       emissivity src=j_nu, and since there are code-technical
     !       reasons to keep the scattering source function in the
     !       sources_get_src_alp() routine, there is no other way than
     !       this ugly way. 
     !
     aligned_alpha_opuv(:) = alp0(inu)
     !
     ! Now we must do a loop over the dust species
     ! 
     do ispec=1,dust_nr_species
        !
        ! Now the absorption opacity
        !
        ! ...First get the angle-averaged absorption opacity
        !
        alpabs = dustdens(ispec) * sources_dustkappa_a(inu,ispec)
        !
        ! ...Then compute the orth and para versions
        !
        orth_alpabs = alpabs *                                    &
             ( (1.d0-epspo)*sources_align_orth(iop,inu,ispec) +   &
                      epspo*sources_align_orth(iop+1,inu,ispec) )
        para_alpabs = alpabs *                                    &
             ( (1.d0-epspo)*sources_align_para(iop,inu,ispec) +   &
                      epspo*sources_align_para(iop+1,inu,ispec) )
        !
        ! ...The grains may not be perfectly aligned. Here aligneff is the
        !    efficiency of alignment. If aligneff==1.0 then the grains are
        !    perfectly aligned. If aligneff=0.0 the grains are not aligned
        !    at all. Partial alignment (0<aligneff<1) is treated as a linear
        !    sum of aligned and non-aligned opacities.
        !
        aligneff1       = 1.d0 - aligneff
        orth_alpabs_eff = aligneff*orth_alpabs + aligneff1*alpabs
        para_alpabs_eff = aligneff*para_alpabs + aligneff1*alpabs
        !
        ! ...Then add
        !
        dum                   = 0.5d0 * ( orth_alpabs_eff + para_alpabs_eff )
        aligned_alpha_opuv(1) = aligned_alpha_opuv(1) + orth_alpabs_eff
        aligned_alpha_opuv(2) = aligned_alpha_opuv(2) + para_alpabs_eff
        aligned_alpha_opuv(3) = aligned_alpha_opuv(3) + dum
        aligned_alpha_opuv(4) = aligned_alpha_opuv(4) + dum
        !
        ! Now add the polarized thermal emission of the aligned
        ! grains
        !
        ! NOTE: The factor 0.5d0*bpl is because I_orth = 0.5*(I+Q)
        !       and I_para = 0.5*(I-Q), so both take care of 50% of
        !       the thermal emission
        !
        bpl = bplanck(dusttemp(ispec),sources_frequencies(inu))
        aligned_jnu_opuv(1)   = aligned_jnu_opuv(1) + 0.5d0*orth_alpabs_eff*bpl
        aligned_jnu_opuv(2)   = aligned_jnu_opuv(2) + 0.5d0*para_alpabs_eff*bpl
     enddo
     !
     ! Calculate the real "source term" S_nu = j_nu / alpha_nu
     !
     aligned_snu_opuv(:)   = aligned_jnu_opuv(:)/(aligned_alpha_opuv(:)+1d-40)
     !
     ! Compute the exp(-tau) and 1-exp(-tau), where for the latter we take
     ! special care of tau << 1.
     !
     xp(:)                 = aligned_alpha_opuv(:)*ds
     exptau_opuv(:)        = exp(-xp(:))
     exptau1_opuv(:)       = 1.d0-exptau_opuv(:)
     if(xp(1).lt.1d-5) exptau1_opuv(1) = xp(1)
     if(xp(2).lt.1d-5) exptau1_opuv(2) = xp(2)
     if(xp(3).lt.1d-5) exptau1_opuv(3) = xp(3)
     if(xp(4).lt.1d-5) exptau1_opuv(4) = xp(4)
     !
     ! Now the first order integration of the RT equation
     !
     aligned_opuv(:)       = exptau_opuv(:)*aligned_opuv(:) +     &
                             exptau1_opuv(:)*aligned_snu_opuv(:)
     !
     ! Recompute the IQUV
     !
     aligned_iquv(1) = aligned_opuv(1) + aligned_opuv(2)
     aligned_iquv(2) = aligned_opuv(1) - aligned_opuv(2)
     aligned_iquv(3) = aligned_opuv(3)
     aligned_iquv(4) = aligned_opuv(4)
     !
     ! Now rotate the Stokes vector back to the original svec
     !
     int_iquv(inu,1) =  aligned_iquv(1)
     int_iquv(inu,2) =  aligned_iquv(2)*cos2a - aligned_iquv(3)*sin2a
     int_iquv(inu,3) =  aligned_iquv(2)*sin2a + aligned_iquv(3)*cos2a
     int_iquv(inu,4) =  aligned_iquv(4)
  enddo
  !
end subroutine pol_integrate_rt_aligned



end module camera_module



