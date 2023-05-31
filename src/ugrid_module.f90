!==========================================================================
!                       Unstructures Grid Module
!
! In RADMC-3D unstructured grids are supported, provided the grid meets
! the following conditions:
!
!   - Each cell is a polyhedron
!   - Each cell must have a convex shape (i.e. convex polyhedron)
!
! The cells at the surface can be either:
!
!   - Open cells going off to infinity. This is the case for e.g.
!     cells generated from a Voronoi tesselation, but it can also
!     be any other kind of tesselation. The advantage is that there
!     exists no "vacuum". Even a starting point of a ray far away
!     will start inside an "open" cell, so it will be straightforward
!     to follow the ray into the object. Note that open cells are
!     treated as "normal" cells for the geometry, but they are
!     treated as "empty" when it comes to density, temperature and
!     such. In the input files for density, temperature etc the
!     values in these "open" cells are included, but they are
!     ignored. The reason for including them nevertheless is that
!     it makes the association between cell index and values easier:
!     If cell[icell] is open, then the value density[icell] is
!     not used. This wastes a bit of memory (and file size), but
!     otherwise it would be confusing which value of density[]
!     belongs to which cell.
!   - Closed cells with one or more cell walls facing the "vacuum"
!     surrounding space. These vacuum-facing cell walls together
!     form the "hull". This is the case for e.g. grids made from
!     a Delaunay triangulation. The advantage is that the grid has
!     a clearly defined outer hull. The disadvantage is that it is
!     not entirely trivial to find the point where a ray enters the
!     grid. If the hull is convex (such as for Delaunay), then
!     it can be done relatively non-expensively using projections.
!     But if the hull is oddly shaped, then the ray may enter, exit
!     and re-enter again, in which case the computation is a bit
!     more involved.
!
! RADMC-3D does not make the grid structure itself. You, the user,
! will have to provide the following information:
!
!   - A list of points that represent the interior (the "center")
!     of each cell. These will be the points where you will specify
!     the values of density, temperature etc.
!   - The volume of each cell, since RADMC-3D has no method to compute it.
!   - A list of cell walls, including support vector, normal vector
!     and the indices of the two adjacent cells. Note that the
!     extent or shape of these cell walls is not specified. Given
!     that each cell is (supposed to be) convex, for radiative
!     transfer we do not need to know.
!   - If your grid has a hull (i.e. cell walls facing "vacuum"), then:
!      - You must assure that the normal vector of each of the
!        cell walls of the hull points outward (into the "vacuum").
!      - Provide a list of indices of the cell walls of that hull.
!        If the hull is convex, then this is enough.
!      - If the hull is not convex, or if it may help with
!        faster ray-tracing, you can specify also:
!        - A list of vertices spanning the hull
!        - For each hull wall, the indices of the vertices spanning
!          it. 
!        - (if possible): for each face, a list of adjacent faces,
!          allowing to speed up the search to the entry point of
!          the next ray.
!
! RADMC-3D will then contruct, for each cell, a list of cell walls
! that enclose each cell. When conducting the radiative transfer,
! it will simply find the nearest crossing of the ray with the cell
! walls, and thus also find which next cell it will go into. It
! does not need to know the cell wall surface areas (which would be
! necessary, for instance, for hydrodynamics). 
!
! RADMC-3D has no way to validate the grid. For instance, if the
! points are not inside the cell walls, it wouldn't know. Also
! if the cell walls are not sufficient to close the cell, again
! it would not know. Also the volume of the cell it cannot verify.
! But the radiative transfer works, as long as you provide the
! self-consistent grid information. 
!==========================================================================
module ugrid_module

  integer :: ugrid_ncells=0              ! Total number of cells
  integer :: ugrid_nwalls=0              ! Total number of cell walls
  integer :: ugrid_nverts=0              ! Total number of vertices (if stored)
  integer :: ugrid_ncells_open=0         ! Number of open (unbounded) cells
  integer :: ugrid_hull_nwalls=0         ! Number of hull walls
  integer :: ugrid_cell_max_nr_walls=0   ! Max number of walls per cell
  integer :: ugrid_cell_max_nr_verts=0   ! Max number of vertices per cell (if stored)
  integer :: ugrid_wall_max_nr_verts=0   ! Max number of vertices per wall
  integer :: ugrid_vert_max_nr_cells=0   ! Max number of cells per vertex (if stored)
  integer :: ugrid_find_cell_start=1     ! Finding cells can be costly. This is the starting cell index, which you can set at will
  integer :: ugrid_max_ray_length=1000000  ! The max nr of steps along a ray
  
  logical :: ugrid_hull_convex=.true.    ! Is the hull (if present) convex?
  
  integer, allocatable :: ugrid_wall_icells(:,:)           ! For each wall: the indices of the two adjacent cells
  integer, allocatable :: ugrid_wall_iverts(:,:)           ! For each wall: the indices of the vertices (if stored)
  integer, allocatable :: ugrid_wall_nverts(:)             ! For each wall: the number of vertices (if stored)
  integer, allocatable :: ugrid_hull_iwalls(:)             ! The indices of the hull walls (if present)
  integer, allocatable :: ugrid_cell_iopen(:)              ! The indices of the open (unbounded) cells (if present)
  integer, allocatable :: ugrid_cell_iwalls(:,:)           ! For each cell: the indices of the cell walls
  integer, allocatable :: ugrid_cell_sgnwalls(:,:)         ! For each cell: the orientation of the cell walls (is n pointing outward = +1 or inward = -1)
  integer, allocatable :: ugrid_cell_iverts(:,:)           ! For each cell: the indices of the vertices (if stored)
  integer, allocatable :: ugrid_cell_nverts(:)             ! For each cell: the number of vertices per cell (if stored)
  integer, allocatable :: ugrid_cell_nwalls(:)             ! For each cell: the number of walls
  integer, allocatable :: ugrid_vert_icells(:,:)           ! For each vert: the indices of the cells (if stored)
  integer, allocatable :: ugrid_vert_ncells(:)             ! For each vert: the number of cells per vertex (if stored)
  double precision, allocatable :: ugrid_cell_volume(:)    ! For each cell: the volume
  double precision, allocatable :: ugrid_cell_size(:)      ! For each cell: the "typical" linear size (useful for estimates)
  double precision, allocatable :: ugrid_wall_s(:,:)       ! For each wall: the support vector
  double precision, allocatable :: ugrid_wall_n(:,:)       ! For each wall: the normal vector
  double precision, allocatable :: ugrid_vertices(:,:)     ! The 3D locations of the vertices (if stored)
  double precision, allocatable :: ugrid_cellcenters(:,:)  ! The 3D locations of the cell centers (if stored)
  double precision, allocatable :: ugrid_bary_matinv(:,:,:)! The inverse matrices used for barycentric interpolation (only for tetraeder grid)

  double precision :: ugrid_reltol = 1d-10      ! The relative tolerance of out of cell / out of wall / negative ds, compared to cell size
  double precision :: ugrid_center_x,ugrid_center_y,ugrid_center_z,ugrid_radius
  
  logical :: ugrid_analyze=.false.              ! Note: The analyze feature works only in serial (not openmp parallel) for now
  integer :: ugrid_cell_wall_cross_checks = 0
  integer :: ugrid_hull_wall_cross_checks = 0
  integer :: ugrid_pic_wall_cross_checks = 0

contains

  !------------------------------------------------------------------
  !                     Geometry subroutines
  !------------------------------------------------------------------
  
  subroutine ugrid_innerp_3d(v1,v2,res)
    implicit none
    double precision :: v1(1:3),v2(1:3),res
    res = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
  end subroutine ugrid_innerp_3d

  subroutine ugrid_crossp_3d(v1,v2,res)
    implicit none
    double precision :: v1(1:3),v2(1:3),res(1:3)
    res(1) = v1(2)*v2(3) - v1(3)*v2(2)
    res(2) = v1(3)*v2(1) - v1(1)*v2(3)
    res(3) = v1(1)*v2(2) - v1(2)*v2(1)
  end subroutine ugrid_crossp_3d

  subroutine ugrid_polygon_get_normal_3d(p,np)
    implicit none
    integer :: np
    double precision :: p(np,1:3)
    double precision :: va(1:3),vb(1:3),vn(1:3),l
    if(np.lt.3) stop 3884
    va(1:3) = p(2,:)-p(1,:)
    vb(1:3) = p(3,:)-p(1,:)
    call ugrid_crossp_3d(va,vb,vn)
    l       = sqrt(vn(1)**2+vn(2)**2+vn(3)**2)
    if(l.eq.0.d0) then
       write(*,*) 'Error: Polygon has equal points'
       stop 3885
    endif
    vn(1:3) = vn(1:3)/l
  end subroutine ugrid_polygon_get_normal_3d

  subroutine ugrid_point_in_convex_polygon_3d(p,np,n,v,res)
    implicit none
    integer :: np,i
    double precision :: p(np,1:3),n(1:3),v(1:3)
    double precision :: dv1(1:3),dv2(1:3),cp(1:3),sp,sc
    logical :: res
    if(np.lt.3) stop 3884
    res      = .true.
    dv1(1:3) = p(np,1:3)-v(1:3)
    dv2(1:3) = p(1,1:3)-v(1:3)
    call ugrid_crossp_3d(dv1,dv2,cp)
    call ugrid_innerp_3d(cp,n,sp)
    do i=1,np-1
       dv1(1:3) = p(i,1:3)-v(1:3)
       dv2(1:3) = p(i+1,1:3)-v(1:3)
       call ugrid_crossp_3d(dv1,dv2,cp)
       call ugrid_innerp_3d(cp,n,sc)
       if(sc*sp<0.d0) res = .false.
    enddo
  end subroutine ugrid_point_in_convex_polygon_3d
  
  subroutine ugrid_polygon_check_convex_3d(p,np,n,res)
    implicit none
    integer :: np,i
    double precision :: p(np,1:3),n(1:3)
    double precision :: dv1(1:3),dv2(1:3),cp(1:3),sp,sc
    logical :: res
    if(np.lt.3) stop 3884
    res = .true.
    if(np.eq.3) return
    dv1(1:3) = p(np,1:3)-p(np-1,1:3)
    dv2(1:3) = p(1,1:3)-p(np,1:3)
    call ugrid_crossp_3d(dv1,dv2,cp)
    call ugrid_innerp_3d(cp,n,sp)
    dv1(1:3) = p(1,1:3)-p(np,1:3)
    dv2(1:3) = p(2,1:3)-p(1,1:3)
    call ugrid_crossp_3d(dv1,dv2,cp)
    call ugrid_innerp_3d(cp,n,sc)
    if(sc*sp<0.d0) res = .false.
    do i=2,np-1
       dv1(1:3) = p(i,1:3)-p(i-1,1:3)
       dv2(1:3) = p(i+1,1:3)-p(i,1:3)
       call ugrid_crossp_3d(dv1,dv2,cp)
       call ugrid_innerp_3d(cp,n,sc)
       if(sc*sp<0.d0) res = .false.
    enddo
  end subroutine ugrid_polygon_check_convex_3d
    
  subroutine ugrid_polygon_check_inplane(p,np,n,reltol,res)
    implicit none
    integer :: np,i
    double precision :: p(np,1:3),n(1:3),reltol
    double precision :: dv(1:3),l,eps
    logical :: res
    if(np.lt.3) stop 3884
    res = .true.
    if(np.eq.3) return
    dv   = p(1,:)-p(np,:)
    l    = sqrt(dv(1)**2+dv(2)**2+dv(3)**2)
    if(l.eq.0.d0) then
       write(*,*) 'Error: Polygon has equal points'
       stop 3886
    endif
    eps  = abs(dv(1)*n(1)+dv(2)*n(2)+dv(3)*n(3))
    if(eps.gt.reltol*l) res=.false.
    do i=1,np-1
       dv   = p(i+1,:)-p(i,:)
       l    = sqrt(dv(1)**2+dv(2)**2+dv(3)**2)
       if(l.eq.0.d0) then
          write(*,*) 'Error: Polygon has equal points'
          stop 3886
       endif
       eps  = abs(dv(1)*n(1)+dv(2)*n(2)+dv(3)*n(3))
       if(eps.gt.reltol*l) res=.false.
    enddo
  end subroutine ugrid_polygon_check_inplane

  subroutine ugrid_ray_cross_plane(v,dir,p0,n,lam,vnew)
    implicit none
    double precision :: v(1:3),dir(1:3),p0(1:3),n(1:3)
    double precision :: ddn,dum(1:3),lam,vnew(1:3)
    ddn  = dir(1)*n(1)+dir(2)*n(2)+dir(3)*n(3)
    if(ddn.eq.0) then
       lam       = -1d99
       vnew(1:3) = -1d99
       return
    endif
    dum(1:3)  = p0(1:3)-v(1:3)
    lam       = (dum(1)*n(1)+dum(2)*n(2)+dum(3)*n(3))/ddn
    vnew(1:3) = v(1:3) + lam*dir(1:3)
  end subroutine ugrid_ray_cross_plane

  subroutine ugrid_ray_cross_convex_polygon(v,dir,p,np,n,lam0,lam,vnew,res)
    implicit none
    integer :: np
    double precision :: v(1:3),dir(1:3),p(np,1:3),n(1:3)
    double precision :: ddn,dum(1:3),lam0,lam,vnew(1:3),p0(1:3)
    logical :: res
    p0(1:3) = p(1,1:3)
    call ugrid_ray_cross_plane(v,dir,p0,n,lam,vnew)
    if(lam.lt.lam0) then
       res = .false.
       return
    endif
    call ugrid_point_in_convex_polygon_3d(p,np,n,vnew,res)
  end subroutine ugrid_ray_cross_convex_polygon

  subroutine ugrid_point_in_wall_plane(iwall,v,tol,res,dist)
    implicit none
    integer :: iwall
    double precision :: n(1:3),v(1:3),dv(1:3),tol,dist
    logical :: res
    n(1:3)  = ugrid_wall_n(iwall,1:3)
    dv(1:3) = ugrid_wall_s(iwall,1:3)-v(1:3)
    dist    = abs(dv(1)*n(1) + dv(2)*n(2) + dv(3)*n(3))
    if(dist.gt.tol) then
       res = .false.
    else
       res = .true.
    endif
  end subroutine ugrid_point_in_wall_plane
  
  subroutine ugrid_point_in_wall(iwall,v,res)
    implicit none
    integer :: i,iwall,nverts,i0,i1
    double precision :: n(1:3),v(1:3)
    double precision :: dv1(1:3),dv2(1:3),cp(1:3),sp,sc
    logical :: res
    if(.not.allocated(ugrid_wall_iverts)) then
       write(*,*) 'Error: Cannot check for point in wall if ugrid_wall_iverts not allocated.'
       stop 7298
    endif
    if(.not.allocated(ugrid_wall_nverts)) then
       write(*,*) 'Error: Cannot check for point in wall if ugrid_wall_nverts not allocated.'
       stop 7298
    endif
    n(1:3)   = ugrid_wall_n(iwall,1:3)
    res      = .true.
    nverts   = ugrid_wall_nverts(iwall)
    i0       = ugrid_wall_iverts(iwall,nverts)
    i1       = ugrid_wall_iverts(iwall,1)
    dv1(1:3) = ugrid_vertices(i0,1:3)-v(1:3)
    dv2(1:3) = ugrid_vertices(i1,1:3)-v(1:3)
    call ugrid_crossp_3d(dv1,dv2,cp)
    call ugrid_innerp_3d(cp,n,sp)
    do i=1,nverts-1
       i0       = ugrid_wall_iverts(iwall,i)
       i1       = ugrid_wall_iverts(iwall,i+1)
       dv1(1:3) = ugrid_vertices(i0,1:3)-v(1:3)
       dv2(1:3) = ugrid_vertices(i1,1:3)-v(1:3)
       call ugrid_crossp_3d(dv1,dv2,cp)
       call ugrid_innerp_3d(cp,n,sc)
       if(sc*sp<0.d0) res = .false.
    enddo
  end subroutine ugrid_point_in_wall
  
  subroutine ugrid_ray_cross_wall(v,dir,iwall,lam0,lam,vnew,res)
    implicit none
    integer :: np,iwall
    double precision :: v(1:3),dir(1:3),p(1:3),n(1:3),s(1:3)
    double precision :: ddn,dum(1:3),lam0,lam,vnew(1:3)
    logical :: res
    s(1:3) = ugrid_wall_s(iwall,1:3)
    n(1:3) = ugrid_wall_n(iwall,1:3)
    call ugrid_ray_cross_plane(v,dir,s,n,lam,vnew)
    if(lam.lt.lam0) then
       res = .false.
       return
    endif
    call ugrid_point_in_wall(iwall,vnew,res)
  end subroutine ugrid_ray_cross_wall

  subroutine ugrid_point_in_cell(v,icell,res,iwallin)
    implicit none
    double precision :: v(1:3)
    integer :: icell,iwall,iiwall,iiwallin
    double precision :: dv(1:3),n(1:3),inp
    logical :: res
    integer, optional :: iwallin
    if(present(iwallin)) then
       iiwallin = iwallin
    else
       iiwallin = -1
    endif
    res = .true.
    do iwall=1,ugrid_cell_nwalls(icell)
       if(ugrid_analyze) ugrid_pic_wall_cross_checks = ugrid_pic_wall_cross_checks + 1
       iiwall  = ugrid_cell_iwalls(icell,iwall)
       if(iiwall.ne.iiwallin) then
          n(1:3)  = ugrid_cell_sgnwalls(icell,iwall)*ugrid_wall_n(iiwall,1:3)
          dv(1:3) = ugrid_wall_s(iiwall,1:3)-v(1:3)
          inp     = dv(1)*n(1) + dv(2)*n(2) + dv(3)*n(3)
          if(inp.lt.0) res=.false.
       endif
    enddo
  end subroutine ugrid_point_in_cell

  subroutine ugrid_closest_cell_wall_distance(v,icell,res)
    implicit none
    double precision :: v(1:3)
    integer :: icell,iwall,iiwall
    double precision :: dv(1:3),n(1:3),inp,res
    res = 1d99
    do iwall=1,ugrid_cell_nwalls(icell)
       iiwall  = ugrid_cell_iwalls(icell,iwall)
       n(1:3)  = ugrid_cell_sgnwalls(icell,iwall)*ugrid_wall_n(iiwall,1:3)
       dv(1:3) = ugrid_wall_s(iiwall,1:3)-v(1:3)
       inp     = dv(1)*n(1) + dv(2)*n(2) + dv(3)*n(3)
       if(inp.lt.0) then
          res=-1.0    ! This should not happen because the point is outside of the cell!
       else
          if(inp.lt.res) then
             res = inp
          endif
       endif
    enddo
  end subroutine ugrid_closest_cell_wall_distance
    
  subroutine ugrid_check_cell_wall_normals()
    implicit none
    integer :: icell
    double precision :: p(1:3)
    logical :: res
    do icell=1,ugrid_ncells
       p(1:3) = ugrid_cellcenters(icell,1:3)
       call ugrid_point_in_cell(p,icell,res)
       if(.not.res) then
          write(*,*) 'The cell ',icell,' has inconsistent cell walls.'
          stop 2214
       endif
    enddo
  end subroutine ugrid_check_cell_wall_normals

  subroutine ugrid_check_cell_wall_supvec()
    implicit none
    integer :: iwall,ivert,icell
    double precision :: va(1:3),vb(1:3),vc(1:3),n(1:3),inp,tol,len
    if(allocated(ugrid_vertices)) then
       do iwall=1,ugrid_nwalls
          icell   = ugrid_wall_icells(iwall,1)
          if(icell.le.0) stop 4329
          if(ugrid_cell_size(icell).gt.0.d0) then
             tol     = ugrid_reltol * ugrid_cell_size(icell)
             ivert   = ugrid_wall_iverts(iwall,1)
             va(1:3) = ugrid_vertices(ivert,1:3)-ugrid_wall_s(iwall,1:3)
             ivert   = ugrid_wall_iverts(iwall,2)
             vb(1:3) = ugrid_vertices(ivert,1:3)-ugrid_wall_s(iwall,1:3)
             ivert   = ugrid_wall_iverts(iwall,3)
             vc(1:3) = ugrid_vertices(ivert,1:3)-ugrid_wall_s(iwall,1:3)
             n(1)    = va(2)*vb(3) - va(3)*vb(2)
             n(2)    = va(3)*vb(1) - va(1)*vb(3)
             n(3)    = va(1)*vb(2) - va(2)*vb(1)
             len     = n(1)**2 + n(2)**2 + n(3)**2
             n(1)    = n(1) / len
             n(2)    = n(2) / len
             n(3)    = n(3) / len
             inp     = abs(n(1)*vc(1)+n(2)*vc(2)+n(3)*vc(3))
             if(inp.gt.tol) then
                write(*,*) 'Error: Cell wall ',iwall,' has support vector that is not in '
                write(*,*) 'the plane spanned by its vertices. Inner product = ',inp
                write(*,*) 'The three vectors va=v_1-s, vb=v_2-s and vc=v_3-s are:'
                write(*,*) va
                write(*,*) vb
                write(*,*) vc
                write(*,*) 'The normal vector of the wall is:'
                write(*,*) n
                stop 5942
             endif
          endif
       enddo
    endif
  end subroutine ugrid_check_cell_wall_supvec

  subroutine ugrid_find_cell_wall_crossing(v,dir,icell,iwallcurr,iwallnext,icellnext,res)
    implicit none
    double precision :: v(1:3),dir(1:3)
    integer :: icell,iwall,iiwall,iwallcurr,iwallnext,icellnext
    double precision :: dum(1:3),n(1:3),res,ddn,dsvn,lam,tol,dist
    logical :: inwall
    res = 1d110
    iwallnext = -1
    icellnext = -1
    if(icell.le.0) stop 1253
    !
    ! If iwallcurr is set, then do some self-consistency checks
    !
    if(iwallcurr.gt.0) then
       !
       ! Check that the wall is indeed among the walls of the cell
       !
       do iwall=1,ugrid_cell_nwalls(icell)
          if(ugrid_cell_iwalls(icell,iwall).eq.iwallcurr) goto 10
       enddo
       write(*,*) 'ERROR in ugrid_find_cell_wall_crossing(): iwallcurr ',iwallcurr
       write(*,*) '      not cell wall of icell ',icell
10     continue
       !
       ! Let's check if we are indeed
       ! in the plane of this wall (only possible if cell size is
       ! given).
       !
       if(ugrid_cell_size(icell).gt.0.d0) then
          tol = ugrid_reltol * ugrid_cell_size(icell)
          call ugrid_point_in_wall_plane(iwallcurr,v,tol,inwall,dist)
          if(.not.inwall) then
             write(*,*) 'ERROR: Photon out of cell ',icell,'! Distance to cell wall ',iwallcurr,' = ',dist
             write(*,*) '       Further information: Cell volume        = ',ugrid_cell_volume(icell)
             write(*,*) '       Further information: Cell size estimate = ',ugrid_cell_size(icell)
             write(*,*) '       Tolerance = ',tol
             stop 8520
          endif
       endif
    endif
    !
    ! Loop over all walls of this cell
    !
    do iwall=1,ugrid_cell_nwalls(icell)
       if(ugrid_analyze) ugrid_cell_wall_cross_checks = ugrid_cell_wall_cross_checks + 1
       iiwall  = ugrid_cell_iwalls(icell,iwall)
       if(iiwall.ne.iwallcurr) then
          n(1:3)  = ugrid_cell_sgnwalls(icell,iwall)*ugrid_wall_n(iiwall,1:3)
          ddn     = dir(1)*n(1)+dir(2)*n(2)+dir(3)*n(3)
          if(ddn.gt.0) then
             dum(1:3) = ugrid_wall_s(iiwall,1:3)-v(1:3)
             dsvn     = dum(1)*n(1)+dum(2)*n(2)+dum(3)*n(3)
             lam      = dsvn/ddn    ! Note: Due to tiny round-off-errors this can be a tiny bit <0
             if(lam.lt.res) then
                res       = lam
                iwallnext = iiwall
             endif
          endif
       endif
    enddo
    !
    ! If a hit found, find the next cell index. If that is -1, then we
    ! will leave the grid. If the grid hull is convex, then this is a
    ! permanent leave. If not, then the ray *might* re-enter the grid.
    !
    if(iwallnext.gt.0) then
       if(ugrid_wall_icells(iwallnext,1).eq.icell) then
          icellnext = ugrid_wall_icells(iwallnext,2)
       else
          icellnext = ugrid_wall_icells(iwallnext,1)
          if(ugrid_wall_icells(iwallnext,2).ne.icell) then
             write(*,*) 'Error: Cell wall inconsistency'
             stop 5472
          endif
       endif
    endif
  end subroutine ugrid_find_cell_wall_crossing
    
  subroutine ugrid_find_hull_wall_crossing(v,dir,iwallcurr,iwallnext,icellnext,res)
    implicit none
    double precision :: v(1:3),dir(1:3)
    integer :: iwall,iiwall,iwallnext,icellnext,icell,iwallcurr
    double precision :: dum(1:3),n(1:3),s(1:3),vnew(1:3),res,ddn,dsvn,lam,lam_front,lam_back
    double precision :: lam_front_max,lam_back_min,tol,dist
    logical :: hit
    logical :: inwall
    res = -1.d0
    lam_front_max = -1
    lam_back_min  = 1d99
    iwallnext = -1
    icellnext = -1
    if(ugrid_hull_nwalls.eq.0) then
       write(*,*) 'Error: hull walls not defined'
       stop 1491
    endif
    !
    ! If iwallcurr is set, then let's check if we are indeed
    ! in the plane of this wall
    !
    if(iwallcurr.gt.0) then
       icell = ugrid_wall_icells(iwallcurr,1)
       if(icell.le.0) stop 2343
       tol = ugrid_reltol * ugrid_cell_size(icell)
       call ugrid_point_in_wall_plane(iwallcurr,v,tol,inwall,dist)
       if(.not.inwall) then
          write(*,*) 'ERROR: Photon out of cell ',icell,'! Distance to hull wall ',iwallcurr,' = ',dist
          write(*,*) '       Cell size    = ',ugrid_cell_size(icell)
          write(*,*) '       Wall s and n = ',ugrid_wall_s(iwallcurr,:),ugrid_wall_n(iwallcurr,:)
          stop 8520
       endif
    endif
    !
    ! Switch between two algorithms
    !
    if(ugrid_hull_convex) then
       !
       ! Hull is convex, so use the convex algorithm
       !
       do iwall=1,ugrid_hull_nwalls
          if(ugrid_analyze) ugrid_hull_wall_cross_checks = ugrid_hull_wall_cross_checks + 1
          iiwall  = ugrid_hull_iwalls(iwall)
          if(iiwall.ne.iwallcurr) then
             s(1)    = ugrid_wall_s(iiwall,1)
             s(2)    = ugrid_wall_s(iiwall,2)
             s(3)    = ugrid_wall_s(iiwall,3)
             n(1)    = ugrid_wall_n(iiwall,1)
             n(2)    = ugrid_wall_n(iiwall,2)
             n(3)    = ugrid_wall_n(iiwall,3)
             call ugrid_ray_cross_plane(v,dir,s,n,lam,vnew)
             ddn     = dir(1)*n(1)+dir(2)*n(2)+dir(3)*n(3)
             if(ddn.lt.0) then
                !
                ! Front of convex hull, as seen from the ray origin
                !
                if(lam.gt.lam_front_max) then
                   lam_front_max = lam
                   iwallnext     = iiwall
                endif
             elseif(ddn.gt.0) then
                !
                ! Back of convex hull, as seen from the ray origin
                !
                if(lam.lt.lam_back_min) then
                   lam_back_min = lam
                endif
             endif
          endif
       enddo
       if(lam_back_min.lt.lam_front_max) then
          !
          ! Miss
          !
          res       = -1
          iwallnext = -1
       else
          res       = lam_front_max
          if(res.le.0.d0) then
             !
             ! Away
             !
             res       = -1
             iwallnext = -1
          endif
       endif
    else
       !
       ! Hull is not necessarily convex, so use the general algorithm
       !
       res = 1d99
       iwallnext = -1
       do iwall=1,ugrid_hull_nwalls
          iiwall  = ugrid_hull_iwalls(iwall)
          if(iiwall.ne.iwallcurr) then
             call ugrid_ray_cross_wall(v,dir,iiwall,0.d0,lam,vnew,hit)
             if(hit) then
                if(lam.lt.res) then
                   res = lam
                   iwallnext = iiwall
                endif
             endif
          endif
       enddo
       if(res.le.0.d0) iwallnext = -1
       if(iwallnext.le.0) res=-1
    endif
    if(iwallnext.gt.0) then
       icellnext = ugrid_wall_icells(iwallnext,1)
       if((icellnext.le.0).or.(ugrid_wall_icells(iwallnext,2).gt.0)) then
          write(*,*) 'Error: Cell wall inconsistency'
          stop 5473
       endif
    else
       icellnext = -1
    endif
  end subroutine ugrid_find_hull_wall_crossing

  subroutine ugrid_find_center_and_radius()
    implicit none
    integer :: icell
    double precision :: r
    ugrid_center_x = 0.d0
    ugrid_center_y = 0.d0
    ugrid_center_z = 0.d0
    ugrid_radius   = 0.d0
    do icell=1,ugrid_ncells
       ugrid_center_x = ugrid_center_x + ugrid_cellcenters(icell,1)
       ugrid_center_y = ugrid_center_y + ugrid_cellcenters(icell,2)
       ugrid_center_z = ugrid_center_z + ugrid_cellcenters(icell,3)
    enddo
    ugrid_center_x = ugrid_center_x / ugrid_ncells
    ugrid_center_y = ugrid_center_y / ugrid_ncells 
    ugrid_center_z = ugrid_center_z / ugrid_ncells 
    do icell=1,ugrid_ncells
       r =  (ugrid_cellcenters(icell,1)-ugrid_center_x)**2 + &
            (ugrid_cellcenters(icell,2)-ugrid_center_y)**2 + &
            (ugrid_cellcenters(icell,3)-ugrid_center_z)**2
       r = sqrt(r)
       if(r.gt.ugrid_radius) then
          ugrid_radius = r
       endif
    enddo
  end subroutine ugrid_find_center_and_radius
  
  !------------------------------------------------------------------
  !                    Ray tracing subroutines
  !------------------------------------------------------------------

  subroutine ugrid_findcell(ray_cart_x,ray_cart_y,ray_cart_z,cellindex)
    implicit none
    doubleprecision :: ray_cart_x,ray_cart_y,ray_cart_z
    integer :: cellindex,index
    double precision :: v(1:3)
    logical :: res,found
    !
    ! NOTE: This can be slow. Better use ugrid_findcell_by_walking()
    !
    ! Put position into vector
    !
    v(1)   = ray_cart_x
    v(2)   = ray_cart_y
    v(3)   = ray_cart_z
    !
    ! Loop over all cells
    !
    found = .false.
    cellindex = -1
    do index=1,ugrid_ncells
       call ugrid_point_in_cell(v,index,res)
       if(res) then
          if(found) then
             write(*,*) 'Error: Found point in multiple cells'
             stop 5711
          endif
          found = .true.
          cellindex = index
       endif
    enddo
  end subroutine ugrid_findcell

  subroutine ugrid_findopencell(ray_cart_x,ray_cart_y,ray_cart_z,cellindex)
    implicit none
    doubleprecision :: ray_cart_x,ray_cart_y,ray_cart_z
    integer :: cellindex,index,iindex
    double precision :: v(1:3)
    logical :: res,found
    !
    ! NOTE: This can be slow. Better use ugrid_findcell_by_walking()
    !
    ! If there are no open cells, then return -1
    !
    if(ugrid_ncells_open.eq.0) then
       cellindex = -1
       return
    endif
    !
    ! Put position into vector
    !
    v(1)   = ray_cart_x
    v(2)   = ray_cart_y
    v(3)   = ray_cart_z
    !
    ! Loop over all open cells
    !
    found = .false.
    cellindex = -1
    do index=1,ugrid_ncells_open
       iindex = ugrid_cell_iopen(index)
       call ugrid_point_in_cell(v,iindex,res)
       if(res) then
          if(found) then
             write(*,*) 'Error: Found point in multiple cells'
             stop 5711
          endif
          found = .true.
          cellindex = iindex
       endif
    enddo
  end subroutine ugrid_findopencell

  subroutine ugrid_findcell_by_walking(ray_cart_x,ray_cart_y,ray_cart_z,cellindex)
    implicit none
    doubleprecision :: ray_cart_x,ray_cart_y,ray_cart_z
    doubleprecision :: x,y,z,dirx,diry,dirz,ds,dsend
    !doubleprecision :: xold,yold,zold  ! For debugging
    double precision :: v(1:3)
    integer :: index,nindex,iwall,niwall,is,cellindex
    logical :: res,arrived
    !
    ! Start from a point where we know the cell index
    !
    if(ugrid_find_cell_start.le.0) ugrid_find_cell_start=1
    x      = ugrid_cellcenters(ugrid_find_cell_start,1)
    y      = ugrid_cellcenters(ugrid_find_cell_start,2)
    z      = ugrid_cellcenters(ugrid_find_cell_start,3)
    dirx   = ray_cart_x - x
    diry   = ray_cart_y - y
    dirz   = ray_cart_z - z
    dsend  = sqrt(dirx**2+diry**2+dirz**2)
    dirx   = dirx / dsend
    diry   = diry / dsend
    dirz   = dirz / dsend
    iwall  = -1
    index  = ugrid_find_cell_start
    !
    ! Walk to the point of interest
    !
    do is=1,ugrid_max_ray_length
       !
       ! Backup current position (for debugging)
       !
       !xold = x
       !yold = y 
       !zold = z
       !
       ! Take a step of the walk
       !
       call ugrid_find_next_location(dsend,x,y,z,dirx,diry,dirz,    &
            index,nindex,ds,iwall,niwall,arrived)
       !
       ! Check that the step is consistent (for debugging)
       !
       !if(index.gt.0) then
       !   v(1)   = 0.5*(x+xold)
       !   v(2)   = 0.5*(y+yold)
       !   v(3)   = 0.5*(z+zold)
       !   call ugrid_point_in_cell(v,index,res)
       !   if(.not.res) then
       !      write(*,*) 'Error in ugrid_findcell_by_walking(): Segment found is inconsistent.'
       !      write(*,*) 'is = ',is,', arrived = ',arrived,', dsend = ',dsend,', dir = ',dirx,diry,dirz
       !      write(*,*) 'index = ',index,', nindex = ',nindex,', iwall = ',iwall,', niwall = ',niwall
       !      write(*,*) 'Cell index = ',index,', Midpoint of segment = ',v(:)
       !      write(*,*) 'Cell center = ',ugrid_cellcenters(index,1:3)
       !      write(*,*) 'Starting point of segment = ',xold,yold,zold
       !      write(*,*) 'Ending point of segment   = ',x,y,z
       !      stop 3789
       !   endif
       !endif
       !
       ! Check that point is indeed on the surface of the next cell, if given,
       ! and if the ray has not stopped before that.
       !
       if((nindex.gt.0).and.(.not.arrived)) then
          v(1)   = x
          v(2)   = y
          v(3)   = z
          call ugrid_point_in_cell(v,nindex,res,niwall)
          if(.not.res) then
             write(*,*) 'Error in ugrid_findcell_by_walking(): Point not in next cell.'
             stop 3790
          endif
       endif
       !
       ! Next
       !
       dsend = sqrt((ray_cart_x-x)**2+(ray_cart_y-y)**2+(ray_cart_z-z)**2)
       if(arrived) goto 10
       index = nindex
       iwall = niwall
    enddo
    write(*,*) 'Error in ugrid_module.ugrid_findcell_by_walking(): max nr of steps along ray exceeded.'
    stop 1000
10  continue
    !
    ! Check if we are indeed in the cell
    !
    if(index.gt.0) then
       v(1)   = ray_cart_x
       v(2)   = ray_cart_y
       v(3)   = ray_cart_z
       call ugrid_point_in_cell(v,index,res)
       if(.not.res) then
          write(*,*) 'Error in ugrid_findcell_by_walking(): Cell found is inconsistent.'
          write(*,*) 'Cell index found = ',index,', Point = ',v(:)
          stop 3788
       endif
    endif
    cellindex = index
  end subroutine ugrid_findcell_by_walking

  subroutine ugrid_find_next_location(ray_dsend,                     &
                ray_cart_x,ray_cart_y,ray_cart_z,                    &
                ray_cart_dirx,ray_cart_diry,ray_cart_dirz,           &
                ray_index,ray_indexnext,ray_ds,ray_curr_iwall,       &
                ray_next_iwall,arrived)
    implicit none
    double precision :: ray_cart_x,ray_cart_y,ray_cart_z
    double precision :: ray_cart_dirx,ray_cart_diry,ray_cart_dirz
    double precision :: ray_dsend,ray_ds
    double precision :: v(1:3),dir(1:3)
    integer :: ray_index,ray_indexnext,ray_curr_iwall,ray_next_iwall
    logical :: arrived
    !
    ! Make sure we only set arrived if we really arrived
    !
    arrived = .false.
    if(ray_dsend.le.0.d0) then
       write(*,*) 'Error in ugrid_find_next_location(): ray_dsend.le.0.d0'
       stop 2572
    endif
    !
    ! Put things into vectors
    !
    v(1)   = ray_cart_x
    v(2)   = ray_cart_y
    v(3)   = ray_cart_z
    dir(1) = ray_cart_dirx
    dir(2) = ray_cart_diry
    dir(3) = ray_cart_dirz
    !
    ! Distinguish between starting outside or inside of grid
    !
    if(ray_index.le.0) then
       !
       ! Ray comes from outside the grid
       !
       call ugrid_find_hull_wall_crossing(v,dir,ray_curr_iwall,  &
            ray_next_iwall,ray_indexnext,ray_ds)
       !
       ! Check if arrived by missing target
       !
       if(ray_ds.lt.0) then
          ray_ds  = 1d99
          arrived = .true.
       endif
    else
       !
       ! Ray starts inside the grid
       !
       call ugrid_find_cell_wall_crossing(v,dir,ray_index,ray_curr_iwall,  &
            ray_next_iwall,ray_indexnext,ray_ds)
       !
       ! Check that if ray_ds is negative (which should actually never happen, but
       ! could happen due to small round-off errors), it is small and not indicative
       ! of a real error
       !
       if(ray_ds.lt.0.d0) then
          if(abs(ray_ds).gt.ugrid_reltol*ray_dsend) then
             write(*,*) 'ERROR in ugrid_find_next_location(): Found negative ds...'
             write(*,*) '  icell = ',ray_index,', ds = ',ray_ds,', dsend = ',ray_dsend
             write(*,*) '  v = ',v,', dir = ',dir,', iwall = ',ray_curr_iwall,', niwall = ',ray_next_iwall
             write(*,*) '  index_next = ',ray_indexnext
             stop 4893
          endif
       endif
       !
       ! Check if arrived
       !
       if(ray_ds.gt.ray_dsend) then
          ray_ds  = ray_dsend
          arrived = .true.
       endif
    endif
    !
    ! Update position
    !
    ray_cart_x = ray_cart_x + ray_ds * ray_cart_dirx
    ray_cart_y = ray_cart_y + ray_ds * ray_cart_diry
    ray_cart_z = ray_cart_z + ray_ds * ray_cart_dirz
    !
  end subroutine ugrid_find_next_location
  
  !------------------------------------------------------------------
  !              Subbox and interpolation subroutines
  !------------------------------------------------------------------

  subroutine ugrid_create_list_cells_per_vertex()
    implicit none
    integer :: icell,ivert,maxncells,iivert
    if((ugrid_nverts.le.0).or.(.not.allocated(ugrid_vertices))) then
       write(*,*) 'ERROR: If you want to create the list of cells per vertex'
       write(*,*) '       you must first read the vertices.'
       stop 4453
    endif
    if(.not.allocated(ugrid_cell_iverts)) then
       write(*,*) 'ERROR: If you want to create the list of cells per vertex'
       write(*,*) '       you must first create the list of vertices per cell.'
       stop 4454
    endif
    allocate(ugrid_vert_ncells(ugrid_nverts))
    !
    ! First find the max nr of cells per vertex
    !
    do icell=1,ugrid_ncells
       do ivert=1,ugrid_cell_nverts(icell)
          iivert = ugrid_cell_iverts(icell,ivert)
          if(iivert.le.0) stop 7337
          ugrid_vert_ncells(iivert) = ugrid_vert_ncells(iivert) + 1
       enddo
    enddo
    maxncells = 0
    do ivert=1,ugrid_nverts
       if(ugrid_vert_ncells(ivert).gt.maxncells) then
          maxncells = ugrid_vert_ncells(ivert)
       endif
    enddo
    if(maxncells.ne.ugrid_vert_max_nr_cells) then
       write(*,*) 'NOTE: Counted ',maxncells,' max nr of cells per vertex. '
       write(*,*) 'File said: ',ugrid_vert_max_nr_cells,'. Taking ',maxncells
       ugrid_vert_max_nr_cells = maxncells
    endif
    !
    ! Now associate the cells to the vertices
    !
    allocate(ugrid_vert_icells(ugrid_nverts,ugrid_vert_max_nr_cells))
    ugrid_vert_ncells(:) = 0
    do icell=1,ugrid_ncells
       do ivert=1,ugrid_cell_nverts(icell)
          iivert = ugrid_cell_iverts(icell,ivert)
          if(iivert.le.0) stop 7337
          ugrid_vert_ncells(iivert) = ugrid_vert_ncells(iivert) + 1
          ugrid_vert_icells(iivert,ugrid_vert_ncells(iivert)) = icell
       enddo
    enddo
  end subroutine ugrid_create_list_cells_per_vertex

  subroutine ugrid_interpolate_from_cells_to_vertices(nv,ncells,nverts,funccell,funcvert)
    implicit none
    integer :: nv,ncells,nverts,icell,ivert,iicell
    double precision :: funccell(nv,ncells),funcvert(nv,nverts),dum(nv)
    if(.not.allocated(ugrid_vert_icells)) then
       write(*,*) 'ERROR: Can only interpolate to the vertices if '
       write(*,*) '       ugrid_vert_icells array is constructed.'
       stop 3347
    endif
    do ivert=1,ugrid_nverts
       dum(:) = 0.d0
       do icell=1,ugrid_vert_ncells(ivert)
          iicell = ugrid_vert_icells(ivert,icell)
          if(iicell.le.0) stop 4438
          dum(:) = dum(:) + funccell(:,iicell)
       enddo
       dum(:) = dum(:) / ugrid_vert_ncells(ivert)
       funcvert(:,ivert) = dum(:)
    enddo
  end subroutine ugrid_interpolate_from_cells_to_vertices
  
  subroutine ugrid_invert_matrix(m)
    implicit none
    double precision :: m(1:3,1:3),c(1:3,1:3),det
    !
    ! Compute the cofactor matrix
    !
    c(1,1) = m(2,2)*m(3,3)-m(2,3)*m(3,2)
    c(1,2) = m(2,3)*m(3,1)-m(2,1)*m(3,3)
    c(1,3) = m(2,1)*m(3,2)-m(2,2)*m(3,1)
    c(2,1) = m(3,2)*m(1,3)-m(3,3)*m(1,2)
    c(2,2) = m(3,3)*m(1,1)-m(3,1)*m(1,3)
    c(2,3) = m(3,1)*m(1,2)-m(3,2)*m(1,1)
    c(3,1) = m(1,2)*m(2,3)-m(1,3)*m(2,2)
    c(3,2) = m(1,3)*m(2,1)-m(1,1)*m(2,3)
    c(3,3) = m(1,1)*m(2,2)-m(1,2)*m(2,1)
    !
    ! Compute the determinant
    !
    det    = m(1,1)*c(1,1)+m(1,2)*c(1,2)+m(1,3)*c(1,3)
    !det    = m(2,1)*c(2,1)+m(2,2)*c(2,2)+m(2,3)*c(2,3)
    !det    = m(3,1)*c(3,1)+m(3,2)*c(3,2)+m(3,3)*c(3,3)
    if(det.eq.0.d0) then
       write(*,*) 'Zero determinant in 3x3 matrix.'
       stop 7452
    endif
    !
    ! Compute the inverse
    !
    m(1,1) = c(1,1)/det
    m(1,2) = c(2,1)/det
    m(1,3) = c(3,1)/det
    m(2,1) = c(1,2)/det
    m(2,2) = c(2,2)/det
    m(2,3) = c(3,2)/det
    m(3,1) = c(1,3)/det
    m(3,2) = c(2,3)/det
    m(3,3) = c(3,3)/det
  end subroutine ugrid_invert_matrix
  
  subroutine ugrid_setup_barycentric_matrices()
    implicit none
    integer :: icell,ivert,icorner
    double precision :: m(1:3,1:3),v4(1:3),d(1:3)
    if(ugrid_cell_max_nr_verts.ne.4) then
       write(*,*) 'ERROR: Barycentric interpolation: Only possible for tetraeder shaped cells!'
       stop 9345
    endif
    allocate(ugrid_bary_matinv(ugrid_ncells,1:3,1:3))
    do icell=1,ugrid_ncells
       ivert    = ugrid_cell_iverts(icell,4)
       if(ivert.le.0) then
          write(*,*) 'ERROR: In interpolation in ugrid: vertex missing.'
          stop 3687
       endif
       v4(1:3)  = ugrid_vertices(ivert,1:3)
       do icorner=1,3
          ivert          = ugrid_cell_iverts(icell,icorner)
          if(ivert.le.0) then
             write(*,*) 'ERROR: In interpolation in ugrid: vertex missing.'
             stop 3687
          endif
          m(1:3,icorner) = ugrid_vertices(ivert,1:3)-v4(1:3)
       enddo
       call ugrid_invert_matrix(m)
       do icorner=1,3
          ugrid_bary_matinv(icell,icorner,1:3) = m(icorner,1:3)
       enddo
    enddo
  end subroutine ugrid_setup_barycentric_matrices
  
  subroutine ugrid_barycentric_coords(icell,v,lam,tol)
    implicit none
    integer :: icell,icorner,ivert
    double precision :: v(1:3),lam(1:4),matrix(1:3,1:3),dum(1:3),tol
    if(icell.lt.1) stop 3443
    if(.not.allocated(ugrid_bary_matinv)) then
       write(*,*) 'Error: For barycentric interpolation, need ugrid_bary_matinv array.'
       stop 8723
    endif
    ivert    = ugrid_cell_iverts(icell,4)
    dum(1:3) = v(1:3) - ugrid_vertices(ivert,1:3)
    do icorner=1,3
       lam(icorner) = ugrid_bary_matinv(icell,icorner,1)*dum(1) + &
                      ugrid_bary_matinv(icell,icorner,2)*dum(2) + &
                      ugrid_bary_matinv(icell,icorner,3)*dum(3)
    enddo
    lam(4) = 1.d0-lam(1)-lam(2)-lam(3)
    do icorner=1,4
       if(lam(icorner).lt.-tol) then
          write(*,*) 'ERROR in ugrid_barycentric_coords(): interpolation out of cell.'
          write(*,*) 'icell = ',icell
          write(*,*) 'lam   = ',lam(1),lam(2),lam(3),lam(4)
          stop 6639
       endif
    enddo
  end subroutine ugrid_barycentric_coords

  subroutine ugrid_interpol_vectorfield_from_vertices(nv,nverts,x,y,z,cellindex,func,res)
    implicit none
    integer :: nv,nverts,iv
    double precision x,y,z
    double precision :: func(nv,nverts),res(nv),dummy,lam(1:4),tol,v(1:3)
    integer :: cellindex,icorner,ivert
    if(nverts.ne.ugrid_nverts) then
       write(*,*) 'ERROR: When using ugrid_interpol_from_vertices() the '
       write(*,*) '       physical values must be given at the vertices.'
       if(nverts.ne.ugrid_ncells) then
          write(*,*) '       Not at the cell centers!'
       endif
       write(*,*) '       The length of the array is unequal to number of vertices.'
       stop 5652
    endif
    if(.not.allocated(ugrid_cell_iverts)) then
       write(*,*) 'ERROR: When using ugrid_interpol_from_vertices() the '
       write(*,*) '       ugrid_cell_iverts array must be allocated.'
       stop 5652
    endif
    if(.not.allocated(ugrid_vertices)) then
       write(*,*) 'ERROR: When using ugrid_interpol_from_vertices() the '
       write(*,*) '       ugrid_vertices array must be allocated.'
       stop 5652
    endif
    res(1:nv) = 0.d0
    if(cellindex.le.0) then
       call ugrid_findcell_by_walking(x,y,z,cellindex)
    endif
    if(cellindex.gt.0) then
       tol = ugrid_reltol * ugrid_cell_size(cellindex)
       if(tol.eq.0.d0) then
          write(*,*) 'In cell ',cellindex,' cell_size is zero.'
          stop 3699
       endif
       v(1) = x
       v(2) = y
       v(3) = z
       call ugrid_barycentric_coords(cellindex,v,lam,tol)
       do iv=1,nv
          dummy = 0.d0
          do icorner=1,4
             ivert = ugrid_cell_iverts(cellindex,icorner)
             if(ivert.le.0) then
                write(*,*) 'ERROR: In interpolation in ugrid: vertex missing.'
                stop 3687
             endif
             dummy = dummy + lam(icorner)*func(iv,ivert)
          enddo
          res(iv) = dummy
       enddo
    endif
  end subroutine ugrid_interpol_vectorfield_from_vertices

  subroutine ugrid_interpol_tensorfield_from_vertices(nv1,nv2,nverts,x,y,z,cellindex,func,res)
    implicit none
    integer :: nv1,nv2,nverts,iv1,iv2
    double precision x,y,z
    double precision :: func(nv1,nv2,nverts),res(nv1,nv2),dummy,lam(1:4),tol,v(1:3)
    integer :: cellindex,icorner,ivert
    if(nverts.ne.ugrid_nverts) then
       write(*,*) 'ERROR: When using ugrid_interpol_from_vertices() the '
       write(*,*) '       physical values must be given at the vertices.'
       if(nverts.ne.ugrid_ncells) then
          write(*,*) '       Not at the cell centers!'
       endif
       write(*,*) '       The length of the array is unequal to number of vertices.'
       stop 5652
    endif
    if(.not.allocated(ugrid_cell_iverts)) then
       write(*,*) 'ERROR: When using ugrid_interpol_from_vertices() the '
       write(*,*) '       ugrid_cell_iverts array must be allocated.'
       stop 5652
    endif
    if(.not.allocated(ugrid_vertices)) then
       write(*,*) 'ERROR: When using ugrid_interpol_from_vertices() the '
       write(*,*) '       ugrid_vertices array must be allocated.'
       stop 5652
    endif
    res(:,:) = 0.d0
    if(cellindex.le.0) then
       call ugrid_findcell_by_walking(x,y,z,cellindex)
    endif
    if(cellindex.gt.0) then
       tol = ugrid_reltol * ugrid_cell_size(cellindex)
       if(tol.eq.0.d0) then
          write(*,*) 'In cell ',cellindex,' cell_size is zero.'
          stop 3699
       endif
       v(1) = x
       v(2) = y
       v(3) = z
       call ugrid_barycentric_coords(cellindex,v,lam,tol)
       do iv1=1,nv1
          do iv2=1,nv2
             dummy = 0.d0
             do icorner=1,4
                ivert = ugrid_cell_iverts(cellindex,icorner)
                if(ivert.le.0) then
                   write(*,*) 'ERROR: In interpolation in ugrid: vertex missing.'
                   stop 3687
                endif
                dummy = dummy + lam(icorner)*func(iv1,iv2,ivert)
             enddo
             res(iv1,iv2) = dummy
          enddo
       enddo
    endif
  end subroutine ugrid_interpol_tensorfield_from_vertices

  subroutine ugrid_interpol_scalar_from_vertices(nverts,x,y,z,cellindex,func,res)
    implicit none
    integer :: nverts
    double precision x,y,z
    double precision :: func(nverts),res,dummy,lam(1:4),tol,v(1:3)
    integer :: cellindex,icorner,ivert
    if(nverts.ne.ugrid_nverts) then
       write(*,*) 'ERROR: When using ugrid_interpol_from_vertices() the '
       write(*,*) '       physical values must be given at the vertices.'
       if(nverts.ne.ugrid_ncells) then
          write(*,*) '       Not at the cell centers!'
       endif
       write(*,*) '       The length of the array is unequal to number of vertices.'
       stop 5652
    endif
    if(.not.allocated(ugrid_cell_iverts)) then
       write(*,*) 'ERROR: When using ugrid_interpol_from_vertices() the '
       write(*,*) '       ugrid_cell_iverts array must be allocated.'
       stop 5652
    endif
    if(.not.allocated(ugrid_vertices)) then
       write(*,*) 'ERROR: When using ugrid_interpol_from_vertices() the '
       write(*,*) '       ugrid_vertices array must be allocated.'
       stop 5652
    endif
    res = 0.d0
    if(cellindex.le.0) then
       call ugrid_findcell_by_walking(x,y,z,cellindex)
    endif
    if(cellindex.gt.0) then
       tol = ugrid_reltol * ugrid_cell_size(cellindex)
       if(tol.eq.0.d0) then
          write(*,*) 'In cell ',cellindex,' cell_size is zero.'
          stop 3699
       endif
       v(1) = x
       v(2) = y
       v(3) = z
       call ugrid_barycentric_coords(cellindex,v,lam,tol)
       dummy = 0.d0
       do icorner=1,4
          ivert = ugrid_cell_iverts(cellindex,icorner)
          if(ivert.le.0) then
             write(*,*) 'ERROR: In interpolation in ugrid: vertex missing.'
             stop 3687
          endif
          dummy = dummy + lam(icorner)*func(ivert)
       enddo
       res = dummy
    endif
  end subroutine ugrid_interpol_scalar_from_vertices

  subroutine ugrid_subbox(nv,nc,nx,ny,nz,x0,x1,y0,y1,z0,z1,  &
                          phi1,theta,phi2,func,funcslice,    &
                          interpol)
    implicit none
    integer :: nv,nc,nx,ny,nz
    double precision :: x0,x1,y0,y1,z0,z1,x,y,z
    double precision :: func(nv,nc),funcslice(nx,ny,nz,nv)
    double precision :: xbk,ybk,zbk
    double precision :: sinphi1,cosphi1,sinphi2,cosphi2,sintheta,costheta
    double precision :: phi1,theta,phi2
    double precision :: r,th,ph,pi,res(1:nv)
    parameter(pi=3.14159265358979323846264338328d0)
    integer :: ix,iy,iz,index,iv
    logical :: interpol
    !
    ! Check if nc is correct
    !
    if(interpol) then
       if(.not.allocated(ugrid_vertices)) then
          write(*,*) 'If second order in ugrid_subbox(): Need to have ugrid_vertices allocated.'
          stop 6385
       endif
       if(nc.ne.ugrid_nverts) then
          write(*,*) 'In ugrid_subbox(): nc is not equal to ugrid_nverts.'
          stop 6386
       endif
       if(ugrid_wall_max_nr_verts.ne.4) then
          write(*,*) 'Sorry: In unstructured grids, interpolation only available for tetrad cells.'
          stop 6387
       endif
    else
       if(nc.ne.ugrid_ncells) then
          write(*,*) 'In ugrid_subbox(): nc is not equal to ugrid_ncells.'
          stop 6386
       endif
    endif
    !
    ! Rotation trigonometry of the box
    !
    sinphi1  = sin(phi1*pi/180.)
    cosphi1  = cos(phi1*pi/180.)
    sintheta = sin(theta*pi/180.)
    costheta = cos(theta*pi/180.)
    sinphi2  = sin(phi2*pi/180.)
    cosphi2  = cos(phi2*pi/180.)
    !
    ! Loop over all subbox grid points
    !
    do iz=1,nz
       do iy=1,ny
          do ix=1,nx
             !
             ! Get the (x,y,z) position
             !
             x = x0 + (x1-x0)*(ix-0.5d0)/(nx*1.d0)
             y = y0 + (y1-y0)*(iy-0.5d0)/(ny*1.d0)
             z = z0 + (z1-z0)*(iz-0.5d0)/(nz*1.d0)
             !
             ! Now make an optional rotation
             !
             xbk = x
             ybk = y
             x   = cosphi1 * xbk + sinphi1 * ybk
             y   =-sinphi1 * xbk + cosphi1 * ybk
             ybk = y
             zbk = z
             y   = costheta * ybk - sintheta * zbk
             z   = sintheta * ybk + costheta * zbk
             xbk = x
             ybk = y
             x   = cosphi2 * xbk + sinphi2 * ybk
             y   =-sinphi2 * xbk + cosphi2 * ybk
             !
             ! Now get the values from the grid
             !
             if(interpol) then
                !
                ! Interpolate in the ugrid
                !
                call ugrid_interpol_vectorfield_from_vertices(nv,ugrid_nverts,x,y,z,-1,func,res)
             else
                !
                ! Find the cell
                !
                call ugrid_findcell_by_walking(x,y,z,index)
                !
                ! If the cell is indeed a cell, then take that for the
                ! next search
                !
                if(index.gt.0) then
                   ugrid_find_cell_start = index
                endif
                !
                ! Take from the cell
                !
                if(index.ge.1) then
                   res(:) = func(:,index)
                else
                   res(:) = 0.d0
                endif
             endif
             !
             ! Put into the slice
             !
             funcslice(ix,iy,iz,:) = res(:)
          enddo
       enddo
    enddo
  end subroutine ugrid_subbox
  
  subroutine ugrid_subbox2(nv1,nv2,nc,nx,ny,nz,x0,x1,y0,y1,z0,z1,  &
                          phi1,theta,phi2,func,funcslice,          &
                          interpol)
    implicit none
    integer :: nv1,nv2,nc,nx,ny,nz
    double precision :: x0,x1,y0,y1,z0,z1,x,y,z
    double precision :: func(nv1,nv2,nc),funcslice(nx,ny,nz,nv1,nv2)
    double precision :: xbk,ybk,zbk
    double precision :: sinphi1,cosphi1,sinphi2,cosphi2,sintheta,costheta
    double precision :: phi1,theta,phi2
    double precision :: r,th,ph,pi,res(1:nv1,1:nv2)
    parameter(pi=3.14159265358979323846264338328d0)
    integer :: ix,iy,iz,index,iv
    logical :: interpol
    !
    ! Check if nc is correct
    !
    if(interpol) then
       if(.not.allocated(ugrid_vertices)) then
          write(*,*) 'If second order in ugrid_subbox(): Need to have ugrid_vertices allocated.'
          stop 6385
       endif
       if(nc.ne.ugrid_nverts) then
          write(*,*) 'In ugrid_subbox(): nc is not equal to ugrid_nverts.'
          stop 6386
       endif
       if(ugrid_wall_max_nr_verts.ne.4) then
          write(*,*) 'Sorry: In unstructured grids, interpolation only available for tetrad cells.'
          stop 6387
       endif
    else
       if(nc.ne.ugrid_ncells) then
          write(*,*) 'In ugrid_subbox(): nc is not equal to ugrid_ncells.'
          stop 6386
       endif
    endif
    !
    ! Rotation trigonometry of the box
    !
    sinphi1  = sin(phi1*pi/180.)
    cosphi1  = cos(phi1*pi/180.)
    sintheta = sin(theta*pi/180.)
    costheta = cos(theta*pi/180.)
    sinphi2  = sin(phi2*pi/180.)
    cosphi2  = cos(phi2*pi/180.)
    !
    ! Loop over all subbox grid points
    !
    do iz=1,nz
       do iy=1,ny
          do ix=1,nx
             !
             ! Get the (x,y,z) position
             !
             x = x0 + (x1-x0)*(ix-0.5d0)/(nx*1.d0)
             y = y0 + (y1-y0)*(iy-0.5d0)/(ny*1.d0)
             z = z0 + (z1-z0)*(iz-0.5d0)/(nz*1.d0)
             !
             ! Now make an optional rotation
             !
             xbk = x
             ybk = y
             x   = cosphi1 * xbk + sinphi1 * ybk
             y   =-sinphi1 * xbk + cosphi1 * ybk
             ybk = y
             zbk = z
             y   = costheta * ybk - sintheta * zbk
             z   = sintheta * ybk + costheta * zbk
             xbk = x
             ybk = y
             x   = cosphi2 * xbk + sinphi2 * ybk
             y   =-sinphi2 * xbk + cosphi2 * ybk
             !
             ! Now get the values from the grid
             !
             if(interpol) then
                !
                ! Interpolate in the ugrid
                !
                call ugrid_interpol_tensorfield_from_vertices(nv1,nv2,ugrid_nverts,x,y,z,-1,func,res)
             else
                !
                ! Find the cell
                !
                call ugrid_findcell_by_walking(x,y,z,index)
                !
                ! Take from the cell
                !
                if(index.ge.1) then
                   res(:,:) = func(:,:,index)
                else
                   res(:,:) = 0.d0
                endif
             endif
             !
             ! Put into the slice
             !
             funcslice(ix,iy,iz,:,:) = res(:,:)
          enddo
       enddo
    enddo
  end subroutine ugrid_subbox2
  
  !------------------------------------------------------------------
  !                     Reading the grid file
  !------------------------------------------------------------------

  subroutine ugrid_read_grid(success)
    implicit none
    integer :: icell,iwall,ivert
    logical :: fex1,fex3,success
    integer :: isavcen,isavvert,isavsize,iconv,isavvol,isavsn
    integer :: idum,iformat,idum1,idum2,counthullw,countopens
    integer :: maxnwalls,iiwall,iivert,maxnverts,ivscan,ivcount
    integer(kind=8) :: iiformat,nn,kk,precis
    double precision :: vec1(1:3),vec2(1:3),vec3(1:3),inp,len,cp(1:3),dp
    integer, allocatable :: iverts(:)
    double precision, allocatable :: vertices(:,:)
    integer :: iv1,iv2,iv3,iv4
    logical :: gotit,res
    call ugrid_cleanup()
    counthullw = 0
    countopens = 0
    maxnwalls  = 0
    maxnverts  = 0
    isavvol    = 1
    isavsn     = 1
    !
    ! Check which file exists
    !
    inquire(file='unstr_grid.inp',exist=fex1)
    inquire(file='unstr_grid.binp',exist=fex3)
    idum=0
    if(fex1) idum=idum+1
    if(fex3) idum=idum+1
    if(idum.gt.1) then
       write(*,*) 'ERROR: Found more than one file unstr_grid.*inp'
       stop
    endif
    if(idum.eq.0) then
       success = .false.
       return
    endif
    !
    ! Now read
    !
    if(fex1) then
       !
       ! The ascii version of the unstr_grid.inp file
       !
       open(unit=1,file='unstr_grid.inp',status='old')
       read(1,*) iformat
       read(1,*) ugrid_ncells
       read(1,*) ugrid_nwalls
       read(1,*) ugrid_nverts
       read(1,*) ugrid_cell_max_nr_walls
       read(1,*) ugrid_cell_max_nr_verts
       read(1,*) ugrid_wall_max_nr_verts
       read(1,*) ugrid_vert_max_nr_cells
       read(1,*) ugrid_ncells_open
       read(1,*) ugrid_hull_nwalls
       if(iformat.ge.2) then
          read(1,*) isavvol
          read(1,*) isavsn
       endif
       read(1,*) isavcen
       read(1,*) isavvert
       read(1,*) isavsize
       read(1,*) iconv
       if(iconv.eq.1) then
          ugrid_hull_convex=.true.
       else
          ugrid_hull_convex=.false.
       endif
       allocate(ugrid_wall_icells(ugrid_nwalls,2))
       allocate(ugrid_cell_volume(ugrid_ncells))
       allocate(ugrid_wall_s(ugrid_nwalls,3))
       allocate(ugrid_wall_n(ugrid_nwalls,3))
       allocate(ugrid_cellcenters(ugrid_ncells,3))
       if(isavvert.ne.0) then
          if(ugrid_wall_max_nr_verts.lt.3) then
             write(*,*) 'Error: When reading in the vertices, ugrid_wall_max_nr_verts must be set.'
             stop 7456
          endif
          allocate(ugrid_wall_nverts(ugrid_nwalls))
          allocate(ugrid_wall_iverts(ugrid_nwalls,ugrid_wall_max_nr_verts))
          allocate(ugrid_vertices(ugrid_nverts,3))
       endif
       if(isavvol.gt.0) then
          do icell=1,ugrid_ncells
             read(1,*) ugrid_cell_volume(icell)
          enddo
       endif
       if(isavsn.gt.0) then
          do iwall=1,ugrid_nwalls
             read(1,*) vec1(1:3),vec2(1:3)
             ugrid_wall_s(iwall,1:3) = vec1(1:3)
             ugrid_wall_n(iwall,1:3) = vec2(1:3)
          enddo
       endif
       do iwall=1,ugrid_nwalls
          read(1,*) idum1,idum2
          if(idum1.le.0) idum1=-1
          if(idum2.le.0) idum2=-1
          ugrid_wall_icells(iwall,1) = idum1
          ugrid_wall_icells(iwall,2) = idum2
       enddo
       if(isavcen.ne.0) then
          do icell=1,ugrid_ncells
             read(1,*) vec1(1:3)
             ugrid_cellcenters(icell,1:3) = vec1(1:3)
          enddo
       endif
       if(isavvert.ne.0) then
          ugrid_wall_nverts(:) = 0
          do iwall=1,ugrid_nwalls
             read(1,*) (ugrid_wall_iverts(iwall,ivert),ivert=1,ugrid_wall_max_nr_verts)
             do ivert=1,ugrid_wall_max_nr_verts
                if(ugrid_wall_iverts(iwall,ivert).gt.0) then
                   ugrid_wall_nverts(iwall) = ivert
                endif
             enddo
          enddo
          do ivert=1,ugrid_nverts
             read(1,*) vec1(1:3)
             ugrid_vertices(ivert,1:3) = vec1(1:3)
          enddo
       endif
       if(isavsize.ne.0) then
          allocate(ugrid_cell_size(ugrid_ncells))
          do icell=1,ugrid_ncells
             read(1,*) ugrid_cell_size(icell)
          enddo
       endif
       close(1)
    else
       !
       ! The binary version of the unstr_grid.binp file
       !
       open(unit=1,file='unstr_grid.binp',status='old',access='stream')
       read(1) iiformat
       read(1) precis
       if(precis.ne.8) stop 9341
       read(1) nn
       ugrid_ncells=nn
       read(1) nn
       ugrid_nwalls=nn
       read(1) nn
       ugrid_nverts=nn
       read(1) nn
       ugrid_cell_max_nr_walls=nn
       read(1) nn
       ugrid_cell_max_nr_verts=nn
       read(1) nn
       ugrid_wall_max_nr_verts=nn
       read(1) nn
       ugrid_vert_max_nr_cells=nn
       read(1) nn
       ugrid_ncells_open=nn
       read(1) nn
       ugrid_hull_nwalls=nn
       if(iiformat.ge.2) then
          read(1) nn
          isavvol=nn
          read(1) nn
          isavsn=nn
       endif
       read(1) nn
       isavcen=nn
       read(1)nn
       isavvert=nn
       read(1)nn
       isavsize=nn
       read(1) nn
       iconv=nn
       if(iconv.eq.1) then
          ugrid_hull_convex=.true.
       else
          ugrid_hull_convex=.false.
       endif
       allocate(ugrid_wall_icells(ugrid_nwalls,2))
       allocate(ugrid_cell_volume(ugrid_ncells))
       allocate(ugrid_wall_s(ugrid_nwalls,3))
       allocate(ugrid_wall_n(ugrid_nwalls,3))
       allocate(ugrid_cellcenters(ugrid_ncells,3))
       if(isavvert.ne.0) then
          if(ugrid_wall_max_nr_verts.lt.3) then
             write(*,*) 'Error: When reading in the vertices, ugrid_wall_max_nr_verts must be set.'
             stop 7456
          endif
          allocate(ugrid_wall_nverts(ugrid_nwalls))
          allocate(ugrid_wall_iverts(ugrid_nwalls,ugrid_wall_max_nr_verts))
          allocate(ugrid_vertices(ugrid_nverts,3))
       endif
       if(isavvol.gt.0) then
          do icell=1,ugrid_ncells
             read(1) ugrid_cell_volume(icell)
          enddo
       endif
       if(isavsn.gt.0) then
          do iwall=1,ugrid_nwalls
             read(1) vec1(1:3),vec2(1:3)
             ugrid_wall_s(iwall,1:3) = vec1(1:3)
             ugrid_wall_n(iwall,1:3) = vec2(1:3)
          enddo
       endif
       do iwall=1,ugrid_nwalls
          read(1) nn,kk
          idum1=nn
          idum2=kk
          if(idum1.le.0) idum1=-1
          if(idum2.le.0) idum2=-1
          ugrid_wall_icells(iwall,1) = idum1
          ugrid_wall_icells(iwall,2) = idum2
       enddo
       if(isavcen.ne.0) then
          do icell=1,ugrid_ncells
             read(1) vec1(1:3)
             ugrid_cellcenters(icell,1:3) = vec1(1:3)
          enddo
       endif
       if(isavvert.ne.0) then
          ugrid_wall_nverts(:) = 0
          do iwall=1,ugrid_nwalls
             do ivert=1,ugrid_wall_max_nr_verts
                read(1) nn
                ugrid_wall_iverts(iwall,ivert)=nn
             enddo
             ugrid_wall_nverts(iwall) = 0
             do ivert=1,ugrid_wall_max_nr_verts
                if(ugrid_wall_iverts(iwall,ivert).gt.0) then
                   if(ugrid_wall_nverts(iwall).ne.ivert-1) then
                      write(*,*) 'Error: In list of vertices for each wall '
                      write(*,*) '       all valid vertices must be listed contiguously.'
                      write(*,*) '       This appears not to be the case for iwall = ',iwall
                      stop 7344
                   endif
                   ugrid_wall_nverts(iwall) = ivert
                endif
             enddo
          enddo
          do ivert=1,ugrid_nverts
             read(1) vec1(1:3)
             ugrid_vertices(ivert,1:3) = vec1(1:3)
          enddo
       endif
       if(isavsize.ne.0) then
          allocate(ugrid_cell_size(ugrid_ncells))
          do icell=1,ugrid_ncells
             read(1) ugrid_cell_size(icell)
          enddo
       endif
       close(1)
    endif
    !
    ! Post-processing
    !
    ! If the s and n vectors are not read in, then try to
    ! construct them from the information we have
    !
    if(isavsn.eq.0) then
       if(isavvert.ne.0) then
          !
          ! We have the vertices of the cell walls, so we can construct
          ! s and n from them
          !
          allocate(vertices(ugrid_wall_max_nr_verts,1:3))
          do iwall=1,ugrid_nwalls
             !
             ! Create the normal vector from the first three vertices
             !
             iv1 = ugrid_wall_iverts(iwall,1)
             iv2 = ugrid_wall_iverts(iwall,2)
             iv3 = ugrid_wall_iverts(iwall,3)
             if((iv1.lt.1).or.(iv1.gt.ugrid_nverts)) stop 8421
             if((iv2.lt.1).or.(iv2.gt.ugrid_nverts)) stop 8422
             if((iv3.lt.1).or.(iv3.gt.ugrid_nverts)) stop 8423
             vec1(1:3) = ugrid_vertices(iv1,:)-ugrid_vertices(iv3,:)
             vec2(1:3) = ugrid_vertices(iv2,:)-ugrid_vertices(iv3,:)
             call ugrid_crossp_3d(vec1,vec2,cp)
             len = sqrt(cp(1)**2+cp(2)**2+cp(3)**2)
             cp(1:3) = cp(1:3) / len
             !
             ! Check if the polygon is really a plane (only necessary for
             ! walls with more than 3 vertices)
             !
             if(ugrid_wall_nverts(iwall).gt.3) then
                do ivert=1,ugrid_wall_nverts(iwall)
                   vertices(ivert,1:3) = ugrid_vertices(ugrid_wall_iverts(iwall,ivert),1:3)
                enddo
                call ugrid_polygon_check_inplane(vertices,ugrid_wall_nverts(iwall),cp,1d-12,res)
                if(.not.res) then
                   write(*,*) 'ERROR: Vertices of wall ',iwall,' do not lie exactly in a plane.'
                   stop 8399
                endif
             endif
             !
             ! Insert the normal vector into the array
             !
             ugrid_wall_n(iwall,1:3) = cp(1:3)
             !
             ! Create the support vector from the mean of the vertex positions
             !
             vec1(1:3) = 0.d0
             do ivert=1,ugrid_wall_nverts(iwall)
                vec1(1:3) = vec1(1:3) + ugrid_vertices(ugrid_wall_iverts(iwall,ivert),1:3)
             enddo
             vec1(1:3) = vec1(1:3) / ugrid_wall_nverts(iwall)
             !
             ! Insert the support vector into the array
             !
             ugrid_wall_s(iwall,1:3) = vec1(1:3)
          enddo
          deallocate(vertices)
       else
          !
          ! We do not have the vertices of the cell walls, so we assume
          ! that the cell walls are like Voronoi cell walls: in the middle
          ! between each pair of points and perpendicular to the line
          ! connecting the points
          !
          do iwall=1,ugrid_nwalls
             idum1 = ugrid_wall_icells(iwall,1)
             idum2 = ugrid_wall_icells(iwall,2)
             if(idum1.le.0) then
                write(*,*) 'ERROR: First cell index of wall ',iwall,' is <0'
                stop 8301
             endif
             if(idum2.le.0) then
                write(*,*) 'ERROR: With Voronoi-type grids each cell wall must'
                write(*,*) '       have a cell on both sides.'
                stop 8302
             endif
             ugrid_wall_s(iwall,1:3) = 0.5d0*(ugrid_cellcenters(idum1,1:3)+ugrid_cellcenters(idum2,1:3))
             ugrid_wall_n(iwall,1:3) = ugrid_cellcenters(idum2,1:3)-ugrid_cellcenters(idum1,1:3)
             len = sqrt(ugrid_wall_n(iwall,1)**2 + ugrid_wall_n(iwall,2)**2 + ugrid_wall_n(iwall,3)**2)
             ugrid_wall_n(iwall,1:3) = ugrid_wall_n(iwall,1:3) / len
          enddo
       endif
    endif
    !
    ! Count cell walls in each cell
    !
    allocate(ugrid_cell_nwalls(ugrid_ncells))
    ugrid_cell_nwalls(:) = 0
    do iwall=1,ugrid_nwalls
       icell = ugrid_wall_icells(iwall,1)
       ugrid_cell_nwalls(icell) = ugrid_cell_nwalls(icell) + 1
       icell = ugrid_wall_icells(iwall,2)
       if(icell.gt.0) then
          ugrid_cell_nwalls(icell) = ugrid_cell_nwalls(icell) + 1
       endif
    enddo
    !
    ! Get the max of that
    !
    maxnwalls = 0
    do icell=1,ugrid_ncells
       if(ugrid_cell_nwalls(icell).gt.maxnwalls) then
          maxnwalls = ugrid_cell_nwalls(icell)
       endif
    enddo
    if(maxnwalls.ne.ugrid_cell_max_nr_walls) then
       write(*,*) 'NOTE: Counted ',maxnwalls,' max nr of walls per cell. File said: ',ugrid_cell_max_nr_walls,'. Taking ',maxnwalls
       ugrid_cell_max_nr_walls = maxnwalls
    endif
    !
    ! Now assign the cell walls to the cells
    !
    ugrid_cell_nwalls(:) = 0
    allocate(ugrid_cell_iwalls(ugrid_ncells,ugrid_cell_max_nr_walls))
    allocate(ugrid_cell_sgnwalls(ugrid_ncells,ugrid_cell_max_nr_walls))
    do iwall=1,ugrid_nwalls
       icell = ugrid_wall_icells(iwall,1)
       if(icell.le.0) then
          write(*,*) 'Error: of each wall, the first cell index must be a cell.'
          write(*,*) ' Only the second can be empty (pointing into vacuum).'
          stop 743
       endif
       ugrid_cell_nwalls(icell) = ugrid_cell_nwalls(icell) + 1
       ugrid_cell_iwalls(icell,ugrid_cell_nwalls(icell))   = iwall
       ugrid_cell_sgnwalls(icell,ugrid_cell_nwalls(icell)) = 1
       icell = ugrid_wall_icells(iwall,2)
       if(icell.gt.0) then
          ugrid_cell_nwalls(icell) = ugrid_cell_nwalls(icell) + 1
          ugrid_cell_iwalls(icell,ugrid_cell_nwalls(icell))   = iwall
          ugrid_cell_sgnwalls(icell,ugrid_cell_nwalls(icell)) = -1
       else
          counthullw = counthullw + 1
       endif
    enddo
    if(counthullw.ne.ugrid_hull_nwalls) then
       write(*,*) 'NOTE: Counted ',counthullw,' hull walls. File said: ',ugrid_hull_nwalls,'. Taking ',counthullw
       ugrid_hull_nwalls = counthullw
    endif
    !
    ! Find all the hull walls
    !
    if(ugrid_hull_nwalls.ne.0) then
       allocate(ugrid_hull_iwalls(ugrid_hull_nwalls))
       counthullw = 0
       do iwall=1,ugrid_hull_nwalls
          icell = ugrid_wall_icells(iwall,2)
          if(icell.le.0) then
             counthullw = counthullw + 1
             ugrid_hull_iwalls(counthullw) = iwall
          endif
       enddo
    endif
    !
    ! If the cell center positions have not been read from the file,
    ! then we make an estimate here. NOTE: This is the reason why
    ! it is advisable to put the support vector of each wall in the
    ! center of each wall, not one of the vertices.
    !
    if(isavcen.eq.0) then
       do icell=1,ugrid_ncells
          ugrid_cellcenters(icell,1:3) = 0.d0
          do iwall=1,ugrid_cell_nwalls(icell)
             iiwall = ugrid_cell_iwalls(icell,iwall)
             ugrid_cellcenters(icell,1:3) = ugrid_cellcenters(icell,1:3) + ugrid_wall_s(iiwall,1:3)
          enddo
          ugrid_cellcenters(icell,1:3) = ugrid_cellcenters(icell,1:3) / ugrid_cell_nwalls(icell)
       enddo
    endif
    !
    ! Check the direction of the normal vectors of the walls, and adjust if necessary
    !
    do iwall=1,ugrid_nwalls
       icell     = ugrid_wall_icells(iwall,1)
       vec1(1:3) = ugrid_cellcenters(icell,1:3)
       icell     = ugrid_wall_icells(iwall,2)
       if(icell.gt.0) then
          vec2(1:3) = ugrid_cellcenters(icell,1:3)   ! Wall inside
       else
          vec2(1:3) = ugrid_wall_s(iwall,1:3)        ! Wall on hull
       endif
       vec3(1:3) = vec2(1:3)-vec1(1:3)
       inp = vec3(1)*ugrid_wall_n(iwall,1) + vec3(2)*ugrid_wall_n(iwall,2) + vec3(3)*ugrid_wall_n(iwall,3)
       if(inp.lt.0) then
          ugrid_wall_n(iwall,1:3) = - ugrid_wall_n(iwall,1:3)
       endif
    enddo
    !
    ! Make sure all normal vectors have length 1
    !
    do iwall=1,ugrid_nwalls
       len = sqrt(ugrid_wall_n(iwall,1)**2 + ugrid_wall_n(iwall,2)**2 + ugrid_wall_n(iwall,3)**2)
       if((len.lt.0.5d0).or.(len.gt.2.d0)) then
          write(*,*) 'Normal vector of wall ',iwall,' is badly normalized: length = ',len
          stop 934
       endif
       ugrid_wall_n(iwall,1:3) = ugrid_wall_n(iwall,1:3) / len
    enddo
    !
    ! If the vertices are read, then link them also to the cells
    !
    if(isavvert.ne.0) then
       allocate(ugrid_cell_nverts(ugrid_ncells))
       !
       ! Count the max number of vertices per cell
       !
       ivcount = ugrid_cell_max_nr_walls*ugrid_wall_max_nr_verts
       if(ivcount.le.0) stop 9031
       allocate(iverts(ivcount))
       do icell=1,ugrid_ncells
          ivcount = 0
          do iwall=1,ugrid_cell_nwalls(icell)
             iiwall = ugrid_cell_iwalls(icell,iwall)
             do ivert=1,ugrid_wall_nverts(iiwall)
                iivert = ugrid_wall_iverts(iiwall,ivert)
                if(iivert.le.0) stop 3301
                gotit  = .false.
                if(ivcount.gt.0) then
                   do ivscan=1,ivcount
                      if(iivert.eq.iverts(ivscan)) gotit=.true.
                   enddo
                endif
                if(.not.gotit) then
                   ivcount = ivcount + 1
                   iverts(ivcount) = iivert
                endif
             enddo
          enddo
          ugrid_cell_nverts(icell) = ivcount
          if(ivcount.gt.maxnverts) then
             maxnverts = ivcount
          endif
       enddo
       if(maxnverts.ne.ugrid_cell_max_nr_verts) then
          write(*,*) 'NOTE: Counted ',maxnverts,' as the max nr of vertices per cell. '
          write(*,*) '      File said: ',ugrid_cell_max_nr_verts,'. Taking ',maxnverts
          ugrid_cell_max_nr_verts = maxnverts
       endif
       !
       ! Now insert the vertex indices into the cells
       !
       allocate(ugrid_cell_iverts(ugrid_ncells,ugrid_cell_max_nr_verts))
       do icell=1,ugrid_ncells
          ivcount = 0
          do iwall=1,ugrid_cell_nwalls(icell)
             iiwall = ugrid_cell_iwalls(icell,iwall)
             do ivert=1,ugrid_wall_nverts(iiwall)
                iivert = ugrid_wall_iverts(iiwall,ivert)
                if(iivert.le.0) stop 3301
                gotit  = .false.
                if(ivcount.gt.0) then
                   do ivscan=1,ivcount
                      if(iivert.eq.iverts(ivscan)) gotit=.true.
                   enddo
                endif
                if(.not.gotit) then
                   ivcount = ivcount + 1
                   iverts(ivcount) = iivert
                endif
             enddo
          enddo
          do ivert=1,ivcount
             ugrid_cell_iverts(icell,ivert) = iverts(ivert)
          enddo
       enddo
       deallocate(iverts)
    endif
    !
    ! If the cell volumes are not read in, then try to compute them
    !
    if(isavvol.eq.0) then
       if(isavvert.eq.0) then
          write(*,*) 'ERROR: Without vertices, we cannot compute the cell volumes'
          stop 3988
       endif
       !
       ! For now we can only compute the cell volumes for tetrahedral cells
       !
       if((ugrid_cell_max_nr_walls.ne.4).or.(ugrid_cell_max_nr_verts.ne.4)) then
          write(*,*) 'ERROR: For the moment we can only compute the cell volumes internally'
          write(*,*) '       for tetrahedral cells. Please provide these volumes in the '
          write(*,*) '       unstr_grid.*inp file.'
          stop 3988
       endif
       !
       ! Now compute the cell volumes of the tetrads
       !
       do icell=1,ugrid_ncells
          iv1 = ugrid_cell_iverts(icell,1)
          iv2 = ugrid_cell_iverts(icell,2)
          iv3 = ugrid_cell_iverts(icell,3)
          iv4 = ugrid_cell_iverts(icell,4)
          if((iv1.lt.1).or.(iv1.gt.ugrid_nverts)) stop 8321
          if((iv2.lt.1).or.(iv2.gt.ugrid_nverts)) stop 8322
          if((iv3.lt.1).or.(iv3.gt.ugrid_nverts)) stop 8323
          if((iv4.lt.1).or.(iv4.gt.ugrid_nverts)) stop 8324
          vec1(1:3) = ugrid_vertices(iv1,:)-ugrid_vertices(iv4,:)
          vec2(1:3) = ugrid_vertices(iv2,:)-ugrid_vertices(iv4,:)
          vec3(1:3) = ugrid_vertices(iv3,:)-ugrid_vertices(iv4,:)
          call ugrid_crossp_3d(vec2,vec3,cp)
          dp = vec1(1)*cp(1) + vec1(2)*cp(2) + vec1(3)*cp(3)
          ugrid_cell_volume(icell) = abs(dp)/6
       enddo
    endif
    !
    ! If the ugrid_cell_size is not read, then create.
    ! This is a typical length scale appropriate for each cell, which can be used to
    ! estimate if the image has to be refined, or if round-off errors are too large.
    !
    if(.not.allocated(ugrid_cell_size)) then
       allocate(ugrid_cell_size(ugrid_ncells))
       do icell=1,ugrid_ncells
          ugrid_cell_size(icell) = ugrid_cell_volume(icell)**0.3333333
       enddo
    endif
    !
    ! Count nr of open cells
    !
    do icell=1,ugrid_ncells
       if(ugrid_cell_volume(icell).le.0.d0) then
          countopens = countopens + 1
       endif
    enddo
    if(countopens.ne.ugrid_ncells_open) then
       write(*,*) 'NOTE: Counted ',countopens,' open (unbounded) cells. File said: ',ugrid_ncells_open,'. Taking ',countopens
       ugrid_ncells_open = countopens
    endif
    !
    ! Find all open cells
    !
    if(ugrid_ncells_open.ne.0) then
       allocate(ugrid_cell_iopen(ugrid_ncells_open))
       countopens = 0
       do icell=1,ugrid_ncells
          if(ugrid_cell_volume(icell).le.0.d0) then
             countopens = countopens + 1
             ugrid_cell_iopen(countopens) = icell
          endif
       enddo
    endif
    !
    ! Check consistency of all cell walls
    !
    call ugrid_check_cell_wall_normals()
    if(allocated(ugrid_vertices)) then
       call ugrid_check_cell_wall_supvec()
    endif
    !
    ! Compute the overal size and mean position
    !
    call ugrid_find_center_and_radius()
    !
    ! Signal
    !
    success = .true.
  end subroutine ugrid_read_grid

  
  subroutine ugrid_cleanup()
    implicit none
    ugrid_ncells=0
    ugrid_nwalls=0
    ugrid_nverts=0
    ugrid_ncells_open=0
    ugrid_hull_nwalls=0
    ugrid_cell_max_nr_walls=0
    ugrid_cell_max_nr_verts=0
    ugrid_wall_max_nr_verts=0
    ugrid_vert_max_nr_cells=0
    ugrid_hull_convex=.true.
    if(allocated(ugrid_wall_icells))  deallocate(ugrid_wall_icells)
    if(allocated(ugrid_wall_iverts))  deallocate(ugrid_wall_iverts)
    if(allocated(ugrid_wall_nverts))  deallocate(ugrid_wall_nverts)
    if(allocated(ugrid_hull_iwalls))  deallocate(ugrid_hull_iwalls)
    if(allocated(ugrid_cell_iopen) )  deallocate(ugrid_cell_iopen) 
    if(allocated(ugrid_cell_iwalls))  deallocate(ugrid_cell_iwalls)
    if(allocated(ugrid_cell_sgnwalls))deallocate(ugrid_cell_sgnwalls)
    if(allocated(ugrid_cell_iverts))  deallocate(ugrid_cell_iverts)
    if(allocated(ugrid_cell_nverts))  deallocate(ugrid_cell_nverts)
    if(allocated(ugrid_cell_nwalls))  deallocate(ugrid_cell_nwalls)
    if(allocated(ugrid_cell_volume))  deallocate(ugrid_cell_volume)
    if(allocated(ugrid_cell_size))    deallocate(ugrid_cell_size)
    if(allocated(ugrid_wall_s)     )  deallocate(ugrid_wall_s)     
    if(allocated(ugrid_wall_n)     )  deallocate(ugrid_wall_n)     
    if(allocated(ugrid_vert_ncells))  deallocate(ugrid_vert_ncells)
    if(allocated(ugrid_vert_icells))  deallocate(ugrid_vert_icells)
    if(allocated(ugrid_vertices)   )  deallocate(ugrid_vertices)   
    if(allocated(ugrid_cellcenters))  deallocate(ugrid_cellcenters)
    if(allocated(ugrid_bary_matinv))  deallocate(ugrid_bary_matinv)
  end subroutine ugrid_cleanup
  
end module ugrid_module

