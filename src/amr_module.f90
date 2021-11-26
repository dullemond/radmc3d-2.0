!==========================================================================
!                    Adaptive Mesh Refinement Module
!
! This is a simple adaptive mesh refinement module. Nothing fancy, but
! functional. It does not support ghost cells in each branch, nor does it do
! zone decomposition for MPI parallel computing. For such complexity we
! refer to e.g. the PARAMESH or CHOMBO library. The objective of the present
! AMR module is to provide a simple mesh routine which allows to refine
! anywhere at will by dividing a cell into 2x2x2 (for 3-D) subcells.  In the
! end this module delivers a base grid in which each grid cell can be either
! a true cell (leaf) or a tree of sub-cells. All the branches are linked via
! pointers (a pointer to the parent and pointers to the children) and
! neighbors are also linked via pointers.  Moreover one can let the module
! produce a 1-D array of pointers to all cells (leafs), so that referencing
! to cells can be done in a do-loop rather than requiring the necessity of
! recursion. A similar thing can be done to all branches (i.e. leafs AND
! nodes). 
!
! The objective of this module is to provide a meshing for e.g.  radiative
! transfer or for diffusion-type equations. It is not quite efficient,
! however, for the explicit solving of hydrodynamics equations, since
! the module is not adapted for parallelization. For that we refer to
! the PARAMESH library or CHOMBO library.
!
!--------------------------------------------------------------------------
! BASE GRID AND TREE STRUCTURED REFINEMENTS
! 
! The base grid (forest) with possibly hierarchically structured 
! refinements (tree). The end-members of the branches are called leafs.
! The leafs are the "true" cells! Each symbol (o or *) is a branch 
! (the amr_branch type). A branch can be a node (consisting of children)
! or a leaf (end-member = AMR grid cell).
!
!  Base grid (level=0)        o----------*----------o----------o----------o
!                                      /||\
!                                    / |  | \
!  First level refinement           *  o  o  o
!                                 /||\   
!                               / |  | \
!  Second level refinement     o  o  o  o
!
! Two types of branches:
!   o = leaf
!   * = node
!
! NOTE: Compared to previous versions of amr_module, the base level is now
!       0 instead of 1. This is in better comparison to the levels of the
!       layers, and is in many respects more consistent.
!
!--------------------------------------------------------------------------
! HOW THE REFEINEMENT WORKS IN 2-D (EXAMPLE)
!
!  |-----------|-----------|-----------------------|
!  |           |           |                       |
!  |           |           |                       |
!  |           |           |                       |
!  |-----|-----|-----------|                       |
!  |     |     |           |                       |
!  |-----|-----|           |                       |
!  |     |     |           |                       |
!  |-----|-----|-----------|-----------------------|
!  |                       |                       |
!  |                       |                       |
!  |                       |                       |
!  |                       |                       |
!  |                       |                       |
!  |                       |                       |
!  |                       |                       |
!  |-----------------------|-----------------------|
!
!--------------------------------------------------------------------------
! Author: C.P. Dullemond
! Date:   07.09.07
!==========================================================================
module amr_module
use constants_module
implicit none
public
integer nvar
parameter(nvar=1)      ! This tells how many data elements stored per cell
!
! ........................................................................
!                     Type definitions for tree
! ........................................................................
!
! Since F90,F95 does not support direct arrays of pointers, one must
! do this indirectly (sorry, things become more messy this way; the
! C-style pointers are more elegant, but alas).
!
type amr_branch_link
   type(amr_branch), pointer :: link
end type amr_branch_link
!
! AMR grid branch 
!
!   This is EITHER a node branch containing (in 3-D) 2x2x2 children OR an
!   endpoint branch (leaf). If it is a node branch, then each of these
!   subcells can be another node branch, or a leaf (=end-point). A leaf is
!   the final true cell we are interested in. The level tells how deep in
!   the refinement this branch is. If level==0 this means that it is directly
!   associated to a `cell' of the base grid. Any level>0 is deeper
!   refinement level. (NOTE: This has been changed in version 0.19 of
!   RADMC-3D, because before that it was level==1 which was the base 
!   grid level)
!
!   The logical 'leaf' tells if this branch is an end-point branch (=leaf)
!   or a node branch which is supposed to be sub-divided into further
!   branches. Of course, if leaf==.true., then the children must be
!   filled, or else the branch is not self-consistent.
!
!   In principle only leaf branches should contain true data. But in the way
!   things are done here, every branch has its own cell data. It is up to
!   the user to decide what to do with them. There is a command to make sure
!   these data represent the average of the current values of the children
!   branches, if you need that.
!
!   NOTE: xc and xi here stand for the x, y and z coordinate: 
!         xc(1) is for x-coordinate, xc(2) for y etc.
!
type amr_branch
   integer :: level,id,branchindex,leafindex
   logical :: leaf
!   doubleprecision :: xc(3),xi(2,3)
   type(amr_branch_link) :: neighbor(2,3)
   type(amr_branch_link), pointer :: child(:,:,:)
   type(amr_branch), pointer :: parent
   integer :: parent_slot(3),ixyzf(3)
end type amr_branch
!
! ........................................................................
!                 Define global module variables
! ........................................................................
!
! An integer specifying the coordinate system type. This is only for 
! read/writing the amr_grid.inp file, as information for the user. It is
! not used in this module. So the user can specify his/her own integer
! coding for coordinate systems. 
!
integer :: amr_coordsystem=-12345
!
! An integer defining the AMR grid style. NOTE: The main oct-tree is
! always the base tree. Rest is extra for the user.
!
!   amr_style  = 0         Simple regular grid (no AMR)
!   amr_style  = 1         Normal oct-tree AMR grid
!   amr_style  = 10        Patchwork tree, somewhat like CHOMBO.
!                          This is the "layers" mode.
!
integer :: amr_style=1
!
! A flag to tell if the AMR tree is allocated or not. For amr_style.ne.0
! this must always be .true.. For amr_style.eq.0 this can be .false.
!
logical :: amr_always_use_tree=.false.
logical :: amr_tree_present=.false.
!
! AMR base grid
!   We start out with a normal regularly spaced grid. Each cell can be
!   a leaf or a branch (the root of a tree) (both are represented by
!   the type amr_branch). In case of a normal simulation without AMR each
!   gridcell is a leaf.  AMR means that some of these gridcells are
!   refined. Each of these refined gridcells spawns its own AMR tree. One could
!   call this grid the 'forest' because it consists of multiple trees. There
!   exists ONE (!) single AMR base grid (i.e. a single forest). Cells in
!   this 'forest' may be void (neither node branch nor leaf) if for instance
!   the grid is not necessary everywhere.
!
integer :: amr_grid_nx,amr_grid_ny,amr_grid_nz
doubleprecision, allocatable :: amr_grid_xi(:,:)
type(amr_branch_link), allocatable :: amr_grid_branch(:,:,:)
!
! 1-D Versions of AMR fine grid 
!
doubleprecision, allocatable :: amr_finegrid_xi(:,:,:)
doubleprecision, allocatable :: amr_finegrid_xc(:,:,:)
integer :: amr_nxyzfmax 
!
! An array of branches and leafs
!   It is useful to set these pointers, so that if you want to do an 
!   operation on all leafs or branches, then a simple loop operation is
!   sufficient, rather than having to dig deep into the trees.
!
!   WARNING: These arrays should use "icell" (or any counter
!            going from 1 to amr_nrleafs) as counter. DO NOT
!            USE "index" to address these arrays. Example:
!              do icell=1,amr_nrleafs
!                 cell=>amr_theleafs(icell)%link
!                 index=amr_theleaf_index(icell)
!                 ...
!              enddo
!            is correct. It allows you to do an operation
!            on all cells. But:
!                 cell=>amr_theleafs(index)%link
!            is wrong! Use instead: 
!                 cell=>amr_index_to_leaf(index)%link
!            which is correct.
!
type(amr_branch_link), allocatable, public :: amr_thebranches(:)
type(amr_branch_link), allocatable, public :: amr_theleafs(:)
integer, allocatable, public :: amr_thebranch_ids(:)
integer, allocatable, public :: amr_thebranch_index(:)
integer, allocatable, public :: amr_theleaf_index(:)
!
! Some AMR global data
!
integer :: amr_levelmax
integer :: amr_leafcount,amr_branchcount
integer :: amr_nrleafs,amr_nrbranches
integer :: amr_last_branch_id
!
! Stuff for the external data storage (the indexing).
! The "index" of a cell is the index you should use
! to access any of your data arrays. Example:
!    index = a%leafindex
!    temperature = mytemperaturearray(index)
!
integer :: amr_nrbranches_max,amr_nrleafs_max
logical :: amr_use_index
integer :: amr_leafindex_freeindex
integer :: amr_branchindex_freeindex
type(amr_branch_link), allocatable, public :: amr_index_to_branch(:)
type(amr_branch_link), allocatable, public :: amr_index_to_leaf(:)
integer, allocatable :: amr_leafindex_holes(:)
integer, allocatable :: amr_branchindex_holes(:)
integer :: amr_leafindex_nrholes=0
integer :: amr_branchindex_nrholes=0
!
! Some general information about the AMR tree
!
!  ...Which dimensions are active?
integer :: amr_xdim,amr_ydim,amr_zdim,amr_xyzdim(1:3),amr_dim,amr_nrchildref
!
!  ...Are there cyclic boundary conditions?
logical :: amr_cyclic_xyz(1:3)
!
! Counter 
!
integer :: amr_octtree_count
!
! Some arrays for layers
!
integer :: amr_layers_nrtot,amr_layers_nrlevels
integer, allocatable :: amr_layers_nr_per_level(:)
integer, allocatable :: amr_layers_level(:)
integer, allocatable :: amr_layers_parent(:)
integer, allocatable :: amr_layers_nxyz(:,:)
integer, allocatable :: amr_layers_nnxyz(:,:)
integer, allocatable :: amr_layers_ixyz(:,:)
type amr_layer_maplink
   integer, pointer :: index(:,:,:)
end type amr_layer_maplink
type(amr_layer_maplink), allocatable :: amr_layers_map(:)
integer :: amr_layers_count_ilayer,amr_layers_count_ix
integer :: amr_layers_count_iy,amr_layers_count_iz
!
! Stuff for the piecewise linear stuff
!
integer :: amr_slope_limiter=1
!
! ........................................................................
!                STUFF FOR CORNER BASED ALGORITHMS
! ........................................................................
!
! A temporary variable set 
!
type(amr_branch_link), public :: amr_loc_corner_cells(1:2,1:2,1:2)
integer :: amr_loc_corner_celllevel(1:2,1:2,1:2)
!integer :: amr_loc_cell_is_corner
logical :: amr_loc_corner(1:2,1:2,1:2)
!
! Main array of vertices
!
integer :: amr_nr_vertices=0,amr_nr_vertices_max=0
type(amr_branch_link), allocatable :: amr_vertex_cells(:,:,:,:)
!integer, allocatable :: amr_vertex_cell_is_corner(:)
!
! Main array for linking cells to their corner vertices
!
integer, allocatable :: amr_cell_corners(:,:,:,:)
!
! Attila Juhasz
! vtk_cell_type(:) VTK cell type of each leaf
! vtk_nr_cells(14) Nr of cells (leafs) of each VTK cell type in the whole mesh
character*5, allocatable :: vtk_cell_type(:)
integer :: vtk_nr_cells(14)
contains
!
!--------------------------------------------------------------------------
!                     AMR Initialization routine
!
! This routine sets up the base grid, defines the arrays of pointers etc.
!
! ARGUMENTS:
!  incl_x,y,z     If true, then use this dimension, if false, then this
!                 dimension is unused. Note that switching off a dimension
!                 cannot be done with setting e.g. nz=1, because upon 
!                 refinement also the z-dimension will be refined. Only
!                 by setting nz=1 AND incl_z=.false. you can truly switch
!                 off a dimension. Note that in that case you must STILL
!                 set the zi(1) and zi(2) because even if a dimension is
!                 not used, it is still assumed to be there passively,
!                 and the zi(1) and zi(2) of each subcell is still put
!                 to these values. If you do not care about this dimension,
!                 then this does not have any effect on your use of the
!                 library, but it could be useful for instance if you
!                 do radiative transfer in an axially symmetric model in
!                 which nothing depends on the phi (=z) dimension, but
!                 photons CAN travel in phi (=z) direction nevertheless.
!                 In that particular example zi(1)=0.d0 and zi(2)=2*pi,
!                 or zi(1)=-pi and zi(2)=pi, whatever you like. Important
!                 is that you do deliver a zi(1:2) array to this routine
!                 even if you do not use the z-direction, because otherwise
!                 the code crashes with a core-dump.
!  nx,ny,nz       The size of the grid in x,y,z direction (always set >=1).
!                 If you set nx=1,ny=1,nz=1 but incl_x=incl_y=incl_z=.true.
!                 then your base grid is basically a single cell, but all
!                 dimensions are active. In AMR this is not strange: it 
!                 simply means that you do not wish to use the option of a 
!                 base grid, but instead you want to define the entire grid 
!                 using the AMR refinement. That is absolutely fine.
!                 NOTE: In the reading and writing subroutines below the
!                       nx.eq.0 is used to disable this dimension. This
!                       works only in the I/O (i.e. in the data files),
!                       but not in the subroutines of this library.
!  xi,yi,zi       The cell interfaces of the base grid 
!                 (always nx+1 resp. ny+1 resp. nz+1 elements; not that
!                 they are ALWAYS arrays with >=2 elements, even if one or
!                 more of the dimensions are switched off.
!  levelmax       The maximum level of refinement you allow (0=no refinement)
!  fillbasegrid   If .eq.1 then the base grid is automatically filled
!                 with leaf branches
!  nrbranchesmax  If nrbranchesmax is set then we allow external data storage.
!                 It means that the branch index is managed in such a
!                 way that free indices are reused. Example: create
!                 4 branches (1,2,3,4), delete branch 2 and 3, then 
!                 create another branch. This branch should have index
!                 2 or 3, not 5. Otherwise you will quickly run out of
!                 computer memory in your externally stored data. Note:
!                 externally stored data means that YOU as user create
!                 an array of your data. The AMR tree will then point
!                 to where in your array the corresponding data to each
!                 branch. This is what the index is for.
!  nrleafsmax     If this is set (must be .le. nrbranchesmax) in addition
!                 to nrbranchesmax, then we allow for less large upper limit
!                 on the nr of leafs than on the nr of branches. This is
!                 useful if the external data storage is very large and
!                 one must not waste memory, and if the user decides to
!                 only store external data for leafs, not for branches.
!                 The user can then allocate its data arrays with nrleafsmax
!                 as max size and be sure that the amr code will complain
!                 if it tends to overflow this.
! x,y,zcyclic     If set, then the neighbors of the cells at the edge
!                 are set to the cells at the other end (cyclic boundary
!                 condition). 
! always_amrtree  If set, then even for a regular grid the amr tree will
!                 be set up. This is how RADMC-3D worked until version
!                 0.27. But from version 0.28 by default the amr tree
!                 will only be set up if at least one cell of the base
!                 grid has child-cells (which, for amr_initialize, means
!                 if levelmax.gt.0). This saves a lot of memory.
!--------------------------------------------------------------------------
subroutine amr_initialize(incl_x,incl_y,incl_z,&
                          nx,ny,nz,xi,yi,zi, &
                          levelmax,fillbasegrid,&
                          nrbranchesmax,nrleafsmax,&
                          xcyclic,ycyclic,zcyclic,&
                          always_amrtree)
implicit none
integer, intent(in) :: nx,ny,nz,fillbasegrid
doubleprecision, intent(in) :: xi(1:nx+1),yi(1:ny+1),zi(1:nz+1)
integer, intent(in) :: levelmax
integer :: ierr,ix,iy,iz,nxyzmax,i,nxmax,nymax,nzmax
integer :: slot(3)
integer :: nnx,nny,nnz,ilevel
integer, optional :: nrbranchesmax,nrleafsmax
logical :: incl_x,incl_y,incl_z
logical, optional :: xcyclic,ycyclic,zcyclic
logical, optional :: always_amrtree
type(amr_branch), pointer :: a,b,c
doubleprecision :: dx,dy,dz
!
! First destroy any possible earlier AMR grid if present
!
if(allocated(amr_grid_branch)) then
   write(stdo,*) 'While initializing AMR: Earlier AMR grid is destroyed'
endif
call amr_cleanup(partial=.true.)
!
! Then decide whether to set up an AMR tree or not
!
if(levelmax.gt.0) then
   !
   ! We have refinement (or at least we wish to allow for refinement)
   ! 
   amr_tree_present = .true.
else
   !
   ! We have no refinement, i.e. we have a regular grid and we do not
   ! intend to refine lateron. We do not need an AMR tree. If, however,
   ! the called insists on having a tree, then so be it.
   !
   amr_tree_present = .false.
   if(present(always_amrtree)) then
      if(always_amrtree) then
         amr_tree_present = .true.
      endif
   endif
endif
!
! Check dimensionality of the problem
!
if(nx.le.0) then
   write(stdo,*) 'ERROR: nx must be >= 1'
   stop
endif
if(ny.le.0) then
   write(stdo,*) 'ERROR: ny must be >= 1'
   stop
endif
if(nz.le.0) then
   write(stdo,*) 'ERROR: nz must be >= 1'
   stop
endif
amr_dim = 0
amr_grid_nx  = nx
amr_grid_ny  = ny
amr_grid_nz  = nz
if(incl_x) then
   amr_xdim = 1
   amr_dim  = amr_dim +1
else
   if(nx.ne.1) then
      write(stdo,*) 'ERROR: If dimension x not used, nx must be 1'
      stop
   endif
   amr_xdim = 0
endif
if(incl_y) then
   amr_ydim = 1
   amr_dim  = amr_dim +1
else
   if(ny.ne.1) then
      write(stdo,*) 'ERROR: If dimension y not used, ny must be 1'
      stop
   endif
   amr_ydim = 0
endif
if(incl_z) then
   amr_zdim = 1
   amr_dim  = amr_dim +1
else
   if(nz.ne.1) then
      write(stdo,*) 'ERROR: If dimension z not used, nz must be 1'
      stop
   endif
   amr_zdim = 0
endif
if(amr_dim==0) then
   write(stdo,*) 'ERROR in AMR Module: Zero-dimensional problem...'
   stop 7120
endif
amr_nrchildref = (1+amr_xdim)*(1+amr_ydim)*(1+amr_zdim)
amr_xyzdim(1)  = amr_xdim
amr_xyzdim(2)  = amr_ydim
amr_xyzdim(3)  = amr_zdim
!
! Cyclic or non-cyclic boundaries
!
amr_cyclic_xyz(1) = .false.
amr_cyclic_xyz(2) = .false.
amr_cyclic_xyz(3) = .false.
if(present(xcyclic)) amr_cyclic_xyz(1) = xcyclic
if(present(ycyclic)) amr_cyclic_xyz(2) = ycyclic
if(present(zcyclic)) amr_cyclic_xyz(3) = zcyclic
!
! Initialize variables
!
amr_leafcount   = 0
amr_branchcount = 0
amr_nrleafs     = 0
amr_nrbranches  = 0
amr_levelmax    = levelmax
amr_last_branch_id = 0
!
! Initialize the indexing stuff (for external storage)
!
if(amr_tree_present) then
   !
   ! The AMR tree is present, so allocate the associated arrays
   !
   if(present(nrbranchesmax)) then
      amr_use_index = .true.
      amr_nrbranches_max = nrbranchesmax
      allocate(amr_index_to_branch(1:nrbranchesmax),STAT=ierr)
      if(ierr.ne.0) then
         write(stdo,*) 'ERROR in AMR Module: Could not allocate amr_index_to_*()'
         stop 271
      endif
      allocate(amr_index_to_leaf(1:nrbranchesmax),STAT=ierr)
      if(ierr.ne.0) then
         write(stdo,*) 'ERROR in AMR Module: Could not allocate amr_index_to_*()'
         stop 271
      endif
      allocate(amr_branchindex_holes(1:nrbranchesmax),STAT=ierr)
      if(ierr.ne.0) then
         write(stdo,*) 'ERROR in AMR Module: Could not allocate amr_branchindex_holes()'
         stop 271
      endif
      allocate(amr_leafindex_holes(1:nrbranchesmax),STAT=ierr)
      if(ierr.ne.0) then
         write(stdo,*) 'ERROR in AMR Module: Could not allocate amr_leafindex_holes()'
         stop 271
      endif
      amr_branchindex_nrholes = 0
      amr_leafindex_nrholes = 0
      do i=1,nrbranchesmax
         nullify(amr_index_to_branch(i)%link)
         nullify(amr_index_to_leaf(i)%link)
      enddo
      if(present(nrleafsmax)) then
         if(nrleafsmax.gt.nrbranchesmax) then
            write(stdo,*) 'ERROR in AMR Module: Cannot have more leafs than branches'
            stop 7308
         endif
         amr_nrleafs_max = nrleafsmax
      else
         amr_nrleafs_max = amr_nrbranches_max
      endif
   else
      if(present(nrleafsmax)) then
         write(stdo,*) 'ERROR in AMR Module: If you specify nrleafsmax, must also specify'
         write(stdo,*) '     nrbranchesmax (which must be .ge. nrleafsmax).'
         stop 7449
      endif
      amr_use_index = .false.
      amr_nrbranches_max = 0
   endif
else
   !
   ! The grid is regular, no AMR tree is present. Only need to set a few numbers
   !
   amr_use_index      = .true.
   amr_levelmax       = 0
   amr_nrleafs_max    = amr_grid_nx*amr_grid_ny*amr_grid_nz
   amr_nrbranches_max = amr_grid_nx*amr_grid_ny*amr_grid_nz
endif
!
! Allocate space for the base grid
! (only if the AMR tree is used)
!
if(amr_tree_present) then
   !
   ! Allocate the array
   ! 
   allocate(amr_grid_branch(1:amr_grid_nx,1:amr_grid_ny,1:amr_grid_nz),STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMR Module: Could not allocate base grid pointers'
      stop 201
   endif
   !
   ! Nullify all pointers
   !
   do iz=1,amr_grid_nz
      do iy=1,amr_grid_ny
         do ix=1,amr_grid_nx
            nullify(amr_grid_branch(ix,iy,iz)%link)
         enddo
      enddo
   enddo
endif
!
! Note that we take the largest of the three dimensions for the size of
! out coordinate grid. What is not used we do not use.
!
nxyzmax = max(amr_grid_nx,amr_grid_ny,amr_grid_nz)
allocate(amr_grid_xi(1:nxyzmax+1,1:3),STAT=ierr)
if(ierr.ne.0) then 
   write(stdo,*) 'ERROR in AMR Module: Could not allocate grid'
   stop 201
endif
!
! Copy base grid to module
!
do ix=1,nx+1
   amr_grid_xi(ix,1) = xi(ix)
enddo
do ix=2,nx+1
   if(xi(ix).le.xi(ix-1)) then
      write(stdo,*) 'ERROR in AMR Module: Non-monotonically increasing grid in xi'
      stop 7301
   endif
enddo
do iy=1,ny+1
   amr_grid_xi(iy,2) = yi(iy)
enddo
do iy=2,ny+1
   if(yi(iy).le.yi(iy-1)) then
      write(stdo,*) 'ERROR in AMR Module: Non-monotonically increasing grid in yi'
      write(stdo,*) iy,yi(iy),yi(iy-1)
      stop 7302
   endif
enddo
do iz=1,nz+1
   amr_grid_xi(iz,3) = zi(iz)
enddo
do iz=2,nz+1
   if(zi(iz).le.zi(iz-1)) then
      write(stdo,*) 'ERROR in AMR Module: Non-monotonically increasing grid in zi'
      stop 7303
   endif
enddo
!
! Allocate the fine grid (which replaces the xi and xc specs in each cell,
! to save memory)
!
nxmax = nx * (2**levelmax)
nymax = ny * (2**levelmax)
nzmax = nz * (2**levelmax)
amr_nxyzfmax = max(nxmax,nymax,nzmax)
allocate(amr_finegrid_xi(1:amr_nxyzfmax+1,1:3,0:levelmax),STAT=ierr) 
if(ierr.ne.0) then
   write(stdo,*) 'ERROR in AMR Module: Could not allocate amr_finegrid_xi...'
   stop
endif
allocate(amr_finegrid_xc(1:amr_nxyzfmax,1:3,0:levelmax),STAT=ierr) 
if(ierr.ne.0) then
   write(stdo,*) 'ERROR in AMR Module: Could not allocate amr_finegrid_xc...'
   stop
endif
!
! Now make this fine grid 
!
! ...First the base grid
!
amr_finegrid_xi(1:nx+1,1,0) = amr_grid_xi(1:nx+1,1)
amr_finegrid_xi(1:ny+1,2,0) = amr_grid_xi(1:ny+1,2)
amr_finegrid_xi(1:nz+1,3,0) = amr_grid_xi(1:nz+1,3)
do ix=1,nx 
   amr_finegrid_xc(ix,1,0)  = 0.5d0 * ( amr_finegrid_xi(ix,1,0) + &
                                        amr_finegrid_xi(ix+1,1,0) )
enddo
do iy=1,ny 
   amr_finegrid_xc(iy,2,0)  = 0.5d0 * ( amr_finegrid_xi(iy,2,0) + &
                                        amr_finegrid_xi(iy+1,2,0) )
enddo
do iz=1,nz 
   amr_finegrid_xc(iz,3,0)  = 0.5d0 * ( amr_finegrid_xi(iz,3,0) + &
                                        amr_finegrid_xi(iz+1,3,0) )
enddo
nnx = nx
nny = ny
nnz = nz
!
! ...Now a loop over all levels
!
if(amr_tree_present) then
   do ilevel=1,levelmax
      !
      ! First the interface grid
      !
      if(incl_x) then
         do ix=1,nnx+1
            amr_finegrid_xi(2*(ix-1)+1,1,ilevel) = amr_finegrid_xi(ix,1,ilevel-1) 
         enddo
         do ix=1,nnx
            dx  = ( amr_finegrid_xi(ix+1,1,ilevel-1) -           &
                    amr_finegrid_xi(ix,1,ilevel-1)  ) * 0.5d0
            amr_finegrid_xi(2*ix,1,ilevel) = amr_finegrid_xi(2*ix-1,1,ilevel) + dx
         enddo
         nnx = nnx*2
      else
         if(nnx.ne.1) stop 6091
         amr_finegrid_xi(1,1,ilevel) = amr_finegrid_xi(1,1,ilevel-1)
         amr_finegrid_xi(2,1,ilevel) = amr_finegrid_xi(2,1,ilevel-1)
      endif
      if(incl_y) then
         do iy=1,nny+1
            amr_finegrid_xi(2*(iy-1)+1,2,ilevel) = amr_finegrid_xi(iy,2,ilevel-1) 
         enddo
         do iy=1,nny
            dy  = ( amr_finegrid_xi(iy+1,2,ilevel-1) -           &
                    amr_finegrid_xi(iy,2,ilevel-1)  ) * 0.5d0
            amr_finegrid_xi(2*iy,2,ilevel) = amr_finegrid_xi(2*iy-1,2,ilevel) + dy
         enddo
         nny = nny*2
      else
         if(nny.ne.1) stop 6092
         amr_finegrid_xi(1,2,ilevel) = amr_finegrid_xi(1,2,ilevel-1)
         amr_finegrid_xi(2,2,ilevel) = amr_finegrid_xi(2,2,ilevel-1)
      endif
      if(incl_z) then
         do iz=1,nnz+1
            amr_finegrid_xi(2*(iz-1)+1,3,ilevel) = amr_finegrid_xi(iz,3,ilevel-1) 
         enddo
         do iz=1,nnz
            dz  = ( amr_finegrid_xi(iz+1,3,ilevel-1) -           &
                    amr_finegrid_xi(iz,3,ilevel-1)  ) * 0.5d0
            amr_finegrid_xi(2*iz,3,ilevel) = amr_finegrid_xi(2*iz-1,3,ilevel) + dz
         enddo
         nnz = nnz*2
      else
         if(nnz.ne.1) stop 6093
         amr_finegrid_xi(1,3,ilevel) = amr_finegrid_xi(1,3,ilevel-1)
         amr_finegrid_xi(2,3,ilevel) = amr_finegrid_xi(2,3,ilevel-1)
      endif
      !
      ! Then the center grid
      !
      do ix=1,nnx 
         amr_finegrid_xc(ix,1,ilevel)  = 0.5d0 * ( amr_finegrid_xi(ix,1,ilevel) + &
                                                   amr_finegrid_xi(ix+1,1,ilevel) )
      enddo
      do iy=1,nny 
         amr_finegrid_xc(iy,2,ilevel)  = 0.5d0 * ( amr_finegrid_xi(iy,2,ilevel) + &
                                                   amr_finegrid_xi(iy+1,2,ilevel) )
      enddo
      do iz=1,nnz 
         amr_finegrid_xc(iz,3,ilevel)  = 0.5d0 * ( amr_finegrid_xi(iz,3,ilevel) + &
                                                   amr_finegrid_xi(iz+1,3,ilevel) )
      enddo
   enddo
endif
!
! If requested we can fill all base grid cells with branches that
! already have an associated cell. In this way we then have already
! a normal regular grid. Then we only need to refine lateron.
!
if(amr_tree_present) then
   !
   ! If the AMR tree is present...
   !
   if(fillbasegrid.eq.1) then
      nullify(b)
      do iz=1,amr_grid_nz
         do iy=1,amr_grid_ny
            do ix=1,amr_grid_nx
               slot(1) = ix
               slot(2) = iy
               slot(3) = iz
               call amr_branch_construct(a,b,slot)
            enddo
         enddo
      enddo
   endif
else
   !
   ! For a regular grid...
   !
   if(fillbasegrid.ne.1) then
      write(stdo,*) 'INTERNAL ERROR: If you use a regular grid and you do not use'
      write(stdo,*) '       an AMR grid, then you must set fillbasegrid=1'
      write(stdo,*) '       in the call to amr_initialize().'
      stop
   endif
   !
   ! We only need to specify the amr_nrleafs and amr_nrbranches
   !
   amr_nrleafs        = amr_nrleafs_max
   amr_nrbranches     = amr_nrbranches_max
endif
!
! Some default value for the piecewise linear stuff (obsolete)
!
amr_slope_limiter = 1  
!
end subroutine amr_initialize


!--------------------------------------------------------------------------
!                          AMR Cleanup
!--------------------------------------------------------------------------
subroutine amr_cleanup(partial)
implicit none
integer :: ierr
integer :: ix,iy,iz,ilayer
type(amr_branch), pointer :: b
logical,optional :: partial
logical :: partial_cleanup
!
if(present(partial)) then
   partial_cleanup = partial
else
   partial_cleanup = .false.
endif
!
! Destroy the grid
!
if(allocated(amr_grid_xi)) then
   deallocate(amr_grid_xi,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMR Module: Could not deallocate grid'
      stop 201
   endif
endif
if(allocated(amr_finegrid_xi)) then
   deallocate(amr_finegrid_xi)
endif
if(allocated(amr_finegrid_xc)) then
   deallocate(amr_finegrid_xc)
endif
!
! If present, destroy the lists of leafs and branches
!
if(allocated(amr_thebranches)) then
   deallocate(amr_thebranches,STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in AMR Module: Could not deallocate list of branches'
      stop 201
   endif
endif
if(allocated(amr_thebranch_ids)) then
   deallocate(amr_thebranch_ids,STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in AMR Module: Could not deallocate list of branch IDs'
      stop 201
   endif
endif
if(allocated(amr_theleafs)) then
   deallocate(amr_theleafs,STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in AMR Module: Could not deallocate list of leafs'
      stop 201
   endif
endif
if(allocated(amr_theleaf_index)) then
   deallocate(amr_theleaf_index,STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in AMR Module: Could not deallocate list of leaf indexes'
      stop 201
   endif
endif
if(allocated(amr_thebranch_index)) then
   deallocate(amr_thebranch_index,STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in AMR Module: Could not deallocate list of branch indexes'
      stop 201
   endif
endif
!
! Now scan the grid for stuff to destroy and finally destroy the grid
!
if(allocated(amr_grid_branch)) then
   !
   ! For each slot in the grid, if there is a tree, then destroy
   !
   do iz=1,amr_grid_nz
      do iy=1,amr_grid_ny
         do ix=1,amr_grid_nx
            if(associated(amr_grid_branch(ix,iy,iz)%link)) then
               b => amr_grid_branch(ix,iy,iz)%link
               call amr_branch_destruct(b)
            endif
         enddo
      enddo
   enddo
   amr_grid_nx=0
   amr_grid_ny=0
   amr_grid_nz=0
   !
   ! Now deallocate the list of trees
   !
   deallocate(amr_grid_branch,STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMR Module: Could not deallocate base grid pointers'
      stop 201
   endif
endif
!
! Deallocate the index lists
!
if(allocated(amr_index_to_branch)) then
   deallocate(amr_index_to_branch,STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in AMR Module: Could not deallocate branchindex list'
      stop 201
   endif
endif
if(allocated(amr_index_to_leaf)) then
   deallocate(amr_index_to_leaf,STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in AMR Module: Could not deallocate leafindex list'
      stop 201
   endif
endif
!
! Set the free index stuff to zero
!
!   ******************************************************************
!   NOTE: These freeindex numbers are just a poor-mans solution to the
!         N^2 problem of constructing branches. Later I should invent
!         a better way to find free indices.
!   ******************************************************************
!
amr_leafindex_freeindex   = 1
amr_branchindex_freeindex = 1
!
! Cleanup the layers
!
if(.not.partial_cleanup) then 
   if(allocated(amr_layers_nr_per_level)) deallocate(amr_layers_nr_per_level)
   if(allocated(amr_layers_level)) deallocate(amr_layers_level)
   if(allocated(amr_layers_parent)) deallocate(amr_layers_parent)
   if(allocated(amr_layers_ixyz)) deallocate(amr_layers_ixyz)
   if(allocated(amr_layers_nxyz)) deallocate(amr_layers_nxyz)
   if(allocated(amr_layers_nnxyz)) deallocate(amr_layers_nnxyz)
   if(allocated(amr_layers_map)) then
      do ilayer=0,amr_layers_nrtot
         if(associated(amr_layers_map(ilayer)%index)) deallocate(amr_layers_map(ilayer)%index)
      enddo
      deallocate(amr_layers_map)
   endif
   amr_layers_nrtot = 0
   amr_layers_nrlevels = 0
endif
!
! Cleanup the corner vertex stuff 
!
amr_nr_vertices = 0
!if(allocated(amr_vertex_cell_is_corner)) deallocate(amr_vertex_cell_is_corner)
if(allocated(amr_vertex_cells)) deallocate(amr_vertex_cells)
if(allocated(amr_cell_corners)) deallocate(amr_cell_corners)
!
end subroutine amr_cleanup


!--------------------------------------------------------------------------
!            Helper routine for listing all branches and leafs
!--------------------------------------------------------------------------
recursive subroutine amr_compute_list_tree(a)
implicit none
type(amr_branch), pointer :: a
integer ix,iy,iz
!
! Update this branch count
!
amr_branchcount = amr_branchcount + 1
if(amr_branchcount.gt.amr_nrbranches) then
   write(stdo,*) 'INTERNAL ERROR in AMR: too many branches'
   stop 346
endif
amr_thebranches(amr_branchcount)%link => a
amr_thebranch_ids(amr_branchcount) = a%id
if(amr_use_index) then
   amr_thebranch_index(amr_branchcount) = a%branchindex
endif
!
! If branch, then call iteratively the children
! Do, incidently, an internal consistency check
!
if(a%leaf) then
   !
   ! It is a leaf, check if indeed no children, as must be.
   ! And if all is OK, then count the leaf
   !
   if(associated(a%child)) then
      write(stdo,*) 'INTERNAL ERROR in AMR Module: Both branch and leaf...'
      write(stdo,*) '    for branch ID=',a%id
      stop
   endif
   amr_leafcount = amr_leafcount + 1
   if(amr_leafcount.gt.amr_nrleafs) then
      write(stdo,*) 'INTERNAL ERROR in AMR Module: too many leafs'
      stop 345
   endif
   amr_theleafs(amr_leafcount)%link => a
   if(amr_use_index) then
      amr_theleaf_index(amr_leafcount) = a%leafindex
   endif
else
   !
   ! It is a branch. Call all children. Check incidently
   ! that indeed all slots are filled, as must be
   !
   if(.not.associated(a%child)) then
      write(stdo,*) 'INTERNAL ERROR in AMR Module: Branch without children'
      write(stdo,*) '    for branch ID=',a%id
      stop
   endif
   do iz=1,1+amr_zdim
      do iy=1,1+amr_ydim
         do ix=1,1+amr_xdim
            if(.not.associated(a%child(ix,iy,iz)%link)) then
               write(stdo,*) 'INTERNAL ERROR in AMR Module: Child of branch not associated...'
               stop 891
            endif
            call amr_compute_list_tree(a%child(ix,iy,iz)%link)
         enddo
      enddo
   enddo
endif
!
end subroutine amr_compute_list_tree


!--------------------------------------------------------------------------
!                     List all branches and leafs
!
! This routine computes a complete list of branches and another complete
! list of leafs. This can be useful if we need to do operations on each
! branch or leaf without having to go recursively.
!--------------------------------------------------------------------------
subroutine amr_compute_list_all
implicit none
integer ix,iy,iz,ierr
!
! Check if the AMR tree exists
!
if(.not.amr_tree_present) then
   write(stdo,*) 'INTERNAL ERROR: Trying to make list of all AMR cells, but'
   write(stdo,*) '      AMR tree is not present (presumably because you use'
   write(stdo,*) '      a regular grid, which, as of version 0.28, does no'
   write(stdo,*) '      longer use the AMR tree).'
   stop
endif
!
! If the lists already exist, then destruct
!
if(allocated(amr_thebranches)) deallocate(amr_thebranches)
if(allocated(amr_theleafs)) deallocate(amr_theleafs)
if(allocated(amr_thebranch_ids)) deallocate(amr_thebranch_ids)
if(allocated(amr_thebranch_index)) deallocate(amr_thebranch_index)
if(allocated(amr_theleaf_index)) deallocate(amr_theleaf_index)
!
! Reset counters
!
amr_leafcount   = 0
amr_branchcount = 0
!
! Allocate the amr_thebranches, amr_theleafs etc arrays
!
if(amr_nrbranches.le.0) stop 5641
if(amr_nrleafs.le.0) stop 5642
allocate(amr_thebranches(1:amr_nrbranches),STAT=ierr)
if(ierr.ne.0) then 
   write(stdo,*) 'ERROR in AMR Module: Could not allocate amr_thebranches'
   stop 201
endif
allocate(amr_thebranch_ids(1:amr_nrbranches),STAT=ierr)
if(ierr.ne.0) then 
   write(stdo,*) 'ERROR in AMR Module: Could not allocate amr_thebranch_ids'
   stop 201
endif
allocate(amr_theleafs(1:amr_nrleafs),STAT=ierr)
if(ierr.ne.0) then 
   write(stdo,*) 'ERROR in AMR Module: Could not allocate amr_theleafs'
   stop 201
endif
!
! If external data storage is used, then also allocate the index arrays
! that allow a direct access to the index without having to go via the
! branch/leaf data structure. 
!
if(amr_use_index) then
   allocate(amr_thebranch_index(1:amr_nrbranches),STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMR Module: Could not allocate amr_thebranch_index'
      stop 201
   endif
   allocate(amr_theleaf_index(1:amr_nrleafs),STAT=ierr)
   if(ierr.ne.0) then 
      write(stdo,*) 'ERROR in AMR Module: Could not allocate amr_theleaf_index'
      stop 201
   endif
endif
!
! Loop over grid
!
do iz=1,amr_grid_nz
   do iy=1,amr_grid_ny
      do ix=1,amr_grid_nx
         if(associated(amr_grid_branch(ix,iy,iz)%link)) then
            call amr_compute_list_tree(amr_grid_branch(ix,iy,iz)%link)
         endif
      enddo
   enddo
enddo
end subroutine amr_compute_list_all


!--------------------------------------------------------------------------
!                Track branch with a certain ID
!
! This subroutine is actually only useable for debugging purposes. It 
! allows you to search for a branch with a given ID number. It is very
! slow, so it is not the preferred way to link branches! Note that this
! routine only works if the list of branches is made.
!--------------------------------------------------------------------------
subroutine amr_track_branch_id(id,a)
implicit none
integer :: i,id
type(amr_branch), pointer :: a
do i=1,amr_nrbranches
   if(amr_thebranch_ids(i).eq.id) then
      a => amr_thebranches(i)%link
      return
   endif
enddo
nullify(a)
write(stdo,*) 'PROBLEM: Could not track branch with ID=',id
end subroutine amr_track_branch_id


!--------------------------------------------------------------------------
!                         Find neighbor
!
! This routine will find the neighbor of branch a in dimension idir
! and direction isgn (+1 = torward right; -1 = toward left). The result
! is returned as a pointer neigh. If no neighbor exists, then a 
! nullified pointer is returned.
!--------------------------------------------------------------------------
subroutine amr_find_neighbor_branch(a,idir,isgn,neigh)
implicit none
type(amr_branch), pointer :: a
type(amr_branch), pointer :: neigh
type(amr_branch), pointer :: b
integer :: idir,isgn,ip,ipo,level,cmax,i
integer :: slotchain(3,amr_levelmax+2),slot(3),chaincount
logical :: repeat
!
! Check
!
if(.not.associated(a)) then
   write(stdo,*) 'ERROR in AMR Module: Pointer to branch not associated'
   stop 3452
endif
!
nullify(neigh)
chaincount = 0
level = a%level
if(isgn.gt.0) then
   ip  = 1
   ipo = 2
else
   ip  = 2
   ipo = 1
endif
!
! Search the tree
!
b => a
repeat = .true.
do while(repeat)
   if(.not.associated(b)) then
      write(stdo,*) 'ERROR in AMR Module: Pointer to branch not associated'
      stop 3453
   endif
   if(associated(b%parent)) then
      !
      ! Parent is still within tree
      !
      if(ip==b%parent_slot(idir)) then
         ! 
         ! Neigbour is within this parent branch
         !
         slot(:) = b%parent_slot(:)
         slot(idir) = ipo
         neigh  => b%parent%child(slot(1),slot(2),slot(3))%link
         repeat = .false.
      else
         !
         ! Neighbor is not within branch: we must go one level more toward base
         ! But memorize the slot, in order to be able to trace back later
         !
         chaincount = chaincount + 1
         if(chaincount.gt.amr_levelmax+1) stop 7209
         slotchain(:,chaincount) = b%parent_slot(:)
         b => b%parent
      endif
   else
      !
      ! Current branch is base branch (= part of base grid)
      !
      select case(idir)
      case(1)
         cmax = amr_grid_nx
      case(2)
         cmax = amr_grid_ny
      case(3)
         cmax = amr_grid_nz
      end select
      slot(:) = b%parent_slot(:)
      slot(idir) = slot(idir) + isgn
      if((slot(idir).ge.1).and.(slot(idir).le.cmax)) then
         if(associated(amr_grid_branch(slot(1),slot(2),slot(3))%link)) then
            neigh => amr_grid_branch(slot(1),slot(2),slot(3))%link
         endif
      elseif(amr_cyclic_xyz(idir)) then
         if(slot(idir).lt.1) then
            slot(idir) = cmax
         elseif(slot(idir).gt.cmax) then
            slot(idir) = 1
         else
            stop 9512
         endif
         if(associated(amr_grid_branch(slot(1),slot(2),slot(3))%link)) then
            neigh => amr_grid_branch(slot(1),slot(2),slot(3))%link
         endif
      endif
      repeat = .false.
   endif
enddo
!
! Now that we in principle found a neighbor, we must go deeper into the
! tree again if the level of the neighbor is more toward base than the 
! level of the original branch where we started
!
if((chaincount.gt.0).and.(associated(neigh))) then
   do i=chaincount,1,-1
      slot(:) = slotchain(:,i)
      slot(idir) = ip
      if(associated(neigh%child)) then
         if(associated(neigh%child(slot(1),slot(2),slot(3))%link)) then
            neigh => neigh%child(slot(1),slot(2),slot(3))%link
         else
            exit
         endif
      else
         exit
      endif
   enddo
endif
end subroutine amr_find_neighbor_branch


!--------------------------------------------------------------------------
!           Make sure that the neighbors of a also point to a
!
! We assume that neigh is the neighbor of a. If the level of neigh is more
! toward the stem of the tree, then apparently neigh is a bigger `cell' than
! a. In that case no linking back is necessary. But if the level is the same
! or more away from the stem, then linking back is necessary, AND also for
! all potential children of neigh. Note that this subroutine assumes
! that a is a branch without further children (as yet).
!
! NOTE: This subroutine is also used to relink neighbors in case a
!       branch is removed. In that case, if the parent of that removed
!       branch is the base grid (pointer == null) then a will be null.
! --------------------------------------------------------------------------
recursive subroutine amr_link_neighbors_back(a,neigh,idir,isgn)
implicit none
type(amr_branch), pointer :: a,neigh
integer :: idir,isgn,ip,ix,iy,iz
!
if(isgn.gt.0) then
   ip  = 1
else
   ip  = 2
endif
!
! Do consistency test
!
if(.not.associated(neigh)) stop 7190
!
! Do linking only if neigh size is equal or smaller than a
!
if(associated(a)) then
   if(a%level>neigh%level) then
      return
   endif
endif
!
! Link neigh back to a
!
neigh%neighbor(ip,idir)%link => a
!
! Now check for children in neigh and if they exist, then also link these
!
select case(idir)
case(1)
   ix = ip
   if(associated(neigh%child)) then
      do iz=1,1+amr_zdim
         do iy=1,1+amr_ydim
            if(associated(neigh%child(ix,iy,iz)%link)) then
               call amr_link_neighbors_back(a,neigh%child(ix,iy,iz)%link,idir,isgn)
            endif
         enddo
      enddo
   endif
case(2)
   iy = ip
   if(associated(neigh%child)) then
      do iz=1,1+amr_zdim
         do ix=1,1+amr_xdim
            if(associated(neigh%child(ix,iy,iz)%link)) then
               call amr_link_neighbors_back(a,neigh%child(ix,iy,iz)%link,idir,isgn)
            endif
         enddo
      enddo
   endif
case(3)
   iz = ip
   if(associated(neigh%child)) then
      do iy=1,1+amr_ydim
         do ix=1,1+amr_xdim
            if(associated(neigh%child(ix,iy,iz)%link)) then
               call amr_link_neighbors_back(a,neigh%child(ix,iy,iz)%link,idir,isgn)
            endif
         enddo
      enddo
   endif
end select
!
end subroutine amr_link_neighbors_back


!--------------------------------------------------------------------------
!              Find and link all neighbors for this branch
! 
! This subroutine will look in all directions to find neighboring
! branches of branch a, and will link the pointers to these branches
! (i.e. it will fill the a%neighbor pointers). It will also make sure
! that the opposite is linked: neighboring branches will link to this
! branch. Note that smaller neighboring branches will also be linked.
! And note that if a neighbor is bigger than the current one, then
! that neighbor will not link to us, because it is already linked to
! the parent or grandparent. This routine is called upon the construction
! of a new branch.
!--------------------------------------------------------------------------
subroutine amr_find_and_link_all_neighbors_branch(a)
implicit none
type(amr_branch), pointer :: a
type(amr_branch), pointer :: neigh
integer ix,iy,iz,idir
!
! Check
!
if(associated(a%child)) then
   write(stdo,*) 'INTERNAL ERROR in AMR Module: Finding and linking'
   write(stdo,*) '    neighbors is only allowed (and is'
   write(stdo,*) '    only necessary) for a fresh branch'
   write(stdo,*) '    that does not yet have children.'
   stop 590
endif
!
! Now do the linking from this branch to the neihbors, and
! also make sure that the linking from these neighbors (and
! possibly their children) goes back to us.
!
do idir=1,3
   if(amr_xyzdim(idir).gt.0) then
      call amr_find_neighbor_branch(a,idir,-1,neigh)
      a%neighbor(1,idir)%link => neigh
      if(associated(neigh)) then
         call amr_link_neighbors_back(a,neigh,idir,-1)
      endif
      call amr_find_neighbor_branch(a,idir,+1,neigh)
      a%neighbor(2,idir)%link => neigh
      if(associated(neigh)) then
         call amr_link_neighbors_back(a,neigh,idir,+1)
      endif
   endif
enddo
!
end subroutine amr_find_and_link_all_neighbors_branch



!!--------------------------------------------------------------------------
!!            Find representative cell values for larger branch
!!
!! Suppose you wish to know the cell values of some branch a which, however,
!! turns out to be refined into a whole tree of smaller branches.  This
!! subroutine will give you the average values and store them in the
!! cell data of the branch.
!! --------------------------------------------------------------------------
!recursive subroutine amr_compute_average_cellvalues_branch(a)
!implicit none
!type(amr_branch), pointer :: a,b
!doubleprecision :: dummy
!integer :: ix,iy,iz
!!
!! Check
!!
!if(.not.associated(a)) then
!   write(stdo,*) 'ERROR in AMR Module: Branch pointer not associated'
!   stop 6230
!endif
!!
!! Only need to do work if branch is a node 
!!
!if(.not.a%leaf) then
!   !
!   ! Reset user data
!   !
!   a%q  = 0.d0
!   !
!   ! Now go to all children, find the values and add them
!   ! to the average over this branch. Recursion will do the
!   ! rest.
!   !
!   dummy = 1.d0/(1.d0*amr_nrchildref)
!   do iz=1,1+amr_zdim
!      do iy=1,1+amr_ydim
!         do ix=1,1+amr_xdim
!            b => a%child(ix,iy,iz)%link
!            call amr_compute_average_cellvalues_branch(b)
!            a%q  = a%q  + dummy * b%q
!         enddo
!      enddo
!   enddo
!endif
!end subroutine amr_compute_average_cellvalues_branch


!--------------------------------------------------------------------------
!                           Refine a branch
!
! Assuming the grid is complete and consistent, but you wish to refine a
! leaf into 2x2x2 (for 3-D) or 2x2 (for 2-D etc). This is the routine for
! doing so. The branch you wish to refine is not allowed to already have
! children, i.e. it must be a leaf. If interpol==0 then these sub-branches
! will simply have the same cell values as the original branch. If
! interpol==1 then first order interpolation with neighboring branches will
! be done.
! --------------------------------------------------------------------------
subroutine amr_branch_refine(a,interpol)
implicit none
type(amr_branch), pointer :: a,b
integer :: interpol,slot(1:3)
integer :: ix,iy,iz,idir,ics,ierr
doubleprecision :: slope(nvar,3)
!
! Check for error
!
if(.not.associated(a)) then
   write(stdo,*) 'ERROR in AMR Module: Pointer to branch to be refined is not associated'
   stop 5901
endif
!
! Only refine if allowed
!
if(a%level.ge.amr_levelmax) then
   write(stdo,*) 'ERROR in AMR Module: I refuse to refine beyond level ',amr_levelmax
   write(stdo,*) '  Branch id = ',a%id
   stop 9301
endif
if(.not.a%leaf) then
   write(stdo,*) 'ERROR in AMR Module: Cannot refine a branch that is not a leaf'
   write(stdo,*) '  Branch id = ',a%id
   stop 9302
endif
!
! Switch off the leaf-status of this branch
!
a%leaf = .false.
!
! If the indexing (for external storage of data) is switched on, then
! we need to free a leaf index
!
if(amr_use_index) then
   call amr_free_leafindex(a%leafindex)
   a%leafindex=0
endif
!
! Decrease the counter for leafs
!
amr_nrleafs = amr_nrleafs - 1
!
! If linear interpolation, then prepare slopes
!
!if(interpol==1) then
!   do idir=1,3
!      ics = 0
!      if(associated(a%neighbor(1,idir)%link)) ics = ics+1
!      if(associated(a%neighbor(2,idir)%link)) ics = ics+2
!      select case(ics)
!      case(0)   ! No neighbors
!         slope(:,idir) = 0.d0
!      case(1)   ! Only left neighbor
!         slope(:,idir) = ( a%q(:) - & 
!              a%neighbor(1,idir)%link%q(:) ) / &
!              ( a%xc(idir) - a%neighbor(1,idir)%link%xc(idir) )
!      case(2)   ! Only right neighbor
!         slope(:,idir) = ( a%neighbor(2,idir)%link%q(:)  - & 
!              a%q(:) ) / &
!              ( a%neighbor(2,idir)%link%xc(idir) - a%xc(idir))
!      case(3)   ! Both left and right neighbor
!         slope(:,idir) = ( a%neighbor(2,idir)%link%q(:)  - & 
!              a%neighbor(1,idir)%link%q(:) ) / &
!              ( a%neighbor(2,idir)%link%xc(idir) - &
!              a%neighbor(1,idir)%link%xc(idir) )
!      end select
!      slope(:,idir) = slope(:,idir) * 0.25 * ( a%xi(2,idir) - a%xi(1,idir) )
!   enddo
!endif
!
! Allocate child-slots
!
allocate(a%child(1:2,1:2,1:2),STAT=ierr)
if(ierr.ne.0) then
   write(stdo,*) 'ERROR in AMR module while refining cell: '
   write(stdo,*) '        Could not allocate child slots to branch.'
   stop
endif
!
! Nullify these pointers first
!
do iz=1,1+amr_zdim
   do iy=1,1+amr_ydim
      do ix=1,1+amr_xdim
         nullify(a%child(ix,iy,iz)%link)
      enddo
   enddo
enddo
!
! Add a branch to each child slot
!
do iz=1,1+amr_zdim
   do iy=1,1+amr_ydim
      do ix=1,1+amr_xdim
         !
         ! Make the refined branch
         !
         slot(1) = ix
         slot(2) = iy
         slot(3) = iz
         call amr_branch_construct(b,a,slot)
         !
         ! Now fill the cell values 
         !
         !select case(interpol) 
         !case(0) 
         !   !
         !   ! No interpolation: just take value from parent 
         !   !
         !   b%q  = a%q
         !case(1)
         !   !
         !   ! First order interpolation
         !   ! *** USE amr_find_average_cellvalues_branch() 
         !   ! *** for the neighbor for the interpolation
         !   !
         !   b%q = a%q + ((ix-1)*2-1)*slope(:,1) + &
         !               ((iy-1)*2-1)*slope(:,2) + &
         !               ((iz-1)*2-1)*slope(:,3) 
         !end select
      enddo
   enddo
enddo
!
end subroutine amr_branch_refine


!--------------------------------------------------------------------------
!                        AMR branch derefine
!
! If you wish to un-refine, or in other words, if you wish to pack the
! entire sub-tree a into a single cell, then this is the subroutine you
! need
!--------------------------------------------------------------------------
subroutine amr_branch_unrefine(a)
implicit none
type(amr_branch), pointer :: a,b
integer :: ix,iy,iz,index
!
! If this branch is already a leaf, then  no derefinement is
! possible. 
!
if(a%leaf) then
   write(stdo,*) 'ERROR in AMR Module: Cannot unrefine a leaf branch, only a node branch'
   write(stdo,*) '  Branch id = ',a%id
   stop 8230
endif
!
! Find the average values for this sub-tree
! 
!call amr_compute_average_cellvalues_branch(a)
!
! Now destruct all children
!
do iz=1,1+amr_zdim
   do iy=1,1+amr_ydim
      do ix=1,1+amr_xdim
         if(.not.associated(a%child(ix,iy,iz)%link)) then
            write(stdo,*) 'INTERNAL ERROR in AMR Module: Cannot unrefine if child missing'
            write(stdo,*) '  Branch id = ',a%id
            stop 1951
         endif
         b => a%child(ix,iy,iz)%link
         call amr_branch_destruct(b)
      enddo
   enddo
enddo
!
! Deallocate the child slots
!
deallocate(a%child)
!
! Now change status to leaf
!
a%leaf = .true.
!
! If the indexing (for external storage of data) is switched on, then
! we need to find a new index slot for this leaf
!
if(amr_use_index) then
   call amr_assign_leafindex(a,index)
   a%leafindex = index
endif
!
! Increase the counter for leafs
!
amr_nrleafs = amr_nrleafs + 1
!
end subroutine amr_branch_unrefine


!--------------------------------------------------------------------------
!               Assign branchindex (for external data storage)
!--------------------------------------------------------------------------
subroutine amr_assign_branchindex(a,index)
implicit none
integer :: index,i
type(amr_branch), pointer :: a
!
! Stupidity check
!
if(.not.amr_use_index) then
   write(stdo,*) 'INTERNAL ERROR in AMR Module: Index mode not active'
   stop 6034
endif
!
! Find possible holes in the indexing, or else find position after last index
!
if(amr_branchindex_nrholes.eq.0) then
   !
   ! As long as the indices are not fragmented, we can simply use the next
   !
   !   ******************************************************************
   !   NOTE: These freeindex numbers are just a poor-mans solution to the
   !         N^2 problem of constructing branches. Later I should invent
   !         a better way to find free indices.
   !   ******************************************************************
   !
   if(amr_branchindex_freeindex.gt.amr_nrbranches_max) then
      write(stdo,*) 'ERROR in AMR module: no branch index free anymore...'
      write(stdo,*) amr_branchindex_freeindex,amr_nrbranches_max
      write(stdo,*) '   Presumably you are refining the AMR grid but you are'
      write(stdo,*) '   running out of memory slots for new cells. This is because'
      write(stdo,*) '   in the AMR module of RADMC-3D you have to estimate'
      write(stdo,*) '   in advance how many cells (including all octs)'
      write(stdo,*) '   you will at least need, so that the appropriate memory'
      write(stdo,*) '   can be allocated. You can increase this a-priori estimate'
      write(stdo,*) '   by increasing the nrbranchesmax argument of the call to'
      write(stdo,*) '   the amr_initialize() subroutine.'
      stop
   endif
   index = amr_branchindex_freeindex
   amr_branchindex_freeindex = amr_branchindex_freeindex + 1
else
   !
   ! Find a free hole from the hole database 
   !
   index = amr_branchindex_holes(amr_branchindex_nrholes)
   amr_branchindex_nrholes = amr_branchindex_nrholes - 1
endif
!
! Check that this branch index is not used yet, else produce error stop
!
if(associated(amr_index_to_branch(index)%link)) stop 3012
!
! Reserve this index
!
amr_index_to_branch(index)%link => a
!
end subroutine amr_assign_branchindex


!--------------------------------------------------------------------------
!               Assign leafindex (for external data storage)
!--------------------------------------------------------------------------
subroutine amr_assign_leafindex(a,index)
implicit none
integer :: index,i
type(amr_branch), pointer :: a
!
! Stupidity check
!
if(.not.amr_use_index) then
   write(stdo,*) 'INTERNAL ERROR in AMR Module: Index mode not active'
   stop 6034
endif
!
! Find possible holes in the indexing, or else find position after last index
!
if(amr_leafindex_nrholes.eq.0) then
   !
   ! As long as the indices are not fragmented, we can simply use the next
   !
   !   ******************************************************************
   !   NOTE: These freeindex numbers are just a poor-mans solution to the
   !         N^2 problem of constructing branches. Later I should invent
   !         a better way to find free indices.
   !   ******************************************************************
   !
   if(amr_leafindex_freeindex.gt.amr_nrleafs_max) then
      write(stdo,*) 'ERROR in AMR module: no leaf index free anymore...'
      write(stdo,*) amr_leafindex_freeindex,amr_nrleafs_max
      write(stdo,*) '   Presumably you are refining the AMR grid but you are'
      write(stdo,*) '   running out of memory slots for new cells. This is because'
      write(stdo,*) '   in the AMR module of RADMC-3D you have to estimate'
      write(stdo,*) '   in advance how many cells (only the leafs, not the octs)'
      write(stdo,*) '   you will at least need, so that the appropriate memory'
      write(stdo,*) '   can be allocated. You can increase this a-priori estimate'
      write(stdo,*) '   by increasing the nrleafsmax argument of the call to'
      write(stdo,*) '   the amr_initialize() subroutine.'
      stop
   endif
   index = amr_leafindex_freeindex
   amr_leafindex_freeindex = amr_leafindex_freeindex + 1
else
   !
   ! Find a free hole from the hole database 
   !
   index = amr_leafindex_holes(amr_leafindex_nrholes)
   amr_leafindex_nrholes = amr_leafindex_nrholes - 1
endif
!
! Check that this leaf index is not used yet, else produce error stop
!
if(associated(amr_index_to_leaf(index)%link)) stop 3012
!
! Reserve this index
!
amr_index_to_leaf(index)%link => a
!
end subroutine amr_assign_leafindex


!--------------------------------------------------------------------------
!               Remove branchindex (for external data storage)
!--------------------------------------------------------------------------
subroutine amr_free_branchindex(index)
implicit none
integer :: index
!
! Stupidity checks
!
if((index.lt.1).or.(index.gt.amr_nrbranches_max)) stop 1205
if(.not.amr_use_index) then
   write(stdo,*) 'INTERNAL ERROR in AMR Module: Index mode not active'
   stop 6036
endif
!
! Check if index was indeed used
!
if(associated(amr_index_to_branch(index)%link)) then
   !
   ! Now free this index
   !
   nullify(amr_index_to_branch(index)%link)
else
   !
   ! Internal error: this index should not have been free
   !
   write(stdo,*) 'INTERNAL ERROR in AMR Module: This index should not have been free:',index
   stop 8342
endif
!
! Decrease freeindex or add a hole in the hole database
!
if(index.eq.amr_branchindex_freeindex-1) then
   amr_branchindex_freeindex = index
else
   if(amr_branchindex_nrholes.ge.amr_nrbranches_max) stop 4922
   amr_branchindex_nrholes = amr_branchindex_nrholes + 1
   amr_branchindex_holes(amr_branchindex_nrholes) = index
endif
!
end subroutine amr_free_branchindex


!--------------------------------------------------------------------------
!               Remove leafindex (for external data storage)
!--------------------------------------------------------------------------
subroutine amr_free_leafindex(index)
implicit none
integer :: index
!
! Stupidity checks
!
if((index.lt.1).or.(index.gt.amr_nrbranches_max)) stop 1205
if(.not.amr_use_index) then
   write(stdo,*) 'INTERNAL ERROR in AMR Module: Index mode not active'
   stop 6036
endif
!
! Check if index was indeed used
!
if(associated(amr_index_to_leaf(index)%link)) then
   !
   ! Now free this index
   !
   nullify(amr_index_to_leaf(index)%link)
else
   !
   ! Internal error: this index should not have been free
   !
   write(stdo,*) 'INTERNAL ERROR in AMR Module: This index should not have been free:',index
   stop 8342
endif
!
! Decrease freeindex or add a hole in the hole database
!
if(index.eq.amr_leafindex_freeindex-1) then
   amr_leafindex_freeindex = index
else
   if(amr_leafindex_nrholes.ge.amr_nrleafs_max) stop 4922
   amr_leafindex_nrholes = amr_leafindex_nrholes + 1
   amr_leafindex_holes(amr_leafindex_nrholes) = index
endif
!
end subroutine amr_free_leafindex


!--------------------------------------------------------------------------
!                      Constructor for a branch 
!
! This routine produces a new branch and (if requested) links it
! automatically to a parent branch or into the base grid.
! 
! ARGUMENTS:
!   a         Pointer to the resulting new branch
!   b         Pointer to parent branch (if it exists). Make sure that
!             b is nullified if no parent branch exists.
!   slot      If parent exists, then this tells where the branch a
!             is inside the branch b. If no parent exists, then the
!             slot gives the position of this tree inside the base grid.
!--------------------------------------------------------------------------
subroutine amr_branch_construct(a,b,slot)
implicit none
type(amr_branch), pointer :: a
type(amr_branch), pointer :: b
integer :: slot(3)
integer :: ierr,ix,iy,iz,i,idir,index
doubleprecision :: dx,dy,dz
!
allocate(a,STAT=ierr)
if(ierr.ne.0) then
   write(stdo,*) 'ERROR in AMR Module: Could not allocate branch'
   stop 101
endif
amr_nrbranches = amr_nrbranches + 1
!
! Associate a unique identification number, which is, however,
! not the same as the position in the amr_thebranches array.
! As opposed to the latter, the ID is a unique and non-varying
! number. It will stay the same for the life time of the branch.
! When a branch is destructed, its ID will never again be used
! by another branch.
!
amr_last_branch_id = amr_last_branch_id + 1
a%id = amr_last_branch_id
!
! Initially a new branch is automatically a 'leaf'
!
a%leaf = .true.
!
! Increase the counter for leafs
!
amr_nrleafs = amr_nrleafs + 1
!
! If the indexing (for external storage of data) is switched on, then
! we need to find a free index slot and reserve this
!
if(amr_use_index) then
   call amr_assign_branchindex(a,index)
   a%branchindex = index
   call amr_assign_leafindex(a,index)
   a%leafindex = index
endif
!
! Make sure to nullify child slots pointer
! In this way this branch is born as a leaf
!
nullify(a%child)
!
! Now link
!
if(associated(b)) then
   !
   ! If parent branch exists, then link to parent branch
   !
   a%parent => b 
   if(.not.associated(b%child)) then
      write(stdo,*) 'ERROR in AMR Module: parent has no child slots'
   endif
   if(associated(b%child(slot(1),slot(2),slot(3))%link)) then
      write(stdo,*) 'ERROR in AMR Module: This child slot is already occupied'
      stop 102
   endif
   b%child(slot(1),slot(2),slot(3))%link => a
   !
   ! Assign level
   !
   a%level = b%level + 1
   !
   ! For the indexing of the amr_finegrid_xi and amr_finegrid_xc arrays
   ! we must calculate the ixyzf
   !
   do idir=1,3
      if(amr_xyzdim(idir).eq.1) then
         a%ixyzf(idir) = 2*(b%ixyzf(idir)-1)+slot(idir)
      else
         a%ixyzf(idir) = 1
      endif
   enddo
   !!
   !! Set coordinates of the branch boundaries and centers
   !!
   !do idir=1,3
   !   if(amr_xyzdim(idir).eq.1) then
   !      dx = ( b%xi(2,idir) - b%xi(1,idir) ) * 0.5d0
   !      if(slot(idir).eq.1) then
   !         a%xi(1,idir) = b%xi(1,idir)
   !         a%xi(2,idir) = b%xi(1,idir) + dx
   !      else
   !         a%xi(1,idir) = b%xi(1,idir) + dx
   !         a%xi(2,idir) = b%xi(2,idir)
   !      endif
   !   else
   !      a%xi(1,idir) = amr_grid_xi(1,idir)
   !      a%xi(2,idir) = amr_grid_xi(2,idir)
   !   endif
   !   a%xc(idir) = 0.5d0 * ( a%xi(1,idir) + a%xi(2,idir) )
   !enddo
   !!
   !! *** FOR TESTING PURPOSES: I CHECK WHETHER THE GRIDS ARE IDENTICAL ***
   !!
   !do idir=1,3
   !   if((a%xi(1,idir).ne.amr_finegrid_xi(a%ixyzf(idir),idir,a%level)).or.    &
   !      (a%xi(2,idir).ne.amr_finegrid_xi(a%ixyzf(idir)+1,idir,a%level)).or.  &
   !      (a%xc(idir).ne.amr_finegrid_xc(a%ixyzf(idir),idir,a%level))) then
   !      write(stdo,*) 'ERROR: AMR REFINED GRID PROBLEM: ILEVEL = ',a%level
   !      write(stdo,*) 'IDIR = ',idir,' ID = ',a%id,' IXYZF = ',a%ixyzf(idir)
   !      write(stdo,*) 'CELL: xi=',a%xi(1:2,idir)
   !      write(stdo,*) 'GRID: xi=',amr_finegrid_xi(a%ixyzf(idir),idir,a%level), &
   !                                amr_finegrid_xi(a%ixyzf(idir)+1,idir,a%level)
   !      write(stdo,*) 'CELL: xc=',a%xc(idir)
   !      write(stdo,*) 'GRID: xc=',amr_finegrid_xc(a%ixyzf(idir),idir,a%level)
   !      stop 1111
   !   endif
   !enddo
   !
else
   !
   ! No parent branch, then link into base grid
   !
   nullify(a%parent)
   !
   ! Assign level
   ! NOTE: This is now 0; used to be 1
   !
   a%level = 0
   !
   ! Safety check
   !
   if(((slot(1).lt.1).or.(slot(1).gt.amr_grid_nx)).or. &
      ((slot(2).lt.1).or.(slot(2).gt.amr_grid_ny)).or. &
      ((slot(3).lt.1).or.(slot(3).gt.amr_grid_nz))) then
      stop 7810
   endif
   !
   if(associated(amr_grid_branch(slot(1),slot(2),slot(3))%link)) then
      write(stdo,*) 'ERROR in AMR Module: Branch already associated...'
      stop 103
   endif
   amr_grid_branch(slot(1),slot(2),slot(3))%link => a
   !
   ! For the indexing of the amr_finegrid_xi and amr_finegrid_xc arrays
   ! we must calculate the ixyzf
   !
   do idir=1,3
      a%ixyzf(idir) = slot(idir)
   enddo
   !!
   !! Set coordinates of the branch boundaries
   !!
   !do idir=1,3
   !   if(amr_xyzdim(idir).eq.1) then
   !      a%xi(1,idir) = amr_grid_xi(slot(idir),idir)
   !      a%xi(2,idir) = amr_grid_xi(slot(idir)+1,idir)
   !   else
   !      a%xi(1,idir) = amr_grid_xi(1,idir)
   !      a%xi(2,idir) = amr_grid_xi(2,idir)
   !   endif
   !   a%xc(idir) = 0.5d0 * ( a%xi(1,idir) + a%xi(2,idir) )
   !enddo
   !!
   !! *** FOR TESTING PURPOSES: I CHECK WHETHER THE GRIDS ARE IDENTICAL ***
   !!
   !do idir=1,3
   !   if((a%xi(1,idir).ne.amr_finegrid_xi(a%ixyzf(idir),idir,a%level)).or.    &
   !      (a%xi(2,idir).ne.amr_finegrid_xi(a%ixyzf(idir)+1,idir,a%level)).or.  &
   !      (a%xc(idir).ne.amr_finegrid_xc(a%ixyzf(idir),idir,a%level))) then
   !      write(stdo,*) 'ERROR: AMR BASE GRID PROBLEM:'
   !      write(stdo,*) 'IDIR = ',idir,' ID = ',a%id,' IXYZF = ',a%ixyzf(idir)
   !      write(stdo,*) 'CELL: xi=',a%xi(1:2,idir)
   !      write(stdo,*) 'GRID: xi=',amr_finegrid_xi(a%ixyzf(idir),idir,a%level), &
   !                                amr_finegrid_xi(a%ixyzf(idir)+1,idir,a%level)
   !      write(stdo,*) 'CELL: xc=',a%xc(idir)
   !      write(stdo,*) 'GRID: xc=',amr_finegrid_xc(a%ixyzf(idir),idir,a%level)
   !      stop 1111
   !   endif
   !enddo
   !
endif
!
! Store where we are (integer coordinates within base grid or
! parent branch)
!
do i=1,3
   a%parent_slot(i) = slot(i)
enddo
!
! Now find and link to neighboring branches, always the deepest
! (most leaflike, i.e. smallest size) branches possible 
!
call amr_find_and_link_all_neighbors_branch(a)
!
! Done constructing. We now have a fully functional branch which is already
! linked to its neighbors (and its neighbors to it).
!
end subroutine amr_branch_construct


!--------------------------------------------------------------------------
!                      Destructor for a branch
!
! OPTIONAL: 
!   nil             If true then force neighbors to point to null after
!                   the destruction, also if a has a parent. Default = 
!                   .false.
!--------------------------------------------------------------------------
recursive subroutine amr_branch_destruct(a,nil)
implicit none
integer :: ix,iy,iz,idir
logical, optional :: nil
logical :: nill
!
type(amr_branch), pointer :: a,b
!
! Check if there is something to destruct
!
if(.not.associated(a)) return
!
! Parse the nil argument
!
if(present(nil)) then
   nill = nil
else
   nill = .false.
endif
!
! If a leaf, then decrease the counter for leafs
!
if(a%leaf) amr_nrleafs = amr_nrleafs - 1
!
! If the indexing (for external storage of data) is switched on, then
! we need to free the index slot 
!
if(amr_use_index) then
   call amr_free_branchindex(a%branchindex)
   if(a%leaf) then
      call amr_free_leafindex(a%leafindex)
   endif
endif
!
! Make sure to also destroy all children
!
if(associated(a%child)) then
   do iz=1,2
      do iy=1,2
         do ix=1,2
            if(associated(a%child(ix,iy,iz)%link)) then
               b => a%child(ix,iy,iz)%link
               call amr_branch_destruct(b)
            endif
         enddo
      enddo
   enddo
   deallocate(a%child)
endif
!
! Go to all neighbors and tell them that we will no longer be there
! NOTE: Base level is now 0; used to be 1.
!
if(a%level==0) then
   !
   ! The to-be-destructed branch is part of the base grid. In this
   ! case all neighbors, if they exist, must point to null when
   ! pointing to this slot
   !
   nullify(b)
   do idir=1,3
      if(associated(a%neighbor(1,idir)%link)) then
         call amr_link_neighbors_back(b,a%neighbor(1,idir)%link,idir,-1)
      endif
      if(associated(a%neighbor(2,idir)%link)) then
         call amr_link_neighbors_back(b,a%neighbor(2,idir)%link,idir,+1)
      endif
   enddo   
else
   ! 
   ! The to-be-destructed branch is part of a tree. All neighbors 
   ! which do not share the same parent will now be linked to the
   ! parent of a. Others will be linked to null.
   !
   nullify(b)
   do idir=1,3
      if(associated(a%neighbor(1,idir)%link)) then
         if((a%parent_slot(idir)==1).and.(.not.nill)) then
            call amr_link_neighbors_back(a%parent,a%neighbor(1,idir)%link,idir,-1)
         else
            call amr_link_neighbors_back(b,a%neighbor(1,idir)%link,idir,-1)
         endif
      endif
      if(associated(a%neighbor(2,idir)%link)) then
         if((a%parent_slot(idir)==2).and.(.not.nill)) then
            call amr_link_neighbors_back(a%parent,a%neighbor(2,idir)%link,idir,+1)
         else
            call amr_link_neighbors_back(b,a%neighbor(2,idir)%link,idir,+1)
         endif
      endif
   enddo
endif
!
! Unlink from parent
!
if(associated(a%parent)) then
   !
   ! There exist a parent branch
   !
   if(associated(a%parent%child)) then
      if(associated(a%parent%child(a%parent_slot(1), &
           a%parent_slot(2),a%parent_slot(3))%link)) then
         nullify(a%parent%child(a%parent_slot(1), &
              a%parent_slot(2),a%parent_slot(3))%link)
      else
         write(stdo,*) 'ERROR in AMR Module: Cannot find parent branch child...'
         stop 1002
      endif
   else
      write(stdo,*) 'ERROR in AMR Module: Cannot find parent branch child...'
      stop 1003
   endif
else
   !
   ! Ths present branch is part of the base grid
   !
   if(associated(amr_grid_branch(a%parent_slot(1), &
          a%parent_slot(2),a%parent_slot(3))%link)) then
      nullify(amr_grid_branch(a%parent_slot(1), &
          a%parent_slot(2),a%parent_slot(3))%link)
   endif
endif
!
! Now deallocate and reduce the number counter of branches
!
deallocate(a)
amr_nrbranches = amr_nrbranches - 1
!
end subroutine amr_branch_destruct


!--------------------------------------------------------------------------
!             Find cell in AMR base grid or in regular grid
! Note: If ix.lt.0 then this means that you are outside of a cell.
!--------------------------------------------------------------------------
subroutine amr_findbasecell(x,y,z,ix,iy,iz)
implicit none
doubleprecision :: x,y,z
integer :: ix,iy,iz
!
! First check out which base cell we are in
!
call amrhunt(amr_grid_xi(:,1),amr_grid_nx+1,x,ix)
call amrhunt(amr_grid_xi(:,2),amr_grid_ny+1,y,iy)
call amrhunt(amr_grid_xi(:,3),amr_grid_nz+1,z,iz)
!
! Check if we are outside of grid
!
if((ix.lt.1).or.(ix.gt.amr_grid_nx)) then
   if(amr_cyclic_xyz(1)) then
      if(ix.lt.1) then
         ix = amr_grid_nx
      else
         ix = 1
      endif
   else
      ix = -1
      iy = -1
      iz = -1
      return
   endif
endif
if((iy.lt.1).or.(iy.gt.amr_grid_ny)) then
   if(amr_cyclic_xyz(2)) then
      if(iy.lt.1) then
         iy = amr_grid_ny
      else
         iy = 1
      endif
   else
      ix = -1
      iy = -1
      iz = -1
      return
   endif
endif
if((iz.lt.1).or.(iz.gt.amr_grid_nz)) then
   if(amr_cyclic_xyz(3)) then
      if(iz.lt.1) then
         iz = amr_grid_nz
      else
         iz = 1
      endif
   else
      ix = -1
      iy = -1
      iz = -1
      return
   endif
endif
!
end subroutine amr_findbasecell


!--------------------------------------------------------------------------
!                       Find cell in AMR grid
!--------------------------------------------------------------------------
subroutine amr_findcell(x,y,z,a)
implicit none
doubleprecision :: x,y,z,dx,dy,dz
type(amr_branch), pointer :: a
integer :: ix,iy,iz,iddr
doubleprecision :: axi(1:2,1:3)
!
! First check out which base cell we are in
!
call amr_findbasecell(x,y,z,ix,iy,iz)
if(ix.lt.0) then
   nullify(a)
   return
endif
!
! Get base cell
!
a => amr_grid_branch(ix,iy,iz)%link
if(.not.associated(a)) then
   write(stdo,*) 'ERROR in AMR Module: Grid cell does not have branch/tree/cell'
   stop 1276
endif
!
! Now dig into tree until we reach leaf
!
do while(.not.a%leaf)
   do iddr=1,3
      axi(1,iddr) = amr_finegrid_xi(a%ixyzf(iddr),iddr,a%level)
      axi(2,iddr) = amr_finegrid_xi(a%ixyzf(iddr)+1,iddr,a%level)
   enddo
   dx = (x-axi(1,1)) / (axi(2,1)-axi(1,1))
   dy = (y-axi(1,2)) / (axi(2,2)-axi(1,2))
   dz = (z-axi(1,3)) / (axi(2,3)-axi(1,3))
   if((dx.gt.0.5d0).and.(amr_xdim.eq.1)) then
      ix = 2
   else
      ix = 1
   endif
   if((dy.gt.0.5d0).and.(amr_ydim.eq.1)) then
      iy = 2
   else
      iy = 1
   endif
   if((dz.gt.0.5d0).and.(amr_zdim.eq.1)) then
      iz = 2
   else
      iz = 1
   endif
   if(.not.associated(a%child)) then
      write(stdo,*) 'ERROR in AMR Module: Branch without children found...'
      stop 1673
   endif
   a => a%child(ix,iy,iz)%link
   if(.not.associated(a)) then
      write(stdo,*) 'ERROR in AMR Module: Tree is not complete...'
      stop 1674
   endif
enddo
return
end subroutine amr_findcell


!-------------------------------------------------------------------------
!                   A MORE GENERAL VERSION OF FINDCELL
!
! This function finds the cell independent of coordinate system and AMR
! refinement or not. Input is the cartesian (!) location x,y,z, and if
! spherical coordinates are used, these are automatically converted into
! r,theta,phi. So it is advisable to always use this subroutine instead
! of the more basic amr_findcell().
!-------------------------------------------------------------------------
subroutine amr_findcellindex_general(x,y,z,cellindex)
  implicit none
  double precision :: x,y,z,r,theta,phi
  type(amr_branch), pointer :: a
  integer :: ix,iy,iz,cellindex
  !
  cellindex=0
  if(amr_coordsystem.lt.100) then
     !
     ! Cartesian coordinates
     !
     if(amr_tree_present) then
        call amr_findcell(x,y,z,a)
        if(associated(a)) then
           cellindex = a%leafindex
        else
           cellindex = 0
        endif
     else
        call amr_findbasecell(x,y,z,ix,iy,iz)
        if(ix.gt.0) then
           cellindex = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
        else
           cellindex = 0
        endif
     endif
  elseif((amr_coordsystem.ge.100).and.(amr_coordsystem.lt.200)) then
     !
     ! Spherical coordinates
     !
     ! Convert x,y,z to spherical coordinates
     !
     call amr_xyz_to_rthphi(x,y,z,r,theta,phi)
     !
     ! Find the cell in which the star resides
     !
     if(theta.eq.pihalf) then
        theta = pihalf - 1d-8   ! Only for cell searching
     endif
     if(amr_tree_present) then
        call amr_findcell(r,theta,phi,a)
        if(associated(a)) then
           cellindex = a%leafindex
        else
           cellindex = 0
        endif
     else
        call amr_findbasecell(r,theta,phi,ix,iy,iz)
        if(ix.gt.0) then
           cellindex = ix+(iy-1)*amr_grid_nx+(iz-1)*amr_grid_nx*amr_grid_ny
        else
           cellindex = 0
        endif
     endif
  else
     write(stdo,*) 'ERROR: Other coordinate not yet implemented'
     stop 6491
  endif
end subroutine amr_findcellindex_general


!--------------------------------------------------------------------------
!                         STRING COMPARING
!--------------------------------------------------------------------------
function stringcompare(string1,string2,len)
implicit none
character*80 string1,string2
integer len
integer i
logical stringcompare
!
if(len.gt.80) stop
if(len.le.0) then
   stringcompare=.false.
   return
endif
stringcompare=.true.
do i=1,len
   if(string1(i:i).ne.string2(i:i)) stringcompare=.false.
enddo
return
end function stringcompare



!--------------------------------------------------------------------------
!             Helper routine for write grid and results
!--------------------------------------------------------------------------
recursive subroutine amr_write_grid_helper(a,formatstyle)
implicit none
type(amr_branch), pointer :: a,b
integer :: formatstyle
integer :: ix,iy,iz,iparent,ileaf,iddr
integer(kind=1) :: byte
doubleprecision :: axi(1:2,1:3)
!
if(formatstyle.eq.1) then 
   !
   ! Ascii output
   !
   if(a%leaf) then
      write(1,*) '0'
   else
      write(1,*) '1'
   endif
elseif(formatstyle.eq.2) then
   !
   ! F77-style unformatted output
   !
   if(a%leaf) then
      ileaf = 0
   else
      ileaf = 1
   endif
   write(1) ileaf
elseif(formatstyle.eq.3) then
   !
   ! C-style binary
   !
   if(a%leaf) then
      ileaf = 0
      byte  = 0
   else
      ileaf = 1
      byte  = 1
   endif
   write(1) byte
else
   write(stdo,*) 'ERROR in AMR module: Do not know format style ',formatstyle
   stop
endif
!
! Now go down the refine tree
!
if(.not.a%leaf) then
   if(.not.associated(a%child)) then
      write(stdo,*) 'ERROR in AMR Module: Branch without children found...'
      stop 7689
   endif
   do iz=1,1+amr_zdim
      do iy=1,1+amr_ydim
         do ix=1,1+amr_xdim
            b => a%child(ix,iy,iz)%link
            if(.not.associated(b)) then
               write(stdo,*) 'ERROR in AMR Module: Not all grid cells are filled'
               stop 6080
            endif
            call amr_write_grid_helper(b,formatstyle)
         enddo
      enddo
   enddo
endif
!
end subroutine amr_write_grid_helper


!--------------------------------------------------------------------------
!             Helper routine for read grid and results
!--------------------------------------------------------------------------
recursive subroutine amr_read_grid_helper(a,formatstyle)
implicit none
type(amr_branch), pointer :: a,b
integer :: formatstyle
integer :: ix,iy,iz,ileaf
doubleprecision xidum(2,3)
integer iparentdum,slotdum(3),iddum
integer(kind=8) :: iiformat,reclen,reclend,nn,kk,idum,iduma(1:7)
integer(kind=1) :: byte
!
! Check association of pointer a
!
if(.not.associated(a)) then
   write(stdo,*) 'ERROR in AMR Module during reading: branch is not yet constructed...'
   stop 7201
endif
!
! Formatted or unformatted?
!
if(formatstyle.eq.1) then
   !
   ! Ascii style
   !
   ! Check if this is a leaf or a node branch
   !
   read(1,*) ileaf
   ileaf = 1-ileaf
elseif(formatstyle.eq.2) then
   !
   ! F77-style unformatted
   !
   read(1) ileaf
   !
   ! Check if this is a leaf or a node branch
   !
   ileaf = 1-ileaf
elseif(formatstyle.eq.3) then
   !
   ! C-style binary I/O
   !
   ! Note: for the ileaf we use kind=1 integer (byte) 
   !
   read(1) byte
   ileaf = byte
   !
   ! Check if this is a leaf or a node branch
   !
   ileaf = 1-ileaf
else
   write(stdo,*) 'ERROR in AMR module: Unknown format style ',formatstyle
   stop
endif
!
! If a node branch, then we need to go deeper into the tree
!
if(ileaf.eq.0) then
   !
   ! Make the subtree 
   !
   call amr_branch_refine(a,0)
   !
   ! Now recursively call the reader for each sub-branch
   !
   if(.not.associated(a%child)) then
      write(stdo,*) 'ERROR in AMR Module: Branch without children found...'
      stop 6689
   endif
   do iz=1,1+amr_zdim
      do iy=1,1+amr_ydim
         do ix=1,1+amr_xdim
            b => a%child(ix,iy,iz)%link
            if(.not.associated(b)) then
               write(stdo,*) 'ERROR in AMR Module in amr_read_grid: child not present.'
               write(stdo,*) '   This is an internal error.'
               stop 5480
            endif
            call amr_read_grid_helper(b,formatstyle)
         enddo
      enddo
   enddo
else
   !
   ! It is a leaf, so check that there are indeed no links
   !
   if(associated(a%child)) then
      write(stdo,*) 'ERROR in AMR Module: Tree of the input file is not consistent'
      write(stdo,*) '   with internally present tree.'
      write(stdo,*) '   Internal tree branch has children that are not expected'
      stop 5481
   endif
endif
!
end subroutine amr_read_grid_helper


!--------------------------------------------------------------------------
!                               Write grid
!--------------------------------------------------------------------------
subroutine amr_write_grid(filename,formatstyle,gridinfo)
implicit none
character*160 :: filename
integer :: formatstyle
integer :: ix,iy,iz,iformat,nx,ny,nz,counter,gridinfo
integer :: incx,incy,incz,ierr,ilayer,ilyr
integer(kind=8) :: iiformat,reclen,reclend,nn,kk,idum,iduma(1:7)
type(amr_branch), pointer :: b
!
! Grid info is obsolete
!
if(gridinfo.ne.0) then
   write(stdo,*) 'ERROR: The extended style amr_grid.inp with gridinfo is'
   write(stdo,*) '       no longer supported.'
   stop
endif
!
! First do a stupidity check
!
if(amr_coordsystem.eq.-12345) then
   write(stdo,*) 'ERROR: The amr_coordsystem was not set before calling amr_write_grid()'
   write(stdo,*) '       This is not important for the AMR module, but could be '
   write(stdo,*) '       confusing for you. Just set amr_coordsystem=0 if you are not'
   write(stdo,*) '       interested in this value.'
   stop
endif
!
! First make sure that the linear lists of gridcells/branches are up to date
!
! IMPORTANT NOTE FOR FUTURE: It must be certain that the precise list order
! of the leafs produced by amr_compute_list_all() is exactly the same as
! the order in which the leafs are written out into the grid file in this
! subroutine. In that way we can simply dump any data as a 1-D array, and
! be sure that each number is correctly associated to its leaf.
!
if(amr_tree_present) then
   call amr_compute_list_all()
endif
!
! Check which dimensions are active and which are not
!
if(amr_xdim.eq.0) then
   nx   = 1
   incx = 0
else
   nx   = amr_grid_nx
   incx = 1
endif
if(amr_ydim.eq.0) then
   ny   = 1
   incy = 0
else
   ny   = amr_grid_ny
   incy = 1
endif
if(amr_zdim.eq.0) then
   nz   = 1
   incz = 0
else
   nz   = amr_grid_nz
   incz = 1
endif
!
! Current file format label
!
iformat = 1
iiformat = iformat
!
! Check if we write formatted or unformatted
!
if(formatstyle.eq.1) then
   !
   ! ASCII output
   !
   open(unit=1,file=filename)
   !
   ! A format number, so that we can always guarantee backward compatibility
   ! if we decide later to have another file format
   !
   write(1,*) iformat
   write(1,*) 
   !
   ! First the style of the AMR grid
   !   0   = Regular grid (no refinement)
   !   1   = Standard oct-tree AMR (but with base grid)
   !   10  = The "layers" mode. This is a kind of patchwork tree.
   !
   write(1,*) amr_style
   !
   ! The coordinate system 
   !
   write(1,*) amr_coordsystem
   !
   ! Gridinfo is obsolete, so it should be 0
   !
   write(1,*) 0
   !
   ! The dimensions of the problem. Note that if e.g. nx=1 then this means
   ! that in the base grid the x-direction only has 1 cell, but AMR refinement
   ! will still take also place in x-direction. Only if nx=0 the x-dimension
   ! is truly disabled. Same is true for y and z.
   !
   write(1,*) 
   write(1,*) incx,incy,incz
   write(1,*) nx,ny,nz
   !
   ! Some more data about the AMR tree
   !
   if((amr_style.ge.1).and.(amr_style.lt.10)) then
      write(1,*) 
      write(1,*) amr_levelmax,amr_nrleafs,amr_nrbranches
   endif
   !
   ! Write the base grid
   !
   write(1,*) 
   do ix=1,amr_grid_nx+1
      write(1,*) amr_grid_xi(ix,1)
   enddo
   write(1,*) 
   do iy=1,amr_grid_ny+1
      write(1,*) amr_grid_xi(iy,2)
   enddo
   write(1,*) 
   do iz=1,amr_grid_nz+1
      write(1,*) amr_grid_xi(iz,3)
   enddo
   write(1,*) 
   write(1,*) 
elseif(formatstyle.eq.2) then
   !
   ! F77-style unformatted output
   !
   open(unit=1,file=filename,form='unformatted')
   !
   ! A format number, so that we can always guarantee backward compatibility
   ! if we decide later to have another file format
   !
   write(1) iformat
   !
   ! Write which AMR style
   !   0   = Regular grid (no refinement)
   !   1   = Standard oct-tree AMR (but with base grid)
   !   10  = Patchwork AMR [NOT YET IMPLEMENTED]
   !
   write(1) amr_style
   !
   ! The coordinate system
   !
   write(1) amr_coordsystem
   !
   ! Gridinfo is obsolete, so it should be 0
   !
   write(1) 0
   !
   ! The dimensions of the problem. Note that if e.g. nx=1 then this means
   ! that in the base grid the x-direction only has 1 cell, but AMR refinement
   ! will still take also place in x-direction. Only if incx=0 it is truly
   ! disabled.
   !
   write(1) incx,incy,incz
   write(1) nx,ny,nz
   !
   ! Some more data about the AMR tree
   !
   if((amr_style.ge.1).and.(amr_style.lt.10)) then
      write(1) amr_levelmax,amr_nrleafs,amr_nrbranches
   endif
   !
   ! Write the base grid
   !
   write(1) (amr_grid_xi(ix,1),ix=1,amr_grid_nx+1)
   write(1) (amr_grid_xi(iy,2),iy=1,amr_grid_ny+1)
   write(1) (amr_grid_xi(iz,3),iz=1,amr_grid_nz+1)
elseif(formatstyle.eq.3) then
   !
   ! C-style binary output
   ! 
   ! Note: All integers are kind=8 integers (8-byte integers) in this mode
   !
   open(unit=1,file=filename,status='replace',access='stream')
   !
   ! A format number, so that we can always guarantee backward compatibility
   ! if we decide later to have another file format
   !
   write(1) iiformat
   !
   ! Write which AMR style
   !   0   = Regular grid (no refinement)
   !   1   = Standard oct-tree AMR (but with base grid)
   !   10  = Patchwork AMR [NOT YET IMPLEMENTED]
   !
   idum = amr_style
   write(1) idum
   !
   ! The coordinate system
   !
   idum = amr_coordsystem
   write(1) idum
   !
   ! Gridinfo is obsolete, so it should be 0
   !
   idum = 0
   write(1) idum
   !
   ! The dimensions of the problem. Note that if e.g. nx=1 then this means
   ! that in the base grid the x-direction only has 1 cell, but AMR refinement
   ! will still take also place in x-direction. Only if incx=0 it is truly
   ! disabled.
   !
   iduma(1) = incx
   iduma(2) = incy
   iduma(3) = incz
   write(1) iduma(1:3)
   iduma(1) = nx
   iduma(2) = ny
   iduma(3) = nz
   write(1) iduma(1:3)
   !
   ! Some more data about the AMR tree
   !
   if((amr_style.ge.1).and.(amr_style.lt.10)) then
      iduma(1) = amr_levelmax
      iduma(2) = amr_nrleafs
      iduma(3) = amr_nrbranches
      write(1) iduma(1:3)
   endif
   !
   ! Write the base grid
   !
   write(1) (amr_grid_xi(ix,1),ix=1,amr_grid_nx+1)
   write(1) (amr_grid_xi(iy,2),iy=1,amr_grid_ny+1)
   write(1) (amr_grid_xi(iz,3),iz=1,amr_grid_nz+1)
else
   write(stdo,*) 'ERROR in AMR module: Unknown format style ',formatstyle
   stop
endif
!
! Write the data of the base grid. Logically those values in the cells
! that are going to be refined lateron are only average values (or they
! are bogus values if the have not bee refreshed using the command
! amr_compute_average_cellvalues_branch() ).
!
if((amr_style.ge.1).and.(amr_style.lt.10)) then 
   !
   ! Oct-tree AMR file structure
   !
   do iz=1,amr_grid_nz
      do iy=1,amr_grid_ny
         do ix=1,amr_grid_nx
            b => amr_grid_branch(ix,iy,iz)%link
            if(.not.associated(b)) then
               write(stdo,*) 'ERROR in AMR Module in amr_write_grid: Not all grid cells are filled'
               stop 6090
            endif
            call amr_write_grid_helper(b,formatstyle)
         enddo
      enddo
   enddo
elseif((amr_style.ge.10).and.(amr_style.lt.20)) then
   !
   ! Layer-style AMR file structure
   !
   if(formatstyle.eq.1) then
      write(1,*) amr_layers_nrlevels,amr_layers_nrtot
   elseif(formatstyle.eq.2) then
      write(1) amr_layers_nrlevels,amr_layers_nrtot
   elseif(formatstyle.eq.3) then
      iduma(1) = amr_layers_nrlevels
      iduma(2) = amr_layers_nrtot
      write(1) iduma(1:2)
   else
      stop
   endif
   do ilayer=1,amr_layers_nrtot
      !
      ! Write the data for this layer
      !
!       BUG: The amr_layers_level() should not be written!
!      if(formatstyle.eq.1) then
!         write(1,*) amr_layers_level(ilayer),amr_layers_parent(ilayer), &
!              amr_layers_ixyz(1:3,ilayer),amr_layers_nxyz(1:3,ilayer)
!      elseif(formatstyle.eq.2) then
!         write(1) amr_layers_level(ilayer),amr_layers_parent(ilayer), &
!              amr_layers_ixyz(1:3,ilayer),amr_layers_nxyz(1:3,ilayer)
!      elseif(formatstyle.eq.3) then
!         iduma(1)   = amr_layers_level(ilayer)
!         iduma(2)   = amr_layers_parent(ilayer)
!         iduma(3:5) = amr_layers_ixyz(1:3,ilayer)
!         iduma(6:8) = amr_layers_nxyz(1:3,ilayer)
!         write(1) iduma(1:8)
!      else
!         stop
!      endif
      if(formatstyle.eq.1) then
         write(1,*) amr_layers_parent(ilayer), &
              amr_layers_ixyz(1:3,ilayer),amr_layers_nxyz(1:3,ilayer)
      elseif(formatstyle.eq.2) then
         write(1) amr_layers_parent(ilayer), &
              amr_layers_ixyz(1:3,ilayer),amr_layers_nxyz(1:3,ilayer)
      elseif(formatstyle.eq.3) then
         iduma(1)   = amr_layers_parent(ilayer)
         iduma(2:4) = amr_layers_ixyz(1:3,ilayer)
         iduma(5:7) = amr_layers_nxyz(1:3,ilayer)
         write(1) iduma(1:7)
      else
         stop
      endif
   enddo
endif
!
! Close file
!
close(1)
!
end subroutine amr_write_grid


!--------------------------------------------------------------------------
!                                  Read grid 
!--------------------------------------------------------------------------
subroutine amr_read_grid(filename,formatstyle,checkspher)
implicit none
character*160 :: filename,string
integer :: formatstyle
logical :: chsp
logical, optional :: checkspher
integer :: ix,iy,iz,nx,ny,nz,counter,l,ilayer,ilyr
integer :: levelmax,nrleafs,nrbranches,gridinfo
integer :: i,iformat
doubleprecision,allocatable :: xi(:),yi(:),zi(:)
integer :: incx,incy,incz,ierr
logical :: incl_x,incl_y,incl_z,zcyclic
integer(kind=8) :: iiformat,reclen,reclend,nn,kk,idum,iduma(1:7)
type(amr_branch), pointer :: b
logical :: midplane
!
! Defaults
!
zcyclic = .false.
!
! Interpret flag checkspher
!
chsp = .false.
if(present(checkspher)) then
   chsp = checkspher
endif
!
! Open file
!
if(formatstyle.eq.1) then
   open(unit=1,file=filename,status='old')
elseif(formatstyle.eq.2) then
   open(unit=1,file=filename,status='old',form='unformatted')
elseif(formatstyle.eq.3) then
   open(unit=1,file=filename,status='old',access='stream')
else
   write(stdo,*) 'INTERNAL ERROR: Do not know format style ',formatstyle
   stop
endif
!
! Which file format?
!
if(formatstyle.eq.1) then
   read(1,*) iformat
elseif(formatstyle.eq.2) then
   read(1) iformat
elseif(formatstyle.eq.3) then
   read(1) iiformat
   iformat=iiformat
else
   write(stdo,*) 'INTERNAL ERROR: Do not know format style ',formatstyle
   stop
endif
!
! Only support iformat 1 and greater
!
if(iformat.lt.1) then
   write(stdo,*) 'ERROR in AMR Module reading ',filename
   write(stdo,*) '  File format number ',iformat,' not supported'
   stop 2018
endif
!
! Check if we support this file format
!
if(iformat.gt.1) then
   write(stdo,*) 'ERROR in AMR Module reading ',filename
   write(stdo,*) '  File format number ',iformat,' not supported'
   stop 2019
endif
!
! Which structure does the AMR grid have?
!
if(formatstyle.eq.1) then
   read(1,*) amr_style
elseif(formatstyle.eq.2) then
   read(1) amr_style
elseif(formatstyle.eq.3) then
   read(1) idum
   amr_style = idum
else
   stop 67
endif
!
! Coordinate ssytem
!
if(formatstyle.eq.1) then
   read(1,*) amr_coordsystem
elseif(formatstyle.eq.2) then
   read(1) amr_coordsystem
elseif(formatstyle.eq.3) then
   read(1) idum
   amr_coordsystem = idum
else
   stop 67
endif
!
! Read whether this file contains redundant grid information 
!
if(formatstyle.eq.1) then
   read(1,*) gridinfo
elseif(formatstyle.eq.2) then
   read(1) gridinfo
elseif(formatstyle.eq.3) then
   read(1) idum
   gridinfo = idum
else
   stop 67
endif
!
! Grid info is obsolete
!
if(gridinfo.ne.0) then
   write(stdo,*) 'ERROR: The extended style amr_grid.inp with gridinfo is'
   write(stdo,*) '       no longer supported.'
   stop
endif
!
! The dimensions of the problem. Note that if e.g. nx=1 then this means
! that in the base grid the x-direction only has 1 cell, but AMR refinement
! will still take also place in x-direction. Only if nx=0 the x-dimension
! is truly disabled. Same is true for y and z.
!
if(formatstyle.eq.1) then
   read(1,*) incx,incy,incz
   read(1,*) nx,ny,nz
elseif(formatstyle.eq.2) then
   read(1) incx,incy,incz
   read(1) nx,ny,nz
elseif(formatstyle.eq.3) then
   read(1) iduma(1:3)
   incx = iduma(1)
   incy = iduma(2)
   incz = iduma(3)
   read(1) iduma(1:3)
   nx = iduma(1)
   ny = iduma(2)
   nz = iduma(3)
else
   stop 67
endif
!
! Check if all dimensions are OK
!
if((nx.lt.1).or.(ny.lt.1).or.(nz.lt.1)) then
   write(stdo,*) 'ERROR in AMR Module: nx,ny or nz < 1'
   stop 5103
endif
!
! Interpret the active/inactive dimensions
!
if(incx.ne.0) then
   incl_x=.true.
else
   incl_x=.false.
   if(nx.ne.1) then 
      write(stdo,*) 'ERROR in AMR Module: If x-dim inactive, must set nx=1'
      stop 1401
   endif
endif
if(incy.ne.0) then
   incl_y=.true.
else
   incl_y=.false.
   if(ny.ne.1) then 
      write(stdo,*) 'ERROR in AMR Module: If y-dim inactive, must set ny=1'
      stop 1402
   endif
endif
if(incz.ne.0) then
   incl_z=.true.
else
   incl_z=.false.
   if(nz.ne.1) then 
      write(stdo,*) 'ERROR in AMR Module: If z-dim inactive, must set nz=1'
      stop 1403
   endif
endif
!
! Some more data about the AMR tree
!
if((amr_style.lt.10).and.(amr_style.ge.0)) then
   if(amr_style.eq.0) then
      !
      ! Regular grid
      !
      levelmax    = 0
      nrleafs     = nx*ny*nz
      nrbranches  = nx*ny*nz
   else
      !
      ! Normal oct-tree for each base grid node
      !
      if(formatstyle.eq.1) then
         read(1,*) levelmax,nrleafs,nrbranches
      elseif(formatstyle.eq.2) then
         read(1) levelmax,nrleafs,nrbranches
      elseif(formatstyle.eq.3) then
         read(1) iduma(1:3)
         levelmax   = iduma(1)
         nrleafs    = iduma(2)
         nrbranches = iduma(3)
      else
         stop 67
      endif
   endif
elseif((amr_style.ge.10).and.(amr_style.lt.20)) then
   !
   ! Layer-style AMR
   !
   if(formatstyle.eq.1) then
      read(1,*) amr_layers_nrlevels,amr_layers_nrtot
   elseif(formatstyle.eq.2) then
      read(1) amr_layers_nrlevels,amr_layers_nrtot
   elseif(formatstyle.eq.3) then
      read(1) iduma(1:2)
      amr_layers_nrlevels = iduma(1)
      amr_layers_nrtot    = iduma(2)
   else
      stop 67
   endif
   !
   ! Do a check
   !
   ix= nx/2
   iy= ny/2
   iz= nz/2
   if(((2*ix.ne.nx).and.incl_x).or. &
      ((2*iy.ne.ny).and.incl_y).or. &
      ((2*iz.ne.nz).and.incl_z)) then
      write(stdo,*) 'ERROR in layers: Can only use layer-style AMR if base grid'
      write(stdo,*) '      sizes in all directions are even (2,4,6,...).'
      stop
   endif
else
   write(stdo,*) 'ERROR in amr: Other AMR types than oct-tree or layer-style not yet implemented'
   stop
endif
!
! Read the base grid
!
allocate(xi(1:nx+1),STAT=ierr)
if(ierr.ne.0) then
   write(stdo,*) 'ERROR in AMR Module: could not allocate xi'
   stop 7027
endif
allocate(yi(1:ny+1),STAT=ierr)
if(ierr.ne.0) then
   write(stdo,*) 'ERROR in AMR Module: could not allocate yi'
   stop 7027
endif
allocate(zi(1:nz+1),STAT=ierr)
if(ierr.ne.0) then
   write(stdo,*) 'ERROR in AMR Module: could not allocate zi'
   stop 7027
endif
if(formatstyle.eq.1) then
   read(1,*) (xi(ix),ix=1,nx+1)
   read(1,*) (yi(iy),iy=1,ny+1)
   read(1,*) (zi(iz),iz=1,nz+1)
elseif(formatstyle.eq.2) then
   read(1) (xi(ix),ix=1,nx+1)
   read(1) (yi(iy),iy=1,ny+1)
   read(1) (zi(iz),iz=1,nz+1)
elseif(formatstyle.eq.3) then
   read(1) (xi(ix),ix=1,nx+1)
   read(1) (yi(iy),iy=1,ny+1)
   read(1) (zi(iz),iz=1,nz+1)
else
   stop 67
endif
!
! If chsp is set, and the coordinates are indeed spherical, then do small
! adjustments to the angular coordinates to avoid problems with round-off
! errors.
!
if(chsp.and.((amr_coordsystem.ge.100).and.(amr_coordsystem.lt.200))) then
   if(.not.incl_x) then
      write(stdo,*) 'ERROR: In spherical coordinates: Cannot switch off r-direction.'
      stop
   endif
   if(.not.incl_y) then
      if(incl_z) then
         write(stdo,*) 'ERROR: In spherical coordinates: Cannot switch off theta-direction but include phi-direction.'
         stop
      endif
      yi(1) = 0.d0
      yi(2) = pi
   else
      if((abs(yi(1)).lt.1d-4*abs(yi(2)-yi(1))).and.(yi(1).ne.0.d0)) then
         write(stdo,*) '   Adjusting theta(1) to exactly 0...'
         yi(1) = 0.d0
      endif
      if((abs(yi(ny+1)-pihalf).lt.1d-4*abs(yi(ny+1)-yi(ny))).and.(yi(ny+1).ne.pihalf)) then
         write(stdo,*) '   Adjusting theta(ny+1) to exactly pi/2...'
         yi(ny+1) = pihalf
      endif
      if((abs(yi(ny+1)-pi).lt.1d-4*abs(yi(ny+1)-yi(ny))).and.(yi(ny+1).ne.pi)) then
         write(stdo,*) '   Adjusting theta(ny+1) to exactly pi...'
         yi(ny+1) = pi
      endif
      do iy=2,ny
         if((abs(yi(iy)-pihalf).lt.0.5d-4*abs(yi(iy+1)-yi(iy-1))).and.(yi(iy).ne.pihalf)) then
            write(string,*) iy
            write(stdo,*) '   Adjusting theta('//trim(adjustl(string))//') to exactly pi/2...'
            yi(iy) = pihalf
         endif
      enddo
      !
      ! For very fine theta grids near the equator with insufficient precision, sometimes
      ! the adjustment to pi/2 is accidently skipped (see 1e-4 factors above). Check for
      ! this, and stop if necessary.
      !
      if(abs(yi(ny+1)-pihalf)<1e-2*abs(yi(ny+1)-yi(1))) then
         !
         ! Should be mirror symmetry case
         !
         if(yi(ny+1).ne.pihalf) then
            write(stdo,*) 'SAFETY STOP: It seems that you want to use mirror symmetry in the'
            write(stdo,*) '   equatorial plane by putting the highest theta_i to pi/2. But'
            write(stdo,*) '   presumably you did not use high enough precision in the '
            write(stdo,*) '   theta_i list in the amr_grid.inp file. The mirror symmetry'
            write(stdo,*) '   is not set. Stopping for safety. If you are absolutely sure'
            write(stdo,*) '   you know what you are doing, you can comment out this stop.'
            stop
         endif
      endif
      !
      ! If theta goes beyond pi/2 then at least ONE theta interface must
      ! be exactly pihalf
      !
      if(yi(ny+1).gt.pihalf) then
         midplane=.false.
         do iy=2,ny
            !if(abs(yi(iy)-pihalf).lt.1d-10) then
            if(yi(iy).eq.pihalf) then
               midplane=.true.
            endif
         enddo
         if(.not.midplane) then
            write(stdo,*) 'ERROR: Theta grid crosses midplane, but there is no'
            write(stdo,*) '       interface at pi/2. Maybe check the precision '
            write(stdo,*) '       of the theta-grid: the equator interface must'
            write(stdo,*) '       be really exactly pi/2 to 12 digits or more.'
            stop
         endif
      endif
   endif
   if(.not.incl_z) then
      zi(1) = 0.d0
      zi(2) = twopi
   else
      if(zi(1).lt.-1d-4) then
         write(stdo,*) 'ERROR: In spherical coordinates, phi must be >=0.'
         stop
      endif
      if(zi(nz+1).gt.twopi+1d-4) then
         write(stdo,*) 'ERROR: In spherical coordinates, phi must be <=2*pi.'
         stop
      endif
      if((abs(zi(1)).lt.1d-2).and.(zi(1).ne.0.d0)) then
         write(stdo,*) '   Adjusting phi(1) to exactly 0...'
         zi(1) = 0.d0
      endif
      if((abs(zi(nz+1)-twopi).lt.1d-2).and.(zi(nz+1).ne.twopi)) then
         write(stdo,*) '   Adjusting phi(nz+1) to exactly 2*pi...'
         zi(nz+1) = twopi
      endif
      if((zi(1).eq.0.d0).and.(zi(nz+1).eq.twopi)) then
         zcyclic = .true.
         if((amr_coordsystem.ge.100).and.(amr_coordsystem.lt.200).and. &
              (nz.eq.1)) then
            write(stdo,*) 'ERROR while reading amr_grid.inp: '
            write(stdo,*) '      When using 1 grid cell in phi-direction,'
            write(stdo,*) '      in spherical coordinates, then you must'
            write(stdo,*) '      switch off the phi-direction.'
            stop
         endif
      else
         write(stdo,*) 'Please note: The phi-grid does not start exactly at 0 and end exactly at 2*pi.'
         write(stdo,*) '             This means that you want to model only a part of 360 degrees.'
         write(stdo,*) '             If this is not your intention, please make sure to set the'
         write(stdo,*) '             phi-grid exactly from 0 to 2*pi in double-precision.'
         write(stdo,*) '             On the other hand, if it IS your intention, then be aware'
         write(stdo,*) '             that RADMC-3D does not handle rays exiting phi=phi_min and'
         write(stdo,*) '             re-entering the grid at phi=phi_max (or vv).'
      endif
   endif
endif
!
! For the layered style AMR, we must now read the data of the layers
!
if((amr_style.ge.10).and.(amr_style.lt.20)) then
   !
   ! Read the layered-style AMR grid data
   ! 
   allocate(amr_layers_nr_per_level(1:amr_layers_nrlevels))
   allocate(amr_layers_level(0:amr_layers_nrtot))     ! Start from 0, i.e. main grid
   allocate(amr_layers_parent(1:amr_layers_nrtot))
   allocate(amr_layers_ixyz(1:3,1:amr_layers_nrtot))
   allocate(amr_layers_nnxyz(1:3,0:amr_layers_nrtot)) ! Start from 0, i.e. main grid
   allocate(amr_layers_nxyz(1:3,0:amr_layers_nrtot))  ! Start from 0, i.e. main grid
   allocate(amr_layers_map(0:amr_layers_nrtot))       ! Start from 0, i.e. main grid
   amr_layers_nr_per_level(:) = 0
   !
   ! First install the layer 0, i.e. the main grid
   !
   amr_layers_level(0)   = 0
   amr_layers_nnxyz(1,0) = nx
   amr_layers_nnxyz(2,0) = ny
   amr_layers_nnxyz(3,0) = nz
   amr_layers_nxyz(1,0)  = max(nx/2,1)
   amr_layers_nxyz(2,0)  = max(ny/2,1)
   amr_layers_nxyz(3,0)  = max(nz/2,1)
   !
   ! Now install all other layers
   !
   do ilayer=1,amr_layers_nrtot
      !
      ! Read the data for this layer
      !
      if(formatstyle.eq.1) then
         read(1,*) amr_layers_parent(ilayer), &
              amr_layers_ixyz(1:3,ilayer),amr_layers_nxyz(1:3,ilayer)
      elseif(formatstyle.eq.2) then
         read(1) amr_layers_parent(ilayer), &
              amr_layers_ixyz(1:3,ilayer),amr_layers_nxyz(1:3,ilayer)
      elseif(formatstyle.eq.3) then
         read(1) iduma(1:7)
         amr_layers_parent(ilayer)   = iduma(1)
         amr_layers_ixyz(1:3,ilayer) = iduma(2:4)
         amr_layers_nxyz(1:3,ilayer) = iduma(5:7)
      else
         stop 67
      endif
      if(incl_x) then
         amr_layers_nnxyz(1,ilayer) = 2*amr_layers_nxyz(1,ilayer)
      else
         amr_layers_nxyz(1,ilayer)  = 1
         amr_layers_nnxyz(1,ilayer) = 1
      endif
      if(incl_y) then
         amr_layers_nnxyz(2,ilayer) = 2*amr_layers_nxyz(2,ilayer)
      else
         amr_layers_nxyz(2,ilayer)  = 1
         amr_layers_nnxyz(2,ilayer) = 1
      endif
      if(incl_z) then
         amr_layers_nnxyz(3,ilayer) = 2*amr_layers_nxyz(3,ilayer)
      else
         amr_layers_nxyz(3,ilayer)  = 1
         amr_layers_nnxyz(3,ilayer) = 1
      endif
      !
      ! Do basic checks
      !
      if((amr_layers_parent(ilayer).lt.0).or.                 &
         (amr_layers_parent(ilayer).ge.ilayer)) then
         write(stdo,*) 'ERROR in layers: id of parent layer out of range...'
         write(stdo,*) '                 Parent id must be < current id.'
         stop
      endif
      !
      ! Compute the level of the current layer
      !
      amr_layers_level(ilayer) = amr_layers_level(amr_layers_parent(ilayer)) + 1
      !
      ! Check that the layer fits in its parent
      !
      if((amr_layers_ixyz(1,ilayer).lt.1).or.                       &
         (amr_layers_ixyz(1,ilayer)+amr_layers_nxyz(1,ilayer)-1.gt. &
          amr_layers_nnxyz(1,amr_layers_parent(ilayer))).or.         &
         (amr_layers_ixyz(2,ilayer).lt.1).or.                       &
         (amr_layers_ixyz(2,ilayer)+amr_layers_nxyz(2,ilayer)-1.gt. &
          amr_layers_nnxyz(2,amr_layers_parent(ilayer))).or.         &
         (amr_layers_ixyz(3,ilayer).lt.1).or.                       &
         (amr_layers_ixyz(3,ilayer)+amr_layers_nxyz(3,ilayer)-1.gt. &
          amr_layers_nnxyz(3,amr_layers_parent(ilayer)))) then
         write(stdo,*) 'ERROR in layers: layer out of bounds'
         write(stdo,*) 'layer number = ',ilayer
         write(stdo,*) 'layer parent = ',amr_layers_parent(ilayer)
         write(stdo,*) 'nnxyz parent = ', &
              amr_layers_nnxyz(1,amr_layers_parent(ilayer)), &
              amr_layers_nnxyz(2,amr_layers_parent(ilayer)), &
              amr_layers_nnxyz(3,amr_layers_parent(ilayer))
         write(stdo,*) 'ixyz         = ', &
              amr_layers_ixyz(1,ilayer),&
              amr_layers_ixyz(2,ilayer),&
              amr_layers_ixyz(3,ilayer)
         write(stdo,*) 'nxyz         = ', &
              amr_layers_nxyz(1,ilayer),&
              amr_layers_nxyz(2,ilayer),&
              amr_layers_nxyz(3,ilayer)
         stop
      endif
      !
      ! Check that the layer does not overlap with another layer
      !
      do ilyr=1,ilayer-1
         if(amr_layers_parent(ilayer).eq.amr_layers_parent(ilyr)) then
            if((amr_layers_ixyz(1,ilayer).le.amr_layers_ixyz(1,ilyr)+amr_layers_nxyz(1,ilyr)-1).and. &
               (amr_layers_ixyz(1,ilayer)+amr_layers_nxyz(1,ilayer)-1.ge.amr_layers_ixyz(1,ilyr)).and. &
               (amr_layers_ixyz(2,ilayer).le.amr_layers_ixyz(2,ilyr)+amr_layers_nxyz(2,ilyr)-1).and. &
               (amr_layers_ixyz(2,ilayer)+amr_layers_nxyz(2,ilayer)-1.ge.amr_layers_ixyz(2,ilyr)).and. &
               (amr_layers_ixyz(3,ilayer).le.amr_layers_ixyz(3,ilyr)+amr_layers_nxyz(3,ilyr)-1).and. &
               (amr_layers_ixyz(3,ilayer)+amr_layers_nxyz(3,ilayer)-1.ge.amr_layers_ixyz(3,ilyr))) then
               write(stdo,*) 'ERROR in layers: Layers overlap'
               stop
            endif
         endif
      enddo
      !
   enddo
   !
   ! Set here some important AMR variables
   !
   ! NOTE: In previous versions there used to be a "+ 1" in the formula
   !       below, but now that the oct-tree AMR levels also start with 0
   !       (contrary to earlier versions, which started with 1), this +1
   !       is no longer necessary
   !
   !!! levelmax   = amr_layers_nrtot  ! BUGFIX Pak Shing Li 03.06.2012
   levelmax   = amr_layers_nrlevels
   nrleafs    = 0
   nrbranches = 0
   do ilayer=0,amr_layers_nrtot
      nrbranches = nrbranches + amr_layers_nnxyz(1,ilayer)* & 
                                amr_layers_nnxyz(2,ilayer)* &
                                amr_layers_nnxyz(3,ilayer)
      nrleafs    = nrleafs    + amr_layers_nnxyz(1,ilayer)* & 
                                amr_layers_nnxyz(2,ilayer)* &
                                amr_layers_nnxyz(3,ilayer)  
      if(ilayer.gt.0) then
         nrleafs    = nrleafs - amr_layers_nxyz(1,ilayer)*  & 
                                amr_layers_nxyz(2,ilayer)*  &
                                amr_layers_nxyz(3,ilayer) 
      endif
   enddo
   !
endif
!
! Check that levelmax is not rediculously large
!
if(levelmax.gt.20) then
   write(stdo,*) 'ERROR in AMR Module: levelmax is too large (>20)...'
   write(stdo,*) '     You can make it larger than 20, in principle, but'
   write(stdo,*) '     since it is highly unlikely (as of Dec 2009) that'
   write(stdo,*) '     such extreme refinement is truly requested, we'
   write(stdo,*) '     assume that this is an error, and therefore stop.'
   write(stdo,*) '     Remove this error message in the code if you are'
   write(stdo,*) '     serious about this high level of refinement.'
   stop
endif
!
! Now initialize the AMR grid
!
if(amr_always_use_tree) then
   call amr_initialize(incl_x,incl_y,incl_z,nx,ny,nz,xi,yi,zi,levelmax,1,   &
                       nrbranches,nrleafs,zcyclic=zcyclic,&
                       always_amrtree=.true.)
else
   call amr_initialize(incl_x,incl_y,incl_z,nx,ny,nz,xi,yi,zi,levelmax,1,   &
                       nrbranches,nrleafs,zcyclic=zcyclic)
endif
!
! Read the data and refine the grid
!
if((amr_style.ge.1).and.(amr_style.lt.10)) then
   !
   ! Oct-tree AMR file structure
   !
   do iz=1,amr_grid_nz
      do iy=1,amr_grid_ny
         do ix=1,amr_grid_nx
            b => amr_grid_branch(ix,iy,iz)%link
            if(.not.associated(b)) then
               write(stdo,*) 'ERROR in AMR Module in amr_read_grid: Not all grid cells are filled'
               stop 2090
            endif
            call amr_read_grid_helper(b,formatstyle)
         enddo
      enddo
   enddo
elseif((amr_style.ge.10).and.(amr_style.lt.20)) then
   !
   ! Layer-style AMR structure.
   ! No need to read any further data. We just have to refine the 
   ! grid. 
   !
   call amr_layers_implement_amr_grid()
   !
endif
!
! Close the file
!
close(1)
!
! Check some things
!
if(amr_nrbranches.ne.amr_nrbranches_max) then
   write(stdo,*) 'AMR Warning: Number of branches smaller than given in header of grid file'
   write(stdo,*) '    In header:      ',amr_nrbranches_max
   write(stdo,*) '    Branches found: ',amr_nrbranches
endif
if(amr_nrleafs.ne.amr_nrleafs_max) then
   write(stdo,*) 'AMR Warning: Number of leafs smaller than given in header of grid file'
   write(stdo,*) '    In header:      ',amr_nrleafs_max
   write(stdo,*) '    Leafs found:    ',amr_nrleafs
endif
!
! Deallocate stuff
!
if(allocated(xi)) deallocate(xi)
if(allocated(yi)) deallocate(yi)
if(allocated(zi)) deallocate(zi)
!
! Now create the list of branches and leafs
!
! IMPORTANT NOTE FOR FUTURE: It must be certain that the precise list order
! of the leafs produced by amr_compute_list_all() is exactly the same as
! the order in which the leafs are read in from the grid file in this
! subroutine. In that way we can simply dump any data as a 1-D array, and
! be sure that each number is correctly associated to its leaf.
!
if(amr_tree_present) then
   call amr_compute_list_all()
endif
!
! Done reading the grid
!
end subroutine amr_read_grid


!--------------------------------------------------------------------------
!             3-D Subbox routine (only for user-convenience)
!                Works only for external (indexed) data
!         Note: 1-D and 2-D slices can be made trivially as well
!
! The goal of this routine is to return a regularly spaced array of a
! variable in the model in a box with dimensions [x0,x1],[y0,y1],[z0,z1]
! and grid cells nx,ny,nz. This is only for analyzing the 3-D data. 
!--------------------------------------------------------------------------
subroutine amr_subbox_3d(nv,nc,nx,ny,nz,x0,x1,y0,y1,z0,z1,            &
                         phi1,theta,phi2,func,funcslice,              &
                         spherical)
  implicit none
  integer :: nv,nc,nx,ny,nz
  double precision :: x0,x1,y0,y1,z0,z1,x,y,z
  double precision :: func(nv,nc),funcslice(nx,ny,nz,nv)
  double precision :: xbk,ybk,zbk
  double precision :: sinphi1,cosphi1,sinphi2,cosphi2,sintheta,costheta
  double precision :: phi1,theta,phi2
  double precision :: r,th,ph
  integer :: ix,iy,iz,index,iv,ixx,iyy,izz
  type(amr_branch), pointer :: a
  logical, optional :: spherical
  logical :: cart,spher
  !
  cart  = .true.
  spher = .false.
  if(present(spherical)) then
     if(spherical) then
        cart  = .false.
        spher = .true.
     endif
  endif
  !
  ! NOTE: As of version 2.0 of RADMC-3D these angles are in degrees, no
  !       longer in radian
  !
  sinphi1  = sin(phi1*pi/180.)
  cosphi1  = cos(phi1*pi/180.)
  sintheta = sin(theta*pi/180.)
  costheta = cos(theta*pi/180.)
  sinphi2  = sin(phi2*pi/180.)
  cosphi2  = cos(phi2*pi/180.)
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
           ! Now find the cell in which we are (NOTE: There is no
           ! interpolation done, so mind the danger for aliasing!)
           !
           if(amr_tree_present) then
              !
              ! Use the AMR tree
              !
              if(cart) then
                 call amr_findcell(x,y,z,a)
              elseif(spher) then
                 call amr_xyz_to_rthphi(x,y,z,r,th,ph)
                 call amr_findcell(r,th,ph,a)
              else
                 stop 6773
              endif
              !
              ! And put the value found in the AMR grid in the 3-D regular
              ! subbox grid
              !
              if(associated(a)) then
                 index = a%leafindex
                 do iv=1,nv
                    funcslice(ix,iy,iz,iv) = func(iv,index)
                 enddo
              else
                 funcslice(ix,iy,iz,1:nv) = 0.d0
              endif
           else
              !
              ! Regular grid, i.e. no AMR tree available
              !
              if(cart) then
                 call amr_findbasecell(x,y,z,ixx,iyy,izz)
              elseif(spher) then
                 call amr_xyz_to_rthphi(x,y,z,r,th,ph)
                 call amr_findbasecell(r,th,ph,ixx,iyy,izz)
              else
                 stop 6873
              endif
              !
              ! And put the value found in the AMR grid in the 3-D regular
              ! subbox grid
              !
              if(ixx.gt.0) then
                 index = ixx+(iyy-1)*amr_grid_nx+(izz-1)*amr_grid_nx*amr_grid_ny
                 do iv=1,nv
                    funcslice(ix,iy,iz,iv) = func(iv,index)
                 enddo
              else
                 funcslice(ix,iy,iz,1:nv) = 0.d0
              endif
           endif
        enddo
     enddo
  enddo
  !
end subroutine amr_subbox_3d


!--------------------------------------------------------------------------
!             3-D Subbox routine (only for user-convenience)
!                Works only for external (indexed) data
!         Note: 1-D and 2-D slices can be made trivially as well
!
! The goal of this routine is to return a regularly spaced array of a
! variable in the model in a box with dimensions [x0,x1],[y0,y1],[z0,z1]
! and grid cells nx,ny,nz. This is only for analyzing the 3-D data. 
!--------------------------------------------------------------------------
subroutine amr_subbox2_3d(nv1,nv2,nc,nx,ny,nz,x0,x1,y0,y1,z0,z1,      &
                         phi1,theta,phi2,func,funcslice,              &
                         spherical)
  implicit none
  integer :: nv1,nv2,nc,nx,ny,nz
  double precision :: x0,x1,y0,y1,z0,z1,x,y,z
  double precision :: func(nv1,nv2,nc),funcslice(nx,ny,nz,nv1,nv2)
  double precision :: xbk,ybk,zbk
  double precision :: sinphi1,cosphi1,sinphi2,cosphi2,sintheta,costheta
  double precision :: phi1,theta,phi2
  double precision :: r,th,ph
  integer :: ix,iy,iz,index,iv1,iv2,ixx,iyy,izz
  type(amr_branch), pointer :: a
  logical, optional :: spherical
  logical :: cart,spher
  !
  cart  = .true.
  spher = .false.
  if(present(spherical)) then
     if(spherical) then
        cart  = .false.
        spher = .true.
     endif
  endif
  !
  ! NOTE: As of version 2.0 of RADMC-3D these angles are in degrees, no
  !       longer in radian
  !
  sinphi1  = sin(phi1*pi/180.)
  cosphi1  = cos(phi1*pi/180.)
  sintheta = sin(theta*pi/180.)
  costheta = cos(theta*pi/180.)
  sinphi2  = sin(phi2*pi/180.)
  cosphi2  = cos(phi2*pi/180.)
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
           ! Now find the cell in which we are (NOTE: There is no
           ! interpolation done, so mind the danger for aliasing!)
           !
           if(amr_tree_present) then
              !
              ! Use the AMR tree
              !
              if(cart) then
                 call amr_findcell(x,y,z,a)
              elseif(spher) then
                 call amr_xyz_to_rthphi(x,y,z,r,th,ph)
                 call amr_findcell(r,th,ph,a)
              else
                 stop 6773
              endif
              !
              ! And put the value found in the AMR grid in the 3-D regular
              ! subbox grid
              !
              if(associated(a)) then
                 index = a%leafindex
                 do iv2=1,nv2
                    do iv1=1,nv1
                       funcslice(ix,iy,iz,iv1,iv2) = func(iv1,iv2,index)
                    enddo
                 enddo
              else
                 funcslice(ix,iy,iz,1:nv1,1:nv2) = 0.d0
              endif
           else
              !
              ! Regular grid, i.e. no AMR tree available
              !
              if(cart) then
                 call amr_findbasecell(x,y,z,ixx,iyy,izz)
              elseif(spher) then
                 call amr_xyz_to_rthphi(x,y,z,r,th,ph)
                 call amr_findbasecell(r,th,ph,ixx,iyy,izz)
              else
                 stop 6873
              endif
              !
              ! And put the value found in the AMR grid in the 3-D regular
              ! subbox grid
              !
              if(ixx.gt.0) then
                 index = ixx+(iyy-1)*amr_grid_nx+(izz-1)*amr_grid_nx*amr_grid_ny
                 do iv2=1,nv2
                    do iv1=1,nv1
                       funcslice(ix,iy,iz,iv1,iv2) = func(iv1,iv2,index)
                    enddo
                 enddo
              else
                 funcslice(ix,iy,iz,1:nv1,1:nv2) = 0.d0
              endif
           endif
        enddo
     enddo
  enddo
  !
end subroutine amr_subbox2_3d


!--------------------------------------------------------------------------
!             3-D sampling routine (only for user-convenience)
!                Works only for external (indexed) data
!
! The goal of this routine is to return the values of a variable in the 
! model at a set of points who's locations are given in a list. These
! points can be completely arbitrary, e.g. [(2.3,4.,-1),(1.7,1.,-1.1),
! (-3,-2,4.2),(0,0,3.2)] as an example of a list of four points.
!--------------------------------------------------------------------------
subroutine amr_sample_3d(nv,nc,npt,xpt,ypt,zpt,    &
                         func,values,              &
                         spherical)
  implicit none
  integer :: nv,nc,npt,ipt
  double precision :: x,y,z,xpt(npt),ypt(npt),zpt(npt)
  double precision :: func(nv,nc),values(npt,nv)
  double precision :: xbk,ybk,zbk
  double precision :: r,th,ph
  integer :: ixx,iyy,izz,index,iv
  type(amr_branch), pointer :: a
  logical, optional :: spherical
  logical :: cart,spher
  !
  cart  = .true.
  spher = .false.
  if(present(spherical)) then
     if(spherical) then
        cart  = .false.
        spher = .true.
     endif
  endif
  !
  do ipt=1,npt
     !
     ! Get the (x,y,z) position
     !
     x = xpt(ipt)
     y = ypt(ipt)
     z = zpt(ipt)
     !
     ! Now find the cell in which we are (NOTE: There is no
     ! interpolation done, so mind the danger for aliasing!)
     !
     if(amr_tree_present) then 
        !
        ! If there is an AMR tree, then we use that
        !
        if(cart) then
           call amr_findcell(x,y,z,a)
        elseif(spher) then
           call amr_xyz_to_rthphi(x,y,z,r,th,ph)
           call amr_findcell(r,th,ph,a)
        else
           stop 6773
        endif
        !
        ! And put the value found in the AMR grid in the values list
        !
        if(associated(a)) then
           index = a%leafindex
           do iv=1,nv
              values(ipt,iv) = func(iv,index)
           enddo
        else
           values(ipt,1:nv) = 0.d0
        endif
     else
        !
        ! If we have a regular grid (no AMR tree) then we use ixx, iyy, izz
        !
        if(cart) then
           call amr_findbasecell(x,y,z,ixx,iyy,izz)
        elseif(spher) then
           call amr_xyz_to_rthphi(x,y,z,r,th,ph)
           call amr_findbasecell(r,th,ph,ixx,iyy,izz)
        else
           stop 6873
        endif
        !
        ! And put the value found in the AMR grid in the values list
        !
        if(ixx.gt.0) then
           index = ixx+(iyy-1)*amr_grid_nx+(izz-1)*amr_grid_nx*amr_grid_ny
           do iv=1,nv
              values(ipt,iv) = func(iv,index)
           enddo
        else
           values(ipt,1:nv) = 0.d0
        endif
     endif
  enddo
  !
end subroutine amr_sample_3d


!--------------------------------------------------------------------------
!             3-D sampling routine (only for user-convenience)
!                Works only for external (indexed) data
!
! The goal of this routine is to return the values of a variable in the 
! model at a set of points who's locations are given in a list. These
! points can be completely arbitrary, e.g. [(2.3,4.,-1),(1.7,1.,-1.1),
! (-3,-2,4.2),(0,0,3.2)] as an example of a list of four points.
!--------------------------------------------------------------------------
subroutine amr_sample2_3d(nv1,nv2,nc,npt,xpt,ypt,zpt,    &
                         func,values,                    &
                         spherical)
  implicit none
  integer :: nv1,nv2,nc,npt,ipt
  double precision :: x,y,z,xpt(npt),ypt(npt),zpt(npt)
  double precision :: func(nv1,nv2,nc),values(npt,nv1,nv2)
  double precision :: xbk,ybk,zbk
  double precision :: r,th,ph
  integer :: ixx,iyy,izz,index,iv1,iv2
  type(amr_branch), pointer :: a
  logical, optional :: spherical
  logical :: cart,spher
  !
  cart  = .true.
  spher = .false.
  if(present(spherical)) then
     if(spherical) then
        cart  = .false.
        spher = .true.
     endif
  endif
  !
  do ipt=1,npt
     !
     ! Get the (x,y,z) position
     !
     x = xpt(ipt)
     y = ypt(ipt)
     z = zpt(ipt)
     !
     ! Now find the cell in which we are (NOTE: There is no
     ! interpolation done, so mind the danger for aliasing!)
     !
     if(amr_tree_present) then 
        !
        ! If there is an AMR tree, then we use that
        !
        if(cart) then
           call amr_findcell(x,y,z,a)
        elseif(spher) then
           call amr_xyz_to_rthphi(x,y,z,r,th,ph)
           call amr_findcell(r,th,ph,a)
        else
           stop 6773
        endif
        !
        ! And put the value found in the AMR grid in the values list
        !
        if(associated(a)) then
           index = a%leafindex
           do iv2=1,nv2
              do iv1=1,nv1
                 values(ipt,iv1,iv2) = func(iv1,iv2,index)
              enddo
           enddo
        else
           values(ipt,1:nv1,1:nv2) = 0.d0
        endif
     else
        !
        ! If we have a regular grid (no AMR tree) then we use ixx, iyy, izz
        !
        if(cart) then
           call amr_findbasecell(x,y,z,ixx,iyy,izz)
        elseif(spher) then
           call amr_xyz_to_rthphi(x,y,z,r,th,ph)
           call amr_findbasecell(r,th,ph,ixx,iyy,izz)
        else
           stop 6873
        endif
        !
        ! And put the value found in the AMR grid in the values list
        !
        if(ixx.gt.0) then
           index = ixx+(iyy-1)*amr_grid_nx+(izz-1)*amr_grid_nx*amr_grid_ny
           do iv2=1,nv2
              do iv1=1,nv1
                 values(ipt,iv1,iv2) = func(iv1,iv2,index)
              enddo
           enddo
        else
           values(ipt,1:nv1,1:nv2) = 0.d0
        endif
     endif
  enddo
  !
end subroutine amr_sample2_3d


!--------------------------------------------------------------------------
!                     FIND R, THETA, PHI FROM X, Y, Z
!--------------------------------------------------------------------------
subroutine amr_xyz_to_rthphi(x,y,z,r,theta,phi)
  implicit none
  double precision :: x,y,z,r,theta,phi
  double precision :: r2
  !
  r2    = x*x + y*y + z*z
  r     = sqrt(r2)
  if(z.eq.0.d0) then
     theta=pihalf
  else
     theta = acos(z/(r+1d-199))
  endif
  if(x.eq.0.d0) then
     if(y.ge.0.d0) then
        phi   = pihalf
     else
        phi   = 3.d0*pihalf 
     endif
  else
   phi   = atan(y/x)
   if(x.gt.0.d0) then
      if(y.lt.0.d0) then
         phi   = phi + twopi
      endif
   else
      phi   = phi + pi
   endif
endif
!################################
if((phi.lt.0.d0).or.(phi.gt.twopi)) stop 1319
!################################
end subroutine amr_xyz_to_rthphi


! !--------------------------------------------------------------------------
! !                     COMPUTE MINIMUM CELL SIZE
! !--------------------------------------------------------------------------
! subroutine amr_min_cellsize(b,dxmin)
! implicit none
! type(amr_branch), pointer :: b
! double precision :: dxmin,x0,x1,r0,t0,t1,sint
! integer :: iddr
! !
! if(amr_coordsystem.lt.100) then
!    !
!    ! Cartesian coordinates
!    !
!    dxmin = 1.d90
!    do iddr=1,3
!       x0 = amr_finegrid_xi(b%ixyzf(iddr),iddr,b%level)
!       x1 = amr_finegrid_xi(b%ixyzf(iddr)+1,iddr,b%level)
!       dxmin = min(dxmin,abs(x1-x0))
!    enddo
! elseif(amr_coordsystem.lt.200) then
!    !
!    ! Spherical coordinates
!    !
!    dxmin = 1.d90
!    r0 = amr_finegrid_xi(b%ixyzf(1),1,b%level)
!    x1 = amr_finegrid_xi(b%ixyzf(1)+1,1,b%level)
!    dxmin = min(dxmin,abs(x1-r0))
!    t0 = amr_finegrid_xi(b%ixyzf(2),2,b%level)
!    t1 = amr_finegrid_xi(b%ixyzf(2)+1,2,b%level)
!    dxmin = min(dxmin,r0*abs(t1-t0))
!    ! Note: We take the max of sin(Theta), because otherwise we get dxmin=0 at the pole
!    x0 = amr_finegrid_xi(b%ixyzf(3),3,b%level)
!    x1 = amr_finegrid_xi(b%ixyzf(3)+1,3,b%level)
!    if(t1.lt.pihalf) then
!       sint = sin(t1)
!    elseif(t0.gt.pihalf) then
!       sint = sin(t0)
!    else
!       sint = 1.d0
!    endif
!    dxmin = min(dxmin,r0*sint*abs(x1-x0))
! else
!    stop
! endif
! end subroutine amr_min_cellsize


!--------------------------------------------------------------------------
!                Hunt routine of Numerical Recipes
!--------------------------------------------------------------------------
SUBROUTINE amrhunt(xx,n,x,jlo)
  INTEGER jlo,n
  doubleprecision x,xx(n)
  INTEGER inc,jhi,jm
  LOGICAL ascnd
  ascnd=xx(n).gt.xx(1)
  if(jlo.le.0.or.jlo.gt.n)then
     jlo=0
     jhi=n+1
     goto 3
  endif
  inc=1
  if(x.ge.xx(jlo).eqv.ascnd)then
1    jhi=jlo+inc
     if(jhi.gt.n)then
        jhi=n+1
     else if(x.ge.xx(jhi).eqv.ascnd)then
        jlo=jhi
        inc=inc+inc
        goto 1
     endif
  else
     jhi=jlo
2    jlo=jhi-inc
     if(jlo.lt.1)then
        jlo=0
     else if(x.lt.xx(jlo).eqv.ascnd)then
        jhi=jlo
        inc=inc+inc
        goto 2
     endif
  endif
3 if(jhi-jlo.eq.1)return
  jm=(jhi+jlo)/2
  if(x.gt.xx(jm).eqv.ascnd)then
     jlo=jm
  else
     jhi=jm
  endif
  goto 3
END SUBROUTINE amrhunt
! C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


!==========================================================================
!                        SUPPORT FOR LAYERS
!
! Layers are a different way of organizing AMR refinement. Instead of
! doing it hierarchically, it is done by layers. But internally these
! layers are translated into the oct-tree hierarchical AMR. It is only
! a useful and easy way to import grids.
!==========================================================================

!--------------------------------------------------------------------------
!                          INSTALL THE LAYERS
!--------------------------------------------------------------------------
subroutine amr_layers_implement_amr_grid()
implicit none
integer :: ierr,ix,iy,iz,ixx,iyy,izz,ilyr,ilayer,indexp
type(amr_branch), pointer :: b
!
! First allocate all the maps
!
allocate(amr_layers_map(0)%index(1:amr_layers_nnxyz(1,0),&
                                 1:amr_layers_nnxyz(2,0),&
                                 1:amr_layers_nnxyz(3,0)),STAT=ierr)
if(ierr.ne.0) then
   write(stdo,*) 'ERROR in layers: could not allocate map'
   stop
endif
do ilayer=1,amr_layers_nrtot
   allocate(amr_layers_map(ilayer)%index(1:amr_layers_nnxyz(1,ilayer),&
                                         1:amr_layers_nnxyz(2,ilayer),&
                                         1:amr_layers_nnxyz(3,ilayer)),STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR in layers: could not allocate map'
      stop
   endif
enddo
!
! Now fill the base map 
!
do iz=1,amr_grid_nz
   do iy=1,amr_grid_ny
      do ix=1,amr_grid_nx
         b => amr_grid_branch(ix,iy,iz)%link
         if(.not.associated(b)) then
            write(stdo,*) 'ERROR in layers: Not all grid cells are filled'
            stop 2190
         endif
         amr_layers_map(0)%index(ix,iy,iz) = b%leafindex
      enddo
   enddo
enddo
!
! Now install each layer
!
do ilayer=1,amr_layers_nrtot
   !
   ! Get the parent
   !
   ilyr = amr_layers_parent(ilayer)
   !
   ! Now go through the patch of parent map that this current layer overlaps
   !
   do iz=0,amr_layers_nxyz(3,ilayer)-1
      do iy=0,amr_layers_nxyz(2,ilayer)-1
         do ix=0,amr_layers_nxyz(1,ilayer)-1
            !
            ! Find the current parent cell
            !
            indexp = amr_layers_map(ilyr)%index(amr_layers_ixyz(1,ilayer)+ix, &
                                                amr_layers_ixyz(2,ilayer)+iy, &
                                                amr_layers_ixyz(3,ilayer)+iz)
            if(indexp.eq.0) then
               write(stdo,*) 'ERROR in layers: cell already has zero index...'
               stop
            endif
            b => amr_index_to_leaf(indexp)%link
            !
            ! Now refine the current parent cell
            !
            call amr_branch_refine(b,0)
            !
            ! Unlink the indexp of the parent
            !
            amr_layers_map(ilyr)%index(amr_layers_ixyz(1,ilayer)+ix, &
                                       amr_layers_ixyz(2,ilayer)+iy, &
                                       amr_layers_ixyz(3,ilayer)+iz) = 0
            !
            ! Link the indices of the child cells
            !
            if(.not.associated(b%child)) then
               write(stdo,*) 'ERROR in AMR Module: branch without children found...'
               stop 7901
            endif
            do izz=0,amr_zdim
               do iyy=0,amr_ydim
                  do ixx=0,amr_xdim
                     amr_layers_map(ilayer)%                         &
                          index(2*ix+ixx+1,2*iy+iyy+1,2*iz+izz+1) =  &
                          b%child(ixx+1,iyy+1,izz+1)%link%leafindex
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
enddo
!
end subroutine amr_layers_implement_amr_grid


!--------------------------------------------------------------------------
!                 RESET NEXT CELL FINDING FOR AMR
! 
! See amr_nextcell() below for more information.
!--------------------------------------------------------------------------
subroutine amr_resetcount()
implicit none
if((amr_style.ge.0).and.(amr_style.lt.10)) then
   !
   ! For normal oct-tree AMR or regular grids
   !
   amr_octtree_count       = 0
   !
elseif((amr_style.ge.10).and.(amr_style.lt.20)) then
   !
   ! For layer-style AMR
   !
   amr_layers_count_ilayer = 0
   amr_layers_count_ix     = 0
   amr_layers_count_iy     = 1
   amr_layers_count_iz     = 1
   !
else
   write(stdo,*) 'ERROR: amr_style ',amr_style,' not known.'
   stop 7209
endif
end subroutine amr_resetcount


!--------------------------------------------------------------------------
!                     FIND NEXT CELL FOR AMR
!
! This subroutine is useful for input and output of data to a file. It
! will automatically find the index of the next cell to read or write to
! file, and does this in such a way that it also goes right when layers
! are used. The reason we need this method is that the layer-style AMR
! is only an I/O thing. Internally RADMC-3D converts this to an oct-tree.
! To easily I/O stuff from the oct-tree-which-is-actually-a-layer-set you
! make a loop:
!   call amr_resetcount()
!   do ileaf=1,amr_nrleafs
!      call amr_nextcell(index)
!      *** I/O stuff to/from cell with index "index" ***
!   enddo
! This will work for regular grids, for oct-tree grids and for layer-style
! grids. 
! NOTE: Do not use amr_nextcell() for internal data processing. It is only
!       meant for I/O, in order to allow for layers (while being also 
!       consistent with oct-tree and regular).
!--------------------------------------------------------------------------
subroutine amr_nextcell(index)
implicit none
integer :: index
!
! Increase counters
!
if((amr_style.ge.0).and.(amr_style.lt.10)) then
   if(amr_tree_present) then
      !
      ! For normal oct-tree AMR
      !
      amr_octtree_count       = amr_octtree_count + 1
      if(amr_octtree_count.le.amr_nrleafs) then
         index                = amr_theleaf_index(amr_octtree_count)
      else
         index = -1
      endif
   else
      !
      ! For regular grid without AMR tree
      !
      amr_octtree_count       = amr_octtree_count + 1
      if(amr_octtree_count.le.amr_nrleafs) then
         index                = amr_octtree_count
      else
         index = -1
      endif
   endif
   !
elseif((amr_style.ge.10).and.(amr_style.lt.20)) then
   !
   ! For layer-style AMR
   !
   amr_layers_count_ix = amr_layers_count_ix + 1
   if(amr_layers_count_ix.gt.amr_layers_nnxyz(1,amr_layers_count_ilayer)) then
      amr_layers_count_iy = amr_layers_count_iy + 1
      amr_layers_count_ix = 1
      if(amr_layers_count_iy.gt.amr_layers_nnxyz(2,amr_layers_count_ilayer)) then
         amr_layers_count_iz = amr_layers_count_iz + 1
         amr_layers_count_iy = 1
         amr_layers_count_ix = 1
         if(amr_layers_count_iz.gt.amr_layers_nnxyz(3,amr_layers_count_ilayer)) then
            amr_layers_count_ilayer = amr_layers_count_ilayer + 1
            amr_layers_count_iz = 1
            amr_layers_count_iy = 1
            amr_layers_count_ix = 1
            if(amr_layers_count_ilayer.gt.amr_layers_nrtot) then
               index = -1
               return
            endif
         endif
      endif
   endif
   index = amr_layers_map(amr_layers_count_ilayer)%index(amr_layers_count_ix, &
              amr_layers_count_iy,amr_layers_count_iz)
   !
else
   write(stdo,*) 'ERROR: amr_style ',amr_style,' not known.'
   stop 7209
endif
end subroutine amr_nextcell

!==========================================================================
!                 ROUTINES FOR CORNER-BASED ALGORITHMS
!==========================================================================


!--------------------------------------------------------------------------
!                       CREATE VERTEX ARRAYS
!
! For the moment the only way to use the vertices (=grid cell corner 
! points) is in a static (= non-adaptive) way. You first have the full
! AMR tree set up, and then you find all corner points and make a large
! array of them. This assigns a unique index number to each vertex: the
! vertexindex. Each cell (=leaf of the AMR tree) will have its own 2x2x2
! array of vertexindex numbers, so that it can quickly find its corner
! vertices. Each vertex will have its own 2x2x2 array of pointers to the
! 8 cells connecting to this vertex. In this way you can quickly find
! neighboring vertices, given one vertex: first go to the 2x2x2 cells,
! and then go to their 2x2x2-1 other vertices.
!
!--------------------------------------------------------------------------
subroutine amr_set_up_vertices(nvmax)
implicit none
integer, optional :: nvmax
integer :: ix,iy,iz,ivt
type(amr_branch), pointer :: b
logical :: abort
!
! Check
!
if(amr_tree_present) then
   if(.not.allocated(amr_grid_branch)) then
      write(stdo,*) 'ERROR: Can only make vertices if grid is set up.'
      stop
   endif
else
   if(amr_grid_nx.le.0) then
      write(stdo,*) 'ERROR: Can only make vertices if grid is set up.'
      stop
   endif
endif
!
! Allocate the amr_cell_corners array
! Note: Only needed if grid uses AMR tree
!
if(amr_tree_present) then
   if(allocated(amr_cell_corners)) deallocate(amr_cell_corners)
   allocate(amr_cell_corners(1:2,1:2,1:2,1:amr_nrleafs_max))
endif
!
! First count how many vertices we will have 
!
! We must make a first guess how many they could maximally be
!
if(present(nvmax)) then
   !
   ! The called knows exactly how many vertices he wants
   !
   amr_nr_vertices_max = nvmax
else
   !
   ! No nr of vertices specified, so we have to make a guess
   !
   amr_nr_vertices_max = max(2*amr_nrleafs_max,(amr_grid_nx+1)*   &
                             (amr_grid_ny+1)*(amr_grid_nz+1))
endif
!
! Entry point for if we have to redo the counting with larger
! amr_nr_vertices_max
!
50 continue
!
! Set up the arrays that connect the corners to the cells and
! vice versa. 
!
! Note: Only necessary if AMR tree used. For regular grids it is trivial
!       to calculate on-the-fly which cell connects to which vertex.
! 
if(amr_tree_present) then
   !
   ! Yes, the AMR tree is present, so...
   !
   ! Clear the cell corner array
   !
   amr_cell_corners(:,:,:,:) = 0
   !
   ! Deallocate vertex arrays
   !
   !if(allocated(amr_vertex_cell_is_corner)) deallocate(amr_vertex_cell_is_corner)
   if(allocated(amr_vertex_cells)) deallocate(amr_vertex_cells)
   !
   ! Allocate them, and clear
   !
   allocate(amr_vertex_cells(1:2,1:2,1:2,1:amr_nr_vertices_max))
   do ivt=1,amr_nr_vertices_max 
      do iz=1,1+amr_zdim
         do iy=1,1+amr_ydim
            do ix=1,1+amr_xdim
               nullify(amr_vertex_cells(ix,iy,iz,ivt)%link)
            enddo
         enddo
      enddo
   enddo
   !allocate(amr_vertex_cell_is_corner(1:amr_nr_vertices_max))
   !amr_vertex_cell_is_corner(:) = 0
   !
   ! Now start installing the vertices
   !
   abort = .false.
   amr_nr_vertices=0
   do iz=1,amr_grid_nz
      do iy=1,amr_grid_ny
         do ix=1,amr_grid_nx
            b => amr_grid_branch(ix,iy,iz)%link
            call amr_install_vertices_of_cell(b,abort)
            if(abort) then
               if(present(nvmax)) then
                  write(stdo,*) 'ERROR: Nr of vertices allocated is too small'
                  stop
               else
                  amr_nr_vertices_max = 2*amr_nr_vertices_max
                  goto 50
               endif
            endif
         enddo
      enddo
   enddo
   !
   ! If nvmax not set, then we redo the entire thing, but now with exactly
   ! the right number of vertices
   !
   ! NOTE: This is the poor-man's version of doing this. Better would be to
   !       have a proper algorithm for computing exactly the number needed.
   !
   if(.not.present(nvmax).and.(amr_nr_vertices_max.gt.amr_nr_vertices*1.3)) then
      !
      ! Set the precise nr of vertices
      !
      amr_nr_vertices_max = amr_nr_vertices
      !
      ! Clear the cell corner array
      !
      amr_cell_corners(:,:,:,:) = 0
      !
      ! Deallocate vertex arrays
      !
!      if(allocated(amr_vertex_cell_is_corner)) deallocate(amr_vertex_cell_is_corner)
      if(allocated(amr_vertex_cells)) deallocate(amr_vertex_cells)
      !
      ! Allocate them, and clear
      !
      allocate(amr_vertex_cells(1:2,1:2,1:2,1:amr_nr_vertices_max))
      do ivt=1,amr_nr_vertices_max 
         do iz=1,1+amr_zdim
            do iy=1,1+amr_ydim
               do ix=1,1+amr_xdim
                  nullify(amr_vertex_cells(ix,iy,iz,ivt)%link)
               enddo
            enddo
         enddo
      enddo
!      allocate(amr_vertex_cell_is_corner(1:amr_nr_vertices_max))
!      amr_vertex_cell_is_corner(:) = 0
      !
      ! Now start installing the vertices
      !
      abort = .false.
      amr_nr_vertices=0
      do iz=1,amr_grid_nz
         do iy=1,amr_grid_ny
            do ix=1,amr_grid_nx
               b => amr_grid_branch(ix,iy,iz)%link
               call amr_install_vertices_of_cell(b,abort)
               if(abort) then
                  write(stdo,*) 'ERROR: Strange: Nr of vertices allocated too small'
                  write(stdo,*) '       But I thought I counted them correctly... Strange.'
                  stop
               endif
            enddo
         enddo
      enddo
   endif
else
   !
   ! If we have a regular grid, then we still must set the amr_nr_vertices and
   ! amr_nr_vertices_max
   !
   ! Note: In case of cyclic boundary conditions we will have vertices at both
   !       boundaries that are, in principle, identical. They are therefore 
   !       double. If you use the AMR tree such identical vertical do not
   !       occur because you can always link them properly. 
   !
   amr_nr_vertices_max = (amr_grid_nx+amr_xdim)*(amr_grid_ny+amr_ydim)*(amr_grid_nz+amr_zdim)
   amr_nr_vertices     = amr_nr_vertices_max
endif
!
end subroutine amr_set_up_vertices

!--------------------------------------------------------------------------
!                     INSTALL VERTICES OF A CELL
!--------------------------------------------------------------------------
recursive subroutine amr_install_vertices_of_cell(b,abort)
implicit none
type(amr_branch), pointer :: b,a
integer :: ix,iy,iz,vertex(1:2,1:2,1:2),level,ilr(1:3)
integer :: vindex
logical :: abort
!
abort = .false.
!
! If we have children, then recursively install their vertices instead
! (which automatically includes the vertices of this branch). Otherwise
! install vertices of this leaf.
!
if(associated(b%child)) then
   !
   ! We have children, so recursively install those
   !
   do iz=1,1+amr_zdim
      do iy=1,1+amr_ydim
         do ix=1,1+amr_xdim
            if(.not.associated(b%child(ix,iy,iz)%link)) then
               write(stdo,*) 'ERROR: Child lost from branch...'
               stop
            endif
            a => b%child(ix,iy,iz)%link
            call amr_install_vertices_of_cell(a,abort)
            if(abort) return
         enddo
      enddo
   enddo
else
   !
   ! We are a leaf, so install all 2x2x2 vertices 
   !
   do iz=1,1+amr_zdim
      do iy=1,1+amr_ydim
         do ix=1,1+amr_xdim
            !
            ! Get the index of this corner; if it is 0, then the corner
            ! vertex does not yet exist
            !
            ilr(1) = ix
            ilr(2) = iy
            ilr(3) = iz
            call amr_get_corner_vertexindex(b,ilr,vindex,.true.)
            if(vindex.eq.0) then
               !
               ! Yes, we have a new vertex
               !
               amr_nr_vertices = amr_nr_vertices + 1
               !
               ! Check that we do not exceed the limit. If we do, then
               ! go back to 
               !
               if(amr_nr_vertices.gt.amr_nr_vertices_max) then
                  abort = .true.
                  return
               endif
               !
               ! Install this new vertex and assure that all cells
               ! that touch this vertex point to it.
               !
               vindex = amr_nr_vertices
               call amr_install_vertex(b,ilr,vindex,.true.)
               !
            endif
            !
            ! Now add this vindex to the amr_cell_corners
            !
            amr_cell_corners(ix,iy,iz,b%leafindex) = vindex
            !
         enddo
      enddo
   enddo
endif
!
end subroutine amr_install_vertices_of_cell

!--------------------------------------------------------------------------
!                       INSTALL A SINGLE VERTEX
!--------------------------------------------------------------------------
subroutine amr_install_vertex(b,ilr,vindex,finddone)
implicit none
type(amr_branch), pointer :: b,a
integer :: ilr(1:3),vindex
integer :: ix,iy,iz
logical :: finddone
!
! First identify all cells cornering on this vertex. This is written
! in the amr_loc_corner_cells() array.
!
if(.not.finddone) then 
   call amr_find_corner_cells(b,ilr,.true.)
endif
!
! Install the cell pointers 
!
do iz=1,1+amr_zdim
   do iy=1,1+amr_ydim
      do ix=1,1+amr_xdim
         amr_vertex_cells(ix,iy,iz,vindex)%link =>  &
              amr_loc_corner_cells(ix,iy,iz)%link
      enddo
   enddo
enddo
!
! Install also the "corner" flags
!
!amr_vertex_cell_is_corner(vindex) = amr_loc_cell_is_corner
!
! Now make sure that all the cells touching this corner actually
! point to it, IFF they have this vertex as a true corner.
!
do iz=1,1+amr_zdim
   do iy=1,1+amr_ydim
      do ix=1,1+amr_xdim
         if(associated(amr_loc_corner_cells(ix,iy,iz)%link)) then
            a => amr_loc_corner_cells(ix,iy,iz)%link
            if(a%leafindex.eq.0) then
               write(stdo,*) 'ERROR while preparing vertex-based grid: '
               write(stdo,*) '      found non-leaf branch on corner...'
               stop
            endif
            if(amr_loc_corner(ix,iy,iz)) then
               amr_cell_corners(2+amr_xdim-ix,2+amr_ydim-iy,&
                                2+amr_zdim-iz,a%leafindex) = vindex
            endif
         endif
      enddo
   enddo
enddo
!
end subroutine amr_install_vertex

!--------------------------------------------------------------------------
!                         FIND CORNER INDEX
!--------------------------------------------------------------------------
subroutine amr_get_corner_vertexindex(b,ilr,vindex,check)
implicit none
type(amr_branch), pointer :: b,a
integer :: ilr(1:3),vindex
integer :: ix,iy,iz,vidx
logical :: check
double precision :: x,y,z
!
! First identify all cells cornering on this vertex. This is written
! in the amr_loc_corner_cells() array.
!
call amr_find_corner_cells(b,ilr,.true.)
!
! Now check all the cells touching this corner
!
vindex=0
do iz=1,1+amr_zdim
   do iy=1,1+amr_ydim
      do ix=1,1+amr_xdim
         if(associated(amr_loc_corner_cells(ix,iy,iz)%link)) then
            a => amr_loc_corner_cells(ix,iy,iz)%link
            if(a%leafindex.eq.0) then
               write(stdo,*) 'ERROR while preparing vertex-based grid: '
               write(stdo,*) '      found non-leaf branch on corner...'
               stop
            endif
            !
            ! Only if this cell touches this vertex as its corner vertex
            !
            if(amr_loc_corner(ix,iy,iz)) then
               vidx = amr_cell_corners(2+amr_xdim-ix,2+amr_ydim-iy,&
                                       2+amr_zdim-iz,a%leafindex)
               !
               ! Check if this corner of this cell has already been obtained
               !
               if(vidx.gt.0) then
                  !
                  ! Check if this is the first non-zero vidx we have
                  !
                  if(vindex.eq.0) then
                     !
                     ! It's the first, so this is THE vindex
                     !
                     vindex = vidx
                     if(.not.check) return
                     !
                     ! For consistency check: get coordinates
                     ! *** Can be removed ***
                     !
                     x = amr_finegrid_xi(a%ixyzf(1)+1+amr_xdim-ix,1,a%level)
                     y = amr_finegrid_xi(a%ixyzf(2)+1+amr_ydim-iy,2,a%level)
                     z = amr_finegrid_xi(a%ixyzf(3)+1+amr_zdim-iz,3,a%level)
                  else
                     !
                     ! It's not the first, so here we can do a consistency check
                     !
                     if(vidx.ne.vindex) then
                        write(stdo,*) 'ERROR: Corner vertexindex mismatch found'
                        write(stdo,*) ix,iy,iz,ilr,b%id,a%id
                        stop
                     endif
                     !
                     ! And another consistency check: the coordinates
                     ! *** Can be removed ***
                     !
                     if((x.ne.amr_finegrid_xi(a%ixyzf(1)+1+amr_xdim-ix,1,a%level)).or.&
                        (y.ne.amr_finegrid_xi(a%ixyzf(2)+1+amr_ydim-iy,2,a%level)).or.&
                        (z.ne.amr_finegrid_xi(a%ixyzf(3)+1+amr_zdim-iz,3,a%level))) then
                        if((.not.amr_cyclic_xyz(1)).and.(.not.amr_cyclic_xyz(2)).and. &
                             (.not.amr_cyclic_xyz(3))) then
                           !
                           ! A single test for all dimensions is OK
                           !
                           write(stdo,*) 'ERROR: Corner coordinate mismatch found:'
                           write(stdo,*) x,y,z
                           write(stdo,*) amr_finegrid_xi(a%ixyzf(1)+1+amr_xdim-ix,1,a%level), &
                                         amr_finegrid_xi(a%ixyzf(2)+1+amr_ydim-iy,2,a%level), &
                                         amr_finegrid_xi(a%ixyzf(3)+1+amr_zdim-iz,3,a%level)
                           stop
                        else
                           !
                           ! Possible cyclic BCs, so test individual dimensions
                           !
                           if(x.ne.amr_finegrid_xi(a%ixyzf(1)+1+amr_xdim-ix,1,a%level)) then
                              if(amr_cyclic_xyz(1)) then
                                 if(abs(abs(x-amr_finegrid_xi(a%ixyzf(1)+1+amr_xdim-ix,1,a%level))/          &
                                      abs(amr_grid_xi(amr_grid_nx+1,1)-amr_grid_xi(1,1))-1.d0).gt. &
                                      1.0d-6) then
                                    write(stdo,*) 'ERROR: Corner coordinate mismatch found:'
                                    write(stdo,*) x,y,z
                                    write(stdo,*) amr_finegrid_xi(a%ixyzf(1)+1+amr_xdim-ix,1,a%level), &
                                         amr_finegrid_xi(a%ixyzf(2)+1+amr_ydim-iy,2,a%level), &
                                         amr_finegrid_xi(a%ixyzf(3)+1+amr_zdim-iz,3,a%level)
                                    stop
                                 endif
                              else
                                 write(stdo,*) 'ERROR: Corner coordinate mismatch found:'
                                 write(stdo,*) x,y,z
                                 write(stdo,*) amr_finegrid_xi(a%ixyzf(1)+1+amr_xdim-ix,1,a%level), &
                                      amr_finegrid_xi(a%ixyzf(2)+1+amr_ydim-iy,2,a%level), &
                                      amr_finegrid_xi(a%ixyzf(3)+1+amr_zdim-iz,3,a%level)
                                 stop
                              endif
                           endif
                           if(y.ne.amr_finegrid_xi(a%ixyzf(2)+1+amr_ydim-iy,2,a%level)) then
                              if(amr_cyclic_xyz(2)) then
                                 if(abs(abs(y-amr_finegrid_xi(a%ixyzf(2)+1+amr_ydim-iy,2,a%level))/          &
                                      abs(amr_grid_xi(amr_grid_ny+1,2)-amr_grid_xi(1,2))-1.d0).gt. &
                                      1.0d-6) then
                                    write(stdo,*) 'ERROR: Corner coordinate mismatch found:'
                                    write(stdo,*) x,y,z
                                    write(stdo,*) amr_finegrid_xi(a%ixyzf(1)+1+amr_xdim-ix,1,a%level), &
                                         amr_finegrid_xi(a%ixyzf(2)+1+amr_ydim-iy,2,a%level), &
                                         amr_finegrid_xi(a%ixyzf(3)+1+amr_zdim-iz,3,a%level)
                                    stop
                                 endif
                              else
                                 write(stdo,*) 'ERROR: Corner coordinate mismatch found:'
                                 write(stdo,*) x,y,z
                                 write(stdo,*) amr_finegrid_xi(a%ixyzf(1)+1+amr_xdim-ix,1,a%level), &
                                      amr_finegrid_xi(a%ixyzf(2)+1+amr_ydim-iy,2,a%level), &
                                      amr_finegrid_xi(a%ixyzf(3)+1+amr_zdim-iz,3,a%level)
                                 stop
                              endif
                           endif
                           if(z.ne.amr_finegrid_xi(a%ixyzf(3)+1+amr_zdim-iz,3,a%level)) then
                              if(amr_cyclic_xyz(3)) then
                                 if(abs(abs(z-amr_finegrid_xi(a%ixyzf(3)+1+amr_zdim-iz,3,a%level))/          &
                                      abs(amr_grid_xi(amr_grid_nz+1,3)-amr_grid_xi(1,3))-1.d0).gt. &
                                      1.0d-6) then
                                    write(stdo,*) 'ERROR: Corner coordinate mismatch found:'
                                    write(stdo,*) x,y,z
                                    write(stdo,*) amr_finegrid_xi(a%ixyzf(1)+1+amr_xdim-ix,1,a%level), &
                                         amr_finegrid_xi(a%ixyzf(2)+1+amr_ydim-iy,2,a%level), &
                                         amr_finegrid_xi(a%ixyzf(3)+1+amr_zdim-iz,3,a%level)
                                    stop
                                 endif
                              else
                                 write(stdo,*) 'ERROR: Corner coordinate mismatch found:'
                                 write(stdo,*) x,y,z
                                 write(stdo,*) amr_finegrid_xi(a%ixyzf(1)+1+amr_xdim-ix,1,a%level), &
                                      amr_finegrid_xi(a%ixyzf(2)+1+amr_ydim-iy,2,a%level), &
                                      amr_finegrid_xi(a%ixyzf(3)+1+amr_zdim-iz,3,a%level)
                                 stop
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
      enddo
   enddo
enddo
!
end subroutine amr_get_corner_vertexindex

!--------------------------------------------------------------------------
!                IDENTIFY ALL 8 CELLS OF A CORNER-VERTEX
!--------------------------------------------------------------------------
subroutine amr_find_corner_cells(b,ilr,reduce)
implicit none
type(amr_branch), pointer :: b
logical :: reduce
integer :: ilr(1:3),inv(1:3)
integer :: ix,iy,iz,levelb,id,ibit,lev,lx,ly,lz,irib(1:3),irdir,iii,irf
!
! Do a check
!
if(.not.b%leaf) then
   write(stdo,*) 'ERROR: Can only find corner cells if b is a leaf.'
   stop
endif
!
! Get the level of cell b
!
levelb = b%level
!
! Find the opposite of ilr
!
inv(1) = 2+amr_xdim-ilr(1)
inv(2) = 2+amr_ydim-ilr(2)
inv(3) = 2+amr_zdim-ilr(3)
!
! The current branch is of course one of them
!
amr_loc_corner_cells(inv(1),inv(2),inv(3))%link => b
amr_loc_corner_celllevel(inv(1),inv(2),inv(3)) = b%level
!
! First identify the neighbors of the current cell/branch in the
! directions ilr
!
if(amr_xdim.eq.1) then
   amr_loc_corner_cells(ilr(1),inv(2),inv(3))%link => b%neighbor(ilr(1),1)%link
   if(associated(amr_loc_corner_cells(ilr(1),inv(2),inv(3))%link)) then
      amr_loc_corner_celllevel(ilr(1),inv(2),inv(3)) = b%neighbor(ilr(1),1)%link%level
   else
      amr_loc_corner_celllevel(ilr(1),inv(2),inv(3)) = 9999
   endif
endif
if(amr_ydim.eq.1) then
   amr_loc_corner_cells(inv(1),ilr(2),inv(3))%link => b%neighbor(ilr(2),2)%link
   if(associated(amr_loc_corner_cells(inv(1),ilr(2),inv(3))%link)) then
      amr_loc_corner_celllevel(inv(1),ilr(2),inv(3)) = b%neighbor(ilr(2),2)%link%level
   else
      amr_loc_corner_celllevel(inv(1),ilr(2),inv(3)) = 9999
   endif
endif
if(amr_zdim.eq.1) then
   amr_loc_corner_cells(inv(1),inv(2),ilr(3))%link => b%neighbor(ilr(3),3)%link
   if(associated(amr_loc_corner_cells(inv(1),inv(2),ilr(3))%link)) then
      amr_loc_corner_celllevel(inv(1),inv(2),ilr(3)) = b%neighbor(ilr(3),3)%link%level
   else
      amr_loc_corner_celllevel(inv(1),inv(2),ilr(3)) = 9999
   endif
endif
!
! So now if one of these 3 cells is not existent, its level is 9999. So
! this makes it easier to test for it.
!
! Now find the 3 cells between the above 3
! If the above 3 cells have a level that is smaller than that of cell b,
! then we must be very careful, because then two of the 2x2x2 slots may
! in fact be the same (larger) cell!
!
! ...First the neighbor of the x and y pair
!
if((amr_xdim.eq.1).and.(amr_ydim.eq.1)) then
   if((amr_loc_corner_celllevel(ilr(1),inv(2),inv(3)).ge.    &
       amr_loc_corner_celllevel(inv(1),ilr(2),inv(3))).and. &
      (amr_loc_corner_celllevel(ilr(1),inv(2),inv(3)).ne.9999)) then
      !
      ! The ilr(1),inv(2),inv(3) cell has the smallest (or equal) 
      ! size of the two. So we use him to find the neighbor.
      !
      amr_loc_corner_cells(ilr(1),ilr(2),inv(3))%link => &
           amr_loc_corner_cells(ilr(1),inv(2),inv(3))%link%neighbor(ilr(2),2)%link
   elseif(amr_loc_corner_celllevel(inv(1),ilr(2),inv(3)).ne.9999) then
      !
      ! We use instead inv(1),ilr(2),inv(3) to find the neighbor.
      !
      amr_loc_corner_cells(ilr(1),ilr(2),inv(3))%link => &
           amr_loc_corner_cells(inv(1),ilr(2),inv(3))%link%neighbor(ilr(1),1)%link
   else
      nullify(amr_loc_corner_cells(ilr(1),ilr(2),inv(3))%link)
   endif
endif
if(associated(amr_loc_corner_cells(ilr(1),ilr(2),inv(3))%link)) then
   amr_loc_corner_celllevel(ilr(1),ilr(2),inv(3)) = & 
        amr_loc_corner_cells(ilr(1),ilr(2),inv(3))%link%level
else
   amr_loc_corner_celllevel(ilr(1),ilr(2),inv(3)) = 9999
endif
!
! ...Then the neighbor of the x and z pair
!
if((amr_xdim.eq.1).and.(amr_zdim.eq.1)) then
   if((amr_loc_corner_celllevel(ilr(1),inv(2),inv(3)).ge.    &
       amr_loc_corner_celllevel(inv(1),inv(2),ilr(3))).and. &
      (amr_loc_corner_celllevel(ilr(1),inv(2),inv(3)).ne.9999)) then
      !
      ! The ilr(1),inv(2),inv(3) cell has the smallest (or equal)
      ! size of the two. So we use him to find the neighbor.
      !
      amr_loc_corner_cells(ilr(1),inv(2),ilr(3))%link => &
           amr_loc_corner_cells(ilr(1),inv(2),inv(3))%link%neighbor(ilr(3),3)%link
   elseif(amr_loc_corner_celllevel(inv(1),inv(2),ilr(3)).ne.9999) then
      !
      ! We use inv(1),inv(2),ilr(3) to find the neighbor.
      !
      amr_loc_corner_cells(ilr(1),inv(2),ilr(3))%link => &
           amr_loc_corner_cells(inv(1),inv(2),ilr(3))%link%neighbor(ilr(1),1)%link
   else
      nullify(amr_loc_corner_cells(ilr(1),inv(2),ilr(3))%link)
   endif
endif
if(associated(amr_loc_corner_cells(ilr(1),inv(2),ilr(3))%link)) then
   amr_loc_corner_celllevel(ilr(1),inv(2),ilr(3)) = & 
        amr_loc_corner_cells(ilr(1),inv(2),ilr(3))%link%level
else
   amr_loc_corner_celllevel(ilr(1),inv(2),ilr(3)) = 9999
endif
!
! ...Finally the neighbor of the y and z pair
!
if((amr_ydim.eq.1).and.(amr_zdim.eq.1)) then
   if((amr_loc_corner_celllevel(inv(1),ilr(2),inv(3)).ge.    &
       amr_loc_corner_celllevel(inv(1),inv(2),ilr(3))).and. & 
      (amr_loc_corner_celllevel(inv(1),ilr(2),inv(3)).ne.9999)) then
      !
      ! The inv(1),ilr(2),inv(3) cell has the smallest (or equal)
      ! size of the two. So we use him to find the neighbor.
      !
      amr_loc_corner_cells(inv(1),ilr(2),ilr(3))%link => &
           amr_loc_corner_cells(inv(1),ilr(2),inv(3))%link%neighbor(ilr(3),3)%link
   elseif(amr_loc_corner_celllevel(inv(1),inv(2),ilr(3)).ne.9999) then
      !
      ! We take inv(1),inv(2),ilr(3) to find the neighbor.
      !
      amr_loc_corner_cells(inv(1),ilr(2),ilr(3))%link => &
           amr_loc_corner_cells(inv(1),inv(2),ilr(3))%link%neighbor(ilr(2),2)%link
   else
      nullify(amr_loc_corner_cells(inv(1),ilr(2),ilr(3))%link)
   endif
endif
if(associated(amr_loc_corner_cells(inv(1),ilr(2),ilr(3))%link)) then
   amr_loc_corner_celllevel(inv(1),ilr(2),ilr(3)) = & 
        amr_loc_corner_cells(inv(1),ilr(2),ilr(3))%link%level
else
   amr_loc_corner_celllevel(inv(1),ilr(2),ilr(3)) = 9999
endif
!
! Now find the one remaining cell diametrically opposite the b cell
! 
if((amr_xdim.eq.1).and.(amr_ydim.eq.1).and.(amr_zdim.eq.1)) then
   !
   ! First find the smallest neigbor of cell b
   !
   lx = amr_loc_corner_celllevel(ilr(1),inv(2),inv(3))
   ly = amr_loc_corner_celllevel(inv(1),ilr(2),inv(3))
   lz = amr_loc_corner_celllevel(inv(1),inv(2),ilr(3))
   if(lx.eq.9999) lx=-9999
   if(ly.eq.9999) ly=-9999
   if(lz.eq.9999) lz=-9999
   if(lx.gt.max(ly,lz)) then
      !
      ! Neighbor in x-direction is smallest of the three
      !
      ! For handling a special case (see later), set some
      ! variable
      !
      irdir = 1
      !
      ! Now find the smallest neighbor of this one 
      !
      ly = amr_loc_corner_celllevel(ilr(1),ilr(2),inv(3))
      lz = amr_loc_corner_celllevel(ilr(1),inv(2),ilr(3))
      if(ly.eq.9999) ly=-9999
      if(lz.eq.9999) lz=-9999
      if(ly.gt.lz) then
         amr_loc_corner_cells(ilr(1),ilr(2),ilr(3))%link => &
            amr_loc_corner_cells(ilr(1),ilr(2),inv(3))%link%neighbor(ilr(3),3)%link
      elseif(lz.ne.-9999) then
         amr_loc_corner_cells(ilr(1),ilr(2),ilr(3))%link => &
            amr_loc_corner_cells(ilr(1),inv(2),ilr(3))%link%neighbor(ilr(2),2)%link
      else
         nullify(amr_loc_corner_cells(ilr(1),ilr(2),ilr(3))%link)
      endif
   elseif(ly.gt.lz) then
      !
      ! Neighbor in y-direction is smallest of the three
      !
      ! For handling a special case (see later), set some
      ! variable
      !
      irdir = 2
      !
      ! Now find the smallest neighbor of this one 
      !
      lx = amr_loc_corner_celllevel(ilr(1),ilr(2),inv(3))
      lz = amr_loc_corner_celllevel(inv(1),ilr(2),ilr(3))
      if(lx.eq.9999) lx=-9999
      if(lz.eq.9999) lz=-9999
      if(lx.gt.lz) then
         amr_loc_corner_cells(ilr(1),ilr(2),ilr(3))%link => &
            amr_loc_corner_cells(ilr(1),ilr(2),inv(3))%link%neighbor(ilr(3),3)%link
      elseif(lz.ne.-9999) then
         amr_loc_corner_cells(ilr(1),ilr(2),ilr(3))%link => &
            amr_loc_corner_cells(inv(1),ilr(2),ilr(3))%link%neighbor(ilr(1),1)%link
      else
         nullify(amr_loc_corner_cells(ilr(1),ilr(2),ilr(3))%link)
      endif
   elseif(lz.ne.-9999) then
      !
      ! Neighbor in z-direction is smallest (or equal) of the three
      !
      ! For handling a special case (see later), set some
      ! variable
      !
      irdir = 3
      !
      ! Now find the smallest neighbor of this one 
      !
      lx = amr_loc_corner_celllevel(ilr(1),inv(2),ilr(3))
      ly = amr_loc_corner_celllevel(inv(1),ilr(2),ilr(3))
      if(lx.eq.9999) lx=-9999
      if(ly.eq.9999) ly=-9999
      if(lx.gt.ly) then
         amr_loc_corner_cells(ilr(1),ilr(2),ilr(3))%link => &
            amr_loc_corner_cells(ilr(1),inv(2),ilr(3))%link%neighbor(ilr(2),2)%link
      elseif(ly.ne.-9999) then
         amr_loc_corner_cells(ilr(1),ilr(2),ilr(3))%link => &
            amr_loc_corner_cells(inv(1),ilr(2),ilr(3))%link%neighbor(ilr(1),1)%link
      else
         nullify(amr_loc_corner_cells(ilr(1),ilr(2),ilr(3))%link)
      endif
   else
      !
      ! No valid neighbors (we must be in the absolute corner of the grid)
      !
      nullify(amr_loc_corner_cells(ilr(1),ilr(2),ilr(3))%link)
   endif
endif
if(associated(amr_loc_corner_cells(ilr(1),ilr(2),ilr(3))%link)) then
   amr_loc_corner_celllevel(ilr(1),ilr(2),ilr(3)) = & 
        amr_loc_corner_cells(ilr(1),ilr(2),ilr(3))%link%level
else
   amr_loc_corner_celllevel(ilr(1),ilr(2),ilr(3)) = 9999
endif
!
! Now find out which of these cells have the present vertex as corner,
! because while they surely touch this vertex, it is not guaranteed that
! this vertex is indeed the cell's corner (in case of AMR refinement).
!
! Set all corner flags to .true. initially
!
amr_loc_corner(:,:,:) = .true.
!
! Surely cell b is a corner cell, so we don't have to do anything
! for this.
!
! Check the neighbors of b; Only do something if a cell is NOT a corner
!
if(amr_xdim.eq.1) then
   lev  = amr_loc_corner_celllevel(ilr(1),inv(2),inv(3))
   if(lev.lt.levelb) then
      id = amr_loc_corner_cells(ilr(1),inv(2),inv(3))%link%id
      if(amr_ydim.eq.1) then
         if(amr_loc_corner_celllevel(ilr(1),ilr(2),inv(3)).ne.9999) then
            if(amr_loc_corner_cells(ilr(1),ilr(2),inv(3))%link%id.eq.id) then
               amr_loc_corner(ilr(1),inv(2),inv(3)) = .false.
               amr_loc_corner(ilr(1),ilr(2),inv(3)) = .false.
            endif
         endif
      endif
      if(amr_zdim.eq.1) then
         if(amr_loc_corner_celllevel(ilr(1),inv(2),ilr(3)).ne.9999) then
            if(amr_loc_corner_cells(ilr(1),inv(2),ilr(3))%link%id.eq.id) then
               amr_loc_corner(ilr(1),inv(2),inv(3)) = .false.
               amr_loc_corner(ilr(1),inv(2),ilr(3)) = .false.
            endif
         endif
      endif
   endif
endif
if(amr_ydim.eq.1) then
   lev  = amr_loc_corner_celllevel(inv(1),ilr(2),inv(3))
   if(lev.lt.levelb) then
      id = amr_loc_corner_cells(inv(1),ilr(2),inv(3))%link%id
      if(amr_xdim.eq.1) then
         if(amr_loc_corner_celllevel(ilr(1),ilr(2),inv(3)).ne.9999) then
            if(amr_loc_corner_cells(ilr(1),ilr(2),inv(3))%link%id.eq.id) then
               amr_loc_corner(inv(1),ilr(2),inv(3)) = .false.
               amr_loc_corner(ilr(1),ilr(2),inv(3)) = .false.
            endif
         endif
      endif
      if(amr_zdim.eq.1) then
         if(amr_loc_corner_celllevel(inv(1),ilr(2),ilr(3)).ne.9999) then
            if(amr_loc_corner_cells(inv(1),ilr(2),ilr(3))%link%id.eq.id) then
               amr_loc_corner(inv(1),ilr(2),inv(3)) = .false.
               amr_loc_corner(inv(1),ilr(2),ilr(3)) = .false.
            endif
         endif
      endif
   endif
endif
if(amr_zdim.eq.1) then
   lev  = amr_loc_corner_celllevel(inv(1),inv(2),ilr(3))
   if(lev.lt.levelb) then
      id = amr_loc_corner_cells(inv(1),inv(2),ilr(3))%link%id
      if(amr_xdim.eq.1) then
         if(amr_loc_corner_celllevel(ilr(1),inv(2),ilr(3)).ne.9999) then
            if(amr_loc_corner_cells(ilr(1),inv(2),ilr(3))%link%id.eq.id) then
               amr_loc_corner(inv(1),inv(2),ilr(3)) = .false.
               amr_loc_corner(ilr(1),inv(2),ilr(3)) = .false.
            endif
         endif
      endif
      if(amr_ydim.eq.1) then
         if(amr_loc_corner_celllevel(inv(1),ilr(2),ilr(3)).ne.9999) then
            if(amr_loc_corner_cells(inv(1),ilr(2),ilr(3))%link%id.eq.id) then
               amr_loc_corner(inv(1),inv(2),ilr(3)) = .false.
               amr_loc_corner(inv(1),ilr(2),ilr(3)) = .false.
            endif
         endif
      endif
   endif
endif
if((amr_xdim.eq.1).and.(amr_ydim.eq.1).and.(amr_zdim.eq.1).and. &
     (amr_loc_corner_celllevel(ilr(1),ilr(2),ilr(3)).ne.9999)) then
   id = amr_loc_corner_cells(ilr(1),ilr(2),ilr(3))%link%id
   if(amr_loc_corner_celllevel(inv(1),ilr(2),ilr(3)).ne.9999) then
      if(amr_loc_corner_cells(inv(1),ilr(2),ilr(3))%link%id.eq.id) then
         amr_loc_corner(ilr(1),ilr(2),ilr(3)) = .false.
         amr_loc_corner(inv(1),ilr(2),ilr(3)) = .false.
      endif
   endif
   if(amr_loc_corner_celllevel(ilr(1),inv(2),ilr(3)).ne.9999) then
      if(amr_loc_corner_cells(ilr(1),inv(2),ilr(3))%link%id.eq.id) then
         amr_loc_corner(ilr(1),ilr(2),ilr(3)) = .false.
         amr_loc_corner(ilr(1),inv(2),ilr(3)) = .false.
      endif
   endif
   if(amr_loc_corner_celllevel(ilr(1),ilr(2),inv(3)).ne.9999) then
      if(amr_loc_corner_cells(ilr(1),ilr(2),inv(3))%link%id.eq.id) then
         amr_loc_corner(ilr(1),ilr(2),ilr(3)) = .false.
         amr_loc_corner(ilr(1),ilr(2),inv(3)) = .false.
      endif
   endif
endif
!
! Now reduce these cells to their smallest children, if they have children
! In this way the resulting branches are always assured to be the leafs.
! But do so only if this is requested by the caller.
!
if(reduce) then
   do iz=1,1+amr_zdim
      do iy=1,1+amr_ydim
         do ix=1,1+amr_xdim
            if(associated(amr_loc_corner_cells(ix,iy,iz)%link)) then
               if(associated(amr_loc_corner_cells(ix,iy,iz)%link%child)) then
                  !
                  ! This cell has children, so we must reduce. 
                  !
                  if(amr_loc_corner(ix,iy,iz)) then
                     !
                     ! In most cases the cell that is to be reduced is 
                     ! a corner cell. This makes the reduction easy, since
                     ! we always know which of the 2x2x2 children to choose.
                     !
                     do while(associated(amr_loc_corner_cells(ix,iy,iz)%link%child)) 
                        amr_loc_corner_cells(ix,iy,iz)%link =>                         &
                             amr_loc_corner_cells(ix,iy,iz)%link%                      &
                             child(2+amr_xdim-ix,2+amr_ydim-iy,2+amr_zdim-iz)%link
                        amr_loc_corner_celllevel(ix,iy,iz) = amr_loc_corner_celllevel(ix,iy,iz) + 1
                     enddo
                  else
                     !
                     ! In special cases we may encounter the situation 
                     ! that the cell-to-be-reduced is NOT a corner cell.
                     ! This is harder. It happens in the following case:
                     !   
                     !   +---+---+-------+     We start with cell A and
                     !   |   |   |       |     vertex "O".
                     !   +---+---+   B   |     Note: 3-D! So "O" is 
                     !   |   | A |       |     one small cell size up,
                     !   +---+---O---+---+     i.e. halfway the vertical     
                     !   |       |c1 |c2 |     ribbon of cell B. The
                     !   |       +---+---+     cell C (split into c1..c4)
                     !   |       |c3 |c4 |     is the cell-to-be-reduced
                     !   +-------+---+---+     into c1 (or c5 above c1)
                     !
                     ! This situation can only occur on a ribbon, not
                     ! on a face nor in a cell itself (if I am not 
                     ! mistaken). The reason is that you start with cell
                     ! "A", and then you go to direct neighbors (sharing a
                     ! face with "A"). These direct neighbors can only be
                     ! non-corner-cells (not having "O" as a corner) if
                     ! they are larger, in which case they have no children,
                     ! by virtue of the neighbor search algorithm. 
                     ! The problem shows up when we go from the neighbor,
                     ! (e.g. from cell "B") to the x-y, x-z and y-z neighbors.
                     ! If these are equally big or smaller than "A", they 
                     ! share only a ribbon with "A". If they are bigger
                     ! than "A", they could also share a face with "A".
                     ! But in that case they would have no children,
                     ! because if they had, then one of these children
                     ! is the one that shares a face with "A", and thus
                     ! the algorithm would choose that pathway to go
                     ! to the x-y, x-z or y-z cell in the first place.
                     ! So: the x-y, x-z or y-z cell that are bigger than 
                     ! "A" but have children (these are the cells that
                     ! can cause troubles) share ONLY a ribbon with "A". 
                     ! In other words: for the x-y, x-z or y-z cells, if
                     ! there is trouble, it can only happen on a ribbon.
                     ! 
                     ! For the final one, the x-y-z cell (opposite "A")
                     ! a similar argument can be made, so also there the
                     ! problematic vertex can only be on a ribbon.
                     !
                     ! So: Here we only have to deal with ribbons. 
                     !
                     ! Which of the 12 ribbons is it? 
                     !
                     ! NOTE: If we are in the diametrically opposite
                     !       cell, it is not clear what irdir is from
                     !       comparing ix,iy,iz with inv(*), but in 
                     !       that case irdir has been already set
                     !       above!
                     !
                     irib(1) = 2+amr_xdim-ix
                     irib(2) = 2+amr_ydim-iy
                     irib(3) = 2+amr_zdim-iz
                     if(ix.eq.inv(1)) then
                        irdir   = 1
                     elseif(iy.eq.inv(2)) then
                        irdir   = 2
                     elseif(iz.eq.inv(3)) then
                        irdir   = 3
                     endif
                     !
                     ! Now recursively reduce, while accounting for
                     ! the fact that the vertex is somewhere along
                     ! the ribbon, but not at the edges of the ribbon.
                     !
                     do while(associated(amr_loc_corner_cells(ix,iy,iz)%link%child).and. &
                          (.not.amr_loc_corner(ix,iy,iz)))
                        !
                        ! This cell is not a corner cell, but has children. 
                        ! This fact alone means that we will refine, and that the
                        ! celllevel of the cell increases.
                        !
                        amr_loc_corner_celllevel(ix,iy,iz) = amr_loc_corner_celllevel(ix,iy,iz) + 1
                        !
                        ! Now find the fine grid x (or y or z)- index of
                        ! this cell. To be more precise: The center of this
                        ! cell will be on an interface of cells of higher
                        ! refinement level. We want to know which interface
                        ! it lies on, at the refinement level of the
                        ! original cell b.
                        !
                        iii = 2*amr_loc_corner_cells(ix,iy,iz)%link%ixyzf(irdir)
                        do irf=2,levelb-amr_loc_corner_cells(ix,iy,iz)%link%level
                           iii = (iii-1)*2+1
                        enddo
                        !
                        ! Now we must find where on the ribbon the vertex
                        ! lies. We distinguish here two scenarios.
                        !
                        if(b%ixyzf(irdir)+ilr(irdir)-1.eq.iii) then
                           !
                           ! Scenario 1: Our vertex is exactly in the middle of
                           ! this ribbon. This means that the child cell will
                           ! have the vertex as its corner. We have to decide
                           ! which of the two children it is.
                           !
                           if(irdir.eq.1) then
                              irib(irdir) = ix
                           elseif(irdir.eq.2) then
                              irib(irdir) = iy
                           elseif(irdir.eq.3) then
                              irib(irdir) = iz
                           else
                              stop
                           endif
                           !
                           ! Now reduce
                           !
                           amr_loc_corner_cells(ix,iy,iz)%link =>                      &
                             amr_loc_corner_cells(ix,iy,iz)%link%                      &
                             child(irib(1),irib(2),irib(3))%link
                           !
                           ! The child has our vertex as its corner
                           !
                           amr_loc_corner(ix,iy,iz)=.true.
                           !
                           ! From this point on we (possibly) continue
                           ! to reduce, but now in the "normal" way
                           !
                        else
                           !
                           ! Scenario 2: Our vertex lies at 0.25, 0.75 or even
                           ! more refined ratios along the ribbon.
                           !
                           if(b%ixyzf(irdir)+ilr(irdir)-1.lt.iii) then
                              !
                              ! The vertex lies below the ribbon center,
                              ! so we take the lower cell.
                              !
                              irib(irdir) = 1
                           else
                              !
                              ! The vertex lies above the ribbon center,
                              ! so we take the upper cell.
                              !
                              irib(irdir) = 2
                           endif
                           !
                           ! Now reduce
                           !
                           amr_loc_corner_cells(ix,iy,iz)%link =>                     &
                                amr_loc_corner_cells(ix,iy,iz)%link%                  &
                                child(irib(1),irib(2),irib(3))%link
                        endif
                     enddo
                     !
                     ! If necessary, continue to reduce in the normal way, 
                     ! assuming that the cell has the vertex as its corner
                     !
                     do while(associated(amr_loc_corner_cells(ix,iy,iz)%link%child)) 
                        amr_loc_corner_cells(ix,iy,iz)%link =>                         &
                             amr_loc_corner_cells(ix,iy,iz)%link%                      &
                             child(2+amr_xdim-ix,2+amr_ydim-iy,2+amr_zdim-iz)%link
                        amr_loc_corner_celllevel(ix,iy,iz) = amr_loc_corner_celllevel(ix,iy,iz) + 1
                     enddo
                  endif
               endif
            endif
         enddo
      enddo
   enddo
endif
!
! Store the flags into a single integer using binary switches
!
!amr_loc_cell_is_corner = 0
!do iz=1,1+amr_zdim
!   do iy=1,1+amr_ydim
!      do ix=1,1+amr_xdim
!         if(amr_loc_corner(ix,iy,iz)) then
!            ibit = ix*(iy**2)*(iz**4)
!            amr_loc_cell_is_corner = amr_loc_cell_is_corner + ibit
!         endif
!      enddo
!   enddo
!enddo
!
end subroutine amr_find_corner_cells


!--------------------------------------------------------------------------
!                 CHECK CONSISTENCY OF VERTEX (SELF-TEST)
!--------------------------------------------------------------------------
subroutine amr_check_vertex(vidx,x,y,z)
  implicit none
  integer :: vidx,ix,iy,iz,ilevel,vidx1
  integer :: icx0,icy0,icz0
  integer :: icx,icy,icz
  double precision, optional :: x,y,z
  type(amr_branch), pointer :: a
  !
  ! Reset coordinate indices
  !
  icx0 = -1
  icy0 = -1
  icz0 = -1
  !
  ! Do a loop over all 8 possible cells that may have this 
  ! vertex as a corner cell
  !
  do iz=1,1+amr_zdim
     do iy=1,1+amr_ydim
        do ix=1,1+amr_xdim
           !
           ! Get this cell
           !
           a => amr_vertex_cells(ix,iy,iz,vidx)%link
           if(associated(a)) then
              !
              ! Yes, this is indeed a cell, so let's go
              !
              ! Check if the cell points back to this vertex.
              ! If not, then this vertex is not a corner point
              ! of the cell, but a point on a ribbon or face.
              !
              vidx1 = amr_cell_corners(2+amr_xdim-ix,2+amr_ydim-iy,&
                                       2+amr_zdim-iz,a%leafindex)
              if(vidx1.eq.vidx) then
                 ! 
                 ! Yes, this vertex is the corner of this
                 ! cell
                 !
                 ! Get the coordinates, but in index form
                 !
                 icx = a%ixyzf(1)+amr_xdim+1-ix
                 icy = a%ixyzf(2)+amr_ydim+1-iy
                 icz = a%ixyzf(3)+amr_zdim+1-iz
                 !
                 ! Refine down to the deepest level
                 !
                 do ilevel=a%level+1,amr_levelmax
                    icx = 2*(icx-1)+1
                    icy = 2*(icy-1)+1
                    icz = 2*(icz-1)+1
                 enddo
                 !
                 ! Compare
                 !
                 if(amr_xdim.ne.0) then
                    if(icx0.le.0) then
                       icx0 =icx
                    else
                       if(icx0.ne.icx) then
                          write(stdo,*) 'ERROR: Inconsistent vertex: vidx = ',vidx
                          write(stdo,*) ' icx=',icx,' icx0=',icx0
                          stop
                       endif
                    endif
                 endif
                 if(amr_ydim.ne.0) then
                    if(icy0.le.0) then
                       icy0 =icy
                    else
                       if(icy0.ne.icy) then
                          write(stdo,*) 'ERROR: Inconsistent vertex: vidx = ',vidx
                          write(stdo,*) ' icy=',icy,' icy0=',icy0
                          stop
                       endif
                    endif
                 endif
                 if(amr_zdim.ne.0) then
                    if(icz0.le.0) then
                       icz0 =icz
                    else
                       if(icz0.ne.icz) then
                          ! Attila Juhasz
                          ! WARNING!!!!!
                          ! Bypassing index consistency check for periodic boundary
                          ! conditions
                          !
                          ! There is a problem in spherical coordinates if
                          ! cyclic periodic boundary condition is used for phi
                          ! then while icz or icz0 is 1 the other one is nz+1
                          ! So I switch that off for periodic boundary
                          ! conditions
                          if (.not.amr_cyclic_xyz(3)) then
                              write(stdo,*) 'ERROR: Inconsistent vertex: vidx = ',vidx
                              write(stdo,*) ' icz=',icz,' icz0=',icz0
                              stop
                          else
                             if ((icz.ne.1).or.(icz0.ne.amr_grid_nz+1)) then
                                write(stdo,*) 'ERROR: Inconsistent vertex: vidx = ',vidx
                                write(stdo,*) ' icz=',icz,' icz0=',icz0
                                stop
                             endif
                          endif 
                       endif
                    endif
                 endif
              endif
           endif
        enddo
     enddo
  enddo
  !
  ! If we have not found even one corner cell, something must be wrong
  !
  if((icx0.lt.0).or.(icy0.lt.0).or.(icz0.lt.0)) then
     write(stdo,*) 'ERROR: Corner cells inconsistent...'
     stop
  endif
  !
  ! If requested, then return the x,y,z coordinates of this point
  !
  if(present(x)) x=amr_finegrid_xi(icx,1,amr_levelmax)
  if(present(y)) y=amr_finegrid_xi(icy,2,amr_levelmax)
  if(present(z)) z=amr_finegrid_xi(icz,3,amr_levelmax)
  !
end subroutine amr_check_vertex


!-------------------------------------------------------------------
!      INTEGER TO STRING, TRIMMED (ONLY FOR POSITIVE INTEGERS)
!-------------------------------------------------------------------
subroutine amr_integer_to_string(int,string)
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
  elseif(int.lt.1000000000) then
     write(string,19) int
19   format(I9)
  else
     write(string,*) int
  endif
end subroutine amr_integer_to_string

!--------------------------------------------------------------------------
!                         OPEN VTK OUTPUT
!
! VTK Stands for Visual Tool Kit, and is a format for visualizing 3-D
! data. The Paraview viewing tool uses it, for instance.
!--------------------------------------------------------------------------
subroutine amr_open_vtk_file(unit,filename,title,data, igrid_coord)
  implicit none
  character*160 :: filename,title
  character*256 :: longstring
  character*80 :: nvtx_string,ncell_string,string,sic(1:8),slen,slen2,form
  integer :: ivtx,icell,idum,index,unit,igrid_coord,nr_wedge_cells, ierr
  double precision :: x,y,z, cart_x, cart_y, cart_z
  logical :: data
  !
  ! Check
  !
  if(.not.amr_tree_present) then
     write(stdo,*) 'ERROR: Need to use the AMR-tree to be able to'
     write(stdo,*) '       create an AMR VTK output file.'
     write(stdo,*) '       Add command line option usetree to radmc3d call.'
     stop
  endif
  if(.not.allocated(amr_thebranches)) then
     call amr_compute_list_all()
  endif
  if(.not.allocated(amr_vertex_cells)) then
     write(stdo,*) 'Setting up the vertices (cell corners) for second order transfer...'
     call amr_set_up_vertices()
  endif
  !
  ! Convert integers to strings
  !
  call amr_integer_to_string(amr_nr_vertices,nvtx_string)
  call amr_integer_to_string(amr_nrleafs,ncell_string)
  !
  ! Allocate the cell type ID array
  !
  allocate(vtk_cell_type(amr_nrleafs), STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in AMR Module: Could not allocate vtk_cell_type()'
     stop 2711
  endif
  vtk_nr_cells(:) = 0
  !
  ! Open file and write header
  !
  open(unit=unit,file=filename)
  !
  ! Write the VTK header
  !
  write(unit,'(A26)') '# vtk DataFile Version 1.0'
  call amr_integer_to_string(len_trim(title),slen)
  form='(A'//trim(slen)//')'
  write(unit,form) trim(title)
  write(unit,'(A5)') 'ASCII'
  write(unit,*) ' '
  !
  ! Write the grid header
  !
  write(unit,'(A25)') 'DATASET UNSTRUCTURED_GRID'
  call amr_integer_to_string(len_trim(nvtx_string),slen)
  form='(A7,A'//trim(slen)//',1X,A6)'
  write(unit,form) 'POINTS ',trim(nvtx_string),'float'
  ! ------------------------------------------------------------------------------
  ! Attila Juhasz
  ! 
  ! Spherical coordinate system
  !
  ! If we have a spherical mesh hexahedron cell type (12) will be used instead of
  ! voxel (11). VTK does not support true spherical grid but we can use this
  ! trick to display data on a spherical mesh with VTK. For high enough resolution 
  ! in phi and theta it will not be really visible that the cell edges are actually 
  ! not curved.  
  !
  ! ------------------------------------------------------------------------------

  if((igrid_coord.ge.100).and.(igrid_coord.lt.200)) then
     do ivtx=1,amr_nr_vertices
        call amr_check_vertex(ivtx,x,y,z)
  !
  ! Transform the spherical coordinates to cartesian system
  !      
        cart_x = x * sin(y) * cos(z)
        cart_y = x * sin(y) * sin(z)
        cart_z = x * cos(y) 
        if (y.eq.pi) then
            cart_x = 0.0
            cart_y = 0.0
        endif
        write(unit,100) cart_x,cart_y,cart_z
        ! NOTE: If you have extreme refinement, you may want to change to E19.12
100 format(3(E13.6,1X))
!100 format(3(E19.12,1X))
     enddo
     
     !
     ! Now go through all cells and check which one is a hexahedron and which
     ! is a wedge
     !

     vtk_nr_cells(:) = 0

     do icell=1,amr_nrleafs
        index = amr_theleaf_index(icell)
        vtk_cell_type(icell) = '12' ! Hexahedron - 12
        call amr_check_vertex(amr_cell_corners(1,1,1,index),x,y,z)
        if ((y.eq.0.).or.(y.eq.pi))  vtk_cell_type(icell) = '13'! VTK_Wedge - 13 
        call amr_check_vertex(amr_cell_corners(2,1,1,index),x,y,z)
        if ((y.eq.0.).or.(y.eq.pi))  vtk_cell_type(icell) = '13'! VTK_Wedge - 13 
        call amr_check_vertex(amr_cell_corners(1,2,1,index),x,y,z)
        if ((y.eq.0.).or.(y.eq.pi))  vtk_cell_type(icell) = '13'! VTK_Wedge - 13 
        call amr_check_vertex(amr_cell_corners(2,2,1,index),x,y,z)
        if ((y.eq.0.).or.(y.eq.pi))  vtk_cell_type(icell) = '13'! VTK_Wedge - 13 
        call amr_check_vertex(amr_cell_corners(1,1,2,index),x,y,z)
        if ((y.eq.0.).or.(y.eq.pi))  vtk_cell_type(icell) = '13'! VTK_Wedge - 13 
        call amr_check_vertex(amr_cell_corners(2,1,2,index),x,y,z)
        if ((y.eq.0.).or.(y.eq.pi))  vtk_cell_type(icell) = '13'! VTK_Wedge - 13 
        call amr_check_vertex(amr_cell_corners(1,2,2,index),x,y,z)
        if ((y.eq.0.).or.(y.eq.pi))  vtk_cell_type(icell) = '13'! VTK_Wedge - 13 
        call amr_check_vertex(amr_cell_corners(2,2,2,index),x,y,z)
        if ((y.eq.0.).or.(y.eq.pi))  vtk_cell_type(icell) = '13'! VTK_Wedge - 13 

        if (vtk_cell_type(icell).eq.'13') then
            vtk_nr_cells(13) = vtk_nr_cells(13) + 1
        else
            vtk_nr_cells(12) = vtk_nr_cells(12) + 1
        endif
     enddo

     ! 
     ! For the moment exclude the wedge cells around the poles and write only
     ! the hexahedron cells
     ! 
     write(unit,*) ' '
     call amr_integer_to_string((vtk_nr_cells(12))*9+vtk_nr_cells(13)*7,string)
     call amr_integer_to_string(len_trim(ncell_string),slen)
     call amr_integer_to_string(len_trim(string),slen2)

     call amr_integer_to_string(vtk_nr_cells(12)+vtk_nr_cells(13),ncell_string)
     form='(A6,A'//trim(slen)//',1X,A'//trim(slen2)//')'
     write(unit,form) 'CELLS ',trim(ncell_string),trim(string)

     !
     ! Write out the vertex indices belonging to a given hexahedron/wedge
     !
     do icell=1,amr_nrleafs
        index = amr_theleaf_index(icell)
        !
        ! If we have a hexahedron cell then go ahead
        !
        if (vtk_cell_type(icell).eq.'12') then
            !
            ! The 3-4 and 7-8 corner index pairs are the switched around 
            !  in a hexahedron cell compared to a voxel
            ! 
            call amr_integer_to_string(amr_cell_corners(1,1,1,index)-1,sic(1))
            call amr_integer_to_string(amr_cell_corners(1,2,1,index)-1,sic(2))
            call amr_integer_to_string(amr_cell_corners(1,2,2,index)-1,sic(3))
            call amr_integer_to_string(amr_cell_corners(1,1,2,index)-1,sic(4))
            call amr_integer_to_string(amr_cell_corners(2,1,1,index)-1,sic(5))
            call amr_integer_to_string(amr_cell_corners(2,2,1,index)-1,sic(6))
            call amr_integer_to_string(amr_cell_corners(2,2,2,index)-1,sic(7))
            call amr_integer_to_string(amr_cell_corners(2,1,2,index)-1,sic(8))

            longstring = '8 '//trim(sic(1))//' '//trim(sic(2))//' '//trim(sic(3))//' '//trim(sic(4))//' '//&
                            trim(sic(5))//' '//trim(sic(6))//' '//trim(sic(7))//' '//trim(sic(8))
            call amr_integer_to_string(len_trim(longstring),slen)
            form='(A'//trim(slen)//')'
            write(unit,form) longstring
        !
        ! If we have a wedge cell
        !
        else
                
            call amr_check_vertex(amr_cell_corners(1,1,1,index),x,y,z)  
            ! 
            ! If we are on the northen hemisphere
            !
            if (y.eq.0.0) then
                call amr_integer_to_string(amr_cell_corners(1,1,1,index)-1,sic(1))
                call amr_integer_to_string(amr_cell_corners(1,2,2,index)-1,sic(2))
                call amr_integer_to_string(amr_cell_corners(1,2,1,index)-1,sic(3))
                call amr_integer_to_string(amr_cell_corners(2,1,1,index)-1,sic(4))
                call amr_integer_to_string(amr_cell_corners(2,2,2,index)-1,sic(5))
                call amr_integer_to_string(amr_cell_corners(2,2,1,index)-1,sic(6))
            !
            ! If we are on the southern hemisphere
            !
            else
                call amr_integer_to_string(amr_cell_corners(2,2,1,index)-1,sic(1))
                call amr_integer_to_string(amr_cell_corners(2,1,2,index)-1,sic(3))
                call amr_integer_to_string(amr_cell_corners(2,1,1,index)-1,sic(2))
                call amr_integer_to_string(amr_cell_corners(1,2,1,index)-1,sic(4))
                call amr_integer_to_string(amr_cell_corners(1,1,2,index)-1,sic(6))
                call amr_integer_to_string(amr_cell_corners(1,1,1,index)-1,sic(5))
            endif

            longstring = '6 '//trim(sic(1))//' '//trim(sic(2))//' '//trim(sic(3))//' '//trim(sic(4))//' '//&
                            trim(sic(5))//' '//trim(sic(6))
            call amr_integer_to_string(len_trim(longstring),slen)
            form='(A'//trim(slen)//')'
            write(unit,form) longstring
        endif
     enddo

     write(unit,*) ' '
     call amr_integer_to_string(len_trim(ncell_string),slen)
     form='(A11,A'//trim(slen)//')'
     write(unit,form) 'CELL_TYPES ',trim(ncell_string)
     do icell=1, amr_nrleafs
        write(unit,'(A2)') vtk_cell_type(icell)   ! 12 = Hexahedron, 13 - VTK_wedge
     enddo
  ! ------------------------------------------------------------------------------
  ! Cartesian coordinate system
  ! ------------------------------------------------------------------------------
  elseif (igrid_coord.lt.100) then
     !
     ! First write out the vertices 
     !
     do ivtx=1,amr_nr_vertices
        call amr_check_vertex(ivtx,x,y,z)
        write(unit,100) x,y,z
        ! NOTE: If you have extreme refinement, you may want to change to E19.12
!100 format(3(E13.6,1X))
     enddo
!100 format(3(E19.12,1X))
     !
     ! Now write out the vertex indices belonging to a given voxel
     !
     write(unit,*) ' '
     call amr_integer_to_string(amr_nrleafs*9,string)
     call amr_integer_to_string(len_trim(ncell_string),slen)
     call amr_integer_to_string(len_trim(string),slen2)
        
     form='(A6,A'//trim(slen)//',1X,A'//trim(slen2)//')'
     write(unit,form) 'CELLS ',trim(ncell_string),trim(string)
     do icell=1,amr_nrleafs
         index = amr_theleaf_index(icell)
         call amr_integer_to_string(amr_cell_corners(1,1,1,index)-1,sic(1))
         call amr_integer_to_string(amr_cell_corners(2,1,1,index)-1,sic(2))
         call amr_integer_to_string(amr_cell_corners(1,2,1,index)-1,sic(3))
         call amr_integer_to_string(amr_cell_corners(2,2,1,index)-1,sic(4))
         call amr_integer_to_string(amr_cell_corners(1,1,2,index)-1,sic(5))
         call amr_integer_to_string(amr_cell_corners(2,1,2,index)-1,sic(6))
         call amr_integer_to_string(amr_cell_corners(1,2,2,index)-1,sic(7))
         call amr_integer_to_string(amr_cell_corners(2,2,2,index)-1,sic(8))
         longstring = '8 '//trim(sic(1))//' '//trim(sic(2))//' '//trim(sic(3))//' '//trim(sic(4))//' '//&
                            trim(sic(5))//' '//trim(sic(6))//' '//trim(sic(7))//' '//trim(sic(8))
         call amr_integer_to_string(len_trim(longstring),slen)
         form='(A'//trim(slen)//')'
         write(unit,form) longstring
    enddo
    write(unit,*) ' '
    call amr_integer_to_string(len_trim(ncell_string),slen)
    form='(A11,A'//trim(slen)//')'
    write(unit,form) 'CELL_TYPES ',trim(ncell_string)
    do icell=1,amr_nrleafs
        write(unit,'(A2)') '11'   ! 11 = Voxel
    enddo
  endif
  if(data) then
     call amr_integer_to_string(amr_nrleafs,ncell_string)
     write(unit,*) ' '
     string = 'CELL_DATA '//trim(ncell_string)
     call amr_integer_to_string(len_trim(string),slen)
     form = '(A'//trim(slen)//')'
     write(unit,form) string
  endif
end subroutine amr_open_vtk_file


!--------------------------------------------------------------------------
!               WRITE A SCALAR FIELD TO THE VTK FILE
!--------------------------------------------------------------------------
subroutine amr_write_vtk_file_scalar(unit,ncl,scalar,log,scalarName)
  implicit none
  character*80 :: string,ncell_string,slen,form,scalarName
  integer :: icell,idum,index,unit,ncl
  double precision :: scalar(1:ncl),dummy
  logical :: log
  !
  ! Convert integers to strings
  !
  ! ----------------------------------------------------------------------
  ! Attila Juhasz 
  ! Add name to the data field header
  !  string = 'SCALARS scalars float'
  ! ----------------------------------------------------------------------
  string = 'SCALARS '//scalarName(1:len_trim(scalarName))//' float'

  call amr_integer_to_string(len_trim(string),slen)
  form = '(A'//trim(slen)//')'
  write(unit,form) string
  string = 'LOOKUP_TABLE default'
  call amr_integer_to_string(len_trim(string),slen)
  form = '(A'//trim(slen)//')'
  write(unit,form) string
  do icell=1,amr_nrleafs
     index = amr_theleaf_index(icell)
     dummy = scalar(index)
     if(log) then
        if(dummy.lt.1d-90) dummy=1d-90
        dummy = log10(dummy)
     endif
     if(abs(dummy).lt.1d-30) dummy=1d-30
     write(unit,200) dummy
200 format(E13.6)
  enddo
end subroutine amr_write_vtk_file_scalar


!--------------------------------------------------------------------------
!               WRITE A SCALAR FIELD TO THE VTK FILE
!--------------------------------------------------------------------------
subroutine amr_write_vtk_file_scalar_alt(unit,nn,ncl,ii,scalar,log,scalarName)
  implicit none
  character*80 :: string,ncell_string,slen,form,scalarName
  integer :: icell,idum,index,unit,ncl,nn,ii
  double precision :: scalar(1:nn,1:ncl),dummy
  logical :: log
  !
  ! Convert integers to strings
  !
  ! ----------------------------------------------------------------------
  ! Attila Juhasz 
  ! Add name to the data field header
  !string = 'SCALARS scalars float'
  ! ----------------------------------------------------------------------
  call amr_integer_to_string(ii,slen)
  string = 'SCALARS '//scalarName(1:len_trim(scalarName))//'_'//slen(1:len_trim(slen))//' float'
  call amr_integer_to_string(len_trim(string),slen)
  form = '(A'//trim(slen)//')'
  write(unit,form) string
  string = 'LOOKUP_TABLE default'
  call amr_integer_to_string(len_trim(string),slen)
  form = '(A'//trim(slen)//')'
  write(unit,form) string
  do icell=1,amr_nrleafs
     index = amr_theleaf_index(icell)
     dummy = scalar(ii,index)
     if(log) then
        if(dummy.lt.1d-90) dummy=1d-90
        dummy = log10(dummy)
     endif
     if(abs(dummy).lt.1d-30) dummy=1d-30
     write(unit,200) dummy
200 format(E13.6)
  enddo
end subroutine amr_write_vtk_file_scalar_alt


!--------------------------------------------------------------------------
!               WRITE A SCALAR FIELD TO THE VTK FILE
!--------------------------------------------------------------------------
subroutine amr_write_vtk_file_scalar_alt2(unit,nn,mm,ncl,ii,jj,scalar,log,scalarName)
  implicit none
  character*80 :: string,ncell_string,slen,form,scalarName,slen2
  integer :: icell,idum,index,unit,ncl,nn,mm,ii,jj
  double precision :: scalar(1:nn,1:mm,1:ncl),dummy
  logical :: log
  !
  ! Convert integers to strings
  !
  ! ----------------------------------------------------------------------
  ! Attila Juhasz 
  ! Add name to the data field header
  !  string = 'SCALARS scalars float'
  ! ----------------------------------------------------------------------
  call amr_integer_to_string(ii,slen)
  call amr_integer_to_string(jj,slen2)
  string = 'SCALARS '//scalarName(1:len_trim(scalarName))//'_'//slen(1:len_trim(slen))//'_'//slen2(1:len_trim(slen2))//' float'

  call amr_integer_to_string(len_trim(string),slen)
  form = '(A'//trim(slen)//')'
  write(unit,form) string
  string = 'LOOKUP_TABLE default'
  call amr_integer_to_string(len_trim(string),slen)
  form = '(A'//trim(slen)//')'
  write(unit,form) string
  do icell=1,amr_nrleafs
     index = amr_theleaf_index(icell)
     dummy = scalar(ii,jj,index)
     if(log) then
        if(dummy.lt.1d-90) dummy=1d-90
        dummy = log10(dummy)
     endif
     if(abs(dummy).lt.1d-30) dummy=1d-30
     write(unit,200) dummy
200 format(E13.6)
  enddo
end subroutine amr_write_vtk_file_scalar_alt2


!--------------------------------------------------------------------------
!               WRITE A VECTOR FIELD TO THE VTK FILE
!--------------------------------------------------------------------------
subroutine amr_write_vtk_file_vector(unit,ncl,vector,igrid_coord,vectorName)
  implicit none
  character*80 :: string,ncell_string,slen,form,vectorName
  integer :: icell,idum,index,unit,ncl,igrid_coord
  double precision :: vector(3,1:ncl),dummy(1:3)
  double precision :: x1, y1, z1, x2, y2, z2, r, phi, theta 


  !
  ! Convert integers to strings
  !
  ! ----------------------------------------------------------------------
  ! Attila Juhasz 
  ! Add name to the data field header
  !  string = 'VECTORS vectors float'
  ! ----------------------------------------------------------------------
  string = 'VECTORS '//vectorName(1:len_trim(vectorName))//' float'
  call amr_integer_to_string(len_trim(string),slen)
  form = '(A'//trim(slen)//')'
  write(unit,form) string
 
  ! ----------------------------------------------------------------------
  ! Attila Juhasz 
  ! 
  ! If we are in spherical coordinate system we have to transform vr,vtheta,vphi
  ! to vx,vy,vz
  ! ----------------------------------------------------------------------
  if((igrid_coord.ge.100).and.(igrid_coord.lt.200)) then   
      do icell=1, amr_nrleafs
        
        !
        ! Calculate cell center coordinates 
        ! This should be done in a more elegant way since the cell centers are
        ! already stored, this is just a quick and dirty solution
        !
        index = amr_theleaf_index(icell)
        call amr_check_vertex(amr_cell_corners(1,1,1,index),x1,y1,z1)
        call amr_check_vertex(amr_cell_corners(2,1,1,index),x2,y2,z2)
        r  = 0.5*(x1+x2)
        call amr_check_vertex(amr_cell_corners(1,1,1,index),x1,y1,z1)   
        call amr_check_vertex(amr_cell_corners(1,2,1,index),x2,y2,z2)   
        theta = 0.5*(y1+y2)
        call amr_check_vertex(amr_cell_corners(1,1,1,index),x1,y1,z1)   
        call amr_check_vertex(amr_cell_corners(1,1,2,index),x2,y2,z2)   
        if ((z2.eq.0.0).and.(amr_cyclic_xyz(3))) z2 = twopi
        phi = 0.5*(z1+z2)
    
        dummy(1) = vector(1,index)*sin(theta)*cos(phi) - vector(3,index)*sin(phi) + vector(2,index)*cos(theta)*cos(phi)    
        dummy(2) = vector(1,index)*sin(theta)*sin(phi) + vector(3,index)*cos(phi) + vector(2,index)*cos(theta)*sin(phi)
        dummy(3) = vector(1,index)*cos(theta) - vector(2,index)*sin(theta) 
       
        if(abs(dummy(1)).lt.1d-30) dummy(1)=1d-30
        if(abs(dummy(2)).lt.1d-30) dummy(2)=1d-30
        if(abs(dummy(3)).lt.1d-30) dummy(3)=1d-30
!        if (vtk_cell_type(icell).eq.'12')  write(unit,200) dummy
        write(unit,200) dummy(1),dummy(2),dummy(3)

      enddo
  else
      do icell=1,amr_nrleafs
         index = amr_theleaf_index(icell)
         dummy(:) = vector(:,index)
         if(abs(dummy(1)).lt.1d-30) dummy(1)=1d-30
         if(abs(dummy(2)).lt.1d-30) dummy(2)=1d-30
         if(abs(dummy(3)).lt.1d-30) dummy(3)=1d-30
!         if (vtk_cell_type(icell).eq.'12')  write(unit,200) dummy
         write(unit,200) dummy(1),dummy(2),dummy(3)
200 format(3(E13.6,1X))
      enddo
  endif
end subroutine amr_write_vtk_file_vector


!==========================================================================
!             ROUTINES FOR REGULAR GRID MANAGEMENT WITHOUT AMR TREE
!==========================================================================

!--------------------------------------------------------------------------
!                 FIND THE IX, IY AND IZ FROM INDEX
!--------------------------------------------------------------------------
subroutine amr_regular_get_ixyz(index,ix,iy,iz)
  implicit none
  integer :: index,ix,iy,iz,nxy,idx1
  nxy  = amr_grid_nx*amr_grid_ny
  iz   = (index-1)/nxy
  idx1 = index-iz*nxy
  iz   = iz+1
  iy   = (idx1-1)/amr_grid_nx
  ix   = idx1-iy*amr_grid_nx
  iy   = iy+1
!###################################
! Safety tests; can be removed later
  if((ix.lt.1).or.(iy.lt.1).or.(iz.lt.1)) then
     write(stdo,*) 'ERROR in computing ix,iy,iz in the regular grid:'
     write(stdo,*) index,ix+(iy-1)*amr_grid_nx+(iz-1)*nxy
     write(stdo,*) ix,iy,iz
     write(stdo,*) amr_grid_nx,amr_grid_ny,amr_grid_nz,nxy
     stop 293
  endif
  if(ix+(iy-1)*amr_grid_nx+(iz-1)*nxy.ne.index) stop 294
!  write(stdo,*) ix,iy,iz
!###################################
end subroutine amr_regular_get_ixyz

!--------------------------------------------------------------------------
!              FIND THE IX, IY AND IZ FROM VERTEX INDEX
!--------------------------------------------------------------------------
subroutine amr_regular_get_vertex_ixyz(index,ix,iy,iz)
  implicit none
  integer :: index,ix,iy,iz,nxy,idx1,nnx,nny
  nnx  = amr_grid_nx+amr_xdim
  nny  = amr_grid_ny+amr_ydim
  nxy  = nnx*nny
  iz   = (index-1)/nxy
  idx1 = index-iz*nxy
  iz   = iz+1
  iy   = (idx1-1)/nnx
  ix   = idx1-iy*nnx
  iy   = iy+1
!###################################
! Safety tests; can be removed later
  if(ix.lt.1) stop 391
  if(iy.lt.1) stop 392
  if(iz.lt.1) stop 393
  if(ix+(iy-1)*nnx+(iz-1)*nxy.ne.index) stop 394
!###################################
end subroutine amr_regular_get_vertex_ixyz


!==========================================================================
!                   ROUTINES FOR PIECEWISE LINEAR STUFF
!==========================================================================

!--------------------------------------------------------------------------
!                FIND AVERAGE SCALAR VALUE IN SOME BRANCH
!--------------------------------------------------------------------------
recursive subroutine amr_find_average_scalarvalue(q,n,b,value)
implicit none
integer :: n,ix,iy,iz
double precision :: q(n),value,dummy
type(amr_branch), pointer :: b
if(associated(b%child)) then
   value = 0.d0
   do iz=1,1+amr_zdim
      do iy=1,1+amr_ydim
         do ix=1,1+amr_xdim
            if(.not.associated(b%child(ix,iy,iz)%link)) then
               write(stdo,*) 'INTERNAL ERROR in AMR Module: Child of branch not associated...'
               stop 791
            endif
            call amr_find_average_scalarvalue(q,n,b%child(ix,iy,iz)%link,dummy)
            value = value + dummy
         enddo
      enddo
   enddo
   value = value / (2**amr_dim)
else
   value = q(b%leafindex)
endif
end subroutine amr_find_average_scalarvalue


!--------------------------------------------------------------------------
!               FIND SLOPE SIGMA OF SOME SCALAR QUANTITY
!--------------------------------------------------------------------------
subroutine amr_find_pl_scalarslope(q,n,b,slope)
implicit none
integer :: n,ilv,idir
double precision :: q(n),slope(1:3),qval(-1:1),dummy
type(amr_branch), pointer :: a,b
double precision :: sa,sb,s1,s2,x0,dxl,dxr
!
! Get the central value q_{i}
!
qval(0) = q(b%leafindex)
!
! Loop over dimensions
!
do idir=1,3
   if(amr_xyzdim(idir).eq.1) then
      ilv = b%level
      x0  = amr_finegrid_xc(b%ixyzf(idir),idir,ilv)
      !
      ! First the left hand
      !
      dxl = x0 - amr_finegrid_xc(b%ixyzf(idir)-1,idir,ilv)
      a   => b%neighbor(1,idir)%link
      if(associated(a)) then
         if(a%leaf) then
            if(a%level.eq.ilv) then
               !
               ! Neighbor cell = leaf of equal size
               !
               qval(-1) = q(a%leafindex)
            else
               !
               ! Neighbor cell = leaf of larger size
               !
               dxl = x0 - amr_finegrid_xc(a%ixyzf(idir),idir,a%level)
               qval(-1) = q(a%leafindex)
            endif
         else
            !
            ! Neighbor cell = branch with children
            !
            call amr_find_average_scalarvalue(q,n,a,dummy)
            qval(-1) = dummy
         endif
      else
         qval(-1) = qval(0)
      endif
      !
      ! Then the right hand
      !
      dxr = amr_finegrid_xc(b%ixyzf(idir)+1,idir,ilv) - x0
      a   => b%neighbor(2,idir)%link
      if(associated(a)) then
         if(a%leaf) then
            if(a%level.eq.ilv) then
               !
               ! Neighbor cell = leaf of equal size
               !
               qval(1) = q(a%leafindex)
            else
               !
               ! Neighbor cell = leaf of larger size
               !
               dxr = amr_finegrid_xc(a%ixyzf(idir),idir,a%level) - x0
               qval(1) = q(a%leafindex)
            endif
         else
            !
            ! Neighbor cell = branch with children
            !
            call amr_find_average_scalarvalue(q,n,a,dummy)
            qval(1) = dummy
         endif
      else
         qval(1) = qval(0)
      endif
      !
      ! Now decide based on slope limiter
      !
      if(amr_slope_limiter.eq.1) then
         !
         ! Use MinMod slope limiter
         !
         sa = ( qval(0) - qval(-1) ) / dxl
         sb = ( qval(1) - qval(0) ) / dxr
         if(sa*sb.lt.0.d0) then
            slope(idir) = 0.d0
         else
            if(abs(sa).lt.abs(sb)) then
               slope(idir) = sa
            else
               slope(idir) = sb
            endif
         endif
      elseif(amr_slope_limiter.eq.2) then
         !
         ! Use superbee slope limiter
         !
         sa = ( qval(0) - qval(-1) ) / dxl
         sb = ( qval(1) - qval(0) ) / dxr
         if(sa*sb.lt.0.d0) then
            slope(idir) = 0.d0
         else
            if(abs(sa).lt.abs(2*sb)) then
               s1 = sa
            else
               s1 = 2*sb
            endif
            if(abs(2*sa).lt.abs(sb)) then
               s2 = 2*sa
            else
               s2 = sb
            endif
            if(s1*s2.lt.0.d0) then
               slope(idir) = 0.d0
            else
               if(abs(s1).gt.abs(s2)) then
                  slope(idir) = s1
               else
                  slope(idir) = s2
               endif
            endif
         endif
      else
         write(stdo,*) 'ERROR in AMR module: slope limiter ',amr_slope_limiter,&
              ' not known. Aborting.'
         stop 
      endif
   else
      slope(idir) = 0.d0
   endif
enddo
end subroutine amr_find_pl_scalarslope


!--------------------------------------------------------------------------
!                  FIND AVERAGE VECTOR VALUE IN SOME BRANCH
!--------------------------------------------------------------------------
recursive subroutine amr_find_average_vectorvalue(q,k,n,b,value)
implicit none
integer :: k,n,ix,iy,iz
double precision :: q(k,n),value(k),dummy(k)
type(amr_branch), pointer :: b
if(associated(b%child)) then
   value(:) = 0.d0
   do iz=1,1+amr_zdim
      do iy=1,1+amr_ydim
         do ix=1,1+amr_xdim
            if(.not.associated(b%child(ix,iy,iz)%link)) then
               write(stdo,*) 'INTERNAL ERROR in AMR Module: Child of branch not associated...'
               stop 791
            endif
            call amr_find_average_vectorvalue(q,k,n,b%child(ix,iy,iz)%link,dummy)
            value(:) = value(:) + dummy(:)
         enddo
      enddo
   enddo
   value(:) = value(:) / (2**amr_dim)
else
   value(:) = q(:,b%leafindex)
endif
end subroutine amr_find_average_vectorvalue


!--------------------------------------------------------------------------
!              FIND SLOPE SIGMA OF SOME VECTORIAL QUANTITY
!--------------------------------------------------------------------------
subroutine amr_find_pl_vectorslope(q,k,n,b,slope)
implicit none
integer :: k,n,ilv,idir,l
double precision :: q(k,n),slope(k,1:3),qval(k,-1:1),dummy(k)
type(amr_branch), pointer :: a,b
double precision :: sa,sb,s1,s2,x0,dxl,dxr
!
! Get the central value q_{i}
!
qval(:,0) = q(:,b%leafindex)
!
! Loop over dimensions
!
do idir=1,3
   if(amr_xyzdim(idir).eq.1) then
      ilv = b%level
      x0  = amr_finegrid_xc(b%ixyzf(idir),idir,ilv)
      !
      ! First the left hand
      !
      dxl = x0 - amr_finegrid_xc(b%ixyzf(idir)-1,idir,ilv)
      a   => b%neighbor(1,idir)%link
      if(associated(a)) then
         if(a%leaf) then
            if(a%level.eq.ilv) then
               !
               ! Neighbor cell = leaf of equal size
               !
               qval(:,-1) = q(:,a%leafindex)
            else
               !
               ! Neighbor cell = leaf of larger size
               !
               dxl = x0 - amr_finegrid_xc(a%ixyzf(idir),idir,a%level)
               qval(:,-1) = q(:,a%leafindex)
            endif
         else
            !
            ! Neighbor cell = branch with children
            !
            call amr_find_average_vectorvalue(q,k,n,a,dummy)
            qval(:,-1) = dummy(:)
         endif
      else
         qval(:,-1) = qval(:,0)
      endif
      !
      ! Then the right hand
      !
      dxr = amr_finegrid_xc(b%ixyzf(idir)+1,idir,ilv) - x0
      a   => b%neighbor(2,idir)%link
      if(associated(a)) then
         if(a%leaf) then
            if(a%level.eq.ilv) then
               !
               ! Neighbor cell = leaf of equal size
               !
               qval(:,1) = q(:,a%leafindex)
            else
               !
               ! Neighbor cell = leaf of larger size
               !
               dxr = amr_finegrid_xc(a%ixyzf(idir),idir,a%level) - x0
               qval(:,1) = q(:,a%leafindex)
            endif
         else
            !
            ! Neighbor cell = branch with children
            !
            call amr_find_average_vectorvalue(q,k,n,a,dummy)
            qval(:,1) = dummy(:)
         endif
      else
         qval(:,1) = qval(:,0)
      endif
      !
      ! Now decide based on slope limiter
      !
      if(amr_slope_limiter.eq.1) then
         !
         ! Use MinMod slope limiter
         !
         do l=1,k
            sa = ( qval(l,0) - qval(l,-1) ) / dxl
            sb = ( qval(l,1) - qval(l,0) ) / dxr
            if(sa*sb.lt.0.d0) then
               slope(l,idir) = 0.d0
            else
               if(abs(sa).lt.abs(sb)) then
                  slope(l,idir) = sa
               else
                  slope(l,idir) = sb
               endif
            endif
         enddo
      elseif(amr_slope_limiter.eq.2) then
         !
         ! Use superbee slope limiter
         !
         do l=1,k
            sa = ( qval(l,0) - qval(l,-1) ) / dxl
            sb = ( qval(l,1) - qval(l,0) ) / dxr
            if(sa*sb.lt.0.d0) then
               slope(l,idir) = 0.d0
            else
               if(abs(sa).lt.abs(2*sb)) then
                  s1 = sa
               else
                  s1 = 2*sb
               endif
               if(abs(2*sa).lt.abs(sb)) then
                  s2 = 2*sa
               else
                  s2 = sb
               endif
               if(s1*s2.lt.0.d0) then
                  slope(l,idir) = 0.d0
               else
                  if(abs(s1).gt.abs(s2)) then
                     slope(l,idir) = s1
                  else
                     slope(l,idir) = s2
                  endif
               endif
            endif
         enddo
      else
         write(stdo,*) 'ERROR in AMR module: slope limiter ',amr_slope_limiter,&
              ' not known. Aborting.'
         stop 
      endif
   else
      slope(:,idir) = 0.d0
   endif
enddo
end subroutine amr_find_pl_vectorslope
!
end module amr_module
