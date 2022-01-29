.. _chap-gridding:

More information about the gridding
***********************************

We already discussed the various types of grids in Section
:ref:`sec-spatial-grid`, and the grid input file structure is described in
Section :ref:`sec-grid-input`. In this chapter let us take a closer look
at the gridding possibilities and things to take special care of.



.. _sec-regular-grid:

Regular grids
=============

A regular grid is called 'grid style 0' in RADMC-3D. It can be used in
Cartesian coordinates as well as in spherical coordinates (Section
:ref:`sec-coord-systems`).

A regular grid, in our definition, is a multi-dimensional grid which is
separable in :math:`x`, :math:`y` and :math:`z` (or in spherical coordinates in
:math:`r`, :math:`\theta` and :math:`\phi`). You specify a 1-D monotonically
increasing array of values :math:`x_1, x_2,\cdots,x_{\mathrm{nx+1}}` which
represent the cell walls in :math:`x-direction`.  You do the same for the other
directions: :math:`y_1, y_2,\cdots,y_{\mathrm{ny+1}}` and :math:`z_1,
z_2,\cdots,z_{\mathrm{nz+1}}`.  The value of, say, :math:`x_2` is the same for
every position in :math:`y` and :math:`z`: this is what we mean with
'separable'.

In Cartesian coordinates RADMC-3D enforces perfectly cubic grid cells (i.e.
linear grids). But that is only to make the image sub-pixeling easier (see
Section :ref:`sec-recursive-subpixeling`). For spherical grids this is not
enforced, and in fact it is strongly encouraged to use non-linear grids in
spherical coordinates. Please read Section :ref:`sec-separable-refinement`
if you use spherical coordinates!

In a regular grid you specify the grids in each direction separately.  For
instance, the x-grid is given by specifying the cell walls in x-direction. If we
have, say, 10 cells in x-direction, we must specify 11 cell wall positions. For
instance: :math:`x_i=\{-5,-4,-3,-2,-1,0,1,2,3,4,5\}`.  For the
:math:`y`-direction and :math:`z`-direction
likewise. Fig. :numref:`fig-regular-grid` shows an example of a 2-D regular grid
of 4x3 cells.

.. _fig-regular-grid:

.. figure:: Figures/base_amr_bare.*
   :width: 50%

   Example of a regular 2-D grid with ``nx=4`` and ``ny=3``.

In Cartesian coordinates we typically define our model in full 3-D.
However, if your problem has translational symmetries, you might also want
to consider the 1-D plane-parallel mode (see Section
:ref:`sec-1d-plane-parallel`). 

In full 3-D Cartesian coordinates the cell sizes *must* be perfectly
cubical, i.e. the spacing in each direction must be the same. If you need a
finer grid in some location, you can use the AMR capabilities discussed
below.

In spherical coordinates you can choose between 1-D spherically symmetric
models, 2-D axisymmetric models or fully 3-D models. In spherical coordinates
you do *not* have restrictions to the cell geometry or grid spacing. 
You can choose any set of numbers :math:`r_1,\cdots,r_{\mathrm{nr}}` as radial
grid, as long as this set of numbers is larger than 0 and monotonically
increasing. The same is true for the :math:`\theta`-grid and the :math:`\phi`-grid.

The precise way how to set up a regular grid using the ``amr_grid.inp`` file is
described in Section :ref:`sec-amr-grid-regular`.  The input of any spatial
variables (such as e.g. the dust density) uses the sequence of grid cells in
the same order as the cells are specified in that ``amr_grid.inp`` file.

For input and output data to file, for stuff on a regular grid, the order of
nested loops over coordinates would be:
::

  do iz=1,amr_grid_nz
     do iy=1,amr_grid_ny
        do ix=1,amr_grid_nx
           << read or write your data >>
        enddo
     enddo
  enddo


For spherical coordinates we have the following association: :math:`x\rightarrow r`,
:math:`y\rightarrow \theta`, :math:`z\rightarrow \phi`.





.. _sec-separable-refinement:

Separable grid refinement in spherical coordinates (important!)
===============================================================

Spherical coordinates are a very powerful way of dealing with
centrally-concentrated problems. For instance, collapsing protostellar cores,
protoplanetary disks, disk galaxies, dust tori around active galactic nuclei,
accretion disks around compact objects, etc. In other words: problems in which a
single central body dominates the problem, and material at all distances from
the central body matters. For example a disk around a young star goes all the
way from 0.01 AU out to 1000 AU, covering 5 orders of magnitude in
radius. Spherical coordinates are the easiest way of dealing with such a huge
radial dynamic range: you simply make a radial grid, where the grid spacing
:math:`r_{i+1}-r_i` scales roughly with :math:`r_i`.

This is called a *logarithmic radial grid*. This is a grid whith a spacing in
which :math:`(r_{i+1}-r_i)/r_i` is constant with :math:`r`. In this way you
assure that you have always the right spatial resolution in :math:`r` at each
radius. In spherical coordinates it is highly recomended to use such a log
spacing. But you can also refine the :math:`r` grid even more (in addition to
the log-spacing). This is also strongly recommended near the inner edge of a
circumstellar shell, for instance.  Or at the inner dust rim of a disk. There
you must refine the :math:`r` grid (by simply making the spacing smaller as you
approach the inner edge from the outside) to assure that the first few cells are
optically thin and that there is a gradual transition from optically thin to
optically thick as you go outward. This is particularly important for, for
instance, the inner rim of a dusty disk.

In spherical coordinates you can vary the spacing in :math:`r`, :math:`\theta`
and :math:`\phi` completely freely. That means: you could have for instance
:math:`r` to be spaced as :math:`1.00, 1.01, 1.03, 1.05, 1.1, 1.2, 1.35,
\cdots`. There is no restriction, as long as the coordinate points are
monotonically increasing. In Figs :numref:`fig-spher-sep-ref-txt` and
:numref:`fig-spher-sep-ref-rays` this is illustrated.

*Note that in addition to separable refinement, also AMR refinement
is possible in spherical coordinates. See Section* :ref:`sec-oct-tree-amr`.

.. _fig-spher-sep-ref-txt:

.. figure:: Figures/spher_grid_ref_txt.*
   :width: 30%

   Example of a spherical 2-D grid in which the radial and :math:`\theta` grids
   are refined in a 'separable' way. In radial direction the inner cells are
   refined ('B' in the right figure) and in :math:`\theta` direction the cells
   near the equatorial plane are refined ('C' in the right figure). This kind of
   grid refinement does not require oct-tree AMR: the grid remains
   separable. For models in which the inner grid edge is also the inner model
   edge (e.g. a simple model of a protoplanetary disk with a sharp inner
   cut-off) this kind of separable grid refinement in :math:`R`-direction may be
   essential to avoid problems with optically thick inner cells (see e.g. Fig.
   :numref:`fig-innerrim-lowres` for an example of what could go wrong if you do
   not do this). Separable grid refinement in :math:`\Theta`-direction is
   typically important for protoplanetary disk models, where the midplane and
   surface layers of the disk need to have sufficient resolution, but any
   possible surrounding spherical nebula may not.

.. _fig-spher-sep-ref-rays:

.. figure:: Figures/spher_grid_ref_rays.*
   :width: 65%
            
   When making an image, RADMC-3D will automatically make 'sub-pixels' to ensure
   that all structure of the model as projected on the sky of the observer are
   spatially resolved.  Extreme grid refinement leads thus to extreme
   sub-pixeling. See Section :ref:`sec-rec-subpixel-spher-coord` for details,
   and ways to prevent excessive sub-pixeling when this is not necessary.

For models of accretion disks it can, for instance, be useful to make sure that
there are more grid points of :math:`\theta` near the equatorial plane
:math:`\theta=\pi/2`. So the grid spacing between :math:`\theta=0.0` and
:math:`\theta=1.0` may be very coarse while between :math:`\theta=1.0` and
:math:`\theta=\pi/2` you may put a finer grid. All of this 'grid refinement' can
be done without the 'AMR' refinement technique: this is the 'separable' grid
refinement, because you can do this separately for :math:`r`, for :math:`\theta`
and for :math:`\phi`.

Sometimes, however, separable refinement may not help you to refine the grid
where necessary. For instance: if you model a disk with a planet in the
disk, then you may need to refine the grid around the planet. You could
refine the grid in principle in a separable way, but you would then have a
large redundancy in cells that are refined by far away from the planet.  Or
if you have a disk with an inner rim that is not exactly at
:math:`r=r_{\mathrm{rim}}`, but is a rounded-off rim. In these cases you need
refinement exactly located at the region of interest. For that you need the
'AMR' refinement (Sections :ref:`sec-oct-tree-amr` and :ref:`sec-layered-amr`).

*Important note:* When using strong refinement in one of the coordinates
:math:`r`, :math:`\theta` or :math:`\phi`, image-rendering and
spectrum-rendering can become very slow, because of the excessive sub-pixeling
this causes. There are ways to limit the sub-pixeling for those cases. See the
Section on sub-pixeling in spherical coordinate: Section
:ref:`sec-rec-subpixel-spher-coord`.



.. _sec-oct-tree-amr:

Oct-tree Adaptive Mesh Refinement
=================================

An oct-tree refinened grid is called 'grid style 1' in RADMC-3D. It can be
used in Cartesian coordinates as well as in spherical coordinates (Section
:ref:`sec-coord-systems`). 

You start from a normal regular base grid (see Section
:ref:`sec-regular-grid`), possibly even with 'separable refinement' (see
Section :ref:`sec-separable-refinement`). You can then split some of the cells
into 2x2x2 subcells (or more precisely: in 1-D 2 subcells, in 2-D 2x2
subcells and in 3-D 2x2x2 subcells). If necessary, each of these 2x2x2
subcells can also be split into further subcells. This can be repeated as
many times as you wish until the desired grid refinement level is reached.
Each refinement step refines the grid by a factor of 2 in linear dimension,
which means in 3-D a factor of 8 in volume. In this way you get, for each
refined cell of the base grid, a tree of refinement. The base grid can have
any size, as long as the number of cells in each direction is an even
number. For instance, you can have a 6x4 base grid in 2-D, and refine cell
(1,2) by one level, so that this cell splits into 2x2 subcells.

Note that it is important to set which dimensions are 'active' and which
are 'non-active'. For instance, if you have a 1-D model with 100 cells and
you tell RADMC-3D (see Section :ref:`sec-amr-grid-oct-tree`) to make a base
grid of 100x1x1 cells, but you still keep all three dimensions 'active'
(see Section :ref:`sec-amr-grid-oct-tree`), then a refinement of cell 1
(which is actually cell (1,1,1)) will split that cell into 2x2x2 subcells,
i.e. it will also refine in y and z direction. Only if you explicitly
switch the y and z dimensions off the AMR will split it into just 2
subcells.

Oct-tree mesh refinement is very powerful, because it allows you to refine
the grid exactly there where you need it. And because we start from a
regular base grid like the grid specified in Section :ref:`sec-regular-grid`,
we can start designing our model on a regular base grid, and then refine
where needed. See Fig. :numref:`fig-oct-tree-amr`.

.. _fig-oct-tree-amr:

.. figure:: Figures/oct_tree_amr_bare.*
   :width: 50%

   Example of a 2-D grid with oct-tree refinement. The base grid has ``nx=4``
   and ``ny=3``. Three levels of refinement are added to this base grid.

The AMR stand for 'Adaptive Mesh Refinement', which may suggest that 
RADMC-3D will refine internally. At the moment this is not yet the case.
The 'adaptive' aspect is left to the user: he/she will have to 'adapt'
the grid such that it is sufficiently refinened where it is needed. In the
future we may allow on-the-fly adaption of the grid, but that is not yet
possible now. 

One problem with oct-tree AMR is that it is difficult to handle such grids in
external plotting programs, or even in programs that set up the grid.  While it
is highly flexible, it is not very user-friendly. Typically you may use this
oct-tree refinement either because you import data from a hydrodynamics code
that works with oct-tree refinement (e.g. FLASH, RAMSES), or when you
internally refine the grid using the ``userdef_module.f90`` (see Chapter
:ref:`chap-internal-setup`). In the former case you are anyway forced to manage
the complexities of AMR, while in the latter case you can make use of the AMR
modules of RADMC-3D internally to handle them. But if you do not need to full
flexibility of oct-tree refinement and want to use a simpler kind of refinement,
then you can use RADMC-3D's alternative refinement mode: the layer-style AMR
described in Section :ref:`sec-layered-amr` below.

The precise way how to set up such an oct-tree grid using the ``amr_grid.inp``
file is described in Section :ref:`sec-amr-grid-oct-tree`.  The input of any
spatial variables (such as e.g. the dust density) uses the sequence of grid
cells in the same order as the cells are specified in that ``amr_grid.inp``
file.



.. _sec-layered-amr:

Layered Adaptive Mesh Refinement
================================

A layer-style refinened grid is called 'grid style 10' in RADMC-3D. It can
be used in Cartesian coordinates as well as in spherical coordinates
(Section :ref:`sec-coord-systems`).

This is an alternative to the full-fledged oct-tree refinement of Section
:ref:`sec-oct-tree-amr`. The main advantage of the layer-style refinement is
that it is far easier to handle by the human brain, and thus easier for 
model setup  and the analysis of the results. 

The idea here is that you start again with a regular grid (like that of
Section :ref:`sec-regular-grid`), but you can now specify a rectangular
region which you want to refine by a factor of 2. The way you do this is by
choosing the starting indices of the rectangle and specifying the size of
the rectangle by setting the number of cells in each direction from that
starting point onward. For instance, setting the starting point at (2,3,1)
and the size at (1,1,1) will simply refine just cell (2,3,1) of the base
grid into a set of 2x2x2 sub-cells. But setting the starting point at
(2,3,1) and the size at (2,2,2) will split cells (2,3,1), (3,3,1), (2,4,1),
(3,4,1), (2,3,2), (3,3,2), (2,4,2) and (3,4,2) each into 2x2x2 subcells.
This in fact is handled as a 4x4x4 regular sub-grid patch. And setting the
starting point at (2,3,1) and the size at (4,6,8) will make an entire
regular sub-grid patch of 8x12x16 cells. Such a sub-grid patch is
called a *layer*.

The nice thing of these layers is that each layer (i.e. subgrid patch) is
handled as a regular sub-grid. The base grid is layer number 0, and the first
layer is layer number 1, etc. Each layer (including the base grid) can contain
multiple sub-layers. The only restriction is that each layer fits entirely
inside its parent layer, and layers with the same parent layer should not
overlap. Each layer can thus have one or more sub-layers, each of which can
again be divided into sub-layers. This builds a tree structure, with the base
layer as the trunk of the tree (this is contrary to the oct-tree structure,
where each base grid *cell* forms the trunk of its own tree). In
Fig. :numref:`fig-twolayer-amr` an example is shown of two layers with the same
parent (= layer 0 = base grid), while in Fig. :numref:`fig-nestedlayer-amr` an
example is shown of two nested layers.

.. _fig-twolayer-amr:

.. figure:: Figures/twolayer_bare.*
   :width: 50%

   Example of a 2-D base grid with ``nx``\ =4 and ``ny``\ =3, with two
   AMR-layers added to it. This example has just one level of refinement, as
   the two layers (brown and green) are on the same level (they have the same
   parent layer = layer 0).

.. _fig-nestedlayer-amr:

.. figure:: Figures/nestedlayer_bare.*
   :width: 50%

   Example of a 2-D base grid with ``nx``\ =4 and ``ny``\ =3, with two nested
   AMR-layers added to it.  This example has two levels of refinement, as layer
   1 (brown) is the parent of layer 2 (green).

If you now want to specify data on this grid, then you simply specify it on
each layer separately, as if each layer is a separate entity. Each layer is
treated as a regular grid, irrespective of whether it contains sub-layers
or not. So if we have a base grid of 4x4x4 grid cells containing two layers:
one starting at (1,1,1) and having (2,2,2) size and another starting at
(3,3,3) and having (1,1,2) size, then we first specify the data on the 
:math:`4^3=64` base grid, then on the :math:`(2\times 2)^3=64` grid cells of the first
layer and then on the 2x2x4=16 cells of the second layer. Each of these
three layers are regular grids, and the data is inputted/outputted in
the same way as if these are normal regular grids (see Section
:ref:`sec-regular-grid`). But instead of just one such regular grid, now
the data file (e.g. ``dust_density.inp``\ ) will contain three
successive lists of numbers, the first for the base grid, the second for
the first layer and the last for the second layer. You may realize at this
point that this will introduce a redundancy. See Subsection
:ref:`sec-layer-amr-redundancy` for a discussion of this redundancy.

The precise way how to set up such an oct-tree grid using the ``amr_grid.inp``
file is described in Section :ref:`sec-amr-grid-layered`.  The input of any
spatial variables (such as e.g. the dust density) uses the sequence of grid
cells in the same order as the cells are specified in that ``amr_grid.inp``
file.


.. _sec-layer-amr-redundancy:

On the 'successively regular' kind of data storage, and its slight redundancy
-------------------------------------------------------------------------------

With the layered grid refinement style there will be *redundant* data in the
data files (such as e.g. the ``dust_density.inp`` file. Each layer is a regular
(sub-)grid and the data will be specified in each of these grid cells of that
regular (sub-)grid.  If then some of these cells are overwritten by a
higher-level layer, these data are then redundant. We could of course have
insistent that only the data in those cells that are not refined by a layer
should be written to (or read from) the data files. But this would require quite
some clever programming on the part of the user to a-priori find out where the
layers are and therefore which cells should be skipped. We have decided that it
is far easier to just insist that each layer (including the base grid, which is
layer number 0) is simply written to the data file as a regular block of
data. The fact that some of this data will be not used (because they reside in
cells that are refined) means that we write more data to file than really exists
in the model. This makes the files larger than strictly necessary, but it makes
the data structure by far easier. Example: suppose you have a base grid of 8x8x8
cells and you replace the inner 4x4x4 cells with a layer of 8x8x8 cells (each
cell being half the size of the original cells).  Then you will have for
instance a ``dust_density.inp`` file containing 1024 values of the density:
:math:`8^3`\ =512 values for the base grid and again :math:`8^3`\ =512 values for
the refinement layer. Of the first :math:`8^3`\ =512 values :math:`4^3`\ =64 values
are ignored (they could have any value as they will not be used). The file is
thus 64 values larger than strictly necessary, which is a redundancy of
64/1024=0.0625. If you would have used the oct-tree refinement style for making
exactly the same grid, you would have only 1024-64=960 values in your file,
making the file 6.25\% smaller. But since 6.25\% is just a very small
difference, we decided that this is not a major problem and the simplicity of
our 'successively regular' kind of data format is more of an advantage than the
6.25\% redundance is a disadvantage.




.. _sec-unstruct-grids:

Unstructured grids (Delaunay, Voronoi, or more general)
=======================================================

RADMC-3D can handle unstructures grids of a variety of types. This works (so
far) only for 3-D cartesian coordinates. The two most well-known are Delaunay
grids and Voronoi grids. But these are just special cases of a more general
unstructured grid capability. RADMC-3D only needs to know, for each grid cell:

#. The volume of the cell
#. The cell walls of the cell

and, conversely, for each cell wall which two cells that the wall separates.  A
cell wall is only defined by two vectors that define the 2D plane in 3D space: a
support vector :math:`{\bf s}` and a normal vector :math:`{\bf n}`. By giving to
RADMC-3D a file ``unstr_grid.inp`` containing all this information, the user determines
the shape of each cell. Each of these cells has an integer index, starting with
1, and increasing by steps of 1 until the number of cells. The density and other
variables in the usual input files such as ``dust_density.inp`` are then associated
to these cells in that order (first value belongs to cell 1, second to cell 2 etc).

This information is enough for RADMC-3D to compute where a ray passes from one
cell to its neighbor, and how the passage through the cell affects the cell's
temperature and scattering source function.

Some of this grid information is redundant, and depends on which information
is contained in ``unstr_grid.inp``. Typically the user provides only the information
about the grid he/she wants to provide, and RADMC-3D will try to complete the rest
if it can (otherwise it will give an error).

For example: For a Voronoi grid you only need to specify, for each cell wall,
which two cells the wall separates. If RADMC-3D lacks further information, it
will by default assume that the cell walls lie exactly in between the two cell
center points and perpendicular to the line connecting them, thus naturally
giving you a Voronoi grid. The support vectors and normal vectors of the walls are then
computed internally in RADMC-3D. Unfortunately, RADMC-3D lacks the capability to
calculate the cell volumes for Voronoi cells, so you will have to provide the
cell volume.

As another example: For a Delaunay grid you only need to specify the vertices
(corner points of the tetraheders), which vertices span a wall, and which two
cells the wall separates (or one cell and the vacuum). The support vectors and
normal vectors of the walls are then computed internally in RADMC-3D. 
For these simple shaped cells (just simplices) RADMC-3D can calculate the volume
itself.

Computing this information for a Delaunay or Voronoi grid can be done using
external software packages such as ``qhull``. In fact, the ``qhull`` package
is built-in into the ``scipy.spatial`` library of Python, and is thefore
easily usable from within Python. To create the necessary ``unstr_grid.inp``
files using these libraries you can use the following tools provided with
the RADMC-3D package:

#. ``python/radmc3d_tools/radmc3d_delaunay_grid.py`` : A tool to generate a 
   Delaunay grid from a given set of vertex points. Note that the cells are
   then created by this tool, and they may not be the same (and almost certainly
   not in the same order) as the cells you may have assumed beforehand. So
   if you obtain a model from someone else that is on a Delaunay grid, and
   the physical values are in the grid cells, you may need to figure out how
   to associate the cells from that other model with the cells generated by
   this package. However, usually the physical values would, in the case of
   a Delaunay grid, be specified at the cell vertices. Then you would need
   to interpolate into the (new) cells, because RADMC-3D needs the physical
   values inside the cells, not at the vertices.
   
#. ``python/radmc3d_tools/radmc3d_voronoi_grid.py`` : A tool to generate
   a Voronoi grid from a given set of cell center points. The physical
   variables are also specified on these points. The vertex points are
   not required as input. Only needed is information about each cell wall:
   which two cells are divided by the wall, and information about each cell:
   it's volume. In voronoi grids there will be 'open cells' which go out
   to infinity. The physical variables in these cells are ignored: These
   cells are assumed to be empty. 

.. _fig-delau:

.. figure:: Figures/delaunay_grid_2d.* 
   :width: 75%

   Example of unstructured grid: a Delaunay grid. Figure is only symbolic, since
   the actual grid is 3D. Blue points are cell centers, orange points are
   vertices.

.. _fig-voronoi:

.. figure:: Figures/voronoi_grid_2d.*
   :width: 75%

   Example of unstructured grid: a Voronoi grid. Figure is only symbolic, since
   the actual grid is 3D. Blue points are cell centers, orange points are
   vertices.

Examples of these two types of grids are shown in
:numref:`fig-delau` and :numref:`fig-voronoi`.

However, you can also generate your own grid designs by specifying the
support vectors and normal vectors and two cells belonging to each wall.
However, this is not an entirely trivial task to do correctly, because
the result has to obey the following conditions:

#. Each cell must be convex in shape.
#. Beware not to forget any necessary wall. A missing wall will likely
   cause very unexpected behavior.
#. There should be no "empty holes" between the cells: inside the grid
   cells should fill space perfectly. This is not a trivial condition
   to fulfill: ill defined grid walls can lead to photons finding themselves
   outside a grid cell, but within the grid. This leads to RADMC-3D
   to crash.
#. If you specify the cell walls with vertices instead of support- and
   normal vectors, the vertices belonging to a cell must lie perfectly
   in same plane.
#. If you have walls that face the vacuum, they are part of what is called
   the hull of the grid (Note: A Delaunay grid has a hull, but a Voronoi
   grid does not: it has instead open cells). It is best (though not
   strictly necessary) to assure that the hull is convex, since that
   makes the code faster.

The way to make your own (non-Delaunay, non-Voronoi) grid is to generate
a file ``unstr_grid.inp`` in the appropriate way. The technical details
how this file is formatted are given in Section :ref:`sec-unstr-grid`.


.. _sec-1d-plane-parallel:

1-D Plane-parallel models
=========================

Sometimes it can be useful to make simple 1-D plane parallel models, for
instance if you want to make a simple 1-D model of a stellar atmosphere.
RADMC-3D is, however, by nature a 3-D code. But as of version 0.31 it features a
genuine 1-D plane-parallel mode as well. This coordinate type has the
number 10. In this mode the :math:`x`- and :math:`y`-coordinates are the
in-plane coordinates, while the :math:`z`-coordinate is the 1-D coordinate.  We
thus have a 1-D grid in the :math:`z`-coordinate, but no grid in :math:`x`- or
:math:`y`-directions.

You can make a 1-D plane-parallel model by setting some settings in the
``amr_grid.inp`` file. Please consult Section :ref:`sec-grid-input` for the
format of this file. The changes/settings you have to do are (see example
below): (1) set the coordinate type number to 10, (2) set the :math:`x` and
:math:`y` dimensions to non-active and (3) setting the cell interfaces in
:math:`x` to -1d90, +1d90, and likewise for :math:`y`. Here is then how it
looks: ::

  1                                     <=== Format number = 1
  0                                     <=== Grid style (0=regular grid)
  10                                    <=== Coordinate type (10=plane-parallel)
  0                                     <=== (obsolete)
  0  0  1                               <=== x and y are non-active, z is active
  1  1  100                             <=== x and y are 1 cell, in z we have 100 cells
  -1e90 1e90                            <=== cell walls in x are at "infinity"
  -1e90 1e90                            <=== cell walls in y are at "infinity"
  zi[1]        zi[2]        zi[3]       ........  zi[nz+1]

The other input files are for the rest as usual, as in the 3-D case.

You can now make your 1-D model as usual. For 1-D plane-parallel problems it is
often useful to put a thermal boundary at the bottom of the model.  For
instance, if the model is a stellar atmosphere, you may want to cap the grid
from below with some given temperature. See Section
:ref:`sec-thermal-boundaries` for details on how to set up thermal boundaries.

In the 1-D plane-parallel mode some things work a bit different than in the
"normal" 3-D mode:

* Images are by default 1x1 pixels, because in a plane-parallel case it
  is useless to have multiple pixels.
* Spectra cannot be made, because "spectrum" is (in RADMC-3D 'language')
  the flux as a function of frequency as seen at a very large distance of
  the object, so that the object is in the "far field". Since the concept of
  "far-field" is no longer meaningful in a plane-parallel case, it is better
  to make frequency-dependent 1x1 pixel images. This gives you the
  frequence-dependent intensity, which is all you should need.
* Stars are not allowed, as they have truly 3-D positions, which is
  inconsistent with the plane-parallel assumption.

But for the rest, most stuff works similarly to the 3-D version. For instance,
you can compute dust temperatures with 
::

  radmc3d mctherm

as usual.



Making a spectrum of the 1-D plane-parallel atmosphere
------------------------------------------------------

As mentioned above, the 'normal' 3-D way of making a spectrum of the 1-D
plane-parallel atmosphere is not possible, because formally the atmosphere is
infinitely extended. Instead you can obtain a spectrum in the form of an
intensity (:math:`\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-2}\,\mathrm{Hz}^{-1}\,\mathrm{ster}^{-1}`)
as a function of wavelength. To do this you ask RADMC-3D to make a
multi-wavelength image of the atmosphere under a certain inclination
(inclination 0 meaning face-on), e.g.: ::

  radmc3d image allwl incl 70

This make an SED at :math:`\lambda=10\,\mu`\ m for the observer seeing the
atmosphere at an inclination of 70 degrees. This produces a file image.out,
described in Section :ref:`sec-image-out`. The image is, in fact, a 1x1 pixel
multi-wavelength image. The ``allwl`` (which stands for 'all wavelengths') means
that the spectral points are the same as those in the ``wavelength_micron.inp``
file (see Section :ref:`sec-wavelengths`).  You can also specify the wavelengths
in a different way, e.g.: ::

  radmc3d image lambdarange 5 20 nlam 10

In fact, see Section :ref:`sec-multi-wavelength-images` and
Section :ref:`sec-set-camera-frequencies` for details.



In 1-D plane-parallel: no star, but incident parallel flux beams
----------------------------------------------------------------

In 1-D plane-parallel geometry it is impossible to include meaningful stars as
sources of photons. This is not a technical issue, but a mathematical truth: a
point in 1-D is in reality a plane in 3-D. As a replacement RADMC-3D offers
(only in 1-D plane-parallel geometry) the possibility of illuminating the 1-D
atmosphere from above with a flux, incident onto the atmosphere in a prescribed
angle. This allows you to model, e.g., the Earth's atmosphere being illuminated
by the sun at a given time of the day.  This is done by providing an ascii file
called ``illum.inp`` which has the following form (similar, but not identical,
to the ``stars.inp`` file, see Section :ref:`sec-stars`): ::

  iformat                           <=== Put this to 2 !
  nillum     nlam
  theta[1]      phi[1]
    .             .        
    .             .        
  theta[nillum] phi[nillum]
  lambda[1]
    .
    .
  lambda[nlam]
  flux[1,illum=1]
    .
    .
  flux[nlam,illum=1]
  flux[1,illum=2]
    .
    .
  flux[nlam,illum=2]
    .
    .
    .
    .
  flux[nlam,illum=nstar]

Here ``nillum`` is the number of illuminating beams you want to
specify. Normally this is 1, unless you have, e.g., a planet around a double
star. The ``theta`` is the angle (in degrees) under which the beam impinges onto
the atmosphere. If you have ``theta``\ =0, then the flux points vertically
downward (sun at zenith). If you have ``theta``\ =89, then the flux points
almost parallel to the atmosphere (sunset). It is not allowed to put ``theta``\
=90.

You can, if you wish, also put the source behind the slab, i.e. ``theta``\
>90. Please note, however, that if you compute the spectrum of the
plane-parallel atmosphere the direct flux from these illumination beams does not
get picked up in the spectrum.



Similarity and difference between 1-D spherical and 1-D plane-parallel
----------------------------------------------------------------------

Note that this 1-D plane-parallel mode is only available in :math:`z`-direction,
and only for cartesian coordinates! For spherical coordinates, a simple switch
to 1-D yields spherically symmetric 1-D radiative transfer, which is, however,
geometrically distinct from 1-D plane-parallel radiative transfer. However, you
can also use a 1-D spherically symmetric setup to 'emulate' 1-D plane parallel
problems: You can make, for instance, a radial grid in which
:math:`r_{\mathrm{nr}}/r_1-1\ll 1`. An example: :math:`r=\{10000.0`,
:math:`10000.1`, :math:`10000.2`, :math:`\cdots,` :math:`10001.0\}`. This is not
perfectly plane-parallel, but sufficiently much so that the difference is
presumably indiscernable.  The spectrum is then automatically that of the entire
large sphere, but by dividing it by the surface area, you can recalculate the
local flux.  In fact, since a plane-parallel model usually is meant to
approximate a tiny part of a large sphere, this mode is presumably even more
realistic than a truly 1-D plane-parallel model.


.. _sec-thermal-boundaries:

Thermal boundaries in Cartesian coordinates
===========================================

By default all boundaries of the computational domain are open, in the sense
that photons can move out freely. The only photons that move into the domain
from the outside are those from the interstellar radiation field (see Section
:ref:`sec-external-source`) and from any stars that are located outside of the
computational domain (see Section :ref:`sec-stars`). For some purposes it might,
however, be useful to have one or more of the six boundaries in 3-D to be
closed. RADMC-3D offers the possibility, in cartesian coordinates, to convert
the boundaries (each of the six separately) to a thermal boundary, i.e. a
blackbody emitter at some user-secified temperature. If you want that the left
X-boundary is a thermal wall at T=100 Kelvin, then you add the following line to
the ``radmc3d.inp`` file: ::

  thermal_boundary_xl = 100

and similarly for xr (right X-boundary), yl, yr, zl and/or zr. You can set this
for each boundary separately, and particularly you can choose to set just one or
just two of the boundaries to thermal boundaries. Note that setting
``thermal_boundary_xl=0`` is equivalent to switching off the thermal boundary.

Note that if you now make an image of the box, the ray-tracer will show you
still the inside of the box, through any possible thermal boundary. In other
words: for the imaging or spectra these thermal boundaries are opaque for
radiation entering the grid, while they are transparent for radiation exiting
the grid. In other words, we see the blackbody emission from the backside walls,
but not of the frontside walls. In this way we can have a look inside the box in
spite of the thermal walls.
