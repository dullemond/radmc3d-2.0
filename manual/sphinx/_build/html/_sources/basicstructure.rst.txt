.. _chap-basic-struct-and-func:

Basic structure and functionality
*********************************

RADMC-3D is a very versatile radiative transfer package with many
possibilities. As a consequence it is a rather complex package. However, we
have tried to keep it still as easy as possible to use as a first-time
user. We tried to do so by keeping many of the sophisticated options
'hidden' and having many default settings already well-chosen. The idea is
that one can already use the code at an entry level, and then gradually work
oneself into the more fancy options.

RADMC-3D is a general-purpose package, so there are no 'built-in' models
inside the ``radmc3d`` executable (Except if you insert one yourself
using the userdef module, see Chapter :ref:`chap-internal-setup`).  For
instance, if you want to model a protoplanetary disk, then you would have to
design the grid and density structure of the disk on this grid yourself. To
make it easier for the user, we have provided several Python-scripts as
examples. Among these examples is indeed a protoplanetary disk model. So
this is as close as we go to 'built-in' models: we provide, for some cases,
already well-developed example models that you, the user, can use
out-of-the-box, or that you can adapt to your needs.

In this chapter we give an overview of the rough functionality of the code
in its simplest form: ignoring all the hidden fancy options and
possibilities. For the details we then refer to the chapters ahead.

.. _sec-dataflow:

Basic dataflow
==============

Let us first clarify the basic philosophy of the code package (details will
be done later). When we talk about RADMC-3D we talk about the
fortran-90 program. The source codes are in the directory ``src/``
and the executable is called ``radmc3d``. This is the code that does
all the main calculations. You can call the code from the bash shell
(in Unix/Linux/MacOSX systems) and you can specify command-line options to
tell RADMC-3D what you want it to do.

The code RADMC-3D is in a way just a dumb computational engine. It has no
physical data (such as opacities or material properties) implemented, nor does
it have any model implemented. It is totally dependent on input files of various
kinds. These input files have filenames that end in ``.inp``, or ``.binp``,
dependent on whether the data in ASCII, or binary form. You, the user, will have
to create these input files. RADMC-3D will simply look if an ``.inp``, or a
``.binp`` file is present, and will switch to ASCII, dependent on which
file-extension it finds.

After you run RADMC-3D (by calling ``radmc3d`` with the appropriate command-line
options) you will see that the code will have produced one or more output files,
with filenames ending in ``.out`` or ``.bout``. Whether RADMC-3D produces ASCII
or binary files, depends on a flag called ``rto_style`` that you can set (see
Chapter :ref:`chap-binary-io`).

*IMPORTANT NOTE: In this manual we will mostly refer to the ASCII form
of input and output files for convenience. But each time we refer to an
\*.inp, \*.dat or \*.out file, we implicitly assume that this could also
be a \*.binp, \*.bdat or \*.bout file.*

This basic dataflow is shown in Fig. :ref:`fig-dataflow-basic`.

.. _fig-dataflow-basic:

.. figure:: Inkscape/dataflow-basic.*
   :width: 50%

   Pictographic representation of the basic dataflow of RADMC-3D. The user
   produces the input files; RADMC-3D reads them, performs the calculation,
   and produces output files. The user can then analyze the output files.

Not always can RADMC-3D produce its output files in one go. Sometimes it has to
use a two-stage procedure: For dust continuum radiative transfer the dust
temperatures are computed first (stage 1), and the images and/or spectra are
rendered after that (stage 2). Between stage 1 and stage 2 an intermediate file
is then produced (with filename ending in ``.dat`` or ``.bdat``),
which in the case of dust continuum radiative transfer is ``dust_temperature.dat``
(or ``*.bdat``).

This basic dataflow is shown in Fig. :ref:`fig-dataflow-twostage`.

.. _fig-dataflow-twostage:

.. figure:: Inkscape/dataflow-twostage.*
   :width: 85%
           
   Pictographic representation of the dataflow of RADMC-3D for the case
   of a 2-stage procedure, such as for dust continuum transfer. An intermediate
   file is produced that will be used by stage 2, but of course the user can
   also analyze the intermediate file itself.

Several of these input files contain large tables, for instance of the density
at each grid point, or the stellar flux at each wavelength bin. It is, of
course, impossible to create these datafiles by hand. The idea is that you
design a program (in any language you like) that creates these datafiles. In
that program you essentially 'program the model'. We have provided a number of
example model setups in the ``examples/`` directory. For these examples models
the setup programs were written in Python (their filenames all start with
``problem_`` and end with ``.py``). For you as the user it is therefore the
easiest to start from one of these examples and modify the Python code to your
needs. However, if you prefer to use another language, you can use the examples
to see how the input files were generated and then program this in another
programming language.

*Note: The Python files called* ``problem_*.py`` *are meant to be edited and
changed by you! They are templates from which you can create your own models.*

For the analysis of the output files created by RADMC-3D you can use your own
favorite plotting or data-analysis software. But also here we provide some tools
in Python. These Python routines are in the ``python/`` directory. Typically you
will create your own program, e.g. ``plot_model.py`` or so, that will use
these subroutines, e.g. by putting in the first line: ``from radmc3dPy import
*``. In this way Python is used also as a post-processing tool. But again: this
can also be done in another language.

This procedure is shown in Fig. :ref:`fig-dataflow-basic-python` for the
single-stage dataflow and in Fig. :ref:`fig-dataflow-twostage-python` for the
two-stage dataflow.

.. _fig-dataflow-basic-python:

.. figure:: Inkscape/dataflow-basic-python.*
   :width: 50%
           
   Pictographic representation of how the Python programs in the example directories
   are used to create the input files of RADMC-3D. 

.. _fig-dataflow-twostage-python:

.. figure:: Inkscape/dataflow-twostage-python.*
   :width: 85%
           
   Pictographic representation of the dataflow of RADMC-3D for the case
   of a 2-stage procedure, such as for dust continuum transfer. An intermediate
   file is produced that will be used by stage 2, but of course the user can
   also analyze the intermediate file itself. 


.. _sec-rad-processes:

Radiative processes
===================

Currently RADMC-3D handles the following radiative processes:

* Dust thermal emission and absorption
  
  RADMC-3D can compute spectra and images in dust continuum. The dust
  temperature must be known in addition to the dust density. In typical
  applications you will know the dust density distribution, but not the dust
  temperature, because the latter is the results of a balance between
  radiative absorption and re-emission. So in order to make spectra and
  images of a dusty object we must first calculate the dust temperature
  consistently. This can be done with RADMC-3D by making it perform a
  'thermal Monte Carlo' simulation (see Chapter :ref:`chap-dust-transfer`).
  This can be a time-consuming computation. But once this is done, RADMC-3D
  writes the resulting dust temperatures out to the file
  ``dust_temperature.dat``, which it can then later use for images and
  spectra. We can then call RADMC-3D again with the command to make an image
  or a spectrum (see Chapter :ref:`chap-dust-transfer`). To summarize: a
  typical dust continuum radiative transfer calculation goes in two stages:

  #. A thermal Monte Carlo simulation with RADMC-3D to compute the dust
     temperatures.
  #. A spectrum or image computation using ray-tracing with RADMC-3D.

* Dust scattering
  
  Dust scattering is automatically included in the thermal Monte Carlo
  simulations described above, as well as in the production of images and
  spectra. For more details, consult Chapter :ref:`chap-dust-transfer`.
  
* Gas atomic/molecular lines
  
  RADMC-3D can compute spectra and images in gas lines (see Chapter
  :ref:`chap-line-transfer`). The images are also known as *channel maps*. To
  compute these, RADMC-3D must know the population densities of the various
  atomic/molecular levels. For now there are the following options how to let
  RADMC-3D know these values:

  * Tell RADMC-3D to assume that the molecules or atoms are in *Local
    Thermodynamic Equilibrium* (LTE), and specify the gas temperature at
    each location to allow RADMC-3D to compute these LTE level populations.
    *Note that in principle one is now faced with the same problem as
    with the dust continuum: we need to know the gas temperature, which we
    typically do not know in advance.* However, computing the gas
    temperature self-consistently is very difficult, because it involves
    many heating and cooling processes, some of which are very complex.
    That is why most line radiative transfer codes assume that the user gives
    the gas temperature as input. We do so as well. If you like, you can
    tell RADMC-3D to use the (previously calculated) dust temperature as the
    gas temperature, for convenience.
    
  * Deliver RADMC-3D an input file with all the level populations
    that you have calculated youself using some method.
    
  * Tell RADMC-3D to compute the level populations according to some
    simple local non-LTE prescription such as the Sobolev approximation
    (*Large Velocity Gradient method*) or the Escape Probability Method.

  Currently RADMC-3D does not have a full non-local non-LTE computation
  method implemented. The reason is that it is very costly, and for many
  applications presumably not worth the computational effort.
  
.. _sec-coord-systems:

Coordinate systems
==================

With RADMC-3D you can specify your density distribution in two coordinate
systems:

* Cartesian coordinates: 3-D
  
  The simplest coordinate system is the Cartesian coordinate system
  :math:`(x,y,z)`. For now each model must be 3-D (i.e. you must specify the
  densities and other quantities as a function of :math:`x`, :math:`y` and :math:`z`).
  
* Cartesian coordinates: 1-D plane-parallel
  
  This is like the normal cartesian coordinates, but now the :math:`x`- and :math:`y`-
  directions are infinitely extended. Only the :math:`z`-direction has
  finite-size cells, and hence the grid is only in :math:`z`-direction.  This mode
  is the usual plane-parallel mode of radiative transfer. See Section
  :ref:`sec-1d-plane-parallel` for more details on this mode.
  
* Cartesian coordinates: 2-D pencil-parallel
  
  This is the intermediate between full 3-D cartesian and 1-D
  plane-parallel.  In this mode only the :math:`x`-direction is infinitely
  extended and a finite grid is in both :math:`y` and :math:`z` directions. This mode is
  only useful in very special cases, and is much less familiar to most - so
  use only when you are confident.

* Spherical coordinates
  
  You can also specify your model in spherical coordinates
  :math:`(r,\theta,\phi)`. These coordinates are related to the cartesian
  ones by:

  .. math::

     \begin{split}
     x &= r \sin\theta \cos\phi \\
     y &= r \sin\theta \sin\phi \\
     z &= r \cos\theta
     \end{split}

  This means that the spatial variables (density, temperature etc) are all
  specified as a function of :math:`(r,\theta,\phi)`. However, the location of the
  stars, the motion and direction of photon packages etc. are still given in
  cartesian coordinates :math:`(x,y,z)`. In other words: any function of space
  :math:`f(\vec x)` will be in spherical coordinates :math:`f(r,\theta,\phi)`, but any
  point-like specification of position :math:`\vec x` will be given as Cartesian
  coordinates :math:`\vec x=(x,y,z)`. This hybrid method allows us to do all
  physics in cartesian coordinates: photon packages or rays are treated
  always in cartesian coordinates, and so is the physics of scattering, line
  emission etc.  Only if RADMC-3D needs to know what the local conditions
  are (dust temperature, gas microturbulence, etc) RADMC-3D looks up which
  coordinates :math:`(r,\theta,\phi)` belong to the current :math:`(x,y,z)` and looks up
  the value of the density, microturbulence etc.\ at that location in the
  :math:`(r,\theta,\phi)` grid. And the same is true if RADMC-3D updates or
  calculates for instance the dust temperature: it will compute the
  :math:`(r,\theta,\phi)` belong to the current :math:`(x,y,z)` and update the
  temperature in the cell belonging to :math:`(r,\theta,\phi)`. For the rest, all
  the physics is done in the Cartesian coordinate system. This has the major
  advantage that we do not need different physics modules for cartesian and
  spherical coordinates. Most parts of the code don't care which coordinate
  system is used: they will do their own work in Cartesian coordinates.
  When using spherical coordinates, please read Section
  :ref:`sec-separable-refinement`.


.. _sec-spatial-grid:

The spatial grid
================

To specify the density or temperature structure (or any other spatial 
variable) as a function of spatial location we must have a grid. There
are two basic types of grids:

The standard gridding is a simple rectangular grid. 

* Cartesian coordinates

  When cartesian coordinates are used, this simply means that each cell is
  defined as :math:`x_l<x<x_r`, :math:`y_l<y<y_r` and :math:`z_l<z<z_r`, where
  :math:`l` and :math:`r` stand for the left and right cell walls respectively.
    
* Spherical coordinates

  When spherical coordinates are used, this simply means that each cell is
  defined as :math:`r_l<r<r_r`, :math:`\theta_l<\theta<\theta_r` and
  :math:`\phi_l<\phi<\phi_r`.  Note therefore that the shape of the cells in
  spherical coordinates is (in real space) curved. For spherical coordinates the
  following four modes are available:
    
  * 1-D Spherical symmetry:

    All spatial functions depend only on :math:`r`.
    
  * 2-D Axial symmetry:

    All spatial functions depend only on :math:`r` and :math:`\theta`.
    
  * 2-D Axial symmetry with mirror symmetry:

    All spatial functions depend only on :math:`r` and :math:`\theta`, where the
    :math:`\theta` grid only covers the part above the :math:`z=0`
    plane. Internally it is in this mode assumed that all quantities below the
    :math:`z=0` plane are equal to those above the plane by mirror symmetry in
    the :math:`z=0` plane.  This saves a factor of two in computational effort
    for Monte Carlo calculations, as well as in memory useage. Note that also
    the resulting output files such as ``dust_temperature.dat`` will only be
    specified for :math:`z>0`.
      
  * 3-D:
    
    All spatial functions depend on all three variables
    :math:`r`, :math:`\theta` and :math:`\phi`.
          
  * 3-D with mirror symmetry:

    All spatial functions depend on all three variables :math:`r`,
    :math:`\theta` and :math:`\phi`, but like in the 2-D case only the upper
    part of the model needs to be specified: the lower part is assumed to be a
    mirror copy.

  When using spherical coordinates, please read Section
  :ref:`sec-separable-refinement`.
         
In all cases these structured grids allow for oct-tree-style grid refinement, or
its simplified version: the layer-style grid refinement. See Section
:ref:`sec-grid-input` and Chapter :ref:`chap-gridding` for more information
about the gridding and the (adaptive) mesh refinement (AMR).

.. _sec-actions:

Computations that RADMC-3D can perform
======================================

The code RADMC-3D (i.e. the executable ``radmc3d``) is *one* code for *many*
actions. Depending on which command-line arguments you give, RADMC-3D can do
various actions. Here is a list:

#. Compute the dust temperature:
   
   With ``radmc3d mctherm`` you call RADMC-3D with the command of performing a
   thermal Monte Carlo simulation to compute the dust temperature under the
   assumption that the dust is in radiative equilibrium with its radiation
   field. This is normally a prerequisite for computing SEDs and images from
   dusty objects (see *computing spectra and images* below).  The output file of
   this computation is ``dust_temperature.dat`` which contains the dust
   temperature everywhere in the model.
   
#. Compute a spectrum or SED:
   
   With ``radmc3d sed`` you call RADMC-3D with the command of performing a
   ray-tracing computation to compute the spectral energy distribution (SED) for
   the model at hand. Typically you first need to have called ``radmc3d
   mctherm`` (see above) beforehand to compute dust temperatures (unless you
   have created the file ``dust_temperature.dat`` yourself because you have a
   special way of computing the dust temperature). With ``radmc3d sed`` the
   spectrum is computed for the wavelengths points given in the file
   ``wavelength_micron.inp``, which is the same wavelength grid that is used for
   ``radmc3d mctherm``. If you want to compute the spectrum at wavelength other
   than those used for the thermal Monte Carlo simulation, you should instead
   call ``radmc3d spectrum``, and you have the full freedom to choose the
   spectral wavelengths points at will, and you can specify these in various
   ways described in Section :ref:`sec-set-camera-frequencies`.  Most easily you
   can create a file called ``camera_wavelength_micron.inp`` (see Section
   :ref:`sec-camera-wavelengths`) and call RADMC-3D using ``radmc3d spectrum
   loadlambda``. In all these cases the vantage point (where is the observer)
   can of course be set as well, see Section :ref:`sec-dust-ray-tracing` and
   Chapter :ref:`chap-images-spectra`.
       
#. Compute an image:
   
   With ``radmc3d image`` you call RADMC-3D with the command of performing a
   ray-tracing computation to compute an image. You must specify the
   wavelength(s) at which you want the image by, for instance, calling RADMC-3D
   as ``radmc3d image lambda 10``, which makes the image at
   :math:`\lambda=10\mu\mathrm{m}`. But there are other ways by which the wavelength(s) can be set, see
   Section :ref:`sec-set-camera-frequencies`.  In all these cases the vantage
   point (where is the observer) can of course be set as well, see Section
   :ref:`sec-dust-ray-tracing` and Chapter :ref:`chap-images-spectra`.
       
#. Compute the local radiation field inside the model:
   
   With ``radmc3d mcmono`` you call RADMC-3D with the command of performing a
   wavelength-by-wavlength monochromatic Monte Carlo simulation (at the
   wavelengths that you specify in the file
   ``mcmono_wavelength_micron.inp``). The output file of this computation is
   ``mean_intensity.out`` which contains the mean intensity :math:`J_\nu` as a
   function of the :math:`(x,y,z)` (cartesian) or :math:`(r,\theta,\phi)`
   (spherical) coordinates at the frequencies :math:`\nu_i\equiv
   10^4c/\lambda_i` where :math:`\lambda_i` are the wavelengths (in
   :math:`\mu`\ m) specified in the file ``mcmono_wavelength_micron.inp``. The
   results of this computation can be interesting for, for instance, models of
   photochemistry. But if you use RADMC-3D only for computing spectra and
   images, then you will not use this.

In addition to the above main methods, you can ask RADMC-3D to do various minor
things as well, which will be described throughout this manual.


How a model is set up and computed: a rough overview
====================================================

A radiative transfer code such as RADMC-3D has the task of computing synthetic
images and spectra of a model that you specify. You tell the code what the dust
and/or gas density distribution in 3-D space is and where the star(s) are, and
the code will then tell you what your cloud looks like in images and/or
spectra. That's basically it. That's the main task of RADMC-3D.

First you have to tell RADMC-3D what 3-D distribution of dust and/or gas you
want it to model. For that you must specify a coordinate system (cartesian or
spherical) and a spatial grid. For cartesian coordinates this grid should be 3-D
(although there are exceptions to this), while for spherical coordinates it can
be 1-D (spherical symmetry), 2-D (axial symmetry) or 3-D (no symmetry). RADMC-3D
is (for most part) a cell-based code, i.e. your grid devides space in cells and
you have to tell RADMC-3D what the average densities of dust and/or gas are in
these cells.

The structure of the grid is specified in a file ``amr_grid.inp`` (see Section
:ref:`sec-grid-input`). All the other data, such as dust density and/or gas
density are specified in other files, but all assume that the grid is given by
``amr_grid.inp``.

We can also specify the locations and properties of one or more stars in the
model. This is done in the ``stars.inp`` (see Section :ref:`sec-stars`) file.

Now suppose we want to compute the appearance of our model in dust continuum. We
will describe this in detail in Chapter :ref:`chap-dust-transfer`, but let us
give a very rough idea here. We write, in addition to the ``amr_grid.inp`` and
``stars.inp`` files, a file ``dust_density.inp`` which specifies the density of
dust in each cell (see Section :ref:`sec-dustdens`).  We also must write the
main input file ``radmc3d.inp`` (see Section :ref:`sec-radmc-inp`), but we can
leave it empty for now. We must give RADMC-3D a dust opacity table in the files
``dustopac.inp`` and for instance ``dustkappa_silicate.inp`` (see Section
:ref:`sec-opacities`). And finally, we have to give RADMC-3D a table of discrete
wavelengths in the file ``wavelength_micron.inp`` that it will use to perform
its calculations on. We then call the ``radmc3d`` code with the keyword
``mctherm`` (see Chapter :ref:`chap-dust-transfer`) to tell it to perform a
Monte Carlo simulation to compute dust temperatures everywhere. RADMC-3D will
write this to the file ``dust_temperature.dat``. If we now want to make a
spectral energy distribution, for instance, we call ``radmc3d sed`` (see Section
:ref:`sec-making-spectra`) and it will write a file called ``spectrum.out``
which is a list of fluxes at the discrete wavelengths we specified in
``wavelength_micron.inp``.  Then we are done: we have computed the spectral
energy distribution of our model. We could also make an image at wavelength 10
:math:`\mu`\ m for instance with ``radmc3d image lambda 10`` (see Section
:ref:`sec-images`). This will write out a file ``image.out`` containing the
image data (see Section :ref:`sec-image-out`).

As you see, RADMC-3D reads all its information from tables in various
files. Since you don't want to make large tables by hand, you will have to write
a little computer program that generates these tables automatically.  You can do
this in any programming language you want. But in the example models (see
Section :ref:`sec-example-models`) we use the programming language Python (see
Section :ref:`sec-requirements`) for this. It is easiest to indeed have a look
at the example models to see how this is (or better: can be) done.

We will explain all these things in much more detail below, and we will discuss
also many other radiative transfer problem types. The above example is really
just meant to give an impression of how RADMC-3D works.


.. _sec-rough-overview-models:

Organization of model directories
=================================

The general philosophy of the RADMC-3D code package is the following. The core
of everything is the fortran code ``radmc3d``. This is the main code which does
the hard work for you: it makes the radiative transfer calculations, makes
images, makes spectra etc. Normally you compile this code just once-and-for-all
(see Chapter :ref:`chap-compilation`), and then simply use the executable
``radmc3d`` for all models. There is an exception to this 'once-and-for-all'
rule described in Section :ref:`sec-special-purpose-compile`, but in the present
chapter we will not use this (see Chapter :ref:`chap-internal-setup` for this
instead). So we will stick here to the philosophy of compiling this code once
and using it for all models.

So how to set up a model? The trick is to present ``radmc3d`` with a set of
input files in which the model is described in all its details. The procedure to
follow is this:

#. The best thing to do (to avoid a mess) is to make a directory for 
   *each model*: one model, one directory. Since ``radmc3d`` reads
   multiple input files, and also outputs a number of files, this is a good
   way to keep organized and we recommend it strongly.  So if we wish to make
   a new model, we make a new directory, or copy an old directory to a new
   name (if we merely want to make small changes to a prior model).
   
#. In this directory we generate the input files according to their required
   format (see Chapter :ref:`chap-input-files`). You can create these input files
   in any way you want. But since many of these input files will/must contain
   huge lists of numbers (for instance, giving the density at each location in
   your model), you will typically want to write some script or program in some
   language (be it either C, C++, Fortran, IDL, GDL, perl, python, you name it)
   that automatically creates these input files. *We recommend using Python,
   because we provide examples and standard subroutines in the programming
   language Python; see below for more details.*  Section
   :ref:`sec-example-models` describes how to use the example Python scripts to
   make these input files with Python.
   
#. When all the input files are created, and we make sure that we are inside the
   model directory, we call ``radmc3d`` with the desired command-line options
   (see Chapter :ref:`chap-command-line-options`). This will do the work for us.
   
#. Once this is done, we can analyze the results by reading the output files
   (see Chapter :ref:`chap-input-files`). To help you reading and analyzing
   these output files you can use a set of Python routines that we created for
   the user (see Chapter :ref:`chap-python-analysis-tools` and Section
   :ref:`sec-install-pythonscripts`). But here again, you are free to use any
   other plotting software and/or data postprocessing packages.


.. _sec-example-models:

Running the example models
==========================

Often the fastest and easiest way to learn a code is simply to analyze and run a
set of example models. They are listed in the ``examples`` directory. Each model
occupies a separate directory. This is also the style we normally recommend:
each model should have its own directory. Of course there are also exceptions to
this rule, and the user is free to organize her/his data in any way he/she
pleases. But in all the examples and throughout this manual each model has its
own directory.

To run an example model, go into the directory of this model, and follow the
directions that are written in the ``README`` file in each of these
directories. *This is under the assumption that you have a full Python
distribution installed on your system, including Numpy and Matplotlib.*

Let us do for instance ``run_simple_1/``::

  cd examples/run_simple_1

Now we must create all the input files for this model. These input files are
all described in chapter :ref:`chap-input-files`, but let us here just
'blindly' follow the example. In this example most (all except one) of the
input files are created using a Python script called ``problem_setup.py``.
To execute this script, this is what you do on the shell::

  python problem_setup.py

This Python script has now created a whole series
of input files, all ending with the extension ``.inp``. To see which
files are created, type the following in the shell::

  ls -l *.inp

There is one file that this example does not create, and that is the file
``dustkappa_silicate.inp``. This is a file that contains the dust opacity in
tabulated form. This is a file that you as the user should provide to the
RADMC-3D code package. The file ``dustkappa_silicate.inp`` is merely an example,
which is an amorphous spherical silicate grain with radius 0.1 micron. But see
Section :ref:`sec-opacities` for more information about the opacities.

Now that the input files are created, we must run ``radmc3d``::

  radmc3d mctherm

This tells RADMC-3D to do the thermal Monte Carlo simulation. This may
take some time. When the model is ready, the prompt of the shell returns. 
To see what files have been created by this run of the code, type::

  ls -l *.dat

You will find the ``dust_temperature.dat`` containing the dust temperature
everywhere in the model. See again chapter :ref:`chap-input-files` for
details of these files. To create a spectral energy distribution (SED)::

  radmc3d sed incl 45.

This will create a file ``spectrum.out``.  To analyze these data you can use the
Python routines delivered with the code (see Chapter
:ref:`chap-python-analysis-tools` and Section :ref:`sec-install-pythonscripts`).

There is a file ``Makefile`` in the directory. This is here only meant to make
it easy to clean the directory. Type ``make cleanmodel`` to clean all the output
from the radmc3d code. Type ``make cleanall`` to clean the directory back to
basics.

Let us now do for instance model ``run_simple_1_userdef/``::

  cd examples/run_simple_1_userdef

This is the same model as above, but now the grid and the dust density are set
up *inside* ``radmc3d``, using the file ``userdef_module.f90`` which is
present in this directory.  See Chapter :ref:`chap-internal-setup` for details
and follow the directions in the ``README`` file. In short: first edit the
variable ``SRC`` in the ``Makefile`` to point to the ``src/`` directory. Then
type ``make``. Then type ``python problem_setup.py`` on the shell command line
(which now only sets up the frequency grid, the star and the ``radmc3d.inp``
file and some small stuff). Now you can run the model.

*Please read the README file in each of the example model directories.
Everything is explained there, including how to make the relevant plots.*
