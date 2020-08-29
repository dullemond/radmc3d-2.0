.. _chap-internal-setup:

Modifying RADMC-3D: Internal setup and user-specified radiative processes
*************************************************************************

It has been mentioned several times before that as an alternative to the
standard `compile once-and-for-all' philosophy, one can also use RADMC-3D by
modifying the code directly so that ``radmc3d`` will have new
functionality that might be of use for you. We refer to Section
:ref:`sec-special-purpose-compile` for an in-depth description of how to
modify the code in a way that is *non-invasive* to the main code. We
urge the reader to read Section :ref:`sec-special-purpose-compile` first
before continuing to read this chapter. In all of the following we assume
that the editings to the fortran files are done in the local way described
in Section :ref:`sec-special-purpose-compile` so that the original source
files in the ``src/`` directory stay unaffected, and only local
copies are edited.


Setting up a model *inside* of RADMC-3D
===========================================

The most common reason for editing the code itself is for setting up the
model *internally* rather than reading in all data via input files. For
a list of advantages and disadvantages of setting models up internally as
opposed to the standard way, see Section :ref:`sec-internalsetup-proscons`
below. Setting up a model within RADMC-3D is done by making a local copy of
the file ``userdef_module.f90`` and editing it (see Section
:ref:`sec-special-purpose-compile`). This file contains a set of standard
subroutines that are called by the main program at special points in the
code. Each subroutine has a special purpose which will be described below.
By keeping a subroutine empty, nothing is done. By filling it with your own
code lines, you can set up the density, temperature or whatever needs to be
set up for the model. In addition to this you can do the following as well:

* Add new variables or arrays in the module header (above the ``contains``
  command), which you can use in the subroutines of the ``userdef_module.f90``
  module. You are completely free to add any new variables you like. A small
  tip: it may be useful (though not required) to start all their names with
  e.g. ``userdef_`` to make sure that no name conflicts with other variables in
  the code happen.
  
* Add new subroutines at will (below the ``contains`` command) which you can
  call from within the standard subroutines.
  
* Introduce your own ``radmc3d`` command-line options (see Section
  :ref:`sec-predef-userdef`).
  
* Introduce your own ``radmc3d.inp`` namelist variables (see Section
  :ref:`sec-predef-userdef`).


.. _fig-dataflow-basic-userdef:

.. figure:: Inkscape/dataflow-basic-userdef.*

   Pictographic representation of the dataflow for the case when you define your
   model *internally* using the ``userdef_module.f90``\ .

Often you still want some of the input data to be still read in in the usual
way, using input files. For instance, you may want to still read the
``dustopac.inp`` and the opacities using the ``dustkappa_xxx.inp`` files. This
is all possible. Typically, you simply keep the files you still want RADMC-3D to
read, and omit the files that contain data that you allocate and set in the
``userdef_module.f90``\ . This is all a bit complicated, so the best way to
learn how to do this is to start from the example directories in which a model
is set up with the ``userdef_module.f90`` method.

In Fig. :numref:`fig-dataflow-basic-userdef` the dataflow for the user-defined
model setup is graphically depicted.


.. _sec-predef-userdef:

The pre-defined subroutines of the userdef_module.f90
======================================================


The idea of the ``userdef_module.f90`` is that it contains a number of standard
pre-defined subroutines that are called from the ``main.f90`` code (and *only*
from there). Just browse through the ``main.f90`` file and search for the
sequence ``calluserdef_`` and you will find all the points where these standard
routines are called. It means that at these points you as the user have
influence on the process of model setup. Here is the list of standard routines
and how they are used. They are ordered roughly in chronological order in which
they are called.

* ``userdef_defaults()``
  
  This subroutine allows you to set the default value of any new parameters you
  may have introduced. If neither on the command line nor in the ``radmc3d.inp``
  file the values of these parameters are set, then they will simply retain this
  default value.
  
* ``userdef_commandline(buffer,numarg,iarg,fromstdi,gotit)``
  
  This subroutine allows you to add your own command-line options for
  ``radmc3d``\ . The routine has a series of standard arguments which you are
  not allowed to change. The ``buffer`` is a string containing the current
  command line option that is parsed. You will check here if it is an option of
  your module, and if yes, activate it.  An example is listed in the code. You
  an also require a second argument, for which also an example is listed in the
  original code.
  
* ``userdef_commandline_postprocessing()``
  
  After the command line options have been read, it can be useful to check if
  the user has not asked for conflicting things. Here you can do such checks.
  
* ``userdef_parse_main_namelist()``
  
  Here you can add your own namelist parameters that read from the
  ``radmc3d.inp`` file. An example is provided in the original code.
  
* ``userdef_main_namelist_postprocessing()``
  
  Also here, after the entire ``radmc3d.inp`` file has been read and
  interpreted, you can do some consistency checks and postprocessing here.
  
* ``userdef_prep_model()``
  
  This routine can be used if you wish to set up the grid not from input files
  but internally. You will have to know how to deal with the ``amr_module.f90``
  module. You can also set your own global frequency grid here. And finally, you
  can set your own stellar sources here. In all cases, if you set these things
  here (which requires you to make the proper memory allocations, or in case of
  the gridding, let the ``amr_module.f90`` do the memory allocations for you)
  the further course of ``radmc3d`` will skip any of its own settings (it will
  simply detect if these arrays are allocated already, and if yes, it will
  simply not read or allocate them anymore).
  
* ``userdef_setup_model()``
  
  This is the place where you can actually make your own model setup.  By the
  time this subroutine is called, all your parameters have been read in, as well
  as all of the other parameters from the original ``radmc3d`` code. So you can
  now set up the dust density, or the gas velocity or you name it. For all of
  these things you will have to allocate the arrays youself (!!!). Once you did
  this, the rest of the ``radmc3d`` code won't read those data anymore, because
  it detects that the corresponding arrays have already been allocated (by
  you). This allows you to completely circumvent the reading of any of the
  following files by making these data yourself here at this location:

    * ``amr_grid.inp`` or in the future the input files for any of the other griding types.
    * ``dust_density.inp``
    * ``dust_temperature.dat``
    * ``gas_density.inp``
    * ``gas_temperature.inp`` 
    * ``gas_velocity.inp``
    * ``microturbulence.inp``
    * ``levelpop_XXX.dat`` 
    * ``numberdens_XXX.inp``

  To learn how to set up a model in this way, we refer you for now to the
  ``ioput_module.f90`` or ``lines_module.f90`` and search for the above file
  names to see how the arrays are allocated and how the data are inserted. I
  apologise for not explaining this in more detail at this point. But examples
  are or will be given in the ``examples/`` directory.
  
* ``userdef_dostuff()``
  
  This routine will be called by the main routine to allow you to do any kind of
  calculation after the main calculation (for instance after the monte carlo
  simulation). This is done within the execution-loop.
  
* ``userdef_compute_levelpop()``
  
  This is a subroutine that can be called by the camera module for
  on-the-fly calculation of level populations according to your own recipe.
  This may be a bit tricky to use, but I hope to be able to provide some
  example(s) in the near future.
  
* ``userdef_srcalp()``
  
  This subroutine allows you to add any emission/absorption process you
  want, even fake ones. For instance, you could use this to create nicely
  volume-rendered images of your 3-D models with fake opacities, which are
  chosen to make the image look nice and/or insight-giving.  You can also
  use this to add physical processes that are not yet implemented in
  RADMC-3D. This subroutine allows you full freedom and flexibility to 
  add emissivity and extinction whereever/however you like. To activate
  it you must set ``incl_userdef_srcalp=1`` in the
  ``radmc3d.inp`` file.
  
* ``userdef_writemodel()``
  
  This allows the user to dump any stuff to file that the user computed
  in this module. You can also use this routine to write out files that would
  have been used normally as input file (like ``amr_grid.inp`` or
  ``dust_density.inp``\ ) so that the Python routines can read them if
  they need. In particular the grid information may be needed by these
  external analysis tools. Here is a list of standard subroutines you can
  call for writing such files:

    * ``write_grid_file()``
    * ``write_dust_density()``
    * ...more to come...
  
For now this is it, more routines will be included in the future.

Note that the ``userdef_compute_levelpop()`` subroutine, in contrast to all the
others, is called not from the ``main.f90`` program but from the
``camera_module.f90`` module. This is why the camera module is the only module
that is higher in compilation ranking than the userdef module (i.e. the userdef
module will be compiled before the camera module). For this reason the userdef
module has no access to the variables of the camera module. For the rest, the
userdef module has access to the variables in all other modules.

Note also that not all input data is meant to be generated in this way. The
following types of data are still supposed to be read from file:

  * Dust opacity data
  * Molecular fundamental data

Please have a look in the ``examples/`` directory for models 
which are set up in this internal way.



.. _sec-internalsetup-proscons:

Some caveats and advantages of internal model setup
===================================================

Setting up the models internally has several advantages as well as
disadvantages compared to the standard way of feeding the models into 
``radmc3d`` via files. The advantages are, among others:

* You can modify the model parameters in ``radmc3d.inp`` and/or in the command
  line options (depending on how you allow the user to set these parameters,
  i.e. in the ``userdef_parse_main_namelist()`` routine and/or in the
  ``userdef_commandline()`` routine. You then do not need to run Python anymore
  (except for setting up the basic files; see examples). Some advantages of
  this:

  * It allows you, for instance, to create a version of the ``radmc3d`` code
    that acts as if it is a special-purpose model. You can specify model
    parameters on the command line (rather than going through the cumbersome
    Python stuff).
    
  * It is faster: even a large model is built up quickly and does not
    require a long read from large input files.
    
* You can make use of the AMR module routines such as the
  ``amr_branch_refine()`` routine, so you can adaptively refine the grid while
  you are setting up the model.

Some of the disadvantages are:

* The model needs to be explicitly written out to file and read into Python or
  any other data plotting package before you can analyze the density structure
  to test if you've done it right. You can explicitly ask ``./radmc3d`` to call
  the ``userdef_writemodel()`` subroutine (which is supposed to be writing out
  all essential data; but that is the user's responsibility) by typing
  ``./radmc3dwritemodel``\ .
  
* Same is true for the grid, and this is potentially even more dangerous if not
  done. You can explicitly ask ``./radmc3d`` to write out the grid file by
  typing ``./radmc3dwritegridfile``\ .  Note that if you call the
  ``write_grid_file()`` subroutine from within ``userdef_writemodel()``\ , then
  you do not have to explicitly type ``./radmc3dwritegridfile`` as well.  Note
  also that ``radmc3d`` will automatically call the ``write_grid_file()``
  subroutine when it writes the results of the thermal Monte Carlo computation,
  if it has its grid from inside (i.e. it has not read the grid from the file
  ``amr_grid.inp``\ .
  
* It requires a bit more knowledge of the internal workings of the ``radmc3d``
  code, as you will need to directly insert code lines in the
  ``userdef_module.f90`` file.




.. _sec-compute-radiation-integrals:

Using the userdef module to compute integrals of :math:`J_\nu`
==============================================================


With the monochromatic Monte Carlo computation (see Section
:ref:`sec-dust-monochromatic-monte-carlo`) we can calculate the mean intensity
:math:`J_\nu` at every location in the model at a user-defined set of
wavelengths. However, as mentioned before, for large models and large numbers of
wavelengths this could easily lead to a data volume that is larger than what the
computer can handle. Since typically the main motivation for computing
:math:`J_\nu` is to compute some integral of the the form:

.. math::

  Q = \int_0^{\infty} J_\nu K_\nu d\nu

where :math:`K_\nu` is some cross section function or so, it may not be
necessary to store the entire function :math:`J` as a function of :math:`nu`.
Instead we would then only by interested in the result of this integral
at each spatial location. 

So it would be useful to allow the user to do this computation internally.  We
should start by initializing :math:`Q(x,y,z)=0` (or :math:`Q(r,\theta,\phi)=0`
if you use spherical coordinates). Then we call the monochromatic Monte Carlo
routine for the first wavelength we want to include, and multiply the resulting
mean intensities with an appropriate :math:`\Delta\nu` and add this to
:math:`Q(x,y,z)`. Then we do the monochromatic Monte Carlo for the next
wavelength and again add to :math:`Q` everywhere. We repeat this until our
integral (at every spatial location on the grid) is finished, and we are
done. This saves a huge amount of memory.

Since this is somewhat hard to explain in this PDF document, we refer to
the example model ``run_example_jnu_integral/``\ .

*STILL IN PROGRESS.*



Some tips and tricks for programming user-defined subroutines
=============================================================

Apart from the standard subroutines that *must* be present in the
``userdef_module.f90`` file (see Section :ref:`sec-predef-userdef`), you are
free to add any subroutines or functions that you want, which you can call from
within the predefined subroutines of Section :ref:`sec-predef-userdef`. You are
completely free to expand this module as you wish. You can add your own
variables, your own arrays, allocate arrays, etc.

Sometimes you may need to know 'where you are' in the grid. For instance, the
subroutine ``userdef_compute_levelpop()`` is called with an argument ``index``\
. This is the index of the current cell from within which the subroutine has
been called. You can now address, for instance, the dust temperature at this
location: ::

  temp = dusttemp(1,index)

(for the case of a single dust species). You may also want to know the
coordinates of the center of the cell. For this, you must first get a pointer to
the AMR-tree structure of this cell. The pointer ``b`` is declared as ::

  type(amr_branch), pointer :: b

Then you can point the pointer to that cell structure
::

  b => amr_index_to_leaf(index)%link

And now you can get the x,y,z-coordinates of the center of the cell:
::

  xc = amr_finegrid_xc(b%ixyzf(1),1,b%level)
  yc = amr_finegrid_xc(b%ixyzf(2),2,b%level)
  zc = amr_finegrid_xc(b%ixyzf(3),3,b%level)

Or the left and right cell walls:
::

  xi_l = amr_finegrid_xi(b%ixyzf(1),1,b%level)
  yi_l = amr_finegrid_xi(b%ixyzf(2),2,b%level)
  zi_l = amr_finegrid_xi(b%ixyzf(3),3,b%level)
  xi_r = amr_finegrid_xi(b%ixyzf(1)+1,1,b%level)
  yi_r = amr_finegrid_xi(b%ixyzf(2)+1,2,b%level)
  zi_r = amr_finegrid_xi(b%ixyzf(3)+1,3,b%level)



Creating your own emission and absorption processes
===================================================

RADMC-3D Allows you to add your own physics to the ray-tracing images and
spectra. At every point during the ray-tracing process, when it computes the
emissivity and extinction coefficients :math:`j_\nu` and :math:`\alpha_\nu` it
calls the ``userdef_srcalp()`` subroutine, giving it the ``index`` in which cell
we are, the frequencies of the different image channels and the ``src`` and
``alp`` arrays which are for resp.\ :math:`j_\nu` and :math:`\alpha_\nu`. You
can *add* any process by ::

  src(:) = src(:) + .....
  alp(:) = alp(:) + .....

where ...... is your formula. You can find the local variables like 
density and temperature using the ``index``\ , e.g.:
::

  rho_g = gasdens(index)

You can be completely free in your choices. If you need some information
that is not usually read into RADMC-3D, you can add read commands in the
``userdef_setup_model()`` subroutine, e.g.:
::

  call read_gas_density(1)


See the example directory ``examples/run_simple_userdefsrc`` for
more ideas.


