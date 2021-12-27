.. _chap-radmc3d-internal-analysis-tools:

Analysis tools inside of radmc3d
********************************

There are also some special purpose features in the Fortran-90 ``radmc3d`` code that can be useful for analyzing complex AMR-gridded models.


.. _sec-subbox:

Making a regularly-spaced datacube ('subbox') of AMR-based models
=================================================================


Because handling AMR-based models in Python or other data analysis packages can
be rather cumbersome, we decided that it would be useful to create the
possibility in ``radmc3d`` to generate 1-D, 2-D or 3-D regularly spaced
'cut-outs' or 'sub-boxes' (whatever you want to call them) of any variable of
the model. 


Creating a subbox
-----------------------------

You can call ``radmc3d`` directly from the shell asking it to make
the subbox. Here is an example::

  ./radmc3d subbox_dust_density subbox_nxyz 64 64 64 subbox_xyz01 -2.e15 2.e15 -2.e15 2.e15 -2.e15 2.e15

which creates a regularly sampled 64x64x64 datacube of the dust density, with :math:`x` grid
between :math:`-2\times 10^{15}\;\mathrm{cm}` and  :math:`+2\times 10^{15}\;\mathrm{cm}` and
likewise for :math:`y` and :math:`z` (note that these box boundaries are the walls of the
regularly spaced cells of the subbox). The file that this creates is called ``dust_density_subbox.out``
(see section :ref:`sec-subbox-file-format` for the format of this file).
For the dust temperature the command is
``./radmc3d subbox_dust_temperature``, in which case the file is called
``dust_temperature_subbox.out``.

You can also rotate the box along three angles: :math:`\phi_1`, :math:`\theta`,
and :math:`\phi_2`, for example::
  
  ./radmc3d subbox_dust_temperature subbox_nxyz 64 64 64 subbox_xyz01 -2.e15 2.e15 -2.e15 2.e15 -2.e15 2.e15 subbox_phi1 30 subbox_theta 60  subbox_phi2 45

(Note that as of version 2.0 of RADMC-3D these angles are in degrees instead of radian).
An example for the level populations would be::

  ./radmc3d subbox_levelpop subbox_nxyz 64 64 64 subbox_xyz01 -2.e15 2.e15 -2.e15 2.e15 -2.e15 2.e15

*Note about subbox for level populations:* By default all level populations will
be written out. However, if you would add the ``subbox_levelpop`` keyword in a
call to RADMC-3D for making an image or spectrum, then it will only write out
the level populations that have been used for that image. Example: ::

  ./radmc3d image lambda 2600 subbox_levelpop subbox_nxyz 64 64 64 subbox_xyz01 -2.e15 2.e15 -2.e15 2.e15 -2.e15 2.e15

would give a much smaller ``'levelpop_co_subbox.out'`` file, because only the
first two levels are included (remember that :math:`\lambda=2600\,\mu`\ m is the
J1-0 line of CO). See Section :ref:`sec-calcstore-levpop` for more information
on how RADMC-3D automatically selects a subset of levels for storage in the
global array (and thus also for writing out to file).



.. _sec-subbox-file-format:

Format of the subbox output files
---------------------------------

All the files produced by the subbox method have the following format:
::

  iformat                                  <=== Typically 2 at present
  nx ny nz nv                              <=== Box of nx*ny*nz cells, each with nv values
  x0 x1 y0 y1 z0 z1                        <=== The x, y and z boundaries of the box
  phi1 theta phi2                          <=== Three rotation angles of the box
                                           <=== Empty line 
  1 2 3 4 ....                             <=== Identifications of the nv values 
                                           <=== Empty line 
  data[ix=1,iy=1,iz=1,iv=1]
  data[ix=2,iy=1,iz=1,iv=1]
  .
  .
  data[ix=nx,iy=1,iz=1,iv=1]
  data[ix=1,iy=2,iz=1,iv=1]
  .
  .
  .
  data[ix=nx,iy=ny,iz=nz,iv=1]
                                           <=== Empty line between components
  data[ix=1,iy=1,iz=1,iv=2]
  .
  .
  data[ix=nx,iy=ny,iz=nz,iv=2]
                                           <=== Empty line between components
  .
  .
  .
                                           <=== Empty line between components
  data[ix=1,iy=1,iz=1,iv=nv]
  .
  .
  data[ix=nx,iy=ny,iz=nz,iv=nv]

and they are always in ascii format. For a subbox of the level populations the
identification numbers are the levels. For instance, if only the populations of
levels 4 and 8 are in this file, then ``nv=2`` and the line with
the identification numbers will be ``48``\ . For all other quantities
(dust density, dust temperature) this line of identification numbers is simply
``123`` etc.

Using the ``radmc3d_tools`` to read the subbox data
---------------------------------------------------

In Section :ref:`sec-simpleread-tools` a set of simple Python tools are
discussed to read a variety of output files from RADMC-3D (as well as input
files to RADMC-3D) for further analysis.

Also for the subbox output there is now a Python function to read those.
Example: First run RADMC-3D:
::
   radmc3d mctherm
   radmc3d subbox_dust_density subbox_nxyz 64 64 64 subbox_xyz01 -2.e14 2.e14 -2.e14 2.e14 -2.e14 2.e14
   radmc3d subbox_dust_temperature subbox_nxyz 64 64 64 subbox_xyz01 -2.e14 2.e14 -2.e14 2.e14 -2.e14 2.e14

Then go into Python and do:
::
   from radmc3d_tools.simpleread import *
   dustdens = read_subbox(name='dust_density')
   dusttemp = read_subbox(name='dust_temperature')
   grid     = dustdens.grid
   import matplotlib.pyplot as plt
   rhodustmin = 1e-18
   plt.figure()
   plt.imshow(np.log10(dustdens.data[:,:,32]+rhodustmin),extent=[grid.x[0],grid.x[-1],grid.y[0],grid.y[-1]])
   plt.figure()
   plt.imshow(dusttemp.data[:,:,32])
   plt.show()

.. _sec-sampling:

Alternative to subbox: arbitrary sampling of AMR-based models
=============================================================

For some purposes it is useful to sample values of various quantities at
arbitrary positions in the grid. The idea is very much like the subbox
method of Section :ref:`sec-subbox`, but instead of a regular subbox grid
the user provides a list of 3-D points where he/she wants to sample the
variables of the model. Here is how to do this. First you must produce
a file containing the list of 3-D positions. The file is called
``sample_points.inp`` and is an ascii file that looks as
follows:
::

  iformat                                  <=== Typically 1 at present
  npt                                      <=== Nr of 3-D sampling points
  xpt[1]  ypt[1]  zpt[1]                   <=== 3-D coordinates of point 1
  xpt[2]  ypt[2]  zpt[2]                   <=== 3-D coordinates of point 2
  xpt[3]  ypt[3]  zpt[3]                   <=== 3-D coordinates of point 3
  ...
  ...

An example for the case in which you want to sample at just one point:
::

  1
  1
  1.49d13   4.02d14   1.03d12

If you want to let RADMC-3D do the sampling of the dust density and
temperature, type (after you have calculated the temperature using
``radmc3dmctherm``\ ):
::

  radmc3d sample-dustdens sample-dusttemp

You can also do the dust temperature calculation and the sampling in one
go:
::

  radmc3d mctherm sample-dustdens sample-dusttemp

You can also do only ``sample-dusttemp`` or only ``sample-dustdens``\ . The
output is written to files ``dust_density_sample.out`` resp.\
``dust_temperature_sample.out``\ . The format of these files is (take dust
density as example): ::

  iformat                                  <=== Typically 2 at present
  npt  nv                                  <=== Nr of point and size of datavector
                                           <=== Empty line
  1 2 3 4 ....                             <=== Identifications of the nv values 
                                           <=== Empty line
  dustdensity[ipt=1,iv=1]
  dustdensity[ipt=2,iv=1]
  ...
  dustdensity[ipt=npt,iv=1]
                                           <=== Empty line between components
  dustdensity[ipt=1,iv=2]
  ...
  dustdensity[ipt=npt,iv=2]
                                           <=== Empty line between components
  ...
                                           <=== Empty line between components
  dustdensity[ipt=npt,iv=nv]

where ``nv`` is in this case the nr of species of dust and 
``iv``\ =``ispecies``\ .

For a sample of the level populations the identification numbers are the
levels. For instance, if only the populations of levels 4 and 8 are in this
file, then ``nv=2`` and the line with the identification numbers
will be ``48``\ . For all other quantities (dust density, dust
temperature) this line of identification numbers is simply ``123``
etc.

Later we will add other possible arrays to sample (at the moment it is only
dust density, dust temperature and level populations). But you can also
implement this yourself. Search in the following files for the following
parts to add your own sampling:

* In ``rtglobal_module.f90``\ : Search for ``do_sample_dustdens`` and add your
  own variable, e.g. ``o_sample_myvariable``\ .
  
* In ``main.f90``\ : Search for ``do_sample_dustdens`` and you will find all
  places where you have to add your own stuff, i.e.  where you will have to add
  statements like ``if(do_sample_myvariable)`` or where you have to set
  ``do_sample_myvariable=.true.`` or reset ``do_sample_myvariable=.false.`` etc.

That should do it.
