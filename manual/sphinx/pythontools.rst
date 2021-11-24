.. _chap-python-analysis-tools:

Python analysis tool set
************************

While the code RADMC-3D is written in fortran-90, there is an extensive set
of tools written in Python that make it easier for the user to set up models
and interpret results. See Section :ref:`sec-install-pythonscripts` for where
they are and how they can be properly installed so that they are easy to
use.

The RADMC-3D package has two support-libraries:

#. ``python/tools/simpleread.py``

   The ``python/tools/simpleread.py`` is a set of functions to read the most
   important data files used by RADMC-3D.  However, the ``simpleread.py`` module
   is very simple, and does not read all RADMC-3D files in all formats. It can
   therefore only be used for certain (simple) models, and is primarily useful
   as a didactical tool.
  
#. ``python/radmc3dPy``

   The ``radmc3dPy`` package is a stand-alone Python package, written by Attila
   Juhasz, meant for the pre- and post-processing of RADMC-3D files.  It has its
   own manual, and has to be installed using e.g.~Python's ``pipinstall``
   method. This is described in the README file in that package.


.. _sec-simpleread-tools:

The simpleread.py library
=========================

For the most rudimentary analysis of the output (or input) files of RADMC-3D you
can use the ``simpleread.py`` file, which you can find in the ``python/tools/``
directory. If everything has been installed correctly, you should be able to
use it within Python like this::

  from radmc3d_tools.simpleread import *

Examples of data files you can read::

  d = read_dustdens()
  d = read_dusttemp()
  d = read_image()
  d = read_spectrum()
  d = read_dustkappa()
  d = read_gastemp()
  d = read_gasvelocity()
  d = read_molnumdens('co')
  d = read_mollevelpop('co')
  d = read_subbox(name='dust_temperature')
  d = read_subbox(name='dust_density')

Of course each one only if the corresponding file is present. Note that 'co' is
just an example molecule. In all these reading functions, except the ones for
images and spectra, the reading function automatically calls::

  grid = read_grid()

which reads the information about the spatial grid. This is then put inside the
``d`` object like this: ``d.grid``.

Here is an example of how you can plot the data (let us take the
``examples/run_simple_1/`` model, after we ran ``radmc3d mctherm``
and ``radmc3d image incl 60 phi 30 lambda 1000``)::

  import matplotlib.pyplot as plt
  from radmc3d_tools.simpleread import *
  import radmc3d_tools.natconst as nc
  tm = read_dusttemp()
  plt.figure()
  plt.plot(tm.grid.x/nc.au,tm.dusttemp[:,16,16])
  plt.xlabel('x [au]')
  plt.ylabel('T [K]')
  im = read_image()
  plt.figure()
  plt.imshow(im.image[:,:,0],vmax=3e-14)
  plt.show()

*Important:* These reading functions are rather basic. At the moment, no binary
file support is included (though this may change), no AMR octree grids can be
read, and several other limitations. For more sophisticated Python tools,
use the radmc3dPy library.

*Note:* For the ``read_subbox()`` function, you need to read the section
on the creation of regular-gridded datacubes of your 3D model, which is
Section :ref:`sec-subbox`. 


The radmc3dPy library
=====================

The ``radmc3dPy`` library is a sophisticated Python library that you
can use for the in-depth analysis of the output (or input) files of RADMC-3D.
It supports most in/output formats of RADMC-3D, including octree grids,
binary file formats etc.

The package is stand-alone, and has its own bitbucket repository:

https://bitbucket.org/at_juhasz/radmc3dpy/

But you can find a copy of this package also inside the RADMC-3D package,
in the directory ``python/radmc3dPy/``\ .

The ``radmc3dPy`` package has its own manual, so we will not reiterate it
here. Instead, please simply open the html manual in that package with a
browser. The entry file of that manual is the ``doc/html/index.html``\ .  On a
Mac you can simply type ``opendoc/html/index.html`` on the command line when you
are in the ``radmc3dPy`` directory. To install ``radmc3dPy`` please consult the
``README`` file in the ``radmc3dPy`` directory.

Once it is installed, you can use ``radmc3dPy`` in Python in the following
way:

#. Make sure to start Python 3 using {\small ipython --matplotlib} if you start
   Python from the command line. If you instead use a Jupyter notebook, make
   sure that as a first line you use ``%matplotlib inline`` to get the plots
   inside the notebook. These are standard Python things, so if you have
   trouble, ask your python friends or system manager.

#. Once you are inside Python you can include ``radmc3dPy`` using a simple
   ``from radmc3dPy import *``\ . This loads a series of radmc3dPy sub-libraries,
   including ``analyze``\ , ``image`` and several others.


We give here a very concise overview of the ``radmc3dPy`` package.
Please refer to the above mentioned stand-alone documentation for more details.


Model creation from within radmc3dPy
====================================

Several of the example models of the RADMC-3D ``examples/`` directory have been
implemented as part of the ``radmc3dPy`` package. This allows you to launch
these models straight from within ``radmc3dPy``\ . But this is merely
optional. You can equally well use the models in the ``examples/`` directory in
the RADMC-3D package, and post-process the results with ``radmc3dPy``\ .

To use one of the ``radmc3dPy``\ -internal models, create a directory
(e.g. ``mymodel``\ ), go into it, and go into iPython. Then type
``from radmc3dPy import *``\ . By typing ``models.getModelNames()`` you get a list
of available models. Suppose we choose the model 'ppdisk', then we would go
about like this (for example): ::

  from radmc3dPy import *
  analyze.writeDefaultParfile('ppdisk')
  setup.problemSetupDust('ppdisk', mdisk='1e-5*ms', gap_rin='[10.0*au]', gap_rout='[40.*au]', gap_drfact='[1e-5]', nz='0')

This example will set up a protoplanetary disk model in 2-D :math:`(r,\theta)`,
with a gap between 10 and 40 au. You can now run RADMC-3D to compute the dust
temperature structure, by calling (on the Linux shell): ::

  radmc3d mctherm

An image can be created with (again on the Lunix shell):
::

  radmc3d image lambda 1000 incl 60

And the image can be displayed (in Python) by
::

  import matplotlib.pyplot as plt
  from matplotlib import cm
  from radmc3dPy import *
  im=image.readImage()
  image.plotImage(im,vmax=3e-3,au=True,cmap=cm.gist_heat)


Diagnostic tools in radmc3dPy
=============================

No matter whether you use the ``radmc3dPy``\ -internal model set, or you create
your own model setup, you can use the extensive tool set inside ``radmc3dPy`` to
analyze the model itself, and the results of RADMC-3D calculations. In
everything below, we assume that you use ``from radmc3dPy import *`` beforehand.


Read the ``amr_grid.inp`` file
------------------------------

Use ``grid=analyze.readGrid()`` to read the information about the
spatial and wavelength grid. 


Read all the spatial data
-------------------------

Using ``data=analyze.readData()`` you read the entire spatial
structure of the model: The dust density, dust temperature, velocity
etc.


Read the ``image.out`` file
---------------------------

Using ``im=image.readImage()`` you read the ``image.out``
file created by RADMC-3D (if you call radmc3d for creating an image).
You can use the ``image.plotImage()`` function to display
the image with the proper axes and color bar.


Read the ``spectrum.out`` file
------------------------------

Any spectrum you create (a file called ``spectrum.out`` can be
read using ``s=analyze.readSpectrum()``\ .
