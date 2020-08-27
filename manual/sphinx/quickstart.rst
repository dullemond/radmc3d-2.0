.. _chap-quick-start:

Quickstarting with RADMC-3D
***************************

In general I recommend reading the manual fully, but it is often useful
to get a quick impression of the package with a quick-start. To make your
first example model, this is what you do:

#. When you read this you have probably already unzipped this package, or
   cloned the git repository.
   You should find, among others, a ``src/`` directory and a
   ``examples/`` directory. Go into the ``src/`` directory.
#. Edit the ``src/Makefile`` file, and make sure to set the
   ``FF`` variable to the Fortran-90 compiler you have installed on
   your system.
#. Type ``make``. If all goes well, this should compile the
   entire code and create an executable called ``radmc3d``.
#. Type ``make install``. If all goes well this should try to
   create a link to ``radmc3d`` in your ``$HOME/bin/``
   directory, where ``$HOME`` is your home directory.
   If this ``$HOME/bin/`` directory does not exist, it will ask to make one.
#. Make sure to have the ``$HOME/bin/`` directory in your path.  If
   you use, for instance, the ``bash`` shell, you do this by setting the
   ``PATH`` variable by adding a line like ``export
   PATH=$HOME/bin:\$PATH`` to your ``$HOME/.bashrc`` file. If you
   change these things you may have to open a new shell to make sure that the
   shell now recognizes the new path.
#. Check if the executable is OK by typing ``radmc3d`` in the
   shell. You should get a small welcoming message by the code.
#. Now enter the directory ``examples/run_simple_1/``. This is
   the simplest example model.
#. Type ``python problem_setup.py`` (Note: you must have a
   working Python distribution on your computer, which is reasonably
   up to date, with ``numpy`` and ``matplotlib`` libraries
   included). This will create a series of input files for RADMC-3D.
#. Type ``radmc3d mctherm``. This should let the code do a Monte
   Carlo run. You should see ``Photon nr 1000``, followed by
   ``Photon nr 2000``, etc until you reach ``Photon nr
   1000000``. The Monte Carlo modeling for the dust temperatures has now
   been done. A file ``dust_temperature.dat`` should have
   been created.
#. Type ``radmc3d image lambda 1000 incl 60 phi 30``. This should
   create an image with the camera at inclination 60 degrees (from pole-on), and
   rotated 30 degrees (along the polar axis, clockwise, i.e.\ the object rotating
   counter-clockwise), at wavelength :math:`\lambda=1000\,\mu\mathrm{m}` (i.e. at 1
   millimeter wavelength). The file that contains the image is ``image.out``.
   It is a text file that can be read with the ``simpleread.py`` tool in the
   directory ``python/radmc3d_tools/``.

If you experience troubles with the above steps, and you cannot fix it,
please read the next chapters for more details. 
