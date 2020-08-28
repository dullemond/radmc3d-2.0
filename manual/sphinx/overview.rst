Overview of the RADMC-3D package
********************************

Introduction
============

The RADMC-3D code is written in fortran-90 and should compile with most f90
compilers without problems. It needs to be compiled only once for each
platform.

The executable is called ``radmc3d`` and it performs all the model
calculations of the RADMC-3D package, for instance the Monte Carlo simulations,
ray-tracing runs (images, spectra), etc. There is also a set of useful
subroutines written in the Python language to use the ``radmc3d`` code,
but ``radmc3d`` can also run without it. In that case the user will
have to write his/her own pre- and post-processing subroutines.

.. _sec-requirements:

Requirements
============

The following pre-installed software is required:

#. *Operating system: Unix-like (e.g. Linux or MacOSX)*

   This package runs under Unix-like environment (e.g. Linux or MacOSX), but has
   not been tested under Windows. There is no particular reason why it should not
   also run under Windows, but it would require different ways of file
   handling. *In this manual we always assume a Unix-like environment,
   in which we will make use of a bash command-line interface (CLI).* We will
   call this the *shell*. 

#. ``make`` or ``gmake``
   
   This is the standard tool for compiling packages on all Unix/Linux/MacOS-based
   systems.
  
#. ``perl``
   
   This is a standard scripting language available on most or all
   Unix/Linux-based systems. If you are in doubt: type ``which perl``
   to find the location of the ``perl`` executable. See http://www.perl.org/
   for details on perl, should you have any
   problems. But on current-day Unix-type operating systems perl is nearly
   always installed in the ``/usr/bin/`` directory. If you do not
   have Perl installed, you can also do without. Its sole use is to copy
   the executables into your home ``$HOME/bin/`` directory for quick
   access from the Unix/Linux/MacOS command line. You can work around that,
   if necessary.

#. *A fortran-90 compiler*
   
   Preferably the ``gfortran`` compiler (which the current
   installation assumes is present on the system). Website: 
   http://gcc.gnu.org/fortran/. Other compilers may work, but have not
   been tested yet.

#. *An OpenMP-fortran-90 compiler (optional)*
   
   Only needed if you want to use the parallelized OpenMP version for the thermal
   Monte Carlo (for faster execution). Preferably the ``GNUOpenMP/GOMP``
   compiler which is an implementation of OpenMP for the Fortran compiler in the
   GNU Compiler Collection. Websites: http://openmp.org/wp and
   http://gcc.gnu.org. Other compilers may work, but we give no
   guarantee.
   
#. *Python version 3 with standard libraries*
   
   The core RADMC-3D code ``radmc3d`` (written in Fortran-90) is the
   raw workhorse code that reads some input files and produces some output
   files. Typically you will not need to worry about the internal workings of the
   RADMC-3D code. All you need to do is produce the proper input files, run
   ``radmc3d``, and read the output files for post-processing (such as
   displaying and analyzing the results). This pre- and post-processing is
   done in Python. This RADMC-3D distribution provides you with the Python tools
   you need, though you will likely want to program your own additional Python code
   to adjust the models to your own needs. To use the Python tools provided in
   this RADMC-3D distribution, you need Python version 3
   (though most things should also work with the depricated Python 2), with a
   set of standard libraries such as ``numpy`` and ``matplotlib``.
   Typically we will assume that Python is used as a Python or iPython command-line
   interface, which we shall call the *Python command line* (as opposed to the
   *shell*). The user can, of course, also use Jupyter Notebooks instead.
   But for the sake of clarity, in this manual we assume the use of Python the
   Python command line.
   

Note that the Monte Carlo code RADMC-3D itself (``radmc3d``) is in Fortran-90. Only the
creation of the input files (and hence the problem definition) and the analysis
of the output files is done in Python. The user is of course welcome to use
other ways to create the input files for RADMC-3D if he/she is not able or
willing to use Python for whatever reason. Therefore Python is not strictly
required for the use of this code. However, all examples and support
infrastructure is provided in Python.


Contents of the RADMC-3D package
================================

RADMC-3D package as a .zip archive
----------------------------------
If you obtain RADMC-3D from its website, it will be
packed in a zip archive called
``radmc-3d_v*.*_dd.mm.yy.zip`` where the ``*.*`` is the version
number and ``dd.mm.yy`` is the date of this version.
To unpack on a linux, unix or Mac OS X machine you type::

  unzip <this archive file>

i.e. for example for ``radmc-3d_v2.0_25.08.20.zip`` you type::

  unzip radmc-3d_v2.0_25.08.20.zip

RADMC-3D package from the github repository
-------------------------------------------
If you obtain RADMC-3D by cloning its github repository, you will get
a copy of the full git repository of RADMC-3D. In principle this is
not much different from unzipping the .zip archive. But it is more
powerful: You can more easily stay up to date with the latest bugfixes,
and you can see the entire development history of this version of the
code. See https://git-scm.com/book/en/v2 for an extensive documentation
of how to use git.

The way to produce a clone of RADMC-3D in the directory where your
shell is currently is, is like this::

  git clone https://github.com/dullemond/radmc3d-2.0.git

This will create the directory ``radmc3d-2.0/``. At any time you can
pull the latest version from the repository like this::

  cd radmc3d-2.0/
  git pull

Keep in mind, however, that while the repository is always the very
latest version, this comes with a (small) risk that some new features
may not have been tested well, or (new) bugs may have been introduced.
Overall, however, we advise to use the github repository instead of
the .zip archive from the website.

Contents of the package
-----------------------

The RADMC-3D package has the following subdirectory
structure::
  
  src/
  python/
  examples/
     run_simple_1/
     run_simple_1_userdef/
     run_simple_1_userdef_refined/
     .
     .
     .
  opac/
  manual/

plus some further directories.

The first directory, ``src/``, contains the fortran-90 source code for
RADMC-3D. The second directory, ``python/``, contains two sets of Python modules
that are useful for model preparation and post-processing. One is a directory
called ``radmc3d_tools/``, which contains some simple Python tools that might be
useful. The other is a directory called ``radmc3dPy/``, which is a high-level
stand-alone Python library developed by Attila Juhasz for RADMC-3D. The third
directory contains a series of example models. The fourth directory,
``opac/`` contains a series of tools and data for creating the opacity
files needed by RADMC-3D (though the example models all have their own
opacity data already included), The fifth directory contains
this manual.

Units: RADMC-3D uses CGS units
==============================

The RADMC-3D package is written such that all units are in CGS (length in
cm, time in sec, frequency in Hz, energy in erg, angle in steradian). There
are exceptions:

* Wavelength is usually written in micron
* Sometimes angles are in degrees (internally in radian, but input as degrees)

