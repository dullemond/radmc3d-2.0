.. _chap-compilation:

Installation of RADMC-3D
************************

Although the RADMC-3D package contains a lot of different software,
the main code is located in the ``src/`` directory, and is
written in Fortran-90. The executable is ``radmc3d``. Here
we explain how to compile the fortran-90 source codes and create 
the executable ``radmc3d``.

.. _sec-makeing:

Compiling the code with 'make'
==============================

To compile the code, enter the ``src/`` directory in your shell. You now
*may* need to edit the ``Makefile`` in this directory using your favorite
text editor and replace the line ::

  FF = gfortran -fopenmp

with a line specifying your own compiler (and possibly OpenMP directive,
if available). If, of course, you use gfortran,
you can keep this line. But if you use, e.g., ifort, then replace the above
line by ::

  FF = ifort -openmp

(note the slightly different OpenMP directive here, too). 
If you save this file, and you are back in the shell, you can compile the
radmc3d code by typing ::

  make

in the shell. If all goes well, you have now created a file called ``radmc3d``
in the ``src/`` directory.

If, for whatever reason, the OpenMP compilation does not work, you can also
compile the code in serial mode. Simply remove the ``-fopenmp`` directive.


The install.perl script
=======================

If instead of typing just ``make`` you type ::

  make install

(or you first type ``make`` and then ``make install``, it is the same), then in
addition to creating the executable, it also automatically executes a perl
script called ``install.perl`` (located also in the ``src/`` directory).  This
PERL script installs the code in such a way that it can be conveniently used in
any directory. What it does is:

* It checks if a ``bin/`` directory is present in your home
  directory (i.e. a ``$HOME/bin/`` directory). If not, it asks if
  you want it to automatically make one.
* It checks if the ``$HOME/bin/`` directory is in the *path* of
  the currently used shell. This is important to allow the computer to look
  for the program ``radmc3d`` in the ``$HOME/bin/`` directory. If you
  use a bash shell, then you can add the following line to your
  ``$HOME/.bashrc``::
    
    export PATH=/myhomedirectory/bin/python:$PATH
    
* It creates a file ``radmc3d`` in this ``$HOME/bin/``
  directory with the correct executable permissions. This file is merely a
  dummy executable, that simply redirects everything to the true ``radmc3d``
  executable located in your current ``src/``
  directory. When you now open a new shell, the path contains the
  ``$HOME/bin/`` directory, and the command ``radmc3d`` is
  recognized. You can also type ``source $HOME/.bashrc`` followed
  by ``rehash``. This also makes sure that your shell recognizes the
  ``radmc3d`` command.
* It checks if a ``python/`` subdirectory exists in the above
  mentioned ``bin/`` directory, i.e.\ a ``$HOME/bin/python/``
  directory. If not, it asks if you want it to automatically create one.
* If yes, then it will copy all the files ending with ``.py`` in
  the ``python/radmc3d_tools/`` directory of the distribution to that
  ``$HOME/bin/python/radmc3d_tools/`` directory. This is useful to
  allow you to make an ``PYTHONPATH`` entry to allow python to find
  these python scripts automatically.

Note that this perl script installs the code only for the user that installs
it. A system-wide installation is not useful, because the code package is not
very big and it should remain in the control of the user which version of the
code he/she uses for each particular problem.

If all went well, then the ``perl.install`` script described here is
called automatically once you type ``make install`` following the
procedure in Section :ref:`sec-makeing`.

Before the installation is recognized by your shell, you must now either
type ``rehash`` in the shell or simply open a new shell. 

How do you know that all went OK? If you type ``radmc3d`` in the
shell the RADMC-3D code should now be executed and give some comments. It
should write::
  
  ================================================================
       WELCOME TO RADMC-3D: A 3-D CONTINUUM AND LINE RT SOLVER    
                                                                  
                           VERSION 2.0                            
                                                                  
                 (c) 2008-2020 Cornelis Dullemond                 
                                                                  
        Please feel free to ask questions. Also please report     
         bugs and/or suspicious behavior without hestitation.     
       The reliability of this code depends on your vigilance!    
                     dullemond@uni-heidelberg.de                  
                                                                  
    To keep up-to-date with bug-alarms and bugfixes, register to  
                      the RADMC-3D forum:                         
             http://radmc3d.ita.uni-heidelberg.de/phpbb/          
                                                                  
               Please visit the RADMC-3D home page at             
   http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/ 
  ================================================================
   
  Nothing to do... Use command line options to generate action:
    mctherm        : Do Monte Carlo simul of thermal radiation
    mcmono         : Do Monte Carlo simul only for computing mean intensity
    spectrum       : Make continuum spectrum
    image          : Make continuum image

on the screen (or for newer versions of RADMC-3D perhaps some more
or different text). This should also work from any other directory.


What to do if this all does not work?
=====================================

In case the above compilation and installation does not work, here is a 
proposed procedure to do problem hunting:

#. First, answer the following questions:

   * Did you type ``make install`` in the ``src/``
     directory? I mean, did you not forget the ``install`` part?
   * Did you put ``$HOME/bin/`` in your path (see above)?
   * If you just added ``$HOME/bin/`` to your path, did you
     follow the rest of the procedure (either closing the current shell and
     opening a new shell or typing the ``source`` and ``rehash`` commands as
     described above)?

   If this does not help, then continue:

#. Close the shell, open a new shell.
#. Go to the RADMC-3D ``src/`` directory.
#. Type ``./radmc3d``. This should give the above message. If
   not, then make sure that the compilation went right in the first place:
#. Type ``rm -f radmc3d``, to make sure that any old executable
   is not still present.
#. Type ``make clean``. This should return the sentence
   ``OBJECT and MODULE files removed.``
#. In case the problem lies with the OpenMP parallellization, you
   could do ``cp Makefile_normal Makefile``, which switches
   off the OpenMP compilation.
#. Then type ``make``. This should produce a set of lines, each
   representing a compilation of a module, e.g. ``gfortran -c -O2
   ./amr_module.f90 -o amr_module.o``, etc. The final line should be
   something like ``gfortran -O2 main.o ..... gascontinuum_module.o -o radmc3d``. If instead there
   is an error message, then do the following:

    * Check if the compiler used (by default ``gfortran``) is
      available on your computer system.
    * If you use an other compiler, check if the compiler options used
      are recognized by your compiler.
    * Check if the executable ``radmc3d`` is now indeed present.
      If it is not present, then something must have gone wrong with the
      compilation. So then please check the compilation and linking stage
      again carefully. 

    If you followed all these procedures, but you still cannot get even the
    executable in the ``src/`` directory to run by typing (in the
    ``src/`` directory) ``./radmc3d`` (don't forget the dot
    slash!), then please contact the author.
#. At this point we assume that the previous point worked. Now go to
   another directory (any one), and type ``radmc3d``.  This should
   also give the above message. If not, but the ``radmc3d`` executable
   was present, then apparently the shell path settings are wrong. Do this:

   * Check if, in the current directory (which is now not ``src/``)
     there is by some accident another copy of the executable
     ``radmc3d``. If yes, please remove it. 
   * Type ``which radmc3d`` to find out if it is recognized at all,
     and if yes, to which location it points. 
   * Did you make sure that the shell path includes the ``$HOME/bin/``
     directory, as it should? Otherwise the shell does not know 
     where to find the ``$HOME/bin/radmc3d`` executable (which is
     a perl link to the ``src/radmc3d`` executable).
   * Does the file ``$HOME/bin/radmc3d`` perl file exist in the
     first place? If no, check why not. 
   * Type ``less $HOME/bin/radmc3d`` and you should
     see a text with first line being ``#!/usr/bin/perl`` and the
     second line being someting like 
     ``system("/Users/user1/radmc-3d/version_2.0/src/radmc3d @ARGV");``
     where the ``/Users/user1`` should of course be the path to
     your home directory, in fact to the directory in which you installed
     RADMC-3D.

If this all brings you no further, please first ask your system administrators
if they can help. If not, then please contact the author.

.. _sec-install-pythonscripts:

Installing the simple Python analysis tools
===========================================

RADMC-3D offers (in addition to the model setup scripts in the ``examples/``
subdirectories) two Python support libraries:

#. ``python/radmc3d_tools/``
   
   This library contains only some bare-bones small Python scripts.
  
#. ``python/radmc3dPy/``

   This library is a sophisticated stand-alone library developed by
   Attila Juhasz, and further maintained together with the RADMC-3D
   main author.

How to install and use the ``python/radmc3d_tools/``
----------------------------------------------------

The installation of the ``python/radmc3d_tools`` should be automatic when you
type ``make install`` in the ``src/`` code directory (see
above). It will copy the files to the ``bin/python/radmc3d_tools/``
directory in your home directory. If this directory does not exist, you
will be asked if you want it to be created. If you confirm (typing 'y'),
then the files from the ``python/radmc3d_tools/`` directory will be
copied into the ``$HOME/bin/python/radmc3d_tools/`` directory.

Now you need to make sure that Python knows that these tools are there.
In Python here are two ways how you can make sure that Python automatically
finds these scripts:

#. Under Unix/Linux/MacOSX you can set the ``PYTHONPATH`` directly in your
   ``.bashrc`` file. For example: in
   ``.bashrc`` (if you use the bash shell) you can write::

     export PYTHONPATH=$HOME/bin/python:$PYTHONPATH

(where ``$HOME`` is your home directory name). 

#. Alternatively you can set the ``PYTHONPATH`` directly from within
   Python with the python command::

     import os
     import sys
     home = os.environ["HOME"]
     sys.path.append(home+'/bin/python')

If all goes well, if you now start Python you should be able to have access to
the basic Python tools of RADMC-3D directly. To test this, try typing ``from
radmc3d_tools.simpleread import *`` in Python. If this gives an error message
that ``simpleread.py`` cannot be found, then please ask your system
administrators how to solve this.

You may ask why first copy these files to ``$HOME/bin/python/radmc3d_tools/``
and not point PYTHONPATH directly to the ``python/radmc3d_tools`` in your RADMC-3D
distribution? The reason is that if you have multiple versions of RADMC-3D on
your computer system, you always are assured that Python finds the python
routines belonging to the latest installation of RADMC-3D (note: only assured if
that latest compilation was done with ``make install``).

Now you should be ready to use the tools. The most important one would be
the ``simpleread.py`` tool, which contains a set of functions for
reading typical RADMC-3D input and output files (though only for regular
model grid, not for octree grids). In a Python command line interface
you can import them by::

  from radmc3d_tools import simpleread

And you can then, for instance, read the dust density file with::

  d = simpleread.read_dustdens()

Here, ``d`` is now an object containing a ``d.grid`` subobject (which contain
information about the grid) and the dust density array ``d.rhodust``. Have a
look at the various functions in ``simpleread``, to see what is available.

How to install and use the ``python/radmc3dPy`` library
----------------------------------------------------------

The installation of the ``python/radmc3dPy`` package is described in the
``python/radmc3dPy/README`` file. In short, by going into the
``python/radmc3dPy/`` directory and typing in the shell::

  python setup.py install --user

it should install itself right into your Python distribution. For instance,
if you have ``anaconda3`` on a Mac, it would copy the files into the
directory ::

  $HOME/.local/lib/python3.7/site-packages/radmc3dPy/

Python knows where to find it there.

Now you should be ready to use ``radmc3dPy``, by importing it::

  import radmc3dPy

``radmc3dPy`` consists of
several sub libraries such as ``radmc3dPy.analyze`` and
``radmc3dPy.image``. For instance, to read the dust density
distribution, you could do this::

  from radmc3dPy import analyze
  d = analyze.readData(ddens=True)

The ``d.rhodust`` array now contains the dust density.

For more information, please consult the ``radmc3dPy`` documentation
in the ``python/radmc3dPy/doc/`` directory.


.. _sec-special-purpose-compile:

Making special-purpose modified versions of RADMC-3D (optional)
===============================================================

For most purposes it should be fine to simply compile the latest version of
RADMC-3D once-and-for-all, and simply use the resulting ``radmc3d``
executable for all models you make. Normally there is no reason to have to
modify the code, because models can be defined quite flexibly by preparing
the various input files for RADMC-3D to your needs. So if you are an 
average user, you can skip to the next subsection without problem.

But sometimes there *is* a good reason to want to modify the code.  For
instance to allow special behavior for a particular model. Or for a model
setup that is simply easier made internally in the code rather than by
preparing large input files. One can imagine some analytic model setup
that might be easier to create internally, so that one can make use of
the full AMR machinery to automatically refine the grid where needed.
Having to do so externally from the code would require you to set up
your own AMR machinery, which would be a waste of time. 

The problem is that if the user would modify the central code for each
special purpose, one would quickly lose track of which modification of the
code is installed right now. 

Here is how this problem is solved in RADMC-3D:

* For most purposes you can achieve your goals by only editing the file
  ``userdef_module.f90``. This is a set of standard subroutines
  that the main code calls at special points in the code, and the user can
  put anything he/she wants into those subroutines. See Chapter :ref:`chap-internal-setup`
  for more information about these standard
  subroutines. This method is the safest way to create special-purpose
  codes. It means (a) that you know that your modification cannot do much
  harm unless you make really big blunders, because these subroutines are
  meant to be modified, and (b) you have all your modifications *only*
  in one single file, leaving the rest of the code untouched.
* You can create a *local* version of the code, without touching
  the main code. Suppose you have a model directory ``run_mymodel`` and for
  this model you want to make a special-purpose version of the code.
  This is what you do:

  #. Copy the Makefile from the ``src/`` directory into ``run_mymodel``.
  #. Copy the ``.f90`` file(s) you want to modify from the ``src/``
     directory into ``run_mymodel``. Usually you only want to modify
     the ``userdef_module.f90`` file, but you can also copy any other 
     file if you want.
  #. In the ``run_mymodel/Makefile`` replace the ``SRC = .`` line with
     ``SRC = XXXXXX``, where ``XXXXXX`` should
     be the *full* path to the ``src/`` directory. An example line
     is given in the Makefile, but is commented out.
  #. In the ``run_mymodel/Makefile`` make sure that all the
     ``.f90`` files that should remain as they are have a ``$(SRC)/``
     in front of the name, and all the ``.f90`` files that
     you want to modify (and which now have a copy in the ``run_mymodel``
     directory) have a ``./`` in front of the name. By
     default all ``.f90`` files have ``$(SRC)/`` in front of
     the name, except the ``userdef_module.f90`` file, which has a
     ``./`` in front of the name because that is the file that is
     usually the one that is going to be edited by you.
  #. Now edit the local ``.f90`` files in the ``run_mymodel`` directory
     in the way you want. See Chapter :ref:`chap-internal-setup` for more details.
  #. Now *inside* the ``run_mymodel`` directory you can now type
     ``make`` and you will create your own local ``radmc3d`` executable.
     NOTE: Do not type ``make install`` in this case, because it should
     remain a local executable, only inside the ``run_mymodel`` directory.
  #. If you want (though this is not required) you can clean up all the
     local ``.o`` and ``.mod`` files by typing ``make
     clean``, so that your ``run_mymodel`` directory is not filled
     with junk.
  #. You can now use this special purpose version of ``radmc3d``
     by simply calling on the command line: ``./radmc3d``, with any
     command-line options you like. Just beware that, depending on the order
     in which you have your paths set (in tcsh or bash) typing just 
     ``radmc3d`` *may* instead use the global version (that you
     may have created in the ``src/`` directory with ``make
     install``). So to be sure to use the *local* version, just put the
     ``./`` in front of the ``radmc3d``.

Note: In chapter :ref:`chap-internal-setup` there is more information on how to
set up models internally in the code using the method described here.

Note: You can use ``make clean`` to remove all the .o and .mod files from your
model directory, because they can be annoying to have hanging around. By typing
``make cleanmodel`` you remove, in addition to the .o and .mod files, also all
model input and output files, with the exception of dust opacity or molecular
data files (because these latter files are usually not created locally by the
``problem_setup.py`` script). By typing ``make cleanall`` you remove everything
*except* the basic files such as the ``Makefile``, any ``.f90`` files, any
``.py`` files, the dust opacity or molecular data files and ``README`` files.
