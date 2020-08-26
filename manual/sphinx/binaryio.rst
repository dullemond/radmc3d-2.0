.. _chap-binary-io:

Binary I/O files
****************


.. _sec-unformatted-overview:

Overview
========


By default all input and output files of RADMC-3D are in ASCII (i.e.\ text)
form. This makes it easier to verify if the files are ok. Also, it is easier
to produce files with the right format and read the output of RADMC-3D. The
disadvantage is that ASCII files are substantially larger than strictly
required to store their information content. For large models, i.e.\ models
with many grid points, this may lead to unpractically large files.

RADMC-3D supports a more compact data format: binary data. In this form, a
double precision variable occupies just 8 bytes, while a single precision
variable occupies just 4 bytes.

Unfortunately, Fortran-90 and Fortran-95 did, for a long time, not support true
binary files. Instead they offered 'f77-unformatted' files, which uses
'records', and is harder to read than true binary files. Recently, however, many
Fortran-90 and Fortran-95 compilers have introduced a true binary format, which
is called 'streaming access'. It is, actually, a Fortran-2003 feature, but has
been retroactively implemented into Fortran-90 and Fortran-95. The gfortran and
g95 compilers have it. Also the ifort compiler has it. Presumably others as
well.

RADMC-3D offers a binary I/O capability. A file containing three double
precision variables will have a length of exactly 24 bytes. Files with this
format will have extensions such as ``.binp``\ , ``.bdat`` or ``.bout``\ .

Here is a (presumably incomplete) list of files that have binary versions:

+------------------------+-----------+------------+
| Name                   | ascii     | binary     |
+------------------------+-----------+------------+
| ``dust_density``       |  ``.inp`` |  ``.binp`` |
+------------------------+-----------+------------+
| ``dust_temperature``   |  ``.inp`` |  ``.binp`` |
+------------------------+-----------+------------+
| ``dust_temperature``   |  ``.dat`` |  ``.bdat`` |
+------------------------+-----------+------------+
| ``gas_density``        |  ``.inp`` |  ``.binp`` |
+------------------------+-----------+------------+
| ``gas_temperature``    |  ``.inp`` |  ``.binp`` |
+------------------------+-----------+------------+
| ``electron_numdens``   |  ``.inp`` |  ``.binp`` |
+------------------------+-----------+------------+
| ``ion_numdens``        |  ``.inp`` |  ``.binp`` |
+------------------------+-----------+------------+
| ``levelpop_***``       |  ``.dat`` |  ``.bdat`` |
+------------------------+-----------+------------+
| ``numberdens_***``     |  ``.inp`` |  ``.binp`` |
+------------------------+-----------+------------+
| ``gas_velocity``       |  ``.inp`` |  ``.binp`` |
+------------------------+-----------+------------+
| ``microturbulence``    |  ``.inp`` |  ``.binp`` |
+------------------------+-----------+------------+
| ``stellarsrc_density`` |  ``.inp`` |  ``.binp`` |
+------------------------+-----------+------------+
| ``mean_intensity``     |  ``.out`` |  ``.bout`` |
+------------------------+-----------+------------+
| ``heatsource``         |  ``.inp`` |  ``.binp`` |
+------------------------+-----------+------------+



.. _sec-switch-to-binary:

How to switch to binary (or back to ascii)
==========================================

Specifying whether RADMC-3D should use ASCII or binary *input* is easy: It will
simply look which extension each input file has, and read it accordingly. If you
present RADMC-3D file input files with extension ``.binp``\ , it will read these
files as binaries.

More tricky is how to tell RADMC-3D to use binary files on *output*. By default,
RADMC-3D will always write ASCII style (``.out`` and ``.dat``\ ). However, if
you add the following line to the ``radmc3d.inp`` file: ::

  rto_style = 3

it will instead use binary output (``.bout`` and ``.bdat``\ ). And, for
completeness (though it is the default anyway), if you set ``rto_style=1``
RADMC-3D will write output in ASCII form. Note that ``rto_style = 2`` is
the old Fortran unformatted data format, which is deprecated.

For the binary form of output you can also tell RADMC-3D to use single-precision
for the main data, to produce smaller output files. This is done by adding the
following line to the ``radmc3d.inp`` file: ::

  rto_single = 1

By default RADMC-3D will always output double precision in the binary format.

*Note:* Images are still outputted in ascii even if you have ``rto_style=3``\
. This is because images are rarely files of huge size, and ascii files are
easier to analyze and check. However, sometimes images can be still quite big
(e.g.\ if you make multi-frequency images). Then it might still be useful to
output binary. If you want to also have the images in binary format, you must
set ::

  writeimage_unformatted = 1

in the ``radmc3d.inp`` file, *or* you add a keyword *imageunform*.


.. _sec-binary-io:

Binary I/O file format of RADMC-3D
==================================

The general format of the files listed in Section
:ref:`sec-unformatted-overview` is similar to the ASCII versions, just binary
this time. There is *one* additional number in the binary version: Right after
the format number comes an integer that gives the precision of the main
data. This number is either 4, meaning that the main data consists of 4-byte
floating point numbers (i.e.\ single precision), or 8, meaning that the main
data consists of 8-byte floating point numbers (i.e.\ double precision). Other
than that additional number, the order of the data is the same.

The following rules apply:

* With the exception of the ``amr_grid.binp`` file (see below),
  all integers are 8-byte integers. 

* Floating point numbers for the main data (i.e.\ the data that
  represents the space-dependent variables) are either 4-byte (single) or
  8-byte (double) precision numbers. Which of the two is specified in the
  second integer of the file (the integer right after the format number,
  see above).

* All other floating point numbers are double precision (i.e.\ 8-byte
  floats).

* For AMR-grids the ``amr_grid.binp`` file contains a huge list
  of 0 or 1 numbers (see Section :ref:`sec-amr-grid-oct-tree`). Since it is
  silly to use 8-byte integers for numbers that are either 0 or 1, the 
  numbers in this list are 1-byte integers (bytes). 


Example: According to Section :ref:`sec-dustdens` the ASCII file
``dust_density.inp`` file has the following format: ::

  iformat                                  <=== Typically 1 at present
  nrcells
  nrspec
  density[1,ispec=1]
  ..
  density[nrcells,ispec=1]
  density[1,ispec=2]
  ..
  ..
  ..
  density[nrcells,ispec=nrspec]

According to the above listed rules the binary file ``dust_density.binp`` file
then has the following format: ::

  <int8:iformat=1>
  <int8:precis=8>
  <int8:nrcells>
  <int8:nrspec>
  <dbl8:density[1,ispec=1]>
  ..
  <dbl8:density[nrcells,ispec=1]>
  <dbl8:density[1,ispec=2]>
  ..
  ..
  ..
  <dbl8:density[nrcells,ispec=nrspec]>

where the ``<int8:precis=8>`` means that this is an 8-byte integer that we call
'precis' (the name is irrelevant here), and it has value 8, and
``<dbl8:density[1,ispec=1]>`` means that this is a double-precision number
(8-byte float). In other words: the first 8 bytes of the file contain the format
number (which is 1 at present). The second 8 bytes contain the number 8, telling
that the main data (i.e.\ the ``density`` data) are double precision
variables. The third set of 8 bytes gives the number of cells, while the fourth
set gives the number of dust species. The data of ``density`` starts as of the
33rd byte of the file. If you want to compress the file even further, and you
are satisfied with single-precision data, then the file would look like: ::

  <int8:iformat=1>
  <int8:precis=4>
  <int8:nrcells>
  <int8:nrspec>
  <flt4:density[1,ispec=1]>
  ..
  <flt4:density[nrcells,ispec=1]>
  <flt4:density[1,ispec=2]>
  ..
  ..
  ..
  <flt4:density[nrcells,ispec=nrspec]>


Another example: According to Section :ref:`sec-dust-monochromatic-monte-carlo`
RADMC-3D can compute the mean intensity of radiation at each grid point at a set
of pre-defined frequencies, and write this out to an ASCII file called
``mean_intensity.out``\ . The contents of this file are: ::

  iformat                                  <=== Typically 2 at present
  nrcells
  nfreq                                    <=== Nr of frequencies 
  freq_1 freq_2 ... freq_nfreq             <=== List of frequencies in Hz
  meanint[1,icell=1]
  meanint[1,icell=2]
  ...
  meanint[1,icell=nrcells]
  meanint[2,icell=1]
  meanint[2,icell=2]
  ...
  meanint[2,icell=nrcells]
  ...
  ...
  ...
  meanint[nfreq,icell=1]
  meanint[nfreq,icell=2]
  ...
  meanint[nfreq,icell=nrcells]

By setting ``rto_style=3`` in the ``radmc3d.inp`` file, however, RADMC-3D will
instead produce a binary file called ``mean_intensity.bout``\ , which has the
contents: ::

  <int8:iformat=2>
  <int8:precis=8>
  <int8:nrcells>
  <int8:nfreq>
  <dbl8:freq_1>
  <dbl8:freq_2>
  ... 
  <dbl8:freq_nfreq>
  <dbl8:meanint[1,icell=1]>
  <dbl8:meanint[1,icell=2]>
  ...
  <dbl8:meanint[1,icell=nrcells]>
  <dbl8:meanint[2,icell=1]>
  <dbl8:meanint[2,icell=2]>
  ...
  <dbl8:meanint[2,icell=nrcells]>
  ...
  ...
  ...
  <dbl8:meanint[nfreq,icell=1]>
  <dbl8:meanint[nfreq,icell=2]>
  ...
  <dbl8:meanint[nfreq,icell=nrcells]>

If you also set ``rto_single=1`` in the ``radmc3d.inp`` file, then you
will get:
::

  <int8:iformat=2>
  <int8:precis=4>
  <int8:nrcells>
  <int8:nfreq>
  <dbl8:freq_1>
  <dbl8:freq_2>
  ... 
  <dbl8:freq_nfreq>
  <flt4:meanint[1,icell=1]>
  <flt4:meanint[1,icell=2]>
  ...
  <flt4:meanint[1,icell=nrcells]>
  <flt4:meanint[2,icell=1]>
  <flt4:meanint[2,icell=2]>
  ...
  <flt4:meanint[2,icell=nrcells]>
  ...
  ...
  ...
  <flt4:meanint[nfreq,icell=1]>
  <flt4:meanint[nfreq,icell=2]>
  ...
  <flt4:meanint[nfreq,icell=nrcells]>

Note that only the mean intensity data (the main data) are single precision
floats.

