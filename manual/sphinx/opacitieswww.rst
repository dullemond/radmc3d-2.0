.. _chap-acquiring-opacities:

Acquiring opacities from the WWW
********************************

Opacities are the basic ingredients necessary for any model with
RADMC-3D. The example models in this package contain example opacities, but
for professional usage of RADMC-3D it may be necessary to get specific
opacity data from the web. These opacity data are usually in a wide variety
of formats. To enable RADMC-3D to read them usually requires a conversion
into RADMC-3D-readable form (see Section :ref:`sec-opacities` for dust
opacities and Section :ref:`sec-molecule-xxx-inp` for gas line opacities).

To make it easier for the user to create RADMC-3D-readable input files
from opacity data downloaded from the web, we now feature a new directory
``opac/`` in the RADMC-3D distribution in which, for several of
the most common WWW databases, we provide Python routines for the conversion.
Please read the ``README_*`` files in this directory and its
subdirectories for details.

Note also that Carsten Dominik made a very nice and easy-to-use tool to
generate dust opacities on the Linux/Mac command line. It is called
``optool``, and can be found on github at
https://github.com/cdominik/optool .
It can produce RADMC-3D-ready dust opacity files.

