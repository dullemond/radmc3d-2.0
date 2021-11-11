Introduction
************

What is RADMC-3D?
=================

RADMC-3D is a software package for astrophysical radiative transfer calculations
in arbitrary 1-D, 2-D or 3-D geometries. It is mainly written for continuum
radiative transfer in dusty media, but also includes modules for gas line
transfer. Typical applications would be protoplanetary disks, pre- and
proto-stellar molecular cloud cores, and similar objects. It does not treat
photoionization of gas, nor does it treat chemistry. It can self-consistently
compute dust temperatures for the radiative transfer, but it is not equipped for
self-consistent gas temperature computations (as this requires detailed coupling
to photochemistry). The main strength of RADMC-3D lies in the flexibility of the
spatial setup of the models: One can create or use parameterized dust and/or gas
density distributions, or one can import these from snapshots of hydrodynamic
simulations.


Capabilities
============

Here is a list of current and planned features. Those features that are now
already working are marked with [+], while those which are not yet (!!) built in
are marked with [-]. Those that are currently being developed are marked with
[.] and those that are ready, but are still in the testing phase are marked with
[t].

#. Coordinate systems:
   
  #. [+] Cartesian coordinates (3-D)
  #. [+] Spherical coordinates (1-D, 2-D and 3-D)

#. Gridding systems (regular and adaptive mesh refinement grids are
   available for cartesian *and* spherical coordinates):

  #. [+] Regular
  #. [+] Adaptive Mesh Refinement: oct-tree style
  #. [+] Adaptive Mesh Refinement: layered ('patch') style
  #. [-] Voronoi gridding *[To be implemented on request]*

#. Radiation mechanisms:

    #. [+] Dust continuum, thermal emission
    #. [+] Dust continuum scattering:

      #. [+] ...in isotropic approximation
      #. [+] ...with full anisotropy
      #. [+] ...with full Stokes and Polarization

    #. [-] Dust quantum heated grains *[To be implemented on request]*
    #. [t] Polarized dust emission by aligned grains *[first test version]*
    #. [+] Gas line transfer (LTE)
    #. [+] Gas line transfer (non-LTE: LVG)
    #. [+] Gas line transfer (non-LTE: LVG + Escape Probability)
    #. [-] Gas line transfer (non-LTE: full transfer)
    #. [+] Gas line transfer with user-defined populations
    #. [+] Gas continuum opacity and emissivity sources

#. Radiation netto sources for continuum:

    #. [+] Discrete stars positioned at will
    #. [t] Continuous 'starlike' source
    #. [t] Continuous 'dissipation' source
    #. [t] External 'interstellar radiation field'

#. Imaging options:

    #. [+] Observer from 'infinite' distance
    #. [+] Zoom-in at will
    #. [+] Flux-conserving imaging, i.e. pixels are recursively refined
    #. [+] A movie-making tool
    #. [+] Multiple wavelengths in a single image
    #. [+] Local observer with perspective view (for PR movies!)

#. Spectrum options:

    #. [+] SED spectrum (spectrum on 'standard' wavelength grid)
    #. [+] Spectrum on any user-specified wavelength grid
    #. [+] Spectrum of user-specified sub-region (pointing)
    #. [t] Specification of size and shape of a primary 'beam' for spectra

#. User flexibility:

    #. [+] Free model specification via tabulated input files
    #. [+] Easy special-purpose compilations of the code (optional)

#. Front-end Python packages:

   #. [+] Python simple tools for RADMC-3D
   #. [+] Python RADMC-3D library {\small\tt radmc3dPy} (author: A. Juhasz)

#.

    #. [+] Stars can be treated as point-sources or as spheres
    #. [+] Option to calculate the mean intensity :math:`J_\nu(\vec x)` in the model
    #. [+] OpenMP parallellization of the Monte Carlo


Version tracker
===============

The RADMC-3D software package in under continuous development. A very
detailed development log-book is found in the git repository.
A more user-friendly overview of the development history can be 
found in this manual, in appendix \ref{chap-development-history}.


Copyright
=========

RADMC-3D was developed from 2007 to 2010/2011 at the Max Planck Institute
for Astronomy in Heidelberg, funded by a Max Planck Research Group grant
from the Max Planck Society. As of 2011 the development continues at the
Institute for Theoretical Astrophysics (ITA) of the Zentrum für Astronomy
(ZAH) at the University of Heidelberg. 

**The use of this software is free of charge. However, it is not allowed
to distribute this package without prior consent of the lead author
(C.P. Dullemond). Please refer any interested user to the web site of this
software where the package is available, which is currently:**

http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d

or the github repository:

https://github.com/dullemond/radmc3d-2.0

The github repository will always have the latest version, but it may
not be always the most stable version (though usually it is). 

Contributing authors
====================

The main author of RADMC-3D is Cornelis P. Dullemond. However, the main
author of the ``radmc3dPy`` Python package is Attila Juhasz.

Numerous people have made contributions to RADMC-3D. Major contributions
are from:

* Michiel Min
* Attila Juhasz
* Adriana Pohl
* Rahul Shetty
* Farzin Sereshti
* Thomas Peters
* Benoit Commercon
* Alexandros Ziampras

The code profited from testing, feedback and bug reports from (incomplete list):
  
* Daniel Harsono
* Rainer Rolffs
* Laszlo Szucs
* Sean Andrews
* Stella Offner
* Chris Beaumont
* Katrin Rosenfeld
* Soren Frimann
* Jon Ramsey
* Seokho Lee
* Blake Hord
* Tilman Birnstiel
* Uma Gorti
  
and others.


Disclaimer
==========

**IMPORTANT NOTICE 1: I/We reject all responsibility for the use of this
package. The package is provided as-is, and we are not responsible for any
damage to hardware or software, nor for incorrect results that may result
from the software. The user is fully responsible for any results from this
code, and we strongly recommend thorough testing of the code before using
its results in any scientific papers.**

**IMPORTANT NOTICE 2: Any publications which involve the use of this
software must mention the name of this software package and cite the
accompanying paper once it is published (Dullemond et al.\ in prep), or
before that the above mentioned web site.**
