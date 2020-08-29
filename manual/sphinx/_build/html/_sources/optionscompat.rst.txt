.. _chap-table-possibilities:

Which options are mutually incompatible?
****************************************

For algorithmic reasons not all options / coordinate systems and all grids
are compatible with each other. Here is an overview of which options/methods
work when. Note that only options/methods for which this is a possible issue
are listed. 


Coordinate systems
==================

Some coordinate systems exclude certain possibilities. Here is a list.

+----------------------------------------------------------+---------+--------+------------------+--------+
| Option/Method:                                           | Cart 3D | Sph 3D | Sph 2D (axisymm) | Sph 1D |
+----------------------------------------------------------+---------+--------+------------------+--------+
| Second order ray-tracing                                 | yes     | yes    | yes              | yes    |
+----------------------------------------------------------+---------+--------+------------------+--------+
| Isotropic scattering                                     | yes     | yes    | yes              | yes    |
+----------------------------------------------------------+---------+--------+------------------+--------+
| An-isotropic scattering for thermal Monte Carlo          | yes     | yes    | yes              | yes    |
+----------------------------------------------------------+---------+--------+------------------+--------+
| An-isotropic scattering for monochromatic Monte Carlo    | yes     | yes    | yes              | yes    |
+----------------------------------------------------------+---------+--------+------------------+--------+
| An-isotropic scattering for images and spectra           | yes     | yes    | yes              | no     |
+----------------------------------------------------------+---------+--------+------------------+--------+
| Full Stokes scattering for thermal Monte Carlo           | yes     | yes    | yes              | yes    |
+----------------------------------------------------------+---------+--------+------------------+--------+
| Full Stokes scattering for monochromatic Monte Carlo     | yes     | yes    | yes              | yes    |
+----------------------------------------------------------+---------+--------+------------------+--------+
| Full Stokes scattering for images and spectra            | yes     | yes    | yes              | no     |
+----------------------------------------------------------+---------+--------+------------------+--------+
| Gas lines                                                | yes     | yes    | yes              | yes    |
+----------------------------------------------------------+---------+--------+------------------+--------+
| Gas lines and Doppler-shift line catching                | yes     | yes    | no               | no     |
+----------------------------------------------------------+---------+--------+------------------+--------+


Scattering off dust grains
==========================

The inclusion of the effect of scattering off dust grains in images and spectra
typically requires a separate Monte Carlo computation for each image. This is
done automatically by RADMC-3D. But it means that there are some technical
limitations.

+------------------------------------------------------+---------------+-------------------------+------------------------------------+
| Option/Method:                                       | No scattering | Isotropic approximation | Full anisotropic/Stokes scattering |
+------------------------------------------------------+---------------+-------------------------+------------------------------------+
| Fast multi-frequency ray tracing for spectra (auto)  | yes           | no                      | no                                 |
+------------------------------------------------------+---------------+-------------------------+------------------------------------+
| Multiple images at different vantage point at once   | yes           | yes                     | yes                                |
+------------------------------------------------------+---------------+-------------------------+------------------------------------+
| Local observer                                       | yes           | yes                     | no                                 |
+------------------------------------------------------+---------------+-------------------------+------------------------------------+



Local observer mode
===================

The local observer mode (Sect. :ref:`sec-local-observer`) is a special mode
for putting the observer in the near-field of the object, or even right in
the middle of the object. It is not meant to be really for science use
(though it can be used for it, to a certain extent), but instead for 
public outreach stuff. However, it is kept relatively basic, because to
make this mode compatible with all the functions of RADMC-3D would require
much more development and that is not worth it at the moment. So here are
the restrictions:

+---------------------------------------------------------+---------------------+
| Option/Method:                                          | Local observer mode |
+---------------------------------------------------------+---------------------+
| Dust isotropic scattering                               | yes                 |
+---------------------------------------------------------+---------------------+
| Dust an-isotropic scattering                            | no                  |
+---------------------------------------------------------+---------------------+
| Multiple images at different vantage point at once      | yes                 |
+---------------------------------------------------------+---------------------+
| Second-order ray-tracing                                | yes                 |
+---------------------------------------------------------+---------------------+
| Doppler-catching of lines                               | no                  |
+---------------------------------------------------------+---------------------+
