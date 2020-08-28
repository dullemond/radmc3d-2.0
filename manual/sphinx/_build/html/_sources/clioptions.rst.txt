.. _chap-command-line-options:

Command-line options
********************

This chapter deals with all the possible command-line options one can 
give when calling the ``radmc3d`` code.


Main commands
=============

In addition to the radmc3d.inp file, which contains many 'steering' parameters,
one can (and even must) give RADMC-3D also command-line options. The most
important (and compulsory) options are the 'command' what RADMC-3D should do.
At the moment you can choose from:

* ``mctherm``\ : Runs RADMC-3D for computing the dust
  temperatures using the Monte Carlo method.
  
* ``spectrum``\ : Runs RADMC-3D for making a spectrum based
  on certain settings. This option requires further command-line
  specifications. See chapter :ref:`chap-images-spectra`.
  
* ``sed``\ : Runs RADMC-3D for making a SED based on certain
  settings. This option requires further command-line specifications.  Note
  that a SED is like a spectrum, but for continuum processes only (no
  lines). See chapter :ref:`chap-images-spectra` for more details.
  
* ``image``\ : Runs RADMC-3D for making an image.  This option
  requires further command-line specifications. See chapter
  :ref:`chap-images-spectra`.
       
* ``movie``\ : Like ``image``\ , but now for a series of different
  vantage points. Useful for making movies in one go, without having to call
  RADMC-3D time and again. *NOTE: This command is still under
  development*. See chapter :ref:`chap-images-spectra`.
  
* ``mcmono``\ : (Only expect use). Runs RADMC-3D for computing
  the local radiation field at each location in the model. This is only
  useful for when you wish to couple RADMC-3D to models of chemistry or so,
  which need the local radiation field. See Section
  :ref:`sec-dust-monochromatic-monte-carlo`.

Example::

  radmc3d mctherm

runs the RADMC-3D code for computing the dust temperatures everywhere using the
Monte Carlo method.

There are also some additional commands that may be useful for diagnostics:

* ``subbox_****``\ : where \*\*\*\* is one of the following:
  ``dust_density``\ , ``dust_temperature``\ . But other
  quantities will follow in later versions. See Section :ref:`sec-subbox`.
  
* ``linelist``\ : Write a list of all the lines included in
  this model.



Additional arguments: general
=============================

Here is a list of command line options, on top of the above listed main
commands (Note: We'll try to be complete, but as the code develops we may
forget to list new options here):

* ``setthreads`` [for MC] The next number sets the number of OpenMP parallel threads
  to be used.
  
* ``npix``\ : [for images] The next number specifies the number of
  pixels in both x and y direction, assuming a square image.

* ``npixx``\ : [for images] The next number specifies the number of
  pixels in x direction only. 

* ``npixy``\ : [for images] The next number specifies the number of
  pixels in y direction only.

* ``nrrefine``\ : [for images and spectra] Specifies a maximum depth of
  refinement of the pixels (see Section :ref:`sec-image-refinement`).

* ``fluxcons``\ : [for images and spectra] Puts nrrefine (see above) to
  a large value to assue flux conservation (see Section :ref:`sec-image-refinement`).

* ``norefine``\ : [for images and spectra] Puts
  nrrefine (see above) to 0 so that each pixel of the image corresponds only
  to 1 ray. This is fast but not reliable and therefore not recommended (see
  Section :ref:`sec-image-refinement`).

* ``nofluxcons``\ : [for images and spectra] As
  ``norefine`` above.

* ``noscat``\ : This option makes RADMC-3D ignore the
  dust scattering process (though not the scattering extinction!) in the
  images, spectra and Monte Carlo simulations. For images and spectra this
  means that no scattering Monte Carlo run has to be performed before each
  image ray tracing (see Section :ref:`sec-scat-monte-carlo`). This can speed
  up the making of images or spectra enormously. This is even more so if you
  make images/spectra of gas lines with LTE, LVG or ESCP methods, because if
  no scattering Monte Carlo needs to be made, ray-tracing can be done
  multi-frequency for each ray, and the populations can be calculated once
  in each cell, and used for all frequencies. That can speed up the line
  rendering enormously -- of course at the cost of not including dust
  scattering. For lines in the infrared and submillimeter, if no large
  grains are present, this is usually OK, because small grains (smaller than
  about 1 micron) have very low scattering albedos in the infrared and
  submillimeter.

* ``ilambda`` or ``inu``\ : [for images] Specify the index of the wavelength
  from the ``wavelength_micron.inp`` file for which a ray-trace image should be
  made.

* ``color``\ : [for images] Allows you to make multiple
  images (each at a different wavelength) in one go. This will make RADMC-3D
  read the file ``color_inus.inp`` (see Section
  :ref:`sec-minor-input-files`) which is a list of indices ``i``
  referring to the ``wavelength_micron.inp`` file for which the
  images should be made. See Section :ref:`sec-set-camera-frequencies` for
  details.

* ``loadcolor``\ : [for images] Same as ``color``\ .

* ``loadlambda``\ : [for images] Allows you to make multiple images (each at a
  different wavelength) in one go. This will make RADMC-3D read the file
  ``camera_wavelength_micron.inp`` to read the precise wavelength points at which you wish
  to make the images. In contrast to ``loadcolor``\ , which only allows you to
  pick from the global set of wavelength used by the Monte Carlo simulation (in
  the file ``wavelength_micron.inp``), with the
  ``camera_wavelength_micron.inp`` files you can
  specify any wavelength you want, and any number of them. See Section
  :ref:`sec-set-camera-frequencies` for details.

* ``sizeau``\ : [for images and spectra] The next number
  specifies the image size in model space in units of AU (=1.496E13
  cm). This image size is measured from the image left to
  right and top to bottom. This gives always square images. This image 
  size in au is observer distance independent. The corresponding image 
  size in arcsec is: image size in arcsec = image size in AU /
  (distance in parsec).

* ``sizepc``\ : [for images and spectra] Same as ``sizeau``\ , but
  now in parsec units.

* ``zoomau``\ : [for images and spectra] The next four numbers set the image
  window precisely by specifying the xleft, xright, ybottom, ytop of the image
  in units of AU. The zero point of the image (the direction of the 2-D image
  point located at (0.0,0.0) in image coordinates) stays the same (i.e. it aims
  toward the 3-D point in model space given by ``pointau`` or ``pointpc``\ ). In
  this way you can move the image window left or with or up or down without
  having to change the ``pointau`` or ``pointpc`` 3-D locations. Also for local
  perspective images it is different if you move the image window in the image
  plane, or if you actually change the direction in which you are looking (for
  images from infinity this is the same). *Note:* If you use this option
  without the ``truepix`` option RADMC-3D will always make square pixels by
  adapting ``npixx`` or ``npixy`` such that together with the ``zoomau`` image
  size you get approximately square pixels. Furthermore, if ``truezoom`` is not
  set, RADMC-3D will alleviate the remaining tiny deviation from square pixel
  shape by slightly (!) adapting the ``zoomau`` window to obtain exactly square
  pixels.

* ``zoompc``\ : [for images and spectra] Same as ``zoomau``\ , but
  now the four numbers are given in units of parsec.

* ``truepix``\ : [for images and spectra] If with ``zoomau`` or ``zoompc`` the
  image window is not square then when specifying ``npix`` one gets non-square
  pixels. Without the ``truepix`` option RADMC-3D will adapt the ``npixx`` or
  ``npixy`` number, and subsequently modify the zoom window a bit such that the
  pixels are square. With the ``truepix`` option RADMC-3D will not change
  ``npixx`` nor ``npixy`` and will allow non-square pixels to form.

* ``truezoom``\ : [for images and spectra] If set, RADMC-3D will always assure
  that the exact zoom window (specified with ``zoomau`` or ``zoompc``\ ) 
  will be used, i.e. if ``truepix`` is *not* set but ``truezoom``
  is set, RADMC-3D will only (!) adapt ``npixx`` or ``npixy`` to get
  *approximately* square pixels.

* ``pointau``\ : [for images and spectra] The subsequent three
  numbers specify a 3-D location in model space toward which the camera is
  pointing for images and spectra.  The (0,0) coordinate in the image plane
  corresponds by definition to a ray going right through this 3-D point.

* ``pointpc``\ : [for images and spectra] Same as ``pointau`` but
  now in units of parsec.

* ``incl``\ : [for images and spectra] For the case when the camera
  is at infinity (i.e. at a large distance so that no local perspective has
  to be taken into account) this inclination specifies the direction toward
  which the camera for images and spectra is positioned. Incl = 0 means
  toward the positive :math:`z`-axis (in cartesian space), incl=90 means toward a
  position in the :math:`x`-:math:`y`-plane and incl=180 means toward the negative
  :math:`z`-axis. The angle is given in degrees.

* ``phi``\ : [for images and spectra] Like ``incl``\ , but now the remaining
  angle, also given in degrees. Examples: ``incl``\ =90 and ``phi``\ =0 means
  that the observer is located at infinity toward the negative :math:`y` axis;
  ``incl``\ =90 and ``phi``\ =90 means that the observer is located at infinity
  toward the negative :math:`x` axis; ``incl``\ =90 and ``phi``\ =180 means that
  the observer is located at infinity toward the positive :math:`y` axis
  (looking back in negative :math:`y` direction). Rotation of the observer
  around the object around the :math:`z`-axis goes in clockwise direction. The
  starting point of this rotation is such that for ``incl``\ =0 and ``phi=0``
  the :math:`(x,y)` in the image plane correspond to the :math:`(x,y)` in the
  3-D space, with :math:`x` pointing toward the right and :math:`y` pointing
  upward. Examples: if we fix the position of the observer at for instance
  ``incl``\ =0 (i.e. we look at the object from the top from the positive
  :math:`z`-axis at infinity downward), then increasing ``phi`` means rotating
  the object counter-clockwise in the image plane.

* ``posang``\ : [for images] This rotates the camera itself around
  the :math:`(0,0)` point in the image plane. 

* ``imageunform``\ : Write out images in binary format

* ``imageformatted``\ : Write out images in text form (default)

* ``tracetau``\ : [for images] If this option is set, then instead
  of ray-tracing a true image, the camera will compute the optical depth
  at the wavelength given by e.g. ``inu`` and puts this into an image
  output as if it were a true image. Can be useful for analysis of models.

* ``tracecolumn``\ : [for images] Like ``tracetau`` but instead of the optical
  depth the simple column depth is computed in
  :math:`\mathrm{g}/\mathrm{cm}^2`. *NOTE: for now only the column depth of the
  dust.*

* ``tracenormal``\ : [for images: Default] Only if you specified 
  ``tracetau`` or ``tracecolumn`` before, and you are in child mode, 
  you may sometimes want to reset to normal imaging mode.

* ``apert`` or ``useapert``\ : [for
  images/spectra] Use the image-plane aperture information from the file
  ``aperture_info.inp``\ .

* ``noapert``\ : [for images/spectra] Do *not* use an image-plane aperture.

* ``nphot_therm``\ : [for MC] The nr of photons for the thermal
  Monte Carlo simulation. But it is better to use the ``radmc3d.inp`` for this
  (see Section :ref:`sec-radmc-inp`), because then you can see afterward with
  which photon statistics the run was done.

* ``nphot_scat``\ : [for MC] The nr of photons for the
  scattering Monte Carlo simulation done before each image (and thus also in
  the spectrum). But it is better to use the ``radmc3d.inp`` for
  this (see Section :ref:`sec-radmc-inp`), because then you can see afterward
  with which photon statistics the run was done.

* ``nphot_mcmono``\ : [for MC] The nr of photons for
  the monochromatic Monte Carlo simulation. But it is better to use the
  ``radmc3d.inp`` for this (see Section :ref:`sec-radmc-inp`),
  because then you can see afterward with which photon statistics the run
  was done.

* ``countwrite``\ : [for MC] The nr of photons between
  'sign of life' outputs in a Monte Carlo run. Default is 1000. That means
  that if you have ``nrphot=10000000`` you will see ten-thousand
  times something like ``Photonnr: 19000`` on your screen. Can
  be annoying. By adding ``countwrite 100000`` to the command
  line, you will only see a message every 100000 photon packages.


Switching on/off of radiation processes
=======================================

You can switch certain radiative processes on or off with the following
command-line options (though often the ``radmc3d.inp`` file also allows this):


* ``inclstar``\ : [for images and spectra] Include stars in
  spectrum or images.

* ``nostar``\ : [for images and spectra] Do *not* include stars
  in spectrum or images. Only the circumstellar / interstellar material is
  imaged as if a perfect coronograph is used.

* ``inclline``\ : Include line emission and extinction
  in the ray tracing (for images and spectra). 

* ``noline``\ : Do not include line emission and extinction
  in the ray tracing (for images and spectra).

* ``incldust``\ : Include dust emission, extinction and
  (unless it is switched off) dust scattering in ray tracing (for images and
  spectra).

* ``nodust``\ : Do not include dust emission, extinction and
  scattering in ray tracing (for images and spectra).

* ``maxnrscat 0``\ : (if dust is included) Do not include
  scattering in the images/spectra created by the camera. With ``maxnrscat 1``
  you limit the scattering in the images/spectra to single-scattering.
  With ``maxnrscat 2`` to double scattering, etc. Can be useful to
  figure out the relative importance of single vs multiple scattering.

