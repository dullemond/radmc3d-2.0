.. _chap-development-history:

Version tracker: Development history
************************************

This version overview is very rough, and has only been started as of version
0.25.


* *Version 0.25*

  * Second order integration, based on a vertex-based grid (as opposed to 
    the usual cell-based grid), implemented. This gives much smoother images,
    and you don't see the blocky cell structure anymore in the images. It 
    requires extra memory, though. See Section :ref:`sec-second-order`.
  * The number of photons for scattering Monte Carlo (i.e. the small MC
    run done before each image, if dust scattering is active) can now be
    chosen to be smaller for when you make a spectrum instead of an
    image. Reason: Since you anyway integrate over the images for a
    spectrum, you do not need the image to 'look nice', i.e. you can
    afford more photon noise. You can set this in ``radmc3d.inp`` by
    setting ``nphot_spec=10000``\ , for instance. See Section
    :ref:`sec-scat-monte-carlo`.

* *Version 0.26*

  * For line transfer: Added the 'doppler catching' method to the
    code. This prevents bad numerical artifacts in images/spectra of regions
    with large velocity gradients, where the doppler-shift between two
    neighboring cells exceeds the intrinsic line width of the material in
    the cell. See Section :ref:`sec-doppler-catching`. 
  * NOTE: Up to, and including, version ``0.26_23.02.11`` this
    method (and for that matter any second order integration of line
    transfer) was not stable when strong shocks or contact discontinuities
    were encountered. This was because interpolation of the source function
    :math:`S_\nu\equiv j_\nu/\alpha_\nu` was done. Experimentation showed that
    interpolation of the emissivity :math:`j_\nu` is much more stable. As of
    version ``0.26_27.02.11`` this is fixed.

* *Version 0.27*
  
  * For line transfer: Implemented the possibility to use a Voigt line
    profile instead of just a Gaussian. This was implemented by Thomas
    Peters, and slightly modified by CPD. It uses the Voigt approximation by
    Humlicek JQSRT 27, 437 (1982) as programmed by Schreier, JQSRT 48, 743
    (1992). It requires a user-defined subroutine ``userdef_compute_lorentz_delta()``
    that sets the value of the Lorentz
    profile delta. This implementation is not yet documented, and may still
    be subject to modification. 
  * Implemented the 'Large Velocity Gradient' (LVG) method (also
    called the Sobolev method) of approximate non-LTE line transfer.
  * Implemented the optically thin populations method.
  * Implemented the possibility of reading linelist molecular data
    instead of full molecular data. *Still needs testing.*
  * Finally implemented the positive ``lines_mode`` modes,
    i.e. in which the level populations are computed and stored globally
    before the ray-tracing. This has been latently in the code somewhat, but
    unfinished. Now it is implemented. The advantage is: it may be under
    some conditions much faster than the on-the-fly computation of the
    populations during the ray-tracing (the negative ``lines_mode``
    modes). Also it allows you to write the populations to file, so that you
    can examine them. Disadvantage: It is memory hungry.
  * The level subset capacilities are now limited to only the storage of
    the levels in the global arrays (for positive ``lines_mode`` modes),
    and to the lines that will appear in
    images/spectra. For the rest, the full set of levels are always used
    from now on.
  * Added a directory '``opac/``\ ' which contains programs
    for generating your own dust opacities using optical constants from the
    web, and for generating your own molecular/atomic input data files using
    data from several web pages. The data from the web are not included, but
    there are README files that point you to the web sites.
  * Tested the 'fisheye' fulldome (OMNIMAX) projection. It seems to
    work! Thanks to Mario Flock.
  * Several small (and bigger) bugfixes
    
    * Fixed bug that showed up when no dust is included.
    * Fixed bug that caused RADMC-3D to crash when using no stars.
    * Fixed bug that caused RADMC-3D to crash when making images at very
      short wavelengths with nearly zero thermal emission.
    * Fixed bug in the AMR module when using second order integration or
      the doppler catching method with certain kinds of AMR-arrangements of
      cells.
    * Fixed many bugs when using a 'piece of a cake' model, i.e.
      using spherical coordinates in 3-D, but having the :math:`\phi`-grid going
      not over the full :math:`0-2\pi` range but e.g. just from 0 to :math:`\pi/4`. It
      is rather rare that one really wants to use such grids (certainly not
      for real physical models, I presume), but for visualization of data it
      might be useful: for instance for visualizing a 3-D disk MHD model,
      which is cut open so you can see also to the midplane. Now it
      works. Thanks to Mario Flock.
    * Fixed bug with the aperture mode for spectra. Thanks to Daniel Harsono.
    * Fixed many bugs in linelist mode; now it works. Thanks to Attila
      Juhasz.
    * Fixed a bug in LVG mode that caused it to fail when AMR was used.
      Thanks to Anika Schmiedeke.
    * Fixed a tiny bug in ``idl/radmc3dfits.pro``\ : filename was
      unused. Thanks to Stella Offner.
    * Retroactive bugfix from version 0.28 (see below): LVG and AMR mode.
    
    For details and for smaller bugfixes, read the ``src/Radmc_3D_LOG.txt`` document. 
  
* *Version 0.28*
  
  * A number of people complained that even without AMR the code
    requires a huge amount of memory. That is because even if no AMR is
    used, the cells are connected via the AMR tree. Since the AMR cells
    contain information about which are the neighboring cells, and each cell
    has 6 neighbors, and slots for 8 child-cells (which are unused in case
    of a regular grid) this wastes a lot of memory space. The first big
    improvement in version 0.28 is that, from now on, the AMR tree is only
    set up and used if the grid indeed has refinement. If RADMC-3D notices
    that the grid is regular, it will not allocate space for the AMR tree,
    and everywhere in the code where the cell-management is done the code
    will switch to regular grid mode. There is now a flag ``amr_tree_present``
    that says whether the AMR tree is present or
    not. Throughout the code there are now if-statements to switch between
    using or not-using the AMR tree. This may make the code a tiny bit
    slower, but this is only a minor reduction of speed. But as a result it
    should now be much easier to load huge regular grid models into memory.
  * A small (but potentially nasty) bug was found and fixed for the case
    when you use LVG mode on a grid with AMR-refinements. For the regular
    grid case (even in version 0.27, when it still used the AMR tree) this
    bug should not have caused problems, but perhaps you might want to check
    nevertheless. Note: This bug is now also retroactively fixed in version
    0.27. See, as always, ``src/Radmc_3D_LOG.txt`` for details.
  * Added the possibility to visualize the location (along the line of
    sight) of the :math:`\tau=1` surface (or any :math:`\tau=\tau_s` surface for that
    matter). See new Section :ref:`sec-tausurf`. This can be very useful
    for getting a 3-D feeling for *where* certain emission comes from.
  
* *Version 0.29*
  
  * The big change in this version is that the whole stuff with the
    global storage of level populations has been improved. In earlier
    versions of RADMC-3D, either the populations of *all* levels of a
    molecule were stored globally (potentially requiring huge amounts of
    memory), or you would have to select a 'subset' of levels to store
    globally. This subset selection had to be done by the user
    ('manually', so to speak). You would have had to think a-priori which
    lines you wish to model, and which levels they connect, and then, in the
    ``lines.inp`` file you would have to select these levels by
    hand. That was cumbersome and prone to error. To avoid having to do this
    you could use 'on-the-fly' calculation of populations (by making the
    ``lines_mode`` negative), but that sometimes caused the code to
    become terribly slow. *Now this is dramatically improved:* From now
    on you can forget about the 'on-the-fly' calculation of populations.
    Just use the 'normal' way by which RADMC-3D first calculates the
    populations and then starts the ray-tracing. The subset-selection is now
    done automatically by RADMC-3D, based on which wavelengths you want to
    make the image(s) or spectra for (see Section
    :ref:`sec-calcstore-levpop`). Now the on-the-fly methods are no longer
    default and should not be used, unless absolutely necessary. Also the
    'manual' subset selection is no longer necessary (though still
    possible if absolutely desired).
  * Added the subbox and sample capabilities to the level populations.
    See Sections :ref:`sec-subbox` and :ref:`sec-sampling`. Note that, in
    order to make it easier to identify which levels were written to file,
    the file formats of ``***_subbox.out`` and ``***_sample.out`` have been
    slightly modified: A list of identification
    numbers is added before the main data. For the dust temperature and dust
    density this list is simply 1 2 3 4 ...  (dust species 1, dust species
    2, dust species 3 ...), which is trivial. For the level populations
    (e.g. the file ``levelpop_co_subbox.out`` and ``levelpop_co_sample.out``
    for the CO molecule) this list is, however,
    essential when not all levels were computed (see Section
    :ref:`sec-calcstore-levpop`). So if only level 4 and level 8 are stored,
    then the identification list is 4 8. 
  * Fixed a bug which caused the code to crash when you put a star
    substantially far outside of the domain and try to make an image or
    spectrum. Thanks, Erika Hamden, for the bug report.
  * Fixed a bug that prevented the ``lines_mode=50`` mode from
    working. Now it works, and we can ask RADMC-3D to read the level
    populations from file (rather than calculating them internally).  Also a
    new section was added to this manual describing this option (Section
    :ref:`sec-nonlte-read-levelpop`).
  * Added VTK output options (see chapter :ref:`chap-vtk-output`) for
    allowing 3-D visualization of your model setups using e.g. Paraview, a
    freely available visualization tool.
  * Fixed a bug that occurred sometimes if a spectrum was made at
    inclination 90 and phi 90. Thanks Stella Offner for reporting this bug.
  
* *Version 0.30*
  
  * Fixed bugs in the Henyey-Greenstein scattering mode.
  * Introduced the new binary I/O feature: No more hassle with
    f77-unformatted records! The new binary mode is much simpler and more
    straightforward. This will help reducing the file sizes for large models.
    See Chapter :ref:`chap-binary-io`. 
  
* *Version 0.31*
  
  * Added the possibility, in cartesian coordinates, to 'close the
    box', in the sense of making the domain boundaries thermal walls.
    Each of the 6 boundaries can be set separately, so you can also have
    just one thermall wall. Also the temperatures can be set separately
    for each of the 6 boundaries. See Section :ref:`sec-thermal-boundaries`.
    
  * Added two new coordinate systems:
    
    * Cartesian 1-D plane-parallel (the only remaining active coordinate
      is :math:`z`). The :math:`x` and :math:`y` dimensions are infinitely extended and 
      have translational symmetry. The photons can, however, travel in 
      full 3-D as always. See Section :ref:`sec-1d-plane-parallel`.
    * Cartesian 2-D pencil-parallel (the two remaining active coordinate
      are :math:`y` and :math:`z`). The :math:`x` dimension is infinitely extended and 
      has translational symmetry. The photons can, however, travel in 
      full 3-D as always.
    * For the 1-D plane-parallel mode it is possible to include parallel
      beams of radiative flux impinging on the 1-D atmosphere.
    * Attila Juhasz has improved the VTK output: Now it also supports
      3-D spherical coordinates. Thanks, Attila!
    
  
* *Version 0.32*
  
  This is an intermediate version in which some stuff for the near-future
  modus of polarization is implemented.
  
* *Version 0.33*
  
  * Some minor technical changes to the doppler-catching integration of
    lines (storing the upper and lower level population instead of the
    jnubase and anubase variables). 
  * Added the classical escape probability to the LVG mode (see Section
    :ref:`sec-lvg` for details).
  * Sped up the filling of the matrix of the statistical equilibrium
    equation.
  * Vastly improved the LVG (and esc prob) method: Instead of the simple
    'lambda iteration style' iteration as it was before, the :math:`A_{ik}` is
    now multiplied with :math:`\beta_{ik}` (the escape probability of the line
    i->k) and the :math:`J_{ik}` is replaced by
    :math:`J_{ik}^{\mathrm{background}}`. This means that the solution is almost
    instant, requiring only 2 or 3 iterations.
  
* *Version 0.34*
  
  Implemented the Modified Random Walk method, based on Min, Dullemond,
  Dominik, de Koter \& Hovenier (2009) A\&A 497, 155, and simplified by
  Robitaille (2010) A\&A 520, 70. But beware: Still in the testing
  phase! By default it is switched off.
  
* *Version 0.35*
  
  * Implemented polarized scattering off randomly oriented
    particles. But beware: Still in the testing phase!
  * Fixed a bug in the modified random walk method (thanks to Daniel
    Harsono for spotting the problem and thanks to Attila Juhasz for
    finding the fix!)
  * Fixed two bugs that made it impossible to use second order
    integration with axially symmetric spherical coordinates and/or
    a finite-size star (thanks to Rolf Kuiper for reporting the bug). 
  * Added the ``sloppy`` command line option to spectrum and
    image making in spherical coordinates. This was necessary because
    RADMC-3D is always trying to make 100\% sure that all cells are picked
    up by the subpixels. In spherical coordinates these cells can be
    extremely non-cubic (they can be extremely flat or needle-like), which
    means that under some projections RADMC-3D feels obliged to do extreme
    sub-pixeling, which can make image- and spectrum-making extremely slow.
    By adding the ``sloppy``  keyword on the command line, RADMC-3D
    will limit it's pubpixeling which could speed up the calculation very
    much (but of course at your own risk!).
  
* *Version 0.38*
  
  * Implemented OpenMP parallellization of the thermal Monte Carlo (by
    Adriana Pohl). Still beta-version.
  * Bugfix in the mean intensity computation (mcmono) mode (thanks to Gwendoline Stephan).
  * Bugfix in the mean intensity computation (mcmono) mode (thanks to Seokho Lee).
  * Major bugfix in aperture mode (thanks to S\o ren Frimann). 
  * Unformatted image format is from now on C-style binary instead of
    F77-style unformatted.
  * The viewimage tool is now ported to Qt by Farzin Sereshti, meaning
    that you can now use viewimage without having an IDL license. Viewimage
    is a very powerful tool to interactively make and view images of your
    model at different wavelengths and viewing angles. It can be found
    in the directory ``viewimage_QT_GUI/``\ .
  * A Python package for RADMC-3D was developed by Attila Juhasz. It is
    included as of RADMC-3D version 0.38 in the directory ``python/``\ .
  
* *Version 0.39*
  
  * Polarization mode is incompatible with mirror mode (in spherical
    coordinates). An error message is now included to catch this.
  * Minor bugfix in ``pick_randomfreq_db()`` (thanks to Seokho Lee).
  * Optimization of the OpenMP parallellization and extension of the
    OpenMP parallellization to the Scattering Monte Carlo computation (both
    by Farzin Sereshti).
  * Bugfix in ``amrray_module.f90``\ : Sometimes one got 'Photon
    outside of cell' error due to a numerical precision round-off
    error. This bug is now (mostly?) fixed.
  * Bugfix in ``sources_module.f90``\ : When using second order
    integration (or doppler catching) for line transfer in spherical
    coordinates, the line doppler shift was not transformed to spherical
    coordinates. This is now fixed.
  * Several bugfixes in the modified random walk method by 
    John Ramsey. The method crashed for extreme optical depth problems
    due to out-of-cell events. Still not 100\% perfect, but better.
  * John Ramsey also proposed two small fixes to the Planck function
    routines so that the events of overflow are caught. Note: This might
    change the results (in a tiny way: at the machine precision level) to
    the extent that a model run by an old version might not yield the same
    values to machine precision, but the differences should not matter in
    any meaningful way.
  
* *Version 0.40*
  
  * The RADMC-3D package is now 'officially' converting from IDL to
    Python wrappers. The Python modules were already there since a long time
    (thanks to Attila Juhasz!). But as of version 0.40 we will no longer
    update/maintain the IDL scripts (though they remain there and should
    remain working), and instead use python as the main setup and analysis
    tools for RADMC-3D. The full conversion will still take some time, but
    should be finished by the end of version 0.40. 
  * Under some circumstances the simple 2x2 pixel plus sub-pixeling
    method for making spectra (default method) can be dangerous. For some
    grid geometries this can lead to under-resolving of the images that are
    integrated to obtain the flux, leading to a too low flux.  So as of now
    15.09.2016 the spectra and SEDs are always by default made with 100x100
    images (and sub-pixeling of course). One can set the number of pixels
    with npix. So if you do ``radmc3dsednostarnpix2`` you get
    the original behavior again. 
  * Bugfix in ``montecarlo_module.f90``\ :
    The internal heat source method (which is still being tested)
    had a bug. The bug manifested itself for optically thin cells with
    non-negligible internal heat production. The energy was not immediately
    added to the cell. It only got added upon re-absorption of that photon
    package. Now this is fixed.
  * I now added some documentation for the heat source method, which
    is useful for e.g. disk viscous accretion heating.
  * Bugfix in ``montecarlo_module.f90``\ : When using mirror
    symmetry in spherical coordinates in the :math:`\theta`-coordinate
    (i.e. modeling only the upper part of the disk and letting RADMC-3D
    assume that the lower part is identical), the distributed source
    luminosity was computed only for the top quadrant, and wasn't multiplied
    by 2. For most applications this does not cause problems, but for the
    heat source (see above), for continuous stellar sources and for the
    thermal origin of the isotropic scattering luminosity (for non-isotropic
    scattering, mirror symmetry was not allowed anyway), this could lead to
    a factor of 2 underestimation (only if mirror symmetry was used, i.e.
    if the :math:`\theta` coordinate was going only up to :math:`\pi/2`). This is now
    fixed. To test if the fix works one can simply make the same model
    again, but now without using mirror symmetry (and thus using twice as
    many cells in :math:`\theta`, to cover both the upper and lower half of the
    object). This should yield (apart from some Monte Carlo noise) the same
    results.
  * Improved the stability of the Modified Random Walk (MRW) method
    a bit further. 
  * Bug fix: scattering mode 3 (tabulated phase function, but not full
    polarization) had a bug which caused images of scattered light to be
    multiplied by some arbitrary number. Reason: as a phase function it
    returned :math:`Z_{11}` instead of :math:`4\pi Z_{11}/\kappa_{\mathrm{scat}}`. Most
    people use either isotropic scattering (scattering mode 1), or
    Henyey-Greenstein (scattering mode 2) or full polarization (scattering
    mode 5), all of which are ok. At any rate: the problem is now fixed,
    so scattering mode 3 should now also work.
  
* *Version 0.41*
  
  * Implemented a first testing version of the aligned grains:
    only polarized thermal emission so far. Still very much a testing
    version.
  * Implemented a method to also allow full Stokes vector polarized
    scattering in the 2-D axisymmetric mode in spherical
    coordinates. Until now the full scattering mode (scattering mode 5)
    was only possible in full 3-D. Note however that anisotropic
    scattering in 2-D axisymmetric models requires scattering mode 5,
    which is the full scattering mode.  It is still not possible to use
    intermediate scattering modes (like henyey-greenstein or any
    scattering mode between 2 and 4) in 2-D axisymmetry. But those
    intermediate modes are anyway more for testing than for real models,
    so that should be ok.
  * Bugfixes to the OpenMP stuff. In particular the OpenMP 
    parallellization of the scattering MC crashed. This is now fixed.
    In general the OpenMP stuff was a bit cleaned up.
  * Bugfix in thermal Monte Carlo with full polarization mode: needed
    to reset the photon package after each thermal absorption/re-emission
    event. Usually the effect is subtle, but had to be fixed.
  * Bugfix in reading the ``scattering_angular_grid.inp``\ :
    the ``theta`` angles should be converted into radian. But
    this file was not officially offered before anyway. 
  * Attila Juhasz has made a large improvement of his python package
    for RADMC-3D. See the ``python/`` directory. This is version
    0.29 of his package. This package now also supports reading and
    writing AMR grids.
  * Bugfix in VTK for 3-D spherical coordinates (thanks Attila Juhasz!).
    Now it should work.
  
* *Version 2.0*

  Version 2.0 is the version after 0.41. We skip version 1.0, because
  version 1.0 could be mistaken for the first version of the code. Version
  2.0 is mostly the same as 0.41, but with a few differences.
  
  * IDL support is removed permanently. From now on, the front-end
    functionality is only in Python. We assume Python 3.
  * Version 0.30.2 of the ``radmc3dPy`` Python package (written by
    Attila Juhasz) has been implemented. It is also being improved, mainly to
    make its use easier (i.e. with more automatic default behavior).
  * A very simple ``simpleread.py`` reading library is provided
    as a 'light version' of radmc3dPy. It contains only some basic reading
    functions, and only for ascii output (no binary files).
  * Some of the standard-output is shortened. You can also call a
    Monte Carlo run with ``radmc3d`` with the command line
    options ``countwrite 100000`` to make RADMC-3D write a
    message only every :math:`10^5` photon packages instead of every thousand.
  * We removed the fortran-unformatted data format from the manual,
    and will remove it from the code in later versions. Use either text
    (ascii) format or binary format. 
  * The manual is now converted to Sphinx, from which the LaTeX version
    and the HTML version can be automatically created.
  * [as of 11.11.2021] BUGFIX: For OpenMP parallel thermal Monte Carlo
    computation of the dust temperatures for multiple grain species or
    sizes, when ``iranfreqmode=1`` (as opposed to the default value of
    ``iranfreqmode=0``), the dust temperatures could acquire errors
    because the ``pick_randomfreq_db()`` subroutine uses the array
    ``db_cumul(:)`` as thread private, but without having it declared
    as such. This led to interference between threads. This is now
    fixed.
