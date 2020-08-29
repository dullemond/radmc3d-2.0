.. _chap-stars:

More information about the treatment of stars
*********************************************

How stars are treated in RADMC-3D is perhaps something that needs some more
background information. This is the structure:

#. *Stars as individual objects:*
   
   The most standard way of injecting stellar light into the model is by putting
   one or more individual stars in the model. A star can be placed anywhere,
   both inside the grid and outside. The main input file specifying their
   location and properties is: ``stars.inp``\ . The stars can be treated in two
   different ways, depending on the setting of the variable ``istar_sphere``
   that can be set to 0 or 1 in the file ``radmc3d.inp`` file.

   * The default is to treat stars as zero-size point sources. This is the way
     it is done if (as is the default) ``istar_sphere=0``\ .  The stars are then
     treated as point sources in spite of the fact that their radius is
     specified as non-zero in the ``stars.inp`` file.  This default mode is the
     easiest and quickest. For most purposes it is perfectly fine. Only if you
     have material very close to a stellar surface it may be important to treat
     the finite size(s) of the star(s).
      
   * If ``istar_sphere=1`` in the ``radmc3d.inp`` file, then all stars are
     treated as spheres, their radii being the radii specified in the
     ``stars.inp`` file. This mode can be tricky, so please read Section
     :ref:`sec-stars-as-spheres`.

#. *Smooth distributions of zillions of stars:*
   
   For modeling galaxies or objects of that size scale, it is of course
   impossible and unnecessary to treat each star individually. So *in addition
   to the individual stars* you can specify spatial distributions of stars,
   assuming that the number of stars is so large that there will always be a
   very large number of them in each cell. Please note that using this
   possibility does *not* exclude the use of individual stars as well. For
   instance, for a galaxy you may want to have distributions of unresolved
   stars, but one single 'star' for the active nucleus and perhaps a few
   individual 'stars' for bright star formation regions or O-star clusters or
   so. The distribution of stars is described in Section
   :ref:`sec-distrib-of-stars`.
  
#. *An external 'interstellar radiation field':*
   
   Often an object is affected not only by the stellar radiation from the stars
   inside the object itself, but also by the diffuse radiation from the many
   near and far stars surrounding the object. This 'Interstellar Radiation
   Field' can be treated by RADMC-3D as well. This is called the 'external
   source' in RADMC-3D. It is described in Section :ref:`sec-external-source`.


.. _sec-stars-as-points:

Stars treated as point sources
==============================

By default the stars are treated as point-sources. Even if the radius is
specified as non-zero in the ``stars.inp`` file, they are still
treated as points. The reason for this is that it is much easier and faster
for the code to treat them as point-sources. Point sources cannot occult
anything in the background, and nothing can partly occult them (they are
only fully or not occulted, of course modulo optical depth of the occulting
object). This approximation is, however, not valid if the spatial scales you
are interested in are not much larger (or even the same or smaller) than the
size of the star. For instance, if we are interested in modeling the
radiative transfer in a disk around a Brown Dwarf, where dust can survive
perhaps even all the way down to the stellar surface, we must take the
non-point-like geometry of the star into account. This is because due to its
size, the star can shine *down* onto the disk, which would not be
possible if the star is treated as a point source. However, for a dust disk
arounda Herbig Ae star, where the dust evaporation radius is at about 0.5
AU, the star can be treated as a point-source without problems.

So if you just use RADMC-3D as-is, or if you explicitly set ``istar_sphere=0``
in the file ``radmc3d.inp``\ , then the stars are all treated as point sources.



.. _sec-stars-as-spheres:

Stars treated as spheres
========================

For problems in which the finite geometrical size of the star (or stars)
is/are important, RADMC-3D has a mode by which the stars are treated as
spheres. This can be necessary for instance if you model a disk around
a Brown Dwarf, where the dusty disk goes all the way down to the stellar
surface. The finite size of the star can thus shine *down* onto the
disk, but only if its finite size is treated as such. In the default
point-source approximation the surface layers of such a disk would be
too cold, because this 'shining down onto the disk' phenomenon is
not treated.

You can switch this mode on by setting ``istar_sphere=1`` in the
file ``radmc3d.inp``\ . Note that no limb darkening or brightening is
included in this mode, and currently RADMC-3D does not have such a mode
available. 

This mode is, however, somewhat complex. A sphere can partly overlap the
grid, while being partly outside the grid. A sphere can also overlap
multiple cells at the same time, engulfing some cells entirely, while only
partly overlapping others. The correct and fast treatment of this makes
the code a bit slower, and required some complex programming. So the user
is at the moment advised to use this mode only if necessary and remain
aware of possible errors for now (as of version 0.17). 

For the Monte Carlo simulations the finite star size means that photon
packages are emitted from the surface of the sphere of the star. It also
means that any photon that re-enters the star during the Monte Carlo
simulation is assumed to be lost.




.. _sec-distrib-of-stars:

Distributions of zillions of stars
==================================

For models of galaxies it is important to be able to have distributed
stellar sources instead of individual stars. The way to implement this
in a model for RADMC-3D is to 

#. Prepare one or more *template stellar spectra*, for instance, one for each
   stellar type you wish to include. These must be specified in the file
   ``stellarsrc_templates.inp`` (see Section
   :ref:`sec-stellarsrc-templates`). Of course the more templates you have, the
   more memory consuming it becomes, which is of particular concern for models
   on large grids. You can of course also take a sum of various stellar types as
   a template. For instance, if we wish to include a 'typical' bulge stellar
   component, then you do not need to treat each stellar type of bulge stars
   separately. You can take the 'average spectrum per gram of average star' as
   the template and thus save memory.
  
#. For each template you must specify the *spatial distribution*,
   i.e. how many stars of each template star are there per unit volume in
   each cell. The stellar density is, in fact, given as gram-of-star/cm\ :math:`^3`
   (i.e. not as number density of stars). The stellar spatial densities
   are specified in the file ``stellarsrc_density.inp`` (see
   Section :ref:`sec-stellarsrc-density`).

Note that if you have a file ``stellarsrc_templates.inp`` in your
model directly, then the stellar sources are automatically switched on.
If you do not want to use them, then you must delete this file. 

The smooth stellar source distributions are nothing else than source
functions for the radiative transfer with the spectral shape of the template
stellar spectra from the ``stellarsrc_templates.inp``\ .  You will
see that if you make a spectrum of your object, then even if the dust
temperature etc is zero everywhere, you still see a spectrum: that of the
stellar template(s). In the Monte Carlo simulations these stellar templates
act as net sources of photons, that subsequently move through the grid in a
Monte Carlo way. 

Note that the smooth stellar source distributions assume that the zillions
of stars that they represent are so small that they do not absorb any
appreciable amount of radiation. They are therefore pure sources, not sinks.



.. _sec-external-source:

The interstellar radiation field: external source of energy
===========================================================

You can include an *isotropic* interstellar radiation field in
RADMC-3D. This will take effect both in the making of spectra and images, as
well as in the Monte Carlo module.

The way to activate this is to make a file ``external_source.inp``
and fill it with the information needed (see Section :ref:`sec-ext-src-inp`).


Role of the external radiation field in Monte Carlo simulations
---------------------------------------------------------------

For the Monte Carlo simulations this means that photons may be launched from
outside inward. The way that this is done is that RADMC-3D will make a
sphere around the entire grid, just large enough to fit in the entire grid
but not larger. Photon packages can freely leave this sphere. But if
necessary, photon packages can be launched from this sphere inward.
RADMC-3D will then calculate the total luminosity of this sphere, which is
:math:`L=4\pi^2 I r_{\mathrm{sphere}}^2` where :math:`I` is the intensity. For
monochromatic Monte Carlo it is simply :math:`I=I_\nu`, while for the thermal
Monte Carlo it is :math:`I=\int_0^\infty I_\nu d\nu`, where :math:`I_\nu` is the
intensity as specified in the file ``external_source.inp``\ .  Note
that if the sphere would have been taken larger, then the luminosity of the
external radiation field would increase. This may seem anti-intuitive. The
trick, however, is that if the sphere is larger, then also more of these
interstellar photons never enter the grid and are lost immediately. That is
why it is so important that RADMC-3D makes the sphere as small as possible,
so that it limits the number of lost photon packages. It also means that you,
the user, would make the grid much larger than the object you are interested
in, then RADMC-3D is forced to make a large sphere, and thus potentially
many photons will get lost: they may enter the outer parts of the grid, but
there they will not get absorbed, nor will they do much. 

In fact, this is a potential difficulty of the use of the external sources:
since the photon packages are lauchned from outside-inward, it may happen that
only few of them will enter in the regions of the model that you, the user, are
interested in. For instance, you are modeling a 3-D molecular cloud complex with
a few dense cold starless cores. Suppose that no stellar sources exist in this
model, only the interstellar radiation field. The temperature in the centers of
these starless cores will be determined by the interstellar radiation field. But
since the cores are very small compared to the total model (e.g. you have used
AMR to refine the grid around/in these cores), the chance of each external
photon package to 'hit' the starless core is small. It means that the larger the
grid or the smaller the starless core, the more photon packages (``nphot``,
see Section :ref:`sec-dust-thermal-monte-carlo`) one must use to make sure that
at least some of them enter the starless cores. If you choose ``phot`` too small
in this case, then the temperature in these cores would remain undetermined
(i.e. they will be zero in the results).



Role of the external radiation field in images and spectra
----------------------------------------------------------

The interstellar radiation field also affects the images and spectra that
you make. Every ray will start at minus-infinity with the intensity given by
the external radiation field, instead of 0 as it would be if no external
radiation field is specified. If you make an image, the background of your
object will then therefore not be black. You can even make silhouette images
like those of the famous silhouette disks in Orion. 

But there is a danger: if you make spectra, then also the background 
radiation is inside the beam, and will thus contribute to the spectrum.
In fact, the larger you make the beam the more you will pick up of the
background. This could thus lead to the spectrum of your source to be
swamped by the background if you do not specify a beam in the spectrum.



.. _sec-internal-source:

Internal heat source
====================

Sometimes the gas and dust inside the object of interest gets heated up by
some internal process such as friction, magnetic reconnection, chemical
reactions, etc. A nice example is the 'viscous heating' inside an
accretion disk. This net heat source can be included in RADMC-3D by creating
a file ``heatsource.inp``\ . The format of the file is described in
Section :ref:`sec-heatsource`. It is the same as for other scalar fields.

With this input file you have to specify in each cell how much energy per
second per cubic centimeter is released in the form of heat. This energy
will then be emitted as radiation by the dust. The way the code does this in
the Bjorkman \& Wood algorithm is that it will launch photon packages from
these cells. The difference with the stellar energy input (see Section
:ref:`sec-distrib-of-stars`) is that the energy is first injected into the
dust of the cell, and then emitted as thermal dust emission. The launching
of the photon package is therefore always a thermal dust emission. In
contrast, in the stellar energy input method of Section
:ref:`sec-distrib-of-stars` the photon package is launched directly, with a
wavelength randomly drawn from the local stellar spectrum shape. The
difference between these two methods will be most apparent for optically
thin models. For very optically thick cases, where the heat source is
released deep inside an optically thick object, both methods will
presumably yield the same result. Nevertheless, it is recommended to
use the heat source method for cases such as chemical or viscous heating
of the gas and dust, even for optically thick cases. 

A note of caution: in spite of the fact that this heat source method allows
you to add additional energy sources, the object of study must still be in
local thermodynamic equilibrium (LTE). If the gas+dust mixture is flowing
and experiences significant adiabatic heating and cooling events, then the
LTE condition is no longer met and RADMC-3D will not be able to give
reliable answers. Sometimes one may be able to fudge this in some clever
way, but one should always be aware that strictly speaking the Bjorkman \&
Wood Monte Carlo method only works if in each cell all energy input (be it
radiative absorption or an internal heat source) is balanced exactly by the
same amount of radiative energy output. The algorithm  computes
the dust temperature on that assumption: it computes how much energy the
cell gains (by the heat source or by absorbing photons) and then it requires
that the temperature of the dust is such that precisely the same amount of
radiative energy is emitted. 



Slow performance of RADMC-3D with heat source
---------------------------------------------

For very optically thick models, such as the inner regions of actively
accreting dusty protoplanetary disks, the use of this heat source can lead
to extremely slow performance. The reason is that all photons originating
from this heat source will start their journey right in the middle of the
most optically thick regions, requiring these photons to make gazillions of
absorption/re-emission events before finally diffusing out. It should in
principle work if the code runs long enough. But one must have some
patience. The use of the Modified Random Walk method (see Section
:ref:`sec-modrandwalk`) would then be useful to speed things up, but still it
can take time.

A few things might be useful to consider. One is that protoplanetary disks
only have such insane optical depths (:math:`\tau\gtrsim 10^5`) if none of the
dust has coagulated to bigger grains. This might be the correct assumption,
especially in the very early phases of protoplanetary disk evolution. But
dust coagulation is known to be quick, so it might equally well be that,
say, 90\% of the small grain dust has already grown to larger grains, which
have less opacity. This is of course just a pure guess. Another thing is
that many MHD models of disk turbulence show that most of the energy is 
not released near the midplane, but instead at one or two scale heights
above the midplane. Both considerations would lower the optical depth for
the energy to get out of the disk, speeding up the calculation. And 
the outcoming spectrum or image it will presumably not be affected that
much, because at the end of the day the effective temperature of the disk
surface must anyway be such that it radiates away the internal heat,
independent of how deep inside the disk this heat is released. 
