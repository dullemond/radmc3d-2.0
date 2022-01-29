# RADMC-3D Version 2.0

## Introduction
RADMC-3D is a tool for astrophysical research. It computes the observational appearance of an astrophysical object on the sky of the observer. It solves the non-local radiative transfer problem of dusty media, including thermal radiative transport and scattering. 

For a very long time (from its inception in 2007 to 2020) it had a version number smaller than 1, because it contained numerous modes and features that were still experimental. The latest (and longest lasting) version was 0.41. But it was high time to make a new version. I decided to skip version 1, because that would give the impression that the code is still new. So we go straight to version 2.

Compared to version 0.41, the new version 2.0 is not very different. The basic Fortran-90 code is still the same. This means that you can use version 2.0 without having to worry that the new version has potential new bugs that were introduced since 0.41. Having said that: it does not mean that verion 2.0 is bug-free. 

What is new in version 2.0 is some of the Python front- and back-end. I also removed old stuff that is no longer necessary. And in version 2.0 IDL front/back-end is no longer supported (if you don't know what IDL is, no need to worry). 

## How to use RADMC-3D
To learn what can be done with RADMC-3D and how to use it, please consult the extensive manual in the directory `manual/`. This manual is written using Sphinx, and therefore has a latex version (`manual/sphinx/_build/latex/radmc3d.pdf`) and an HTML version (`manual/sphinx/_build/html/index.html`).

Simple example models are found in the `examples/` directory. The README files in these directories explain how to run the models.

## Branches
Updates to RADMC-3D are done carefully. A major new features is typically first offered in a branch (git branch, in github, see the button that says 'master'; if you click on it, you will see other branches).

A branch of particular interest is `camera_openmp`, which is the OpenMP parallellization of the images and spectra, which was implemented by Patrick Sheehan. You are strongly encouraged to try out this branch and test it. If I do not receive bug reports within the coming months, I will merge this into the master branch.

Also a branch of particular interest is `unstruct`, in which unstructured grids are implemented (Delaunay, Voronoi or home-made). If you try this out, please let me know if it works correctly and well, and does not break earlier features of RADMC-3D (backward compatibility). There are likely to be still numerous bugs, so be careful. Once I am confident of this new feature, I will merge it into the master branch.

## News
Bugfix [11.11.2021] in OpenMP parallellization of thermal Monte Carlo affecting models with multiple dust species and `iranfreqmode=1` in the `radmc3d.inp` file when using `setthreads X` with X>1. Please upgrade to version of 11.11.2021 of the master branch.

## Website
For more information, please consult the website http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d
