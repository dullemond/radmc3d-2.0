"""
SIMPLEREAD: A simple Python reading tool for RADMC-3D

By C.P. Dullemond and A. Juhasz (2020)

This Python script contains a few very basic functions for reading the input
and output files of RADMC-3D. They are mostly meant for pedagogical purpose: to
get a feeling for the structure of these files, and for those who just started
with RADMC-3D. These functions only read the ascii versions of these files, not
the binary versions. Also, they are very incomplete, and work only for some very
simple models (e.g. no octree grids, only simple dust opacities, no line data
files).

Usage: In a model directory, in Python (use Python 3):

  from radmc3d_tools.simpleread import *

You can then read, for instance, the dust temperature computed by RADMC-3D
using the radmc3d mctherm command:

  d = read_dusttemp()

Then d contains d.grid and d.dusttemp, the latter being the array of dust
temperatures. There are various read_***() functions: have a look in this
script file.

For real scientific use, it is recommended to use the radmc3dPy package, which is
also included in the RADMC-3D distribution. It is located here (starting from the
root of the RADMC-3D package):

  python/radmc3dPy/

"""
from __future__ import print_function

try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use the python module of RADMC-3D you need to install Numpy')
    print(traceback.format_exc())

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None
    print('Warning')
    print('matplotlib.pyplot cannot be imported')
    print('Without matplotlib you can use the python module to set up a model but you will not be able to plot things')
    print('or display images')

import glob

class simplereaddataobject(object):
    """
    Generic data object for the RADMC-3D simpleread.py functions.
    """
    def __init__(self,datatype):
        self.datatype = datatype

def read_grid():
    """
    Reading the amr_grid.inp file, but only for regular grids (not for octree ones). 

    ARGUMENTS:
      None

    RETURNS:
      Data object containing:

        .nx          Nr of cells in x direction (in spherical coordinates this is the r-direction)
        .ny          Nr of cells in y direction (in spherical coordinates this is the theta-direction)
        .nz          Nr of cells in z direction (in spherical coordinates this is the phi-direction)
        .nxi etc     The same as .nx etc, but now +1, giving the number of cell interfaces
        .crd_sys     String indicating the coordinate system ('car' or 'sph')
        .x,.y,.z     The x, y and z grid cell center locations
        .xi,.yi,.zi  The x, y and z grid cell interface locations (each array 1 longer than the x, y, z ones)

    """
    grid  = simplereaddataobject('grid')
    fname = 'amr_grid.inp'
    print('Reading '+fname)
    data  = np.fromfile(fname, count=-1, sep=" ", dtype=np.float64)

    # Read the header
    hdr   = np.array(data[:10], dtype=np.int)
    data  = data[10:]

    # Check the file format
    if hdr[0] != 1:
        msg = 'Unknown format number in amr_grid.inp'
        raise RuntimeError(msg)
    if hdr[2] < 100:
        grid.crd_sys = 'car'
    elif (hdr[2] >= 100) & (hdr[2] < 200):
        grid.crd_sys = 'sph'
    else:
        raise ValueError('Unsupported coordinate system identification in the ' + fname + ' file.')

    # Get active dimensions
    grid.act_dim = hdr[4:7]

    # Get the number of cells in each dimensions
    grid.nx = hdr[7]
    grid.ny = hdr[8]
    grid.nz = hdr[9]
    grid.nxi, grid.nyi, grid.nzi = grid.nx + 1, grid.ny + 1, grid.nz + 1

    # Get the cell interfaces
    grid.xi = data[:grid.nxi]
    data = data[grid.nxi:]
    grid.yi = data[:grid.nyi]
    data = data[grid.nyi:]
    grid.zi = data[:grid.nzi]

    # Compute the cell centers
    if grid.crd_sys == 'car':
        grid.x = (grid.xi[0:grid.nx] + grid.xi[1:grid.nx + 1]) * 0.5
        grid.y = (grid.yi[0:grid.ny] + grid.yi[1:grid.ny + 1]) * 0.5
        grid.z = (grid.zi[0:grid.nz] + grid.zi[1:grid.nz + 1]) * 0.5
    else:
        grid.x = np.sqrt(grid.xi[0:grid.nx] * grid.xi[1:grid.nx + 1])
        grid.y = (grid.yi[0:grid.ny] + grid.yi[1:grid.ny + 1]) * 0.5
        grid.z = (grid.zi[0:grid.nz] + grid.zi[1:grid.nz + 1]) * 0.5

    # Now return the grid object
    return grid
    
def read_dustdens(indexorder='fortran'):
    """
    Reading the dust_density.inp file, but only for regular grids, and only
    for text data format (not binary).

    ARGUMENTS:
      indexorder        If 'fortran' then converting array to fortran 
                        index order (default). Else use Python/C order.

    RETURNS:
      Data object containing:

        .grid           A grid object (see read_grid())
        .rhodust        An array with the dust density values (in g/cm^3)
    """
    grid      = read_grid()
    dustdens  = simplereaddataobject('dust_density')
    dustdens.grid = grid
    fname     = 'dust_density.inp'
    print('Reading '+fname)
    data      = np.fromfile(fname, count=-1, sep=" ", dtype=np.float64)
    
    # Read the header
    hdr       = np.array(data[:3], dtype=np.int)
    data      = data[3:]

    # Check the file format
    if hdr[0] != 1:
        msg = 'Unknown format number in dust_density.inp'
        raise RuntimeError(msg)

    # Get the number of cells, and check against
    nrcells = grid.nx*grid.ny*grid.nz
    if(hdr[1]!=nrcells):
        msg = 'Number of grid cells in dust_density.inp inconsistent with amr_grid.inp'
        raise RuntimeError(msg)

    # Get the number of dust species
    dustdens.nrspec = hdr[2]

    # Convert the rest of the data to the proper shape
    data = np.reshape(data, [dustdens.nrspec, grid.nz, grid.ny, grid.nx])

    # If indexorder is set to 'fortran', then the inner index of the array
    # should be left (even though in Python the inner index is right). This
    # is to assure that the index order in the Python arrays is the same as
    # in the RADMC-3D code. But by setting indexorder to anything else, you
    # can keep Python natural order (which is equal to C index order), in
    # which the inner index is the rightmost index.
    if indexorder=='fortran':
        data = np.swapaxes(data, 0, 3)
        data = np.swapaxes(data, 1, 2)

    # Now add this to the object
    dustdens.rhodust = data

    # Return the dustdens object
    return dustdens

def read_dusttemp(indexorder='fortran'):
    """
    Reading the dust_temperature.dat file, but only for regular grids, and only
    for text data format (not binary).

    ARGUMENTS:
      indexorder        If 'fortran' then converting array to fortran 
                        index order (default). Else use Python/C order.

    RETURNS:
      Data object containing:

        .grid           A grid object (see read_grid())
        .dusttemp       An array with the dust temperature values (in K)
    """
    grid      = read_grid()
    dusttemp  = simplereaddataobject('dust_temperature')
    dusttemp.grid = grid
    fname     = 'dust_temperature.dat'
    print('Reading '+fname)
    data      = np.fromfile(fname, count=-1, sep=" ", dtype=np.float64)
    
    # Read the header
    hdr       = np.array(data[:3], dtype=np.int)
    data      = data[3:]

    # Check the file format
    if hdr[0] != 1:
        msg = 'Unknown format number in dust_temperature.inp'
        raise RuntimeError(msg)

    # Get the number of cells, and check against
    nrcells = grid.nx*grid.ny*grid.nz
    if(hdr[1]!=nrcells):
        msg = 'Number of grid cells in dust_temperature.inp inconsistent with amr_grid.inp'
        raise RuntimeError(msg)

    # Get the number of dust species
    dusttemp.nrspec = hdr[2]

    # Convert the rest of the data to the proper shape
    data = np.reshape(data, [dusttemp.nrspec, grid.nz, grid.ny, grid.nx])

    # If indexorder is set to 'fortran', then the inner index of the array
    # should be left (even though in Python the inner index is right). This
    # is to assure that the index order in the Python arrays is the same as
    # in the RADMC-3D code. But by setting indexorder to anything else, you
    # can keep Python natural order (which is equal to C index order), in
    # which the inner index is the rightmost index.
    if indexorder=='fortran':
        data = np.swapaxes(data, 0, 3)
        data = np.swapaxes(data, 1, 2)

    # Now add this to the object
    dusttemp.dusttemp = data

    # Return the dusttemp object
    return dusttemp

def read_image(indexorder='fortran'):
    """
    Reading the image.out file.

    ARGUMENTS:
      indexorder        If 'fortran' then converting array to fortran 
                        index order (default). Else use Python/C order.

    RETURNS:
      Data object containing:

        .freq           Frequency at which the image is taken
        .image          An array with the image intensity in erg/(s.cm^2.Hz.ster)
    """
    pc        = 3.08572e18     # Parsec                  [cm]
    cc        = 2.99792458e10  # Light speed             [cm/s]
    image     = simplereaddataobject('image')
    fname     = 'image.out'
    print('Reading '+ fname)
    with open(fname, 'r') as rfile:
        dum = ''

        # Format number
        iformat = int(rfile.readline())

        # Nr of pixels
        dum = rfile.readline()
        dum = dum.split()
        image.nx = int(dum[0])
        image.ny = int(dum[1])

        # Nr of frequencies
        image.nfreq = int(rfile.readline())
        image.nwav = image.nfreq

        # Pixel sizes
        dum = rfile.readline()
        dum = dum.split()
        image.sizepix_x = float(dum[0])
        image.sizepix_y = float(dum[1])

        # Wavelength of the image
        image.wav = np.zeros(image.nwav, dtype=np.float64)
        for iwav in range(image.nwav):
            image.wav[iwav] = float(rfile.readline())
        image.freq = cc / image.wav * 1e4

        # Read the rest of the data
        data = np.fromfile(rfile, count=-1, sep=" ", dtype=np.float64)
        
    # Convert the rest of the data to the proper shape
    if iformat == 1:
        # We have a normal total intensity image
        image.stokes = False
        data = np.reshape(data, [image.nfreq, image.ny, image.nx])
        if indexorder=='fortran':
            data = np.swapaxes(data, 0, 2)

    elif iformat == 3:
        # We have the full stokes image
        image.stokes = True
        data = np.reshape(data, [image.nfreq, 4, image.ny, image.nx])
        if indexorder=='fortran':
            data = np.swapaxes(data, 0, 3)
            data = np.swapaxes(data, 1, 2)

    else:
        msg = 'Unknown format number in image.out'
        raise ValueError(msg)

    # Add this to the object
    image.image = data
    
    # Conversion from erg/s/cm/cm/Hz/ster to Jy/pixel
    conv = image.sizepix_x * image.sizepix_y / pc**2. * 1e23
    image.imageJyppix = image.image * conv
    
    # Create the x and y axes in units of cm
    image.x = ((np.arange(image.nx, dtype=np.float64) + 0.5) - image.nx / 2) * image.sizepix_x
    image.y = ((np.arange(image.ny, dtype=np.float64) + 0.5) - image.ny / 2) * image.sizepix_y
    
    # Return object
    return image

def read_spectrum(dpc=1.):
    """
    Reading the spectrum.out file.

    ARGUMENTS:
      dpc               Distance of observer in parsec (default=1)

    RETURNS:
      Data object containing:

        .wav            Wavelength array of the spectrum in micron
        .freq           Frequency array of the spectrum in Hertz
        .fnu            An array with the flux F_nu at dpc parsec in erg/(s.cm^2.Hz)
    """
    cc        = 2.99792458e10  # Light speed             [cm/s]
    spectrum  = simplereaddataobject('spectrum')
    fname     = 'spectrum.out'
    print('Reading '+ fname)
    with open(fname, 'r') as rfile:
        # Read the format number
        iformat = int(rfile.readline())

        # Read the number of wavelengths
        spectrum.nwav = int(rfile.readline())
        spectrum.freq = spectrum.nwav

        # Read a blank line
        dum = rfile.readline()

        # Read the rest of the data
        data = np.fromfile(rfile, count=-1, sep=" ", dtype=np.float64)

    # Reshape the spectrum
    data          = np.reshape(data, [spectrum.nwav, 2])
    spectrum.wav  = data[:,0]
    spectrum.freq = 1e4 * cc / spectrum.wav
    spectrum.fnu  = data[:,1]

    # Rescale spectrum to other distance
    spectrum.fnu /= dpc**2

    # Return the spectrum
    return spectrum

def read_dustkappa(species=None):
    """
    Reading a dust opacity file (but only the basic dustkappa_*.inp type,
    not the dustkapscatmat_*.inp type).

    ARGUMENTS:
      species           The dust species: Reading dustkappa_<species>.inp
                        If unspecified, read_dustkappa will search for such a file.
                        If it finds a single one, it will read that. Otherwise it
                        will request you to specify species.

    RETURNS:
      Data object containing:

        .wav            Wavelength array of the spectrum in micron
        .freq           Frequency array of the spectrum in Hertz
        .kappa_abs      Absorption opacity in cm^2/gram-of-dust
        .kappa_sca      Scattering opacity in cm^2/gram-of-dust
        .kappa_g        The g-coefficient (between -1 and 1) for non-isotropic scattering
    """
    cc        = 2.99792458e10  # Light speed             [cm/s]

    # Find which dust opacity to read
    if species is None:
        fnames = glob.glob('dustkappa_*.inp')
        if len(fnames)==0:
            msg = 'No file of type dustkappa_*.inp is found in this directory.'
            raise RuntimeError(msg)
        if len(fnames)>1:
            msg = 'More than one file of type dustkappa_*.inp is found in this directory. Please specify the name of the dust species as keyword species.'
            raise RuntimeError(msg)
        species = fnames[0]
    if species[0:10]=='dustkappa_': species = species[10:]
    if species[-4:]=='.inp': species = species[:-4]

    # Read that dust opacity
    dustkappa = simplereaddataobject('dustkappa')
    fname     = 'dustkappa_'+species+'.inp'
    print('Reading '+ fname)
    with open(fname, 'r') as rfile:
        # Check the file format (skipping comments)
        iformat_str = rfile.readline()
        while iformat_str[0]=='#':
            iformat_str = rfile.readline()
        iformat = int(iformat_str)

        # Read the number of wavelength points
        dustkappa.nwav = int(rfile.readline())
        dustkappa.freq = dustkappa.nwav

        # Now read the rest of the data
        data = np.fromfile(rfile, count=-1, sep=" ", dtype=np.float64)

    # Reshape the data
    if iformat==1:
        data = np.reshape(data, [dustkappa.nwav, 2])
    elif iformat==2:
        data = np.reshape(data, [dustkappa.nwav, 3])
    elif iformat==3:
        data = np.reshape(data, [dustkappa.nwav, 4])
    else:
        msg = 'Format number of kappa file not known'
        raise RuntimeError(msg)

    # Extract the information
    dustkappa.wav        = data[:,0]
    dustkappa.freq       = 1e4 * cc / dustkappa.wav
    dustkappa.kappa_abs  = data[:,1]
    if iformat>1:
        dustkappa.kappa_sca = data[:,2]
    if iformat>2:
        dustkappa.kappa_g = data[:,3]

    # Return dustkappa
    return dustkappa

def read_gastemp(indexorder='fortran'):
    """
    Reading the gas_temperature.inp file, but only for regular grids, and only
    for text data format (not binary).

    ARGUMENTS:
      indexorder        If 'fortran' then converting array to fortran
                        index order (default). Else use Python/C order.

    RETURNS:
      Data object containing:

        .grid           A grid object (see read_grid())
        .gastemp        An array with the gas temperature values (in K)
    """
    grid      = read_grid()
    gastemp   = simplereaddataobject('gas_temperature')
    gastemp.grid = grid
    fname     = 'gas_temperature.inp'
    print('Reading '+fname)
    data      = np.fromfile(fname, count=-1, sep=" ", dtype=np.float64)

    # Read the header
    hdr       = np.array(data[:2], dtype=np.int)
    data      = data[2:]

    # Check the file format
    if hdr[0] != 1:
        msg = 'Unknown format number in gas_temperature.inp'
        raise RuntimeError(msg)

    # Get the number of cells, and check against
    nrcells = grid.nx*grid.ny*grid.nz
    if(hdr[1]!=nrcells):
        msg = 'Number of grid cells in gas_temperature.inp inconsistent with amr_grid.inp'
        raise RuntimeError(msg)

    # Convert the rest of the data to the proper shape
    data = np.reshape(data, [grid.nz, grid.ny, grid.nx])

    # If indexorder is set to 'fortran', then the inner index of the array
    # should be left (even though in Python the inner index is right). This
    # is to assure that the index order in the Python arrays is the same as
    # in the RADMC-3D code. But by setting indexorder to anything else, you
    # can keep Python natural order (which is equal to C index order), in
    # which the inner index is the rightmost index.
    if indexorder=='fortran':
        data = np.swapaxes(data, 0, 2)

    # Now add this to the object
    gastemp.gastemp = data

    # Return the gastemp object
    return gastemp

def read_gasvelocity(indexorder='fortran'):
    """
    Reading the gas_velocity.inp file, but only for regular grids, and only
    for text data format (not binary).

    ARGUMENTS:
      indexorder        If 'fortran' then converting array to fortran
                        index order (default). Else use Python/C order.

    RETURNS:
      Data object containing:

        .grid           A grid object (see read_grid())
        .gasvelo        An array with the gas velocity values (in cm/s)
    """
    grid      = read_grid()
    gasvelo   = simplereaddataobject('gas_velocity')
    gasvelo.grid = grid
    fname     = 'gas_velocity.inp'
    print('Reading '+fname)
    data      = np.fromfile(fname, count=-1, sep=" ", dtype=np.float64)

    # Read the header
    hdr       = np.array(data[:2], dtype=np.int)
    data      = data[2:]

    # Check the file format
    if hdr[0] != 1:
        msg = 'Unknown format number in gas_velocity.inp'
        raise RuntimeError(msg)

    # Get the number of cells, and check against
    nrcells = grid.nx*grid.ny*grid.nz
    if(hdr[1]!=nrcells):
        msg = 'Number of grid cells in gas_velocity.inp inconsistent with amr_grid.inp'
        raise RuntimeError(msg)

    # Convert the rest of the data to the proper shape
    data = np.reshape(data, [grid.nz, grid.ny, grid.nx, 3])

    # If indexorder is set to 'fortran', then the inner index of the array
    # should be left (even though in Python the inner index is right). This
    # is to assure that the index order in the Python arrays is the same as
    # in the RADMC-3D code. But by setting indexorder to anything else, you
    # can keep Python natural order (which is equal to C index order), in
    # which the inner index is the rightmost index.
    if indexorder=='fortran':
        data = np.swapaxes(data, 0, 3)
        data = np.swapaxes(data, 1, 2)

    # Now add this to the object
    gasvelo.velocity = data

    # Return the gasvelo object
    return gasvelo

def read_molnumdens(molecule,indexorder='fortran'):
    """
    Reading a numberdens_xxxx.inp file, but only for regular grids, and only
    for text data format (not binary).

    ARGUMENTS:
      molecule          Then name of the molecule (e.g. 'co')
      indexorder        If 'fortran' then converting array to fortran
                        index order (default). Else use Python/C order.

    RETURNS:
      Data object containing:

        .grid           A grid object (see read_grid())
        .numdens        An array with the number density values of that molecule (in 1/cm^3)
    """
    grid      = read_grid()
    molnumdens   = simplereaddataobject('molecule_number_density')
    molnumdens.grid = grid
    fname     = 'numberdens_'+molecule+'.inp'
    print('Reading '+fname)
    data      = np.fromfile(fname, count=-1, sep=" ", dtype=np.float64)

    # Read the header
    hdr       = np.array(data[:2], dtype=np.int)
    data      = data[2:]

    # Check the file format
    if hdr[0] != 1:
        msg = 'Unknown format number in '+fname
        raise RuntimeError(msg)

    # Get the number of cells, and check against
    nrcells = grid.nx*grid.ny*grid.nz
    if(hdr[1]!=nrcells):
        msg = 'Number of grid cells in '+fname+' inconsistent with amr_grid.inp'
        raise RuntimeError(msg)

    # Convert the rest of the data to the proper shape
    data = np.reshape(data, [grid.nz, grid.ny, grid.nx])

    # If indexorder is set to 'fortran', then the inner index of the array
    # should be left (even though in Python the inner index is right). This
    # is to assure that the index order in the Python arrays is the same as
    # in the RADMC-3D code. But by setting indexorder to anything else, you
    # can keep Python natural order (which is equal to C index order), in
    # which the inner index is the rightmost index.
    if indexorder=='fortran':
        data = np.swapaxes(data, 0, 2)

    # Now add this to the object
    molnumdens.numdens = data

    # Return the molnumdens object
    return molnumdens

def read_mollevelpop(molecule,indexorder='fortran'):
    """
    Reading a levelpop_xxxx.dat file, but only for regular grids, and only
    for text data format (not binary).

    ARGUMENTS:
      molecule          Then name of the molecule (e.g. 'co')
      indexorder        If 'fortran' then converting array to fortran
                        index order (default). Else use Python/C order.

    RETURNS:
      Data object containing:

        .grid           A grid object (see read_grid())
        .pop            An array with the level population values of that molecule (in units of number density 1/cm^3)
    """
    grid      = read_grid()
    mollevelpop   = simplereaddataobject('molecule_level_populations')
    mollevelpop.grid = grid
    fname     = 'levelpop_'+molecule+'.dat'
    print('Reading '+fname)
    data      = np.fromfile(fname, count=-1, sep=" ", dtype=np.float64)

    # Read the header
    hdr       = np.array(data[:3], dtype=np.int)
    data      = data[3:]

    # Check the file format
    if hdr[0] != 1:
        msg = 'Unknown format number in '+fname
        raise RuntimeError(msg)

    # Get the number of cells, and check against
    nrcells = grid.nx*grid.ny*grid.nz
    if(hdr[1]!=nrcells):
        msg = 'Number of grid cells in '+fname+' inconsistent with amr_grid.inp'
        raise RuntimeError(msg)

    # Number of levels
    mollevelpop.nlevels = hdr[2]

    # Identities of the levels (for associating them to the levels in the molecule_xxxx.inp file)
    hdr       = np.array(data[:mollevelpop.nlevels], dtype=np.int)
    data      = data[mollevelpop.nlevels:]
    mollevelpop.levels = hdr

    # Convert the rest of the data to the proper shape
    data = np.reshape(data, [grid.nz, grid.ny, grid.nx, mollevelpop.nlevels])

    # If indexorder is set to 'fortran', then the inner index of the array
    # should be left (even though in Python the inner index is right). This
    # is to assure that the index order in the Python arrays is the same as
    # in the RADMC-3D code. But by setting indexorder to anything else, you
    # can keep Python natural order (which is equal to C index order), in
    # which the inner index is the rightmost index.
    if indexorder=='fortran':
        data = np.swapaxes(data, 0, 3)
        data = np.swapaxes(data, 1, 2)

    # Now add this to the object
    mollevelpop.pop = data

    # For convenience, also compute the relative level populations
    # (which sum up to 1 in each cell)
    mollevelpop.relpop = np.zeros_like(mollevelpop.pop)
    if indexorder=='fortran':
        dum = mollevelpop.pop.sum(axis=0)
        for i in range(mollevelpop.nlevels):
            mollevelpop.relpop[i,:,:,:] = mollevelpop.pop[i,:,:,:]/dum[:,:,:]
    else:
        dum = mollevelpop.pop.sum(axis=3)
        for i in range(mollevelpop.nlevels):
            mollevelpop.relpop[:,:,:,i] = mollevelpop.pop[:,:,:,i]/dum[:,:,:]

    # Return the mollevelpop object
    return mollevelpop

def read_spy_image(indexorder='fortran',usenan=True):
    """
    Reading the spy_image.out file. A spy image is a 3D data block
    where for each pixel of an image the intermediate intensities 
    along the entire ray belonging to that pixel are stored as well.
    The ray length typically varies, because each point is a crossing
    of the ray with a grid wall. That varies from ray to ray. The
    data block has a maximum ray length that should (hopefully) be
    long enough so that all rays have a length (in number of points)
    less than that. That also means that the spy image has a lot of
    empty points (because the full array, including the non-used
    points, is stored). These are marked with an intensity value 
    of -1e90. To make it easier to handle using numpy as matplotlib,
    the read_spy_image() function will put all values of these 
    empty/not-used points to np.nan ('not a number').

    ARGUMENTS:
      indexorder        If 'fortran' then converting array to fortran 
                        index order (default). Else use Python/C order.
      usenan            If set to True, then all not-used datapoints are
                        set to np.nan. If False, then the not-used datapoints
                        are recognized by a value of -1e90 for the intensity.

    RETURNS:
      Data object containing:

        .freq           Frequency at which the image is taken
        .intensity      An array with the intensity at each point along
                        the ray of each pixel in erg/(s.cm^2.Hz.ster)
        .x              The x location of each point at which the .intensity is given
        .y              The y location of each point at which the .intensity is given
        .z              The z location of each point at which the .intensity is given

    EXAMPLE:

      radmc3d image lambda 10 incl 70 phi 30 nofluxcons spymode
      ipython --matplotlib
      from simpleread import *
      a = read_spy_image()
      ix_pixel = a.nx//2   # The middle of the image
      iy_pixel = a.ny//2   # The middle of the image
      plt.figure()
      plt.semilogy(a.z[1:,ix_pixel,iy_pixel],a.intensity[1:,ix_pixel,iy_pixel],'.-')
      plt.ylim(ymin=1e-6*a.intensity_max)
      plt.xlabel('z coordinate along ray [cm]')
      plt.ylabel(r'$I_\nu\;[\mathrm{erg}\,\mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{Hz}^{-1}\mathrm{ster}^{-1}]$')
      plt.show()

    """
    pc        = 3.08572e18     # Parsec                  [cm]
    cc        = 2.99792458e10  # Light speed             [cm/s]
    image     = simplereaddataobject('spy_image')
    fname     = 'spy_image.out'
    print('Reading '+ fname)
    with open(fname, 'r') as rfile:
        dum = ''

        # Format number
        iformat = int(rfile.readline())

        # Nr of pixels
        dum = rfile.readline()
        dum = dum.split()
        image.nx = int(dum[0])
        image.ny = int(dum[1])

        # Max nr of points along each ray
        image.raylen = int(rfile.readline())

        # Nr of wavelength is always 1 for spy images
        image.nfreq = 1
        image.nwav = 1

        # Pixel sizes
        dum = rfile.readline()
        dum = dum.split()
        image.sizepix_x = float(dum[0])
        image.sizepix_y = float(dum[1])

        # Wavelength of the image
        image.wav = float(rfile.readline())
        image.freq = cc / image.wav * 1e4

        # Read the rest of the data
        data = np.fromfile(rfile, count=-1, sep=" ", dtype=np.float64)
        
    # Convert the rest of the data to the proper shape
    if iformat == 1:
        # We have a normal total intensity image
        image.stokes = False
        data = np.reshape(data, [image.ny, image.nx, image.raylen, 4])
        image.x = data[:,:,:,0].squeeze()
        image.y = data[:,:,:,1].squeeze()
        image.z = data[:,:,:,2].squeeze()
        image.intensity = data[:,:,:,3].squeeze()
        if indexorder=='fortran':
            image.x = np.swapaxes(image.x, 0, 2)
            image.y = np.swapaxes(image.y, 0, 2)
            image.z = np.swapaxes(image.z, 0, 2)
            image.intensity = np.swapaxes(image.intensity, 0, 2)

    elif iformat == 3:
        # We have the full stokes image
        image.stokes = True
        data = np.reshape(data, [image.ny, image.nx, image.raylen, 7])
        image.x = data[:,:,:,0].squeeze()
        image.y = data[:,:,:,1].squeeze()
        image.z = data[:,:,:,2].squeeze()
        image.intensity = data[:,:,:,3:6]
        if indexorder=='fortran':
            image.x = np.swapaxes(image.x, 0, 3)
            image.y = np.swapaxes(image.y, 0, 3)
            image.z = np.swapaxes(image.z, 0, 3)
            image.x = np.swapaxes(image.x, 1, 2)
            image.y = np.swapaxes(image.y, 1, 2)
            image.z = np.swapaxes(image.z, 1, 2)
            image.intensity = np.swapaxes(image.intensity, 0, 3)
            image.intensity = np.swapaxes(image.intensity, 1, 2)
    else:
        msg = 'Unknown format number in spy_image.out'
        raise ValueError(msg)

    # Determine the maximum value
    image.intensity_max = image.intensity.max()

    # Set non-used points to np.nan
    if usenan:
        mask=image.intensity<0
        image.intensity[mask]=np.nan
        image.x[mask]=np.nan
        image.y[mask]=np.nan
        image.z[mask]=np.nan
    
    # Create the image plane x and y axes in units of cm
    image.image_x = ((np.arange(image.nx, dtype=np.float64) + 0.5) - image.nx / 2) * image.sizepix_x
    image.image_y = ((np.arange(image.ny, dtype=np.float64) + 0.5) - image.ny / 2) * image.sizepix_y
    
    # Return object
    return image

