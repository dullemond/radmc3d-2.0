import matplotlib.pyplot as plt
import numpy as np
from natconst import *
from viewarr import *
from radmc3d_tools.simpleread import *
import os

os.system('radmc3d subbox_dust_density subbox_xyz01 -1.5e14 1.5e14 -1.5e14 1.5e14 -1.5e14 1.5e14 subbox_nxyz 100 100 100')
q = read_subbox(name='dust_density',indexorder='fortran')

#os.system('radmc3d subbox_dust_temperature subbox_xyz01 -1.5e14 1.5e14 -1.5e14 1.5e14 -1.5e14 1.5e14 subbox_nxyz 100 100 100')
#q = read_subbox(name='dust_temperature',indexorder='fortran')

slicearr(q.data)
