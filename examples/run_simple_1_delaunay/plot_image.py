import matplotlib.pyplot as plt
import numpy as np
from natconst import *
from viewarr import *
from radmc3dPy.image import *
import os

os.system('radmc3d image lambda 1000 incl 60 phi 30 nostar nofluxcons')

im = readImage()

plotImage(im)
