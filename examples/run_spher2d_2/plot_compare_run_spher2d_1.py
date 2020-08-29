from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import math
import os

from radmc3dPy.analyze import *  
from radmc3dPy.image import *  
from radmc3dPy.natconst import * 

os.chdir('../run_spher2d_1/')
sed_basic   = readSpectrum()
os.chdir('../run_spher2d_2/')
sed_refined = readSpectrum()

lam = sed_basic[:,0]
nu  = 1e4*cc/lam
fnu_basic   = sed_basic[:,1]
fnu_refined = sed_refined[:,1]

plt.figure()
plt.loglog(lam,nu*fnu_basic,label='Basic grid')
plt.loglog(lam,nu*fnu_refined,label='Refined grid')
plt.title('Effect of grid refinement in the SED')
plt.ylim((1e-9,1e-4))
plt.xlabel(r'$\lambda\,[\mu\mathrm{m}]$')
plt.ylabel(r'$\nu F_\nu\,[\mathrm{erg}\,\mathrm{cm}^{-1}\,\mathrm{s}^{-1}]$')
plt.legend()
