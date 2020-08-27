import numpy as np
import matplotlib.pyplot as plt

data_orig = np.loadtxt('callindex.out_silD03.txt')
data_ld   = np.loadtxt('../astrosil_laordraine/astrosilicate_laordraine93.lnk')

plt.figure()
plt.plot(data_orig[:,0],data_orig[:,3]+1)
plt.plot(data_ld[:,0],data_ld[:,1])
plt.xscale('log')

plt.figure()
plt.plot(data_orig[:,0],data_orig[:,4])
plt.plot(data_ld[:,0],data_ld[:,2])
plt.xscale('log')
plt.yscale('log')
plt.show()

with open('astrosilicate_draine03.lnk','w') as f:
    f.write('# Optical constants for "astronomical silicate" from Draine (2003)\n')
    f.write('# ApJ 598, 1017. Material is MgFeSiO4.\n')
    f.write('# The material density is 3.3 g/cm^3.\n')
    f.write('# When you use this file, please cite the above paper.\n')
    f.write('# Columns are: lambda [micron], n, k\n')
    for inu in range(len(data_orig[:,0])):
        f.write('{0:13.6e}  {1:13.6e}  {2:13.6e}\n'.format(data_orig[inu,0],data_orig[inu,3]+1,data_orig[inu,4]))

