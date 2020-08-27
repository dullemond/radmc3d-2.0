#
# This Python script demonstrates the use of the simpleread and
# simpleplot tools. They are bare-bones tools for analyzing RADMC-3D
# input/output files. For more sophisticated stuff, use radmc3dPy.
#
# Make sure to have done the following beforhand:
#
#  First compile RADMC-3D
#  Then run:
#   python problem_setup.py
#   radmc3d mctherm
#

from radmc3d_tools.simpleread import *
from radmc3d_tools.simpleplot import *
import os

#
# View a 2-D slice of the 3-D array of the setup
#
d = read_dustdens()
surface(d.rhodust[:,:,16])

#
# View a 2-D slice of the 3-D array of the dust temperature
#
d = read_dusttemp()
surface(d.dusttemp[:,:,16])

#
# Make and plot an example image
#
os.system("radmc3d image incl 60 phi 30 lambda 1000")
im = read_image()
image(im.image[:,:,0].T,range=[0,3e-14])

plt.show()
