#!/usr/bin/env python

try:
    from numpy.distutils.core import setup, Extension

    ext = Extension(name = 'radmc3dPy._bhmie',
                    sources = ['radmc3dPy/bhmie.f90'])

except:
    msg = "numpy.distutils.core could not be imported. It is required to build the fast, fortran version of the mie scattering code, bhmie.\n"
    msg = "Falling back to the slower python version"
    print(msg)
    ext = None
    from distutils.core import setup

import os, sys
from subprocess import Popen, PIPE

def findFiles(src_dir, *wildcards):

    src_dir = src_dir.strip()
    while src_dir[-1]=='/':
        src_dir = src_dir[:-1]


    # Find all directory names
    dirList = Popen(['find '+src_dir+' -name "*"'], shell=True, \
            stdout=PIPE).communicate()[0].split()
   

    foundList = []
    for i in range(len(dirList)):
        if os.path.isdir(dirList[i]):
            # Find the appropriate files within each directory
            fileList = []
            for wc in wildcards:
                com = 'ls -1 '+dirList[i].decode('utf-8')+'/'+wc
                dum = Popen([com], shell=True, 
                        stdout=PIPE, stderr=PIPE).communicate()[0].split()

                if len(dum)>0:
                    flist = [s.decode('utf-8') for s in dum]
                    fileList.extend(flist)

            if len(fileList)>0:
                foundList.append((dirList[i], fileList))

    return foundList    
   
fileList = findFiles('./radmc3dPy', './doc/html/', '*.py')

python_files = []
for i in range(len(fileList)):
    for j in range(len(fileList[i][1])):
        python_files.append(fileList[i][1][j])

moduleNames = []
packageNames = []
for i in range(len(python_files)):

    ind1 = python_files[i].strip()[::-1].find('/')
    dum  = python_files[i].strip()[-ind1:-3]

    if dum.strip()!='__init__':
        moduleNames.append('radmc3dPy.'+dum)
    else:
        sdum = python_files[i].split('/')[1:-1]

        txt = sdum[0]
        if len(sdum)>1:
            for imod in range(1,len(sdum)):
                txt += '.'+sdum[imod]

        packageNames.append(str(txt))


if ext is not None:
    setup(name='radmc3dPy',
          version='0.30.2',
          requires=['numpy', 'scipy', 'matplotlib', 'astropy'], 
          description='Python module for RADMC3D',
          author='Attila Juhasz',
          author_email='juhasz@ast.cam.ac.uk',
          packages=packageNames,
          ext_modules=[ext])
else:
    setup(name='radmc3dPy',
          version='0.30.2',
          requires=['numpy', 'scipy', 'matplotlib', 'astropy'], 
          description='Python module for RADMC3D',
          author='Attila Juhasz',
          author_email='juhasz@ast.cam.ac.uk',
          packages=packageNames)


