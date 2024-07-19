import numpy as np
import os

class subBox(object):
    """
    Class for RADMC-3D's subbox functionality.

    A subbox is simply a regular-gridded sample of some quantity (e.g. the
    dust density) on the true model grid. If the model grid is regularly
    gridded, there is no real point in using a subbox. But if your grid uses
    AMR (octree or patch-based), then it is not easy to handle these data.
    Using subbox you can ask RADMC-3D to map the data on a regular grid 
    (the subbox). The word "subbox" reflects the idea that you can easily
    "zoom in" to any structure by choosing the box size to be much smaller
    than the full grid. 
    """

    def __init__(self):
        self.datatype = ""

    def makeSubbox(self,datatype,bounds,nxyz,phi1=0.,theta=0.,phi2=0.):
        if((datatype=='dust_density') or (datatype=='rhodust')):
            self.datatype       = 'dust_density'
            radmc3d_subbox_type = 'subbox_dust_density'
        elif((datatype=='dust_temperature') or (datatype=='dusttemp') or (datatype=='tdust')):
            self.datatype       = 'dust_temperature'
            radmc3d_subbox_type = 'subbox_dust_temperature'
        elif((datatype=='lines_levelpop') or (datatype=='levelpop')):
            self.datatype       = 'lines_levelpop'
            radmc3d_subbox_type = 'subbox_levelpop'
        else:
            errmsg = "Subbox data type not recognized: "+datatype
            raise ValueError(errmsg)
        if(os.path.isfile('./radmc3d')):
            com = './radmc3d'
        else:
            com = 'radmc3d'
        com += ' '+radmc3d_subbox_type
        com += ' subbox_xyz01 {0} {1} {2} {3} {4} {5}'.format(bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5])
        com += ' subbox_nxyz {0} {1} {2}'.format(nxyz[0],nxyz[1],nxyz[2])
        com += ' subbox_phi1 {0}'.format(phi1)
        com += ' subbox_theta {0}'.format(theta)
        com += ' subbox_phi2 {0}'.format(phi2)
        print(com)
        os.system(com)

    def readSubbox(self,datatype=None,fname=None,spec=None,indexorder='fortran'):
        if datatype is not None:
            if((datatype=='dust_density') or (datatype=='rhodust')):
                self.datatype       = 'dust_density'
            elif((datatype=='dust_temperature') or (datatype=='dusttemp') or (datatype=='tdust')):
                self.datatype       = 'dust_temperature'
            elif((datatype=='lines_levelpop') or (datatype=='levelpop')):
                self.datatype       = 'lines_levelpop'
            else:
                errmsg = "Subbox data type not recognized: "+datatype
                raise ValueError(errmsg)
        if(self.datatype=='dust_density'):
            if fname is None: fname = 'dust_density_subbox.out'
        elif(self.datatype=='dust_temperature'):
            if fname is None: fname = 'dust_temperature_subbox.out'
        elif(self.datatype=='lines_levelpop'):
            if fname is None:
                assert spec is not None, "For level populations subbox: specify spec."
                fname = 'levelpop_'+spec+'_subbox.out'
        else:
            errmsg = "Subbox data type not recognized: "+self.datatype
            raise ValueError(errmsg)
        with open(fname,'r') as f:
            line = f.readline()
            iformat = int(line)
            assert iformat==2, "Subbox file format inconsistent"
            print('Reading '+fname)
            data = np.fromfile(f, count=-1, sep=" ", dtype=np.float64)
        hdr  = np.array(data[:4], dtype=np.int_)
        self.nx = hdr[0]
        self.ny = hdr[1]
        self.nz = hdr[2]
        self.nv = hdr[3]
        data = data[4:]
        hdr  = np.array(data[:6], dtype=np.float64)
        self.x0 = hdr[0]
        self.x1 = hdr[1]
        self.y0 = hdr[2]
        self.y1 = hdr[3]
        self.z0 = hdr[4]
        self.z1 = hdr[5]
        self.xi = np.linspace(self.x0,self.x1,self.nx+1)
        self.yi = np.linspace(self.y0,self.y1,self.ny+1)
        self.zi = np.linspace(self.z0,self.z1,self.nz+1)
        self.x  = 0.5*(self.xi[1:]+self.xi[:-1])
        self.y  = 0.5*(self.yi[1:]+self.yi[:-1])
        self.z  = 0.5*(self.zi[1:]+self.zi[:-1])
        data = data[6:]
        hdr  = np.array(data[:3], dtype=np.float64)
        self.phi1  = hdr[0]
        self.theta = hdr[1]
        self.phi2  = hdr[2]
        data = data[3:]
        hdr  = np.array(data[:self.nv], dtype=np.int_)
        data = data[self.nv:]
        if(self.datatype=='lines_levelpop'):
            self.level_subset = hdr
        if((self.datatype=='dust_density')or(self.datatype=='dust_temperature')):
            self.dust_ispec = hdr
        data = data.reshape((self.nv,self.nz,self.ny,self.nx))
        if indexorder=='fortran':
            data = np.swapaxes(data, 0, 3)
            data = np.swapaxes(data, 1, 2)
        self.data = data
        if(self.datatype=='dust_density'): self.rhodust  = self.data
        if(self.datatype=='dust_temperature'): self.dusttemp = self.data
        if(self.datatype=='lines_levelpop'): self.levelpop = self.data
