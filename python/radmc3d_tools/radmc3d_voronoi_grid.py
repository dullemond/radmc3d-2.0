import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay,Voronoi,ConvexHull
import struct
from tqdm import tqdm    # For a nice progress bar

class Voronoigrid(object):
    def __init__(self,points,bbox=None):
        """
        Create a Voronoi grid for RADMC-3D from a set of 3D points.

        Arguments:
          points      The points given as a numpy array points[npnts,3].
        
        Options:
          bbox               If set to [[xmin,xmax],[ymin,ymax],[zmin,zmax]] this sets
                             the bounding box: only cells that fit entirely inside this
                             box are considered actual cells (with finite volume). The
                             cells that have >=1 point outside the box are treated as
                             "open regions", see below.

        Note on "open regions": These are "cells" that are not closed, and thus 
        have an infinite volume. These regions are still useful for guiding the
        ray of RADMC-3D to the right cells, but they are considered to be empty.
        If you use bbox, then, mathematically speaking, one could close these
        cells (as well as the cells that do not fully fit into the box) at the
        box walls. But the additional complexity is not worth it, so instead, 
        these cells are considered empty.
        """
        self.save_cellcenters = True    # Must be!
        self.save_size        = False   # For now
        #
        # The Qhull algorithms are sensitive to very large numbers, so let's
        # scale everything. First compute a useful scale factor
        #
        scale            = points.max(axis=0)-points.min(axis=0)
        scale            = scale.max()
        self.scale       = scale
        #
        # Create the Voronoi tesselation
        #
        print('Running scipy.spatial.Voronoi() to create Voronoi tesselation...')
        self.vor         = Voronoi(points/scale)
        self.cell_points = self.vor.points*scale
        self.vertices    = self.vor.vertices*scale
        self.npnts       = points.shape[0]
        self.nverts      = 0       # For now
        #
        # Compute the cell volumes. This also automatically notices
        # which points correspond to "real" cells (with a finite volume)
        # and which points correspond to "open" regions going out to
        # infinity (with an infinite volume). These "open" regions are
        # marked by a volume = 0.0
        #
        print('Computing the cell volumes...')
        self.compute_cell_volumes()
        #
        # If bbox is given, then reject all cells that have Voronoi vertices
        # that are outside of the box. 
        #
        self.bbox = bbox
        if bbox is not None:
            print('Implementing bounding box by rejecting cells that are not entirely in the bbox...')
            self.bbox_reject_cells()
        #
        # Now create the cell walls
        #
        print('Creating the cell walls...')
        self.create_cellwalls()
        print('Linking the cell walls back to the cells...')
        self.link_walls_to_cells()
        #
        # Compute some diagnostics
        #
        self.compute_diagnostics()

    def compute_cell_volumes(self):
        """
        Compute the volume of the Voronoi regions.
        From https://stackoverflow.com/questions/19634993/volume-of-voronoi-cell-python
        Answer by ybeltukov
        """
        scale = self.scale
        self.cell_volumes = np.zeros(self.vor.npoints)
        #for i in range(self.npnts):  # Without progress bar
        for i in tqdm(range(self.npnts)):
            indices = self.vor.regions[self.vor.point_region[i]]
            if -1 in indices: # some regions can be open
                self.cell_volumes[i] = 0.0
            else:
                self.cell_volumes[i] = scale**3 * ConvexHull(self.vertices[indices]/scale).volume

    def bbox_reject_cells(self):
        """
        Reject cells that have one or more Voronoi vertices outside of the bbox.
        This is done by setting its cell_volume to 0.0.
        NOTE: Must call this after compute_cell_volumes, because that overwrites
        these.
        """
        vor   = self.vor
        bbox  = self.bbox
        #
        # Find the Voronoi vertices that are outside of the bbox
        #
        maskx = np.logical_or(self.vertices[:,0]<bbox[0][0],self.vertices[:,0]>bbox[0][1])
        masky = np.logical_or(self.vertices[:,1]<bbox[1][0],self.vertices[:,1]>bbox[1][1])
        maskz = np.logical_or(self.vertices[:,2]<bbox[2][0],self.vertices[:,2]>bbox[2][1])
        mask  = np.logical_or(maskx,np.logical_or(masky,maskz))
        ivout = set(np.where(mask)[0])
        #
        # Mark all regions that contain one of the above vertices as volume 0.0
        #
        #for i in tqdm(range(self.npnts)):
        for i in range(self.npnts):
            indices = set(self.vor.regions[self.vor.point_region[i]])
            if len(indices&ivout)>0:
                self.cell_volumes[i] = 0.0
        
    def create_cellwalls(self):
        points           = self.cell_points
        self.nwalls      = self.vor.ridge_points.shape[0]
        #
        # The support vectors
        #
        self.wall_s      = 0.5*(points[self.vor.ridge_points[:,1],:] + points[self.vor.ridge_points[:,0],:])
        #
        # The normals
        #
        self.wall_n      = points[self.vor.ridge_points[:,1],:] - points[self.vor.ridge_points[:,0],:]
        l = (self.wall_n**2).sum(axis=1)**0.5
        self.wall_n     /= l[:,None]
        self.wall_icells = self.vor.ridge_points

    def link_walls_to_cells(self):
        self.cell_iwalls = [ [] for _ in range(self.npnts) ]
        self.cell_wsign  = [ [] for _ in range(self.npnts) ]
        for i in range(self.nwalls):
            icell = self.wall_icells[i,0]
            self.cell_iwalls[icell].append(i)
            self.cell_wsign[icell].append(1)
            icell = self.wall_icells[i,1]
            self.cell_iwalls[icell].append(i)
            self.cell_wsign[icell].append(-1)

    def compute_diagnostics(self):
        self.ncells        = self.npnts
        self.cell_iopen    = np.where((self.cell_volumes<=0.0))[0]
        self.ncells_open   = len(self.cell_iopen)
        self.ncells_closed = len(np.where((self.cell_volumes>0.0))[0])
        assert self.ncells_open+self.ncells_closed==self.ncells, 'Internal error.'
        self.volume_total   = self.cell_volumes.sum()
        l=np.array([len(self.cell_iwalls[i]) for i in range(self.npnts)])
        self.cell_max_nr_walls = l.max()
        self.cell_max_nr_verts = 0
        self.wall_max_nr_verts = 0
        self.vert_max_nr_cells = 0
        self.hull_nwalls       = 0

    def visualize_cells(self,icells=None,alpha=0.9,colors="C1",bbox=None):
        import mpl_toolkits.mplot3d as a3
        if icells is None: icells=np.arange(grid.ncells)
        vor    = self.vor
        if np.isscalar(icells):
            icells=np.array([icells])
        ax     = a3.Axes3D(plt.figure())
        polys  = []
        for icell in icells:
            w      = self.cell_iwalls[icell]
            if(len(icells)==1): 
                colors = list(map("C{}".format, range(len(w))))
            for iwall in w:
                ivert = vor.ridge_vertices[iwall].copy()
                ivert.pop(-1)
                polys.append(self.vertices[ivert])
        pc = a3.art3d.Poly3DCollection(polys, facecolor=colors, edgecolor="k", alpha=alpha)
        ax.add_collection3d(pc)
        if bbox is None:
            bbox = self.bbox
        if bbox is not None:
            ax.set_xlim([bbox[0][0],bbox[0][1]])
            ax.set_ylim([bbox[1][0],bbox[1][1]])
            ax.set_zlim([bbox[2][0],bbox[2][1]])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()
        return ax

    def visualize_points(self,icells=None,colors="C1"):
        import mpl_toolkits.mplot3d
        if icells is None: icells=np.arange(grid.ncells)
        fig = plt.figure()
        ax  = fig.add_subplot(projection='3d')
        ax.scatter(self.cell_points[icells,0],self.cell_points[icells,1],self.cell_points[icells,2])
        plt.show()
        
    def write_radmc3d_unstr_grid(self,bin=False):
        self.compute_diagnostics()
        iformat=1
        isave_cellcenters = int(self.save_cellcenters)
        isave_vertices    = 0
        isave_size        = int(self.save_size)
        if bin:
            ihdr   = np.array([iformat,                        # Format number
                               8,                              # 8 means double precision
                               self.ncells,                    # Nr of cells (both "closed" and "open" ones)
                               self.nwalls,                    # Nr of cell walls in total
                               self.nverts,                    # Nr of vertices in total (0 means: not relevant)
                               self.cell_max_nr_walls,         # Max number of walls per cell
                               self.cell_max_nr_verts,         # Max number of vertices per cell (0 means: not relevant)
                               self.wall_max_nr_verts,         # Max number of vertices per wall (required if saving vertices)
                               self.vert_max_nr_cells,         # Max number of cells per vertex (0 means: not relevant)
                               self.ncells_open,               # Nr of "open" cells
                               self.hull_nwalls,               # Nr of cell walls at the surface
                               isave_cellcenters,              # Include a list of cell center positions?
                               isave_vertices,                 # Include the list of vertices?
                               isave_size,                     # Include the list of estimated linear cell sizes?
                               0],                             # Are the surface cell walls convex?
                              dtype='int64')                   # Create an array of integers for the header
            strh          = str(len(ihdr))+'Q'                 # Create a format string for struct for the header
            strp          = str(3*self.ncells)+'d'             # Create a format string for struct for the points
            strpv         = str(3*self.nverts)+'d'             # Create a format string for struct for the points
            strv          = str(self.ncells)+'d'               # Create a format string for struct for the volumes
            strw          = str(6*self.nwalls)+'d'             # Create a format string for struct for the walls
            stri          = str(2*self.nwalls)+'Q'             # Create a format string for struct for the indices
            stric         = str(self.ncells_open)+'Q'          # Create a format string for struct for the open cell indices
            striw         = str(self.hull_nwalls)+'Q'          # Create a format string for struct for the hull wall indices
            striwv        = str(self.nwalls*self.wall_max_nr_verts)+'Q' # Create a format string for struct for the wall vertex indices
            data_cellvols = self.cell_volumes.ravel()
            data_walls    = np.hstack((self.wall_s,self.wall_n)).ravel()
            indices       = (self.wall_icells+1).ravel()
            sdat          = struct.pack(strh,*ihdr) \
                            +struct.pack(strv,*data_cellvols) \
                            +struct.pack(strw,*data_walls) \
                            +struct.pack(stri,*indices)        # Create a binary image of the data
            if(self.save_cellcenters):
                data_points   = self.cell_points.ravel()
                sdat     += struct.pack(strp,*data_points)
            with open('unstr_grid.binp','w+b') as f:
                f.write(sdat)                                  # Write data in binary form
        else:
            with open('unstr_grid.inp','w') as f:
                f.write('{}\n'.format(iformat))                # Format number
                f.write('{}\n'.format(self.ncells))            # Nr of cells (both "closed" and "open" ones)
                f.write('{}\n'.format(self.nwalls))            # Nr of walls
                f.write('{}\n'.format(self.nverts))            # Nr of vertices
                f.write('{}\n'.format(self.cell_max_nr_walls)) # Max number of walls per cell
                f.write('{}\n'.format(self.cell_max_nr_verts)) # Max number of vertices per cell (0 means: not relevant)
                f.write('{}\n'.format(self.wall_max_nr_verts)) # Max number of vertices per wall (required if saving vertices)
                f.write('{}\n'.format(self.vert_max_nr_cells)) # Max number of cells per vertex (0 means: not relevant)
                f.write('{}\n'.format(self.ncells_open))       # Nr of "open" cells
                f.write('{}\n'.format(self.hull_nwalls))       # Nr of cell walls at the surface
                f.write('{}\n'.format(isave_cellcenters))      # Include a list of cell center positions?
                f.write('{}\n'.format(isave_vertices))         # Include the list of vertices?
                f.write('{}\n'.format(isave_size))             # Include the list of cell sizes?
                f.write('0\n')                                 # Are the surface cell walls convex?
                np.savetxt(f,self.cell_volumes)                # The cell volumes (0="open" cell)
                data=np.hstack((self.wall_s,self.wall_n))      # Glue the support and direction vectors
                np.savetxt(f,data)                             # Write the wall support and direction vectors
                np.savetxt(f,self.wall_icells+1,fmt='%d')      # Indices of cells are on each side of the wall (starting with 1, fortran style!)
                if(self.save_cellcenters):
                    data = self.cell_points
                    np.savetxt(f,data)                         # Write the cell centers

#
# EXAMPLE (uncomment to try)
# 
# npt    = 300
# npttry = npt*4
# radius = 1.0
# x      = (2*np.random.random(npttry)-1)*radius
# y      = (2*np.random.random(npttry)-1)*radius
# z      = (2*np.random.random(npttry)-1)*radius
# r      = np.sqrt(x**2+y**2+z**2)
# mask   = r<=radius
# x = x[mask]
# y = y[mask]
# z = z[mask]
# assert len(x)>=npt, 'Increase npttry'
# x = x[:npt]
# y = y[:npt]
# z = z[:npt]
# 
# points = np.vstack([x,y,z]).T
# 
# bbox   = [[-radius,radius],[-radius,radius],[-radius,radius]]
# #bbox   = None
# 
# grid   = Voronoigrid(points,bbox=bbox)
# 
# icells = np.where(grid.cell_volumes>0)[0]
# #grid.visualize_cells(icells[1],alpha=0.8)
# #grid.visualize_cells(icells[10:20],alpha=0.8)
# #grid.visualize_cells(icells,alpha=0.8) #,bbox=[[-1.5,2.5],[-1.5,2.5],[-1.5,2.5]])
# 
# #grid.visualize_points(icells)
# 
# grid.write_radmc3d_unstr_grid(bin=False)
