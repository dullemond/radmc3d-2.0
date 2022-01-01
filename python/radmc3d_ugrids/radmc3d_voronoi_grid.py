import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay,Voronoi,ConvexHull
from tqdm import tqdm    # For a nice progress bar

class Voronoigrid(object):
    def __init__(self,points,bbox=None):
        """
        Create a Voronoi grid for RADMC-3D from a set of 3D points.

        Arguments:
          points      The points given as a numpy array points[npnts,3].
        
        Options:
          bbox        If set to [[xmin,xmax],[ymin,ymax],[zmin,zmax]] this sets
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
        #
        # Create the Voronoi tesselation
        #
        print('Running scipy.spatial.Voronoi() to create Voronoi tesselation...')
        self.vor    = Voronoi(points)
        self.points = self.vor.points
        self.npnts  = points.shape[0]
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
        self.cell_volumes = np.zeros(self.vor.npoints)
        #for i in range(self.npnts):  # Without progress bar
        for i in tqdm(range(self.npnts)):
            indices = self.vor.regions[self.vor.point_region[i]]
            if -1 in indices: # some regions can be open
                self.cell_volumes[i] = 0.0
            else:
                self.cell_volumes[i] = ConvexHull(self.vor.vertices[indices]).volume

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
        maskx = np.logical_or(vor.vertices[:,0]<bbox[0][0],vor.vertices[:,0]>bbox[0][1])
        masky = np.logical_or(vor.vertices[:,1]<bbox[1][0],vor.vertices[:,1]>bbox[1][1])
        maskz = np.logical_or(vor.vertices[:,2]<bbox[2][0],vor.vertices[:,2]>bbox[2][1])
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
        self.nwalls      = self.vor.ridge_points.shape[0]
        self.wall_v      = 0.5*(points[self.vor.ridge_points[:,1],:] + points[self.vor.ridge_points[:,0],:])
        self.wall_n      = points[self.vor.ridge_points[:,1],:] - points[self.vor.ridge_points[:,0],:]
        l = (self.wall_n**2).sum(axis=1)**0.5
        self.wall_n     /= l[:,None]
        self.wall_cells  = self.vor.ridge_points

    def link_walls_to_cells(self):
        self.cell_walls  = [ [] for _ in range(self.npnts) ]
        self.cell_wsign  = [ [] for _ in range(self.npnts) ]
        for i in range(self.nwalls):
            icell = self.wall_cells[i,0]
            self.cell_walls[icell].append(i)
            self.cell_wsign[icell].append(1)
            icell = self.wall_cells[i,1]
            self.cell_walls[icell].append(i)
            self.cell_wsign[icell].append(-1)

    def compute_diagnostics(self):
        self.n_cells_open   = len(np.where((self.cell_volumes<=0.0))[0])
        self.n_cells_closed = len(np.where((self.cell_volumes>0.0))[0])
        self.volume_total   = self.cell_volumes.sum()

    def visualize_cells(self,icells,alpha=0.9,colors="C1",bbox=None):
        import mpl_toolkits.mplot3d as a3
        vor    = self.vor
        if np.isscalar(icells):
            icells=np.array([icells])
        ax     = a3.Axes3D(plt.figure())
        polys  = []
        for icell in icells:
            w      = self.cell_walls[icell]
            if(len(icells)==1): 
                colors = list(map("C{}".format, range(len(w))))
            for iwall in w:
                polys.append(vor.vertices[vor.ridge_vertices[iwall]])
        pc = a3.art3d.Poly3DCollection(polys, facecolor=colors, edgecolor="k", alpha=alpha)
        ax.add_collection3d(pc)
        ax.dist=10
        ax.azim=30
        ax.elev=10
        if bbox is None:
            bbox = self.bbox
        if bbox is not None:
            ax.set_xlim([bbox[0][0],bbox[0][1]])
            ax.set_ylim([bbox[1][0],bbox[1][1]])
            ax.set_zlim([bbox[2][0],bbox[2][1]])
        plt.show()
        
npt    = 300
x      = np.random.random(npt)
y      = np.random.random(npt)
z      = np.random.random(npt)

points = np.vstack([x,y,z]).T

#bbox   = [[0,1],[0,1],[0,1]]
bbox   = None

grid   = Voronoigrid(points,bbox=bbox)

icells = np.where(grid.cell_volumes>0)[0]
#grid.visualize_cells(icells[1],alpha=0.8)
#grid.visualize_cells(icells[10:20],alpha=0.8)
grid.visualize_cells(icells,alpha=0.8,bbox=[[-1.5,2.5],[-1.5,2.5],[-1.5,2.5]])
