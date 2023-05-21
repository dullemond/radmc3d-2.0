import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay,Voronoi,ConvexHull
import struct
from tqdm import tqdm    # For a nice progress bar

class Voronoigrid(object):
    def __init__(self,points,bbox=None,qhull_options=None,compvolume=True):
        """
        Create a Voronoi grid for RADMC-3D from a set of 3D points.

        Arguments:
          points             The points given as a numpy array points[npnts,3].
                             If set to None, the present unstr_grid.*inp file will be 
                             read, but only the points. The grid will be re-constructed
                             from these points.
        Options:
          bbox               If set to [[xmin,xmax],[ymin,ymax],[zmin,zmax]] this sets
                             the bounding box: only cells that fit entirely inside this
                             box are considered actual cells (with finite volume). The
                             cells that have >=1 point outside the box are treated as
                             "open regions", see below.
          qhull_options      For passing to the Voronoi generator of SciPy
          compvolume         If False, then skip the computation of the cell volumes.
                             Can be useful if you import these values externally.
                             Default is True.

        Note on "open regions": These are "cells" that are not closed, and thus 
        have an infinite volume. These regions are still useful for guiding the
        ray of RADMC-3D to the right cells, but they are considered to be empty.
        If you use bbox, then, mathematically speaking, one could close these
        cells (as well as the cells that do not fully fit into the box) at the
        box walls. But the additional complexity is not worth it, so instead, 
        these cells are considered empty.
        """
        #
        # Set flags
        #
        self.save_cellcenters = True    # Must be!
        self.save_vertices    = False   # For now
        self.save_volumes     = True    # As long as RADMC-3D cannot do this internally
        self.save_sn_vectors  = False   # Not necessary
        self.save_size        = False   # For now
        self.qhull_options    = qhull_options
        #
        # If points is None, then read points from a pre-existing file
        #
        if(points is None):
            print('Reading cell center points from grid file...')
            self.read_radmc3d_unstr_grid()
            points = self.cell_points
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
        self.vor         = Voronoi(points/scale,qhull_options=qhull_options)
        self.cell_points = self.vor.points*scale
        self.vertices    = self.vor.vertices*scale
        self.npnts       = points.shape[0]
        self.nverts      = 0       # For now
        #
        # Find all open cells
        #
        self.find_open_cells()
        #
        # Now create the cell walls
        #
        print('Creating the cell walls...')
        self.create_cellwalls()
        #
        # Create the Delaunay triangulation for completing the cell wall list
        #
        print('Running scipy.spatial.Delaunay() to create Delaunay triangulation (for wall completeness check)...')
        self.tri = Delaunay(points/scale,qhull_options=qhull_options)
        #
        # Link cell walls back to cells
        #
        print('Linking the cell walls back to the cells...')
        self.link_walls_to_cells()
        #
        # Complete the cell wall list
        #
        print('Completing the cell walls of the open cells...')
        self.complete_cellwalls()
        #
        # Link cell walls back to cells again
        #
        print('Linking the cell walls back to the cells again...')
        self.link_walls_to_cells()
        #
        # Link vertices to cells
        #
        print('Linking the vertices back to the cells...')
        self.link_vertices_to_cells(remove_nonverts=True)
        #
        # Compute the cell volumes. This also automatically notices
        # which points correspond to "real" cells (with a finite volume)
        # and which points correspond to "open" regions going out to
        # infinity (with an infinite volume). These "open" regions are
        # marked by a volume = 0.0
        #
        if compvolume:
            print('Computing the cell volumes...')
            self.compute_cell_volumes()
        #
        # If bbox is given, then reject all cells that have Voronoi vertices
        # that are outside of the box.
        #
        self.bbox = bbox
        if compvolume and bbox is not None:
            print('Implementing bounding box by rejecting cells that are not entirely in the bbox...')
            self.bbox_reject_cells()
        #
        # Compute some diagnostics
        #
        self.compute_diagnostics()

    def find_open_cells(self):
        self.cell_iopen = []
        for i in range(self.npnts):
            indices = self.vor.regions[self.vor.point_region[i]]
            if -1 in indices: self.cell_iopen.append(i)
    
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

    def link_vertices_to_cells(self,remove_nonverts=True):
        self.cell_iverts = [ [] for _ in range(self.npnts) ]
        for i in range(self.npnts):
            iverts = self.vor.regions[self.vor.point_region[i]]
            iverts = np.array(iverts)
            if remove_nonverts:
                iverts = iverts[iverts>=0]
            self.cell_iverts[i] = iverts
                
    def crossp(self,v1,v2):
        v1 = np.array(v1)
        v2 = np.array(v2)
        if(len(v1.shape)==1):
            nv = 1
            assert(len(v2.shape)==1), 'v1 and v2 must has same shape'
        else:
            nv = v1.shape[0]
            assert(len(v2.shape)==2), 'v1 and v2 must has same shape'
            assert(v2.shape[0]==nv), 'v1 and v2 must has same shape'
        dm = v1.shape[-1]
        v1 = v1.reshape((nv,dm))
        v2 = v2.reshape((nv,dm))
        if dm==2:
            c    = np.zeros(nv)
            c[:] = v1[:,0]*v2[:,1] - v1[:,1]*v2[:,0]
        else:
            c      = np.zeros((nv,dm))
            c[:,0] = v1[:,1]*v2[:,2] - v1[:,2]*v2[:,1]
            c[:,1] = v1[:,2]*v2[:,0] - v1[:,0]*v2[:,2]
            c[:,2] = v1[:,0]*v2[:,1] - v1[:,1]*v2[:,0]
        return c

    def create_cellwalls(self,remove_nonverts=True):
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
        self.wall_iverts = [ [] for _ in range(self.nwalls) ]
        for i in range(self.nwalls):
            iverts = self.vor.ridge_vertices[i]
            iverts = np.array(iverts)
            if remove_nonverts:
                iverts = iverts[iverts>=0]
            self.wall_iverts[i] = iverts

    def complete_cellwalls(self):
        points          = self.cell_points
        indptr, indices = self.tri.vertex_neighbor_vertices
        add_wall_s      = []
        add_wall_n      = []
        add_wall_icells = []
        add_wall_iverts = []
        #for icell in self.cell_iopen:
        for icell in range(self.npnts):
            neighbors_voronoi  = list(self.find_neighboring_cells(icell))
            lneighvor = len(neighbors_voronoi)
            neighbors_delaunay = list(indices[indptr[icell]:indptr[icell+1]])
            lneighdel = len(neighbors_delaunay)
            assert lneighdel>=lneighvor, 'ERROR: Nr of Ridges of Voronoi larger than nr of neighbors in Delaunay'
            if lneighdel>lneighvor:
                icell_missing = list(set(neighbors_delaunay)-set(neighbors_voronoi))
                ss = ""
                for i in icell_missing:
                    ss += " {}".format(i)
                #print('In open cell {} missing neighbors are: '.format(icell)+ss)
                for i in icell_missing:
                    if i>icell:  # Avoid double counting
                        add_wall_s.append(0.5*(points[i,:] + points[icell,:]))
                        n = points[i,:] - points[icell,:]
                        l = ((n**2).sum())**0.5
                        n /= l
                        add_wall_n.append(n)
                        add_wall_icells.append(np.array([icell,i]))
                        add_wall_iverts.append([-1])
        addnwall = len(add_wall_n)
        add_wall_s = np.array(add_wall_s)
        add_wall_n = np.array(add_wall_n)
        add_wall_icells = np.array(add_wall_icells)
        if(addnwall>0):
            print('Found {} missing walls and adding them.'.format(addnwall))
            self.wall_s = np.vstack((self.wall_s,add_wall_s))
            self.wall_n = np.vstack((self.wall_n,add_wall_n))
            self.wall_icells = np.vstack((self.wall_icells,add_wall_icells))
            self.wall_iverts = self.wall_iverts + add_wall_iverts
            self.nwalls += addnwall
        assert self.nwalls==len(self.tri.vertex_neighbor_vertices[1])//2,'ERROR: Something went wrong with the completion of the cell walls'

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
        #self.cell_iopen    = np.where((self.cell_volumes<=0.0))[0]
        self.ncells_open   = len(self.cell_iopen)
        if hasattr(self,'cell_volumes'):
            if self.bbox is None:
                self.ncells_closed = len(np.where((self.cell_volumes>0.0))[0])
                assert self.ncells_open+self.ncells_closed==self.ncells, 'Internal error.'
            self.volume_total   = self.cell_volumes.sum()
        else:
            if self.bbox is None:
                self.ncells_closed = self.ncells - self.ncells_open
        l=np.array([len(self.cell_iwalls[i]) for i in range(self.npnts)])
        self.cell_max_nr_walls = l.max()
        self.cell_max_nr_verts = 0
        self.wall_max_nr_verts = 0
        self.vert_max_nr_cells = 0
        self.hull_nwalls       = 0

    def check_point_in_wall_plane(self,p,iwall,getdiag=False):
        verts = self.wall_get_vertices(iwall,close=True)
        if(verts.shape[0]<4): return False
        dv    = verts[1:,:]-verts[:-1,:]
        dvp   = p-verts[:-1,:]
        scale = np.sqrt((dv**2).sum(axis=1)).max()
        dv   /= scale
        dvp  /= scale
        cp    = self.crossp(dvp[:-1],dvp[1:])
        err   = np.abs((cp[0:-2]*dv[2:-1]).sum(axis=1)).max()
        tol   = 1e-12
        res   = err<tol
        if(not getdiag): return res
        return res,verts,dv,dvp,cp,err
        
    def check_point_in_cell(self,p,icell,iwallin=None,getdiag=False):
        iwalls = np.array(self.cell_iwalls[icell])
        sign   = np.array(self.cell_wsign[icell])
        if iwallin is not None:
            # Exclude wall on which the point is currently
            mask   = iwalls!=iwallin
            iwalls = iwalls[mask]
            sign   = sign[mask]
            # But then check if the point is indeed on the plane of this wall
            onwall = self.check_point_in_wall_plane(p,iwallin)
        s      = self.wall_s[iwalls,:]
        n      = sign[:,None]*self.wall_n[iwalls,:]
        inp    = ((s-p)*n).sum(axis=1)
        inpmin = inp.min()
        res    = (inpmin>=0)
        if iwallin is not None:
            res = res and onwall
        if(not getdiag): return res
        return res,iwalls,sign,s,n,inp

    def find_neighboring_cells(self,icell):
        ineigh = self.wall_icells[self.cell_iwalls[icell]].flatten()
        ineigh = ineigh[ineigh!=icell]
        return ineigh

    def wall_get_vertices(self,iwall,close=True):
        iverts = np.array(self.wall_iverts[iwall])
        iverts = iverts[iverts>=0]
        if close: iverts = np.hstack((iverts,iverts[0]))
        verts  = self.vertices[iverts]
        return verts
        
    def visualize_cells(self,icells=None,alpha=0.9,colors="C1",bbox=None,wireframe=False,ax=None):
        if icells is None: icells=np.arange(self.ncells)
        vor    = self.vor
        if np.isscalar(icells):
            icells=np.array([icells])
        iwalls = []
        for icell in icells:
            w      = self.cell_iwalls[icell]
            iwalls += list(w)
        iwalls = list(set(iwalls))
        if(len(icells)==1): 
            colors = list(map("C{}".format, range(len(iwalls))))
        polys  = []
        opens  = []
        for iwall in iwalls:
            ivert = self.wall_iverts[iwall] # vor.ridge_vertices[iwall].copy()
            o = -1 in vor.ridge_vertices[iwall]
            polys.append(self.vertices[ivert])
            opens.append(o)
        if wireframe:
            if ax is None:
                import mpl_toolkits.mplot3d
                fig = plt.figure()
                ax  = fig.add_subplot(projection='3d')
            for i in range(len(polys)):
                if np.isscalar(colors):
                    color = colors
                else:
                    color = colors[i]
                p = polys[i]
                p = np.vstack((p,p[0]))
                if opens[i]:
                    ax.plot(p[:,0],p[:,1],p[:,2],':',color=color)
                else:
                    ax.plot(p[:,0],p[:,1],p[:,2],color=color)
        else:
            if ax is None:
                import mpl_toolkits.mplot3d as a3
                ax = a3.Axes3D(plt.figure())
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

    def visualize_points(self,icells=None,colors="C1",ax=None):
        import mpl_toolkits.mplot3d
        if icells is None: icells=np.arange(self.ncells)
        if ax is None:
            fig = plt.figure()
            ax  = fig.add_subplot(projection='3d')
        ax.scatter(self.cell_points[icells,0],self.cell_points[icells,1],self.cell_points[icells,2])
        plt.show()
        return ax
        
    def write_radmc3d_unstr_grid(self,bin=False):
        self.compute_diagnostics()
        iformat           = 2
        isave_cellcenters = int(self.save_cellcenters)
        isave_vertices    = int(self.save_vertices)
        isave_volumes     = int(self.save_volumes)
        isave_sn_vectors  = int(self.save_sn_vectors)
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
                               isave_volumes,                  # Include a list of cell volumes?
                               isave_sn_vectors,               # Include a list of s and n vectors?
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
            #data_walls    = np.hstack((self.wall_s,self.wall_n)).ravel()
            indices       = (self.wall_icells+1).ravel()
            sdat          = struct.pack(strh,*ihdr) \
                            +struct.pack(strv,*data_cellvols) \
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
                f.write('{}\n'.format(isave_volumes))          # Include a list of cell volumes?
                f.write('{}\n'.format(isave_sn_vectors))       # Include a list of s and n vectors?
                f.write('{}\n'.format(isave_cellcenters))      # Include a list of cell center positions?
                f.write('{}\n'.format(isave_vertices))         # Include the list of vertices?
                f.write('{}\n'.format(isave_size))             # Include the list of cell sizes?
                f.write('0\n')                                 # Are the surface cell walls convex?
                np.savetxt(f,self.cell_volumes)                # The cell volumes (0="open" cell)
                #data=np.hstack((self.wall_s,self.wall_n))      # Glue the support and direction vectors
                #np.savetxt(f,data)                             # Write the wall support and direction vectors
                np.savetxt(f,self.wall_icells+1,fmt='%d')      # Indices of cells are on each side of the wall (starting with 1, fortran style!)
                if(self.save_cellcenters):
                    data = self.cell_points
                    np.savetxt(f,data)                         # Write the cell centers

    def read_radmc3d_unstr_grid(self):
        import os.path
        txt = os.path.isfile('unstr_grid.inp')
        bin = os.path.isfile('unstr_grid.binp')
        assert (bin and not txt) or (txt and not bin), 'Need either unstr_grid.inp or unstr_grid.binp to read. But not both.'
        nihead = 14
        if(bin):
            with open('unstr_grid.binp','r+b') as f:
                header            = np.fromfile(f,dtype=int,count=nihead+1)
                self.ncells       = header[2]
                self.nwalls       = header[3]
                isave_cellcenters = header[11]
                assert isave_cellcenters==1, 'Error: Cannot read unstr_grid.binp as a Voronoi because cell centers not given.'
                # We skip the rest, except for the cell centers
                np.fromfile(f,dtype=float,count=self.ncells)   # Skip the cell volumes
                np.fromfile(f,dtype=float,count=self.nwalls*6) # Skip the cell walls s and n vectors
                np.fromfile(f,dtype=int,count=self.nwalls*2)   # Skip the cell walls cell indices
                self.cell_points  = np.fromfile(f,dtype=float,count=self.ncells*3) # Read the cell centers
                self.cell_points  = self.cell_points.reshape((self.ncells,3))
        else:
            with open('unstr_grid.inp','r+') as f:
                iformat       = int(f.readline())
                self.ncells   = int(f.readline())
                self.nwalls   = int(f.readline())
                self.nverts   = int(f.readline())
                idum          = int(f.readline())
                idum          = int(f.readline())
                idum          = int(f.readline())
                idum          = int(f.readline())
                idum          = int(f.readline())
                idum          = int(f.readline())
                isave_cellcenters = int(f.readline())
                idum          = int(f.readline())
                idum          = int(f.readline())
                idum          = int(f.readline())
                dum           = np.loadtxt(f,float,max_rows=self.ncells)    # Skip the cell volumes
                dum           = np.loadtxt(f,float,max_rows=self.nwalls)    # Skip the wall s and n vectors
                dum           = np.loadtxt(f,int,max_rows=self.nwalls)      # Skip the wall cell indices
                self.cell_points = np.loadtxt(f,float,max_rows=self.ncells) # Read the cell centers
                
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
