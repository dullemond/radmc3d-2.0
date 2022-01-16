import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay,ConvexHull
import struct
from tqdm import tqdm    # For a nice progress bar

class Delaunaygrid(object):
    def __init__(self,vertices,save_cellcenters=False,save_vertices=True, \
                 save_volumes=False,save_sn_vectors=False,qhull_options=None):
        """
        Create a Delaunay grid for RADMC-3D from a set of 3D vertices. Note that
        this means that the vertices are not representative of the cells. The
        values of density, temperature etc inside each cell must be either
        computed as a mean of the cell corners (=vertices) or as the linear 
        interpolation between them. In the latter case the cell does not contain 
        a single value, but instead a linear function. This allows for higher-order
        integration of the radiative transfer equation.

        Arguments:
          vertices           The vertices given as a numpy array vertices[npnts,3].
                             If set to None, the present unstr_grid.*inp file will be 
                             read, but only the vertices. The grid will be re-constructed
                             from these vertices.
          qhull_options      For passing to the Delaunay generator of SciPy

        Options:
          save_vertices:     If True, then the 3D locations of the vertices are 
                             written to the unstr_grid.inp file. This is not strictly
                             necessary for RADMC-3D, since the cell shapes are 
                             already defined by the cell walls alone. However, if
                             you want to make use to second-order radiative transfer
                             integration in RADMC-3D, then the vertices have to be
                             specified by setting save_vertices=True.
          save_cellcenters:  If True, then the 3D locations of the centers of the
                             cells are written to the unstr_grid.inp file. This is
                             not strictly necessary, as RADMC-3D can compute these
                             positions from the mean of the support vector positions
                             of the cell walls.
          save_volumes:      If True, then write the cell volumes to the file. This is
                             not really necessary, because RADMC-3D can calculate the
                             volumes of the tetraheders internally.
          save_sn_vectors:   If True, then write the support- and normal vectors to
                             the file. This is not really necessary, because RADMC-3D 
                             can calculate these internally. But if save_sn_vectors=False,
                             then save_vertices must be True. 
        """
        #
        # Set flags
        #
        self.save_cellcenters = save_cellcenters
        self.save_vertices    = save_vertices
        self.save_volumes     = save_volumes
        self.save_sn_vectors  = save_sn_vectors
        self.save_size        = False   # For now
        self.qhull_options    = qhull_options
        if (not self.save_sn_vectors) or (not self.save_volumes):
            assert self.save_vertices, 'ERROR: If you want the s and n vectors and/or the volumes to be calculated internally by RADMC-3D, you have to write out the vertices.'
        #
        # If vertices is None, then read vertices from a pre-existing file
        #
        if(vertices is None):
            print('Reading vertices from grid file...')
            self.read_radmc3d_unstr_grid()
            vertices = self.vertices
        #
        # The Qhull algorithms are sensitive to very large numbers, so let's
        # scale everything. First compute a useful scale factor
        #
        scale            = vertices.max(axis=0)-vertices.min(axis=0)
        scale            = scale.max()
        self.scale       = scale
        #
        # Create the Delaunay Triangulation
        #
        print('Running scipy.spatial.Delaunay() to create Delaunay triangulation...')
        self.tri         = Delaunay(vertices/scale,qhull_options=qhull_options)
        self.vertices    = self.tri.points*scale
        self.nverts      = vertices.shape[0]
        self.ncells      = self.tri.nsimplex
        self.ncells_open = 0
        self.icells_open = np.array([])
        #
        # Create cell center points
        #
        self.create_cell_centers()
        #
        # Now create the cell walls
        #
        print('Creating the cell walls...')
        self.create_cellwalls()
        print('Linking the cell walls back to the cells...')
        self.link_walls_to_cells()
        #
        # Link vertices to cells
        #
        print('Linking the vertices back to the cells...')
        self.link_vertices_to_cells()
        #
        # Compute the cell volumes
        #
        print('Computing the cell volumes...')
        self.compute_cell_volumes()
        #
        # Compute some diagnostics
        #
        self.compute_diagnostics()

    def compute_cell_volumes(self,useCH=False):
        """
        Compute the volume of the Delaunay cells
        """
        scale = self.scale
        self.cell_volumes = np.zeros(self.ncells)
        if useCH:
            #for i in range(self.tri.nsimplex):  # Without progress bar
            for i in tqdm(range(self.ncells)):
                indices = self.tri.simplices[i]
                self.cell_volumes[i] = scale**3 * ConvexHull(self.vertices[indices]/scale).volume
        else:
            v  = self.vertices[self.cell_iverts]
            cp = self.crossp((v[:,1,:]-v[:,-1,:]),(v[:,2,:]-v[:,-1,:]))
            self.cell_volumes = np.abs(((v[:,0,:]-v[:,-1,:])*cp).sum(axis=1))/6

    def create_cell_centers(self):
        self.cell_points = self.vertices[self.tri.simplices,:].mean(axis=1)
        
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

    def create_cellwalls(self):
        tri = self.tri
        #
        # First seek all walls facing vacuum (parts of the hull)
        #
        list_icell = []
        list_ivert = []
        for l in range(4):
            i  = np.where(tri.neighbors[:,l]<0)[0]
            k  = [0,1,2,3]
            k.pop(l)
            k  = np.array(k)
            vi = tri.simplices[i[:,None],k[None,:]]
            list_icell = list_icell + list(i)
            list_ivert = list_ivert + list(vi)
        list_icell = np.array(list_icell)
        list_ivert = np.array(list_ivert)
        #
        # Compute their support vectors, pointing to the
        # center of their faces
        #
        sv = self.vertices[list_ivert,:].mean(axis=1)
        #
        # Create the normal vectors
        #
        va = self.vertices[list_ivert[:,1]]-self.vertices[list_ivert[:,0]]
        vb = self.vertices[list_ivert[:,2]]-self.vertices[list_ivert[:,0]]
        n  = self.crossp(va,vb)
        l  = (n*n).sum(axis=1)
        assert(l.min()!=0), 'Error: One of the faces (cell walls) of the grid has zero surface area.'
        l  = np.sqrt(l)
        n /= l[:,None]
        #
        # Make sure the normals point outward
        #
        pc = self.cell_points[list_icell]
        dv = sv-pc
        sn = np.sign((dv*n).sum(axis=1))
        n *= sn[:,None]
        #
        # Store
        #
        hull_icell  = list_icell
        hull_ivert  = list_ivert
        hull_svec   = sv
        hull_n      = n
        hull_nwalls = len(hull_icell)
        hull_iwalls = np.arange(hull_nwalls)
        #
        # Now seek all walls inside. Avoid double-counting. 
        #
        ncells = self.ncells
        neigh  = tri.neighbors.copy()
        list_icell = []
        list_ivert = []
        for icell in range(ncells):
            for ineigh in range(4):
                icell_other = neigh[icell,ineigh]
                if icell_other>=0:
                    #
                    # Create the cell wall
                    #
                    k  = [0,1,2,3]
                    k.pop(ineigh)
                    k  = np.array(k)
                    list_ivert.append(tri.simplices[icell,k])
                    #
                    # Link both cells
                    #
                    list_icell.append(np.array([icell,icell_other]))
                    #
                    # Avoid double-counting
                    #
                    j = np.where(neigh[icell_other]==icell)[0]
                    assert len(j)==1, 'Something went wrong with Delaunay grid.'
                    ineigh_other = j[0]
                    neigh[icell,ineigh] = -1
                    neigh[icell_other,ineigh_other] = -1
        #
        # Convert to arrays
        #
        list_icell = np.array(list_icell)
        list_ivert = np.array(list_ivert)
        #
        # Compute the support vectors for all these inner walls, pointing to the
        # center of this face
        #
        sv = self.vertices[list_ivert,:].mean(axis=1)
        #
        # Create the normal vectors
        #
        va = self.vertices[list_ivert[:,1]]-self.vertices[list_ivert[:,0]]
        vb = self.vertices[list_ivert[:,2]]-self.vertices[list_ivert[:,0]]
        n  = self.crossp(va,vb)
        l  = (n*n).sum(axis=1)
        assert(l.min()!=0), 'Error: One of the faces (cell walls) of the grid has zero surface area.'
        l  = np.sqrt(l)
        n /= l[:,None]
        #
        # Make sure they point from cell 0 to cell 1 (the order in icell)
        #
        dv = self.cell_points[list_icell[:,1]]-self.cell_points[list_icell[:,0]]
        sn = np.sign((dv*n).sum(axis=1))
        n *= sn[:,None]
        #
        # Now put into the standard format
        #
        self.hull_nwalls = hull_nwalls
        self.nwalls      = hull_nwalls + len(list_icell)
        self.wall_s      = np.vstack((hull_svec,sv))
        self.wall_n      = np.vstack((hull_n,n))
        hicell           = np.zeros((hull_nwalls,2),dtype=np.int32)-1
        hicell[:,0]      = hull_icell[:]
        self.wall_icells = np.vstack((hicell,list_icell))
        self.wall_iverts = np.vstack((hull_ivert,list_ivert))
        self.hull_iwalls = hull_iwalls
        
    def link_walls_to_cells(self):
        self.cell_iwalls = [ [] for _ in range(self.ncells) ]
        self.cell_wsign  = [ [] for _ in range(self.ncells) ]
        for i in range(self.nwalls):
            icell = self.wall_icells[i,0]
            self.cell_iwalls[icell].append(i)
            self.cell_wsign[icell].append(1)
            icell = self.wall_icells[i,1]
            if icell>=0:
                self.cell_iwalls[icell].append(i)
                self.cell_wsign[icell].append(-1)

    def link_vertices_to_cells(self):
        self.cell_iverts = [ [] for _ in range(self.ncells) ]
        for iwall in range(self.nwalls):
            icell = self.wall_icells[iwall,0]
            self.cell_iverts[icell] += list(self.wall_iverts[iwall])
            icell = self.wall_icells[iwall,1]
            if icell>=0:
                self.cell_iverts[icell] += list(self.wall_iverts[iwall])
        for icell in range(self.ncells):
            self.cell_iverts[icell] = np.array(list(set(self.cell_iverts[icell])))
        self.cell_iverts = np.array(self.cell_iverts)

    def compute_diagnostics(self):
        self.ncells_closed = self.ncells
        self.volume_total   = self.cell_volumes.sum()
        self.cell_max_nr_walls = 4
        self.cell_max_nr_verts = 4
        self.wall_max_nr_verts = 3
        self.vert_max_nr_cells = 0
        
    def check_point_in_cell(self,p,icell,getdiag=False):
        iwalls = np.array(self.cell_iwalls[icell])
        sign   = np.array(self.cell_wsign[icell])
        s      = self.wall_s[iwalls,:]
        n      = sign[:,None]*self.wall_n[iwalls,:]
        inp    = ((s-p)*n).sum(axis=1)
        inpmin = inp.min()
        res    = (inpmin>=0)
        if(not getdiag): return res
        return res,iwalls,sign,s,n,inp

    def find_neighboring_cells(self,icell):
        ineigh = self.wall_icells[grid.cell_iwalls[icell]].flatten()
        ineigh = ineigh[ineigh!=icell]
        return ineigh
    
    def visualize_cells(self,icells=None,alpha=0.9,colors="C1",bbox=None):
        import mpl_toolkits.mplot3d as a3
        if icells is None: icells=np.arange(self.ncells)
        tri    = self.tri
        if np.isscalar(icells):
            icells=np.array([icells])
        ax     = a3.Axes3D(plt.figure())
        polys  = []
        for icell in icells:
            w      = self.cell_iwalls[icell]
            if(len(icells)==1): 
                colors = list(map("C{}".format, range(len(w))))
            for iwall in w:
                polys.append(self.vertices[self.wall_iverts[iwall,:]])
        pc = a3.art3d.Poly3DCollection(polys, facecolor=colors, edgecolor="k", alpha=alpha)
        ax.add_collection3d(pc)
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
        if icells is None: icells=np.arange(self.ncells)
        fig = plt.figure()
        ax  = fig.add_subplot(projection='3d')
        ax.scatter(self.cell_points[icells,0],self.cell_points[icells,1],self.cell_points[icells,2])
        plt.show()
        
    def visualize_vertices(self,ivert=None,colors="C1"):
        import mpl_toolkits.mplot3d
        if ivert is None: ivert=np.arange(self.nverts)
        fig = plt.figure()
        ax  = fig.add_subplot(projection='3d')
        ax.scatter(self.vertices[ivert,0],self.vertices[ivert,1],self.vertices[ivert,2])
        plt.show()
        
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
                               1],                             # Are the surface cell walls convex?
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
            sdat          = struct.pack(strh,*ihdr)            # Create a binary image of the data
            if(self.save_volumes):
                sdat     += struct.pack(strv,*data_cellvols)
            if(self.save_sn_vectors):
                sdat     += struct.pack(strw,*data_walls)
            sdat         += struct.pack(stri,*indices)
            if(self.save_cellcenters):
                data_points   = self.cell_points.ravel()
                sdat     += struct.pack(strp,*data_points)
            if(self.save_vertices):
                assert self.wall_iverts.shape[1]==self.wall_max_nr_verts, 'Max nr of vertices per wall is inconsistent'
                data_wall_iverts = (self.wall_iverts+1).ravel()
                sdat     += struct.pack(striwv,*data_wall_iverts)
                data_verts    = self.vertices.ravel()
                sdat     += struct.pack(strpv,*data_verts)
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
                f.write('1\n')                                 # Are the surface cell walls convex?
                if(self.save_volumes):
                    np.savetxt(f,self.cell_volumes)            # The cell volumes (0="open" cell)
                if(self.save_sn_vectors):
                    data=np.hstack((self.wall_s,self.wall_n))  # Glue the support and direction vectors
                    np.savetxt(f,data)                         # Write the wall support and direction vectors
                np.savetxt(f,self.wall_icells+1,fmt='%d')      # Indices of cells are on each side of the wall (starting with 1, fortran style!)
                if(self.save_cellcenters):
                    data = self.cell_points
                    np.savetxt(f,data)                         # Write the cell centers
                if(self.save_vertices):
                    assert self.wall_iverts.shape[1]==self.wall_max_nr_verts, 'Max nr of vertices per wall is inconsistent'
                    np.savetxt(f,self.wall_iverts+1,fmt='%d')  # Indices of vertices of each wall (starting with 1, fortran style!)
                    data = self.vertices
                    np.savetxt(f,data)                         # Write the vertices

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
                self.nverts       = header[4]
                self.wall_max_nr_verts = header[7]
                assert self.wall_max_nr_verts>0, 'Error: Cannot read unstr_grid.inp as a Delaunay because max nr of vertices per wall not given.'
                isave_cellcenters = header[11]
                isave_vertices    = header[12]
                assert isave_vertices==1, 'Error: Cannot read unstr_grid.binp as a Delaunay because vertices not given.'
                # We skip the rest, except for the vertices
                np.fromfile(f,dtype=float,count=self.ncells)   # Skip the cell volumes
                np.fromfile(f,dtype=float,count=self.nwalls*6) # Skip the cell walls s and n vectors
                np.fromfile(f,dtype=int,count=self.nwalls*2)   # Skip the cell walls cell indices
                if isave_cellcenters==1: dum = np.fromfile(f,dtype=float,count=self.ncells*3) # Skip the cell centers
                np.fromfile(f,dtype=int,count=self.nwalls*self.wall_max_nr_verts)   # Skip the cell walls vert indices
                self.vertices = np.fromfile(f,dtype=float,count=self.nverts*3) # Read the vertices
                self.vertices  = self.vertices.reshape((self.nverts,3))
        else:
            with open('unstr_grid.inp','r+') as f:
                iformat       = int(f.readline())
                self.ncells   = int(f.readline())
                self.nwalls   = int(f.readline())
                self.nverts   = int(f.readline())
                idum          = int(f.readline())
                idum          = int(f.readline())
                self.wall_max_nr_verts = int(f.readline())
                assert self.wall_max_nr_verts>0, 'Error: Cannot read unstr_grid.inp as a Delaunay because max nr of vertices per wall not given.'
                idum          = int(f.readline())
                idum          = int(f.readline())
                idum          = int(f.readline())
                isave_cellcenters = int(f.readline())
                isave_vertices    = int(f.readline())
                assert isave_vertices==1, 'Error: Cannot read unstr_grid.binp as a Delaunay because vertices not given.'
                idum          = int(f.readline())
                idum          = int(f.readline())
                dum           = np.loadtxt(f,float,max_rows=self.ncells)    # Skip the cell volumes
                dum           = np.loadtxt(f,float,max_rows=self.nwalls)    # Skip the wall s and n vectors
                dum           = np.loadtxt(f,int,max_rows=self.nwalls)      # Skip the wall cell indices
                if isave_cellcenters==1: dum = np.loadtxt(f,float,max_rows=self.ncells) # Skip the cell centers
                dum           = np.loadtxt(f,float,max_rows=self.nwalls)    # Skip the vertex links to walls
                self.vertices = np.loadtxt(f,float,max_rows=self.nverts)    # Read the vertices

    def interpolate_to_cell_centers(self,q_at_vertices):
        shape         = list(q_at_vertices.shape)
        shape[0]      = self.ncells
        q_in_cells    = np.zeros(shape)
        for icorner in range(4):
            iverts         = self.cell_iverts[:,icorner]
            q_in_cells[:] += q_at_vertices[iverts]
        q_in_cells /= 4
        return q_in_cells

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
# grid   = Delaunaygrid(points,save_cellcenters=True,save_vertices=True)
# 
# # icells = np.where(grid.cell_volumes>0)[0]
# # #grid.visualize_cells(icells[1],alpha=0.8)
# # #grid.visualize_cells(icells[10:20],alpha=0.8)
# # grid.visualize_cells(icells,alpha=0.8) #,bbox=[[-1.5,2.5],[-1.5,2.5],[-1.5,2.5]])
# # 
# # grid.visualize_points(icells)
# # 
# # #grid.write_radmc3d_unstr_grid(bin=True)
# 
# grid.write_radmc3d_unstr_grid(bin=False)
