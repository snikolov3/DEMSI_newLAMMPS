from __future__ import print_function, absolute_import, division

import numpy as np
import scipy.spatial

nc_avail = False
try:
    import netCDF4 as nc4
    nc_avail=True
except:
    print("Warning: netCDF4 python package not installed, netcdf functionality will not be available")

class Bond:
    def __init__(self, bondid, id1, id2, type = 1, *args, **kwargs):
        """
        Create a bond
        
        Parameters:
        -----------
        bondid : int
          id of bond
        i1 : int
          id of particle 1
        i2 : int
          id of particle 2
        type : int
          Bond type
        **kwargs : any additional per-bond properties (e.g. 'id')
        """
        self.bondid = bondid
        self.id1 = id1
        self.id2 = id2
        self.type = type
        for k,v in kwargs.items():
            setattr(self, k, v)

class Particles:
    """
    Python class to hold data structure for a collection of particles
    
    Parameters
    ----------
    x : 1D array of float
        x-coordinates of particles
    y : 1D array of float
        y-coordinates of particles
    radius : 1D array of float
        radii of particles
    iceFraction : 1D array of float
        Fraction of each particle covered in sea ice. Dimension is nParticles.
        If unspecified, set to [1.0] for all particles.
    iceThickness : 1D array of float
        Thickness of initial sea ice within element, in meters. Dimension
        is nParticles If unspecified, set to [0.1] for all particles.
    density : 1D array of float
        Area densities of particles. Assumed 0.1 kg/m^2 if not specified
    type : 1D array of int
        Particle types, assumed 1 for all particles if not specified
    id : 1D array of int
        Global ids of all particles, assumed 1..nParticles if not specified
    xlo, ylo, xhi, yhi, zlo, zhi : float
        Domain bounds. Assumed 
    xperiodic, yperiodic : bool
        Indicate periodicity, default is False
        
    Examples
    --------
    
    To create a collection of monosized particles on a lattice 
    (lattice dimensions of 100 X 200, particle radius = 5000)
    
    >>>import pydemsi
    >>>mylattice = pydemsi.Particles.from_lattice(100,20,5000)
    
    Optionally, form bonds between lattice neighbors:
    >>>mylattice.create_bonds()
    
    To create a collection of particles given x, y coordinates and radii:
    (pydemsi fills in the remaining parameters with defaults, as noted above)
    >>>mystructure = pydemsi.Particles(x, y, radii)
    
    In either case, the structure (`mylattice` or `mystructure`) can be written to a variety of formats, e.g.:
    >>>mylattice.to_particle_netcdf("mylatice.nc")
    >>>mylattice.to_lammps_data("mylattice.data")
    >>>mylattice.to_vtk("mylattice.vtk") #for ParaView viz.
    
    If matplotlib is installed, you can view the structure (optionally, with bonds drawn in):
    >>>mylattice.plot(draw_bonds=True)
    
    """
    def __init__(self, x, y, radius, bonds = [],
                 iceFraction = None, iceThickness = None,
                 density = None, type = None, id = None,
                 xlo = None, ylo = None, xhi = None, yhi = None,
                 zlo = None, zhi = None,
                 xperiodic = False, yperiodic = False, **kwargs):
        """
        Build a list of particles from numpy arrays.

        Parameters
        ---------- 
        x, y, radius : 1-d array-like of float
           Lists or 1-d numpy arrays of particle x, y coordinates and radii. Must 
           all be same length.
        density : 1-d array-like of float, or single float
            
        **kwargs : any additional per-particle or global properties.
           Any keyword/value pair can be passed to define additional attributes.
        """
        
        defaultDensity = 900 #kg/m^2; defaults to this value unless otherwise indicated
        defaultThickness = 1.0 # iceThickness and iceFraction were swapped before
        defaultFraction  = 0.1 # iceThickness and iceFraction were swapped before
        
        if np.ndim(x) != 1:
            print("Error: x must be 1-d list or array")
            return

        self.num_particles = len(x)
        self.bonds = bonds
        self.num_bonds = len(bonds)

        if len(y) != self.num_particles:
            print("Error: number of y entries does not match number of particles")
            return
        if len(radius) != self.num_particles:
            print("Error: number of radius entries does not match number of particles")

        if np.any(np.array(x)-np.array(radius) < 0):
            print("Warning: some particles extend below x = 0, this may be problematic with DEMSI grid restrictions")
        if np.any(np.array(y)-np.array(radius) < 0):
            print("Warning: some particles extend below y = 0, this may be problematic with DEMSI grid restrictions")
            
        if density is not None:
            if (len(density) == 1):
                density = [density]*self.num_particles
            else:
                if (len(density) != self.num_particles):
                    print("Error: number of density entries does not match number of particles")
                    return
        else:
            density = np.array([defaultDensity]*self.num_particles) #assume 900 kg/m^2 by default
                
        self.x = np.array(x)
        self.y = np.array(y)
        self.radius = np.array(radius)
        self.density = np.array(density)
        
        #Need additional error checking here
        if iceFraction is not None:
            self.iceFraction = np.array(iceFraction)
            
        if iceThickness is not None:
            self.iceThickness = np.array(iceThickness)

        for k,v in kwargs.items():
            setattr(self, k, np.array(v))
            
        #Set defaults if not specifically set

        #Masses, density and iceThickness/iceFraction definitions
        # Define mass based on iceThickness/iceFraction if available, 
        # otherwise assume iceThickness = 0.1m, iceFraction = 1
                
        if iceThickness is not None:
            if iceFraction is not None:
                self.mass = iceFraction*iceThickness*self.density*np.pi*self.radius**2
            else:
                print("Error: iceThickness is defined, but iceFraction is not")
                return
        elif iceThickness is None and iceFraction is None:
            self.mass = self.density*np.pi*self.radius**2*defaultThickness
            self.iceThickness = np.array([defaultThickness]*self.num_particles)
            self.iceFraction = np.array([defaultFraction]*self.num_particles)
        
        self.type = np.array([1]*self.num_particles) if type is None else np.array(type)
        self.id = np.arange(1,self.num_particles+1) if id is None else np.array(id)

        if len(self.type) != self.num_particles or len(self.id) != self.num_particles:
            print("Error: number of type and/or id entries does not match number of particles")
            return
        
        maxrad = np.max(self.radius)
        self.xlo = np.min(self.x) - maxrad if xlo is None else xlo 
        self.ylo = np.min(self.y) - maxrad if ylo is None else ylo
        self.xhi = np.max(self.x) + maxrad if xhi is None else xhi 
        self.yhi = np.max(self.y) + maxrad if yhi is None else yhi
        self.zlo = -maxrad if zlo is None else zlo
        self.zhi = maxrad if zhi is None else zhi
        
        self.xperiodic = xperiodic
        self.yperiodic = yperiodic
        self.zperiodic = True
         
        self.bonds = bonds

    @classmethod
    def from_lammps_data(cls, filename):
        """
        Create Particles object from LAMMPS data file
        
        Parameters
        ----------
        filename : string
            Name of lammps data file
            
        """
        bonds = []
        with open(filename,"r") as f:
            for line in f:
                words = line.split()
                if len(words) > 1:
                    if words[1] == "atoms":
                        natoms = int(words[0])
                    elif words[1] == "bonds":
                        nbonds = int(words[0])
                    elif words[0] == "Atoms":
                        f.readline()
                        id = []
                        x = []
                        y = []
                        radius = []
                        type = []
                        density = []
                        for _ in range(natoms):
                            words = f.readline().split()
                            id.append(int(words[0]))
                            type.append(int(words[1]))
                            radius.append(float(words[2])*0.5)
                            density.append(float(words[3]))
                            x.append(float(words[4]))
                            y.append(float(words[5]))
                    elif words[0] == "Bonds":
                        f.readline()                                              
                        for _ in range(nbonds):
                            id, type, id1, id2 = map(int, f.readline().split())
                            bonds.append(Bond(type, id1, id2, id = id))
                                
                if len(words) > 2:
                    if words[1] == "atom" and words[2] == "types":
                        n_atom_types = int(words[0])
                    if words[1] == "bond" and words[2] == "types":
                        n_bond_types = int(words[1])
                if len(words) > 3:
                    if words[2] == "xlo" and words[3] == "xhi":
                        xlo, xhi = map(float, words[0:2])
                    if words[2] == "ylo" and words[3] == "yhi":
                        ylo, yhi = map(float, words[0:2])
                    if words[2] == "zlo" and words[3] == "zhi":
                        zlo, zhi = map(float, words[0:2])
                    
           
            return cls(x, y, radius, id = id, type = type, density = density,
                        xlo = xlo, ylo = ylo, zlo = zlo,
                        xhi = xhi, yhi = yhi, zhi = zhi, bonds = bonds)
                                          
    @classmethod
    def from_particle_netcdf(cls, filename):
        """
        Create Particles object from particles netcdf file
        
        Parameters
        ----------
        filename : string
            Name of netcdf file
            
        """
        if not nc_avail:
            print("Function requires netcdf python package, exiting")
            return
        
        fileIn = nc4.Dataset(filename,"r")

        nParticles = len(fileIn.dimensions["nParticles"])

        id = fileIn.variables["globalID"]
        try:
            type = fileIn.variables["type"]
        except:
            type = np.ones(nParticles)
        x = fileIn.variables["x"][:,0]
        y = fileIn.variables["x"][:,1]
        radius = fileIn.variables["radius"]
        try:
            iceThickness = fileIn.variables["iceThickness"]
        except:
            iceThickness = np.ones(nParticles)
        try:
            iceFraction = fileIn.variables["iceFraction"]
        except:
            iceFraction = np.ones(nParticles)

        # Does not set up the bod correctly it seems
        #nBonds = len(fileIn.dimensions["nBonds"])
        #if nBonds > 0:
        #    bonds = fileIn.variables["bonds"]
        #else:
        #    bonds = []

        nBonds = 0
        bonds = []

        return cls(x, y, radius, id = id, type = type,
                   iceThickness = iceThickness, iceFraction = iceFraction,
                   bonds = bonds)

    
    @classmethod
    def from_lammps_dump(cls, filename, requested_timestep):
        """
        Create Particles object from lammps dump file. Dump file
        must contain fields as follows:
        id, type, mass, diameter, x, y, z
        
        Parameters
        ----------
        filename : string
            Name of lammps dump file
        requested_timestep : int
            Time step in the dump file to extract particles from.
            
        """
        id = []
        type = []
        density = []
        radius = []
        x = []
        y = []
        with open(filename, "r") as input_file:
            for line in input_file:
                if "ITEM: TIMESTEP" in line:
                    timestep = int(input_file.readline())
                    if timestep == requested_timestep:
                        if "ITEM: NUMBER OF ATOMS" not in input_file.readline():
                            print("Missing expected 'ITEM: NUMBER OF ATOMS' line in dump file\n")
                            return -1
                        num_atoms = int(input_file.readline())
                        if "ITEM: BOX BOUNDS" not in input_file.readline():
                            print("Missing expected 'ITEM: BOX BOUNDS' line in dump file\n")
                            return -1
                        xlo, xhi = map(float, input_file.readline().split())
                        ylo, yhi = map(float, input_file.readline().split())
                        zlo, zhi = map(float, input_file.readline().split())
                        if "ITEM: ATOMS" not in input_file.readline():
                            print("Missing expected 'ITEM: ATOMS' line in dump file\n")
                            return -1
                        else:
                            for i in range(0, num_atoms):
                                words = input_file.readline().split()
                                id.append(int(words[0]))
                                type.append(int(words[1]))
                                mass = float(words[2])
                                rad = 0.5*float(words[3])
                                density.append(mass/(np.pi*rad*rad))
                                radius.append(rad)
                                x.append(float(words[4]))
                                y.append(float(words[5]))
                            return cls(x, y, radius, id = id, type = type, 
                                       density = density, xlo = xlo, ylo = ylo, zlo = zlo, 
                                          xhi = xhi, yhi = yhi, zhi = zhi)
                    
        print("Could not find requested time step "+str(requested_timestep)+" in dump file")
        
    @classmethod
    def from_lattice(cls, nx, ny, radius):
        """
        Create monosized particles on a regular square lattice
        
        All other particle settings are set to default values (see `Particles` doc string),
        but can easily be changed subsequently.
        
        Parameters
        ----------
        nx : int
            number of lattice points (particles) in x direction
        ny : int
            number of lattice points (particles) in y direction
        radius : float
            radius of all particles
        
        """
        x = np.linspace(radius, (nx*2-1)*radius, nx)
        y = np.linspace(radius, (ny*2-1)*radius, ny)
        xx, yy = np.meshgrid(x, y)
        return cls(xx.flatten(), yy.flatten(), [radius]*xx.size)

    def to_vmd(self, filename):
        with open(filename+".pdb", "w") as pdbfile:
            pdbfile.write("CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf\n"
                      %(self.xhi-self.xlo, self.yhi-self.ylo, self.zhi-self.zlo, 90.0, 90.0, 90.0))
            rads = self.radius
            minrad = np.min(rads)
            maxrad = np.max(rads)
            drad = {}
            for id, type, x, y, rad in zip(self.id, self.type, self.x, self.y, self.radius):
                current_id = id
                if current_id >= 100000: 
                    current_id = current_id % 99999
                if maxrad == minrad: 
                    beta = 1.0            
                else:
                    beta = int((rad - minrad)/(maxrad - minrad)*100)/100;
                
                pdbfile.write("ATOM  %5d %4s %3s  %4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n"
                              %(current_id, "COL", "COL", type, x, y, 0.0, 0.0, beta))
                drad[beta] = rad  
            pdbfile.write("END\n")

        with open(filename+".vmd", "w") as vmdfile:
            vmdfile.write("mol new "+str(filename)+".pdb\n")
            vmdfile.write("mol delrep 0 top\n")
            for j,(k, v) in enumerate(drad.items()):
                vmdfile.write("set sel%d [atomselect top {beta %g}]\n" %(j,k))
                vmdfile.write("$sel%d set radius 1\n" %(j))
            for j,(k, v) in enumerate(drad.items()):
                vmdfile.write("mol representation VDW %g 10.0\n" %(v))
                vmdfile.write("mol selection [$sel%d text]\n" %(j))
                vmdfile.write("mol color ColorID 1\n")
                vmdfile.write("mol addrep top\n")
                
    def to_vtk(self, filename):
        with open(filename, "w") as f:
            f.write("# vtk DataFile Version 2.0\n")
            f.write("Time step 0\n")
            f.write("ASCII\n")
            f.write("DATASET UNSTRUCTURED_GRID\n")
            f.write("POINTS "+str(len(self.x))+" float\n")
            for x,y in zip(self.x, self.y):
                f.write(str(x)+" "+str(y)+" 0\n")
            f.write("POINT_DATA "+str(len(self.x))+"\n")
            f.write("SCALARS Diameter float 1\n")
            f.write("LOOKUP_TABLE default\n")
            for rad in self.radius:
                f.write(str(rad*2)+" ")
            f.write("\nSCALARS Type int 1\n")
            f.write("LOOKUP_TABLE default\n")
            for type in self.type:
                f.write(str(type)+" ")

    def to_particle_netcdf(self, filename):
        """
        Write `Particles` object to netcdf file
        
        Parameters
        ----------
        filename : string
            Name of netcdf file to write to. If a file by that name already exists,
            it will not be overwritten, and the function will return having written nothing.
        """
        if not nc_avail:
            print("Function requires netcdf python package, exiting")
            return
        
        fileOut = nc4.Dataset(filename,"w",format="NETCDF3_CLASSIC")
        
        nBonds = len(self.bonds)
        
        fileOut.nTypes = np.max(self.type)
        fileOut.maxRadius = np.max(self.radius)

        fileOut.createDimension("nParticles",self.num_particles)
        fileOut.createDimension("TWO",2)

        var = fileOut.createVariable("globalID","i",dimensions=["nParticles"])
        var[:] = self.id
        var.units = "-"
        
        var = fileOut.createVariable("type","i",dimensions=["nParticles"])
        var[:] = self.type
        var.units = "-"

        var = fileOut.createVariable("x","d",dimensions=["nParticles","TWO"])
        var[:,0] = self.x
        var[:,1] = self.y
        var.units = "m"

        var = fileOut.createVariable("radius","d",dimensions=["nParticles"])
        var[:] = self.radius
        var.units = "m"
        
        var = fileOut.createVariable("iceFraction","d",dimensions=["nParticles"])
        var[:] = self.iceFraction
        var.units = "-"

        var = fileOut.createVariable("iceThickness","d",dimensions=["nParticles"])
        var[:] = self.iceThickness
        var.units = "m"

        if nBonds > 0:            
            bondsNetcdf = np.array([[bond.id1, bond.id2] for bond in self.bonds]).astype(np.int)
            fileOut.createDimension("nBonds",len(self.bonds))            
            var = fileOut.createVariable("bonds","i",dimensions=["nBonds","TWO"])
            var[:] = bondsNetcdf
            var.units = "-"

        fileOut.close()


    def to_lammps_data(self, filename, comment="LAMMPS data file created with pydemsi\n"):
        """
        Write `Particles` object to LAMMPS data file
        
        Parameters
        ----------
        filename : string
            Name of data file to write to. If a file by that name already exists,
            it will be overwritten. 
        comment : string, optional
            Comment to append to top of lammps data file.
        """
        with open(filename, "w") as f:
            f.write(comment)
            f.write(str(self.num_particles)+" atoms\n")
            if len(self.bonds) > 0:
                f.write(str(len(self.bonds))+" bonds\n")
            f.write("\n")
            f.write(str(np.max(self.type))+" atom types\n")
            if len(self.bonds) > 0:
                f.write("1 bond types\n")
            f.write("\n")
            f.write(str(self.xlo)+" "+str(self.xhi)+" xlo xhi\n")
            f.write(str(self.ylo)+" "+str(self.yhi)+" ylo yhi\n")
            f.write(str(self.zlo)+" "+str(self.zhi)+" zlo zhi\n")
            f.write("\nAtoms\n\n")
            for id, type, rad, dens, x, y in zip(self.id, self.type, self.radius, self.density, self.x, self.y):
                f.write(str(id)+" "+str(type)+" "+str(rad)+" "+str(dens)+" "+
                        str(x)+" "+str(y)+" 0\n")
            if len(self.bonds) > 0:
                f.write("\nBonds\n\n")
                for b in self.bonds:
                    f.write(str(b.bondid)+" "+str(b.type)+" "+str(b.id1)+" "+str(b.id2)+"\n")

    def create_bonds(self, min_overlap = 0):
        """
        Create bonds between any particles with overlap >= min_overlap.

        Overlap is defined as :math:`R_1 + R_2 - \|\vec{r}_1 - \vec{r}_2\|`, where
        :math: `\vec{r}` is the particle position, :math:`R` is the particle radius.
        
        Typical use is min_overlap = 0, meaning any particles in contact will be bonded.

        Note, you can also set min_overlap < 0, in which case any particles separated by a distance 
        less than abs(min_overlap) will be bonded.
        
        Parameters:
        -----------
        min_overlap : float
            Minimum overlap, see above.

        """
        try:
            import scipy.spatial
        except:
            print("Error: bond search requires scipy.spatial, exiting")
            return

        points = np.vstack((self.x, self.y)).T 
        maxrad = np.max(self.radius)
        if min_overlap < 0:
            maxrad += np.abs(min_overlap)

        tree = scipy.spatial.KDTree(points)
        distance_matrix = tree.sparse_distance_matrix(tree, 2*maxrad)
        bonds = {(self.id[k[0]], self.id[k[1]]) : v for k,v in distance_matrix.items() 
                 if self.radius[k[0]] + self.radius[k[1]] - v >= min_overlap and k[0] < k[1] and self.type[k[0]] != 0 and self.type[k[1]] != 0}
        self.bonds = [Bond(i+1, k[0], k[1]) for i,k in enumerate(bonds.keys())]

    def get_bond_element_indices(self, min_overlap = 0):
        """
        Determine element types for bonds between any particles with overlap >= min_overlap.

        Overlap is defined as :math:`R_1 + R_2 - \|\vec{r}_1 - \vec{r}_2\|`, where
        :math: `\vec{r}` is the particle position, :math:`R` is the particle radius.

        Typical use is min_overlap = 0, meaning any particles in contact will be bonded.

        Note, you can also set min_overlap < 0, in which case any particles separated by a distance
        less than abs(min_overlap) will be bonded.

        Parameters:
        -----------
        min_overlap : float
            Minimum overlap, see above.

        """
        try:
            import scipy.spatial
        except:
            print("Error: bond search requires scipy.spatial, exiting")
            return

        points = np.vstack((self.x, self.y)).T
        maxrad = np.max(self.radius)
        if min_overlap < 0:
            maxrad += np.abs(min_overlap)

        tree = scipy.spatial.KDTree(points)
        distance_matrix = tree.sparse_distance_matrix(tree, 2*maxrad)
        bonds = {(k[0], k[1]) : v for k,v in distance_matrix.items()
                 if self.radius[k[0]] + self.radius[k[1]] - v >= min_overlap and k[0] < k[1]}
        return [Bond(i+1, k[0], k[1]) for i,k in enumerate(bonds.keys())]

    @classmethod
    def _new_from_indices(cls, obj, keep_indices):
        #Only keep bonds that have particles that are both in keep_indices
        newbonds = obj.bonds
        if len(keep_indices) == 0:
            print("Warning, all particles have been removed!")
            return
        
        keep_ids = obj.id[keep_indices]
        if len(obj.bonds) > 0:
            newbonds = [b for b in obj.bonds if b.id1 in keep_ids and b.id2 in keep_ids]
    
        return cls(obj.x[keep_indices], obj.y[keep_indices], obj.radius[keep_indices],
                     iceFraction = obj.iceFraction[keep_indices], iceThickness = obj.iceThickness[keep_indices],
                     density = obj.density[keep_indices], type = obj.type[keep_indices], id = obj.id[keep_indices],
                     xlo = obj.xlo, ylo = obj.ylo, xhi = obj.xhi, yhi = obj.yhi, zlo = obj.zlo, zhi = obj.zhi,
                     bonds = newbonds)
    
    def reorder_ids(self):
        """
        Set particle ID's to run from 1, N; change id's in bond entries as needed
        """
        new_ids = np.arange(1,len(self.id)+1)
        map_old_new = dict(zip(self.id, new_ids))
        for b in self.bonds:
            b.id1 = map_old_new[b.id1]
            b.id2 = map_old_new[b.id2]
        self.id = new_ids   
        
    def shift_particles(self, offset_x=0, offset_y=0):
        """
        Shift all particles such that minimum x extent is set to offset_x, 
        minimum y extent is set to offset_y, and box bounds correspond to the 
        min/max extent of any particle.
        """
        self.x = self.x - np.min(self.x - self.radius) + offset_x;
        self.y = self.y - np.min(self.y - self.radius) + offset_y;
        
        self.xlo = np.min(self.x - self.radius)
        self.ylo = np.min(self.y - self.radius)
        self.xhi = np.max(self.x + self.radius)
        self.yhi = np.max(self.y + self.radius)
    
    def types_from_image(self, imagefile, dx=1, dy=1, types=[], pixelvalues = []):
        """
        Modify particle types based on a grayscale image.

        The `types` and `pixelvalues` arrays are used to indicate which pixelvalues to
        map to which particle types; use a value of type <=0 to indicate a pixel
        value for which particles will be removed. In all cases, particles are mapped to pixel values
        based on the pixel that their center falls in.
        
        To align the image to the lammps data file, the bottom left corner of
        the image is assumed to be at xlo, ylo, and the horizontal axis corresponds to the
        x-direction. The input pixel size as given by `dx and  `dy`
        is then used to scale the image. Any particles outside the image are also removed.
        
        """

        try:
            img = scipy.ndimage.imread(imagefile).T.flipud() #Check this!
        except:
            print("Could not read image file "+imagefile)
            return
    
        xc = (self.x - self.xlo)//dx
        yc = (self.y - self.ylo)//dy
        keep_indices = [(xc > 0) & (yc > 0) & (xc <=img.shape[0]) & (yc <= img.shape[1])]
        
        
        values_at_centers = img[xc[keep_indices], yc[keep_indices]]
    
        for type, value in zip(types, values):
            indices = np.where(values_at_centers == value)[0].tolist()
            if type <= 0:
                keep_indices = keep_indices | np.logical_not(indices)
            else:
                self.type[indices] = type
        
        return self._new_from_indices(self, keep_indices)
        
        

    def remove_particles(self, remove_indices = [], remove_ids = []):
        """
        Remove particles with indices contained in `remove_indices`, or ids
        contained in `remove_ids`.
        
        Parameters
        ----------
        remove_indices : array of int
            Indices of particles to remove. These are any numbers of integers
            between 0..nParticles, corresponding to the indices of particle
            properties in the object (i.e. not the particle id's)
        remove_ids : array of int
            Ids of particles to remove, as defined in the `id` attribute. Both remove_indices
            and remove_ids can be specified, in which case particles that fit either 
            criterion will be removed
            
        
        """
        remove_indices = set(remove_indices).union(set(np.where(np.isin(self.id, remove_ids))[0]))
        keep_indices = set(np.arange(0,len(self.id))) - remove_indices
        inds = np.array(list(keep_indices)) #Convert to correct type for fancy indexing below
        return self._new_from_indices(self, inds), inds
        
    def plot(self, draw_bonds = False):
        """
        Plot particle structure using matplotlib.
        
        Parameters
        ----------
        draw_bonds : bool, optional
            If `True`, bonds will be drawn as lines between bonded particles.
            
        """
        import matplotlib.pyplot as plt
        import matplotlib.collections
        import matplotlib.patches
        fig, ax = plt.subplots()
        patches = [plt.Circle((x,y), 1.0*rad) for x, y, rad in zip(self.x, self.y, self.radius)]
        coll = matplotlib.collections.PatchCollection(patches, facecolors='black')
        ax.add_collection(coll)
        
        if draw_bonds:
            id_to_index = {id:index for index,id in enumerate(self.id)}
            xycoords = [((self.x[id_to_index[b.id1]], self.y[id_to_index[b.id1]]),
                         (self.x[id_to_index[b.id2]], self.y[id_to_index[b.id2]]))
                         for b in self.bonds]
            [ax.plot([xy[0][0], xy[1][0]],[xy[0][1], xy[1][1]],'-r', linewidth=3) for xy in xycoords]
            #bond_patches = [matplotlib.patches.ConnectionPatch(xy[0], xy[1], "data") for xy in xycoords]
            #print(bond_patches)
            #bond_coll = matplotlib.collections.PatchCollection(bond_patches)
            #ax.add_collection(bond_coll)
            
            
        ax.set_aspect('equal')
        ax.set_xlim([self.xlo, self.xhi])
        ax.set_ylim([self.ylo, self.yhi])
        plt.show()
        
    def __add__(self, other):
        """
        Combine two Particles objects into one. Also, any additional parameters added to each
        object, e.g. via the 'kwargs' argument to the constructor, will not be preserved in
        the combined object.
        """
        
        iceFraction = np.hstack((self.iceFraction, other.iceFraction))
        iceThickness = np.hstack((self.iceThickness, other.iceThickness))
        
        x = np.hstack((self.x, other.x))
        y = np.hstack((self.y, other.y))
        radius = np.hstack((self.radius, other.radius))
        type = np.hstack((self.type, other.type))
        density = np.hstack((self.density, other.density))
        
        #Offset id's
        offset = np.max(self.id)
        id = np.hstack((self.id, other.id + offset))
        bond_offset = 0
        if len(self.bonds) > 0:
            bond_offset = np.max([b.bondid for b in self.bonds])
        otherbonds = [Bond(b.bondid + bond_offset, b.id1 + offset, b.id2 + offset, type = b.type) for b in other.bonds]
        bonds = self.bonds + otherbonds
        
        return Particles(x, y, radius, bonds=bonds,
                         iceFraction = iceFraction, iceThickness = iceThickness,
                         density = density, type = type, id = id,
                         xlo = np.min((self.xlo, other.xlo)), xhi = np.max((self.xhi, other.xhi)),
                         ylo = np.min((self.ylo, other.ylo)), yhi = np.max((self.yhi, other.yhi)),
                         zlo = np.min((self.zlo, other.zlo)), zhi = np.max((self.zhi, other.zhi)),
                         xperiodic = self.xperiodic & other.xperiodic,
                         yperiodic = self.yperiodic & other.yperiodic)
