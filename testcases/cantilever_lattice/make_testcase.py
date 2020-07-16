import numpy as np
from netCDF4 import Dataset
import netCDF4
import sys

sys.path.append("../../utils/initialize_particles/")
import pydemsi

sys.path.append("../../utils/testcases/")
from testcase_utils import write_out_grid

# grid
nx=101
ny=101

lx=200000
ly=200000

dx, dy = lx/(nx-1), ly/(ny-1)

lat = np.zeros((nx,ny))
lon = np.zeros((nx,ny))

write_out_grid("grid.nc", nx, ny, dx, dy, lat, lon)


# particles
lattice = pydemsi.Particles.from_lattice(30,3,2500)
lattice.create_bonds()

#Define different types for particles near edges, to flag in lammps script for 'pulling'
#(could also define these in lammps based on regions)
lattice.type[np.where(lattice.x <= np.min(lattice.x)+3000)[0]] = 4
lattice.type[np.where(lattice.x >= np.max(lattice.x)-3000)[0]] = 5

#Shift positions in both directions
lattice.shift_particles(offset_x = 25000, offset_y = 92500)
#lattice.plot(draw_bonds = True)

lattice.to_particle_netcdf("particles_in.nc")
