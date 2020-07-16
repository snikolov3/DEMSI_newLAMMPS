import sys
sys.path.append("../../utils/initialize_particles/")
import pydemsi
import numpy as np
from netCDF4 import Dataset
import netCDF4
sys.path.append("../../utils/testcases/")
from testcase_utils import write_out_grid


# grid definition
lxGrid = 30000
lyGrid = 30000


# particles
particlesRadius = 2500.0
particlesNx = 2
particlesNy = 1


# bonded particles
particlesBonded = pydemsi.Particles.from_lattice(particlesNx,particlesNy,particlesRadius)
particlesBonded.create_bonds()

particlesBonded.type[0] = 1
particlesBonded.type[1] = 2

xOffset = 0.5 * (lxGrid - 4.0 * particlesRadius)
yOffset = 0.5 * lyGrid - particlesRadius

particlesBonded.shift_particles(offset_x = xOffset, offset_y = yOffset)

particlesBonded.to_particle_netcdf("particles_in_bonded.nc")

# unbonded particles
particlesUnbonded = pydemsi.Particles.from_lattice(particlesNx,particlesNy,particlesRadius)

particlesUnbonded.type[0] = 1
particlesUnbonded.type[1] = 2

xOffset = 0.5 * (lxGrid - 4.0 * particlesRadius)
yOffset = 0.5 * lyGrid - particlesRadius

particlesUnbonded.shift_particles(offset_x = xOffset, offset_y = yOffset)

particlesUnbonded.to_particle_netcdf("particles_in_unbonded.nc")


# grid
nxGrid = 2
nyGrid = 2

dxGrid, dyGrid = lxGrid/(nxGrid-1), lyGrid/(nyGrid-1)

latGrid = np.zeros((nxGrid,nyGrid))
lonGrid = np.zeros((nxGrid,nyGrid))

write_out_grid("grid.nc", nxGrid, nyGrid, dxGrid, dyGrid, latGrid, lonGrid)
