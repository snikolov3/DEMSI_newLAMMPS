import sys
sys.path.append("../../utils/initialize_particles/")
import pydemsi
import numpy as np
from netCDF4 import Dataset
import math
import os
sys.path.append("../../utils/testcases/")
from testcase_utils import plot_particles_debug, write_out_grid

import matplotlib.pyplot as plt

# compress simulation stats
initialWidth = 1000000

xMiddleCompress = initialWidth * 0.5
yMiddleCompress = initialWidth * 0.5

# load particles
particles = pydemsi.Particles.from_particle_netcdf("particles_in_expand.nc")


# shift to origin in middle
particles.x[:] = particles.x[:] - xMiddleCompress
particles.y[:] = particles.y[:] - yMiddleCompress


# remove elements exterior to desired domain
print("Remove elements outside desired domain")
xMin = -200000
xMax =  200000
yMin = -200000
yMax =  200000

particlesToRemove = []
for iParticle in range(0,particles.num_particles):

    if (particles.x[iParticle]-particles.radius[iParticle] <= xMin or particles.x[iParticle]+particles.radius[iParticle] >= xMax or
        particles.y[iParticle]-particles.radius[iParticle] <= yMin or particles.y[iParticle]+particles.radius[iParticle] >= yMax):

        particlesToRemove.append(particles.id[iParticle])




print("Particles to remove: %i" %(len(particlesToRemove)))
particles, keptIndices = particles.remove_particles(remove_ids = particlesToRemove)




#output
particles.x[:] = particles.x[:] + 250000.0
particles.y[:] = particles.y[:] + 250000.0
particles.type[:] = 1

particles.create_bonds(min_overlap = -1000.0)
particles.to_particle_netcdf("particles_in.nc")




# grid
nx = 11
ny = 11
dx = 50000
dy = 50000
lat = np.zeros((nx,ny))
lon = np.zeros((nx,ny))


filenameOut = "grid.nc"
write_out_grid(filenameOut, nx, ny, dx, dy, lat, lon)
