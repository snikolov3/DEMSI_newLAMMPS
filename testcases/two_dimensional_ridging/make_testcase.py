from __future__ import print_function

import sys
sys.path.append("../../utils/initialize_particles/")
import pydemsi
sys.path.append("../../utils/testcases/")
from testcase_utils import write_out_grid, write_out_fixed_forcing, add_init_polygons, add_init_connectivity
import numpy as np
from netCDF4 import Dataset
import argparse
import os
from radical_voronoi_tessellation import radical_voronoi_tessellation

# arguments
parser = argparse.ArgumentParser(description='Create Arctic Basin test case initial particle file.')

parser.add_argument('-i', dest='filenameIn', help='Input particle distribution filename.', required=True)

args = parser.parse_args()


# particle domain size
filein = Dataset(args.filenameIn,"r")
xRange = filein.xRange
yRange = filein.yRange
filein.close()


    
# grid
Lx = xRange
Ly = yRange
    
nx = 101
ny = 101
dx = Lx / float(nx-1)
dy = Ly / float(ny-1)

latGRID = np.zeros((nx,ny))
lonGRID = np.zeros((nx,ny))

filenameOut = "grid.nc"
write_out_grid(filenameOut, nx, ny, dx, dy, latGRID, lonGRID)





# load particles
particles = pydemsi.Particles.from_particle_netcdf(args.filenameIn)

# set element types
particles.type[:] = 1

# set coastal particles
for iParticle in range(0,len(particles.type)):
    if (particles.x[iParticle] > 1000000.0):
        particles.type[iParticle] = 2

# remove excess
#removeIndices = []
#for iParticle in range(0,len(particles.type)):
#    if (particles.x[iParticle] > 1000000.0+30000.0):
#        removeIndices.append(iParticle)
#print("Particles to cull: ", len(removeIndices))
#particles.remove_particles(remove_indices = removeIndices)


# write out
particles.to_particle_netcdf("particles_in_init.nc")


# perform radical voronoi tessellation of particle distribution
voroExecutable = os.environ['DEMSI_VORO_EXE']
radical_voronoi_tessellation("particles_in_init.nc", "radical_voronoi.tmp.nc", voroExecutable, 0.0, Lx, 0.0, Ly)

# load the tessellation
filePolygons = Dataset("radical_voronoi.tmp.nc","r")

nParticlesPolygon = len(filePolygons.dimensions["nParticles"])
maxVertices = len(filePolygons.dimensions["maxVertices"])

globalIDPolygon = filePolygons.variables["globalID"][:]
nVerticesOnCell = filePolygons.variables["nVerticesOnCell"][:]
xVertex = filePolygons.variables["xVertex"][:]
yVertex = filePolygons.variables["yVertex"][:]

filePolygons.close()


# add initial particles
maxCellSize = 0.0
add_init_polygons("particles_in_init.nc", maxCellSize, maxVertices, nVerticesOnCell, xVertex, yVertex)

# add init connectivity info
add_init_connectivity("particles_in_init.nc")


# forcing
ugx = np.zeros((nx,ny))
ugy = np.zeros((nx,ny))
ugx[:] = 10.0
ugy[:] =  0.0
write_out_fixed_forcing("forcing_fixed.nc", nx, ny, ugx, ugy)
