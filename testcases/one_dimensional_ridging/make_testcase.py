import sys
sys.path.append("../../utils/initialize_particles/")
import pydemsi
import numpy as np
from netCDF4 import Dataset
import netCDF4
import math
sys.path.append("../../utils/testcases/")
from testcase_utils import write_out_grid, write_out_fixed_forcing

# grid
nx=102
ny=2

lx=1010000
ly=10000

dx, dy = lx/(nx-1), ly/(ny-1)

lat = np.zeros((nx,ny))
lon = np.zeros((nx,ny))

write_out_grid("grid.nc", nx, ny, dx, dy, lat, lon)

# forcing
ugx = np.zeros((nx,ny))
ugy = np.zeros((nx,ny))
ugx[:] = 10.0

write_out_fixed_forcing("forcing_fixed.nc", nx, ny, ugx, ugy)

# initial particles
radius = 5000.0
nx = 101
ny = 1
nParticles = nx * ny
lattice = pydemsi.Particles.from_lattice(nx,ny,radius)
lattice.type[:] = 1
lattice.type[-1] = 2 # final coastal point
lattice.create_bonds()

lattice.to_particle_netcdf("particles_in.nc")

# add polygon to output file
maxCellSize = math.sqrt(2.0) * 2.0 * radius

maxVertices = 4
nVerticesOnCell = np.zeros(nParticles,dtype="i")
xVertex = np.zeros((nParticles,maxVertices))
yVertex = np.zeros((nParticles,maxVertices))

nVerticesOnCell[:] = 4

for iParticle in range(0,nParticles):
    xVertex[iParticle,0] = lattice.x[iParticle] - radius
    yVertex[iParticle,0] = lattice.y[iParticle] - radius
    xVertex[iParticle,1] = lattice.x[iParticle] + radius
    yVertex[iParticle,1] = lattice.y[iParticle] - radius
    xVertex[iParticle,2] = lattice.x[iParticle] + radius
    yVertex[iParticle,2] = lattice.y[iParticle] + radius
    xVertex[iParticle,3] = lattice.x[iParticle] - radius
    yVertex[iParticle,3] = lattice.y[iParticle] + radius

nEdges = nParticles - 1
xEdges = np.zeros(nEdges)
yEdges = np.zeros(nEdges)
cellsOnEdge = np.zeros((nEdges,2),dtype="i")
edgesOnCell = np.zeros((nParticles,maxVertices),dtype="i")
cellsOnCell = np.zeros((nParticles,maxVertices),dtype="i")

cellsOnEdge[:] = -1
edgesOnCell[:] = -1
cellsOnCell[:] = -1

cellsOnCell[0,1] = 1
for iParticle in range(1,nParticles-1):
    cellsOnCell[iParticle,0] = iParticle-1
    cellsOnCell[iParticle,1] = iParticle+1
cellsOnCell[nParticles-1,0] = nParticles-2

edgesOnCell[0,1] = 0
for iParticle in range(1,nParticles-1):
    edgesOnCell[iParticle,0] = iParticle-1
    edgesOnCell[iParticle,1] = iParticle
edgesOnCell[nParticles-1,0] = nParticles-2

fileParticlesOut = Dataset("particles_in.nc","a")

fileParticlesOut.maxCellSize = maxCellSize

fileParticlesOut.createDimension("maxVertices",maxVertices)
fileParticlesOut.createDimension("nEdges",nEdges)
#fileParticlesOut.createDimension("TWO",2)

var = fileParticlesOut.createVariable("nVerticesOnCell","i",dimensions=["nParticles"])
var[:] = nVerticesOnCell[:]
var.units = "-"

var = fileParticlesOut.createVariable("xVertex","d",dimensions=["nParticles","maxVertices"])
var[:] = xVertex[:]
var.units = "m"

var = fileParticlesOut.createVariable("yVertex","d",dimensions=["nParticles","maxVertices"])
var[:] = yVertex[:]
var.units = "m"

var = fileParticlesOut.createVariable("xEdge","d",dimensions=["nEdges"])
var[:] = xEdges[:]
var.units = "m"

var = fileParticlesOut.createVariable("yEdge","d",dimensions=["nEdges"])
var[:] = yEdges[:]
var.units = "m"

var = fileParticlesOut.createVariable("cellsOnEdge","i",dimensions=["nEdges","TWO"])
var[:] = cellsOnEdge[:]
var.units = "-"

var = fileParticlesOut.createVariable("edgesOnCell","i",dimensions=["nParticles","maxVertices"])
var[:] = edgesOnCell[:]
var.units = "-"

var = fileParticlesOut.createVariable("cellsOnCell","i",dimensions=["nParticles","maxVertices"])
var[:] = cellsOnCell[:]
var.units = "-"

fileParticlesOut.close()
