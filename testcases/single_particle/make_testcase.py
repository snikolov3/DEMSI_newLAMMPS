import sys
import numpy as np
from netCDF4 import Dataset
import netCDF4
import matplotlib.pyplot as plt
sys.path.append("../../utils/testcases/")
from testcase_utils import write_out_grid, write_particle_file

#-------------------------------------------------------------------------------

def write_out_grid(filenameOut, nx, ny, dx, dy, lat, lon):

    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.createDimension("TWO",2)
    fileOut.createDimension("nx",nx)
    fileOut.createDimension("ny",ny)

    resolutionVar = fileOut.createVariable("resolution","d",dimensions=["TWO"])
    resolutionVar[:] = [dx,dy]

    dataVar = fileOut.createVariable("latitude","d",dimensions=["nx","ny"])
    dataVar[:] = lat[:]

    dataVar = fileOut.createVariable("longitude","d",dimensions=["nx","ny"])
    dataVar[:] = lon[:]

    fileOut.close()

#-------------------------------------------------------------------------------

def write_particle_file(filenameOut, nParticles, d, nTypes, t, x, y, r):

    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.maxRadius = np.amax(r)
    fileOut.nTypes = nTypes

    fileOut.createDimension("nParticles",nParticles)
    fileOut.createDimension("TWO",2)

    var = fileOut.createVariable("globalID","i",dimensions=["nParticles"])
    var.units = "-"
    var[:] = d[:]

    var = fileOut.createVariable("type","i",dimensions=["nParticles"])
    var.units = "-"
    var[:] = t[:]

    var = fileOut.createVariable("x","d",dimensions=["nParticles","TWO"])
    var.units = "m"
    var[:,0] = x[:]
    var[:,1] = y[:]

    var = fileOut.createVariable("radius","d",dimensions=["nParticles"])
    var.units = "m"
    var[:] = r[:]

    fileOut.close()

#-------------------------------------------------------------------------------

# grid
nx=11
ny=11

lx=1000000
ly=1000000

dx, dy = lx/(nx-1), ly/(ny-1)

lat = np.zeros((nx,ny))
lon = np.zeros((nx,ny))

write_out_grid("grid.nc", nx, ny, dx, dy, lat, lon)


# particles
nParticles = 1

nTypes = 1

x = np.zeros(nParticles)
y = np.zeros(nParticles)
r = np.zeros(nParticles)
d = np.zeros(nParticles, dtype="i")
t = np.ones(nParticles, dtype="i")

x[0] = 10000.0
y[0] = 500000.0
r[0] = 2500.0
d[0] = 1
t[0] = 1

filenameOut = "particles_in.nc"
write_particle_file(filenameOut, nParticles, d, nTypes, t, x, y, r)
