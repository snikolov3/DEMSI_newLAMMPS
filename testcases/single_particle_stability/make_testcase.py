import sys
import numpy as np
from netCDF4 import Dataset
import netCDF4
import matplotlib.pyplot as plt
import math
sys.path.append("../../utils/testcases/")
from testcase_utils import write_out_grid, write_out_fixed_forcing, write_particle_file

#-------------------------------------------------------------------------------

def make_particles(filenameOut, radius, thickness):

    # particles
    nParticles = 1

    nTypes = 1

    x = np.zeros(nParticles)
    y = np.zeros(nParticles)
    r = np.zeros(nParticles)
    d = np.zeros(nParticles, dtype="i")
    t = np.ones(nParticles, dtype="i")
    f = np.zeros(nParticles)
    h = np.zeros(nParticles)

    x[0] = 2.0 * radius
    y[0] = 2.0 * radius

    r[0] = radius

    f[0] = 1.0
    h[0] = thickness

    d[0] = 1

    write_particle_file(filenameOut, nParticles, d, nTypes, t, x, y, r, f=f, h=h)

#-------------------------------------------------------------------------------

# grid and forcing
area = 1e5
radius = math.sqrt(area / math.pi)

nx=101
ny=2


ly = 4.0*radius
lx = (ly / (ny-1)) * (nx-1)


dx, dy = lx/(nx-1), ly/(ny-1)


ugx = np.ones((nx,ny)) * 0.1
ugy = np.zeros((nx,ny)) * 0.1


lat = np.zeros((nx,ny))
lon = np.zeros((nx,ny))

write_out_grid("grid.nc", nx, ny, dx, dy, lat, lon)

write_out_fixed_forcing("forcing_fixed.nc", nx, ny, ugx, ugy)

make_particles("particles_in_thick.nc",radius,1.0)
make_particles("particles_in_thin.nc",radius,0.001)
