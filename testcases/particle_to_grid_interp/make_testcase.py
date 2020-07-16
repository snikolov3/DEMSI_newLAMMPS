#!/usr/bin/env python
import sys
import numpy as np
from netCDF4 import Dataset
import netCDF4
sys.path.append("../../utils/testcases/")
from testcase_utils import write_out_grid, write_particle_file

# grid
nx=11
ny=11

lx=500
ly=500

dx, dy = lx/(nx-1), ly/(ny-1)

lat = np.zeros((nx,ny))
lon = np.zeros((nx,ny))

write_out_grid("grid.nc", nx, ny, dx, dy, lat, lon)

# particles
nx = 20
ny = 20

nParticles = nx * ny

nTypes = 1

x = np.zeros(nParticles)
y = np.zeros(nParticles)
r = np.zeros(nParticles)
h = np.zeros(nParticles)
a = np.zeros(nParticles)
d = np.zeros(nParticles, dtype="i")
t = np.ones(nParticles, dtype="i")

radius = 12.5

ij = 0
for i in range(0,nx):
    for j in range(0,ny):

        x[ij] = (float(i ) + 0.5) * radius * 2.0
        y[ij] = (float(j ) + 0.5) * radius * 2.0

        r[ij] = radius
        h[ij] = 1.0
        a[ij] = 1.0

        d[ij] = ij + 1

        ij = ij + 1

filenameOut = "particles_in.nc"
write_particle_file(filenameOut, nParticles, d, nTypes, t, x, y, r, f=a, h=h)
