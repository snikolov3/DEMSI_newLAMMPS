import numpy as np
from netCDF4 import Dataset
import netCDF4
import matplotlib.pyplot as plt
from numpy.random import RandomState
import math
import sys
sys.path.append("../../utils/testcases/")
from testcase_utils import write_out_grid, write_out_fixed_forcing, write_out_time_varying_forcing, write_particle_file

#-------------------------------------------------------------------------------

def make_particles(filenameOut, random=False, initDistribution=False, hmin=0.8, hmax=1.0):

    # particles
    nx = 100
    ny = 50
    radius = 2500.0

    nParticles = nx * ny
    if (initDistribution): nParticles = nParticles * 2

    nTypes = 1

    x = np.zeros(nParticles)
    y = np.zeros(nParticles)
    r = np.zeros(nParticles)
    d = np.zeros(nParticles, dtype="i")
    t = np.zeros(nParticles, dtype="i")
    f = np.zeros(nParticles)
    h = np.zeros(nParticles)

    if (initDistribution):

        maxVertices = 4

        nv = np.zeros((nParticles), dtype="i")
        nv[:] = 4
        xv = np.zeros((nParticles,maxVertices), dtype="i")
        yv = np.zeros((nParticles,maxVertices), dtype="i")



    rndState = RandomState(1234)

    # upper real particles
    ij = 0
    for i in range(0,nx):
        for j in range(0,ny):

            x[ij] = (float(i     ) + 0.5) * radius * 2.0
            y[ij] = (float(j + ny) + 0.5) * radius * 2.0

            r[ij] = radius

            f[ij] = 1.0
            if (random):
                h[ij] = rndState.uniform(hmin,hmax)
            else:
                h[ij] = 1.0

            t[ij] = 1

            d[ij] = ij + 1

            if (initDistribution):

                xv[ij,0] = x[ij] - radius
                xv[ij,1] = x[ij] + radius
                xv[ij,2] = x[ij] + radius
                xv[ij,3] = x[ij] - radius

                yv[ij,0] = y[ij] - radius
                yv[ij,1] = y[ij] - radius
                yv[ij,2] = y[ij] + radius
                yv[ij,3] = y[ij] + radius

                # boundary polygons need to be larger for complete conservation in remapping
                if (i == 0):
                    xv[ij,0] = xv[ij,0] - 2500.0
                    xv[ij,3] = xv[ij,3] - 2500.0

                if (i == nx-1):
                    xv[ij,1] = xv[ij,1] + 2500.0
                    xv[ij,2] = xv[ij,2] + 2500.0

                if (j == ny-1):
                    yv[ij,2] = yv[ij,2] + 2500.0;
                    yv[ij,3] = yv[ij,3] + 2500.0;

            ij = ij + 1

    # lower init particles
    if (initDistribution):

        for i in range(0,nx):
            for j in range(0,ny):

                x[ij] = (float(i) + 0.5) * radius * 2.0
                y[ij] = (float(j) + 0.5) * radius * 2.0

                r[ij] = radius

                f[ij] = 0.0
                h[ij] = 0.0

                d[ij] = 0.0

                t[ij] = 0

                xv[ij,0] = x[ij] - radius
                xv[ij,1] = x[ij] + radius
                xv[ij,2] = x[ij] + radius
                xv[ij,3] = x[ij] - radius

                yv[ij,0] = y[ij] - radius
                yv[ij,1] = y[ij] - radius
                yv[ij,2] = y[ij] + radius
                yv[ij,3] = y[ij] + radius

                # boundary polygons need to be larger for complete conservation in remapping
                if (i == 0):
                    xv[ij,0] = xv[ij,0] - 2500.0
                    xv[ij,3] = xv[ij,3] - 2500.0

                if (i == nx-1):
                    xv[ij,1] = xv[ij,1] + 2500.0
                    xv[ij,2] = xv[ij,2] + 2500.0

                if (j == 0):
                    yv[ij,0] = yv[ij,0] - 2500.0;
                    yv[ij,1] = yv[ij,1] - 2500.0;

                ij = ij + 1

    if (initDistribution):
        write_particle_file(filenameOut, nParticles, d, nTypes, t, x, y, r, f, h, maxVertices, nv, xv, yv)
    else:
        write_particle_file(filenameOut, nParticles, d, nTypes, t, x, y, r, f, h)

#-------------------------------------------------------------------------------

# grid and forcing
nx=101
ny=101

lx=500000
ly=500000

dx, dy = lx/(nx-1), ly/(ny-1)

omega = 0.5e-3
lda = 8e5

xc, yc = 0.5*lx, 0.5*ly
yc = yc - 50000.0

ix, iy = np.indices((nx,ny))

rmag = np.sqrt((ix*dx - xc)**2 + (iy*dy-yc)**2)
rmag[rmag == 0] = 1e-10;

rxn = (ix*dx-xc)/rmag
ryn = (iy*dx-yc)/rmag
kx = -ryn
ky = rxn
prefac = np.minimum(omega*rmag, lda/rmag)
ugx = prefac*kx
ugy = prefac*ky

lat = np.zeros((nx,ny))
lon = np.zeros((nx,ny))

write_out_grid("grid.nc", nx, ny, dx, dy, lat, lon)

write_out_fixed_forcing("forcing_fixed.nc", nx, ny, ugx, ugy)

write_out_time_varying_forcing("forcing_varying.0001.nc", nx, ny, ugx, ugy)

make_particles("particles_in.nc")
make_particles("particles_in_random.nc", random=True)
make_particles("particles_in_random_init.nc", random=True, initDistribution=True)
make_particles("particles_in_stability.nc", random=True, initDistribution=False, hmin=0.0, hmax=1.0)
