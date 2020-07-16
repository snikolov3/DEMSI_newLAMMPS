import sys
sys.path.append("../../utils/initialize_particles/")
import pydemsi
import numpy as np
from netCDF4 import Dataset
import math
import os

sys.path.append("../../utils/testcases/")
from testcase_utils import write_out_grid, plot_particles_debug

sys.path.append("../../utils/stereoprojection/")
from demsi_xy2geodetic import demsi_xy2geodetic

sys.path.append("../../utils/visualization/")
from make_particle_plot import plot_particles

import matplotlib.pyplot as plt

# data input location
dataDir = os.environ['DEMSI_DATA_DIR']

# compress simulation stats
initialWidth = 24000000
compressRatio = 0.4
totalTimesteps = 100000
selectedTimestep = 40000

# desired domain and grid
Lx = 7e6
Ly = 11e6

x1 = Lx * 0.5 + 300000.0
y1 = Ly * 0.5

theta = -25.0

nx = 100

# infer corrected domain values so pole is midway between grid nodes
dx = Lx / float(nx-1)
nx1 = int(math.floor(x1 / dx))
x1 = (float(nx1-1) + 0.5) * dx

dy = dx
ny1 = int(math.floor(y1 / dy))
y1 = (float(ny1-1) + 0.5) * dy
ny2 = int(math.ceil((Ly - y1) / dy))
ny = ny1 + ny2
Ly = dx * float(ny-1)




# calculate DEMSI grid file
xGRID = np.zeros((nx,ny))
yGRID = np.zeros((nx,ny))
for ix in range(0,nx):
    for iy in range(0,ny):
        xGRID[ix,iy] = ix * dx - x1
        yGRID[ix,iy] = iy * dy - y1

thetaRadians = np.radians(theta)
c, s = np.cos(thetaRadians), np.sin(thetaRadians)

xRotate = np.zeros((nx,ny))
yRotate = np.zeros((nx,ny))
for ix in range(0,nx):
    for iy in range(0,ny):
        xRotate[ix,iy] = c * xGRID[ix,iy] - s * yGRID[ix,iy]
        yRotate[ix,iy] = s * xGRID[ix,iy] + c * yGRID[ix,iy]

latGRID, lonGRID = demsi_xy2geodetic(xRotate, yRotate, 1)

filenameOut = "grid.nc"
write_out_grid(filenameOut, nx, ny, dx, dy, latGRID, lonGRID)



# final compress element distribution size
CompressSize = initialWidth * (1.0 - (2.0 * (1.0 - compressRatio) * selectedTimestep) / totalTimesteps)
print("CompressSize: ", CompressSize)

xMiddleCompress = initialWidth * 0.5
yMiddleCompress = initialWidth * 0.5



# load particles
particles = pydemsi.Particles.from_particle_netcdf("particles_in_expand.nc")


# shift to origin in middle
particles.x[:] = particles.x[:] - xMiddleCompress
particles.y[:] = particles.y[:] - yMiddleCompress



# remove elements exterior to desired domain
print("Remove elements outside desired domain")
xMin = -x1
xMax = Lx - x1
yMin = -y1
yMax = Ly - y1

particlesToRemove = []
for iParticle in range(0,particles.num_particles):

    if (particles.x[iParticle]-particles.radius[iParticle] <= xMin or particles.x[iParticle]+particles.radius[iParticle] >= xMax or
        particles.y[iParticle]-particles.radius[iParticle] <= yMin or particles.y[iParticle]+particles.radius[iParticle] >= yMax):

        particlesToRemove.append(particles.id[iParticle])

print("Particles to remove: %i" %(len(particlesToRemove)))

particles, keptIndices = particles.remove_particles(remove_ids = particlesToRemove)



# lat lon DEMSI
print("Get DEMSI lat lon")

thetaRadians = np.radians(theta)
c, s = np.cos(thetaRadians), np.sin(thetaRadians)

xRotate = np.zeros(particles.num_particles)
yRotate = np.zeros(particles.num_particles)
for iParticle in range(0,particles.num_particles):

    xRotate[iParticle] = c * particles.x[iParticle] - s * particles.y[iParticle]
    yRotate[iParticle] = s * particles.x[iParticle] + c * particles.y[iParticle]

latDEMSI, lonDEMSI = demsi_xy2geodetic(xRotate, yRotate, 1)







# determine whether on land or not

# load ETOPO data
print("Load ETOPO1 data")
filenameETOPO1 = dataDir + "/etopo1/ETOPO1_Ice_g_gmt4.grd"
print(filenameETOPO1)

fileETOPO1 = Dataset(filenameETOPO1, "r")

nxETOPO1 = len(fileETOPO1.dimensions["x"])
nyETOPO1 = len(fileETOPO1.dimensions["y"])

lonETOPO1 = fileETOPO1.variables["x"][:]
latETOPO1 = fileETOPO1.variables["y"][:]
zETOPO1 = fileETOPO1.variables["z"][:]

lonETOPO1[:] = np.radians(lonETOPO1[:])
latETOPO1[:] = np.radians(latETOPO1[:])

fileETOPO1.close()







# interpolate ETOPO1
print("Interpolate ETOPO1")
dLonLat = 0.00029088820867

for iParticle in range(0,particles.num_particles):

    iLon1 = int(math.floor((lonDEMSI[iParticle] + math.pi) / dLonLat))
    iLon2 = iLon1 + 1

    wLon1 = (lonETOPO1[iLon2]    - lonDEMSI[iParticle] ) / (lonETOPO1[iLon2] - lonETOPO1[iLon1])
    wLon2 = (lonDEMSI[iParticle] - lonETOPO1[iLon1]    ) / (lonETOPO1[iLon2] - lonETOPO1[iLon1])

    iLat1 = int(math.floor((latDEMSI[iParticle] + 0.5*math.pi) / dLonLat))
    iLat2 = iLat1 + 1

    wLat1 = (latETOPO1[iLat2]    - latDEMSI[iParticle] ) / (latETOPO1[iLat2] - latETOPO1[iLat1])
    wLat2 = (latDEMSI[iParticle] - latETOPO1[iLat1]    ) / (latETOPO1[iLat2] - latETOPO1[iLat1])

    zInterp = zETOPO1[iLat1,iLon1] * wLat1 * wLon1 + \
              zETOPO1[iLat1,iLon2] * wLat1 * wLon2 + \
              zETOPO1[iLat2,iLon1] * wLat2 * wLon1 + \
              zETOPO1[iLat2,iLon2] * wLat2 * wLon2

    if (zInterp > 0.0):
        particles.type[iParticle] = 2
    else:
        particles.type[iParticle] = 1



# MPAS ice concentration and thickness
print("Load MPAS grid")
filenameMPAS = dataDir + "/mpas-seaice/QU60km_polar/seaice_QU_60km_polar.nc"
fileMPAS = Dataset(filenameMPAS, "r")

nCells = len(fileMPAS.dimensions["nCells"])

latCell = fileMPAS.variables["latCell"][:]
lonCell = fileMPAS.variables["lonCell"][:]

fileMPAS.close()

xCell = np.zeros(nCells)
yCell = np.zeros(nCells)
zCell = np.zeros(nCells)
for iCell in range(0,nCells):
    xCell[iCell] = math.cos(latCell[iCell]) * math.sin(lonCell[iCell])
    yCell[iCell] = math.cos(latCell[iCell]) * math.cos(lonCell[iCell])
    zCell[iCell] = math.sin(latCell[iCell])



print("Load MPAS data")
filenameMPAS = dataDir + "/mpas-seaice/QU60km_polar/timeSeriesStatsMonthly.2000-01.nc"
fileMPAS = Dataset(filenameMPAS, "r")

iceAreaCell   = fileMPAS.variables["timeMonthly_avg_iceAreaCell"][0,:]
iceVolumeCell = fileMPAS.variables["timeMonthly_avg_iceVolumeCell"][0,:]

fileMPAS.close()

xDEMSI = np.zeros(nCells)
yDEMSI = np.zeros(nCells)
zDEMSI = np.zeros(nCells)
for iParticle in range(0,particles.num_particles):
    xDEMSI[iParticle] = math.cos(latDEMSI[iParticle]) * math.sin(lonDEMSI[iParticle])
    yDEMSI[iParticle] = math.cos(latDEMSI[iParticle]) * math.cos(lonDEMSI[iParticle])
    zDEMSI[iParticle] = math.sin(latDEMSI[iParticle])


print("Find closest MPAS points")
for iParticle in range(0,particles.num_particles):

    if (particles.type[iParticle] == 1):

        dx = xCell[:]-xDEMSI[iParticle]
        dy = yCell[:]-yDEMSI[iParticle]
        dz = zCell[:]-zDEMSI[iParticle]

        distances = np.add(np.multiply(dx, dx), np.add(np.multiply(dy, dy), np.multiply(dz, dz)))
        iCellMin = np.argmin(distances)

        particles.iceFraction [iParticle] = iceAreaCell[iCellMin]
        particles.iceThickness[iParticle] = iceVolumeCell[iCellMin]

    else:

        particles.iceFraction [iParticle] = 1.0
        particles.iceThickness[iParticle] = 100.0







# get bond connectivity
print("Get bond connectivity info")
bondElementIndices = particles.get_bond_element_indices()
bondElementIndices1 = np.zeros((len(bondElementIndices)),dtype="i")
bondElementIndices2 = np.zeros((len(bondElementIndices)),dtype="i")
bondElementTypes1 = np.zeros((len(bondElementIndices)),dtype="i")
bondElementTypes2 = np.zeros((len(bondElementIndices)),dtype="i")

for i in range(0,len(bondElementIndices)):
    bondElementIndices1[i] = bondElementIndices[i].id1
    bondElementIndices2[i] = bondElementIndices[i].id2



# Reduce radius so no overlap - only needed temporarily
print("Get max bond overlaps")
maxOverlap = np.zeros(particles.num_particles)
for bondElementIndex1, bondElementIndex2 in zip(bondElementIndices1, bondElementIndices2):
    separation = math.sqrt(math.pow(particles.x[bondElementIndex1]-particles.x[bondElementIndex2],2) + \
                           math.pow(particles.y[bondElementIndex1]-particles.y[bondElementIndex2],2))
    overlap = particles.radius[bondElementIndex1] + particles.radius[bondElementIndex2] - separation
    maxOverlap[bondElementIndex1] = max(maxOverlap[bondElementIndex1], overlap)
    maxOverlap[bondElementIndex2] = max(maxOverlap[bondElementIndex2], overlap)

print("Reduce radius by half max overlap")
for iParticle in range(0,particles.num_particles):
    particles.radius[iParticle] = particles.radius[iParticle] - 0.6 * maxOverlap[iParticle]



# remove seas
print("find closest element to north pole")
distance = np.add(np.multiply(particles.x,particles.x),np.multiply(particles.y,particles.y))
iParticleStart = np.argmin(distance)


print("Flood fill Arctic ocean")
nextParticles = np.array([iParticleStart])

while (len(nextParticles) > 0):

    particles.type[nextParticles] = -1
    bondElementTypes1 = particles.type[bondElementIndices1]
    bondElementTypes2 = particles.type[bondElementIndices2]

    nextParticles1 = np.where((np.isin(bondElementIndices1,nextParticles)) & (bondElementTypes2 == 1))[0]
    nextParticles1 = bondElementIndices2[nextParticles1]

    nextParticles2 = np.where((np.isin(bondElementIndices2,nextParticles)) & (bondElementTypes1 == 1))[0]
    nextParticles2 = bondElementIndices1[nextParticles2]

    nextParticles = np.array(list(set(np.concatenate((nextParticles1,nextParticles2)))))




print("Discard non-contiguous ocean points")
particles.type[np.where(particles.type ==  1)[0]] = 2
particles.type[np.where(particles.type == -1)[0]] = 1




print("Set coastal elements")
nCoastLayers = 2

type1 = 1
type2 = 2

for iCoastLayer in range(0,nCoastLayers):

    bondElementTypes1 = particles.type[bondElementIndices1]
    bondElementTypes2 = particles.type[bondElementIndices2]

    coastIndices1 = np.where((bondElementTypes1 == type1) & (bondElementTypes2 == type2))
    coastIndices1 = bondElementIndices2[coastIndices1]
    coastIndices2 = np.where((bondElementTypes1 == type2) & (bondElementTypes2 == type1))
    coastIndices2 = bondElementIndices1[coastIndices2]
    coastIndices = np.concatenate((coastIndices1,coastIndices2))

    particles.type[coastIndices] = 3 + iCoastLayer

    type1 = 3 + iCoastLayer





# remove interior land
print("Remove interior land")
indicesInteriorLand = np.where(particles.type == 2)[0]
particles, keptIndices = particles.remove_particles(remove_indices = indicesInteriorLand)
particles.type = np.where(particles.type > 1, 2, 1)

# cull elements with no sea ice area
print("Cull elements with no sea ice area")
elementsToCull = np.where((particles.iceFraction == 0.0) & (particles.type == 1))[0]
particles, keptIndices = particles.remove_particles(remove_indices = elementsToCull)





# final adjustment of particle positions for all positive DEMSI domain
particles.x[:] = particles.x[:] + x1
particles.y[:] = particles.y[:] + y1


particles.iceFraction[:] = 1.0
particles.iceThickness[:] = 1.0




# output
print("Write out particle file")
particles.bonds = []
particles.to_particle_netcdf("particles_in.nc")

print("Make particle file plots")
xMin = np.amin(np.subtract(particles.x,particles.radius))
xMax = np.amax(np.add(particles.x,particles.radius))
yMin = np.amin(np.subtract(particles.y,particles.radius))
yMax = np.amax(np.add(particles.y,particles.radius))

argsPlot = {"filenameIn":"particles_in_init.nc", "filenameOut":"particles_type.png", "xplotmin":xMin, "xplotmax":xMax, "yplotmin":yMin, "yplotmax":yMax, "gridFilename":None, "plotType":True, "useLegend":True, "removeTicks":False, "plotBonds":False, "varname":None, "cmin":None, "cmax":None, "coastline":False}
plot_particles(argsPlot)

argsPlot = {"filenameIn":"particles_in_init.nc", "filenameOut":"particles_iceFraction.png", "xplotmin":xMin, "xplotmax":xMax, "yplotmin":yMin, "yplotmax":yMax, "gridFilename":None, "plotType":False, "useLegend":False, "removeTicks":False, "plotBonds":False, "varname":"iceFraction", "cmin":None, "cmax":None, "coastline":False}
plot_particles(argsPlot)

argsPlot = {"filenameIn":"particles_in_init.nc", "filenameOut":"particles_iceThickness.png", "xplotmin":xMin, "xplotmax":xMax, "yplotmin":yMin, "yplotmax":yMax, "gridFilename":None, "plotType":False, "useLegend":False, "removeTicks":False, "plotBonds":False, "varname":"iceThickness", "cmin":None, "cmax":None, "coastline":False}
plot_particles(argsPlot)

argsPlot = {"filenameIn":"particles_in_init.nc", "filenameOut":"particles_iceFraction2.png", "xplotmin":None, "xplotmax":None, "yplotmin":None, "yplotmax":None, "gridFilename":"grid.nc", "plotType":False, "useLegend":False, "removeTicks":False, "plotBonds":False, "varname":"iceFraction", "cmin":None, "cmax":None, "coastline":False}
plot_particles(argsPlot)
