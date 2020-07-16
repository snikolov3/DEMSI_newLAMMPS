import sys
sys.path.append("../../utils/initialize_particles/")
import pydemsi
from radical_voronoi_tessellation import radical_voronoi_tessellation
import numpy as np
from netCDF4 import Dataset
import math
import os
import argparse

sys.path.append("../../utils/testcases/")
from testcase_utils import write_out_grid, plot_particles_debug, add_init_connectivity, add_init_polygons

sys.path.append("../../utils/stereoprojection/")
from demsi_xy2geodetic import demsi_xy2geodetic

sys.path.append("../../utils/visualization/")
from make_particle_plot import plot_particles

import matplotlib.pyplot as plt


# arguments
parser = argparse.ArgumentParser(description='Create Arctic Basin test case initial particle file.')

parser.add_argument('-i', dest='filenameIn', help='Input particle distribution filename.', required=True)
parser.add_argument('-w', dest='initialWidth', help='Initial width of the input particle distribution.', required=True, type=float)

args = parser.parse_args()


# Compressed state
# python make_testcase_init.py -i mypack_compress_step40000.nc -w 24000000.0
#filenameIn = "mypack_compress_step40000.nc"
#initialWidth = 24000000.0

# Lloyds algorithm
# python make_testcase_init.py -i particles_in_replicate.nc -w 18000000.0
#filenameIn = "particles_in_replicate.nc"
#initialWidth = 18000000.0




# data input location
dataDir = os.environ['DEMSI_DATA_DIR']


# middle of the input distribution
xMiddleCompress = args.initialWidth * 0.5
yMiddleCompress = args.initialWidth * 0.5




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







# load particles
particles = pydemsi.Particles.from_particle_netcdf(args.filenameIn)


# perform radical voronoi tessellation of particle distribution
print("Get radical voronoi tessellation")
voroExecutable = os.environ['DEMSI_VORO_EXE']
radical_voronoi_tessellation(args.filenameIn, "radical_voronoi.tmp.nc", voroExecutable, 0.0, args.initialWidth, 0.0, args.initialWidth)

# load the tessellation
filePolygons = Dataset("radical_voronoi.tmp.nc","r")

nParticlesPolygon = len(filePolygons.dimensions["nParticles"])
maxVertices = len(filePolygons.dimensions["maxVertices"])

globalIDPolygon = filePolygons.variables["globalID"][:]
nVerticesOnCell = filePolygons.variables["nVerticesOnCell"][:]
xVertex = filePolygons.variables["xVertex"][:]
yVertex = filePolygons.variables["yVertex"][:]

filePolygons.close()





# shift to origin in middle
particles.x[:] = particles.x[:] - xMiddleCompress
particles.y[:] = particles.y[:] - yMiddleCompress

xVertex[:] = xVertex[:] - xMiddleCompress
yVertex[:] = yVertex[:] - yMiddleCompress



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

particles, keptIndices1 = particles.remove_particles(remove_ids = particlesToRemove)



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

xDEMSI = np.zeros(particles.num_particles)
yDEMSI = np.zeros(particles.num_particles)
zDEMSI = np.zeros(particles.num_particles)
for iParticle in range(0,particles.num_particles):
    xDEMSI[iParticle] = math.cos(latDEMSI[iParticle]) * math.sin(lonDEMSI[iParticle])
    yDEMSI[iParticle] = math.cos(latDEMSI[iParticle]) * math.cos(lonDEMSI[iParticle])
    zDEMSI[iParticle] = math.sin(latDEMSI[iParticle])


print("Find closest MPAS points")
for iParticle in range(0,particles.num_particles):

    dx = xCell[:]-xDEMSI[iParticle]
    dy = yCell[:]-yDEMSI[iParticle]
    dz = zCell[:]-zDEMSI[iParticle]

    distances = np.add(np.multiply(dx, dx), np.add(np.multiply(dy, dy), np.multiply(dz, dz)))
    iCellMin = np.argmin(distances)

    particles.iceFraction [iParticle] = iceAreaCell[iCellMin]
    particles.iceThickness[iParticle] = min(iceVolumeCell[iCellMin],5.0)







# get bond connectivity
print("Get bond connectivity info")
bondElementIndices = particles.get_bond_element_indices(min_overlap = -10000.0)
bondElementIndices1 = np.zeros((len(bondElementIndices)),dtype="i")
bondElementIndices2 = np.zeros((len(bondElementIndices)),dtype="i")
bondElementTypes1 = np.zeros((len(bondElementIndices)),dtype="i")
bondElementTypes2 = np.zeros((len(bondElementIndices)),dtype="i")

for i in range(0,len(bondElementIndices)):
    bondElementIndices1[i] = bondElementIndices[i].id1
    bondElementIndices2[i] = bondElementIndices[i].id2



# Reduce radius so no overlap - only needed temporarily
#print("Get max bond overlaps")
#maxOverlap = np.zeros(particles.num_particles)
#for bondElementIndex1, bondElementIndex2 in zip(bondElementIndices1, bondElementIndices2):
#    separation = math.sqrt(math.pow(particles.x[bondElementIndex1]-particles.x[bondElementIndex2],2) + \
#                           math.pow(particles.y[bondElementIndex1]-particles.y[bondElementIndex2],2))
#    overlap = particles.radius[bondElementIndex1] + particles.radius[bondElementIndex2] - separation
#    maxOverlap[bondElementIndex1] = max(maxOverlap[bondElementIndex1], overlap)
#    maxOverlap[bondElementIndex2] = max(maxOverlap[bondElementIndex2], overlap)

#print("Reduce radius by half max overlap")
#for iParticle in range(0,particles.num_particles):
#    particles.radius[iParticle] = particles.radius[iParticle] - 0.6 * maxOverlap[iParticle]



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
nCoastLayers = 3

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
particles, keptIndices2 = particles.remove_particles(remove_indices = indicesInteriorLand)
particles.type = np.where(particles.type > 1, 2, 1)

# set elements with no sea ice to type 0
print("Set elements with no sea ice to type 0")
elementsToSetType = np.where((particles.iceFraction == 0.0) & (particles.type == 1))[0]
particles.type[elementsToSetType] = 0



# set coastal points thickness and fraction
for iParticle in range(0,particles.num_particles):

    if (particles.type[iParticle] == 2):

        particles.iceFraction [iParticle] = 1.0
        particles.iceThickness[iParticle] = 100.0





# final adjustment of particle positions for all positive DEMSI domain
particles.x[:] = particles.x[:] + x1
particles.y[:] = particles.y[:] + y1

xVertex[:] = xVertex[:] + x1
yVertex[:] = yVertex[:] + y1


#particles.iceFraction[:] = 1.0
#particles.iceThickness[:] = 1.0




# get polygons for reduced particle set
print("Get mapping between polygon particles and culled particles")
polygonMapping = {}
for iParticle in range(0,len(particles.id)):
    polygonMapping[iParticle] = keptIndices1[keptIndices2[iParticle]]



print("Extract polygon data for culled particles")
nVerticesOnCellCulled = np.zeros(len(particles.id),dtype="i")
xVertexCulled = np.zeros((len(particles.id), maxVertices))
yVertexCulled = np.zeros((len(particles.id), maxVertices))
for iParticle in range(0,len(particles.id)):
    nVerticesOnCellCulled[iParticle] = nVerticesOnCell[polygonMapping[iParticle]]
    for iVertex in range(0,maxVertices):
        xVertexCulled[iParticle,iVertex] = xVertex[polygonMapping[iParticle],iVertex]
        yVertexCulled[iParticle,iVertex] = yVertex[polygonMapping[iParticle],iVertex]


# find maximum cell size
print("Find max cell size")
maxCellSize = 0.0
for iParticle in range(0,len(particles.id)):
    for iVertex in range(1,nVerticesOnCellCulled[iParticle]):
        distance = math.sqrt(math.pow(xVertexCulled[iParticle,0]-xVertexCulled[iParticle,iVertex],2) +
                             math.pow(yVertexCulled[iParticle,0]-yVertexCulled[iParticle,iVertex],2))
        maxCellSize = max(maxCellSize,distance)
print("maxCellSize: ", maxCellSize)



# output
print("Write out particle file")
#particles.bonds = []
particles.create_bonds(min_overlap = -10000.0)
particles.to_particle_netcdf("particles_in_init.nc")

# add polygon to output file
add_init_polygons("particles_in_init.nc", maxCellSize, maxVertices, nVerticesOnCellCulled, xVertexCulled, yVertexCulled)

# add init connectivity info
add_init_connectivity("particles_in_init.nc")

# empty output
print("Write out no ice particle file")
for iParticle in range(0,particles.num_particles):
    if (particles.type[iParticle] == 1):
        particles.iceFraction [iParticle] = 0.0
        particles.iceThickness[iParticle] = 0.0
        particles.type[iParticle] = 0

#particles.bonds = []
particles.create_bonds(min_overlap = -10000.0)
particles.to_particle_netcdf("particles_in_init_no_ice.nc")

# add polygon to output file
add_init_polygons("particles_in_init_no_ice.nc", maxCellSize, maxVertices, nVerticesOnCellCulled, xVertexCulled, yVertexCulled)

# add init connectivity info
add_init_connectivity("particles_in_init_no_ice.nc")



print("Min/Max particle radius: ", np.amin(particles.radius), np.amax(particles.radius))



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
