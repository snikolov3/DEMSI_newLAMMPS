from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
import numpy as np
import scipy, sys
import scipy.sparse, scipy.spatial
import random


#filename = "./output/particles_out.0001-01-01_00:00:00.nc"
filename = "./output/particles_out.0001-01-01_00:02:10.nc"

# load data
filein = Dataset(filename,"r")

nParticles = len(filein.dimensions["nParticles"])

x = filein.variables["x"][:]
y = filein.variables["y"][:]
r = filein.variables["radius"][:]
globalID = filein.variables["globalID"][:]

nBonds = len(filein.dimensions["nBonds"])

bondCrackFraction = filein.variables["bondCrackFraction"][:,:]
bondGlobalIDs = filein.variables["bondGlobalIDs"][:,:]

b1p1 = filein.variables["bondEndPoint1Particle1"][:]
b2p1 = filein.variables["bondEndPoint2Particle1"][:]
b1p2 = filein.variables["bondEndPoint1Particle2"][:]
b2p2 = filein.variables["bondEndPoint2Particle2"][:]

filein.close()


# derived quantities
bondFraction = np.zeros(nBonds)
for iBond in range(0,nBonds):
    bondFraction[iBond] = bondCrackFraction[iBond,1] - bondCrackFraction[iBond,0]

globalIDToIndex = {}
for iParticle in range(0,nParticles):
    globalIDToIndex[globalID[iParticle]] = iParticle

    
# determine connected regions
row = []
col = []
data = []

for iBond in range(0,nBonds):
    if (bondFraction[iBond] > 0.0):
        row.append(globalIDToIndex[bondGlobalIDs[iBond,0]])
        col.append(globalIDToIndex[bondGlobalIDs[iBond,1]])
        data.append(1.0)
        row.append(globalIDToIndex[bondGlobalIDs[iBond,1]])
        col.append(globalIDToIndex[bondGlobalIDs[iBond,0]])
        data.append(1.0)

mtx = scipy.sparse.csr_matrix((data, (row, col)), shape=(nParticles, nParticles))
n_components, labels = scipy.sparse.csgraph.connected_components(mtx, directed=False, return_labels=True)


# randomize the regions
randomizer = np.zeros(n_components)
for i in range(0,n_components):
    randomizer[i] = i
random.shuffle(randomizer)

for i in range(0,nParticles):
    labels[i] = randomizer[labels[i]]



# bonds
segments = []
bondFractionsLines = []

for iBond in range(0,nBonds):

    bondFraction = bondCrackFraction[iBond,1] - bondCrackFraction[iBond,0]

    xEndPoint1 = 0.5 * (b1p1[iBond,0] + b1p2[iBond,0])
    yEndPoint1 = 0.5 * (b1p1[iBond,1] + b1p2[iBond,1])

    xEndPoint2 = 0.5 * (b2p1[iBond,0] + b2p2[iBond,0])
    yEndPoint2 = 0.5 * (b2p1[iBond,1] + b2p2[iBond,1])

    segments.append([(xEndPoint1,yEndPoint1),(xEndPoint2,yEndPoint2)])
    bondFractionsLines.append(bondFraction)

lc = LineCollection(segments,cmap=plt.get_cmap("viridis"))
lc.set_array(np.array(bondFractionsLines))



# connected regions
patches = []
colors = []

for iParticle in range(0, nParticles):
    patches.append(Circle((x[iParticle], y[iParticle]), r[iParticle]))
    colors.append(labels[iParticle])

pc = PatchCollection(patches,cmap=plt.get_cmap("jet"))
pc.set_array(np.array(colors))




# plot connected regions
fig, axis = plt.subplots()

axis.add_collection(pc)
#axis.add_collection(lc)

axis.set_xlim((50000,450000))
axis.set_ylim((50000,450000))
plt.savefig("regions.png")


points = np.zeros((nParticles,2))
for iPoint in range(0,nParticles):
    points[iPoint,0] = x[iPoint]
    points[iPoint,1] = y[iPoint]
vor = scipy.spatial.Voronoi(points)


segments = []

nRidges = vor.ridge_points.shape[0]
for iRidge in range(0,nRidges):
    iPoint1 = vor.ridge_points[iRidge,0]
    iPoint2 = vor.ridge_points[iRidge,1]
    if (labels[iPoint1] != labels[iPoint2]):
       iVertex1 = vor.ridge_vertices[iRidge][0]
       iVertex2 = vor.ridge_vertices[iRidge][1]
       if (iVertex1 > -1 and iVertex2 > -1):
           x1 = vor.vertices[iVertex1,0]
           y1 = vor.vertices[iVertex1,1]
           x2 = vor.vertices[iVertex2,0]
           y2 = vor.vertices[iVertex2,1]
           segments.append([[x1,y1],[x2,y2]])
       
lc = LineCollection(segments)

# plot connected regions
fig, axis = plt.subplots()

axis.add_collection(lc)

axis.set_xlim((50000,450000))
axis.set_ylim((50000,450000))
plt.savefig("region_boundaries.png")
