from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

#filename = "./output/particles_out.0001-01-01_00:00:00.nc"
filename = "./output/particles_out.0001-01-01_00:03:30.nc"

filein = Dataset(filename,"r")

nParticles = len(filein.dimensions["nParticles"])

nBonds = len(filein.dimensions["nBonds"])

bondCrackFraction = filein.variables["bondCrackFraction"][:,:]

b1p1 = filein.variables["bondEndPoint1Particle1"][:]
b2p1 = filein.variables["bondEndPoint2Particle1"][:]
b1p2 = filein.variables["bondEndPoint1Particle2"][:]
b2p2 = filein.variables["bondEndPoint2Particle2"][:]


bondFractions = []
segments = []
bondFractionsLines = []

for iBond in range(0,nBonds):

    bondFraction = bondCrackFraction[iBond,1] - bondCrackFraction[iBond,0]
    bondFractions.append(bondFraction)

    xEndPoint1 = 0.5 * (b1p1[iBond,0] + b1p2[iBond,0])
    yEndPoint1 = 0.5 * (b1p1[iBond,1] + b1p2[iBond,1])

    xEndPoint2 = 0.5 * (b2p1[iBond,0] + b2p2[iBond,0])
    yEndPoint2 = 0.5 * (b2p1[iBond,1] + b2p2[iBond,1])

    #if (bondFraction < 1.0):
    segments.append([(xEndPoint1,yEndPoint1),(xEndPoint2,yEndPoint2)])
    bondFractionsLines.append(bondFraction)
    

filein.close()

lineCollection = LineCollection(segments,cmap=plt.get_cmap("viridis"))
lineCollection.set_array(np.array(bondFractionsLines))

# plot
fig, axis = plt.subplots()
axis.add_collection(lineCollection)
axis.set_xlim((0,500000))
axis.set_ylim((0,500000))
plt.savefig("bonds.png")
