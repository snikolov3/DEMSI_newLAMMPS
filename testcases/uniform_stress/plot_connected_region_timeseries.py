from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt
import numpy as np
import scipy, sys
import scipy.sparse

filenames = sorted(glob.glob("./output/*"))

nComponents = []

for filename in filenames:

    print(filename)

    filein = Dataset(filename,"r")

    nParticles = len(filein.dimensions["nParticles"])

    globalID = filein.variables["globalID"][:]

    globalIDToIndex = {}
    for iParticle in range(0,nParticles):
        globalIDToIndex[globalID[iParticle]] = iParticle

    try:
        
        nBonds = len(filein.dimensions["nBonds"])

        bondGlobalIDs = filein.variables["bondGlobalIDs"][:,:]
        bondCrackFraction = filein.variables["bondCrackFraction"][:,:]

        bondFraction = bondCrackFraction[:,1] - bondCrackFraction[:,0]

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
        nComponents.append(scipy.sparse.csgraph.connected_components(mtx, directed=False, return_labels=False))

    except:
        
        nComponents.append(nParticles)


# plot
fig, axis = plt.subplots()
axis.semilogy(nComponents)
axis.set_xlabel("Time")
axis.set_ylabel("Number connected regions")
plt.savefig("connected_region_number_timeseries.png")
