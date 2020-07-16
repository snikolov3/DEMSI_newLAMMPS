import sys
from netCDF4 import Dataset
import math
import numpy as np
import argparse

#-------------------------------------------------------------------------------

def expand_particle_distribution(filenameIn, filenameOut, xDomainSize, yDomainSize):

    # read in distribution
    fileIn = Dataset(filenameIn,"r")

    periodic = fileIn.periodic
    if (not periodic):
        print("distribution must be periodic")
        sys.exit()

    xRangeIn = fileIn.xRange
    yRangeIn = fileIn.yRange

    nParticlesIn = len(fileIn.dimensions["nParticles"])

    xIn = fileIn.variables["x"][:]
    yIn = fileIn.variables["y"][:]
    radiusIn = fileIn.variables["radius"][:]

    fileIn.close()

    # decide how many replicas of domain to create
    nReplicateX = int(math.ceil(xDomainSize / xRangeIn))
    nReplicateY = int(math.ceil(yDomainSize / yRangeIn))
    nReplicate = nReplicateX * nReplicateY

    nParticlesOut = nParticlesIn * nReplicate
    xRangeOut = xRangeIn * nReplicateX
    yRangeOut = yRangeIn * nReplicateY

    # create output variables
    xOut = np.zeros(nParticlesOut)
    yOut = np.zeros(nParticlesOut)
    radiusOut = np.zeros(nParticlesOut)

    iParticleOut = 0

    for i in range(0,nReplicateX):
        for j in range(0,nReplicateY):

            for iParticleIn in range(0,nParticlesIn):

                xOut[iParticleOut] = xIn[iParticleIn] + i * xRangeIn
                yOut[iParticleOut] = yIn[iParticleIn] + j * yRangeIn
                radiusOut[iParticleOut] = radiusIn[iParticleIn]

                iParticleOut += 1

    # output
    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.periodic = "YES"
    fileOut.xRange = xRangeOut
    fileOut.yRange = yRangeOut

    fileOut.createDimension("nParticles",nParticlesOut)

    var = fileOut.createVariable("globalID","i",dimensions=["nParticles"])
    var[:] = np.arange(1, nParticlesOut+1, dtype="i")

    var = fileOut.createVariable("x","d",dimensions=["nParticles"])
    var[:] = xOut[:]

    var = fileOut.createVariable("y","d",dimensions=["nParticles"])
    var[:] = yOut[:]

    var = fileOut.createVariable("radius","d",dimensions=["nParticles"])
    var[:] = radiusOut[:]

    fileOut.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Replicate a periodic particle distribution to a desired size.')
    parser.add_argument('-i', dest='filenameIn', required=True, help='Input particle distribution.')
    parser.add_argument('-o', dest='filenameOut', required=True, help='Output particle distribution.')
    parser.add_argument('--dx', dest='xDomainSize', required=True, type=float, help='Minimum desired x domain size.')
    parser.add_argument('--dy', dest='yDomainSize', required=True, type=float, help='Minimum desired x domain size.')

    args = parser.parse_args()

    expand_particle_distribution(args.filenameIn, args.filenameOut, args.xDomainSize, args.yDomainSize)
