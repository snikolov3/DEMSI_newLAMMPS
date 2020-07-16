from netCDF4 import Dataset
import math
import glob

#-------------------------------------------------------------------------------

def vector_mag(x1, y1, x2, y2):

    return math.sqrt(math.pow(x2-x1,2) + math.pow(y2-y1,2))

#-------------------------------------------------------------------------------

filenameExpression = "./output/*"

filenames = sorted(glob.glob(filenameExpression))

for filename in filenames:

    fileIn = Dataset(filename,"r")

    nBonds = len(fileIn.dimensions["nBonds"])

    bondEndPoint1Particle1 = fileIn.variables["bondEndPoint1Particle1"][:]
    bondEndPoint2Particle1 = fileIn.variables["bondEndPoint2Particle1"][:]
    bondEndPoint1Particle2 = fileIn.variables["bondEndPoint1Particle2"][:]
    bondEndPoint2Particle2 = fileIn.variables["bondEndPoint2Particle2"][:]
    bondCrackFraction = fileIn.variables["bondCrackFraction"][:]
    bondGlobalIDs = fileIn.variables["bondGlobalIDs"][:]

    fileIn.close()

    for iBond in range(0,nBonds):

        distance1 = vector_mag(\
            bondEndPoint1Particle1[iBond,0],bondEndPoint1Particle1[iBond,1],\
            bondEndPoint1Particle2[iBond,0],bondEndPoint1Particle2[iBond,1])

        distance2 = vector_mag(\
            bondEndPoint2Particle1[iBond,0],bondEndPoint2Particle1[iBond,1],\
            bondEndPoint2Particle2[iBond,0],bondEndPoint2Particle2[iBond,1])

        # check if bond exists
        if (bondCrackFraction[iBond,0] < bondCrackFraction[iBond,1]):

            print(iBond, distance1, distance2, True)

        else:

            print(iBond, distance1, distance2, False)
