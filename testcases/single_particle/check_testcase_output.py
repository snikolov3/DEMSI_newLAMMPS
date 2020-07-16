from netCDF4 import Dataset
import glob, sys, math

outputDirectory = sys.argv[1]

xExpected = 10010.0
yExpected = 500000.0
epsilon = 1e-8

filenames = sorted(glob.glob("%s/*" %(outputDirectory)))

fileOut = Dataset(filenames[-1],"r")

x = fileOut.variables["x"][:,0]
y = fileOut.variables["x"][:,1]

returnCode = 0

# check x position
if (math.fabs(x[0] - xExpected) < epsilon):
    print("x position is correct at %f" %(xExpected))
else:
    print("x position of particle incorrect at %f instead of %f" %(x[0],xExpected))
    returnCode = -1


# check y position
if (y[0] == yExpected):
    print("y position is correct at %f" %(yExpected))
else:
    print("y position of particle incorrect at %f instead of %f" %(y[0],yExpected))
    returnCode = -1

fileOut.close()

exit(returnCode)
