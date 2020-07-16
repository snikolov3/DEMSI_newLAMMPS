import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt

filenames = sorted(glob.glob("./output/particles_out*.nc"))

fileIn = Dataset(filenames[0],"r")

nParticles = len(fileIn.dimensions["nParticles"])

x0 = fileIn.variables["x"][:]

fileIn.close()

for iParticle in range(0,nParticles):

    print(iParticle)

    x = []

    for filename in filenames:

        #print("  ", filename)

        fileIn = Dataset(filename,"r")

        x.append(fileIn.variables["x"][iParticle]-x0[iParticle])

        fileIn.close()

    plt.plot(x)


plt.show()
